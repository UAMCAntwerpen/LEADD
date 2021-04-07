#include "Fragmentation.hpp"

void AssignImmutableAtomIndices(RDKit::ROMol& mol) {
  unsigned immutable_idx = 0u;
  for (RDKit::Atom* atom : mol.atoms()) {
    atom->setProp<unsigned>("ImmutableIdx", immutable_idx);
    ++immutable_idx;
  };
};

void AssignMMFFAtomTypes(RDKit::ROMol& mol) {
  RDKit::MMFF::MMFFMolProperties mmff_properties = RDKit::MMFF::MMFFMolProperties(mol);
  for (RDKit::Atom* atom : mol.atoms()) {
    unsigned atom_idx = atom->getIdx();
    unsigned mmff_atom_type = mmff_properties.getMMFFAtomType(atom_idx);
    atom->setProp<unsigned>("MMFF_atom_type", mmff_atom_type);
  };
};

void FlagRingAtoms(RDKit::ROMol& mol) {
  const RDKit::RingInfo* ring_info = mol.getRingInfo();
  for (RDKit::Atom* atom : mol.atoms()) {
    if (ring_info->numAtomRings(atom->getIdx()) > 0) {
      atom->setProp<bool>("WasInRing", true);
    } else {
      atom->setProp<bool>("WasInRing", false);
    };
  };
};

void FlagAtoms(RDKit::ROMol& mol) {
  unsigned immutable_idx = 0u;
  RDKit::MMFF::MMFFMolProperties mmff_properties = RDKit::MMFF::MMFFMolProperties(mol);
  const RDKit::RingInfo* ring_info = mol.getRingInfo();
  for (RDKit::Atom* atom : mol.atoms()) {
    unsigned atom_idx = atom->getIdx();
    unsigned mmff_atom_type = mmff_properties.getMMFFAtomType(atom_idx);
    atom->setProp<unsigned>("ImmutableIdx", immutable_idx);
    atom->setProp<unsigned>("MMFF_atom_type", mmff_atom_type);
    if (ring_info->numAtomRings(atom_idx) > 0) {
      atom->setProp<bool>("WasInRing", true);
    } else {
      atom->setProp<bool>("WasInRing", false);
    };
    ++immutable_idx;
  };
};

void CopyFlags(const RDKit::Atom* source, RDKit::Atom* destination) {
  destination->setProp<unsigned>("ImmutableIdx", source->getProp<unsigned>("ImmutableIdx"));
  destination->setProp<unsigned>("MMFF_atom_type", source->getProp<unsigned>("MMFF_atom_type"));
  destination->setProp<bool>("WasInRing", source->getProp<bool>("WasInRing"));
};

void ExpandRingSystem(std::vector<int>& ring_system, std::list<std::vector<int>>& rings) {
  bool expanded_ring_system = false;
  std::list<std::vector<int>>::iterator it = rings.begin();
  while (it != rings.end()) {
    bool erased_ring = false;
    for (int atom_idx : *it) {
      if (std::find(ring_system.begin(), ring_system.end(), atom_idx) != ring_system.end()) {
        ring_system.insert(ring_system.end(), it->begin(), it->end());
        it = rings.erase(it);
        erased_ring = true;
        expanded_ring_system = true;
        break;
      };
    };
    if (!erased_ring) {
      ++it;
    };
  };
  if (expanded_ring_system) {
    ExpandRingSystem(ring_system, rings);
  };
};

std::vector<std::vector<int>> DefineRingSystems(const RDKit::ROMol& mol) {
  std::vector<std::vector<int>> ring_systems;
  RDKit::RingInfo* ring_info = mol.getRingInfo();
  std::vector<std::vector<int>> sssr(ring_info->atomRings());
  std::list<std::vector<int>> remaining_rings(sssr.begin(), sssr.end());
  for (const auto& ring : sssr) {
    if (std::find(remaining_rings.begin(), remaining_rings.end(), ring) != remaining_rings.end()) {
      std::vector<int> ring_system = ring;
      ExpandRingSystem(ring_system, remaining_rings);
      std::sort(ring_system.begin(), ring_system.end());
      ring_system.erase(std::unique(ring_system.begin(), ring_system.end()), ring_system.end());
      ring_systems.push_back(ring_system);
    };
  };
  return ring_systems;
};

// Convert a dummy molecule (i.e. a RDKit::ROMol with Connection-encoding pseudoatoms) into a Pseudofragment.
Pseudofragment DummyMolToPseudofragment(const RDKit::ROMol& mol, boost::format& formatter) {
  // Initialize the pseudomol of the Pseudofragment.
  RDKit::RWMol pseudomol = RDKit::RWMol(mol);

  // Initialize the Pseudofragment's connection table.
  ConnectionsTable connections;

  // Initialize a RingInfo object.
  bool has_ring = false;
  RDKit::RingInfo* ring_info = mol.getRingInfo();

  // Initialize a vector to store the indices of the molecule's pseudoatoms.
  std::vector<unsigned> pseudoatom_indices;

  // Loop over the molecule's atoms.
  for (RDKit::Atom* atom : pseudomol.atoms()) {
    // For all encountered pseudoatoms (i.e. atoms with an atomic number of 0),
    // store their index and add the Connection they encode to the connections table.
    if (atom->getAtomicNum() == 0) {
      std::uint8_t atom_type, neighbor_type;
      std::tie(neighbor_type, atom_type) = DecodeConnection(atom->getIsotope());
      unsigned atom_idx = atom->getIdx();
      pseudoatom_indices.push_back(atom_idx);
      // The pseudoatom has only a single neighbor so no true iteration is necessary.
      RDKit::ROMol::ADJ_ITER nbr_it, nbr_end_it;
      std::tie(nbr_it, nbr_end_it) = mol.getAtomNeighbors(atom);
      unsigned nbr_idx = *nbr_it;
      RDKit::Bond::BondType bond_type = mol.getBondBetweenAtoms(atom_idx, nbr_idx)->getBondType();
      unsigned bond_type_int = bond_type_int_table.at(bond_type);
      Connection connection(neighbor_type, atom_type, bond_type_int, formatter);
      connections.AddConnection(connection, nbr_idx);
      // For all regular atoms, label them according to whether they are part of a
      // ring system or not.
    } else {
      if (!has_ring && ring_info->numAtomRings(atom->getIdx()) > 0) {
        has_ring = true;
      };
    };
  };

  // Delete the pseudoatoms in sorted descending order (to avoid triggering atom
  // reindexing).
  std::sort(pseudoatom_indices.begin(), pseudoatom_indices.end(), std::greater<unsigned>());
  for (unsigned idx : pseudoatom_indices) {
    pseudomol.removeAtom(idx);
  };

  // Encode the dummy molecule as a SMILES string and generate the Pseudofragment.
  std::string smiles = RDKit::MolToSmiles(mol);
  Pseudofragment pseudofragment(pseudomol, connections, smiles, has_ring, false);
  return pseudofragment;
};

// Function to convert an RDKit::Atom into a Pseudofragment.
Pseudofragment AtomToPseudofragment(const RDKit::ROMol& mol, const RDKit::Atom* atom, boost::format& formatter) {
  // Initialize the Pseudofragment's pseudomol and connections table.
  RDKit::RWMol pseudomol;
  ConnectionsTable connections;

  // Initialize the dummy molecule used to create the Connection-encoding SMILES.
  RDKit::RWMol dummymol;

  // Add the base Atom to the pseudomol and dummy molecule. Apparently its
  // important to create a copy of the atom before adding it, otherwise it
  // adds an empty pseudoatom.
  RDKit::Atom atom_copy = *atom;
  unsigned pseudomol_atom_idx = pseudomol.addAtom(&atom_copy);
  unsigned dummymol_atom_idx = dummymol.addAtom(&atom_copy);

  // Loop over the neighboring Atoms of the base Atom in the owning Mol.
  RDKit::ROMol::ADJ_ITER nbr_it, nbr_end_it;
  std::tie(nbr_it, nbr_end_it) = mol.getAtomNeighbors(atom);
  while (nbr_it != nbr_end_it) {
    unsigned nbr_idx = *nbr_it;
    const RDKit::Atom* neighbor = mol.getAtomWithIdx(nbr_idx);

    // If the neighboring atom is a pseudoatom, decode its isotope into atom types.
    // If not, get the atom types from the Atom properties.
    std::uint8_t atom_type, neighbor_type;
    if (neighbor->getAtomicNum() == 0) {
      std::tie(atom_type, neighbor_type) = DecodeConnection(neighbor->getIsotope());
    }
    else {
      atom_type = atom_copy.getProp<unsigned>("MMFF_atom_type");
      neighbor_type = neighbor->getProp<unsigned>("MMFF_atom_type");
    };

    // Create a Connection object symbolizing the bond between the Atom and the Neighbor.
    RDKit::Bond::BondType bond_type = mol.getBondBetweenAtoms(atom->getIdx(), nbr_idx)->getBondType();
    unsigned bond_type_int = bond_type_int_table.at(bond_type);
    Connection connection(atom_type, neighbor_type, bond_type_int, formatter);

    // Add the Connection to the connections table.
    connections.AddConnection(connection, pseudomol_atom_idx);

    // Create a pseudoatom with an isotope encoding the start- and end atom types
    // of the Connection point and add it to the dummy molecule.
    RDKit::Atom pseudoatom(0);
    pseudoatom.setIsotope(connection.GetEncoded());
    unsigned pseudoatom_idx = dummymol.addAtom(&pseudoatom);
    dummymol.addBond(dummymol_atom_idx, pseudoatom_idx, bond_type);

    ++nbr_it;
  };

  // Encode the dummy molecule as a SMILES string and generate the Pseudofragment.
  std::string smiles = RDKit::MolToSmiles(dummymol);
  Pseudofragment pseudofragment(pseudomol, connections, smiles, false, false);
  return pseudofragment;
};

std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> SplitAtRings(const RDKit::ROMol& mol) {
  // Initialize containers to store the data needed to split the molecule informaively.
  std::vector<unsigned> split_bond_indices;
  std::vector<RDKit::Bond::BondType> split_bond_types;
  std::vector<std::pair<unsigned, unsigned>> split_bond_labels;

  // Define the ring systems of the molecule.
  std::vector<std::vector<int>> ring_systems = DefineRingSystems(mol);
  // Find the bonds that connect ring systems to the acyclic parts of the molecule.
  for (const auto& ring_system : ring_systems) {
    for (int atom_idx : ring_system) {
      const RDKit::Atom* atom = mol.getAtomWithIdx(atom_idx);
      RDKit::ROMol::ADJ_ITER neighbors_it, neighbors_end_it;
      std::tie(neighbors_it, neighbors_end_it) = mol.getAtomNeighbors(atom);
      while (neighbors_it != neighbors_end_it) {
        int neighbor_idx = *neighbors_it;
        if (std::find(ring_system.begin(), ring_system.end(), neighbor_idx) == ring_system.end()) {
          const RDKit::Bond* bond = mol.getBondBetweenAtoms(atom_idx, neighbor_idx);
          unsigned bond_idx = bond->getIdx();
          // Flag the identified bond for cleavage only if it hasn't been
          // flagged before. This could happen if two rings are connected
          // directly through a bond.
          if (std::find(split_bond_indices.begin(), split_bond_indices.end(), bond_idx) == split_bond_indices.end()) {
            RDKit::Bond::BondType bond_type = bond->getBondType();
            split_bond_indices.push_back(bond_idx);
            split_bond_types.push_back(bond_type);
            // Use as post-fragmentation labels Connection encoding integers.
            Connection connection(bond->getBeginAtom()->getProp<unsigned>("MMFF_atom_type"), bond->getEndAtom()->getProp<unsigned>("MMFF_atom_type"), bond_type_int_table[bond_type]);
            Connection mirrored_connection = connection.Mirror();
            std::pair<unsigned, unsigned> p(mirrored_connection.GetEncoded(), connection.GetEncoded());
            split_bond_labels.push_back(p);
          };
        };
        ++neighbors_it;
      };
    };
  };

  // Fragment the molecule by cutting the specified bonds.
  RDKit::ROMol* fragmented_mol = RDKit::MolFragmenter::fragmentOnBonds(mol, split_bond_indices, true, &split_bond_labels, &split_bond_types);
  // Extract the resulting fragments as separate RDKit::ROMols.
  std::vector<int> fragment_map;
  std::vector<std::vector<int>> fragment_atom_map;
  std::vector<RDKit::ROMOL_SPTR> fragments = RDKit::MolOps::getMolFrags(*fragmented_mol, true, &fragment_map, &fragment_atom_map);
  delete fragmented_mol;

  // Classify the fragments into acyclic and cyclic fragments, and copy the
  // RDKit::Atom properties from the original molecule to the corresponding fragments.
  // This is necessary since RDKit::MolOps::getMolFrags creates new graphs.
  std::vector<RDKit::ROMOL_SPTR> acyclic_fragments, ring_fragments;
  for (unsigned fragment_idx = 0; fragment_idx < fragments.size(); ++fragment_idx) {
    bool has_ring = false;
    RDKit::ROMOL_SPTR fragment = fragments[fragment_idx];
    const std::vector<int>& atom_map = fragment_atom_map[fragment_idx];
    for (unsigned mapping_idx = 0; mapping_idx < atom_map.size(); ++mapping_idx) {
      unsigned new_atom_idx = mapping_idx;
      RDKit::Atom* new_atom = fragment->getAtomWithIdx(new_atom_idx);
      if (new_atom->getAtomicNum() != 0) {
        unsigned original_atom_idx = atom_map[mapping_idx];
        const RDKit::Atom* original_atom = mol.getAtomWithIdx(original_atom_idx);
        CopyFlags(original_atom, new_atom);
        if (!has_ring && new_atom->getProp<bool>("WasInRing")) {
          has_ring = true;
        };
      };
    };
    if (has_ring) {
      ring_fragments.push_back(fragment);
    } else {
      acyclic_fragments.push_back(fragment);
    };
  };

  // Return the resulting fragments.
  std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> output(std::move(acyclic_fragments), std::move(ring_fragments));
  return output;
};
