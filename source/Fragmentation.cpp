#include "Fragmentation.hpp"

void ExpandRingSystem(std::set<unsigned>& ring_system, std::vector<std::set<unsigned>>& rings) {
  // Iterate over the rings that may be annexed by the ring system.
  bool expanded = false;
  std::vector<std::set<unsigned>>::iterator ring_it = rings.begin();
  while (ring_it != rings.end()) {
    // Check if the ring overlaps with the ring system.
    const std::set<unsigned>& ring = *ring_it;
    std::vector<unsigned> shared_atoms;
    std::set_intersection(ring_system.begin(), ring_system.end(), ring.begin(), ring.end(), std::back_inserter(shared_atoms));
    // If it does, annex it.
    if (!shared_atoms.empty()) {
      ring_system.insert(ring.begin(), ring.end());
      ring_it = rings.erase(ring_it);
      expanded = true;
    } else {
      ++ring_it;
    };
  };
  // Recursively annex rings.
  if (expanded) {
    ExpandRingSystem(ring_system, rings);
  };
};

std::vector<std::set<unsigned>> DefineRingSystems(const RDKit::ROMol& molecule) {
  std::vector<std::set<unsigned>> ring_systems;
  // Find the molecule's SSSR.
  RDKit::RingInfo* ring_info = molecule.getRingInfo();
  std::vector<std::vector<int>> tmp_sssr = ring_info->atomRings();
  std::vector<std::set<unsigned>> sssr;
  sssr.reserve(tmp_sssr.size());
  for (const std::vector<int>& ring : tmp_sssr) {
    sssr.emplace_back(ring.begin(), ring.end());
  };
  // Iterate over the rings in the SSSR.
  std::vector<std::set<unsigned>> remaining_rings = sssr;
  for (const std::set<unsigned>& ring : sssr) {
    // Check if the ring was already assigned to a ring system.
    std::vector<std::set<unsigned>>::const_iterator ring_it = std::find(remaining_rings.begin(), remaining_rings.end(), ring);
    // If it wasn't, expand it into a ring system.
    if (ring_it != remaining_rings.end()) {
      std::set<unsigned> ring_system = ring;
      remaining_rings.erase(ring_it);
      ExpandRingSystem(ring_system, remaining_rings);
      ring_systems.push_back(ring_system);
    };
  };
  return ring_systems;
};

std::vector<std::uint32_t> GetDummyAtomTypes(const RDKit::ROMol& molecule) {
  std::vector<std::uint32_t> atom_types (molecule.getNumAtoms());
  std::fill(atom_types.begin(), atom_types.end(), 0);
  return atom_types;
};

std::vector<std::uint32_t> GetAtomicNumberAtomTypes(const RDKit::ROMol& molecule) {
  std::vector<std::uint32_t> atom_types (molecule.getNumAtoms());
  for (size_t atom_idx = 0; atom_idx < molecule.getNumAtoms(); ++atom_idx) {
    const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
    atom_types[atom_idx] = atom->getAtomicNum();
  };
  return atom_types;
};

std::vector<std::uint32_t> GetMMFFAtomTypes(const RDKit::ROMol& molecule) {
  // For some reason the MMFFMolProperties constructor modifies the input RDKit::ROMol.
  // To avoid different outputs we make a copy.
  RDKit::ROMol tmp_molecule = molecule;
  RDKit::MMFF::MMFFMolProperties mmff_properties = RDKit::MMFF::MMFFMolProperties(tmp_molecule);
  std::vector<std::uint32_t> atom_types (tmp_molecule.getNumAtoms());
  for (std::uint32_t atom_idx = 0; atom_idx < tmp_molecule.getNumAtoms(); ++atom_idx) {
    atom_types[atom_idx] = mmff_properties.getMMFFAtomType(atom_idx);
  };
  return atom_types;
};

std::vector<std::uint32_t> GetMorganAtomTypes(const RDKit::ROMol& molecule, unsigned radius, bool consider_chirality) {
  RDKit::MorganFingerprints::BitInfoMap feature_atom_map;
  RDKit::SparseIntVect<std::uint32_t>* fingerprint = RDKit::MorganFingerprints::getFingerprint(molecule, radius, nullptr, nullptr, consider_chirality, true, false, false, &feature_atom_map, true);
  delete fingerprint;
  std::vector<std::uint32_t> atom_types (molecule.getNumAtoms());
  for (std::uint32_t atom_idx = 0; atom_idx < molecule.getNumAtoms(); ++atom_idx) {
    bool done = false;
    for (const auto& [feature, atom_radius_pairs] : feature_atom_map) {
      for (const auto& [atom_id, feature_radius] : atom_radius_pairs) {
        if (atom_id == atom_idx && radius == feature_radius) {
          atom_types[atom_idx] = feature;
          done = true;
          break;
        };
      };
      if (done) {
        break;
      };
    };
  };
  return atom_types;
};

std::vector<std::uint32_t> GetHashedMorganAtomTypes(const RDKit::ROMol& molecule, unsigned radius, unsigned n_bits, bool consider_chirality) {
  RDKit::MorganFingerprints::BitInfoMap feature_atom_map;
  RDKit::SparseIntVect<std::uint32_t>* fingerprint = RDKit::MorganFingerprints::getHashedFingerprint(molecule, radius, n_bits, nullptr, nullptr, consider_chirality, true, false, &feature_atom_map, true);
  delete fingerprint;
  std::vector<std::uint32_t> atom_types (molecule.getNumAtoms());
  for (std::uint32_t atom_idx = 0; atom_idx < molecule.getNumAtoms(); ++atom_idx) {
    bool done = false;
    for (const auto& [feature, atom_radius_pairs] : feature_atom_map) {
      for (const auto& [atom_id, feature_radius] : atom_radius_pairs) {
        if (atom_id == atom_idx && radius == feature_radius) {
          atom_types[atom_idx] = feature;
          done = true;
          break;
        };
      };
      if (done) {
        break;
      };
    };
  };
  return atom_types;
};

std::vector<std::uint32_t> GetAtomTypes(const RDKit::ROMol& molecule, const FragmentationSettings& settings) {
  FragmentationSettings::AtomTyping atom_typing = settings.GetAtomTyping();
  if (atom_typing == FragmentationSettings::AtomTyping::DUMMY) {
    return GetDummyAtomTypes(molecule);
  };
  if (atom_typing == FragmentationSettings::AtomTyping::ATOMIC_NUMBER) {
    return GetAtomicNumberAtomTypes(molecule);
  };
  if (atom_typing == FragmentationSettings::AtomTyping::MMFF) {
    return GetMMFFAtomTypes(molecule);
  };
  if (atom_typing == FragmentationSettings::AtomTyping::MORGAN) {
    return GetMorganAtomTypes(molecule, settings.GetMorganRadius(), settings.MorganConsiderChirality());
  };
  if (atom_typing == FragmentationSettings::AtomTyping::HASHED_MORGAN) {
    return GetHashedMorganAtomTypes(molecule, settings.GetMorganRadius(), settings.GetHashedMorganNBits(), settings.MorganConsiderChirality());
  };
};

void AssignOldIdxAtomProperty(RDKit::ROMol& molecule) {
  for (RDKit::Atom* atom : molecule.atoms()) {
    atom->setProp<unsigned>("OldIdx", atom->getIdx());
  };
};

void RemoveBonds(RDKit::RWMol& molecule, const std::vector<unsigned>& bond_indices) {
  std::vector<RDKit::Bond*> bonds (bond_indices.size());
  for (size_t i = 0; i < bond_indices.size(); ++i) {
    bonds[i] = molecule.getBondWithIdx(bond_indices[i]);
  };
  for (RDKit::Bond* bond : bonds) {
    molecule.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
  };
};

std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> SplitAtRings(const RDKit::ROMol& molecule, bool copy_old_idx_property) {
  // Define the ring systems of the molecule.
  std::vector<std::set<unsigned>> ring_systems = DefineRingSystems(molecule);
  // Find the bonds that delimit ring systems.
  std::vector<unsigned> ring_delimiting_bond_indices;
  // To do so, iterate over the ring systems.
  for (const std::set<unsigned>& ring_system : ring_systems) {
    // Iterate over the atoms in the ring system.
    for (unsigned atom_idx : ring_system) {
      const RDKit::Atom* atom = molecule.getAtomWithIdx(atom_idx);
      // Iterate over the neighboring atoms.
      RDKit::ROMol::ADJ_ITER neighbors_it, neighbors_end_it;
      std::tie(neighbors_it, neighbors_end_it) = molecule.getAtomNeighbors(atom);
      while (neighbors_it != neighbors_end_it) {
        unsigned neighbor_idx = *neighbors_it;
        // If the neighboring atom doesn't pertain to the ring system the bond
        // delimits a ring.
        if (ring_system.find(neighbor_idx) == ring_system.end()) {
          const RDKit::Bond* bond = molecule.getBondBetweenAtoms(atom_idx, neighbor_idx);
          ring_delimiting_bond_indices.push_back(bond->getIdx());
        };
        ++neighbors_it;
      };
    };
  };
  // Remove duplicate bond indices.
  std::sort(ring_delimiting_bond_indices.begin(), ring_delimiting_bond_indices.end());
  ring_delimiting_bond_indices.erase(std::unique(ring_delimiting_bond_indices.begin(), ring_delimiting_bond_indices.end()), ring_delimiting_bond_indices.end());
  // Fragment the molecule by splitting the aforementioned bonds. This shouldn't
  // alter atom indexing.
  RDKit::RWMol fragmented_molecule = molecule;
  RemoveBonds(fragmented_molecule, ring_delimiting_bond_indices);
  // If necessary, assign the OldIdx atom property. If we don't do this it's
  // assumed to have been present in the input molecule.
  if (!copy_old_idx_property) {
    AssignOldIdxAtomProperty(fragmented_molecule);
  };
  // Convert the fragments into separate molecules. Atom properties shouldn't be
  // lost in this step.
  std::vector<RDKit::ROMOL_SPTR> fragments = RDKit::MolOps::getMolFrags(fragmented_molecule, false);
  // Determine which of the resulting fragments are cyclic.
  std::vector<RDKit::ROMOL_SPTR> cyclic_fragments, acyclic_fragments;
  for (const RDKit::ROMOL_SPTR fragment : fragments) {
    bool is_cycle = false;
    for (const RDKit::Atom* atom : fragment->atoms()) {
      unsigned atom_idx = atom->getProp<unsigned>("OldIdx");
      for (const std::set<unsigned>& ring_system : ring_systems) {
        if (ring_system.find(atom_idx) != ring_system.end()) {
          is_cycle = true;
          break;
        };
      };
      if (is_cycle) {
        break;
      };
    };
    if (is_cycle) {
      cyclic_fragments.push_back(fragment);
    } else {
      acyclic_fragments.push_back(fragment);
    };
  };
  return {cyclic_fragments, acyclic_fragments};
};

RDKit::ROMOL_SPTR PathToFragment(const RDKit::ROMol& molecule, const std::vector<int>& path, bool copy_old_idx_property) {
  std::map<int, int> atom_mapping;
  RDKit::ROMol* fragment = RDKit::Subgraphs::pathToSubmol(molecule, path, false, atom_mapping);
  for (auto[old_atom_idx, new_atom_idx] : atom_mapping) {
    const RDKit::Atom* old_atom = molecule.getAtomWithIdx(old_atom_idx);
    RDKit::Atom* new_atom = fragment->getAtomWithIdx(new_atom_idx);
    // The call to pathToSubmol erases atom properties, which is why we need to
    // either copy over the OldIdx or assign a new one.
    if (copy_old_idx_property) {
      new_atom->setProp<unsigned>("OldIdx", old_atom->getProp<unsigned>("OldIdx"));
    } else {
      new_atom->setProp<unsigned>("OldIdx", new_atom_idx);
    };
  };
  return RDKit::ROMOL_SPTR(fragment);
};

std::vector<RDKit::ROMOL_SPTR> AtomsToFragments(RDKit::ROMol& molecule, bool copy_old_idx_property) {
  std::vector<RDKit::ROMOL_SPTR> fragments;
  for (RDKit::Atom* atom : molecule.atoms()) {
    RDKit::RWMol* fragment = new RDKit::RWMol;
    unsigned fragment_atom_idx = fragment->addAtom(atom);
    if (!copy_old_idx_property) {
      RDKit::Atom* fragment_atom = fragment->getAtomWithIdx(fragment_atom_idx);
      fragment_atom->setProp<unsigned>("OldIdx", atom->getIdx());
    };
    fragments.push_back(RDKit::ROMOL_SPTR(fragment));
  };
  return fragments;
};

std::vector<RDKit::ROMOL_SPTR> NoFragments(RDKit::ROMol& molecule) {
  std::vector<RDKit::ROMOL_SPTR> fragments;
  fragments.emplace_back(new RDKit::ROMol(molecule));
  return fragments;
};

std::vector<RDKit::ROMOL_SPTR> LinearPathsToFragments(RDKit::ROMol& molecule, unsigned minimum_n_bonds, unsigned maximum_n_bonds, bool copy_old_idx_property) {
  std::vector<RDKit::ROMOL_SPTR> fragments;
  if (minimum_n_bonds == 0) {
    fragments = AtomsToFragments(molecule, copy_old_idx_property);
    minimum_n_bonds = 1;
  };
  if (maximum_n_bonds > 0){
    std::map<int, std::list<std::vector<int>>> paths_by_length = RDKit::findAllPathsOfLengthsMtoN(molecule, minimum_n_bonds, maximum_n_bonds);
    for (const auto& [length, paths] : paths_by_length) {
      for (const std::vector<int>& path : paths) {
        RDKit::ROMOL_SPTR fragment = PathToFragment(molecule, path, copy_old_idx_property);
        fragments.push_back(fragment);
      };
    };
  };
  return fragments;
};

std::vector<RDKit::ROMOL_SPTR> SubgraphsToFragments(RDKit::ROMol& molecule, unsigned minimum_n_bonds, unsigned maximum_n_bonds, bool copy_old_idx_property) {
  std::vector<RDKit::ROMOL_SPTR> fragments;
  if (minimum_n_bonds == 0) {
    fragments = AtomsToFragments(molecule, copy_old_idx_property);
    minimum_n_bonds = 1;
  };
  if (maximum_n_bonds > 0) {
    std::map<int, std::list<std::vector<int>>> subgraphs_by_size = RDKit::findAllSubgraphsOfLengthsMtoN(molecule, minimum_n_bonds, maximum_n_bonds);
    for (const auto& [size, paths] : subgraphs_by_size) {
      for (const std::vector<int>& path : paths) {
        RDKit::ROMOL_SPTR fragment = PathToFragment(molecule, path, copy_old_idx_property);
        fragments.push_back(fragment);
      };
    };
  };
  return fragments;
};

std::vector<RDKit::ROMOL_SPTR> CircularEnvironmentsToFragments(RDKit::ROMol& molecule, unsigned minimum_radius, unsigned maximum_radius, bool copy_old_idx_property) {
  std::vector<RDKit::ROMOL_SPTR> fragments;
  if (minimum_radius == 0) {
    fragments = AtomsToFragments(molecule, copy_old_idx_property);
    minimum_radius = 1;
  };
  if (maximum_radius > 0) {
    for (unsigned radius = minimum_radius; radius <= maximum_radius; ++radius) {
      for (unsigned atom_idx = 0; atom_idx < molecule.getNumAtoms(); ++atom_idx) {
        std::vector<int> path = RDKit::findAtomEnvironmentOfRadiusN(molecule, radius, atom_idx);
        RDKit::ROMOL_SPTR fragment = PathToFragment(molecule, path, copy_old_idx_property);
        fragments.push_back(fragment);
      };
    };
  };
  return fragments;
};

std::vector<RDKit::ROMOL_SPTR> FragmentSystematically(RDKit::ROMol& molecule, const FragmentationSettings& settings, bool copy_old_idx_property) {
  FragmentationSettings::SystematicFragmentation systematic_fragmentation = settings.GetSystematicFragmentation();
  if (systematic_fragmentation == FragmentationSettings::SystematicFragmentation::NONE) {
    return NoFragments(molecule);
  };
  if (systematic_fragmentation == FragmentationSettings::SystematicFragmentation::LINEAR) {
    return LinearPathsToFragments(molecule, settings.GetMinFragmentSize(), settings.GetMaxFragmentSize(), copy_old_idx_property);
  };
  if (systematic_fragmentation == FragmentationSettings::SystematicFragmentation::SUBGRAPH) {
    return SubgraphsToFragments(molecule, settings.GetMinFragmentSize(), settings.GetMaxFragmentSize(), copy_old_idx_property);
  };
  if (systematic_fragmentation == FragmentationSettings::SystematicFragmentation::CIRCULAR) {
    return CircularEnvironmentsToFragments(molecule, settings.GetMinFragmentSize(), settings.GetMaxFragmentSize(), copy_old_idx_property);
  };
};

std::vector<bool> RingFragmentMask(const std::vector<RDKit::ROMOL_SPTR>& fragments) {
  std::vector<bool> mask (fragments.size());
  for (size_t fragment_idx = 0; fragment_idx < fragments.size(); ++fragment_idx) {
    const RDKit::ROMOL_SPTR fragment = fragments[fragment_idx];
    std::vector<std::vector<int>> sssr = fragment->getRingInfo()->atomRings();
    bool has_ring = false;
    if (!sssr.empty()) {
      has_ring = true;
    };
    mask[fragment_idx] = has_ring;
  };
  return mask;
};

std::vector<bool> RingPartFragmentMask(const RDKit::ROMol& molecule, const std::vector<RDKit::ROMOL_SPTR>& fragments) {
  // Determine which atoms form part of rings in the source molecule.
  std::vector<std::vector<int>> sssr = molecule.getRingInfo()->atomRings();
  std::set<unsigned> ring_atoms;
  for (const std::vector<int>& ring : sssr) {
    ring_atoms.insert(ring.begin(), ring.end());
  };
  std::set<unsigned>::const_iterator end_it = ring_atoms.end();
  // Create the fragments mask.
  std::vector<bool> mask (fragments.size());
  // Iterate over the fragments' atoms.
  for (size_t fragment_idx = 0; fragment_idx < fragments.size(); ++fragment_idx) {
    const RDKit::ROMOL_SPTR fragment = fragments[fragment_idx];
    bool is_ring_part = false;
    for (const RDKit::Atom* atom : fragment->atoms()) {
      // If one of the atoms used to be part of a ring, label the fragment as a "ring part".
      unsigned old_atom_idx = atom->getProp<unsigned>("OldIdx");
      if (ring_atoms.find(old_atom_idx) != end_it) {
        is_ring_part = true;
        break;
      };
    };
    mask[fragment_idx] = is_ring_part;
  };
  return mask;
};

Pseudofragment FragmentToPseudofragment(const RDKit::ROMol& molecule, const RDKit::ROMOL_SPTR fragment, std::vector<std::uint32_t> molecule_atom_types, bool has_ring, bool is_ring_part) {
  if (molecule.getNumAtoms() != molecule_atom_types.size()) {
    throw std::invalid_argument("Atom-atom type number mismatch");
  };
  // Collect the old indices of the fragment's atoms.
  std::set<unsigned> fragment_atom_indices;
  for (unsigned atom_idx = 0; atom_idx < fragment->getNumAtoms(); ++atom_idx) {
    const RDKit::Atom* atom = fragment->getAtomWithIdx(atom_idx);
    fragment_atom_indices.insert(atom->getProp<unsigned>("OldIdx"));
  };
  // Initialize the Pseudofragment's ConnectionsTable.
  ConnectionsTable connections_table;
  // Iterate over the molecule's atoms that correspond to the fragment.
  for (const RDKit::Atom* fragment_atom : fragment->atoms()) {
    unsigned fragment_atom_idx = fragment_atom->getIdx();
    unsigned molecule_atom_idx = fragment_atom->getProp<unsigned>("OldIdx");
    const RDKit::Atom* molecule_atom = molecule.getAtomWithIdx(molecule_atom_idx);
    // Iterate over the neighboring atoms.
    RDKit::ROMol::ADJ_ITER neighbors_it, neighbors_end_it;
    std::tie(neighbors_it, neighbors_end_it) = molecule.getAtomNeighbors(molecule_atom);
    while (neighbors_it != neighbors_end_it) {
      unsigned molecule_neighbor_idx = *neighbors_it;
      // If the neighboring atom doesn't pertain to the fragment, we define a ConnectionPoint.
      if (fragment_atom_indices.find(molecule_neighbor_idx) == fragment_atom_indices.end()) {
        const RDKit::Bond* bond = molecule.getBondBetweenAtoms(molecule_atom_idx, molecule_neighbor_idx);
        RDKit::Bond::BondType bond_type = bond->getBondType();
        std::uint32_t start_atom_type = molecule_atom_types[molecule_atom_idx];
        std::uint32_t end_atom_type = molecule_atom_types[molecule_neighbor_idx];
        Connection connection (start_atom_type, end_atom_type, bond_type_int_table.at(bond_type));
        connections_table.AddConnection(connection, fragment_atom_idx);
      };
      ++neighbors_it;
    };
  };
  std::string smiles = MakePseudomolSMILES(*fragment, connections_table);
  return Pseudofragment(*fragment, connections_table, smiles, has_ring, is_ring_part);
};

std::vector<Pseudofragment> MakePseudofragments(RDKit::ROMol& molecule, const FragmentationSettings& settings) {
  // Create a container to store the fragments.
  std::vector<RDKit::ROMOL_SPTR> fragments;

  // Calculate the molecule's atom types.
  std::vector<std::uint32_t> atom_types = GetAtomTypes(molecule, settings);

  // Flag atoms with their current indices.
  AssignOldIdxAtomProperty(molecule);

  size_t n_cyclic_fragments = 0;
  std::vector<bool> ring_part_fragment_mask;
  if (settings.FragmentRings()) {
    // If rings ought to be fragmented, fragment the entire molecule systematically.
    fragments = FragmentSystematically(molecule, settings, true);
    // Identify the fragments that are parts of rings.
    ring_part_fragment_mask = RingPartFragmentMask(molecule, fragments);
  } else {
    // Otherwise, first split the molecule into cyclic and acyclic parts.
    auto[cyclic_structures, acyclic_structures] = SplitAtRings(molecule, true);
    // The cyclic parts are converted as-is into fragments.
    fragments.insert(fragments.end(), cyclic_structures.begin(), cyclic_structures.end());
    n_cyclic_fragments = cyclic_structures.size();
    // The acyclic structures are subjected to further systematic fragmentation.
    for (const RDKit::ROMOL_SPTR acyclic_structure : acyclic_structures) {
      std::vector<RDKit::ROMOL_SPTR> acyclic_fragments = FragmentSystematically(*acyclic_structure, settings, true);
      fragments.insert(fragments.end(), acyclic_fragments.begin(), acyclic_fragments.end());
    };
    // If ring structures are treated as fragments, no fragment can be a ring part.
    ring_part_fragment_mask.resize(fragments.size(), false);
  };

  // Convert the fragments into Pseudofragment objects.
  std::vector<Pseudofragment> pseudofragments;
  pseudofragments.reserve(fragments.size());
  for (size_t fragment_idx = 0; fragment_idx < fragments.size(); ++fragment_idx) {
    const RDKit::ROMOL_SPTR fragment = fragments[fragment_idx];
    bool has_ring = false;
    if (fragment_idx < n_cyclic_fragments) {
      has_ring = true;
    };
    bool is_ring_part = ring_part_fragment_mask[fragment_idx];
    Pseudofragment pseudofragment = FragmentToPseudofragment(molecule, fragment, atom_types, has_ring, is_ring_part);
    pseudofragments.push_back(pseudofragment);
  };
  return pseudofragments;
};
//
// void AssignImmutableAtomIndices(RDKit::ROMol& mol) {
//   unsigned immutable_idx = 0u;
//   for (RDKit::Atom* atom : mol.atoms()) {
//     atom->setProp<unsigned>("ImmutableIdx", immutable_idx);
//     ++immutable_idx;
//   };
// };
//
// void AssignMMFFAtomTypes(RDKit::ROMol& mol) {
//   RDKit::MMFF::MMFFMolProperties mmff_properties = RDKit::MMFF::MMFFMolProperties(mol);
//   for (RDKit::Atom* atom : mol.atoms()) {
//     unsigned atom_idx = atom->getIdx();
//     unsigned mmff_atom_type = mmff_properties.getMMFFAtomType(atom_idx);
//     atom->setProp<unsigned>("MMFF_atom_type", mmff_atom_type);
//   };
// };
//
// void FlagRingAtoms(RDKit::ROMol& mol) {
//   const RDKit::RingInfo* ring_info = mol.getRingInfo();
//   for (RDKit::Atom* atom : mol.atoms()) {
//     if (ring_info->numAtomRings(atom->getIdx()) > 0) {
//       atom->setProp<bool>("WasInRing", true);
//     } else {
//       atom->setProp<bool>("WasInRing", false);
//     };
//   };
// };
//
// void FlagAtoms(RDKit::ROMol& mol) {
//   unsigned immutable_idx = 0u;
//   RDKit::MMFF::MMFFMolProperties mmff_properties = RDKit::MMFF::MMFFMolProperties(mol);
//   const RDKit::RingInfo* ring_info = mol.getRingInfo();
//   for (RDKit::Atom* atom : mol.atoms()) {
//     unsigned atom_idx = atom->getIdx();
//     unsigned mmff_atom_type = mmff_properties.getMMFFAtomType(atom_idx);
//     atom->setProp<unsigned>("ImmutableIdx", immutable_idx);
//     atom->setProp<unsigned>("MMFF_atom_type", mmff_atom_type);
//     if (ring_info->numAtomRings(atom_idx) > 0) {
//       atom->setProp<bool>("WasInRing", true);
//     } else {
//       atom->setProp<bool>("WasInRing", false);
//     };
//     ++immutable_idx;
//   };
// };
//
// void CopyFlags(const RDKit::Atom* source, RDKit::Atom* destination) {
//   destination->setProp<unsigned>("ImmutableIdx", source->getProp<unsigned>("ImmutableIdx"));
//   destination->setProp<unsigned>("MMFF_atom_type", source->getProp<unsigned>("MMFF_atom_type"));
//   destination->setProp<bool>("WasInRing", source->getProp<bool>("WasInRing"));
// };
//
// // Convert a dummy molecule (i.e. a RDKit::ROMol with Connection-encoding pseudoatoms) into a Pseudofragment.
// Pseudofragment DummyMolToPseudofragment(const RDKit::ROMol& mol, boost::format& formatter) {
//   // Initialize the pseudomol of the Pseudofragment.
//   RDKit::RWMol pseudomol = RDKit::RWMol(mol);
//
//   // Initialize the Pseudofragment's connection table.
//   ConnectionsTable connections;
//
//   // Initialize a RingInfo object.
//   bool has_ring = false;
//   RDKit::RingInfo* ring_info = mol.getRingInfo();
//
//   // Initialize a vector to store the indices of the molecule's pseudoatoms.
//   std::vector<unsigned> pseudoatom_indices;
//
//   // Loop over the molecule's atoms.
//   for (RDKit::Atom* atom : pseudomol.atoms()) {
//     // For all encountered pseudoatoms (i.e. atoms with an atomic number of 0),
//     // store their index and add the Connection they encode to the connections table.
//     if (atom->getAtomicNum() == 0) {
//       std::uint8_t atom_type, neighbor_type;
//       std::tie(neighbor_type, atom_type) = DecodeConnection(atom->getIsotope());
//       unsigned atom_idx = atom->getIdx();
//       pseudoatom_indices.push_back(atom_idx);
//       // The pseudoatom has only a single neighbor so no true iteration is necessary.
//       RDKit::ROMol::ADJ_ITER nbr_it, nbr_end_it;
//       std::tie(nbr_it, nbr_end_it) = mol.getAtomNeighbors(atom);
//       unsigned nbr_idx = *nbr_it;
//       RDKit::Bond::BondType bond_type = mol.getBondBetweenAtoms(atom_idx, nbr_idx)->getBondType();
//       unsigned bond_type_int = bond_type_int_table.at(bond_type);
//       Connection connection(neighbor_type, atom_type, bond_type_int, formatter);
//       connections.AddConnection(connection, nbr_idx);
//       // For all regular atoms, label them according to whether they are part of a
//       // ring system or not.
//     } else {
//       if (!has_ring && ring_info->numAtomRings(atom->getIdx()) > 0) {
//         has_ring = true;
//       };
//     };
//   };
//
//   // Delete the pseudoatoms in sorted descending order (to avoid triggering atom
//   // reindexing).
//   std::sort(pseudoatom_indices.begin(), pseudoatom_indices.end(), std::greater<unsigned>());
//   for (unsigned idx : pseudoatom_indices) {
//     pseudomol.removeAtom(idx);
//   };
//
//   // Encode the dummy molecule as a SMILES string and generate the Pseudofragment.
//   std::string smiles = RDKit::MolToSmiles(mol);
//   Pseudofragment pseudofragment(pseudomol, connections, smiles, has_ring, false);
//   return pseudofragment;
// };
//
// // Function to convert an RDKit::Atom into a Pseudofragment.
// Pseudofragment AtomToPseudofragment(const RDKit::ROMol& mol, const RDKit::Atom* atom, boost::format& formatter) {
//   // Initialize the Pseudofragment's pseudomol and connections table.
//   RDKit::RWMol pseudomol;
//   ConnectionsTable connections;
//
//   // Initialize the dummy molecule used to create the Connection-encoding SMILES.
//   RDKit::RWMol dummymol;
//
//   // Add the base Atom to the pseudomol and dummy molecule. Apparently its
//   // important to create a copy of the atom before adding it, otherwise it
//   // adds an empty pseudoatom.
//   RDKit::Atom atom_copy = *atom;
//   unsigned pseudomol_atom_idx = pseudomol.addAtom(&atom_copy);
//   unsigned dummymol_atom_idx = dummymol.addAtom(&atom_copy);
//
//   // Loop over the neighboring Atoms of the base Atom in the owning Mol.
//   RDKit::ROMol::ADJ_ITER nbr_it, nbr_end_it;
//   std::tie(nbr_it, nbr_end_it) = mol.getAtomNeighbors(atom);
//   while (nbr_it != nbr_end_it) {
//     unsigned nbr_idx = *nbr_it;
//     const RDKit::Atom* neighbor = mol.getAtomWithIdx(nbr_idx);
//
//     // If the neighboring atom is a pseudoatom, decode its isotope into atom types.
//     // If not, get the atom types from the Atom properties.
//     std::uint8_t atom_type, neighbor_type;
//     if (neighbor->getAtomicNum() == 0) {
//       std::tie(atom_type, neighbor_type) = DecodeConnection(neighbor->getIsotope());
//     }
//     else {
//       atom_type = atom_copy.getProp<unsigned>("MMFF_atom_type");
//       neighbor_type = neighbor->getProp<unsigned>("MMFF_atom_type");
//     };
//
//     // Create a Connection object symbolizing the bond between the Atom and the Neighbor.
//     RDKit::Bond::BondType bond_type = mol.getBondBetweenAtoms(atom->getIdx(), nbr_idx)->getBondType();
//     unsigned bond_type_int = bond_type_int_table.at(bond_type);
//     Connection connection(atom_type, neighbor_type, bond_type_int, formatter);
//
//     // Add the Connection to the connections table.
//     connections.AddConnection(connection, pseudomol_atom_idx);
//
//     // Create a pseudoatom with an isotope encoding the start- and end atom types
//     // of the Connection point and add it to the dummy molecule.
//     RDKit::Atom pseudoatom(0);
//     pseudoatom.setIsotope(connection.GetEncoded());
//     unsigned pseudoatom_idx = dummymol.addAtom(&pseudoatom);
//     dummymol.addBond(dummymol_atom_idx, pseudoatom_idx, bond_type);
//
//     ++nbr_it;
//   };
//
//   // Encode the dummy molecule as a SMILES string and generate the Pseudofragment.
//   std::string smiles = RDKit::MolToSmiles(dummymol);
//   Pseudofragment pseudofragment(pseudomol, connections, smiles, false, false);
//   return pseudofragment;
// };
