#include <fstream>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <omp.h>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include "PseudofragmentDB.hpp"

// Function to divide a range of size N into equally sized chunks.
std::vector<std::pair<unsigned, unsigned>> EquallySizedChunks(unsigned size, unsigned n_chunks, bool one_based_index = false) {
  std::vector<std::pair<unsigned, unsigned>> chunks;
  unsigned chunk_size = size / n_chunks;
  // Define up to the last chunk.
  unsigned begin_idx = 0, end_idx = 0;
  if (one_based_index) {
    ++begin_idx;
  };
  for (unsigned n = 0; n < n_chunks - 1; ++n) {
    end_idx = begin_idx + chunk_size - 1;
    chunks.push_back(std::make_pair(begin_idx, end_idx));
    begin_idx = end_idx + 1;
  };
  // Define the last chunk.
  end_idx = size;
  if (!one_based_index) {
    --end_idx;
  };
  chunks.push_back(std::make_pair(begin_idx, end_idx));
  return chunks;
};

// Function that concatenates rings (as detected through SSSR) sharing atleast
// one atom into a single ring system.
std::vector<std::vector<int>> DefineRingSystems(const std::vector<std::vector<int>>& rings) {
  // Create a vector containing the ring indices of the rings detected through
  // SSSR ring perception.
  std::vector<int> ring_indices(rings.size());
  std::iota(ring_indices.begin(), ring_indices.end(), 0);

  // Initialize a vector to store the fused-ring ring systems.
  std::vector<std::vector<int>> ring_systems;

  // Loop over each (remaining) ring in "ring_indices" and check whether there
  // is any overlap between its atom indices and the atom indices of the other
  // (remaining) rings. If so, fuse both rings and repeat the process until
  // the ring system can no longer be expanded further. At said point, the
  // indices of the fused rings are removed from the "ring_indices" vector,
  // and if any rings remain unclassified the process is repeated with a new
  // ring as starting point.

  // Select a ring starting point for the ring system.
next_ring_system:
  for (int ring_idx1 : ring_indices) {
    if (ring_indices.size() > 0) {
      std::vector<int> ring_system (rings[ring_idx1]);
      std::vector<int> ring_indices_system { ring_idx1 };
      // Remove the chosen ring from the rings to evaluate (i.e. "remaining rings").
      ring_indices.erase(std::remove(ring_indices.begin(), ring_indices.end(), ring_idx1), ring_indices.end());

    // Loop over the remaining rings again to find more ring system members.
    check_further_expansion:
      bool expanded = false;
      for (int ring_idx2 : ring_indices) {
        const std::vector<int>& ring2 = rings[ring_idx2];
        // Check whether the chosen pair of rings have any atoms in common.
        // If so, fuse them.
        for (int atom_idx : ring_system) {
          if (std::find(ring2.begin(), ring2.end(), atom_idx) != ring2.end()) {
            ring_system.insert(ring_system.end(), ring2.begin(), ring2.end());
            ring_indices_system.push_back(ring_idx2);
            ring_indices.erase(std::remove(ring_indices.begin(), ring_indices.end(), ring_idx2), ring_indices.end());
            expanded = true;
            break;
          };
        };
      };

      // If the ring system was expanded during this iteration, attempt
      // another iteration to see if the ring system can be expanded further.
      if (expanded) {
        goto check_further_expansion;
      // If the ring system wasn't expanded it has reached its maximum size and
      // is stored in the container vector after removing duplicate atom indices.
      } else {
        std::sort(ring_system.begin(), ring_system.end());
        ring_system.erase(std::unique(ring_system.begin(), ring_system.end()), ring_system.end());
        ring_systems.push_back(ring_system);
        goto next_ring_system;
      };
    } else {
      break;
    };
  };
  return ring_systems;
};

// Class to store a Murcko decomposition of a molecule.
class MurckoDecomposition {
public:
  // Member variables that store the indices of the decomposition.
  std::vector<int> linkers, side_chains;
  std::vector<std::vector<int>> rings;

  // Class constructor
  MurckoDecomposition(const RDKit::ROMOL_SPTR mol, const RDKit::RingInfo* ring_info) {
    // Extract the Murcko Scaffold as a RDKit::Mol and convert it into a
    // vector of indices.
    RDKit::ROMol* framework = RDKit::MurckoDecompose(*mol);
    std::vector<unsigned> framework_indices;
    for (RDKit::ROMol::AtomIterator ai = framework->beginAtoms(); ai != framework->endAtoms(); ++ai) {
      framework_indices.push_back((*ai)->getProp<unsigned>("ImmutableIdx"));
    };
    delete framework;

    // Create a vector of integer vectors containing the SSSR.
    std::vector<std::vector<int>> sssr (ring_info->atomRings());

    // Combine the SSSR into the largest possible ring systems.
    rings = DefineRingSystems(sssr);

    // Define the atom indices that form part of linkers.
    for (unsigned atom_idx : framework_indices) {
      bool is_linker = true;
      for (const auto& ring_system : rings) {
        if (std::find(ring_system.begin(), ring_system.end(), atom_idx) != ring_system.end()) {
          is_linker = false;
          break;
        };
      };
      if (is_linker) {
        linkers.push_back(atom_idx);
      };
    };

    // Define the atom indices that form part of side chains.
    for (RDKit::ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms(); ++ai) {
      int atom_idx = (*ai)->getIdx();
      if (std::find(framework_indices.begin(), framework_indices.end(), atom_idx) == framework_indices.end()) {
        side_chains.push_back(atom_idx);
      };
    };
  };
};

// Function to write an entry into the log.
void Log(std::ofstream& logger, const std::string& severity, const std::string& message, std::string smiles = "") {
  auto time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  if (smiles == "") {
    logger << "[" << std::strtok(std::ctime(&time), "\n") << "] " << severity << ": " << message << "\n";
  } else {
    logger << "[" << std::strtok(std::ctime(&time), "\n") << "] [" << smiles << "] " << severity << ": " << message << "\n";
  };
};

// Function to calculate the MMFF94 atom types of a RDKit::ROMol's atoms and store
// them as a property.
int AssignMMFFAtomTypes(RDKit::ROMOL_SPTR mol, bool assign_immutable_indices = false) {
  RDKit::MMFF::MMFFMolProperties MMFF_properties = RDKit::MMFF::MMFFMolProperties(*mol);
  for (unsigned idx = 0, max_idx = mol->getNumAtoms(); idx < max_idx; ++idx) {
    unsigned MMFF_atom_type = MMFF_properties.getMMFFAtomType(idx);
    RDKit::Atom* atom = mol->getAtomWithIdx(idx);
    atom->setProp("MMFF_atom_type", MMFF_atom_type);
    if (assign_immutable_indices) {
      atom->setProp("ImmutableIdx", atom->getIdx());
    };
  };
  return 0;
};


// Function to split a molecule by cutting the bonds connecting ring systems to
// the rest of the molecule.
std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> SplitAtRings(const RDKit::ROMOL_SPTR mol, const MurckoDecomposition& murcko) {
  // Initialize containers to store the data needed to split the molecule informaively.
  std::vector<unsigned> split_bond_indices;
  std::vector<RDKit::Bond::BondType> split_bond_types;
  std::vector<std::pair<unsigned, unsigned>> split_bond_labels;
  // Loop over the ring systems according to the Murcko decomposition.
  for (const auto& ring_system : murcko.rings) {
    // Loop over the atoms in the ring system.
    for (int atom_idx : ring_system) {
      const RDKit::Atom* atom = mol->getAtomWithIdx(atom_idx);
      // Loop over the neighboring atoms.
      RDKit::ROMol::ADJ_ITER nbr_idx, end_nbr_idx;
      boost::tie(nbr_idx, end_nbr_idx) = mol->getAtomNeighbors(atom);
      while (nbr_idx != end_nbr_idx) {
        // If the neighbor isn't part of a ring system, the bond between
        // both ought to be split.
        if (std::find(ring_system.begin(), ring_system.end(), *nbr_idx) == ring_system.end()) {
          const RDKit::Bond* bond = mol->getBondBetweenAtoms(atom_idx, *nbr_idx);
          unsigned bond_idx = bond->getIdx();
          // Flag the identified bond for cleavage only if it hasn't been
          // flagged before. This could happen if two rings are connected
          // directly through a bond.
          if (std::find(split_bond_indices.begin(), split_bond_indices.end(), bond_idx) == split_bond_indices.end()) {
            split_bond_indices.push_back(bond->getIdx());
            split_bond_types.push_back(bond->getBondType());
            // Use as post-fragmentation labels Connection encoding integers.
            Connection connection(bond->getBeginAtom()->getProp<unsigned>("MMFF_atom_type"), bond->getEndAtom()->getProp<unsigned>("MMFF_atom_type"), 0);
            Connection mirrored_connection = connection.Mirror();
            std::pair<unsigned, unsigned> p (mirrored_connection.GetEncoded(), connection.GetEncoded());
            split_bond_labels.push_back(p);
          };
        };
        ++nbr_idx;
      };
    };
  };

  // Fragment the molecule by cutting the specified bonds.
  RDKit::ROMol* fragmented = RDKit::MolFragmenter::fragmentOnBonds(*mol, split_bond_indices, true, &split_bond_labels, &split_bond_types);

  // Extract the resulting fragments as separate RDKit::ROMols.
  std::vector<int> frag_map;
  std::vector<std::vector<int>> atom_map;
  std::vector<RDKit::ROMOL_SPTR> blocks = RDKit::MolOps::getMolFrags(*fragmented, true, &frag_map, &atom_map);
  delete fragmented;

  // Determine which of the fragments are rings.
  int i = 0;
  std::vector <int> ring_block_indices;
  for (const auto& block_atom_map : atom_map) {
    for (int atom_idx : block_atom_map) {
      for (const auto& ring_system : murcko.rings) {
        if (std::find(ring_system.begin(), ring_system.end(), atom_idx) != ring_system.end()) {
          ring_block_indices.push_back(i);
          goto next_block;
        };
      };
    };
  next_block:
    ++i;
  };

  // Separate the fragments into ring and acyclic fragments.
  std::sort(ring_block_indices.begin(), ring_block_indices.end(), std::greater<int>());
  std::vector<RDKit::ROMOL_SPTR> ring_blocks;
  for (int idx : ring_block_indices) {
    ring_blocks.push_back(std::move(blocks[idx]));
    blocks.erase(blocks.begin() + idx);
  };

  // Return the resulting fragments.
  std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> output (blocks, ring_blocks);
  return output;
};

// Convert a dummy molecule (i.e. a RDKit::ROMol with Connection-encoding pseudoatoms) into a Pseudofragment.
Pseudofragment DummyMolToPseudofragment(const RDKit::ROMOL_SPTR mol, const RDKit::RingInfo* ring_info, bool has_ring, boost::format& formatter) {
  // Initialize the pseudomol of the Pseudofragment.
  RDKit::RWMol pseudomol = RDKit::RWMol(*mol);

  // Initialize the Pseudofragment's connection table.
  ConnectionsTable connections;

  // Initialize a flag to signal whether the dummymol was part of a ring system.
  bool ring_part = false;

  // Initialize a vector to store the indices of the molecule's pseudoatoms.
  std::vector<unsigned> pseudoatom_indices;

  // Loop over the molecule's atoms.
  for (RDKit::RWMol::AtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    // For all encountered pseudoatoms (i.e. atoms with an atomic number of 0),
    // store their index and add the Connection they encode to the connections table.
    if ((*ai)->getAtomicNum() == 0) {
      std::uint8_t atom_type, neighbor_type;
      std::tie(neighbor_type, atom_type) = DecodeConnection((*ai)->getIsotope());
      unsigned atom_idx = (*ai)->getIdx();
      pseudoatom_indices.push_back(atom_idx);
      // The pseudoatom has only a single neighbor so no true iteration is necessary.
      RDKit::ROMol::ADJ_ITER nbr_idx, end_nbr_idx;
      boost::tie(nbr_idx, end_nbr_idx) = mol->getAtomNeighbors(*ai);
      RDKit::Bond::BondType bond_type = mol->getBondBetweenAtoms(atom_idx, *nbr_idx)->getBondType();
      unsigned bond_type_int = bond_type_int_table.at(bond_type);
      Connection connection(neighbor_type, atom_type, bond_type_int, formatter);
      connections.AddConnection(connection, *nbr_idx);
      // For all regular atoms, label them according to whether they are part of a
      // ring system or not.
    }
    else {
      if (ring_info->numAtomRings((*ai)->getIdx()) > 0) {
        (*ai)->setProp<bool>("WasInRing", true);
        if (!has_ring) {
          ring_part = true;
        };
      } else {
        (*ai)->setProp<bool>("WasInRing", false);
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
  std::string smiles = RDKit::MolToSmiles(*mol);
  Pseudofragment pseudofragment(pseudomol, connections, smiles, has_ring, ring_part);
  return pseudofragment;
};

// Function that converts a "submol" (i.e. a molecule created by extracting a
// subgraph from a RDKit::ROMol into a Pseudofragment.
Pseudofragment SubmolToPseudofragment(const RDKit::ROMol* submol, const RDKit::ROMOL_SPTR sourcemol, const std::map<int, int>& atom_map, const RDKit::RingInfo* ring_info, bool has_ring, boost::format& formatter) {
  // Initialize a dummy molecule to represent the Pseudofragment. This dummymol
  // is used to generate the connection-encoding canonical SMILES.
  RDKit::RWMol dummymol = RDKit::RWMol(*submol);

  // Initialize the Pseudofragment's connection table.
  ConnectionsTable connections;

  // Initialize a flag to signal whether the submol is part of a ring system in the sourcemol.
  bool ring_part = false;

  // Loop over the entries in the atom map and extract the sourcemol and submol Atoms.
  for (const auto& map_entry : atom_map) {
    unsigned source_atom_idx = map_entry.first;
    unsigned submol_atom_idx = map_entry.second;

    const RDKit::Atom* submol_atom = submol->getAtomWithIdx(submol_atom_idx);
    const RDKit::Atom* source_atom = sourcemol->getAtomWithIdx(source_atom_idx);

    // Loop over the neighbors of the sourcemol's Atom.
    RDKit::ROMol::ADJ_ITER nbr_iter, end_nbr_iter;
    std::tie(nbr_iter, end_nbr_iter) = sourcemol->getAtomNeighbors(source_atom);
    while (nbr_iter != end_nbr_iter) {
      // If the neighbor isn't part of the submol (i.e. isn't in the atom map of
      // the submol) a Connection to the neighbor must be created.
      unsigned nbr_idx = *nbr_iter;
      if (atom_map.count(nbr_idx) == 0) {

        // If the sourcemol Atom was in a ring, label the submol atom as a former ring member.
        if (ring_info->numAtomRings(source_atom_idx) > 0) {
          submol_atom->setProp<bool>("WasInRing", true);
          if (!has_ring) {
            ring_part = true;
          };
        } else {
          submol_atom->setProp<bool>("WasInRing", false);
        };

        const RDKit::Atom* neighbor = (*sourcemol)[nbr_idx];

        // Create the Connection object and add it to the connections table.
        std::uint8_t atom_type, neighbor_type;
        RDKit::Bond::BondType bond_type = sourcemol->getBondBetweenAtoms(source_atom_idx, nbr_idx)->getBondType();
        if (neighbor->getAtomicNum() == 0) {
          std::tie(atom_type, neighbor_type) = DecodeConnection(neighbor->getIsotope());
        } else {
          atom_type = source_atom->getProp<unsigned>("MMFF_atom_type");
          neighbor_type = neighbor->getProp<unsigned>("MMFF_atom_type");
        };
        unsigned bond_type_int = bond_type_int_table.at(bond_type);
        Connection connection(atom_type, neighbor_type, bond_type_int, formatter);
        connections.AddConnection(connection, submol_atom_idx);

        // Add a Connection-encoding pseudoatom to the dummymol.
        RDKit::Atom pseudoatom(0);
        pseudoatom.setIsotope(connection.GetEncoded());
        unsigned pseudoatom_idx = dummymol.addAtom(&pseudoatom);
        dummymol.addBond(submol_atom_idx, pseudoatom_idx, bond_type);
      };
      ++nbr_iter;
    };
  };

  // Encode the dummy molecule as a SMILES string and generate the Pseudofragment.
  std::string smiles = RDKit::MolToSmiles(dummymol);
  Pseudofragment pseudofragment(*submol, connections, smiles, has_ring, ring_part);
  return pseudofragment;
};

// Function to convert an RDKit::Atom into a Pseudofragment.
Pseudofragment AtomToPseudofragment(const RDKit::ROMOL_SPTR mol, const RDKit::Atom* atom, const RDKit::RingInfo* ring_info, boost::format& formatter) {
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

  // Set an Atom property flag to signal that whether the atom is part of ring or
  // not. This is done to be coherent with the multi-atom Pseudofragments.
  bool ring_part = false;
  if (ring_info->numAtomRings(atom->getIdx()) > 0) {
    pseudomol.getAtomWithIdx(pseudomol_atom_idx)->setProp<bool>("WasInRing", true);
    ring_part = true;
  } else {
    pseudomol.getAtomWithIdx(pseudomol_atom_idx)->setProp<bool>("WasInRing", false);
  };

  // Loop over the neighboring Atoms of the base Atom in the owning Mol.
  RDKit::ROMol::ADJ_ITER nbr_iter, end_nbr_iter;
  std::tie(nbr_iter, end_nbr_iter) = mol->getAtomNeighbors(atom);
  while (nbr_iter != end_nbr_iter) {
    unsigned nbr_idx = *nbr_iter;
    RDKit::Atom* neighbor = mol->getAtomWithIdx(nbr_idx);

    // If the neighboring atom is a pseudoatom, decode its isotope into atom types.
    // If not, get the atom types from the Atom properties.
    std::uint8_t atom_type, neighbor_type;
    if (neighbor->getAtomicNum() == 0) {
      std::tie(atom_type, neighbor_type) = DecodeConnection(neighbor->getIsotope());
    } else {
      atom_type = atom_copy.getProp<unsigned>("MMFF_atom_type");
      neighbor_type = neighbor->getProp<unsigned>("MMFF_atom_type");
    };

    // Create a Connection object symbolizing the bond between the Atom and the Neighbor.
    RDKit::Bond::BondType bond_type = mol->getBondBetweenAtoms(atom->getIdx(), neighbor->getIdx())->getBondType();
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

    ++nbr_iter;
  };

  // Encode the dummy molecule as a SMILES string and generate the Pseudofragment.
  std::string smiles = RDKit::MolToSmiles(dummymol);
  Pseudofragment pseudofragment(pseudomol, connections, smiles, false, ring_part);
  return pseudofragment;
};

// Function to fragment a RDKit::ROMol and convert it into Pseudofragments.
std::vector<Pseudofragment> FragmentMol(RDKit::ROMOL_SPTR mol, unsigned minsize, unsigned maxsize, bool ringfrag, bool purge, boost::format& formatter) {
  // Verify that the molecule is large enough to be fragmented according to the
  // user specified criteria. Otherwise the molecule is skipped.
  std::vector<Pseudofragment> pseudofragments;
  unsigned nbonds = mol->getNumBonds();
  if (nbonds < minsize) {
    return pseudofragments;
  };

  // Initialize a container to store the RDKit::ROMol objects that will be
  // subjected to fragmentation.
  std::vector<RDKit::ROMOL_SPTR> blocks_to_fragment;

  // Perceive the rings of the molecule.
  RDKit::RingInfo* ring_info = mol->getRingInfo();
  int nrings = ring_info->numRings();

  // Assign MMFF94 atom types. These are the atom types used to define Connections.
  AssignMMFFAtomTypes(mol, true);

  // If the molecule contains rings and ring fragmentation is enabled, fragment
  // the molecule systematically.
  if (nrings > 0 && ringfrag) {
    // Generate 1-Atom Pseudofragments (if specified).
    if (minsize == 0) {
      for (RDKit::ROMol::AtomIterator ai = mol->beginAtoms(); ai != mol->endAtoms(); ++ai) {
        RDKit::Atom* atom = *ai;
        Pseudofragment pseudofragment = AtomToPseudofragment(mol, atom, ring_info, formatter);
        pseudofragments.push_back(pseudofragment);
      };
    };
    // Generate subgraph-based Pseudofragments of the specified sizes. Note that
    // when the "minsize" is set to 0 it defaults to 1 during function execution.
    std::map<int, std::list<std::vector<int>>> subgraph_paths = RDKit::findAllSubgraphsOfLengthsMtoN(*mol, minsize, maxsize);
    for (const auto& paths : subgraph_paths) {
      for (const auto& path : paths.second) {
        std::map<int, int> atom_map;
        RDKit::ROMol* submol = RDKit::Subgraphs::pathToSubmol(*mol, path, false, atom_map);
        Pseudofragment pseudofragment = SubmolToPseudofragment(submol, mol, atom_map, ring_info, false, formatter);
        pseudofragments.push_back(pseudofragment);
        delete submol;
      };
    };

  // If the molecule contains rings but ring fragmentation is disabled, attempt
  // to isolate the rings from the linkers and side chains by performing a Murcko
  // decomposition.
  } else if (nrings > 0 && !ringfrag) {
    const MurckoDecomposition murcko(mol, ring_info);
    // If the molecule is made up of exclusively rings and ring fragmentation is
    // disabled, the molecule is skipped.
    if (murcko.linkers.empty() && murcko.side_chains.empty()) {
      return pseudofragments;
      // If not, cut the bonds connecting the rings with the rest of the molecule.
    } else {
      std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> blocks = SplitAtRings(mol, murcko);
      // Separate the resulting "blocks" into rings and others.
      blocks_to_fragment = blocks.first;
      const std::vector<RDKit::ROMOL_SPTR>& ring_blocks = blocks.second;
      // Ring blocks are transformed directly into Pseudofragments. Acyclic "blocks"
      // are sent to systematic fragmentation.
      for (RDKit::ROMOL_SPTR ring_block : ring_blocks) {
        ring_info = ring_block->getRingInfo();
        Pseudofragment pseudofragment = DummyMolToPseudofragment(ring_block, ring_info, true, formatter);
        pseudofragments.push_back(pseudofragment);
      };
    };

    // If the molecule lacks rings, it can be fragmented systematically in its entirety.
  } else if (nrings == 0) {
    blocks_to_fragment.push_back(mol);
  };

  // Loop over the acyclic "blocks" and subject them to systematic fragmentation.
  for (RDKit::ROMOL_SPTR open_block : blocks_to_fragment) {
    ring_info = open_block->getRingInfo();
    // Only blocks with a size greater than 1 can be fragmented.
    if (open_block->getNumHeavyAtoms() > 1) {
      // Generate 1-Atom Pseudofragments (if specified).
      if (minsize == 0) {
        for (RDKit::ROMol::AtomIterator ai = open_block->beginAtoms(); ai != open_block->endAtoms(); ++ai) {
          RDKit::Atom* atom = *ai;
          // If the Atom isn't a pseudoatom, convert it into a Pseudofragment.
          if (atom->getAtomicNum() != 0) {
            Pseudofragment pseudofragment = AtomToPseudofragment(open_block, atom, ring_info, formatter);
            pseudofragments.push_back(pseudofragment);
          };
        };
      };
      // Create a vector to store the atom indices of the pseudoatoms in the "block".
      // This is necessary in order to filter out the subgraphs that contain them later on.
      std::vector<int> pseudoatom_bond_indices;
      for (RDKit::ROMol::BondIterator bi = open_block->beginBonds(); bi != open_block->endBonds(); ++bi) {
        if (!(*bi)->getBeginAtom()->getAtomicNum() || !(*bi)->getEndAtom()->getAtomicNum()) {
          pseudoatom_bond_indices.push_back((*bi)->getIdx());
        };
      };
      // Generate subgraph-based Pseudofragments of the specified sizes. Note that
      // when the "minsize" is set to 0 it defaults to 1 during function execution.
      std::map<int, std::list<std::vector<int>>> subgraph_paths = RDKit::findAllSubgraphsOfLengthsMtoN(*open_block, minsize, maxsize);
      for (const auto& paths : subgraph_paths) {
        for (const auto& path : paths.second) {
          // Check if the path is valid. If a bond involving a pseudoatom is part
          // of it the path is invalid.
          bool valid_path = true;
          for (int i : path) {
            if (std::find(pseudoatom_bond_indices.begin(), pseudoatom_bond_indices.end(), i) != pseudoatom_bond_indices.end()) {
              valid_path = false;
              break;
            };
          };
          // If the path is valid, convert it into a Pseudofragment.
          if (valid_path) {
            std::map<int, int> atom_map;
            RDKit::ROMol* submol = RDKit::Subgraphs::pathToSubmol(*open_block, path, false, atom_map);
            Pseudofragment pseudofragment = SubmolToPseudofragment(submol, open_block, atom_map, ring_info, false, formatter);
            pseudofragments.push_back(pseudofragment);
            delete submol;
          };
        };
      };
    // If the size of the block is 1, whether it is added or not to the database
    // depends on the user specified "purge" argument.
    } else {
      if (!purge) {
        Pseudofragment pseudofragment = DummyMolToPseudofragment(open_block, ring_info, false, formatter);
        pseudofragments.push_back(pseudofragment);
      };
    };
  };

  // Return the generated Pseudofragments.
  return pseudofragments;
};

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser and parse the command line arguments.
  std::string input, output;
  unsigned minsize, maxsize;
  bool ringfrag, purge;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to the input file containing molecules in SMILES or SDF format.")
    ("output,o", boost::program_options::value<std::string>(&output)->required(),
      "Path to the output SQLite3 database containing the fragments.")
    ("minsize,n", boost::program_options::value<unsigned>(&minsize)->default_value(0),
      "Minimum size of the generated fragments (in number of bonds). If set to 0 individual atoms are considered fragments.")
    ("maxsize,x", boost::program_options::value<unsigned>(&maxsize)->default_value(4),
      "Maximum size of the generated fragments (in number of bonds).")
    ("ringfrag,r", boost::program_options::bool_switch(&ringfrag)->default_value(false),
      "Flag to enable ring fragmentation. If not provided ring systems will be isolated and treated as a single fragment.")
    ("purge,p", boost::program_options::bool_switch(&purge)->default_value(false),
      "Flag to remove single-atom fragments resulting from ring isolation (regardless of the --minsize value).");

  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("input", 1);
  positionals_description.add("output", 1);
  boost::program_options::variables_map vm;
  boost::program_options::command_line_parser parser(argc, argv);
  parser.options(description);
  parser.positional(positionals_description);
  boost::program_options::parsed_options parsed_options = parser.run();
  boost::program_options::store(parsed_options, vm);
  if (vm.count("help")) {
    std::cout << description;
    return 1;
  };
  boost::program_options::notify(vm);

  if (minsize > maxsize) {
    throw std::runtime_error(std::string("--minsize can not be greater than --maxsize."));
  };

  // Determine the number of OpenMP threads to use.
  int n_threads;
  char* omp_num_threads = std::getenv("OMP_NUM_THREADS");
  if (omp_num_threads == nullptr) {
    n_threads = omp_get_max_threads();
  } else {
    n_threads = std::stoi(omp_num_threads);
  };
  omp_set_num_threads(n_threads);

  // Set up a logger.
  std::string log_file = output + ".log";
  std::ofstream logger(log_file);
  Log(logger, "INFO", "Application started.");
  std::stringstream args_message;
  args_message << "Command line arguments: --input=" << input << ", --output="
    << output << ", --minsize=" << minsize << ", --maxsize=" << maxsize
    << ", --ringfrag=" << ringfrag << ", --purge=" << purge;
  Log(logger, "INFO", args_message.str());

  // Open a connection to the output database.
  sqlite3* database;
  // If the output database already exists, reopen it.
  // NOTE: std::filesystem requires standard C++17.
  if (std::filesystem::exists(output)) {
    sqlite3_open(output.c_str(), &database);
    Log(logger, "INFO", "Existing database was opened.");
  // If it doesn't exist yet, initialize it.
  } else {
    database = InitializePseudofragmentDB(output);
    Log(logger, "INFO", "New database was initialized.");
  };

  // Start a transaction to insert all the information.
  sqlite3_exec(database, "BEGIN TRANSACTION", NULL, NULL, NULL);
  Log(logger, "INFO", "Started database transaction.");

  // Create a molecule input stream based on the extension of the input file.
  std::filesystem::path extension = std::filesystem::path(input).extension();

  // If the input file is a SMILES file:
  if (extension == ".smi" || extension == ".ism" || extension == ".smiles") {
    RDKit::SmilesMolSupplier supplier(input);
    Log(logger, "INFO", "Created molecule input stream.");

    // Initialize a progress bar.
    unsigned n_mols = supplier.length();
    supplier.reset();
    boost::progress_display progress(n_mols);

    // Divide the input file into N equally sized chunks, where N is the number of
    // OpenMP threads.
    std::vector<std::pair<unsigned, unsigned>> chunks = EquallySizedChunks(n_mols, n_threads);

    #pragma omp parallel for
    for (int thrid = 0; thrid < n_threads; ++thrid) {
      // Prepare the statements required to work with the database.
      sqlite3_stmt* insert_source;
      sqlite3_stmt* insert_pseudofragment;
      sqlite3_stmt* insert_connection;
      sqlite3_stmt* insert_source_pseudofragment_relationship;
      sqlite3_stmt* insert_pseudofragment_connection_relationship;
      sqlite3_stmt* select_pseudofragment_id;
      sqlite3_stmt* select_connection_id;
      sqlite3_prepare_v2(database, "INSERT INTO sources (smiles) VALUES (?)", -1, &insert_source, NULL);
      sqlite3_prepare_v2(database, "INSERT INTO pseudofragments (smiles, ring_part, has_ring, level, pickle, pickle_size, frequency) VALUES (?, ?, ?, ?, ?, ?, 1) ON CONFLICT (smiles, ring_part) DO UPDATE SET frequency=frequency+1", -1, &insert_pseudofragment, NULL);
      sqlite3_prepare_v2(database, "INSERT INTO connections (start_atom_type, end_atom_type, bond_type, string, frequency) VALUES (?, ?, ?, ?, 1) ON CONFLICT (start_atom_type, end_atom_type, bond_type) DO UPDATE SET frequency=frequency+1", -1, &insert_connection, NULL);
      sqlite3_prepare_v2(database, "INSERT INTO sources_pseudofragments (source_id, pseudofragment_id, frequency) VALUES (?, ?, 1) ON CONFLICT (source_id, pseudofragment_id) DO UPDATE SET frequency=frequency+1", -1, &insert_source_pseudofragment_relationship, NULL);
      sqlite3_prepare_v2(database, "INSERT OR IGNORE INTO pseudofragments_connections (pseudofragment_id, connection_id, frequency) VALUES (?, ?, ?)", -1, &insert_pseudofragment_connection_relationship, NULL);
      sqlite3_prepare_v2(database, "SELECT id FROM pseudofragments WHERE smiles = ? AND ring_part = ?", -1, &select_pseudofragment_id, NULL);
      sqlite3_prepare_v2(database, "SELECT id FROM connections WHERE start_atom_type = ? AND end_atom_type = ? AND bond_type = ?", -1, &select_connection_id, NULL);

      const std::pair<unsigned, unsigned> chunk = chunks[thrid];
      unsigned chunk_size = chunk.second - chunk.first;

      // Set up a formatter to convert a Connection into a tuple-like string.
      boost::format connection_formatter("(%d, %d, %d)");

      // Create a local supplier that begins at the start of the file chunk.
      RDKit::SmilesMolSupplier local_supplier(input);
      local_supplier.moveTo(chunk.first);

      // Loop over the molecules in the file chunk.
      for (unsigned i = 0; i <= chunk_size; ++i) {
        RDKit::ROMOL_SPTR mol(local_supplier.next());

        // If the molecule was succesfully read, fragment it.
        if (mol) {
          std::string smiles = MolToSmiles(*mol);

          std::vector<Pseudofragment> pseudofragments;
          try {
            // Fragment the molecule into pseudofragments.
            pseudofragments = FragmentMol(mol, minsize, maxsize, ringfrag, purge, connection_formatter);
          } catch (...) {
            Log(logger, "ERROR", "Fragmentation protocol failed.", smiles);
            ++progress;
            continue;
          };

          #pragma omp critical
          {
            // If the fragmentation didn't yield any fragments skip to the next input molecule.
            if (pseudofragments.empty()) {
              Log(logger, "WARNING", "Fragmentation didn't yield any fragments.", smiles);
            // Otherwise, insert the data into the database.
            } else {
              // (1) Insert the source molecule.
              sqlite3_int64 source_id = InsertSourceMolInDB(smiles, database, insert_source);
              // (2) Insert the Pseudofragments.
              for (auto& pseudofragment : pseudofragments) {
                InsertPseudofragmentInDB(pseudofragment, insert_pseudofragment);
                // (3) Add the source-pseudofragment relationship to the joining table.
                sqlite3_int64 pseudofragment_id = InsertSourcePseudofragmentRelationshipInDB(source_id, pseudofragment, select_pseudofragment_id, insert_source_pseudofragment_relationship);
                // (4) Insert the Pseudofragments Connections and add the pseudofragment-
                //     connection relationships to the joining table.
                InsertConnectionsInDB(pseudofragment, pseudofragment_id, insert_connection, select_connection_id, insert_pseudofragment_connection_relationship);
              };
            };
            // Advance the progress bar.
            ++progress;
          };
        // If the molecule couldn't be read, report the error.
        } else {
          Log(logger, "ERROR", "Failed to read next compound.");
        };
      };
    };

  // If the input file is a SDF file:
  } else if (extension == ".sdf") {
    RDKit::SDMolSupplier supplier(input);
    Log(logger, "INFO", "Created molecule input stream.");

    // Initialize a progress bar.
    unsigned n_mols = supplier.length();
    supplier.reset();
    boost::progress_display progress(n_mols);

    // Divide the input file into N equally sized chunks, where N is the number of
    // OpenMP threads.
    std::vector<std::pair<unsigned, unsigned>> chunks = EquallySizedChunks(n_mols, n_threads);

    #pragma omp parallel for
    for (int thrid = 0; thrid < n_threads; ++thrid) {
      // Prepare the statements required to work with the database.
      sqlite3_stmt* insert_source;
      sqlite3_stmt* insert_pseudofragment;
      sqlite3_stmt* insert_connection;
      sqlite3_stmt* insert_source_pseudofragment_relationship;
      sqlite3_stmt* insert_pseudofragment_connection_relationship;
      sqlite3_stmt* select_pseudofragment_id;
      sqlite3_stmt* select_connection_id;
      sqlite3_prepare_v2(database, "INSERT INTO sources (smiles) VALUES (?)", -1, &insert_source, NULL);
      sqlite3_prepare_v2(database, "INSERT INTO pseudofragments (smiles, ring_part, has_ring, level, pickle, pickle_size, frequency) VALUES (?, ?, ?, ?, ?, ?, 1) ON CONFLICT (smiles, ring_part) DO UPDATE SET frequency=frequency+1", -1, &insert_pseudofragment, NULL);
      sqlite3_prepare_v2(database, "INSERT INTO connections (start_atom_type, end_atom_type, bond_type, string, frequency) VALUES (?, ?, ?, ?, 1) ON CONFLICT (start_atom_type, end_atom_type, bond_type) DO UPDATE SET frequency=frequency+1", -1, &insert_connection, NULL);
      sqlite3_prepare_v2(database, "INSERT INTO sources_pseudofragments (source_id, pseudofragment_id, frequency) VALUES (?, ?, 1) ON CONFLICT (source_id, pseudofragment_id) DO UPDATE SET frequency=frequency+1", -1, &insert_source_pseudofragment_relationship, NULL);
      sqlite3_prepare_v2(database, "INSERT OR IGNORE INTO pseudofragments_connections (pseudofragment_id, connection_id, frequency) VALUES (?, ?, ?)", -1, &insert_pseudofragment_connection_relationship, NULL);
      sqlite3_prepare_v2(database, "SELECT id FROM pseudofragments WHERE smiles = ? AND ring_part = ?", -1, &select_pseudofragment_id, NULL);
      sqlite3_prepare_v2(database, "SELECT id FROM connections WHERE start_atom_type = ? AND end_atom_type = ? AND bond_type = ?", -1, &select_connection_id, NULL);

      const std::pair<unsigned, unsigned> chunk = chunks[thrid];
      unsigned chunk_size = chunk.second - chunk.first;

      // Set up a formatter to convert a Connection into a tuple-like string.
      boost::format connection_formatter("(%d, %d, %d)");

      // Create a local supplier that begins at the start of the file chunk.
      RDKit::SDMolSupplier local_supplier(input);
      local_supplier.moveTo(chunk.first);

      // Loop over the molecules in the file chunk.
      for (unsigned i = 0; i <= chunk_size; ++i) {
        RDKit::ROMOL_SPTR mol (local_supplier.next());

        // If the molecule was succesfully read, fragment it.
        if (mol) {
          std::string smiles = MolToSmiles(*mol);

          std::vector<Pseudofragment> pseudofragments;
          try {
            // Fragment the molecule into pseudofragments.
            pseudofragments = FragmentMol(mol, minsize, maxsize, ringfrag, purge, connection_formatter);
          } catch (...) {
            Log(logger, "ERROR", "Fragmentation protocol failed.", smiles);
            ++progress;
            continue;
          };

          #pragma omp critical
          {
            // If the fragmentation didn't yield any fragments skip to the next input molecule.
            if (pseudofragments.empty()) {
              Log(logger, "WARNING", "Fragmentation didn't yield any fragments.", smiles);
            // Otherwise, insert the data into the database.
            } else {
              // (1) Insert the source molecule.
              sqlite3_int64 source_id = InsertSourceMolInDB(smiles, database, insert_source);
              // (2) Insert the Pseudofragments.
              for (auto& pseudofragment : pseudofragments) {
                InsertPseudofragmentInDB(pseudofragment, insert_pseudofragment);
                // (3) Add the source-pseudofragment relationship to the joining table.
                sqlite3_int64 pseudofragment_id = InsertSourcePseudofragmentRelationshipInDB(source_id, pseudofragment, select_pseudofragment_id, insert_source_pseudofragment_relationship);
                // (4) Insert the Pseudofragments Connections and add the pseudofragment-
                //     connection relationships to the joining table.
                InsertConnectionsInDB(pseudofragment, pseudofragment_id, insert_connection, select_connection_id, insert_pseudofragment_connection_relationship);
              };
            };
            // Advance the progress bar.
            ++progress;
          };
        // If the molecule couldn't be read, report the error.
        } else {
          Log(logger, "ERROR", "Failed to read next compound.");
        };
      };
    };

  // If the molecule input format is something else:
  } else {
    Log(logger, "FATAL", "Unrecognized input format.");
    throw std::runtime_error(std::string("Unrecognized input format. Input should be either SMILES or SDF."));
  };
  Log(logger, "INFO", "Finished input molecules fragmentation.");

  // Commit the transaction. The size of the transaction could potentially be increased.
  Log(logger, "INFO", "Committing database transaction.");
  sqlite3_exec(database, "COMMIT TRANSACTION", NULL, NULL, NULL);
  Log(logger, "INFO", "Successfully committed database transaction.");

  // Gather and store statistics on future query time to aid the query planner.
  Log(logger, "INFO", "Started SQLite3 optimizer protocol.");
  sqlite3_exec(database, "PRAGMA optimize", NULL, NULL, NULL);
  Log(logger, "INFO", "SQLite3 query planner optimization ended.");

  // Vacuum the database to combat memory fragmentation and decrease size overhead.
  Log(logger, "INFO", "Cleaning up the database.");
  sqlite3_exec(database, "VACUUM", NULL, NULL, NULL);
  Log(logger, "INFO", "Succesfully cleaned the database.");

  // Close the database connection.
  sqlite3_close(database);
  Log(logger, "INFO", "Closed database connection.");

  // Close the log.
  Log(logger, "INFO", "Application terminated normally.");
  logger.close();

  // Signal success.
  return 0;
};
