#include "Reconstruction.hpp"

// Definition of the version integer of the ReconstructedMol & associated classes.
extern const unsigned reconstruction_version = 20200610;

// Definition of a variable that keeps track of whether the RDKit::MolPickler was configured.
bool molpickler_configured = false;

// Definition of a function to configure the RDKit::MolPickler to pickle RDKit::Atom properties.
void ConfigureMolpickler() {
  if (molpickler_configured) {
    return;
  };
  RDKit::MolPickler::setDefaultPickleProperties(RDKit::PicklerOps::AtomProps);
  molpickler_configured = true;
};

// Class SizeController
SizeController::SizeController() = default;
SizeController::SizeController(unsigned max_size, float mean, float stdev, float peak) : max_size(max_size), mean(mean), stdev(stdev), peak(peak) {
  assert(peak > 0.0f && peak < 1.0f);
  // Define a normal distribution using the user provided mean and standard
  // deviation. This distribution describes the probability of not changing
  // the object's size.
  mean_squared = mean * mean;
  two_mean = 2.0f * mean;
  two_stdev_squared = 2.0f * stdev * stdev;
  denominator = stdev * sqrt(2.0f * M_PI);
  // Calculate the probability densities of the normal distribution for all
  // integers in the range from 0 to the maximum size.
  float probability_density_sum = 0.0f;
  std::unordered_map<unsigned, float> probability_densities;
  for (unsigned size = 0u; size <= max_size; ++size) {
    float density = NormalProbabilityDensity(size);
    probability_densities.insert({ size, density });
    probability_density_sum += density;
  };
  // Fit a discrete function to the shape of the normal distribution that
  // transforms the probability densities into probabilities of no size change.
  for (const auto& sd : probability_densities) {
    change_probabilities.insert({ sd.first, sd.second / probability_density_sum });
  };
  // Transform the discrete function so that its maximum probability matches the
  // user provided peak probability and it describes the probability of changing size
  float scaling_factor = peak / change_probabilities[static_cast<unsigned>(mean)];
  for (auto& sp : change_probabilities) {
    sp.second = 1.0f - sp.second * scaling_factor;
  };
  // Define a logistic function using as the inflection point the mean
  // and as growth rate the steepness of the normal distribution. This
  // logistic function describes the probability of growth at a given size.
  k = 0.6826895f / 2.0f * stdev;
  // Calculate the growth probabilities at all sizes up to the maximum size.
  // The probability of growing at size 0 is hardcoded to be 1 (i.e. shrinking
  // is not possible).
  growth_probabilities.insert({ 0u, 1.0f });
  for (unsigned size = 1; size <= max_size; ++size) {
    growth_probabilities.insert({size, LogisticValue(size)});
  };
  // Create discrete distributions from the calculated probabilities.
  for (const auto& sp : change_probabilities) {
    unsigned size = sp.first;
    float change_probability = sp.second;
    change_distributions.insert({ size, std::discrete_distribution<unsigned>({ 1.0f - change_probability, change_probability }) });
  };
  for (const auto& sp : growth_probabilities) {
    unsigned size = sp.first;
    float growth_probability = sp.second;
    growth_distributions.insert({ size, std::discrete_distribution<unsigned>({ 1.0f - growth_probability, growth_probability }) });
  };
};

int SizeController::Decide(unsigned size, std::mt19937& prng) {
  if (size > max_size) {
    return -1;
  };
  if (DecideChange(size, prng)) {
    if (DecideGrowth(size, prng)) {
      return 1;
    } else {
      return -1;
    };
  } else {
    return 0;
  };
};

bool SizeController::DecideChange(unsigned size, std::mt19937& prng) {
  if (size > max_size) {
    return true;
  };
  return change_distributions.at(size)(prng);
};

bool SizeController::DecideGrowth(unsigned size, std::mt19937& prng) {
  if (size >= max_size) {
    return false;
  };
  return growth_distributions.at(size)(prng);
};

float SizeController::NormalProbabilityDensity(unsigned size) const {
  float e_exp = -(size * size - size * two_mean + mean_squared) / two_stdev_squared;
  return exp(e_exp) / denominator;
};

float SizeController::LogisticValue(unsigned size) const {
  return 1.0f / (1.0f + exp(k * (size - mean)));
};


// Class EvolutionReport
EvolutionReport::EvolutionReport() = default;
EvolutionReport::EvolutionReport(const std::string& output_directory_path) {
  // Define the paths to the output CSV files.
  std::filesystem::path path (output_directory_path);
  std::filesystem::path scores_csv_path = path;
  std::filesystem::path n_ring_atoms_csv_path = path;
  std::filesystem::path operation_frequencies_csv_path = path;
  scores_csv_path.append("reconstructions_scores.csv");
  n_ring_atoms_csv_path.append("reconstructions_nrings.csv");
  operation_frequencies_csv_path.append("operation_frequencies.csv");
  // Initialize the output streams for the CSV files.
  scores_csv.open(scores_csv_path);
  n_ring_atoms_csv.open(n_ring_atoms_csv_path);
  operation_frequencies_csv.open(operation_frequencies_csv_path);
  // Write the headers of the CSV files.
  scores_csv << "Generation,NScoringCalls,Scores\n";
  n_ring_atoms_csv << "Generation, NRings\n";
  operation_frequencies_csv << ",PeripheralExpansion,InternalExpansion,PeripheralDeletion,InternalDeletion,";
  operation_frequencies_csv << "PeripheralSubstitution,InternalSubstitution,PeripheralCrossover,InternalCrossover,Translation,StereoFlip\n";
};

void EvolutionReport::WriteScores(unsigned generation, unsigned n_scoring_calls, const std::list<ReconstructedMol>& reconstructions) {
  scores_csv << generation << "," << n_scoring_calls;
  for (const ReconstructedMol& reconstruction : reconstructions) {
    scores_csv << "," << reconstruction.GetScore();
  };
  scores_csv << "\n";
};

void EvolutionReport::WriteNRings(unsigned generation, const std::list<ReconstructedMol>& reconstructions) {
  n_ring_atoms_csv << generation;
  for (const ReconstructedMol& reconstruction : reconstructions) {
    n_ring_atoms_csv << "," << reconstruction.GetNRingAtoms();
  };
  n_ring_atoms_csv << "\n";
};

void EvolutionReport::WriteOperationFrequencies() {
  operation_frequencies_csv << "NAttempts,";
  operation_frequencies_csv << attempt_frequencies["peripheral_expansion"] << ",";
  operation_frequencies_csv << attempt_frequencies["internal_expansion"] << ",";
  operation_frequencies_csv << attempt_frequencies["peripheral_deletion"] << ",";
  operation_frequencies_csv << attempt_frequencies["internal_deletion"] << ",";
  operation_frequencies_csv << attempt_frequencies["peripheral_substitution"] << ",";
  operation_frequencies_csv << attempt_frequencies["internal_substitution"] << ",";
  operation_frequencies_csv << attempt_frequencies["peripheral_crossover"] << ",";
  operation_frequencies_csv << attempt_frequencies["internal_crossover"] << ",";
  operation_frequencies_csv << attempt_frequencies["translation"] << ",";
  operation_frequencies_csv << attempt_frequencies["stereo_flip"] << "\n";
  operation_frequencies_csv << "NSuccess,";
  operation_frequencies_csv << success_frequencies["peripheral_expansion"] << ",";
  operation_frequencies_csv << success_frequencies["internal_expansion"] << ",";
  operation_frequencies_csv << success_frequencies["peripheral_deletion"] << ",";
  operation_frequencies_csv << success_frequencies["internal_deletion"] << ",";
  operation_frequencies_csv << success_frequencies["peripheral_substitution"] << ",";
  operation_frequencies_csv << success_frequencies["internal_substitution"] << ",";
  operation_frequencies_csv << success_frequencies["peripheral_crossover"] << ",";
  operation_frequencies_csv << success_frequencies["internal_crossover"] << ",";
  operation_frequencies_csv << success_frequencies["translation"] << ",";
  operation_frequencies_csv << success_frequencies["stereo_flip"] << "\n";
};

void EvolutionReport::WriteWeights(const ReconstructedMol& reconstruction, const std::string& file_name) {
  // Create a CSV file in the output directory.
  std::filesystem::path weights_csv_path = directory;
  weights_csv_path.append(file_name);
  std::ofstream weights_csv (weights_csv_path);

  // Create the header of the CSV file.
  weights_csv << "Connection,Atom_idx\n";

  // Write the IDs and Weights of the Pseudofragments compatible with each ConnectionPoint.
  for (const auto& c : reconstruction.GetConnections()) {
    const Connection& connection = c.first;
    const CONNECTION_POINT_VECTOR& cpoints = c.second;
    for (const ConnectionPoint& cpoint : cpoints) {
      const std::vector<unsigned>& acyclic_compatible_ids = cpoint.GetIDs().GetAcyclicIDs();
      const std::vector<float>& acyclic_compatible_weights = cpoint.GetWeights().GetAcyclicWeights();
      const std::vector<unsigned>& ring_compatible_ids = cpoint.GetIDs().GetRingIDs();
      const std::vector<float>& ring_compatible_weights = cpoint.GetWeights().GetRingWeights();
      assert(acyclic_compatible_ids.size() == acyclic_compatible_weights.size());
      assert(ring_compatible_ids.size() == ring_compatible_weights.size());
      weights_csv << "(" << connection.GetStartAtomType() << ":" << connection.GetEndAtomType() << ":" << connection.GetBondType() << ")";
      weights_csv << "," << cpoint.GetAtomIdx() << ",AcyclicIDs";
      for (unsigned id : acyclic_compatible_ids) {
        weights_csv << "," << id;
      };
      weights_csv << "\n,,AcyclicWeights";
      for (float weight : acyclic_compatible_weights) {
        weights_csv << "," << weight;
      };
      weights_csv << "\n,,RingIDs";
      for (unsigned id : ring_compatible_ids) {
        weights_csv << "," << id;
      };
      weights_csv << "\n,,RingWeights";
      for (float weight : ring_compatible_weights) {
        weights_csv << "," << weight;
      };
      weights_csv << "\n";
    };
  };
  for (const auto& c : reconstruction.GetInternalConnections()) {
    const Connection& connection = c.first;
    const CONNECTION_POINT_VECTOR& cpoints = c.second;
    for (const ConnectionPoint& cpoint : cpoints) {
      const std::vector<unsigned>& acyclic_compatible_ids = cpoint.GetIDs().GetAcyclicIDs();
      const std::vector<float>& acyclic_compatible_weights = cpoint.GetWeights().GetAcyclicWeights();
      const std::vector<unsigned>& ring_compatible_ids = cpoint.GetIDs().GetRingIDs();
      const std::vector<float>& ring_compatible_weights = cpoint.GetWeights().GetRingWeights();
      assert(acyclic_compatible_ids.size() == acyclic_compatible_weights.size());
      assert(ring_compatible_ids.size() == ring_compatible_weights.size());
      weights_csv << "(" << connection.GetStartAtomType() << ":" << connection.GetEndAtomType() << ":" << connection.GetBondType() << ")";
      weights_csv << "," << cpoint.GetAtomIdx() << ",AcyclicIDs";
      for (unsigned id : acyclic_compatible_ids) {
        weights_csv << "," << id;
      };
      weights_csv << "\n,,AcyclicWeights";
      for (float weight : acyclic_compatible_weights) {
        weights_csv << "," << weight;
      };
      weights_csv << "\n,,RingIDs";
      for (unsigned id : ring_compatible_ids) {
        weights_csv << "," << id;
      };
      weights_csv << "\n,,RingWeights";
      for (float weight : ring_compatible_weights) {
        weights_csv << "," << weight;
      };
      weights_csv << "\n";
    };
  };

  // Close the CSV file.
  weights_csv.close();
};

void EvolutionReport::Close() {
  scores_csv.close();
  n_ring_atoms_csv.close();
  operation_frequencies_csv.close();
};


// Class MolBrick
// Default constructor that initializes empty variables.
MolBrick::MolBrick() = default;

// MolBrick constructor that doesn't assign new immutable indices to the Atoms
// in the pseudomol, but instead keeps the ones of the Pseudofragment. Requires
// that the Pseudofragment be flagged before, but the ConnectionsTable should still
// be based on the "original" mutable indices assigned upon Peudofragment construction.
MolBrick::MolBrick(const Pseudofragment& pseudofragment, unsigned pseudofragment_id, unsigned schidx, float wgt, ReconstructedMol* owner) :
  Pseudofragment(pseudofragment), pseudofragment_id(pseudofragment_id), schematic_idx(schidx), weight(wgt), owner(owner) {
  // Check that the Pseudofragment's Atoms have immutable atom indices.
  for (const RDKit::Atom* atom : pseudomol.atoms()) {
    if (!atom->hasProp("ImmutableIdx")) {
      throw std::runtime_error("The indexless MolBrick constructor expects immutable indices to be pre-assigned to the Pseudofragment.");
    } else {
      unsigned immutable_atom_idx = atom->getProp<unsigned>("ImmutableIdx");
      atom_indices.push_back(immutable_atom_idx);
      if (immutable_atom_idx > max_atom_idx) {
        max_atom_idx = immutable_atom_idx;
      };
    };
  };
  // Update the connection table using the immutable atom indices.
  for (auto& c : connections) {
    for (ConnectionPoint& cpoint : c.second) {
      cpoint.atom_idx = pseudomol.getAtomWithIdx(cpoint.atom_idx)->getProp<unsigned>("ImmutableIdx");
    };
  };
};

// Copy-based MolBrick constructor. The Pseudofragment copy constructor
// is called in place of the Pseudofragment base constructor in order
// to copy the Pseudofragments attributes to the MolBrick. This leaves
// a valid copy of the template Pseudofragment.
MolBrick::MolBrick(const Pseudofragment & pseudofragment, unsigned pseudofragment_id, unsigned start_atom_idx, unsigned schidx, float wgt, ReconstructedMol* owner) :
  Pseudofragment(pseudofragment), pseudofragment_id(pseudofragment_id), schematic_idx(schidx), weight(wgt), owner(owner) {
  // Assign immutable atom indices to the Pseudofragment's atoms.
  unsigned immutable_atom_idx = start_atom_idx;
  for (RDKit::ROMol::AtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    (*ai)->setProp<unsigned>("ImmutableIdx", immutable_atom_idx);
    atom_indices.push_back(immutable_atom_idx);
    ++immutable_atom_idx;
  };
  max_atom_idx = immutable_atom_idx - 1;
  // Update the connection table using the immutable atom indices.
  for (auto& c : connections) {
    for (ConnectionPoint& cpoint : c.second) {
      cpoint.atom_idx = pseudomol.getAtomWithIdx(cpoint.atom_idx)->getProp<unsigned>("ImmutableIdx");
    };
  };
};

// Move-based MolBrick constructor. The Pseudofragment's attributes
// are moved to the MolBrick object instead of being copied. This
// is more performant but leaves the Pseudofragment in an undefined
// state.
MolBrick::MolBrick(Pseudofragment & pseudofragment, unsigned pseudofragment_id, unsigned start_atom_idx, unsigned schidx, float wgt, ReconstructedMol * owner) :
  pseudofragment_id(pseudofragment_id), schematic_idx(schidx), weight(wgt), owner(owner) {
  // Move the Pseudofragment's attributes to the MolBrick.
  pseudomol = std::move(pseudofragment.pseudomol);
  connections = std::move(pseudofragment.connections);
  smiles = std::move(pseudofragment.smiles);
  has_ring = pseudofragment.has_ring;
  ring_part = pseudofragment.ring_part;
  level = pseudofragment.level;
  // Assign immutable atom indices to the Pseudofragment's atoms.
  unsigned immutable_atom_idx = start_atom_idx;
  for (RDKit::ROMol::AtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    (*ai)->setProp<unsigned>("ImmutableIdx", immutable_atom_idx);
    atom_indices.push_back(immutable_atom_idx);
    ++immutable_atom_idx;
  };
  max_atom_idx = immutable_atom_idx - 1;
  // Update the connection table using the immutable atom indices.
  for (auto& c : connections) {
    for (ConnectionPoint& cpoint : c.second) {
      cpoint.atom_idx = pseudomol.getAtomWithIdx(cpoint.atom_idx)->getProp<unsigned>("ImmutableIdx");
    };
  };
};

bool MolBrick::operator== (const MolBrick & other) const {
  return (smiles == other.smiles) && (ring_part == other.ring_part) && (owner == other.owner);
};

bool MolBrick::IsPeripheral() const {
  unsigned n_neighbors = 0;
  for (unsigned adjacency_bool : owner->schematic.adjacency[schematic_idx]) {
    if (adjacency_bool == 1) {
      ++n_neighbors;
      if (n_neighbors > 1) {
        return false;
      };
    };
  };
  return true;
};

bool MolBrick::IsInternal() const {
  unsigned n_neighbors = 0;
  for (unsigned adjacency_bool : owner->schematic.adjacency[schematic_idx]) {
    if (adjacency_bool == 1) {
      ++n_neighbors;
      if (n_neighbors > 1) {
        return true;
      };
    };
  };
  return false;
};

std::vector<MolBrick*> MolBrick::GetNeighbors() const {
  std::vector<MolBrick*> neighbors;
  unsigned neighbor_idx = 0;
  for (unsigned adjacency_bool : owner->schematic.adjacency[schematic_idx]) {
    if (adjacency_bool == 1) {
      neighbors.push_back(&owner->bricks[neighbor_idx]);
    };
    ++neighbor_idx;
  };
  return neighbors;
};

ConnectionsTable MolBrick::GetAvailableConnections(std::mt19937& prng, const MolBrick* ghost_brick, bool transfer_weights) const {
  ConnectionsTable available_connections;
  const ConnectionsTable& owner_connections = owner->connections;
  CONNECTIONS_MAP::const_iterator cmap_it, cmap_end_it = owner_connections.end();
  std::vector<const ConnectionPoint*> available_cpoints;
  CONNECTION_POINT_VECTOR::const_iterator cp_it, cp_begin_it, cp_end_it;
  std::vector<unsigned> available_cpoint_indices;
  std::unordered_map<unsigned, unsigned> cpoint_frequency;
  // Loop over the MolBrick's Connections to determine which Connections are free.
  for (const auto& c : connections) {
    const Connection& connection = c.first;
    const std::vector<ConnectionPoint>& brick_cpoints = c.second;
    cmap_it = owner_connections.find(connection);
    // If the Connection is present in the owner ReconstructedMol's peripheral
    // Connections, it's an available Connection.
    if (cmap_it != cmap_end_it) {
      const std::vector<ConnectionPoint>& owner_cpoints = cmap_it->second;
      cp_begin_it = owner_cpoints.begin();
      cp_end_it = owner_cpoints.end();
      available_cpoints.reserve(owner_cpoints.size());
      available_cpoint_indices.reserve(owner_cpoints.size());
      // For each of the MolBrick's Connection's ConnectionPoints, determine if
      // an equivalent ConnectionPoint was found previously. Two ConnectionPoints
      // are considered equivalent when they share the same atom index.
      cpoint_frequency.clear();
      for (const ConnectionPoint& cpoint : brick_cpoints) {
        // If it wasn't found previously count its number of occurrences in the
        // owner ReconstructedMol's Connections and store it.
        if (cpoint_frequency.find(cpoint.atom_idx) == cpoint_frequency.end()) {
          unsigned cnt = std::count(cp_begin_it, cp_end_it, cpoint);
          cpoint_frequency.insert({ cpoint.atom_idx, cnt });
        };
        // While the ConnectionPoint's frequency is positive, incoming equivalent
        // ConnectionPoints are deemed to be available.
        unsigned& frequency = cpoint_frequency[cpoint.atom_idx];
        if (frequency > 0) {
          available_cpoints.push_back(&cpoint);
          available_cpoint_indices.push_back(cpoint.atom_idx);
          --frequency;
        };
      };
      // Add the available ConnectionPoints to the result.
      if (!available_cpoints.empty()) {
        if (transfer_weights) {
          available_connections.AddConnections(connection, available_cpoints);
        } else {
          available_connections.AddConnections(connection, available_cpoint_indices);
          available_cpoint_indices.clear();
        };
        available_cpoints.clear();
      };
    };
  };
  // If a "ghost brick" was provided, the Connection to said brick is considered
  // to be available.
  if (ghost_brick != nullptr) {
    // Assert that the subject and ghost MolBrick are part of the same owner.
    assert(ghost_brick->owner == owner);
    const ReconstructionSchematic& schematic = owner->schematic;
    unsigned ghost_idx = ghost_brick->schematic_idx;
    // Check whether the subject and ghost MolBrick are adjacent.
    if (schematic.adjacency[schematic_idx][ghost_idx] == 1) {
      // If they are retrieve the Connection and add it to the available connections.
      unsigned interactor_idx = schematic.start_atom_idx[schematic_idx][ghost_idx];
      Connection connection(schematic.start_atom_type[schematic_idx][ghost_idx], schematic.end_atom_type[schematic_idx][ghost_idx], schematic.bond_type[schematic_idx][ghost_idx]);
      if (transfer_weights) {
        // If more than one ConnectionPoint share the same atom index choose a random one.
        const ConnectionPoint* interactor_cpoint = owner->internal_connections.GetConnectionPointWithIdx(connection, interactor_idx, prng);
        available_connections.AddConnection(connection, interactor_cpoint);
      } else {
        available_connections.AddConnection(connection, interactor_idx);
      };
    };
  };
  return available_connections;
};

ConnectionsTable MolBrick::GetAvailableConnections(std::mt19937& prng, const std::vector<MolBrick*>& ghost_bricks, bool transfer_weights) const {
  ConnectionsTable available_connections;
  const ConnectionsTable& owner_connections = owner->connections;
  CONNECTIONS_MAP::const_iterator cmap_it, cmap_end_it = owner_connections.end();
  std::vector<const ConnectionPoint*> available_cpoints;
  CONNECTION_POINT_VECTOR::const_iterator cp_it, cp_begin_it, cp_end_it;
  std::vector<unsigned> available_cpoint_indices;
  std::unordered_map<unsigned, unsigned> cpoint_frequency;
  // Loop over the MolBrick's Connections to determine which Connections are free.
  for (const auto& c : connections) {
    const Connection& connection = c.first;
    const std::vector<ConnectionPoint>& brick_cpoints = c.second;
    cmap_it = owner_connections.find(connection);
    // If the Connection is present in the owner ReconstructedMol's peripheral
    // Connections, it's an available Connection.
    if (cmap_it != cmap_end_it) {
      const std::vector<ConnectionPoint>& owner_cpoints = cmap_it->second;
      cp_begin_it = owner_cpoints.begin();
      cp_end_it = owner_cpoints.end();
      available_cpoints.reserve(owner_cpoints.size());
      available_cpoint_indices.reserve(owner_cpoints.size());
      // For each of the MolBrick's Connection's ConnectionPoints, determine if
      // an equivalent ConnectionPoint was found previously. Two ConnectionPoints
      // are considered equivalent when they share the same atom index.
      cpoint_frequency.clear();
      for (const ConnectionPoint& cpoint : brick_cpoints) {
        // If it wasn't found previously count its number of occurrences in the
        // owner ReconstructedMol's Connections and store it.
        if (cpoint_frequency.find(cpoint.atom_idx) == cpoint_frequency.end()) {
          unsigned cnt = std::count(cp_begin_it, cp_end_it, cpoint);
          cpoint_frequency.insert({ cpoint.atom_idx, cnt });
        };
        // While the ConnectionPoint's frequency is positive, incoming equivalent
        // ConnectionPoints are deemed to be available.
        unsigned& frequency = cpoint_frequency[cpoint.atom_idx];
        if (frequency > 0) {
          available_cpoints.push_back(&cpoint);
          available_cpoint_indices.push_back(cpoint.atom_idx);
          --frequency;
        };
      };
      // Add the available ConnectionPoints to the result.
      if (!available_cpoints.empty()) {
        if (transfer_weights) {
          available_connections.AddConnections(connection, available_cpoints);
        } else {
          available_connections.AddConnections(connection, available_cpoint_indices);
          available_cpoint_indices.clear();
        };
        available_cpoints.clear();
      };
    };
  };
  // If a vector of "ghost brick" was provided, the Connection to said bricks are
  // considered to be available.
  unsigned ghost_idx, interactor_idx;
  const ConnectionPoint* interactor_cpoint;
  for (MolBrick* ghost_brick : ghost_bricks) {
    // Assert that the subject and ghost MolBrick are part of the same owner.
    assert(ghost_brick->owner == owner);
    const ReconstructionSchematic& schematic = owner->schematic;
    ghost_idx = ghost_brick->schematic_idx;
    // Check whether the subject and ghost MolBrick are adjacent.
    if (schematic.adjacency[schematic_idx][ghost_idx] == 1) {
      // If they are retrieve the Connection and add it to the available connections.
      interactor_idx = schematic.start_atom_idx[schematic_idx][ghost_idx];
      Connection connection(schematic.start_atom_type[schematic_idx][ghost_idx], schematic.end_atom_type[schematic_idx][ghost_idx], schematic.bond_type[schematic_idx][ghost_idx]);
      if (transfer_weights) {
        // If more than one ConnectionPoint share the same atom index choose a random one.
        interactor_cpoint = owner->internal_connections.GetConnectionPointWithIdx(connection, interactor_idx, prng);
        available_connections.AddConnection(connection, interactor_cpoint);
      } else {
        available_connections.AddConnection(connection, interactor_idx);
      };
    };
  };
  return available_connections;
};

void MolBrick::UpdateImmutableIndices(unsigned start_atom_idx) {
  std::unordered_map<unsigned, unsigned> atom_indices_map;
  atom_indices.clear();
  unsigned immutable_atom_idx = start_atom_idx;
  for (RDKit::ROMol::AtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    atom_indices_map.insert({ (*ai)->getProp<unsigned>("ImmutableIdx"), immutable_atom_idx });
    (*ai)->setProp<unsigned>("ImmutableIdx", immutable_atom_idx);
    atom_indices.push_back(immutable_atom_idx);
    ++immutable_atom_idx;
  };
  max_atom_idx = immutable_atom_idx - 1;

  // Update the connection table using the immutable atom indices.
  // This approach might also be better for the MolBrick constructor.
  for (auto& c : connections) {
    for (ConnectionPoint& cpoint : c.second) {
      cpoint.atom_idx = atom_indices_map[cpoint.atom_idx];
    };
  };
};

const ReconstructedMol* MolBrick::GetOwner() const {
  return owner;
};

unsigned MolBrick::GetFragmentID() const {
  return pseudofragment_id;
};

unsigned MolBrick::GetMaxAtomIdx() const {
  return max_atom_idx;
};

unsigned MolBrick::GetSchematicIdx() const {
  return schematic_idx;
};

const std::vector<unsigned>& MolBrick::GetImmutableAtomIndices() const {
  return atom_indices;
};

float MolBrick::GetWeight() const {
  return weight;
};


// Class ReconstructionSchematic
ReconstructionSchematic::ReconstructionSchematic() = default;
ReconstructionSchematic::ReconstructionSchematic(unsigned size) {
  this->Expand(size);
};

void ReconstructionSchematic::Expand(unsigned bump) {
  unsigned new_dimension = dimension + bump;
  adjacency.resize(new_dimension);
  for (std::vector<unsigned>& column : adjacency) {
    column.resize(new_dimension);
  };
  start_atom_idx.resize(new_dimension);
  for (std::vector<unsigned>& column : start_atom_idx) {
    column.resize(new_dimension);
  };
  end_atom_idx.resize(new_dimension);
  for (std::vector<unsigned>& column : end_atom_idx) {
    column.resize(new_dimension);
  };
  start_atom_type.resize(new_dimension);
  for (std::vector<unsigned>& column : start_atom_type) {
    column.resize(new_dimension);
  };
  end_atom_type.resize(new_dimension);
  for (std::vector<unsigned>& column : end_atom_type) {
    column.resize(new_dimension);
  };
  bond_type.resize(new_dimension);
  for (std::vector<unsigned>& column : bond_type) {
    column.resize(new_dimension);
  };
  dimension += bump;
};

void ReconstructionSchematic::UnitSchematic() {
  assert(dimension == 0);
  this->Expand();
  adjacency[0][0] = 0;
  start_atom_idx[0][0] = 0;
  end_atom_idx[0][0] = 0;
  start_atom_type[0][0] = 0;
  end_atom_type[0][0] = 0;
  bond_type[0][0] = 0;
  ++dimension;
};

void ReconstructionSchematic::ClearCell(unsigned row_idx, unsigned column_idx) {
  adjacency[row_idx][column_idx] = 0;
  start_atom_idx[row_idx][column_idx] = 0;
  end_atom_idx[row_idx][column_idx] = 0;
  start_atom_type[row_idx][column_idx] = 0;
  end_atom_type[row_idx][column_idx] = 0;
  bond_type[row_idx][column_idx] = 0;
};

void ReconstructionSchematic::ClearRow(unsigned row_idx) {
  for (unsigned column_idx = 0; column_idx < adjacency[row_idx].size(); ++column_idx) {
    this->ClearCell(row_idx, column_idx);
  };
};

void ReconstructionSchematic::ClearColumn(unsigned column_idx) {
  for (unsigned row_idx = 0; row_idx < adjacency[column_idx].size(); ++row_idx) {
    this->ClearCell(row_idx, column_idx);
  };
};

unsigned ReconstructionSchematic::GetDimension() const {
  return dimension;
};

const std::vector<std::vector<unsigned>>& ReconstructionSchematic::GetAdjacencyMatrix() const {
  return adjacency;
};

const std::vector<std::vector<unsigned>>& ReconstructionSchematic::GetStartAtomIdxMatrix() const {
  return start_atom_idx;
};

const std::vector<std::vector<unsigned>>& ReconstructionSchematic::GetEndAtomIdxMatrix() const {
  return end_atom_idx;
};

const std::vector<std::vector<unsigned>>& ReconstructionSchematic::GetStartAtomTypeMatrix() const {
  return start_atom_type;
};

const std::vector<std::vector<unsigned>>& ReconstructionSchematic::GetEndAtomTypeMatrix() const {
  return end_atom_type;
};

const std::vector<std::vector<unsigned>>& ReconstructionSchematic::GetBondTypeMatrix() const {
  return bond_type;
};

void ReconstructionSchematic::PrintAdjacencyMatrix() const {
  for (const std::vector<unsigned>& column : adjacency) {
    for (unsigned cell : column) {
      std::cout << cell << " ";
    };
    std::cout << std::endl;
  };
};

void ReconstructionSchematic::PrintStartAtomIdxMatrix() const {
  for (const std::vector<unsigned>& column : start_atom_idx) {
    for (unsigned cell : column) {
      std::cout << cell << " ";
    };
    std::cout << std::endl;
  };
};

void ReconstructionSchematic::PrintEndAtomIdxMatrix() const {
  for (const std::vector<unsigned>& column : end_atom_idx) {
    for (unsigned cell : column) {
      std::cout << cell << " ";
    };
    std::cout << std::endl;
  };
};

void ReconstructionSchematic::PrintStartAtomTypeMatrix() const {
  for (const std::vector<unsigned>& column : start_atom_type) {
    for (unsigned cell : column) {
      std::cout << cell << " ";
    };
    std::cout << std::endl;
  };
};

void ReconstructionSchematic::PrintEndAtomTypeMatrix() const {
  for (const std::vector<unsigned>& column : end_atom_type) {
    for (unsigned cell : column) {
      std::cout << cell << " ";
    };
    std::cout << std::endl;
  };
};

void ReconstructionSchematic::PrintBondTypeMatrix() const {
  for (const std::vector<unsigned>& column : bond_type) {
    for (unsigned cell : column) {
      std::cout << cell << " ";
    };
    std::cout << std::endl;
  };
};

// Class ReconstructionsInventory
ReconstructionsInventory::ReconstructionsInventory() = default;
ReconstructionsInventory::ReconstructionsInventory(float gamma) : gamma(gamma) {};
ReconstructionsInventory::ReconstructionsInventory(const std::list<ReconstructedMol> & reconstructions, float gamma) : gamma(gamma) {
  for (const ReconstructedMol& reconstruction : reconstructions) {
    Add(reconstruction);
  };
};

void ReconstructionsInventory::Add(const MolBrick & brick) {
  if (brick.has_ring) {
    ring.push_back(brick);
  } else if (brick.ring_part) {
    ring_part.push_back(brick);
  } else {
    acyclic.push_back(brick);
  };
};

void ReconstructionsInventory::Add(const ReconstructedMol & reconstruction) {
  for (const auto& brick : reconstruction.bricks) {
    Add(brick.second);
  };
};

void ReconstructionsInventory::Add(const ReconstructedMol * reconstruction) {
  for (const auto& brick : reconstruction->bricks) {
    Add(brick.second);
  };
};

void ReconstructionsInventory::Remove(const ReconstructedMol & reconstruction) {
  for (const auto& brick : reconstruction.bricks) {
    Remove(brick.second);
  };
};

void ReconstructionsInventory::Remove(const ReconstructedMol * reconstruction) {
  for (const auto& brick : reconstruction->bricks) {
    Remove(brick.second);
  };
};

void ReconstructionsInventory::Remove(const MolBrick & brick) {
  std::vector<MolBrick>::iterator it, end_it;
  if (brick.has_ring) {
    end_it = ring.end();
    it = std::find(ring.begin(), end_it, brick);
    assert(it != end_it);
    if (it != end_it - 1) {
      *it = std::move(ring.back());
    };
    ring.pop_back();
  } else if (brick.ring_part) {
    end_it = ring_part.end();
    it = std::find(ring_part.begin(), end_it, brick);
    assert(it != end_it);
    if (it != end_it - 1) {
      *it = std::move(ring_part.back());
    };
    ring_part.pop_back();
  } else {
    end_it = acyclic.end();
    it = std::find(acyclic.begin(), end_it, brick);
    assert(it != end_it);
    if (it != end_it - 1) {
      *it = std::move(acyclic.back());
    };
    acyclic.pop_back();
  };
};

void ReconstructionsInventory::SubmitToAdditionQueue(const MolBrick * brick) {
  addition_queue.push_back(*brick);
};

void ReconstructionsInventory::SubmitToAdditionQueue(const MolBrick & brick) {
  addition_queue.push_back(brick);
};

void ReconstructionsInventory::SubmitToDeletionQueue(const MolBrick * brick) {
  deletion_queue.push_back(*brick);
};

void ReconstructionsInventory::SubmitToDeletionQueue(const MolBrick & brick) {
  deletion_queue.push_back(brick);
};

void ReconstructionsInventory::ClearQueues() {
  addition_queue.clear();
  deletion_queue.clear();
};

void ReconstructionsInventory::CommitQueuedOperations() {
  // Add all the molecules in the addition queue.
  for (const MolBrick& brick : addition_queue) {
    Add(brick);
  };
  // Remove all the molecules in the deletion queue.
  for (const MolBrick& brick : deletion_queue) {
    Remove(brick);
  };
  ClearQueues();
  assert(addition_queue.empty());
  assert(deletion_queue.empty());
};

void ReconstructionsInventory::CalcAcyclicWeights() {
  acyclic_weights.clear();
  acyclic_weights.reserve(acyclic.size());
  for (const MolBrick& brick : acyclic) {
    acyclic_weights.push_back(brick.weight * std::pow(brick.owner->score, gamma));
  };
};

void ReconstructionsInventory::CalcRingWeights() {
  ring_weights.clear();
  ring_weights.reserve(ring.size());
  for (const MolBrick& brick : ring) {
    ring_weights.push_back(brick.weight * std::pow(brick.owner->score, gamma));
  };
};

void ReconstructionsInventory::CalcRingPartWeights() {
  ring_part_weights.clear();
  ring_part_weights.reserve(ring_part.size());
  for (const MolBrick& brick : ring_part) {
    ring_part_weights.push_back(brick.weight * std::pow(brick.owner->score, gamma));
  };
};

void ReconstructionsInventory::CalcWeights() {
  CalcRingWeights();
  CalcRingPartWeights();
  CalcAcyclicWeights();
};

void ReconstructionsInventory::SetGamma(float new_gamma) {
  gamma = new_gamma;
};

void ReconstructionsInventory::Clear() {
  acyclic.clear();
  ring.clear();
  ring_part.clear();
  acyclic_weights.clear();
  ring_weights.clear();
  ring_part_weights.clear();
  ClearQueues();
};

const std::vector<MolBrick>& ReconstructionsInventory::GetRingBricks() const {
  return ring;
};

const std::vector<MolBrick>& ReconstructionsInventory::GetRingPartBricks() const {
  return ring_part;
};

const std::vector<MolBrick>& ReconstructionsInventory::GetAcyclicBricks() const {
  return acyclic;
};

const std::vector<float>& ReconstructionsInventory::GetRingBrickWeights() const {
  return ring_weights;
};

const std::vector<float>& ReconstructionsInventory::GetRingPartBrickWeights() const {
  return ring_part_weights;
};

const std::vector<float>& ReconstructionsInventory::GetAcyclicBrickWeights() const {
  return acyclic_weights;
};


// Class DeletionPoint
DeletionPoint::DeletionPoint(MolBrick * deletable) :
  deletable(deletable),
  has_ring(deletable->has_ring),
  ring_part(deletable->ring_part),
  weight(deletable->weight) {};

void DeletionPoint::PerceiveNeighborsAvailableConnections(std::mt19937& prng) {
  assert(neighbors_available_connections.empty());
  neighbors = deletable->GetNeighbors();
  assert(neighbors.size() > 0);
  for (MolBrick* neighbor : neighbors) {
    neighbors_available_connections.insert({ neighbor, neighbor->GetAvailableConnections(prng, deletable) });
  };
  n_neighbors = neighbors.size();
  n_other_neighbors = n_neighbors - 1;
  evaluated_neighbors.reserve(n_neighbors);
};

void DeletionPoint::Evaluate(ConnectionCompatibilities& compatibilities, std::mt19937& prng) {
  assert(!evaluated);
  assert(!generated_all_matchings);
  if (neighbors_available_connections.empty()) {
    PerceiveNeighborsAvailableConnections(prng);
  };
  // Loop over the neighboring MolBricks. In each iteration the current MolBrick
  // is considered the new "core" that will move into the position of the deleted
  // MolBrick.
  std::vector<MolBrick*> other_neighbors;
  std::vector<MolBrick*>::iterator it;
  std::map<int, const MolBrick*> vtx_ids_neighbors;
  std::map<const MolBrick*, int> neighbors_vtx_ids;
  std::map<int, const Connection*> vtx_ids_connections;
  std::vector<std::pair<int, int>> edges;
  BipartiteGraph graph;
  PATH matching;
  for (MolBrick* core : neighbors) {
    // Mark the core MolBrick as evaluated.
    evaluated_neighbors.push_back(core);
    // Create a copy of the neighbors list excluding the core MolBrick.
    other_neighbors = neighbors;
    it = std::find(other_neighbors.begin(), other_neighbors.end(), core);
    *it = std::move(other_neighbors.back());
    other_neighbors.pop_back();
    // Setup the BipartiteGraph to run a Maximum Bipartite Matching to see if the
    // core MolBrick can connect to all other neighboring MolBricks.
    vtx_ids_neighbors.clear();
    neighbors_vtx_ids.clear();
    vtx_ids_connections.clear();
    edges.clear();
    int id = 1;
    for (MolBrick* other_neighbor : other_neighbors) {
      vtx_ids_neighbors.insert({ id, other_neighbor });
      neighbors_vtx_ids.insert({ other_neighbor, id });
      ++id;
    };
    for (const auto& c : neighbors_available_connections[core]) {
      for (const auto& cp : c.second) {
        vtx_ids_connections.insert({ id, &c.first });
        for (MolBrick* other_neighbor : other_neighbors) {
          if (neighbors_available_connections[other_neighbor].IsCompatibleWith(c.first, compatibilities)) {
            edges.push_back(std::make_pair(neighbors_vtx_ids[other_neighbor], id));
          };
        };
        ++id;
      };
    };
    // Run the Hopcroft-Karp algorithm to find a random solution to the Maximum
    // Bipartite Matching problem.
    if (vtx_ids_connections.size() < n_other_neighbors || edges.size() < n_other_neighbors) {
      continue;
    };
    graph = BipartiteGraph(n_other_neighbors, vtx_ids_connections.size(), edges);
    matching = graph.HopcroftKarp(prng);
    // If the cardinality of the matching is equal to the number of involved
    // neighbors, the DeletionPoint's MolBrick is deletable by connecting the
    // core neighboring MolBrick to the rest of neighbors. Store the matching.
    if (matching.size() == n_other_neighbors) {
      other_neighbors_matchings.insert({ core, std::map<const MolBrick*, const Connection*>{} });
      std::map<const MolBrick*, const Connection*>& core_matching = other_neighbors_matchings[core];
      for (Edge* edge : matching) {
        const MolBrick* match_brick = vtx_ids_neighbors[edge->GetStart()->GetID()];
        const Connection* match_connection = vtx_ids_connections[edge->GetEnd()->GetID()];
        core_matching.insert({ match_brick, match_connection });
      };
      is_deletable = true;
      break;
    };
  };
  evaluated = true;
  if (evaluated_neighbors.size() == n_neighbors) {
    generated_all_matchings = true;
  };
};

void DeletionPoint::GenerateMatchings(ConnectionCompatibilities& compatibilities, std::mt19937& prng) {
  if (generated_all_matchings) {
    return;
  };
  if (neighbors_available_connections.empty()) {
    PerceiveNeighborsAvailableConnections(prng);
  };
  // Loop over the neighboring MolBricks that haven't been evaluated yet.
  std::vector<MolBrick*> other_neighbors;
  std::vector<MolBrick*>::iterator it, begin_eval_it = evaluated_neighbors.begin(), end_eval_it = evaluated_neighbors.end();
  std::map<int, const MolBrick*> vtx_ids_neighbors;
  std::map<const MolBrick*, int> neighbors_vtx_ids;
  std::map<int, const Connection*> vtx_ids_connections;
  std::vector<std::pair<int, int>> edges;
  BipartiteGraph graph;
  PATH matching;
  for (MolBrick* core : neighbors) {
    if (!evaluated || (std::find(begin_eval_it, end_eval_it, core) == end_eval_it)) {
      // Mark the core MolBrick as evaluated.
      evaluated_neighbors.push_back(core);
      end_eval_it = evaluated_neighbors.end();
      // Create a copy of the neighbors list excluding the core MolBrick.
      other_neighbors = neighbors;
      it = std::find(other_neighbors.begin(), other_neighbors.end(), core);
      *it = std::move(other_neighbors.back());
      other_neighbors.pop_back();
      // Setup the BipartiteGraph to run a Maximum Bipartite Matching to see if the
      // core MolBrick can connect to all other neighboring MolBricks.
      vtx_ids_neighbors.clear();
      neighbors_vtx_ids.clear();
      vtx_ids_connections.clear();
      edges.clear();
      int id = 1;
      for (MolBrick* other_neighbor : other_neighbors) {
        vtx_ids_neighbors.insert({ id, other_neighbor });
        neighbors_vtx_ids.insert({ other_neighbor, id });
        ++id;
      };
      for (const auto& c : neighbors_available_connections[core]) {
        for (const auto& cp : c.second) {
          vtx_ids_connections.insert({ id, &c.first });
          for (MolBrick* other_neighbor : other_neighbors) {
            if (neighbors_available_connections[other_neighbor].IsCompatibleWith(c.first, compatibilities)) {
              edges.push_back(std::make_pair(neighbors_vtx_ids[other_neighbor], id));
            };
          };
          ++id;
        };
      };
      // Run the Hopcroft-Karp algorithm to find a random solution to the Maximum
      // Bipartite Matching problem.
      if (vtx_ids_connections.size() < n_other_neighbors || edges.size() < n_other_neighbors) {
        continue;
      };
      graph = BipartiteGraph(n_other_neighbors, vtx_ids_connections.size(), edges);
      matching = graph.HopcroftKarp(prng);
      // If the cardinality of the matching is equal to the number of involved
      // neighbors, the DeletionPoint's MolBrick is deletable by connecting the
      // core neighboring MolBrick to the rest of neighbors. Store the matching.
      if (matching.size() == n_other_neighbors) {
        // Create an empty matching entry for the core MolBrick.
        other_neighbors_matchings.insert({ core, std::map<const MolBrick*, const Connection*>{} });
        std::map<const MolBrick*, const Connection*>& core_matching = other_neighbors_matchings[core];
        for (Edge* edge : matching) {
          const MolBrick* match_brick = vtx_ids_neighbors[edge->GetStart()->GetID()];
          const Connection* match_connection = vtx_ids_connections[edge->GetEnd()->GetID()];
          core_matching.insert({ match_brick, match_connection });
        };
      };
    };
  };
  generated_all_matchings = true;
};

// Class SubstitutionPoint
SubstitutionPoint::SubstitutionPoint(MolBrick* substitutable, bool guided) :
  substitutable(substitutable),
  has_ring(substitutable->has_ring),
  ring_part(substitutable->ring_part),
  guided(guided) {};

void SubstitutionPoint::PerceiveNeighborsAvailableConnections(std::mt19937& prng) {
  assert(neighbors_available_connections.empty());
  std::vector<MolBrick*> neighbors = substitutable->GetNeighbors();
  for (MolBrick* neighbor : neighbors) {
    neighbors_available_connections.insert({ neighbor, neighbor->GetAvailableConnections(prng, substitutable) });
  };
  n_neighbors = neighbors_available_connections.size();
};

void SubstitutionPoint::GenerateCombinations(std::mt19937& prng) {
  assert(combinations.empty());
  // If the neighboring MolBrick's available Connections weren't
  // perceived yet, do it now.
  if (neighbors_available_connections.empty()) {
    PerceiveNeighborsAvailableConnections(prng);
  };
  size_t n = neighbors_available_connections.size();
  assert(n >= 1);
  neighbors_combination.reserve(n);
  // Create a set of iterators to traverse the available
  // connections of the neighboring MolBricks.
  std::vector<CONNECTIONS_MAP::const_iterator> its, begin_its, end_its;
  its.reserve(n);
  begin_its.reserve(n);
  end_its.reserve(n);
  for (const auto& nac : neighbors_available_connections) {
    neighbors_combination.push_back(nac.first);
    const ConnectionsTable& available_connections = nac.second;
    its.push_back(available_connections.begin());
    begin_its.push_back(available_connections.begin());
    end_its.push_back(available_connections.end());
  };
  // Generate all unique neighbor available Connection
  // combinations using an "odometer-like" approach.
  CONNECTIONS_COMBINATION combination;
  while (its[0] != end_its[0]) {
    // Retrieve the current combination of the odometer.
    for (auto& it : its) {
      combination.push_back(it->first);
    };
    // Store the combination.
    combinations.insert(combination);
    combination.clear();
    // Increment the odometer by 1 (i.e. increment the iterators).
    ++its[n - 1];
    for (size_t i = n - 1; (i > 0) && (its[i] == end_its[i]); --i) {
      its[i] = begin_its[i];
      ++its[i - 1];
    };
  };
};

void SubstitutionPoint::PrepareHopcroftKarp(std::mt19937& prng) {
  assert(max_vtx_id == 0);
  assert(vtx_ids_neighbors.empty());
  assert(neighbors_vtx_ids.empty());
  // Assign IDs to the Vertices that will represent the neighboring
  // MolBricks in the BipartiteGraph.
  if (neighbors_available_connections.empty()) {
    PerceiveNeighborsAvailableConnections(prng);
  };
  max_vtx_id = 1;
  for (const auto& nac : neighbors_available_connections) {
    vtx_ids_neighbors.insert({ max_vtx_id, nac.first });
    neighbors_vtx_ids.insert({ nac.first, max_vtx_id });
    ++max_vtx_id;
  };
  // Flag the SubstitutionPoint.
  hopcroft_ready = true;
};

PATH SubstitutionPoint::HopcroftKarp(const MolBrick & substitute, ConnectionCompatibilities & compatibilities, std::mt19937 & prng) {
  // IMPROVE: Better BipartiteGraph construction.
  // Assign IDs to the Vertices that will represent the neighbor's
  // compatible Connections in the BipartiteGraph.
  assert(hopcroft_ready);
  vtx_ids_connections.clear();
  edges.clear();
  int vtx_id = max_vtx_id;
  PATH matching;
  // IMPROVE: Might not be worth to track this
  std::unordered_set<MolBrick*> used_bricks;
  for (const auto& c : substitute.connections) {
    for (const auto& cp : c.second) {
      vtx_ids_connections.insert({ vtx_id, &c.first });
      for (const auto& nac : neighbors_available_connections) {
        if (nac.second.IsCompatibleWith(c.first, compatibilities)) {
          used_bricks.insert(nac.first);
          edges.push_back(std::make_pair(neighbors_vtx_ids[nac.first], vtx_id));
        };
      };
      ++vtx_id;
    };
  };

  // If one of the neighboring MolBricks is not compatible with the
  // putative substitute MolBrick it isn't a valid substitute.
  if (used_bricks.size() < n_neighbors) {
    return matching;
  };

  // Run the Hopcroft-Karp algorithm to find a solution to the Maximum
  // Bipartite Matching problem.
  graph = BipartiteGraph(vtx_ids_neighbors.size(), vtx_ids_connections.size(), edges);
  matching = graph.HopcroftKarp(prng);
  return matching;
};

void SubstitutionPoint::EvaluateSubstitutability(ConnectionQueryResults& query_results, std::mt19937& prng) {
  // Evaluate how a Substitution operation could modulate the number of rings.
  assert(!evaluated);

  // If the Connection combinations weren't generated yet, do it now.
  if (combinations.empty()) {
    GenerateCombinations(prng);
  };

  // Retrieve the acyclic and ring substitute Pseudofragments of the
  // MolBrick taking into account all possible Connection combinations.
  assert(!retrieved_acyclic_substitutes);
  assert(!retrieved_ring_substitutes);
  std::pair<std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::iterator, bool> insert_p;
  // For each Connections combination.
  unsigned n_evaluated_combinations = 0;
  for (const CONNECTIONS_COMBINATION& combination : combinations) {
    ++n_evaluated_combinations;
    // Attempt to retrieve acyclic substitutes.
    insert_p = acyclic_substitutes.insert({ combination, query_results.SlidingIntersection(combination, false, false) });
    // Check if acyclic substitutes were found.
    if (!insert_p.first->second.first.empty()) {
      has_substitutes = true;
      has_acyclic_substitutes = true;
    };
    // Attempt to retrieve ring substitutes.
    insert_p = ring_substitutes.insert({ combination, query_results.SlidingIntersection(combination, true, false) });
    // Check if ring substitutes were found.
    if (!insert_p.first->second.first.empty()) {
      has_substitutes = true;
      has_ring_substitutes = true;
    };
    if (has_acyclic_substitutes && has_ring_substitutes) {
      break;
    };
  };
  if (combinations.size() == n_evaluated_combinations) {
    retrieved_acyclic_substitutes = true;
    retrieved_ring_substitutes = true;
  };

  // Check how the number of rings could be modulated during
  // a Substitution.
  if (has_substitutes) {
    if (substitutable->has_ring) {
      if (has_acyclic_substitutes) {
        rings_can_be_removed = true;
      } else {
        rings_can_be_kept_constant = true;
      };
    } else {
      if (has_ring_substitutes) {
        rings_can_be_added = true;
      } else {
        rings_can_be_kept_constant = true;
      };
    };
  };
  // Flag the SubstitutionPoint as evaluated.
  evaluated = true;
};

void SubstitutionPoint::EvaluateSubstitutability(ReconstructionsInventory & inventory, ConnectionCompatibilities & compatibilities, std::mt19937 & prng) {
  assert(!evaluated);
  assert(hopcroft_ready);
  PATH matching;
  // If the SubstitutionPoint pertains to a ring MolBrick loop over all the acyclic
  // MolBricks in the ReconstructionsInventory and check if any of them are a valid
  // substitute.
  if (has_ring) {
    for (const MolBrick& brick : inventory.acyclic) {
      matching = HopcroftKarp(brick, compatibilities, prng);
      if (matching.size() == n_neighbors) {
        rings_can_be_removed = true;
        break;
      };
    };
  // If the SubstitutionPoint pertains to an acyclic MolBrick loop over all the ring
  // MolBricks in the ReconstructionsInventory and check if any of them are a valid
  // substitute.
  } else {
    for (const MolBrick& brick : inventory.ring) {
      matching = HopcroftKarp(brick, compatibilities, prng);
      if (matching.size() == n_neighbors) {
        rings_can_be_added = true;
        break;
      };
    };
  };
  // Flag the SubstitutionPoint.
  evaluated = true;
};

void SubstitutionPoint::RetrieveAcyclicSubstitutes(ConnectionQueryResults& query_results, std::mt19937& prng) {
  if (retrieved_acyclic_substitutes) {
    return;
  };
  // If the Connection combinations weren't generated yet, do it now.
  if (combinations.empty()) {
    GenerateCombinations(prng);
  };
  // Loop over the Connection combinations.
  std::pair<std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::iterator, bool> insert_p;
  for (const CONNECTIONS_COMBINATION& combination : combinations) {
    // If substitutes weren't retrieved for this combination,
    // attempt to retrieve them.
    if (acyclic_substitutes.find(combination) == acyclic_substitutes.end()) {
      insert_p = acyclic_substitutes.insert({ combination, query_results.SlidingIntersection(combination, false, false) });
      // Check if acyclic substitutes were found.
      if (!insert_p.first->second.first.empty()) {
        has_substitutes = true;
        has_acyclic_substitutes = true;
      };
    };
  };
  // Flag the SubstitutionPoint.
  retrieved_acyclic_substitutes = true;
};

void SubstitutionPoint::RetrieveRingSubstitutes(ConnectionQueryResults& query_results, std::mt19937& prng) {
  if (retrieved_ring_substitutes) {
    return;
  };
  // If the Connection combinations weren't generated yet, do it now.
  if (combinations.empty()) {
    GenerateCombinations(prng);
  };
  // Loop over the Connection combinations.
  std::pair<std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::iterator, bool> insert_p;
  for (const CONNECTIONS_COMBINATION& combination : combinations) {
    // If substitutes weren't retrieved for this combination,
    // attempt to retrieve them.
    if (ring_substitutes.find(combination) == ring_substitutes.end()) {
      insert_p = ring_substitutes.insert({ combination, query_results.SlidingIntersection(combination, true, false) });
      // Check if ring substitutes were found.
      if (!insert_p.first->second.first.empty()) {
        has_substitutes = true;
        has_ring_substitutes = true;
      };
    };
  };
  // Flag the SubstitutionPoint.
  retrieved_ring_substitutes = true;
};

void SubstitutionPoint::RetrieveAcyclicSubstitutes(ReconstructionsInventory& inventory, ConnectionCompatibilities& compatibilities, std::mt19937& prng) {
  assert(hopcroft_ready);
  assert(!retrieved_acyclic_substitutes);
  PATH matching;
  inventory_acyclic_substitutes.first.reserve(inventory.acyclic.size());
  inventory_acyclic_substitutes.second.reserve(inventory.acyclic.size());
  MolBrick* brick;
  for (unsigned i = 0; i < inventory.acyclic.size(); ++i) {
    brick = &inventory.acyclic[i];
    matching = HopcroftKarp(*brick, compatibilities, prng);
    if (matching.size() == n_neighbors) {
      inventory_acyclic_substitutes.first.push_back(brick);
      inventory_acyclic_substitutes.second.push_back(inventory.acyclic_weights[i]);
    };
  }
  retrieved_acyclic_substitutes = true;
  inventory_acyclic_substitutes.first.shrink_to_fit();
  inventory_acyclic_substitutes.second.shrink_to_fit();
  if (!inventory_acyclic_substitutes.first.empty()) {
    has_acyclic_substitutes = true;
  };
};

void SubstitutionPoint::RetrieveRingSubstitutes(ReconstructionsInventory& inventory, ConnectionCompatibilities& compatibilities, std::mt19937& prng) {
  assert(hopcroft_ready);
  assert(!retrieved_ring_substitutes);
  PATH matching;
  inventory_ring_substitutes.first.reserve(inventory.ring.size());
  inventory_ring_substitutes.second.reserve(inventory.ring.size());
  MolBrick* brick;
  for (unsigned i = 0; i < inventory.ring.size(); ++i) {
    brick = &inventory.ring[i];
    matching = HopcroftKarp(*brick, compatibilities, prng);
    if (matching.size() == n_neighbors) {
      inventory_ring_substitutes.first.push_back(brick);
      inventory_ring_substitutes.second.push_back(inventory.ring_weights[i]);
    };
  }
  retrieved_ring_substitutes = true;
  inventory_ring_substitutes.first.shrink_to_fit();
  inventory_ring_substitutes.second.shrink_to_fit();
  if (!inventory_ring_substitutes.first.empty()) {
    has_ring_substitutes = true;
  };
};

void SubstitutionPoint::CalcAcyclicBasedWeight() {
  assert(!calculated_weights);
  assert(retrieved_acyclic_substitutes);
  combination_weights.reserve(combinations.size());

  // If guided evolution is enabled the weights of the retrieved query
  // results have to be modified using the weights of the ConnectionPoints
  // involved in each Connections combination.
  if (guided) {
    // Create references to the owner ReconstructedMol's ConnectionTables.
    ConnectionsTable& peripheral_connections = substitutable->owner->connections;
    ConnectionsTable& internal_connections = substitutable->owner->internal_connections;
    std::vector<ConnectionPoint*> cpoints, all_cpoints;
    // Loop over the Connection combinations and the resulting query results.
    for (auto& s : acyclic_substitutes) {
      const CONNECTIONS_COMBINATION& combination = s.first;
      QUERY_RESULT& qr = s.second;
      const std::vector<unsigned>& qr_ids = qr.first;
      std::vector<float>& qr_weights = qr.second;
      // Retrieve the ConnectionPoints potentially involved in the combination.
      for (unsigned cidx = 0; cidx < n_neighbors; ++cidx) {
        cpoints = peripheral_connections.GetConnectionPoints(combination[cidx], neighbors_combination[cidx]->connections[combination[cidx]]);
        if (!cpoints.empty()) {
          all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
        };
        cpoints = internal_connections.GetConnectionPoints(combination[cidx], neighbors_combination[cidx]->connections[combination[cidx]]);
        if (!cpoints.empty()) {
          all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
        };
      };
      assert(!all_cpoints.empty());
      // Reset the query result weights to 0.
      float combination_weight = 0;
      std::fill(qr_weights.begin(), qr_weights.end(), 0);
      // Loop over the potentially involved ConnectionPoints.
      for (ConnectionPoint* cpoint : all_cpoints) {
        const std::vector<unsigned>& cp_ids = cpoint->GetIDs().GetAcyclicIDs();
        const std::vector<float>& cp_weights = cpoint->GetWeights().GetAcyclicWeights();
        std::vector<unsigned>::const_iterator it, cp_begin_it = cp_ids.begin(), cp_end_it = cp_ids.end();
        // For each Pseudofragment ID in the original query result, find the
        // ConnectionPoint's weight of the same Pseudofragment.
        for (unsigned qr_idx = 0; qr_idx < qr_ids.size(); ++qr_idx) {
          unsigned qr_id = qr_ids[qr_idx];
          it = std::lower_bound(cp_begin_it, cp_end_it, qr_id);
          assert(it != cp_end_it);
          unsigned cp_idx = it - cp_begin_it;
          float cp_weight = cp_weights[cp_idx];
          qr_weights[qr_idx] += cp_weight;
          combination_weight += cp_weight;
        };
      };
      all_cpoints.clear();
      // Store the weights.
      combination_weights.push_back(combination_weight);
      // Update the global InsertionPoint weight.
      weight += combination_weight;
    };

  // If guided evolution isn't enabled it suffices to sum the weights of the
  // query results.
  } else {
    for (const auto& s : acyclic_substitutes) {
      float combination_weight = 0;
      for (float w : s.second.second) {
        combination_weight += w;
      };
      combination_weights.push_back(combination_weight);
      weight += combination_weight;
    };
  };
  // Flag the SubstitutionPoint.
  calculated_weights = true;
};

void SubstitutionPoint::CalcRingBasedWeight() {
  assert(!calculated_weights);
  assert(retrieved_ring_substitutes);
  combination_weights.reserve(combinations.size());

  // If guided evolution is enabled the weights of the retrieved query
  // results have to be modified using the weights of the ConnectionPoints
  // involved in each Connections combination.
  if (guided) {
    // Create references to the owner ReconstructedMol's ConnectionTables.
    ConnectionsTable& peripheral_connections = substitutable->owner->connections;
    ConnectionsTable& internal_connections = substitutable->owner->internal_connections;
    std::vector<ConnectionPoint*> cpoints, all_cpoints;
    // Loop over the Connection combinations and the resulting query results.
    for (auto& s : ring_substitutes) {
      const CONNECTIONS_COMBINATION& combination = s.first;
      QUERY_RESULT& qr = s.second;
      const std::vector<unsigned>& qr_ids = qr.first;
      std::vector<float>& qr_weights = qr.second;
      // Retrieve the ConnectionPoints potentially involved in the combination.
      for (unsigned cidx = 0; cidx < n_neighbors; ++cidx) {
        cpoints = peripheral_connections.GetConnectionPoints(combination[cidx], neighbors_combination[cidx]->connections[combination[cidx]]);
        if (!cpoints.empty()) {
          all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
        };
        cpoints = internal_connections.GetConnectionPoints(combination[cidx], neighbors_combination[cidx]->connections[combination[cidx]]);
        if (!cpoints.empty()) {
          all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
        };
      };
      // Reset the query result weights to 0.
      float combination_weight = 0;
      std::fill(qr_weights.begin(), qr_weights.end(), 0);
      // Loop over the potentially involved ConnectionPoints.
      for (ConnectionPoint* cpoint : all_cpoints) {
        const std::vector<unsigned>& cp_ids = cpoint->GetIDs().GetRingIDs();
        const std::vector<float>& cp_weights = cpoint->GetWeights().GetRingWeights();
        std::vector<unsigned>::const_iterator it, cp_begin_it = cp_ids.begin(), cp_end_it = cp_ids.end();
        // For each Pseudofragment ID in the original query result, find the
        // ConnectionPoint's weight of the same Pseudofragment.
        for (unsigned qr_idx = 0; qr_idx < qr_ids.size(); ++qr_idx) {
          unsigned qr_id = qr_ids[qr_idx];
          it = std::lower_bound(cp_begin_it, cp_end_it, qr_id);
          assert(it != cp_end_it);
          unsigned cp_idx = it - cp_begin_it;
          float cp_weight = cp_weights[cp_idx];
          qr_weights[qr_idx] += cp_weight;
          combination_weight += cp_weight;
        };
      };
      all_cpoints.clear();
      // Store the weights.
      combination_weights.push_back(combination_weight);
      // Update the global InsertionPoint weight.
      weight += combination_weight;
    };

  // If guided evolution isn't enabled it suffices to sum the weights of the
  // query results.
  } else {
    for (const auto& s : ring_substitutes) {
      float combination_weight = 0;
      for (float w : s.second.second) {
        combination_weight += w;
      };
      combination_weights.push_back(combination_weight);
      weight += combination_weight;
    };
  };
  // Flag the SubstitutionPoint.
  calculated_weights = true;
};

void SubstitutionPoint::CalcInventoryAcyclicBasedWeight() {
  weight = 0;
  if (inventory_acyclic_substitutes.first.empty()) {
    return;
  };
  for (float w : inventory_acyclic_substitutes.second) {
    weight += w;
  };
};

void SubstitutionPoint::CalcInventoryRingBasedWeight() {
  weight = 0;
  if (inventory_ring_substitutes.first.empty()) {
    return;
  };
  for (float w : inventory_ring_substitutes.second) {
    weight += w;
  };
};

std::tuple<CONNECTIONS_COMBINATION, unsigned, float> SubstitutionPoint::SampleAcyclicSubstitute(std::mt19937& prng) {
  assert(has_acyclic_substitutes);
  assert(calculated_weights);
  // Choose a Connections combination through weighted selection.
  std::discrete_distribution<unsigned> distribution(combination_weights.begin(), combination_weights.end());
  unsigned idx = distribution(prng);
  // Move an iterator to the position of the chosen Connections combination.
  std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::const_iterator it = acyclic_substitutes.begin();
  std::advance(it, idx);
  const CONNECTIONS_COMBINATION& connections_combination = it->first;
  // Choose an insert through weighted selection.
  const QUERY_RESULT& query_result = it->second;
  assert(!query_result.first.empty());
  distribution = std::discrete_distribution<unsigned>(query_result.second.begin(), query_result.second.end());
  idx = distribution(prng);
  unsigned pseudofragment_id = query_result.first[idx];
  float pseudofragment_weight = query_result.second[idx];
  return std::make_tuple(connections_combination, pseudofragment_id, pseudofragment_weight);
};

std::tuple<CONNECTIONS_COMBINATION, unsigned, float> SubstitutionPoint::SampleRingSubstitute(std::mt19937& prng) {
  assert(has_ring_substitutes);
  assert(calculated_weights);
  // Choose a Connections combination through weighted selection.
  std::discrete_distribution<unsigned> distribution(combination_weights.begin(), combination_weights.end());
  unsigned idx = distribution(prng);
  // Move an iterator to the position of the chosen Connections combination.
  std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::const_iterator it = ring_substitutes.begin();
  std::advance(it, idx);
  const CONNECTIONS_COMBINATION& connections_combination = it->first;
  // Choose an insert through weighted selection.
  const QUERY_RESULT& query_result = it->second;
  assert(!query_result.first.empty());
  distribution = std::discrete_distribution<unsigned>(query_result.second.begin(), query_result.second.end());
  idx = distribution(prng);
  unsigned pseudofragment_id = query_result.first[idx];
  float pseudofragment_weight = query_result.second[idx];
  return std::make_tuple(connections_combination, pseudofragment_id, pseudofragment_weight);
};


// Class InsertionPoint
InsertionPoint::InsertionPoint(MolBrick* owner, std::mt19937& prng, bool guided) : owner(owner), guided(guided) {
  // Retrieve the neighboring MolBricks.
  std::vector<MolBrick*> neighbors = owner->GetNeighbors();
  n_neighbors = neighbors.size();
  // Perceive the available Connections of the neighbors.
  for (MolBrick* neighbor : neighbors) {
    neighbors_available_connections.insert({ neighbor, neighbor->GetAvailableConnections(prng, owner) });
  };
  // Generate all possible neighbor combinations of any size.
  std::vector<std::vector<MolBrick*>> neighbors_combinations;
  for (unsigned n = 1; n <= n_neighbors; ++n) {
    Combinations(neighbors, n, neighbors_combinations);
  };
  n_neighbors_combinations = neighbors_combinations.size();
  // Declare a set of iterators to traverse the available
  // Connections of the neighboring MolBricks.
  std::vector<CONNECTIONS_MAP::const_iterator> its, begin_its, end_its;
  its.reserve(n_neighbors + 1u);
  begin_its.reserve(n_neighbors + 1u);
  end_its.reserve(n_neighbors + 1u);
  // Initialize some other variables for performant loops.
  ConnectionsTable owner_available_connections;
  CONNECTIONS_COMBINATION connections_combination;
  connections_combination.reserve(n_neighbors + 1u);
  std::pair<std::map<MOLBRICKS_COMBINATION, std::set<CONNECTIONS_COMBINATION>>::iterator, bool> emplace_p;
  std::pair<std::set<CONNECTIONS_COMBINATION>::iterator, bool> insert_p;
  MOLBRICKS_COMBINATION owner_plus_neighbors;
  // For each neighbor combination:
  for (const auto& neighbors_combination : neighbors_combinations) {
    owner_plus_neighbors = neighbors_combination;
    // Add the owner MolBrick to the "neighbors" combination.
    owner_plus_neighbors.push_back(owner);
    unsigned n = owner_plus_neighbors.size();
    // Create an empty entry for the neighbors combination.
    emplace_p = combinations.emplace(owner_plus_neighbors, std::set<CONNECTIONS_COMBINATION>());
    std::set<CONNECTIONS_COMBINATION>& connections_combinations = emplace_p.first->second;
    // Retrieve the owner MolBrick's available Connections (i.e.
    // combination ghost bricks).
    owner_available_connections = owner->GetAvailableConnections(prng, neighbors_combination);
    // Define the iterators for both the owner's and neighbors
    // available Connections.
    for (MolBrick* neighbor : neighbors_combination) {
      const ConnectionsTable& neighbor_available_connections = neighbors_available_connections[neighbor];
      its.push_back(neighbor_available_connections.begin());
      begin_its.push_back(neighbor_available_connections.begin());
      end_its.push_back(neighbor_available_connections.end());
    };
    its.push_back(owner_available_connections.begin());
    begin_its.push_back(owner_available_connections.begin());
    end_its.push_back(owner_available_connections.end());
    // Use the available Connections of the owner and neighbor
    // MolBricks to generate the Connection combinations using an
    // "odometer-like" approach.
    while (its[0] != end_its[0]) {
      // Retrieve the current combination of the odometer.
      for (auto& it : its) {
        connections_combination.push_back(it->first);
      };
      // Store the Connections combination alongisde its MolBrick combination:
      //  - If guided evolution will be used it's important to not lose
      //    track of from which MolBrick each Connection was extracted,
      //    since the ConnectionPoint weights may be different for each MolBrick.
      //  - Otherwise, we are free to sort the Connections combination.
      //    This reduces the number of unique combinations and hence the
      //    number of multiple set intersections to be carried out.
      if (!guided) {
        std::sort(connections_combination.begin(), connections_combination.end());
      };
      insert_p = connections_combinations.insert(connections_combination);
      if (insert_p.second) {
        ++n_connections_combinations;
      };
      connections_combination.clear();
      // Increment the odometer by 1 (i.e. increment the iterators).
      ++its[n - 1];
      for (unsigned i = n - 1; (i > 0) && (its[i] == end_its[i]); --i) {
        its[i] = begin_its[i];
        ++its[i - 1];
      };
    };
    // Reset all the iterators.
    its.clear();
    begin_its.clear();
    end_its.clear();
  };
};

void InsertionPoint::RetrieveAcyclicInserts(ConnectionQueryResults& query_results) {
  // Retrieve the InsertionPoint's insertable acylic Pseudofragments taking
  // into account all possible Connection combinations.
  assert(acyclic_inserts.empty());
  std::pair<std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::iterator, bool> it_b_p;
  // Loop over the Connection combinations.
  for (const auto& combination : combinations) {
    for (const CONNECTIONS_COMBINATION& connection_combination : combination.second) {
      // If the query result for the Connections combination hasn't been
      // calculated yet, calculate it.
      if (acyclic_inserts.find(connection_combination) == acyclic_inserts.end()) {
        // Create an empty query result entry for the combination.
        it_b_p = acyclic_inserts.emplace(connection_combination, QUERY_RESULT());
        QUERY_RESULT& qr = it_b_p.first->second;
        // Calculate the multiple set intersection of the Pseudofragments
        // compatible with the combination.
        qr = query_results.SlidingIntersection(connection_combination, false, false);
        if (!qr.first.empty()) {
          has_acyclic_inserts = true;
        };
      };
    };
  };
  // Flag the SubstitutionPoint.
  retrieved_acyclic_inserts = true;
};

void InsertionPoint::RetrieveRingInserts(ConnectionQueryResults& query_results) {
  // Retrieve the InsertionPoint's insertable acylic Pseudofragments taking
  // into account all possible Connection combinations.
  assert(ring_inserts.empty());
  std::pair<std::map<CONNECTIONS_COMBINATION, QUERY_RESULT>::iterator, bool> it_b_p;
  // Loop over the Connection combinations.
  for (const auto& combination : combinations) {
    for (const CONNECTIONS_COMBINATION& connection_combination : combination.second) {
      // If the query result for the Connections combination hasn't been
      // calculated yet, calculate it.
      if (ring_inserts.find(connection_combination) == ring_inserts.end()) {
        // Create an empty query result entry for the combination.
        it_b_p = ring_inserts.emplace(connection_combination, QUERY_RESULT());
        QUERY_RESULT& qr = it_b_p.first->second;
        // Calculate the multiple set intersection of the Pseudofragments
        // compatible with the combination.
        qr = query_results.SlidingIntersection(connection_combination, true, false);
        if (!qr.first.empty()) {
          has_ring_inserts = true;
        };
      };
    };
  };
  // Flag the SubstitutionPoint.
  retrieved_ring_inserts = true;
};

void InsertionPoint::CalcAcyclicBasedWeight() {
  assert(has_acyclic_inserts);
  assert(!calculated_weights);
  combination_weights.reserve(n_connections_combinations);

  // If the evolution is being guided, the weights of the Connection combinations
  // must be calculated based on the ConnectionPoints involved in the combination.
  if (guided) {
    // Create references to the owner ReconstructedMol's ConnectionTables.
    ConnectionsTable& peripheral_connections = owner->owner->connections;
    ConnectionsTable& internal_connections = owner->owner->internal_connections;
    std::vector<ConnectionPoint*> cpoints, all_cpoints;
    // Loop over the neighbor MolBrick combinations.
    for (const auto& combination : combinations) {
      const MOLBRICKS_COMBINATION& neighbors_combination = combination.first;
      // Loop over the Connections combinations.
      for (const CONNECTIONS_COMBINATION& connections_combination : combination.second) {
        // Retrieve the raw query result associated with the Connections combination.
        QUERY_RESULT& query_result = acyclic_inserts[connections_combination];
        const std::vector<unsigned>& qr_ids = query_result.first;
        // Retrieve all the neighbors' ConnectionPoints that could be involved in the
        // Connections combination. Since the ConnectionPoint may be considered
        // available either because it's peripheral or connected to another neighbor
        // of the neighbors combination, both the peripheral and internal ConnectionTables
        // ought to be searched.
        assert(neighbors_combination.size() == connections_combination.size());
        for (unsigned cidx = 0; cidx < connections_combination.size(); ++cidx) {
          cpoints = peripheral_connections.GetConnectionPoints(connections_combination[cidx], neighbors_combination[cidx]->connections[connections_combination[cidx]]);
          if (!cpoints.empty()) {
            all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
          };
          cpoints = internal_connections.GetConnectionPoints(connections_combination[cidx], neighbors_combination[cidx]->connections[connections_combination[cidx]]);
          if (!cpoints.empty()) {
            all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
          };
        };
        // Initialize a vector to hold the adjusted weights of the query result.
        float combination_weight = 0;
        std::vector<float> adjusted_qr_weights (qr_ids.size());
        // Loop over the potentially involved ConnectionPoints.
        for (ConnectionPoint* cpoint : all_cpoints) {
          const std::vector<unsigned>& cp_ids = cpoint->GetIDs().GetAcyclicIDs();
          const std::vector<float>& cp_weights = cpoint->GetWeights().GetAcyclicWeights();
          std::vector<unsigned>::const_iterator it, cp_begin_it = cp_ids.begin(), cp_end_it = cp_ids.end();
          // For each Pseudofragment ID in the original query result, find the
          // ConnectionPoint's weight of the same Pseudofragment.
          for (unsigned qr_idx = 0; qr_idx < qr_ids.size(); ++qr_idx) {
            unsigned qr_id = qr_ids[qr_idx];
            it = std::lower_bound(cp_begin_it, cp_end_it, qr_id);
            assert(it != cp_end_it);
            unsigned cp_idx = it - cp_begin_it;
            float cp_weight = cp_weights[cp_idx];
            adjusted_qr_weights[qr_idx] += cp_weight;
            combination_weight += cp_weight;
          };
        };
        all_cpoints.clear();
        // Store the weights.
        insert_weights.insert({ {neighbors_combination, connections_combination}, std::move(adjusted_qr_weights) });
        combination_weights.push_back(combination_weight);
        // Update the global InsertionPoint weight.
        weight += combination_weight;
      };
    };

  // If the evolution isn't guided, the weights of each Connections combination
  // can be calculated directly from the query results.
  } else {
    // Loop over the neighbor MolBrick combinations.
    for (const auto& combination : combinations) {
      const MOLBRICKS_COMBINATION& neighbors_combination = combination.first;
      // Loop over the Connections combinations.
      for (const CONNECTIONS_COMBINATION& connections_combination : combination.second) {
        // Retrieve the raw query result associated with the Connections combination.
        QUERY_RESULT& query_result = acyclic_inserts[connections_combination];
        const std::vector<float>& qr_weights = query_result.second;
        // Calculate the Connections combination's weight as the sum of its
        // Pseudofragment weights.
        float combination_weight = 0;
        for (float w : qr_weights) {
          combination_weight += w;
        };
        // Store the weights.
        insert_weights.insert({ {neighbors_combination, connections_combination}, std::vector<float>() });
        combination_weights.push_back(combination_weight);
        // Update the global InsertionPoint weight.
        weight += combination_weight;
      };
    };
  };
  calculated_weights = true;
};

void InsertionPoint::CalcRingBasedWeight() {
  assert(has_ring_inserts);
  assert(!calculated_weights);
  combination_weights.reserve(n_connections_combinations);

  // If the evolution is being guided, the weights of the Connection combinations
  // must be calculated based on the ConnectionPoints involved in the combination.
  if (guided) {
    // Create references to the owner ReconstructedMol's ConnectionTables.
    ConnectionsTable& peripheral_connections = owner->owner->connections;
    ConnectionsTable& internal_connections = owner->owner->internal_connections;
    std::vector<ConnectionPoint*> cpoints, all_cpoints;
    // Loop over the neighbor MolBrick combinations.
    for (const auto& combination : combinations) {
      const MOLBRICKS_COMBINATION& neighbors_combination = combination.first;
      // Loop over the Connections combinations.
      for (const CONNECTIONS_COMBINATION& connections_combination : combination.second) {
        // Retrieve the raw query result associated with the Connections combination.
        QUERY_RESULT& query_result = ring_inserts[connections_combination];
        const std::vector<unsigned>& qr_ids = query_result.first;
        // Retrieve all the neighbors' ConnectionPoints that could be involved in the
        // Connections combination. Since the ConnectionPoint may be considered
        // available either because it's peripheral or connected to another neighbor
        // of the neighbors combination, both the peripheral and internal ConnectionTables
        // ought to be searched.
        assert(neighbors_combination.size() == connections_combination.size());
        for (unsigned cidx = 0; cidx < connections_combination.size(); ++cidx) {
          cpoints = peripheral_connections.GetConnectionPoints(connections_combination[cidx], neighbors_combination[cidx]->connections[connections_combination[cidx]]);
          if (!cpoints.empty()) {
            all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
          };
          cpoints = internal_connections.GetConnectionPoints(connections_combination[cidx], neighbors_combination[cidx]->connections[connections_combination[cidx]]);
          if (!cpoints.empty()) {
            all_cpoints.insert(all_cpoints.end(), cpoints.begin(), cpoints.end());
          };
        };
        // Initialize a vector to hold the adjusted weights of the query result.
        float combination_weight = 0;
        std::vector<float> adjusted_qr_weights(qr_ids.size());
        // Loop over the potentially involved ConnectionPoints.
        for (ConnectionPoint* cpoint : all_cpoints) {
          const std::vector<unsigned>& cp_ids = cpoint->GetIDs().GetRingIDs();
          const std::vector<float>& cp_weights = cpoint->GetWeights().GetRingWeights();
          std::vector<unsigned>::const_iterator it, cp_begin_it = cp_ids.begin(), cp_end_it = cp_ids.end();
          // For each Pseudofragment ID in the original query result, find the
          // ConnectionPoint's weight of the same Pseudofragment.
          for (unsigned qr_idx = 0; qr_idx < qr_ids.size(); ++qr_idx) {
            unsigned qr_id = qr_ids[qr_idx];
            it = std::lower_bound(cp_begin_it, cp_end_it, qr_id);
            assert(it != cp_end_it);
            unsigned cp_idx = it - cp_begin_it;
            float cp_weight = cp_weights[cp_idx];
            adjusted_qr_weights[qr_idx] += cp_weight;
            combination_weight += cp_weight;
          };
        };
        all_cpoints.clear();
        // Store the weights.
        insert_weights.insert({ {neighbors_combination, connections_combination}, std::move(adjusted_qr_weights) });
        combination_weights.push_back(combination_weight);
        // Update the global InsertionPoint weight.
        weight += combination_weight;
      };
    };

  // If the evolution isn't guided, the weights of each Connections combination
  // can be calculated directly from the query results.
  } else {
    // Loop over the neighbor MolBrick combinations.
    for (const auto& combination : combinations) {
      const MOLBRICKS_COMBINATION& neighbors_combination = combination.first;
      // Loop over the Connections combinations.
      for (const CONNECTIONS_COMBINATION& connections_combination : combination.second) {
        // Retrieve the raw query result associated with the Connections combination.
        QUERY_RESULT& query_result = ring_inserts[connections_combination];
        const std::vector<float>& qr_weights = query_result.second;
        // Calculate the Connections combination's weight as the sum of its
        // Pseudofragment weights.
        float combination_weight = 0;
        for (float w : qr_weights) {
          combination_weight += w;
        };
        // Store the weights.
        insert_weights.insert({ {neighbors_combination, connections_combination}, std::vector<float>() });
        combination_weights.push_back(combination_weight);
        // Update the global InsertionPoint weight.
        weight += combination_weight;
      };
    };
  };
  calculated_weights = true;
};

void InsertionPoint::GenerateMatchingsForBrick(const MolBrick* brick, ConnectionCompatibilities& compatibilities, std::mt19937& prng) {
  // Generates all possible matchings between the InsertionPoint's neighbors
  // combinations and the putative insert MolBrick's Connections.
  // NOTE: an InsertionPoint can only store matchings for a single insert MolBrick at a time.
  assert(matchings.empty());

  // Loop over the neighboring MolBrick combinations.
  for (const auto& combination : combinations) {
    MOLBRICKS_COMBINATION neighbors_combination = combination.first;
    // Remove the owner MolBrick from the combination so that only the real
    // neighbors remain.
    neighbors_combination.pop_back();

    // Get the InsertionPoint's owner MolBrick's available Connections.
    ConnectionsTable owner_available_connections = owner->GetAvailableConnections(prng, neighbors_combination);

    // Determine how the MolBrick insert should connect to the neighboring
    // MolBricks of the InsertionPoint by finding a solution to the Maximum
    // Bipartite Matching problem.
    std::vector<std::pair<int, int>> edges;
    std::map<int, const MolBrick*> vtx_ids_neighbors;
    std::map<const MolBrick*, int> neighbors_vtx_ids;
    std::map<int, const Connection*> vtx_ids_connections;

    int id = 1;
    vtx_ids_neighbors.insert({ id, owner });
    neighbors_vtx_ids.insert({ owner, id });
    ++id;
    for (MolBrick* neighbor : neighbors_combination) {
      vtx_ids_neighbors.insert({ id, neighbor });
      neighbors_vtx_ids.insert({ neighbor, id });
      ++id;
    };
    for (const auto& c : brick->GetConnections()) {
      for (const auto& cp : c.second) {
        vtx_ids_connections.insert({ id, &c.first });
        if (owner_available_connections.IsCompatibleWith(c.first, compatibilities)) {
          edges.push_back(std::make_pair(1, id));
        };
        for (MolBrick* neighbor : neighbors_combination) {
          if (neighbors_available_connections[neighbor].IsCompatibleWith(c.first, compatibilities)) {
            edges.push_back(std::make_pair(neighbors_vtx_ids[neighbor], id));
          };
        };
        ++id;
      };
    };

    BipartiteGraph graph(vtx_ids_neighbors.size(), vtx_ids_connections.size(), edges);
    PATH matching = graph.HopcroftKarp(prng);

    // If a perfect matching was found, store it.
    if (matching.size() == neighbors_combination.size() + 1) {
      std::pair<std::map<MOLBRICKS_COMBINATION, std::map<const MolBrick*, const Connection*>>::iterator, bool> it_b_p;
      it_b_p = matchings.emplace(neighbors_combination, std::map<const MolBrick*, const Connection*>());
      std::map<const MolBrick*, const Connection*>& neighbors_matching = it_b_p.first->second;
      for (Edge* edge : matching) {
        const MolBrick* match_brick = vtx_ids_neighbors[edge->GetStart()->GetID()];
        const Connection* match_connection = vtx_ids_connections[edge->GetEnd()->GetID()];
        neighbors_matching.insert({ match_brick, match_connection });
      };
      ++n_matchings;
    };
  };
};

std::tuple<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION, unsigned, float> InsertionPoint::SampleAcyclicInsert(std::mt19937& prng) {
  assert(has_acyclic_inserts);
  assert(calculated_weights);
  // Choose a neighbors combination through weighted selection.
  std::discrete_distribution<unsigned> distribution (combination_weights.begin(), combination_weights.end());
  unsigned idx = distribution(prng);
  // Move an iterator to the position of the chosen neighbors combination.
  std::map<std::pair<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION>, std::vector<float>>::const_iterator it = insert_weights.begin();
  std::advance(it, idx);
  const MOLBRICKS_COMBINATION& neighbors_combination = it->first.first;
  const CONNECTIONS_COMBINATION& connections_combination = it->first.second;
  // Choose an insert through weighted selection.
  QUERY_RESULT& query_result = acyclic_inserts[it->first.second];
  if (guided) {
    const std::vector<float>& adjusted_weights = it->second;
    assert(!it->second.empty());
    distribution = std::discrete_distribution<unsigned>(adjusted_weights.begin(), adjusted_weights.end());
  } else {
    distribution = std::discrete_distribution<unsigned>(query_result.second.begin(), query_result.second.end());
  };
  idx = distribution(prng);
  unsigned pseudofragment_id = query_result.first[idx];
  float pseudofragment_weight = query_result.second[idx];
  return std::make_tuple(neighbors_combination, connections_combination, pseudofragment_id, pseudofragment_weight);
};

std::tuple<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION, unsigned, float> InsertionPoint::SampleRingInsert(std::mt19937& prng) {
  assert(has_ring_inserts);
  assert(calculated_weights);
  // Choose a neighbors combination through weighted selection.
  std::discrete_distribution<unsigned> distribution(combination_weights.begin(), combination_weights.end());
  unsigned idx = distribution(prng);
  // Move an iterator to the position of the chosen neighbors combination.
  std::map<std::pair<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION>, std::vector<float>>::const_iterator it = insert_weights.begin();
  std::advance(it, idx);
  const MOLBRICKS_COMBINATION& neighbors_combination = it->first.first;
  const CONNECTIONS_COMBINATION& connections_combination = it->first.second;
  // Choose an insert through weighted selection.
  QUERY_RESULT& query_result = ring_inserts[it->first.second];
  if (guided) {
    const std::vector<float>& adjusted_weights = it->second;
    assert(!it->second.empty());
    distribution = std::discrete_distribution<unsigned>(adjusted_weights.begin(), adjusted_weights.end());
  } else {
    distribution = std::discrete_distribution<unsigned>(query_result.second.begin(), query_result.second.end());
  };
  idx = distribution(prng);
  unsigned pseudofragment_id = query_result.first[idx];
  float pseudofragment_weight = query_result.second[idx];
  return std::make_tuple(neighbors_combination, connections_combination, pseudofragment_id, pseudofragment_weight);
};


// Class GeneticLogEntry
GeneticLogEntry::GeneticLogEntry(GeneticLogEntry::Type type, const Connection& connection, unsigned atom_idx, unsigned pseudofragment_id) :
  type(type), connection(connection), atom_idx(atom_idx), pseudofragment_id(pseudofragment_id) {};

GeneticLogEntry::Type GeneticLogEntry::GetType() const {
  return type;
};

const Connection& GeneticLogEntry::GetConnection() const {
  return connection;
};

unsigned GeneticLogEntry::GetAtomIdx() const {
  return atom_idx;
};

unsigned GeneticLogEntry::GetPseudofragmentID() const {
  return pseudofragment_id;
};


// Class GeneticLog
GeneticLog::GeneticLog() = default;

bool GeneticLog::SetImportantPseudofragmentID(unsigned pseudofragment_id) {
  if (pseudofragment_id1 == pseudofragment_id || pseudofragment_id2 == pseudofragment_id) {
    return true;
  };
  if (pseudofragment_id1 == 0u) {
    pseudofragment_id1 = pseudofragment_id;
    return true;
  };
  if (pseudofragment_id2 == 0u) {
    pseudofragment_id2 = pseudofragment_id;
    return true;
  };
  return false;
};

void GeneticLog::AddEntry(GeneticLogEntry::Type type, const Connection& connection, unsigned atom_idx, unsigned pseudofragment_id) {
  std::map<unsigned, std::vector<GeneticLogEntry>>::iterator it = entries.find(pseudofragment_id);
  if (it != entries.end()) {
    it->second.emplace_back(type, connection, atom_idx, pseudofragment_id);
  } else {
    entries.insert({pseudofragment_id, std::vector<GeneticLogEntry>{GeneticLogEntry(type, connection, atom_idx, pseudofragment_id)}});
  };
};

void GeneticLog::Clear() {
  pseudofragment_id1 = 0u;
  pseudofragment_id2 = 0u;
  entries.clear();
};

unsigned GeneticLog::GetImportantPseudofragmentID1() const {
  return pseudofragment_id1;
};

unsigned GeneticLog::GetImportantPseudofragmentID2() const {
  return pseudofragment_id2;
};

const std::map<unsigned, std::vector<GeneticLogEntry>>& GeneticLog::GetEntries() const {
  return entries;
};


// Class ReconstructedMol
ReconstructedMol::ReconstructedMol() = default;
ReconstructedMol::ReconstructedMol(const MolBrick & brick) :
  pseudomol(brick.pseudomol),
  connections(brick.connections),
  level(brick.level),
  max_atom_idx(brick.max_atom_idx),
  max_schematic_idx(brick.schematic_idx) {
  assert(brick.schematic_idx == 0);
  if (brick.HasRing()) {
    n_ring_atoms = brick.level;
  } else {
    n_ring_atoms = 0;
  };
  bricks.insert({ brick.schematic_idx, brick });
  schematic.UnitSchematic();
};

int ReconstructedMol::EvaluatePeripheralCompatibility(ConnectionQueryResults & query_results) const {
  // Evaluates what kind of Pseudofragments can be annealed in a Peripheral
  // Expansion by consulting the ConnectionQueryResults.
  // Possible return values:
  //  -1 = Can only be expanded with acyclic Pseudofragments.
  //   0 = Can be expanded with both acyclic and ring Pseudofragments (including ring_parts).
  //   1 = Can only be expanded with cyclic Pseudofragments.
  bool ring_compatible = false, acyclic_compatible = false;

  // Loop over the peripheral Connections. If it is still unknown if the
  // ReconstructedMol can bind peripherally Pseudofragments of a specific type,
  // check if the Connection supports said type of fragment.
  for (const auto& c : connections) {
    const Connection& connection = c.first;
    if (!ring_compatible) {
      if (query_results.ring[connection][1].first.size() > 0) {
        ring_compatible = true;
      };
    };
    if (!acyclic_compatible) {
      if (query_results.acyclic[connection][1].first.size() > 0) {
        acyclic_compatible = true;
      };
    };
    // If it has been determined that the ReconstructedMol can bind both
    // cylic and acyclic Pseudofragments, return this evaluation and shortcircuit.
    if (ring_compatible && acyclic_compatible) {
      return 0;
    };
  };

  // If the ReconstructedMol can't bind both types of Pseudofragments it can
  // only bind one. Return with the corresponding value.
  if (acyclic_compatible) {
    return -1;
  } else if (ring_compatible) {
    return 1;
  } else {
    // In theory, during fragmentation, when a bond is broken and a Connection
    // is defined there ought to be atleast one compatible Pseudofragment.
    // However, when defining fragments through a subgraph search, bonds
    // connecting to peripheral atoms are only broken into Connections on the
    // internal side of the molecule, whereas on the peripheral side no path of
    // 2 atoms length can be defined. While peripheral atoms will still be
    // represented in the database, if the compatibility definition used when
    // precalculating compatible Pseudofragments is too strict, some connections
    // won't be deemed compatible with any Pseudofragment. This might occur to
    // exclusively peripheral functional groups (e.g. nitro) when using a high
    // stringency level (e.g. 4).
    return -2;
  };
};

int ReconstructedMol::EvaluateBrickSet(const std::vector<MolBrick*>& b) const {
  // Evaluates the types of MolBricks that form a set.
  // Possible return values:
  //  -1 = Consists of only acyclic MolBricks.
  //   0 = Consists of both acyclic and cyclic MolBricks (including ring_parts).
  //   1 = Consists of only cyclic MolBricks.
  bool contains_rings = false;
  bool contains_acyclics = false;
  for (const MolBrick* brick : b) {
    if (brick->has_ring || brick->ring_part) {
      contains_rings = true;
    } else {
      contains_acyclics = true;
    };
    if (contains_rings && contains_acyclics) {
      return 0;
    };
  };
  if (contains_acyclics) {
    return -1;
  } else if (contains_rings) {
    return 1;
  };
  throw std::runtime_error(std::string("Brick set evaluation failed\n"));
};

std::pair<const Connection*, const ConnectionPoint*> ReconstructedMol::ChoosePeripheralConnectionPoint(ConnectionQueryResults& query_results, bool has_ring, bool ring_part, std::mt19937& prng) const {
  // Initialize vectors to contain the candidate Connections and their weights.
  std::vector<const Connection*> candidates;
  std::vector<float> weights;
  candidates.reserve(connections.n_connection_points);
  weights.reserve(connections.n_connection_points);

  // If the weight should be calculated based on ring Pseudofragments:
  if (has_ring) {
    // For each peripheral Connection
    for (const auto& c : connections) {
      // Retrieve the compatible Pseudofragments
      QUERY_RESULT& qr = query_results.ring[c.first][1];
      // Sum their weights and use it as the Connection weight.
      float weights_sum = 0;
      for (float weight : qr.second) {
        weights_sum += weight;
      };
      // Store the Connection and its weight.
      candidates.push_back(&c.first);
      weights.push_back(weights_sum);
    };
  // If the weight should be calculated based on ring_part Pseudofragments:
  } else if (ring_part) {
    for (const auto& c : connections) {
      QUERY_RESULT& qr = query_results.ring_part[c.first][1];
      float weights_sum = 0;
      for (float weight : qr.second) {
        weights_sum += weight;
      };
      candidates.push_back(&c.first);
      weights.push_back(weights_sum);
    };
  // If the weight should be calculated based on acyclic Pseudofragments:
  } else {
    for (const auto& c : connections) {
      QUERY_RESULT& qr = query_results.acyclic[c.first][1];
      float weights_sum = 0;
      for (float weight : qr.second) {
        weights_sum += weight;
      };
      candidates.push_back(&c.first);
      weights.push_back(weights_sum);
    };
  };

  // Perform a weighted selection of the Connection.
  std::discrete_distribution<unsigned> distribution(weights.begin(), weights.end());
  const Connection* chosen_connection = candidates[distribution(prng)];
  // Choose a random ConnectionPoint of said Connection.
  unsigned idx = 0;
  const std::vector<ConnectionPoint>& cpoints = connections.at(*chosen_connection);
  unsigned n_cpoints = cpoints.size();
  if (n_cpoints > 1) {
    if (n_cpoints <= max_size_uniform_int_distributions) {
      idx = uniform_int_distributions[n_cpoints](prng);
    } else {
      std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_cpoints - 1);
      idx = cpoint_distribution(prng);
    };
  };
  const ConnectionPoint* chosen_cpoint = &cpoints[idx];
  return std::make_pair(chosen_connection, chosen_cpoint);
};

std::pair<const Connection*, const ConnectionPoint*> ReconstructedMol::ChoosePeripheralConnectionPoint(bool has_ring, bool ring_part, std::mt19937& prng) const {
  // Initialize vectors to contain the candidate ConnectionPoints and their weights.
  std::vector<std::pair<const Connection*, const ConnectionPoint*>> candidates;
  std::vector<float> weights;
  candidates.reserve(connections.n_connection_points);
  weights.reserve(connections.n_connection_points);

  // If the weight should be calculated based on ring Pseudofragments:
  if (has_ring) {
    // For each ConnectionPoint
    for (const auto& c : connections) {
      for (const auto& cp : c.second) {
        // Sum up the current weights of the compatible ring Pseudofragments.
        float weights_sum = 0;
        for (float weight : cp.weights.GetRingWeights()) {
          weights_sum += weight;
        };
        // Store the ConnectionPoint and its weight.
        candidates.push_back({ &c.first, &cp });
        weights.push_back(weights_sum);
      };
    };
  // If the weight should be calculated based on ring_part Pseudofragments:
  } else if (ring_part) {
    for (const auto& c : connections) {
      for (const auto& cp : c.second) {
        float weights_sum = 0;
        for (float weight : cp.weights.GetRingPartWeights()) {
          weights_sum += weight;
        };
        candidates.push_back({ &c.first, &cp });
        weights.push_back(weights_sum);
      };
    };
  // If the weight should be calculated based on acyclic Pseudofragments:
  } else {
    for (const auto& c : connections) {
      for (const auto& cp : c.second) {
        float weights_sum = 0;
        for (float weight : cp.weights.GetAcyclicWeights()) {
          weights_sum += weight;
        };
        candidates.push_back({ &c.first, &cp });
        weights.push_back(weights_sum);
      };
    };
  };

  // Perform the weighted selection.
  std::discrete_distribution<unsigned> distribution(weights.begin(), weights.end());
  return candidates[distribution(prng)];
};

MolBrick* ReconstructedMol::ChooseBrickFromSet(const std::vector<MolBrick*>& b, bool has_ring, bool ring_part, std::mt19937 & prng) const {
  unsigned n_bricks = b.size();
  std::vector<MolBrick*> candidates;
  std::vector<float> weights;
  candidates.reserve(n_bricks);
  weights.reserve(n_bricks);
  if (has_ring) {
    for (MolBrick* brick : b) {
      if (brick->has_ring) {
        candidates.push_back(brick);
        weights.push_back(brick->weight);
      };
    };
  } else if (ring_part) {
    for (MolBrick* brick : b) {
      if (brick->ring_part) {
        candidates.push_back(brick);
        weights.push_back(brick->weight);
      };
    };
  } else {
    for (MolBrick* brick : b) {
      if (!brick->has_ring && !brick->ring_part) {
        candidates.push_back(brick);
        weights.push_back(brick->weight);
      };
    };
  };
  std::discrete_distribution<unsigned> distribution(weights.begin(), weights.end());
  return candidates[distribution(prng)];
};

const MolBrick* ReconstructedMol::GetAtomSourceBrick(unsigned immutable_atom_idx) const {
  std::vector<unsigned>::const_iterator it, end_it;
  for (auto& brick : bricks) {
    const std::vector<unsigned>& atom_indices = brick.second.atom_indices;
    end_it = atom_indices.end();
    it = std::find(atom_indices.begin(), end_it, immutable_atom_idx);
    if (it != end_it) {
      return &brick.second;
    };
  };
  throw std::runtime_error(std::string("The provided immutable atom index isn't part of any of the ReconstructedMol's MolBricks. Sanity has been lost."));
};

std::vector<MolBrick*> ReconstructedMol::GetPeripheralBricks() {
  std::vector<MolBrick*> peripheral_bricks;
  unsigned schematic_idx = 0;
  unsigned n_adjacent_bricks = 0;
  for (const auto& row : schematic.adjacency) {
    n_adjacent_bricks = 0;
    for (const unsigned adjacency_bool : row) {
      if (adjacency_bool == 1) {
        ++n_adjacent_bricks;
        // For whatever reason I can't substitute the dual if
        // statement with a goto on MSVC17.
        if (n_adjacent_bricks > 1) {
          break;
        };
      };
    };
    if (n_adjacent_bricks == 1) {
      peripheral_bricks.push_back(&(bricks[schematic_idx]));
    };
    ++schematic_idx;
  };
  return peripheral_bricks;
};

std::vector<MolBrick*> ReconstructedMol::GetInternalBricks() {
  std::vector<MolBrick*> internal_bricks;
  unsigned schematic_idx = 0;
  unsigned n_adjacent_bricks;
  for (const auto& row : schematic.adjacency) {
    n_adjacent_bricks = 0;
    for (const unsigned adjacency_bool : row) {
      if (adjacency_bool == 1) {
        ++n_adjacent_bricks;
        if (n_adjacent_bricks > 1) {
          internal_bricks.push_back(&bricks[schematic_idx]);
          break;
        };
      };
    };
    ++schematic_idx;
  };
  return internal_bricks;
};


// ReconstructedMol's genetic algorithm operator functions.
bool ReconstructedMol::PeripheralExpansion(sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, bool guided, bool update_stereo) {
  // In a peripheral expansion a new Pseudofragment is annealed to one of the
  // free connection points of the ReconstructedMol.

  // If no peripheral Connections are available, signal function failure.
  if (connections.Empty()) {
    return false;
  };

  // Evaluate with which kind of Pseudofragments the ReconstructedMol can be expanded:
  //  -1 = Can only be expanded with acyclic Pseudofragments.
  //   0 = Can be expanded with both acyclic and ring Pseudofragments (including ring_parts).
  //   1 = Can only be expanded with cyclic Pseudofragments.
  int compatibility = EvaluatePeripheralCompatibility(query_results);
  if (compatibility == -2) {
    std::cout << "WARNING: No compatible Pseudofragments were found for a Connection." << std::endl;
    return false;
  };

  // If both kinds of Pseudofragments can be annealed, decide which one to use.
  // Otherwise, use the only type that can be used.
  bool add_ring = false;
  if (compatibility == 0) {
    if (n_ring_atoms < controller.max_size) {
      add_ring = controller.DecideChange(n_ring_atoms, prng);
    };
  } else if (compatibility == 1) {
    // If the peripheral Connections are only compatible with ring Pseudofragments
    // but the maximum number of rings was already reached, signal failure.
    if (n_ring_atoms >= controller.max_size) {
      return false;
    };
    add_ring = true;
  };

  bool has_ring = false, ring_part = false;
  if (add_ring) {
    has_ring = true;
  };

  // Perform a weighted selection of a peripheral ConnectionPoint that is compatible
  // with the specified type of Pseudofragment, as well as a compatible Pseudofragment.
  std::pair<const Connection*, const ConnectionPoint*> choice;
  unsigned pseudofragment_id;
  float pseudofragment_weight;
  if (guided) {
    // If guided evolution is enabled the weights used to select the ConnectionPoint
    // are those of the ConnectionPoints themselves.
    choice = ChoosePeripheralConnectionPoint(has_ring, ring_part, prng);
    // Get the IDs and weights of the Pseudofragments that are compatible with the
    // chosen Connection.
    const std::vector<unsigned>* ids;
    const std::vector<float>* weights;
    if (has_ring) {
      ids = &choice.second->ids.GetRingIDs();
      weights = &choice.second->weights.GetRingWeights();
    } else if (ring_part) {
      ids = &choice.second->ids.GetRingPartIDs();
      weights = &choice.second->weights.GetRingPartWeights();
    } else {
      ids = &choice.second->ids.GetAcyclicIDs();
      weights = &choice.second->weights.GetAcyclicWeights();
    };
    // Check if any compatible fragments exist. If they do, sample the ID of
    // one of them randomly. If not, signal function failure. The latter should
    // almost never be the case (see EvaluatePeripheralCompatibility).
    if (ids->empty()) {
      return false;
    };
    // Perform a weighted selection of a random Pseudofragment.
    std::discrete_distribution<unsigned> distribution(weights->begin(), weights->end());
    unsigned idx = distribution(prng);
    pseudofragment_id = (*ids)[idx];
    pseudofragment_weight = (*weights)[idx];
  } else {
    // If the evolution isn't guided the ConnectionQueryResults weights are used to
    // select a Connection and thereafter a random ConnectionPoint is chosen.
    choice = ChoosePeripheralConnectionPoint(query_results, has_ring, ring_part, prng);
    const QUERY_RESULT* qr;
    if (has_ring) {
      qr = &query_results.ring[*choice.first][1];
    } else if (ring_part) {
      qr = &query_results.ring_part[*choice.first][1];
    } else {
      qr = &query_results.acyclic[*choice.first][1];
    };
    if (qr->first.empty()) {
      return false;
    };
    std::discrete_distribution<unsigned> distribution(qr->second.begin(), qr->second.end());
    unsigned idx = distribution(prng);
    pseudofragment_id = qr->first[idx];
    pseudofragment_weight = qr->second[idx];
  };
  // Unpack the chosen Connection and ConnectionPoint's atom index.
  Connection chosen_connection = *choice.first;
  unsigned immutable_reconstruction_interactor_idx = choice.second->atom_idx;
  RDKit::Bond::BondType bond_type = int_bond_type_table[chosen_connection.bond_type];

  // Retrieve the Pseudofragment with the specified ID.
  Pseudofragment pseudofragment = GetPseudofragmentByIDFromDB(pseudofragment_id, select_statement);

  // Determine which ReconstructionSchematic index should be assigned to the
  // new MolBrick. If any indices became available due to a deletion, take one
  // of them. If not, assign a new one.
  unsigned schematic_idx;
  if (available_schematic_idxs.empty()) {
    schematic.Expand();
    schematic_idx = max_schematic_idx + 1;
    ++max_schematic_idx;
  } else {
    schematic_idx = available_schematic_idxs.back();
    available_schematic_idxs.pop_back();
  };

  // Convert the Pseudofragment to a MolBrick.
  MolBrick brick(pseudofragment, pseudofragment_id, max_atom_idx + 1, schematic_idx, pseudofragment_weight, this);

  // For the chosen ReconstructedMol's ConnectionPoint, randomly choose a
  // compatible interaction atom on the MolBrick.
  Connection compatible_connection = brick.connections.GetRandomCompatibleConnection(chosen_connection, query_results.compatibilities, prng);
  const CONNECTION_POINT_VECTOR& brick_cpoint_vector = brick.connections[compatible_connection];
  unsigned n_connections = brick_cpoint_vector.size();
  unsigned connection_idx = 0;
  if (n_connections > 1) {
    if (n_connections <= max_size_uniform_int_distributions) {
      connection_idx = uniform_int_distributions[n_connections](prng);
    } else {
      std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
      connection_idx = cpoint_distribution(prng);
    };
  };
  unsigned immutable_brick_interactor_idx = brick_cpoint_vector[connection_idx].atom_idx;

  // Update the pseudomol by combining the atoms of the ReconstructedMol's
  // pseudomol and the MolBrick's pseudomol and by creating a bond between
  // the specified interaction atoms.
  pseudomol.insertMol(brick.pseudomol);
  unsigned mutable_reconstruction_interactor_idx = GetMutableAtomIdx(pseudomol, immutable_reconstruction_interactor_idx);
  unsigned mutable_brick_interactor_idx = GetMutableAtomIdx(pseudomol, immutable_brick_interactor_idx);
  pseudomol.addBond(mutable_reconstruction_interactor_idx, mutable_brick_interactor_idx, bond_type);

  // Get the schematic index of the interacting brick.
  const MolBrick* interactor = GetAtomSourceBrick(immutable_reconstruction_interactor_idx);
  unsigned interactor_schematic_idx = interactor->schematic_idx;

  // Update the peripheral and internal connection tables. This entails:
  //  (1) Adding all the Connections from the new brick to the peripheral
  //      ConnectionTable.
  //  (2) Moving the two Connections involved in connecting the original
  //      ReconstructedMol to the new brick from the peripheral to the internal
  //      ConnectionTable.
  if (guided) {
    std::vector<ConnectionPoint*> cpoints;
    for (auto& c : brick.connections) {
      cpoints = connections.AddConnections(c.first, c.second);
      // Assign fresh weights to the insert's ConnectionPoints.
      query_results.AssignFreshWeights(c.first, cpoints);
    };
    assert(AllConnectionsHaveWeights());
  } else {
    for (auto& c : brick.connections) {
      connections.AddConnections(c.first, c.second);
    };
  };
  connections.MoveConnectionTo(chosen_connection, immutable_reconstruction_interactor_idx, internal_connections);
  connections.MoveConnectionTo(compatible_connection, immutable_brick_interactor_idx, internal_connections);

  // Reset the GeneticLog and record the new neighbor-to-insert link.
  if (guided) {
    genetic_log.Clear();
    genetic_log.SetImportantPseudofragmentID(pseudofragment_id);
    genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, chosen_connection, immutable_reconstruction_interactor_idx, pseudofragment_id);
  };

  // Update the ReconstructionSchematic.
  schematic.adjacency[interactor_schematic_idx][brick.schematic_idx] = 1;
  schematic.adjacency[brick.schematic_idx][interactor_schematic_idx] = 1;
  schematic.start_atom_type[interactor_schematic_idx][brick.schematic_idx] = chosen_connection.start_atom_type;
  schematic.start_atom_type[brick.schematic_idx][interactor_schematic_idx] = compatible_connection.start_atom_type;
  schematic.end_atom_type[interactor_schematic_idx][brick.schematic_idx] = chosen_connection.end_atom_type;
  schematic.end_atom_type[brick.schematic_idx][interactor_schematic_idx] = compatible_connection.end_atom_type;
  schematic.bond_type[interactor_schematic_idx][brick.schematic_idx] = chosen_connection.bond_type;
  schematic.bond_type[brick.schematic_idx][interactor_schematic_idx] = compatible_connection.bond_type;
  schematic.start_atom_idx[interactor_schematic_idx][brick.schematic_idx] = immutable_reconstruction_interactor_idx;
  schematic.start_atom_idx[brick.schematic_idx][interactor_schematic_idx] = immutable_brick_interactor_idx;
  schematic.end_atom_idx[interactor_schematic_idx][brick.schematic_idx] = immutable_brick_interactor_idx;
  schematic.end_atom_idx[brick.schematic_idx][interactor_schematic_idx] = immutable_reconstruction_interactor_idx;

  // Update the rest of the ReconstructedMol's attributes.
  level += brick.level;
  if (add_ring) {
    n_ring_atoms += brick.level;
  };
  max_atom_idx = brick.max_atom_idx;
  bricks.insert({ brick.schematic_idx, brick });
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  // Signal function success.
  return true;
};

bool ReconstructedMol::InternalExpansion(sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, bool guided, bool update_stereo) {
  // In an internalExpansion a Pseudofragment is inserted between one of the
  // ReconstructedMol's MolBricks and a combination of the MolBrick's neighbors.

  // Retrieve the internal MolBricks.
  std::vector<MolBrick*> internal_bricks = GetInternalBricks();

  // If no internal bricks exist signal function failure.
  if (internal_bricks.empty()) {
    return false;
  };

  // Use each internal MolBrick to create an InsertionPoint.
  std::vector<InsertionPoint> insertion_points;
  insertion_points.reserve(internal_bricks.size());
  for (MolBrick* brick : internal_bricks) {
    // IMPROVE: Check to see if the number of combinations is
    // below the max multiple intersections.
    insertion_points.push_back(InsertionPoint(brick, prng, guided));
  };

  // Evaluate if the operation could increase the number of rings. Note that
  // in case it couldn't it doesn't mean that acyclic inserts exist either.
  bool rings_can_be_added = false;
  if (n_ring_atoms < controller.max_size) {
    for (InsertionPoint& ip : insertion_points) {
      ip.RetrieveRingInserts(query_results);
      if (ip.has_ring_inserts) {
        rings_can_be_added = true;
        break;
      };
    };
  };

  // Decide whether the number of rings should be increased or not.
  bool increase_rings = false;
  if (rings_can_be_added) {
    if (n_ring_atoms < controller.max_size) {
      increase_rings = controller.DecideChange(n_ring_atoms, prng);
    };
  };

  // Retrieve the pertinent type of Pseudofragment inserts.
  std::vector<InsertionPoint*> candidates;
  std::vector<float> candidates_weights;
  candidates.reserve(insertion_points.size());
  candidates_weights.reserve(insertion_points.size());
  if (increase_rings) {
    for (InsertionPoint& ip : insertion_points) {
      if (!ip.retrieved_ring_inserts) {
        ip.RetrieveRingInserts(query_results);
      };
      if (ip.has_ring_inserts) {
        ip.CalcRingBasedWeight();
        candidates.push_back(&ip);
        candidates_weights.push_back(ip.weight);
      };
    };
  } else {
    for (InsertionPoint& ip : insertion_points) {
      if (!ip.retrieved_acyclic_inserts) {
        ip.RetrieveAcyclicInserts(query_results);
      };
      if (ip.has_acyclic_inserts) {
        ip.CalcAcyclicBasedWeight();
        candidates.push_back(&ip);
        candidates_weights.push_back(ip.weight);
      };
    };
  };

  // If no insertable Pseudofragments were found signal function failure.
  if (candidates.empty()) {
    assert(!increase_rings);
    return false;
  };

  // Choose an InsertionPoint through weighted selection, and sample a
  // Pseudofragment insert from it.
  std::discrete_distribution<unsigned> distribution(candidates_weights.begin(), candidates_weights.end());
  InsertionPoint* insertion_point = candidates[distribution(prng)];
  MolBrick* owner = insertion_point->owner;
  unsigned owner_schematic_idx = owner->schematic_idx;
  std::tuple<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION, unsigned, float> sample;
  if (increase_rings) {
    sample = insertion_point->SampleRingInsert(prng);
  } else {
    sample = insertion_point->SampleAcyclicInsert(prng);
  };
  MOLBRICKS_COMBINATION& neighbors_combination = std::get<0>(sample);
  neighbors_combination.pop_back(); // Remove the last "neighbor", which is really the owner.
  const ConnectionsTable& owner_available_connections = owner->GetAvailableConnections(prng, neighbors_combination);

  // Retrieve the Pseudofragment with the specified ID.
  unsigned pseudofragment_id = std::get<2>(sample);
  Pseudofragment pseudofragment = GetPseudofragmentByIDFromDB(pseudofragment_id, select_statement);

  // Determine which ReconstructionSchematic index should be assigned to the
  // new MolBrick. If any indices became available due to a deletion, take one
  // of them. If not, assign a new one.
  unsigned insert_schematic_idx;
  if (available_schematic_idxs.empty()) {
    schematic.Expand();
    insert_schematic_idx = max_schematic_idx + 1;
    ++max_schematic_idx;
  } else {
    insert_schematic_idx = available_schematic_idxs.back();
    available_schematic_idxs.pop_back();
  };

  // Convert the Pseudofragment to a MolBrick.
  MolBrick insert(pseudofragment, pseudofragment_id, max_atom_idx + 1, insert_schematic_idx, std::get<3>(sample), this);

  // Reset the GeneticLog.
  if (guided) {
    genetic_log.Clear();
    genetic_log.SetImportantPseudofragmentID(insert.pseudofragment_id);
  };

  // Update the pseudomol by combining the atoms of the ReconstructedMol's
  // pseudomol and the MolBrick's pseudomol.
  pseudomol.insertMol(insert.pseudomol);

  // Determine how the MolBrick insert should connect to the MolBricks of the
  // InsertionPoint by finding a solution to the Maximum Bipartite Matching problem.

  // First, create maps correlating MolBricks/ConnectionPoints to integer IDs.
  // These will be the IDs of the BipartiteGraph's Vertices. Secondly, define
  // the insert MolBrick - Neighboring MolBrick Connections relationships.
  // These will be the Edges. Note that herein the InsertionPoint's owner MolBrick
  // is also considered a neighbor.
  // IMPROVE: Better BipartiteGraph construction.
  std::vector<std::pair<int, int>> edges;
  std::map<int, MolBrick*> vtx_ids_neighbors;
  std::map<MolBrick*, int> neighbors_vtx_ids;
  std::map<int, const Connection*> vtx_ids_connections;

  int id = 1;
  vtx_ids_neighbors.insert({ id, owner });
  neighbors_vtx_ids.insert({ owner, id });
  ++id;
  for (MolBrick* neighbor : neighbors_combination) {
    vtx_ids_neighbors.insert({ id, neighbor });
    neighbors_vtx_ids.insert({ neighbor, id });
    ++id;
  };
  for (const auto& c : insert.connections) {
    for (const auto& cp : c.second) {
      vtx_ids_connections.insert({ id, &c.first });
      if (owner_available_connections.IsCompatibleWith(c.first, query_results.compatibilities)) {
        edges.push_back(std::make_pair(1, id));
      };
      for (MolBrick* neighbor : neighbors_combination) {
        if (insertion_point->neighbors_available_connections[neighbor].IsCompatibleWith(c.first, query_results.compatibilities)) {
          edges.push_back(std::make_pair(neighbors_vtx_ids[neighbor], id));
        };
      };
      ++id;
    };
  };

  // Run the Hopcroft-Karp algorithm to find a solution to the Maximum
  // Bipartite Matching problem.
  BipartiteGraph graph(vtx_ids_neighbors.size(), vtx_ids_connections.size(), edges);
  PATH matching = graph.HopcroftKarp(prng);

  // Sanity check. Make sure that the cardinality of the matching
  // corresponds to the number of neighboring MolBricks. Otherwise
  // the substitute can't possibly be valid.
  assert(matching.size() == neighbors_combination.size() + 1);

  // Connect the MolBrick insert to the InsertionPoint's owner MolBrick's neighbors
  // according to the aforementioned matching.
  Connection original_forward_connection, original_backward_connection, new_forward_connection, new_backward_connection;
  ConnectionsTable insert_available_connections = insert.connections;
  for (MolBrick* neighbor : neighbors_combination) {
    unsigned neighbor_schematic_idx = neighbor->schematic_idx;
    ConnectionsTable& neighbor_available_connections = insertion_point->neighbors_available_connections[neighbor];

    // Remove the original Connection between the owner and the neighbor.
    // and break the bond between the owner and the neighbor.
    original_forward_connection = Connection(schematic.start_atom_type[owner_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[owner_schematic_idx][neighbor_schematic_idx], schematic.bond_type[owner_schematic_idx][neighbor_schematic_idx]);
    original_backward_connection = Connection(schematic.start_atom_type[neighbor_schematic_idx][owner_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][owner_schematic_idx], schematic.bond_type[neighbor_schematic_idx][owner_schematic_idx]);
    unsigned owner_interactor_idx = schematic.start_atom_idx[owner_schematic_idx][neighbor_schematic_idx];
    unsigned original_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][owner_schematic_idx];
    internal_connections.MoveConnectionTo(original_forward_connection, owner_interactor_idx, connections);
    internal_connections.MoveConnectionTo(original_backward_connection, original_neighbor_interactor_idx, connections);
    pseudomol.removeBond(GetMutableAtomIdx(pseudomol, owner_interactor_idx), GetMutableAtomIdx(pseudomol, original_neighbor_interactor_idx));

    // Record the deleted owner-to-neighbor and neighbor-to-owner in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_forward_connection, owner_interactor_idx, neighbor->pseudofragment_id);
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, original_neighbor_interactor_idx, owner->pseudofragment_id);
    };

    // Add the Connections between the MolBrick insert and the neighboring MolBrick.
    int vtx_id = neighbors_vtx_ids[neighbor];
    PATH::const_iterator it = std::find_if(matching.begin(), matching.end(),
      [=](const Edge* edge) {
        return edge->GetStart()->GetID() == vtx_id;
      });
    assert(it != matching.end());
    new_forward_connection = *vtx_ids_connections[(*it)->GetEnd()->GetID()];
    new_backward_connection = neighbor_available_connections.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

    const CONNECTION_POINT_VECTOR& insert_cpoint_vector = insert_available_connections[new_forward_connection];
    unsigned n_connections = insert_cpoint_vector.size();
    assert(n_connections > 0);
    unsigned connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
        connection_idx = cpoint_distribution(prng);
      };
    };
    unsigned insert_interactor_idx = insert_cpoint_vector[connection_idx].atom_idx;

    const CONNECTION_POINT_VECTOR& neighbor_cpoint_vector = neighbor_available_connections[new_backward_connection];
    n_connections = neighbor_cpoint_vector.size();
    assert(n_connections > 0);
    connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
        connection_idx = cpoint_distribution(prng);
      };
    };
    unsigned new_neighbor_interactor_idx = neighbor_cpoint_vector[connection_idx].atom_idx;

    ConnectionPoint* cpoint = internal_connections.AddConnection(new_forward_connection, insert_interactor_idx);
    if (guided) {
      query_results.AssignFreshWeights(new_forward_connection, cpoint);
    };
    insert_available_connections.RemoveConnection(new_forward_connection, insert_interactor_idx);
    connections.MoveConnectionTo(new_backward_connection, new_neighbor_interactor_idx, internal_connections);

    // Record the new neighbor-to-insert link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, new_neighbor_interactor_idx, insert.pseudofragment_id);
    };

    // Clear the ReconstructionsSchematic's entries corresponding to
    // the InsertionPoint's owner MolBrick and the neighboring MolBrick.
    schematic.ClearCell(owner_schematic_idx, neighbor_schematic_idx);
    schematic.ClearCell(neighbor_schematic_idx, owner_schematic_idx);

    // Update the ReconstructionsSchematic's entries corresponding to
    // the the insert-neighbor Connection.
    schematic.adjacency[insert_schematic_idx][neighbor_schematic_idx] = 1;
    schematic.adjacency[neighbor_schematic_idx][insert_schematic_idx] = 1;
    schematic.start_atom_type[insert_schematic_idx][neighbor_schematic_idx] = new_forward_connection.start_atom_type;
    schematic.start_atom_type[neighbor_schematic_idx][insert_schematic_idx] = new_backward_connection.start_atom_type;
    schematic.end_atom_type[insert_schematic_idx][neighbor_schematic_idx] = new_forward_connection.end_atom_type;
    schematic.end_atom_type[neighbor_schematic_idx][insert_schematic_idx] = new_backward_connection.end_atom_type;
    schematic.bond_type[insert_schematic_idx][neighbor_schematic_idx] = new_forward_connection.bond_type;
    schematic.bond_type[neighbor_schematic_idx][insert_schematic_idx] = new_backward_connection.bond_type;
    schematic.start_atom_idx[insert_schematic_idx][neighbor_schematic_idx] = insert_interactor_idx;
    schematic.start_atom_idx[neighbor_schematic_idx][insert_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[insert_schematic_idx][neighbor_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[neighbor_schematic_idx][insert_schematic_idx] = insert_interactor_idx;

    // Create the Bond between the MolBrick insert and the neighboring MolBrick.
    RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
    pseudomol.addBond(GetMutableAtomIdx(pseudomol, insert_interactor_idx), GetMutableAtomIdx(pseudomol, new_neighbor_interactor_idx), bond_type);
  };

  // For the InsertionPoint's owner MolBrick only adding Connections and
  // updating the ReconstructionsSchematic is necessary. Updating the owner
  // last is important to free up Connections first.
  PATH::const_iterator it = std::find_if(matching.begin(), matching.end(),
    [=](const Edge* edge) {
      return edge->GetStart()->GetID() == 1;
    });
  assert(it != matching.end());
  new_forward_connection = *vtx_ids_connections[(*it)->GetEnd()->GetID()];
  new_backward_connection = owner_available_connections.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

  const CONNECTION_POINT_VECTOR& insert_cpoint_vector = insert_available_connections[new_forward_connection];
  unsigned n_connections = insert_cpoint_vector.size();
  assert(n_connections > 0);
  unsigned connection_idx = 0;
  if (n_connections > 1) {
    if (n_connections <= max_size_uniform_int_distributions) {
      connection_idx = uniform_int_distributions[n_connections](prng);
    } else {
      std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
      connection_idx = cpoint_distribution(prng);
    };
  };
  unsigned insert_interactor_idx = insert_cpoint_vector[connection_idx].atom_idx;

  const CONNECTION_POINT_VECTOR& owner_cpoint_vector = owner_available_connections.at(new_backward_connection);
  n_connections = owner_cpoint_vector.size();
  assert(n_connections > 0);
  connection_idx = 0;
  if (n_connections > 1) {
    if (n_connections <= max_size_uniform_int_distributions) {
      connection_idx = uniform_int_distributions[n_connections](prng);
    } else {
      std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
      connection_idx = cpoint_distribution(prng);
    };
  };
  unsigned owner_interactor_idx = owner_cpoint_vector[connection_idx].atom_idx;

  ConnectionPoint* cpoint = internal_connections.AddConnection(new_forward_connection, insert_interactor_idx);
  if (guided) {
    query_results.AssignFreshWeights(new_forward_connection, cpoint);
  };
  insert_available_connections.RemoveConnection(new_forward_connection, insert_interactor_idx);
  connections.MoveConnectionTo(new_backward_connection, owner_interactor_idx, internal_connections);

  // Record the new owner-to-insert link in the GeneticLog.
  if (guided) {
    genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, owner_interactor_idx, insert.pseudofragment_id);
  };

  // Update the ReconstructionsSchematic's entries corresponding to
  // the the insert-neighbor Connection.
  schematic.adjacency[insert_schematic_idx][owner_schematic_idx] = 1;
  schematic.adjacency[owner_schematic_idx][insert_schematic_idx] = 1;
  schematic.start_atom_type[insert_schematic_idx][owner_schematic_idx] = new_forward_connection.start_atom_type;
  schematic.start_atom_type[owner_schematic_idx][insert_schematic_idx] = new_backward_connection.start_atom_type;
  schematic.end_atom_type[insert_schematic_idx][owner_schematic_idx] = new_forward_connection.end_atom_type;
  schematic.end_atom_type[owner_schematic_idx][insert_schematic_idx] = new_backward_connection.end_atom_type;
  schematic.bond_type[insert_schematic_idx][owner_schematic_idx] = new_forward_connection.bond_type;
  schematic.bond_type[owner_schematic_idx][insert_schematic_idx] = new_backward_connection.bond_type;
  schematic.start_atom_idx[insert_schematic_idx][owner_schematic_idx] = insert_interactor_idx;
  schematic.start_atom_idx[owner_schematic_idx][insert_schematic_idx] = owner_interactor_idx;
  schematic.end_atom_idx[insert_schematic_idx][owner_schematic_idx] = owner_interactor_idx;
  schematic.end_atom_idx[owner_schematic_idx][insert_schematic_idx] = insert_interactor_idx;

  // Create the Bond between the MolBrick insert and the owner MolBrick.
  RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
  pseudomol.addBond(GetMutableAtomIdx(pseudomol, insert_interactor_idx), GetMutableAtomIdx(pseudomol, owner_interactor_idx), bond_type);

  // Add the remainining insert MolBrick's Connections as peripheral connections.
  if (guided) {
    std::vector<ConnectionPoint*> cpoints;
    for (auto& c : insert_available_connections) {
      cpoints = connections.AddConnections(c.first, c.second);
      query_results.AssignFreshWeights(c.first, cpoints);
    };
    assert(AllConnectionsHaveWeights());
  } else {
    for (auto& c : insert_available_connections) {
      connections.AddConnections(c.first, c.second);
    };
  };

  // Update the rest of the ReconstructedMol's attributes.
  level += insert.level;
  if (increase_rings) {
    n_ring_atoms += insert.level;
  };
  max_atom_idx = insert.max_atom_idx;
  bricks.insert({ insert_schematic_idx, insert });
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  // Signal function success.
  return true;
};


bool ReconstructedMol::PeripheralDeletion(SizeController & controller, std::mt19937 & prng, bool guided, bool update_stereo) {
  // In a peripheral deletion one of the peripheral MolBricks of the
  // ReconstructedMol (that is, a MolBrick bound through a single Connection
  // to the rest of the ReconstructedMol) is deleted.

  // If the ReconstructedMol is made up of a single MolBrick a deletion isn't
  // allowed, and function failure is signalled.
  if (bricks.size() == 1) {
    return false;
  };

  // Retrieve the peripheral MolBricks that are candidates to be deleted.
  std::vector<MolBrick*> peripheral_bricks = GetPeripheralBricks();

  // Evaluate the type of the peripheral MolBricks.
  // Possible return values:
  //  -1 = Consists of only acyclic MolBricks.
  //   0 = Consists of both acyclic and cyclic MolBricks (including ring_parts).
  //   1 = Consists of only cyclic MolBricks.
  int evaluation = EvaluateBrickSet(peripheral_bricks);

  // If the peripheral bricks include both cyclic and acyclic instances,
  // decide if the number of rings should be reduced. Otherwise, remove
  // whatever type of brick is available.
  bool delete_ring = false;
  if (evaluation == -1) {
    delete_ring = false;
  } else if (evaluation == 0) {
    delete_ring = controller.DecideChange(n_ring_atoms, prng);
  } else if (evaluation == 1) {
    delete_ring = true;
  };

  bool has_ring, ring_part;
  if (delete_ring) {
    has_ring = true;
    ring_part = false;
  } else {
    has_ring = false;
    ring_part = false;
  };

  // Perform a weighted selection of the desired type of MolBrick. The selection
  // is weighted for "fairness" reasons, since the expansion selections are
  // weighted. If  deletions weren't weighted there would be an influx/outflux
  // imbalance for less frequent MolBricks.
  MolBrick* brick = ChooseBrickFromSet(peripheral_bricks, has_ring, ring_part, prng);

  // Get the mutable atom indices of the MolBrick to be deleted.
  std::vector<unsigned> brick_mutable_atom_indices = GetMutableAtomIndices(pseudomol, brick->atom_indices);

  // Remove the atoms in descending sorted order to avoid atom index renumbering.
  std::sort(brick_mutable_atom_indices.begin(), brick_mutable_atom_indices.end(), std::greater<unsigned>());
  for (auto idx : brick_mutable_atom_indices) {
    pseudomol.removeAtom(idx);
  };

  // Determine which MolBrick is the neighbor to the one to be deleted. Since
  // the MolBrick is peripheral it can only have a single neighbor.
  unsigned neighbor_schematic_idx = 0;
  for (unsigned adjacency_bool : schematic.adjacency[brick->schematic_idx]) {
    if (adjacency_bool == 1) {
      break;
    };
    ++neighbor_schematic_idx;
  };

  // Update the ReconstructedMol's connection tables. This entails:
  //  (1) Moving the Connections involved in the deleted MolBrick-neighbor
  //      interaction from the internal to the peripheral ConnectionTable.
  //  (2) Removing the MolBrick's Connections from the peripheral ConnectionTable.
  Connection forward_connection(schematic.start_atom_type[neighbor_schematic_idx][brick->schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][brick->schematic_idx], schematic.bond_type[neighbor_schematic_idx][brick->schematic_idx]);
  Connection backward_connection(schematic.start_atom_type[brick->schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[brick->schematic_idx][neighbor_schematic_idx], schematic.bond_type[brick->schematic_idx][neighbor_schematic_idx]);
  unsigned immutable_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][brick->schematic_idx];
  unsigned immutable_brick_interactor_idx = schematic.start_atom_idx[brick->schematic_idx][neighbor_schematic_idx];

  internal_connections.MoveConnectionTo(forward_connection, immutable_neighbor_interactor_idx, connections);
  internal_connections.MoveConnectionTo(backward_connection, immutable_brick_interactor_idx, connections);
  for (const auto& c : brick->connections) {
    connections.RemoveConnections(c.first, c.second);
  };

  // Reset the GeneticLog and record the deleted link in it.
  if (guided) {
    genetic_log.Clear();
    genetic_log.SetImportantPseudofragmentID(brick->pseudofragment_id);
    genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, forward_connection, immutable_neighbor_interactor_idx, brick->pseudofragment_id);
  };

  // Update the ReconstructionSchematic. In a deletion this means clearing
  // the values of the ReconstructionSchematic's row and column corresponding
  // to the MolBrick in question. Note that the matrix itself isn't resized,
  // meaning that the blank row and column are made available for future MolBrick's.
  schematic.ClearRow(brick->schematic_idx);
  schematic.ClearColumn(brick->schematic_idx);
  // Store the deleted MolBrick's schematic index for future use by new MolBricks.
  available_schematic_idxs.push_back(brick->schematic_idx);

  // Update the rest of the ReconstructedMol's attributes.
  level -= brick->level;
  if (delete_ring) {
    n_ring_atoms -= brick->level;
  };
  bricks.erase(brick->schematic_idx);
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  // Signal function success.
  return true;
};


bool ReconstructedMol::InternalDeletion(ConnectionQueryResults & query_results, SizeController & controller, std::mt19937 & prng, bool guided, bool update_stereo) {
  // In an internal deletion one of the internal MolBricks of the ReconstructedMol
  // (that is, a MolBrick bound to multiple neighboring MolBricks) is deleted.
  // To do so one of its neighbors must be able to connect to all of the other
  // neighbors of the deleted MolBrick.

  // Retrieve the internal MolBricks. If none exist signal function failure.
  std::vector<MolBrick*> internal_bricks = GetInternalBricks();
  if (internal_bricks.empty()) {
    return false;
  };

  // Loop over the internal MolBricks and convert them into DeletionPoints.
  std::vector<DeletionPoint> deletion_points;
  deletion_points.reserve(internal_bricks.size());
  for (MolBrick* brick : internal_bricks) {
    DeletionPoint deletion_point(brick);
    // Evaluate whether the DeletionPoint pertains to a deletable MolBrick.
    // If it does, store it.
    deletion_point.Evaluate(query_results.compatibilities, prng);
    if (deletion_point.is_deletable) {
      deletion_points.push_back(std::move(deletion_point));
    };
  };

  // If none of the internal MolBricks are deletable signal function failure.
  if (deletion_points.empty()) {
    return false;
  };

  // Evaluate how the InternalDeletion could modulate the number of rings.
  bool rings_can_be_kept = false, rings_can_be_removed = false;
  for (DeletionPoint& dp : deletion_points) {
    if (dp.has_ring) {
      rings_can_be_removed = true;
    } else {
      rings_can_be_kept = true;
    };
  };

  // Decide how the InternalDeletion should modulate the number of rings.
  bool decrease_rings = false;
  if (rings_can_be_kept && rings_can_be_removed) {
    if (controller.DecideChange(n_ring_atoms, prng)) {
      decrease_rings = true;
    };
  } else if (rings_can_be_removed) {
    decrease_rings = true;
  } else {
    decrease_rings = false;
  };

  // Choose a DeletionPoint (i.e. MolBrick to be deleted)
  // through weighted selection.
  std::vector<DeletionPoint*> deletion_point_candidates;
  std::vector<float> deletion_point_weights;
  if (decrease_rings) {
    for (DeletionPoint& dp : deletion_points) {
      if (dp.has_ring) {
        deletion_point_candidates.push_back(&dp);
        deletion_point_weights.push_back(dp.weight);
      };
    };
  } else {
    for (DeletionPoint& dp : deletion_points) {
      if (!dp.has_ring) {
        deletion_point_candidates.push_back(&dp);
        deletion_point_weights.push_back(dp.weight);
      };
    };
  };
  std::discrete_distribution<unsigned> dp_distribution(deletion_point_weights.begin(), deletion_point_weights.end());
  DeletionPoint* deletion_point = deletion_point_candidates[dp_distribution(prng)];
  MolBrick* deleted = deletion_point->deletable;
  unsigned deleted_schematic_idx = deleted->schematic_idx;

  // Reset the GeneticLog.
  if (guided) {
    genetic_log.Clear();
    genetic_log.SetImportantPseudofragmentID(deleted->pseudofragment_id);
  };

  // Generate the remaining matchings for the chosen DeletionPoint and randomly
  // choose one of them.
  deletion_point->GenerateMatchings(query_results.compatibilities, prng);
  unsigned idx = 0, n_matchings = deletion_point->other_neighbors_matchings.size();
  if (n_matchings > 1) {
    if (n_matchings <= max_size_uniform_int_distributions) {
      idx = uniform_int_distributions[n_matchings](prng);
    } else {
      std::uniform_int_distribution<unsigned> matchings_distribution(0, n_matchings - 1);
      idx = matchings_distribution(prng);
    };
  };
  std::map<const MolBrick*, std::map<const MolBrick*, const Connection*>>::iterator it = deletion_point->other_neighbors_matchings.begin();
  std::advance(it, idx);

  // Retrieve the MolBrick that will become the "core" of the new Connections
  // (i.e. the one that connects to all the remaining neighboring MolBrick).
  const MolBrick* core = it->first;
  unsigned core_schematic_idx = core->schematic_idx;
  assert(deletion_point->neighbors_available_connections[core].n_connection_points >= deletion_point->n_other_neighbors);

  // Remove the atoms of the substituted MolBrick in sorted descending order
  // to avoid atom reindexing.
  std::vector<unsigned> deleted_mutable_atom_indices = GetMutableAtomIndices(pseudomol, deleted->atom_indices);
  std::sort(deleted_mutable_atom_indices.begin(), deleted_mutable_atom_indices.end(), std::greater<unsigned>());
  for (unsigned atom_idx : deleted_mutable_atom_indices) {
    pseudomol.removeAtom(atom_idx);
  };

  // Update the ReconstructedMol's Connections and ReconstructionSchematic.
  const MolBrick* neighbor;
  unsigned neighbor_schematic_idx, n_connections;
  unsigned deleted_interactor_idx, original_neighbor_interactor_idx, core_interactor_idx, new_neighbor_interactor_idx;
  Connection original_forward_connection, original_backward_connection, new_forward_connection, new_backward_connection;
  ConnectionsTable core_available_connections = deletion_point->neighbors_available_connections[core];

  // First remove the Connection between the deleted MolBrick and the new core MolBrick.
  deleted_interactor_idx = schematic.start_atom_idx[deleted_schematic_idx][core_schematic_idx];
  core_interactor_idx = schematic.start_atom_idx[core_schematic_idx][deleted_schematic_idx];
  original_forward_connection = Connection(schematic.start_atom_type[deleted_schematic_idx][core_schematic_idx], schematic.end_atom_type[deleted_schematic_idx][core_schematic_idx], schematic.bond_type[deleted_schematic_idx][core_schematic_idx]);
  original_backward_connection = Connection(schematic.start_atom_type[core_schematic_idx][deleted_schematic_idx], schematic.end_atom_type[core_schematic_idx][deleted_schematic_idx], schematic.bond_type[core_schematic_idx][deleted_schematic_idx]);
  internal_connections.RemoveConnection(original_forward_connection, deleted_interactor_idx);
  internal_connections.MoveConnectionTo(original_backward_connection, core_interactor_idx, connections);

  // Record the deleted core-to-deleted link in the GeneticLog.
  if (guided) {
    genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, core_interactor_idx, deleted->pseudofragment_id);
  };

  // Thereafter do the same for all the other neighbors of the deleted MolBrick.
  for (const auto& match : it->second) {
    neighbor = match.first;
    neighbor_schematic_idx = neighbor->schematic_idx;
    ConnectionsTable& neighbor_available_connections = deletion_point->neighbors_available_connections[neighbor];

    // Retrieve the original Connection/ConnectionPoints involved.
    deleted_interactor_idx = schematic.start_atom_idx[deleted_schematic_idx][neighbor_schematic_idx];
    original_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][deleted_schematic_idx];
    original_forward_connection = Connection(schematic.start_atom_type[deleted_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[deleted_schematic_idx][neighbor_schematic_idx], schematic.bond_type[deleted_schematic_idx][neighbor_schematic_idx]);
    original_backward_connection = Connection(schematic.start_atom_type[neighbor_schematic_idx][deleted_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][deleted_schematic_idx], schematic.bond_type[neighbor_schematic_idx][deleted_schematic_idx]);

    // Remove the original Connections.
    internal_connections.RemoveConnection(original_forward_connection, deleted_interactor_idx);
    internal_connections.MoveConnectionTo(original_backward_connection, original_neighbor_interactor_idx, connections);

    // Record the deleted neighbor-to-deleted link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, original_neighbor_interactor_idx, deleted->pseudofragment_id);
    };

    // Select random ConnectionPoints for the Maximum Bipartite Matching solution.
    new_forward_connection = *match.second;
    new_backward_connection = neighbor_available_connections.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

    const CONNECTION_POINT_VECTOR& core_cpoint_vector = core_available_connections[new_forward_connection];
    n_connections = core_cpoint_vector.size();
    assert(n_connections > 0);
    idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
        idx = cpoint_distribution(prng);
      };
    };
    core_interactor_idx = core_cpoint_vector[idx].atom_idx;

    const CONNECTION_POINT_VECTOR& neighbor_cpoint_vector = neighbor_available_connections[new_backward_connection];
    n_connections = neighbor_cpoint_vector.size();
    assert(n_connections > 0);
    idx = 0;
    if (n_connections > 1) {
      idx = uniform_int_distributions[n_connections](prng);
      if (n_connections <= max_size_uniform_int_distributions) {
        idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
        idx = cpoint_distribution(prng);
      };
    };
    new_neighbor_interactor_idx = neighbor_cpoint_vector[idx].atom_idx;

    // Add the new Connections.
    connections.MoveConnectionTo(new_forward_connection, core_interactor_idx, internal_connections);
    connections.MoveConnectionTo(new_backward_connection, new_neighbor_interactor_idx, internal_connections);
    core_available_connections.RemoveConnection(new_forward_connection, core_interactor_idx);

    // Record the new core-to-neighbor and neighbor-to-core links in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_forward_connection, core_interactor_idx, neighbor->pseudofragment_id);
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, new_neighbor_interactor_idx, core->pseudofragment_id);
    };

    // Update the ReconstructionsSchematic.
    schematic.adjacency[core_schematic_idx][neighbor_schematic_idx] = 1;
    schematic.adjacency[neighbor_schematic_idx][core_schematic_idx] = 1;
    schematic.start_atom_type[core_schematic_idx][neighbor_schematic_idx] = new_forward_connection.start_atom_type;
    schematic.start_atom_type[neighbor_schematic_idx][core_schematic_idx] = new_backward_connection.start_atom_type;
    schematic.end_atom_type[core_schematic_idx][neighbor_schematic_idx] = new_forward_connection.end_atom_type;
    schematic.end_atom_type[neighbor_schematic_idx][core_schematic_idx] = new_backward_connection.end_atom_type;
    schematic.bond_type[core_schematic_idx][neighbor_schematic_idx] = new_forward_connection.bond_type;
    schematic.bond_type[neighbor_schematic_idx][core_schematic_idx] = new_backward_connection.bond_type;
    schematic.start_atom_idx[core_schematic_idx][neighbor_schematic_idx] = core_interactor_idx;
    schematic.start_atom_idx[neighbor_schematic_idx][core_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[core_schematic_idx][neighbor_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[neighbor_schematic_idx][core_schematic_idx] = core_interactor_idx;

    // Create the new bonds.
    unsigned mutable_core_interactor_idx = GetMutableAtomIdx(pseudomol, core_interactor_idx);
    unsigned mutable_neighbor_interactor_idx = GetMutableAtomIdx(pseudomol, new_neighbor_interactor_idx);
    RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
    pseudomol.addBond(mutable_core_interactor_idx, mutable_neighbor_interactor_idx, bond_type);
  };

  // Finish updating the ConnectionsTables by adjusting the peripheral connections.
  for (const auto& c : deleted->GetAvailableConnections(prng)) {
    connections.RemoveConnections(c.first, c.second);
  };

  // Update the rest of the ReconstructionSchematic by clearing the data
  // pertaining to the deleted MolBrick.
  schematic.ClearRow(deleted_schematic_idx);
  schematic.ClearColumn(deleted_schematic_idx);
  available_schematic_idxs.push_back(deleted_schematic_idx);

  // Update the rest of the ReconstructedMol's attributes.
  level -= deleted->level;
  if (decrease_rings) {
    n_ring_atoms -= deleted->level;
  };
  bricks.erase(deleted_schematic_idx);
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  // Signal function success.
  return true;
};


bool ReconstructedMol::Substitution(const std::string & location, sqlite3_stmt * select_statement, ConnectionQueryResults & query_results, SizeController & controller, std::mt19937 & prng, bool guided, bool update_stereo) {
  // In a substitution a peripheral or internal MolBrick is substituted by
  // a Pseudofragment that has atleast the Connections involved in binding
  // the substituted MolBrick to its neighbors.

  // Currently Substitutions aren't allowed when the ReconstructedMol consists
  // of a single MolBrick. While the Substitution of the only MolBrick for
  // another random one would be possible it would imply completely disregarding
  // the starting base ReconstructedMol, which doesn't truely fit the definition
  // of an evolutive process.
  if (bricks.size() <= 1) {
    return false;
  };

  // Retrieve the ReconstructedMol's MolBricks of the type to be substituted
  // (i.e. peripheral or internal).
  std::vector<MolBrick*> bricks_to_evaluate;
  unsigned n_bricks;
  if (location == "peripheral") {
    bricks_to_evaluate = GetPeripheralBricks();
    n_bricks = bricks_to_evaluate.size();
  } else if (location == "internal") {
    bricks_to_evaluate = GetInternalBricks();
    n_bricks = bricks_to_evaluate.size();
    // If there are no internal MolBricks in the ReconstructedMol signal function
    // failure. I could default to a peripheral substitution in this case, but
    // it might cloud operation probabilities during development.
    if (n_bricks == 0) {
      return false;
    };
  } else {
    throw std::runtime_error("Invalid location option for Substitution.");
  };

  // Convert the subject MolBricks into SubstitutionPoints.
  std::vector<SubstitutionPoint> substitution_points;
  substitution_points.reserve(n_bricks);
  for (MolBrick* brick : bricks_to_evaluate) {
    substitution_points.push_back(SubstitutionPoint(brick, guided));
  };

  // Loop over the identified SubstitutionPoints and evaluate how
  // the number of rings can be modulated.
  bool substitutes_exist = false, rings_can_be_added = false, rings_can_be_removed = false, rings_can_be_kept_constant = false;
  for (SubstitutionPoint& sp : substitution_points) {
    sp.EvaluateSubstitutability(query_results, prng);
    if (sp.has_substitutes) {
      if (sp.rings_can_be_added) {
        if (n_ring_atoms < controller.max_size) {
          rings_can_be_added = true;
          substitutes_exist = true;
        };
      };
      if (sp.rings_can_be_removed) {
        rings_can_be_removed = true;
        substitutes_exist = true;
      };
      if (sp.rings_can_be_kept_constant) {
        rings_can_be_kept_constant = true;
        substitutes_exist = true;
      };
      if (rings_can_be_added && rings_can_be_removed && rings_can_be_kept_constant) {
        break;
      };
    };
  };

  // If no substitute Pseudofragments exist, signal function failure. This
  // may happen since due to the mismatch in the compatibility definition
  // stringency between genetic operators (e.g. PeripheralExpansion uses
  // a non-strict definition while Substitution uses a strict definition).
  if (!substitutes_exist) {
    return false;
  };

  // If both kinds of Pseudofragments can be used as substitutes,
  // decide which one to use. Otherwise, use the only type that can be used.
  bool constant_rings = false, increase_rings = false, decrease_rings = false;
  if (rings_can_be_kept_constant) {
    // If the number of rings can be kept constant, increased or decreased:
    if (rings_can_be_added && rings_can_be_removed) {
      int decision = controller.Decide(n_ring_atoms, prng);
      if (decision == 0) {
        constant_rings = true;
      } else if (decision == 1) {
        increase_rings = true;
      } else {
        decrease_rings = true;
      };
      // If the number of rings can be kept constant or increased:
    } else if (rings_can_be_added) {
      bool change = false;
      if (n_ring_atoms < controller.max_size) {
        change = controller.DecideChange(n_ring_atoms, prng);
      };
      if (change) {
        increase_rings = true;
      } else {
        constant_rings = true;
      };
      // If the number of rings can be kept constant or decreased:
    } else if (rings_can_be_removed) {
      bool change = controller.DecideChange(n_ring_atoms, prng);
      if (change) {
        decrease_rings = true;
      } else {
        constant_rings = true;
      };
      // If the number of rings can only be kept constant:
    } else {
      constant_rings = true;
    };
  } else if (rings_can_be_added) {
    // If the number of rings can be increased or decreased:
    if (rings_can_be_removed) {
      bool grow = controller.DecideGrowth(n_ring_atoms, prng);
      if (grow) {
        increase_rings = true;
      } else {
        decrease_rings = true;
      };
      // If the number of rings can only be increased:
    } else {
      increase_rings = true;
    };
    // If the number of rings can only be decreased:
  } else {
    decrease_rings = true;
  };

  // Loop over the subset of SubstitutionPoints and attempt to retrieve
  // substitute Pseudofragments for them that modify the ReconstructedMol's
  // number of rings in the desired way. If substitutes are found, calculate
  // the weight of the SubstitutionPoint.
  unsigned n_substitution_points = substitution_points.size();
  std::vector<SubstitutionPoint*> substitition_point_candidates;
  std::vector<float> substitution_point_weights;
  substitition_point_candidates.reserve(n_substitution_points);
  substitution_point_weights.reserve(n_substitution_points);
  // If the number of rings should be kept constant
  if (constant_rings) {
    for (SubstitutionPoint& sp : substitution_points) {
      // and the SubstitutionPoint pertains to a ring MolBrick
      if (sp.has_ring) {
        // retrieve ring substitutes.
        if (!sp.retrieved_ring_substitutes) {
          sp.RetrieveRingSubstitutes(query_results, prng);
        };
        if (sp.has_ring_substitutes) {
          sp.CalcRingBasedWeight();
          substitition_point_candidates.push_back(&sp);
          substitution_point_weights.push_back(sp.weight);
        };
      // and the SubstitutionPoint pertains to an acyclic MolBrick
      } else {
        // retrieve acyclic substitutes.
        if (!sp.retrieved_acyclic_substitutes) {
          sp.RetrieveAcyclicSubstitutes(query_results, prng);
        };
        if (sp.has_acyclic_substitutes) {
          sp.CalcAcyclicBasedWeight();
          substitition_point_candidates.push_back(&sp);
          substitution_point_weights.push_back(sp.weight);
        };
      };
    };
  // If the number of rings should be increased
  } else if (increase_rings) {
    for (SubstitutionPoint& sp : substitution_points) {
      // and the SubstitutionPoint pertains to an acyclic MolBrick
      if (!sp.has_ring) {
        // retrieve ring substitutes.
        if (!sp.retrieved_ring_substitutes) {
          sp.RetrieveRingSubstitutes(query_results, prng);
        };
        if (sp.has_ring_substitutes) {
          sp.CalcRingBasedWeight();
          substitition_point_candidates.push_back(&sp);
          substitution_point_weights.push_back(sp.weight);
        };
      };
    };
  // If the number of rings should be decreased
  } else if (decrease_rings) {
    for (SubstitutionPoint& sp : substitution_points) {
      // and the SubstitutionPoint pertains to a ring MolBrick
      if (sp.has_ring) {
        // retrieve acyclic substitutes.
        if (!sp.retrieved_acyclic_substitutes) {
          sp.RetrieveAcyclicSubstitutes(query_results, prng);
        };
        if (sp.has_acyclic_substitutes) {
          sp.CalcAcyclicBasedWeight();
          substitition_point_candidates.push_back(&sp);
          substitution_point_weights.push_back(sp.weight);
        };
      };
    };
  };

  // Choose a SubstitutionPoint (i.e. MolBrick to be substituted)
  // through weighted selection.
  assert(!substitition_point_candidates.empty());
  std::discrete_distribution<unsigned> sp_distribution(substitution_point_weights.begin(), substitution_point_weights.end());
  unsigned substitution_point_idx = sp_distribution(prng);
  SubstitutionPoint* substitution_point = substitition_point_candidates[substitution_point_idx];
  MolBrick* substituted = substitution_point->substitutable;
  unsigned brick_schematic_idx = substituted->schematic_idx;

  // Choose a MolBrick compatible with the aforementioned SubstitutionPoint
  // through weighted selection.
  std::tuple<CONNECTIONS_COMBINATION, unsigned, float> sample;
  if (constant_rings) {
    if (substitution_point->has_ring) {
      sample = substitution_point->SampleRingSubstitute(prng);
    } else {
      sample = substitution_point->SampleAcyclicSubstitute(prng);
    };
  } else if (increase_rings) {
    sample = substitution_point->SampleRingSubstitute(prng);
  } else if (decrease_rings) {
    sample = substitution_point->SampleAcyclicSubstitute(prng);
  };
  unsigned pseudofragment_id = std::get<1>(sample);
  float pseudofragment_weight = std::get<2>(sample);

  // Retrieve the Pseudofragment with the specified ID and convert it to a MolBrick.
  Pseudofragment pseudofragment = GetPseudofragmentByIDFromDB(pseudofragment_id, select_statement);
  MolBrick substitute(pseudofragment, pseudofragment_id, max_atom_idx + 1, brick_schematic_idx, pseudofragment_weight, this);

  // Reset the GeneticLog.
  if (guided) {
    genetic_log.Clear();
    genetic_log.SetImportantPseudofragmentID(substituted->pseudofragment_id);
    genetic_log.SetImportantPseudofragmentID(substitute.pseudofragment_id);
  };

  // Update the pseudomol. This entails removing the atoms of the MolBrick
  // to be substituted, adding the atoms of the new substitute MolBrick and
  // re-establishing the bonds between them.

  // Remove the atoms of the substituted MolBrick in sorted descending order
  // to avoid atom reindexing.
  std::vector<unsigned> substituted_mutable_atom_indices = GetMutableAtomIndices(pseudomol, substituted->atom_indices);
  std::sort(substituted_mutable_atom_indices.begin(), substituted_mutable_atom_indices.end(), std::greater<unsigned>());
  for (unsigned idx : substituted_mutable_atom_indices) {
    pseudomol.removeAtom(idx);
  };

  // Add the atoms of the substitute MolBrick.
  pseudomol.insertMol(substitute.pseudomol);

  // Determine how the substitute MolBrick should connect to the neighboring
  // MolBricks by finding a solution to the Maximum Bipartite Matching problem.

  // First, create maps correlating MolBricks/ConnectionPoints to integer IDs.
  // These will be the IDs of the BipartiteGraph's Vertices.
  // Secondly, define the Substitute MolBrick - Neighboring MolBrick
  // Connections relationships. These will be the Edges.
  // IMPROVE: Better BipartiteGraph construction.
  std::vector<std::pair<int, int>> edges;
  std::map<int, const MolBrick*> vtx_ids_neighbors;
  std::map<const MolBrick*, int> neighbors_vtx_ids;
  std::map<int, const Connection*> vtx_ids_connections;

  int id = 1;
  for (const auto& nac : substitution_point->neighbors_available_connections) {
    vtx_ids_neighbors.insert({ id, nac.first });
    neighbors_vtx_ids.insert({ nac.first, id });
    ++id;
  };
  for (const auto& c : substitute.connections) {
    for (const auto& cp : c.second) {
      vtx_ids_connections.insert({ id, &c.first });
      for (const auto& nac : substitution_point->neighbors_available_connections) {
        if (nac.second.IsCompatibleWith(c.first, query_results.compatibilities)) {
          edges.push_back(std::make_pair(neighbors_vtx_ids[nac.first], id));
        };
      };
      ++id;
    };
  };

  // Run the Hopcroft-Karp algorithm to find a solution to the Maximum
  // Bipartite Matching problem.
  BipartiteGraph graph(vtx_ids_neighbors.size(), vtx_ids_connections.size(), edges);
  PATH matching = graph.HopcroftKarp(prng);

  // Sanity check. Make sure that the cardinality of the matching
  // corresponds to the number of neighboring MolBricks. Otherwise
  // the substitute can't possibly be valid.
  assert(matching.size() == substitution_point->n_neighbors);

  // Reconnect the atoms, store the Connections to carry out the connection
  // map update and update the schematic.
  PATH::const_iterator m_it, m_begin_it = matching.begin(), m_end_it = matching.end();
  unsigned neighbor_schematic_idx, n_connections, connection_idx;
  unsigned original_brick_interactor_idx, original_neighbor_interactor_idx, new_brick_interactor_idx, new_neighbor_interactor_idx;
  Connection original_forward_connection, original_backward_connection, new_forward_connection, new_backward_connection;
  ConnectionsTable substituted_available_connections = substituted->connections;
  ConnectionsTable substitute_available_connections = substitute.connections;
  // Loop over the neighbors of the substituted MolBrick and:
  for (auto& nac : substitution_point->neighbors_available_connections) {
    MolBrick* neighbor = nac.first;
    neighbor_schematic_idx = neighbor->schematic_idx;

    // Retrieve the original Connection/ConnectionPoints involved.
    original_brick_interactor_idx = schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx];
    original_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx];
    original_forward_connection = Connection(schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx]);
    original_backward_connection = Connection(schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx]);

    // Remove the original Connections.
    internal_connections.RemoveConnection(original_forward_connection, original_brick_interactor_idx);
    internal_connections.MoveConnectionTo(original_backward_connection, original_neighbor_interactor_idx, connections);
    substituted_available_connections.RemoveConnection(original_forward_connection, original_brick_interactor_idx);

    // Record the deleted neighbor-to-substituted link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, original_neighbor_interactor_idx, substituted->pseudofragment_id);
    };

    // Select random ConnectionPoints for the Maximum Bipartite Matching solution.
    int neighbor_vtx_id = neighbors_vtx_ids[nac.first];
    m_it = std::find_if(m_begin_it, m_end_it,
      [=](const Edge* edge) {
        return edge->GetStart()->GetID() == neighbor_vtx_id;
      });
    int connection_vtx_id = (*m_it)->GetEnd()->GetID();
    new_forward_connection = *vtx_ids_connections[connection_vtx_id];
    new_backward_connection = nac.second.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

    const CONNECTION_POINT_VECTOR& brick_cpoint_vector = substitute_available_connections[new_forward_connection];
    n_connections = brick_cpoint_vector.size();
    connection_idx = 0;
    if (n_connections > 1) {
      connection_idx = uniform_int_distributions[n_connections](prng);
    };
    new_brick_interactor_idx = brick_cpoint_vector[connection_idx].atom_idx;

    const CONNECTION_POINT_VECTOR& neighbor_cpoint_vector = nac.second[new_backward_connection];
    n_connections = neighbor_cpoint_vector.size();
    connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
        connection_idx = cpoint_distribution(prng);
      };
    };
    new_neighbor_interactor_idx = neighbor_cpoint_vector[connection_idx].atom_idx;

    // Add the new Connections.
    ConnectionPoint* cpoint = internal_connections.AddConnection(new_forward_connection, new_brick_interactor_idx);
    if (guided) {
      query_results.AssignFreshWeights(new_forward_connection, cpoint);
    };
    connections.MoveConnectionTo(new_backward_connection, new_neighbor_interactor_idx, internal_connections);
    substitute_available_connections.RemoveConnection(new_forward_connection, new_brick_interactor_idx);

    // Record the new neighbor-to-substitute link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, new_neighbor_interactor_idx, substitute.pseudofragment_id);
    };

    // Update the ReconstructionsSchematic.
    schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.start_atom_type;
    schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.start_atom_type;
    schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.end_atom_type;
    schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.end_atom_type;
    schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.bond_type;
    schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.bond_type;
    schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx] = new_brick_interactor_idx;
    schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[brick_schematic_idx][neighbor_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[neighbor_schematic_idx][brick_schematic_idx] = new_brick_interactor_idx;

    // Create the new bonds.
    unsigned mutable_neighbor_interactor_idx = GetMutableAtomIdx(pseudomol, new_neighbor_interactor_idx);
    unsigned mutable_brick_interactor_idx = GetMutableAtomIdx(pseudomol, new_brick_interactor_idx);
    RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
    pseudomol.addBond(mutable_neighbor_interactor_idx, mutable_brick_interactor_idx, bond_type);
  };

  // Finish updating the ConnectionsTables by adjusting the peripheral connections.
  if (guided) {
    std::vector<ConnectionPoint*> cpoints;
    for (auto& c : substitute_available_connections) {
      cpoints = connections.AddConnections(c.first, c.second);
      query_results.AssignFreshWeights(c.first, cpoints);
    };
    assert(AllConnectionsHaveWeights());
  } else {
    for (auto& c : substitute_available_connections) {
      connections.AddConnections(c.first, c.second);
    };
  };
  for (const auto& c : substituted_available_connections) {
    connections.RemoveConnections(c.first, c.second);
  };

  // Update the rest of the ReconstructedMol's attributes.
  level += (substitute.level - substituted->level);
  // IMPROVE: Check that the ring number control still works
  if (increase_rings) {
    n_ring_atoms += substitute.level;
  } else if (decrease_rings) {
    n_ring_atoms -= substituted->level;
  } else {
    if (substitute.has_ring || substituted->has_ring) {
      assert(substitute.has_ring && substituted->has_ring);
      n_ring_atoms += (substitute.level - substituted->level);
    };
  };
  max_atom_idx = substitute.max_atom_idx;
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // Move operator to insert substitute in bricks?
  bricks[brick_schematic_idx] = substitute;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  // Signal function success.
  return true;
};

bool ReconstructedMol::Crossover(const std::string & location, ConnectionQueryResults & query_results, SizeController & controller, std::mt19937 & prng, ReconstructionsInventory & inventory, bool guided, bool update_stereo) {
  // In a crossover operation a MolBrick of the subject ReconstructedMol is
  // substituted by one of the MolBricks forming part of the current pool of
  // ReconstructedMols.

  // Note that this operation isn't a true genetic crossover since only the
  // acceptor ReconstructedMol is modified, while the donor ReconstructedMol
  // remains unchanged. The rationale behind this decision is as follows:
  //  (1) All ReconstructedMols evolve every generation. This is done to optimize
  //      the efficiency of each generation by minimizing the number of calls
  //      to external scoring functions. If the Crossover operation affected
  //      both parties it's likely that some ReconstructedMols would effectively
  //      undergo multiple mutations every generation, which would be confounding.
  //  (2) While the decision is score guided, all ReconstructedMols are allowed
  //      to Crossover with all others. As such, if a true Crossover were to
  //      take place between two ReconstructedMols with a large disparity in
  //      their scores the best ReconstructedMol would incur a high risk of
  //      becoming worse. In silico, since the "genetic material" can be
  //      duplicated at will, we can afford a "greedy" evolution approach
  //      where only one molecule benefits of the Crossover event.
  //  (3) Related to the above, if a Crossover mutation were to modify two
  //      ReconstructedMols, with one of them raising their score and the other
  //      one lowering it, it would become unclear whether the operation was
  //      beneficial and whether the changes should be kept or reverted.

  // Currently Crossovers aren't allowed when the ReconstructedMol consists
  // of a single MolBrick. While the substitution of the only MolBrick for
  // another random one would be possible it would imply completely disregarding
  // the starting base ReconstructedMol, which doesn't truely fit the definition
  // of an evolutive process.
  if (bricks.size() <= 1) {
    return false;
  };

  // Retrieve the ReconstructedMol's MolBricks of the type to be substituted
  // (i.e. peripheral or internal).
  std::vector<MolBrick*> bricks_to_evaluate;
  unsigned n_bricks;
  if (location == "peripheral") {
    bricks_to_evaluate = GetPeripheralBricks();
    n_bricks = bricks_to_evaluate.size();
  } else if (location == "internal") {
    bricks_to_evaluate = GetInternalBricks();
    n_bricks = bricks_to_evaluate.size();
    // If there are no internal MolBricks in the ReconstructedMol signal function
    // failure. I could default to a peripheral substitution in this case, but
    // it might cloud operation probabilities during development.
    if (n_bricks == 0) {
      return false;
    };
  } else {
    throw std::runtime_error(std::string("Invalid location option for Crossover."));
  };

  // Convert the subject MolBricks into SubstitutionPoints.
  std::vector<SubstitutionPoint> substitution_points;
  substitution_points.reserve(n_bricks);
  for (MolBrick* brick : bricks_to_evaluate) {
    substitution_points.push_back(SubstitutionPoint(brick, guided));
  };

  // Loop over the SubstititutionPoints and do the necessary preparations
  // for efficient Maximum Bipartite Matching solving.
  for (SubstitutionPoint& sp : substitution_points) {
    sp.PrepareHopcroftKarp(prng);
  };

  // Loop over the identified SubstitutionPoints and evaluate how
  // the number of rings can be modulated.
  // NOTE: If the ReconstructionsInventory has been kept up to date it's guaranteed to
  // contain atleast one valid substitute MolBrick: the substitutable MolBrick of the
  // SubstitutionPoint. Hence, it's not crucial to assert that this is the case.
  bool rings_can_be_added = false, rings_can_be_removed = false;
  for (SubstitutionPoint& sp : substitution_points) {
    sp.EvaluateSubstitutability(inventory, query_results.compatibilities, prng);
    if (sp.rings_can_be_added) {
      if (n_ring_atoms < controller.max_size) {
        rings_can_be_added = true;
      };
    };
    if (sp.rings_can_be_removed) {
      rings_can_be_removed = true;
    };
    if (rings_can_be_added && rings_can_be_removed) {
      break;
    };
  };

  // If both kinds of Pseudofragments can be used as substitutes,
  // decide which one to use. Otherwise, use the only type that can be used.
  bool constant_rings = false, increase_rings = false, decrease_rings = false;
  if (rings_can_be_added && rings_can_be_removed) {
    int decision = controller.Decide(n_ring_atoms, prng);
    if (decision == 0) {
      constant_rings = true;
    } else if (decision == 1) {
      increase_rings = true;
    } else if (decision == -1) {
      decrease_rings = true;
    };
  } else if (rings_can_be_added) {
    bool change = false;
    if (n_ring_atoms < controller.max_size) {
      change = controller.DecideChange(n_ring_atoms, prng);
    };
    if (change) {
      increase_rings = true;
    } else {
      constant_rings = true;
    };
  } else if (rings_can_be_removed) {
    bool change = controller.DecideChange(n_ring_atoms, prng);
    if (change) {
      decrease_rings = true;
    } else {
      constant_rings = true;
    };
  } else {
    constant_rings = true;
  };

  // Loop over the subset of SubstitutionPoints and attempt to retrieve
  // substitute Pseudofragments for them that modify the ReconstructedMol's
  // number of rings in the desired way. If substitutes are found, calculate
  // the weight of the SubstitutionPoint.
  unsigned n_substitution_points = substitution_points.size();
  std::vector<SubstitutionPoint*> substitition_point_candidates;
  std::vector<float> substitution_point_weights;
  substitition_point_candidates.reserve(n_substitution_points);
  substitution_point_weights.reserve(n_substitution_points);
  bool substitutes_exist = false;
  // If the number of rings should be kept constant
  if (constant_rings) {
    // Calculate the weights of both acyclic and ring MolBricks in the inventory.
    inventory.CalcWeights();
    for (SubstitutionPoint& sp : substitution_points) {
      // and the SubstitutionPoint pertains to a ring MolBrick
      if (sp.has_ring) {
        sp.RetrieveRingSubstitutes(inventory, query_results.compatibilities, prng);
        assert(sp.has_ring_substitutes);
        sp.CalcInventoryRingBasedWeight();
        substitition_point_candidates.push_back(&sp);
        substitution_point_weights.push_back(sp.weight);
      // and the SubstitutionPoint pertains to an acyclic MolBrick
      } else {
        // retrieve acyclic substitutes.
        sp.RetrieveAcyclicSubstitutes(inventory, query_results.compatibilities, prng);
        assert(sp.has_acyclic_substitutes);
        sp.CalcInventoryAcyclicBasedWeight();
        substitition_point_candidates.push_back(&sp);
        substitution_point_weights.push_back(sp.weight);
      };
    };
  // If the number of rings should be increased
  } else if (increase_rings) {
    // Calculate the weights of the ring MolBricks in the inventory.
    inventory.CalcRingWeights();
    for (SubstitutionPoint& sp : substitution_points) {
      // and the SubstitutionPoint pertains to an acyclic MolBrick
      if (!sp.has_ring) {
        // retrieve ring substitutes.
        sp.RetrieveRingSubstitutes(inventory, query_results.compatibilities, prng);
        if (sp.has_ring_substitutes) {
          substitutes_exist = true;
        };
        sp.CalcInventoryRingBasedWeight();
        substitition_point_candidates.push_back(&sp);
        substitution_point_weights.push_back(sp.weight);
      };
    };
    assert(substitutes_exist);
  // If the number of rings should be decreased
  } else if (decrease_rings) {
    // Calculate the weights of the acyclic MolBricks in the inventory.
    inventory.CalcAcyclicWeights();
    for (SubstitutionPoint& sp : substitution_points) {
      // and the SubstitutionPoint pertains to a ring MolBrick
      if (sp.has_ring) {
        // retrieve acyclic substitutes.
        sp.RetrieveAcyclicSubstitutes(inventory, query_results.compatibilities, prng);
        if (sp.has_acyclic_substitutes) {
          substitutes_exist = true;
        };
        sp.CalcInventoryAcyclicBasedWeight();
        substitition_point_candidates.push_back(&sp);
        substitution_point_weights.push_back(sp.weight);
      };
    };
    assert(substitutes_exist);
  };

  // Choose a SubstitutionPoint (i.e. MolBrick to be substituted)
  // through weighted selection.
  std::discrete_distribution<unsigned> sp_distribution(substitution_point_weights.begin(), substitution_point_weights.end());
  SubstitutionPoint* substitution_point = substitition_point_candidates[sp_distribution(prng)];
  MolBrick* substituted = substitution_point->substitutable;
  unsigned brick_schematic_idx = substituted->schematic_idx;

  // Choose a MolBrick compatible with the aforementioned SubstitutionPoint
  // through weighted selection.
  INVENTORY_RESULT* inventory_result;
  // If the number of rings should be kept constant
  if (constant_rings) {
    // and the SubstitutionPoint pertains to a ring MolBrick
    if (substitution_point->has_ring) {
      inventory_result = &substitution_point->inventory_ring_substitutes;
      // and the SubstitutionPoint pertains to an acyclic MolBrick
    } else {
      inventory_result = &substitution_point->inventory_acyclic_substitutes;
    };
  // If the number of rings should be increased
  } else if (increase_rings) {
    inventory_result = &substitution_point->inventory_ring_substitutes;
  // If the number of rings should be decreased
  } else {
    inventory_result = &substitution_point->inventory_acyclic_substitutes;
  };
  assert(!inventory_result->first.empty());
  unsigned idx = 0;
  if (inventory_result->first.size() > 1) {
    std::discrete_distribution<unsigned> brick_distribution(inventory_result->second.begin(), inventory_result->second.end());
    idx = brick_distribution(prng);
  };
  MolBrick* inventory_substitute = inventory_result->first[idx];

  // Create a copy of the substitute to update its attributes without affecting
  // the original substitute MolBrick (which is still part of the ReconstructionsInventory).
  MolBrick substitute = *inventory_substitute;
  substitute.UpdateImmutableIndices(max_atom_idx + 1);
  substitute.schematic_idx = brick_schematic_idx;
  substitute.owner = this;

  // Reset the GeneticLog.
  if (guided) {
    genetic_log.Clear();
    genetic_log.SetImportantPseudofragmentID(substituted->pseudofragment_id);
    genetic_log.SetImportantPseudofragmentID(substitute.pseudofragment_id);
  };

  // Update the pseudomol. This entails removing the atoms of the MolBrick
  // to be substituted, adding the atoms of the new substitute MolBrick and
  // reestablishing the bonds between them.

  // Remove the atoms of the substituted MolBrick in sorted descending order
  // to avoid atom reindexing.
  std::vector<unsigned> substituted_mutable_atom_indices = GetMutableAtomIndices(pseudomol, substituted->atom_indices);
  std::sort(substituted_mutable_atom_indices.begin(), substituted_mutable_atom_indices.end(), std::greater<unsigned>());
  for (unsigned atom_idx : substituted_mutable_atom_indices) {
    pseudomol.removeAtom(atom_idx);
  };

  // Add the atoms of the substitute MolBrick.
  pseudomol.insertMol(substitute.pseudomol);

  // Determine how the substitute MolBrick should connect to the neighboring
  // MolBricks by finding a solution to the Maximum Bipartite Matching problem.
  PATH matching = substitution_point->HopcroftKarp(substitute, query_results.compatibilities, prng);

  // Sanity check. Make sure that the cardinality of the matching
  // corresponds to the number of neighboring MolBricks. Otherwise
  // the substitute can't possibly be valid.
  assert(matching.size() == substitution_point->n_neighbors);

  // Reconnect the atoms, store the Connections to carry out the connection
  // map update and update the schematic.
  PATH::const_iterator m_it, m_begin_it = matching.begin(), m_end_it = matching.end();
  unsigned neighbor_schematic_idx, n_connections, connection_idx;
  unsigned original_brick_interactor_idx, original_neighbor_interactor_idx, new_brick_interactor_idx, new_neighbor_interactor_idx;
  Connection original_forward_connection, original_backward_connection, new_forward_connection, new_backward_connection;
  ConnectionsTable substituted_available_connections = substituted->connections;
  ConnectionsTable substitute_available_connections = substitute.connections;
  // Loop over the neighbors of the substituted MolBrick and:
  for (auto& nac : substitution_point->neighbors_available_connections) {
    MolBrick* neighbor = nac.first;
    neighbor_schematic_idx = neighbor->schematic_idx;

    // Retrieve the original Connection/ConnectionPoints involved.
    original_brick_interactor_idx = schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx];
    original_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx];
    original_forward_connection = Connection(schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx]);
    original_backward_connection = Connection(schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx]);

    // Remove the original Connections.
    internal_connections.RemoveConnection(original_forward_connection, original_brick_interactor_idx);
    internal_connections.MoveConnectionTo(original_backward_connection, original_neighbor_interactor_idx, connections);
    substituted_available_connections.RemoveConnection(original_forward_connection, original_brick_interactor_idx);

    // Record the deleted neighbor-to-substituted link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, original_neighbor_interactor_idx, substituted->pseudofragment_id);
    };

    // Select random ConnectionPoints for the Maximum Bipartite Matching solution.
    int neighbor_vtx_id = substitution_point->neighbors_vtx_ids[nac.first];
    m_it = std::find_if(m_begin_it, m_end_it,
      [=](const Edge* edge) {
        return edge->GetStart()->GetID() == neighbor_vtx_id;
      });
    int connection_vtx_id = (*m_it)->GetEnd()->GetID();
    new_forward_connection = *(substitution_point->vtx_ids_connections[connection_vtx_id]);
    new_backward_connection = nac.second.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

    const CONNECTION_POINT_VECTOR& brick_cpoint_vector = substitute_available_connections[new_forward_connection];
    n_connections = brick_cpoint_vector.size();
    connection_idx = 0;
    if (n_connections > 1) {
      connection_idx = uniform_int_distributions[n_connections](prng);
    };
    new_brick_interactor_idx = brick_cpoint_vector[connection_idx].atom_idx;

    const CONNECTION_POINT_VECTOR& neighbor_cpoint_vector = nac.second[new_backward_connection];
    n_connections = neighbor_cpoint_vector.size();
    connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> cpoint_distribution(0, n_connections - 1);
        connection_idx = cpoint_distribution(prng);
      };
    };
    new_neighbor_interactor_idx = neighbor_cpoint_vector[connection_idx].atom_idx;

    // Add the new Connections.
    ConnectionPoint* cpoint = internal_connections.AddConnection(new_forward_connection, new_brick_interactor_idx);
    if (guided) {
      query_results.AssignFreshWeights(new_forward_connection, cpoint);
    };
    connections.MoveConnectionTo(new_backward_connection, new_neighbor_interactor_idx, internal_connections);
    substitute_available_connections.RemoveConnection(new_forward_connection, new_brick_interactor_idx);

    // Record the new neighbor-to-substitute link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, new_neighbor_interactor_idx, substitute.pseudofragment_id);
    };

    // Update the ReconstructionsSchematic.
    schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.start_atom_type;
    schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.start_atom_type;
    schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.end_atom_type;
    schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.end_atom_type;
    schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.bond_type;
    schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.bond_type;
    schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx] = new_brick_interactor_idx;
    schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[brick_schematic_idx][neighbor_schematic_idx] = new_neighbor_interactor_idx;
    schematic.end_atom_idx[neighbor_schematic_idx][brick_schematic_idx] = new_brick_interactor_idx;

    // Create the new bonds.
    unsigned mutable_neighbor_interactor_idx = GetMutableAtomIdx(pseudomol, new_neighbor_interactor_idx);
    unsigned mutable_brick_interactor_idx = GetMutableAtomIdx(pseudomol, new_brick_interactor_idx);
    RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
    pseudomol.addBond(mutable_neighbor_interactor_idx, mutable_brick_interactor_idx, bond_type);
  };

  // Finish updating the ConnectionsTables by adjusting the peripheral connections.
  if (guided) {
    std::vector<ConnectionPoint*> cpoints;
    for (auto& c : substitute_available_connections) {
      cpoints = connections.AddConnections(c.first, c.second);
      query_results.AssignFreshWeights(c.first, cpoints);
    };
    assert(AllConnectionsHaveWeights());
  } else {
    for (auto& c : substitute_available_connections) {
      connections.AddConnections(c.first, c.second);
    };
  };
  for (const auto& c : substituted_available_connections) {
    connections.RemoveConnections(c.first, c.second);
  };

  // Update the rest of the ReconstructedMol's attributes.
  level += (substitute.level - substituted->level);
  if (increase_rings) {
    n_ring_atoms += substitute.level;
  } else if (decrease_rings) {
    n_ring_atoms -= substituted->level;
  } else {
    if (substitute.has_ring || substituted->has_ring) {
      assert(substitute.has_ring && substituted->has_ring);
      n_ring_atoms += (substitute.level - substituted->level);
    };
  };
  max_atom_idx = substitute.max_atom_idx;
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // IMPROVE: Move operator to insert substitute in bricks?
  bricks[brick_schematic_idx] = substitute;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  // Signal function success.
  return true;
};

bool ReconstructedMol::Translation(ConnectionQueryResults& query_results, std::mt19937& prng, bool guided, bool update_stereo) {
  // In a Translation operation one of the ReconstructedMol's MolBricks is
  // moved from its original position to a different one. In practice this
  // means performing a deletion and an expansion in tandem. Note that:
  // (1) A Translation operation can also be a rotation in place.
  // (2) The translated MolBrick may be inserted peripherally or internally,
  //     regardless of its source position.

  // If the ReconstructedMol consists of a single MolBrick the latter can't be
  // moved. Atleast two MolBricks must be present for a potential Rotation.
  if (bricks.size() == 1) {
    return false;
  };

  // ###################### Translation MolBrick choice #######################
  // Select a random deletable MolBrick.
  std::vector<MolBrick*> candidates;
  // Loop over the ReconstructedMol's MolBricks.
  for (auto& b : bricks) {
    MolBrick* brick = &b.second;
    // If the MolBrick is peripheral it's always deletable.
    if (brick->IsPeripheral()) {
      candidates.push_back(brick);
      // If the MolBrick is internal it has to be converted into a DeletionPoint
      // and its deletability evaluated exhaustively.
    } else {
      DeletionPoint deletion_point(brick);
      deletion_point.Evaluate(query_results.compatibilities, prng);
      if (deletion_point.is_deletable) {
        candidates.push_back(brick);
      };
    };
  };
  // Any molecule with more than one MolBrick is guaranteed to have atleast two deletable MolBricks.
  unsigned n_candidates = candidates.size();
  assert(n_candidates >= 2);
  unsigned idx;
  if (n_candidates <= max_size_uniform_int_distributions) {
    idx = uniform_int_distributions[n_candidates](prng);
  } else {
    std::uniform_int_distribution<unsigned> distribution(0, n_candidates - 1);
    idx = distribution(prng);
  };
  MolBrick* brick = candidates[idx];
  unsigned brick_schematic_idx = brick->schematic_idx;

  // Reset the GeneticLog.
  genetic_log.Clear();
  genetic_log.SetImportantPseudofragmentID(brick->pseudofragment_id);

  // ########################### MolBrick deletion ############################
  // Determine whether the chosen MolBrick is peripheral or internal and carry
  // out the relevant steps to delete it.
  // If the MolBrick is peripheral:
  if (brick->IsPeripheral()) {
    // Determine which MolBrick is the neighbor to the one to be deleted. Since
    // the MolBrick is peripheral it can only have a single neighbor.
    unsigned neighbor_schematic_idx = 0;
    for (unsigned adjacency_bool : schematic.adjacency[brick->schematic_idx]) {
      if (adjacency_bool == 1) {
        break;
      };
      ++neighbor_schematic_idx;
    };

    // Update the ReconstructedMol's connection tables. This entails:
    //  (1) Moving the Connections involved in the deleted MolBrick-neighbor
    //      interaction from the internal to the peripheral ConnectionTable.
    //  (2) Removing the MolBrick's Connections from the peripheral ConnectionTable.
    Connection forward_connection(schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx]);
    Connection backward_connection(schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx]);
    unsigned immutable_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx];
    unsigned immutable_brick_interactor_idx = schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx];

    internal_connections.MoveConnectionTo(forward_connection, immutable_neighbor_interactor_idx, connections);
    internal_connections.MoveConnectionTo(backward_connection, immutable_brick_interactor_idx, connections);
    for (const auto& c : brick->connections) {
      connections.RemoveConnections(c.first, c.second);
    };

    // Record the deleted neighbor-to-deleted link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, forward_connection, immutable_neighbor_interactor_idx, brick->pseudofragment_id);
    };

    // Remove the Bond between the subject MolBrick and its neighbor MolBrick.
    pseudomol.removeBond(GetMutableAtomIdx(pseudomol, immutable_brick_interactor_idx), GetMutableAtomIdx(pseudomol, immutable_neighbor_interactor_idx));

    // Update the ReconstructionSchematic. In a deletion this means clearing
    // the values of the ReconstructionSchematic's row and column corresponding
    // to the MolBrick in question.
    schematic.ClearRow(brick_schematic_idx);
    schematic.ClearColumn(brick_schematic_idx);

  // If the MolBrick is internal:
  } else {
    // Convert it into a DeletionPoint.
    DeletionPoint deletion_point(brick);

    // Generate the neighbor matchings for the chosen DeletionPoint and randomly
    // choose one of them.
    deletion_point.GenerateMatchings(query_results.compatibilities, prng);
    unsigned matching_idx = 0, n_matchings = deletion_point.other_neighbors_matchings.size();
    assert(n_matchings > 0);
    if (n_matchings > 1) {
      if (n_matchings <= max_size_uniform_int_distributions) {
        matching_idx = uniform_int_distributions[n_matchings](prng);
      } else {
        std::uniform_int_distribution<unsigned> matchings_distribution(0, n_matchings - 1);
        matching_idx = matchings_distribution(prng);
      };
    };
    std::map<const MolBrick*, std::map<const MolBrick*, const Connection*>>::iterator it = deletion_point.other_neighbors_matchings.begin();
    std::advance(it, matching_idx);

    // Retrieve the MolBrick that will become the "core" of the new Connections
    // (i.e. the one that connects to all the remaining neighboring MolBrick).
    const MolBrick* core = it->first;
    unsigned core_schematic_idx = core->schematic_idx;
    assert(deletion_point.neighbors_available_connections[core].n_connection_points >= deletion_point.n_other_neighbors);

    // Update the ReconstructedMol's Connections and ReconstructionSchematic.
    Connection new_forward_connection, new_backward_connection;
    ConnectionsTable core_available_connections = deletion_point.neighbors_available_connections[core];

    // First remove the Connection and Bond between the deleted MolBrick and the
    // new core MolBrick.
    unsigned brick_interactor_idx = schematic.start_atom_idx[brick_schematic_idx][core_schematic_idx];
    unsigned core_interactor_idx = schematic.start_atom_idx[core_schematic_idx][brick_schematic_idx];
    Connection original_forward_connection(schematic.start_atom_type[brick_schematic_idx][core_schematic_idx], schematic.end_atom_type[brick_schematic_idx][core_schematic_idx], schematic.bond_type[brick_schematic_idx][core_schematic_idx]);
    Connection original_backward_connection(schematic.start_atom_type[core_schematic_idx][brick_schematic_idx], schematic.end_atom_type[core_schematic_idx][brick_schematic_idx], schematic.bond_type[core_schematic_idx][brick_schematic_idx]);
    internal_connections.RemoveConnection(original_forward_connection, brick_interactor_idx);
    internal_connections.MoveConnectionTo(original_backward_connection, core_interactor_idx, connections);
    pseudomol.removeBond(GetMutableAtomIdx(pseudomol, brick_interactor_idx), GetMutableAtomIdx(pseudomol, core_interactor_idx));

    // Record the deleted core-to-deleted link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, core_interactor_idx, brick->pseudofragment_id);
    };

    // Thereafter do the same for all the other neighbors of the deleted MolBrick.
    for (const auto& match : it->second) {
      const MolBrick* neighbor = match.first;
      unsigned neighbor_schematic_idx = neighbor->schematic_idx;
      ConnectionsTable& neighbor_available_connections = deletion_point.neighbors_available_connections[neighbor];

      // Retrieve the original Connection/ConnectionPoints involved.
      brick_interactor_idx = schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx];
      unsigned original_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx];
      original_forward_connection = Connection(schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx], schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx]);
      original_backward_connection = Connection(schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx], schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx]);

      // Remove the original Connections and Bond.
      internal_connections.RemoveConnection(original_forward_connection, brick_interactor_idx);
      internal_connections.MoveConnectionTo(original_backward_connection, original_neighbor_interactor_idx, connections);
      pseudomol.removeBond(GetMutableAtomIdx(pseudomol, brick_interactor_idx), GetMutableAtomIdx(pseudomol, original_neighbor_interactor_idx));

      // Record the deleted neighbor-to-deleted link in the GeneticLog.
      if (guided) {
        genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, original_neighbor_interactor_idx, brick->pseudofragment_id);
      };

      // Select random ConnectionPoints for the Maximum Bipartite Matching solution.
      new_forward_connection = *match.second;
      new_backward_connection = neighbor_available_connections.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

      const CONNECTION_POINT_VECTOR& core_cpoint_vector = core_available_connections[new_forward_connection];
      unsigned n_connections = core_cpoint_vector.size();
      assert(n_connections > 0);
      unsigned cpoint_idx = 0;
      if (n_connections > 1) {
        if (n_connections <= max_size_uniform_int_distributions) {
          cpoint_idx = uniform_int_distributions[n_connections](prng);
        } else {
          std::uniform_int_distribution<unsigned> cpoints_distribution(0, n_connections - 1);
          cpoint_idx = cpoints_distribution(prng);
        };
      };
      core_interactor_idx = core_cpoint_vector[cpoint_idx].atom_idx;

      const CONNECTION_POINT_VECTOR& neighbor_cpoint_vector = neighbor_available_connections[new_backward_connection];
      n_connections = neighbor_cpoint_vector.size();
      assert(n_connections > 0);
      cpoint_idx = 0;
      if (n_connections > 1) {
        if (n_connections <= max_size_uniform_int_distributions) {
          cpoint_idx = uniform_int_distributions[n_connections](prng);
        } else {
          std::uniform_int_distribution<unsigned> cpoints_distribution(0, n_connections - 1);
          cpoint_idx = cpoints_distribution(prng);
        };
      };
      unsigned new_neighbor_interactor_idx = neighbor_cpoint_vector[cpoint_idx].atom_idx;

      // Add the new Connections.
      connections.MoveConnectionTo(new_forward_connection, core_interactor_idx, internal_connections);
      connections.MoveConnectionTo(new_backward_connection, new_neighbor_interactor_idx, internal_connections);
      core_available_connections.RemoveConnection(new_forward_connection, core_interactor_idx);

      // Record the new core-to-neighbor and neighbor-to-core links in the GeneticLog.
      if (guided) {
        genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_forward_connection, core_interactor_idx, neighbor->pseudofragment_id);
        genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, new_neighbor_interactor_idx, core->pseudofragment_id);
      };

      // Update the ReconstructionsSchematic.
      schematic.adjacency[core_schematic_idx][neighbor_schematic_idx] = 1;
      schematic.adjacency[neighbor_schematic_idx][core_schematic_idx] = 1;
      schematic.start_atom_type[core_schematic_idx][neighbor_schematic_idx] = new_forward_connection.start_atom_type;
      schematic.start_atom_type[neighbor_schematic_idx][core_schematic_idx] = new_backward_connection.start_atom_type;
      schematic.end_atom_type[core_schematic_idx][neighbor_schematic_idx] = new_forward_connection.end_atom_type;
      schematic.end_atom_type[neighbor_schematic_idx][core_schematic_idx] = new_backward_connection.end_atom_type;
      schematic.bond_type[core_schematic_idx][neighbor_schematic_idx] = new_forward_connection.bond_type;
      schematic.bond_type[neighbor_schematic_idx][core_schematic_idx] = new_backward_connection.bond_type;
      schematic.start_atom_idx[core_schematic_idx][neighbor_schematic_idx] = core_interactor_idx;
      schematic.start_atom_idx[neighbor_schematic_idx][core_schematic_idx] = new_neighbor_interactor_idx;
      schematic.end_atom_idx[core_schematic_idx][neighbor_schematic_idx] = new_neighbor_interactor_idx;
      schematic.end_atom_idx[neighbor_schematic_idx][core_schematic_idx] = core_interactor_idx;

      // Create the new bonds.
      unsigned mutable_core_interactor_idx = GetMutableAtomIdx(pseudomol, core_interactor_idx);
      unsigned mutable_neighbor_interactor_idx = GetMutableAtomIdx(pseudomol, new_neighbor_interactor_idx);
      RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
      pseudomol.addBond(mutable_core_interactor_idx, mutable_neighbor_interactor_idx, bond_type);
    };

    // Finish updating the ConnectionsTables by adjusting the peripheral connections.
    for (const auto& c : brick->GetAvailableConnections(prng)) {
      connections.RemoveConnections(c.first, c.second);
    };

    // Update the rest of the ReconstructionSchematic by clearing the data
    // pertaining to the deleted MolBrick.
    schematic.ClearRow(brick_schematic_idx);
    schematic.ClearColumn(brick_schematic_idx);
  };

  // ######################## Reinsertion point choice #########################
  // Enumerate all points at which the deleted MolBrick could be added again.
  // This includes both peripheral ConnectionPoints and internal InsertionPoints.
  std::vector<const Connection*> peripheral_points;
  peripheral_points.reserve(connections.n_connection_points);
  for (const auto& c : connections) {
    if (brick->connections.IsCompatibleWith(c.first, query_results.compatibilities)) {
      for (const ConnectionPoint& cpoint : c.second) {
        peripheral_points.push_back(&c.first);
      };
    };
  };

  std::vector<InsertionPoint> internal_points;
  internal_points.reserve(bricks.size());
  for (auto& b : bricks) {
    if (&b.second != brick) {
      InsertionPoint insertion_point(&b.second, prng);
      insertion_point.GenerateMatchingsForBrick(brick, query_results.compatibilities, prng);
      if (insertion_point.n_matchings > 0) {
        internal_points.push_back(std::move(insertion_point));
      };
    };
  };

  // Randomly choose one of the points to which it could be moved.
  unsigned n_peripheral_points = peripheral_points.size();
  unsigned n_internal_points = internal_points.size();
  unsigned n_points = n_peripheral_points + n_internal_points;
  assert(n_points > 0);
  unsigned point_idx = 0;
  if (n_points > 1) {
    if (n_points <= max_size_uniform_int_distributions) {
      point_idx = uniform_int_distributions[n_points](prng);
    } else {
      std::uniform_int_distribution<unsigned> distribution(0, n_points - 1);
      point_idx = distribution(prng);
    };
  };

  // ########################## MolBrick reinsertion ###########################
  // If the chosen point is a peripheral ConnectionPoint perform the operations
  // analogous to a PeripheralExpansion.
  if ((n_peripheral_points > 0) && (point_idx < n_peripheral_points)) {
    // For the chosen type of Connection, randomly choose two compatible
    // interaction atoms (one on the ReconstructedMol and one on the MolBrick).
    Connection chosen_connection = *peripheral_points[point_idx];
    RDKit::Bond::BondType bond_type = int_bond_type_table[chosen_connection.bond_type];
    const CONNECTION_POINT_VECTOR& reconstruction_cpoint_vector = connections[chosen_connection];
    unsigned n_connections = reconstruction_cpoint_vector.size();
    unsigned connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> distribution(0, n_connections - 1);
        connection_idx = distribution(prng);
      };
    };
    unsigned immutable_reconstruction_interactor_idx = reconstruction_cpoint_vector[connection_idx].atom_idx;

    Connection compatible_connection = brick->connections.GetRandomCompatibleConnection(chosen_connection, query_results.compatibilities, prng);
    const CONNECTION_POINT_VECTOR& brick_cpoint_vector = brick->connections[compatible_connection];
    n_connections = brick_cpoint_vector.size();
    connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> distribution(0, n_connections - 1);
        connection_idx = distribution(prng);
      };
    };
    unsigned immutable_brick_interactor_idx = brick_cpoint_vector[connection_idx].atom_idx;

    // Update the pseudomol by creating a bond between the specified interaction atoms.
    unsigned mutable_reconstruction_interactor_idx = GetMutableAtomIdx(pseudomol, immutable_reconstruction_interactor_idx);
    unsigned mutable_brick_interactor_idx = GetMutableAtomIdx(pseudomol, immutable_brick_interactor_idx);
    pseudomol.addBond(mutable_reconstruction_interactor_idx, mutable_brick_interactor_idx, bond_type);

    // Get the schematic index of the interacting brick.
    const MolBrick* interactor = GetAtomSourceBrick(immutable_reconstruction_interactor_idx);
    unsigned interactor_schematic_idx = interactor->schematic_idx;

    // Update the peripheral and internal connection tables. This entails:
    //  (1) Adding all the Connections from the new brick to the peripheral
    //      ConnectionTable.
    //  (2) Moving the two Connections involved in connecting the original
    //      ReconstructedMol to the new brick from the peripheral to the internal
    //      ConnectionTable.
    if (guided) {
      std::vector<ConnectionPoint*> cpoints;
      for (auto& c : brick->connections) {
        cpoints = connections.AddConnections(c.first, c.second);
        query_results.AssignFreshWeights(c.first, cpoints);
      };
    } else {
      for (auto& c : brick->connections) {
        connections.AddConnections(c.first, c.second);
      };
    };
    connections.MoveConnectionTo(chosen_connection, immutable_reconstruction_interactor_idx, internal_connections);
    connections.MoveConnectionTo(compatible_connection, immutable_brick_interactor_idx, internal_connections);

    // Record the new neighbor-to-insert link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, chosen_connection, immutable_reconstruction_interactor_idx, brick->pseudofragment_id);
    };

    // Update the ReconstructionSchematic.
    schematic.adjacency[interactor_schematic_idx][brick_schematic_idx] = 1;
    schematic.adjacency[brick_schematic_idx][interactor_schematic_idx] = 1;
    schematic.start_atom_type[interactor_schematic_idx][brick_schematic_idx] = chosen_connection.start_atom_type;
    schematic.start_atom_type[brick_schematic_idx][interactor_schematic_idx] = compatible_connection.start_atom_type;
    schematic.end_atom_type[interactor_schematic_idx][brick_schematic_idx] = chosen_connection.end_atom_type;
    schematic.end_atom_type[brick_schematic_idx][interactor_schematic_idx] = compatible_connection.end_atom_type;
    schematic.bond_type[interactor_schematic_idx][brick_schematic_idx] = chosen_connection.bond_type;
    schematic.bond_type[brick_schematic_idx][interactor_schematic_idx] = compatible_connection.bond_type;
    schematic.start_atom_idx[interactor_schematic_idx][brick_schematic_idx] = immutable_reconstruction_interactor_idx;
    schematic.start_atom_idx[brick_schematic_idx][interactor_schematic_idx] = immutable_brick_interactor_idx;
    schematic.end_atom_idx[interactor_schematic_idx][brick_schematic_idx] = immutable_brick_interactor_idx;
    schematic.end_atom_idx[brick_schematic_idx][interactor_schematic_idx] = immutable_reconstruction_interactor_idx;

  // If the chosen point is an internal InsertionPoint:
  } else {
    // Retrieve the InsertionPoint and its owner MolBrick.
    assert(!internal_points.empty());
    point_idx -= n_peripheral_points;
    InsertionPoint& insertion_point = internal_points[point_idx];
    MolBrick* owner = insertion_point.owner;
    unsigned owner_schematic_idx = owner->schematic_idx;
    assert(insertion_point.n_matchings > 0);

    // Choose a random matching from the InsertionPoint.
    std::map<MOLBRICKS_COMBINATION, std::map<const MolBrick*, const Connection*>>::iterator it = insertion_point.matchings.begin();
    unsigned ip_idx = 0;
    if (insertion_point.n_matchings > 1) {
      if (insertion_point.n_matchings <= max_size_uniform_int_distributions) {
        ip_idx = uniform_int_distributions[insertion_point.n_matchings](prng);
      } else {
        std::uniform_int_distribution<unsigned> distribution(0, insertion_point.n_matchings - 1);
        ip_idx = distribution(prng);
      };
    };
    std::advance(it, ip_idx);
    const MOLBRICKS_COMBINATION& neighbors_combination = it->first;
    std::map<const MolBrick*, const Connection*>& matching = it->second;
    ConnectionsTable owner_available_connections = owner->GetAvailableConnections(prng, neighbors_combination);

    // Connect the MolBrick insert to the InsertionPoint's owner
    // MolBrick's neighbors according to the aforementioned matching.
    Connection original_forward_connection, original_backward_connection, new_forward_connection, new_backward_connection;
    ConnectionsTable brick_available_connections = brick->connections;
    for (MolBrick* neighbor : neighbors_combination) {
      unsigned neighbor_schematic_idx = neighbor->schematic_idx;
      ConnectionsTable& neighbor_available_connections = insertion_point.neighbors_available_connections[neighbor];

      // Remove the original Connection between the owner and the neighbor.
      // and break the bond between the owner and the neighbor.
      original_forward_connection = Connection(schematic.start_atom_type[owner_schematic_idx][neighbor_schematic_idx], schematic.end_atom_type[owner_schematic_idx][neighbor_schematic_idx], schematic.bond_type[owner_schematic_idx][neighbor_schematic_idx]);
      original_backward_connection = Connection(schematic.start_atom_type[neighbor_schematic_idx][owner_schematic_idx], schematic.end_atom_type[neighbor_schematic_idx][owner_schematic_idx], schematic.bond_type[neighbor_schematic_idx][owner_schematic_idx]);
      unsigned owner_interactor_idx = schematic.start_atom_idx[owner_schematic_idx][neighbor_schematic_idx];
      unsigned original_neighbor_interactor_idx = schematic.start_atom_idx[neighbor_schematic_idx][owner_schematic_idx];

      internal_connections.MoveConnectionTo(original_forward_connection, owner_interactor_idx, connections);
      internal_connections.MoveConnectionTo(original_backward_connection, original_neighbor_interactor_idx, connections);

      pseudomol.removeBond(GetMutableAtomIdx(pseudomol, owner_interactor_idx), GetMutableAtomIdx(pseudomol, original_neighbor_interactor_idx));

      // Record the deleted owner-to-neighbor and neighbor-to-owner in the GeneticLog.
      if (guided) {
        genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_forward_connection, owner_interactor_idx, neighbor->pseudofragment_id);
        genetic_log.AddEntry(GeneticLogEntry::Type::DELETION, original_backward_connection, original_neighbor_interactor_idx, owner->pseudofragment_id);
      };

      // Add the Connections between the MolBrick insert and the neighboring MolBrick.
      new_forward_connection = *matching[neighbor];
      new_backward_connection = neighbor_available_connections.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

      const CONNECTION_POINT_VECTOR& brick_cpoint_vector = brick_available_connections[new_forward_connection];
      unsigned n_connections = brick_cpoint_vector.size();
      assert(n_connections > 0);
      unsigned cpoint_idx = 0;
      if (n_connections > 1) {
        if (n_connections <= max_size_uniform_int_distributions) {
          cpoint_idx = uniform_int_distributions[n_connections](prng);
        } else {
          std::uniform_int_distribution<unsigned> distribution(0, n_connections - 1);
          cpoint_idx = distribution(prng);
        };
      };
      unsigned brick_interactor_idx = brick_cpoint_vector[cpoint_idx].atom_idx;

      const CONNECTION_POINT_VECTOR& neighbor_cpoint_vector = neighbor_available_connections[new_backward_connection];
      n_connections = neighbor_cpoint_vector.size();
      assert(n_connections > 0);
      cpoint_idx = 0;
      if (n_connections > 1) {
        if (n_connections <= max_size_uniform_int_distributions) {
          cpoint_idx = uniform_int_distributions[n_connections](prng);
        } else {
          std::uniform_int_distribution<unsigned> distribution(0, n_connections - 1);
          cpoint_idx = distribution(prng);
        };
      };
      unsigned new_neighbor_interactor_idx = neighbor_cpoint_vector[cpoint_idx].atom_idx;

      ConnectionPoint* cpoint = internal_connections.AddConnection(new_forward_connection, brick_interactor_idx);
      if (guided) {
        query_results.AssignFreshWeights(new_forward_connection, cpoint);
      };
      brick_available_connections.RemoveConnection(new_forward_connection, brick_interactor_idx);
      connections.MoveConnectionTo(new_backward_connection, new_neighbor_interactor_idx, internal_connections);

      // Record the new neighbor-to-insert link in the GeneticLog.
      if (guided) {
        genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, new_neighbor_interactor_idx, brick->pseudofragment_id);
      };

      // Clear the ReconstructionsSchematic's entries corresponding to
      // the InsertionPoint's owner MolBrick and the neighboring MolBrick.
      schematic.ClearCell(owner_schematic_idx, neighbor_schematic_idx);
      schematic.ClearCell(neighbor_schematic_idx, owner_schematic_idx);

      // Update the ReconstructionsSchematic's entries corresponding to
      // the the insert-neighbor Connection.
      schematic.adjacency[brick_schematic_idx][neighbor_schematic_idx] = 1;
      schematic.adjacency[neighbor_schematic_idx][brick_schematic_idx] = 1;
      schematic.start_atom_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.start_atom_type;
      schematic.start_atom_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.start_atom_type;
      schematic.end_atom_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.end_atom_type;
      schematic.end_atom_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.end_atom_type;
      schematic.bond_type[brick_schematic_idx][neighbor_schematic_idx] = new_forward_connection.bond_type;
      schematic.bond_type[neighbor_schematic_idx][brick_schematic_idx] = new_backward_connection.bond_type;
      schematic.start_atom_idx[brick_schematic_idx][neighbor_schematic_idx] = brick_interactor_idx;
      schematic.start_atom_idx[neighbor_schematic_idx][brick_schematic_idx] = new_neighbor_interactor_idx;
      schematic.end_atom_idx[brick_schematic_idx][neighbor_schematic_idx] = new_neighbor_interactor_idx;
      schematic.end_atom_idx[neighbor_schematic_idx][brick_schematic_idx] = brick_interactor_idx;

      // Create the Bond between the MolBrick insert and the neighboring MolBrick.
      RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
      pseudomol.addBond(GetMutableAtomIdx(pseudomol, brick_interactor_idx), GetMutableAtomIdx(pseudomol, new_neighbor_interactor_idx), bond_type);
    };

    // For the InsertionPoint's owner MolBrick only adding Connections and
    // updating the ReconstructionsSchematic is necessary. Updating the owner
    // last is important to free up Connections first.
    new_forward_connection = *matching[owner];
    new_backward_connection = owner_available_connections.GetRandomCompatibleConnection(new_forward_connection, query_results.compatibilities, prng);

    const CONNECTION_POINT_VECTOR& brick_cpoint_vector = brick_available_connections[new_forward_connection];
    unsigned n_connections = brick_cpoint_vector.size();
    assert(n_connections > 0);
    unsigned connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> distribution(0, n_connections - 1);
        connection_idx = distribution(prng);
      };
    };
    unsigned brick_interactor_idx = brick_cpoint_vector[connection_idx].atom_idx;

    const CONNECTION_POINT_VECTOR& owner_cpoint_vector = owner_available_connections.at(new_backward_connection);
    n_connections = owner_cpoint_vector.size();
    assert(n_connections > 0);
    connection_idx = 0;
    if (n_connections > 1) {
      if (n_connections <= max_size_uniform_int_distributions) {
        connection_idx = uniform_int_distributions[n_connections](prng);
      } else {
        std::uniform_int_distribution<unsigned> distribution(0, n_connections - 1);
        connection_idx = distribution(prng);
      };
    };
    unsigned owner_interactor_idx = owner_cpoint_vector[connection_idx].atom_idx;

    ConnectionPoint* cpoint = internal_connections.AddConnection(new_forward_connection, brick_interactor_idx);
    if (guided) {
      query_results.AssignFreshWeights(new_forward_connection, cpoint);
    };
    brick_available_connections.RemoveConnection(new_forward_connection, brick_interactor_idx);
    connections.MoveConnectionTo(new_backward_connection, owner_interactor_idx, internal_connections);

    // Record the new owner-to-insert link in the GeneticLog.
    if (guided) {
      genetic_log.AddEntry(GeneticLogEntry::Type::CREATION, new_backward_connection, owner_interactor_idx, brick->pseudofragment_id);
    };

    // Update the ReconstructionsSchematic's entries corresponding to
    // the the insert-neighbor Connection.
    schematic.adjacency[brick_schematic_idx][owner_schematic_idx] = 1;
    schematic.adjacency[owner_schematic_idx][brick_schematic_idx] = 1;
    schematic.start_atom_type[brick_schematic_idx][owner_schematic_idx] = new_forward_connection.start_atom_type;
    schematic.start_atom_type[owner_schematic_idx][brick_schematic_idx] = new_backward_connection.start_atom_type;
    schematic.end_atom_type[brick_schematic_idx][owner_schematic_idx] = new_forward_connection.end_atom_type;
    schematic.end_atom_type[owner_schematic_idx][brick_schematic_idx] = new_backward_connection.end_atom_type;
    schematic.bond_type[brick_schematic_idx][owner_schematic_idx] = new_forward_connection.bond_type;
    schematic.bond_type[owner_schematic_idx][brick_schematic_idx] = new_backward_connection.bond_type;
    schematic.start_atom_idx[brick_schematic_idx][owner_schematic_idx] = brick_interactor_idx;
    schematic.start_atom_idx[owner_schematic_idx][brick_schematic_idx] = owner_interactor_idx;
    schematic.end_atom_idx[brick_schematic_idx][owner_schematic_idx] = owner_interactor_idx;
    schematic.end_atom_idx[owner_schematic_idx][brick_schematic_idx] = brick_interactor_idx;

    // Create the Bond between the MolBrick insert and the owner MolBrick.
    RDKit::Bond::BondType bond_type = int_bond_type_table[new_forward_connection.bond_type];
    pseudomol.addBond(GetMutableAtomIdx(pseudomol, brick_interactor_idx), GetMutableAtomIdx(pseudomol, owner_interactor_idx), bond_type);

    // Add the remainining insert MolBrick's Connections as peripheral connections.
    if (guided) {
      std::vector<ConnectionPoint*> cpoints;
      for (auto& c : brick_available_connections) {
        cpoints = connections.AddConnections(c.first, c.second);
        query_results.AssignFreshWeights(c.first, cpoints);
      };
    } else {
      for (auto& c : brick_available_connections) {
        connections.AddConnections(c.first, c.second);
      };
    };
  };

  // Flag the ReconstructedMol as modified.
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // If any chiral centers or stereochemical double bonds have unspecified
  // stereochemistry, assign a random one.
  if (update_stereo) {
    AssignUnspecifiedStereochemistry(prng);
  };

  if (guided) {
    assert(AllConnectionsHaveWeights());
  };

  // Signal function success.
  return true;
};

bool ReconstructedMol::StereoFlip(std::mt19937& prng) {
  // NOTE: This function will only work if the stereocenter and
  // stereo bond indices are up to date. This means that
  // AssignUnspecifiedStereochemistry must be called before.
  if (!stereo_updated) {
    AssignUnspecifiedStereochemistry(prng);
  };
  unsigned n_chiral_centers = chiral_center_atom_indices.size();
  unsigned n_stereo_bonds = stereo_bond_indices.size();
  unsigned n = n_chiral_centers + n_stereo_bonds;
  // If the molecule has no chiral centers or stereo bonds signal
  // function failure.
  if (n == 0) {
    return false;
  };
  // If it does, pick a random chiral center/stereo bond whose
  // stereochemistry to invert.
  unsigned stereo_idx = 0;
  if (n > 1) {
    std::uniform_int_distribution<unsigned> distribution(0, n - 1);
    stereo_idx = distribution(prng);
  };
  // If a chiral center was chosen:
  if (stereo_idx < n_chiral_centers) {
    unsigned atom_idx = chiral_center_atom_indices[stereo_idx];
    RDKit::Atom* chiral_center = pseudomol.getAtomWithIdx(atom_idx);
    RDKit::Atom::ChiralType chiral_tag = chiral_center->getChiralTag();
    assert(chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CW || chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CCW);
    if (chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CW) {
      chiral_center->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
    } else {
      chiral_center->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
    };
  // If a stereo double bond was chosen:
  } else {
    unsigned bond_idx = stereo_bond_indices[stereo_idx - n_chiral_centers];
    RDKit::Bond* stereo_bond = pseudomol.getBondWithIdx(bond_idx);
    RDKit::Bond::BondStereo stereo_tag = stereo_bond->getStereo();
    assert(stereo_tag == RDKit::Bond::BondStereo::STEREOE || stereo_tag == RDKit::Bond::BondStereo::STEREOZ);
    if (stereo_tag == RDKit::Bond::BondStereo::STEREOE) {
      stereo_bond->setStereo(RDKit::Bond::BondStereo::STEREOZ);
    } else {
      stereo_bond->setStereo(RDKit::Bond::BondStereo::STEREOE);
    };
  };

  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;

  // Signal function success.
  return true;
};

bool ReconstructedMol::Evolve(ReconstructionSettings& settings, sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, ReconstructionsInventory& inventory, EvolutionReport& report, unsigned n_failures) {
  // Determine which operations are possible. Checking whether a peripheral
  // deletion or expansion is possible is trivial and should be done. The
  // others I can afford to let fail.

  // Choose an evolutionary operator based on the probabilities defined in
  // the ReconstructionSettings.
  const std::string& operation = settings.ChooseOperation(prng);

  // std::cout << "Attempting operation: " << operation << std::endl;

  ++report.attempt_frequencies.at(operation);

  bool success = false;
  if (operation == "peripheral_expansion") {
    success = PeripheralExpansion(select_statement, query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "internal_expansion") {
    success = InternalExpansion(select_statement, query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "peripheral_deletion") {
    success = PeripheralDeletion(controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "internal_deletion") {
    success = InternalDeletion(query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "peripheral_substitution") {
    success = Substitution("peripheral", select_statement, query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "internal_substitution") {
    success = Substitution("internal", select_statement, query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "peripheral_crossover") {
    success = Crossover("peripheral", query_results, controller, prng, inventory, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "internal_crossover") {
    success = Crossover("peripheral", query_results, controller, prng, inventory, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "translation") {
    success = Translation(query_results, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
  } else if (operation == "stereo_flip") {
    success = StereoFlip(prng);
  };

  if (success) {
    // std::cout << "Operation succeeded!\n";
    ++report.success_frequencies.at(operation);
  };

  while (!success) {
    inventory.ClearQueues();
    ++n_failures;
    // std::cout << "Operation " << operation << " failed! Retrying..." << std::endl;
    if (n_failures >= settings.GetMaxAttemptsPerGeneration()) {
      return false;
    };
    success = Evolve(settings, select_statement, query_results, controller, prng, inventory, report, n_failures);
  };
  assert(success);
  return true;
};

void ReconstructedMol::AddLabelledPseudofragment(const Pseudofragment& pseudofragment) {
  // Determine which schematic index to use for the new MolBrick.
  unsigned brick_schematic_idx;
  if (available_schematic_idxs.empty()) {
    schematic.Expand();
    if (bricks.empty()) {
      brick_schematic_idx = 0u;
    } else {
      brick_schematic_idx = ++max_schematic_idx;
    };
  } else {
    brick_schematic_idx = available_schematic_idxs.back();
    available_schematic_idxs.pop_back();
  };

  // Construct the new MolBrick.
  MolBrick brick(pseudofragment, 0u, brick_schematic_idx, 1.0f, this);

  // Initialize iterators to check if a given Atom pertains to the new MolBrick.
  std::vector<unsigned>::const_iterator new_brick_it, new_brick_begin_it = brick.atom_indices.begin(), new_brick_end_it = brick.atom_indices.end();

  // Check which MolBricks of the ReconstructedMol neighbor the new MolBrick and
  // record the Connections between them.
  // Loop over the Atoms in the MolBrick.
  for (const RDKit::Atom* brick_atom : brick.pseudomol.atoms()) {
    // Get the corresponding Atom in the ReconstructedMol.
    unsigned immutable_atom_idx = brick_atom->getProp<unsigned>("ImmutableIdx");
    unsigned atom_type = brick_atom->getProp<unsigned>("MMFF_atom_type");
    const RDKit::Atom* reconstruction_atom = GetAtomWithImmutableIdx(pseudomol, immutable_atom_idx);
    // Loop over the neighboring atoms.
    RDKit::ROMol::ADJ_ITER nbr_it, nbr_end_it;
    std::tie(nbr_it, nbr_end_it) = pseudomol.getAtomNeighbors(reconstruction_atom);
    while (nbr_it != nbr_end_it) {
      unsigned neighbor_atom_idx = *nbr_it;
      const RDKit::Atom* neighbor_atom = pseudomol.getAtomWithIdx(neighbor_atom_idx);
      unsigned immutable_neighbor_atom_idx = neighbor_atom->getProp<unsigned>("ImmutableIdx");
      // If the neighboring atom isn't part of the same MolBrick, check to which
      // neighbor it pertains.
      new_brick_it = std::find(new_brick_begin_it, new_brick_end_it, immutable_neighbor_atom_idx);
      if (new_brick_it == new_brick_end_it) {
        // Iterate over the MolBricks in the ReconstructedMol.
        for (const auto& b : bricks) {
          unsigned neighbor_brick_schematic_idx = b.first;
          const MolBrick& neighbor_brick = b.second;
          // If the MolBrick contains the atom in question record it as a neighbor.
          std::vector<unsigned>::const_iterator it, end_it = neighbor_brick.atom_indices.end();
          it = std::find(neighbor_brick.atom_indices.begin(), end_it, immutable_neighbor_atom_idx);
          if (it != end_it) {
            // Retrieve the necessary attributes of the neighboring Atom to
            // update the ReconstructedMol.
            unsigned neighbor_atom_type = neighbor_atom->getProp<unsigned>("MMFF_atom_type");
            const RDKit::Bond* bond = pseudomol.getBondBetweenAtoms(reconstruction_atom->getIdx(), neighbor_atom_idx);
            unsigned bond_type_int = bond_type_int_table[bond->getBondType()];
            Connection connection (atom_type, neighbor_atom_type, bond_type_int);
            // Update the ReconstructionSchematic.
            schematic.adjacency[brick_schematic_idx][neighbor_brick_schematic_idx] = 1u;
            schematic.adjacency[neighbor_brick_schematic_idx][brick_schematic_idx] = 1u;
            schematic.start_atom_idx[brick_schematic_idx][neighbor_brick_schematic_idx] = immutable_atom_idx;
            schematic.start_atom_idx[neighbor_brick_schematic_idx][brick_schematic_idx] = immutable_neighbor_atom_idx;
            schematic.end_atom_idx[brick_schematic_idx][neighbor_brick_schematic_idx] = immutable_neighbor_atom_idx;
            schematic.end_atom_idx[neighbor_brick_schematic_idx][brick_schematic_idx] = immutable_atom_idx;
            schematic.start_atom_type[brick_schematic_idx][neighbor_brick_schematic_idx] = atom_type;
            schematic.start_atom_type[neighbor_brick_schematic_idx][brick_schematic_idx] = neighbor_atom_type;
            schematic.end_atom_type[brick_schematic_idx][neighbor_brick_schematic_idx] = neighbor_atom_type;
            schematic.end_atom_type[neighbor_brick_schematic_idx][brick_schematic_idx] = atom_type;
            schematic.bond_type[brick_schematic_idx][neighbor_brick_schematic_idx] = bond_type_int;
            schematic.bond_type[neighbor_brick_schematic_idx][brick_schematic_idx] = bond_type_int;
            // Update the ReconstructedMol's ConnectionsTables.
            internal_connections.AddConnection(connection, immutable_atom_idx);
            internal_connections.AddConnection(connection.Mirror(), immutable_neighbor_atom_idx);
            break;
          };
        };
      };
      ++nbr_it;
    };
  };

  // Update the rest of the ReconstructedMol's attributes.
  level += brick.level;
  if (brick.has_ring) {
    n_ring_atoms += brick.level;
  };
  if (brick.max_atom_idx > max_atom_idx) {
    max_atom_idx = brick.max_atom_idx;
  };
  sanitized_mol_updated = false;
  smiles_updated = false;
  sanitized_smiles_updated = false;
  fingerprint_updated = false;
  stereo_updated = false;

  // Store the new MolBrick.
  bricks.insert({brick_schematic_idx, std::move(brick)});
};

void ReconstructedMol::SetID(unsigned new_id) {
  id = new_id;
};

void ReconstructedMol::TakeBricksOwnership() {
  for (auto& brick : bricks) {
    brick.second.owner = this;
  };
};

std::pair<Connection, bool> ReconstructedMol::HasKnownConnections(const ConnectionQueryResults& query_results) const {
  const ConnectionCompatibilities& compatibilities = query_results.GetConnectionCompatibilities();
  for (const auto& c : connections) {
    if (!compatibilities.HasConnection(c.first)) {
      return std::make_pair(c.first, false);
    };
  };
  for (const auto& c : internal_connections) {
    if (!compatibilities.HasConnection(c.first)) {
      return std::make_pair(c.first, false);
    };
  };
  return std::make_pair(Connection(), true);
};

void ReconstructedMol::AssignWeightsToConnections(ConnectionQueryResults & query_results) {
  for (auto& c : connections) {
    query_results.AssignFreshWeights(c.first, c.second);
  };
  for (auto& c : internal_connections) {
    query_results.AssignFreshWeights(c.first, c.second);
  };
};

bool ReconstructedMol::AllConnectionsHaveWeights() const {
  for (const auto& c : connections) {
    for (const auto& cpoint : c.second) {
      if (cpoint.GetIDs().GetAcyclicIDs().empty() && cpoint.GetIDs().GetRingIDs().empty()) {
        return false;
      };
    };
  };
  for (const auto& c : internal_connections) {
    for (const auto& cpoint : c.second) {
      if (cpoint.GetIDs().GetAcyclicIDs().empty() && cpoint.GetIDs().GetRingIDs().empty()) {
        return false;
      };
    };
  };
  return true;
};

void ReconstructedMol::ClearConnectionWeights() {
  connections.ClearWeights();
  internal_connections.ClearWeights();
};

void ReconstructedMol::SetScore(float new_score) {
  previous_score = score;
  score = new_score;
};

void ReconstructedMol::SetChildFlag(bool flag) {
  is_child = flag;
};

void ReconstructedMol::SetParentID(unsigned new_parent_id) {
  parent_id = new_parent_id;
};

void ReconstructedMol::Sanitize() {
  if (sanitized_mol_updated) {
    return;
  };
  try {
    // Reset the sanitized molecule as a copy of the unsanitized pseudomol.
    sanitized_mol = pseudomol;
    // Add explicit hydrogens to replace the molecule's pseudoatoms.
    // (i.e. saturate the ReconstructedMol's Connections with hydrogens).
    for (const auto& c : connections) {
      const Connection& connection = c.first;
      const std::vector<ConnectionPoint>& cpoints = c.second;
      unsigned bond_order = connection.bond_type;
      assert(bond_order == 1 || bond_order == 2 || bond_order == 3);
      for (const ConnectionPoint& cpoint : cpoints) {
        RDKit::Atom* atom = GetAtomWithImmutableIdx(sanitized_mol, cpoint.atom_idx);
        atom->setNumExplicitHs(atom->getNumExplicitHs() + bond_order);
      };
    };
    // Loop over the molecule in search of Sulfur and Phosphorus atoms and adjust
    // their number of explicit hydrogens to match their lowest possible valence.
    // NOTE: If working with other elements for which simple hydrogen saturation
    // isn't guaranteed to yield a legal chemotype this section ought to be expanded.
    for (RDKit::RWMol::AtomIterator ai = sanitized_mol.beginAtoms(); ai != sanitized_mol.endAtoms(); ++ai) {
      // If the atom is a phosphorus:
      if ((*ai)->getAtomicNum() == 15) {
        // Calculate the "heavy atom valence" as the number of bonds involved in
        // bonding to heavy atoms.
        unsigned valence = (*ai)->getTotalValence();
        unsigned n_hydrogens = (*ai)->getTotalNumHs();
        unsigned heavy_atom_valence = valence - n_hydrogens;
        // assert(valence >= 0 && valence <= 5);
        assert(valence >= 0 && valence <= 6);
        if (valence == 3 || valence == 5) {
          // If this valence is 0 (i.e. the phosphorus is detached) convert it to phosphine.
          if (heavy_atom_valence == 0) {
            (*ai)->setNumExplicitHs(3);
          // If the valence is odd it's already valid and the excess hydrogens are removed.
          } else if (heavy_atom_valence % 2) {
            (*ai)->setNumExplicitHs(0);
          // If the valence is even push it up to the nearest odd valence.
          } else {
            (*ai)->setNumExplicitHs(1);
          };
        };
      // If the atom is a sulfur:
      } else if ((*ai)->getAtomicNum() == 16) {
        // Calculate the "heavy atom valence" as the number of bonds involved in
        // bonding to heavy atoms.
        unsigned valence = (*ai)->getTotalValence();
        unsigned n_hydrogens = (*ai)->getTotalNumHs();
        unsigned heavy_atom_valence = valence - n_hydrogens;
        assert(valence >= 0 && valence <= 6);
        if (valence == 2 || valence == 4 || valence == 6) {
          // If this valence is 0 (i.e. the sulfur is detached) convert it to hydrogen sulfide.
          if (heavy_atom_valence == 0) {
            (*ai)->setNumExplicitHs(2);
          // If the valence is odd push it up to the nearest even valence.
          } else if (heavy_atom_valence % 2) {
            (*ai)->setNumExplicitHs(1);
          // If the valence is even it's already valid and the excess hydrogens are removed.
          } else {
            (*ai)->setNumExplicitHs(0);
          };
        };
      };
    };
    // Wrap up the molecule sanitization and remove the hydrogens.
    // This is important since explicit hydrogens affect the fingerprint.
    RDKit::MolOps::sanitizeMol(sanitized_mol);
    RDKit::MolOps::removeHs(sanitized_mol);
    sanitized_mol_updated = true;
  } catch (const std::exception& e) {
    std::cout << "ERROR: Sanitization failed for ReconstructedMol " << id << ": " << GetSMILES() << std::endl;
    throw;
  };
};

void ReconstructedMol::AssignUnspecifiedStereochemistry(std::mt19937& prng) {
  // Stereochemistry can only be inferred from sanitized molecules. If
  // the molecule isn't sanitized yet, do it now.
  Sanitize();
  // Create a Bernoulli distribution to select a random stereoisomer.
  // For chirality:
  //  0 = CW = R
  //  1 = CCW = S
  // For bond configurations:
  //  0 = Z = CIS
  //  1 = E = TRANS
  std::bernoulli_distribution distribution(0.5);
  // Initialize/reset the output vectors
  std::vector<unsigned> cip_ranks;
  std::vector<std::pair<int, int>> neighbors;
  chiral_center_atom_indices.clear();
  stereo_bond_indices.clear();
  // Assign CIP ranks to the sanitized molcule's atoms.
  RDKit::Chirality::assignAtomCIPRanks(sanitized_mol, cip_ranks);
  // Search for chiral centers in the sanitized molecule.
  for (RDKit::RWMol::AtomIterator ai = sanitized_mol.beginAtoms(); ai != sanitized_mol.endAtoms(); ++ai) {
    // Get the equivalent atom in the non-sanitized molecule.
    RDKit::Atom* sanitized_atom = *ai;
    unsigned atom_idx = sanitized_atom->getIdx();
    RDKit::Atom* atom = pseudomol.getAtomWithIdx(atom_idx);
    bool legal_center = false, has_dupes = true;
    std::tie(legal_center, has_dupes) = RDKit::Chirality::isAtomPotentialChiralCenter(*ai, sanitized_mol, cip_ranks, neighbors);
    neighbors.clear();
    // If the atom is a chiral center without specified chirality set one randomly.
    if (legal_center && !has_dupes) {
      if (sanitized_atom->getChiralTag() == RDKit::Atom::CHI_UNSPECIFIED) {
        if (distribution(prng)) {
          sanitized_atom->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
          atom->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CCW);
        } else {
          sanitized_atom->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
          atom->setChiralTag(RDKit::Atom::CHI_TETRAHEDRAL_CW);
        };
      };
      // Store the chiral atom's index.
      chiral_center_atom_indices.push_back(atom_idx);
    // If the atom isn't a chiral center make sure its chiral tag is set to "unspecified".
    } else {
      sanitized_atom->setChiralTag(RDKit::Atom::CHI_UNSPECIFIED);
      atom->setChiralTag(RDKit::Atom::CHI_UNSPECIFIED);
    };
  };
  // Find stereochemistry-definining double bonds and randomly assign a bond
  // configuration to the unspecified ones. Those previously flagged as stereo
  // bonds and that no longer fulfill the criteria are cleared of their flag.
  RDKit::MolOps::findPotentialStereoBonds(sanitized_mol, true);
  // Loop over the bonds in the sanitized molecule.
  for (RDKit::RWMol::BondIterator bi = sanitized_mol.beginBonds(); bi != sanitized_mol.endBonds(); ++bi) {
    // Get the equivalent bond in the non-sanitized molecule.
    RDKit::Bond* sanitized_bond = *bi;
    unsigned bond_idx = sanitized_bond->getIdx();
    RDKit::Bond* bond = pseudomol.getBondWithIdx(bond_idx);
    RDKit::Bond::BondStereo bond_stereo = sanitized_bond->getStereo();
    // If the bond isn't a stereochemical bond skip it.
    if (bond_stereo == RDKit::Bond::BondStereo::STEREONONE) {
      continue;
    };
    // If the bond is a stereochemical bond without defined stereochemistry, assign a random one.
    if (bond_stereo == RDKit::Bond::BondStereo::STEREOANY) {
      if (distribution(prng)) {
        sanitized_bond->setStereo(RDKit::Bond::BondStereo::STEREOE);
        bond->setStereo(RDKit::Bond::BondStereo::STEREOE);
      } else {
        sanitized_bond->setStereo(RDKit::Bond::BondStereo::STEREOZ);
        bond->setStereo(RDKit::Bond::BondStereo::STEREOE);
      };
    };
    stereo_bond_indices.push_back(bond_idx);
  };
  stereo_updated = true;
};

void ReconstructedMol::GenerateSMILES() {
  if (smiles_updated) {
    return;
  };
  RDKit::RWMol dummymol(pseudomol);
  RDKit::Bond::BondType bond_type;
  RDKit::Atom pseudoatom(0);
  unsigned start_atom_idx, pseudoatom_idx;
  for (const auto& c : connections) {
    bond_type = int_bond_type_table[c.first.bond_type];
    pseudoatom.setIsotope(c.first.encoded);
    for (const ConnectionPoint& cpoint : c.second) {
      start_atom_idx = GetMutableAtomIdx(dummymol, cpoint.atom_idx);
      pseudoatom_idx = dummymol.addAtom(&pseudoatom);
      dummymol.addBond(start_atom_idx, pseudoatom_idx, bond_type);
    };
  };
  smiles = RDKit::MolToSmiles(dummymol);
  smiles_updated = true;
};

void ReconstructedMol::GenerateSanitizedSMILES() {
  if (sanitized_smiles_updated) {
    return;
  };
  Sanitize();
  sanitized_smiles = RDKit::MolToSmiles(sanitized_mol);
  sanitized_smiles_updated = true;
};

void ReconstructedMol::GenerateFingerprint() {
  if (fingerprint_updated) {
    return;
  };
  Sanitize();
  // WARNING: At the time of writing (RDKit version 2020.09.4) the
  // RDKit::SparseIntVect copy constructor  and assignment operators were bugged.
  // By the time LEADD is released this  should be fixed, but beware if you are
  // using older versions of the RDKit.
  const RDKit::SparseIntVect<std::uint32_t>* tmp = RDKit::MorganFingerprints::getFingerprint(sanitized_mol, 2);
  fingerprint = *tmp;
  delete tmp;
  fingerprint_updated = true;
};

float ReconstructedMol::Score(const RDKit::SparseIntVect<std::uint32_t>* reference) {
  GenerateFingerprint();
  previous_score = score;
  score = RDKit::TanimotoSimilarity(fingerprint, *reference);
  return score;
};

unsigned ReconstructedMol::GetID() const {
  return id;
};

double ReconstructedMol::GetSimilarity(ReconstructedMol& reconstruction) {
  GenerateFingerprint();
  reconstruction.GenerateFingerprint();
  return RDKit::TanimotoSimilarity(fingerprint, reconstruction.fingerprint);
};

bool ReconstructedMol::StereoIsUpdated() const {
  return stereo_updated;
};

const RDKit::RWMol& ReconstructedMol::GetPseudomol() const {
  return pseudomol;
};

const RDKit::RWMol& ReconstructedMol::GetSanitizedMol() {
  Sanitize();
  return sanitized_mol;
};

const ConnectionsTable& ReconstructedMol::GetConnections() const {
  return connections;
};

const ConnectionsTable& ReconstructedMol::GetInternalConnections() const {
  return internal_connections;
};

const std::unordered_map<unsigned, MolBrick>& ReconstructedMol::GetBricks() const {
  return bricks;
};

const ReconstructionSchematic& ReconstructedMol::GetSchematic() const {
  return schematic;
};

unsigned ReconstructedMol::GetNBricks() const {
  return bricks.size();
};

unsigned ReconstructedMol::GetLevel() const {
  return level;
};

unsigned ReconstructedMol::GetNRingAtoms() const {
  return n_ring_atoms;
};

unsigned ReconstructedMol::GetMaxAtomIdx() const {
  return max_atom_idx;
};

unsigned ReconstructedMol::GetMaxSchematicIdx() const {
  return max_schematic_idx;
};

float ReconstructedMol::GetScore() const {
  return score;
};

const std::string& ReconstructedMol::GetSMILES() {
  GenerateSMILES();
  return smiles;
};

const std::string& ReconstructedMol::GetSanitizedSMILES() {
  GenerateSanitizedSMILES();
  return sanitized_smiles;
};

const RDKit::SparseIntVect<std::uint32_t>& ReconstructedMol::GetFingerprint() {
  GenerateFingerprint();
  return fingerprint;
};

bool ReconstructedMol::IsChild() const {
  return is_child;
};

unsigned ReconstructedMol::GetParentID() const {
  return parent_id;
};

void ReconstructedMol::Draw(const std::string& file_path) {
  Sanitize();
  RDDepict::compute2DCoords(sanitized_mol);
  for (RDKit::RWMol::AtomIterator ai = sanitized_mol.beginAtoms(); ai != sanitized_mol.endAtoms(); ++ai) {
    (*ai)->setProp<unsigned>(RDKit::common_properties::atomNote, (*ai)->getProp<unsigned>("ImmutableIdx"));
  };
  std::ofstream output_stream (file_path);
  RDKit::MolDraw2DSVG drawer (300, 300, output_stream);
  drawer.drawMolecule(sanitized_mol);
  drawer.finishDrawing();
  output_stream.close();
};


// Class EvolutionGuide
EvolutionGuide::EvolutionGuide() = default;
EvolutionGuide::EvolutionGuide(sqlite3* database, const ReconstructionSettings& settings) {
  // Set the learning reinforcements and similarity thresholds.
  acyclic_threshold = settings.GetAcyclicLearningSimilarityThreshold();
  ring_threshold = settings.GetRingLearningSimilarityThreshold();
  acyclic_positive_reinforcement = settings.GetAcyclicPositiveReinforcement();
  acyclic_negative_reinforcement = settings.GetAcyclicNegativeReinforcement();
  ring_positive_reinforcement = settings.GetRingPositiveReinforcement();
  ring_negative_reinforcement = settings.GetRingNegativeReinforcement();
  assert(acyclic_threshold >= 0.0f && acyclic_threshold <= 1.0f);
  assert(acyclic_positive_reinforcement >= 0.0f);
  assert(acyclic_negative_reinforcement >= 0.0f && acyclic_negative_reinforcement <= 1.0f);
  assert(ring_positive_reinforcement >= 0.0f);
  assert(ring_negative_reinforcement >= 0.0f && ring_negative_reinforcement <= 1.0f);
  // Count the number of Pseudofragments in the database.
  unsigned n_fragments = GetPseudofragmentCount(database);
  // Set up the necessary handles to work with the HDF5 similarity file.
  simatrix_file = H5::H5File(settings.GetSimilarityMatrixFile(), H5F_ACC_RDONLY);
  simatrix_dataset = simatrix_file.openDataSet("simatrix");
  simatrix_dataspace = simatrix_dataset.getSpace();
  float_type = simatrix_dataset.getFloatType();
  // Configure the size of the selection hyperslab.
  hyperslab_size[0] = 1;
  hyperslab_size[1] = n_fragments;
  // Set up buffers to store an entire row of read values from the HDF5 file.
  buffer1 = new float[n_fragments];
  buffer2 = new float[n_fragments];
  buffer3 = new float[n_fragments];
  buffer_size[0] = n_fragments;
  buffer_dataspace = H5::DataSpace(1, buffer_size);
};

void EvolutionGuide::AdjustWeights(ConnectionPoint* cpoint, float* buffer, float sign) {
  // Set the reinforcement rate based on the sign of the weight change.
  float acyclic_reinforcement = acyclic_positive_reinforcement;
  float ring_reinforcement = ring_positive_reinforcement;
  if (sign < 0.0) {
    acyclic_reinforcement = acyclic_negative_reinforcement;
    ring_reinforcement = ring_negative_reinforcement;
  };
  if (acyclic_reinforcement > 0.0f) {
    for (size_t idx = 0; idx < cpoint->ids.acyclic.size(); ++idx) {
      unsigned id = cpoint->ids.acyclic[idx];
      float similarity = buffer[id - 1];
      if (similarity >= acyclic_threshold) {
        cpoint->weights.acyclic[idx] *= (1 + sign * acyclic_reinforcement * similarity);
      };
    };
  };
  if (ring_reinforcement > 0.0f) {
    for (size_t idx = 0; idx < cpoint->ids.ring.size(); ++idx) {
      unsigned id = cpoint->ids.ring[idx];
      float similarity = buffer[id - 1];
      if (similarity >= ring_threshold) {
        cpoint->weights.ring[idx] *= (1 + sign * ring_reinforcement * similarity);
      };
    };
  };
};

void EvolutionGuide::AdjustWeights(std::list<ReconstructedMol>& reconstructions, const std::set<unsigned>& survivor_ids) {
  // Loop over the child ReconstructedMols.
  std::set<unsigned>::const_iterator survivors_end_it = survivor_ids.end();
  for (ReconstructedMol& child : reconstructions) {
    if (child.IsChild()) {
      // If either the child or its parent survived, the child's GeneticLog is
      // used to update the ConnectionPoint weights of the child and/or parent
      // respectively.
      bool child_survived = false;
      bool parent_survived = false;
      if (survivor_ids.find(child.GetID()) != survivors_end_it) {
        child_survived = true;
      };
      if (survivor_ids.find(child.GetParentID()) != survivors_end_it) {
        parent_survived = true;
      };

      if (child_survived || parent_survived) {
        // Get the child's parent.
        unsigned parent_id = child.GetParentID();
        std::list<ReconstructedMol>::iterator it, end_it = reconstructions.end();
        it = std::find_if(reconstructions.begin(), end_it,
          [=](const ReconstructedMol& reconstruction) {
            return reconstruction.GetID() == parent_id;
          });
        if (it == end_it) {
          throw std::runtime_error("Parent ReconstructedMol of a child ReconstructedMol couldn't be found.");
        };
        ReconstructedMol& parent = *it;

        // Determine how the last operation affected the score.
        bool score_increased = false;
        float sign = 1.0;
        // If the score wasn't affected the weight of the involved fragment is
        // always decreased to push the algorithm towards innovative solutions.
        // This helps escape local minima.
        if (child.score == child.previous_score) {
          sign = -1.0;
        } else if (child.score > child.previous_score) {
          score_increased = true;
        };

        // If the child's GeneticLog defines some important Pseudofragments,
        // read their similarity data into their private buffers.
        // Note: -1 because SQLite3 is 1-indexed while HDF5 is 0-indexed.
        const GeneticLog& log = child.genetic_log;
        unsigned pseudofragment_id1 = log.GetImportantPseudofragmentID1();
        unsigned pseudofragment_id2 = log.GetImportantPseudofragmentID2();
        if (pseudofragment_id1 > 0u && pseudofragment_id1 != buffer1_id) {
          hyperslab_offset[0] = pseudofragment_id1 - 1;
          simatrix_dataspace.selectHyperslab(H5S_SELECT_SET, hyperslab_size, hyperslab_offset);
          simatrix_dataset.read(buffer1, float_type, buffer_dataspace, simatrix_dataspace);
          buffer1_id = pseudofragment_id1;
        };
        if (pseudofragment_id2 > 0u  && pseudofragment_id2 != buffer2_id) {
          hyperslab_offset[0] = pseudofragment_id2 - 1;
          simatrix_dataspace.selectHyperslab(H5S_SELECT_SET, hyperslab_size, hyperslab_offset);
          simatrix_dataset.read(buffer2, float_type, buffer_dataspace, simatrix_dataspace);
          buffer2_id = pseudofragment_id2;
        };

        // Iterate over the GeneticLogEntries, stratified according to the involved
        // Pseudofragment's ID.
        for (const auto& e : log.GetEntries()) {
          unsigned pseudofragment_id = e.first;
          const std::vector<GeneticLogEntry>& entries = e.second;

          // Choose the appropiate buffer for the Pseudofragment.
          float* buffer = nullptr;
          // If the involved Pseudofragment already has an associated pre-loaded buffer,
          // select it.
          if (pseudofragment_id == buffer1_id) {
            buffer = buffer1;
          } else if (pseudofragment_id == buffer2_id) {
            buffer = buffer2;
          } else if (pseudofragment_id == buffer3_id) {
            buffer = buffer3;
          // Otherwise, load the similarity data into buffer 3.
          } else {
            hyperslab_offset[0] = pseudofragment_id - 1;
            simatrix_dataspace.selectHyperslab(H5S_SELECT_SET, hyperslab_size, hyperslab_offset);
            simatrix_dataset.read(buffer3, float_type, buffer_dataspace, simatrix_dataspace);
            buffer3_id = pseudofragment_id;
            buffer = buffer3;
          };

          for (const GeneticLogEntry& entry : entries) {
            GeneticLogEntry::Type type = entry.GetType();
            const Connection& connection = entry.GetConnection();
            unsigned atom_idx = entry.GetAtomIdx();

            // Determine the sign of the weight change based on the genetic
            // operation outcome and type.
            if (score_increased) {
              if (type == GeneticLogEntry::Type::DELETION) {
                sign = -1.0;
              };
            } else {
              if (type == GeneticLogEntry::Type::CREATION) {
                sign = -1.0;
              };
            };

            // If the parent survived, update its weights.
            // Retrieve a pointer to the ConnectionPoints referenced by the entry.
            // If multiple ConnectionPoints share an atom index this will be > 1.
            if (parent_survived) {
              std::vector<ConnectionPoint*> cpoints;
              if (parent.connections.HasConnection(connection, atom_idx)) {
                cpoints = parent.connections.GetConnectionPointsWithIdx(connection, atom_idx);
              } else if (parent.internal_connections.HasConnection(connection, atom_idx)) {
                cpoints = parent.internal_connections.GetConnectionPointsWithIdx(connection, atom_idx);
              } else {
                throw std::runtime_error("ConnectionPoint referenced by GeneticLog doesn't exist.");
              };
              for (ConnectionPoint* cpoint : cpoints) {
                AdjustWeights(cpoint, buffer, sign);
              };
            };

            // If the child survived, update its weights.
            if (child_survived) {
              std::vector<ConnectionPoint*> cpoints;
              if (child.connections.HasConnection(connection, atom_idx)) {
                cpoints = child.connections.GetConnectionPointsWithIdx(connection, atom_idx);
              } else if (child.internal_connections.HasConnection(connection, atom_idx)) {
                cpoints = child.internal_connections.GetConnectionPointsWithIdx(connection, atom_idx);
              } else {
                throw std::runtime_error("ConnectionPoint referenced by GeneticLog doesn't exist.");
              };
              for (ConnectionPoint* cpoint : cpoints) {
                AdjustWeights(cpoint, buffer, sign);
              };
            };
          };
        };
      };
    };
  };
};

void EvolutionGuide::Cleanup() {
  float_type.close();
  simatrix_dataspace.close();
  simatrix_dataset.close();
  simatrix_file.close();
  buffer_dataspace.close();
  delete buffer1;
  delete buffer2;
  delete buffer3;
};


// Standalone function to retrieve a pointer to the Atom with the specified immutable index.
RDKit::Atom* GetAtomWithImmutableIdx(RDKit::RWMol& pseudomol, unsigned immutable_atom_idx) {
  for (RDKit::RWMol::AtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    if ((*ai)->getProp<unsigned>("ImmutableIdx") == immutable_atom_idx) {
      return *ai;
    };
  };
  throw std::runtime_error("Atom retrieval failed.\n");
};

// Standalone functions to convert pseudomol immutable atom indices into their mutable counterparts.
unsigned GetMutableAtomIdx(const RDKit::RWMol & pseudomol, unsigned immutable_atom_idx) {
  for (RDKit::ROMol::ConstAtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    if ((*ai)->getProp<unsigned>("ImmutableIdx") == immutable_atom_idx) {
      return (*ai)->getIdx();
    };
  };
  throw std::runtime_error("Atom idx retrieval failed.\n");
};

std::vector<unsigned> GetMutableAtomIndices(const RDKit::RWMol & pseudomol, const std::vector<unsigned> & immutable_atom_indices) {
  std::vector<unsigned> mutable_atom_indices;
  mutable_atom_indices.reserve(immutable_atom_indices.size());
  std::vector<unsigned>::const_iterator it, begin_it = immutable_atom_indices.begin(), end_it = immutable_atom_indices.end();
  unsigned immutable_atom_idx;
  for (RDKit::ROMol::ConstAtomIterator ai = pseudomol.beginAtoms(); ai != pseudomol.endAtoms(); ++ai) {
    immutable_atom_idx = (*ai)->getProp<unsigned>("ImmutableIdx");
    it = begin_it;
    while ((it = std::find(it, end_it, immutable_atom_idx)) != end_it) {
      mutable_atom_indices.push_back((*ai)->getIdx());
      ++it;
    };
  };
  return mutable_atom_indices;
};

ReconstructedMol ConvertToReconstructedMol(const RDKit::ROMol& mol, bool fragment_rings, bool names_as_scores, boost::format& formatter) {
  // Initialize an empty ReconstructedMol.
  ReconstructedMol reconstruction;
  // Set its molecule objects to the input RDKit::ROMol.
  reconstruction.pseudomol = mol;
  FlagAtoms(reconstruction.pseudomol);
  reconstruction.sanitized_mol = reconstruction.pseudomol;
  // If the ring systems ought to be fragmented:
  if (fragment_rings) {
    // Loop over the atoms in the molecule, convert them into Pseudofragments
    // and add them to the ReconstructedMol as MolBricks.
    for (const RDKit::Atom* atom : reconstruction.pseudomol.atoms()) {
      Pseudofragment pseudofragment = AtomToPseudofragment(reconstruction.pseudomol, atom, formatter);
      reconstruction.AddLabelledPseudofragment(pseudofragment);
    };
  // If ring systems should be treated as separate Pseudofragments.
  } else {
    // Separate the molecule into acyclic and cyclic fragments.
    std::vector<RDKit::ROMOL_SPTR> acyclic_fragments, ring_fragments;
    std::tie(acyclic_fragments, ring_fragments) = SplitAtRings(reconstruction.pseudomol);
    // Convert the cyclic fragments into Pseudofragments and add them to the
    // ReconstructedMol as MolBricks.
    for (RDKit::ROMOL_SPTR ring_fragment : ring_fragments) {
      Pseudofragment pseudofragment = DummyMolToPseudofragment(*ring_fragment, formatter);
      reconstruction.AddLabelledPseudofragment(pseudofragment);
    };
    // Loop over the acyclic atoms in the molecule, convert them into Pseudofragments
    // and add them to the ReconstructedMol as MolBricks.
    for (const RDKit::Atom* atom : reconstruction.pseudomol.atoms()) {
      if (!atom->getProp<bool>("WasInRing")) {
        Pseudofragment pseudofragment = AtomToPseudofragment(reconstruction.pseudomol, atom, formatter);
        reconstruction.AddLabelledPseudofragment(pseudofragment);
      };
    };
  };
  // If required, use the molecule's original name as its score.
  if (names_as_scores) {
    float score = std::stof(mol.getProp<std::string>("_Name"));
    reconstruction.SetScore(score);
  };

  return reconstruction;
};
