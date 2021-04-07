#pragma once
#ifndef _RECONSTRUCTION_HPP_
#define _RECONSTRUCTION_HPP_

#include <fstream>
#include <tuple>
#include <filesystem>
#include <assert.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include "H5Cpp.h"
#include "Graph.hpp"
#include "PseudofragmentDB.hpp"
#include "ConnectionQueryResults.hpp"
#include "ReconstructionSettings.hpp"
#include "Fragmentation.hpp"

// Declaration of a RDKit::Chirality function. This function is implemented in Chirality.cpp,
// but for some reason isn't made available to the user in Chirality.h.
namespace RDKit {
  namespace Chirality {
    typedef std::pair<int, int> INT_PAIR;
    typedef std::vector<INT_PAIR> INT_PAIR_VECT;
    std::pair<bool, bool> isAtomPotentialChiralCenter(const Atom* atom, const ROMol& mol, const UINT_VECT& ranks, INT_PAIR_VECT& nbrs);
  };
};

class MolBrick;
class ReconstructedMol;

// Type to store MolBrick pointers and their score-adjusted weights. Intended
// for use when retrieving substitutes from the ReconstructionsInventory.
typedef std::pair<std::vector<MolBrick*>, std::vector<float>> INVENTORY_RESULT;

typedef std::vector<MolBrick*> MOLBRICKS_COMBINATION;

// Declaration of the version integer of the ReconstructedMol & associated classes.
extern const unsigned reconstruction_version;

// Declaration of a variable that keeps track of whether the RDKit::MolPickler was configured.
extern bool molpickler_configured;

// Declaration of a function to configure the RDKit::MolPickler to pickle RDKit::Atom properties.
void ConfigureMolpickler();

class SizeController {
    // Maximum value covered by the discrete distribution
    unsigned max_size;
    // Normal distribution parameters
    float mean, stdev;
    float mean_squared, two_mean, two_stdev_squared, denominator;
    // Peak discrete distribution probability of no change
    float peak;
    // Logistic function growth rate
    float k;
    // Event probabilities
    std::unordered_map<unsigned, float> change_probabilities;
    std::unordered_map<unsigned, float> growth_probabilities;
    // Event decision distributions
    std::unordered_map<unsigned, std::discrete_distribution<unsigned>> change_distributions;
    std::unordered_map<unsigned, std::discrete_distribution<unsigned>> growth_distributions;

  public:
    // Don't use the default constructor! It's provided only to be able to
    // initialize the class without using initializer lists.
    SizeController();
    SizeController(unsigned max_size, float mean, float stdev, float peak);
    int Decide(unsigned size, std::mt19937& prng);
    bool DecideChange(unsigned size, std::mt19937& prng);
    bool DecideGrowth(unsigned size, std::mt19937& prng);

  private:
    float NormalProbabilityDensity(unsigned size) const;
    float LogisticValue(unsigned size) const;

  friend class ReconstructedMol;
};

class EvolutionReport {
    std::filesystem::path directory;
    std::ofstream scores_csv;
    std::ofstream n_ring_atoms_csv;
    std::ofstream operation_frequencies_csv;

    std::map<std::string, unsigned> attempt_frequencies{
      {"peripheral_expansion", 0},    {"internal_expansion", 0},
      {"peripheral_deletion", 0},     {"internal_deletion", 0},
      {"peripheral_substitution", 0}, {"internal_substitution", 0},
      {"peripheral_crossover", 0},    {"internal_crossover", 0},
      {"translation", 0},             {"stereo_flip", 0}
    };
    std::map<std::string, unsigned> success_frequencies{
      {"peripheral_expansion", 0},    {"internal_expansion", 0},
      {"peripheral_deletion", 0},     {"internal_deletion", 0},
      {"peripheral_substitution", 0}, {"internal_substitution", 0},
      {"peripheral_crossover", 0},    {"internal_crossover", 0},
      {"translation", 0},             {"stereo_flip", 0}
    };

  public:
    EvolutionReport();
    EvolutionReport(const std::string& output_directory_path);

    void WriteScores(unsigned generation, const std::list<ReconstructedMol>& reconstructions);
    void WriteNRings(unsigned generation, const std::list<ReconstructedMol>& reconstructions);
    void WriteOperationFrequencies();
    void WriteWeights(const ReconstructedMol& reconstruction, const std::string& file_name);

    void Close();

  friend class ReconstructedMol;
};

class MolBrick : public Pseudofragment {
  unsigned pseudofragment_id = 0u;
  unsigned max_atom_idx = 0u, schematic_idx = 0u;
  std::vector<unsigned> atom_indices;
  // The weight of the MolBrick is the weight of the source
  // Pseudofragment based on its frequency in the database
  // and its size.
  float weight = 0.0f;
  ReconstructedMol* owner = nullptr;

public:
  MolBrick();
  MolBrick(const Pseudofragment& pseudofragment, unsigned pseudofragment_id, unsigned schidx = 0, float wgt = 0, ReconstructedMol* owner = nullptr);
  MolBrick(const Pseudofragment& pseudofragment, unsigned pseudofragment_id, unsigned start_atom_idx = 0, unsigned schidx = 0, float wgt = 0, ReconstructedMol* owner = nullptr);
  MolBrick(Pseudofragment& pseudofragment, unsigned pseudofragment_id, unsigned start_atom_idx = 0, unsigned schidx = 0, float wgt = 0, ReconstructedMol* owner = nullptr);

  bool operator== (const MolBrick& other) const;

  bool IsPeripheral() const;
  bool IsInternal() const;

  std::vector<MolBrick*> GetNeighbors() const;
  ConnectionsTable GetAvailableConnections(std::mt19937& prng, const MolBrick* ghost_brick = nullptr, bool transfer_weights = false) const;
  ConnectionsTable GetAvailableConnections(std::mt19937& prng, const std::vector<MolBrick*>& ghost_bricks, bool transfer_weights = false) const;

  void UpdateImmutableIndices(unsigned start_atom_idx);

  const ReconstructedMol* GetOwner() const;
  unsigned GetFragmentID() const;
  unsigned GetMaxAtomIdx() const;
  unsigned GetSchematicIdx() const;
  const std::vector<unsigned>& GetImmutableAtomIndices() const;
  float GetWeight() const;

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = reconstruction_version) {
    ar & boost::serialization::base_object<Pseudofragment>(*this);
    ar & pseudofragment_id & max_atom_idx & schematic_idx & atom_indices & weight;
  };

  friend class InsertionPoint;
  friend class DeletionPoint;
  friend class SubstitutionPoint;
  friend class ReconstructedMol;
  friend class ReconstructionsInventory;
  friend class boost::serialization::access;
};

class ReconstructionSchematic {
  std::vector<std::vector<unsigned>> adjacency, start_atom_idx, end_atom_idx, start_atom_type, end_atom_type, bond_type;
  unsigned dimension = 0;

public:
  ReconstructionSchematic();
  ReconstructionSchematic(unsigned size);

  unsigned GetDimension() const;
  const std::vector<std::vector<unsigned>>& GetAdjacencyMatrix() const;
  const std::vector<std::vector<unsigned>>& GetStartAtomIdxMatrix() const;
  const std::vector<std::vector<unsigned>>& GetEndAtomIdxMatrix() const;
  const std::vector<std::vector<unsigned>>& GetStartAtomTypeMatrix() const;
  const std::vector<std::vector<unsigned>>& GetEndAtomTypeMatrix() const;
  const std::vector<std::vector<unsigned>>& GetBondTypeMatrix() const;

  void PrintAdjacencyMatrix() const;
  void PrintStartAtomIdxMatrix() const;
  void PrintEndAtomIdxMatrix() const;
  void PrintStartAtomTypeMatrix() const;
  void PrintEndAtomTypeMatrix() const;
  void PrintBondTypeMatrix() const;

private:
  void Expand(unsigned bump = 1);
  void UnitSchematic();
  void ClearCell(unsigned row_idx, unsigned column_idx);
  void ClearRow(unsigned row_idx);
  void ClearColumn(unsigned column_idx);

  template <class Archive>
  void serialize(Archive& ar, const unsigned version = reconstruction_version) {
    ar & adjacency & start_atom_idx & end_atom_idx & start_atom_type & end_atom_type & bond_type & dimension;
  };

  friend class MolBrick;
  friend class ReconstructedMol;
  friend class boost::serialization::access;
};

class ReconstructionsInventory {
  std::vector<MolBrick> acyclic, ring, ring_part;
  std::vector<float> acyclic_weights, ring_weights, ring_part_weights;
  std::vector<MolBrick> addition_queue;
  std::vector<MolBrick> deletion_queue;
  float gamma = 1;

public:
  ReconstructionsInventory();
  ReconstructionsInventory(float gamma);
  ReconstructionsInventory(const std::list<ReconstructedMol>& reconstructions, float gamma);

  void Add(const MolBrick& brick);
  void Add(const ReconstructedMol& reconstruction);
  void Add(const ReconstructedMol* reconstruction);
  void Remove(const MolBrick& brick);
  void Remove(const ReconstructedMol& reconstruction);
  void Remove(const ReconstructedMol* reconstruction);

  void SubmitToAdditionQueue(const MolBrick* brick);
  void SubmitToAdditionQueue(const MolBrick& brick);
  void SubmitToDeletionQueue(const MolBrick* brick);
  void SubmitToDeletionQueue(const MolBrick& brick);
  void ClearQueues();
  void CommitQueuedOperations();

  void CalcRingWeights();
  void CalcRingPartWeights();
  void CalcAcyclicWeights();
  void CalcWeights();

  void SetGamma(float new_gamma);
  void Clear();

  const std::vector<MolBrick>& GetRingBricks() const;
  const std::vector<MolBrick>& GetRingPartBricks() const;
  const std::vector<MolBrick>& GetAcyclicBricks() const;
  const std::vector<float>& GetRingBrickWeights() const;
  const std::vector<float>& GetRingPartBrickWeights() const;
  const std::vector<float>& GetAcyclicBrickWeights() const;

  friend class SubstitutionPoint;
  friend class ReconstructedMol;
};

class DeletionPoint {
  MolBrick* deletable;
  bool has_ring, ring_part;
  float weight;

  std::vector<MolBrick*> neighbors, evaluated_neighbors;
  std::map<const MolBrick*, ConnectionsTable> neighbors_available_connections;
  std::map<const MolBrick*, std::map<const MolBrick*, const Connection*>> other_neighbors_matchings;

  unsigned n_neighbors = 0, n_other_neighbors = 0;
  bool evaluated = false, is_deletable = false, generated_all_matchings = false;

public:
  DeletionPoint(MolBrick* deletable);

  void PerceiveNeighborsAvailableConnections(std::mt19937& prng);
  void Evaluate(ConnectionCompatibilities& compatibilities, std::mt19937& prng);
  void GenerateMatchings(ConnectionCompatibilities& compatibilities, std::mt19937& prng);

  friend class ReconstructedMol;
};

class SubstitutionPoint {
  MolBrick* substitutable;
  bool has_ring, ring_part;

  unsigned n_neighbors = 0;
  MOLBRICKS_COMBINATION neighbors_combination;
  std::map<MolBrick*, ConnectionsTable> neighbors_available_connections;
  std::set<CONNECTIONS_COMBINATION> combinations;

  bool retrieved_acyclic_substitutes = false, retrieved_ring_substitutes = false;
  bool has_substitutes = false, has_ring_substitutes = false, has_acyclic_substitutes = false;
  bool rings_can_be_added = false, rings_can_be_removed = false, rings_can_be_kept_constant = false;

  std::map<CONNECTIONS_COMBINATION, QUERY_RESULT> acyclic_substitutes, ring_substitutes;
  INVENTORY_RESULT inventory_acyclic_substitutes, inventory_ring_substitutes;

  int max_vtx_id = 0;
  std::map<int, const MolBrick*> vtx_ids_neighbors;
  std::map<const MolBrick*, int> neighbors_vtx_ids;
  std::map<int, const Connection*> vtx_ids_connections;
  std::vector<std::pair<int, int>> edges;
  BipartiteGraph graph;

  bool calculated_weights = false;
  float weight = 0;
  std::vector<float> combination_weights;

  bool guided = false, evaluated = false, hopcroft_ready = false;

public:
  SubstitutionPoint(MolBrick* substitutable, bool guided = false);

  void PerceiveNeighborsAvailableConnections(std::mt19937& prng);
  void GenerateCombinations(std::mt19937& prng);

  void PrepareHopcroftKarp(std::mt19937& prng);
  PATH HopcroftKarp(const MolBrick& substitute, ConnectionCompatibilities& compatibilities, std::mt19937& prng);

  void EvaluateSubstitutability(ConnectionQueryResults& query_results, std::mt19937& prng);
  void EvaluateSubstitutability(ReconstructionsInventory& inventory, ConnectionCompatibilities& compatibilities, std::mt19937& prng);

  void RetrieveAcyclicSubstitutes(ConnectionQueryResults& query_results, std::mt19937& prng);
  void RetrieveRingSubstitutes(ConnectionQueryResults& query_results, std::mt19937& prng);
  void RetrieveAcyclicSubstitutes(ReconstructionsInventory& inventory, ConnectionCompatibilities& compatibilities, std::mt19937& prng);
  void RetrieveRingSubstitutes(ReconstructionsInventory& inventory, ConnectionCompatibilities& compatibilities, std::mt19937& prng);

  void CalcAcyclicBasedWeight();
  void CalcRingBasedWeight();
  void CalcInventoryAcyclicBasedWeight();
  void CalcInventoryRingBasedWeight();

  std::tuple<CONNECTIONS_COMBINATION, unsigned, float> SampleAcyclicSubstitute(std::mt19937& prng);
  std::tuple<CONNECTIONS_COMBINATION, unsigned, float> SampleRingSubstitute(std::mt19937& prng);

  friend class ReconstructedMol;
};


class InsertionPoint {
  MolBrick* owner;
  std::unordered_map<MolBrick*, ConnectionsTable> neighbors_available_connections;
  std::map<MOLBRICKS_COMBINATION, std::set<CONNECTIONS_COMBINATION>> combinations;
  unsigned n_neighbors, n_neighbors_combinations, n_connections_combinations = 0;
  bool guided = false;

  bool retrieved_acyclic_inserts = false, retrieved_ring_inserts = false;
  bool has_acyclic_inserts = false, has_ring_inserts = false;
  std::map<CONNECTIONS_COMBINATION, QUERY_RESULT> acyclic_inserts, ring_inserts;

  bool calculated_weights = false;
  float weight = 0;
  std::vector<float> combination_weights;
  std::map<std::pair<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION>, std::vector<float>> insert_weights;

  std::map<MOLBRICKS_COMBINATION, std::map<const MolBrick*, const Connection*>> matchings;
  unsigned n_matchings = 0;

public:
  InsertionPoint(MolBrick* brick, std::mt19937& prng, bool guided = false);

  void RetrieveAcyclicInserts(ConnectionQueryResults& query_results);
  void RetrieveRingInserts(ConnectionQueryResults& query_results);

  void CalcAcyclicBasedWeight();
  void CalcRingBasedWeight();

  std::tuple<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION, unsigned, float> SampleAcyclicInsert(std::mt19937& prng);
  std::tuple<MOLBRICKS_COMBINATION, CONNECTIONS_COMBINATION, unsigned, float> SampleRingInsert(std::mt19937& prng);

  void GenerateMatchingsForBrick(const MolBrick* brick, ConnectionCompatibilities& compatibilities, std::mt19937& prng);

  friend class MolBrick;
  friend class ReconstructedMol;
};


// During a genetic operation the ReconstructionSchematic (i.e. the chromosomal
// meta-graph) is altered, creating and/or destroying links between MolBricks
// (i.e. genes). A GeneticLogEntry keeps track of one of those changes.
class GeneticLogEntry {
public:
  enum class Type {CREATION, DELETION};
private:
  Type type; // Type of meta-graph edge modification
  Connection connection; // Involved Connection
  unsigned atom_idx; // Atom index of the involved ConnectionPoint
  unsigned pseudofragment_id; // Database ID of the inserted or deleted MolBrick

public:
  GeneticLogEntry(GeneticLogEntry::Type type, const Connection& connection, unsigned atom_idx, unsigned pseudofragment_id);
  GeneticLogEntry::Type GetType() const;
  const Connection& GetConnection() const;
  unsigned GetAtomIdx() const;
  unsigned GetPseudofragmentID() const;
};


class GeneticLog {
private:
  // The GeneticLog stores the IDs of up to 2 Pseudofragments that had special
  // relevance during the logged operation. For efficiency reasons each of these
  // Pseudofragments is allocated its own buffer.
  unsigned pseudofragment_id1 = 0u, pseudofragment_id2 = 0u;
  // All modifications that took place during the logged operation, stratified
  // according to the involved Pseudofragment ID. This allows for more efficient
  // usage of the file system by the EvolutionGuide.
  std::map<unsigned, std::vector<GeneticLogEntry>> entries;

public:
  GeneticLog();
  bool SetImportantPseudofragmentID(unsigned pseudofragment_id);
  void AddEntry(GeneticLogEntry::Type type, const Connection& connection, unsigned atom_idx, unsigned pseudofragment_id);
  void Clear();
  unsigned GetImportantPseudofragmentID1() const;
  unsigned GetImportantPseudofragmentID2() const;
  const std::map<unsigned, std::vector<GeneticLogEntry>>& GetEntries() const;
};


class ReconstructedMol {
  unsigned id = 0;
  RDKit::RWMol pseudomol;
  RDKit::RWMol sanitized_mol;

  ConnectionsTable connections;
  ConnectionsTable internal_connections;
  std::unordered_map<unsigned, MolBrick> bricks;
  ReconstructionSchematic schematic;

  unsigned level = 0, n_ring_atoms = 0, max_atom_idx = 0, max_schematic_idx = 0;
  std::vector<unsigned> available_schematic_idxs;
  std::vector<unsigned> chiral_center_atom_indices, stereo_bond_indices;

  bool sanitized_mol_updated = false, smiles_updated = false, sanitized_smiles_updated = false, fingerprint_updated = false, stereo_updated = false;
  std::string smiles, sanitized_smiles;
  RDKit::SparseIntVect<std::uint32_t> fingerprint;
  float score = 0, previous_score = 0;
  bool is_child = false;
  unsigned parent_id = 0;
  GeneticLog genetic_log;

public:
  ReconstructedMol();
  ReconstructedMol(const MolBrick& brick);

  bool PeripheralExpansion(sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, bool guided = false, bool update_stereo = false);
  bool InternalExpansion(sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, bool guided = false, bool update_stereo = false);
  bool PeripheralDeletion(SizeController& controller, std::mt19937& prng, bool guided = false, bool update_stereo = false);
  bool InternalDeletion(ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, bool guided = false, bool update_stereo = false);
  bool Substitution(const std::string& location, sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, bool guided = false, bool update_stereo = false);
  bool Crossover(const std::string& location, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, ReconstructionsInventory& inventory, bool guided = false, bool update_stereo = false);
  bool Translation(ConnectionQueryResults& query_results, std::mt19937& prng, bool guided = false, bool update_stereo = false);
  bool StereoFlip(std::mt19937& prng);
  bool Evolve(ReconstructionSettings& settings, sqlite3_stmt* select_statement, ConnectionQueryResults& query_results, SizeController& controller, std::mt19937& prng, ReconstructionsInventory& inventory, EvolutionReport& report, unsigned n_failures = 0);

  void SetID(unsigned new_id);
  void TakeBricksOwnership();
  std::pair<Connection, bool> HasKnownConnections(const ConnectionQueryResults& query_results) const;
  void AssignWeightsToConnections(ConnectionQueryResults& query_results);
  bool AllConnectionsHaveWeights() const;
  void ClearConnectionWeights();
  void SetScore(float new_score);
  void SetChildFlag(bool flag);
  void SetParentID(unsigned new_parent_id);

  void Sanitize();
  void AssignUnspecifiedStereochemistry(std::mt19937& prng);
  void GenerateSMILES();
  void GenerateSanitizedSMILES();
  void GenerateFingerprint();
  float Score(const RDKit::SparseIntVect<std::uint32_t>* reference);
  double GetSimilarity(ReconstructedMol& reconstruction);
  bool StereoIsUpdated() const;

  unsigned GetID() const;
  const RDKit::RWMol& GetPseudomol() const;
  const RDKit::RWMol& GetSanitizedMol();
  const ConnectionsTable& GetConnections() const;
  const ConnectionsTable& GetInternalConnections() const;
  const std::unordered_map<unsigned, MolBrick>& GetBricks() const;
  const ReconstructionSchematic& GetSchematic() const;
  unsigned GetNBricks() const;
  unsigned GetLevel() const;
  unsigned GetNRingAtoms() const;
  unsigned GetMaxAtomIdx() const;
  unsigned GetMaxSchematicIdx() const;
  float GetScore() const;
  const std::string& GetSMILES();
  const std::string& GetSanitizedSMILES();
  const RDKit::SparseIntVect<std::uint32_t>& GetFingerprint();
  bool IsChild() const;
  unsigned GetParentID() const;
  void Draw(const std::string& file_path);

  // Serialization protocol. Note that the protocol is implemented with the
  // boost::serialization library but built around the RDKit's MolPickler
  // serializer. Therefore the behaviour of the serializer and deserializer
  // are different (i.e. split).
  template <class Archive>
  void save(Archive& ar, const unsigned version = reconstruction_version) const {
    if (!molpickler_configured) {
      ConfigureMolpickler();
    };
    ar << id;
    std::string pickle;
    RDKit::MolPickler::pickleMol(pseudomol, pickle);
    ar << pickle;
    pickle.clear();
    RDKit::MolPickler::pickleMol(sanitized_mol, pickle);
    ar << pickle;
    ar << connections << internal_connections;
    ar << bricks;
    ar << schematic;
    ar << level << n_ring_atoms << max_atom_idx << max_schematic_idx;
    ar << available_schematic_idxs << chiral_center_atom_indices << stereo_bond_indices;
    ar << sanitized_mol_updated << smiles_updated << sanitized_smiles_updated << stereo_updated;
    ar << smiles << sanitized_smiles;
    ar << score << previous_score;
    ar << is_child;
    ar << parent_id;
  };

  template <class Archive>
  void load(Archive& ar, const unsigned version = reconstruction_version) {
    if (!molpickler_configured) {
      ConfigureMolpickler();
    };
    ar >> id;
    std::string pickle;
    ar >> pickle;
    RDKit::MolPickler::molFromPickle(pickle, pseudomol);
    pickle.clear();
    ar >> pickle;
    RDKit::MolPickler::molFromPickle(pickle, sanitized_mol);
    ar >> connections >> internal_connections;
    ar >> bricks;
    ar >> schematic;
    ar >> level >> n_ring_atoms >> max_atom_idx >> max_schematic_idx;
    ar >> available_schematic_idxs >> chiral_center_atom_indices >> stereo_bond_indices;
    ar >> sanitized_mol_updated >> smiles_updated >> sanitized_smiles_updated >> stereo_updated;
    ar >> smiles >> sanitized_smiles;
    ar >> score >> previous_score;
    ar >> is_child;
    ar >> parent_id;
    fingerprint_updated = false;
    GenerateFingerprint();
  };

private:
  void AddLabelledPseudofragment(const Pseudofragment& pseudofragment);

  int EvaluatePeripheralCompatibility(ConnectionQueryResults& query_results) const;
  int EvaluateBrickSet(const std::vector<MolBrick*>& bricks) const;

  std::pair<const Connection*, const ConnectionPoint*> ChoosePeripheralConnectionPoint(ConnectionQueryResults& query_results, bool has_ring, bool ring_part, std::mt19937& prng) const;
  std::pair<const Connection*, const ConnectionPoint*> ChoosePeripheralConnectionPoint(bool has_ring, bool ring_part, std::mt19937& prng) const;
  MolBrick* ChooseBrickFromSet(const std::vector<MolBrick*>& bricks, bool has_ring, bool ring_part, std::mt19937& prng) const;

  const MolBrick* GetAtomSourceBrick(unsigned immutable_atom_idx) const;
  std::vector<MolBrick*> GetPeripheralBricks();
  std::vector<MolBrick*> GetInternalBricks();

  template <class Archive>
  void serialize(Archive& ar, const unsigned version = reconstruction_version) {
    boost::serialization::split_member(ar, *this, reconstruction_version);
  };

  friend class MolBrick;
  friend class ReconstructionsInventory;
  friend class InsertionPoint;
  friend class SubstitutionPoint;
  friend class EvolutionGuide;
  friend class LEADD;
  friend class boost::serialization::access;
  friend ReconstructedMol ConvertToReconstructedMol(const RDKit::ROMol& mol, bool fragment_rings, bool names_as_scores, boost::format& formatter);
};

class EvolutionGuide {
  // 3 buffers to store similarity data in memory are defined. Buffers 1 and 2
  // are dedicated to the "important" Pseudofragments, as defined by the GeneticLog.
  // Buffer 3 is used for all the "unimportant" Pseudofragments, and is potentially
  // reloaded multiple times.
  float* buffer1 = nullptr;
  float* buffer2 = nullptr;
  float* buffer3 = nullptr;
  unsigned buffer1_id = 0, buffer2_id = 0, buffer3_id = 0;
  hsize_t hyperslab_size[2] = { 0, 0 };
  hsize_t hyperslab_offset[2] = { 0, 0 };
  hsize_t buffer_size[1] = {0};
  H5::DataSpace buffer_dataspace;
  H5::H5File simatrix_file;
  H5::DataSet simatrix_dataset;
  H5::DataSpace simatrix_dataspace;
  H5::FloatType float_type;
  float acyclic_threshold = 0.0f, ring_threshold = 0.0f;
  float acyclic_positive_reinforcement = 0.0f, acyclic_negative_reinforcement = 0.0f;
  float ring_positive_reinforcement = 0.0f, ring_negative_reinforcement = 0.0f;

public:
  EvolutionGuide();
  EvolutionGuide(sqlite3* database, const ReconstructionSettings& settings);
  void AdjustWeights(ConnectionPoint* cpoint, float* buffer, float sign);
  void AdjustWeights(std::list<ReconstructedMol>& reconstructions, const std::set<unsigned>& survivor_ids);
  void Cleanup();
};

template <class T>
void AdvanceCombinations(const std::vector<T>& elements, std::vector<T>& combination, size_t start, size_t end, size_t index, size_t k, std::vector<std::vector<T>>& combinations) {
  if (index == k) {
    combinations.push_back(combination);
    return;
  };
  for (size_t i = start; (i <= end) && (end - i + 1 >= k - index); ++i) {
    combination[index] = elements[i];
    AdvanceCombinations(elements, combination, i + 1, end, index + 1, k, combinations);
  };
};

template <class T>
void Combinations(const std::vector<T>& elements, size_t k, std::vector<std::vector<T>>& combinations) {
  std::vector<T> combination(k);
  AdvanceCombinations(elements, combination, 0, elements.size() - 1, 0, k, combinations);
};

RDKit::Atom* GetAtomWithImmutableIdx(RDKit::RWMol& pseudomol, unsigned immutable_atom_idx);
unsigned GetMutableAtomIdx(const RDKit::RWMol& pseudomol, unsigned immutable_atom_idx);
std::vector<unsigned> GetMutableAtomIndices(const RDKit::RWMol& pseudomol, const std::vector<unsigned>& immutable_atom_indices);

// Function to convert a RDKit::ROMol to a ReconstructedMol.
ReconstructedMol ConvertToReconstructedMol(const RDKit::ROMol& mol, bool fragment_rings, bool names_as_scores, boost::format& formatter);

#endif
