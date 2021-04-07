#include "ReconstructionSettings.hpp"

// Class ReconstructionSettings
ReconstructionSettings::ReconstructionSettings() {
  GenerateOperationDistributions();
};

ReconstructionSettings::ReconstructionSettings(const std::string& settings_file, bool check_paths) {
  // Iterate over the lines in the input file.
  std::ifstream input_stream(settings_file);
  std::string line, key, value;
  while (std::getline(input_stream, line)) {
    if (line[0] == '#') {
      continue;
    };
    // Parse the line into key and value.
    std::stringstream line_stream (line);
    line_stream >> key;
    line_stream.get();
    std::getline(line_stream, value);
    // If the value is blank skip the line.
    if (StringIsBlank(value)) {
      continue;
    };
    // Assign the value to the right variable according to the key.
    if (key == "FRAGMENT_DATABASE_FILE") {
      fragment_database_file = value;
      if (check_paths){
        if (!FileIsAccessible(fragment_database_file)) {
          throw std::runtime_error("ERROR: FRAGMENT_DATABASE_FILE doesn't reference an existing file.");
        };
      };
    } else if (key == "CONNECTION_QUERY_RESULTS_FILE") {
      connection_query_results_file = value;
      if (check_paths) {
        if (!FileIsAccessible(connection_query_results_file)) {
          throw std::runtime_error("ERROR: CONNECTION_QUERY_RESULTS_FILE doesn't reference an existing file.");
        };
      };

    } else if (key == "ACYCLIC_FREQUENCY_GAMMA") {
      gammas["acyclic_frequency"] = std::stod(value);
    } else if (key == "ACYCLIC_LEVEL_GAMMA") {
      gammas["acyclic_level"] = std::stod(value);
    } else if (key == "RING_FREQUENCY_GAMMA") {
      gammas["ring_frequency"] = std::stod(value);
    } else if (key == "RING_LEVEL_GAMMA") {
      gammas["ring_level"] = std::stod(value);
    } else if (key == "SCORE_GAMMA") {
      gammas["score"] = std::stod(value);

    } else if (key == "PERIPHERAL_EXPANSION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["peripheral_expansion"] = probability;
    } else if (key == "PERIPHERAL_DELETION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["peripheral_deletion"] = probability;
    } else if (key == "PERIPHERAL_SUBSTITUTION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["peripheral_substitution"] = probability;
    } else if (key == "PERIPHERAL_CROSSOVER_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["peripheral_crossover"] = probability;
    } else if (key == "INTERNAL_EXPANSION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["internal_expansion"] = probability;
    } else if (key == "INTERNAL_DELETION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["internal_deletion"] = probability;
    } else if (key == "INTERNAL_SUBSTITUTION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["internal_substitution"] = probability;
    } else if (key == "INTERNAL_CROSSOVER_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["internal_crossover"] = probability;
    } else if (key == "TRANSLATION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["translation"] = probability;
    } else if (key == "STEREO_FLIP_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["stereo_flip"] = probability;
    } else if (key == "RING_EXPANSION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["ring_expansion"] = probability;
    } else if (key == "RING_DELETION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["ring_deletion"] = probability;
    } else if (key == "RING_SUBSTITUTION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["ring_substitution"] = probability;
    } else if (key == "RING_CROSSOVER_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["ring_crossover"] = probability;
    } else if (key == "RING_TRANSLATION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["ring_translation"] = probability;
    } else if (key == "CYCLIZATION_PROBABILITY") {
      unsigned probability = std::stoi(value);
      assert(probability >= 0u);
      probabilities["cyclization"] = probability;

    } else if (key == "PRNG_SEED") {
      prng_seed = std::stoi(value);
      assert(prng_seed >= 0u);
      if (prng_seed > 0u) {
        seed_prng = true;
      } else {
        seed_prng = false;
      };

    } else if (key == "N_RING_ATOMS_MEAN") {
      n_ring_atoms_mean = std::stod(value);
      assert(n_ring_atoms_mean >= 0.0f);
    } else if (key == "N_RING_ATOMS_STDEV") {
      n_ring_atoms_stdev = std::stod(value);
      assert(n_ring_atoms_stdev >= 0.0f);
    } else if (key == "MAX_N_RING_ATOMS") {
      max_n_ring_atoms = std::stoi(value);
      assert(max_n_ring_atoms >= 0u);
    } else if (key == "PEAK_N_RING_ATOMS_PROBABILITY") {
      peak_n_ring_atoms_probability = std::stod(value);
      assert(peak_n_ring_atoms_probability >= 0.0f);

    } else if (key == "RING_SIZE_MEAN") {
      ring_size_mean = std::stod(value);
      assert(ring_size_mean >= 0.0f);
    } else if (key == "RING_SIZE_STDEV") {
      ring_size_stdev = std::stod(value);
      assert(ring_size_stdev >= 0.0f);
    } else if (key == "MAX_RING_SIZE") {
      max_ring_size = std::stoi(value);
      assert(max_ring_size >= 0u);
    } else if (key == "PEAK_RING_SIZE_PROBABILITY") {
      peak_ring_size_probability = std::stod(value);
      assert(peak_ring_size_probability >= 0.0f);

    } else if (key == "MIN_SEED_SIZE") {
      min_seed_size = std::stoi(value);
      assert(min_seed_size >= 0u);
    } else if (key == "ASSIGN_UNSPECIFIED_STEREO") {
      assign_unspecified_stereo = std::stoi(value);

    } else if (key == "N_SEEDS") {
      n_seeds = std::stoi(value);
      assert(n_seeds >= 0u);
    } else if (key == "N_GENERATIONS") {
      n_generations = std::stoi(value);
      assert(n_generations >= 0u);
    } else if (key == "N_CHILDREN_PER_GENERATION") {
      n_children_per_generation = std::stoi(value);
      assert(n_children_per_generation >= 0u);
    } else if (key == "N_SURVIVORS_PER_GENERATION") {
      n_survivors_per_generation = std::stoi(value);
      assert(n_survivors_per_generation >= 0u);
    } else if (key == "MAX_CHILD_SIMILARITY") {
      max_child_similarity = std::stod(value);
      assert(max_child_similarity >= 0.0f && max_child_similarity <= 1.0f);
      if (max_child_similarity > 0.0f) {
        similarity_filter = true;
      } else {
        similarity_filter = false;
      };

    } else if (key == "TERMINATION_SCORE") {
      termination_score = std::stod(value);
      assert(termination_score >= 0.0f);
    } else if (key == "AVERAGE_SCORE_AS_TERMINATION_CRITERION") {
      average_score_as_termination_criterion = std::stoi(value);
      assert(max_attempts_per_generation >= 0u);
    } else if (key == "MAX_GENERATIONS_STUCK") {
      max_generations_stuck = std::stoi(value);
      assert(max_generations_stuck >= 0u);
    } else if (key == "MAX_ATTEMPTS_PER_GENERATION") {
      max_attempts_per_generation = std::stoi(value);
      assert(max_attempts_per_generation >= 0u);
    } else if (key == "MAX_MULTIPLE_INTERSECTIONS") {
      max_multiple_intersections = std::stoi(value);
      assert(max_multiple_intersections >= 0u);
    } else if (key == "MAX_MAXIMUM_BIPARTITE_MATCHINGS") {
      max_maximum_bipartite_matchings = std::stoi(value);
      assert(max_maximum_bipartite_matchings >= 0u);

    } else if (key == "FEATURE_LIBRARY_FILE") {
      feature_library_file = value;
      if (check_paths) {
        if (!FileIsAccessible(feature_library_file)) {
          throw std::runtime_error("ERROR: FEATURE_LIBRARY_FILE doesn't reference an existing file.");
        };
      };
    } else if (key == "MAX_SASCORE") {
      max_sascore = std::stod(value);
      assert(max_sascore >= 0.0 && max_sascore <= 10.0);
    } else if (key == "SASCORE_HEURISTIC_MU") {
      sascore_heuristic_mu = std::stod(value);
    } else if (key == "SASCORE_HEURISTIC_SIGMA") {
      sascore_heuristic_sigma = std::stod(value);
    } else if (key == "USE_SASCORE_HEURISTIC") {
      sascore_heuristic = std::stoi(value);

    } else if (key == "SIMILARITY_MATRIX_FILE") {
      similarity_matrix_file = value;
      if (check_paths) {
        if (!FileIsAccessible(similarity_matrix_file)) {
          throw std::runtime_error("ERROR: SIMILARITY_MATRIX_FILE doesn't reference an existing file.");
        };
      };
    } else if (key == "ACYCLIC_LEARNING_SIMILARITY_THRESHOLD") {
      acyclic_learning_similarity_threshold = std::stod(value);
      assert(acyclic_learning_similarity_threshold >= 0.0f && acyclic_learning_similarity_threshold <= 1.0f);
    } else if (key == "RING_LEARNING_SIMILARITY_THRESHOLD") {
      ring_learning_similarity_threshold = std::stod(value);
      assert(ring_learning_similarity_threshold >= 0.0f && ring_learning_similarity_threshold <= 1.0f);
    } else if (key == "ACYCLIC_POSITIVE_REINFORCEMENT") {
      acyclic_positive_reinforcement = std::stod(value);
      assert(acyclic_positive_reinforcement >= 0.0f);
    } else if (key == "ACYCLIC_NEGATIVE_REINFORCEMENT") {
      acyclic_negative_reinforcement = std::stod(value);
      assert(acyclic_negative_reinforcement >= 0.0f && acyclic_negative_reinforcement <= 1.0f);
    } else if (key == "RING_POSITIVE_REINFORCEMENT") {
      ring_positive_reinforcement = std::stod(value);
      assert(ring_positive_reinforcement >= 0.0f);
    } else if (key == "RING_NEGATIVE_REINFORCEMENT") {
      ring_negative_reinforcement = std::stod(value);
      assert(ring_negative_reinforcement >= 0.0f && ring_negative_reinforcement <= 1.0f);

    } else if (key == "RESTART_INPUT_FILE") {
      restart_input_file = value;
      restart = true;
      if (check_paths) {
        if (!FileIsAccessible(restart_input_file)) {
          throw std::runtime_error("ERROR: RESTART_INPUT_FILE doesn't reference an existing file.");
        };
      };
    } else if (key == "RESET_WEIGHTS_ON_RESTART") {
      reset_weights_on_restart = std::stoi(value);
    } else if (key == "RESTART_OUTPUT_FILE") {
      restart_output_file = value;
      save_progress = true;
      if (check_paths) {
        if (!FileRootDirectoryIsAccessible(restart_output_file)) {
          throw std::runtime_error("ERROR: RESTART_OUTPUT_FILE isn't accessible.");
        };
      };
    } else if (key == "N_GENERATIONS_PER_SAVE") {
      n_generations_per_save = std::stoi(value);
      assert(n_generations_per_save >= 0u);

    } else if (key == "SCORING_FUNCTION_INPUT_FILE") {
      scoring_function_input_file = value;
      if (check_paths) {
        if (!FileRootDirectoryIsAccessible(scoring_function_input_file)) {
          throw std::runtime_error("ERROR: SCORING_FUNCTION_INPUT_FILE isn't accessible.");
        };
      };
    } else if (key == "SCORING_FUNCTION_CALL") {
      scoring_function_call = value;
    } else if (key == "SCORES_OUTPUT_FILE") {
      scores_output_file = value;
      if (check_paths) {
        if (!FileRootDirectoryIsAccessible(scores_output_file)) {
          throw std::runtime_error("ERROR: SCORES_OUTPUT_FILE isn't accessible.");
        };
      };
    } else if (key == "TEMPLATE_MOL_SMILES") {
      template_mol_smiles = value;
    } else if (key == "SCORE_FIRST_POPULATION") {
      score_first_population = std::stoi(value);

    } else if (key == "WRITE_REPORT") {
      write_report = std::stoi(value);
    } else if (key == "MONITOR_BEST_MOLECULE") {
      monitor_best_molecule = std::stoi(value);

    } else {
      std::stringstream ss;
      ss << "ERROR: Invalid key: '" << key << "'." << std::endl;
      throw std::runtime_error(ss.str());
    };
    value.clear();
  };
  input_stream.close();

  // Determine if a SAScore filter is being used.
  if (!feature_library_file.empty() && max_sascore < 10.0) {
    sascore_filter = true;
  } else {
    sascore_filter = false;
  };
  // If the SAScore heuristic is requested, verify that the feature library is specified.
  if (feature_library_file.empty() && sascore_heuristic) {
    throw std::runtime_error("ERROR: A FeatureLibrary must be provided to use the SAScore heuristic.");
  };

  // Determine if guided evolution is being used.
  if (!similarity_matrix_file.empty()) {
    if (acyclic_positive_reinforcement > 0.0f || acyclic_negative_reinforcement > 0.0f || ring_positive_reinforcement > 0.0f || ring_negative_reinforcement > 0.0f) {
      guided_evolution = true;
    };
  } else {
    guided_evolution = false;
  };

  // Determine if an external scoring function is being used.
  if (!scoring_function_input_file.empty() && !scoring_function_call.empty() && !scores_output_file.empty()) {
    external_scoring_function = true;
  } else {
    external_scoring_function = false;
  };

  // Create a probability distribution using the provided operation probabilities.
  GenerateOperationDistributions();
};

void ReconstructionSettings::GenerateOperationDistributions() {
  operations.clear();
  operation_weights.clear();
  for (const auto& p : probabilities) {
    operations.push_back(p.first);
    operation_weights.push_back(p.second);
  };
  operations_distribution = std::discrete_distribution<unsigned>(operation_weights.begin(), operation_weights.end());
};

const std::string& ReconstructionSettings::ChooseOperation(std::mt19937& prng) {
  return operations[operations_distribution(prng)];
};

void ReconstructionSettings::CheckRequiredInputFiles() const {
  if (fragment_database_file.empty()) {
    throw std::runtime_error("ERROR: FRAGMENT_DATABASE_FILE must be specified.");
  };
  if (connection_query_results_file.empty()) {
    throw std::runtime_error("ERROR: CONNECTION_QUERY_RESULTS_FILE must be specified.");
  };
};

void ReconstructionSettings::CheckScoringFunctionValidity() const {
  // If an external scoring function isn't being used, make sure that a template
  // molecule's SMILES was provided.
  if (!external_scoring_function) {
    if (template_mol_smiles.empty()) {
      throw std::runtime_error("ERROR: Either SCORING_FUNCTION_CALL or TEMPLATE_MOL_SMILES must be specified.");
    };
  };
};

const std::string& ReconstructionSettings::GetFragmentDatabaseFile() const {
  return fragment_database_file;
};

const std::string& ReconstructionSettings::GetConnectionQueryResultsFile() const {
  return connection_query_results_file;
};

float ReconstructionSettings::GetAcyclicFrequencyGamma() const {
  return gammas.at("acyclic_frequency");
};

float ReconstructionSettings::GetAcyclicLevelGamma() const {
  return gammas.at("acyclic_level");
};

float ReconstructionSettings::GetRingFrequencyGamma() const {
  return gammas.at("ring_frequency");
};

float ReconstructionSettings::GetRingLevelGamma() const {
  return gammas.at("ring_level");
};

float ReconstructionSettings::GetScoreGamma() const {
  return gammas.at("score");
};

unsigned ReconstructionSettings::GetOperationProbability(const std::string& operation) const {
  return probabilities.at(operation);
};

void ReconstructionSettings::SetPRNGSeed(unsigned new_prng_seed) {
  prng_seed = new_prng_seed;
};

unsigned ReconstructionSettings::GetPRNGSeed() const {
  return prng_seed;
};

bool ReconstructionSettings::SeedPRNG() const {
  return seed_prng;
};

void ReconstructionSettings::SetNRingAtomsMean(float new_n_ring_atoms_mean) {
  n_ring_atoms_mean = new_n_ring_atoms_mean;
};

void ReconstructionSettings::SetNRingAtomsSTDEV(float new_n_ring_atoms_stdev) {
  n_ring_atoms_stdev = new_n_ring_atoms_stdev;
};

void ReconstructionSettings::SetMaxNRingAtoms(unsigned new_max_n_ring_atoms) {
  max_n_ring_atoms = new_max_n_ring_atoms;
};

void ReconstructionSettings::SetPeakNRingAtomsProbability(float new_peak_n_ring_atoms_probability) {
  peak_n_ring_atoms_probability = new_peak_n_ring_atoms_probability;
};

float ReconstructionSettings::GetNRingAtomsMean() const {
  return n_ring_atoms_mean;
};

float ReconstructionSettings::GetNRingAtomsSTDEV() const {
  return n_ring_atoms_stdev;
};

unsigned ReconstructionSettings::GetMaxNRingAtoms() const {
  return max_n_ring_atoms;
};

float ReconstructionSettings::GetPeakNRingAtomsProbability() const {
  return peak_n_ring_atoms_probability;
};

float ReconstructionSettings::GetRingSizeMean() const {
  return ring_size_mean;
};

float ReconstructionSettings::GetRingSizeSTDEV() const {
  return ring_size_stdev;
};

unsigned ReconstructionSettings::GetMaxRingSize() const {
  return max_ring_size;
};

float ReconstructionSettings::GetPeakRingSizeProbability() const {
  return peak_ring_size_probability;
};

unsigned ReconstructionSettings::GetMinSeedSize() const {
  return min_seed_size;
};

bool ReconstructionSettings::AssignUnspecifiedStereo() const {
  return assign_unspecified_stereo;
};

void ReconstructionSettings::SetNSeeds(unsigned new_n_seeds) {
  n_seeds = new_n_seeds;
};

void ReconstructionSettings::SetNGenerations(unsigned new_n_generations) {
  n_generations = new_n_generations;
};

void ReconstructionSettings::SetNChildrenPerGeneration(unsigned new_n_children_per_generation) {
  n_children_per_generation = new_n_children_per_generation;
};

void ReconstructionSettings::SetNSurvivorsPerGeneration(unsigned new_n_survivors_per_generation) {
  n_survivors_per_generation = new_n_survivors_per_generation;
};

unsigned ReconstructionSettings::GetNSeeds() const {
  return n_seeds;
};

unsigned ReconstructionSettings::GetNGenerations() const {
  return n_generations;
};

unsigned ReconstructionSettings::GetNChildrenPerGeneration() const {
  return n_children_per_generation;
};

unsigned ReconstructionSettings::GetNSurvivorsPerGeneration() const {
  return n_survivors_per_generation;
};

float ReconstructionSettings::GetMaxChildSimilarity() const {
  return max_child_similarity;
};

bool ReconstructionSettings::UsingSimilarityFilter() const {
  return similarity_filter;
};

float ReconstructionSettings::GetTerminationScore() const {
  return termination_score;
};

void ReconstructionSettings::SetAverageScoreAsTerminationCriterion(bool new_average_score_as_termination_criterion) {
  average_score_as_termination_criterion = new_average_score_as_termination_criterion;
};

bool ReconstructionSettings::AverageScoreAsTerminationCriterion() const {
  return average_score_as_termination_criterion;
};

unsigned ReconstructionSettings::GetMaxGenerationsStuck() const {
  return max_generations_stuck;
};

unsigned ReconstructionSettings::GetMaxAttemptsPerGeneration() const {
  return max_attempts_per_generation;
};

unsigned ReconstructionSettings::GetMaxMultipleIntersections() const {
  return max_multiple_intersections;
};

unsigned ReconstructionSettings::GetMaxMaximumBipartiteMatchings() const {
  return max_maximum_bipartite_matchings;
};

const std::string& ReconstructionSettings::GetFeatureLibraryFile() const {
  return feature_library_file;
};

double ReconstructionSettings::GetMaxSAScore() const {
  return max_sascore;
};

double ReconstructionSettings::GetSAScoreHeuristicMu() const {
  return sascore_heuristic_mu;
};

double ReconstructionSettings::GetSAScoreHeuristicSigma() const {
  return sascore_heuristic_sigma;
};

bool ReconstructionSettings::UsingSAScoreFilter() const {
  return sascore_filter;
};

bool ReconstructionSettings::UsingSAScoreHeuristic() const {
  return sascore_heuristic;
};

const std::string& ReconstructionSettings::GetSimilarityMatrixFile() const {
  return similarity_matrix_file;
};

float ReconstructionSettings::GetAcyclicLearningSimilarityThreshold() const {
  return acyclic_learning_similarity_threshold;
};

float ReconstructionSettings::GetRingLearningSimilarityThreshold() const {
  return ring_learning_similarity_threshold;
};

float ReconstructionSettings::GetAcyclicPositiveReinforcement() const {
  return acyclic_positive_reinforcement;
};

float ReconstructionSettings::GetAcyclicNegativeReinforcement() const {
  return acyclic_negative_reinforcement;
};

float ReconstructionSettings::GetRingPositiveReinforcement() const {
  return ring_positive_reinforcement;
};

float ReconstructionSettings::GetRingNegativeReinforcement() const {
  return ring_negative_reinforcement;
};

void ReconstructionSettings::SetUsingGuidedEvolution(bool new_guided_evolution) {
  guided_evolution = new_guided_evolution;
};

bool ReconstructionSettings::UsingGuidedEvolution() const {
  return guided_evolution;
};

const std::string& ReconstructionSettings::GetRestartInputFile() const {
  return restart_input_file;
};

const std::string& ReconstructionSettings::GetRestartOutputFile() const {
  return restart_output_file;
};

unsigned ReconstructionSettings::GetNGenerationsPerSave() const {
  return n_generations_per_save;
};

bool ReconstructionSettings::SaveProgress() const {
  return save_progress;
};

bool ReconstructionSettings::Restart() const {
  return restart;
};

bool ReconstructionSettings::ResetWeightsOnRestart() const {
  return reset_weights_on_restart;
};

const std::string& ReconstructionSettings::GetScoringFunctionInputFile() const {
  return scoring_function_input_file;
};

const std::string& ReconstructionSettings::GetScoringFunctionCall() const {
  return scoring_function_call;
};

const std::string& ReconstructionSettings::GetScoresOutputFile() const {
  return scores_output_file;
};

const std::string& ReconstructionSettings::GetTemplateMolSMILES() const {
  return template_mol_smiles;
};

bool ReconstructionSettings::UsingExternalScoringFunction() const {
  return external_scoring_function;
};

bool ReconstructionSettings::ScoreFirstPopulation() const {
  return score_first_population;
};

bool ReconstructionSettings::WriteReport() const {
  return write_report;
};

bool ReconstructionSettings::MonitorBestMolecule() const {
  return monitor_best_molecule;
};

void ReconstructionSettings::Print() const {
  std::cout << "FRAGMENT_DATABASE_FILE: '" << fragment_database_file << "'" << std::endl;
  std::cout << "CONNECTION_QUERY_RESULTS_FILE: '" << connection_query_results_file << "'" << std::endl;
  std::cout << "ACYCLIC_FREQUENCY_GAMMA: " << gammas.at("acyclic_frequency") << std::endl;
  std::cout << "ACYCLIC_LEVEL_GAMMA: " << gammas.at("acyclic_level") << std::endl;
  std::cout << "RING_FREQUENCY_GAMMA: " << gammas.at("ring_frequency") << std::endl;
  std::cout << "RING_LEVEL_GAMMA: " << gammas.at("ring_level") << std::endl;
  std::cout << "SCORE_GAMMA: " << gammas.at("score") << std::endl;
  std::cout << "PERIPHERAL_EXPANSION_PROBABILITY: " << probabilities.at("peripheral_expansion") << std::endl;
  std::cout << "PERIPHERAL_DELETION_PROBABILITY: " << probabilities.at("peripheral_deletion") << std::endl;
  std::cout << "PERIPHERAL_SUBSTITUTION_PROBABILITY: " << probabilities.at("peripheral_substitution") << std::endl;
  std::cout << "PERIPHERAL_CROSSOVER_PROBABILITY: " << probabilities.at("peripheral_crossover") << std::endl;
  std::cout << "INTERNAL_EXPANSION_PROBABILITY: " << probabilities.at("internal_expansion") << std::endl;
  std::cout << "INTERNAL_DELETION_PROBABILITY: " << probabilities.at("internal_deletion") << std::endl;
  std::cout << "INTERNAL_SUBSTITUTION_PROBABILITY: " << probabilities.at("internal_substitution") << std::endl;
  std::cout << "INTERNAL_CROSSOVER_PROBABILITY: " << probabilities.at("internal_crossover") << std::endl;
  std::cout << "TRANSLATION_PROBABILITY: " << probabilities.at("translation") << std::endl;
  std::cout << "STEREO_FLIP_PROBABILITY: " << probabilities.at("stereo_flip") << std::endl;
  std::cout << "PRNG_SEED: " << prng_seed << std::endl;
  std::cout << "SEED_PRNG: " << seed_prng << std::endl;
  std::cout << "N_RING_ATOMS_MEAN: " << n_ring_atoms_mean << std::endl;
  std::cout << "N_RING_ATOMS_STDEV: " << n_ring_atoms_stdev << std::endl;
  std::cout << "MAX_N_RING_ATOMS: " << max_n_ring_atoms << std::endl;
  std::cout << "PEAK_N_RING_ATOMS_PROBABILITY: " << peak_n_ring_atoms_probability << std::endl;
  std::cout << "MIN_SEED_SIZE: " << min_seed_size << std::endl;
  std::cout << "ASSIGN_UNSPECIFIED_STEREO: " << assign_unspecified_stereo << std::endl;
  std::cout << "N_SEEDS: " << n_seeds << std::endl;
  std::cout << "N_GENERATIONS: " << n_generations << std::endl;
  std::cout << "N_CHILDREN_PER_GENERATION: " << n_children_per_generation << std::endl;
  std::cout << "N_SURVIVORS_PER_GENERATION: " << n_survivors_per_generation << std::endl;
  std::cout << "MAX_CHILD_SIMILARITY: " << max_child_similarity << std::endl;
  std::cout << "USING_SIMILARITY_FILTER: " << similarity_filter << std::endl;
  std::cout << "TERMINATION_SCORE: " << termination_score << std::endl;
  std::cout << "AVERAGE_SCORE_AS_TERMINATION_CRITERION: " << average_score_as_termination_criterion << std::endl;
  std::cout << "MAX_GENERATIONS_STUCK: " << max_generations_stuck << std::endl;
  std::cout << "MAX_ATTEMPTS_PER_GENERATION: " << max_attempts_per_generation << std::endl;
  std::cout << "FEATURE_LIBRARY_FILE: '" << feature_library_file << "'" << std::endl;
  std::cout << "MAX_SASCORE: " << max_sascore << std::endl;
  std::cout << "SASCORE_HEURISTIC_MU: " << sascore_heuristic_mu << std::endl;
  std::cout << "SASCORE_HEURISTIC_SIGMA: " << sascore_heuristic_sigma << std::endl;
  std::cout << "USING_SASCORE_FILTER: " << sascore_filter << std::endl;
  std::cout << "USING_SASCORE_HEURISTIC: " << sascore_heuristic << std::endl;
  std::cout << "SIMILARITY_MATRIX_FILE: '" << similarity_matrix_file << "'" << std::endl;
  std::cout << "ACYCLIC_LEARNING_SIMILARITY_THRESHOLD: " << acyclic_learning_similarity_threshold << std::endl;
  std::cout << "RING_LEARNING_SIMILARITY_THRESHOLD: " << ring_learning_similarity_threshold << std::endl;
  std::cout << "ACYCLIC_POSITIVE_REINFORCEMENT: " << acyclic_positive_reinforcement << std::endl;
  std::cout << "ACYCLIC_NEGATIVE_REINFORCEMENT: " << acyclic_negative_reinforcement << std::endl;
  std::cout << "RING_POSITIVE_REINFORCEMENT: " << ring_positive_reinforcement << std::endl;
  std::cout << "RING_NEGATIVE_REINFORCEMENT: " << ring_negative_reinforcement << std::endl;
  std::cout << "USING_GUIDED_EVOLUTION: " << guided_evolution << std::endl;
  std::cout << "RESTART_INPUT_FILE: '" << restart_input_file << "'" << std::endl;
  std::cout << "RESTART_OUTPUT_FILE: '" << restart_output_file << "'" << std::endl;
  std::cout << "N_GENERATIONS_PER_SAVE: " << n_generations_per_save << std::endl;
  std::cout << "SAVE_PROGRESS: " << save_progress << std::endl;
  std::cout << "RESTART: " << restart << std::endl;
  std::cout << "RESET_WEIGHTS_ON_RESTART: " << reset_weights_on_restart << std::endl;
  std::cout << "SCORING_FUNCTION_INPUT_FILE: '" << scoring_function_input_file << "'" << std::endl;
  std::cout << "SCORING_FUNCTION_CALL: '" << scoring_function_call << "'" << std::endl;
  std::cout << "SCORES_OUTPUT_FILE: '" << scores_output_file << "'" << std::endl;
  std::cout << "TEMPLATE_MOL_SMILES: '" << template_mol_smiles << "'" << std::endl;
  std::cout << "USING_EXTERNAL_SCORING_FUNCTION: " << external_scoring_function << std::endl;
  std::cout << "SCORE_FIRST_POPULATION: " << score_first_population << std::endl;
  std::cout << "WRITE_REPORT: " << write_report << std::endl;
  std::cout << "MONITOR_BEST_MOLECULE: " << monitor_best_molecule << std::endl;
};

bool StringIsBlank(const std::string& str) {
  return str.find_first_not_of(" \t\n\v\f\r") == str.npos;
};

bool FileIsAccessible(const std::string& path_str) {
  std::filesystem::file_status status = std::filesystem::status(path_str);
  if (std::filesystem::exists(status)) {
    return true;
  };
  return false;
};

bool FileRootDirectoryIsAccessible(const std::string& path_str) {
  std::filesystem::path path(path_str);
  if (std::filesystem::exists(path.parent_path())) {
    return true;
  };
  return false;
};
