#include "LEADD.hpp"

LEADD::LEADD(const ReconstructionSettings& settings, const std::string& output_directory_path) :
  settings(settings) {
  // Open a connection to the SQLite3 database containing the fragments.
  sqlite3_open(settings.GetFragmentDatabaseFile().c_str(), &database);
  // Prepare a SQLite3 statement to retrieve Pseudofragments by ID.
  sqlite3_prepare_v2(database, "SELECT pickle, pickle_size FROM pseudofragments WHERE id = ?", -1, &select_fragment_by_id, NULL);

  // Deserialize the file containing precalculated ConnectionQueryResults.
  std::ifstream cqr_stream(settings.GetConnectionQueryResultsFile(), std::ifstream::binary);
  boost::archive::binary_iarchive cqr_archive (cqr_stream);
  cqr_archive >> query_results;
  cqr_stream.close();

  // If the SAScore filter or heuristic are enabled, deserialize the FeatureLibrary.
  if (settings.UsingSAScoreFilter() || settings.UsingSAScoreHeuristic()) {
    std::ifstream feature_library_stream(settings.GetFeatureLibraryFile(), std::ifstream::binary);
    boost::archive::binary_iarchive feature_library_archive(feature_library_stream);
    feature_library_archive >> feature_library;
    feature_library_stream.close();
  };

  // Initialize the SAScore heuristic with the appropiate parameters.
  if (settings.UsingSAScoreHeuristic()) {
    sascore_heuristic = SAScoreHeuristic(settings.GetSAScoreHeuristicMu(), settings.GetSAScoreHeuristicSigma());
  };

  // If guided evolution is enabled, create a EvolutionGuide with access to
  // the fragments similarity matrix.
  if (settings.UsingGuidedEvolution()) {
    guide = EvolutionGuide(database, settings);
  };

  // Initialize a SizeController to decide whether to decrease, increase or keep
  // constant the number of ring atoms in a molecule.
  controller = SizeController(settings.GetMaxNRingAtoms(), settings.GetNRingAtomsMean(), settings.GetNRingAtomsSTDEV(), settings.GetPeakNRingAtomsProbability());

  // Seed the Mersenne Twister pseudorandom number generator.
  if (settings.SeedPRNG()) {
    prng.seed(settings.GetPRNGSeed());
  } else {
    std::random_device rd;
    prng.seed(rd());
  };

  // Either initialize a random population or load a pre-generated one from a file.
  if (settings.Restart()) {
    LoadPopulation();
  } else {
    InitializeRandomPopulation();
  };

  // Initialize an inventory object to keep track of the MolBricks that make up
  // the current pool of ReconstructedMols.
  inventory = ReconstructionsInventory(population, settings.GetScoreGamma());

  // If the internal topological similarity scoring function is being used,
  // create the fingerprint of the molecule that the algorithm will try to replicate.
  if (!settings.UsingExternalScoringFunction()) {
    RDKit::ROMol* reference_mol = RDKit::SmilesToMol(settings.GetTemplateMolSMILES());
    RDKit::ROMol* reference_mol_noh = RDKit::MolOps::removeHs(*reference_mol);
    reference = RDKit::MorganFingerprints::getFingerprint(*reference_mol_noh, 2);
    delete reference_mol;
    delete reference_mol_noh;
  };

  // If requested, create an EvolutionReport to write some statistics to CSV files.
  if (settings.WriteReport()) {
    report = EvolutionReport(output_directory_path);
  };

  // Allocate memory for vectors that are reused frequently during evolution.
  parent_weights.reserve(population.size());
  parent_indices.reserve(settings.GetNChildrenPerGeneration());

  // Store the highest ID assigned to a ReconstructedMol.
  max_individual_id = population.back().GetID();
};

void LEADD::MakePopulationTakeBrickOwnership() {
  for (ReconstructedMol& reconstruction : population) {
    reconstruction.TakeBricksOwnership();
  };
};

void LEADD::AssignConnectorWeightsToPopulation() {
  for (ReconstructedMol& reconstruction : population) {
    reconstruction.AssignWeightsToConnections(query_results);
  };
};

void LEADD::SeedPopulation() {
  // Verify that the population is empty.
  assert(population.empty());
  // Prepare the general SQLite3 statement.
  sqlite3_stmt* select_fragment_info;
  sqlite3_prepare_v2(database, "SELECT id, level, frequency FROM pseudofragments WHERE has_ring = ? AND ring_part = ?", -1, &select_fragment_info, NULL);
  // Initialize the objects to store the Pseudofragment IDs and weights.
  PseudofragmentIDs ids;
  PseudofragmentWeights weights;
  // Extract the relevant settings from the ReconstructionSettings.
  float acyclic_freq_gamma = settings.GetAcyclicFrequencyGamma();
  float acyclic_level_gamma = settings.GetAcyclicLevelGamma();
  float ring_freq_gamma = settings.GetRingFrequencyGamma();
  float ring_level_gamma = settings.GetRingLevelGamma();
  // Select the acyclic Pseudofragments.
  float level = 0, freq = 0;
  sqlite3_bind_int(select_fragment_info, 1, 0);
  sqlite3_bind_int(select_fragment_info, 2, 0);
  int rc = sqlite3_step(select_fragment_info);
  while (rc == SQLITE_ROW) {
    ids.AddAcyclicID(sqlite3_column_int64(select_fragment_info, 0));
    level = static_cast<float>(sqlite3_column_int(select_fragment_info, 1));
    freq = static_cast<float>(sqlite3_column_int(select_fragment_info, 2));
    weights.AddAcyclicWeight(std::pow(freq, acyclic_freq_gamma) * std::pow(level, acyclic_level_gamma));
    rc = sqlite3_step(select_fragment_info);
  };
  sqlite3_clear_bindings(select_fragment_info);
  sqlite3_reset(select_fragment_info);
  // Select the Pseudofragments containing rings.
  sqlite3_bind_int(select_fragment_info, 1, 1);
  sqlite3_bind_int(select_fragment_info, 2, 0);
  rc = sqlite3_step(select_fragment_info);
  while (rc == SQLITE_ROW) {
    ids.AddRingID(sqlite3_column_int64(select_fragment_info, 0));
    level = static_cast<float>(sqlite3_column_int(select_fragment_info, 1));
    freq = static_cast<float>(sqlite3_column_int(select_fragment_info, 2));
    weights.AddRingWeight(std::pow(freq, ring_freq_gamma) * std::pow(level, ring_level_gamma));
    rc = sqlite3_step(select_fragment_info);
  };
  sqlite3_clear_bindings(select_fragment_info);
  sqlite3_reset(select_fragment_info);
  sqlite3_finalize(select_fragment_info);
  // Create the probability distributions for the retrieved Pseudofragment IDs
  // based on their weights.
  std::discrete_distribution<unsigned> acyclic_distribution(weights.GetAcyclicWeights().begin(), weights.GetAcyclicWeights().end());
  std::discrete_distribution<unsigned> ring_distribution(weights.GetRingWeights().begin(), weights.GetRingWeights().end());
  // Sample a ReconstuctedMol seed.
  unsigned n_seeds = settings.GetNSeeds();
  unsigned idx = 0, pseudofragment_id = 0;
  float weight = 0;
  for (unsigned seed_id = 1; seed_id <= n_seeds; ++seed_id) {
    bool ring_seed = controller.DecideChange(0, prng);
    if (ring_seed) {
      idx = ring_distribution(prng);
      pseudofragment_id = ids.GetRingIDs()[idx];
      weight = weights.GetRingWeights()[idx];
    } else {
      idx = acyclic_distribution(prng);
      pseudofragment_id = ids.GetAcyclicIDs()[idx];
      weight = weights.GetAcyclicWeights()[idx];
    };
    Pseudofragment pseudofragment = GetPseudofragmentByIDFromDB(pseudofragment_id, select_fragment_by_id);
    MolBrick brick(pseudofragment, pseudofragment_id, 0, 0, weight, nullptr);
    ReconstructedMol reconstruction(brick);
    reconstruction.SetID(seed_id);
    population.push_back(std::move(reconstruction));
  };
  // Allow the population's ReconstructedMols to take ownership of their MolBricks.
  MakePopulationTakeBrickOwnership();
  // If necessary, assign new IDs and weights arrays to the ReconstructedMol's
  // ConnectionPoints.
  if (settings.UsingGuidedEvolution()) {
    AssignConnectorWeightsToPopulation();
  };
};

void LEADD::GrowIndividuals() {
  for (ReconstructedMol& reconstruction : population) {
    while (reconstruction.GetLevel() < settings.GetMinSeedSize()) {
      if (!reconstruction.GetConnections().Empty()) {
        reconstruction.PeripheralExpansion(select_fragment_by_id, query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
      } else {
        break;
      };
    };
  };
};

void LEADD::InitializeRandomPopulation() {
  // Erase the current population.
  population.clear();
  // Seed new individuals with a single random MolBrick.
  SeedPopulation();
  // Grow each individual up to the size specified by the ReconstructionSettings.
  GrowIndividuals();
};

void LEADD::SetPopulation(const std::list<ReconstructedMol>& new_population, bool reset_weights) {
  boost::format formatter("(%d, %d, %d)");
  // When a population comes from outside sources we can't be certain that the
  // Connections within the ReconstructedMols are also recorded in the ConnectionQueryResults.
  // The ReconstructedMols with unknown Connections are skipped, since our operators
  // can't handle them.
  unsigned id = 1;
  std::list<ReconstructedMol> valid_population;
  for (const ReconstructedMol& individual : new_population) {
    std::pair<Connection, bool> connection_known = individual.HasKnownConnections(query_results);
    if (connection_known.second) {
      if (id > max_individual_id) {
        max_individual_id = id;
      };
      valid_population.push_back(individual);
      valid_population.back().SetID(id++);
    } else {
      connection_known.first.GenerateString(formatter);
      std::cout << "WARNING: " << RDKit::MolToSmiles(individual.GetPseudomol()) << " contains an unknown Connection " << connection_known.first.GetString() << " and won't be added to the population." << std::endl;
    };
  };
  if (valid_population.empty()) {
    throw std::runtime_error("The valid population is empty");
  };
  population = valid_population;
  MakePopulationTakeBrickOwnership(); // This must be done BEFORE creating the inventory!
  inventory = ReconstructionsInventory(population, settings.GetScoreGamma());
  if (settings.UsingGuidedEvolution() && reset_weights) {
    AssignConnectorWeightsToPopulation();
  };
};

void LEADD::SetPopulationFromSMILES(const std::vector<std::string>& new_population_smiles) {
  std::list<ReconstructedMol> new_population = MakePopulationFromSMILES(new_population_smiles);
  SetPopulation(new_population, true);
};

void LEADD::LoadPopulation() {
  // Load the population contained in the serialized file.
  std::list<ReconstructedMol> new_population;
  std::ifstream input_stream(settings.GetRestartInputFile(), std::ifstream::binary);
  boost::archive::binary_iarchive archive(input_stream);
  archive >> new_population;
  input_stream.close();
  // Replace the current population with the new one.
  SetPopulation(new_population, settings.ResetWeightsOnRestart());
};

void LEADD::SavePopulation() const {
  std::ofstream output_stream(settings.GetRestartOutputFile(), std::ofstream::binary);
  boost::archive::binary_oarchive archive(output_stream);
  archive << population;
  output_stream.close();
};

void LEADD::WritePopulationSMILES(const std::string& file_path) {
  std::ofstream output_stream(file_path);
  for (ReconstructedMol& individual : population) {
    output_stream << individual.GetSanitizedSMILES() << " " << individual.GetScore() << "\n";
  };
  output_stream.close();
};

void LEADD::SampleParents() {
  // Clear the previously stored parent weights and indices.
  parent_weights.clear();
  parent_indices.clear();
  // Calculate the sampling weights of the population according to their scores.
  float cumulative_score = 0.0f;
  for (const ReconstructedMol& individual : population) {
    float score = individual.GetScore();
    cumulative_score += score;
    parent_weights.push_back(std::pow(score, settings.GetScoreGamma()));
  };
  // Sample some individuals to be parents for the next generation.
  // If all individuals have a weight of 0 the choice is purely random.
  if (cumulative_score == 0.0f) {
    std::uniform_int_distribution<unsigned> distribution (0, settings.GetNChildrenPerGeneration() - 1);
    for (unsigned child_n = 0; child_n < settings.GetNChildrenPerGeneration(); ++child_n) {
      parent_indices.push_back(distribution(prng));
    };
  // Otherwise the selection is fitness-proportionate.
  } else {
    std::discrete_distribution<unsigned> distribution(parent_weights.begin(), parent_weights.end());
    for (unsigned child_n = 0; child_n < settings.GetNChildrenPerGeneration(); ++child_n) {
      parent_indices.push_back(distribution(prng));
    };
  };
};

unsigned LEADD::GenerateChildren() {
  // Sample some individuals to be the parents of the upcoming children.
  SampleParents();

  // Iterate over the previously sampled parent ReconstructedMols.
  unsigned n_generated_children = 0;
  std::list<ReconstructedMol>::const_iterator it, begin_it = population.begin();
  for (unsigned parent_idx : parent_indices) {
    it = begin_it;
    std::advance(it, parent_idx);
    const ReconstructedMol& parent = *it;
    // Create a copy of the parent which will become the child ReconstructedMol.
    ReconstructedMol child = parent;
    child.TakeBricksOwnership();
    // Attempt to evolve the child.
    bool evolved = child.Evolve(settings, select_fragment_by_id, query_results, controller, prng, inventory, report);
    if (!evolved) {
      continue;
    };
    // Check if the child fulfills the criteria to be included in the population.
    // The child must be sufficiently distinct from other members of the population
    // (including both parents and children).
    bool too_similar = false;
    if (settings.UsingSimilarityFilter()) {
      for (ReconstructedMol& individual : population) {
        float similarity = child.GetSimilarity(individual);
        if (similarity > settings.GetMaxChildSimilarity()) {
          too_similar = true;
          break;
        };
      };
    };
    if (too_similar) {
      continue;
    };
    // The child must be deemed to be synthetically feasible.
    if (settings.UsingSAScoreFilter()) {
      double sascore = 0.0;
      if (child.StereoIsUpdated()) {
        sascore = SAScore(child.GetSanitizedMol(), child.GetFingerprint(), feature_library, false);
      } else {
        sascore = SAScore(child.GetSanitizedMol(), child.GetFingerprint(), feature_library, true);
      };
      if (sascore > settings.GetMaxSAScore()) {
        continue;
      };
    };
    // If the criteria are fulfilled, add it to the population and label it as
    // a child.
    child.SetChildFlag(true);
    child.SetID(++max_individual_id);
    child.SetParentID(parent.GetID());
    population.push_back(std::move(child));
    population.back().TakeBricksOwnership();
    ++n_generated_children;
    // Record the child's constituent MolBricks in the ReconstructionsInventory.
    inventory.Add(population.back());
  };
  return n_generated_children;
};

void LEADD::FinishGeneration() {
  // Increment the generation counter.
  ++generation_n;

  // If required, write the statistics of this generation to the report..
  if (settings.WriteReport()) {
    UpdateReport();
  };

  // Determine the score of the present population.
  float score = 0;
  if (settings.AverageScoreAsTerminationCriterion()) {
    for (const ReconstructedMol& individual : population) {
      score += individual.GetScore();
    };
    score /= population.size();
  } else {
    score = population.front().GetScore();
  };

  // If the score of the population increased during the last generation reset
  // the counter of unfruitful generations.
  if (score > best_score) {
    best_score = score;
    n_generations_stuck = 0;
    // If requested, draw the best molecule and write out its ConnectionPoint weights.
    if (settings.WriteReport() && settings.MonitorBestMolecule()) {
      ReportOnBestMolecule();
    };
  // Otherwise increment the counter of unfruitful generations.
  } else {
    ++n_generations_stuck;
  };

  // If requested, save the serialized population to a file.
  if ((settings.GetNGenerationsPerSave() > 0) && settings.SaveProgress()) {
    if (generation_n % settings.GetNGenerationsPerSave() == 0) {
      SavePopulation();
    };
  };
};

unsigned LEADD::SelectivePressure() {
  // If required, modify the individuals scores using the SAScore heuristic.
  if (settings.UsingSAScoreHeuristic()) {
    for (ReconstructedMol& individual : population) {
      double sascore = 0.0;
      if (individual.StereoIsUpdated()) {
        sascore = SAScore(individual.GetSanitizedMol(), individual.GetFingerprint(), feature_library, false);
      } else {
        sascore = SAScore(individual.GetSanitizedMol(), individual.GetFingerprint(), feature_library, true);
      };
      individual.score = sascore_heuristic(individual.score, sascore);
    };
  };

  // Sort the population according to their scores in descending order,
  // and record which individuals will survive this generation.
  population.sort(
    [](const ReconstructedMol& r1, const ReconstructedMol& r2) {
      return r2.GetScore() < r1.GetScore();
    });
  std::list<ReconstructedMol>::const_iterator it = population.begin(), end_it = population.end();
  unsigned n_survivors = 0, max_n_survivors = settings.GetNSurvivorsPerGeneration();
  std::set<unsigned> survivor_ids;
  while (it != end_it && n_survivors < max_n_survivors) {
    survivor_ids.insert(it->GetID());
    ++n_survivors;
    ++it;
  };

  // Adjust the ConnectionPoint weights of the surviving individuals.
  // This may include both children and parents that reproduced. For the parents
  // the child's GeneticLog is used.
  if (settings.UsingGuidedEvolution()) {
    guide.AdjustWeights(population, survivor_ids);
  };

  // Retain the individuals with the highest scores.
  while (it != population.end()) {
    inventory.Remove(*it);
    population.erase(it++);
  };

  // De-flag all individuals as children.
  for (ReconstructedMol& individual : population) {
    individual.SetChildFlag(false);
  };

  // Wrap up this generation.
  FinishGeneration();
  return n_survivors;
};

bool LEADD::TerminationCriteriaMet() {
  if (generation_n >= settings.GetNGenerations()) {
    return true;
  };
  if (n_generations_stuck >= settings.GetMaxGenerationsStuck()) {
    return true;
  };
  if (best_score >= settings.GetTerminationScore()) {
    return true;
  };
  return false;
};

void LEADD::ScoreChildrenBySimilarity() {
  assert(reference != nullptr);
  for (ReconstructedMol& individual : population) {
    if (individual.IsChild()) {
      individual.Score(reference);
    };
  };
};

unsigned LEADD::ScoreChildrenExternally() {
  // Generate the SMILES of the sanitized ReconstructedMols and write them to
  // a file alongside the ReconstructedMol's ID.
  std::ofstream output_stream (settings.GetScoringFunctionInputFile());
  for (ReconstructedMol& individual : population) {
    // Only those ReconstructedMols flagged as children are scored.
    if (individual.IsChild()) {
      assert(individual.GetID() != 0u);
      output_stream << individual.GetSanitizedSMILES() << " " << individual.GetID() << "\n";
    };
  };
  output_stream.close();
  // Call the external scoring function using the aforementioned file as input.
  int result = std::system(settings.GetScoringFunctionCall().c_str());
  assert(result == 0);
  std::filesystem::path scores_path(settings.GetScoresOutputFile());
  if (!std::filesystem::exists(scores_path)) {
    throw std::runtime_error("Scores output file couldn't be read.");
  };
  // Read the ReconstructedMol IDs and their calculated scores from the scoring
  // function's output file.
  std::ifstream input_stream(settings.GetScoresOutputFile());
  std::string line, id, score;
  std::map<unsigned, float> ids_scores;
  while (std::getline(input_stream, line)) {
    std::stringstream line_stream(line);
    line_stream >> id >> score;
    ids_scores.insert({std::stoi(id) , std::stof(score)});
  };
  input_stream.close();
  // Update the ReconstructedMols' scores.
  unsigned n_scored_molecules = 0u;
  std::map<unsigned, float>::const_iterator it, end_it = ids_scores.end();
  for (ReconstructedMol& individual : population) {
    if (individual.IsChild()) {
      // If the ReconstructedMol was successfully scored, assign it its new score.
      it = ids_scores.find(individual.GetID());
      if (it != end_it) {
        individual.SetScore(it->second);
        ++n_scored_molecules;
      // If it failed to be scored, set the score to 0.
      } else {
        individual.SetScore(0.0f);
      };
    };
  };
  // Return the number of molecules that were scored.
  return n_scored_molecules;
};

void LEADD::ScoreChildren() {
  // If an external scoring function was specified call it.
  if (settings.UsingExternalScoringFunction()) {
    ScoreChildrenExternally();
  // Otherwise use the similarity to the template molecule as a score.
  } else {
    ScoreChildrenBySimilarity();
  };
};

void LEADD::UpdateReport() {
  report.WriteScores(generation_n, population);
  report.WriteNRings(generation_n, population);
};

void LEADD::WriteOperationFrequenciesToReport() {
  report.WriteOperationFrequencies();
};

void LEADD::ReportOnBestMolecule() {
  ReconstructedMol& rank1 = population.front();
  std::string score_text = std::to_string(rank1.GetScore());
  rank1.Draw("rank1_" + score_text + ".svg");
  report.WriteWeights(rank1, "rank1_weights_" + score_text + ".csv");
};

void LEADD::Cleanup() {
  // Close the EvolutionReport.
  if (settings.WriteReport()) {
    report.Close();
  };

  // Close the connection to the database.
  sqlite3_close(database);

  // Close the connection to the similarity matrix.
  if (settings.UsingGuidedEvolution()) {
    guide.Cleanup();
  };
};

unsigned LEADD::GetGenerationNumber() const {
  return generation_n;
};

float LEADD::GetBestScore() const {
  return best_score;
};

ReconstructedMol& LEADD::GetBestIndividual() {
  return population.front();
};

std::list<ReconstructedMol>& LEADD::GetPopulation() {
  return population;
};

const ReconstructionSettings& LEADD::GetSettings() const {
  return settings;
};

EvolutionReport& LEADD::GetReport() {
  return report;
};

std::list<ReconstructedMol> MakePopulationFromSMILES(const std::vector<std::string>& population_smiles, bool fragment_rings) {
  boost::format formatter("(%d, %d, %d)");
  std::list<ReconstructedMol> population;
  for (const std::string& smiles : population_smiles) {
    RDKit::ROMol* molecule = RDKit::SmilesToMol(smiles);
    population.push_back(ConvertToReconstructedMol(*molecule, fragment_rings, false, formatter));
    delete molecule;
  };
  return population;
};
