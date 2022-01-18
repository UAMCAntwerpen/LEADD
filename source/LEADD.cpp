#include "LEADD.hpp"

LEADD::LEADD(const ReconstructionSettings& settings, const std::string& output_directory_path) :
  settings(settings),
  output_directory_path(output_directory_path),
  database(settings.GetFragmentDatabaseFile()) {
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
  unsigned max_id = 0;
  for (const ReconstructedMol& individual : population) {
    unsigned individual_id = individual.GetID();
    if (individual_id > max_id) {
      max_id = individual_id;
    };
  };
  max_individual_id = max_id;
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
  population.clear();
  // Create a ring Pseudofragment sampling distribution.
  std::map<unsigned, float> candidate_pseudofragments;
  for (const auto& [connection, query_results_by_connection_frequency] : query_results.GetStrictRingResults()) {
    for (const auto& [connection_frequency, query_result] : query_results_by_connection_frequency) {
      const auto& [pseudofragment_ids, pseudofragment_weights] = query_result;
      for (size_t i = 0; i < pseudofragment_ids.size(); ++i) {
        unsigned pseudofragment_id = pseudofragment_ids[i];
        float pseudofragment_weight = pseudofragment_weights[i];
        candidate_pseudofragments.insert({pseudofragment_id, pseudofragment_weight});
      };
    };
  };
  std::vector<float> candidate_pseudofragment_weights;
  candidate_pseudofragment_weights.reserve(candidate_pseudofragments.size());
  for (auto [pseudofragment_id, pseudofragment_weight] : candidate_pseudofragments) {
    candidate_pseudofragment_weights.push_back(pseudofragment_weight);
  };
  std::discrete_distribution<size_t> distribution (candidate_pseudofragment_weights.begin(), candidate_pseudofragment_weights.end());
  // Sample N random Pseudofragments and convert them to ReconstructedMols.
  for (unsigned seed_id = 1; seed_id <= settings.GetNSeeds(); ++seed_id) {
    // Choose a ring Pseudofragment.
    size_t pseudofragment_idx = distribution(prng);
    std::map<unsigned, float>::const_iterator pseudofragment_it = candidate_pseudofragments.begin();
    std::advance(pseudofragment_it, pseudofragment_idx);
    unsigned pseudofragment_id = pseudofragment_it->first;
    float pseudofragment_weight = pseudofragment_it->second;
    // Fetch the chosen Pseudofragment from the database.
    auto [pseudofragment, pseudofragment_frequency] = database.SelectPseudofragmentWithID(pseudofragment_id);
    // Convert it to a ReconstructedMol.
    MolBrick brick (pseudofragment, pseudofragment_id, 0, 0, pseudofragment_weight, nullptr);
    population.emplace_back(brick); // Construct the ReconstructedMol in place.
    ReconstructedMol& reconstruction = population.back();
    reconstruction.SetID(seed_id);
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
    while (reconstruction.GetSize() < settings.GetMinSeedSize()) {
      if (!reconstruction.GetConnections().Empty()) {
        reconstruction.PeripheralExpansion(database, query_results, controller, prng, settings.UsingGuidedEvolution(), settings.AssignUnspecifiedStereo());
      } else {
        break;
      };
    };
  };
};

void LEADD::InitializeRandomPopulation() {
  // Seed new individuals with a single random MolBrick.
  SeedPopulation();
  // Grow each individual up to the size specified by the ReconstructionSettings.
  GrowIndividuals();
};

void LEADD::SetPopulation(const std::list<ReconstructedMol>& new_population, bool reset_weights) {
  unsigned id = 1;
  std::list<ReconstructedMol> valid_population;
  for (const ReconstructedMol& individual : new_population) {
    // Verify that the individuals have connections. Otherwise they can't evolve.
    if (individual.GetConnections().Empty() && individual.GetInternalConnections().Empty()) {
      std::cout << "WARNING: " << RDKit::MolToSmiles(individual.GetPseudomol()) << " has no Connections and won't be added to the population." << std::endl;
      continue;
    };
    // When a population comes from outside sources we can't be certain that the
    // Connections within the ReconstructedMols are also recorded in the
    // ConnectionQueryResults. The ReconstructedMols with unknown Connections are
    // skipped, since our operators can't handle them.
    std::pair<Connection, bool> connection_known = individual.HasKnownConnections(query_results);
    if (!connection_known.second) {
      std::cout << "WARNING: " << RDKit::MolToSmiles(individual.GetPseudomol()) << " contains an unknown Connection " << connection_known.first.GetString() << " and won't be added to the population." << std::endl;
      continue;
    }
    if (id > max_individual_id) {
      max_individual_id = id;
    };
    valid_population.push_back(individual);
    valid_population.back().SetID(id++);
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
  SavePopulation(settings.GetRestartOutputFile());
};

void LEADD::SavePopulation(const std::string& file_path) const {
  std::ofstream output_stream(file_path, std::ofstream::binary);
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
    std::uniform_int_distribution<unsigned> distribution (0, population.size() - 1);
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
    bool evolved = child.Evolve(settings, database, query_results, controller, prng, inventory, report);
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
  float score = 0.0;
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

unsigned LEADD::ScoreChildrenBySimilarity() {
  assert(reference != nullptr);
  unsigned n_scored_molecules = 0;
  for (ReconstructedMol& individual : population) {
    if (individual.IsChild()) {
      individual.Score(reference);
      ++n_scored_molecules;
    };
  };
  return n_scored_molecules;
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

unsigned LEADD::ScoreChildren() {
  unsigned n_scored_molecules = 0;
  // If an external scoring function was specified call it.
  if (settings.UsingExternalScoringFunction()) {
    n_scored_molecules = ScoreChildrenExternally();
  // Otherwise use the similarity to the template molecule as a score.
  } else {
    n_scored_molecules = ScoreChildrenBySimilarity();
  };
  // Update the amount of scored molecules.
  n_scoring_calls += n_scored_molecules;
  return n_scored_molecules;
};

void LEADD::IncreaseNScoringCalls(unsigned n) {
  n_scoring_calls += n;
};

void LEADD::UpdateReport() {
  report.WriteScores(generation_n, n_scoring_calls, population);
  report.WriteNRings(generation_n, population);
};

void LEADD::WriteOperationFrequenciesToReport() {
  report.WriteOperationFrequencies();
};

void LEADD::ReportOnBestMolecule() {
  ReconstructedMol& rank1 = population.front();
  std::string score_text = std::to_string(rank1.GetScore());
  std::filesystem::path base_path (output_directory_path);
  std::filesystem::path drawing_path = base_path.append("rank1_" + score_text + ".svg");
  std::filesystem::path weights_path = base_path.append("rank1_weights_" + score_text + ".csv");
  rank1.Draw(drawing_path);
  report.WriteWeights(rank1, weights_path);
};

void LEADD::Cleanup() {
  // Close the EvolutionReport.
  if (settings.WriteReport()) {
    report.Close();
  };

  // Close the connection to the database.
  database.Close();

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

std::list<ReconstructedMol> MakePopulationFromSMILES(const std::vector<std::string>& population_smiles, const FragmentationSettings& settings) {
  std::list<ReconstructedMol> population;
  for (const std::string& smiles : population_smiles) {
    RDKit::ROMol* molecule = RDKit::SmilesToMol(smiles);
    population.emplace_back(*molecule, settings);
    delete molecule;
  };
  return population;
};
