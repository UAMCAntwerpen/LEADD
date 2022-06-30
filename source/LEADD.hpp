#pragma once
#ifndef _LEADD_HPP_
#define _LEADD_HPP_

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/serialization/library_version_type.hpp>
#include <boost/serialization/list.hpp>
#include "Reconstruction.hpp"
#include "SAScore.hpp"

class LEADD {
private:
  std::string output_directory_path;
  ReconstructionSettings settings;
  PseudofragmentDB database;
  ConnectionQueryResults query_results;
  FeatureLibrary feature_library;
  SAScoreHeuristic sascore_heuristic;
  EvolutionGuide guide;

  std::list<ReconstructedMol> population;
  ReconstructionsInventory inventory;
  EvolutionReport report;

  SizeController controller;
  std::mt19937 prng;

  unsigned generation_n = 0;
  unsigned n_scoring_calls = 0;
  float best_score = 0.0;
  unsigned n_generations_stuck = 0;
  unsigned max_individual_id = 0;

  std::vector<float> parent_weights;
  std::vector<unsigned> parent_indices;

  const RDKit::SparseIntVect<std::uint32_t>* reference = nullptr;

private:
  void MakePopulationTakeBrickOwnership();
  void AssignConnectorWeightsToPopulation();
  void SampleParents();
  unsigned ScoreChildrenBySimilarity();
  unsigned ScoreChildrenExternally();
  void FinishGeneration();

public:
  LEADD(const ReconstructionSettings& settings, const std::string& output_directory_path);

  void SeedPopulation();
  void GrowIndividuals();
  void InitializeRandomPopulation();
  void SetPopulation(const std::list<ReconstructedMol>& new_population, bool reset_weights = true);
  void LoadPopulation();
  void SavePopulation() const;
  void SavePopulation(const std::string& file_path) const;
  void WritePopulationSMILES(const std::string& file_path);

  unsigned GenerateChildren();
  unsigned ScoreChildren();
  unsigned SelectivePressure();
  bool TerminationCriteriaMet();

  void IncreaseNScoringCalls(unsigned n);
  void UpdateReport();
  void WriteOperationFrequenciesToReport();
  void ReportOnBestMolecule();

  void Cleanup();

  unsigned GetGenerationNumber() const;
  float GetBestScore() const;
  ReconstructedMol& GetBestIndividual();
  std::list<ReconstructedMol>& GetPopulation();
  const ReconstructionSettings& GetSettings() const;
  EvolutionReport& GetReport();
};

std::list<ReconstructedMol> MakePopulationFromSMILES(const std::vector<std::string>& population_smiles, const FragmentationSettings& settings);

#endif
