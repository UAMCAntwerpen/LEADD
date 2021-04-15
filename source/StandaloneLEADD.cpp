#include <boost/program_options.hpp>
#include "LEADD.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser.
  std::string settings_file, output_directory;
  bool verbose = false;
  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("settings,s", boost::program_options::value<std::string>(&settings_file)->required(),
      "Path to the settings file containing parameters for the run.")
    ("output,o", boost::program_options::value<std::string>(&output_directory)->required(),
      "Path to the output directory in which the output files are stored.")
    ("verbose,v", boost::program_options::bool_switch(&verbose)->default_value(false),
      "Flag to print the generation number and its score to standard output.");
  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("settings", 1);
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

  // Read the settings file and verify that the settings are valid.
  ReconstructionSettings settings;
  try {
    settings = ReconstructionSettings(settings_file);
    settings.CheckScoringFunctionValidity();
    settings.Print();
  } catch (const std::runtime_error& error) {
    std::cout << error.what() << std::endl;
    return 1;
  };

  // Initialize a LEADD object containing a population of ReconstructedMols.
  LEADD leadd (settings, output_directory);

  // Assign preliminary scores to the population before starting the evolution.
  if (settings.ScoreFirstPopulation()) {
    // Label all individuals as "children" since they are the only ones
    // that are scored by the scoring function.
    for (ReconstructedMol& individual : leadd.GetPopulation()) {
      individual.SetChildFlag(true);
    };
    // Score the children.
    leadd.ScoreChildren();
    // Remove the child label
    for (ReconstructedMol& individual : leadd.GetPopulation()) {
      individual.SetChildFlag(false);
    };
  };

  // If required, write the initial statistics to the report.
  if (settings.WriteReport()) {
    leadd.UpdateReport();
  };

  // Main evolution loop. While the termination criteria aren't met:
  while (!leadd.TerminationCriteriaMet()) {
    // Sample some parent individuals and reproduce them asexually to expand the
    // population with children.
    leadd.GenerateChildren();
    // Score the children.
    leadd.ScoreChildren();
    // Retain the best individuals.
    leadd.SelectivePressure();
    // If required, print the generation number and its score to standard output.
    if (verbose) {
      std::cout << leadd.GetGenerationNumber() << " " << leadd.GetBestScore() << std::endl;
    };
  };

  // If requested, draw the best molecule and write out its ConnectionPoint weights.
  if (settings.WriteReport()) {
    leadd.WriteOperationFrequenciesToReport();
    if (settings.MonitorBestMolecule()) {
      leadd.ReportOnBestMolecule();
    };
  };

  // Write out the SMILES of the all individuals within the population alongside
  // their corresponding scores in the second column.
  std::filesystem::path output_molecules_path (output_directory);
  output_molecules_path.append("designed_molecules.smi");
  leadd.WritePopulationSMILES(output_molecules_path);

  // If requested, save the final serialized population in a file.
  if (settings.SaveProgress()) {
    leadd.SavePopulation();
  };

  // Close all open files.
  leadd.Cleanup();

  return 0;
};
