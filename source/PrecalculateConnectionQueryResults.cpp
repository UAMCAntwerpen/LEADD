#include <filesystem>
#include <boost/program_options.hpp>
#include "ReconstructionSettings.hpp"
#include "ConnectionQueryResults.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser.
  std::string input, output, settings_file;
  unsigned stringency;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to the input SQLite3 fragments database.")
    ("output,o", boost::program_options::value<std::string>(&output)->required(),
      "Path to the output binary file containing the precalculated connection query results.")
    ("settings,t", boost::program_options::value<std::string>(&settings_file),
      "Path to the settings file containing parameters to calculate the weights of the fragments.")
    ("stringency,s", boost::program_options::value<unsigned>(&stringency)->default_value(1),
      "Integer between 0 and 2 symbolizing the stringency of the lax connection compatibility definition (0 = valence, 1 = lax, 2 = strict).");

  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("input", 1);
  positionals_description.add("output", 1);
  positionals_description.add("stringency", 1);
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

  // Validate that the provided arguments are valid.
  if (stringency > 2) {
    throw std::invalid_argument("--stringency must be 0 <= s <= 2");
  };
  std::filesystem::path input_path(input);
  if (!std::filesystem::exists(input_path)) {
    throw std::invalid_argument("Input database path is invalid.");
  };

  // Initialize a ReconstructionSettings object.
  ReconstructionSettings settings = ReconstructionSettings();
  if (!settings_file.empty()) {
    std::filesystem::path settings_path(settings_file);
    if (!std::filesystem::exists(settings_path)) {
      throw std::invalid_argument("Settings file path is invalid.");
    } else {
      settings = ReconstructionSettings(settings_file, false);
    };
  };

  // Open a connection to the input database.
  PseudofragmentDB database (input);

  // Compute which Pseudofragments are compatible with each Connection.
  ConnectionQueryResults query_results (database, stringency, settings.GetAcyclicFrequencyGamma(), settings.GetAcyclicSizeGamma(), settings.GetRingFrequencyGamma(), settings.GetRingSizeGamma());

  // Store the data in a binary file.
  std::ofstream outstream(output, std::ofstream::binary);
  boost::archive::binary_oarchive binar(outstream);
  binar << query_results;
  outstream.close();

  // Close the database.
  database.Close();

  // Signal success.
  return 0;
};
