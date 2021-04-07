#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/program_options.hpp>
#include "ConnectionQueryResults.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser.
  std::string input;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to a serialized ConnectionQueryResults object (.cqr file).");

  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("input", 1);
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

  // Validate that the provided input CQR file exists.
  std::filesystem::path input_path(input);
  if (!std::filesystem::exists(input_path)) {
    throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
  };

  // Deserialize the ConnectionQueryResults file.
  ConnectionQueryResults query_results;
  std::ifstream instream (input, std::ifstream::binary);
  boost::archive::binary_iarchive binar (instream);
  binar >> query_results;

  // Print the parameters to the standard output.
  std::cout << "Acyclic fragment frequency gamma: " << query_results.GetAcyclicFrequencyGamma() << std::endl;
  std::cout << "Acyclic fragment level gamma: " << query_results.GetAcyclicLevelGamma() << std::endl;
  std::cout << "Ring fragment frequency gamma: " << query_results.GetRingFrequencyGamma() << std::endl;
  std::cout << "Ring fragment level gamma: " << query_results.GetRingLevelGamma() << std::endl;

  // Signal success.
  return 0;
};
