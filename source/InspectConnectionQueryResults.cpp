#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/program_options.hpp>
#include "ConnectionQueryResults.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser.
  std::string input;
  bool verbose = false, very_verbose = false;
  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to a serialized ConnectionQueryResults object (.cqr file).")
    ("verbose,v", boost::program_options::bool_switch(&verbose)->default_value(false),
      "Flag to print the ConnectionCompatibilities to standard output.")
    ("very_verbose,y", boost::program_options::bool_switch(&very_verbose)->default_value(false),
      "Flag to print the ConnectionQueryResults to standard output.");
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
    throw std::invalid_argument("Input CQR path is invalid.");
  };

  // Deserialize the ConnectionQueryResults file.
  ConnectionQueryResults query_results;
  std::ifstream instream (input, std::ifstream::binary);
  boost::archive::binary_iarchive binar (instream);
  binar >> query_results;
  instream.close();

  // Print the parameters to the standard output.
  std::cout << "Stringency: " << query_results.GetConnectionCompatibilities().GetStringency() << std::endl;
  std::cout << "Acyclic fragment frequency gamma: " << query_results.GetAcyclicFrequencyGamma() << std::endl;
  std::cout << "Acyclic fragment size gamma: " << query_results.GetAcyclicSizeGamma() << std::endl;
  std::cout << "Ring fragment frequency gamma: " << query_results.GetRingFrequencyGamma() << std::endl;
  std::cout << "Ring fragment size gamma: " << query_results.GetRingSizeGamma() << std::endl;

  // Calculate and print some statistics to the standard output.
  const ConnectionCompatibilities& compatibilities = query_results.GetConnectionCompatibilities();
  const COMPATIBILITY_TABLE& compatibility_table = compatibilities.GetCompatibilityTable();
  size_t n_connections = compatibility_table.size();
  double average_n_compatible_connections = 0.0;
  for (const auto& [connection, compatible_connections] : compatibility_table) {
    average_n_compatible_connections += compatible_connections.size();
  };
  average_n_compatible_connections /= n_connections;

  double average_n_compatible_fragments_strict = 0.0;
  double average_n_compatible_fragments_lax = 0.0;
  for (const auto& [connection, qrbn] : query_results.GetStrictAcyclicResults()) {
    average_n_compatible_fragments_strict += qrbn.at(1).first.size();
  };
  for (const auto& [connection, qrbn] : query_results.GetStrictRingResults()) {
    average_n_compatible_fragments_strict += qrbn.at(1).first.size();
  };
  for (const auto& [connection, qrbn] : query_results.GetAcyclicResults()) {
    average_n_compatible_fragments_lax += qrbn.at(1).first.size();
  };
  for (const auto& [connection, qrbn] : query_results.GetRingResults()) {
    average_n_compatible_fragments_lax += qrbn.at(1).first.size();
  };
  average_n_compatible_fragments_strict /= n_connections;
  average_n_compatible_fragments_lax /= n_connections;

  std::cout << "# connections: " << n_connections << std::endl;
  std::cout << "Average number of compatible connections per connection: " << average_n_compatible_connections << std::endl;
  std::cout << "Average number of compatible fragments per connection (strict): " << average_n_compatible_fragments_strict << std::endl;
  std::cout << "Average number of compatible fragments per connection (lax): " << average_n_compatible_fragments_lax << std::endl;

  // If requested, print the ConnectionCompatibilities to the standard output.
  if (verbose) {
    compatibilities.Print();
  };

  // If requested, print the ConnectionQueryResults to the standard output.
  if (very_verbose) {
    std::cout << "Strict ConnectionQueryResults" << std::endl;
    query_results.PrintStrict();
    std::cout << "Lax ConnectionQueryResults" << std::endl;
    query_results.Print();
  };

  // Signal success.
  return 0;
};
