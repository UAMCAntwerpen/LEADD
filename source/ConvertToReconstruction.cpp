#include <iostream>
#include <boost/serialization/list.hpp>
#include <boost/program_options.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include "Reconstruction.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser and parse the command line arguments.
  std::string input, output, settings_file;
  bool names_as_scores;

  boost::program_options::options_description description("Options");
  description.add_options()
   ("help,h",
     "Display this help message.")
   ("input,i", boost::program_options::value<std::string>(&input)->required(),
     "Path to the input file containing molecules in SMILES format.")
   ("output,o", boost::program_options::value<std::string>(&output)->required(),
     "Path to the output file containing serialized reconstruction versions of the molecules.")
   ("settings,s", boost::program_options::value<std::string>(&settings_file),
     "Path to the settings file specifying how to fragment molecules.")
   ("nascores,n", boost::program_options::bool_switch(&names_as_scores)->default_value(false),
     "Flag to use the molecules' names as their scores.");

  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("input", 1);
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

  // Try to read the settings file.
  FragmentationSettings settings;
  try {
    settings = FragmentationSettings(settings_file);
    settings.Print();
  } catch (const std::runtime_error& error) {
    std::cout << error.what() << std::endl;
    return 1;
  };

  // Initialize the list to store the ReconstructedMols.
  std::list<ReconstructedMol> reconstructions;

  // Iterate over the input molecules.
  RDKit::SmilesMolSupplier supplier(input);
  while (!supplier.atEnd()) {
    RDKit::ROMOL_SPTR molecule (supplier.next());
    // Convert the molecule to ReconstructedMol.
    reconstructions.emplace_back(*molecule, settings, names_as_scores);
  };

  // Assign unique IDs to the ReconstructedMols.
  unsigned id = 1;
  for (ReconstructedMol& reconstruction : reconstructions) {
    reconstruction.SetID(id++);
  };

  // Store the serialized ReconstructedMols in the output file.
  std::ofstream output_stream(output, std::ofstream::binary);
  boost::archive::binary_oarchive archive(output_stream);
  archive << reconstructions;
  output_stream.close();

  // Signal success.
  return 0;
};
