#include <fstream>
#include <filesystem>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include "SAScore.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser and parse the command line arguments.
  std::string input_path, library_path, output_path;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input_path)->required(),
      "Path to the input file containing molecules in SMILES or SDF format.")
    ("library,l", boost::program_options::value<std::string>(&library_path)->required(),
      "Path to the serialized feature library containing the feature scores.")
    ("output,o", boost::program_options::value<std::string>(&output_path)->required(),
      "Path to the output text file containing the SAScores.");

  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("input", 1);
  positionals_description.add("library", 1);
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

  // Deserialize the FeatureLibrary.
  FeatureLibrary feature_library;
  std::ifstream library_stream (library_path, std::ifstream::binary);
  boost::archive::binary_iarchive binary_archive (library_stream);
  binary_archive >> feature_library;

  // Create an output stream for the SAScores.
  std::ofstream output_stream (output_path);

  // Create a molecule input stream based on the extension of the input file.
  std::filesystem::path extension = std::filesystem::path(input_path).extension();

  // If the input file is a SMILES file:
  if (extension == ".smi" || extension == ".ism") {
    RDKit::SmilesMolSupplier supplier(input_path);

    // Initialize a progress bar.
    unsigned n_mols = supplier.length();
    supplier.reset();
    boost::progress_display progress(n_mols);

    // Loop over the input molecules.
    std::unique_ptr<RDKit::ROMol> mol;
    while (!supplier.atEnd()) {
      mol.reset(supplier.next());
      // Verify that the molecule was successfully read.
      if (mol) {
        // Calculate its SAScore and write it to the output file.
        double sascore = SAScore(*mol, feature_library);
      } else {
        std::cout << "WARNING: Molecule couldn't be read." << std::endl;
      };
      ++progress;
    };

  // If the input file is a SDF file:
  } else if (extension == ".sdf") {
    RDKit::SDMolSupplier supplier(input_path);

    // Initialize a progress bar.
    unsigned n_mols = supplier.length();
    supplier.reset();
    boost::progress_display progress(n_mols);

    // Loop over the input molecules.
    std::unique_ptr<RDKit::ROMol> mol;
    while (!supplier.atEnd()) {
      mol.reset(supplier.next());
      // Verify that the molecule was successfully read.
      if (mol) {
        // Calculate its SAScore and write it to the output file.
        double sascore = SAScore(*mol, feature_library);
      } else {
        std::cout << "WARNING: Molecule couldn't be read." << std::endl;
      };
      ++progress;
    };

  // If the molecule input format is something else, throw an error.
  } else {
    throw std::runtime_error("Unrecognized input format. Input should be either SMILES or SDF.");
  };

  // Close the output stream.
  output_stream.close();

  // Signal success.
  return 0;
};
