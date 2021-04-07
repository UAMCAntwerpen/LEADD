#include <fstream>
#include <filesystem>
#include <boost/progress.hpp>
#include <boost/program_options.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include "SAScore.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser and parse the command line arguments.
  std::string input, output;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to the input file containing molecules in SMILES or SDF format.")
    ("output,o", boost::program_options::value<std::string>(&output)->required(),
      "Path to the output file containing the serialized feature library.");

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

  // Initialize an empty FeatureLibrary.
  FeatureLibrary feature_library;

  // Create a molecule input stream based on the extension of the input file.
  std::filesystem::path extension = std::filesystem::path(input).extension();

  // If the input file is a SMILES file:
  if (extension == ".smi" || extension == ".ism") {
    RDKit::SmilesMolSupplier supplier(input);

    // Initialize a progress bar.
    unsigned n_mols = supplier.length();
    supplier.reset();
    boost::progress_display progress(n_mols);

    // Add the features of the input molecules to the FeatureLibrary.
    std::unique_ptr<RDKit::ROMol> mol;
    while (!supplier.atEnd()) {
      mol.reset(supplier.next());
      // Verify that the molecule was successfully read.
      if (mol) {
        feature_library.AddMolecule(*mol);
      };
      ++progress;
    };

  // If the input file is a SDF file:
  } else if (extension == ".sdf") {
    RDKit::SDMolSupplier supplier(input);

    // Initialize a progress bar.
    unsigned n_mols = supplier.length();
    supplier.reset();
    boost::progress_display progress(n_mols);

    // Add the features of the input molecules to the FeatureLibrary.
    std::unique_ptr<RDKit::ROMol> mol;
    while (!supplier.atEnd()) {
      mol.reset(supplier.next());
      // Verify that the molecule was successfully read.
      if (mol) {
        feature_library.AddMolecule(*mol);
      };
      ++progress;
    };

  // If the molecule input format is something else, throw an error.
  } else {
    throw std::runtime_error("Unrecognized input format. Input should be either SMILES or SDF.");
  };

  // Calculate the scores of the features in the FeatureLibrary.
  feature_library.CalcFeatureScores();

  // Write out the FeatureLibrary to the output file.
  std::ofstream output_stream (output, std::ofstream::binary);
  boost::archive::binary_oarchive binary_archive (output_stream);
  binary_archive << feature_library;

  // Signal success.
	return 0;
};
