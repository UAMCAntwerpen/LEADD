#include <iostream>
#include <boost/program_options.hpp>
#include "H5Cpp.h"
#include "PseudofragmentDB.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser.
  std::string input, output;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to the input SQLite3 fragments database.")
    ("output,o", boost::program_options::value<std::string>(&output)->required(),
      "Path to the output HDF5 similarity matrix.");

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

  // Determine the number of OpenMP threads to use.
  int n_threads = DetermineOMPNThreads();
  omp_set_num_threads(n_threads);

  // Open a connection to the input database.
  PseudofragmentDB database (input);
  sqlite3_int64 n_pseudofragments = database.GetNPseudofragments();

  // Divide the database into N equally sized chunks, where N is the number of
  // OpenMP threads.
  std::vector<std::pair<unsigned, unsigned>> database_chunks = EquallySizedChunks(n_pseudofragments, n_threads, true);

  // Generate fingerprints for the Pseudofragments in the database.
  std::cout << "Generating fingerprints..." << std::endl;
  std::vector<RDKit::SparseIntVect<std::uint32_t>*> fingerprints (n_pseudofragments);
  #pragma omp parallel for
  for (int thrid = 0; thrid < n_threads; ++thrid) {
    const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];
    PseudofragmentDB::PseudofragmentIterator it (database, chunk.first, chunk.second);
    while (!it.AtEnd()) {
      // Generate the Pseudofragment's fingerprint.
      RDKit::SparseIntVect<std::uint32_t>* fingerprint = it.GetPseudofragment().GetFingerprint();
      // Store the fingerprint.
      sqlite3_int64 pseudofragment_idx = it.GetPseudofragmentID() - 1; // -1 because SQLite3 uses 1-based indexing
      fingerprints[pseudofragment_idx] = fingerprint;
      ++it;
    };
    // Finalize the SQLite3 iteration statement.
    it.Finalize();
  };

  // Close the database.
  database.Close();

  // Create the output HDF5 file.
  H5::H5File matrix_file(output, H5F_ACC_EXCL);

  // Define the dataspaces of the memory buffer that holds the similarity values
  // and the file's similarity matrix (i.e. define their shapes and dimensions).
  hsize_t memory_dimensions[1] = { n_pseudofragments };
  hsize_t file_dimensions[2] = { n_pseudofragments, n_pseudofragments };
  H5::DataSpace memory_dataspace(1, memory_dimensions);
  H5::DataSpace file_dataspace(2, file_dimensions);

  // Create the similarity matrix dataset.
  H5::FloatType datatype(H5::PredType::IEEE_F32LE);
  H5::DataSet similarity_matrix = matrix_file.createDataSet("simatrix", datatype, file_dataspace);

  // Define the dimensions of a row in the similarity matrix.
  hsize_t hyperslab_dimensions[2] = {1, n_pseudofragments};

  // Iterate over all pairs of fingerprints and calculate their topological similarities.
  std::cout << "Calculating similarity values..." << std::endl;
  boost::progress_display progress (n_pseudofragments);
  #pragma omp parallel for
  for (int thrid = 0; thrid < n_threads; ++thrid) {
    // Create a buffer to store similarity values.
    std::vector<float> buffer (n_pseudofragments);
    const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];
    // Loop over the fingerprints in the chunk.
    for (unsigned fragment_idx1 = chunk.first - 1; fragment_idx1 <= chunk.second - 1; ++fragment_idx1) {
      const RDKit::SparseIntVect<std::uint32_t>* fingerprint1 = fingerprints[fragment_idx1];
      // Loop over all the fingerprints.
      for (unsigned fragment_idx2 = 0; fragment_idx2 < n_pseudofragments; ++fragment_idx2) {
        const RDKit::SparseIntVect<std::uint32_t>* fingerprint2 = fingerprints[fragment_idx2];
        // Compute and store the similarity value.
        buffer[fragment_idx2] = RDKit::TanimotoSimilarity(*fingerprint1, *fingerprint2);
      };
      // Write the similarity values as a row in the matrix.
      hsize_t hyperslab_offset[2] = {fragment_idx1, 0};
      #pragma omp critical
      {
        file_dataspace.selectHyperslab(H5S_SELECT_SET, hyperslab_dimensions, hyperslab_offset);
        similarity_matrix.write(buffer.data(), datatype, memory_dataspace, file_dataspace);
        ++progress;
      };
    };
  };

  // Destroy all fingerprints.
  for (RDKit::SparseIntVect<std::uint32_t>* fingerprint : fingerprints) {
    delete fingerprint;
  };

  // Close the HDF5 file.
  similarity_matrix.close();
  matrix_file.close();

  // Signal success.
  return 0;
};
