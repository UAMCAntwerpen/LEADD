#include <iostream>
#include <omp.h>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include "H5Cpp.h"
#include "PseudofragmentDB.hpp"

std::vector<std::pair<unsigned, unsigned>> EquallySizedChunks(unsigned size, unsigned n_chunks, bool one_based_index = false) {
  std::vector<std::pair<unsigned, unsigned>> chunks;
  unsigned chunk_size = size / n_chunks;
  // Define up to the last chunk.
  unsigned begin_idx = 0, end_idx = 0;
  if (one_based_index) {
    ++begin_idx;
  };
  for (unsigned n = 0; n < n_chunks - 1; ++n) {
    end_idx = begin_idx + chunk_size - 1;
    chunks.push_back(std::make_pair(begin_idx, end_idx));
    begin_idx = end_idx + 1;
  };
  // Define the last chunk.
  end_idx = size;
  if (!one_based_index) {
    --end_idx;
  };
  chunks.push_back(std::make_pair(begin_idx, end_idx));
  return chunks;
};

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

  // Open a connection to the input database.
  sqlite3* database;
  sqlite3_open(input.c_str(), &database);

  // Assert that SQLite3 was compiled with multithreading support.
  assert(sqlite3_threadsafe());

  // Initialize a vector to store a row of similarity values.
  unsigned n_pseudofragments = GetPseudofragmentCount(database);
  assert(n_pseudofragments > 0);
  std::vector<float> similarity_array(n_pseudofragments);

  // Create the output HDF5 file.
  H5::H5File matrix_file(output, H5F_ACC_EXCL);

  // Define the datatype of the matrix's elements.
  H5::FloatType datatype(H5::PredType::NATIVE_FLOAT);
  datatype.setOrder(H5T_ORDER_LE);

  // Define the dataspaces of the memory buffer that holds the similarity values
  // and the file's similarity matrix (i.e. define their shapes and dimensions).
  hsize_t memory_dimensions[1] = { n_pseudofragments };
  hsize_t file_dimensions[2] = { n_pseudofragments, n_pseudofragments };
  H5::DataSpace memory_dataspace(1, memory_dimensions);
  H5::DataSpace file_dataspace(2, file_dimensions);

  // Create the similarity matrix dataset.
  H5::DataSet similarity_matrix = matrix_file.createDataSet("simatrix", datatype, file_dataspace);

  // Define the dimensions of a row in the similarity matrix.
  hsize_t size[2] = {1, n_pseudofragments};

  // Determine the number of OpenMP threads to use.
  int n_threads;
  char* omp_num_threads = std::getenv("OMP_NUM_THREADS");
  if (omp_num_threads == nullptr) {
    n_threads = omp_get_max_threads();
  } else {
    n_threads = std::stoi(omp_num_threads);
  };
  omp_set_num_threads(n_threads);

  // Divide the database into N equally sized chunks, where N is the number of
  // OpenMP threads.
  std::vector<std::pair<unsigned, unsigned>> database_chunks = EquallySizedChunks(n_pseudofragments, n_threads, true);

  // Initialize a progress bar.
  boost::progress_display progress (n_pseudofragments);

  #pragma omp parallel for default(shared)
  for (int thrid = 0; thrid < n_threads; ++thrid) {
    const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];

    // Create a buffer to store similarity values between a reference fragment
    // and all the other fragments in the database.
    std::vector<float> buffer;
    buffer.reserve(n_pseudofragments);

    // Create a pair of iterators to traverse the database, one to traverse the
    // specified chunk and one to traverse the entire database for each fragment
    // of the chunk, defining all possible pairs of Pseudofragments.
    PseudofragmentIterator chunk_it (database, chunk.first, chunk.second);
    PseudofragmentIterator db_it (database);
    assert(chunk_it.GetResultCode() == SQLITE_ROW);
    assert(db_it.GetResultCode() == SQLITE_ROW);

    // Define the offset of the HDF5 row selection hyperslab (i.e. row idx).
    hsize_t offset[2] = {chunk_it.GetID() - 1, 0};

    // Loop over the corresponding pairs of Pseudofragments in the database. For each
    // pair calculate the similarity coefficient and store it in the similarities
    // array. Once the array is full (e.g. all pairs with a single Pseudofragment
    // have been examined) transfer it to the HDF5 matrix.
    while (chunk_it.GetResultCode() == SQLITE_ROW) {
      // Get the fingerprint for the reference Pseudofragment.
      const RDKit::SparseIntVect<std::uint32_t>* fp1 = chunk_it->GetFingerprint();
      // Calculate the similarity between the reference and all the other Pseudofragments.
      while (db_it.GetResultCode() == SQLITE_ROW) {
        const RDKit::SparseIntVect<std::uint32_t>* fp2 = db_it->GetFingerprint();
        buffer.push_back(RDKit::TanimotoSimilarity(*fp1, *fp2));
        delete fp2;
        ++db_it;
      };
      // Store the similarity values as a row in the matrix.
      #pragma omp critical
      {
        file_dataspace.selectHyperslab(H5S_SELECT_SET, size, offset);
        similarity_matrix.write(buffer.data(), datatype, memory_dataspace, file_dataspace);
        ++progress;
      };
      // Reset the nested loop and move on to the next row.
      buffer.clear();
      db_it.Reset();
      delete fp1;
      ++offset[0];
      ++chunk_it;
    };
  };

  // Close the database and HDF5 file.
  similarity_matrix.close();
  matrix_file.close();
  sqlite3_close(database);

  // Signal success.
  return 0;
};
