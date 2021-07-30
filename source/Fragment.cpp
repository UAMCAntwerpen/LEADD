#include <boost/program_options.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include "Fragmentation.hpp"
#include "PseudofragmentDB.hpp"

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser and parse the command line arguments.
  std::string input, settings_file, output;
  unsigned purge_fragment_size_threshold;

  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("input,i", boost::program_options::value<std::string>(&input)->required(),
      "Path to the input file containing molecules in SMILES format.")
    ("output,o", boost::program_options::value<std::string>(&output)->required(),
      "Path to the output SQLite3 database containing the fragments.")
    ("settings,s", boost::program_options::value<std::string>(&settings_file),
      "Path to the settings file specifying how to fragment molecules.")
    ("purge_fragment_size_threshold,p", boost::program_options::value<unsigned>(&purge_fragment_size_threshold)->default_value(0),
      "Size threshold (in number of atoms) below which fragments are deleted.");

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

  // Assert that SQLite3 was compiled with multithreading support.
  assert(sqlite3_threadsafe());

  // Determine the number of OpenMP threads to use.
  int n_threads = DetermineOMPNThreads();
  omp_set_num_threads(n_threads);

  // Open a connection to the output database and, if necessary, initialize the tables.
  PseudofragmentDB database (output);
  // Start a database transaction.
  database.BeginTransaction();

  // Create an input molecules stream.
  RDKit::SmilesMolSupplier molecule_supplier(input);

  // Divide the input file into N equally sized chunks, where N is the number of
  // OpenMP threads.
  unsigned n_molecules = molecule_supplier.length();
  molecule_supplier.reset();
  std::vector<std::pair<unsigned, unsigned>> chunks = EquallySizedChunks(n_molecules, n_threads);

  // Initialize a progress bar.
  boost::progress_display progress(n_molecules);

  // Main loop.
  #pragma omp parallel for shared(settings, database, progress)
  for (int thrid = 0; thrid < n_threads; ++thrid) {
    const std::pair<unsigned, unsigned> chunk = chunks[thrid];
    unsigned chunk_size = chunk.second - chunk.first;

    // Create a local supplier that begins at the start of the file chunk.
    RDKit::SmilesMolSupplier local_molecule_supplier(input);
    local_molecule_supplier.moveTo(chunk.first);

    // Loop over the molecules in the file chunk.
    for (unsigned i = 0; i <= chunk_size; ++i) {
      RDKit::ROMOL_SPTR molecule (local_molecule_supplier.next());

      // If the molecule wasn't succesfully read, skip it.
      if (!molecule) {
        continue;
      };

      // Fragment the molecule into Pseudofragments.
      std::vector<Pseudofragment> pseudofragments = MakePseudofragments(*molecule, settings);

      // Identify Pseudofragments without ConnectionPoints.
      std::set<size_t> pseudofragments_to_remove_indices;
      for (size_t pseudofragment_idx = 0; pseudofragment_idx < pseudofragments.size(); ++pseudofragment_idx) {
        const Pseudofragment& pseudofragment = pseudofragments[pseudofragment_idx];
        if (pseudofragment.GetConnections().Empty()) {
          pseudofragments_to_remove_indices.insert(pseudofragment_idx);
        };
      };

      // Identify Pseudofragments that are too small.
      if (purge_fragment_size_threshold > 1) {
        for (size_t pseudofragment_idx = 0; pseudofragment_idx < pseudofragments.size(); ++pseudofragment_idx) {
          const Pseudofragment& pseudofragment = pseudofragments[pseudofragment_idx];
          if (pseudofragment.GetSize() < purge_fragment_size_threshold) {
            pseudofragments_to_remove_indices.insert(pseudofragment_idx);
          };
        };
      };

      // Remove the aforementioned Pseudofragments.
      for (std::set<size_t>::reverse_iterator it = pseudofragments_to_remove_indices.rbegin(); it != pseudofragments_to_remove_indices.rend(); ++it) {
        pseudofragments.erase(pseudofragments.begin() + *it);
      };

      if (pseudofragments.empty()) {
        continue;
      };

      // Insert the remaining fragments.
      // SQLite3 is thread-safe, but only one write transaction may be active
      // concurrently. This, coupled to the need of executing multiple SQL
      // statements to insert a molecule forces us to use a mutex.
      #pragma omp critical
      {
        // Insert the source molecule in the database.
        sqlite3_int64 source_molecule_id = database.InsertSourceMolecule(*molecule);
        // Insert the Pseudofragments in the database.
        for (const Pseudofragment& pseudofragment : pseudofragments) {
          sqlite3_int64 pseudofragment_id = database.InsertPseudofragment(pseudofragment);
          // Record the source molecule-fragment relationship.
          database.InsertSourceMoleculePseudofragmentRelationship(source_molecule_id, pseudofragment_id);
          // Insert the Pseudofragment's Connections.
          for (const auto& [connection, connection_points] : pseudofragment.GetConnections()) {
            sqlite3_int64 connection_id = database.InsertConnection(connection, connection_points.size());
            // Record the fragment-connection relationship.
            database.InsertPseudofragmentConnectionRelationship(pseudofragment_id, connection_id, connection_points.size());
          };
        };
        // Advance the progress bar.
        ++progress;
      };
    };
  };

  // Commit the database transaction.
  database.CommitTransaction();
  // Optimize the database's layout and query planner.
  database.Optimize();
  // Close the database.
  database.Close();

  // Signal success.
  return 0;
};
