#include <filesystem>
#include <boost/program_options.hpp>
#include "ReconstructionSettings.hpp"
#include "ConnectionQueryResults.hpp"

// Function that retrieves a list of all unique Connections in the database.
CONNECTIONS_SET GetConnectionsFromDB(sqlite3* database) {
  // Set up a Connection string formatter.
  boost::format formatter("(%d, %d, %d)");

  // Initialize the vector to store the Connections.
  CONNECTIONS_SET connections;

  // Prepare the SQLite3 statement to retrieve the Connections from the database.
  sqlite3_stmt* select_statement;
  sqlite3_prepare_v2(database, "SELECT start_atom_type, end_atom_type, bond_type FROM connections", -1, &select_statement, NULL);

  // Create Connection objects using the data stored in the database.
  std::uint8_t start_atom_type, end_atom_type, bond_type;
  int rc = sqlite3_step(select_statement);
  while (rc == SQLITE_ROW) {
    start_atom_type = sqlite3_column_int(select_statement, 0);
    end_atom_type = sqlite3_column_int(select_statement, 1);
    bond_type = sqlite3_column_int(select_statement, 2);
    Connection connection(start_atom_type, end_atom_type, bond_type, formatter);
    connections.insert(connection);
    rc = sqlite3_step(select_statement);
  };
  return connections;
};

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
      "Integer between 0 and 2 symbolizing the stringency of the connection compatibility definition.");

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
    throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
  };
  std::filesystem::path input_path(input);
  if (!std::filesystem::exists(input_path)) {
    throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
  };

  // Initialize a ReconstructionSettings object.
  ReconstructionSettings settings = ReconstructionSettings();
  if (!settings_file.empty()) {
    std::filesystem::path settings_path(settings_file);
    if (!std::filesystem::exists(settings_path)) {
      throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
    } else {
      settings = ReconstructionSettings(settings_file, false);
    };
  };

  // Open a connection to the input database.
  sqlite3* database;
  sqlite3_open(input.c_str(), &database);

  // Retrieve the connections from the database and store them in a set.
  CONNECTIONS_SET connections = GetConnectionsFromDB(database);

  // Determine pairwise Connection compatibilities according to the user-specified
  // stringency criterion.
  ConnectionCompatibilities compatibilities(connections, stringency);

  // Compute which Pseudofragments are compatible with each Connection.
  ConnectionQueryResults query_results (compatibilities, database, settings.GetAcyclicFrequencyGamma(), settings.GetAcyclicLevelGamma(), settings.GetRingFrequencyGamma(), settings.GetRingLevelGamma());

  // Store the data in a binary file.
  std::ofstream outstream(output, std::ofstream::binary);
  boost::archive::binary_oarchive binar(outstream);
  binar << query_results;

  // Signal success.
  return 0;
};
