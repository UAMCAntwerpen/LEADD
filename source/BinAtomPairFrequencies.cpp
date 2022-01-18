#include <iostream>
#include <fstream>
#include <filesystem>
#include <boost/program_options.hpp>
#include "ConnectionQueryResults.hpp"

template <typename T>
class Bins {
  std::vector<T> bin_edges;
  typename std::vector<T>::const_iterator it, begin_it, end_it;
  size_t n_bins = 0;
  size_t max_bin_idx = 0;
  T lower_bin_edge_bound;
  T upper_bin_edge_bound;
  std::vector<std::pair<T, T>> bin_bounds;
  std::map<size_t, unsigned> bins;

public:
  Bins(const std::vector<T>& bin_edges) :
    bin_edges(bin_edges),
    begin_it(bin_edges.begin()),
    end_it(bin_edges.end()),
    n_bins(bin_edges.size() - 1),
    max_bin_idx(n_bins -1),
    lower_bin_edge_bound(bin_edges[0]),
    upper_bin_edge_bound(bin_edges[n_bins]) {
    bin_bounds.resize(n_bins);
    for (size_t bin_idx = 0; bin_idx < n_bins; ++bin_idx) {
      double bin_lower_bound = bin_edges[bin_idx];
      double bin_upper_bound = bin_edges[bin_idx + 1];
      bin_bounds[bin_idx] = {bin_lower_bound, bin_upper_bound};
      bins.insert({bin_idx, 0});
    };
  };

  std::map<size_t, unsigned>::const_iterator begin() const {
    return bins.begin();
  };

  std::map<size_t, unsigned>::const_iterator end() const {
    return bins.end();
  };

  bool Add(const T& value) {
    if (value < lower_bin_edge_bound || value > upper_bin_edge_bound) {
      return false;
    };
    it = std::upper_bound(begin_it, end_it, value);
    size_t bin_idx = max_bin_idx;
    if (it != end_it) {
      bin_idx = std::distance(begin_it, it) - 1;
    };
    bins[bin_idx] += 1;
    return true;
  };

  std::pair<T, T> BinBounds(size_t bin_idx) const {
    return bin_bounds[bin_idx];
  };

  std::string BinString(size_t bin_idx) const {
    auto[bin_lower_bound, bin_upper_bound] = bin_bounds.at(bin_idx);
    std::stringstream bin_string;
    bin_string << "[" << bin_lower_bound << ", " << bin_upper_bound;
    if (bin_idx == max_bin_idx) {
      bin_string << "]";
    } else {
      bin_string << ")";
    };
    return bin_string.str();
  };

  void Print() const {
    for (auto[bin_idx, bin_count] : bins) {
      std::cout << BinString(bin_idx) << ": " << bin_count << std::endl;
    };
  };
};

int main(int argc, const char* argv[]) {
  // Set up a command line argument parser.
  std::string database_path, cqr_path, output;
  boost::program_options::options_description description("Options");
  description.add_options()
    ("help,h",
      "Display this help message.")
    ("database,d", boost::program_options::value<std::string>(&database_path)->required(),
      "Path to a input SQLite3 fragments database.")
    ("cqr,c", boost::program_options::value<std::string>(&cqr_path)->required(),
      "Path to a input serialized ConnectionQueryResults object (.cqr file).")
    ("output,o", boost::program_options::value<std::string>(&output)->required(),
      "Path to an output CSV file containing the pollution data.");
  boost::program_options::positional_options_description positionals_description;
  positionals_description.add("database", 1);
  positionals_description.add("cqr", 1);
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

  // Validate that the provided input files exist.
  if (!std::filesystem::exists(std::filesystem::path(database_path))) {
    throw std::invalid_argument("Input database path is invalid.");
  };
  if (!std::filesystem::exists(std::filesystem::path(cqr_path))) {
    throw std::invalid_argument("Input CQR path is invalid.");
  };

  // Open a connection to the input database.
  PseudofragmentDB database (database_path);

  // Deserialize the ConnectionQueryResults file.
  ConnectionQueryResults query_results;
  std::ifstream instream (cqr_path, std::ifstream::binary);
  boost::archive::binary_iarchive binar (instream);
  binar >> query_results;
  instream.close();

  // Retrieve all unique bond types.
  std::vector<std::uint32_t> bond_types;
  sqlite3_stmt* select_bond_types;
  int sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT DISTINCT bond_type FROM connections", -1, &select_bond_types, nullptr);
  sqlite3_result_code = sqlite3_step(select_bond_types);
  while (sqlite3_result_code == SQLITE_ROW) {
    std::uint32_t bond_type = sqlite3_column_int64(select_bond_types, 0);
    bond_types.push_back(bond_type);
    sqlite3_result_code = sqlite3_step(select_bond_types);
  };
  sqlite3_result_code = sqlite3_finalize(select_bond_types);

  // Retrieve the connection frequencies.
  std::unordered_map<Connection, double> connection_frequencies;
  PseudofragmentDB::ConnectionIterator db_it (database);
  while (!db_it.AtEnd()) {
    const Connection& connection = db_it.GetConnection();
    double connection_frequency = db_it.GetConnectionFrequency();
    connection_frequencies.insert({connection, connection_frequency});
    ++db_it;
  };
  db_it.Finalize();

  // Initialize bins to count the occurrences of connection frequency distances.
  std::vector<double> bin_edges {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35,
                                 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,
                                 0.80, 0.85, 0.90, 0.95, 1.00};
  Bins bins (bin_edges);

  // Iterate over the Connections.
  const ConnectionCompatibilities& compatibilities = query_results.GetConnectionCompatibilities();
  for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
    std::uint32_t start_atom_type = connection.GetStartAtomType();
    // Retrieve the frequencies of the atom pairs responsible for Connection compatibility.
    std::vector<double> atom_pair_frequencies;
    atom_pair_frequencies.reserve(compatible_connections.size());
    for (const Connection& compatible_connection : compatible_connections) {
      std::uint32_t compatible_start_atom_type = compatible_connection.GetStartAtomType();
      // Get the frequency of the atom pair.
      std::unordered_map<Connection, double>::const_iterator cfreq_it, cfreq_end_it = connection_frequencies.end();
      double atom_pair_frequency = 0.0;
      for (std::uint32_t bond_type : bond_types) {
        Connection tmp_connection (start_atom_type, compatible_start_atom_type, bond_type);
        cfreq_it = connection_frequencies.find(tmp_connection);
        if (cfreq_it != cfreq_end_it) {
          atom_pair_frequency += cfreq_it->second;
        };
      };
      atom_pair_frequencies.push_back(atom_pair_frequency);
    };
    std::sort(atom_pair_frequencies.begin(), atom_pair_frequencies.end(), std::greater<double>());
    // Calculate the frequency difference between the most prevalent atom pair and the rest.
    double max_atom_pair_frequency = atom_pair_frequencies[0];
    for (size_t i = 1; i < atom_pair_frequencies.size(); ++i) {
      double atom_pair_frequency_difference = 1.0 - (atom_pair_frequencies[i] / max_atom_pair_frequency);
      bins.Add(atom_pair_frequency_difference);
    };
  };

  // Print the bin counts to standard output.
  bins.Print();

  // Write the bin counts to a CSV file.
  std::ofstream outstream (output);
  for (auto[bin_idx, bin_count] : bins) {
    auto[bin_lower_bound, bin_upper_bound] = bins.BinBounds(bin_idx);
    outstream << bin_lower_bound << "," << bin_upper_bound << "," << bin_count << "\n";
  };
  outstream.close();

  // Signal success.
  return 0;
};
