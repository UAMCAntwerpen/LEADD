#include "ConnectionCompatibilities.hpp"

// Definition of the version integer of the ConnectionQueryResults class.
extern const unsigned connection_compatibilities_version = 20210615;

ConnectionCompatibilities::ConnectionCompatibilities() = default;
ConnectionCompatibilities::ConnectionCompatibilities(const PseudofragmentDB& database, unsigned stringency) : stringency(stringency) {
  // Initialize a progress bar.
  std::cout << "Calculating ConnectionCompatibilities..." << std::endl;
  sqlite3_int64 n_connections = database.GetNConnections();
  boost::progress_display progress (n_connections);
  // If the stringency level is 0 all Connections with the same bond type are compatible.
  // This corresponds mostly to a valence model.
  if (stringency == 0) {
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
    // Retrieve the Connections with each bond type.
    std::map<std::uint32_t, std::vector<Connection>> connections_with_bond_type;
    sqlite3_stmt* select_connections_with_bond_type;
    sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT start_atom_type, end_atom_type FROM connections WHERE bond_type = ?", -1, &select_connections_with_bond_type, nullptr);
    for (std::uint32_t bond_type : bond_types) {
      auto [it, inserted] = connections_with_bond_type.insert({bond_type, {}});
      std::vector<Connection>& connections = it->second;
      sqlite3_result_code = sqlite3_bind_int64(select_connections_with_bond_type, 1, bond_type);
      sqlite3_result_code = sqlite3_step(select_connections_with_bond_type);
      while (sqlite3_result_code == SQLITE_ROW) {
        std::uint32_t start_atom_type = sqlite3_column_int64(select_connections_with_bond_type, 0);
        std::uint32_t end_atom_type = sqlite3_column_int64(select_connections_with_bond_type, 1);
        connections.emplace_back(start_atom_type, end_atom_type, bond_type);
        sqlite3_result_code = sqlite3_step(select_connections_with_bond_type);
      };
      sqlite3_clear_bindings(select_connections_with_bond_type);
      sqlite3_reset(select_connections_with_bond_type);
    };
    sqlite3_result_code = sqlite3_finalize(select_connections_with_bond_type);
    // Define the compatibility table.
    for (const auto& [bond_type, connections] : connections_with_bond_type) {
      for (const Connection& connection : connections) {
        compatibility_table.insert({connection, CONNECTIONS_SET(connections.begin(), connections.end())});
        ++progress;
      };
    };
  // If the stringency level is 1 a Connection is considered compatible if:
  //  (1) The Connection's starting atom type has been previously observed
  //      paired with the query Connection's starting atom type.
  //  (2) It shares the query Connection's BondType.
  } else if (stringency == 1) {
    // Record the observed atom pairings.
    std::unordered_map<std::uint32_t, std::unordered_set<std::uint32_t>> observed_atom_pairs;
    std::unordered_map<std::uint32_t, std::unordered_set<std::uint32_t>>::iterator atom_pairs_it;
    // Iterate over the Connections in the database.
    PseudofragmentDB::ConnectionIterator connection_it (database);
    while (!connection_it.AtEnd()) {
      const Connection& connection = connection_it.GetConnection();
      // Insert a entry in the compatibility table for the Connection.
      compatibility_table.insert({connection, {}});
      // Store the start-end atom type pair.
      std::uint32_t start_atom_type = connection.GetStartAtomType();
      std::uint32_t end_atom_type = connection.GetEndAtomType();
      atom_pairs_it = observed_atom_pairs.find(start_atom_type);
      if (atom_pairs_it == observed_atom_pairs.end()) {
        observed_atom_pairs.insert({start_atom_type, {end_atom_type}});
      } else {
        atom_pairs_it->second.insert(end_atom_type);
      };
      // Store the end-start atom type pair.
      atom_pairs_it = observed_atom_pairs.find(end_atom_type);
      if (atom_pairs_it == observed_atom_pairs.end()) {
        observed_atom_pairs.insert({end_atom_type, {start_atom_type}});
      } else {
        atom_pairs_it->second.insert(start_atom_type);
      };
      ++connection_it;
    };
    connection_it.Finalize();
    // Determine the number of OpenMP threads to use.
    int n_threads = DetermineOMPNThreads();
    omp_set_num_threads(n_threads);
    // Divide the database into N equally sized chunks, where N is the number of
    // OpenMP threads.
    std::vector<std::pair<unsigned, unsigned>> database_chunks = EquallySizedChunks(n_connections, n_threads, true);
    // Have each thread loop over its chunk of the database to define the ConnectionCompatibilities.
    #pragma omp parallel for
    for (int thrid = 0; thrid < n_threads; ++thrid) {
      sqlite3_stmt* select_compatible_atom_types;
      int sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT end_atom_type FROM connections WHERE start_atom_type = ? AND bond_type = ?", -1, &select_compatible_atom_types, nullptr);
      const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];
      PseudofragmentDB::ConnectionIterator chunk_it (database, chunk.first, chunk.second);
      while (!chunk_it.AtEnd()) {
        const Connection& connection = chunk_it.GetConnection();
        CONNECTIONS_SET& compatible_connections = compatibility_table[connection];
        std::uint32_t start_atom_type = connection.GetStartAtomType();
        std::uint32_t bond_type = connection.GetBondType();
        // For each compatible starting atom type, fetch all matching Connections.
        for (std::uint32_t compatible_start_atom_type : observed_atom_pairs[start_atom_type]) {
          sqlite3_bind_int64(select_compatible_atom_types, 1, compatible_start_atom_type);
          sqlite3_bind_int64(select_compatible_atom_types, 2, bond_type);
          sqlite3_result_code = sqlite3_step(select_compatible_atom_types);
          while (sqlite3_result_code == SQLITE_ROW) {
            std::uint32_t compatible_end_atom_type = sqlite3_column_int64(select_compatible_atom_types, 0);
            // Record the Connection compatibility.
            compatible_connections.emplace(compatible_start_atom_type, compatible_end_atom_type, bond_type);
            sqlite3_result_code = sqlite3_step(select_compatible_atom_types);
          };
          sqlite3_clear_bindings(select_compatible_atom_types);
          sqlite3_reset(select_compatible_atom_types);
        };
        ++chunk_it;
        #pragma omp critical
        {
          ++progress;
        };
      };
      chunk_it.Finalize();
      sqlite3_finalize(select_compatible_atom_types);
    };
    // // Define the ConnectionCompatibilities.
    // sqlite3_stmt* select_compatible_atom_types;
    // int sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT end_atom_type FROM connections WHERE start_atom_type = ? AND bond_type = ?", -1, &select_compatible_atom_types, nullptr);
    // // Iterate over the database Connections again.
    // connection_it = PseudofragmentDB::ConnectionIterator(database);
    // while (!connection_it.AtEnd()) {
    //   const Connection& connection = connection_it.GetConnection();
    //   auto [it, inserted] = compatibility_table.insert({connection, {}});
    //   CONNECTIONS_SET& compatible_connections = it->second;
    //   std::uint32_t start_atom_type = connection.GetStartAtomType();
    //   std::uint32_t bond_type = connection.GetBondType();
    //   // For each compatible starting atom type, fetch all matching Connections.
    //   for (std::uint32_t compatible_start_atom_type : observed_atom_pairs[start_atom_type]) {
    //     sqlite3_bind_int64(select_compatible_atom_types, 1, compatible_start_atom_type);
    //     sqlite3_bind_int64(select_compatible_atom_types, 2, bond_type);
    //     sqlite3_result_code = sqlite3_step(select_compatible_atom_types);
    //     while (sqlite3_result_code == SQLITE_ROW) {
    //       std::uint32_t compatible_end_atom_type = sqlite3_column_int64(select_compatible_atom_types, 0);
    //       // Record the Connection compatibility.
    //       compatible_connections.emplace(compatible_start_atom_type, compatible_end_atom_type, bond_type);
    //       sqlite3_result_code = sqlite3_step(select_compatible_atom_types);
    //     };
    //     sqlite3_clear_bindings(select_compatible_atom_types);
    //     sqlite3_reset(select_compatible_atom_types);
    //   };
    //   ++progress;
    //   ++connection_it;
    // };
    // connection_it.Finalize();
    // sqlite3_finalize(select_compatible_atom_types);
  // If the stringency level is 2 only the mirrored Connection is considered
  // to be compatible.
  } else if (stringency == 2) {
    PseudofragmentDB::ConnectionIterator it (database);
    while (!it.AtEnd()) {
      const Connection& connection = it.GetConnection();
      compatibility_table.insert({connection, {connection.Mirror()}});
      ++progress;
      ++it;
    };
    it.Finalize();
  // Verify that the provided stringency level is valid.
  } else {
    throw std::invalid_argument("Stringency criterion has to be one of 0, 1 or 2");
  };
};

CONNECTIONS_SET& ConnectionCompatibilities::operator[](const Connection& connection) {
  return compatibility_table[connection];
};

const CONNECTIONS_SET& ConnectionCompatibilities::at(const Connection & connection) const {
  return compatibility_table.at(connection);
};

bool ConnectionCompatibilities::AreCompatible(const Connection& connection1, const Connection& connection2) const {
  // Check if Connection 1 is compatible with Connection 2.
  COMPATIBILITY_TABLE::const_iterator table_it = compatibility_table.find(connection1);
  if (table_it == compatibility_table.end()) {
    throw std::runtime_error("Connection wasn't found in ConnectionCompatibilities.");
  };
  const CONNECTIONS_SET& compatible_connections1 = table_it->second;
  CONNECTIONS_SET::const_iterator compatible_it = compatible_connections1.find(connection2);
  if (compatible_it == compatible_connections1.end()) {
    return false;
  };
  // Check if Connection 2 is compatible with Connection 1.
  table_it = compatibility_table.find(connection2);
  if (table_it == compatibility_table.end()) {
    throw std::runtime_error("Connection wasn't found in ConnectionCompatibilities.");
  };
  const CONNECTIONS_SET& compatible_connections2 = table_it->second;
  compatible_it = compatible_connections2.find(connection1);
  if (compatible_it == compatible_connections2.end()) {
    return false;
  };
  return true;
};

bool ConnectionCompatibilities::AreCompatible(const ConnectionsTable& connections_table, const Connection& connection) const {
  const CONNECTIONS_SET& compatible_connections = compatibility_table.at(connection);
  CONNECTIONS_SET::const_iterator it, end_it = compatible_connections.end();
  for (const auto& [c, cpoints] : connections_table) {
    it = compatible_connections.find(c);
    if (it != end_it) {
      return true;
    };
  };
  return false;
};

Connection ConnectionCompatibilities::GetRandomCompatibleConnection(const ConnectionsTable& connections_table, const Connection& connection, std::mt19937& prng) const {
  CONNECTIONS_VECTOR connection_candidates;
  const CONNECTIONS_SET& compatible_connections = compatibility_table.at(connection);
  CONNECTIONS_SET::const_iterator it, end_it = compatible_connections.end();
  for (const auto& [c, cpoints] : connections_table) {
    it = compatible_connections.find(c);
    if (it != end_it) {
      connection_candidates.push_back(c);
    };
  };
  // Since the ConnectionsTable is a std::unordered_map its iteration order isn't
  // well defined. We sort the vector to ensure we get the same results with the
  // same PRNG seed.
  std::sort(connection_candidates.begin(), connection_candidates.end());
  unsigned idx = 0;
  unsigned n_compatible_connections = connection_candidates.size();
  assert(n_compatible_connections > 0);
  if (n_compatible_connections > 1) {
    if (n_compatible_connections <= max_size_uniform_int_distributions) {
      idx = uniform_int_distributions[n_compatible_connections](prng);
    } else {
      std::uniform_int_distribution<unsigned> distribution(0, n_compatible_connections - 1);
      idx = distribution(prng);
    };
  };
  return connection_candidates[idx];
};

unsigned ConnectionCompatibilities::GetStringency() const {
  return stringency;
};

bool ConnectionCompatibilities::HasConnection(const Connection& connection) const {
  if (compatibility_table.find(connection) == compatibility_table.end()) {
    return false;
  };
  return true;
};

const COMPATIBILITY_TABLE& ConnectionCompatibilities::GetCompatibilityTable() const {
  return compatibility_table;
};

void ConnectionCompatibilities::Print() const {
  for (const auto& c : compatibility_table) {
    std::cout << c.first.GetString() << ": ";
    for (const Connection& connection : c.second) {
      std::cout << connection.GetString() << " ";
    };
    std::cout << std::endl;
  };
};
