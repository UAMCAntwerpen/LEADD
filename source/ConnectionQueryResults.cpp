#include "ConnectionQueryResults.hpp"

// Definition of the version integer of the ConnectionQueryResults class.
extern const unsigned connection_query_results_version = 20210615;

// Class ConnectionQueryResults
ConnectionQueryResults::ConnectionQueryResults() = default;
ConnectionQueryResults::ConnectionQueryResults(const PseudofragmentDB& database, unsigned stringency, float acyclic_frequency_gamma, float acyclic_size_gamma, float ring_frequency_gamma, float ring_size_gamma) :
  compatibilities(database, stringency),
  acyclic_frequency_gamma(acyclic_frequency_gamma),
  acyclic_size_gamma(acyclic_size_gamma),
  ring_frequency_gamma(ring_frequency_gamma),
  ring_size_gamma(ring_size_gamma) {
  CalculateStrictQueryResults(database);
  CalculateQueryResults(database);
};

void ConnectionQueryResults::CalculateStrictQueryResults(const PseudofragmentDB& database) {
  // Ensure the containers are empty.
  strict_acyclic.clear();
  strict_ring.clear();
  strict_ring_part.clear();
  // Initialize the containers with the appropiate structure.
  for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
    strict_acyclic.insert({ connection, QUERY_RESULTS_BY_NUM{} });
  };
  strict_ring = strict_acyclic;
  strict_ring_part = strict_acyclic;
  // Determine the number of OpenMP threads to use.
  int n_threads = DetermineOMPNThreads();
  omp_set_num_threads(n_threads);
  // Divide the input file into N equally sized chunks, where N is the number of
  // OpenMP threads.
  sqlite3_int64 n_connections = database.GetNConnections();
  std::vector<std::pair<unsigned, unsigned>> database_chunks = EquallySizedChunks(n_connections, n_threads, true);
  // Initialize a progress bar.
  std::cout << "Calculating ConnectionQueryResults according to the strict compatibility definition..." << std::endl;
  boost::progress_display progress(n_connections);
  #pragma omp parallel for
  for (int thrid = 0; thrid < n_threads; ++thrid) {
    // Prepare the necessary SQLite3 statements.
    sqlite3_stmt* select_connection_id;
    sqlite3_stmt* select_max_connection_frequency;
    sqlite3_stmt* select_pseudofragments_with_connection;
    int sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT id FROM connections WHERE start_atom_type = ? AND end_atom_type = ? AND bond_type = ?", -1, &select_connection_id, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT MAX(frequency) FROM pseudofragments_connections WHERE connection_id = ?", -1, &select_max_connection_frequency, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT id, ring_part, has_ring, size, frequency FROM pseudofragments WHERE id IN (SELECT pseudofragment_id FROM pseudofragments_connections WHERE connection_id = ? AND frequency >= ?)", -1, &select_pseudofragments_with_connection, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    // Iterate over the Connections in the database chunk.
    const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];
    PseudofragmentDB::ConnectionIterator chunk_it (database, chunk.first, chunk.second);
    while (!chunk_it.AtEnd()) {
      const Connection& connection = chunk_it.GetConnection();
      // According to the strict compatibility definition the mirrored Connection
      // is the only one considered compatible.
      Connection mirrored_connection = connection.Mirror();
      sqlite3_bind_int64(select_connection_id, 1, mirrored_connection.GetStartAtomType());
      sqlite3_bind_int64(select_connection_id, 2, mirrored_connection.GetEndAtomType());
      sqlite3_bind_int64(select_connection_id, 3, mirrored_connection.GetBondType());
      sqlite3_result_code = sqlite3_step(select_connection_id);
      assert(sqlite3_result_code == SQLITE_ROW);
      sqlite3_int64 mirrored_connection_id = sqlite3_column_int64(select_connection_id, 0);
      sqlite3_clear_bindings(select_connection_id);
      sqlite3_reset(select_connection_id);
      // Determine the maximum number of instances of the mirrored Connection on a
      // single Pseudofragment.
      sqlite3_bind_int64(select_max_connection_frequency, 1, mirrored_connection_id);
      sqlite3_result_code = sqlite3_step(select_max_connection_frequency);
      assert(sqlite3_result_code == SQLITE_ROW);
      unsigned max_connection_frequency = sqlite3_column_int64(select_max_connection_frequency, 0);
      sqlite3_clear_bindings(select_max_connection_frequency);
      sqlite3_reset(select_max_connection_frequency);
      // Iterate over the number of mirrored Connection instances.
      QUERY_RESULTS_BY_NUM& acyclic_qrbn = strict_acyclic[connection];
      QUERY_RESULTS_BY_NUM& ring_qrbn = strict_ring[connection];
      QUERY_RESULTS_BY_NUM& ring_part_qrbn = strict_ring_part[connection];
      for (unsigned connection_frequency = 1; connection_frequency <= max_connection_frequency; ++connection_frequency) {
        // Create the containers to store the IDs and weights of Pseudofragments with
        // N instances of the mirrored connection.
        QUERY_RESULTS_BY_NUM::iterator qrbn_it;
        bool inserted = false;
        std::tie(qrbn_it, inserted) = acyclic_qrbn.insert({connection_frequency, QUERY_RESULT()});
        QUERY_RESULT& acyclic_qr = qrbn_it->second;
        std::tie(qrbn_it, inserted) = ring_qrbn.insert({connection_frequency, QUERY_RESULT()});
        QUERY_RESULT& ring_qr = qrbn_it->second;
        std::tie(qrbn_it, inserted) = ring_part_qrbn.insert({connection_frequency, QUERY_RESULT()});
        QUERY_RESULT& ring_part_qr= qrbn_it->second;
        // Iterate over the Pseudofragments that have at least the specified number
        // of mirrored Connection instances.
        sqlite3_bind_int64(select_pseudofragments_with_connection, 1, mirrored_connection_id);
        sqlite3_bind_int64(select_pseudofragments_with_connection, 2, connection_frequency);
        sqlite3_result_code = sqlite3_step(select_pseudofragments_with_connection);
        while (sqlite3_result_code == SQLITE_ROW) {
          sqlite3_int64 pseudofragment_id = sqlite3_column_int64(select_pseudofragments_with_connection, 0);
          bool is_ring_part = sqlite3_column_int(select_pseudofragments_with_connection, 1);
          bool has_ring = sqlite3_column_int(select_pseudofragments_with_connection, 2);
          unsigned pseudofragment_size = sqlite3_column_int(select_pseudofragments_with_connection, 3);
          unsigned pseudofragment_frequency = sqlite3_column_int64(select_pseudofragments_with_connection, 4);
          // Store the Pseudofragment's ID and weight in the appropiate containers.
          if (!has_ring && !is_ring_part) {
            float pseudofragment_weight = std::pow(pseudofragment_frequency, acyclic_frequency_gamma) * std::pow(pseudofragment_size, acyclic_size_gamma);
            acyclic_qr.first.push_back(pseudofragment_id);
            acyclic_qr.second.push_back(pseudofragment_weight);
          };
          if (has_ring) {
            float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
            ring_qr.first.push_back(pseudofragment_id);
            ring_qr.second.push_back(pseudofragment_weight);
          };
          if (is_ring_part) {
            float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
            ring_part_qr.first.push_back(pseudofragment_id);
            ring_part_qr.second.push_back(pseudofragment_weight);
          };
          sqlite3_result_code = sqlite3_step(select_pseudofragments_with_connection);
        };
        sqlite3_clear_bindings(select_pseudofragments_with_connection);
        sqlite3_reset(select_pseudofragments_with_connection);
      };
      ++chunk_it;
      #pragma omp critical
      {
        ++progress;
      };
    };
    // Finalize all SQLite3 statements.
    chunk_it.Finalize();
    sqlite3_result_code = sqlite3_finalize(select_connection_id);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_result_code = sqlite3_finalize(select_max_connection_frequency);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_result_code = sqlite3_finalize(select_pseudofragments_with_connection);
    assert(sqlite3_result_code == SQLITE_OK);
  };
  // // Prepare the necessary SQLite3 statements.
  // sqlite3_stmt* select_max_connection_frequency;
  // sqlite3_stmt* select_pseudofragments_with_connection;
  // int sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT MAX(frequency) FROM pseudofragments_connections WHERE connection_id = ?", -1, &select_max_connection_frequency, nullptr);
  // assert(sqlite3_result_code == SQLITE_OK);
  // sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT id, ring_part, has_ring, size, frequency FROM pseudofragments WHERE id IN (SELECT pseudofragment_id FROM pseudofragments_connections WHERE connection_id = ? AND frequency >= ?)", -1, &select_pseudofragments_with_connection, nullptr);
  // assert(sqlite3_result_code == SQLITE_OK);
  // // Iterate over the Connections in the database.
  // std::cout << "Calculating ConnectionQueryResults according to the strict compatibility definition..." << std::endl;
  // for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
  //   CONNECTION_QUERY_RESULTS::iterator cqr_it;
  //   bool inserted = false;
  //   // Create the containers to store the IDs and weights of compatible Pseudofragments.
  //   std::tie(cqr_it, inserted) = strict_acyclic.insert({connection, QUERY_RESULTS_BY_NUM()});
  //   QUERY_RESULTS_BY_NUM& acyclic_qresults = cqr_it->second;
  //   std::tie(cqr_it, inserted) = strict_ring.insert({connection, QUERY_RESULTS_BY_NUM()});
  //   QUERY_RESULTS_BY_NUM& ring_qresults = cqr_it->second;
  //   std::tie(cqr_it, inserted) = strict_ring_part.insert({connection, QUERY_RESULTS_BY_NUM()});
  //   QUERY_RESULTS_BY_NUM& ring_part_qresults = cqr_it->second;
  //   // According to the strict compatibility definition the mirrored Connection
  //   // is the only one considered compatible.
  //   Connection mirrored_connection = connection.Mirror();
  //   sqlite3_int64 mirrored_connection_id = database.SelectConnectionID(mirrored_connection);
  //   // Determine the maximum number of instances of the mirrored Connection on a
  //   // single Pseudofragment.
  //   sqlite3_bind_int64(select_max_connection_frequency, 1, mirrored_connection_id);
  //   sqlite3_result_code = sqlite3_step(select_max_connection_frequency);
  //   assert(sqlite3_result_code == SQLITE_ROW);
  //   unsigned max_connection_frequency = sqlite3_column_int(select_max_connection_frequency, 0);
  //   assert(max_connection_frequency > 0);
  //   sqlite3_clear_bindings(select_max_connection_frequency);
  //   sqlite3_reset(select_max_connection_frequency);
  //   // Iterate over the number of mirrored Connection instances.
  //   for (unsigned connection_frequency = 1; connection_frequency <= max_connection_frequency; ++connection_frequency) {
  //     // Create the containers to store the IDs and weights of Pseudofragments with
  //     // N instances of the mirrored connection.
  //     QUERY_RESULTS_BY_NUM::iterator qrbn_it;
  //     std::tie(qrbn_it, inserted) = acyclic_qresults.insert({connection_frequency, QUERY_RESULT()});
  //     QUERY_RESULT& acyclic_qresult = qrbn_it->second;
  //     std::tie(qrbn_it, inserted) = ring_qresults.insert({connection_frequency, QUERY_RESULT()});
  //     QUERY_RESULT& ring_qresult = qrbn_it->second;
  //     std::tie(qrbn_it, inserted) = ring_part_qresults.insert({connection_frequency, QUERY_RESULT()});
  //     QUERY_RESULT& ring_part_qresult = qrbn_it->second;
  //     // Iterate over the Pseudofragments that have at least the specified number
  //     // of mirrored Connection instances.
  //     sqlite3_bind_int64(select_pseudofragments_with_connection, 1, mirrored_connection_id);
  //     sqlite3_bind_int64(select_pseudofragments_with_connection, 2, connection_frequency);
  //     sqlite3_result_code = sqlite3_step(select_pseudofragments_with_connection);
  //     while (sqlite3_result_code == SQLITE_ROW) {
  //       sqlite3_int64 pseudofragment_id = sqlite3_column_int64(select_pseudofragments_with_connection, 0);
  //       bool is_ring_part = sqlite3_column_int(select_pseudofragments_with_connection, 1);
  //       bool has_ring = sqlite3_column_int(select_pseudofragments_with_connection, 2);
  //       unsigned pseudofragment_size = sqlite3_column_int(select_pseudofragments_with_connection, 3);
  //       unsigned pseudofragment_frequency = sqlite3_column_int(select_pseudofragments_with_connection, 4);
  //       // Store the Pseudofragment's ID and weight in the appropiate containers.
  //       if (!has_ring && !is_ring_part) {
  //         float pseudofragment_weight = std::pow(pseudofragment_frequency, acyclic_frequency_gamma) * std::pow(pseudofragment_size, acyclic_size_gamma);
  //         acyclic_qresult.first.push_back(pseudofragment_id);
  //         acyclic_qresult.second.push_back(pseudofragment_weight);
  //       };
  //       if (has_ring) {
  //         float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
  //         ring_qresult.first.push_back(pseudofragment_id);
  //         ring_qresult.second.push_back(pseudofragment_weight);
  //       };
  //       if (is_ring_part) {
  //         float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
  //         ring_part_qresult.first.push_back(pseudofragment_id);
  //         ring_part_qresult.second.push_back(pseudofragment_weight);
  //       };
  //       sqlite3_result_code = sqlite3_step(select_pseudofragments_with_connection);
  //     };
  //     sqlite3_clear_bindings(select_pseudofragments_with_connection);
  //     sqlite3_reset(select_pseudofragments_with_connection);
  //   };
  // };
  // // Finalize the SQLite3 statements.
  // sqlite3_result_code = sqlite3_finalize(select_max_connection_frequency);
  // assert(sqlite3_result_code == SQLITE_OK);
  // sqlite3_result_code = sqlite3_finalize(select_pseudofragments_with_connection);
  // assert(sqlite3_result_code == SQLITE_OK);
};

// void ConnectionQueryResults::CalculateQueryResults(const PseudofragmentDB& database) {
//   // Ensure the containers are empty.
//   acyclic.clear();
//   ring.clear();
//   ring_part.clear();
//   // Initialize the containers with the appropiate structure.
//   for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
//     acyclic.insert({ connection, QUERY_RESULTS_BY_NUM{{1, QUERY_RESULT{}}} });
//   };
//   ring = acyclic;
//   ring_part = acyclic;
//   // Determine the number of OpenMP threads to use.
//   int n_threads = DetermineOMPNThreads();
//   omp_set_num_threads(n_threads);
//   // Divide the input file into N equally sized chunks, where N is the number of
//   // OpenMP threads.
//   sqlite3_int64 n_connections = database.GetNConnections();
//   std::vector<std::pair<unsigned, unsigned>> database_chunks = EquallySizedChunks(n_connections, n_threads, true);
//   // Initialize a progress bar.
//   std::cout << "Calculating ConnectionQueryResults according to the lax compatibility definition..." << std::endl;
//   boost::progress_display progress(n_connections);
//   #pragma omp parallel for
//   for (int thrid = 0; thrid < n_threads; ++thrid) {
//     // Prepare the necessary SQLite3 statements.
//     sqlite3_stmt* select_pseudofragments_with_connection;
//     int sqlite3_result_code = sqlite3_prepare_v2(database.GetDatabaseConnection(), "SELECT id, ring_part, has_ring, size, frequency FROM pseudofragments WHERE id IN (SELECT pseudofragment_id FROM pseudofragments_connections WHERE start_atom_type = ? AND end_atom_type = ? AND bond_type = ?)", -1, &select_pseudofragments_with_connection, nullptr);
//     assert(sqlite3_result_code == SQLITE_OK);
//
//     // Iterate over the Connections in the database chunk.
//     const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];
//     PseudofragmentDB::ConnectionIterator chunk_it (database, chunk.first, chunk.second);
//     while (!chunk_it.AtEnd()) {
//       const Connection& connection = chunk_it.GetConnection();
//       // Initialize containers to store the compatible Pseudofragment IDs and weights.
//       QUERY_RESULT acyclic_qr;
//       QUERY_RESULT ring_qr;
//       QUERY_RESULT ring_part_qr;
//       // Iterate over the Connection's compatible Connections.
//       const CONNECTIONS_SET& compatible_connections = compatibilities[connection];
//       for (const Connection& compatible_connection : compatible_connections) {
//         // Iterate over the Pseudofragments that have the compatible Connection.
//         sqlite3_bind_int64(select_pseudofragments_with_connection, 1, compatible_connection.GetStartAtomType());
//         sqlite3_bind_int64(select_pseudofragments_with_connection, 2, compatible_connection.GetEndAtomType());
//         sqlite3_bind_int64(select_pseudofragments_with_connection, 3, compatible_connection.GetBondType());
//         sqlite3_result_code = sqlite3_step(select_pseudofragments_with_connection);
//         while (sqlite3_result_code == SQLITE_ROW) {
//           sqlite3_int64 pseudofragment_id = sqlite3_column_int64(select_pseudofragments_with_connection, 0);
//           bool is_ring_part = sqlite3_column_int(select_pseudofragments_with_connection, 1);
//           bool has_ring = sqlite3_column_int(select_pseudofragments_with_connection, 2);
//           unsigned pseudofragment_size = sqlite3_column_int(select_pseudofragments_with_connection, 3);
//           unsigned pseudofragment_frequency = sqlite3_column_int64(select_pseudofragments_with_connection, 4);
//           // Store the Pseudofragment's ID and weight in the appropiate containers.
//           if (!has_ring && !is_ring_part) {
//             float pseudofragment_weight = std::pow(pseudofragment_frequency, acyclic_frequency_gamma) * std::pow(pseudofragment_size, acyclic_size_gamma);
//             acyclic_qr.first.push_back(pseudofragment_id);
//             acyclic_qr.second.push_back(pseudofragment_weight);
//           };
//           if (has_ring) {
//             float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
//             ring_qr.first.push_back(pseudofragment_id);
//             ring_qr.second.push_back(pseudofragment_weight);
//           };
//           if (is_ring_part) {
//             float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
//             ring_part_qr.first.push_back(pseudofragment_id);
//             ring_part_qr.second.push_back(pseudofragment_weight);
//           };
//
//           // Sort the QUERY_RESULTs
//           // Unique the QUERY_RESULTs
//           // Insert the QUERY_RESULTs
//
//           sqlite3_result_code = sqlite3_step(select_pseudofragments_with_connection);
//         };
//         sqlite3_clear_bindings(select_pseudofragments_with_connection);
//         sqlite3_reset(select_pseudofragments_with_connection);
//       };
//       ++chunk_it;
//     };
//     chunk_it.Finalize();
//     sqlite3_finalize(select_pseudofragments_with_connection);
//   };
// };

void ConnectionQueryResults::CalculateQueryResults(const PseudofragmentDB& database) {
  // Ensure the containers are empty.
  acyclic.clear();
  ring.clear();
  ring_part.clear();
  // Initialize the containers with the appropiate structure.
  for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
    acyclic.insert({ connection, QUERY_RESULTS_BY_NUM{{1, QUERY_RESULT{}}} });
  };
  ring = acyclic;
  ring_part = acyclic;
  // Determine the number of OpenMP threads to use.
  int n_threads = DetermineOMPNThreads();
  omp_set_num_threads(n_threads);
  // Divide the input file into N equally sized chunks, where N is the number of
  // OpenMP threads.
  sqlite3_int64 n_pseudofragments = database.GetNPseudofragments();
  std::vector<std::pair<unsigned, unsigned>> database_chunks = EquallySizedChunks(n_pseudofragments, n_threads, true);
  // Initialize a progress bar.
  std::cout << "Calculating ConnectionQueryResults according to the lax compatibility definition..." << std::endl;
  boost::progress_display progress(n_pseudofragments);
  #pragma omp parallel for
  for (int thrid = 0; thrid < n_threads; ++thrid) {
    // Iterate over the Pseudofragments in the database chunk.
    const std::pair<unsigned, unsigned>& chunk = database_chunks[thrid];
    PseudofragmentDB::PseudofragmentIterator chunk_it (database, chunk.first, chunk.second);
    while (!chunk_it.AtEnd()) {
      sqlite3_int64 pseudofragment_id = chunk_it.GetPseudofragmentID();
      unsigned pseudofragment_frequency = chunk_it.GetPseudofragmentFrequency();
      const Pseudofragment& pseudofragment = chunk_it.GetPseudofragment();
      unsigned pseudofragment_size = pseudofragment.GetSize();
      const ConnectionsTable& pseudofragment_connections = pseudofragment.GetConnections();
      // Iterate over the Connections in the database.
      for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
        // If the Pseudofragment is compatible with the Connection record the
        // Pseudofragment's ID and weight in the appropiate containers.
        if (compatibilities.AreCompatible(pseudofragment_connections, connection)) {
          if (!pseudofragment.HasRing() && !pseudofragment.IsRingPart()) {
            QUERY_RESULT& query_result = acyclic[connection][1];
            float pseudofragment_weight = std::pow(pseudofragment_frequency, acyclic_frequency_gamma) * std::pow(pseudofragment_size, acyclic_size_gamma);
            #pragma omp critical
            {
              query_result.first.push_back(pseudofragment_id);
              query_result.second.push_back(pseudofragment_weight);
            };
          };
          if (pseudofragment.HasRing()) {
            QUERY_RESULT& query_result = ring[connection][1];
            float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
            #pragma omp critical
            {
              query_result.first.push_back(pseudofragment_id);
              query_result.second.push_back(pseudofragment_weight);
            };
          };
          if (pseudofragment.IsRingPart()) {
            QUERY_RESULT& query_result = ring_part[connection][1];
            float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
            #pragma omp critical
            {
              query_result.first.push_back(pseudofragment_id);
              query_result.second.push_back(pseudofragment_weight);
            };
          };
        };
      };
      ++chunk_it;
      #pragma omp critical
      {
        ++progress;
      };
    };
    chunk_it.Finalize();
  };
  // Since the threads may write to the query results asynchronously, sort the
  // query results.
  std::cout << "Sorting query results..." << std::endl;
  for (auto& [connection, qrbn] : acyclic) {
    QUERY_RESULT& query_result = qrbn[1];
    SortQueryResult(query_result);
  };
  for (auto& [connection, qrbn] : ring) {
    QUERY_RESULT& query_result = qrbn[1];
    SortQueryResult(query_result);
  };
  for (auto& [connection, qrbn] : ring_part) {
    QUERY_RESULT& query_result = qrbn[1];
    SortQueryResult(query_result);
  };
};

// void ConnectionQueryResults::CalculateQueryResults(const PseudofragmentDB& database) {
//   // Ensure the containers are empty.
//   acyclic.clear();
//   ring.clear();
//   ring_part.clear();
//   // Initialize a ConnectionQueryResults template.
//   CONNECTION_QUERY_RESULTS template_query_results;
//   for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
//     template_query_results.insert({ connection, QUERY_RESULTS_BY_NUM{{1, QUERY_RESULT{}}} });
//   };
//   // Initialize the member variables using the template.
//   acyclic = template_query_results;
//   ring = template_query_results;
//   ring_part = template_query_results;
//   // Initialize a progress bar.
//   std::cout << "Calculating ConnectionQueryResults according to the provided ConnectionCompatibilities (lax compatibility definition)..." << std::endl;
//   sqlite3_int64 n_pseudofragments = database.GetNPseudofragments();
//   boost::progress_display progress(n_pseudofragments);
//   // Iterate over the Pseudofragments in the database.
//   PseudofragmentDB::PseudofragmentIterator it (database);
//   while (!it.AtEnd()) {
//     sqlite3_int64 pseudofragment_id = it.GetPseudofragmentID();
//     unsigned pseudofragment_frequency = it.GetPseudofragmentFrequency();
//     const Pseudofragment& pseudofragment = it.GetPseudofragment();
//     unsigned pseudofragment_size = pseudofragment.GetSize();
//     const ConnectionsTable& pseudofragment_connections = pseudofragment.GetConnections();
//     // Iterate over the Connections in the database.
//     for (const auto& [connection, compatible_connections] : compatibilities.GetCompatibilityTable()) {
//       // If the Pseudofragment is compatible with the Connection record the
//       // Pseudofragment's ID and weight in the appropiate containers.
//       if (compatibilities.AreCompatible(pseudofragment_connections, connection)) {
//         if (!pseudofragment.HasRing() && !pseudofragment.IsRingPart()) {
//           QUERY_RESULT& query_result = acyclic[connection][1];
//           float pseudofragment_weight = std::pow(pseudofragment_frequency, acyclic_frequency_gamma) * std::pow(pseudofragment_size, acyclic_size_gamma);
//           query_result.first.push_back(pseudofragment_id);
//           query_result.second.push_back(pseudofragment_weight);
//         };
//         if (pseudofragment.HasRing()) {
//           QUERY_RESULT& query_result = ring[connection][1];
//           float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
//           query_result.first.push_back(pseudofragment_id);
//           query_result.second.push_back(pseudofragment_weight);
//         };
//         if (pseudofragment.IsRingPart()) {
//           QUERY_RESULT& query_result = ring_part[connection][1];
//           float pseudofragment_weight = std::pow(pseudofragment_frequency, ring_frequency_gamma) * std::pow(pseudofragment_size, ring_size_gamma);
//           query_result.first.push_back(pseudofragment_id);
//           query_result.second.push_back(pseudofragment_weight);
//         };
//       };
//     };
//     ++it;
//     ++progress;
//   };
//   it.Finalize();
// };

void ConnectionQueryResults::QueryResultIntersection(const QUERY_RESULT& qr1, const QUERY_RESULT& qr2, QUERY_RESULT& intersection, bool ignore_weights) {
  std::vector<unsigned>::const_iterator first_ids1 = qr1.first.begin();
  std::vector<unsigned>::const_iterator first_ids2 = qr2.first.begin();
  std::vector<unsigned>::const_iterator last_ids1 = qr1.first.end();
  std::vector<unsigned>::const_iterator last_ids2 = qr2.first.end();
  std::back_insert_iterator<std::vector<unsigned>> ids_intersection = std::back_inserter(intersection.first);

  unsigned qr1_size = qr1.first.size();
  unsigned qr2_size = qr2.first.size();
  unsigned max_intersection_size;
  if (qr1_size > qr2_size) {
    max_intersection_size = qr2_size;
  } else {
    max_intersection_size = qr1_size;
  };

  if (ignore_weights) {
    intersection.first.reserve(max_intersection_size);
    while (first_ids1 != last_ids1 && first_ids2 != last_ids2) {
      if (*first_ids1 < *first_ids2) {
        ++first_ids1;
      } else if (*first_ids2 < *first_ids1) {
        ++first_ids2;
      } else {
        *ids_intersection = *first_ids1;
        ++first_ids1;
        ++first_ids2;
        ++ids_intersection;
      };
    };
  } else {
    std::vector<float>::const_iterator first_weights1 = qr1.second.begin();
    std::vector<float>::const_iterator first_weights2 = qr2.second.begin();
    std::back_insert_iterator<std::vector<float>> weights_intersection = std::back_inserter(intersection.second);
    intersection.first.reserve(max_intersection_size);
    intersection.second.reserve(max_intersection_size);
    while (first_ids1 != last_ids1 && first_ids2 != last_ids2) {
      if (*first_ids1 < *first_ids2) {
        ++first_ids1;
        ++first_weights1;
      } else if (*first_ids2 < *first_ids1) {
        ++first_ids2;
        ++first_weights2;
      } else {
        *ids_intersection = *first_ids1;
        *weights_intersection = *first_weights1;
        ++first_ids1;
        ++first_weights1;
        ++first_ids2;
        ++first_weights2;
        ++ids_intersection;
        ++weights_intersection;
      };
    };
  };
};

QUERY_RESULT ConnectionQueryResults::SlidingIntersection(const CONNECTIONS_COMBINATION& combination, bool has_ring, bool is_ring_part, bool ignore_weights) {
  // Count the frequency of each Connection in the CONNECTIONS_COMBINATION.
  std::vector<Connection> connections;
  std::map<Connection, unsigned> connection_frequencies;
  connections.reserve(combination.size());
  CONNECTIONS_COMBINATION::const_iterator begin_it = combination.begin(), end_it = combination.end();
  for (const Connection& connection : combination) {
    connections.push_back(connection);
    unsigned frequency = std::count(begin_it, end_it, connection);
    connection_frequencies.insert({ connection, frequency });
  };

  // If there is only a single unique Connection retrieve the result immediately.
  int n_connections = connections.size();
  if (n_connections == 1) {
    if (has_ring) {
      return strict_ring[connections[0]][connection_frequencies[connections[0]]];
    } else if (is_ring_part) {
      return strict_ring_part[connections[0]][connection_frequencies[connections[0]]];
    } else {
      return strict_acyclic[connections[0]][connection_frequencies[connections[0]]];
    };
  };
  // Otherwise set up the multiple set intersection.
  std::vector<QUERY_RESULT*> qr_pointers;
  qr_pointers.reserve(n_connections);

  // Sort the Connections according to the number of the Pseudofragments
  // within the ConnectionQueryResults that are compatible with the
  // Connection atleast N times. Sorting them this way reduces the number
  // of required iterations by QueryResultIntersection.
  // Thereafter, retrieve the QUERY_RESULTs of the sorted Connections.
  if (has_ring) {
    std::sort(connections.begin(), connections.end(),
      [&](const Connection& c1, const Connection& c2) {
        return strict_ring[c1][connection_frequencies[c1]].first.size() < strict_ring[c2][connection_frequencies[c2]].first.size();
      });
    for (const Connection& connection : connections) {
      qr_pointers.push_back(&strict_ring[connection][connection_frequencies[connection]]);
    };
  } else if (is_ring_part) {
    std::sort(connections.begin(), connections.end(),
      [&](const Connection& c1, const Connection& c2) {
        return strict_ring_part[c1][connection_frequencies[c1]].first.size() < strict_ring_part[c2][connection_frequencies[c2]].first.size();
      });
    for (const Connection& connection : connections) {
      qr_pointers.push_back(&strict_ring_part[connection][connection_frequencies[connection]]);
    };
  } else {
    std::sort(connections.begin(), connections.end(),
      [&](const Connection& c1, const Connection& c2) {
        return strict_acyclic[c1][connection_frequencies[c1]].first.size() < strict_acyclic[c2][connection_frequencies[c2]].first.size();
      });
    for (const Connection& connection : connections) {
      qr_pointers.push_back(&strict_acyclic[connection][connection_frequencies[connection]]);
    };
  };

  // Calculate the multiple set intersection.
  QUERY_RESULT last_intersection, current_intersection;
  // Do a quick check to ensure that there is an intersection to be calculated.
  // Since the query results are sorted it would suffice with checking the size
  // of the first element of the vector.
  if (qr_pointers[0]->first.size() == 0) {
    return last_intersection;
  };
  QueryResultIntersection(*qr_pointers[0], *qr_pointers[1], last_intersection, ignore_weights);
  for (int i = 2; i < n_connections; ++i) {
    QueryResultIntersection(last_intersection, *qr_pointers[i], current_intersection, ignore_weights);
    last_intersection = std::move(current_intersection);
  };
  return last_intersection;
};

void ConnectionQueryResults::AssignFreshWeights(const Connection& connection, ConnectionPoint& cpoint) {
  // Assign the acyclic IDs and weights.
  QUERY_RESULT* query_result = &acyclic[connection][1];
  cpoint.ids.acyclic = query_result->first;
  cpoint.weights.acyclic = query_result->second;
  // Assign the cyclic IDs and weights.
  query_result = &ring[connection][1];
  cpoint.ids.ring = query_result->first;
  cpoint.weights.ring = query_result->second;
  // Assign the ring part IDs and weights.
  query_result = &ring_part[connection][1];
  cpoint.ids.ring_part = query_result->first;
  cpoint.weights.ring_part = query_result->second;
};

void ConnectionQueryResults::AssignFreshWeights(const Connection& connection, std::vector<ConnectionPoint>& cpoints) {
  QUERY_RESULT* acyclic_qresult = &acyclic[connection][1];
  QUERY_RESULT* ring_qresult = &ring[connection][1];
  QUERY_RESULT* ring_part_qresult = &ring_part[connection][1];
  for (ConnectionPoint& cpoint : cpoints) {
    // Assign the acyclic IDs and weights.
    cpoint.ids.acyclic = acyclic_qresult->first;
    cpoint.weights.acyclic = acyclic_qresult->second;
    // Assign the cyclic IDs and weights.
    cpoint.ids.ring = ring_qresult->first;
    cpoint.weights.ring = ring_qresult->second;
    // Assign the ring part IDs and weights.
    cpoint.ids.ring_part = ring_part_qresult->first;
    cpoint.weights.ring_part = ring_part_qresult->second;
  };
};

void ConnectionQueryResults::AssignFreshWeights(const Connection& connection, const std::vector<ConnectionPoint*>& cpoints) {
  QUERY_RESULT* acyclic_query_result = &acyclic[connection][1];
  QUERY_RESULT* ring_query_result = &ring[connection][1];
  QUERY_RESULT* ring_part_query_result = &ring_part[connection][1];
  for (ConnectionPoint* cpoint : cpoints) {
    // Assign the acyclic IDs and weights.
    cpoint->ids.acyclic = acyclic_query_result->first;
    cpoint->weights.acyclic = acyclic_query_result->second;
    // Assign the cyclic IDs and weights.
    cpoint->ids.ring = ring_query_result->first;
    cpoint->weights.ring = ring_query_result->second;
    // Assign the ring part IDs and weights.
    cpoint->ids.ring_part = ring_part_query_result->first;
    cpoint->weights.ring_part = ring_part_query_result->second;
  };
};

const CONNECTION_QUERY_RESULTS& ConnectionQueryResults::GetAcyclicResults() const {
  return acyclic;
};

const CONNECTION_QUERY_RESULTS& ConnectionQueryResults::GetRingResults() const {
  return ring;
};

const CONNECTION_QUERY_RESULTS& ConnectionQueryResults::GetRingPartResults() const {
  return ring_part;
};

const CONNECTION_QUERY_RESULTS& ConnectionQueryResults::GetStrictAcyclicResults() const {
  return strict_acyclic;
};

const CONNECTION_QUERY_RESULTS& ConnectionQueryResults::GetStrictRingResults() const {
  return strict_ring;
}

const CONNECTION_QUERY_RESULTS& ConnectionQueryResults::GetStrictRingPartResults() const {
  return strict_ring_part;
};

const ConnectionCompatibilities& ConnectionQueryResults::GetConnectionCompatibilities() const {
  return compatibilities;
};

float ConnectionQueryResults::GetAcyclicFrequencyGamma() const {
  return acyclic_frequency_gamma;
};

float ConnectionQueryResults::GetAcyclicSizeGamma() const {
  return acyclic_size_gamma;
};

float ConnectionQueryResults::GetRingFrequencyGamma() const {
  return ring_frequency_gamma;
};

float ConnectionQueryResults::GetRingSizeGamma() const {
  return ring_size_gamma;
};

void ConnectionQueryResults::Print() const {
  std::cout << "Acyclic Pseudofragments:" << std::endl;
  for (const auto& qrbn : acyclic) {
    std::cout << qrbn.first.GetString() << std::endl;
    for (const auto& qr : qrbn.second) {
      std::cout << qr.first << ": ";
      for (const auto& id : qr.second.first) {
        std::cout << id << " ";
      };
      std::cout << std::endl;
    };
  };
  std::cout << "Ring Pseudofragments:" << std::endl;
  for (const auto& qrbn : ring) {
    std::cout << qrbn.first.GetString() << std::endl;
    for (const auto& qr : qrbn.second) {
      std::cout << qr.first << ": ";
      for (const auto& id : qr.second.first) {
        std::cout << id << " ";
      };
      std::cout << std::endl;
    };
  };
  std::cout << "Ring part Pseudofragments:" << std::endl;
  for (const auto& qrbn : ring_part) {
    std::cout << qrbn.first.GetString() << std::endl;
    for (const auto& qr : qrbn.second) {
      std::cout << qr.first << ": ";
      for (const auto& id : qr.second.first) {
        std::cout << id << " ";
      };
      std::cout << std::endl;
    };
  };
};

void ConnectionQueryResults::PrintStrict() const {
  std::cout << "Acyclic Pseudofragments:" << std::endl;
  for (const auto& qrbn : strict_acyclic) {
    std::cout << qrbn.first.GetString() << std::endl;
    for (const auto& qr : qrbn.second) {
      std::cout << qr.first << ": ";
      for (const auto& id : qr.second.first) {
        std::cout << id << " ";
      };
      std::cout << std::endl;
    };
  };
  std::cout << "Ring Pseudofragments:" << std::endl;
  for (const auto& qrbn : strict_ring) {
    std::cout << qrbn.first.GetString() << std::endl;
    for (const auto& qr : qrbn.second) {
      std::cout << qr.first << ": ";
      for (const auto& id : qr.second.first) {
        std::cout << id << " ";
      };
      std::cout << std::endl;
    };
  };
  std::cout << "Ring part Pseudofragments:" << std::endl;
  for (const auto& qrbn : strict_ring_part) {
    std::cout << qrbn.first.GetString() << std::endl;
    for (const auto& qr : qrbn.second) {
      std::cout << qr.first << ": ";
      for (const auto& id : qr.second.first) {
        std::cout << id << " ";
      };
      std::cout << std::endl;
    };
  };
};

// Standalone functions to work with QUERY_RESULTs.
void ConcatenateQueryResults(QUERY_RESULT& recipient, const QUERY_RESULT& sender, bool ignore_weights) {
  if (ignore_weights) {
    recipient.first.insert(recipient.first.end(), sender.first.begin(), sender.first.end());
  } else {
    recipient.first.insert(recipient.first.end(), sender.first.begin(), sender.first.end());
    recipient.second.insert(recipient.second.end(), sender.second.begin(), sender.second.end());
  };
};

void SortQueryResult(QUERY_RESULT& query_result, bool ignore_weights) {
  std::vector<size_t> permutation = GetSortPermutation(query_result.first);
  ApplySortPermutation(query_result.first, permutation);
  if (!ignore_weights) {
    ApplySortPermutation(query_result.second, permutation);
  };
};

void UniqueQueryResult(QUERY_RESULT& query_result, bool ignore_weights) {
  std::vector<unsigned>::iterator first_ids = query_result.first.begin(), last_ids = query_result.first.end(), result_ids = query_result.first.begin();
  if (first_ids == last_ids) {
    return;
  };
  if (ignore_weights) {
    while (first_ids != last_ids) {
      if (*result_ids != *(++first_ids)) {
        *(++result_ids) = *first_ids;
      };
    };
    query_result.first.resize(std::distance(query_result.first.begin(), result_ids));
  } else {
    std::vector<float>::iterator first_weights = query_result.second.begin(), result_weights = query_result.second.begin();
    while (first_ids != last_ids) {
      ++first_weights;
      if (*result_ids != *(++first_ids)) {
        *(++result_ids) = *first_ids;
        *(++result_weights) = *first_weights;
      };
    };
    unsigned n_elements = std::distance(query_result.first.begin(), result_ids);
    query_result.first.resize(n_elements);
    query_result.second.resize(n_elements);
  };
};
