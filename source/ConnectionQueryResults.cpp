#include "ConnectionQueryResults.hpp"

// Definition of the version integer of the ConnectionQueryResults class.
extern const unsigned connection_query_results_version = 20200306;

// Class ConnectionQueryResults
ConnectionQueryResults::ConnectionQueryResults() = default;
ConnectionQueryResults::ConnectionQueryResults(const ConnectionCompatibilities & compatibilities, sqlite3 * database, float acyclic_frequency_gamma, float acyclic_level_gamma, float ring_frequency_gamma, float ring_level_gamma) :
  compatibilities(compatibilities),
  acyclic_frequency_gamma(acyclic_frequency_gamma),
  acyclic_level_gamma(acyclic_level_gamma),
  ring_frequency_gamma(ring_frequency_gamma),
  ring_level_gamma(ring_level_gamma) {
  CalculateStrictQueryResults(database);
  CalculateQueryResults(database);
};

void ConnectionQueryResults::CalculateQueryResults(sqlite3* database) {
  assert(acyclic.empty());
  assert(ring.empty());
  assert(ring_part.empty());
  // Prepare the required SQLite3 statements.
  sqlite3_stmt* select_n_pseudofragments;
  sqlite3_stmt* select_pseudofragments;
  sqlite3_prepare_v2(database, "SELECT id, ring_part, has_ring, level, pickle, pickle_size, frequency FROM pseudofragments", -1, &select_pseudofragments, NULL);
  sqlite3_prepare_v2(database, "SELECT COUNT(*) FROM pseudofragments", -1, &select_n_pseudofragments, NULL);
  // Determine the number of Pseudofragments in the database.
  sqlite3_step(select_n_pseudofragments);
  unsigned n_pseudofragments = sqlite3_column_int64(select_n_pseudofragments, 0);
  // Initialize a progress bar.
  boost::progress_display progress(n_pseudofragments);
  // Initialize a ConnectionQueryResults template.
  CONNECTION_QUERY_RESULTS template_query_results;
  for (const auto& c : compatibilities.GetCompatibilityTable()) {
    const Connection& connection = c.first;
    template_query_results.insert({ connection, QUERY_RESULTS_BY_NUM{{1, QUERY_RESULT{}}} });
  };
  // Assign to the member variables the structure of the template.
  acyclic = template_query_results;
  ring = template_query_results;
  ring_part = template_query_results;
  // Loop over the Pseudofragments in the database and retrieve their data.
  std::stringstream blob;
  int rc = sqlite3_step(select_pseudofragments);
  while (rc == SQLITE_ROW) {
    int id = sqlite3_column_int64(select_pseudofragments, 0);
    int rpart = sqlite3_column_int(select_pseudofragments, 1);
    int has_ring = sqlite3_column_int(select_pseudofragments, 2);
    int level = sqlite3_column_int(select_pseudofragments, 3);
    int blob_size = sqlite3_column_int(select_pseudofragments, 5);
    int frequency = sqlite3_column_int(select_pseudofragments, 6);
    // Deserialize the Pseudofragment.
    Pseudofragment pseudofragment;
    blob.write((const char*)sqlite3_column_blob(select_pseudofragments, 4), blob_size);
    boost::archive::binary_iarchive archive(blob);
    pseudofragment.load(archive);
    // Loop over the Connections in the database and check if the Pseudofragment
    // is compatible with them.
    const ConnectionsTable& pseudofragment_connections = pseudofragment.GetConnections();
    for (const auto& c : template_query_results) {
      const Connection& connection = c.first;
      if (pseudofragment_connections.IsCompatibleWith(connection, compatibilities)) {
        if (rpart) {
          QUERY_RESULT& query_result = ring_part[connection][1];
          query_result.first.push_back(id);
          query_result.second.push_back(std::pow(frequency, acyclic_frequency_gamma) * std::pow(level, acyclic_level_gamma));
        } else if (has_ring) {
          QUERY_RESULT& query_result = ring[connection][1];
          query_result.first.push_back(id);
          query_result.second.push_back(std::pow(frequency, ring_frequency_gamma) * std::pow(level, ring_level_gamma));
        } else {
          QUERY_RESULT& query_result = acyclic[connection][1];
          query_result.first.push_back(id);
          query_result.second.push_back(std::pow(frequency, acyclic_frequency_gamma) * std::pow(level, acyclic_level_gamma));
        };
      };
    };
    // Clear the blob and Pseudofragment and move on to the next fragment.
    blob.str(std::string());
    blob.clear();
    ++progress;
    rc = sqlite3_step(select_pseudofragments);
  };
  // Destroy the SQLite3 statements to avoid memory leaks.
  sqlite3_finalize(select_n_pseudofragments);
  sqlite3_finalize(select_pseudofragments);
};

void ConnectionQueryResults::CalculateStrictQueryResults(sqlite3* database) {
  assert(strict_acyclic.empty());
  assert(strict_ring.empty());
  assert(strict_ring_part.empty());
  // Set up a Connection string formatter.
  boost::format formatter("(%d, %d, %d)");
  // Prepare the required SQLite3 statements.
  sqlite3_stmt* select_connection_id;
  sqlite3_stmt* select_max_connection_freq;
  sqlite3_stmt* select_pseudofragments_with_connection;
  sqlite3_prepare_v2(database, "SELECT id FROM connections WHERE start_atom_type = ? AND end_atom_type = ? AND bond_type = ?", -1, &select_connection_id, NULL);
  sqlite3_prepare_v2(database, "SELECT MAX(frequency) FROM pseudofragments_connections WHERE connection_id = ?", -1, &select_max_connection_freq, NULL);
  sqlite3_prepare_v2(database, "SELECT id, ring_part, has_ring, level, frequency FROM pseudofragments WHERE id IN (SELECT pseudofragment_id FROM pseudofragments_connections WHERE connection_id = ? AND frequency >= ?)", -1, &select_pseudofragments_with_connection, NULL);
  // Initialize empty containers to store the query results.
  QUERY_RESULT acyclic_qresult, ring_qresult, ring_part_qresult;
  std::pair<CONNECTION_QUERY_RESULTS::iterator, bool> it_bool_pair;
  // Loop over the databases Connections.
  for (const auto& c : compatibilities.GetCompatibilityTable()) {
    const Connection& connection = c.first;
    // Mirror the Connection. The mirrored Connection is the only
    // one considered compatible with the original connection.
    Connection mirrored_connection = connection.Mirror(formatter);
    // Create empty entries in the ConnectionQueryResults corresponding
    // to the current Connection.
    it_bool_pair = strict_acyclic.emplace(connection, QUERY_RESULTS_BY_NUM());
    QUERY_RESULTS_BY_NUM& acyclic_qresults = it_bool_pair.first->second;
    it_bool_pair = strict_ring.emplace(connection, QUERY_RESULTS_BY_NUM());
    QUERY_RESULTS_BY_NUM& ring_qresults = it_bool_pair.first->second;
    it_bool_pair = strict_ring_part.emplace(connection, QUERY_RESULTS_BY_NUM());
    QUERY_RESULTS_BY_NUM& ring_part_qresults = it_bool_pair.first->second;
    // Retrieve the IDs of the mirrored Connection from the database.
    sqlite3_bind_int(select_connection_id, 1, mirrored_connection.GetStartAtomType());
    sqlite3_bind_int(select_connection_id, 2, mirrored_connection.GetEndAtomType());
    sqlite3_bind_int(select_connection_id, 3, mirrored_connection.GetBondType());
    sqlite3_step(select_connection_id);
    unsigned connection_id = sqlite3_column_int64(select_connection_id, 0);
    // Retrieve the maximum number of occurrences of the mirrored
    // Connection in a single Pseudofragment.
    sqlite3_bind_int64(select_max_connection_freq, 1, connection_id);
    sqlite3_step(select_max_connection_freq);
    unsigned max_connection_freq = sqlite3_column_int(select_max_connection_freq, 0);
    // Loop over the range of the number of mirrored Connection occurrences.
    for (unsigned freq = 1; freq <= max_connection_freq; ++freq) {
      // Retrieve the Pseudofragments that have the mirrored Connection
      // atleast N number of times.
      sqlite3_bind_int64(select_pseudofragments_with_connection, 1, connection_id);
      sqlite3_bind_int(select_pseudofragments_with_connection, 2, freq);
      int rc = sqlite3_step(select_pseudofragments_with_connection);
      while (rc == SQLITE_ROW) {
        unsigned id = sqlite3_column_int64(select_pseudofragments_with_connection, 0);
        bool is_ring_part = sqlite3_column_int(select_pseudofragments_with_connection, 1);
        bool has_ring = sqlite3_column_int(select_pseudofragments_with_connection, 2);
        unsigned level = sqlite3_column_int(select_pseudofragments_with_connection, 3);
        unsigned frequency = sqlite3_column_int(select_pseudofragments_with_connection, 4);
        // Clasify the Pseudofragment depending on whether it's cyclic or not.
        if (has_ring) {
          ring_qresult.first.push_back(id);
          ring_qresult.second.push_back(std::pow(frequency, ring_frequency_gamma) * std::pow(level, ring_level_gamma));
        } else if (is_ring_part) {
          ring_part_qresult.first.push_back(id);
          ring_part_qresult.second.push_back(std::pow(frequency, acyclic_frequency_gamma) * std::pow(level, acyclic_level_gamma));
        } else {
          acyclic_qresult.first.push_back(id);
          acyclic_qresult.second.push_back(std::pow(frequency, acyclic_frequency_gamma) * std::pow(level, acyclic_level_gamma));
        };
        rc = sqlite3_step(select_pseudofragments_with_connection);
      };
      // Store the query results.
      acyclic_qresults.insert({ freq, std::move(acyclic_qresult) });
      ring_qresults.insert({ freq, std::move(ring_qresult) });
      ring_part_qresults.insert({ freq, std::move(ring_part_qresult) });
      // Clear the containers to make space for the next query results.
      acyclic_qresult.first.clear();
      acyclic_qresult.second.clear();
      ring_qresult.first.clear();
      ring_qresult.second.clear();
      ring_part_qresult.first.clear();
      ring_part_qresult.second.clear();
      // Reset the Pseudofragment select statement.
      sqlite3_clear_bindings(select_pseudofragments_with_connection);
      sqlite3_reset(select_pseudofragments_with_connection);
    };
    // Reset the Connection select statements.
    sqlite3_clear_bindings(select_connection_id);
    sqlite3_clear_bindings(select_max_connection_freq);
    sqlite3_reset(select_connection_id);
    sqlite3_reset(select_max_connection_freq);
  };
  // Destruct the SQLite3 statements to avoid memory leaks.
  sqlite3_finalize(select_connection_id);
  sqlite3_finalize(select_max_connection_freq);
  sqlite3_finalize(select_pseudofragments_with_connection);
};

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

void ConnectionQueryResults::AssignFreshWeights(const Connection& connection, ConnectionPoint* cpoint) {
  // Assign the acyclic IDs and weights.
  QUERY_RESULT* query_result = &acyclic[connection][1];
  cpoint->ids.acyclic = query_result->first;
  cpoint->weights.acyclic = query_result->second;
  // Assign the cyclic IDs and weights.
  query_result = &ring[connection][1];
  cpoint->ids.ring = query_result->first;
  cpoint->weights.ring = query_result->second;
  // Assign the ring part IDs and weights.
  query_result = &ring_part[connection][1];
  cpoint->ids.ring_part = query_result->first;
  cpoint->weights.ring_part = query_result->second;
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

float ConnectionQueryResults::GetAcyclicLevelGamma() const {
  return acyclic_level_gamma;
};

float ConnectionQueryResults::GetRingFrequencyGamma() const {
  return ring_frequency_gamma;
};

float ConnectionQueryResults::GetRingLevelGamma() const {
  return ring_level_gamma;
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

