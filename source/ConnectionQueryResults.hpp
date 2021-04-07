// To be compiled with "-ffast-math". This boosts exponentiation performance.

#pragma once
#ifndef _CONNECTION_QUERY_RESULTS_HPP_
#define _CONNECTION_QUERY_RESULTS_HPP_

#include <vector>
#include <math.h>
#include <sqlite3.h>
#include <boost/serialization/vector.hpp>
#include <boost/progress.hpp>
#include "Pseudofragment.hpp"

// Declaration of the version integer of the ConnectionQueryResults class.
extern const unsigned connection_query_results_version;

// Definitions of STL containers to store query results.
typedef std::pair<std::vector<unsigned>, std::vector<float>> QUERY_RESULT;
typedef std::unordered_map<unsigned, QUERY_RESULT> QUERY_RESULTS_BY_NUM;
typedef std::unordered_map<Connection, QUERY_RESULTS_BY_NUM> CONNECTION_QUERY_RESULTS;

// Definition of a type to store an (ordered) combination of Connections.
typedef std::vector<Connection> CONNECTIONS_COMBINATION;

// Class to calculate and store the results of SQLite3 connection queries. These
// queries retrieve the IDs and weights of the Pseudofragments compatible with
// each Connection.
class ConnectionQueryResults {
  // Containers to store the results. The results are calculated twice according to
  // a strict and lax compatibility definition.
  CONNECTION_QUERY_RESULTS acyclic, ring, ring_part, strict_acyclic, strict_ring, strict_ring_part;
  // Associative array describing which Connections are compatible with each other.
  ConnectionCompatibilities compatibilities;
  // The weight exponents used to calculate the Pseudofragment weights.
  float acyclic_frequency_gamma = 0.0f, acyclic_level_gamma = 0.0f, ring_frequency_gamma = 0.0f, ring_level_gamma = 0.0f;

public:
  ConnectionQueryResults();
  ConnectionQueryResults(const ConnectionCompatibilities& compatibilities, sqlite3* database, float acyclic_frequency_gamma = 1, float acyclic_level_gamma = 1, float ring_frequency_gamma = 0.5, float ring_level_gamma = 0);

  // Function to calculate the query results according to an arbitrary compatibility
  // definition (intended to be a lax one).
  void CalculateQueryResults(sqlite3* database);
  // Function to calculate the query results according to the strictest definition.
  void CalculateStrictQueryResults(sqlite3* database);

  // Function to calculate the multiple query results intersection associated with
  // a given combination of Connections.
  QUERY_RESULT SlidingIntersection(const CONNECTIONS_COMBINATION& combination, bool has_ring, bool is_ring_part, bool ignore_weights = false);

  // Functions to assign a copy of the original query results as IDs/weights to a
  // ConnectionPoint.
  void AssignFreshWeights(const Connection& connection, ConnectionPoint& cpoint);
  void AssignFreshWeights(const Connection& connection, ConnectionPoint* cpoint);
  void AssignFreshWeights(const Connection& connection, const std::vector<ConnectionPoint*>& cpoints);
  void AssignFreshWeights(const Connection& connection, std::vector<ConnectionPoint>& cpoints);

  const CONNECTION_QUERY_RESULTS& GetAcyclicResults() const;
  const CONNECTION_QUERY_RESULTS& GetRingResults() const;
  const CONNECTION_QUERY_RESULTS& GetRingPartResults() const;
  const CONNECTION_QUERY_RESULTS& GetStrictAcyclicResults() const;
  const CONNECTION_QUERY_RESULTS& GetStrictRingResults() const;
  const CONNECTION_QUERY_RESULTS& GetStrictRingPartResults() const;
  const ConnectionCompatibilities& GetConnectionCompatibilities() const;
  float GetAcyclicFrequencyGamma() const;
  float GetAcyclicLevelGamma() const;
  float GetRingFrequencyGamma() const;
  float GetRingLevelGamma() const;
  void Print() const;
  void PrintStrict() const;

private:
  // Function to calculate the intersection of two query results.
  void QueryResultIntersection(const QUERY_RESULT& qr1, const QUERY_RESULT& qr2, QUERY_RESULT& intersection, bool ignore_weights = false);

  template <class Archive>
  void serialize(Archive& ar, const unsigned version = connection_query_results_version) {
    ar & acyclic & ring & ring_part & strict_acyclic & strict_ring & strict_ring_part & compatibilities & acyclic_frequency_gamma & acyclic_level_gamma & ring_frequency_gamma & ring_level_gamma;
  };

  friend class boost::serialization::access;
  friend class ConnectionPoint;
  friend class MolBrick;
  friend class ReconstructedMol;
};

// Function to concatenate two query results.
void ConcatenateQueryResults(QUERY_RESULT& recipient, const QUERY_RESULT& sender, bool ignore_weights = false);

// Function to remove duplicate elements from a query result.
void UniqueQueryResult(QUERY_RESULT& query_result, bool ignore_weights = false);

#endif
