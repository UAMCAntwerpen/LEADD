#pragma once
#ifndef _PSEUDOFRAGMENT_HPP_
#define _PSEUDOFRAGMENT_HPP_

#include <vector>
#include <random>
#include <algorithm>
#include <GraphMol/GraphMol.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/MolPickler.h>
#include <boost/serialization/vector.hpp>
// #include <boost/serialization/map.hpp>
#include "Connection.hpp"

class ConnectionPoint;

// Definition & forward declaration of associative arrays to convert
// RDKit::BondTypes to equivalent integers. Ideally one would use a
// getter to retrieve this information from RDKit::Bond::BondType,
// but this doesn't seem to be implemented.
typedef std::map<unsigned, RDKit::Bond::BondType> INT_BOND_TYPE_TABLE;
typedef std::map<RDKit::Bond::BondType, unsigned> BOND_TYPE_INT_TABLE;
extern INT_BOND_TYPE_TABLE int_bond_type_table;
extern BOND_TYPE_INT_TABLE bond_type_int_table;

// Definitions of STL containers of Connections. Provided for convenience.
typedef std::vector<ConnectionPoint> CONNECTION_POINT_VECTOR;
typedef std::vector<const ConnectionPoint*> CONNECTION_POINT_PTR_VECTOR;
// typedef std::map<Connection, CONNECTION_POINT_VECTOR> CONNECTIONS_MAP;
typedef std::unordered_map<Connection, CONNECTION_POINT_VECTOR> CONNECTIONS_MAP;
// Definition & declation of an associative array containing precalculated uniform
// distributions of given sizes.
typedef std::unordered_map<unsigned, std::uniform_int_distribution<unsigned>> INT_DISTRIBUTION_MAP;
extern INT_DISTRIBUTION_MAP uniform_int_distributions;
extern const unsigned max_size_uniform_int_distributions;

// Version integer of the Pseudofragment & associated classes.
extern const unsigned pseudofragment_version;

// Class to store the (database) IDs of Pseudofragments, separated into the
// different categories of Pseudofragments.
class PseudofragmentIDs {
  std::vector<unsigned> acyclic, ring, ring_part;

public:
  PseudofragmentIDs();

  void AddAcyclicID(unsigned id);
  void AddRingID(unsigned id);
  void AddRingPartID(unsigned id);
  void ClearAcyclicIDs();
  void ClearRingIDs();
  void ClearRingPartIDs();
  void Clear();

  const std::vector<unsigned>& GetAcyclicIDs() const;
  const std::vector<unsigned>& GetRingIDs() const;
  const std::vector<unsigned>& GetRingPartIDs() const;

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = pseudofragment_version) {
    ar & acyclic & ring & ring_part;
  };

  friend class boost::serialization::access;
  friend class MolBrick;
  friend class ConnectionQueryResults;
  friend class EvolutionGuide;
};

// Class to store the weights of Pseudofragments, separated into the
// different categories of Pseudofragments. This container is meant
// to be associated with an equal-dimensions PseudofragmentIDs object.
class PseudofragmentWeights {
  std::vector<float> acyclic, ring, ring_part;

public:
  PseudofragmentWeights();

  void AddAcyclicWeight(float weight);
  void AddRingWeight(float weight);
  void AddRingPartWeight(float weight);
  void ClearAcyclicWeights();
  void ClearRingWeights();
  void ClearRingPartWeights();
  void Clear();

  const std::vector<float>& GetAcyclicWeights() const;
  const std::vector<float>& GetRingWeights() const;
  const std::vector<float>& GetRingPartWeights() const;

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = pseudofragment_version) {
    ar & acyclic & ring & ring_part;
  };

  friend class boost::serialization::access;
  friend class MolBrick;
  friend class ConnectionQueryResults;
  friend class EvolutionGuide;
};

// Class to store a Connection object centered on a specific atom
// of a molecule. Hence, it's a specific instance of a Connection.
class ConnectionPoint {
  // Immutable index of the atom on which it's centered.
  unsigned atom_idx = 0;
  // Containers to store the IDs and corresponding weights of
  // Pseudofragments that are compatible with the Connection.
  PseudofragmentIDs ids;
  PseudofragmentWeights weights;

public:
  ConnectionPoint();
  ConnectionPoint(unsigned atom_idx);

  bool operator==(const ConnectionPoint& other) const;

  unsigned GetAtomIdx() const;
  const PseudofragmentIDs& GetIDs() const;
  const PseudofragmentWeights& GetWeights() const;
  void Print() const;

private:
  void ClearWeights();

  template <class Archive>
  void serialize(Archive& ar, const unsigned version = pseudofragment_version) {
    ar & atom_idx & ids & weights;
  };

  friend class ConnectionsTable;
  friend class Pseudofragment;
  friend class ConnectionQueryResults;
  friend class MolBrick;
  friend class ReconstructedMol;
  friend class EvolutionGuide;
  friend class boost::serialization::access;
};

// Class to store the ConnectionPoints of a Pseudofragment.
class ConnectionsTable {
  // Associative array storing all ConnectionPoint instances of a
  // specific Connection.
  CONNECTIONS_MAP cmap;
  unsigned n_connections = 0, n_connection_points = 0;

public:
  ConnectionsTable();

  CONNECTIONS_MAP::iterator begin();
  CONNECTIONS_MAP::const_iterator begin() const;
  CONNECTIONS_MAP::iterator end();
  CONNECTIONS_MAP::const_iterator end() const;
  CONNECTIONS_MAP::iterator find(const Connection& connection);
  CONNECTIONS_MAP::const_iterator find(const Connection& connection) const;

  // Operators to access elements in the associative array.
  CONNECTION_POINT_VECTOR& operator[](const Connection& connection);
  const CONNECTION_POINT_VECTOR& at(const Connection& connection) const;

  // Functions to add/remove ConnectionPoints from the ConnectionsTable.
  ConnectionPoint* AddConnection(const Connection& connection, const ConnectionPoint& cpoint);
  ConnectionPoint* AddConnection(const Connection& connection, const ConnectionPoint* cpoint);
  ConnectionPoint* AddConnection(const Connection& connection, const unsigned atom_idx);
  std::vector<ConnectionPoint*> AddConnections(const Connection& connection, CONNECTION_POINT_VECTOR& cpoints);
  std::vector<ConnectionPoint*> AddConnections(const Connection& connection, const CONNECTION_POINT_PTR_VECTOR& cpoints);
  std::vector<ConnectionPoint*> AddConnections(const Connection& connection, const std::vector<unsigned>& atom_indices);
  bool RemoveConnection(const Connection& connection, const ConnectionPoint& cpoint);
  bool RemoveConnection(const Connection& connection, const ConnectionPoint* cpoint);
  bool RemoveConnection(const Connection& connection, const unsigned atom_idx);
  bool RemoveConnection(const Connection& connection);
  bool RemoveConnections(const Connection& connection, const CONNECTION_POINT_VECTOR& cpoints);
  bool RemoveConnections(const Connection& connection, const CONNECTION_POINT_PTR_VECTOR& cpoints);
  bool RemoveConnections(const Connection& connection, const std::vector<unsigned>& atom_indices);
  // Function to transfer a ConnectionPoint from one ConnectionsTable to another.
  ConnectionPoint* MoveConnectionTo(const Connection& connection, const unsigned atom_idx, ConnectionsTable& recipient);

  // Functions to check if the ConnectionsTable contains a Connection/specific ConnectionPoint.
  bool HasConnection(const Connection& connection) const;
  bool HasConnection(const Connection& connection, const unsigned atom_idx) const;

  // Functions to retrieve pointers to specific ConnectionPoints.
  ConnectionPoint* GetConnectionPointWithIdx(const Connection& connection, unsigned atom_idx, std::mt19937& prng);
  std::vector<ConnectionPoint*> GetConnectionPointsWithIdx(const Connection& connection, unsigned atom_idx);
  std::vector<ConnectionPoint*> GetConnectionPointsWithIndices(const Connection& connection, std::vector<unsigned> atom_indices);
  std::vector<ConnectionPoint*> GetConnectionPoints(const Connection& connection, const std::vector<ConnectionPoint>& cpoints);

  // Function to check if any of the contained Connections are compatible with a query Connection.
  bool IsCompatibleWith(const Connection& connection, ConnectionCompatibilities& compatibilities) const;
  // Function to retrieve a random contained Connection compatible with a query Connection.
  Connection GetRandomCompatibleConnection(const Connection& connection, ConnectionCompatibilities& compatibilities, std::mt19937& prng) const;

  bool Empty() const;
  unsigned Size() const;
  void Print() const;
  void ClearWeights();

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = pseudofragment_version) {
    ar & cmap & n_connections & n_connection_points;
  };

  friend class boost::serialization::access;
  friend class ReconstructedMol;
};

// Class to store a molecular subgraph alongside its ConnectionPoints as a fragment.
class Pseudofragment {
protected:
  // Molecule object containing the molecular graph of the fragment.
  // Note that a Pseudofragment is a subgraph of a real molecular graph and
  // therefore is not guaranteed to respect valence rules.
  RDKit::RWMol pseudomol;
  // Fragment's ConnectionsTable.
  ConnectionsTable connections;
  // SMILES string encoding the molecular graph alongside the ConnectionsTable
  // as pseudoatoms.
  std::string smiles;
  // Boolean indicating whether the fragment contains a ring system.
  bool has_ring = false;
  // Boolean indicating whether the fragment was extracted from a ring system.
  bool ring_part = false;
  // Integer expressing the size of the fragment (in number of heavy atoms)
  unsigned level = 0;

public:
  Pseudofragment();
  Pseudofragment(const RDKit::ROMol& pmol, const ConnectionsTable& cs, const std::string& sml, bool hr, bool rp);
  Pseudofragment(const RDKit::RWMol& pmol, const ConnectionsTable& cs, const std::string& sml, bool hr, bool rp);

  bool operator== (const Pseudofragment& other) const;

  // Function to create a chemically-sane equivalent of the Pseudofragment's
  // molecular graph by fixing valence issues through explicit hydrogen addition.
  RDKit::RWMol SanitizeCopy() const;

  // Function to generate a ECFP fingerprint of the (sanitized) molecular graph.
  const RDKit::SparseIntVect<std::uint32_t>* GetFingerprint() const;

  const RDKit::ROMol& GetPseudomol() const;
  const ConnectionsTable& GetConnections() const;
  const std::string& GetSMILES() const;
  bool HasRing() const;
  bool IsRingPart() const;
  unsigned GetLevel() const;

  // Serialization protocol. Note that the protocol is implemented with the
  // boost::serialization library but built around the RDKit's MolPickler
  // serializer. Therefore the behaviour of the serializer and deserializer
  // are different (i.e. split).
  template<class Archive>
  void save(Archive& ar, const unsigned version = pseudofragment_version) const {
    std::string pickle;
    RDKit::MolPickler::pickleMol(pseudomol, pickle);
    ar << pickle;
    ar << connections;
    ar << smiles;
    ar << has_ring;
    ar << ring_part;
    ar << level;
  };

  template<class Archive>
  void load(Archive& ar, const unsigned version = pseudofragment_version) {
    std::string pickle;
    ar >> pickle;
    RDKit::MolPickler::molFromPickle(pickle, pseudomol);
    ar >> connections;
    ar >> smiles;
    ar >> has_ring;
    ar >> ring_part;
    ar >> level;
  };

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = pseudofragment_version) {
    boost::serialization::split_member(ar, *this, pseudofragment_version);
  };

  friend class boost::serialization::access;
  friend class MolBrick;
};

#endif
