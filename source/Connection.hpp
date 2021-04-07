#pragma once
#ifndef _CONNECTION_HPP_
#define _CONNECTION_HPP_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <boost/format.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/unordered_map.hpp>

class Connection;

// Declaration of the version integer of the Connection & associated classes.
extern const unsigned connection_version;

// Definitions of STL containers of Connections. Provided for convenience.
typedef std::vector<Connection> CONNECTIONS_VECTOR;
typedef std::unordered_set<Connection> CONNECTIONS_SET;
typedef std::unordered_map<Connection, CONNECTIONS_SET> COMPATIBILITY_TABLE;

// Class to represent the link between two atoms in a molecule.
// A connection is essentially 3 integers bundled together.
class Connection {
  // Three (8-bit) integers encoding:
  //  (1) start_atom_type: The atom type of the atom on which the
  //      Connections is centered.
  //  (2) end_atom_type: The atom type of the atom to which the above
  //      used to be connected in the source molecule from which the
  //      Connection was extracted.
  unsigned start_atom_type = 0, end_atom_type = 0, bond_type = 0;
  // A 16-bit integer encoding the start_atom_type and end_atom_type
  // as two separate 8-bit integers.
  std::uint16_t encoded = 0;
  // A string representation of the the triplet of integers.
  std::string string;

public:
  Connection();
  Connection(unsigned sat, unsigned eat, unsigned bt);
  Connection(unsigned sat, unsigned eat, unsigned bt, boost::format& fmt);

  bool operator==(const Connection& other) const;
  bool operator<(const Connection& other) const;

  // Generates a formatted string describing the Connection. Useful for
  // visualization.
  void GenerateString(boost::format& fmt);

  // Function to convert the start_atom_type and end_atom_type into
  // a single 16-bit integer.
  void Encode();

  // Function to decode the atom types encoded as a 16-bit integer
  // into two separate 8-bit integers.
  std::pair<std::uint8_t, std::uint8_t> Decode() const;

  unsigned GetStartAtomType() const;
  unsigned GetEndAtomType() const;
  unsigned GetBondType() const;
  unsigned GetEncoded() const;
  const std::string& GetString() const;

  // Functions to generate a copy of the Connection in which the
  // start_atom_type and end_atom_type are swapped.
  Connection Mirror() const;
  Connection Mirror(boost::format& fmt) const;

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = connection_version) {
    ar & start_atom_type & end_atom_type & bond_type & encoded & string;
  };

  friend struct std::hash<Connection>;
  friend class Pseudofragment;
  friend class MolBrick;
  friend class ReconstructedMol;
  friend class boost::serialization::access;
};

// Standard library hash function specialization for Connection objects.
// Necessary to use Connections as keys in "unordered" containers.
template<>
struct std::hash<Connection> {
  std::size_t operator() (const Connection& c) const {
    return ((std::hash<std::uint8_t>()(c.start_atom_type) ^
      (std::hash<std::uint8_t>()(c.end_atom_type) << 1)) >> 1) ^
      (std::hash<std::uint8_t>()(c.bond_type) << 1);
  };
};

// Standalone function to decode atom type encoding 16-bit  integers into two
// 8-bit integers.
std::pair<std::uint8_t, std::uint8_t> DecodeConnection(unsigned encoded);

// Class to store which Connections are compatible with each other,
// according to a compatibility definition.
class ConnectionCompatibilities {
  // Associate container to store Connection-Connection compatibility
  // relationships.
  COMPATIBILITY_TABLE compatibility_table;
  // Integer symbolizing the compatibility definition.
  unsigned stringency = 0;

public:
  ConnectionCompatibilities();
  ConnectionCompatibilities(const CONNECTIONS_SET& connections, unsigned stringency = 1);

  // Operators to access elements in the compatibility table.
  CONNECTIONS_SET& operator[](const Connection& connection);
  const CONNECTIONS_SET& at(const Connection& connection) const;

  bool HasConnection(const Connection& connection) const;
  const COMPATIBILITY_TABLE& GetCompatibilityTable() const;
  void Print() const;

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = connection_version) {
    ar & compatibility_table & stringency;
  };

  friend class boost::serialization::access;
  friend class ConnectionsTable;
};

#endif
