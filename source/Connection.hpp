#pragma once
#ifndef _CONNECTION_HPP_
#define _CONNECTION_HPP_

#include <iostream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// Declaration of the version integer of the Connection & associated classes.
extern const unsigned connection_version;

std::uint64_t EncodeAtomTypes(std::uint32_t start_atom_type, std::uint32_t end_atom_type);
std::pair<std::uint32_t, std::uint32_t> DecodeAtomTypes(std::uint64_t encoded_atom_types);

// Class to represent the link between two atoms in a molecule.
class Connection {
  // A connection is essentially 3 integers bundled together.
  //  (1) start_atom_type: The atom type of the atom on which the
  //      Connection is centered.
  //  (2) end_atom_type: The atom type of the atom to which the above
  //      used to be connected in the source molecule from which the
  //      Connection was extracted.
  //  (3) bond_type: The bond type of the bond connecting the aforementioned atoms.
  std::uint32_t start_atom_type = 0, end_atom_type = 0, bond_type = 0;

public:
  Connection();
  Connection(std::uint32_t sat, std::uint32_t eat, std::uint32_t bt);

  bool operator==(const Connection& other) const;
  bool operator<(const Connection& other) const;

  // Function to generate a copy of the Connection in which the
  // start_atom_type and end_atom_type are swapped.
  Connection Mirror() const;

  std::uint32_t GetStartAtomType() const;
  std::uint32_t GetEndAtomType() const;
  std::uint32_t GetBondType() const;
  std::uint64_t GetEncodedAtomTypes() const;
  std::string GetString() const;

  void Print() const;

private:
  template <class Archive>
  void serialize(Archive& ar, const unsigned version = connection_version) {
    ar & start_atom_type & end_atom_type & bond_type;
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
    std::size_t seed = 0;
    boost::hash_combine(seed, c.start_atom_type);
    boost::hash_combine(seed, c.end_atom_type);
    boost::hash_combine(seed, c.bond_type);
    return seed;
  };
};

#endif
