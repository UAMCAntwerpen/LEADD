#include "Connection.hpp"

// Definition of the version integer of the Connection & associated classes.
extern const unsigned connection_version = 20210615;

std::uint64_t EncodeAtomTypes(std::uint32_t start_atom_type, std::uint32_t end_atom_type) {
  std::uint64_t encoded_atom_types = ((std::uint64_t) start_atom_type) << 32 | end_atom_type;
  return encoded_atom_types;
};

std::pair<std::uint32_t, std::uint32_t> DecodeAtomTypes(std::uint64_t encoded_atom_types) {
  std::uint32_t start_atom_type = (std::uint32_t) (encoded_atom_types >> 32);
  std::uint32_t end_atom_type = (std::uint32_t) encoded_atom_types;
  return std::pair<std::uint32_t, std::uint32_t> (start_atom_type, end_atom_type);
};

// Class Connection
Connection::Connection() = default;
Connection::Connection(std::uint32_t sat, std::uint32_t eat, std::uint32_t bt) :
  start_atom_type(sat),
  end_atom_type(eat),
  bond_type(bt) {};

bool Connection::operator==(const Connection& other) const {
  return (start_atom_type == other.start_atom_type &&
    end_atom_type == other.end_atom_type &&
    bond_type == other.bond_type);
};

bool Connection::operator<(const Connection& other) const {
  if (start_atom_type == other.start_atom_type) {
    if (end_atom_type == other.end_atom_type) {
      return bond_type < other.bond_type;
    } else {
      return end_atom_type < other.end_atom_type;
    };
  } else {
    return start_atom_type < other.start_atom_type;
  };
};

Connection Connection::Mirror() const {
  return Connection(end_atom_type, start_atom_type, bond_type);
};

std::uint32_t Connection::GetStartAtomType() const {
  return start_atom_type;
};

std::uint32_t Connection::GetEndAtomType() const {
  return end_atom_type;
};

std::uint32_t Connection::GetBondType() const {
  return bond_type;
};

std::uint64_t Connection::GetEncodedAtomTypes() const {
  return EncodeAtomTypes(start_atom_type, end_atom_type);
};

std::string Connection::GetString() const {
  std::stringstream ss;
  ss << "(" << start_atom_type << "," << end_atom_type << "," << bond_type << ")";
  return ss.str();
};

void Connection::Print() const {
  std::cout << GetString() << std::endl;
};
