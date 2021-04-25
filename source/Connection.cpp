#include "Connection.hpp"

// Definition of the version integer of the Connection & associated classes.
extern const unsigned connection_version = 20200306;

// Class Connection
Connection::Connection() = default;
Connection::Connection(unsigned sat, unsigned eat, unsigned bt) :
  start_atom_type(sat),
  end_atom_type(eat),
  bond_type(bt) {
  Encode();
};
Connection::Connection(unsigned sat, unsigned eat, unsigned bt, boost::format & fmt) :
  start_atom_type(sat),
  end_atom_type(eat),
  bond_type(bt),
  string(str(fmt % sat % eat % bt)) {
  Encode();
};

bool Connection::operator==(const Connection & other) const {
  return (start_atom_type == other.start_atom_type &&
    end_atom_type == other.end_atom_type &&
    bond_type == other.bond_type);
};

bool Connection::operator<(const Connection & other) const {
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

void Connection::GenerateString(boost::format & fmt) {
  string = str(fmt % start_atom_type % end_atom_type % bond_type);
};

void Connection::Encode() {
  encoded = static_cast<std::uint8_t>(start_atom_type) | (static_cast<std::uint8_t>(end_atom_type) << 8);
};

std::pair<std::uint8_t, std::uint8_t> Connection::Decode() const {
  std::pair<std::uint8_t, std::uint8_t> p(encoded & 0xFF, (encoded >> 8) & 0xFF);
  return p;
};

unsigned Connection::GetStartAtomType() const {
  return start_atom_type;
};

unsigned Connection::GetEndAtomType() const {
  return end_atom_type;
};

unsigned Connection::GetBondType() const {
  return bond_type;
};

unsigned Connection::GetEncoded() const {
  return static_cast<unsigned>(encoded);
};

const std::string& Connection::GetString() const {
  return string;
};

Connection Connection::Mirror() const {
  return Connection(end_atom_type, start_atom_type, bond_type);
};

Connection Connection::Mirror(boost::format & fmt) const {
  return Connection(end_atom_type, start_atom_type, bond_type, fmt);
};


// Connection decoding function to be used on standalone integers.
std::pair<std::uint8_t, std::uint8_t> DecodeConnection(unsigned encoded) {
  std::pair<std::uint8_t, std::uint8_t> p(encoded & 0xFF, (encoded >> 8) & 0xFF);
  return p;
};


// Class ConnectionCompatibilities
ConnectionCompatibilities::ConnectionCompatibilities() = default;
ConnectionCompatibilities::ConnectionCompatibilities(const CONNECTIONS_SET & connections, unsigned stringency) : stringency(stringency) {
  CONNECTIONS_SET::const_iterator it, begin_it = connections.begin(), end_it = connections.end();
  unsigned start_atom_type, end_atom_type, bond_type;

  // If the stringency level is 0 all Connections with the same bond type are compatible.
  // This corresponds mostly to a valence model.
  if (stringency == 0) {
    for (const Connection& connection : connections) {
      CONNECTIONS_SET compatible_connections;
      for (const Connection& query_connection : connections) {
        if (connection.GetBondType() == query_connection.GetBondType()) {
          compatible_connections.insert(query_connection);
        };
      };
      compatibility_table.insert({ connection, compatible_connections });
    };

  // If the stringency level is 1 a Connection is considered compatible if:
  //  (1) The Connection's starting atom type has been previously observed
  //      paired with the query Connection's starting atom type.
  //  (2) It shares the query Connection's BondType.
  } else if (stringency == 1) {
    std::unordered_map<unsigned, std::unordered_set<unsigned>> observed_pairings;
    std::unordered_map<unsigned, std::unordered_set<unsigned>>::const_iterator p_it;
    // Store the observed atom pairings.
    for (const Connection& connection : connections) {
      start_atom_type = connection.GetStartAtomType();
      end_atom_type = connection.GetEndAtomType();
      p_it = observed_pairings.find(start_atom_type);
      if (p_it == observed_pairings.end()) {
        observed_pairings.insert({ start_atom_type, std::unordered_set<unsigned> {end_atom_type} });
      } else {
        observed_pairings[start_atom_type].insert(end_atom_type);
      };
      p_it = observed_pairings.find(end_atom_type);
      if (p_it == observed_pairings.end()) {
        observed_pairings.insert({ end_atom_type, std::unordered_set<unsigned> {start_atom_type} });
      } else {
        observed_pairings[end_atom_type].insert(start_atom_type);
      };
    };
    // Define a functor that evaluates if a Connection is compatible with the
    // last-defined Connection in the scope based on the latter's start_atom_type
    // and bond_type.
    std::function<bool(const Connection&)> IsCompatible =
      [&](const Connection& c) {
      bool valid_pairing = false, same_bond_type = false;
      const std::unordered_set<unsigned>& partners = observed_pairings[start_atom_type];
      if (partners.find(c.GetStartAtomType()) != partners.end()) {
        valid_pairing = true;
      };
      if (bond_type == c.GetBondType()) {
        same_bond_type = true;
      };
      return valid_pairing && same_bond_type;
    };
    // Use the observed atom pairings and the query Connection's BondType to
    // determine which Connections are compatible with it.
    for (const Connection& connection : connections) {
      CONNECTIONS_SET compatible_connections;
      start_atom_type = connection.GetStartAtomType();
      bond_type = connection.GetBondType();
      it = std::find_if(begin_it, end_it, IsCompatible);
      while (it != end_it) {
        compatible_connections.insert(*it);
        ++it;
        it = std::find_if(it, end_it, IsCompatible);
      };
      compatibility_table.insert({ connection, compatible_connections });
    };

  // If the stringency level is 3 only the mirrored Connection is considered
  // to be compatible.
} else if (stringency == 2) {
    boost::format formatter("(%d, %d, %d)");
    for (const Connection& connection : connections) {
      CONNECTIONS_SET compatible_connections{ connection.Mirror(formatter) };
      compatibility_table.insert({ connection, compatible_connections });
    };

  // Verify that the provided stringency level is valid.
  } else {
    throw std::runtime_error("Stringency criterion has to be one of 0, 1 or 2");
  };
};

CONNECTIONS_SET& ConnectionCompatibilities::operator[](const Connection& connection) {
  return compatibility_table[connection];
};

const CONNECTIONS_SET& ConnectionCompatibilities::at(const Connection & connection) const {
  return compatibility_table.at(connection);
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
