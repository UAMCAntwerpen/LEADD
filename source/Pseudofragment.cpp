#include "Pseudofragment.hpp"

extern const unsigned pseudofragment_version = 20210615;

// Definition of the integer-RDKit::BondType associative arrays.
INT_BOND_TYPE_TABLE int_bond_type_table({
  { 0,  RDKit::Bond::UNSPECIFIED}, { 1,        RDKit::Bond::SINGLE},
  { 2,       RDKit::Bond::DOUBLE}, { 3,        RDKit::Bond::TRIPLE},
  { 4,    RDKit::Bond::QUADRUPLE}, { 5,     RDKit::Bond::QUINTUPLE},
  { 6,     RDKit::Bond::HEXTUPLE}, { 7,   RDKit::Bond::ONEANDAHALF},
  { 8,  RDKit::Bond::TWOANDAHALF}, { 9, RDKit::Bond::THREEANDAHALF},
  {10, RDKit::Bond::FOURANDAHALF}, {11,  RDKit::Bond::FIVEANDAHALF},
  {12,     RDKit::Bond::AROMATIC}, {13,         RDKit::Bond::IONIC},
  {14,     RDKit::Bond::HYDROGEN}, {15,   RDKit::Bond::THREECENTER},
  {16,    RDKit::Bond::DATIVEONE}, {17,        RDKit::Bond::DATIVE},
  {18,      RDKit::Bond::DATIVEL}, {19,       RDKit::Bond::DATIVER},
  {20,        RDKit::Bond::OTHER}, {21,          RDKit::Bond::ZERO}
});

BOND_TYPE_INT_TABLE bond_type_int_table({
  {RDKit::Bond::UNSPECIFIED,   0}, {RDKit::Bond::SINGLE,        1},
  {RDKit::Bond::DOUBLE,        2}, {RDKit::Bond::TRIPLE,        3},
  {RDKit::Bond::QUADRUPLE,     4}, {RDKit::Bond::QUINTUPLE,     5},
  {RDKit::Bond::HEXTUPLE,      6}, {RDKit::Bond::ONEANDAHALF,   7},
  {RDKit::Bond::TWOANDAHALF,   8}, {RDKit::Bond::THREEANDAHALF, 9},
  {RDKit::Bond::FOURANDAHALF, 10}, {RDKit::Bond::FIVEANDAHALF, 11},
  {RDKit::Bond::AROMATIC,     12}, {RDKit::Bond::IONIC,        13},
  {RDKit::Bond::HYDROGEN,     14}, {RDKit::Bond::THREECENTER,  15},
  {RDKit::Bond::DATIVEONE,    16}, {RDKit::Bond::DATIVE,       17},
  {RDKit::Bond::DATIVEL,      18}, {RDKit::Bond::DATIVER,      19},
  {RDKit::Bond::OTHER,        20}, {RDKit::Bond::ZERO,         21}
});

BOND_TYPE_INT_TABLE bond_type_order_table ({
  {RDKit::Bond::UNSPECIFIED,   1}, {RDKit::Bond::SINGLE,        1},
  {RDKit::Bond::DOUBLE,        2}, {RDKit::Bond::TRIPLE,        3},
  {RDKit::Bond::QUADRUPLE,     4}, {RDKit::Bond::QUINTUPLE,     5},
  {RDKit::Bond::HEXTUPLE,      6}, {RDKit::Bond::ONEANDAHALF,   1},
  {RDKit::Bond::TWOANDAHALF,   2}, {RDKit::Bond::THREEANDAHALF, 3},
  {RDKit::Bond::FOURANDAHALF,  4}, {RDKit::Bond::FIVEANDAHALF,  5},
  {RDKit::Bond::AROMATIC,      1}, {RDKit::Bond::IONIC,         1},
  {RDKit::Bond::HYDROGEN,      1}, {RDKit::Bond::THREECENTER,   1},
  {RDKit::Bond::DATIVEONE,     1}, {RDKit::Bond::DATIVE,        1},
  {RDKit::Bond::DATIVEL,       1}, {RDKit::Bond::DATIVER,       1},
  {RDKit::Bond::OTHER,         1}, {RDKit::Bond::ZERO,          0}
});

// Pre-construction of uniform distributions of given sizes.
INT_DISTRIBUTION_MAP uniform_int_distributions{
  {2,  std::uniform_int_distribution<unsigned>(0, 1)},
  {3,  std::uniform_int_distribution<unsigned>(0, 2)},
  {4,  std::uniform_int_distribution<unsigned>(0, 3)},
  {5,  std::uniform_int_distribution<unsigned>(0, 4)},
  {6,  std::uniform_int_distribution<unsigned>(0, 5)},
  {7,  std::uniform_int_distribution<unsigned>(0, 6)},
  {8,  std::uniform_int_distribution<unsigned>(0, 7)},
  {9,  std::uniform_int_distribution<unsigned>(0, 8)},
  {10, std::uniform_int_distribution<unsigned>(0, 9)},
  {11, std::uniform_int_distribution<unsigned>(0, 10)},
  {12, std::uniform_int_distribution<unsigned>(0, 11)},
  {13, std::uniform_int_distribution<unsigned>(0, 12)},
  {14, std::uniform_int_distribution<unsigned>(0, 13)},
  {15, std::uniform_int_distribution<unsigned>(0, 14)}
};

extern const unsigned max_size_uniform_int_distributions = 15;


// Class PseudofragmentIDs
const std::vector<unsigned>& PseudofragmentIDs::GetAcyclicIDs() const {
  return acyclic;
};

const std::vector<unsigned>& PseudofragmentIDs::GetRingIDs() const {
  return ring;
};

const std::vector<unsigned>& PseudofragmentIDs::GetRingPartIDs() const {
  return ring_part;
};

void PseudofragmentIDs::AddAcyclicID(unsigned id) {
  acyclic.push_back(id);
};

void PseudofragmentIDs::AddRingID(unsigned id) {
  ring.push_back(id);
};

void PseudofragmentIDs::AddRingPartID(unsigned id) {
  ring_part.push_back(id);
};

void PseudofragmentIDs::ClearAcyclicIDs() {
  acyclic.clear();
};
void PseudofragmentIDs::ClearRingIDs() {
  ring.clear();
};

void PseudofragmentIDs::ClearRingPartIDs() {
  ring_part.clear();
};

void PseudofragmentIDs::Clear() {
  ClearAcyclicIDs();
  ClearRingIDs();
  ClearRingPartIDs();
};


// Class PseudofragmentWeights
const std::vector<float>& PseudofragmentWeights::GetAcyclicWeights() const {
  return acyclic;
};

const std::vector<float>& PseudofragmentWeights::GetRingWeights() const {
  return ring;
};

const std::vector<float>& PseudofragmentWeights::GetRingPartWeights() const {
  return ring_part;
};

void PseudofragmentWeights::AddAcyclicWeight(float weight) {
  acyclic.push_back(weight);
};

void PseudofragmentWeights::AddRingWeight(float weight) {
  ring.push_back(weight);
};

void PseudofragmentWeights::AddRingPartWeight(float weight) {
  ring_part.push_back(weight);
};

void PseudofragmentWeights::ClearAcyclicWeights() {
  acyclic.clear();
};

void PseudofragmentWeights::ClearRingWeights() {
  ring.clear();
};

void PseudofragmentWeights::ClearRingPartWeights() {
  ring_part.clear();
};

void PseudofragmentWeights::Clear() {
  ClearAcyclicWeights();
  ClearRingWeights();
  ClearRingPartWeights();
};


// Class ConnectionPoint
ConnectionPoint::ConnectionPoint() = default;
ConnectionPoint::ConnectionPoint(unsigned atom_idx) : atom_idx(atom_idx) {};

bool ConnectionPoint::operator==(const ConnectionPoint& other) const {
  return atom_idx == other.atom_idx;
};

unsigned ConnectionPoint::GetAtomIdx() const {
  return atom_idx;
};

const PseudofragmentIDs& ConnectionPoint::GetIDs() const {
  return ids;
};

const PseudofragmentWeights& ConnectionPoint::GetWeights() const {
  return weights;
};

void ConnectionPoint::Print() const {
  std::cout << "Atom idx: " << atom_idx << std::endl;
  std::cout << "Acyclic IDs: ";
  for (unsigned id : ids.GetAcyclicIDs()) {
    std::cout << id << " ";
  };
  std::cout << std::endl << "Acyclic weights: ";
  for (float weight : weights.GetAcyclicWeights()) {
    std::cout << weight << " ";
  };
  std::cout << std::endl << "Ring IDs: ";
  for (unsigned id : ids.GetRingIDs()) {
    std::cout << id << " ";
  };
  std::cout << std::endl << "Ring weights: ";
  for (float weight : weights.GetRingWeights()) {
    std::cout << weight << " ";
  };
  std::cout << std::endl << "Ring part IDs: ";
  for (unsigned id : ids.GetRingPartIDs()) {
    std::cout << id << " ";
  };
  std::cout << std::endl << "Ring part weights: ";
  for (float weight : weights.GetRingPartWeights()) {
    std::cout << weight << " ";
  };
  std::cout << std::endl;
};

void ConnectionPoint::ClearWeights() {
  ids.Clear();
  weights.Clear();
};


// Class ConnectionsTable
CONNECTIONS_MAP::iterator ConnectionsTable::begin() {
  return cmap.begin();
};

CONNECTIONS_MAP::const_iterator ConnectionsTable::begin() const {
  return cmap.begin();
};

CONNECTIONS_MAP::iterator ConnectionsTable::end() {
  return cmap.end();
};

CONNECTIONS_MAP::const_iterator ConnectionsTable::end() const {
  return cmap.end();
};

CONNECTIONS_MAP::iterator ConnectionsTable::find(const Connection& connection) {
  return cmap.find(connection);
};

CONNECTIONS_MAP::const_iterator ConnectionsTable::find(const Connection& connection) const {
  return cmap.find(connection);
};

CONNECTION_POINT_VECTOR& ConnectionsTable::operator[](const Connection& connection) {
  return cmap[connection];
};

const CONNECTION_POINT_VECTOR& ConnectionsTable::at(const Connection& connection) const {
  return cmap.at(connection);
};

ConnectionPoint* ConnectionsTable::AddConnection(const Connection& connection, const ConnectionPoint& cpoint) {
  CONNECTIONS_MAP::iterator it = cmap.find(connection);
  if (it != cmap.end()) {
    CONNECTION_POINT_VECTOR& cpoints = it->second;
    cpoints.push_back(cpoint);
    ++n_connection_points;
    return &cpoints.back();
  } else {
    auto[it, inserted] = cmap.insert({connection, {cpoint}});
    CONNECTION_POINT_VECTOR& cpoints = it->second;
    ++n_connections;
    ++n_connection_points;
    return &cpoints.back();
  };
};

ConnectionPoint* ConnectionsTable::AddConnection(const Connection& connection, const unsigned atom_idx) {
  ConnectionPoint cpoint (atom_idx);
  CONNECTIONS_MAP::iterator it = cmap.find(connection);
  if (it != cmap.end()) {
    CONNECTION_POINT_VECTOR& cpoints = it->second;
    cpoints.push_back(cpoint);
    ++n_connection_points;
    return &cpoints.back();
  } else {
    auto[it, inserted] = cmap.insert({connection, {cpoint}});
    CONNECTION_POINT_VECTOR& cpoints = it->second;
    ++n_connections;
    ++n_connection_points;
    return &cpoints.back();
  };
};

std::vector<ConnectionPoint*> ConnectionsTable::AddConnections(const Connection& connection, CONNECTION_POINT_VECTOR& cpoints) {
  std::vector<ConnectionPoint*> inserted_cpoints;
  unsigned n_cpoints = cpoints.size();
  inserted_cpoints.reserve(n_cpoints);
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it != cmap.end()) {
    CONNECTION_POINT_VECTOR& cpoints_vec = cmap_it->second;
    cpoints_vec.reserve(cpoints_vec.size() + n_cpoints);
    for (ConnectionPoint& cpoint : cpoints) {
      cpoints_vec.push_back(cpoint);
      inserted_cpoints.push_back(&cpoints_vec.back());
    };
    n_connection_points += n_cpoints;
  } else {
    std::pair<CONNECTIONS_MAP::iterator, bool> it_b_p = cmap.insert(std::make_pair(connection, cpoints));
    ++n_connections;
    n_connection_points += n_cpoints;
    for (ConnectionPoint& cpoint : it_b_p.first->second) {
      inserted_cpoints.push_back(&cpoint);
    };
  };
  return inserted_cpoints;
};

std::vector<ConnectionPoint*> ConnectionsTable::AddConnections(const Connection& connection, const CONNECTION_POINT_PTR_VECTOR& cpoints) {
  std::vector<ConnectionPoint*> inserted_cpoints;
  unsigned n_cpoints = cpoints.size();
  inserted_cpoints.reserve(n_cpoints);
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it != cmap.end()) {
    CONNECTION_POINT_VECTOR& cpoint_vec = cmap_it->second;
    cpoint_vec.reserve(cpoint_vec.size() + n_cpoints);
    for (const ConnectionPoint* cpoint : cpoints) {
      cpoint_vec.push_back(*cpoint);
      inserted_cpoints.push_back(&cpoint_vec.back());
    };
    n_connection_points += n_cpoints;
  } else {
    std::pair<CONNECTIONS_MAP::iterator, bool> it_b_p = cmap.emplace(connection, CONNECTION_POINT_VECTOR());
    CONNECTION_POINT_VECTOR& cpoint_vec = it_b_p.first->second;
    for (const ConnectionPoint* cpoint : cpoints) {
      cpoint_vec.push_back(*cpoint);
      inserted_cpoints.push_back(&cpoint_vec.back());
    };
    ++n_connections;
    n_connection_points += n_cpoints;
  };
  return inserted_cpoints;
};

std::vector<ConnectionPoint*> ConnectionsTable::AddConnections(const Connection& connection, const std::vector<unsigned>& atom_indices) {
  std::vector<ConnectionPoint*> inserted_cpoints;
  unsigned n = atom_indices.size();
  inserted_cpoints.reserve(n);
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it != cmap.end()) {
    CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
    cpoints.reserve(cpoints.size() + n);
    for (unsigned atom_idx : atom_indices) {
      cpoints.push_back(ConnectionPoint(atom_idx));
      inserted_cpoints.push_back(&cpoints.back());
    };
    n_connection_points += n;
  } else {
    std::pair<CONNECTIONS_MAP::iterator, bool> it_b_p = cmap.emplace(connection, CONNECTION_POINT_VECTOR());
    CONNECTION_POINT_VECTOR& cpoints = it_b_p.first->second;
    cpoints.reserve(n);
    for (unsigned atom_idx : atom_indices) {
      cpoints.push_back(ConnectionPoint(atom_idx));
      inserted_cpoints.push_back(&cpoints.back());
    };
    ++n_connections;
    n_connection_points += n;
  };
  return inserted_cpoints;
};

bool ConnectionsTable::RemoveConnection(const Connection& connection, const ConnectionPoint& cpoint) {
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    throw std::runtime_error("Target Connection for deletion wasn't found.");
    return false;
  } else {
    CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
    CONNECTION_POINT_VECTOR::iterator vec_it, vec_end_it = cpoints.end();
    vec_it = std::find(cpoints.begin(), vec_end_it, cpoint);
    if (vec_it == vec_end_it) {
      throw std::runtime_error("Target ConnectionPoint for deletion wasn't found.");
      return false;
    } else {
      if (cpoints.size() > 1) {
        if (vec_it == vec_end_it - 1) {
          cpoints.pop_back();
        } else {
          *vec_it = std::move(cpoints.back());
          cpoints.pop_back();
        };
      } else {
        cmap.erase(cmap_it);
        --n_connections;
      };
      --n_connection_points;
      return true;
    };
  };
};

bool ConnectionsTable::RemoveConnection(const Connection& connection, const unsigned atom_idx) {
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    throw std::runtime_error("Target Connection for deletion wasn't found.");
    return false;
  } else {
    CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
    CONNECTION_POINT_VECTOR::iterator vec_it, vec_end_it = cpoints.end();
    vec_it = std::find_if(cpoints.begin(), vec_end_it,
      [=](const ConnectionPoint& cpoint) {
        return cpoint.atom_idx == atom_idx;
      });
    if (vec_it == vec_end_it) {
      throw std::runtime_error("Target ConnectionPoint for deletion wasn't found.");
      return false;
    } else {
      if (cpoints.size() > 1) {
        if (vec_it == vec_end_it - 1) {
          cpoints.pop_back();
        } else {
          *vec_it = std::move(cpoints.back());
          cpoints.pop_back();
        };
      } else {
        cmap.erase(cmap_it);
        --n_connections;
      };
      --n_connection_points;
      return true;
    };
  };
};

bool ConnectionsTable::RemoveConnection(const Connection& connection) {
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    throw std::runtime_error("Target Connection for deletion wasn't found.");
    return false;
  } else {
    unsigned n_connections_vec = cmap_it->second.size();
    cmap.erase(cmap_it);
    --n_connections;
    n_connection_points -= n_connections_vec;
    return true;
  };
};

bool ConnectionsTable::RemoveConnections(const Connection& connection, const CONNECTION_POINT_VECTOR& cpoints) {
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    throw std::runtime_error("Target Connection for deletion wasn't found.");
    return false;
  } else {
    CONNECTION_POINT_VECTOR& cpoint_vec = cmap_it->second;
    CONNECTION_POINT_VECTOR::iterator vec_it, vec_end_it, vec_begin_it = cpoint_vec.begin();
    for (const ConnectionPoint& cpoint : cpoints) {
      vec_end_it = cpoint_vec.end();
      vec_it = std::find(vec_begin_it, vec_end_it, cpoint);
      if (vec_it == vec_end_it) {
        throw std::runtime_error("Target ConnectionPoint for deletion wasn't found.");
        return false;
      } else {
        if (cpoint_vec.size() > 1) {
          if (vec_it == vec_end_it - 1) {
            cpoint_vec.pop_back();
          } else {
            *vec_it = std::move(cpoint_vec.back());
            cpoint_vec.pop_back();
          };
        } else {
          cmap.erase(cmap_it);
          --n_connections;
        };
        --n_connection_points;
      };
    };
    return true;
  };
};

bool ConnectionsTable::RemoveConnections(const Connection& connection, const std::vector<unsigned>& atom_indices) {
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    throw std::runtime_error("Target Connection for deletion wasn't found.");
    return false;
  } else {
    CONNECTION_POINT_VECTOR& cpoint_vec = cmap_it->second;
    CONNECTION_POINT_VECTOR::iterator vec_it, vec_end_it, vec_begin_it = cpoint_vec.begin();
    for (const unsigned atom_idx : atom_indices) {
      vec_end_it = cpoint_vec.end();
      vec_it = std::find_if(vec_begin_it, vec_end_it,
        [=](const ConnectionPoint& cpoint) {
          return cpoint.atom_idx == atom_idx;
        });
      if (vec_it == vec_end_it) {
        throw std::runtime_error("Target ConnectionPoint for deletion wasn't found.");
        return false;
      } else {
        if (cpoint_vec.size() > 1) {
          if (vec_it == vec_end_it - 1) {
            cpoint_vec.pop_back();
          } else {
            *vec_it = std::move(cpoint_vec.back());
            cpoint_vec.pop_back();
          };
        } else {
          cmap.erase(cmap_it);
          --n_connections;
        };
        --n_connection_points;
      };
    };
    return true;
  };
};

ConnectionPoint* ConnectionsTable::MoveConnectionTo(const Connection& connection, const unsigned atom_idx, ConnectionsTable& recipient) {
  // Search for the Connection in both the source and recipient ConnectionTables.
  CONNECTIONS_MAP::iterator sender_it = cmap.find(connection), recipient_it = recipient.cmap.find(connection);
  // If the Connection isn't present in the source table signal function failure.
  if (sender_it == cmap.end()) {
    throw std::runtime_error("Target Connection to move wasn't found.");
  };
  // Search for the ConnectionPoint in the source ConnectionTable's Connection entry.
  CONNECTION_POINT_VECTOR& sender_cpoint_vec = sender_it->second;
  CONNECTION_POINT_VECTOR::iterator vec_it, vec_begin_it = sender_cpoint_vec.begin(), vec_end_it = sender_cpoint_vec.end();
  vec_it = std::find_if(vec_begin_it, vec_end_it,
    [=](const ConnectionPoint& cpoint) {
      return cpoint.atom_idx == atom_idx;
    });
  // If the ConnectionPoint isn't present in the Connection entry signal function failure.
  if (vec_it == vec_end_it) {
    throw std::runtime_error("Target ConnectionPoint to move wasn't found.");
  };
  // Add the Connection to the recipient ConnectionTable.
  ConnectionPoint* inserted_cpoint;
  if (recipient_it == recipient.cmap.end()) {
    std::pair<CONNECTIONS_MAP::iterator, bool> it_b_p = recipient.cmap.emplace(connection, CONNECTION_POINT_VECTOR());
    CONNECTION_POINT_VECTOR& cpoints = it_b_p.first->second;
    cpoints.push_back(std::move(*vec_it));
    inserted_cpoint = &cpoints.back();
    ++recipient.n_connections;
  } else {
    CONNECTION_POINT_VECTOR& cpoints = recipient_it->second;
    cpoints.push_back(std::move(*vec_it));
    inserted_cpoint = &cpoints.back();
  };
  ++recipient.n_connection_points;
  // Remove the Connection from the recipient ConnectionTable.
  if (sender_cpoint_vec.size() > 1) {
    if (vec_it == vec_end_it - 1) {
      sender_cpoint_vec.pop_back();
    } else {
      *vec_it = std::move(sender_cpoint_vec.back());
      sender_cpoint_vec.pop_back();
    };
  } else {
    cmap.erase(sender_it);
    --n_connections;
  };
  --n_connection_points;
  return inserted_cpoint;
};

bool ConnectionsTable::HasConnection(const Connection& connection) const {
  CONNECTIONS_MAP::const_iterator cmap_it = cmap.find(connection);
  return cmap_it != cmap.end();
};

bool ConnectionsTable::HasConnection(const Connection& connection, const unsigned atom_idx) const {
  CONNECTIONS_MAP::const_iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    return false;
  };
  const CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
  for (const ConnectionPoint& cpoint : cpoints) {
    if (cpoint.atom_idx == atom_idx) {
      return true;
    };
  };
  return false;
};

ConnectionPoint* ConnectionsTable::GetConnectionPointWithIdx(const Connection& connection, unsigned atom_idx, std::mt19937& prng) {
  // Retrieve all matching ConnectionPoints.
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    throw std::runtime_error("Connection couldn't be found.");
  };
  CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
  std::vector<ConnectionPoint*> matches;
  for (ConnectionPoint& cpoint : cpoints) {
    if (cpoint.atom_idx == atom_idx) {
      matches.push_back(&cpoint);
    };
  };
  // Select a random matching ConnectionPoint.
  unsigned n_cpoints = matches.size();
  if (n_cpoints == 0) {
    throw std::runtime_error("ConnectionPoint couldn't be found.");
  };
  if (n_cpoints > 1) {
    // It isn't necessary to check that the distribution exists since we will
    // never have more than 15 ConnectionPoints on the same atom.
    return matches[uniform_int_distributions[matches.size()](prng)];
  } else {
    return matches[0];
  };
};

std::vector<ConnectionPoint*> ConnectionsTable::GetConnectionPointsWithIdx(const Connection& connection, unsigned atom_idx) {
  std::vector<ConnectionPoint*> matches;
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    return matches;
  };
  CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
  for (ConnectionPoint& cpoint : cpoints) {
    if (cpoint.atom_idx == atom_idx) {
      matches.push_back(&cpoint);
    };
  };
  return matches;
};

std::vector<ConnectionPoint*> ConnectionsTable::GetConnectionPointsWithIndices(const Connection& connection, std::vector<unsigned> atom_indices) {
  std::vector<ConnectionPoint*> matches;
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    return matches;
  };
  CONNECTION_POINT_VECTOR& cpoints = cmap_it->second;
  std::sort(atom_indices.begin(), atom_indices.end());
  atom_indices.erase(std::unique(atom_indices.begin(), atom_indices.end()), atom_indices.end());
  for (unsigned atom_idx : atom_indices) {
    for (ConnectionPoint& cpoint : cpoints) {
      if (cpoint.atom_idx == atom_idx) {
        matches.push_back(&cpoint);
      };
    };
  };
  return matches;
};

std::vector<ConnectionPoint*> ConnectionsTable::GetConnectionPoints(const Connection& connection, const std::vector<ConnectionPoint>& cpoints) {
  std::vector<ConnectionPoint*> matches;
  CONNECTIONS_MAP::iterator cmap_it = cmap.find(connection);
  if (cmap_it == cmap.end()) {
    return matches;
  };
  CONNECTION_POINT_VECTOR& map_cpoints = cmap_it->second;
  std::vector<unsigned> atom_indices;
  atom_indices.reserve(cpoints.size());
  for (const ConnectionPoint& cpoint : cpoints) {
    atom_indices.push_back(cpoint.atom_idx);
  };
  std::sort(atom_indices.begin(), atom_indices.end());
  atom_indices.erase(std::unique(atom_indices.begin(), atom_indices.end()), atom_indices.end());
  for (unsigned atom_idx : atom_indices) {
    for (ConnectionPoint& map_cpoint : map_cpoints) {
      if (map_cpoint.atom_idx == atom_idx) {
        matches.push_back(&map_cpoint);
      };
    };
  };
  return matches;
};

// Function that returns a re-numbered copy of the ConnectionsTable.
ConnectionsTable ConnectionsTable::RenumberAtoms(const std::vector<unsigned>& atom_order) const {
  ConnectionsTable renumbered_table = *this;
  for (auto& [connection, cpoints] : renumbered_table) {
    for (ConnectionPoint& cpoint : cpoints) {
      unsigned old_atom_idx = cpoint.atom_idx;
      unsigned new_atom_idx = atom_order.at(old_atom_idx);
      cpoint.atom_idx = new_atom_idx;
    };
  };
  return renumbered_table;
};

INVERTED_CONNECTIONS_MAP ConnectionsTable::GetInvertedTable() const {
  std::map<unsigned, std::multiset<Connection>> inverted_table;
  std::map<unsigned, std::multiset<Connection>>::iterator it;
  for (const auto& [connection, cpoints] : cmap) {
    for (const ConnectionPoint& cpoint : cpoints) {
      unsigned atom_idx = cpoint.GetAtomIdx();
      it = inverted_table.find(atom_idx);
      if (it == inverted_table.end()) {
        inverted_table.insert({atom_idx, {connection}});
      } else {
        it->second.insert(connection);
      };
    };
  };
  return inverted_table;
};

unsigned ConnectionsTable::Size() const {
  return n_connections;
};

bool ConnectionsTable::Empty() const {
  if (n_connections == 0) {
    return true;
  };
  return false;
};

void ConnectionsTable::Print() const {
  for (const auto& c : cmap) {
    std::cout << c.first.GetString() << ": ";
    for (const ConnectionPoint& cpoint : c.second) {
      std::cout << cpoint.atom_idx << " ";
    };
    std::cout << std::endl;
  };
};

void ConnectionsTable::ClearWeights() {
  for (auto& c : cmap) {
    CONNECTION_POINT_VECTOR& cpoints = c.second;
    for (ConnectionPoint& cpoint : cpoints) {
      cpoint.ClearWeights();
    };
  };
};


void AssignImmutableAtomIndices(RDKit::RWMol& molecule) {
  unsigned immutable_idx = 0;
  for (RDKit::Atom* atom : molecule.atoms()) {
    atom->setProp<unsigned>("ImmutableIdx", immutable_idx);
    ++immutable_idx;
  };
};

unsigned GetImmutableAtomIdx(const RDKit::Atom* atom) {
  return atom->getProp<unsigned>("ImmutableIdx");
};

RDKit::Atom* GetAtomWithImmutableIdx(RDKit::RWMol& molecule, unsigned immutable_atom_idx) {
  for (RDKit::Atom* atom : molecule.atoms()) {
    if (atom->getProp<unsigned>("ImmutableIdx") == immutable_atom_idx) {
      return atom;
    };
  };
  throw std::runtime_error("Atom with specified ImmutableIdx couldn't be found.");
};


RDKit::RWMol SanitizePseudomol(const RDKit::RWMol& pseudomol, const ConnectionsTable& connections_table, bool use_immutable_indices) {
  try {
    // Create a copy of the unsanitized pseudomol.
    RDKit::RWMol sanitized_mol = pseudomol;
    // Replace the ConnectionPoints with explicit hydrogens.
    for (const auto& [connection, cpoints] : connections_table) {
      RDKit::Bond::BondType bond_type = int_bond_type_table[connection.GetBondType()];
      unsigned bond_order = bond_type_order_table[bond_type];
      for (const ConnectionPoint& cpoint : cpoints) {
        RDKit::Atom* atom = nullptr;
        if (use_immutable_indices) {
          atom = GetAtomWithImmutableIdx(sanitized_mol, cpoint.GetAtomIdx());
        } else {
          atom = sanitized_mol.getAtomWithIdx(cpoint.GetAtomIdx());
        };
        atom->setNumExplicitHs(atom->getNumExplicitHs() + bond_order);
      };
    };
    // Loop over the molecule in search of elements that aren't guaranteed to
    // form a legal chemotype after hydrogen saturation. Note that if working
    // with other elements for which simple hydrogen saturation  isn't guaranteed
    // to yield sensible chemotypes this section ought to be expanded.
    for (RDKit::Atom* atom : sanitized_mol.atoms()) {
      // If the atom is a phosphorus:
      if (atom->getAtomicNum() == 15) {
        unsigned valence = atom->getTotalValence();
        // We only treat cases where the valence is conventional.
        if (valence == 3 || valence == 5) {
          // Calculate the "heavy atom valence" as the number of bonds involved in
          // bonding to heavy atoms.
          unsigned n_hydrogens = atom->getTotalNumHs();
          unsigned heavy_atom_valence = valence - n_hydrogens;
          // We set the number of explicit hydrogens to the number required to
          // reach the lowest valid valence.
          unsigned n_explicit_hydrogens_to_sanitize = 0;
          if (heavy_atom_valence <= 5) {
            n_explicit_hydrogens_to_sanitize = 5 - heavy_atom_valence;
          };
          if (heavy_atom_valence <= 3) {
            n_explicit_hydrogens_to_sanitize = 3 - heavy_atom_valence;
          };
          atom->setNumExplicitHs(n_explicit_hydrogens_to_sanitize);
        };
      // If the atom is a sulfur:
      } else if (atom->getAtomicNum() == 16) {
        unsigned valence = atom->getTotalValence();
        if (valence == 2 || valence == 4 || valence == 6) {
          unsigned n_hydrogens = atom->getTotalNumHs();
          unsigned heavy_atom_valence = valence - n_hydrogens;
          unsigned n_explicit_hydrogens_to_sanitize = 0;
          if (heavy_atom_valence <= 6) {
            n_explicit_hydrogens_to_sanitize = 6 - heavy_atom_valence;
          };
          if (heavy_atom_valence <= 4) {
            n_explicit_hydrogens_to_sanitize = 4 - heavy_atom_valence;
          };
          if (heavy_atom_valence <= 2) {
            n_explicit_hydrogens_to_sanitize = 2 - heavy_atom_valence;
          };
          atom->setNumExplicitHs(n_explicit_hydrogens_to_sanitize);
        };
      };
    };
    // Remove the hydrogens and implicitly wrap up molecule sanitization.
    // This is important since explicit hydrogens affect the fingerprint.
    RDKit::MolOps::removeHs(sanitized_mol);
    return sanitized_mol;
  } catch (const std::exception& e) {
    std::cout << "ERROR: Sanitization failed molecule: " << RDKit::MolToSmiles(pseudomol) << std::endl;
    throw;
  };
};


// Class Pseudofragment
Pseudofragment::Pseudofragment() = default;
// Pseudofragment::Pseudofragment(const RDKit::ROMol& pmol, const ConnectionsTable& cs,
//   const std::string& sml, bool hr, bool rp) :
//   pseudomol(pmol),
//   connections(cs),
//   smiles(sml),
//   has_ring(hr),
//   ring_part(rp) {};
Pseudofragment::Pseudofragment(const RDKit::RWMol& pmol, const ConnectionsTable& cs,
  const std::string& sml, bool hr, bool rp) :
  pseudomol(pmol),
  connections(cs),
  smiles(sml),
  has_ring(hr),
  ring_part(rp) {};

bool Pseudofragment::operator== (const Pseudofragment & other) const {
  return (smiles == other.smiles) && (ring_part == other.ring_part);
};

const RDKit::ROMol& Pseudofragment::GetPseudomol() const {
  return pseudomol;
};

const ConnectionsTable& Pseudofragment::GetConnections() const {
  return connections;
};

const std::string& Pseudofragment::GetConnectionEncodingSMILES() const {
  return smiles;
};

std::string Pseudofragment::GetSanitizedSMILES() const {
  RDKit::RWMol sanitized_fragment = SanitizePseudomol(pseudomol, connections, false);
  return RDKit::MolToSmiles(sanitized_fragment);
};

RDKit::SparseIntVect<std::uint32_t>* Pseudofragment::GetFingerprint() const {
  RDKit::RWMol sanitized_fragment = SanitizePseudomol(pseudomol, connections, false);
  return RDKit::MorganFingerprints::getFingerprint(sanitized_fragment, 2);
};

unsigned Pseudofragment::GetSize() const {
  return pseudomol.getNumAtoms();
};

bool Pseudofragment::HasRing() const {
  return has_ring;
};

bool Pseudofragment::IsRingPart() const {
  return ring_part;
};
