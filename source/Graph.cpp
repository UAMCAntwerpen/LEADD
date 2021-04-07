#include "Graph.hpp"

// Tables containing ready-made arrays that can be shuffled to create the indices
// for random container iteration. These are convenient to avoid iterator
// invalidation by shuffling the containers in question.
std::unordered_map<int, std::vector<int>> vertex_shuffle_vectors{
  {1,  {0}},
  {2,  {0, 1}},
  {3,  {0, 1, 2}},
  {4,  {0, 1, 2, 3}},
  {5,  {0, 1, 2, 3, 4}},
  {6,  {0, 1, 2, 3, 4, 5}},
  {7,  {0, 1, 2, 3, 4, 5, 6}},
  {8,  {0, 1, 2, 3, 4, 5, 6, 7}},
  {9,  {0, 1, 2, 3, 4, 5, 6, 7, 8}},
  {10, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}},
  {11, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}},
  {12, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}},
  {13, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}},
  {14, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}},
  {15, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}},
};

std::unordered_map<int, std::vector<int>> edge_shuffle_vectors{
  {1,  {0}},
  {2,  {0, 1}},
  {3,  {0, 1, 2}},
  {4,  {0, 1, 2, 3}},
  {5,  {0, 1, 2, 3, 4}},
  {6,  {0, 1, 2, 3, 4, 5}},
  {7,  {0, 1, 2, 3, 4, 5, 6}},
  {8,  {0, 1, 2, 3, 4, 5, 6, 7}},
  {9,  {0, 1, 2, 3, 4, 5, 6, 7, 8}},
  {10, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}},
  {11, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}},
  {12, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}},
  {13, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}},
  {14, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}},
  {15, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}},
};

int max_vertex_shuffle_size = vertex_shuffle_vectors.size();
int max_edge_shuffle_size = edge_shuffle_vectors.size();

// Standalone function to augment the current matching in a HopcroftKarp iteration.
void AugmentMatching(PATH& matching, const std::vector<PATH>& shortest_augmenting_paths, VERTEX_SET& free_U, VERTEX_SET& free_V) {
  PATH::iterator m_it, m_begin_it, m_end_it;
  std::vector<Edge*> edges_to_add;
  edges_to_add.reserve(shortest_augmenting_paths[0].size());
  // Loop over the shortest augmenting Paths and use them to augment the matching.
  for (const PATH& path : shortest_augmenting_paths) {
    for (Edge* edge : path) {
      m_begin_it = matching.begin();
      m_end_it = matching.end();
      m_it = matching.find(edge);
      // In order to not lose track of which Vertices are free the matched Edges
      // are used first (e.g. Vertices are freed up).
      if (m_it != m_end_it) {
        matching.erase(m_it);
        free_U.insert(*edge->GetStart());
        free_V.insert(*edge->GetEnd());
      }
      else {
        edges_to_add.push_back(edge);
      };
    };
    // And thereafter the free Edges are used (e.g. Vertices are matched).
    for (Edge* edge : edges_to_add) {
      matching.insert(edge);
      free_U.erase(*edge->GetStart());
      free_V.erase(*edge->GetEnd());
    }
    edges_to_add.clear();
  };
};

// ############################### Class Vertex ################################
// Public class members
Vertex::Vertex(int id, bool flag) : id(id), flag(flag) {};

bool Vertex::operator== (const Vertex& other) const {
  return id == other.id;
};

bool Vertex::operator!= (const Vertex& other) const {
  return id != other.id;
};

int Vertex::GetID() const {
  return id;
};

const std::list<Edge*>& Vertex::GetEdges() const {
  return edges;
};

bool Vertex::GetFlag() const {
  return flag;
};

// Private class members
Vertex* Vertex::FollowEdge(const Edge& edge) const {
  assert(std::find(edges.begin(), edges.end(), &edge) != edges.end());
  if (edge.directed) {
    assert(this == edge.start);
    return edge.end;
  }
  else {
    return edge.GetOther(this);
  };
};

Vertex* Vertex::FollowEdge(Edge* edge) const {
  assert(std::find(edges.begin(), edges.end(), edge) != edges.end());
  if (edge->directed) {
    assert(this == edge->start);
    return edge->end;
  }
  else {
    return edge->GetOther(this);
  };
};

Vertex* Vertex::FollowEdge(int edge_id) const {
  CONST_EDGE_PTR_IT it;
  it = std::find_if(edges.begin(), edges.end(),
    [=](Edge* edge) {
      return edge->id == edge_id;
    });
  assert(it != edges.end());
  if ((*it)->directed) {
    assert(this == (*it)->start);
    return (*it)->end;
  }
  else {
    return (*it)->GetOther(this);
  };
};

std::vector<Vertex*> Vertex::GetAdjacentVertices() const {
  std::vector<Vertex*> adjacent_vertices;
  adjacent_vertices.reserve(edges.size());
  for (Edge* edge : edges) {
    adjacent_vertices.push_back(edge->GetOther(this));
  };
  adjacent_vertices.resize(adjacent_vertices.size());
  return adjacent_vertices;
};

Edge* Vertex::AddEdge(Edge& edge) {
  edges.push_back(&edge);
  return edges.back();
};

Edge* Vertex::AddEdge(Edge* edge) {
  edges.push_back(edge);
  return edges.back();
};

EDGE_PTR_IT Vertex::RemoveEdge(const Edge& edge) {
  EDGE_PTR_IT it = std::find(edges.begin(), edges.end(), &edge);
  assert(it != edges.end());
  return edges.erase(it);
};

EDGE_PTR_IT Vertex::RemoveEdge(Edge* edge) {
  EDGE_PTR_IT it = std::find(edges.begin(), edges.end(), edge);
  assert(it != edges.end());
  return edges.erase(it);
};

EDGE_PTR_IT Vertex::RemoveEdge(int edge_id) {
  EDGE_PTR_IT it = std::find_if(edges.begin(), edges.end(),
    [=](Edge* edge) {
      return edge->id == edge_id;
    });
  assert(it != edges.end());
  return edges.erase(it);
};
// #############################################################################


// ################################ Class Edge #################################
// Public class members
// const qualify constructor arguments?
Edge::Edge(int id, Vertex& start, Vertex& end, bool directed, int weight, bool flag) :
  id(id), start(&start), end(&end), directed(directed), weight(weight), flag(flag) {};

Edge::Edge(int id, Vertex* start, Vertex* end, bool directed, int weight, bool flag) :
  id(id), start(start), end(end), directed(directed), weight(weight), flag(flag) {};

bool Edge::operator== (const Edge& other) const {
  return id == other.id;
};

bool Edge::operator!= (const Edge& other) const {
  return id != other.id;
};

int Edge::GetID() const {
  return id;
};

Vertex* Edge::GetStart() const {
  return start;
};

Vertex* Edge::GetEnd() const {
  return end;
};

int Edge::GetWeight() const {
  return weight;
};

bool Edge::IsDirected() const {
  return directed;
};

bool Edge::GetFlag() const {
  return flag;
};

// Private class members
Vertex* Edge::GetOther(const Vertex& vertex) const {
  assert(start == &vertex || end == &vertex);
  if (start == &vertex) {
    return end;
  }
  else if (end == &vertex) {
    return start;
  };
  throw std::runtime_error(std::string("Couldn't retrieve other vertex"));
};

Vertex* Edge::GetOther(const Vertex* vertex) const {
  assert(start == vertex || end == vertex);
  if (start == vertex) {
    return end;
  }
  else if (end == vertex) {
    return start;
  };
  throw std::runtime_error(std::string("Couldn't retrieve other vertex"));
};

Vertex* Edge::GetOther(int vertex_id) const {
  assert(start->id == vertex_id || end->id == vertex_id);
  if (start->id == vertex_id) {
    return end;
  }
  else if (end->id == vertex_id) {
    return start;
  };
  throw std::runtime_error(std::string("Couldn't retrieve other vertex"));
};
// #############################################################################


// ################################ Class Graph ################################
// Public class members
Graph::Graph() : max_vertex_id(0), max_edge_id(0) {};

const std::list<Edge>& Graph::GetEdges() const {
  return edges;
};

const std::list<Vertex>& Graph::GetVertices() const {
  return vertices;
};

// Protected class members
Edge* Graph::GetEdgeWithID(int edge_id) {
  EDGE_IT it;
  it = std::find_if(edges.begin(), edges.end(),
    [=](const Edge& edge) {
      return edge.id == edge_id;
    });
  assert(it != edges.end());
  return &(*it);
};

// When adding Edges one should verify that both the start and end Vertices are
// part of the same Graph.
Edge* Graph::AddEdge(Vertex& start, Vertex& end, bool directed, int weight, int id, bool flag) {
  int edge_id;
  if (id) {
    assert(std::all_of(edges.begin(), edges.end(),
      [=](const Edge& edge) {
        return edge.id != id;
      }));
    edge_id = id;
    if (id > max_edge_id) {
      max_edge_id = id;
    };
  }
  else {
    edge_id = max_edge_id + 1;
    ++max_edge_id;
  };
  Edge edge(edge_id, start, end, directed, weight, flag);
  edges.push_back(edge);
  Edge* vec_edge = &edges.back();
  start.AddEdge(vec_edge);
  end.AddEdge(vec_edge);
  return vec_edge;
};

Edge* Graph::AddEdge(Vertex* start, Vertex* end, bool directed, int weight, int id, bool flag) {
  int edge_id;
  if (id) {
    assert(std::all_of(edges.begin(), edges.end(),
      [=](const Edge& edge) {
        return edge.id != id;
      }));
    edge_id = id;
    if (id > max_edge_id) {
      max_edge_id = id;
    };
  }
  else {
    edge_id = max_edge_id + 1;
    ++max_edge_id;
  };
  Edge edge(edge_id, start, end, directed, weight, flag);
  edges.push_back(edge);
  Edge* vec_edge = &edges.back();
  start->AddEdge(vec_edge);
  end->AddEdge(vec_edge);
  return vec_edge;
};

EDGE_IT Graph::RemoveEdge(const Edge& edge) {
  EDGE_IT it = std::find(edges.begin(), edges.end(), edge);
  assert(it != edges.end());
  edge.start->RemoveEdge(edge);
  edge.end->RemoveEdge(edge);
  return edges.erase(it);
};

EDGE_IT Graph::RemoveEdge(Edge* edge) {
  EDGE_IT it = std::find(edges.begin(), edges.end(), *edge);
  assert(it != edges.end());
  edge->start->RemoveEdge(edge);
  edge->end->RemoveEdge(edge);
  return edges.erase(it);
};

EDGE_IT Graph::RemoveEdge(int edge_id) {
  EDGE_IT it = std::find_if(edges.begin(), edges.end(),
    [=](const Edge& edge) {
      return edge.id == edge_id;
    });
  assert(it != edges.end());
  it->start->RemoveEdge(*it);
  it->end->RemoveEdge(*it);
  return edges.erase(it);
};

void Graph::SetEdgeFlags(bool flag) {
  for (Edge& edge : edges) {
    edge.flag = flag;
  };
};

bool Graph::HasEdge(const Edge& edge) const {
  CONST_EDGE_IT it, end_it = edges.end();
  it = std::find(edges.begin(), end_it, edge);
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

bool Graph::HasEdge(const Edge* edge) const {
  CONST_EDGE_IT it, end_it = edges.end();
  it = std::find(edges.begin(), end_it, *edge);
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

bool Graph::HasEdge(int edge_id) const {
  CONST_EDGE_IT it, end_it = edges.end();
  it = std::find_if(edges.begin(), end_it,
    [=](const Edge& edge) {
      return edge.id == edge_id;
    });
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

// Private class members
// To be moved into protected if there are no classes inheriting from Graph.
Edge* Graph::AddEdge(int start_id, int end_id, bool directed, int weight, int id, bool flag) {
  int edge_id;
  if (id) {
    assert(std::all_of(edges.begin(), edges.end(),
      [=](const Edge& edge) {
        return edge.id != id;
      }));
    edge_id = id;
    if (id > max_edge_id) {
      max_edge_id = id;
    };
  }
  else {
    edge_id = max_edge_id + 1;
    ++max_edge_id;
  };
  Vertex* start = GetVertexWithID(start_id);
  Vertex* end = GetVertexWithID(end_id);
  Edge edge(edge_id, start, end, directed, weight, flag);
  edges.push_back(edge);
  Edge* vec_edge = &edges.back();
  start->AddEdge(vec_edge);
  end->AddEdge(vec_edge);
  return vec_edge;
};

Vertex* Graph::GetVertexWithID(int vertex_id) {
  VERTEX_IT it;
  it = std::find_if(vertices.begin(), vertices.end(),
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  assert(it != vertices.end());
  return &(*it);

};

Vertex* Graph::AddVertex(int id, bool flag) {
  int vertex_id;
  if (id) {
    assert(std::all_of(vertices.begin(), vertices.end(),
      [=](const Vertex& vertex) {
        return vertex.id != id;
      }));
    vertex_id = id;
    if (id > max_vertex_id) {
      max_vertex_id = id;
    };
  }
  else {
    vertex_id = max_vertex_id + 1;
    ++max_vertex_id;
  };
  vertices.push_back(Vertex(vertex_id, flag));
  return &vertices.back();
};

VERTEX_IT Graph::RemoveVertex(const Vertex& vertex) {
  VERTEX_IT it = std::find(vertices.begin(), vertices.end(), vertex);
  assert(it != vertices.end());
  for (Edge* edge : vertex.edges) {
    RemoveEdge(edge);
  };
  return vertices.erase(it);
};

VERTEX_IT Graph::RemoveVertex(Vertex* vertex) {
  VERTEX_IT it = std::find(vertices.begin(), vertices.end(), *vertex);
  assert(it != vertices.end());
  for (Edge* edge : vertex->edges) {
    RemoveEdge(edge);
  };
  return vertices.erase(it);
};

VERTEX_IT Graph::RemoveVertex(int vertex_id) {
  VERTEX_IT it = std::find_if(vertices.begin(), vertices.end(),
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  assert(it != vertices.end());
  for (Edge* edge : it->edges) {
    RemoveEdge(edge);
  };
  return vertices.erase(it);
};

void Graph::SetVertexFlags(bool flag) {
  for (Vertex& vertex : vertices) {
    vertex.flag = flag;
  };
};

bool Graph::HasVertex(const Vertex& vertex) const {
  CONST_VERTEX_IT it, end_it = vertices.end();
  it = std::find(vertices.begin(), end_it, vertex);
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

bool Graph::HasVertex(const Vertex* vertex) const {
  CONST_VERTEX_IT it, end_it = vertices.end();
  it = std::find(vertices.begin(), end_it, *vertex);
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

bool Graph::HasVertex(int vertex_id) const {
  CONST_VERTEX_IT it, end_it = vertices.end();
  it = std::find_if(vertices.begin(), end_it,
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};
// #############################################################################


// ############################ Class LayeredGraph #############################
// Public class members
LayeredGraph::LayeredGraph() = default;

const std::unordered_map<int, std::list<Vertex>>& LayeredGraph::GetLayeredVertices() const {
  return layered_vertices;
};

void LayeredGraph::Print() {
  for (const auto& layer : layered_vertices) {
    std::cout << "Layer " << layer.first << ": ";
    for (const Vertex& vertex : layer.second) {
      std::cout << vertex.id << " ";
    };
    std::cout << std::endl;
  };
  std::cout << "Edges:\n";
  for (const Edge& edge : edges) {
    std::cout << edge.id << ": " << edge.start->id << "->" << edge.end->id << std::endl;
  };
};

// Private class members
std::list<Vertex>& LayeredGraph::operator[] (int layer_idx) {
  return layered_vertices[layer_idx];
};

void LayeredGraph::AddLayer() {
  layered_vertices.insert({ max_layer_id + 1, std::list<Vertex>{} });
  ++max_layer_id;
};

void LayeredGraph::AddLayer(int layer_id) {
  assert(layered_vertices.find(layer_id) == layered_vertices.end());
  layered_vertices.insert({ layer_id, std::list<Vertex>{} });
  ++max_layer_id;
};

Vertex* LayeredGraph::GetVertexWithID(int vertex_id, int layer_id) {
  std::list<Vertex>& layer = layered_vertices[layer_id];
  VERTEX_IT it = std::find_if(layer.begin(), layer.end(),
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  assert(it != layer.end());
  return &(*it);
};

Vertex* LayeredGraph::GetVertexWithID(int vertex_id) {
  for (auto& layer : layered_vertices) {
    for (Vertex& vertex : layer.second) {
      if (vertex.id == vertex_id) {
        return &vertex;
      };
    };
  };
  throw std::runtime_error(std::string("Couldn't get Vertex with ID"));
};

Vertex* LayeredGraph::AddVertex(int layer_id, int id, bool flag) {
  bool has_layer;
  int vertex_id;
  std::unordered_map<int, std::list<Vertex>>::const_iterator it = layered_vertices.find(layer_id);
  if (it != layered_vertices.end()) {
    has_layer = true;
  }
  else {
    has_layer = false;
  };
  if (id) {
    for (const auto& layer : layered_vertices) {
      assert(std::all_of(layer.second.begin(), layer.second.end(),
        [=](const Vertex& vertex) {
          return vertex.id != id;
        }));
    };
    vertex_id = id;
    if (id > max_vertex_id) {
      max_vertex_id = id;
    };
  }
  else {
    vertex_id = max_vertex_id + 1;
    ++max_vertex_id;
  };
  if (!has_layer) {
    AddLayer(layer_id);
  };
  layered_vertices[layer_id].push_back(Vertex(vertex_id, flag));
  return &layered_vertices[layer_id].back();
};

VERTEX_IT LayeredGraph::RemoveVertex(const Vertex & vertex) {
  VERTEX_IT it, begin_it, end_it;
  for (auto& layer : layered_vertices) {
    begin_it = layer.second.begin();
    end_it = layer.second.end();
    it = std::find(begin_it, end_it, vertex);
    if (it != end_it) {
      for (Edge* edge : it->edges) {
        RemoveEdge(edge);
      };
      return layer.second.erase(it);
    };
  };
  throw std::runtime_error(std::string("Couldn't remove Vertex"));
};

VERTEX_IT LayeredGraph::RemoveVertex(Vertex * vertex) {
  VERTEX_IT it, begin_it, end_it;
  for (auto& layer : layered_vertices) {
    begin_it = layer.second.begin();
    end_it = layer.second.end();
    it = std::find(begin_it, end_it, *vertex);
    if (it != end_it) {
      for (Edge* edge : it->edges) {
        RemoveEdge(edge);
      };
      return layer.second.erase(it);
    };
  };
  throw std::runtime_error(std::string("Couldn't remove Vertex"));
};

VERTEX_IT LayeredGraph::RemoveVertex(int vertex_id) {
  VERTEX_IT it, begin_it, end_it;
  for (auto& layer : layered_vertices) {
    begin_it = layer.second.begin();
    end_it = layer.second.end();
    it = std::find_if(begin_it, end_it,
      [=](const Vertex& vertex) {
        return vertex.id == vertex_id;
      });
    if (it != end_it) {
      for (Edge* edge : it->edges) {
        RemoveEdge(edge);
      };
      return layer.second.erase(it);
    };
  };
  throw std::runtime_error(std::string("Couldn't remove Vertex"));
};

void LayeredGraph::SetVertexFlags(bool flag) {
  for (auto& layer : layered_vertices) {
    for (Vertex& vertex : layer.second) {
      vertex.flag = flag;
    };
  };
};

bool LayeredGraph::HasVertex(const Vertex & vertex) const {
  CONST_VERTEX_IT it, end_it;
  for (const auto& layer : layered_vertices) {
    end_it = layer.second.end();
    it = std::find(layer.second.begin(), end_it, vertex);
    if (it != end_it) {
      return true;
    };
  };
  return false;
};

bool LayeredGraph::HasVertex(const Vertex * vertex) const {
  CONST_VERTEX_IT it, end_it;
  for (const auto& layer : layered_vertices) {
    end_it = layer.second.end();
    it = std::find(layer.second.begin(), end_it, *vertex);
    if (it != end_it) {
      return true;
    };
  };
  return false;
};

bool LayeredGraph::HasVertex(int vertex_id) const {
  CONST_VERTEX_IT it, end_it;
  for (const auto& layer : layered_vertices) {
    end_it = layer.second.end();
    it = std::find_if(layer.second.begin(), end_it,
      [=](const Vertex& vertex) {
        return vertex.id == vertex_id;
      });
    if (it != end_it) {
      return true;
    };
  };
  return false;
};
// #############################################################################


// ########################### Class BipartiteGraph ############################
// Public class members
BipartiteGraph::BipartiteGraph() = default;

BipartiteGraph::BipartiteGraph(int u_size, int v_size, const std::vector<std::pair<int, int>> & connections) {
  for (int i = 1; i <= u_size; ++i) {
    AddVertexToU(i);
  };
  for (int i = u_size + 1; i <= u_size + v_size; ++i) {
    AddVertexToV(i);
  };
  for (const auto& pair : connections) {
    Vertex* start = GetVertexWithID(pair.first);
    Vertex* end = GetVertexWithID(pair.second);
    AddEdge(start, end);
  };
};

const std::list<Vertex>& BipartiteGraph::GetU() const {
  return U;
};

const std::list<Vertex>& BipartiteGraph::GetV() const {
  return V;
};

void BipartiteGraph::Print() {
  std::cout << "Vertices in U: ";
  for (const Vertex& vertex : U) {
    std::cout << vertex.id << " ";
  };
  std::cout << std::endl;
  std::cout << "Vertices in V: ";
  for (const Vertex& vertex : V) {
    std::cout << vertex.id << " ";
  };
  std::cout << std::endl;
  std::cout << "Edges:\n";
  for (const Edge& edge : edges) {
    std::cout << edge.id << ": " << edge.start->id << "->" << edge.end->id << std::endl;
  };
};

// Private class members
Vertex* BipartiteGraph::GetVertexWithID(int vertex_id) {
  VERTEX_IT it;
  it = std::find(U.begin(), U.end(), vertex_id);
  if (it != U.end()) {
    return &(*it);
  };
  it = std::find(V.begin(), V.end(), vertex_id);
  assert(it != V.end());
  return &(*it);
};

Vertex* BipartiteGraph::AddVertexToU(int id, bool flag) {
  int vertex_id;
  if (id) {
    assert(std::all_of(U.begin(), U.end(),
      [=](const Vertex& vertex) {
        return vertex.id != id;
      }));
    assert(std::all_of(V.begin(), V.end(),
      [=](const Vertex& vertex) {
        return vertex.id != id;
      }));
    vertex_id = id;
    if (id > max_vertex_id) {
      max_vertex_id = id;
    };
  }
  else {
    vertex_id = max_vertex_id + 1;
    ++max_vertex_id;
  };
  U.push_back(Vertex(vertex_id, flag));
  return &U.back();
};

Vertex* BipartiteGraph::AddVertexToV(int id, bool flag) {
  int vertex_id;
  if (id) {
    assert(std::all_of(U.begin(), U.end(),
      [=](const Vertex& vertex) {
        return vertex.id != id;
      }));
    assert(std::all_of(V.begin(), V.end(),
      [=](const Vertex& vertex) {
        return vertex.id != id;
      }));
    vertex_id = id;
    if (id > max_vertex_id) {
      max_vertex_id = id;
    };
  }
  else {
    vertex_id = max_vertex_id + 1;
    ++max_vertex_id;
  };
  V.push_back(Vertex(vertex_id, flag));
  return &V.back();
};

VERTEX_IT BipartiteGraph::RemoveVertex(const Vertex & vertex) {
  VERTEX_IT end_it = U.end();
  VERTEX_IT it = std::find(U.begin(), end_it, vertex);
  if (it != end_it) {
    for (Edge* edge : it->edges) {
      RemoveEdge(edge);
    };
    return U.erase(it);
  }
  else {
    end_it = V.end();
    it = std::find(V.begin(), end_it, vertex);
    assert(it != end_it);
    for (Edge* edge : it->edges) {
      RemoveEdge(edge);
    };
    return V.erase(it);
  };
};

VERTEX_IT BipartiteGraph::RemoveVertex(Vertex * vertex) {
  VERTEX_IT end_it = U.end();
  VERTEX_IT it = std::find(U.begin(), end_it, *vertex);
  if (it != end_it) {
    for (Edge* edge : it->edges) {
      RemoveEdge(edge);
    };
    return U.erase(it);
  }
  else {
    end_it = V.end();
    it = std::find(V.begin(), end_it, *vertex);
    assert(it != end_it);
    for (Edge* edge : it->edges) {
      RemoveEdge(edge);
    };
    return V.erase(it);
  };
};

VERTEX_IT BipartiteGraph::RemoveVertex(int vertex_id) {
  VERTEX_IT end_it = U.end();
  VERTEX_IT it = std::find_if(U.begin(), end_it,
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  if (it != end_it) {
    for (Edge* edge : it->edges) {
      RemoveEdge(edge);
    };
    return U.erase(it);
  }
  else {
    end_it = V.end();
    std::list<Vertex>::iterator it = std::find_if(V.begin(), end_it,
      [=](const Vertex& vertex) {
        return vertex.id == vertex_id;
      });
    assert(it != end_it);
    for (Edge* edge : it->edges) {
      RemoveEdge(edge);
    };
    return V.erase(it);
  };
};

void BipartiteGraph::SetVertexFlags(bool flag) {
  for (Vertex& vertex : U) {
    vertex.flag = flag;
  };
  for (Vertex& vertex : V) {
    vertex.flag = flag;
  };
};

bool BipartiteGraph::HasVertex(const Vertex & vertex) const {
  CONST_VERTEX_IT it, end_it = U.end();
  it = std::find(U.begin(), end_it, vertex);
  if (it != end_it) {
    return true;
  };
  end_it = V.end();
  it = std::find(V.begin(), end_it, vertex);
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

bool BipartiteGraph::HasVertex(const Vertex * vertex) const {
  CONST_VERTEX_IT it, end_it = U.end();
  it = std::find(U.begin(), end_it, *vertex);
  if (it != end_it) {
    return true;
  };
  end_it = V.end();
  it = std::find(V.begin(), end_it, *vertex);
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

bool BipartiteGraph::HasVertex(int vertex_id) const {
  CONST_VERTEX_IT it, end_it = U.end();
  it = std::find_if(U.begin(), end_it,
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  if (it != end_it) {
    return true;
  };
  end_it = V.end();
  it = std::find_if(V.begin(), end_it,
    [=](const Vertex& vertex) {
      return vertex.id == vertex_id;
    });
  if (it != end_it) {
    return true;
  }
  else {
    return false;
  };
};

VERTEX_SET BipartiteGraph::GetUSet() const {
  return std::unordered_set<Vertex>(U.begin(), U.end());
};

VERTEX_SET BipartiteGraph::GetVSet() const {
  return std::unordered_set<Vertex>(V.begin(), V.end());
};

std::vector<PATH> BipartiteGraph::ShortestAugmentingPaths(const PATH & matching, const VERTEX_SET & free_U, const VERTEX_SET & free_V, std::mt19937 & prng) {
  // Initialize a vector to store the shortest augmenting Paths.
  std::vector<PATH> shortest_augmenting_paths;

  // If there are no free Vertices on both sides of the BipartiteGraph there
  // can't be an augmenting Path and the Maximum Cardinality Matching has
  // been found.
  if (free_U.size() == 0 || free_V.size() == 0) {
    return shortest_augmenting_paths;
  };

  // Initialize an iterator to check if an Edge is in the current matching.
  PATH::const_iterator m_end_it = matching.end();

  // Initialize iterators to check if a given Vertex is free.
  VERTEX_SET::const_iterator v_end_it = free_V.end();

  // Initialize a LayeredGraph to store the results of a Breadth-First Search that:
  //  (1) Connects layers alternating between the use of matching and free Edges.
  //      this means that the search is limited to potential augmenting Paths.
  //  (2) Terminates when encountering an augmenting Path (i.e. free Vertices
  //      in any layer other than the first one). With this termination criterion
  //      we ensure that the LayeredGraph contains all shortest augmenting Paths.
  LayeredGraph layered_graph;

  // Starting with an arbitrary set of Vertices from the BipartiteGraph (by
  // convention U), add the free Vertices of the set as the first layer of the
  // LayeredGraph.
  int layer_id = 1;
  for (const Vertex& vertex : free_U) {
    layered_graph.AddVertex(layer_id, vertex.id);
  };

  // Perform the Breadth-First Search by following the Edges associated with the
  // Vertices of the LayeredGraph's last layer. If they lead to a free Vertex
  // the process terminates.
  SetVertexFlags(false);
  bool has_augmenting_path = false;
  bool discovered_vertex = false;
  int other_vertex_id;
  Vertex* other_vertex;
  Vertex* bipartite_vertex;
  Vertex* other_bipartite_vertex;
  while (!has_augmenting_path) {
    // If the layer number is odd it means that in order to continue looking for
    // augmenting Paths the next BipartiteGraph's Edge in the Path:
    //  (1) Ought to be a free Edge.
    //  (2) Terminates in a Vertex that hasn't been visited yet.
    //  (3) Terminates in a Vertex from set V.
    discovered_vertex = false;
    if (layer_id % 2) {
      // Loop over the Vertices that haven't been visited yet.
      for (Vertex& vertex : layered_graph[layer_id]) {
        bipartite_vertex = GetVertexWithID(vertex.id);
        if (!bipartite_vertex->flag) {
          // Loop over the Edges of the Vertex.
          for (Edge* edge : bipartite_vertex->edges) {
            // Continue only if the Edge is free.
            if (matching.find(edge) == m_end_it) {
              other_bipartite_vertex = edge->GetOther(bipartite_vertex);
              if (!other_bipartite_vertex->flag) {
                discovered_vertex = true;
                other_vertex_id = other_bipartite_vertex->id;
                if (layered_graph.HasVertex(other_vertex_id)) {
                  other_vertex = layered_graph.GetVertexWithID(other_vertex_id);
                }
                else {
                  other_vertex = layered_graph.AddVertex(layer_id + 1, other_vertex_id);
                };
                layered_graph.AddEdge(vertex, *other_vertex, false, 1, edge->id);
                // If the termination Vertex is free, an augmenting Path has been
                // found and this is the last iteration of the loop.
                if (free_V.find(*other_vertex) != v_end_it) {
                  has_augmenting_path = true;
                };
              };
            };
          };
          // Label the Graph's Vertex as visited. This is important to avoid falling
          // into cycles during the search.
          bipartite_vertex->flag = true;
        };
      };
      // Conversely, if the layer number is even the next BipartiteGraph's Edge:
      //  (1) Ought to be a matched Edge.
      //  (2) Terminates in a Vertex that hasn't been visited yet.
      //  (3) Terminates in a Vertex from set U.
    }
    else {
      for (Vertex& vertex : layered_graph[layer_id]) {
        bipartite_vertex = GetVertexWithID(vertex.id);
        if (!bipartite_vertex->flag) {
          // Loop over the Edges of the Vertex.
          for (Edge* edge : bipartite_vertex->edges) {
            // Continue only if the Edge is matched.
            if (matching.find(edge) != m_end_it) {
              other_bipartite_vertex = edge->GetOther(bipartite_vertex);
              if (!other_bipartite_vertex->flag) {
                discovered_vertex = true;
                other_vertex_id = other_bipartite_vertex->id;
                if (layered_graph.HasVertex(other_vertex_id)) {
                  other_vertex = layered_graph.GetVertexWithID(other_vertex_id);
                }
                else {
                  other_vertex = layered_graph.AddVertex(layer_id + 1, other_vertex_id);
                };
                layered_graph.AddEdge(vertex, *other_vertex, false, 1, edge->id);
              };
            };
          };
          // Label the graph's vertex as visited. this is important to avoid
          // falling into cycles during the search.
          bipartite_vertex->flag = true;
        };
      };
    };
    ++layer_id;
    if (!discovered_vertex) {
      if (layer_id % 2) {
        break;
      }
      else {
        return shortest_augmenting_paths;
      };
    };
  };

  // Find the free Vertices of the last layer of the LayeredGraph and trace back
  // to the first layer of the LayeredGraph. Due to the constraints imposed
  // during the Breadth-First Search we know that all Paths starting from a free
  // Vertex in the last layer are augmenting Paths. Therefore, it holds true that
  // the last layer is made out of V-vertices.
  // The augmenting Paths have to be maximally disjoint. Hence, during the search
  // the identified augmenting Paths are flagged as intransitable. Combined
  // these steps are analogous to a Depth-First Search.

  // Loop over the free Vertices of the LayeredGraph's last layer that haven't
  // been used yet in a maximally disjoint Path (in random order).
  std::vector<int>* vertex_shuffle;
  std::vector<int>* edge_shuffle;
  std::vector<int> large_vertex_shuffle, large_edge_shuffle;

  std::list<Vertex>& last_layer = layered_graph[layer_id];
  int shuffle_size = last_layer.size();

  if (shuffle_size <= max_vertex_shuffle_size) {
    vertex_shuffle = &vertex_shuffle_vectors[shuffle_size];
  }
  else {
    large_vertex_shuffle.resize(shuffle_size);
    std::iota(large_vertex_shuffle.begin(), large_vertex_shuffle.end(), 0);
    vertex_shuffle = &large_vertex_shuffle;
  };
  std::shuffle(vertex_shuffle->begin(), vertex_shuffle->end(), prng);

  // assert(shuffle_size < max_vertex_shuffle_size);
  // std::vector<int>& vertex_shuffle = vertex_shuffle_vectors[shuffle_size];
  // std::shuffle(vertex_shuffle.begin(), vertex_shuffle.end(), prng);

  Vertex* current_vertex;
  Vertex* next_vertex;
  Vertex* interrogated_vertex;
  VERTEX_IT last_it, last_begin_it = last_layer.begin(), above_it, above_end_it;
  EDGE_PTR_IT e_it, e_begin_it;
  bool valid_edge_found;
  PATH path;
  std::vector<Vertex*> vertices_in_path;
  for (int v_idx : *vertex_shuffle) {
    last_it = last_begin_it;
    std::advance(last_it, v_idx);
    // Check if the Vertex is free.
    if (free_V.find(*last_it) != v_end_it) {
      path.clear();
      vertices_in_path.clear();
      current_vertex = &(*last_it);
      current_vertex->flag = true;
      // Trace back from the chosen free Vertex of the LayeredGraph's last layer
      // to a free Vertex of the LayeredGraph's first layer. Note that all
      // Vertices in the first layer are guaranteed to be free.
      for (int above_layer_id = layer_id - 1; above_layer_id > 0; --above_layer_id) {
        std::list<Vertex>& above_layer = layered_graph[above_layer_id];
        above_end_it = above_layer.end();
        valid_edge_found = false;
        // Pick a random valid Vertex's Edge. Valid edges are those that when
        // followed lead to a Vertex in the upper layer of the LayeredGraph
        // that hasn't been visited yet.

        e_begin_it = current_vertex->edges.begin();
        shuffle_size = current_vertex->edges.size();

        if (shuffle_size <= max_edge_shuffle_size) {
          edge_shuffle = &edge_shuffle_vectors[shuffle_size];
        }
        else {
          large_edge_shuffle.resize(shuffle_size);
          std::iota(large_edge_shuffle.begin(), large_edge_shuffle.end(), 0);
          edge_shuffle = &large_edge_shuffle;
        };
        std::shuffle(edge_shuffle->begin(), edge_shuffle->end(), prng);

        // assert(shuffle_size < max_edge_shuffle_size);
        // std::vector<int>& edge_shuffle = edge_shuffle_vectors[shuffle_size];
        // std::shuffle(edge_shuffle.begin(), edge_shuffle.end(), prng);

        for (int e_idx : *edge_shuffle) {
          e_it = e_begin_it;
          std::advance(e_it, e_idx);
          interrogated_vertex = (*e_it)->GetOther(current_vertex);
          above_it = std::find(above_layer.begin(), above_end_it, *interrogated_vertex);
          if (!interrogated_vertex->flag && above_it != above_end_it) {
            next_vertex = interrogated_vertex;
            valid_edge_found = true;
            break;
          };
        };
        // If the Vertex has no valid Edges move on to the next free Vertex of
        // the last layer.
        if (!valid_edge_found) {
          goto next_iteration;
        };
        // Add the chosen Edge to the Path and move up to the next layer.
        path.insert(GetEdgeWithID((*e_it)->id));
        vertices_in_path.push_back(next_vertex);
        current_vertex = next_vertex;
      };
      // If the Path was fully defined (i.e. no skips to the next Vertex), add
      // it to the return container and flag the Vertices involved as visited.
      shortest_augmenting_paths.push_back(path);
      for (Vertex* vertex : vertices_in_path) {
        vertex->flag = true;
      };
    };
  next_iteration:
    continue;
  };

  // Return the shortest augmenting Paths.
  return shortest_augmenting_paths;
};

PATH BipartiteGraph::HopcroftKarp(std::mt19937 & prng) {
  // Initialize an empty matching.
  PATH matching;

  // Get the free Vertices of the BipartiteGraph. In the first iteration
  // all vertices are free.
  VERTEX_SET free_U = GetUSet();
  VERTEX_SET free_V = GetVSet();

  // Define the perfect matching cardinality to be achieved.
  int n_u_vertices = free_U.size();
  int n_v_vertices = free_V.size();
  if (n_u_vertices == 0 || n_v_vertices == 0) {
    return matching;
  };
  int cardinality = 0, max_cardinality;
  if (n_u_vertices >= n_v_vertices) {
    max_cardinality = n_v_vertices;
  }
  else {
    max_cardinality = n_u_vertices;
  };

  // Find the shortest augmenting Paths given the current matching.
  std::vector<PATH> shortest_augmenting_paths = ShortestAugmentingPaths(matching, free_U, free_V, prng);
  if (shortest_augmenting_paths.empty()) {
    return matching;
  };

  // Augment the matching with the aforementioned Paths and update which
  // Vertices are free.
  AugmentMatching(matching, shortest_augmenting_paths, free_U, free_V);
  cardinality = matching.size();

  // Repeat the process until either a perfect matching has been found or no
  // further augmenting Paths exist, in which case a maximum matching was found.
  while (cardinality < max_cardinality) {
    shortest_augmenting_paths = ShortestAugmentingPaths(matching, free_U, free_V, prng);
    if (shortest_augmenting_paths.empty()) {
      break;
    };
    AugmentMatching(matching, shortest_augmenting_paths, free_U, free_V);
    cardinality = matching.size();
  };

  // Return the solution to the maximum bipartite matching problem.
  return matching;
};

// #############################################################################

