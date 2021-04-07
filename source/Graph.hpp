#pragma once
#ifndef _GRAPH_HPP_
#define _GRAPH_HPP_

#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <list>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <assert.h>

class Vertex;
class Edge;

typedef std::list<Vertex>::iterator VERTEX_IT;
typedef std::list<Vertex>::const_iterator CONST_VERTEX_IT;
typedef std::list<Edge>::iterator EDGE_IT;
typedef std::list<Edge>::const_iterator CONST_EDGE_IT;
typedef std::list<Edge*>::iterator EDGE_PTR_IT;
typedef std::list<Edge*>::const_iterator CONST_EDGE_PTR_IT;
typedef std::unordered_set<Edge*> PATH;
typedef std::unordered_set<Vertex> VERTEX_SET;

class Vertex {
private:
  int id;
  std::list<Edge*> edges;
  bool flag;

public:
  Vertex(int id, bool flag = false);
  bool operator== (const Vertex& other) const;
  bool operator!= (const Vertex& other) const;

  int GetID() const;
  const std::list<Edge*>& GetEdges() const;
  bool GetFlag() const;

private:
  Vertex* FollowEdge(const Edge& edge) const;
  Vertex* FollowEdge(Edge* edge) const;
  Vertex* FollowEdge(int edge_id) const;
  std::vector<Vertex*> GetAdjacentVertices() const;

  Edge* AddEdge(Edge& edge);
  Edge* AddEdge(Edge* edge);
  EDGE_PTR_IT RemoveEdge(const Edge& edge);
  EDGE_PTR_IT RemoveEdge(Edge* edge);
  EDGE_PTR_IT RemoveEdge(int edge_id);

  friend struct std::hash<Vertex>;
  friend class Edge;
  friend class Graph;
  friend class LayeredGraph;
  friend class BipartiteGraph;
};

template<>
struct std::hash<Vertex> {
  std::size_t operator() (const Vertex& vertex) const {
    return std::hash<int>()(vertex.id);
  };
};

class Edge {
private:
  int id;
  Vertex* start;
  Vertex* end;
  bool directed;
  int weight;
  bool flag;

public:
  Edge(int id, Vertex& start, Vertex& end, bool directed = false, int weight = 1, bool flag = false);
  Edge(int id, Vertex* start, Vertex* end, bool directed = false, int weight = 1, bool flag = false);
  bool operator== (const Edge& other) const;
  bool operator!= (const Edge& other) const;

  int GetID() const;
  Vertex* GetStart() const;
  Vertex* GetEnd() const;
  int GetWeight() const;
  bool IsDirected() const;
  bool GetFlag() const;

private:
  Vertex* GetOther(const Vertex& vertex) const;
  Vertex* GetOther(const Vertex* vertex) const;
  Vertex* GetOther(int vertex_id) const;

  friend struct std::hash<Edge>;
  friend class Vertex;
  friend class Graph;
  friend class LayeredGraph;
  friend class BipartiteGraph;
};

template<>
struct std::hash<Edge> {
  std::size_t operator() (const Edge& edge) const {
    return std::hash<int>()(edge.id);
  };
};

class Graph {
private:
  std::list<Vertex> vertices;

protected:
  std::list<Edge> edges;
  int max_vertex_id;
  int max_edge_id;

public:
  Graph();
  const std::list<Vertex>& GetVertices() const;
  const std::list<Edge>& GetEdges() const;

protected:
  Edge* GetEdgeWithID(int edge_id);
  Edge* AddEdge(Vertex& start, Vertex& end, bool directed = false, int weight = 1, int id = 0, bool flag = false);
  Edge* AddEdge(Vertex* start, Vertex* end, bool directed = false, int weight = 1, int id = 0, bool flag = false);
  EDGE_IT RemoveEdge(const Edge& edge);
  EDGE_IT RemoveEdge(Edge* edge);
  EDGE_IT RemoveEdge(int edge_id);
  void SetEdgeFlags(bool flag);
  bool HasEdge(const Edge& edge) const;
  bool HasEdge(const Edge* edge) const;
  bool HasEdge(int edge_id) const;

private:
  Edge* AddEdge(int start_id, int end_id, bool directed = false, int weight = 1, int id = 0, bool flag = false);
  Vertex* GetVertexWithID(int vertex_id);
  Vertex* AddVertex(int id = 0, bool flag = false);
  VERTEX_IT RemoveVertex(const Vertex& vertex);
  VERTEX_IT RemoveVertex(Vertex* vertex);
  VERTEX_IT RemoveVertex(int vertex_id);
  void SetVertexFlags(bool flag);
  bool HasVertex(const Vertex& vertex) const;
  bool HasVertex(const Vertex* vertex) const;
  bool HasVertex(int vertex_id) const;
};

class LayeredGraph : private Graph {
  std::unordered_map<int, std::list<Vertex>> layered_vertices;
  int max_layer_id = 0;

public:
  LayeredGraph();
  const std::unordered_map<int, std::list<Vertex>>& GetLayeredVertices() const;
  void Print();

private:
  std::list<Vertex>& operator[] (int layer_idx);
  void AddLayer();
  void AddLayer(int layer_id);
  Vertex* GetVertexWithID(int vertex_id, int layer_id);
  Vertex* GetVertexWithID(int vertex_id);
  Vertex* AddVertex(int layer_id, int id = 0, bool flag = false);
  VERTEX_IT RemoveVertex(const Vertex& vertex);
  VERTEX_IT RemoveVertex(Vertex* vertex);
  VERTEX_IT RemoveVertex(int vertex_id);
  void SetVertexFlags(bool flag);
  bool HasVertex(const Vertex& vertex) const;
  bool HasVertex(const Vertex* vertex) const;
  bool HasVertex(int vertex_id) const;

  friend class Graph;
  friend class BipartiteGraph;
};

class BipartiteGraph : private Graph {
  std::list<Vertex> U, V;

public:
  BipartiteGraph();
  BipartiteGraph(int u_size, int v_size, const std::vector<std::pair<int, int>>& connections);
  const std::list<Vertex>& GetU() const;
  const std::list<Vertex>& GetV() const;
  void Print();
  PATH HopcroftKarp(std::mt19937& prng);

private:
  Vertex* GetVertexWithID(int vertex_id);
  Vertex* AddVertexToU(int id = 0, bool flag = false);
  Vertex* AddVertexToV(int id = 0, bool flag = false);
  VERTEX_IT RemoveVertex(const Vertex& vertex);
  VERTEX_IT RemoveVertex(Vertex* vertex);
  VERTEX_IT RemoveVertex(int vertex_id);
  void SetVertexFlags(bool flag);
  bool HasVertex(const Vertex& vertex) const;
  bool HasVertex(const Vertex* vertex) const;
  bool HasVertex(int vertex_id) const;

  VERTEX_SET GetUSet() const;
  VERTEX_SET GetVSet() const;
  std::vector<PATH> ShortestAugmentingPaths(const PATH& matching, const VERTEX_SET& free_U, const VERTEX_SET& free_V, std::mt19937& prng);
};

void AugmentMatching(PATH& matching, const std::vector<PATH>& shortest_augmenting_paths, VERTEX_SET& free_U, VERTEX_SET& free_V);

#endif
