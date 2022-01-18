#pragma once
#ifndef _PSEUDOFRAGMENT_DB_HPP_
#define _PSEUDOFRAGMENT_DB_HPP_

#include <omp.h>
#include <sqlite3.h>
#include <boost/progress.hpp>
#include "Pseudofragment.hpp"

class PseudofragmentDB {
  std::string database_path;
  sqlite3* database;
  sqlite3_stmt* select_n_pseudofragments;
  sqlite3_stmt* select_n_connections;
  sqlite3_stmt* select_source_molecule_id;
  sqlite3_stmt* select_pseudofragment_id;
  sqlite3_stmt* select_connection_id;
  sqlite3_stmt* select_connection;
  sqlite3_stmt* select_pseudofragment;
  sqlite3_stmt* insert_source_molecule;
  sqlite3_stmt* insert_pseudofragment;
  sqlite3_stmt* insert_connection;
  sqlite3_stmt* insert_source_pseudofragment_relationship;
  sqlite3_stmt* insert_pseudofragment_connection_relationship;

public:
  PseudofragmentDB(const std::string& database_path);

  void BeginTransaction() const;
  void CommitTransaction() const;
  void Optimize() const;
  void Close() const;

  sqlite3_int64 SelectSourceMoleculeID(const std::string& smiles) const;
  sqlite3_int64 SelectSourceMoleculeID(const RDKit::ROMol& molecule) const;
  sqlite3_int64 SelectPseudofragmentID(const std::string& smiles, bool ring_part) const;
  sqlite3_int64 SelectPseudofragmentID(const Pseudofragment& pseudofragment) const;
  sqlite3_int64 SelectConnectionID(std::uint32_t start_atom_type, std::uint32_t end_atom_type, std::uint32_t bond_type) const;
  sqlite3_int64 SelectConnectionID(const Connection& connection) const;
  std::pair<Connection, unsigned> SelectConnectionWithID(sqlite3_int64 connection_id) const;
  std::pair<Pseudofragment, unsigned> SelectPseudofragmentWithID(sqlite3_int64 pseudofragment_id) const;
  sqlite3_int64 GetNPseudofragments() const;
  sqlite3_int64 GetNConnections() const;

  sqlite3_int64 InsertSourceMolecule(const std::string& smiles) const;
  sqlite3_int64 InsertSourceMolecule(const RDKit::ROMol& molecule) const;
  sqlite3_int64 InsertPseudofragment(const Pseudofragment& pseudofragment) const;
  sqlite3_int64 InsertConnection(const Connection& connection, unsigned n) const;
  void InsertSourceMoleculePseudofragmentRelationship(sqlite3_int64 source_molecule_id, sqlite3_int64 pseudofragment_id) const;
  void InsertPseudofragmentConnectionRelationship(sqlite3_int64 pseudofragment_id, sqlite3_int64 connection_id, unsigned frequency) const;

  sqlite3* GetDatabaseConnection() const;

  class PseudofragmentIterator {
	private:
		sqlite3_stmt* select_pseudofragments;
    int sqlite3_result_code = 0;
		sqlite3_int64 pseudofragment_id = 0;
    Pseudofragment pseudofragment;
    unsigned pseudofragment_frequency = 0;

	public:
    PseudofragmentIterator(const PseudofragmentDB& pseudofragment_db, sqlite3_int64 begin_id = 1, sqlite3_int64 end_id = 0);
		const Pseudofragment& operator*() const;
		PseudofragmentIterator& operator++();
		PseudofragmentIterator operator++(int);
		bool AtEnd() const;
    sqlite3_int64 GetPseudofragmentID() const;
    const Pseudofragment& GetPseudofragment() const;
    unsigned GetPseudofragmentFrequency() const;
		void Finalize() const;
	};

  class ConnectionIterator {
  private:
    sqlite3_stmt* select_connections;
    int sqlite3_result_code = 0;
    sqlite3_int64 connection_id = 0;
    Connection connection;
    unsigned connection_frequency = 0;

  public:
    ConnectionIterator(const PseudofragmentDB& pseudofragment_db, sqlite3_int64 begin_id = 1, sqlite3_int64 end_id = 0);
    const Connection& operator*() const;
    ConnectionIterator& operator++();
    ConnectionIterator operator++(int);
    bool AtEnd() const;
    sqlite3_int64 GetConnectionID() const;
    const Connection& GetConnection() const;
    unsigned GetConnectionFrequency() const;
    void Finalize() const;
  };

private:
  void Open();
  bool TablesExist() const;
  void InitializeTables() const;
  void PrepareStatements();
  void FinalizeStatements() const;
};

// Function to determine the number of OpenMP threads to use.
int DetermineOMPNThreads();

// Function to divide a range of size N into equally sized chunks.
std::vector<std::pair<unsigned, unsigned>> EquallySizedChunks(unsigned size, unsigned n_chunks, bool one_based_index = false);

#endif
