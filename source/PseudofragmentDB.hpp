#pragma once
#ifndef _PSEUDOFRAGMENT_DB_HPP_
#define _PSEUDOFRAGMENT_DB_HPP_

#include <sqlite3.h>
#include "Pseudofragment.hpp"

// Function to initialize a Pseudofragment database with the correct schema.
sqlite3* InitializePseudofragmentDB(const std::string& output);

// Function to insert a source molecule (as a SMILES string) into the database.
// Returns the database ID of the inserted source molecule.
sqlite3_int64 InsertSourceMolInDB(std::string smiles, sqlite3* database, sqlite3_stmt* statement);

// Function to insert a Pseudofragment into the database.
void InsertPseudofragmentInDB(const Pseudofragment& pseudofragment, sqlite3_stmt* statement);

// Function to insert a Connection into the database.
void InsertConnectionsInDB(const Pseudofragment& pseudofragment, sqlite3_int64 pseudofragment_id, sqlite3_stmt* insert_connection, sqlite3_stmt* select_connection_id, sqlite3_stmt* insert_pseudofragment_connection_relationship);

// Function to insert a source molecule-Pseudofragment relationship into the database.
// Returns the database ID of the Pseudofragment in question.
sqlite3_int64 InsertSourcePseudofragmentRelationshipInDB(sqlite3_int64 source_id, const Pseudofragment& pseudofragment, sqlite3_stmt* select_statement, sqlite3_stmt* insert_statement);

// Function to retrieve a Pseudofragment from the database by its database ID.
Pseudofragment GetPseudofragmentByIDFromDB(unsigned id, sqlite3_stmt* select_statement);

// Function to retrieve the total number of Pseudofragments in the database.
unsigned GetPseudofragmentCount(sqlite3* database);

// Class to iterate over a Pseudofragment database and deserialize the
// Pseudofragments.
class PseudofragmentIterator {
  // SQLite3 statement used to create the iterator.
  sqlite3_stmt* statement;
  // Binary string to store the serialized Pseudofragment.
  std::stringstream blob;
  // Size of the binary string.
  int blob_size = 0;
  // Last accessed Pseudofragment
  Pseudofragment pseudofragment;
  // ID of the last accessed Pseudofragment.
  unsigned id = 0;
  // Integer recording the last SQLite3 result code.
  int result_code;

public:
  PseudofragmentIterator(sqlite3* database);
  PseudofragmentIterator(sqlite3* database, unsigned start_idx, unsigned end_idx);

  // Iterator dereferencing operators.
  const Pseudofragment& operator*() const;
  const Pseudofragment* operator->() const;

  // Operator to advance the iterator by 1.
  void operator++();

  unsigned GetID() const;
  int GetResultCode() const;

  // Function to reset the iterator to its starting position.
  void Reset();
};

#endif
