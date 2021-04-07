#include "PseudofragmentDB.hpp"

sqlite3* InitializePseudofragmentDB(const std::string& output) {
  // Open a connection to the database.
  sqlite3* database;
  sqlite3_open(output.c_str(), &database);

  // Initialize a vector to store the SQL statements. Used purely for cleaner code.
  std::vector<const char*> statements{
    "PRAGMA foreign_keys = ON;",

    "CREATE TABLE sources ("
    "id INTEGER PRIMARY KEY,"
    "smiles TEXT NOT NULL UNIQUE);",

    "CREATE TABLE pseudofragments ("
    "id INTEGER PRIMARY KEY,"
    "smiles TEXT NOT NULL,"
    "ring_part INTEGER NOT NULL,"
    "has_ring INTEGER NOT NULL,"
    "level INTEGER NOT NULL,"
    "pickle BLOB NOT NULL,"
    "pickle_size INTEGER NOT NULL,"
    "frequency INTEGER NOT NULL,"
    "UNIQUE (smiles, ring_part));",
    "CREATE UNIQUE INDEX covering_index_pseudofragments_bools ON pseudofragments(ring_part, has_ring, id, level, frequency);",

    "CREATE TABLE connections ("
    "id INTEGER PRIMARY KEY,"
    "start_atom_type INTEGER NOT NULL,"
    "end_atom_type INTEGER NOT NULL,"
    "bond_type INTEGER NOT NULL,"
    "string TEXT NOT NULL UNIQUE,"
    "frequency INTEGER NOT NULL,"
    "UNIQUE (start_atom_type, end_atom_type, bond_type));",

    "CREATE TABLE sources_pseudofragments ("
    "source_id INTEGER NOT NULL,"
    "pseudofragment_id INTEGER NOT NULL,"
    "frequency INTEGER NOT NULL,"
    "FOREIGN KEY (source_id) REFERENCES sources(id),"
    "FOREIGN KEY (pseudofragment_id) REFERENCES pseudofragments(id),"
    "UNIQUE (source_id, pseudofragment_id));",

    "CREATE TABLE pseudofragments_connections ("
    "pseudofragment_id INTEGER NOT NULL,"
    "connection_id INTEGER NOT NULL,"
    "frequency INTEGER NOT NULL,"
    "FOREIGN KEY (pseudofragment_id) REFERENCES pseudofragments(id),"
    "FOREIGN KEY (connection_id) REFERENCES connections(id),"
    "UNIQUE (connection_id, pseudofragment_id, frequency));"
  };

  // Execute all the stored statements.
  for (auto statement : statements) {
    sqlite3_exec(database, statement, NULL, NULL, NULL);
  };

  return database;
};

sqlite3_int64 InsertSourceMolInDB(std::string smiles, sqlite3* database, sqlite3_stmt* statement) {
  sqlite3_bind_text(statement, 1, smiles.c_str(), -1, NULL);
  sqlite3_step(statement);
  sqlite3_int64 source_id = sqlite3_last_insert_rowid(database);
  sqlite3_clear_bindings(statement);
  sqlite3_reset(statement);
  return source_id;
};

void InsertPseudofragmentInDB(const Pseudofragment& pseudofragment, sqlite3_stmt* statement) {
  // Serialize the pseudofragment as a binary string and calculate its size. The
  // size is stored together with the string itself in the database to make it
  // easier to read it back in.
  std::stringstream blob_stream;
  boost::archive::binary_oarchive archive(blob_stream);
  pseudofragment.save(archive);
  std::string blob = blob_stream.str();
  int blob_size = blob.size();
  sqlite3_bind_text(statement, 1, pseudofragment.GetSMILES().c_str(), -1, NULL);
  sqlite3_bind_int(statement, 2, pseudofragment.IsRingPart());
  sqlite3_bind_int(statement, 3, pseudofragment.HasRing());
  sqlite3_bind_int(statement, 4, pseudofragment.GetLevel());
  sqlite3_bind_blob(statement, 5, blob.c_str(), blob_size, NULL);
  sqlite3_bind_int(statement, 6, blob_size);
  sqlite3_step(statement);
  sqlite3_clear_bindings(statement);
  sqlite3_reset(statement);
};

void InsertConnectionsInDB(const Pseudofragment& pseudofragment, sqlite3_int64 pseudofragment_id, sqlite3_stmt* insert_connection, sqlite3_stmt* select_connection_id, sqlite3_stmt* insert_pseudofragment_connection_relationship) {
  for (auto& connection : pseudofragment.GetConnections()) {
    // Insert the Connection in the database.
    sqlite3_bind_int(insert_connection, 1, connection.first.GetStartAtomType());
    sqlite3_bind_int(insert_connection, 2, connection.first.GetEndAtomType());
    sqlite3_bind_int(insert_connection, 3, connection.first.GetBondType());
    sqlite3_bind_text(insert_connection, 4, connection.first.GetString().c_str(), -1, NULL);
    sqlite3_step(insert_connection);
    sqlite3_clear_bindings(insert_connection);
    sqlite3_reset(insert_connection);
    // Retrieve the Connection's ID in the database. This is necessary since
    // UPSERT statements don't update the "last_rowid" attribute when updating.
    sqlite3_bind_int(select_connection_id, 1, connection.first.GetStartAtomType());
    sqlite3_bind_int(select_connection_id, 2, connection.first.GetEndAtomType());
    sqlite3_bind_int(select_connection_id, 3, connection.first.GetBondType());
    sqlite3_step(select_connection_id);
    sqlite3_int64 connection_id = sqlite3_column_int64(select_connection_id, 0);
    sqlite3_clear_bindings(select_connection_id);
    sqlite3_reset(select_connection_id);
    // Insert the Pseudofragment-Connection relationship in the database.
    sqlite3_bind_int64(insert_pseudofragment_connection_relationship, 1, pseudofragment_id);
    sqlite3_bind_int64(insert_pseudofragment_connection_relationship, 2, connection_id);
    sqlite3_bind_int(insert_pseudofragment_connection_relationship, 3, connection.second.size());
    sqlite3_step(insert_pseudofragment_connection_relationship);
    sqlite3_clear_bindings(insert_pseudofragment_connection_relationship);
    sqlite3_reset(insert_pseudofragment_connection_relationship);
  };
};

sqlite3_int64 InsertSourcePseudofragmentRelationshipInDB(sqlite3_int64 source_id, const Pseudofragment& pseudofragment, sqlite3_stmt* select_statement, sqlite3_stmt* insert_statement) {
  // Determine the pseudofragment ID by running a SELECT query.
  sqlite3_bind_text(select_statement, 1, pseudofragment.GetSMILES().c_str(), -1, NULL);
  sqlite3_bind_int(select_statement, 2, pseudofragment.IsRingPart());
  sqlite3_step(select_statement);
  sqlite3_int64 pseudofragment_id = sqlite3_column_int64(select_statement, 0);
  sqlite3_clear_bindings(select_statement);
  sqlite3_reset(select_statement);
  // Insert the relationship in the joining table.
  sqlite3_bind_int64(insert_statement, 1, source_id);
  sqlite3_bind_int64(insert_statement, 2, pseudofragment_id);
  sqlite3_step(insert_statement);
  sqlite3_clear_bindings(insert_statement);
  sqlite3_reset(insert_statement);
  // Return the ID of the involved Pseudofragment.
  return pseudofragment_id;
};

Pseudofragment GetPseudofragmentByIDFromDB(unsigned id, sqlite3_stmt* select_statement) {
  // Run the SQLite3 query to locate the Pseudofragment with the specified ID.
  sqlite3_bind_int(select_statement, 1, id);
  sqlite3_step(select_statement);
  // Deserialize the binary string to recreate the Pseudofragment object.
  int blob_size = sqlite3_column_int(select_statement, 1);
  std::stringstream blob;
  blob.write((const char*)sqlite3_column_blob(select_statement, 0), blob_size);
  sqlite3_clear_bindings(select_statement);
  sqlite3_reset(select_statement);
  Pseudofragment pseudofragment;
  boost::archive::binary_iarchive archive(blob);
  pseudofragment.load(archive);
  return pseudofragment;
};

unsigned GetPseudofragmentCount(sqlite3* database) {
  sqlite3_stmt* statement;
  sqlite3_prepare_v2(database, "SELECT COUNT(*) FROM pseudofragments", -1, &statement, NULL);
  sqlite3_step(statement);
  return sqlite3_column_int64(statement, 0);
};


// Class PseudofragmentIterator
PseudofragmentIterator::PseudofragmentIterator(sqlite3* database) {
  sqlite3_prepare_v2(database, "SELECT id, pickle, pickle_size FROM pseudofragments", -1, &statement, NULL);
  result_code = sqlite3_step(statement);
  if (result_code == SQLITE_DONE) {
    return;
  } else {
    id = sqlite3_column_int64(statement, 0);
    blob_size = sqlite3_column_int(statement, 2);
    std::stringstream blob;
    // IMPROVE: static cast instead of C cast
    blob.write((const char*)sqlite3_column_blob(statement, 1), blob_size);
    boost::archive::binary_iarchive archive(blob);
    pseudofragment.load(archive);
  };
};

PseudofragmentIterator::PseudofragmentIterator(sqlite3* database, unsigned start_idx, unsigned end_idx) {
  sqlite3_prepare_v2(database, "SELECT id, pickle, pickle_size FROM pseudofragments WHERE id BETWEEN ? AND ?", -1, &statement, NULL);
  sqlite3_bind_int64(statement, 1, start_idx);
  sqlite3_bind_int64(statement, 2, end_idx);
  result_code = sqlite3_step(statement);
  if (result_code == SQLITE_DONE) {
    return;
  } else {
    id = sqlite3_column_int64(statement, 0);
    blob_size = sqlite3_column_int(statement, 2);
    std::stringstream blob;
    // IMPROVE: static cast instead of C cast
    blob.write((const char*)sqlite3_column_blob(statement, 1), blob_size);
    boost::archive::binary_iarchive archive(blob);
    pseudofragment.load(archive);
  };
};

const Pseudofragment& PseudofragmentIterator::operator*() const {
  assert(result_code != SQLITE_DONE);
  return pseudofragment;
};

const Pseudofragment* PseudofragmentIterator::operator->() const {
  assert(result_code != SQLITE_DONE);
  return &pseudofragment;
};

void PseudofragmentIterator::operator++() {
  assert(result_code == SQLITE_ROW);
  result_code = sqlite3_step(statement);
  if (result_code == SQLITE_DONE) {
    return;
  };
  id = sqlite3_column_int64(statement, 0);
  blob_size = sqlite3_column_int(statement, 2);
  blob.str(std::string());
  blob.write((const char*)sqlite3_column_blob(statement, 1), blob_size);
  pseudofragment = Pseudofragment();
  boost::archive::binary_iarchive archive(blob);
  pseudofragment.load(archive);
};

unsigned PseudofragmentIterator::GetID() const {
  return id;
};

int PseudofragmentIterator::GetResultCode() const {
  return result_code;
};

void PseudofragmentIterator::Reset() {
  sqlite3_reset(statement);
  result_code = SQLITE_ROW;
  operator++();
};
