#include "PseudofragmentDB.hpp"

PseudofragmentDB::PseudofragmentDB(const std::string& database_path) :
  database_path(database_path) {
  // Open a connection to the SQLite3 database.
  Open();
  // If necessary, initialize the tables.
  if (!TablesExist()) {
    InitializeTables();
  };
  // Prepare all the required SQLite3 statements.
  PrepareStatements();
};

void PseudofragmentDB::Open() {
  int sqlite3_result_code = sqlite3_open(database_path.c_str(), &database);
  assert(sqlite3_result_code == SQLITE_OK);
};

bool PseudofragmentDB::TablesExist() const {
  // Prepare the statement to check if a table exists.
  sqlite3_stmt* statement;
  int sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name= ?", -1, &statement, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  // Check if the tables exist.
  std::vector<std::string> table_names {"sources", "pseudofragments", "connections", "sources_pseudofragments", "pseudofragments_connections"};
  for (const std::string& table_name : table_names) {
    sqlite3_bind_text(statement, 1, table_name.c_str(), -1, nullptr);
    sqlite3_result_code = sqlite3_step(statement);
    assert(sqlite3_result_code == SQLITE_ROW);
    sqlite3_int64 count = sqlite3_column_int64(statement, 0);
    if (count == 0) {
      sqlite3_finalize(statement);
      return false;
    };
    sqlite3_clear_bindings(statement);
    sqlite3_reset(statement);
  };
  sqlite3_finalize(statement);
  return true;
};

void PseudofragmentDB::InitializeTables() const {
  // Initialize a vector to store the SQL statements. Used purely for cleaner code.
  std::vector<const char*> statements{
    "PRAGMA foreign_keys = ON;",

    "CREATE TABLE sources ("
    "id INTEGER PRIMARY KEY,"
    "smiles TEXT NOT NULL,"
    "frequency INTEGER NOT NULL,"
    "UNIQUE (smiles));",

    "CREATE TABLE pseudofragments ("
    "id INTEGER PRIMARY KEY,"
    "smiles TEXT NOT NULL,"
    "ring_part INTEGER NOT NULL,"
    "has_ring INTEGER NOT NULL,"
    "size INTEGER NOT NULL,"
    "pickle BLOB NOT NULL,"
    "frequency INTEGER NOT NULL,"
    "UNIQUE (smiles, ring_part));",
    "CREATE UNIQUE INDEX covering_index_pseudofragments_bools ON pseudofragments(ring_part, has_ring, id, size, frequency);",

    "CREATE TABLE connections ("
    "id INTEGER PRIMARY KEY,"
    "start_atom_type INTEGER NOT NULL,"
    "end_atom_type INTEGER NOT NULL,"
    "bond_type INTEGER NOT NULL,"
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
  for (const char* statement : statements) {
    int sqlite3_result_code = sqlite3_exec(database, statement, nullptr, nullptr, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
  };
};

void PseudofragmentDB::PrepareStatements() {
  int sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT COUNT(*) FROM pseudofragments", -1, &select_n_pseudofragments, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT COUNT(*) FROM connections", -1, &select_n_connections, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT id FROM sources WHERE smiles = ?", -1, &select_source_molecule_id, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT id FROM pseudofragments WHERE smiles = ? AND ring_part = ?", -1, &select_pseudofragment_id, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT id FROM connections WHERE start_atom_type = ? AND end_atom_type = ? AND bond_type = ?", -1, &select_connection_id, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT start_atom_type, end_atom_type, bond_type, frequency FROM connections WHERE id = ?", -1, &select_connection, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "SELECT pickle, frequency FROM pseudofragments WHERE id = ?", -1, &select_pseudofragment, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "INSERT INTO sources (smiles, frequency) VALUES (?, 1) ON CONFLICT (smiles) DO UPDATE SET frequency=frequency+1", -1, &insert_source_molecule, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "INSERT INTO pseudofragments (smiles, ring_part, has_ring, size, pickle, frequency) VALUES (?, ?, ?, ?, ?, 1) ON CONFLICT (smiles, ring_part) DO UPDATE SET frequency=frequency+1", -1, &insert_pseudofragment, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "INSERT INTO connections (start_atom_type, end_atom_type, bond_type, frequency) VALUES (?, ?, ?, 1) ON CONFLICT (start_atom_type, end_atom_type, bond_type) DO UPDATE SET frequency=frequency+?", -1, &insert_connection, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "INSERT INTO sources_pseudofragments (source_id, pseudofragment_id, frequency) VALUES (?, ?, 1) ON CONFLICT (source_id, pseudofragment_id) DO UPDATE SET frequency=frequency+1", -1, &insert_source_pseudofragment_relationship, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_prepare_v2(database, "INSERT OR IGNORE INTO pseudofragments_connections (pseudofragment_id, connection_id, frequency) VALUES (?, ?, ?)", -1, &insert_pseudofragment_connection_relationship, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
};

void PseudofragmentDB::FinalizeStatements() const {
  int sqlite3_result_code = sqlite3_finalize(select_n_pseudofragments);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(select_n_connections);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(select_source_molecule_id);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(select_pseudofragment_id);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(select_connection_id);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(select_connection);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(select_pseudofragment);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(insert_source_molecule);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(insert_pseudofragment);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(insert_connection);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(insert_source_pseudofragment_relationship);
  assert(sqlite3_result_code == SQLITE_OK);
  sqlite3_result_code = sqlite3_finalize(insert_pseudofragment_connection_relationship);
  assert(sqlite3_result_code == SQLITE_OK);
};

void PseudofragmentDB::BeginTransaction() const {
  int sqlite3_result_code = sqlite3_exec(database, "BEGIN TRANSACTION", nullptr, nullptr, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
};

void PseudofragmentDB::CommitTransaction() const {
  int sqlite3_result_code = sqlite3_exec(database, "COMMIT TRANSACTION", nullptr, nullptr, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
};

void PseudofragmentDB::Optimize() const {
  // Gather and store statistics on future query time to aid the query planner.
  int sqlite3_result_code = sqlite3_exec(database, "PRAGMA optimize", nullptr, nullptr, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
  // Vacuum the database to combat memory fragmentation and decrease size overhead.
  sqlite3_result_code = sqlite3_exec(database, "VACUUM", nullptr, nullptr, nullptr);
  assert(sqlite3_result_code == SQLITE_OK);
};

void PseudofragmentDB::Close() const {
  // Destroy the prepared statements.
  FinalizeStatements();
  // Close the database.
  int sqlite3_result_code = sqlite3_close(database);
  assert(sqlite3_result_code == SQLITE_OK);
};

sqlite3_int64 PseudofragmentDB::SelectSourceMoleculeID(const std::string& smiles) const {
  sqlite3_bind_text(select_source_molecule_id, 1, smiles.c_str(), -1, nullptr);
  int sqlite3_result_code = sqlite3_step(select_source_molecule_id);
  assert(sqlite3_result_code == SQLITE_ROW);
  sqlite3_int64 source_molecule_id = sqlite3_column_int64(select_source_molecule_id, 0);
  sqlite3_clear_bindings(select_source_molecule_id);
  sqlite3_reset(select_source_molecule_id);
  return source_molecule_id;
};

sqlite3_int64 PseudofragmentDB::SelectSourceMoleculeID(const RDKit::ROMol& molecule) const {
  return SelectSourceMoleculeID(RDKit::MolToSmiles(molecule));
};

sqlite3_int64 PseudofragmentDB::SelectPseudofragmentID(const std::string& smiles, bool ring_part) const {
  sqlite3_bind_text(select_pseudofragment_id, 1, smiles.c_str(), -1, nullptr);
  sqlite3_bind_int(select_pseudofragment_id, 2, ring_part);
  int sqlite3_result_code = sqlite3_step(select_pseudofragment_id);
  assert(sqlite3_result_code == SQLITE_ROW);
  sqlite3_int64 pseudofragment_id = sqlite3_column_int64(select_pseudofragment_id, 0);
  sqlite3_clear_bindings(select_pseudofragment_id);
  sqlite3_reset(select_pseudofragment_id);
  return pseudofragment_id;
};

sqlite3_int64 PseudofragmentDB::SelectPseudofragmentID(const Pseudofragment& pseudofragment) const {
  return SelectPseudofragmentID(pseudofragment.GetConnectionEncodingSMILES(), pseudofragment.IsRingPart());
};

sqlite3_int64 PseudofragmentDB::SelectConnectionID(std::uint32_t start_atom_type, std::uint32_t end_atom_type, std::uint32_t bond_type) const {
  sqlite3_bind_int64(select_connection_id, 1, start_atom_type);
  sqlite3_bind_int64(select_connection_id, 2, end_atom_type);
  sqlite3_bind_int64(select_connection_id, 3, bond_type);
  int sqlite3_result_code = sqlite3_step(select_connection_id);
  assert(sqlite3_result_code == SQLITE_ROW);
  sqlite3_int64 connection_id = sqlite3_column_int64(select_connection_id, 0);
  sqlite3_clear_bindings(select_connection_id);
  sqlite3_reset(select_connection_id);
  return connection_id;
};

sqlite3_int64 PseudofragmentDB::SelectConnectionID(const Connection& connection) const {
  return SelectConnectionID(connection.GetStartAtomType(), connection.GetEndAtomType(), connection.GetBondType());
};

std::pair<Pseudofragment, unsigned> PseudofragmentDB::SelectPseudofragmentWithID(sqlite3_int64 pseudofragment_id) const {
  sqlite3_bind_int64(select_pseudofragment, 1, pseudofragment_id);
  int sqlite3_result_code = sqlite3_step(select_pseudofragment);
  assert(sqlite3_result_code == SQLITE_ROW);
  // De-serialize the Pseudofragment.
  std::stringstream blob;
  const char* blob_buffer = (const char*) sqlite3_column_blob(select_pseudofragment, 0);
  int blob_size = sqlite3_column_bytes(select_pseudofragment, 0);
  blob.write(blob_buffer, blob_size);
  boost::archive::binary_iarchive archive(blob);
  Pseudofragment pseudofragment;
  pseudofragment.load(archive);
  // Retrieve its frequency.
  unsigned frequency = sqlite3_column_int64(select_pseudofragment, 1);
  sqlite3_clear_bindings(select_pseudofragment);
  sqlite3_reset(select_pseudofragment);
  // Return both values.
  return {pseudofragment, frequency};
};

std::pair<Connection, unsigned> PseudofragmentDB::SelectConnectionWithID(sqlite3_int64 connection_id) const {
  sqlite3_bind_int64(select_connection, 1, connection_id);
  int sqlite3_result_code = sqlite3_step(select_connection);
  assert(sqlite3_result_code == SQLITE_ROW);
  std::uint32_t start_atom_type = sqlite3_column_int64(select_connection, 0);
  std::uint32_t end_atom_type = sqlite3_column_int64(select_connection, 1);
  std::uint32_t bond_type = sqlite3_column_int64(select_connection, 2);
  unsigned frequency = sqlite3_column_int64(select_connection, 3);
  Connection connection (start_atom_type, end_atom_type, bond_type);
  sqlite3_clear_bindings(select_connection);
  sqlite3_reset(select_connection);
  return {connection, frequency};
};

sqlite3_int64 PseudofragmentDB::GetNPseudofragments() const {
  int sqlite3_result_code = sqlite3_step(select_n_pseudofragments);
  assert(sqlite3_result_code == SQLITE_ROW);
  sqlite3_int64 n_pseudofragments = sqlite3_column_int64(select_n_pseudofragments, 0);
  sqlite3_reset(select_n_pseudofragments);
  return n_pseudofragments;
};

sqlite3_int64 PseudofragmentDB::GetNConnections() const {
  int sqlite3_result_code = sqlite3_step(select_n_connections);
  assert(sqlite3_result_code == SQLITE_ROW);
  sqlite3_int64 n_connections = sqlite3_column_int64(select_n_connections, 0);
  sqlite3_reset(select_n_connections);
  return n_connections;
};

sqlite3_int64 PseudofragmentDB::InsertSourceMolecule(const std::string& smiles) const {
  // Insert the molecule in the database.
  sqlite3_bind_text(insert_source_molecule, 1, smiles.c_str(), -1, nullptr);
  int sqlite3_result_code = sqlite3_step(insert_source_molecule);
  assert(sqlite3_result_code == SQLITE_DONE);
  sqlite3_clear_bindings(insert_source_molecule);
  sqlite3_reset(insert_source_molecule);
  // Retrieve the ID of the inserted molecule.
  sqlite3_int64 source_molecule_id = SelectSourceMoleculeID(smiles);
  return source_molecule_id;
};

sqlite3_int64 PseudofragmentDB::InsertSourceMolecule(const RDKit::ROMol& molecule) const {
  return InsertSourceMolecule(RDKit::MolToSmiles(molecule));
};

sqlite3_int64 PseudofragmentDB::InsertPseudofragment(const Pseudofragment& pseudofragment) const {
  // Serialize the Pseudofragment.
  std::stringstream blob_stream;
  boost::archive::binary_oarchive archive(blob_stream);
  pseudofragment.save(archive);
  std::string blob = blob_stream.str();
  // Insert the Pseudofragment in the database.
  std::string smiles = pseudofragment.GetConnectionEncodingSMILES();
  sqlite3_bind_text(insert_pseudofragment, 1, smiles.c_str(), -1, nullptr);
  sqlite3_bind_int(insert_pseudofragment, 2, pseudofragment.IsRingPart());
  sqlite3_bind_int(insert_pseudofragment, 3, pseudofragment.HasRing());
  sqlite3_bind_int(insert_pseudofragment, 4, pseudofragment.GetSize());
  sqlite3_bind_blob(insert_pseudofragment, 5, blob.c_str(), blob.size(), nullptr);
  int sqlite3_result_code = sqlite3_step(insert_pseudofragment);
  assert(sqlite3_result_code == SQLITE_DONE);
  sqlite3_clear_bindings(insert_pseudofragment);
  sqlite3_reset(insert_pseudofragment);
  // Retrieve the Pseudofragment's ID.
  sqlite3_int64 pseudofragment_id = SelectPseudofragmentID(smiles, pseudofragment.IsRingPart());
  return pseudofragment_id;
};

sqlite3_int64 PseudofragmentDB::InsertConnection(const Connection& connection, unsigned n) const {
  // Insert the Connection in the database.
  sqlite3_bind_int64(insert_connection, 1, connection.GetStartAtomType());
  sqlite3_bind_int64(insert_connection, 2, connection.GetEndAtomType());
  sqlite3_bind_int64(insert_connection, 3, connection.GetBondType());
  sqlite3_bind_int64(insert_connection, 4, n);
  int sqlite3_result_code = sqlite3_step(insert_connection);
  assert(sqlite3_result_code == SQLITE_DONE);
  sqlite3_clear_bindings(insert_connection);
  sqlite3_reset(insert_connection);
  // Retrieve the Connection's ID.
  sqlite3_int64 connection_id = SelectConnectionID(connection.GetStartAtomType(), connection.GetEndAtomType(), connection.GetBondType());
  return connection_id;
};

void PseudofragmentDB::InsertSourceMoleculePseudofragmentRelationship(sqlite3_int64 source_molecule_id, sqlite3_int64 pseudofragment_id) const {
  sqlite3_bind_int64(insert_source_pseudofragment_relationship, 1, source_molecule_id);
  sqlite3_bind_int64(insert_source_pseudofragment_relationship, 2, pseudofragment_id);
  int sqlite3_result_code = sqlite3_step(insert_source_pseudofragment_relationship);
  assert(sqlite3_result_code == SQLITE_DONE);
  sqlite3_clear_bindings(insert_source_pseudofragment_relationship);
  sqlite3_reset(insert_source_pseudofragment_relationship);
};

void PseudofragmentDB::InsertPseudofragmentConnectionRelationship(sqlite3_int64 pseudofragment_id, sqlite3_int64 connection_id, unsigned frequency) const {
  sqlite3_bind_int64(insert_pseudofragment_connection_relationship, 1, pseudofragment_id);
  sqlite3_bind_int64(insert_pseudofragment_connection_relationship, 2, connection_id);
  sqlite3_bind_int64(insert_pseudofragment_connection_relationship, 3, frequency);
  int sqlite3_result_code = sqlite3_step(insert_pseudofragment_connection_relationship);
  assert(sqlite3_result_code == SQLITE_DONE);
  sqlite3_clear_bindings(insert_pseudofragment_connection_relationship);
  sqlite3_reset(insert_pseudofragment_connection_relationship);
};

sqlite3* PseudofragmentDB::GetDatabaseConnection() const {
  return database;
};

PseudofragmentDB::PseudofragmentIterator::PseudofragmentIterator(const PseudofragmentDB& pseudofragment_db, sqlite3_int64 begin_id, sqlite3_int64 end_id) {
  if (end_id > 0) {
    assert(end_id >= begin_id);
    sqlite3_result_code = sqlite3_prepare_v2(pseudofragment_db.database, "SELECT id, pickle, frequency FROM pseudofragments WHERE id BETWEEN ? AND ?", -1, &select_pseudofragments, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_bind_int64(select_pseudofragments, 1, begin_id);
    sqlite3_bind_int64(select_pseudofragments, 2, end_id);
  } else {
    sqlite3_result_code = sqlite3_prepare_v2(pseudofragment_db.database, "SELECT id, pickle, frequency FROM pseudofragments WHERE id >= ?", -1, &select_pseudofragments, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_bind_int64(select_pseudofragments, 1, begin_id);
  };
	++(*this);
};

const Pseudofragment& PseudofragmentDB::PseudofragmentIterator::operator*() const {
  return pseudofragment;
};

PseudofragmentDB::PseudofragmentIterator& PseudofragmentDB::PseudofragmentIterator::operator++() {
	sqlite3_result_code = sqlite3_step(select_pseudofragments);
  if (!AtEnd()) {
    // Get the Pseudofragment's database properties.
    pseudofragment_id = sqlite3_column_int64(select_pseudofragments, 0);
    pseudofragment_frequency = sqlite3_column_int64(select_pseudofragments, 2);
    // De-serialize the Pseudofragment.
    std::stringstream blob;
    const char* blob_buffer = (const char*) sqlite3_column_blob(select_pseudofragments, 1);
    int blob_size = sqlite3_column_bytes(select_pseudofragments, 1);
    blob.write(blob_buffer, blob_size);
    boost::archive::binary_iarchive archive(blob);
    Pseudofragment tmp_pseudofragment;
    tmp_pseudofragment.load(archive);
    pseudofragment = tmp_pseudofragment;
  };
	return *this;
};

PseudofragmentDB::PseudofragmentIterator PseudofragmentDB::PseudofragmentIterator::operator++(int) {
	PseudofragmentIterator previous_state = *this;
	++(*this);
	return previous_state;
};

bool PseudofragmentDB::PseudofragmentIterator::AtEnd() const {
	return sqlite3_result_code == SQLITE_DONE;
};

sqlite3_int64 PseudofragmentDB::PseudofragmentIterator::GetPseudofragmentID() const {
  return pseudofragment_id;
};

const Pseudofragment& PseudofragmentDB::PseudofragmentIterator::GetPseudofragment() const {
  return pseudofragment;
};

unsigned PseudofragmentDB::PseudofragmentIterator::GetPseudofragmentFrequency() const {
  return pseudofragment_frequency;
};

void PseudofragmentDB::PseudofragmentIterator::Finalize() const {
	int sqlite3_result_code = sqlite3_finalize(select_pseudofragments);
  assert(sqlite3_result_code == SQLITE_OK);
};

PseudofragmentDB::ConnectionIterator::ConnectionIterator(const PseudofragmentDB& pseudofragment_db, sqlite3_int64 begin_id, sqlite3_int64 end_id) {
  if (end_id > 0) {
    assert(end_id >= begin_id);
    sqlite3_result_code = sqlite3_prepare_v2(pseudofragment_db.database, "SELECT id, start_atom_type, end_atom_type, bond_type, frequency FROM connections WHERE id BETWEEN ? and ?", -1, &select_connections, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_bind_int64(select_connections, 1, begin_id);
    sqlite3_bind_int64(select_connections, 2, end_id);
  } else {
    sqlite3_result_code = sqlite3_prepare_v2(pseudofragment_db.database, "SELECT id, start_atom_type, end_atom_type, bond_type, frequency FROM connections WHERE id >= ?", -1, &select_connections, nullptr);
    assert(sqlite3_result_code == SQLITE_OK);
    sqlite3_bind_int64(select_connections, 1, begin_id);
  };
  ++(*this);
};

const Connection& PseudofragmentDB::ConnectionIterator::operator*() const {
  return connection;
};

PseudofragmentDB::ConnectionIterator& PseudofragmentDB::ConnectionIterator::operator++() {
	sqlite3_result_code = sqlite3_step(select_connections);
  if (!AtEnd()) {
    connection_id = sqlite3_column_int64(select_connections, 0);
    connection_frequency = sqlite3_column_int64(select_connections, 4);
    std::uint32_t start_atom_type = sqlite3_column_int64(select_connections, 1);
    std::uint32_t end_atom_type = sqlite3_column_int64(select_connections, 2);
    std::uint32_t bond_type = sqlite3_column_int64(select_connections, 3);
    connection = Connection(start_atom_type, end_atom_type, bond_type);
  };
	return *this;
};

PseudofragmentDB::ConnectionIterator PseudofragmentDB::ConnectionIterator::operator++(int) {
	ConnectionIterator previous_state = *this;
	++(*this);
	return previous_state;
};

bool PseudofragmentDB::ConnectionIterator::AtEnd() const {
	return sqlite3_result_code == SQLITE_DONE;
};

sqlite3_int64 PseudofragmentDB::ConnectionIterator::GetConnectionID() const {
  return connection_id;
};

const Connection& PseudofragmentDB::ConnectionIterator::GetConnection() const {
  return connection;
};

unsigned PseudofragmentDB::ConnectionIterator::GetConnectionFrequency() const {
  return connection_frequency;
};

void PseudofragmentDB::ConnectionIterator::Finalize() const {
	int sqlite3_result_code = sqlite3_finalize(select_connections);
  assert(sqlite3_result_code == SQLITE_OK);
};

int DetermineOMPNThreads() {
  int n_threads = 0;
  char* omp_num_threads = std::getenv("OMP_NUM_THREADS");
  if (omp_num_threads == nullptr) {
    n_threads = omp_get_max_threads();
  } else {
    n_threads = std::stoi(omp_num_threads);
  };
  return n_threads;
};

std::vector<std::pair<unsigned, unsigned>> EquallySizedChunks(unsigned size, unsigned n_chunks, bool one_based_index) {
  std::vector<std::pair<unsigned, unsigned>> chunks;
  unsigned chunk_size = size / n_chunks;
  // Define up to the last chunk.
  unsigned begin_idx = 0, end_idx = 0;
  if (one_based_index) {
    ++begin_idx;
  };
  for (unsigned n = 0; n < n_chunks - 1; ++n) {
    end_idx = begin_idx + chunk_size - 1;
    chunks.push_back(std::make_pair(begin_idx, end_idx));
    begin_idx = end_idx + 1;
  };
  // Define the last chunk.
  end_idx = size;
  if (!one_based_index) {
    --end_idx;
  };
  chunks.push_back({begin_idx, end_idx});
  return chunks;
};
