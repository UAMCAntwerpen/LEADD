#pragma once
#ifndef _FRAGMENTATION_HPP_
#define _FRAGMENTATION_HPP_

#include "Pseudofragment.hpp"
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

// Function to flag the RDKit::ROMol's atoms with an immutable index property.
void AssignImmutableAtomIndices(RDKit::ROMol& mol);

// Function to calculate the MMFF94 atom types of a RDKit::ROMol's atoms and store
// them as a property.
void AssignMMFFAtomTypes(RDKit::ROMol& mol);

// Function to flag whether the RDKit::ROMol's atoms are in a ring or not.
void FlagRingAtoms(RDKit::ROMol& mol);

// Function to flag the RDKit::ROMol's atoms with all of the above properties.
void FlagAtoms(RDKit::ROMol& mol);

void CopyFlags(const RDKit::Atom* source, RDKit::Atom* destination);

// Function to recursively expand a ring system by annexing smallest rings that share an atom.
void ExpandRingSystem(std::vector<int>& ring_system, std::list<std::vector<int>>& rings);

// Function to retrieve the ring systems of a molecule, defined as those rings of the SSSR
// that have atleast one overlapping atom.
std::vector<std::vector<int>> DefineRingSystems(const RDKit::ROMol& mol);

// Function to convert a dummy molecule (i.e. a RDKit::ROMol with Connection-encoding pseudoatoms) into a Pseudofragment.
Pseudofragment DummyMolToPseudofragment(const RDKit::ROMol& mol, boost::format& formatter);

// Function to convert an RDKit::Atom into a Pseudofragment.
Pseudofragment AtomToPseudofragment(const RDKit::ROMol& mol, const RDKit::Atom* atom, boost::format& formatter);

// Function to separate the cyclic parts of the molecule from the acyclic ones.
std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> SplitAtRings(const RDKit::ROMol& mol);

#endif
