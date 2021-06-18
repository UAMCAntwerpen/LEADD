#pragma once
#ifndef _FRAGMENTATION_HPP_
#define _FRAGMENTATION_HPP_

#include "Pseudofragment.hpp"
#include "FragmentationSettings.hpp"
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/new_canon.h>

// Recursively expand a ring system by annexing smallest rings that share an atom.
void ExpandRingSystem(std::set<unsigned>& ring_system, std::vector<std::set<unsigned>>& rings);

// Retrieve the ring systems of a molecule, defined as those rings of the SSSR
// that have atleast one overlapping atom.
std::vector<std::set<unsigned>> DefineRingSystems(const RDKit::ROMol& molecule);

// Assigns dummy atom types (i.e. 0) to every atom. During molecule design this
// will correspond to a valence model.
std::vector<std::uint32_t> GetDummyAtomTypes(const RDKit::ROMol& molecule);

// Uses atomic numbers as atom types.
std::vector<std::uint32_t> GetAtomicNumberAtomTypes(const RDKit::ROMol& molecule);

// Calculates the MMFF94 atom types of a molecule.
std::vector<std::uint32_t> GetMMFFAtomTypes(const RDKit::ROMol& molecule);

// Calculates the Morgan atom types of a molecule. For each atom only the one
// corresponding with the specified radius is kept.
std::vector<std::uint32_t> GetMorganAtomTypes(const RDKit::ROMol& molecule, unsigned radius = 2, bool consider_chirality = false);

// Calculates the aforementioned Morgan atom types of a molecule but hashes them
// to a smaller number between 0 and n_bits.
std::vector<std::uint32_t> GetHashedMorganAtomTypes(const RDKit::ROMol& molecule, unsigned radius = 2, unsigned n_bits = 2048, bool consider_chirality = false);

// Calculates the atom types defined by the settings.
std::vector<std::uint32_t> GetAtomTypes(const RDKit::ROMol& molecule, const FragmentationSettings& settings);

// Assign the current atom index as the atom's OldIdx property.
// All functions that create fragments either set or carry over this property.
void AssignOldIdxAtomProperty(RDKit::ROMol& molecule);

// Removes a set of bonds from a molecule without any sort of corrections. We
// use this as a substitute to RDKit::MolFragmenter::fragmentOnBonds because
// the latter corrects the hydrogen count on aromatic heteroatoms.
void RemoveBonds(RDKit::RWMol& molecule, const std::vector<unsigned>& bond_indices);

// Separate the cyclic parts of the molecule from the acyclic ones.
std::pair<std::vector<RDKit::ROMOL_SPTR>, std::vector<RDKit::ROMOL_SPTR>> SplitAtRings(const RDKit::ROMol& molecule, bool copy_old_idx_property = false);

// Extract a substructure from a molecule as another molecule (fragment).
RDKit::ROMOL_SPTR PathToFragment(const RDKit::ROMol& molecule, const std::vector<int>& path, bool copy_old_idx_property = false);

// Extract individual atoms as fragments.
std::vector<RDKit::ROMOL_SPTR> AtomsToFragments(RDKit::ROMol& molecule, bool copy_old_idx_property = false);

// Dummy function. Doesn't do anything but return the input molecule.
std::vector<RDKit::ROMOL_SPTR> NoFragments(RDKit::ROMol& molecule);

// Defines all possible linear paths of a given size range within a molecule and
// converts each one into a fragment.
std::vector<RDKit::ROMOL_SPTR> LinearPathsToFragments(RDKit::ROMol& molecule, unsigned minimum_n_bonds, unsigned maximum_n_bonds, bool copy_old_idx_property = false);

// Defines all possible subgraphs of a given size range within a molecule and
// converts each one into a fragment.
std::vector<RDKit::ROMOL_SPTR> SubgraphsToFragments(RDKit::ROMol& molecule, unsigned minimum_n_bonds, unsigned maximum_n_bonds, bool copy_old_idx_property = false);

// Defines all possible circular environmetns of a given size range within a
// molecule and converts each one into a fragment.
std::vector<RDKit::ROMOL_SPTR> CircularEnvironmentsToFragments(RDKit::ROMol& molecule, unsigned minimum_radius, unsigned maximum_radius, bool copy_old_idx_property = false);

// Systematically defines molecular fragments according to the specified settings.
std::vector<RDKit::ROMOL_SPTR> FragmentSystematically(RDKit::ROMol& molecule, const FragmentationSettings& settings, bool copy_old_idx_property = false);

// Determines which fragments contain rings.
std::vector<bool> RingFragmentMask(const std::vector<RDKit::ROMOL_SPTR>& fragments);

// Determines which fragments contain parts of ring systems.
std::vector<bool> RingPartFragmentMask(const RDKit::ROMol& molecule, const std::vector<RDKit::ROMOL_SPTR>& fragments);

// Determine the order in which atoms will be written to canonical SMILES.
std::vector<unsigned> GetCanonicalSMILESAtomOrder(const RDKit::ROMol& molecule);

// Determine the canonical atom order of a molecule.
std::vector<unsigned> GetCanonicalAtomOrder(const RDKit::ROMol& molecule);

// Renumber a molecule's atoms to follow the canonical atom order.
RDKit::ROMOL_SPTR CanonicalizeAtomOrder(const RDKit::ROMol& molecule);

// Erases all molecular and atomic properties.
void EraseProperties(RDKit::RWMol& molecule);

// Generates the canonical connectivity-encoding CXSMILES of a fragment.
std::string MakePseudomolSMILES(const RDKit::ROMol& pseudomol, const ConnectionsTable& connections);

// Converts a fragment into a Pseudofragment.
Pseudofragment FragmentToPseudofragment(const RDKit::ROMol& molecule, const RDKit::ROMOL_SPTR fragment, std::vector<std::uint32_t> molecule_atom_types, bool has_ring, bool is_ring_part);

std::vector<Pseudofragment> MakePseudofragments(RDKit::ROMol& molecule, const FragmentationSettings& settings);
//
// // Function to flag the RDKit::ROMol's atoms with an immutable index property.
// void AssignImmutableAtomIndices(RDKit::ROMol& mol);
//
// // Function to calculate the MMFF94 atom types of a RDKit::ROMol's atoms and store
// // them as a property.
// void AssignMMFFAtomTypes(RDKit::ROMol& mol);
//
// // Function to flag whether the RDKit::ROMol's atoms are in a ring or not.
// void FlagRingAtoms(RDKit::ROMol& mol);
//
// // Function to flag the RDKit::ROMol's atoms with all of the above properties.
// void FlagAtoms(RDKit::ROMol& mol);
//
// void CopyFlags(const RDKit::Atom* source, RDKit::Atom* destination);
//
// // Function to convert a dummy molecule (i.e. a RDKit::ROMol with Connection-encoding pseudoatoms) into a Pseudofragment.
// Pseudofragment DummyMolToPseudofragment(const RDKit::ROMol& mol, boost::format& formatter);
//
// // Function to convert an RDKit::Atom into a Pseudofragment.
// Pseudofragment AtomToPseudofragment(const RDKit::ROMol& mol, const RDKit::Atom* atom, boost::format& formatter);

#endif
