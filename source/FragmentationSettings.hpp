#pragma once
#ifndef _FRAGMENTATION_SETTINGS_HPP_
#define _FRAGMENTATION_SETTINGS_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>

class FragmentationSettings {
public:
  enum class AtomTyping {DUMMY, ATOMIC_NUMBER, MMFF, MORGAN, HASHED_MORGAN};
  enum class SystematicFragmentation {NONE, LINEAR, SUBGRAPH, CIRCULAR};

private:
  AtomTyping atom_typing = AtomTyping::MORGAN;
  SystematicFragmentation systematic_fragmentation = SystematicFragmentation::NONE;

  unsigned min_fragment_size = 0, max_fragment_size = 0;

  bool fragment_rings = false;

  unsigned morgan_radius = 1;
  bool morgan_consider_chirality = false;
  unsigned hashed_morgan_n_bits = 2048;

public:
  FragmentationSettings();
  FragmentationSettings(const std::string& settings_file);

  void SetAtomTyping(AtomTyping new_atom_typing);
  void SetSystematicFragmentation(SystematicFragmentation new_systematic_fragmentation);
  void SetMinFragmentSize(unsigned new_min_fragment_size);
  void SetMaxFragmentSize(unsigned new_max_fragment_size);
  void SetFragmentRings(bool new_fragment_rings);
  void SetMorganRadius(unsigned new_morgan_radius);
  void SetMorganConsiderChirality(bool new_morgan_consider_chirality);
  void SetHashedMorganNBits(unsigned new_hashed_morgan_n_bits);

  AtomTyping GetAtomTyping() const;
  SystematicFragmentation GetSystematicFragmentation() const;
  unsigned GetMinFragmentSize() const;
  unsigned GetMaxFragmentSize() const;
  bool FragmentRings() const;
  unsigned GetMorganRadius() const;
  bool MorganConsiderChirality() const;
  unsigned GetHashedMorganNBits() const;

  void Print() const;
};

bool StringIsBlank(const std::string& str);

#endif
