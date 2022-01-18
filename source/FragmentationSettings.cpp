#include "FragmentationSettings.hpp"

FragmentationSettings::FragmentationSettings() = default;
FragmentationSettings::FragmentationSettings(const std::string& settings_file) {
  // Iterate over the lines in the input file.
  std::ifstream input_stream(settings_file);
  std::string line, key, value;
  while (std::getline(input_stream, line)) {
    if (line[0] == '#') {
      continue;
    };
    // Parse the line into key and value.
    std::stringstream line_stream (line);
    line_stream >> key;
    line_stream.get();
    std::getline(line_stream, value);
    // If the value is blank skip the line.
    if (StringIsBlank(value)) {
      continue;
    };
    // Assign the value to the right variable according to the key.
    if (key == "ATOM_TYPING_SCHEME") {
      if (value == "DUMMY") {
        atom_typing = AtomTyping::DUMMY;
      } else if (value == "ATOMIC_NUMBER") {
        atom_typing = AtomTyping::ATOMIC_NUMBER;
      } else if (value == "MMFF") {
        atom_typing = AtomTyping::MMFF;
      } else if (value == "MORGAN") {
        atom_typing = AtomTyping::MORGAN;
      } else if (value == "HASHED_MORGAN") {
        atom_typing = AtomTyping::HASHED_MORGAN;
      } else {
        throw std::runtime_error("ERROR: Invalid ATOM_TYPING_SCHEME.");
      };
    } else if (key == "SYSTEMATIC_FRAGMENTATION_SCHEME") {
      if (value == "NONE") {
        systematic_fragmentation = SystematicFragmentation::NONE;
      } else if (value == "LINEAR") {
        systematic_fragmentation = SystematicFragmentation::LINEAR;
      } else if (value == "SUBGRAPH") {
        systematic_fragmentation = SystematicFragmentation::SUBGRAPH;
      } else if (value == "CIRCULAR") {
        systematic_fragmentation = SystematicFragmentation::CIRCULAR;
      } else {
        throw std::runtime_error("ERROR: Invalid SYSTEMATIC_FRAGMENTATION_SCHEME.");
      };
    } else if (key == "MIN_FRAGMENT_SIZE") {
      min_fragment_size = std::stoi(value);
    } else if (key == "MAX_FRAGMENT_SIZE") {
      max_fragment_size = std::stoi(value);
      assert(max_fragment_size >= min_fragment_size);
    } else if (key == "FRAGMENT_RINGS") {
      fragment_rings = std::stoi(value);
    } else if (key == "MORGAN_RADIUS") {
      morgan_radius = std::stoi(value);
    } else if (key == "MORGAN_CONSIDER_CHIRALITY") {
      morgan_consider_chirality = std::stoi(value);
    } else if (key == "HASHED_MORGAN_N_BITS") {
      hashed_morgan_n_bits = std::stoi(value);
      assert(hashed_morgan_n_bits > 0);
    } else {
      std::stringstream ss;
      ss << "ERROR: Invalid key: '" << key << "'." << std::endl;
      throw std::runtime_error(ss.str());
    };
    line.clear();
    key.clear();
    value.clear();
  };
  input_stream.close();
};

void FragmentationSettings::SetAtomTyping(AtomTyping new_atom_typing) {
  atom_typing = new_atom_typing;
};

void FragmentationSettings::SetSystematicFragmentation(SystematicFragmentation new_systematic_fragmentation) {
  systematic_fragmentation = new_systematic_fragmentation;
};

void FragmentationSettings::SetMinFragmentSize(unsigned new_min_fragment_size) {
  min_fragment_size = new_min_fragment_size;
};

void FragmentationSettings::SetMaxFragmentSize(unsigned new_max_fragment_size) {
  max_fragment_size = new_max_fragment_size;
};

void FragmentationSettings::SetFragmentRings(bool new_fragment_rings) {
  fragment_rings = new_fragment_rings;
};

void FragmentationSettings::SetMorganRadius(unsigned new_morgan_radius) {
  morgan_radius = new_morgan_radius;
};

void FragmentationSettings::SetMorganConsiderChirality(bool new_morgan_consider_chirality) {
  morgan_consider_chirality = new_morgan_consider_chirality;
};

void FragmentationSettings::SetHashedMorganNBits(unsigned new_hashed_morgan_n_bits) {
  hashed_morgan_n_bits = new_hashed_morgan_n_bits;
};

FragmentationSettings::AtomTyping FragmentationSettings::GetAtomTyping() const {
  return atom_typing;
};

FragmentationSettings::SystematicFragmentation FragmentationSettings::GetSystematicFragmentation() const {
  return systematic_fragmentation;
};

unsigned FragmentationSettings::GetMinFragmentSize() const {
  return min_fragment_size;
};

unsigned FragmentationSettings::GetMaxFragmentSize() const {
  return max_fragment_size;
};

bool FragmentationSettings::FragmentRings() const {
  return fragment_rings;
};

unsigned FragmentationSettings::GetMorganRadius() const {
  return morgan_radius;
};

bool FragmentationSettings::MorganConsiderChirality() const {
  return morgan_consider_chirality;
};

unsigned FragmentationSettings::GetHashedMorganNBits() const {
  return hashed_morgan_n_bits;
};

void FragmentationSettings::Print() const {
  std::cout << "ATOM_TYPING_SCHEME: " << static_cast<std::underlying_type<AtomTyping>::type>(atom_typing) << std::endl;
  std::cout << "SYSTEMATIC_FRAGMENTATION_SCHEME: " << static_cast<std::underlying_type<SystematicFragmentation>::type>(systematic_fragmentation) << std::endl;
  std::cout << "MIN_FRAGMENT_SIZE: " << min_fragment_size << std::endl;
  std::cout << "MAX_FRAGMENT_SIZE: " << max_fragment_size << std::endl;
  std::cout << "FRAGMENT_RINGS: " << fragment_rings << std::endl;
  std::cout << "MORGAN_RADIUS: " << morgan_radius << std::endl;
  std::cout << "MORGAN_CONSIDER_CHIRALITY: " << morgan_consider_chirality << std::endl;
  std::cout << "HASHED_MORGAN_N_BITS: " << hashed_morgan_n_bits << std::endl;
};

bool StringIsBlank(const std::string& str) {
  return str.find_first_not_of(" \t\n\v\f\r") == str.npos;
};
