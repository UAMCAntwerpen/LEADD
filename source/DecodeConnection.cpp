#include <iostream>
#include <tuple>
#include <bitset>

std::pair<std::uint32_t, std::uint32_t> DecodeAtomTypes(std::uint64_t encoded_atom_types) {
  std::uint32_t start_atom_type = (std::uint32_t) (encoded_atom_types >> 32);
  std::uint32_t end_atom_type = (std::uint32_t) encoded_atom_types;
  return std::pair<std::uint32_t, std::uint32_t> (start_atom_type, end_atom_type);
};

int main(int argc, char** argv) {
  std::uint64_t encoded_atom_types = std::stoull(argv[1]);
  std::uint32_t start_atom_type, end_atom_type;
  std::tie(start_atom_type, end_atom_type) = DecodeAtomTypes(encoded_atom_types);
  std::bitset<32> start_atom_type_bits (start_atom_type);
  std::bitset<32> end_atom_type_bits (end_atom_type);
  std::bitset<64> encoded_atom_types_bits (encoded_atom_types);

  std::cout << "Encoded atom types:" << std::endl;
  std::cout << encoded_atom_types << std::endl;
  std::cout << encoded_atom_types_bits << std::endl;
  std::cout << "Start atom type:" << std::endl;
  std::cout << start_atom_type << std::endl;
  std::cout << start_atom_type_bits << std::endl;
  std::cout << "End atom type:" << std::endl;
  std::cout << end_atom_type << std::endl;
  std::cout << end_atom_type_bits << std::endl;

  return 0;
};
