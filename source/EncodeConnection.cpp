#include <iostream>
#include <bitset>

std::uint64_t EncodeAtomTypes(std::uint32_t start_atom_type, std::uint32_t end_atom_type) {
  std::uint64_t encoded_atom_types = ((std::uint64_t) start_atom_type) << 32 | end_atom_type;
  return encoded_atom_types;
};

int main(int argc, char** argv) {
  std::uint32_t start_atom_type = std::atoi(argv[1]);
  std::uint32_t end_atom_type = std::atoi(argv[2]);
  std::uint64_t encoded_atom_types = EncodeAtomTypes(start_atom_type, end_atom_type);
  std::bitset<32> start_atom_type_bits (start_atom_type);
  std::bitset<32> end_atom_type_bits (end_atom_type);
  std::bitset<64> encoded_atom_types_bits (encoded_atom_types);

  std::cout << "Start atom type:" << std::endl;
  std::cout << start_atom_type << std::endl;
  std::cout << start_atom_type_bits << std::endl;
  std::cout << "End atom type:" << std::endl;
  std::cout << end_atom_type << std::endl;
  std::cout << end_atom_type_bits << std::endl;
  std::cout << "Encoded atom types:" << std::endl;
  std::cout << encoded_atom_types << std::endl;
  std::cout << encoded_atom_types_bits << std::endl;

  return 0;
};
