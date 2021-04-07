#include <iostream>
#include <array>

std::array<std::uint8_t, 2> Decode(unsigned encoded) {
  return std::array<std::uint8_t, 2> {encoded & 0xFF, (encoded >> 8) & 0xFF};
};

int main(int argc, char** argv) {
  std::array<std::uint8_t, 2> decoded = Decode(std::atoi(argv[1]));

  for (std::uint8_t i : decoded) {
    std::cout << (int)i << " ";
  };
  std::cout << "\n";

  return 0;
};
