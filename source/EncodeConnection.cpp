#include <iostream>

unsigned Encode(std::uint8_t a, std::uint8_t b) {
  return a | (b << 8);
};

int main(int argc, char** argv) {

  unsigned encoded = Encode(std::atoi(argv[1]), std::atoi(argv[2]));
  std::cout << encoded << "\n";

  return 0;
};
