#include <bitset>
#include <iostream>
#include <sstream>

int main() {
  int a, b, c, x, u;

  u = 0x3FFFFFFF;
  a = 0x386;
  b = 0x440190;
  c = 0x1E03;
  x = 0;

  x = (a & ~c) & ~b;

  std::cout << std::bitset<26>(x) << std::endl;

  return 0;
}