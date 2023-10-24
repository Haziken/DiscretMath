#pragma once

#include <bitset>

namespace dml {
template <typename T, size_t size> class Array {
public:
  Array();

private:
  std::bitset<size> m_bits;
};
} // namespace dml