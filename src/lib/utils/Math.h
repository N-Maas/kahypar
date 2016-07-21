/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_UTILS_MATH_H_
#define SRC_LIB_UTILS_MATH_H_

namespace utils {
// See: Warren, Henry S. Hacker's Delight (Second Edition), p.61
template <typename T>
static inline T nextPowerOfTwoCeiled(T x) {
  return 0x80000000 >> (__builtin_clz(x - 1) - 1);
}

template <typename T>
static inline T _rol(const T x, const uint r) {
  return (x << r) | (x >> (sizeof(T) * 8 - r));
}

template <typename T>
static inline T _rol(const T x) {
  return _rol(x, x & (sizeof(T) * 8 - 1));
}
}  // namespace utils
#endif  // SRC_LIB_UTILS_MATH_H_
