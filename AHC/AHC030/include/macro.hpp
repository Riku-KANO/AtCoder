#pragma once

using ll = long long;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using byte = unsigned char;
#define FOR(i, s, n) for (int i = s; i < n; i++)
#define Vec std::vector

typedef int Response;
typedef uint64_t Hash;

template <typename T>
concept Integer = std::is_same_v<T, i16> ||
    std::is_same_v<T, i32> ||
    std::is_same_v<T, i64> ||
    std::is_same_v<T, u16> ||
    std::is_same_v<T, u32> ||
    std::is_same_v<T, u64>;