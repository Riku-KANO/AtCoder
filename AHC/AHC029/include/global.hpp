#pragma once
#include <iostream>
#include <random>
#include <chrono>
#include <utility>
#include <cmath>
#include <string>

#define pii std::pair<int, int>
#define Vec std::vector
#define FOR(i, n) for(int i = 0; i < (int)(n); i++)
#define MONEY long long
#define LOG_INFO(message, ...) fprintf(stderr, "[INFO] " message "\n", ##__VA_ARGS__)
#define LOG_ERROR(message, ...) fprintf(stderr, "[ERROR] " message "\n", ##__VA_ARGS__)
#define LOG_WARNING(message, ...) fprintf(stderr, "[WARNING] " message "\n", ##__VA_ARGS__)

using ll = long long;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 = int16_t;
using i32 = int16_t;
using i64 = int16_t;
using byte = unsigned char;

// --------------------- constants ----------------------
constexpr int INF = 1 << 30;
constexpr long long LINF = 1LL << 60;
constexpr int MAX_CARD_TYPE = 5;
constexpr double DEFAULT_TL = 1.90;

// --------------------- global variables ----------------------
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
std::uniform_int_distribution<> rand01(0, 1);
std::uniform_real_distribution<> randReal(0, 1);
clock_t start_time;

inline double get_time(clock_t startTime)
{
  return (double)(clock() - startTime) / CLOCKS_PER_SEC;
}

template <typename T>
const T& clamp(const T& value, const T& low, const T& high) {
    return std::min(std::max(value, low), high);
}


inline double gauss(double mu = 0.0, double sigma = 1.0)
{
  std::normal_distribution<double> dist(mu, sigma);
  return dist(mt);
}