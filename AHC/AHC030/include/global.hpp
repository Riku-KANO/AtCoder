#pragma once
/**
 * @file global.hpp
 * @author Rick (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-12
 * 
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <istream>
#include <string>
#include <array>
#include <vector>
#include <set>
#include <queue>
#include <algorithm>
#include <thread>
#include <memory>
#include <utility>
#include <chrono>
#include <random>
#include <numeric>
#include <type_traits>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>

// boost multiprecision
#include <boost/multiprecision>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>
namespace mp = boost::multiprecision;
using FloatType = mp::number<mp::cpp_dec_float<50>>;
#include <macro.hpp>
#include <constant.hpp>