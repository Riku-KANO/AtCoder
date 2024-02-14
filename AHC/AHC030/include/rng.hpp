#pragma once
/**
 * @file rng.hpp
 * @author Rick (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-13
 * 
 */
#include <global.hpp>

struct Rng
{
    u64 seed;
    std::mt19937 mt;

    Rng()
    {
        std::random_device seed_gen();
        seed = (u64)seed_gen();
        mt = std::mt(seed);
    }

    Rng(u64 _seed): seed(_seed), mt(_seed){}

    /**
     * @brief random value of [0, m)
     * 
     * @param m 
     * @return template<typename T> 
     */
    template<typename T> inline T generate_randint(T m)
    {
        return (T)mt() % m;
    }

    /**
     * @brief random value of closed form [l, r]
     * 
     * @param l 
     * @param r 
     * @return template<typename T> 
     */
    template<typename T> inline T generage_randint(T l, T r)
    {
        assert(r >= l);
        T d = r - l + 1;
        return l + (T)mt() % d;
    }

    /**
     * @brief generate real number random value inside [l, r).
     * 
     * @tparam T 
     * @param l 
     * @param r 
     * @return T 
     */
    template<typename T> inline T generate_randreal(T l, T r)
    {
        assert(r >= l);
        constexpr u64 MAX_D = (u64)1e18;
        T  d = r - l;
        return l + (T)(mt() % MAX_D) / MAX_D * d;
    }

    /**
     * @brief pick one random element from std::vector<T>.
     * 
     * @tparam T 
     * @param v 
     * @return T 
     */
    template<typename T> inline T random_pick(const Vec<T>& v)
    {
        const int NUM_ELEMENT = v.size();
        assert(NUM_ELEMENT > 0);
        int idx = (u32)mt() % NUM_ELEMENT;
        return v[idx];
    }
}