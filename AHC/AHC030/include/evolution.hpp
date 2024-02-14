#pragma once
/**
 * @file evolution.hpp
 * @author Rick (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-12
 * 
 */
#include <global.hpp>
#include <util.hpp>
#include <io.hpp>

struct Environment
{
    int N;
    std::array<int, MAX_NUM_CELL> fixed_v;
    Vec<Communication> conditions;

    Environment()
    {
        fixed_v.fill(-1);
    }
    Environment(int n): N(n)
    {
        fixed_v.fill(-1);
    }

    /**
     * @brief update environment from query and response.
     * When the size of query position is 1, fill the value.
     * 
     * @param output 
     * @param response 
     */
    void update(const Output& output, Response response)
    {
        conditions.push_back(Communication{output, response});

        if(output.n_pos == 1)
        {
            auto[i, j] = output.positions.front();
            fixed_v[get_index(i, j, N)] = response;
        }
    }

    double calc_fitness(const Specie& specie)
    {
        double fitness = 0.0;
        
    }
};

/**
 * @brief 
 * @todo fitnessをmultiprecisionにするかどうか
 */
struct Specie
{
    Hash hash;
    Vec<Point> points;
    double fitness;


    inline bool operator==(const Specie& other)
    {
        return this->hash == other.hash;
    }

    inline bool operator!=(const Specie& other)
    {
        return this->hash != other.hash;
    }

    inline bool operator>(const Specie& other)
    {
        return this->fitness > other.fitness;
    }

    inline bool operator>=(const Specie& other)
    {
        return this->fitness >= other.fitness;
    }

    inline bool operator<(const Specie& other)
    {
        return this->fitness < other.fitness;
    }

    inline bool operator<=(const Specie& other)
    {
        return this->fitness <= other.fitness;
    }


    static constexpr u64 MOD = (u64)1e9 + 7;
    static constexpr u64 base[40] = {
       820608610,  92754042,  22394908, 670808574, 312706920, 100352652,
       264053809, 456513636, 574762218,  67573902, 680304321, 336976067,
       749554964, 142004569, 388132101, 633100979, 252163558, 436545860,
       150060349, 129563184, 299160421, 748783071,  69924266, 688828624,
       283103665,  70889374,  72092303, 676334153, 772871052, 879973420,
       526226928, 462891622, 293849627, 941041968,  91063185, 285279698,
       281664676, 371649357, 635775890, 840404143
    };

    static Hash calc_hash(const Vec<Point>& _points);

    static Hash calc_hash(const Vec<Point>& _points) {
        Hash ret = 0;
        int m = _points.size();

        for(int i = 0; i < m; i++)
        {
            ret += _points[i].i * base[2 * i] + _points[i].j * base[2 * i + 1];
        }

        return ret % MOD;
    }

};