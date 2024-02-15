#pragma once
/**
 * @file Solver.hpp
 * @author Rick (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-12
 * 
 */
#include <global.hpp>
#include <io.hpp>
#include <Property.hpp>
#include <util.hpp>
#include <rng.hpp>
#include <evolution.hpp>

class Solver
{
public:
    Solver();
    Solver(const Input& _input, const Property& _props) = noexcept;

    void init() = noexcept;
    Output query(int turn);
    void update(const Output& output, Response res);
private:
    Input input;
    Property props;
    Rng rng;
    std::array<int, MAX_NUM_CELL> exp_v;
    Vec<OilField> oils;

    std::queue<Output> initial_query_queue;

    Environment env;
    Vec<Specie> species;
    std::set<Specie> extinct;
    Specie ans_specie;

    /**
     * @brief 
     * 
     * @return Vec<Specie> 
     */
    Vec<Specie> select();

    /**
     * 必要性があるか？
     * @brief 
     * 
     * @param s 
     * @param t 
     * @return Specie 
     */
    Specie crossover(const Specie& s, const Specie& t);

    /**
     * @brief 
     * 
     * @param s 
     * @return Specie 
     */
    Specie mutate(const Specie& s);

    /**
     * @brief 
     * 
     * @param specie 
     * @return double 
     */
    double calc_fitness(const Specie& specie);
};