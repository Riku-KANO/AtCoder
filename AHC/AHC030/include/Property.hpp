#pragma once
/**
 * @file Property.hpp
 * @author Rick (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-02-12
 * 
 */
#include <global.hpp>


struct Property
{
    bool verbose;
    u64 seed;
    double time_limit;
    int num_thread;

    Property();

    static Property read_property(int argc, char *argv[]);
    
    static constexpr double default_time_limit = 2.9;
    static constexpr bool default_verbose = false;
    static constexpr int default_num_thread = 2;
};