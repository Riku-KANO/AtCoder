#pragma once
/**
 * @file App.hpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */

#include <Property.hpp>
#include <io.hpp>
#include <Solver.hpp>
#include <global.hpp>

#ifdef RICK
#include <logger.hpp>
#endif 

class App
{
public:
    App(const Input& _input, const Property& _props);
    void init() noexcept;
    void run();

    void summary();

private:
    const Input input;
    const Property props;
    Solver solver;
    Judge judge;
    clock_t start_time;
#ifdef RICK
    const Logger logger;
#endif
};