/**
 * @file Property.cpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */

#include <Property.hpp>

Property::Property()
{
    verbose = default_verbose;
    time_limit = default_time_limit;

    std::random_device seed_gen;
    seed = (u64)seed_gen();
}

Property Property::read_property(int argc, char *argv[])
{
    Property ret;

    for (int i = 0; i < argc; i++)
    {
        /**
         * verbose
         * 
         */
        if (std::strcmp(argv[i], "--verbose") == 0)
        {
            if (i + 1 > argc)
            {
                std::cerr << "No arguments\n";
            }

            bool verbose = std::stoull(argv[i + 1]);
            ret.verbose = verbose;
        }

        /**
         * seed 
         * 
         */
        if (std::strcmp(argv[i], "--seed") == 0)
        {
            if (i + 1 > argc)
            {
                std::cerr << "No arguments\n";
            }

            u64 seed = std::stoull(argv[i + 1]);
            ret.seed = seed;
        }

        /**
         * time limit 
         * 
         */
        if (std::strcmp(argv[i], "--TL") == 0)
        {
            if (i + 1 > argc)
            {
                std::cerr << "No arguments\n";
            }
            double _TL = std::stod(argv[i + 1]);
            ret.time_limit = _TL;
        }
    }

    return ret;
}
