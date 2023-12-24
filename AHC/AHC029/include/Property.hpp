#pragma once
#include <string>
#include <cstring>
#include "global.hpp"

struct Property
{
    double TL = DEFAULT_TL;
    unsigned int seed;
    Property():seed(seed_gen()){}
};

Property arg_parse(int argc, char *argv[])
{
  Property param;
  for (int i = 0; i < argc; i++)
  {
    if (strcmp(argv[i], "--seed") == 0)
    {
      if (i + 1 > argc)
      {
        LOG_ERROR("No arguments.");
      }
      int _seed = std::stoi(argv[i + 1]);
      param.seed = _seed;
    }
    if (strcmp(argv[i], "--TL") == 0)
    {
      if (i + 1 > argc)
      {
        LOG_ERROR("No arguments.");
      }
      double _TL = std::stod(argv[i + 1]);
      param.TL = _TL;
    }
  }
  return param;
}