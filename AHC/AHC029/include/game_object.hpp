#pragma once

#include "global.hpp"

struct Card
{
    int type;
    long long work;
    Card() {}
    Card(int t, long long w) : type(t), work(w) {}
};

struct Project
{
    long long remain;
    long long initial_remain;
    MONEY value;
    Project() {}
    Project(long long h, long long v) : remain(h), initial_remain(h), value(v) {}
};

struct Action
{
    int card_index;
    int target_project;
    Action()
        : card_index(0),
          target_project(0) {}
    Action(int _card_index, int _target_project)
        : card_index(_card_index),
          target_project(_target_project) {}
};