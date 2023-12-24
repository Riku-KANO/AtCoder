#pragma once

#include "global.hpp"
#include "game_object.hpp"

struct Input 
{
    int n_card, n_project, n_given, n_turn;

    static Input get_input(){
        Input ret;

        std::cin >> ret.n_card >> ret.n_project >> ret.n_given >> ret.n_turn;
        
        return ret;
    }
};