#pragma once
#include "game_object.hpp"
#include "global.hpp"

class IOServer
{
public:
    inline void init(int n, int m, int k, int t);
    inline void read_action(const Action &action);
    inline void read_next_card(int next_card_index);
    inline Vec<Project> serve_initial_projects() const;
    inline Vec<Card> serve_initial_cards() const;
    inline Vec<Project> serve_projects() const;
    inline MONEY serve_money() const;
    inline Vec<std::pair<Card, MONEY>> serve_cards() const;
    inline void increment_turn();
private:
    int cur_turn;
    int n_card;
    int n_project;
    int n_given;
    int n_turn;
    int n_invest;
    MONEY money;

    Vec<Project> initial_projects;
    Vec<Project> projects;

    Vec<Card> initial_cards;
    Vec<Vec<std::pair<Card, MONEY>>> cards; // [turn][k] -> Card,MONEY pair

    inline void read_input();
};