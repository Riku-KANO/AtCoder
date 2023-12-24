#pragma once

#include <queue>
#include "input.hpp"
#include "Property.hpp"
#include "Player.hpp"

class Game
{
public:
    Game();
    void init(const Input& input, const Property& props);
    void run();

private:
    const int MAX_INVEST = 20;
    int n_card;    // 2 <= N <= 7
    int n_given;   // 2 <= K <= 5
    int n_project; // 2 <= M <= 8
    int n_turn;    // T = 1000
    int n_invest;
    MONEY money;
    Vec<Card> hands;
    Vec<Project> projects;
    Vec<bool> disappeared;
    player::Player player;
    
    inline void read_initial_projects();
    inline void read_initial_cards();
    inline void process_action(const Action& action);
    inline void read_projects();
    inline void read_money();
    inline Vec<std::pair<Card, MONEY>> read_cards(int turn);
    inline void process_choosed_card(int next_card_index, const Action& action, int turn);
    inline void update(const Action& action, int next_card_index);

#ifdef LOCAL
    std::queue<Project> queued_projects;

    Vec<Vec<std::pair<Card, MONEY>>> next_cards; // [turn][k] -> Card,MONEY pair

    inline void local_init();
    inline void local_input();
    inline void update_projects_and_money(const Action& action);
    inline void update_hands(int next_card_index, const Action& action, int turn);
#endif
};