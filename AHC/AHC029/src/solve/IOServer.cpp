#include "IOServer.hpp"

inline void IOServer::init(int n, int m, int k, int t)
{
    this->cur_turn = 0;
    this->n_card = n;
    this->n_project = m;
    this->n_given = k;
    this->n_turn = t;
    this->n_invest = 0;
    this->money = 0;

    this->read_input();

}

inline void IOServer::read_action(const Action &action, Vec<Card> *hands)
{

}

inline void IOServer::read_next_card(int next_card_index)
{
}

inline Vec<Project> IOServer::serve_initial_projects() const
{
    return this->initial_projects;
}

inline Vec<Card> IOServer::serve_initial_cards() const
{
    return this->initial_cards;
}

inline Vec<Project> IOServer::serve_projects() const
{
    return this->cur_projects;
}

inline MONEY IOServer::serve_money() const
{
    return this->money;
}

inline Vec<std::pair<Card, MONEY>> IOServer::serve_cards() const
{
    return this->cards[this->cur_turn];
}

inline void IOServer::read_input()
{

    initial_projects.resize(n_project);
    projects.resize(n_turn * n_project);
    initial_cards.resize(n_card);
    cards.resize(n_turn, Vec<std::pair<Card, MONEY>>(n_given));

    for (Project &project : initial_projects)
    {
        std::cin >> project.remain >> project.value;
    }

    for (Project &project : projects)
    {
        std::cin >> project.remain >> project.value;
    }

    for (Card &card : initial_cards)
    {
        std::cin >> card.type >> card.w;
    }

    for (Vec<std::pair<Card, MONEY>> &turn_cards : cards)
    {
        for (auto &[card, cost] : turn_cards)
        {
            std::cin >> card.type >> card.w >> cost;
        }
    }
}

inline void IOServer::increment_turn() {
    this->cur_turn++;
}