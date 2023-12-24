#pragma once

#include "game_object.hpp"
#include "Property.hpp"
#include "input.hpp"

namespace player
{
    class Generator
    {
    public:
        /**
         * @brief
         *
         * @param _type_cnts
         * @param _n_invest
         * @param _n_project
         */
        inline void init(
            const Vec<int> &_type_cnts,
            int _n_invest,
            int _n_project);

        inline void update(int _n_invest);

        /**
         * @brief
         *
         * @return Project
         */
        inline Project generate_next_project();

        /**
         * @brief
         *
         * @param n_project
         * @return Vec<Project>
         */
        inline Vec<Project> generate_projects(int n_gen);

        /**
         * @brief
         *
         * @param n_given
         * @return Vec<std::pair<Card, MONEY>>
         */
        inline Vec<std::pair<Card, MONEY>> generate_cards(int n_gen);

    private:
        int n_invest;
        int n_project;
        Vec<int> type_cnts;
        Vec<double> accum;
        std::mt19937 rng;

        /**
         * @brief
         *
         * @return Project
         */
        inline Project _generate_next_project();

        /**
         * @brief
         *
         * @return std::pair<Card, MONEY>
         */
        inline std::pair<Card, MONEY> _generate_next_card();
    };


    class GameState
    {
    public:
        Vec<Card> hands;
        Vec<Project> projects;
        Vec<std::pair<Card, MONEY>> next_cards;
        Action pre_action;
        MONEY money;
        int n_invest;
        /**
         * @brief
         *
         * @param _hands
         * @param _projects
         * @param _turn
         */
        inline void init(
            const Vec<Card> &_hands,
            const Vec<Project> &_projects,
            MONEY _money,
            int _turn,
            int _n_given,
            int _n_turn,
            int _n_invest,
            const Vec<int> &card_type_cnts);

        /**
         * @brief
         *
         * @param limit
         * @return true
         * @return false
         */
        inline bool is_over(int limit);

        /**
         * @brief 
         * 
         * @param action 
         */
        inline void request_action(const Action &action);

        /**
         * @brief 
         * 
         */
        inline void prepare_next_cards();

        /**
         * @brief 
         * 
         * @param card_index 
         */
        inline void request_card(int card_index);

        inline void increment_turn();

        /**
         * @brief
         *
         * @return double
         */
        inline double evaluate();

    private:
        int start_turn;
        int cur_turn;
        int n_card;
        int n_project;
        int n_given;
        int n_turn;
        Generator generator;
    };

    class Player
    {
    public:
        /**
         * @brief
         *
         * @param _cards
         * @param _projects
         * @param props
         */
        inline void init(
            const Vec<Card> &_cards,
            const Vec<Project> &_projects,
            const Input &input,
            const Property &props);

        /**
         * @brief
         *
         * @param _cards
         * @param _projects
         * @param _n_invest
         * @param _turn
         */
        inline void read_status(
            const Vec<Card> &_cards,
            const Vec<Project> &_projects,
            int _n_invest,
            int _turn);

        /**
         * @brief
         *
         * @param _money
         */
        inline void read_money(
            MONEY _money);

        /**
         * @brief
         *
         * @param _next_cards
         */
        inline void read_next_cards(
            const Vec<std::pair<Card, MONEY>> &_next_cards);

        /**
         * @brief
         *
         * @return Action
         */
        inline Action action();

        /**
         * @brief
         *
         * @param next_cards
         * @return int
         */
        inline int choose_card(
            const Vec<std::pair<Card, MONEY>> &next_cards);

        /**
         * @brief Get the card type cnts object
         *
         * @return Vec<int>
         */
        inline Vec<int> get_card_type_cnts() const;

    private:
        int n_card;
        int n_project;
        int n_given;
        int n_turn;
        int cur_turn;
        int n_invest;
        Vec<Card> hands;
        Vec<Project> projects;
        Vec<int> card_type_cnts; //
        MONEY money;
        std::mt19937 rng;

        inline double random_playout(
            const Action &first_action,
            int limit);

        /**
         * @brief
         *
         */
        inline void update_hiddens();

        /**
         * @brief
         *
         * @return Action
         */
        inline Action rulebase_action();

        /**
         * @brief
         *
         * @param next_cards
         * @return int
         */
        inline int rulebase_card_choice(
            const Vec<std::pair<Card, MONEY>> &next_cards);

        /**
         * @brief Referencing its hands
         *
         * @return Action
         */
        inline Action random_action();

        /**
         * @brief Referemcing game state's hands
         *
         * @param gamestate
         * @return Action
         */
        inline Action random_action(const GameState &gamestate);

        /**
         * @brief
         *
         * @param next_cards
         * @return int
         */
        inline int choose_next_card_at_random(
            const Vec<std::pair<Card, MONEY>> &next_cards);

        /**
         * @brief
         *
         * @param gamestate
         * @return int
         */
        inline int choose_next_card_at_random(const GameState &gamestate);
    };

    /**
     * @brief
     *
     * @param next_cards
     * @param money
     * @return int
     */
    inline int choose_next_card_at_random(
        const Vec<std::pair<Card, MONEY>> &next_cards,
        MONEY money);

    /**
     * @brief 買うことができるカードのインデックスを返す
     *
     * @param next_cards
     * @param money
     * @return Vec<int>
     */
    inline Vec<int> filter_available_cards(
        const Vec<std::pair<Card, MONEY>> &next_cards,
        MONEY money,
        int n_invest);

    /**
     * @brief Get the next actions object
     *
     * @param hands
     * @param projects
     * @param money
     * @return Vec<Action>
     */
    inline Vec<Action> get_next_actions(
        const Vec<Card> &hands,
        const Vec<Project> &projects,
        MONEY money);

}