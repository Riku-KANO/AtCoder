#include "Player.hpp"

#include <iostream>

namespace player
{
    inline void Player::init(
        const Vec<Card> &_cards,
        const Vec<Project> &_projects,
        const Input &input,
        const Property &props)
    {
        hands = _cards;
        projects = _projects;

        n_card = input.n_card;
        n_project = input.n_project;
        n_given = input.n_given;
        n_turn = input.n_turn;
        n_invest = 0;

        money = 0;
        rng = std::mt19937(seed_gen());
        this->card_type_cnts.resize(MAX_CARD_TYPE, 0);

        // count the number of types
        for (Card card : _cards)
        {
            card_type_cnts[card.type]++;
        }
    }

    inline void Player::read_status(
        const Vec<Card> &_cards,
        const Vec<Project> &_projects,
        int _n_invest,
        int _turn)
    {
        hands = _cards;
        projects = _projects;
        n_invest = _n_invest;
        cur_turn = _turn;
    }

    inline void Player::read_money(MONEY _money)
    {
        this->money = _money;
    }

    inline void Player::read_next_cards(
        const Vec<std::pair<Card, MONEY>> &next_cards)
    {
        int n_given = next_cards.size();
        FOR(i, n_given)
        {
            if (i == 0)
            {
                continue;
            }
            this->card_type_cnts[next_cards[i].first.type]++;
        }
    }

    inline Action Player::action()
    {
        update_hiddens();
        Vec<Action> next_actions = get_next_actions(hands, projects, money);
        Vec<double> scores;
        int limit = 50;
#ifdef DEBUG
        limit = 5;
#endif
        for (Action action : next_actions)
        {
            double score = random_playout(action, limit);
            scores.emplace_back(score);
        }

        // arg max
        double highest = -1e19;
        int target_action_idx = -1;
        FOR(idx, next_actions.size())
        {
            if (scores[idx] > highest)
            {
                highest = scores[idx];
                target_action_idx = idx;
            }
        }

#ifdef DEBUG
        if (target_action_idx == -1)
        {
            LOG_ERROR("Player cannot find next action index at turn %3d, at line %d",
                      cur_turn,
                      __LINE__);
        }
#endif

        return next_actions[target_action_idx];
    }

    inline int Player::choose_card(
        const Vec<std::pair<Card, MONEY>> &next_cards)
    {
        // todo
        return 0;
    }

    inline Vec<int> Player::get_card_type_cnts() const
    {
        return this->card_type_cnts;
    }

    inline double Player::random_playout(
        const Action &first_action,
        int limit = 50)
    {
        GameState gamestate;
        gamestate.init(
            hands,
            projects,
            money,
            cur_turn,
            n_given,
            n_turn,
            n_invest,
            get_card_type_cnts());

        Action next_action = first_action;
        while (!gamestate.is_over(limit))
        {
            gamestate.request_action(next_action);
            gamestate.prepare_next_cards();
            int card_index = choose_next_card_at_random(gamestate);
            gamestate.request_card(card_index);

            gamestate.increment_turn();

            next_action = random_action(gamestate);
        }

        return gamestate.evaluate();
    }

    inline void Player::update_hiddens()
    {
        // todo
    }

    inline Action Player::rulebase_action()
    {
    }

    inline int Player::rulebase_card_choice(const Vec<std::pair<Card, MONEY>> &next_cards)
    {
    }

    inline Action Player::random_action()
    {
        Vec<Action> next_actions = get_next_actions(hands, projects, money);
        int idx = (u32)rng() % next_actions.size();
        return next_actions[idx];
    }

    inline Action Player::random_action(const GameState &gamestate)
    {
        Vec<Action> next_actions = get_next_actions(
            gamestate.hands,
            gamestate.projects,
            gamestate.money);
        int idx = (u32)rng() % next_actions.size();
        return next_actions[idx];
    }

    inline void GameState::init(
        const Vec<Card> &_hands,
        const Vec<Project> &_projects,
        MONEY _money,
        int _turn,
        int _n_given,
        int _n_turn,
        int _n_invest,
        const Vec<int> &card_type_cnts)
    {
        start_turn = _turn;
        cur_turn = _turn;
        n_given = _n_given;
        hands = _hands;
        projects = _projects;
        money = _money;

        n_card = hands.size();
        n_project = projects.size();
        n_turn = _n_turn;
        n_invest = _n_invest;

        generator.init(card_type_cnts, n_invest, n_project);
    }

    inline int Player::choose_next_card_at_random(
        const Vec<std::pair<Card, MONEY>> &next_cards)
    {
        Vec<int> available_indices = filter_available_cards(next_cards, money, n_invest);
        int idx = (u32)mt() % available_indices.size();
        return idx;
    }

    inline int Player::choose_next_card_at_random(
        const GameState &gamestate)
    {
        Vec<int> available_indices = filter_available_cards(gamestate.next_cards, gamestate.money, gamestate.n_invest);
        int idx = (u32)mt() % available_indices.size();
        return idx;
    }

    inline bool GameState::is_over(int limit)
    {
        return cur_turn - start_turn == limit || cur_turn == n_turn;
    }

    inline void GameState::request_action(const Action &action)
    {
        int type = hands[action.card_index].type;
        long long work = hands[action.card_index].work;
        int target_project = action.target_project;

        auto process_work = [&](int target_project, long long amount)
        {
            projects[target_project].remain -= amount;
            if (projects[target_project].remain <= 0)
            {
                this->money += projects[target_project].value;
                Project next_project = generator.generate_next_project();
                projects[target_project] = next_project;
            }
        };

        auto process_cancel = [&](int target_project)
        {
            Project next_project = generator.generate_next_project();
            projects[target_project] = next_project;
        };

        if (type == 0)
        {
            process_work(target_project, work);
        }
        else if (type == 1)
        {
            FOR(idx, n_project)
            {
                process_work(target_project, work);
            }
        }
        else if (type == 2)
        {
            process_cancel(target_project);
        }
        else if (type == 3)
        {
            FOR(idx, n_project)
            {
                process_cancel(target_project);
            }
        }
        else if (type == 4)
        {
            this->n_invest++;
            generator.update(this->n_invest);
        }

        this->pre_action = action;
    }

    inline void GameState::prepare_next_cards()
    {
        this->next_cards = this->generator.generate_cards(n_given);
    }

    inline void GameState::request_card(int card_index)
    {
        auto [next_card, cost] = next_cards[card_index];
        hands[pre_action.card_index] = next_card;
        this->money -= cost;
    }

    inline void GameState::increment_turn()
    {
        this->cur_turn++;
    }

    inline double GameState::evaluate()
    {
        double score = 0.0;

        score += (double)this->money * cur_turn / n_turn;
        for (const Project &project : projects)
        {
            score += (double)project.value / project.remain * (1LL << n_invest);
        }

        for (const Card &card : hands)
        {
            if(card.type == 0)
            {
                score += (double)card.work;
            }
            else if(card.type == 1)
            {
                score += (double)card.work * n_project;
            }
            else if(card.type == 4)
            {
                score += (double)10 * (1 << n_invest);
            }
        }
        return score;
    }

    inline void Generator::init(
        const Vec<int> &_type_cnts,
        int _n_invest,
        int _n_project)
    {
        this->n_invest = _n_invest;
        this->n_project = _n_project;
        this->type_cnts = _type_cnts;
        this->rng = std::mt19937(seed_gen());

        accum.resize(MAX_CARD_TYPE);

        // create accumulated array
        bool enough_type = true;
        int total_cnt = 0;
        FOR(type, MAX_CARD_TYPE)
        {
            total_cnt += type_cnts[type];
            if (type_cnts[type] == 0)
            {
                enough_type = false;
            }
        }

        Vec<double> type_ratio(MAX_CARD_TYPE);

        if (enough_type)
        {
            type_ratio[0] = (double)type_cnts[0] / total_cnt;
            type_ratio[1] = (double)type_cnts[1] / total_cnt;
            type_ratio[2] = (double)type_cnts[2] / total_cnt;
            type_ratio[3] = (double)type_cnts[3] / total_cnt;
            type_ratio[4] = (double)type_cnts[4] / total_cnt;
        }
        else
        {
            int x0 = (u32)rng() % 20 + 1;
            int x1 = (u32)rng() % 10 + 1;
            int x2 = (u32)rng() % 10 + 1;
            int x3 = (u32)rng() % 5 + 1;
            int x4 = (u32)rng() % 3 + 1;
            int total_x = x0 + x1 + x2 + x3 + x4;
            type_ratio[0] = (double)x0 / total_x;
            type_ratio[1] = (double)x1 / total_x;
            type_ratio[2] = (double)x2 / total_x;
            type_ratio[3] = (double)x3 / total_x;
            type_ratio[4] = (double)x4 / total_x;
        }

        this->accum[0] = type_ratio[0];
        this->accum[1] = this->accum[0] + type_ratio[1];
        this->accum[2] = this->accum[1] + type_ratio[2];
        this->accum[3] = this->accum[2] + type_ratio[3];
        this->accum[4] = this->accum[3] + type_ratio[4]; // must be 1.0
    }

    inline void Generator::update(int _n_invest)
    {
        this->n_invest = _n_invest;
    }

    inline Project Generator::generate_next_project()
    {
        return this->_generate_next_project();
    }

    inline Vec<Project> Generator::generate_projects(int n_gen)
    {
        Vec<Project> projects;
        FOR(_, n_gen)
        {
            projects.push_back(this->_generate_next_project());
        }
        return projects;
    }

    inline Vec<std::pair<Card, MONEY>> Generator::generate_cards(int n_gen)
    {
        Vec<std::pair<Card, MONEY>> next_cards;

        next_cards.push_back(std::pair<Card, MONEY>{Card(0, 1LL << n_invest), 0LL});

        FOR(_, n_gen - 1)
        {
            next_cards.push_back(this->_generate_next_card());
        }

        return next_cards;
    }

    inline Project Generator::_generate_next_project()
    {
        const double MIN_B = 2.0;
        const double MAX_B = 8.0;
        std::uniform_real_distribution<> dist_b(MIN_B, MAX_B);
        double b = dist_b(rng);
        double h = std::round(std::pow(2.0, b)) * (1 << n_invest);
        double v = std::round(std::pow(2.0, clamp(gauss(b, 0.5), 0.0, 10.0))) * (1 << n_invest);
        return Project(h, v);
    }

    inline std::pair<Card, MONEY> Generator::_generate_next_card()
    {
        int next_type;
        long long next_work = 0;
        MONEY next_cost;

        double random_value = randReal(rng);

        next_type = std::lower_bound(
                        this->accum.begin(),
                        this->accum.end(),
                        random_value) -
                    this->accum.begin();

#ifdef DEBUG
        if (next_type < 0 || next_type >= MAX_CARD_TYPE)
        {
            LOG_ERROR("Generated card type must be within 0 and 4 at line %d", __LINE__);
        }
#endif
        if (next_type == 0 || next_type == 1)
        {
            next_work = ((u64)rng() % 50 + 1) * 1LL << n_invest;
        }

        // generate cost
        if (next_type == 0)
        {
            double w_tmp = (double)next_work / (1 << n_invest);
            next_cost = (long long)clamp((int)gauss(w_tmp, w_tmp / 3), 1, 10'000) * (1 << n_invest);
        }
        else if (next_type == 1)
        {
            double w_tmp = (double)next_work / (1 << n_invest);
            next_cost = (long long)clamp((int)gauss(w_tmp * n_project, w_tmp * n_project / 3), 1, 10'000) * (1 << n_invest);
        }
        else if (next_type == 2)
        {
            next_cost = ((u32)rng() % 11) * (1 << n_invest);
        }
        else if (next_type == 3)
        {
            next_cost = ((u32)rng() % 11) * (1 << n_invest);
        }
        else if (next_type == 4)
        {
            next_cost = ((u32)rng() % 801 + 200) * (1 << n_invest);
        }

        return {Card(next_type, next_work), next_cost};
    }

    inline int choose_next_card_at_random(
        const Vec<std::pair<Card, MONEY>> &next_cards,
        MONEY money,
        int n_invest)
    {
        Vec<int> available_indices = filter_available_cards(next_cards, money, n_invest);
        int idx = (u32)mt() % available_indices.size();
        return idx;
    }

    inline Vec<int> filter_available_cards(const Vec<std::pair<Card, MONEY>> &next_cards, MONEY money, int n_invest)
    {
        Vec<int> available_indices;
        FOR(idx, next_cards.size())
        {
            auto [next_card, cost] = next_cards[idx];
            if (next_card.type == 4 && n_invest == 20)
            {
                continue;
            }
            if (cost <= money)
            {
                available_indices.push_back(idx);
            }
        }
        return available_indices;
    }

    inline Vec<Action> get_next_actions(const Vec<Card> &hands, const Vec<Project> &projects, MONEY money)
    {
        Vec<Action> next_actions;
        int n_card = hands.size();
        int n_project = projects.size();
        FOR(card_index, n_card)
        {
            Card card = hands[card_index];
            Action next_action;
            next_action.card_index = card_index;
            if (card.type == 0 || card.type == 2)
            {
                FOR(target_project, n_project)
                {
                    next_action.target_project = target_project;
                    next_actions.push_back(next_action);
                }
            }
            else
            {
                next_actions.push_back(next_action);
            }
        }
        return next_actions;
    }
}