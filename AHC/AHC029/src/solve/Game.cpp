#include "Game.hpp"

Game::Game() {}

void Game::init(const Input &input, const Property &props)
{
    n_card = input.n_card;
    n_project = input.n_project;
    n_given = input.n_given;
    n_turn = input.n_turn;
    n_invest = 0;

#ifdef LOCAL

    local_init();

#else

    read_initial_projects();
    read_initial_cards();
    this->disappeared.resize(n_project, false);
#endif

    player.init(this->hands, this->projects, input, props);

#ifdef DEBUG
    LOG_INFO("Game Initialization done at %6.4f[s]", get_time(start_time));
#endif
}

void Game::run()
{
    const int N_TURN = n_turn;
    for (int turn = 0; turn < N_TURN; turn++)
    {

        // action phase
        player.read_status(hands, projects, n_invest, turn);
        Action action = player.action();

        this->process_action(action);

        this->read_projects();
        this->read_money();
        player.read_money(money);
        Vec<std::pair<Card, MONEY>> given_cards = this->read_cards(turn);

        // card choice phase
        int next_card_index = player.choose_card(given_cards);

        this->process_choosed_card(next_card_index, action, turn);

// for visualizer comment
#ifdef VISUALIZER

#endif

#ifdef DEBUG
        LOG_INFO("TURN %03d money: %10d, at %6.4f[s]", get_time(start_time));
#endif

        // turn end update
        this->update(action, next_card_index);
    } // game loop

#ifdef DEBUG
    LOG_INFO("Game Ended at %6.4f[s]", get_time(start_time));
#endif
}

inline void Game::process_action(const Action &action)
{
#ifdef LOCAL
    this->update_projects_and_money(action);
#else
    // output
    std::cout << action.card_index << " " << action.target_project << std::endl;

    // check disappear or not
    int type = hands[action.card_index].type;
    if (type == 0)
    {
        int work = hands[action.card_index].work;
        if (work > projects[action.target_project].remain)
        {
            disappeared[action.target_project] = true;
        }
    }
    else if (type == 1)
    {
        int work = hands[action.card_index].work;
        FOR(i, n_project)
        {
            if (work > projects[i].remain)
            {
                disappeared[i] = true;
            }
        }
    }
    else if (type == 2)
    {
        disappeared[action.target_project] = true;
    }
    else if (type == 3)
    {
        FOR(i, n_project)
        {
            disappeared[i] = true;
        }
    }
#endif
}

inline void Game::read_initial_projects()
{
    this->projects.clear();
    this->projects.resize(n_project);
    for (Project &project : this->projects)
    {
        std::cin >> project.remain >> project.value;
        project.initial_remain = project.remain;
    }
}

inline void Game::read_initial_cards()
{
    this->hands.clear();
    this->hands.resize(n_card);
    for (Card &card : this->hands)
    {
        std::cin >> card.type >> card.work;
    }
}

inline void Game::read_projects()
{
#ifdef LOCAL
#else
    int idx = 0;
    for (Project &project : this->projects)
    {
        std::cin >> project.remain >> project.value;
        if (disappeared[idx])
        {
            project.initial_remain = project.remain;
            disappeared[idx] = false;
        }
        idx++;
    }
#endif
}

inline void Game::read_money()
{
#ifdef LOCAL
#else
    std::cin >> this->money;
#endif
}

inline Vec<std::pair<Card, MONEY>> Game::read_cards(int turn)
{
    Vec<std::pair<Card, MONEY>> cards;
#ifdef LOCAL
    cards = next_cards[turn];
#else
    cards.resize(n_given);
    for (auto &[card, cost] : cards)
    {
        std::cin >> card.type >> card.work >> cost;
    }
#endif
    return cards;
}

inline void Game::process_choosed_card(int next_card_index, const Action &action, int turn)
{
#ifdef LOCAL
    this->update_hands(next_card_index, action, turn);
#else
    std::cout << next_card_index << std::endl;
#endif
}

inline void Game::update(const Action &action, int next_card_index)
{
}

#ifdef LOCAL
inline void Game::local_init()
{

    local_input();
}

inline void Game::local_input()
{

    projects.resize(n_project);
    hands.resize(n_card);
    next_cards.resize(n_turn, Vec<std::pair<Card, MONEY>>(n_given));

    for (Project &project : projects)
    {
        std::cin >> project.remain >> project.value;
    }
    FOR(i, n_turn * n_project)
    {
        Project project;
        std::cin >> project.remain >> project.value;
        queued_projects.push(project);
    }

    for (Card &card : hands)
    {
        std::cin >> card.type >> card.work;
    }

    for (Vec<std::pair<Card, MONEY>> &turn_cards : next_cards)
    {
        for (auto &[card, cost] : turn_cards)
        {
            std::cin >> card.type >> card.work >> cost;
        }
    }
}

inline void Game::update_projects_and_money(const Action &action)
{
    Card card = hands[action.card_index];
    int type = card.type;
    long long work = card.work;
    int target_project = action.target_project;

    auto process_work = [&](int target_project, int amount)
    {
        projects[target_project].remain -= amount;
        if (projects[target_project].remain <= 0)
        {
            this->money += projects[target_project].value;
            Project next_project = queued_projects.front();
            queued_projects.pop();
            next_project.remain *= 1 << n_invest;
            next_project.value *= 1 << n_invest;
            projects[target_project] = next_project;
            projects[target_project].initial_remain = next_project.remain;
        }
    };

    auto process_cancel = [&](int target_project)
    {
        Project next_project = queued_projects.front();
        queued_projects.pop();
        next_project.remain *= 1 << n_invest;
        next_project.value *= 1 << n_invest;
        projects[target_project] = next_project;
        projects[target_project].initial_remain = next_project.remain;
    };

    if (type == 0)
    {
        process_work(target_project, work);
    }
    else if (type == 1)
    {
        FOR(i, n_project)
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
        FOR(i, n_project)
        {
            process_cancel(i);
        }
    }
    else if (type == 4)
    {
        this->n_invest++;
    }
    else
    {
        LOG_ERROR("INVALID action type");
        std::exit(1);
    }
}

inline void Game::update_hands(int next_card_index, const Action &action, int turn)
{
    auto [next_card, cost] = next_cards[turn][next_card_index];
    hands[action.card_index] = next_card;
    money -= cost;
}
#endif