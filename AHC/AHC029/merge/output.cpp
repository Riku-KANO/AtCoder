#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <utility>

#define pii std::pair<int, int>
#define Vec std::vector
#define FOR(i, n) for(int i = 0; i < (int)(n); i++)
#define MONEY long long
#define LOG_INFO(message, ...) fprintf(stderr, "[INFO] " message "\n", ##__VA_ARGS__)
#define LOG_ERROR(message, ...) fprintf(stderr, "[ERROR] " message "\n", ##__VA_ARGS__)
#define LOG_WARNING(message, ...) fprintf(stderr, "[WARNING] " message "\n", ##__VA_ARGS__)

using ll = long long;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 = int16_t;
using i32 = int16_t;
using i64 = int16_t;
using byte = unsigned char;

// --------------------- constants ----------------------
constexpr int INF = 1 << 30;
constexpr long long LINF = 1LL << 60;
constexpr int MAX_CARD_TYPE = 5;
constexpr double DEFAULT_TL = 1.90;

// --------------------- global variables ----------------------
std::random_device seed_gen;
std::mt19937 mt(seed_gen());
std::uniform_int_distribution<> rand01(0, 1);
std::uniform_real_distribution<> randReal(0, 1);
clock_t start_time;

inline double get_time(clock_t startTime)
{
  return (double)(clock() - startTime) / CLOCKS_PER_SEC;
}

template <typename T>
const T& clamp(const T& value, const T& low, const T& high) {
    return std::min(std::max(value, low), high);
}


inline double gauss(double mu = 0.0, double sigma = 1.0)
{
  std::normal_distribution<double> dist(mu, sigma);
  return dist(mt);
}


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


struct Input 
{
    int n_card, n_project, n_given, n_turn;

    static Input get_input(){
        Input ret;

        std::cin >> ret.n_card >> ret.n_project >> ret.n_given >> ret.n_turn;
        
        return ret;
    }
};



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


class App
{
public:
    App(const Property& _props);
    void init();
    void run();
    void summary();

private:
    const Property props;
    Game game;
};

// ******************************* main ************************************
int main(int argc, char *argv[])
{
    Property props;
    start_time = clock();

#ifdef LOCAL
    props = arg_parse(argc, argv);
#endif

    App app(props);

    app.init();
    app.run();
#ifdef LOCAL
    app.summary();
#endif
    return 0;
}
// *************************************************************************


App::App(const Property& _props): props(_props){}

void App::init() {
    
    Input input = Input::get_input();
    
    game.init(input, props);

    LOG_INFO("Application Initialization done at %6.4f[s]", get_time(start_time));
}

void App::run() {
    
    game.run();

    LOG_INFO("Application Run Successfuly Ended at %6.4f[s]", get_time(start_time));
}

void App::summary() {

}


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
        for (Action action : next_actions)
        {
            double score = random_playout(action, 50);
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