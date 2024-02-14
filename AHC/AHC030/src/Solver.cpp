/**
 * @file Solver.cpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */
#include <Solver.hpp>

Solver::Solver() {}

Solver::Solver(const Input &_input, const Property &_props)
{
    this->input = _input;
    this->oils = _input.oil_fields;
    this->props = props;
    this->rng = Rng(_props.seed);
    this->exp_v.fill(0);
}

void Solver::init()
{
    this->exp_v.fill(0);
    this->species.clear();
    this->extinct.clear();

    const int MAX_SPECIE = 1000;
    for (int i = 0; i < MAX_SPECIE; i++)
    {
        Specie specie;
        for (int j = 0; j < input.M; j++)
        {
            int imin = rng.generage_randint(0, oils[j].height - 1);
            int jmin = rng.generage_randint(0, oils[j].width - 1);
            specie.points.push_back(Point{imin, jmin});
        }
        specie.hash = Specie::calc_hash(specie.points);
        species.push_back(specie);
    }
}

inline void Solver::update(const Output &output, Response res)
{
    const int N = input.N;
    const int M = input.M;

    env.update(output, res);

    if (output.type == QueryType::ANSWER && res == 0)
    {
        extinct.insert(ans_specie);
    }

    if (output.type == QueryType::QUERY && output.n_pos == 1)
    {
        auto [i, j] = output.positions.front();
        exp_v[get_index(i, j, N)] = res;
    }
}

inline Output Solver::query(int turn)
{
    const int N = input.N;
    const int M = input.M;
    const int NUM_CELL = N * N;
    const double eps = input.eps;
    Output ret;

    // 初手はテキトーに探索箇所をいくつか生成してqueueに入れておく。。
    if (turn == 0)
    {
        const int NUM_INITIAL_DATA = 2 * N * N;
        const int NUM_PICK = NUM_CELL / 2;
        Vec<int> position_indices(NUM_CELL);
        std::iota(position_indices.begin(), position_indices.end(), 0);
        for (int i = 0; i < NUM_INITIAL_DATA; i++)
        {
            Output initial_output;
            Vec<Point> &positions = initial_output.positions;
            initial_output.type == QueryType::QUERY;
            initial_output.n_pos = NUM_PICK;
            positions.resize(NUM_PICK);
            std::shuffle(position_indices.begin(), position_indices.end(), rng.mt);
            std::sort(position_indices.begin(), position_indices.begin() + NUM_PICK);
            for (int j = 0; j < NUM_PICK; j++)
            {
                int r = position_indices[j] / N;
                int c = position_indices[j] % N;
                positions[j] = Point{r, c};
            }

            initial_query_queue.push(initial_output);
        }
    }

    // キューにoutputが入ってたらそれを出す。
    if (!initial_query_queue.empty())
    {
        ret = initial_query_queue.front();
        initial_query_queue.pop();

        return ret;
    }

    const int NUM_REMAIN = 10;
    const int NUM_NEXT_GENERATION = 100;
    for (Specie &specie : species)
    {
        specie.fitness = env.calc_fitness(specie);
    }

    std::sort(species.begin(), species.end(), [&species](const Specie &lhs, const Specie &rhs)
              { return lhs.fitness > rhs.specie; });

    while(species.size() > NUM_REMAIN)
    {
        species.pop_back();
    }

    Vec<Specie> next_generation(NUM_REMAIN); // 後で増えていく。とりあえずNUM_REMAINだけ確保
    Vec<Specie> remain_species = selection();



    Vec<double> cum_fitness(NUM_REMAIN, 0.0);
    cum_fitness[0] = species[0].fitness;
    for (int i = 1; i < NUM_REMAIN; i++)
    {
        cum_fitness[i] += cum_fitness[i - 1] + species[i].fitness;
    }

    for (; (int)next_generation.size() < NUM_NEXT_GENERATION;)
    {
        double rand_value = rng.generate_randreal(0, cum_fitness.back());
        int target_specie_idx = std::lower_bound(cum_fitness.begin(), cum_fitness.end(), rand_value) - cum_fitness.begin();
        Specie next_specie = remain_species[i].mutate();
        next_specie.hash = Specie::calc_hash(next_specie.points);
        species.push_back(next_specie);
    }

    // remove duplicates
    std::sort(species.begin(), species.end(), [](const Specie& lhs, const Specie& rhs){
        return lhs.hash < rhs.hash;
    });
    species.erase(std::unique(species.begin(), species.end()), species.end());



    std::array<int, MAX_NUM_CELL> count;
    Vec<Vec<int>> data(NUM_CELL, Vec<int>(NUM_NEXT_GENERATION));
    for (const Specie &specie : next_generation)
    {
        count.fill(0);
        for (int oid = 0; oid < M; oid++)
        {
            auto [imin, jmin] = specie.points[oid];
            for (auto [i, j] : oils[oid].positions)
            {
                count[get_index(imin + i, jmin + j, N)] += 1;
            }
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                int idx = get_index(i, j, N);
                data[idx].push_back(count[idx]);
            }
        }
    }

    const int NUM_TARGET_POSITION = NUM_CELL / 2;
    // specieによる予測値の分散が大きいところを調べる。
    // TODO: 計算量チェック
    Vec<double> sigma2s(NUM_CELL);
    for (int idx = 0; idx < NUM_CELL; idx++)
    {
        double sigma2 = math::var(data[idx]);
        sigma2s.emplace_back(sigma2);
    }

    // cell id のsort
    Vec<int> cell_indices(NUM_CELL);
    std::iota(cell_indices.begin(), cell_indices.end(), 0);
    std::sort(
        cell_indices.begin(),
        cell_indices.end(),
        [&sigma2s](int i, int j)
        {
            return sigma2s[i] > sigma2s[j];
        });

    // assign
    ret.positions.clear();
    ret.positions.resize(NUM_TARGET_POSITION);
    for (int idx = 0; idx < NUM_TARGET_POSITION; idx++)
    {
        ret.positions[idx] = Point{
            cell_indices[idx] / N,
            cell_indices[idx] % N};
    }

    return ret;
}

inline Vec<Specie> Solver::selection(int num_remain)
{
}

inline Specie Solver::crossover(const Specie &s, const Specie &t)
{
}

/**
 * @brief
 *
 * @param s
 * @return Specie
 */
inline Specie Solver::mutate(const Specie &s)
{
    const int N = input.N;
    const int M = input.M;
    const double eps = input.eps;
    Specie new_specie = s;

    // swap positions
    int rand_value;
    Vec<int> oil_indices(M);
    std::iota(oil_indices.begin(), oil_indices.end(), 0);
    std::shuffle(oil_indices.begin(), oil_indices.end(), rng.mt);
    for (int from = 0; from < M; from++)
    {
        int to = oil_indices[from];
        new_specie.points[from] = Point
        {
            std::min(s.points[to].i, oils[from].height - 1),
                std::min(s.points[to].j, oils[from].width - 1)
        }
    }

    // position move
    for (int i = 0; i < M; i++)
    {
        rand_value = rng.generage_randint(0, 2);
        if (rand_value == 1)
        {
            Vec<std::pair<int, int>> dxys;
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    if (dx == 0 && dy == 0)
                    {
                        continue;
                    }
                    int nx = new_specie.points[i].i + dy;
                    int ny = new_specie.points[i].j + dx;
                    if (is_inside(ny, nx, N))
                    {
                        dxys.push_back({dx, dy});
                    }
                }
            }

            auto [dx, dy] = rng.random_pick(dxys);

            new_specie.points[i] = Point{
                new_specie.points[i].i + dy,
                new_specie.points[i].j + dx};
        }
    }

    new_specie.hash = Specie::calc_hash(new_specie.points);

    return new_specie;
}