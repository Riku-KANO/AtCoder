#pragma once
/**
 * @file io.hpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */
#include <global.hpp>
#include <OilField.hpp>

struct Input
{
    int N;      // grid size
    int M;      // the number of oil fields
    double eps; // error parameter

    Vec<OilField> oil_fields;

    Vec<int> imin;
    Vec<int> jmin;
    Vec<Vec<int>> v;
    Vec<double> epsilons;

    static Input read_input()
    {
        Input ret;
        std::cin >> ret.N >> ret.M >> ret.eps;

        ret.oil_fields.resize(ret.M);

        for (OilField &oil : ret.oil_fields)
        {
            int d;
            std::cin >> d;
            oil.positions.resize(d);
            for (Point &point : oil.positions)
            {
                std::cin >> point.i >> point.j;
                oil.height = std::max(oil.height, point.i + 1);
                oil.width = std::max(oil.width, point.j + 1);
            }
        }

#ifdef LOCAL
        // input for local judge
        imin.resize(this->M);
        jmin.resize(this->M);
        v.resize(this->N, Vec<int>(this->N));
        epsilons.resize(this->N * this->N * 2);

        for (int i = 0; i < this->M; i++)
        {
            std::cin >> imin[i] >> jmin[i];
        }

        for (int i = 0; i < this->N; i++)
        {
            for (int j = 0; j < this->N; j++)
            {
                std::cin >> v[i][j];
            }
        }

        for (int i = 0; i < this->N * this->N * 2; i++)
        {
            std::cin >> epsilons[i];
        }

#endif

        return ret;
    }
};

enum class QueryType
{
    QUERY,
    ANSWER
};

struct Output
{
    QueryType type;
    int n_pos;
    Vec<Point> positions;
    Output(): n_pos(0){};

    void show()
    {
        std::cout << "qa"[type == QueryType::ANSWER] << " ";
        std::cout << n_pos;
        for (const Point &position : this->positions)
        {
            std::cout << " " << position.i << " " << position.j;
        }

        std::cout << "\n";
        std::cout.flush();
    }
};

struct Communication
{
    Output query;
    Response response;
};

class Judge
{
public:
    Judge() {}

    void init(const Input &input)
    {
        N = input.N;
        M = input.M;
        eps = input.eps;
        imin = input.imin;
        jmin = input.jmin;
        v = input.v;

        for(double e: input.epsilons)
        {
            epsilons.push(e);
        }

        cost = 0.0;
        score = 1e9;
    }

    int respond(const Output &output)
    {
        assert(output.n_pos == (int)output.positions.size());
        int res;

#ifdef LOCAL
        switch (output.type)
        {
        case QueryType::ANSWER:
            res = judge_answer(output);
            break;

        case QueryType::QUERY:
            res = judge_query(output);
            break;

        default:
            break;
        }
#else
        std::cin >> res;
#endif

        if (output.type == QueryType::QUERY)
        {
            this->add_cost(output.n_pos);
        }
        return res;
    }

    int calc_score(const Output& output)
    {
        int res = this->judge_answer(output);
        if (res == 0) {
            return 1e9;
        }
        return std::round(1e6 * cost);
    }

private:
    int N;
    int M;
    double eps;
    Vec<int> imin;
    Vec<int> jmin;
    Vec<Vec<int>> v;
    std::queue<double> epsilons;

    double cost;
    double score;

    int judge_answer(const Output& output)
    {
        int res = 1;

        Vec<Vec<int>> expect(this->N, Vec<int>(this->N, 0));
        for(const Point& point: output.positions)
        {
            auto[i, j] = point;
            expect[i][j] = 1;
        }

        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                res &= (expect[i][j] ^ (v[i][j]>0));
            }
        }

        return res;
    }


    int judge_query(const Output& output)
    {
        int res;
        int k = output.n_pos;
        int vs = 0;
        double e = epsilons.front(); epsilons.pop();

        for(const Point& point: output.positions) {
            auto[i, j] = point;
            vs += this->v[i][j];
        }
        double mu = (k - vs) * eps + vs * (1-eps);
        double diff = k * eps * (1 - eps) * e;
        res = std::max(0, (int)std::round(mu + diff));
        return res;
    }

    void add_cost(int n_pos) {
        this->cost += 1.0 / std::sqrt(n_pos);
    }
};