/**
 * @file App.cpp
 * @author Rick (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-02-12
 *
 */
#include <App.hpp>
App::App(const Input& _input, const Property& _props) : input(_input), props(_props)
{
    solver = Solver(input, props);
    fprintf(
        stderr,
        "N: %2d, M: %2d, eps: %.2f\n",
        this->input.N, this->input.M, this->input.eps);
}

void App::init()
{
    solver.init();
    judge.init(this->input);
}

void App::run()
{
    const int MAX_TURN = 2 * input.N * input.N;
    Output output;
    int res;
    for(int turn = 0; turn < MAX_TURN; turn ++)
    {
        output = solver.query(turn);
        output.show();

        res = judge.respond(output);

        if(output.type == QueryType::ANSWER && res == 1)
        {
            break;
        }

        solver.update(output, res);
    }

    int final_score = judge.calc_score(output);

    std::cerr << "elapsed: " << get_time(this->start_time) << "\n";
    std::cerr << "score: " << final_score << "\n";
}

