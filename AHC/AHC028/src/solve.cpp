#pragma GCC target("avx2")
#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <tuple>
#include <bitset>
#include <queue>
#include <deque>
#include <complex>
#include <memory>
#include <algorithm>
#include <numeric>
#include <utility>
#include <functional>
#include <random>
#include <iomanip>
#include <limits>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
using namespace std;
using ll = long long;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;

#define PII pair<int, int>
#define FOR(i, s, e) for (int i = s; i < e; i++)
#define RFOR(i, s, e) for (int i = s, i >= e; i--)
#define Vec vector

template <typename T>
void chmin(T &x, T y) { x = min(x, y); }
template <typename T>
void chmax(T &x, T y) { x = max(x, y); }

random_device seed_gen;
mt19937 mt(seed_gen());
clock_t start_time;
uniform_real_distribution<> rand01(0, 1);

constexpr int INF = 1 << 30;
constexpr long long LINF = 1LL << 60;
constexpr double TL = 1.98;
const int DI[4] = {0, 1, 0, -1};
const int DJ[4] = {1, 0, -1, 0};
const double PI = std::acos(-1);
constexpr int MAX_N = 15;
constexpr int MAX_M = 200;
constexpr int TURN_LIMIT = 5000;

double get_time(clock_t start_time)
{
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

//-- globals --
int N;
int M;
int si, sj;
string A[MAX_N];
string T[MAX_M];

//-------------

struct Word
{
    string prefixs[5];
    string suffixs[5];

    Word() {}
    Word(string s)
    {
        FOR(i, 0, 5)
        {
            string pre = s.substr(0, i + 1);
            prefixs[i] = pre;
            string suf = s.substr(4 - i);
            suffixs[i] = suf;
        }
    }
};

struct State
{
    vector<int> p;
    vector<int> path;
    int score;
    int eval;

    int calc_score()
    {
        int ci = si;
        int cj = sj;
        int total_cost = 0;
        for (int idx : path)
        {
            int ni = idx / N;
            int nj = idx % N;
            int di = abs(ni - ci);
            int dj = abs(nj - cj);
            total_cost += (di + dj + 1);
            ci = ni;
            cj = nj;
        }
        score = 10000 - total_cost;
        return total_cost;
    }

    int eval_score()
    {
        int ci = si;
        int cj = sj;
        int total_cost = 0;
        for (int idx : path)
        {
            int ni = idx / N;
            int nj = idx % N;
            int di = abs(ni - ci);
            int dj = abs(nj - cj);
            total_cost += (di + dj + 1) * (di+dj+1);
            ci = ni;
            cj = nj;
        }
        eval = total_cost;
        return total_cost;
    }

};

class Solver
{
public:
    vector<vector<vector<int>>> G;
    vector<vector<int>> letters; // [26] -> (positions)
    vector<Word> words;
    vector<vector<int>> match; // i-th suffix to j-th prefix mathing degree
    Solver() {}
    void init()
    {
        words.resize(M);
        FOR(i, 0, M)
        {
            words[i] = Word(T[i]);
        }

        cerr << "word done\n";

        match.resize(M, vector<int>(M));
        FOR(i, 0, M)
        {
            FOR(j, 0, M)
            {
                if (i == j)
                    continue;
                int how_match = 0;
                Word w1 = words[i];
                Word w2 = words[j];
                FOR(m, 0, 5)
                {
                    if (w1.suffixs[m] != w2.prefixs[m])
                    {
                        break;
                    }
                    how_match++;
                }
                match[i][j] = how_match;
            }
        }

        letters.resize(26);
        FOR(i, 0, N)
        FOR(j, 0, N)
        {
            int idx = i * N + j;
            int c = A[i][j] - 'A';
            letters[c].emplace_back(idx);
        }
        cerr << "letter done\n";

        G.resize(N * N, vector<vector<int>>(26));
        FOR(i, 0, N)
        {
            FOR(j, 0, N)
            {
                int from = i * N + j;
                FOR(ii, 0, N)
                {
                    FOR(jj, 0, N)
                    {
                        int to = ii * N + jj;
                        int c = A[ii][jj] - 'A';
                        G[from][c].emplace_back(to);
                    }
                }

                FOR(c, 0, 26)
                {
                    sort(G[from][c].begin(), G[from][c].end(), [&](int i, int j)
                         {
                        int i1 = i / N;
                        int j1 = i % N;
                        int i2 = j / N;
                        int j2 = j % N;
                        int i0 = from / N;
                        int j0 = from % N;
                        int di10 = i1-i0;
                        int dj10 = j1-j0;
                        int di20 = i2-i0;
                        int dj20 = j2-j0;
                        return di10*di10 + dj10*dj10 < di20*di20 + dj20*dj20; });
                }
            }
        }
        cerr << "graph done\n";
    }

    void solve()
    {
        _solve();
    }

private:
    void _solve()
    {
        const int MAX_STATE = 50;
        vector<State> states(MAX_STATE);
        vector<int> indices(M);
        iota(indices.begin(), indices.end(), 0);
        cerr << "state init done\n";

        State cur_state;
        cur_state.p = indices;
        shuffle(cur_state.p.begin(), cur_state.p.end(), mt);
        // cur_state.p = generate_initial();
        int ci = si;
        int cj = sj;
        for (int i : cur_state.p)
        {
            for (char ch : T[i])
            {
                int cidx = ci * N + cj;
                int c = ch - 'A';
                int nidx = G[cidx][c].front();
                int ni = nidx / N;
                int nj = nidx % N;
                ci = ni;
                cj = nj;
                cur_state.path.emplace_back(nidx);
            }
        }

        cur_state.calc_score();
        cur_state.eval_score();
        State best_state = cur_state;
        constexpr double start_temp = 100000000;
        constexpr double end_temp = 10000;
        int iter = 0;
        while (1)
        {
            iter++;
            double cur_time = get_time(start_time);
            if(cur_time > TL) {
                break;
            }

            int pre_score = cur_state.score;
            int pre_eval = cur_state.eval;

            int i1, i2;
            i1 = (u32)mt() % M;
            i2 = i1;
            
            while (i1 == i2 || abs(i1 - i2) > 100)
            {
                i2 = (u32)mt() % M;
            }
            swap(cur_state.p[i1], cur_state.p[i2]);
            cur_state.path.clear();
            string next_path = concatenate(cur_state.p);
            int ci = si;
            int cj = sj;
            for(int i = 0; i < next_path.size(); i++)
            {
                char ch = next_path[i];
                int cidx = ci * N + cj;
                int c = ch - 'A';
                
                int nidx = G[cidx][c].front();
                int ni = nidx / N;
                int nj = nidx % N;
                ci = ni;
                cj = nj;
                cur_state.path.emplace_back(nidx);
            }

            cur_state.calc_score();
            cur_state.eval_score();
            if (cur_state.eval < best_state.eval)
            {
                best_state = cur_state;
#ifdef LOCAL
                cerr << "best updated!\n";
#endif
            } else {
                double cur_temp = ((TL-cur_time) * start_temp + cur_time * end_temp) / TL;
                double diff = cur_state.eval - pre_eval;
                double thre = rand01(mt);

                if(exp(-diff / cur_temp) > thre) {
                    swap(cur_state.p[i1], cur_state.p[i2]);
                }
            }

           
        }

        for (int idx : best_state.path)
        {
            int ci = idx / N;
            int cj = idx % N;
            cout << ci << " " << cj << "\n";
        }
        cout.flush();
        cerr << "best score: " << best_state.score << "\n";
        cerr << "num iter: " << iter << "\n";
    }

  

    string concatenate(const vector<int>& perm) {
        string ret = "";

        for(int i = 0; i < M; i++) {
            int idx = perm[i];
            if(i == 0) {
                ret += T[idx];
            } else {
                int size = ret.size();
                string suf1 = ret.substr(size - 1);
                string suf2 = ret.substr(size - 2);
                string suf3 = ret.substr(size - 3);
                string suf4 = ret.substr(size - 4);
                string suf5 = ret.substr(size - 5);
                string pre1 = T[idx].substr(0, 1);
                string pre2 = T[idx].substr(0, 2);
                string pre3 = T[idx].substr(0, 3);
                string pre4 = T[idx].substr(0, 4);
                string pre5 = T[idx].substr(0, 5);
                if(suf5 == pre5) {
                    continue;
                } else if(suf4 == pre4) {
                    ret += T[idx].substr(4);
                } else if(suf3 == pre3) {
                    ret += T[idx].substr(3);
                } else if(suf2 == pre2) {
                    ret += T[idx].substr(2);
                } else if(suf1 == pre1) {
                    ret += T[idx].substr(1);
                } else {
                    ret += T[idx];
                }
            }
        }

        return ret;
    }

    void greedy()
    {
        vector<PII> ans;
        int ci = si;
        int cj = sj;

        FOR(i, 0, M)
        {
            if (ans.size() == TURN_LIMIT)
                break;
            for (char ch : T[i])
            {
                int idx = ci * N + cj;
                int c = ch - 'A';
                int nidx = G[idx][c].front();
                int ni = nidx / N;
                int nj = nidx % N;
                ans.emplace_back(ni, nj);
                ci = ni;
                cj = nj;
            }
        }

        for (PII a : ans)
        {
            cout << a.first << " " << a.second << "\n";
        }
        cout.flush();
    }
};

int main()
{
    start_time = clock();
    cin >> N >> M;
    cin >> si >> sj;
    cerr << N << " " << M << " " << si << " " << sj << "\n";
    FOR(i, 0, N)
        cin >> A[i];
    FOR(i, 0, M)
        cin >> T[i];
    cerr << "Input done\n";
    Solver solver;
    solver.init();
    solver.solve();
    cerr << "Elapsed time: " << get_time(start_time) << "\n";
    return 0;
}