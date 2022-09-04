#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include <set>
#include <tuple>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <iterator>
#include <functional>
#include <utility>

#include <boost/math/distributions.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#define DEBUG
// #define TEST

#define rep(i, n) for (int i = 0; i < (n); i++)
#define ll long long
#define pii pair<int, int>
#define EncodedPat std::vector<std::pair<char, char>>
#define EncodedChange std::vector<std::pair<std::pair<int, int>, char>>
#define dmatrix boost::numeric::ublas::matrix<double>
#define dvector boost::numeric::ublas::vector<double>
#define Patterns tuple<std::vector<std::string>, EncodedPat, EncodedChange>

using namespace std;
namespace bmath = boost::math;
namespace ublas = boost::numeric::ublas;

struct Operation
{
    std::vector<ll> costA = std::vector<ll>(10);
    std::vector<ll> costB = std::vector<ll>(10);
};

const int dx[4] = {1, 0, -1, 0};
const int dy[4] = {0, 1, 0, -1};
const int dx8[8] = {-1, -1, -1, 0, 1, 1, 1, 0};
const int dy8[8] = {-1, 0, 1, -1, -1, 0, 1, 1};
const int INF = 1 << 30;
const ll LINF = 1LL << 60;
const double epsilon = 1e-6;

// global variable//////
double TL;
double TIME_LIMIT_90;
double DECAY; // スコアが返ってきて０だった時、想定スコアにこのDECAYをかける。そこから名目的なスコアを作成。
double SEARCH_TIME;
bool USE_LOG;
int X, M, C, E;
std::vector<Operation> productLine(21);

clock_t start_time;
clock_t cur_time;
double best_eval_score = -1e20;
std::vector<std::string> best_eval_pattern;
std::set<std::vector<std::string>> pattern_set;
vector<double> score_table(301, 0.0);
vector<vector<bool>> fixed_op;
int num_use_change = 0;

const int HOUR = 86'400;
const vector<int> operation_time = {0, 3, 5, 8, 10, 12, 14, 16, 18};
vector<pair<int, ll>> heavy_cost_order; // 機械に応じてコストの差がありすぎるのでここで管理する。
////////////////////////

class Config
{
public:
    double tl = 5.0;
    double decay = 0.90;
    bool use_log = false;
    std::string file_name = "../output/config.txt";
    Config() {}
    void init()
    {
        TL = this->tl;
        TIME_LIMIT_90 = this->tl * 0.95;
        DECAY = this->decay;
        USE_LOG = this->use_log;
        return;
    }

    void save()
    {
        ofstream ofs;
        ofs.open(this->file_name, ios::out);
        cout << "N: " << X << "\n";
        cout << "M: " << M << "\n";
        cout << "C: " << C << "\n";
        cout << "E: " << E << "\n";
        cout << "TL: " << this->tl << "\n";
        cout << "DECAY: " << this->decay << "\n";
        cout << "USE LOG: " << this->use_log << endl;
        ofs.close();
        return;
    }
};

random_device seed_gen;
mt19937 mt(seed_gen());
bmath::normal_distribution<> distribution(0, 1);

bool out_bound(int ni, int nj, int h, int w)
{
    if (ni < 0 || nj < 0 || ni >= h || nj >= w)
        return true;
    return false;
}

void debug()
{
    std::cerr << "### DEBUG ###\n";
    return;
}

void test()
{
    std::cerr << "### TEST ###\n";
    return;
}

void input(int &X, int &M, int &C, int &E, std::vector<Operation> &productLine)
{
    cin >> X >> M >> C >> E;
    productLine.resize(M);
    for (int i = 0; i < M; i++)
    {
        for (int p = 1; p <= 9; p++)
        {
            cin >> productLine[i].costA[p] >> productLine[i].costB[p];
        }
    }
    return;
}

void init(std::vector<ll> &scores, std::vector<double> &eval_scores, std::vector<ll> &exp_scores)
{
    std::cerr << "\n####### INIT START ######" << endl;
    Config config;
    config.init();

    // initial score for "11111...111" size = 2*X * M
    rep(i, M) exp_scores[0] += productLine[i].costA[1] + productLine[i].costB[1];

    scores[0] = 0;
    eval_scores[0] = 0; //評価値も0にしておく。
    best_eval_pattern = std::vector<std::string>(M, std::string(2 * X, '1'));
    SEARCH_TIME = TIME_LIMIT_90 / (double)E;
    vector<vector<bool>> fixed_tmp(M, vector<bool>(2 * X, false));
    fixed_op = fixed_tmp;
    //ルーレットのスコアテーブル
    score_table.resize(E + 1);

    // 重たいコストを管理して９番目の稼働のコスト値でソートする。
    rep(i, M) {
        heavy_cost_order.push_back({i, productLine[i].costA[9]});
    }
    sort(heavy_cost_order.begin(), heavy_cost_order.end(), [](pair<int,ll>a, pair<int,ll>b){return a.second > b.second;});

#ifdef TEST
    std::cerr << "Parameter saving..." << std::endl;
    config.save();
#endif
    std::cerr << "\n#####################" << std::endl;
    std::cerr << "X: " << X << ", M: " << M << ", C: " << C << ", E: " << E << "\n";
    std::cerr << "TL: " << TL << "\n";
    std::cerr << "DECAY: " << DECAY << "\n";
    std::cerr << "SEARCH TIME:" << SEARCH_TIME << "\n";
    std::cerr << "#####################\n\n";
    return;
}

// メモ
// 全部でX週, 設備はM個, 注文N件
// 注文は納期d, 工程の数J, 工程を処理する設備rと所要時間t
// dは秒単位
// 計画プログラムのVは稼働変更の上限を上回った回数、Dは納期遅れの注文の数

// 稼働パターン1: 稼働させない
// 稼働パターン2: 9:00-12:00
// 稼働パターン3: 13:00-18:00
// 稼働パターン4: 9:00-12:00 + 13:00-18:00
// 稼働パターン5: 9:00-12:00 + 13:00-20:00
// 稼働パターン6: 9:00-12:00 + 13:00-22:00
// 稼働パターン7: 9:00-12:00 + 13:00-24:00
// 稼働パターン8: 9:00-12:00 + 13:00-26:00
// 稼働パターン9: 9:00-12:00 + 13:00-28:00
// １週間は平日5日と休日2日からなり平日と休日でオペレーションコストが異なる。
// 稼働パターンはMこの文字列で表される。P_1, .., P_M
// それぞれのパターンは長さ2Xの文字列(平日と休日で2つ × X週)
// 稼働パターンが変化しすぎ作業員が手配できにくくなるので上限がある。設備ごとにC回まで。
// ある週の平日と翌週の平日で稼働パターンが異なっていること、およびある週の休日と翌週の休日で稼働パターンが異なっていることを指します。
// 調整プログラムは、以上の条件の下で適切な稼働パターンを探すために、この後で説明する計画プログラムと以下の流れでやりとりを行います。
// 1. 調整プログラムが稼働パターンを出力する
// 2. 計画プログラムが調整プログラムの出力と注文データを使って計画を立てる
// 3. 計画プログラムが計画の評価指標を計算する
// 4. 調整プログラムに評価結果が返される
// このインタラクティブ操作をE回繰り返す。

// 計画プログラムは以下の二つの条件が満たされるように計画を立てます。
// j>1 のとき、作業 (i,j) が開始できるのは作業 (i,j−1) が開始した時刻以降である
// j>1 のとき、作業 (i,j) が終了できるのは作業 (i,j−1) が終了した時刻以降である

// 最終テストでは
// E = 50, 100, 300 ぞれぞれ1000個のテストケースで判断。

// 制約
// 10 <= M <= 20
// 8 <= X <= 16
// 2 <= C <= 8
// Eは50, 100, 300のいずれか

// アイデア
// ・ブラックボックス最適化・ベイズ最適化？　+ 貪欲 <=　最有力
// ・焼きなまし??
// ・遺伝的アルゴリズム???
//

//　その他
// C++の E=300の時のIOの所要時間は0.1~0.2s程度。ちなみにPythonのIOは4.5秒程度。遅すぎるのでnumpyの行列演算ライブラリは意味をなさない。
// 時間的な制約の関係でループ時間が短すぎると次の候補が見つからない。初期点サンプリングを多めにやってもいいかもしれない。<- 嘘かもしれない。時間のせいではない
// もしくは次の候補を貪欲的に選ぶのもありかもしれない。
// ベイズ最適化を実行してみたがE=300で1ループ１０秒以上かかる時がある。計算量がきつい。
// E = 300の時は実行してみてわかったがルーレットをやっているだけで3Gを超える時がある。ベイズ最適化してる暇は最初にしかない。
// 0000.txtにおいて、ルーレット+ベイズ最適化 で2.65G, ルーレットだけは2.4G。多少なりの効果はあり。
// 0000.txtで対数を使わないベイズの方が2.6Gを上回る傾向が強い。

// LU分解で求める。
ublas::matrix<double> calc_invmat(ublas::matrix<double> A)
{
    ublas::matrix<double> B = ublas::identity_matrix<double>(A.size1());
    ublas::permutation_matrix<> pm(A.size1());
    ublas::lu_factorize(A, pm);
    ublas::lu_substitute(A, pm, B);
    return B;
}

// ２つのデータ間の距離関数（オリジナルカーネル）
// カーネルの性質(正定性・対称性・非退化性)は全て満たすようになっている。
// RBFカーネルを利用したりできそう
double my_kernel(const std::vector<string> &a, const std::vector<string> &b)
{
    double ret = 0.;
    const double d = 10.0 * M * X;
    for (int m = 0; m < M; m++)
    {
        for (int x = 0; x < 2 * X; x++)
        {
            int aa = a[m][x] - '0';
            int bb = b[m][x] - '0';
            if (x % 2 == 0)
                ret += 5. * (double)(aa - bb) * (aa - bb); // TODO: 係数の調整
            else
                ret += 2. * (double)(aa - bb) * (aa - bb); // TODO: 係数の調整
        }
    }
    ret = exp(-ret / d);
    return ret;
}

// next pointsが1個の場合の関数
pair<double, double> gaussian_process(const std::vector<std::vector<std::string>> &pat, const dvector &eval_scores, const std::vector<std::string> &next_point, int epoch)
{
    pair<double, double> ret;
    int ker_size = epoch;
    dmatrix x(ker_size, ker_size);
    dmatrix x2 = ublas::identity_matrix<double>(ker_size); // x2は正則行列を足したもの。
    for (int i = 0; i < ker_size; i++)
    {
        for (int j = 0; j < ker_size; j++)
        {
            x(i, j) = my_kernel(pat[i], pat[j]);
            x2(i, j) += x(i, j);
        }
        // cerr << i << " OK\n";
    }
    dvector aux(ker_size), aux2(ker_size), aux3(ker_size);
    rep(i, ker_size) aux(i) = my_kernel(pat[i], next_point);

    dmatrix ker_inv = calc_invmat(x);
    dmatrix ker_inv2 = calc_invmat(x2);
    aux2 = ublas::prod(ker_inv, eval_scores);
    // rep(i, ker_size) rep(j, ker_size) aux2(i) += log10(eval_scores[j] + 1) * ker_inv(i, j);
    double mu = ublas::inner_prod(aux, aux2);
    // rep(i, ker_size)
    // {
    //     mu += aux[i] * aux2[i];
    //     cerr << "aux: " << aux[i] << " " << aux2[i] << "\n";
    // }
    aux3 = ublas::prod(ker_inv2, aux);
    // rep(i, ker_size) rep(j, ker_size) aux3[i] += ker_inv2(i, j) * aux[j]; // TODO: debug不足  要チェック
    double var = my_kernel(next_point, next_point) - ublas::inner_prod(aux, aux3);
    // rep(i, ker_size) var -= aux3[i] * aux[i];
    ret = {mu, var};
    // cerr << mu << " " << var << endl;
    return ret;
}

std::vector<double> EI(const std::vector<double> &mu, const std::vector<double> &var, double best)
{
    std::vector<double> lamb;
    rep(i, (int)mu.size()) lamb.push_back((mu[i] - best) / (var[i] * 1.0));
    std::vector<double> z(mu.size());
    rep(i, (int)mu.size()) z[i] = (mu[i] - best) * bmath::cdf(distribution, lamb[i]) + var[i] * bmath::pdf(distribution, lamb[i]);
    return z;
}

double EI(double mu, double var, double best)
{
    double lamb = (mu - best) / (var * 1.0);
    double z = (mu - best - 0.01) * bmath::cdf(distribution, lamb) + var * bmath::pdf(distribution, lamb);
    return z;
}

// O(N)の探索
int argmax(std::vector<double> &x)
{
    int ret = 0;
    double M = -1e10;
    rep(i, (int)x.size())
    {
        if (x[i] > M)
        {
            M = x[i];
            ret = i;
        }
    }
    return ret;
}


// sampling v3 コストが最も高いものを早めに打ち切るようなサンプリング
Patterns first_sampling_ver3(int ep, const std::vector<std::vector<std::string>> &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_changes, const std::vector<ll> &scores, auto &RD)
{
    std::vector<std::string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_change_pat;
    if (ep == 1)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }
    }
    if(ep == 2) {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }

        bool flag=true;
        int h_id=0;
        while(flag) {
            int m_id = heavy_cost_order[h_id].first;
            int change_week; 
            for(int x = X-1; x >= 0; x--) {
                if(RD[1][m_id][x] > epsilon) {
                    change_week=x+1;
                    break;
                }
            }
            ret_change_pat.push_back({{m_id, change_week * 2}, '1'});
            num_use_change++;
            if(num_use_change < C) {
                ret_change_pat.push_back({{m_id, change_week * 2 + 1}, '1'});
                num_use_change++;
            }
            h_id++;
            if(h_id == M or num_use_change == C) break;
        }
        // decode
        for(auto p: ret_change_pat) {
            int m = p.first.first;
            for(int x = p.first.second; x < 2 * X; x+=2) {
                ret_pat[m][x] = p.second;
            }
        }
    }
    return {ret_pat, ret_enc_pat, ret_change_pat};
}

Patterns first_sampling_ver2(int ep, const std::vector<std::vector<std::string>> &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_changes, const std::vector<ll> &scores, auto &RD)
{

    std::vector<std::string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_change_pat;
    if (ep == 1)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }
    }
    else if (ep - 2 < C and num_use_change < C)
    {

        bool found = false;
        ret_change_pat = n_changes[ep - 1];
        for (int x = 1; x < X / 2 + 1 and !found; ++x)
        {
            for (int m = 0; m < M and !found; m++)
            {
                cerr << RD[1][m][x] << endl;
                if (RD[1][m][x] < epsilon and !fixed_op[m][2 * x])
                {
                    cerr << "m, x: " << m << " " << x << endl;
                    fixed_op[m][2 * x] = true;
                    found = true;
                    num_use_change++;
                    rep(i, M)
                    {
                        ret_pat[i] = string(2 * X, '9');
                        ret_enc_pat[i] = {'9', '9'};
                    }
                    ret_change_pat.push_back({{m, 2 * x}, '1'});
                    if (num_use_change < C)
                    {
                        ret_change_pat.push_back({{m, 2 * x + 1}, '1'});
                        fixed_op[m][2 * x + 1] = true;
                        num_use_change++;
                    }
                    for (auto p : ret_change_pat)
                    {
                        int pos_m = p.first.first;
                        for (int pos_x = p.first.second; pos_x < 2 * X; pos_x += 2)
                        {
                            ret_pat[pos_m][pos_x] = p.second;
                        }
                    }
                }
            }
        }
    }
    else
    {
        uniform_int_distribution<> randM(0, M - 1);
        uniform_int_distribution<> rand8(1, 8);
        uniform_int_distribution<> rand1(0, 1);

        int m = randM(mt);
        int wh = rand1(mt);
        int nx = rand8(mt);
        ret_pat = patterns[ep - 1];
        ret_enc_pat = basic_pat[ep - 1];
        if (wh)
            ret_enc_pat[m].first = '0' + nx;
        else
            ret_enc_pat[m].second = '0' + nx;
        ret_change_pat = n_changes[ep - 1];
    }
    return {ret_pat, ret_enc_pat, ret_change_pat};
}

// ver1
Patterns first_sampling(int ep, const std::vector<std::vector<std::string>> &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_changes, const std::vector<ll> &scores)
{
    std::vector<std::string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_change_pat;
    if (ep == 1)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '2');
            ret_enc_pat[i] = {'2', '2'};
        }
    }

    else if (ep == 2)
    {

        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '3');
            rep(j, X)
            {
                ret_pat[i][j * 2 + 1] = '4';
            }
            ret_enc_pat[i] = {'3', '4'};
        }
    }
    else if (ep == 3)
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '5');
            rep(j, X)
            {
                ret_pat[i][j * 2 + 1] = '6';
            }
            ret_enc_pat[i] = {'5', '6'};
        }
    else if (ep == 4)
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '6');
            rep(j, X)
            {
                    ret_pat[i][j * 2 + 1] = '7';
            }
                ret_enc_pat[i] = {'6', '7'};
        }

    else if (ep == 5)
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '7');
            rep(j, X)
            {
                    ret_pat[i][j * 2 + 1] = '8';
            }
                ret_enc_pat[i] = {'7', '8'};
        }
    else if (ep == 6)
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '8');
            rep(j, X)
            {
                ret_pat[i][j * 2 + 1] = '9';
            }
            ret_enc_pat[i] = {'8', '9'};
        }
    else if (ep == 7)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }
    }
    else if (ep == 8)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '1');
            rep(j, X) ret_pat[i][j * 2] = '9';
            ret_enc_pat[i] = {'9', '1'};
        }
    }
    else if (ep == 9)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            rep(j, X) ret_pat[i][j * 2] = '1';
            ret_enc_pat[i] = {'1', '9'};
        }
    }
    return {ret_pat, ret_enc_pat, ret_change_pat};
}

Patterns binary_search_sampling()
{
    Patterns ret;
    // TO DO
    return ret;
}

// temporary
void latin_hyper_cube_samplsng()
{
    return;
}

// temporary
std::vector<string> thompson_sampling()
{
    std::vector<std::string> ret;
    return ret;
}

bool check_pattern(const std::vector<std::string> &pat)
{
    int c = 0;
    rep(m, M)
    {
        char ap = pat[m][0], bp = pat[m][1];
        rep(x, X)
        {
            if (ap != pat[m][2 * x])
            {
                c++;
                ap = pat[m][2 * x];
            }
            if (bp != pat[m][2 * x + 1])
            {
                c++;
                bp = pat[m][2 * x + 1];
            }
            if (c > C)
                return false;
        }
    }
    return true;
}

int chtoi(char ch)
{
    return ch - '0';
}

// 探索点の候補 メモリが消えるバグはここになし。
Patterns modify_candidate(int ep, const EncodedPat &basic_pat, const EncodedChange &n_change, auto &V, auto &D, auto &RO, auto &RD, auto &cum_RO, auto &cum_RD)
{
    // TODO: implement
    EncodedPat pat = basic_pat;
    EncodedChange changes = n_change;
    uniform_int_distribution<> rand3(0, 2);
    uniform_int_distribution<> rand9(1, 9);
    uniform_int_distribution<> randM(0, M - 1);
    uniform_int_distribution<> randX(2, X - 1);
    uniform_int_distribution<> AorB(0, 1);
    uniform_int_distribution<> rand4(1, 5);
    // cerr << changes.size() << endl;
    assert(changes.size() <= C and changes.size() >= 0);
    // b=0: basicのupdate, b=1: changeのアップデート

    for (int b = 0; b < 2; b++)
    {
        int judge = rand3(mt); // 0: 1つ上のレベルの稼働にする, 1: 1つ下のレペルの稼働にする。 2: 何もしない。
        int m = randM(mt);
        int x = randX(mt);
        int wh = AorB(mt);
        int nx_val = rand9(mt);
        int range = rand4(mt);
        if (b == 0)
        {
            if (RD[ep][m][x] > 0.2) // TODO: cumで置き換えてみる
            {
                if (wh)
                    pat[m].first = '0' + min(chtoi(pat[m].first) + range, 9);
                else
                    pat[m].second = '0' + min(chtoi(pat[m].first) + range, 9);
            }
            else if (RD[ep][m][x] <= 0.2)
            {
                if (wh)
                    pat[m].first = '0' + max(chtoi(pat[m].first) - range, 1);
                else
                    pat[m].second = '0' + max(chtoi(pat[m].first) - range, 1);
            }
            // if (judge == 0)
            // {
            //     if (wh == 0)
            //     {
            //         pat[m].first = '0' + min(chtoi(pat[m].first) + 1 + rand3(mt), 9);
            //     }
            //     else
            //     {
            //         pat[m].second = '0' + min(chtoi(pat[m].second) + 1 + rand3(mt), 9);
            //     }
            // }
            // else if (judge == 1)
            // {
            //     if (wh == 0)
            //     {
            //         pat[m].first = '0' + max(chtoi(pat[m].first) - 1 - rand3(mt), 1);
            //     }
            //     else
            //     {
            //         pat[m].second = '0' + max(chtoi(pat[m].second) - 1 - rand3(mt), 1);
            //     }
            // }
        }
        else
        {
            bool flag = false;
            if (judge == 0)
            { //稼働変化をたす
                if (changes.size() == C)
                    continue;
                pair<pair<int, int>, int> to_change;
                if (RD[ep][m][x] > 0.2)
                {
                    if (wh)
                    {
                        to_change.first = (pii){m, 2 * x};
                        to_change.second = chtoi(basic_pat[m].first) - range;
                    }
                    else
                    {
                        to_change.first = (pii){m, 2 * x + 1};
                        to_change.second = chtoi(basic_pat[m].second) - range;
                    }
                    if (to_change.second <= 0 || to_change.second >= 10)
                        continue;
                    for (auto &change : changes)
                    {
                        if (change.first == to_change.first)
                        {
                            int mx = max({to_change.second, min(9, chtoi(change.second) + 1)});
                            change.second = '0' + mx;
                            flag = true;
                            break;
                        }
                    }
                }
                else
                {
                    if (wh)
                    {
                        to_change.first = (pii){m, 2 * x};
                        to_change.second = chtoi(basic_pat[m].first) - range;
                    }
                    else
                    {
                        to_change.first = (pii){m, 2 * x + 1};
                        to_change.second = chtoi(basic_pat[m].second) - 1;
                    }

                    if (to_change.second <= 0 || to_change.second >= 10)
                        continue;
                    for (auto &change : changes)
                    {
                        if (change.first == to_change.first)
                        {
                            int mn = min({to_change.second, max(1, chtoi(change.second) - 1)});
                            change.second = '0' + mn;
                            flag = true;
                            break;
                        }
                    }
                }
                if (flag)
                    break;

                else
                    changes.push_back({to_change.first, '0' + to_change.second});
            }
            else if (judge == 1)
            {
                if (changes.size() == 0)
                    continue;
                uniform_int_distribution<> randDelete(0, changes.size() - 1);
                int d_loc = randDelete(mt);
                // cerr << d_loc << " " << changes[d_loc].second << endl;
                changes.erase(changes.begin() + d_loc);
            }
        }
    }

    sort(changes.begin(), changes.end());

    // decode
    std::vector<string> ret_pat(M);
    rep(i, M)
    {
        string sa(X, pat[i].first), sb(X, pat[i].second);
        string tmp = "";
        rep(j, X)
        {
            tmp += sa[j];
            tmp += sb[j];
        }
        // cerr << tmp << endl;
        ret_pat[i] = tmp;
    }
    for (auto change : changes)
    {
        int m = change.first.first;
        int x = change.first.second;
        char val = change.second;
        // cerr << "m, x, val: " << m << " " << x << " " << val << endl;
        for (int i = x; i < 2 * X; i += 2)
            ret_pat[m][i] = val;
    }

    tuple<std::vector<string>, EncodedPat, EncodedChange> ret = {ret_pat, pat, changes};
    // cerr << ret[0].size() << endl;
    // for (auto s : ret)
    // {
    //     cerr << s << "\n";
    // }
    return ret;
}

Patterns roulette_select()
{
    Patterns ret;

    return ret;
}

// ベイズ最適化をやってみる解法
Patterns suggest(int ep, auto &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_change, const vector<ll> &scores, const std::vector<double> &eval_scores, const std::vector<ll> &exp_scores, auto &V, auto &D, auto &RO, auto &RD)
{
    if (ep <= 9)
    {
        Patterns ret = first_sampling(ep, patterns, basic_pat, n_change, scores);
        // Patterns ret = first_sampling_ver2(ep, patterns, basic_pat, n_change, scores, RD);
        // Patterns ret = first_sampling_ver3(ep, patterns, basic_pat, n_change, scores, RD);
        if (get<0>(ret)[0].size() == 2 * X)
            return ret;
    }
    clock_t suggest_start_time = clock();
    // vector<double> mu_points, var_points;
    vector<vector<string>> x_points(patterns.cbegin(), patterns.cbegin() + ep);
    dvector y_points(ep);

    rep(i, ep)
    {
        if (USE_LOG)
            y_points[i] = log10(1 + eval_scores[i]);
        else
            y_points[i] = eval_scores[i];
    }

    double y_max = (USE_LOG) ? log10(1 + best_eval_score) : best_eval_score;

    // 評価スコアの累積和とルーレット
    vector<double> score_table(ep + 1, 0.);
    rep(e, ep) score_table[e + 1] = score_table[e] + y_points[e];
    uniform_real_distribution<> roulette(0, score_table[ep]);

    // ROとRDの累積話 TODO: 後々改善して計算量を落とす.添字に注意
    vector<ll> cum_RO(M, 0);
    vector<double> cum_RD(M, 0.0);
    
    // > --------------------------- search loop -------------------------------<
    cur_time = clock();
    double best_acq_score = -1e20;
    vector<string> next_pattern;
    EncodedPat best_enc_pat;
    EncodedChange best_change;
    // cerr << "TIME: " << suggest_start_time << " " << cur_time << endl;
    // cerr << (double)(cur_time - suggest_start_time) / CLOCKS_PER_SEC << endl;
    int loop = 0;
    while ((double)(cur_time - suggest_start_time) / CLOCKS_PER_SEC < SEARCH_TIME)
    {
        loop++;
        double val = roulette(mt);
        int idx = lower_bound(score_table.begin(), score_table.end(), val) - score_table.begin();
        idx--;
        EncodedPat b_pat = basic_pat[idx];
        EncodedChange c_pat = n_change[idx];

        // cerr << "epoch: " << ep << ", id: " << idx << ", val: " << val << ", max: " << score_table[ep] << " " << patterns[idx][0][0] << " " << basic_pat[idx][0].first << endl;
        Patterns next_candidate = modify_candidate(idx, b_pat, c_pat, V, D, RO, RD, cum_RO, cum_RD); // TODO: 次の候補点の選び方が大切。ここの実装が大事。重点的に取り組む
        // for(auto np: get<0>(next_candidate)) cerr << np << "\n";
        if (pattern_set.find(get<0>(next_candidate)) != pattern_set.end())
            continue;
        pair<double, double> gp = gaussian_process(x_points, y_points, get<0>(next_candidate), ep);
        // mu_points.push_back(gp.first);
        // var_points.push_back(gp.second);

        // cerr << "acq: " << endl;
        double acq_score = EI(gp.first, gp.second, y_max);
        // cerr << "done:" << endl;
        // double acq_score = 1;
        if (acq_score > -1e10)
        {
            best_acq_score = acq_score;
            next_pattern = get<0>(next_candidate);
            best_enc_pat = get<1>(next_candidate);
            best_change = get<2>(next_candidate);
        }
        cur_time = clock();
        break;
        // cerr << "loop end " << endl;
    }
    return {next_pattern, best_enc_pat, best_change};
}

void calc_expected_score(int epoch, const std::vector<std::vector<std::string>> &patterns, std::vector<ll> &exp_scores)
{
    ll cost = 0;
    for (int i = 0; i < M; i++)
    {
        for (int x = 0; x < X; x++)
        {
            int a = patterns[epoch][i][2 * x] - '0';
            int b = patterns[epoch][i][2 * x + 1] - '0';
            cost += productLine[i].costA[a] + productLine[i].costB[b];
        }
    }
    exp_scores[epoch] = round(1e9 * (10. - log10((double)(cost / X))));
}

void get_result(int epoch, std::vector<ll> &scores, std::vector<ll> &V, std::vector<ll> &D, auto &RO, auto &RD, auto &eval_scores, auto &exp_scores)
{
    ll score, v, d;
    cin >> score >> v >> d;
    // cerr << score << " " << v << " " << d << "\n";
    scores[epoch] = score;
    eval_scores[epoch] = (double)exp_scores[epoch] * pow(DECAY, d);
    // if (epoch == 5)
    // eval_scores[epoch] *= 0.5;
    V[epoch] = v;
    D[epoch] = d;
    for (int m = 0; m < M; ++m)
    {
        for (int x = 0; x < X; ++x)
        {
            double rd, ro;
            cin >> rd >> ro;
            RD[epoch][m][x] = rd;
            RO[epoch][m][x] = ro;
            // cerr << "m, x: " << RD[epoch][m][x] << " " << RO[epoch][m][x] << endl;
        }
    }
    return;
}

void show_result(int epoch, const std::vector<std::vector<std::string>> &pat)
{
    for (int m = 0; m < M; m++)
    {
        cout << pat[epoch][m] << "\n";
    }
    return;
}

// コードが汚いけれど...とりあえず Done is better than perfect
int main()
{
    start_time = clock();
    input(X, M, C, E, productLine);
    std::vector<ll> V(E + 1, 0), D(E + 1, 0);
    ll best_score = -1;
    ll score = -1;
    std::vector<std::vector<std::vector<double>>> RD(E + 1, vector<vector<double>>(M, vector<double>(X, 0.0)));
    std::vector<std::vector<std::vector<ll>>> RO(E + 1, vector<vector<ll>>(M, vector<ll>(X, 0)));
    string ini_pat(2 * X, '1');
    std::vector<std::vector<std::string>> patterns(E + 1, std::vector<std::string>(M, ini_pat)); // 稼働パターン
    vector<string>best_patterns;
    std::vector<EncodedPat> basic_patterns(E + 1, EncodedPat(M, {'1', '1'}));
    std::vector<EncodedChange> change_patterns(E + 1); // あるepochの稼働パターンの時の稼動変更位置(m, x)
    std::vector<ll> scores(E + 1, 0);                  // 計画プログラムから帰ってくるスコア
    std::vector<double> eval_scores(E + 1, 0.);        // 評価に用いるスコア(ベイズ最適化とか) 納期が遅れまくりだと小さくなる。
    std::vector<ll> expected_scores(E + 1, 0);         // 納期が全て間に合うと仮定した時の期待されるスコア

    std::cerr << "***************************************** START **************************************\n";
    // Init
    init(scores, eval_scores, expected_scores);
    // >-------------------------- reaction loop ----------------------------<
    for (int epoch = 1; epoch <= E; epoch++)
    {
        cerr << "########## EPOCH: " << setfill('0') << setw(3) << epoch << " ##########\n";
        cur_time = clock();
        cerr << "ELAPSED: " << (double)(cur_time - start_time) / CLOCKS_PER_SEC << "\n";
        Patterns next_pattern = suggest(epoch, patterns, basic_patterns, change_patterns, scores, eval_scores, expected_scores, V, D, RO, RD);
        for (auto np : get<0>(next_pattern))
            cerr << np << "\n";
        patterns[epoch] = get<0>(next_pattern);
        pattern_set.insert(get<0>(next_pattern));
        basic_patterns[epoch] = get<1>(next_pattern);
        change_patterns[epoch] = get<2>(next_pattern);
        calc_expected_score(epoch, patterns, expected_scores);
        show_result(epoch, patterns);
        get_result(epoch, scores, V, D, RO, RD, eval_scores, expected_scores);

        // update
        if (eval_scores[epoch] > best_eval_score)
        {
            best_eval_score = eval_scores[epoch];
            best_eval_pattern = patterns[epoch];
        }
        cerr << "SCORE: " << scores[epoch] << ", EVAL: " << eval_scores[epoch] << ", EXP: " << expected_scores[epoch] << "\n";

#ifdef TEST
#endif
    }
    // >-------------------------- reaction loop ----------------------------<

#ifdef DEBUG
    for(int x = 0; x < X; x++) {
        cerr << "WEEK: " << x+1 << endl;
        for(int m = 0; m < M; m++) {
            cerr << "MACHINE: " << m << ", RD: "<< RD[1][x][m] << ", RO: " << RO[1][x][m] << endl;
        }
    }
    cerr << "cost list" << endl;
    for(auto p: heavy_cost_order) {
        cerr << "(machine, cost): " << p.first + 1 << " " << p.second << endl; 
    }

    cerr << "cost" << endl;
    for(int i = 0; i < M; i++) {
        cerr << "### MACHINE: " << i + 1 << endl;
        for(int j = 1; j <= 9; j++){
            cerr << "A,B: " << productLine[i].costA[j] << " " << productLine[i].costB[j] << endl;
        }
    }

    for(auto s: best_eval_pattern) {
        cerr << s << endl;
    } 
#endif
    return 0;
}