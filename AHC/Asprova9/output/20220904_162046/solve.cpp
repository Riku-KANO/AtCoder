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

// #define DEBUG
// #define TEST
// #define MAIN1
#define MAIN2

#define SHOW_EACH_PATTERN


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
    std::vector<ll> cost_order;
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
bool USE_QUADRATIC;
int X, M, C, E;
std::vector<Operation> productLine(21);
std::vector<ll> V(301, 0), D(301, 0);
std::vector<std::vector<std::vector<double>>> RD(301, std::vector<std::vector<double>>(20, std::vector<double>(16, 0.0)));
std::vector<std::vector<std::vector<ll>>> RO(301, std::vector<std::vector<ll>>(20, std::vector<ll>(16, 0)));

clock_t start_time;
clock_t cur_time;
ll best_score = -10;
double best_eval_score = -1e20;
ll best_basic_opt_score = 0;
int best_score_epoch = -1;
int best_eval_epoch = -1;
std::vector<std::string> best_score_pattern;
std::vector<std::string> best_eval_pattern;
std::vector<std::string> best_basic_opt_pattern;
std::set<std::vector<std::string>> pattern_set;
std::vector<double> score_table(301, 0.0);
std::vector<std::vector<bool>> fixed_op;
int num_use_change = 0;

const int HOUR = 86'400;
const std::vector<int> operation_time = {0, 3, 5, 8, 10, 12, 14, 16, 18};
const std::vector<double> operation_time_ratio = {0.0, 1. / 6, 5. / 18 / 4. / 9, 5. / 9, 6. / 9, 7. / 9, 8. / 9, 1.};
std::vector<pair<double, pair<char, char>>> operation_time_ratio_double;

std::vector<std::pair<int, ll>> heavy_cost_order; // 機械に応じてコストの差がありすぎるのでここで管理する。
std::vector<std::pair<std::pair<int,char>, ll>> heavy_cost_order_AB; 

std::vector<std::set<int>> machine_relation(20);
////////////////////////
class Config
{
public:
    double tl = 5.0;
    double decay = 0.6;
    bool use_log = false;
    bool use_quadratic = true;
    std::string file_name = "../output/config.txt";
    Config() {}
    void init()
    {
        TL = this->tl;
        TIME_LIMIT_90 = this->tl * 0.95;
        DECAY = this->decay;
        USE_LOG = this->use_log;
        USE_QUADRATIC = this->use_quadratic;
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
        cout << "USE LOG: " << this->use_log << std::endl;
        cout << "USE QUADRATIC: " << this->use_quadratic << std::endl;
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

int chtoi(char ch)
{
    return ch - '0';
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
    std::cerr << "\n####### INIT START ######" << std::endl;
    Config config;
    config.init();

    // initial score for "11111...111" size = 2*X * M
    rep(i, M) exp_scores[0] += productLine[i].costA[1] + productLine[i].costB[1];

    scores[0] = 0;
    eval_scores[0] = 0; //評価値も0にしておく。
    best_eval_pattern = std::vector<std::string>(M, std::string(2 * X, '1'));
    SEARCH_TIME = TIME_LIMIT_90 / (double)E;
    std::vector<std::vector<bool>> fixed_tmp(M, std::vector<bool>(2 * X, false));
    fixed_op = fixed_tmp;
    //ルーレットのスコアテーブル
    score_table.resize(E + 1);

    // 重たいコストを管理して９番目の稼働のコスト値でソートする。
    rep(i, M)
    {
        heavy_cost_order.push_back({i, productLine[i].costA[9]});
        heavy_cost_order_AB.push_back({{i, 'A'}, productLine[i].costA[9]});
        heavy_cost_order_AB.push_back({{i, 'B'}, productLine[i].costB[9]});
    }
    sort(heavy_cost_order.begin(), heavy_cost_order.end(), [](pair<int, ll> a, pair<int, ll> b)
         { return a.second > b.second; });

    sort(heavy_cost_order_AB.begin(), heavy_cost_order_AB.end(), [](pair<pair<int, char>,ll>&a, pair<pair<int,char>, ll>&b)
         { return a.second > b.second;});

    // operation_time_ratio_double
    for (int i = 1; i <= 9; i++)
    {
        for (int j = 1; j <= 9; j++)
        {
            double ratio = ((double)operation_time[i - 1] * 5 + (double)operation_time[j - 1] * 2) / 126.0;
            operation_time_ratio_double.push_back({ratio, {'0' + i, '0' + j}});
        }
    }
    sort(operation_time_ratio_double.begin(), operation_time_ratio_double.end());

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

// #制約
// 10 <= M <= 20
// 8 <= X <= 16
// 2 <= C <= 8
// Eは50, 100, 300のいずれか
// RDが負荷率 ROが遅れ作業の個数
//
// アイデア
// ・ブラックボックス最適化・ベイズ最適化？+ 貪欲<=最有力
//

// その他
// C++の E=300の時のIOの所要時間は0.1~0.2s程度。ちなみにPythonのIOは4.5秒程度。遅すぎるのでnumpyの行列演算ライブラリは意味をなさない。
// 時間的な制約の関係でループ時間が短すぎると次の候補が見つからない。初期点サンプリングを多めにやってもいいかもしれない。<- 嘘かもしれない。時間のせいではない
// もしくは次の候補を貪欲的に選ぶのもありかもしれない。
// ベイズ最適化を実行してみたがE=300で1ループ１０秒以上かかる時がある。計算量がきつい。
// E = 300の時は実行してみてわかったがルーレットをやっているだけで3Gを超える時がある。ベイズ最適化してる暇は最初にしかない。
// 0000.txtにおいて、ルーレット+ベイズ最適化 で2.65G, ルーレットだけは2.4G。多少なりの効果はあり。
// 0000.txtで対数を使わないベイズの方が2.6Gを上回る傾向が強い。
// 色々やってみた結果ベイズ最適化はあまり使えなさそう。。

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
        // std::cerr << i << " OK\n";
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
    //     std::cerr << "aux: " << aux[i] << " " << aux2[i] << "\n";
    // }
    aux3 = ublas::prod(ker_inv2, aux);
    // rep(i, ker_size) rep(j, ker_size) aux3[i] += ker_inv2(i, j) * aux[j]; // TODO: debug不足  要チェック
    double var = my_kernel(next_point, next_point) - ublas::inner_prod(aux, aux3);
    // rep(i, ker_size) var -= aux3[i] * aux[i];
    ret = {mu, var};
    // std::cerr << mu << " " << var << std::endl;
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
Patterns first_sampling_ver3(int ep, const std::vector<std::vector<std::string>> &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_changes, const std::vector<ll> &scores)
{
    std::vector<Patterns> ret;
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
    if (ep == 2)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }

        bool flag = true;
        int h_id = 0;
        while (flag)
        {
            int m_id = heavy_cost_order[h_id].first;
            int change_week;
            for (int x = X - 1; x >= 0; x--)
            {
                if (RD[1][m_id][x] > epsilon)
                {
                    change_week = x + 1;
                    break;
                }
            }
            ret_change_pat.push_back({{m_id, change_week * 2}, '1'});
            num_use_change++;
            if (num_use_change < C)
            {
                ret_change_pat.push_back({{m_id, change_week * 2 + 1}, '1'});
                num_use_change++;
            }
            h_id++;
            if (h_id == M or num_use_change == C)
                break;
        }
        // decode
        for (auto p : ret_change_pat)
        {
            int m = p.first.first;
            for (int x = p.first.second; x < 2 * X; x += 2)
            {
                ret_pat[m][x] = p.second;
            }
        }
    }
    return {ret_pat, ret_enc_pat, ret_change_pat};
}

Patterns first_sampling_ver2(int ep, const std::vector<std::vector<std::string>> &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_changes, const std::vector<ll> &scores)
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
                std::cerr << RD[1][m][x] << std::endl;
                if (RD[1][m][x] < epsilon and !fixed_op[m][2 * x])
                {
                    std::cerr << "m, x: " << m << " " << x << std::endl;
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
std::vector<Patterns> first_sampling(int ep, const std::vector<std::vector<std::string>> &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_changes, const std::vector<ll> &scores)
{
    std::vector<Patterns> ret;
    std::vector<std::string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_change_pat;
    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '2');
        ret_enc_pat[i] = {'2', '2'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '3');
        rep(j, X)
        {
            ret_pat[i][j * 2 + 1] = '4';
        }
        ret_enc_pat[i] = {'3', '4'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '5');
        rep(j, X)
        {
            ret_pat[i][j * 2 + 1] = '6';
        }
        ret_enc_pat[i] = {'5', '6'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '6');
        rep(j, X)
        {
            ret_pat[i][j * 2 + 1] = '7';
        }
        ret_enc_pat[i] = {'6', '7'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '7');
        rep(j, X)
        {
            ret_pat[i][j * 2 + 1] = '8';
        }
        ret_enc_pat[i] = {'7', '8'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '8');
        rep(j, X)
        {
            ret_pat[i][j * 2 + 1] = '9';
        }
        ret_enc_pat[i] = {'8', '9'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = std::string(2 * X, '9');
        ret_enc_pat[i] = {'9', '9'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = std::string(2 * X, '1');
        rep(j, X) ret_pat[i][j * 2] = '9';
        ret_enc_pat[i] = {'9', '1'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

    rep(i, M)
    {
        ret_pat[i] = string(2 * X, '9');
        rep(j, X) ret_pat[i][j * 2] = '1';
        ret_enc_pat[i] = {'1', '9'};
    }
    ret.push_back({ret_pat, ret_enc_pat, ret_change_pat});

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

// 変え方の２パターン目。２つのマシーンを選んで変更をさせる。蜜柑星
Patterns modify_candidate2(int ep, const EncodedPat &basic_pat, const EncodedChange &change_pattern)
{
    EncodedPat pat = basic_pat;
    EncodedChange change_pat = change_pattern;
    uniform_int_distribution<> rand3(0, 2);
    uniform_int_distribution<> rand9(1, 9);
    uniform_int_distribution<> randM(0, M - 1);
    uniform_int_distribution<> randX(2, X - 1);
    uniform_int_distribution<> AorB(0, 1);
    uniform_int_distribution<> rand4(1, 5);
    assert(change_pat.size() <= C and change_pat.size() >= 0);

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
        }
        else
        {
            bool flag = false;
            if (judge == 0)
            { //稼働変化をたす
                if (change_pat.size() == C)
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
                    for (auto &change : change_pat)
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
                    for (auto &change : change_pat)
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
                    change_pat.push_back({to_change.first, '0' + to_change.second});
            }
            else if (judge == 1)
            {
                if (change_pat.size() == 0)
                    continue;
                uniform_int_distribution<> randDelete(0, change_pat.size() - 1);
                int d_loc = randDelete(mt);
                // std::cerr << d_loc << " " << change_pat[d_loc].second << std::endl;
                change_pat.erase(change_pat.begin() + d_loc);
            }
        }
    }

    sort(change_pat.begin(), change_pat.end());

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
        // std::cerr << tmp << std::endl;
        ret_pat[i] = tmp;
    }
    for (auto change : change_pat)
    {
        int m = change.first.first;
        int x = change.first.second;
        char val = change.second;
        // std::cerr << "m, x, val: " << m << " " << x << " " << val << std::endl;
        for (int i = x; i < 2 * X; i += 2)
            ret_pat[m][i] = val;
    }

    tuple<std::vector<string>, EncodedPat, EncodedChange> ret = {ret_pat, pat, change_pat};
    // std::cerr << ret[0].size() << std::endl;
    // for (auto s : ret)
    // {
    //     std::cerr << s << "\n";
    // }
    return ret;
}

// 探索点の候補 メモリが消えるバグはここになし。
Patterns modify_candidate(int ep, const EncodedPat &basic_pat, const EncodedChange &change_pattern, auto &V, auto &D)
{
    // TODO: implement
    EncodedPat pat = basic_pat;
    EncodedChange change_pat = change_pattern;
    uniform_int_distribution<> rand3(0, 2);
    uniform_int_distribution<> rand9(1, 9);
    uniform_int_distribution<> randM(0, M - 1);
    uniform_int_distribution<> randX(2, X - 1);
    uniform_int_distribution<> AorB(0, 1);
    uniform_int_distribution<> rand4(1, 5);
    // std::cerr << changes.size() << std::endl;
    assert(change_pat.size() <= C and change_pat.size() >= 0);
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
        }
        else
        {
            bool flag = false;
            if (judge == 0)
            { //稼働変化をたす
                if (change_pat.size() == C)
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
                    for (auto &change : change_pat)
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
                    for (auto &change : change_pat)
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
                    change_pat.push_back({to_change.first, '0' + to_change.second});
            }
            else if (judge == 1)
            {
                if (change_pat.size() == 0)
                    continue;
                uniform_int_distribution<> randDelete(0, change_pat.size() - 1);
                int d_loc = randDelete(mt);
                // std::cerr << d_loc << " " << change_pat[d_loc].second << std::endl;
                change_pat.erase(change_pat.begin() + d_loc);
            }
        }
    }

    sort(change_pat.begin(), change_pat.end());

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
        // std::cerr << tmp << std::endl;
        ret_pat[i] = tmp;
    }
    for (auto change : change_pat)
    {
        int m = change.first.first;
        int x = change.first.second;
        char val = change.second;
        // std::cerr << "m, x, val: " << m << " " << x << " " << val << std::endl;
        for (int i = x; i < 2 * X; i += 2)
            ret_pat[m][i] = val;
    }

    tuple<std::vector<string>, EncodedPat, EncodedChange> ret = {ret_pat, pat, change_pat};
    // std::cerr << ret[0].size() << std::endl;
    // for (auto s : ret)
    // {
    //     std::cerr << s << "\n";
    // }
    return ret;
}

Patterns roulette_select()
{
    Patterns ret;

    return ret;
}

// ベイズ最適化をやってみる解法
Patterns suggest(int ep, auto &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &n_change, const std::vector<ll> &scores, const std::vector<double> &eval_scores, const std::vector<ll> &exp_scores, auto &V, auto &D)
{
    if (ep <= 9)
    {
        // Patterns ret = first_sampling(ep, patterns, basic_pat, n_change, scores);
        Patterns ret = first_sampling_ver2(ep, patterns, basic_pat, n_change, scores);
        // Patterns ret = first_sampling_ver3(ep, patterns, basic_pat, n_change, scores, RD);
        // if (get<0>(ret)[0].size() == 2 * X)
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
    // std::cerr << "TIME: " << suggest_start_time << " " << cur_time << std::endl;
    // std::cerr << (double)(cur_time - suggest_start_time) / CLOCKS_PER_SEC << std::endl;
    int loop = 0;
    while ((double)(cur_time - suggest_start_time) / CLOCKS_PER_SEC < SEARCH_TIME)
    {
        loop++;
        double val = roulette(mt);
        int idx = lower_bound(score_table.begin(), score_table.end(), val) - score_table.begin();
        idx--;
        EncodedPat b_pat = basic_pat[idx];
        EncodedChange c_pat = n_change[idx];

        // std::cerr << "epoch: " << ep << ", id: " << idx << ", val: " << val << ", max: " << score_table[ep] << " " << patterns[idx][0][0] << " " << basic_pat[idx][0].first << std::endl;
        Patterns next_candidate = modify_candidate(idx, b_pat, c_pat, V, D, cum_RO, cum_RD); // TODO: 次の候補点の選び方が大切。ここの実装が大事。重点的に取り組む
        // for(auto np: get<0>(next_candidate)) std::cerr << np << "\n";
        if (pattern_set.find(get<0>(next_candidate)) != pattern_set.end())
            continue;
        pair<double, double> gp = gaussian_process(x_points, y_points, get<0>(next_candidate), ep);
        // mu_points.push_back(gp.first);
        // var_points.push_back(gp.second);

        // std::cerr << "acq: " << std::endl;
        double acq_score = EI(gp.first, gp.second, y_max);
        // std::cerr << "done:" << std::endl;
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
        // std::cerr << "loop end " << std::endl;
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

void get_result(int epoch, std::vector<ll> &scores, std::vector<ll> &V, std::vector<ll> &D, auto &eval_scores, auto &exp_scores)
{
    ll score, v, d;
    cin >> score >> v >> d;
    // std::cerr << score << " " << v << " " << d << "\n";
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
            // std::cerr << "m, x: " << RD[epoch][m][x] << " " << RO[epoch][m][x] << std::endl;
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

#ifdef MAIN1
int main()
{
    start_time = clock();
    input(X, M, C, E, productLine);
    std::vector<ll> V(E + 1, 0), D(E + 1, 0);
    ll best_score = -1;
    ll score = -1;

    string ini_pat(2 * X, '1');
    std::vector<std::vector<std::string>> patterns(E + 1, std::vector<std::string>(M, ini_pat)); // 稼働パターン
    std::vector<std::string> best_patterns;
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
        std::cerr << "########## EPOCH: " << setfill('0') << setw(3) << epoch << " ##########\n";
        cur_time = clock();
        std::cerr << "ELAPSED: " << (double)(cur_time - start_time) / CLOCKS_PER_SEC << "\n";
        Patterns next_pattern = suggest(epoch, patterns, basic_patterns, change_patterns, scores, eval_scores, expected_scores, V, D);
        for (auto np : get<0>(next_pattern))
            std::cerr << np << "\n";
        patterns[epoch] = std::get<0>(next_pattern);
        pattern_set.insert(std::get<0>(next_pattern));
        basic_patterns[epoch] = std::get<1>(next_pattern);
        change_patterns[epoch] = std::get<2>(next_pattern);
        calc_expected_score(epoch, patterns, expected_scores);
        show_result(epoch, patterns);
        get_result(epoch, scores, V, D, eval_scores, expected_scores);

        // update
        if (eval_scores[epoch] > best_eval_score)
        {
            best_eval_score = eval_scores[epoch];
            best_eval_pattern = patterns[epoch];
        }
        std::cerr << "SCORE: " << scores[epoch] << ", EVAL: " << eval_scores[epoch] << ", EXP: " << expected_scores[epoch] << "\n";
    }
    // >-------------------------- reaction loop ----------------------------<

#ifdef DEBUG
    for (int x = 0; x < X; x++)
    {
        std::cerr << "WEEK: " << x + 1 << std::endl;
        for (int m = 0; m < M; m++)
        {
            std::cerr << "MACHINE: " << m << ", RD: " << RD[1][x][m] << ", RO: " << RO[1][x][m] << std::endl;
        }
    }
    std::cerr << "cost list" << std::endl;
    for (auto p : heavy_cost_order)
    {
        std::cerr << "(machine, cost): " << p.first + 1 << " " << p.second << std::endl;
    }

    std::cerr << "cost" << std::endl;
    for (int i = 0; i < M; i++)
    {
        std::cerr << "### MACHINE: " << i + 1 << std::endl;
        for (int j = 1; j <= 9; j++)
        {
            std::cerr << "A,B: " << productLine[i].costA[j] << " " << productLine[i].costB[j] << std::endl;
        }
    }

    for (auto s : best_eval_pattern)
    {
        std::cerr << s << std::endl;
    }
#endif
    return 0;
}
#endif

// ガウス過程かいきのsuggest。重たいので使わない。
std::vector<Patterns> suggest2(int ep, auto &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &change_pat, const std::vector<ll> &scores, const std::vector<double> &eval_scores, const std::vector<ll> &exp_scores, auto &V, auto &D)
{
    if (ep == 1)
    {
        std::vector<Patterns> ret = first_sampling(ep, patterns, basic_pat, change_pat, scores);
        // std::vector<Patterns> ret = first_sampling_ver2(ep, patterns, basic_pat, change_pat, scores, RD);
        // std::vector<Patterns> ret = first_sampling_ver3(ep, patterns, basic_pat, change_pat, scores, RD);
        // if (get<0>(ret)[0].size() == 2 * X)
        return ret;
    }
    std::vector<Patterns> ret;
    clock_t suggest_start_time = clock();
    // vector<double> mu_points, var_points;
    std::vector<std::vector<std::string>> x_points(patterns.cbegin(), patterns.cbegin() + ep);
    dvector y_points(ep);

    rep(i, ep)
    {
        if (USE_LOG)
            y_points[i] = log10(1 + eval_scores[i]);
        else if (USE_QUADRATIC)
            y_points[i] = eval_scores[i] * eval_scores[i];
        else
            y_points[i] = eval_scores[i];
    }

    double y_max = (USE_LOG) ? log10(1 + best_eval_score) : best_eval_score;

    // 評価スコアの累積和とルーレット
    std::vector<double> score_table(ep + 1, 0.);
    rep(e, ep) score_table[e + 1] = score_table[e] + y_points[e];
    uniform_real_distribution<> roulette(0, score_table[ep]);

    // ROとRDの累積話 TODO: 後々改善して計算量を落とす.添字に注意
    std::vector<ll> cum_RO(M, 0);
    std::vector<double> cum_RD(M, 0.0);

    // > --------------------------- search loop -------------------------------<
    cur_time = clock();
    double best_acq_score = -1e20;
    std::vector<std::string> next_pattern;
    EncodedPat best_enc_pat;
    EncodedChange best_change;
    // std::cerr << "TIME: " << suggest_start_time << " " << cur_time << std::endl;
    // std::cerr << (double)(cur_time - suggest_start_time) / CLOCKS_PER_SEC << std::endl;
    int loop = 0;
    while ((double)(cur_time - suggest_start_time) / CLOCKS_PER_SEC < SEARCH_TIME)
    {
        loop++;
        double val = roulette(mt);
        int idx = lower_bound(score_table.begin(), score_table.end(), val) - score_table.begin();
        idx--;
        EncodedPat b_pat = basic_pat[idx];
        EncodedChange c_pat = change_pat[idx];

        Patterns next_candidate = modify_candidate(idx, b_pat, c_pat, V, D, cum_RO, cum_RD); // TODO: 次の候補点の選び方が大切。ここの実装が大事。重点的に取り組む
        if (pattern_set.find(get<0>(next_candidate)) != pattern_set.end())
            continue;
        pair<double, double> gp = gaussian_process(x_points, y_points, get<0>(next_candidate), ep);

        double acq_score = EI(gp.first, gp.second, y_max);
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
        // std::cerr << "loop end " << std::endl;
    }
    Patterns tmp = {next_pattern, best_enc_pat, best_change};
    ret.push_back(tmp);
    return ret;
}

#ifdef MAIN2

Patterns greedy_change(auto &patterns, auto &opt_pat)
{
    std::vector<std::pair<ll, std::pair<int, int>>> exp_improve_score;
    std::cerr << "BEST epoch: " << best_score_epoch << std::endl;
    for (int m = 0; m < M; m++)
    {
        for (int x = X - 1; x >= 0; x--)
        {
            if (RD[best_score_epoch][m][x] > 0.1)
            {
                if (x == X - 1)
                {
                    exp_improve_score.push_back({0, {m, 2 * X}});
                    exp_improve_score.push_back({0, {m, 2 * X + 1}});
                    break;
                }
                ll base_a=productLine[m].costA[opt_pat[m].first-'0'];
                ll base_b=productLine[m].costB[opt_pat[m].second-'0'];

                exp_improve_score.push_back({(X - x - 1) * (base_a - productLine[m].costA[1]), {m, 2 * x + 2}});
                exp_improve_score.push_back({(X - x - 1) * (base_b - productLine[m].costB[1]), {m, 2 * x + 3}});

                // upper_bound  todo
                int time_ida = opt_pat[m].first-'1';
                int time_idb = opt_pat[m].second-'1';
                double total_hour=5.*operation_time[time_ida]+ 2.*operation_time[time_idb];
                int cur_pat_a = opt_pat[m].first-'0';
                int cur_pat_b = opt_pat[m].second-'0';
                bool oka=false;
                bool okb=false;
                // 今はaとbが同じと考えているので同じループで処理
                for(int start_pat = cur_pat_a; start_pat >= 1; start_pat--) {
                    double total_hour_a = 5.0*operation_time[start_pat-1] + 2.0*operation_time[time_idb];
                    double total_hour_b = 5.0*operation_time[time_ida] + 2.0*operation_time[start_pat-1];
                    //std::cerr << "pat: " << start_pat << " " << total_hour_a << " " << total_hour_b << " " << total_hour << endl;
                    if(total_hour_a/total_hour < RD[best_score_epoch][m][x] and !oka) {
                        int pat_a = min(start_pat+2, cur_pat_a);
                        //std::cerr << "machine: " << m+1 << std::endl;
                        //std::cerr << "THIS is: " << (ll)(X-x)*(base_a-productLine[m].costA[pat_a]) << std::endl;
                        //std::cerr << "Versus: " << (ll)(X-x-1)*(base_a-productLine[m].costA[1]) << std::endl;
                        exp_improve_score.push_back({(ll)(X-x)*(base_a-productLine[m].costA[pat_a]), {m, 2*x}});
                        oka=true;
                    }

                    if(total_hour_b/total_hour < RD[best_score_epoch][m][x] and !okb) {
                        int pat_b = min(start_pat+2, cur_pat_b);
                        exp_improve_score.push_back({(ll)(X-x)*(base_b-productLine[m].costB[pat_b]), {m, 2*x+1}});
                        okb=true;
                    }
                }

                break;
            }
        }
    }
    sort(exp_improve_score.rbegin(), exp_improve_score.rend());
    std::vector<std::string> pat = best_score_pattern;
    EncodedPat enc_pat = opt_pat;
    EncodedChange enc_change;
    for (int c = 0; c < C; c++)
    {
        int m = exp_improve_score[c].second.first;
        int x = exp_improve_score[c].second.second;
        std::cerr << "next change: " << m << " " << x << " , score: " << exp_improve_score[c].first << std::endl;
        enc_change.push_back({{m, x}, '1'});
    }

    // decoding
    for (auto change : enc_change)
    {
        int m = change.first.first;
        int start = change.first.second;
        char ch = change.second;
        for (int x = start; x < 2 * X; x += 2)
        {
            pat[m][x] = ch;
        }
        // std::cerr << pat[m] << std::endl;
    }

    return {pat, enc_pat, enc_change};
}

// 重要そうなミスのパターン
// ある稼働を変化させたときに、それがどのように周りに変化を起こすか。
// 1. 他に影響を及ぼすパターン
// 作業が遅すぎることを示唆。特に1週目のdelayが大事そう。as
// 2. 影響を及ぼさないパターン
// 自身のレーンで問題が発生しているパターン->そのレーンだけをうまいこと最適化
void analyze_failure(int ep, auto &patterns, auto &basic_pat, auto &opt_pat, auto &scores, int pre_target, pii pre_val)
{
    for (int m = 0; m < M; m++)
    {
        for (int x = 0; x < X; x++)
        {
            if (RD[ep - 1][m][x] > 0)
            {
                machine_relation[pre_target].insert(m);
            }
        }
    }
}

// 稼働率の最大値を返す
double get_max_operation_rate(int ep, int m_id)
{
    double operation_rate_max = -1.0;
    for (int x = 0; x < X; x++)
    {
        operation_rate_max = max(RD[ep - 1][m_id][x], operation_rate_max);
    }
    return operation_rate_max;
}

Patterns binary_suggest(int ep, auto &patterns, auto &basic_pat, auto &opt_pat, auto &scores, int &pre_target, pii &pre_val, std::vector<bool> &visited)
{
    std::vector<string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_enc_change;
    if (ep == 1)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '8');
            ret_enc_pat[i] = {'8', '8'};
        }
        return {ret_pat, ret_enc_pat, ret_enc_change};
    }
    else if (ep == 2)
    {
        // 88888..88でもうまくいかなかった場合
        if (scores[ep - 1] == 0)
        {
            rep(i, M)
            {
                ret_pat[i] = string(2 * X, '9');
                ret_enc_pat[i] = {'9', '9'};
                opt_pat[i] = {'9', '9'};
            }
            pre_val = {9, 9};
            pre_target = heavy_cost_order[0].first;
            return {ret_pat, ret_enc_pat, ret_enc_change};
        }
        else
        {
            int m_start = 0;
            int m_id = heavy_cost_order[m_start].first;
            pair<int, int> nx_val = {pre_val.first / 2, pre_val.second / 2};
            visited[nx_val.first] = true;
            rep(i, M)
            {
                if (i == m_id)
                {
                    char nx = '0' + nx_val.first;
                    ret_pat[i] = std::string(2 * X, nx);
                    ret_enc_pat[i] = {nx, nx};
                }
                else
                {
                    ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                    ret_enc_pat[i] = opt_pat[i];
                }
            }
            pre_target = m_id;
            pre_val = {pre_val.first / 2, pre_val.second / 2};
            return {ret_pat, ret_enc_pat, ret_enc_change};
        }
    }
    else
    {
        bool go_next = false;
        if (scores[ep - 1] == 0)
        {
            pii nx_val = {pre_val.first + 1, pre_val.second + 1};
            if (nx_val >= opt_pat[pre_target].first - '0')
            {
                go_next = true;
            }
            else
            {
                rep(i, M)
                {
                    if (i == pre_target)
                    {
                        char nx = '0' + nx_val.first;
                        ret_pat[i] = std::string(2 * X, nx);
                        ret_enc_pat[i] = {nx, nx};
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                pre_val = nx_val;
                visited[nx_val.first] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
        else
        {
            char opt_val = '0' + pre_val.first;
            opt_pat[pre_target] = {opt_val, opt_val};
            pii nx_val = {pre_val.first - 1, pre_val.second - 1};
            if (visited[nx_val.first])
            {
                go_next = true;
            }
            else
            {
                rep(i, M)
                {
                    if (i == pre_target)
                    {
                        char nx = '0' + nx_val.first;
                        ret_pat[i] = std::string(2 * X, nx);
                        ret_enc_pat[i] = {nx, nx};
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                pre_val = nx_val;
                visited[nx_val.first] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
        if (go_next)
        {
            int nx_m_id = -1;
            for (int i = 0; i < M - 1; ++i)
            {
                if (heavy_cost_order[i].first == pre_target)
                {
                    nx_m_id = heavy_cost_order[i + 1].first;
                    break;
                }
            }
            // 調べる機械がもうないとき、空データを送る。
            if (nx_m_id == -1)
            {
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
            else
            {
                // 次のマシーンの探索。
                double max_operation_rate = get_max_operation_rate(ep, nx_m_id);
                // 100%の時

                // 55..555を調べてもらう。
                pre_target = nx_m_id;
                pre_val = {5, 5};
                rep(i, M)
                {
                    if (i == nx_m_id)
                    {
                        ret_pat[i] = std::string(2 * X, '5');
                        ret_enc_pat[i] = {'5', '5'};
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                for (int i = 1; i <= 9; i++)
                    visited[i] = false;
                visited[5] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
    }
}

// startを全て9にするパターン //
Patterns binary_suggest2(int ep, auto &patterns, std::vector<EncodedPat> &basic_pat, EncodedPat &opt_pat, auto &scores, int &pre_target, pii &pre_val, std::vector<bool> &visited, pii &nx_step)
{
    std::vector<string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_enc_change;
    // std::cerr << pre_val.first << " " << pre_val.second << std::endl;
    assert(pre_val.first >= 1 and pre_val.first <= 9 and pre_val.second >= 1 and pre_val.second <= 9);
    if (ep == 1)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }
        pre_val = {9, 9};
        nx_step = {4, 4};
        return {ret_pat, ret_enc_pat, ret_enc_change};
    }
    else
    {

        bool go_next = false;
        if (scores[ep - 1] == 0)
        {
            analyze_failure(ep, patterns, basic_pat, opt_pat, scores, pre_target, pre_val);

            pii nx_val = {min(pre_val.first + nx_step.first, 9), min(pre_val.second + nx_step.second, 9)};
            if (nx_step.first == 0 and nx_step.second == 0)
            {
                go_next = true;
            }
            else
            {
                rep(i, M)
                {
                    if (i == pre_target)
                    {
                        char nxa = '0' + nx_val.first;
                        char nxb = '0' + nx_val.second;
                        ret_pat[i] = std::string(2 * X, nxa);
                        rep(x, X)
                        {
                            ret_pat[i][2 * x + 1] = nxb;
                        }
                        ret_enc_pat[i] = {nxa, nxb};
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        rep(x, X)
                        {
                            ret_pat[i][2 * x + 1] = opt_pat[i].second;
                        }
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                pre_val = nx_val;
                nx_step = {nx_step.first / 2, nx_step.second / 2};
                // visited[nx_val] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
        else
        {
            char opt_vala = '0' + pre_val.first;
            char opt_valb = '0' + pre_val.second;
            opt_pat[pre_target] = {opt_vala, opt_valb};
            pii nx_val = {max(pre_val.first - nx_step.first, 1), max(pre_val.second - nx_step.second, 1)};
            if (nx_step.first == 0 and nx_step.second == 0)
            {
                go_next = true;
            }
            else
            {
                rep(i, M)
                {
                    if (i == pre_target)
                    {
                        char nxa = '0' + nx_val.first;
                        char nxb = '0' + nx_val.second;
                        ret_pat[i] = std::string(2 * X, nxa);
                        if (nxa != nxb)
                        {
                            for (int x = 0; x < X; x++)
                            {
                                ret_pat[i][2 * x + 1] = nxb;
                            }
                        }
                        ret_enc_pat[i] = {nxa, nxb};
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        if (opt_pat[i].first != opt_pat[i].second)
                        {
                            for (int x = 0; x < X; x++)
                            {
                                ret_pat[i][2 * x + 1] = opt_pat[i].second;
                            }
                        }
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                pre_val = nx_val;
                nx_step = {nx_step.first / 2, nx_step.second / 2};
                // visited[nx_val] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
        if (go_next)
        {
            int nx_m_id = -1;
            for (int i = 0; i < M - 1; ++i)
            {
                if (heavy_cost_order[i].first == pre_target)
                {
                    nx_m_id = heavy_cost_order[i + 1].first;
                    break;
                }
            }
            // 調べる機械がもうないとき、空データを送る。
            if (nx_m_id == -1)
            {
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
            else
            {
                // 次のマシーンの探索。
                double max_operation_rate = get_max_operation_rate(ep, nx_m_id);
                // 100%近くの時
                if (max_operation_rate > 0.9)
                {
                    // 55..555を調べてもらう。
                    pre_target = nx_m_id;
                    pre_val = {5, 5};
                    nx_step = {2, 2};
                    rep(i, M)
                    {
                        if (i == nx_m_id)
                        {
                            ret_pat[i] = std::string(2 * X, '5');
                            ret_enc_pat[i] = {'5', '5'};
                        }
                        else
                        {
                            ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                            ret_enc_pat[i] = opt_pat[i];
                        }
                    }
                    for (int i = 1; i <= 9; i++)
                        visited[i] = false;
                    // visited[5]=true;
                    return {ret_pat, ret_enc_pat, ret_enc_change};
                }
                else
                {
                    // 稼働率に余裕がある時
                    // std::cerr << "FREEEEE!!!" << std::endl;

                    auto nx_pat = upper_bound(operation_time_ratio.begin(), operation_time_ratio.end(), max_operation_rate);

                    // pair<double, pair<char, char>> val = {max_operation_rate, {'.', ','}};
                    // auto nx_pat = upper_bound(operation_time_ratio_double.begin(), operation_time_ratio_double.end(), val);
                    // nx_pat++;

                    int val = nx_pat - operation_time_ratio.begin();
                    pre_target = nx_m_id;
                    pii nx_val = {val + 1, val + 1};
                    // pii nx_val = {(*nx_pat).second.first-'0', (*nx_pat).second.second-'0'};
                    std::cerr << "TIME RATIO DOUBLE FOUND: " << max_operation_rate << std::endl;
                    // std::cerr << (*nx_pat).first << " " << nx_val.first << " " << nx_val.second << std::endl;
                    pre_val = nx_val;
                    nx_step = {nx_val.first / 2, nx_val.second / 2};
                    rep(i, M)
                    {
                        if (i == nx_m_id)
                        {
                            ret_pat[i] = std::string(2 * X, '0' + nx_val.first);
                            for (int x = 0; x < X; ++x)
                            {
                                ret_pat[i][2 * x + 1] = '0' + nx_val.second;
                            }
                            ret_enc_pat[i] = {'0' + nx_val.first, '0' + nx_val.second};
                        }
                        else
                        {
                            ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                            rep(x, X)
                            {
                                ret_pat[i][2 * x + 1] = opt_pat[i].second;
                            }
                            ret_enc_pat[i] = opt_pat[i];
                        }
                    }
                    return {ret_pat, ret_enc_pat, ret_enc_change};
                }
            }
        }
    }
}

// startを全て9にするパターンA, Bを分けて最適化するパターン //
Patterns binary_suggest3(int ep, auto &patterns, std::vector<EncodedPat> &basic_pat, EncodedPat &opt_pat, auto &scores, std::pair<int, char> &pre_target, int &pre_val, int &nx_step)
{
    std::vector<std::string> ret_pat(M);
    EncodedPat ret_enc_pat(M);
    EncodedChange ret_enc_change;
    // std::cerr << pre_val.first << " " << pre_val.second << std::endl;
    std::cerr << "pre_target: (" << pre_target.first << ", " << pre_target.second << ")" << std::endl;
    assert(pre_val >= 1 and pre_val <= 9);
    if (ep == 1)
    {
        rep(i, M)
        {
            ret_pat[i] = string(2 * X, '9');
            ret_enc_pat[i] = {'9', '9'};
        }
        pre_val = 9;
        nx_step = 4;
        return {ret_pat, ret_enc_pat, ret_enc_change};
    }
    else
    {

        bool go_next = false;
        if (scores[ep - 1] == 0)
        {
            // analyze_failure(ep, patterns, basic_pat, opt_pat, scores, pre_target, pre_val);

            int nx_val = min(pre_val + nx_step, 9);
            if (nx_step==0)
            {
                go_next = true;
            }
            else
            {
                rep(i, M)
                {
                    if (i == pre_target.first)
                    {
                        char nx = '0' + nx_val;
                        ret_pat[i] = std::string(2 * X, nx);
                        rep(x, X)
                        {
                            if(pre_target.second=='A') ret_pat[i][2 * x+1] = opt_pat[i].second;
                            else if(pre_target.second=='B') ret_pat[i][2 * x] = opt_pat[i].first;
                        }
                        if(pre_target.second=='A') ret_enc_pat[i] = {nx, opt_pat[i].second};
                        else ret_enc_pat[i]={opt_pat[i].first, nx};
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        rep(x, X)
                        {
                            ret_pat[i][2 * x + 1] = opt_pat[i].second;
                        }
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                pre_val = nx_val;
                nx_step = nx_step / 2;
                // visited[nx_val] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
        else
        {
            char opt_val = '0' + pre_val;
            if(pre_target.second=='A') opt_pat[pre_target.first].first = opt_val;
            else if(pre_target.second=='B') opt_pat[pre_target.first].second = opt_val;
            int nx_val = max(pre_val - nx_step, 1);
            if (nx_step == 0)
            {
                go_next = true;
            }
            else
            {
                rep(i, M)
                {
                    if (i == pre_target.first)
                    {
                        char nx = '0' + nx_val;
                        int m_id=pre_target.first;
                        ret_pat[i] = std::string(2 * X, nx);
                        rep(x,X) {
                            if(pre_target.second=='A')ret_pat[i][2*x+1] = opt_pat[m_id].second;
                            else ret_pat[i][2*x]=opt_pat[m_id].first;
                        }
                        //ターゲット以外のもの
                        ret_enc_pat[i]=opt_pat[i];
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        if (opt_pat[i].first != opt_pat[i].second)
                        {
                            for (int x = 0; x < X; x++)
                            {
                                ret_pat[i][2 * x + 1] = opt_pat[i].second;
                            }
                        }
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                pre_val = nx_val;
                nx_step = nx_step / 2;
                // visited[nx_val] = true;
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
        }
        if (go_next)
        {
            int nx_m_id = -1;
            char nx_ab;
            for (int i = 0; i < (int)heavy_cost_order_AB.size(); ++i)
            {
                if (heavy_cost_order_AB[i].first.first == pre_target.first and heavy_cost_order_AB[i].first.second == pre_target.second)
                {
                    nx_m_id = heavy_cost_order_AB[i + 1].first.first;
                    nx_ab=heavy_cost_order_AB[i+1].first.second;
                    break;
                }
            }
            // 調べる機械がもうないとき、空データを送る。
            if (nx_m_id == -1)
            {
                return {ret_pat, ret_enc_pat, ret_enc_change};
            }
            else
            {
                // 次のマシーンの探索。
                
                // 55..555を調べてもらう。
                pre_target = {nx_m_id, nx_ab};
                pre_val = 5;
                nx_step = 2;
                rep(i, M)
                {
                    if (i == nx_m_id)
                    {
                        ret_pat[i] = std::string(2 * X, '5');
                        for(int x = 0; x < X; x++) {
                            if(nx_ab != 'A') ret_pat[i][x*2] = opt_pat[i].first;
                            else ret_pat[i][x*2+1] = opt_pat[i].second;
                        }
                        ret_enc_pat[i] = {'5', '5'};
                        if(nx_ab!='A') ret_enc_pat[i].first=opt_pat[i].first;
                        else ret_enc_pat[i].second=opt_pat[i].second;
                    }
                    else
                    {
                        ret_pat[i] = std::string(2 * X, opt_pat[i].first);
                        rep(x,X) {
                            ret_pat[i][2*x+1]=opt_pat[i].second;
                        }
                        ret_enc_pat[i] = opt_pat[i];
                    }
                }
                return {ret_pat, ret_enc_pat, ret_enc_change};
                
                
            }
        }
    }
}

// binary_suggest2で得られる最終状態の稼働パターンは平日休日の両方を下げることができない。
// コストの高いものからの片方ずつ下げていく、もしくは途中で稼働パターンを下げることにより小さいものが見つかる。DFS
// 両方下げるということはもうできない。他が9999..999のパターンの時でさえうまくいっていなかったDFS
// Patterns fine_opt_suggest(){

// }

Patterns random_suggest(int ep, auto &patterns, const std::vector<EncodedPat> &basic_pat, const std::vector<EncodedChange> &change_pat, const std::vector<ll> &scores, const std::vector<double> &eval_scores, const std::vector<ll> &exp_scores, auto &V, auto &D)
{
    clock_t suggest_start_time = clock();
    vector<vector<string>> x_points(patterns.cbegin(), patterns.cbegin() + ep);
    dvector y_points(ep);
    rep(i, ep)
    {
        if (USE_LOG)
            y_points[i] = log10(1 + eval_scores[i]);
        else if (USE_QUADRATIC)
            y_points[i] = eval_scores[i] * eval_scores[i];
        else
            y_points[i] = eval_scores[i];
    }

    // 評価スコアの累積和とルーレット
    vector<double> score_table(ep + 1, 0.);
    rep(e, ep) score_table[e + 1] = score_table[e] + y_points[e];
    uniform_real_distribution<> roulette(0, score_table[ep]);

    // > --------------------------- search loop -------------------------------<
    cur_time = clock();
    vector<string> next_pattern;
    EncodedPat best_enc_pat;
    EncodedChange best_change;
    Patterns next_candidate;

    while (1)
    {
        std::cerr << "search loop" << std::endl;
        double val = roulette(mt);
        int idx = lower_bound(score_table.begin(), score_table.end(), val) - score_table.begin();
        idx--; //最初に0を入れてるのでindexを１つ下げる。
        EncodedPat b_pat = basic_pat[idx];
        EncodedChange c_pat = change_pat[idx];

        next_candidate = modify_candidate(idx, b_pat, c_pat, V, D); // TODO: 次の候補点の選び方が大切。ここの実装が大事。重点的に取り組む
        // for(auto np: get<0>(next_candidate)) std::cerr << np << "\n";
        if (pattern_set.find(get<0>(next_candidate)) != pattern_set.end())
            continue;
        else
            break;
    }
    return next_candidate;
}

void update_heavy_cost_order(auto &opt_pat)
{
    heavy_cost_order.clear();
    for (int i = 0; i < M; i++)
    {
        int aid = opt_pat[i].first - '0';
        int bid = opt_pat[i].second - '0';
        ll cost = 5 * productLine[i].costA[aid] + 2 * productLine[i].costB[bid];
        heavy_cost_order.push_back({i, cost});
    }
    sort(heavy_cost_order.begin(), heavy_cost_order.end(), [](std::pair<int, ll> &a, std::pair<int, ll> &b)
         { return a.second > b.second; });
}

Patterns fine_opt_suggest(int ep, char &before_target, int &start_id, std::vector<ll> &scores, EncodedPat &opt_pat)
{
    std::vector<std::string> next_pattern(M);
    EncodedPat next_enc_pat(M);
    next_enc_pat = opt_pat;
    EncodedChange next_change;
    heavy_cost_order[1].first;

    if (start_id == M)
    {
        return {next_pattern, next_enc_pat, next_change};
    }
    if (before_target == 'Z')
    {
        next_pattern = best_basic_opt_pattern;

        // start_id=heavy_cost_order[0].first;
        int m_id = heavy_cost_order[start_id].first;
        cerr << "m_id: " << m_id << endl;
        char nx = next_enc_pat[m_id].first - 1;
        cerr << "nx: " << nx << endl;
        if (nx == '0')
            nx = '1';
        rep(i, X) next_pattern[m_id][2 * i] = nx;
        next_enc_pat[m_id].first = nx;
        before_target = 'A';
        return {next_pattern, next_enc_pat, next_change};
    }
    else
    {
        if (before_target == 'A')
        {
            before_target = 'B';
            if (scores[ep - 1] > 0)
            {
                return fine_opt_suggest(ep, before_target, ++start_id, scores, opt_pat);
            }
            else
            {
                next_pattern = best_eval_pattern;

                int m_id = heavy_cost_order[start_id].first;
                cerr << m_id << endl;
                // start_id=heavy_cost_order[0].first;
                char nx = next_enc_pat[m_id].second - 1;
                if (nx == '0')
                    nx = '1';
                cerr << nx << endl;
                next_enc_pat[m_id].second = nx;
                for (int m = 0; m < M; m++)
                {
                    for (int x = 0; x < X; x++)
                    {
                        if (m == m_id)
                        {
                            next_pattern[m][2 * x] = opt_pat[m].first;
                            next_pattern[m][2 * x + 1] = nx;
                        }
                        else
                        {
                            next_pattern[m][2 * x] = opt_pat[m].first;
                            next_pattern[m][2 * x + 1] = opt_pat[m].second;
                        }
                    }
                }
                next_enc_pat[m_id].second = nx;
                start_id++;
                return {next_pattern, next_enc_pat, next_change};
            }
        }
        else if (before_target == 'B')
        {
            if (start_id == M)
            {
                cerr << "get" << endl;
                return (Patterns){next_pattern, next_enc_pat, next_change};
            }
            next_pattern = best_eval_pattern;
            before_target = 'A';
            int m_id = heavy_cost_order[start_id].first;
            // start_id=heavy_cost_order[0].first;
            char nx = next_enc_pat[m_id].first - 1;
            if (nx == '0')
                nx = '1';
            // rep(i, X) next_pattern[m_id][2*i] = nx;
            next_enc_pat[start_id].first = nx;
            for (int m = 0; m < M; m++)
            {
                for (int x = 0; x < X; x++)
                {
                    if (m == m_id)
                    {
                        next_pattern[m][2 * x] = nx;
                        next_pattern[m][2 * x + 1] = opt_pat[m].second;
                    }
                    else
                    {
                        next_pattern[m][2 * x] = opt_pat[m].first;
                        next_pattern[m][2 * x + 1] = opt_pat[m].second;
                    }
                }
            }
            next_enc_pat[m_id].first = nx;
            return {next_pattern, next_enc_pat, next_change};
        }
    }
}

int main()
{
    start_time = clock();
    input(X, M, C, E, productLine);
    ll best_score = -1;
    ll score = -1;
    string ini_pat(2 * X, '1');
    std::vector<std::vector<std::string>> patterns(E + 1, std::vector<std::string>(M, ini_pat)); // 稼働パターン

    // std::vector<std::string> best_pattern;
    std::vector<EncodedPat> basic_patterns(E + 1, EncodedPat(M, {'9', '9'}));
    std::vector<EncodedChange> change_patterns(E + 1); // あるepochの稼働パターンの時の稼動変更位置(m, x)

    std::vector<ll> scores(E + 1, 0);           // 計画プログラムから帰ってくるスコア
    std::vector<double> eval_scores(E + 1, 0.); // 評価に用いるスコア(ベイズ最適化とか) 納期が遅れまくりだと小さくなる。
    std::vector<ll> expected_scores(E + 1, 0);  // 納期が全て間に合うと仮定した時の期待されるスコア
    std::vector<bool> visited(10, false);
    EncodedPat optimized_pat(M, {'9', '9'});
    std::pair<int, int> pre_val = {9, 9};
    std::pair<int, int> nx_step = {4, 4};
    int pre_val_ab = 9;
    int nx_step_ab=4;

    bool basic_opt = false;
    bool fine_opt = false;
    bool fine_start = false;
    bool change_opt = false;
    bool greedy_opt = false;
    bool greedy_opt2 = false;

    char before_target = 'Z';
    int heavy_id = 0;
    std::cerr << "***************************************** START **************************************\n";
    // Init
    init(scores, eval_scores, expected_scores);
    int pre_target = heavy_cost_order[0].first;
    std::pair<int,char> pre_target_ab = heavy_cost_order_AB[0].first;
    // >-------------------------- reaction loop ----------------------------<
    for (int epoch = 1; epoch <= E; epoch++)
    {
        std::cerr << "########## EPOCH: " << setfill('0') << setw(3) << epoch << " ##########\n";
        cur_time = clock();
        std::cerr << "ELAPSED: " << (double)(cur_time - start_time) / CLOCKS_PER_SEC << "\n";

        // std::vector<Patterns> next_patterns = suggest2(epoch, patterns, basic_patterns, change_patterns, scores, eval_scores, expected_scores, V, D);
        // Patterns next_pattern = binary_suggest(epoch, patterns, basic_patterns, optimized_pat, scores, pre_target, pre_val, visited);

        Patterns next_pattern;
        std::vector<Patterns> next_patterns;
        bool null_ans = true;
        if (!basic_opt)
        {
            // next_pattern = binary_suggest2(epoch, patterns, basic_patterns, optimized_pat, scores, pre_target, pre_val, visited, nx_step);
            next_pattern=binary_suggest3(epoch, patterns, basic_patterns, optimized_pat, scores, pre_target_ab, pre_val_ab, nx_step_ab);
            next_patterns = {next_pattern};
            if (std::get<0>(next_pattern)[0].size() != 2 * X)
            {
                std::cerr << std::get<0>(next_pattern)[0] << std::endl;
                basic_opt = true;
                best_basic_opt_score = best_score;
                best_basic_opt_pattern = best_score_pattern;
            }
        }
        if (basic_opt and !greedy_opt) //一時的にnullansを外す
        {
            next_pattern = greedy_change(patterns, optimized_pat);
            next_patterns = {next_pattern};
            null_ans = false;
            greedy_opt = true;
            update_heavy_cost_order(optimized_pat);
        }
        if (greedy_opt and !fine_opt and null_ans) //一時的にnullansを
        {
            cerr << heavy_id << " " << before_target << endl;
            next_pattern = fine_opt_suggest(epoch, before_target, heavy_id, scores, optimized_pat);
            cerr << std::get<0>(next_pattern)[0].size() << endl;
            cerr << std::get<0>(next_pattern)[0] << endl;
            next_patterns = {next_pattern};
            if (std::get<0>(next_pattern)[0].size() != 2 * X)
            {
                fine_opt = true;
                cerr << "GO" << endl;
            }
            fine_start = true;
        }
        if (fine_opt and !greedy_opt2)
        {
            next_pattern = greedy_change(patterns, optimized_pat);
            next_patterns = {next_pattern};
            null_ans = false;
            greedy_opt2 = true;
            update_heavy_cost_order(optimized_pat);
        }
        if (greedy_opt2 and null_ans)
        {
            update_heavy_cost_order(optimized_pat);
            next_pattern = random_suggest(epoch, patterns, basic_patterns, change_patterns, scores, eval_scores, expected_scores, V, D);
            // std::vector<std::string> tmp = best_basic_opt_pattern;
            // int heavy_id = heavy_cost_order[0].first;
            // char now = best_basic_opt_pattern[heavy_id][0];
            // rep(i, X - 1) tmp[heavy_id][2 * i + 1 + 2] = now - 1;
            // get<0>(next_pattern) = tmp;
            next_patterns.push_back(next_pattern);
            
            // next_pattern = ;
        }
        for (int iter = 0; iter < (int)next_patterns.size(); iter++)
        {

            if (iter + epoch > E)
                break;
            Patterns np = next_patterns[iter];
            std::vector<std::string> nx_pat = std::get<0>(np);
            EncodedPat nx_bs_pat = std::get<1>(np);
            EncodedChange nx_ch_pat = std::get<2>(np);
            assert(nx_pat[0].size() == 2 * X);
            if (nx_pat[0].size() != 2 * X)
            {
                nx_pat = patterns[epoch + iter - 1];
                nx_bs_pat = basic_patterns[epoch + iter - 1];
                nx_ch_pat = change_patterns[epoch + iter - 1];
            }
            else
            {
                pattern_set.insert(nx_pat);
                rep(m, M)
                {
                    patterns[epoch + iter][m] = nx_pat[m];
                    basic_patterns[epoch + iter][m] = nx_bs_pat[m];
                    change_patterns[epoch + iter] = nx_ch_pat;
                }
            }
            calc_expected_score(epoch + iter, patterns, expected_scores);
            show_result(epoch + iter, patterns);
            get_result(epoch + iter, scores, V, D, eval_scores, expected_scores);

            // update
            if (eval_scores[epoch + iter] > best_eval_score)
            {
                best_eval_epoch = epoch + iter;
                best_eval_score = eval_scores[epoch + iter];
                best_eval_pattern = patterns[epoch + iter];
            }
            if (scores[epoch + iter] > best_score)
            {
                best_score_epoch = epoch + iter;
                best_score = scores[epoch + iter];
                best_score_pattern = patterns[epoch + iter];
            }
            if (fine_start and scores[epoch + iter] > best_basic_opt_score)
            {
                cerr << scores[epoch + iter] << " " << best_basic_opt_score << endl;
                best_basic_opt_score = scores[epoch + iter];
                rep(j, M) optimized_pat[j] = basic_patterns[epoch + iter][j];
            }

#ifdef SHOW_EACH_PATTERN
            for (auto np : std::get<0>(np))
            {
                cerr << np << std::endl;
            }
#endif

            std::cerr << "SCORE: " << scores[epoch + iter] << ", EVAL: " << eval_scores[epoch + iter] << ", EXP: " << expected_scores[epoch + iter] << "\n";
        }
        epoch += (next_patterns.size() - 1);
    }
    // >-------------------------- reaction loop ----------------------------<

#ifdef DEBUG
    for (int x = 0; x < X; x++)
    {
        std::cerr << "WEEK: " << x + 1 << std::endl;
        for (int m = 0; m < M; m++)
        {
            std::cerr << "MACHINE: " << m << ", RD: " << RD[1][x][m] << ", RO: " << RO[1][x][m] << std::endl;
        }
    }
    std::cerr << "cost list" << std::endl;
    for (auto p : heavy_cost_order)
    {
        std::cerr << "(machine, cost): " << p.first + 1 << " " << p.second << std::endl;
    }

    std::cerr << "cost" << std::endl;
    for (int i = 0; i < M; i++)
    {
        std::cerr << "### MACHINE: " << i + 1 << std::endl;
        for (int j = 1; j <= 9; j++)
        {
            std::cerr << "A,B: " << productLine[i].costA[j] << " " << productLine[i].costB[j] << std::endl;
        }
    }

    for (auto s : best_score_pattern)
    {
        std::cerr << s << std::endl;
    }
#endif
    return 0;
}
#endif