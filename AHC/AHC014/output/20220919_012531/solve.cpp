#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include <stack>
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

// #include <boost/math/distributions.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/vector.hpp>
// #include <boost/numeric/ublas/lu.hpp>
// #include <boost/numeric/ublas/triangular.hpp>

// #define DEBUG
// #define TEST


#define rep(i, n) for (int i = 0; i < (n); i++)
#define ll long long
#define pii pair<int, int>
// #define dmatrix boost::numeric::ublas::matrix<double>
// #define dvector boost::numeric::ublas::vector<double>

using namespace std;
// namespace bmath = boost::math;
// namespace ublas = boost::numeric::ublas;

const int dx[4] = {1, 0, -1, 0};
const int dy[4] = {0, 1, 0, -1};
const int dx8[8] = {-1, -1, -1, 0, 1, 1, 1, 0}; //clockwise
const int dy8[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
map<pair<int,int>, int> direction_map = {
        {{-1, -1}, 0},
        {{-1,  0}, 1},
        {{-1,  1}, 2},
        {{ 0,  1}, 3},
        {{ 1,  1}, 4},
        {{ 1,  0}, 5},
        {{ 1, -1}, 6},
        {{ 0, -1}, 7},
    };
const int INF = 1 << 30;
const ll LINF = 1LL << 60;
const double epsilon = 1e-6;

// global variable//////
struct Vector2 {
    int x;
    int y;
    Vector2(){}
    Vector2(int _x, int _y): x(_x), y(_y){}
};

struct Point {
    int x, y;
    bool exist = false;
    std::vector<bool> used;
    Point(){
        used.resize(8, false);
    }
    Point(int _x, int _y): x(_x), y(_y) {
        used.resize(8, false);
    }
    Point(int _x, int _y, bool _exist): x(_x), y(_y), exist(_exist){
        used.resize(8, false);
    }
};

struct Parameters {
    double TL = 5.00;
    double TIME_LIMIT_90 = 4.50;
    double DECAY;
    std::string file_name="../output/config.txt";

    int N, M;
    std::vector<Point> init_points;

    Parameters(){}
    Parameters(int n, int m, std::vector<Point>& _init_points): N(n), M(m) {
        init_points.resize(m);
        init_points = _init_points;
    }

    void save() {
        ofstream ofs;
        ofs.open(file_name, ios::out);
        cout << "TL: " << TL << "\n";
        cout << "DECAY: " << DECAY << "\n";
        ofs.close();
        return;
    } 
};



struct Rectangle {
    std::vector<Point> vert;
    Rectangle(){
        vert.resize(4);
    }
    Rectangle(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4){
        vert.resize(4);
        vert[0] = Point{x1, y1};
        vert[1] = Point{x2, y2};
        vert[2] = Point{x3, y3};
        vert[3] = Point{x4, y4};
    };
    Rectangle(Point a, Point b, Point c, Point d){
        vert.resize(4);
        vert[0] = a;
        vert[1] = b;
        vert[2] = c;
        vert[3] = d;
        for(int i = 0; i < 4; i++) {
            assert(vert[i].x == vert[(i+1)%4].x || vert[i].y == vert[(i+1)%4].y);
        }
    };
    Rectangle(Point a[]) {
        vert.resize(4);
        for(int i = 0; i < 4; i++) vert[i] = a[i];
    };
    Rectangle(std::vector<Point> &a) {
        vert.resize(4);
        for(int i = 0; i < 4; i++) vert[i] = a[i];
    }
};

clock_t start_time;
clock_t cur_time;
////////////////////////

random_device seed_gen;
mt19937 mt(seed_gen());
// bmath::normal_distribution<> distribution(0, 1);

bool out_bound(int nx, int ny, int n)
{
    return (nx < 0 || ny < 0 || nx >= n || ny >= n) ? true : false;
}

template<class T> T pow2(T x) {
    return x * x;
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

void test(){
    std::cerr << "### TEST ###\n";
    return;
}

// # memo
// 左下が(0, 0), 右上が(N-1, N-1)
// 外縁に近い方がスコアが良い。
// 内部の探索をするメリットがほとんどない。<- 状況によってはアリかもしれない。
// 作れる格子点を洗い出すのも一つの手
// 初期状態で内側の四角形領域に格子点が埋まっているがここから出ようと思うと45度傾いた長方形を作る必要がある。
// この四角形を最初に探索してみるのもありかもしれない。
// 再帰

// # 仮説
// 実験してみた結果、角に格子点を設けるのは理論的に無理。諦める。角から1, 2マス空いた場所は狙おうと思えば狙える。
//



class State {
    public:
        int n;
        int m;
        int center;
        ll S = 0;
        std::vector<std::vector<Point>> points;
        std::vector<Point> stack_points;
        std::vector<Rectangle> rectangles;
        State(){}
        State(int _n, int _m, const Parameters &param): n(_n), m(_m) {
            center = (_n-1)/2;
            points.resize(n);
            for(std::vector<Point> &v: points) {
                v.resize(n);
            }
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    points[i][j].x = j;
                    points[i][j].y = i;
                    S += pow2(i-center) + pow2(j-center) + 1;
                }
            }
            for(Point _p: param.init_points) {
                points[_p.y][_p.x] = _p;
                stack_points.push_back(_p);
            }
        }

        std::vector<Rectangle> search_rectangles(int x, int y) {
            std::vector<Rectangle> ret;

            // clockwise-search
            for(int dir = 0; dir < 8; dir ++ ) {
                Rectangle candidate;
                candidate.vert[0] = Point{x, y};
                int cx = x;
                int cy = y;
                for(int n_rot = 0; n_rot < 4; n_rot++) {
                    Point nx_point;
                    int rec_dir = (8 + dir + 2 * n_rot) % 8;
                    if(go_till_find(cx, cy, rec_dir, nx_point, x, y)) {
                        cx = nx_point.x;
                        cy = nx_point.y;
                        // std::cerr << "FOUND AT (" << nx_point.x << "," << nx_point.y << "), from[" << x << "," << y << "]" << std::endl;
                        if(n_rot != 3) {
                            candidate.vert[n_rot + 1] = nx_point;
                        } else if(n_rot == 3 and nx_point.x == x and nx_point.y == y){
                            ret.push_back(candidate);
                        }
                    } else {
                        break;
                    }
                }
            }

            // counter-clockwise-search
            // for(int dir = 0; dir < 8; dir ++ ) {
            //     Rectangle candidate;
            //     candidate.vert[0] = Point{x, y};
            //     int cx = x;
            //     int cy = y;
            //     for(int n_rot = 0; n_rot < 4; n_rot++) {
            //         Point nx_point;
            //         int rec_dir = (8 + dir + 2 * n_rot) % 8;
            //         if(go_till_find(cx, cy, rec_dir, nx_point, x, y)) {
            //             cx = nx_point.x;
            //             cy = nx_point.y;
            //             // std::cerr << "FOUND AT (" << nx_point.x << "," << nx_point.y << "), from[" << x << "," << y << "]" << std::endl;
            //             if(n_rot != 3) {
            //                 candidate.vert[n_rot + 1] = nx_point;
            //             } else if(n_rot == 3 and nx_point.x == x and nx_point.y == y){
            //                 ret.push_back(candidate);
            //             }
            //         } else {
            //             break;
            //         }
            //     }
            // }

            // for(int dir = 0; dir < 8; dir ++ ) {
            //     Rectangle candidate;
            //     candidate.vert[0] = Point{x, y};
            //     for(int n_rot = 0; n_rot < 4; n_rot++) {
            //         Point nx_point;
            //         int rec_dir = (8 + dir + 2 * n_rot) % 8; // this line is different from above codes.
            //         if(go_till_find(x, y, rec_dir, nx_point)) {
            //             if(n_rot != 3) {
            //                 candidate.vert[n_rot + 1] = nx_point;
            //             } else if(n_rot == 3 and nx_point.x == x and nx_point.y == y){
            //                 ret.push_back(candidate);
            //             }
            //         } else {
            //             break;
            //         }
            //     }
            // }
            return ret;
        }

        bool go_till_find(int x, int y, int dir, Point& nx_point, int sx, int sy) {
            bool found = false;
            int cx = x;
            int cy = y;
            while(!found) {
                if(points[cy][cx].used[dir] || points[cy][cx].used[(dir+4)%4]) {
                    // if the current position has an edge of rectangl along the direction, stop the search.
                    // std::cerr << "INVALID!!!!" << std::endl;
                    return false;
                }
                int nx = cx + dx8[dir];
                int ny = cy + dy8[dir];
                if(out_bound(nx, ny, n)) {
                    return false;
                }
                if(points[ny][nx].exist or (nx == sx and ny == sy)) {
                    nx_point = points[ny][nx];
                    return true;
                }
                cx = nx;
                cy = ny;
            }
            return found ? true : false;
        }

        void add_rectangle(Rectangle rect) {
            rectangles.push_back(rect);
            this->_update_points_on_rectangles(rect);
        }
        
    private:
        void _update_points_on_rectangles(Rectangle &rect) {
            for(int i = 0; i < 4; i++) {
                Point from = rect.vert[i];
                Point to = rect.vert[(i+1)%4];
                int dx = to.x-from.x;
                int dy = to.y-from.y;
                // std::cerr << "From:(" << from.x << "," << from.y <<"), TO: (" << to.x << "," << to.y << ")" << std::endl;
                int d = max(abs(dx), abs(dy));
                int ex = dx/d;
                int ey = dy/d;
                int dir = direction_map[{ex, ey}];
                // std::cerr << "direction: " << dir << " " << ex << " " << ey << std::endl;
                for(int j = 0; j <= d; j++) {
                    int cx = from.x + j * dx8[dir];
                    int cy = from.y + j * dy8[dir];
                    if(j < d){
                        points[cy][cx].used[dir] = true;
                    }
                    if(j > 0) {
                        points[cy][cx].used[(dir+4)%4] = true;
                    }
                }
            }
        }
};

struct Result {
    int K;
    std::vector<Rectangle> rects;
    std::vector<Point> points;
    long long score;
    Result(){}
    Result(int k, std::vector<Rectangle> r): K(k), rects(r) {}
    Result(int k, std::vector<Rectangle> r, std::vector<Point> p): K(k), rects(r), points(p) {}
    Result(int k, std::vector<Rectangle> r, long long s): K(k), rects(r), score(s){}  
    Result(int k, std::vector<Rectangle> r, std::vector<Point> p, long long s): K(k), rects(r), points(p), score(s){}  
};


class Solver {
    public:
        Parameters params;
        State state;
        std::vector<Result> result_list;

        Solver() {
            std::cerr << "## ENTERING INITIALIZATION " << std::endl;
            int n, m; cin >> n >> m;
            std::vector<Point> init_points(m);
            for(int i = 0; i < m; i++) {
                std::cin >> init_points[i].x >> init_points[i].y;
                init_points[i].exist = true;
            }
            params = Parameters(n, m, init_points);
            state = State(n, m, params);
            std::cerr << "N: " << n << std::endl;
            std::cerr << "M: " << m <<std::endl;

            std::cerr << "## FINISHED INIT ##" << std::endl;
        }


        void show_result() {
            ll score_max=0;
            Result best;
            for(Result res: result_list) {
                if(res.score > score_max) {
                    score_max = res.score;
                    best = res;
                }
            }
            int K = best.K;
            cout << K << endl;
            for(Rectangle r: best.rects) {
                for(int j = 0; j < 4; j ++ ) {
                    cout << r.vert[j].x << " " << r.vert[j].y << " ";
                }
                cout << "\n";
            }
            return ;
        }

        double calc_score() {
            double ret;
            long long w = 0;
            for(Point p: state.stack_points) {
                w += pow2(p.x-state.center) + pow2(p.y-state.center) + 1;
            }
            ret = round((double)1e6 * pow2(params.N) * w / params.M / state.S);
            return ret;
        }

        double calc_result_score(int idx) {
            double ret;
            long long w = 0;
            for(Point p: result_list[idx].points) {
                w += pow2(p.x-state.center) + pow2(p.y-state.center) + 1;
            }
            ret = round((double)1e6 * pow2(params.N) * w / params.M / state.S);
            return ret;
        }

        // greedy solution
        void baseline() {
            std::cerr << "## ENTERING baseline codes" << std::endl;
            cur_time = clock();
            bool time_limit_flag = false;

            while(!time_limit_flag) {
                for(int x = 0; x < params.N and !time_limit_flag; x++) {
                    for(int y = 0; y < params.N and !time_limit_flag; y++) {
                        if(state.points[y][x].exist) continue;
                        std::vector<Rectangle> rec_candidates = state.search_rectangles(x, y);
                        // std::err << rec_candidates.size() << std::endl;
                        if(rec_candidates.size() > 0) {
                            int cand_sz = rec_candidates.size();
                            uniform_int_distribution<> rand(0, cand_sz - 1);
                            int rec_id = rand(mt);
                            // std::cerr << "rectangle_id: " << rec_id << std::endl;
                            Rectangle nx_rec = rec_candidates[rec_id];
                            state.points[y][x].exist = true;
                            state.add_rectangle(nx_rec);
                            state.stack_points.push_back(Point{x, y});
                            // state.rectangles.push_back(nx_rec);
                        }
                        time_limit_flag = this->time_flag_check();
                        // std::cerr << (double)cur_time / CLOCKS_PER_SEC<< std::endl;
                    }
                }
            }
            long long score = this->calc_score();
            result_list.push_back(Result{(int)state.rectangles.size(), state.rectangles, score});
        }

        //start search from center point. if found point, resart from center
        void baseline_ver2() {
            std::cerr << "## ENTERING baseline codes" << std::endl;
            cur_time = clock();
            bool time_limit_flag = false;
            std::vector<std::vector<bool>> visited(params.N, std::vector<bool>(params.N, false));
            std::vector<Vector2> search_order(pow2(params.N));
            int cur_x = state.center, cur_y = state.center;
            int cur_dir = 0;
            for(int i = 0; i < (int)search_order.size(); i++) {
                search_order[i] = Vector2{cur_x, cur_y};
                visited[cur_y][cur_x] = true;
                int next_dir = (cur_dir + 1)%4;
                int next_x = cur_x + dx[next_dir];
                int next_y = cur_y + dy[next_dir];
                if(visited[next_y][next_x]) {
                    cur_x = cur_x + dx[cur_dir];
                    cur_y = cur_y + dy[cur_dir];
                } else {
                    cur_x = next_x;
                    cur_y = next_y;
                    cur_dir = next_dir;
                }
            }
            for(int i = 0; i < (int)search_order.size(); i++) {
                std::cerr << search_order[i].x << " " << search_order[i].y << std::endl;
            }
            std::cerr << search_order.size() << std::endl;

            bool restart_flag = false;
            while(!time_limit_flag) {
                restart_flag = false;
                for(int i = 0; i < (int)search_order.size() and !time_limit_flag and !restart_flag; i++) {
                    int x = search_order[i].x;
                    int y = search_order[i].y;
                    if(state.points[y][x].exist) continue;
                    std::vector<Rectangle> rec_candidates = state.search_rectangles(x, y);
                    // std::err << rec_candidates.size() << std::endl;
                    if(rec_candidates.size() > 0) {
                        int cand_sz = rec_candidates.size();
                        uniform_int_distribution<> rand(0, cand_sz - 1);
                        int rec_id = rand(mt);
                        // std::cerr << "rectangle_id: " << rec_id << std::endl;
                        Rectangle nx_rec = rec_candidates[rec_id];
                        state.points[y][x].exist = true;
                        state.add_rectangle(nx_rec);
                        restart_flag=true;
                        state.stack_points.push_back(Point{x, y});
                        // state.rectangles.push_back(nx_rec);
                    }
                    time_limit_flag = this->time_flag_check();
                    // std::cerr << (double)cur_time / CLOCKS_PER_SEC<< std::endl;
                    
                }
            }
        }

        void run() {
            // >--------------------------- codes ---------------------------<
            this->baseline();
            // >-------------------------------------------------------------< 
            this->show_result();
            this->summary();
            #ifdef DEBUG
            this->debug();
            #endif
        }

        bool time_flag_check() {
            cur_time = clock();
            return ((double)cur_time - start_time) / CLOCKS_PER_SEC < this->params.TIME_LIMIT_90 ? false : true;
        }

        void summary() {
            std::cerr << "###################### SUMMARY ######################" << std::endl;
            std::cerr << "SCORE: " << this->calc_score() << std::endl;

            std::cerr << "#####################################################" << std::endl;

        }

        void debug() {
            for(int x = 0; x < params.N; x++) {
                for(int y = 0; y < params.N; y++) {
                    std::cerr << "(LOCATION):" << "(" << x << "," << y << ")" << std::endl;
                    Point p = state.points[y][x];
                    std::cerr << "(LINE) ";
                    for(int i = 0; i < 8; i++) {
                        std::cerr << p.used[i] << " ";
                    }
                    std::cerr << std:: endl;
                    std::cerr << "(EXIST)" << p.exist << std::endl;
                    std::cerr << std::endl;
                }
            }
        }
};


int main() {
    start_time = clock();
    Solver solver;
    solver.run();
    return 0;
}