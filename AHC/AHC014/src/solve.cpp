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
#include <numeric>
#include <bitset>
#include <string>
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

#define LOCAL
#ifndef _DEBUG
#define _DEBUG 0
#endif
#ifndef _TEST
#define _TEST 0
#endif
#ifndef _OUTPUT_FILE
#define _OUTPUT_FILE "../output/config.txt"
#endif


#define rep(i, n) for (int i = 0; i < (n); i++)
#define ll long long
#define pii pair<int, int>
#define Grid std::vector<std::vector<Point>>
// #define dmatrix boost::numeric::ublas::matrix<double>
// #define dvector boost::numeric::ublas::vector<double>

using namespace std;
// namespace bmath = boost::math;
// namespace ublas = boost::numeric::ublas;

const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, 1, 0, -1};
const int dcx[4] = {-1, -1, 1, 1}; //cross
const int dcy[4] = {-1, 1, 1, -1};
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

struct Vector2 {
    int x;
    int y;
    Vector2(){};
    Vector2(int _x, int _y): x(_x), y(_y){};
    Vector2 operator+(const Vector2& rhs);
    Vector2 operator-(const Vector2& rhs);
    int get_direction() {
        int m = max(x, y);
        int ex = x / m;
        int ey = y / m;
        return direction_map[{ex, ey}];
    }
};

Vector2 Vector2::operator+(const Vector2& rhs) {
    return Vector2{this->x + rhs.x, this->y + rhs.y};
}

Vector2 Vector2::operator-(const Vector2& rhs) {
    return Vector2{this->x - rhs.x, this->y - rhs.y};
}


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
    void show() {
        std::cerr << "(" << x << "," << y << ")" << std::endl;
    }
};

struct Parameters {
    double TL = 5.00;
    double TIME_LIMIT_90 = 4.80;
    double DECAY;
    std::string file_name = _OUTPUT_FILE;

    int N;
    int M;
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
    };
    Rectangle(Point a[]) {
        vert.resize(4);
        for(int i = 0; i < 4; i++) vert[i] = a[i];
    };
    Rectangle(std::vector<Point> &a) {
        vert.resize(4);
        for(int i = 0; i < 4; i++) vert[i] = a[i];
    }

    void show() {
        for(int i = 0; i < 4; i++) {
            std::cerr << i << ":(" << vert[i].x << "," << vert[i].y << "), ";
        }
        std::cerr << std::endl;
    }
};

struct RectangleTree {
    Rectangle rect;
    std::vector<RectangleTree*> tree;
    RectangleTree(){
        tree.resize(4, nullptr);
        tree[0] = this;
    }
    RectangleTree(const Rectangle& _rect): rect(_rect){
        tree.resize(4, nullptr);
        tree[0] = this;
    }
    void add_subtree(int idx) {
        this->tree[idx] = new RectangleTree;
        this->tree[idx]->tree[0] = this;
    }
};

clock_t start_time;
clock_t cur_time;

random_device seed_gen;
mt19937 mt(seed_gen());
// bmath::normal_distribution<> distribution(0, 1);

bool out_bound(int nx, int ny, int n)
{
    return (nx < 0 || ny < 0 || nx >= n || ny >= n);
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
        Grid points;
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
        State(const Parameters &param) {
            n = param.N;
            m = param.M;
            center = (param.N-1)/2;
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
                        #if _DEBUG
                        // std::cerr << "FOUND AT (" << nx_point.x << "," << nx_point.y << "), from[" << x << "," << y << "]" << std::endl;
                        #endif
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

            return ret;
        }

        bool go_till_find(int x, int y, int dir, Point& nx_point, int sx, int sy) {
            bool found = false;
            int cx = x;
            int cy = y;
            while(!found) {
                int nx = cx + dx8[dir];
                int ny = cy + dy8[dir];
                if(out_bound(nx, ny, n)) {
                    return false;
                }
                if(points[ny][nx].used[(dir+4)%8]) {
                    // if the current position has an edge of rectangl along the direction, stop the search.
                    // std::cerr << "INVALID!!!!" << std::endl;
                    return false;
                }
                if(points[ny][nx].exist or (nx == sx and ny == sy)) {
                    nx_point = points[ny][nx];
                    return true;
                }
                cx = nx;
                cy = ny;
            }
            return found;
        }

        void add_rectangle(Rectangle rect) {
            rectangles.push_back(rect);
            this->_update_points_on_rectangles(rect);
        }

        std::vector<Rectangle> recursive_point_search(Point start_point, int depth){
            std::vector<Rectangle> now;
            Grid grid = points;
            bool flag = true;
            std::cerr << "START RECURSIVE SEARCH" << std::endl;
            this->_recursive_search_ver3(start_point, start_point, depth, grid, now, flag);
            std::cerr << "now size: " << now.size() << std::endl;
            if(!now.empty()) {
                for(Rectangle rect: now) {
                    std::cerr << "OK?" << std::endl;
                    this->_update_points_on_rectangles(rect);
                }
            }
            return now;
        }
        
    private:
        void _update_points_on_edge(const Point &from, const Point &to, Grid& grid) {
            int dx = to.x-from.x;
            int dy = to.y-from.y;
            int d = max(abs(dx), abs(dy));
            std::cerr << d << " " << dx << " " << dy << std::endl;
            int ex = dx/d;
            int ey = dy/d;
            int dir = direction_map[{ex, ey}];
            for(int j = 0; j <= d; j++) {
                int cx = from.x + j * dx8[dir];
                int cy = from.y + j * dy8[dir];
                if(j < d){
                    grid[cy][cx].used[dir] = true;
                }
                if(j > 0) {
                    grid[cy][cx].used[(dir+4)%8] = true;
                }
                // std::cerr << "used: " << points[cy][cx].used[dir] << " " << points[cy][cx].used[(dir+4)%8] << std::endl;
                // for(bool b: points[cy][cx].used) std::cerr << b << " ";
                // std::cerr << std::endl;
            }
        }

        void _update_points_on_rectangles(Rectangle &rect) {
            rect.show();
            for(int i = 0; i < 4; i++) {
                Point from = rect.vert[i];
                Point to = rect.vert[(i+1)%4];
                this->_update_points_on_edge(from, to, this->points);
            }
        }

        Grid _update_points_on_rectangles(Grid grid, Rectangle &rect) {
            Grid ret = grid;
            for(int i = 0; i < 4; i++) {
                Point from = rect.vert[i];
                Point to = rect.vert[(i+1)%4];
                this->_update_points_on_edge(from, to, ret);
            }
            return ret;
        }

        // keep track of the state of grid.
        // TODO implement direction argument.
        void _recursive_search_ver3(Point start_point, Point pre_point, int depth, Grid &grid, std::vector<Rectangle>& cur_rec_list, bool &construct_flag) {
            pre_point = start_point;
            for(int dir = 0; dir < 8; dir++) {
                std::cerr << "DIR: " << dir << ", start: (" << start_point.x << "," << start_point.y << ")" << std::endl;
                int e1_dir = dir;
                int e2_dir = (dir + 2) % 8;
                if(grid[start_point.y][start_point.x].used[e1_dir] or grid[start_point.y][start_point.x].used[e2_dir]) {
                    continue;
                }
                int ce1p_x = start_point.x + dx8[e1_dir];
                int ce1p_y = start_point.y + dy8[e1_dir];
                int ce2p_x = start_point.x + dx8[e2_dir];
                int ce2p_y = start_point.y + dy8[e2_dir];
                int crossp_x = start_point.x + dx8[e1_dir] + dx8[e2_dir];
                int crossp_y = start_point.y + dy8[e1_dir] + dy8[e2_dir];

                while(!out_bound(ce1p_x, ce1p_y, n)) {
                    ce2p_x = start_point.x + dx8[e2_dir];
                    ce2p_y = start_point.y + dy8[e2_dir];
                    crossp_x = ce1p_x + dx8[e2_dir];
                    crossp_y = ce1p_y + dy8[e2_dir];
                    std::cerr << "POINTS: (" << ce1p_x << ", " << ce1p_y << "), ("<< ce2p_x << "," << ce2p_y << "), (" << crossp_x << "," << crossp_y << ")"<< std::endl;
                    Point p2 = grid[ce1p_y][ce1p_x];
                    while(!out_bound(ce2p_x, ce2p_y, n) and !out_bound(crossp_x, crossp_y, n)) {
                        Point p3 = grid[ce2p_y][ce2p_x];
                        Point p4 = grid[crossp_y][crossp_x];
                        // std::cerr << "OK!" << std::endl;
                    std::cerr << "POINTS: (" << ce1p_x << ", " << ce1p_y << "), ("<< ce2p_x << "," << ce2p_y << "), (" << crossp_x << "," << crossp_y << ")"<< std::endl;
                        if(depth == 0 and p2.exist and p3.exist and p4.exist) {
                            // std::cerr << "depth is zero!" << std::endl;
                            cur_rec_list.push_back(Rectangle{start_point, p2, p4, p3});
                            return;
                        } else if(depth == 0) {
                            ce2p_x += dx8[e2_dir];
                            ce2p_y += dy8[e2_dir];
                            crossp_x += dx8[e2_dir];
                            crossp_y += dy8[e2_dir];
                            continue;
                        }
                        std::vector<Rectangle> rec_list = cur_rec_list;


                        bool valid_rect = _check_rect(start_point, p2, p3, p4, grid, {e1_dir, e2_dir});
                        // std::cerr << "RECT is " << valid_rect << std::endl;
                        if(valid_rect) {
                            Rectangle rect = Rectangle{start_point, p2, p4, p3};
                            Grid cur_points_state = grid;
                            bool can_construct = true;

                            if(!grid[p2.y][p2.x].exist) {
                                std::cerr << "P2 SEARCH FROM (" << p2.x << "," << p2.y << ")" << std::endl;
                                cur_points_state[p2.y][p2.x].exist = true;
                                this->_recursive_search_ver3(p2, pre_point, depth-1, cur_points_state, rec_list, can_construct);
                            }
                            if(!grid[p3.y][p3.x].exist) {
                                std::cerr << "P3 SEARCH FROM (" << p3.x << "," << p3.y << ")" << std::endl;
                                cur_points_state[p3.y][p3.x].exist = true;
                                this->_recursive_search_ver3(p3, pre_point, depth-1, cur_points_state, rec_list, can_construct);
                            }
                            if(!grid[p4.y][p4.x].exist) {
                                std::cerr << "P4 SEARCH FROM (" << p4.x << "," << p4.y << ")" << std::endl;
                                cur_points_state[p4.y][p4.x].exist = true;
                                this->_recursive_search_ver3(p4, pre_point, depth-1, cur_points_state, rec_list, can_construct);
                            }
                            if(!can_construct) {
                                continue;
                            }
                            grid = cur_points_state;
                            cur_rec_list = rec_list;
                            return;
                        }
                        if(grid[p3.y][p3.x].exist) {
                            break;
                        }
                        ce2p_x += dx8[e2_dir];
                        ce2p_y += dy8[e2_dir];
                        crossp_x += dx8[e2_dir];
                        crossp_y += dy8[e2_dir];
                    }
                    // std::cerr << "POINTS: (" << ce1p_x << ", " << ce1p_y << "), ("<< ce2p_x << "," << ce2p_y << "), (" << crossp_x << "," << crossp_y << ")"<< std::endl;
                    if(grid[p2.y][p2.x].exist) {
                        break;
                    }
                    ce1p_x += dx8[e1_dir];
                    ce1p_y += dy8[e1_dir];
                }
            }
            construct_flag = false;
        } // ver_3

        // [start points]: start point of search on grid.
        // [remains]: The number of remain search depth
        // [grid]: keep track of the current situation
        void _recursive_search_ver2(std::vector<Point> start_points, std::vector<Point> pre_points, std::vector<int> remains, Grid grid, std::vector<Rectangle> cur_rec_list, std::vector<std::vector<Rectangle>> &res) {
            if(start_points.empty()) {
                std::cerr << "start_points does not contain any value!" << std::endl;
                std::cerr << "function: void _recursive_search_ver2" << std::endl;
                exit(0);
            }
            // First situation. All depth is equall
            if(start_points.size() == 1) {
                Point start_point = start_points[0];
                pre_points = std::vector<Point>(3, start_point);
                for(int dir = 0; dir < 8; dir++) {
                    int e1_dir = dir;
                    int e2_dir = (dir + 2) % 8;
                    if(grid[start_point.y][start_point.x].used[e1_dir] or grid[start_point.y][start_point.x].used[e2_dir]) {
                        continue;
                    }
                    int ce1p_x = start_point.x + dx8[e1_dir];
                    int ce1p_y = start_point.y + dy8[e1_dir];
                    int ce2p_x = start_point.x + dx8[e2_dir];
                    int ce2p_y = start_point.y + dy8[e2_dir];
                    int crossp_x = start_point.x + dx8[e1_dir] + dx8[e2_dir];
                    int crossp_y = start_point.y + dy8[e1_dir] + dy8[e2_dir];

                    while(!out_bound(ce1p_x, ce1p_y, n)) {
                        
                        ce2p_x = start_point.x + dx8[e2_dir];
                        ce2p_y = start_point.y + dy8[e2_dir];
                        crossp_x = ce1p_x + dx8[e2_dir];
                        crossp_y = ce1p_y + dy8[e2_dir];
                        Point &p2 = grid[ce1p_y][ce1p_x];
                        // Point p2 = Point{ce1p_x, ce1p_y, cur_points_state[ce1p_y][ce1p_x].exist};

                        while(!out_bound(ce2p_x, ce2p_y, n)) {
                            Point &p3 = grid[ce2p_y][ce2p_x];
                            // Point p3 = Point{ce2p_x, ce2p_y, cur_points_state[ce2p_y][ce2p_x].exist};
                            Point &p4 = grid[crossp_y][crossp_x];
                            // Point p4 = Point{crossp_x, crossp_y, cur_points_state[crossp_y][crossp_x].exist};
                            pre_points = start_points;
                            bool valid_rect = _check_rect(start_point, p2, p3, p4, grid, {e1_dir, e2_dir});
                            if(valid_rect) {
                                Rectangle rect = Rectangle{start_point, p2, p3, p4};
                                Grid cur_points_state = grid;
                                std::vector<int> nx_remains = remains;
                                if(grid[p2.y][p2.x].exist) {
                                    nx_remains[0] = 0;
                                }
                                if(grid[p3.y][p3.x].exist) {
                                    nx_remains[1] = 0;
                                }
                                if(grid[p4.y][p4.x].exist) {
                                    nx_remains[2] = 0;
                                }
                                cur_rec_list.push_back(rect);
                                
                                // cur_points_state = this->_update_points_on_rectangles(grid, rect);
                                std::vector<Point> nx_points = {p2, p3, p4};
                                this->_recursive_search_ver2(nx_points, pre_points, nx_remains, cur_points_state, cur_rec_list, res);

                                if(!grid[p2.y][p2.x].exist and nx_remains[0] != 0) {
                                    cur_points_state[p2.y][p2.x].exist = true;
                                    std::vector<int> nx_r = nx_remains;
                                    nx_r[0]--;
                                    this->_recursive_search_ver2(nx_points, pre_points, nx_r, cur_points_state, cur_rec_list, res);
                                }
                                if(!grid[p3.y][p3.x].exist and nx_remains[1] != 0) {
                                    cur_points_state[p3.y][p3.x].exist = true;
                                    std::vector<int> nx_r = nx_remains;
                                    nx_r[1]--;
                                    this->_recursive_search_ver2(nx_points, pre_points, nx_r, cur_points_state, cur_rec_list, res);
                                }
                                if(!grid[p4.y][p4.x].exist and nx_remains[2] != 0) {
                                    cur_points_state[p4.y][p4.x].exist = true;
                                    std::vector<int> nx_r = nx_remains;
                                    nx_r[2]--;
                                    this->_recursive_search_ver2(nx_points, pre_points, nx_r, cur_points_state, cur_rec_list, res);
                                }
                            
                            }
                            if(grid[p3.y][p3.x].exist) {
                                break;
                            }
                            ce2p_x += dx8[e2_dir];
                            ce2p_y += dy8[e2_dir];
                            crossp_x += dx8[e2_dir];
                            crossp_y += dy8[e2_dir];
                        }
                        if(grid[p2.y][p2.x].exist) {
                            break;
                        }
                        ce1p_x += dx8[e1_dir];
                        ce1p_y += dy8[e1_dir];
                    }
                }
            } else { // second situation. In the recursive search 
                pre_points = start_points;
                for(int p_idx = 0; p_idx < 3; p_idx++) {
                    Point start_point = start_points[p_idx];
                    for(int dir = 0; dir < 8; dir++) {
                        int e1_dir = dir;
                        int e2_dir = (dir + 2) % 8;
                        if(grid[start_point.y][start_point.x].used[e1_dir] or grid[start_point.y][start_point.x].used[e2_dir]) {
                            continue;
                        }
                        int ce1p_x = start_point.x + dx8[e1_dir];
                        int ce1p_y = start_point.y + dy8[e1_dir];
                        int ce2p_x = start_point.x + dx8[e2_dir];
                        int ce2p_y = start_point.y + dy8[e2_dir];
                        int crossp_x = start_point.x + dx8[e1_dir] + dx8[e2_dir];
                        int crossp_y = start_point.y + dy8[e1_dir] + dy8[e2_dir];

                        while(!out_bound(ce1p_x, ce1p_y, n)) {
                            
                            ce2p_x = start_point.x + dx8[e2_dir];
                            ce2p_y = start_point.y + dy8[e2_dir];
                            crossp_x = ce1p_x + dx8[e2_dir];
                            crossp_y = ce1p_y + dy8[e2_dir];
                            Point &p2 = grid[ce1p_y][ce1p_x];
                            // Point p2 = Point{ce1p_x, ce1p_y, cur_points_state[ce1p_y][ce1p_x].exist};

                            while(!out_bound(ce2p_x, ce2p_y, n)) {
                                Point &p3 = grid[ce2p_y][ce2p_x];
                                // Point p3 = Point{ce2p_x, ce2p_y, cur_points_state[ce2p_y][ce2p_x].exist};
                                Point &p4 = grid[crossp_y][crossp_x];
                                // Point p4 = Point{crossp_x, crossp_y, cur_points_state[crossp_y][crossp_x].exist};
                                pre_points = start_points;
                                bool valid_rect = _check_rect(start_point, p2, p3, p4, grid, {e1_dir, e2_dir});
                                if(valid_rect) {
                                    Rectangle rect = Rectangle{start_point, p2, p3, p4};
                                    Grid cur_points_state = grid;
                                    std::vector<int> nx_remains = remains;
                                    if(grid[p2.y][p2.x].exist) {
                                        nx_remains[0] = 0;
                                    }
                                    if(grid[p3.y][p3.x].exist) {
                                        nx_remains[1] = 0;
                                    }
                                    if(grid[p4.y][p4.x].exist) {
                                        nx_remains[2] = 0;
                                    }
                                    cur_rec_list.push_back(rect);
                                    
                                    std::vector<Point> nx_points = pre_points;
                                    // nx_points[p_idx] = 

                                    if(!grid[p2.y][p2.x].exist and nx_remains[0] != 0) {
                                        cur_points_state[p2.y][p2.x].exist = true;
                                        
                                        std::vector<int> nx_r = nx_remains;
                                        nx_r[0]--;
                                        this->_recursive_search_ver2(nx_points, pre_points, nx_r, cur_points_state, cur_rec_list, res);
                                    }
                                }
                                if(grid[p3.y][p3.x].exist) {
                                    break;
                                }
                                ce2p_x += dx8[e2_dir];
                                ce2p_y += dy8[e2_dir];
                                crossp_x += dx8[e2_dir];
                                crossp_y += dy8[e2_dir];
                            }
                            if(grid[p2.y][p2.x].exist) {
                                break;
                            }
                            ce1p_x += dx8[e1_dir];
                            ce1p_y += dy8[e1_dir];
                        }
                    }
                    
                }
            }

        }
        // TODO implement this codes    
        void _recursive_search(Point start_point, int d, Grid grid, std::vector<Rectangle> cur_rec_list, std::vector<std::vector<Rectangle>> &res) {            
            for(int dir = 0; dir < 8; dir++) {
                int edge1_dir = dir;
                int edge2_dir = (dir + 2) % 8;
                if(start_point.used[edge1_dir] || start_point.used[edge2_dir]) {
                    continue;
                }
                int ce1p_x = start_point.x + dx8[edge1_dir];
                int ce1p_y = start_point.y + dy8[edge1_dir];
                int ce2p_x = start_point.x + dx8[edge2_dir];
                int ce2p_y = start_point.y + dy8[edge2_dir];
                int crossp_x = start_point.x + dx8[edge1_dir] + dx8[edge2_dir];
                int crossp_y = start_point.y + dy8[edge1_dir] + dy8[edge2_dir];
                // Grid cur_points_state = grid;
                // Starting from start_point, searching the rectangle's vertices
                while(!out_bound(ce1p_x, ce1p_y, n)) {
                    
                    ce2p_x = start_point.x + dx8[edge2_dir];
                    ce2p_y = start_point.y + dy8[edge2_dir];
                    crossp_x = ce1p_x + dx8[edge2_dir];
                    crossp_y = ce1p_y + dy8[edge2_dir];
                    Point &p2 = grid[ce1p_y][ce1p_x];
                    // Point p2 = Point{ce1p_x, ce1p_y, cur_points_state[ce1p_y][ce1p_x].exist};

                    while(!out_bound(ce2p_x, ce2p_y, n)) {
                        Point &p3 = grid[ce2p_y][ce2p_x];
                        // Point p3 = Point{ce2p_x, ce2p_y, cur_points_state[ce2p_y][ce2p_x].exist};
                        Point &p4 = grid[crossp_y][crossp_x];
                        // Point p4 = Point{crossp_x, crossp_y, cur_points_state[crossp_y][crossp_x].exist};
                        bool valid_rect = _check_rect(start_point, p2, p3, p4, grid, {edge1_dir, edge2_dir});
                        std::cerr << valid_rect << std::endl;
                        if(valid_rect) {
                            Rectangle rect = Rectangle{start_point, p2, p4, p3};
                            rect.show();
                            Grid cur_points_state = grid;
                            cur_rec_list.push_back(rect);
                            if(d == 0 and p2.exist and p3.exist and p4.exist) {
                                res.push_back(cur_rec_list); // add rectangles list in result list
                            } else if (d == 0) {
                                return ;
                            } else {    
                                cur_points_state = this->_update_points_on_rectangles(grid, rect);
                                cur_points_state[p2.y][p2.x].exist = true;
                                cur_points_state[p3.y][p3.x].exist = true;
                                cur_points_state[p4.y][p4.x].exist = true;
                                if(!grid[p2.y][p2.x].exist) {
                                    this->_recursive_search(p2, d-1, cur_points_state, cur_rec_list, res);
                                }
                                if(!grid[p3.y][p3.x].exist) {
                                    this->_recursive_search(p3, d-1, cur_points_state, cur_rec_list, res);
                                }
                                if(!grid[p4.y][p4.x].exist) {
                                    this->_recursive_search(p4, d-1, cur_points_state, cur_rec_list, res);
                                }
                            }
                        }
                        if(grid[p3.y][p3.x].exist) {
                            break;
                        }
                        ce2p_x += dx8[edge2_dir];
                        ce2p_y += dy8[edge2_dir];
                        crossp_x += dx8[edge2_dir];
                        crossp_y += dy8[edge2_dir];
                    }
                    if(grid[p2.y][p2.x].exist) {
                        break;
                    }
                    ce1p_x += dx8[edge1_dir];
                    ce1p_y += dy8[edge1_dir];
                }
            }
            return ;
        }

        bool _check_rect(const Point &p1, const Point &p2, const Point &p3, const Point &p4, const Grid &grid, pair<int,int> dir) {
            // dir.first // e12 // e34
            // dir.second // e13 // e24
            // p2->p4
            int e24x = p2.x + dx8[dir.second];
            int e24y = p2.y + dy8[dir.second];
            int e34x = p3.x + dx8[dir.first];
            int e34y = p3.y + dy8[dir.first];
            std::cerr << "CHECK RECTANGLE" << std::endl;
            while(true) {
                if(e24x == p4.x and e24y == p4.y) {
                    break;
                }
                if(grid[e24y][e24x].used[(dir.second+4)%8] || grid[e24y][e24x].exist) {
                    // std::cerr << "OUT!" << std::endl;
                    return false;
                }
                e24x += dx8[dir.second];
                e24y += dy8[dir.second];
                // std::cerr << e24x << " " << e24y << std::endl;
            }
            // p3->p4
            while(true) {
                if(e34x == p4.x and e34y == p4.y) {
                    break;
                }
                if(grid[e34y][e34x].used[(dir.first+4)%8] || grid[e34y][e34x].exist) {
                    return false;
                }
                e34x += dx8[dir.first];
                e34y += dy8[dir.first];
                // std::cerr << e34x << " " << e34y << std::endl;
            }
            return true;
        }
};

struct Result {
    int K;
    std::vector<Rectangle> rects;
    std::vector<Point> points;
    long long score;
    Result(){}
    Result(std::vector<Rectangle> r): K((int)r.size()), rects(r){}
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
        Result best_result;

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
            best_result.score = -1;
            std::cerr << "N: " << n << std::endl;
            std::cerr << "M: " << m <<std::endl;

            std::cerr << "## FINISHED INIT ##" << std::endl;
        }


        void show_result() {
            int K = best_result.K;
            cout << K << endl;
            for(Rectangle r: best_result.rects) {
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

        // third solution
        Result search(bool& run_flag) {
            Result ret;
            while(!run_flag) {
                std::vector<Rectangle> rec_list;
                Point start;
                #if _TEST
                uniform_int_distribution<> randN(0, params.N-1);
                start = Point{randN(mt), randN(mt)};
                if(state.points[start.y][start.x].exist) {
                    continue;
                }
                #else 
                
                #endif
                std::cerr << "ENTERING RECURSIVE SEARCH" << std::endl;
                rec_list = this->state.recursive_point_search(start, 1);
                if(rec_list.empty()) {
                    std::cerr << "NOT FOUND!!!" << std::endl; 
                    continue;
                }
                // reverse(rec_list.begin(), rec_list.end());
                for(Rectangle rec: rec_list) {
                    ret.rects.push_back(rec);
                }
                run_flag = this->time_flag_check();
            }
            ret.K = ret.rects.size();
            return ret;
        }

        // greedy solution
        Result baseline(bool& run_flag) {
            std::cerr << "## ENTERING baseline codes" << std::endl;
            cur_time = clock();
            bool time_limit_flag = false;
            bool found_flag = true;
            while(!time_limit_flag and found_flag) {
                found_flag = false;
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
                            found_flag = true;
                            // state.rectangles.push_back(nx_rec);
                        }
                        time_limit_flag = this->time_flag_check();
                        // std::cerr << (double)cur_time / CLOCKS_PER_SEC<< std::endl;
                    }
                }
                if(!found_flag) {
                    std::cerr << "NOT FOUND!!!! TIME: " << (double)cur_time/CLOCKS_PER_SEC << std::endl;
                    
                }
            }
            run_flag = time_limit_flag;
            long long score = this->calc_score();
            return Result{(int)state.rectangles.size(), state.rectangles, score};
        }

        //start search from center point. if found point, resart from center
        Result baseline_ver2(bool &run_flag) {
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
            bool found_flag = true;
            bool restart_flag = true;
            while(time_limit_flag and found_flag) {
                found_flag = false;
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
                        found_flag=true;
                        // state.rectangles.push_back(nx_rec);
                    }
                    time_limit_flag = this->time_flag_check();
                    // std::cerr << (double)cur_time / CLOCKS_PER_SEC<< std::endl;
                    
                }
            }
            long long score = this->calc_score();
            run_flag = time_limit_flag;
            return Result{(int)state.rectangles.size(), state.rectangles, score};
            
        }

        void run() {
            // >--------------------------- codes ---------------------------<
            bool run_flag = false;
            while(!run_flag) {
                Result res = this->search(run_flag);
                #if _DEBUG
                if(res.rects.empty()) {
                    std::cerr << "THE RESULT DOES NOT CONTAIN VALUE" << std::endl;
                    exit(0);
                }
                #endif
                if(res.score > best_result.score) {
                    best_result = res;
                }
                result_list.push_back(res);
                if(run_flag)break;
                this->init();
            }
            // >-------------------------------------------------------------< 
            this->show_result();
            this->summary();
            #if _DEBUG
            this->debug();
            #endif
        }

        bool time_flag_check() {
            cur_time = clock();
            return ((double)cur_time - start_time) / CLOCKS_PER_SEC > this -> params.TIME_LIMIT_90;
        }

        void summary() {
            std::cerr << "###################### SUMMARY ######################" << std::endl;
            std::cerr << "ELAPSED TIME: " << ((double)clock()-start_time) / CLOCKS_PER_SEC << " s" << std::endl;
            std::cerr << "BEST SCORE: " << best_result.score << std::endl;

            std::cerr << "#####################################################" << std::endl;
        }

        void init() {
            state.points.clear();
            state.rectangles.clear();
            state.stack_points.clear();
            state = State(params);
            return ;
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