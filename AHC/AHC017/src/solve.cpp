#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
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
#include <limits>

// #include <boost/math/distributions.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/lu.hpp>
// #include <boost/numeric/ublas/triangular.hpp>
// #include <boost/numeric/ublas/vector.hpp>

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
#define Vec std::vector<int>
// #define dmatrix boost::numeric::ublas::matrix<float>
// #define dvector boost::numeric::ublas::vector<double>

using namespace std;

template <typename T>
T relu(T x)
{
    return std::max(x, 0);
}

template <typename T>
T sigmoid(T x)
{
    return (double)1 / (1 + std::exp(x));
}

template <typename T>
T silu(T x)
{
    return x * sigmoid(x);
}

bool operator==(const pii &lhs, const pii &rhs)
{
    return ((lhs.first == rhs.first and lhs.second == rhs.second) or (lhs.first == rhs.second and lhs.second == rhs.first));
}

// dmatrix scatter_add(dmatrix src, std::vector<int> indices, int dim_size) {
//     int r = src.size1();
//     int c = src.size2();
//     dmatrix ret(dim_size, c);
//     for(int i = 0; i < r; i++) {
//         int t_idx = indices[i];
//         for(int j = 0; j < c; j++) {
//             ret(t_idx, j) = src(i, j);
//         }
//     }
//     return ret;
// }

// // to do
// dmatrix load_weights() {
//     dmatrix ret;
//     return
// }

// // input for Net
// struct Input {
//     dmatrix edge; // [edge_size, 2]
//     bool adj[1001][1001];
//     std::vector<std::vector<int>> graph;
//     std::vector<std::vector<int>> dist;
//     std::vector<bool> mask;
//     int time;
// };

// ?
// class Linear {
// public:

// private:
// };

// class Encoder {
// public:
//     std::vector<std::vector<int>> edge_list; // edge ids for any node i

//     Encoder() {}
// private:
// };

// class Net {
// public:
//     Net() {}
// private:
// };

// global variable
namespace glb
{
    int N, M, D, K;
    int U[3001], V[3001], W[3001];
    int X[1001], Y[1001];
    std::vector<std::vector<pair<int, int>>> G; // (u (v,w))
    std::vector<std::vector<int>> dist;
    int edge_id_matrix[1001][1001]; // (u,v) -> eid
    const int R = 500;
};

const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, 1, 0, -1};
const int dcx[4] = {-1, -1, 1, 1}; // cross
const int dcy[4] = {-1, 1, 1, -1};
const int dx8[8] = {-1, -1, -1, 0, 1, 1, 1, 0}; // clockwise
const int dy8[8] = {-1, 0, 1, 1, 1, 0, -1, -1};

const int INF = 1 << 30;
const ll LINF = 1LL << 60;
const double epsilon = 1e-6;
const double PI = std::acos(-1);

Vec dijkstra(int start, const std::vector<std::vector<bool>> &deleted)
{
    typedef pair<int, pii> Edge;

    int n = glb::N;
    std::vector<int> dist(n + 1, INF);
    std::vector<bool> visited(n + 1, false);
    dist[start] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, start}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;

        for (pii ne : glb::G[from])
        {
            int to = ne.first;
            if (deleted[from][to])
                continue;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return dist;
}

std::tuple<Vec, Vec> dijkstra(int start)
{
    typedef pair<int, pii> Edge;

    int n = glb::N;
    Vec dist(n + 1, INF);
    Vec pre(n + 1, -1);
    std::vector<bool> visited(n + 1, false);
    dist[start] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, start}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;

        for (pii ne : glb::G[from])
        {
            int to = ne.first;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist, pre};
}

// blocked three edges
std::tuple<int, std::vector<int>> dijkstra(int s, int t, int e1, int e2, int e3)
{
    typedef pair<int, pii> Edge;
    pii e1p = {glb::U[e1], glb::V[e1]};
    pii e2p = {glb::U[e2], glb::V[e2]};
    pii e3p = {glb::U[e3], glb::V[e3]};
    int n = glb::N;
    std::vector<int> dist(n + 1, INF);
    std::vector<int> pre(n + 1, -1);

    dist[s] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, s}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;
        int cost = edge.first;
        if (cost > dist[t])
            break;
        if (cost > dist[from])
            continue;
        for (pii ne : glb::G[from])
        {
            int to = ne.first;
            if (pii{from, to} == e1p or pii{from, to} == e2p or pii{from, to} == e3p)
                continue;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist[t], pre};
}

// blocked two edges
std::tuple<int, std::vector<int>> dijkstra(int s, int t, int e1, int e2)
{
    typedef pair<int, pii> Edge;
    pii e1p = {glb::U[e1], glb::V[e1]};
    pii e2p = {glb::U[e2], glb::V[e2]};
    int n = glb::N;
    std::vector<int> dist(n + 1, INF);
    std::vector<int> pre(n + 1, -1);

    dist[s] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, s}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;
        int cost = edge.first;
        if (cost > dist[t])
            break;
        if (cost > dist[from])
            continue;
        for (pii ne : glb::G[from])
        {
            int to = ne.first;
            if (pii{from, to} == e1p or pii{from, to} == e2p)
                continue;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist[t], pre};
}

// blocked one edge
std::tuple<int, std::vector<int>> dijkstra(int s, int t, int e1)
{
    typedef pair<int, pii> Edge;
    pii e1p = {glb::U[e1], glb::V[e1]};
    int n = glb::N;
    std::vector<int> dist(n + 1, INF);
    std::vector<int> pre(n + 1, -1);
    dist[s] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, s}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;
        int cost = edge.first;
        if (cost > dist[t])
            break;
        if (cost > dist[from])
            continue;

        for (pii ne : glb::G[from])
        {
            int to = ne.first;
            if (pii{from, to} == e1p)
                continue;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist[t], pre};
}

std::tuple<int, std::vector<int>> dijkstra(int s, int t)
{
    typedef pair<int, pii> Edge;
    int n = glb::N;
    std::vector<int> dist(n + 1, INF);
    std::vector<int> pre(n + 1, -1);
    dist[s] = 0;
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> pq;
    pq.push(Edge{0, {-1, s}});
    while (!pq.empty())
    {
        auto edge = pq.top();
        pq.pop();
        int from = edge.second.second;
        int cost = edge.first;
        if (cost > dist[t])
            break;
        if (cost > dist[from])
            continue;

        for (pii ne : glb::G[from])
        {
            int to = ne.first;
            if (dist[to] > dist[from] + ne.second)
            {
                dist[to] = dist[from] + ne.second;
                pre[to] = from;
                pq.push(Edge{dist[to], {from, to}});
            }
        }
    }
    return {dist[t], pre};
}

struct Parameters
{

    double TL = 6.00;
    double TIME_LIMIT_90 = 5.40;
    double TIME_LIMIT_95 = 5.70;
    double DECAY;
    std::string file_name = _OUTPUT_FILE;

    int N;
    int M;

    Parameters()
    {
#ifdef _LOCAL
        this->TL = 40;
        this->TIME_LIMIT_90 = 36.;
        this->TIME_LIMIT_95 = 38.;
#endif
    }

    void save()
    {
        ofstream ofs;
        ofs.open(file_name, ios::out);
        cout << "TL: " << TL << "\n";
        cout << "DECAY: " << DECAY << "\n";
        ofs.close();
        return;
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

template <class T>
T pow2(T x)
{
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

void test()
{
    std::cerr << "### TEST ###\n";
    return;
}

double get_time()
{
    cur_time = clock();
    return (double)(cur_time - start_time) / CLOCKS_PER_SEC;
}

// # memo
// nihongo nyuuryoku dekinai
//
// # Idea
// find all closed cycle. calculate Surface are by Green's theoorem.
// Algorithm.
// 1. Find big cycle.
// 2. Enumerate edges on the cycle and find which one is very important.
//  (By cutting the edges, how much areas are added to the original cycle.)

struct Edge
{
    int id;
    int u, v;
    int w;
};

class Node
{
public:
    int id;
    int parent;
    Vec children;
    Vec edges;

    Node()
    {
        this->parent = -1;
    }
    Node(int _id) : id(_id)
    {
        this->parent = -1;
    }

private:
};

struct Triangle
{
    int s, t, u;
    std::vector<int> est, etu, eus;
    int length;
    Triangle(int _s, int _t, int _u, Vec _est, Vec _etu, Vec _eus, int _length)
    {
        this->s = _s;
        this->t = _t;
        this->u = _u;
        this->est = _est;
        this->etu = _etu;
        this->eus = _eus;
        this->length = _length;
    }
};

bool operator>(const Triangle &lhs, const Triangle &rhs)
{
    return lhs.length > rhs.length;
}

bool operator<(const Triangle &lhs, const Triangle &rhs)
{
    return lhs.length < rhs.length;
}

class Network
{
public:
    Vec nodes;
    std::map<int, Node> node_map;
    std::map<int, std::map<int, Vec>> edges;
    std::map<int, std::vector<pii>> graph; // u -> (v, eid)

    int city_cnt;

    Network()
    {
        for (int i = 0; i < glb::N; i++)
        {
            this->nodes.push_back(i);
            this->node_map[i] = Node(i);
        }
        for (int i = 0; i < glb::M; i++)
        {
            int u = glb::U[i];
            int v = glb::V[i];
            this->edges[u][v].push_back(i);
            this->edges[v][u].push_back(i);
            this->graph[u].push_back({v, i});
            this->graph[v].push_back({u, i});
        }
        this->city_cnt = glb::N;
    }

    void step(int episode)
    {
        std::vector<Triangle> triangles = this->find_triangle(episode);
        if (episode == 1)
            std::sort(triangles.begin(), triangles.end());
        this->merge_cities2(triangles);
    }

    int get_leader(int child)
    {
        int ret = child;
        while (this->node_map[ret].parent != -1)
        {
            ret = this->node_map[ret].parent;
        }
        return ret;
    }

    // Output function for python visualization.
    void show_nodes_pos()
    {
        auto dfs = [&](auto dfs, int nid, int &cnt) -> pii
        {
            pii ret = {0, 0};
            if (nid < glb::N)
            {
                cnt++;
                return {glb::X[nid], glb::Y[nid]};
            }
            else
            {
                for (int chid : this->node_map[nid].children)
                {
                    pii tmp = dfs(dfs, chid, cnt);
                    ret.first += tmp.first;
                    ret.second += tmp.second;
                }
                return ret;
            }
        };

        std::cout << this->nodes.size() << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        for (int nid : this->nodes)
        {

            int num = 0;
            std::cout << "NID: " << nid << std::endl;
            pii res = dfs(dfs, nid, num);
            std::cout << (double)res.first / num << " " << (double)res.second / num << std::endl;
        }
    }

    void show_edge_conn()
    {
        std::vector<pii> show_list;
        for (int s : this->nodes)
        {
            for (auto p : this->graph[s])
            {
                int t = p.first;
                if (s > t)
                    continue;
                show_list.push_back({s, t});
            }
        }
        std::sort(show_list.begin(), show_list.end());

        std::cout << show_list.size() << std::endl;
        for (pii p : show_list)
        {
            std::cout << p.first << " " << p.second << std::endl;
        }
    }

    std::vector<Triangle> find_triangle(int episode)
    {
        std::vector<Triangle> ret;
        int num_node = this->nodes.size();
        std::map<int, bool> used;

        for (int i = 0; i < num_node; i++)
        {
            int s = this->nodes[i];
            used[s] = true;
            std::map<int, std::map<int, bool>> used_edge1;
            for (auto e1 : this->graph[s])
            {
                int t = e1.first;
                int est = e1.second;
                int wst = glb::W[est];
                if (used[t] or used_edge1[s][t])
                    continue;
                used_edge1[s][t] = used_edge1[t][s] = true;
                std::map<int, std::map<int, bool>> used_edge2;
                for (auto e2 : this->graph[t])
                {
                    int u = e2.first;
                    int etu = e2.second;
                    int wtu = glb::W[etu];
                    if (used[u] or used_edge2[t][u])
                        continue;
                    used_edge2[t][u] = used_edge2[u][t] = true;
                    std::map<int, std::map<int, bool>> used_edge3;

                    for (auto e3 : this->graph[u])
                    {
                        int v = e3.first;
                        int euv = e3.second;
                        int wuv = glb::W[euv];
                        if (used_edge3[u][v])
                            continue;
                        used_edge3[v][u] = used_edge3[u][v] = true;
                        if (s == v)
                        {
                            std::vector<int> ests, etus, euss;
                            for (int eid : this->edges[s][t])
                                ests.push_back(eid);
                            for (int eid : this->edges[t][u])
                                etus.push_back(eid);
                            for (int eid : this->edges[u][s])
                                euss.push_back(eid);
                            if (episode == 1)
                                ret.push_back({s, t, u, ests, etus, euss, wst + wtu + wuv});
                            else
                                ret.push_back({s, t, u, ests, etus, euss, 0});
                        }
                    }
                }
            }
        }
        return ret;
    }

private:
    //
    // another algorithm
    std::vector<Triangle> find_triangle2(int episode)
    {
        std::vector<Triangle> ret;
        int num_node = this->nodes.size();
        std::map<int, bool> used;

        for (int i = 0; i < num_node; i++)
        {
            int s = this->nodes[i];
            used[s] = true;
            std::map<int, std::map<int, bool>> used_edge1;
            for (auto e1 : this->graph[s])
            {
                int t = e1.first;
                int est = e1.second;
                int wst = glb::W[est];
                if (used[t] or used_edge1[s][t])
                    continue;
                used_edge1[s][t] = used_edge1[t][s] = true;
                std::map<int, std::map<int, bool>> used_edge2;
                for (auto e2 : this->graph[t])
                {
                    int u = e2.first;
                    int etu = e2.second;
                    int wtu = glb::W[etu];
                    if (used[u] or used_edge2[t][u])
                        continue;
                    used_edge2[t][u] = used_edge2[u][t] = true;
                    int nsu = this->edges[s][u].size();
                    int ntu = this->edges[t][u].size();
                    if (nsu > 0 and ntu > 0)
                    {
                        int euv = this->edges[u][s][0];
                        int wuv = glb::W[euv];
                        std::vector<int> ests, etus, euss;
                        for (int eid : this->edges[s][t])
                            ests.push_back(eid);
                        for (int eid : this->edges[t][u])
                            etus.push_back(eid);
                        for (int eid : this->edges[u][s])
                            euss.push_back(eid);
                        if (episode == 1)
                            ret.push_back({s, t, u, ests, etus, euss, wst + wtu + wuv});
                        else
                            ret.push_back({s, t, u, ests, etus, euss, 0});
                    }
                }
            }
        }
        return ret;
    }

    void merge_cities(const std::vector<Triangle> &triangles)
    {
        std::map<int, bool> merged;

        for (Triangle tri : triangles)
        {
            int s, t, u;
            Vec est, etu, eus;
            s = tri.s;
            t = tri.t;
            u = tri.u;
            est = tri.est;
            etu = tri.etu;
            eus = tri.eus;
            if (merged[s] or merged[t] or merged[u])
                continue;
            std::map<int, bool> in_group;
            in_group[s] = in_group[t] = in_group[u] = true;
            std::vector<int> edge_id_list;
            edge_id_list.push_back(est[0]); // enough to see only first element.
            edge_id_list.push_back(etu[0]);
            edge_id_list.push_back(eus[0]);
            Node new_node(this->city_cnt);
            new_node.children.push_back(s);
            new_node.children.push_back(t);
            new_node.children.push_back(u);
            for (int eid : est)
                new_node.edges.push_back(eid);
            for (int eid : etu)
                new_node.edges.push_back(eid);
            for (int eid : eus)
                new_node.edges.push_back(eid);

            std::queue<int> q;
            q.push(est[0]);
            q.push(etu[0]);
            q.push(eus[0]);
            while (!q.empty())
            {
                int start_eid = q.front();
                q.pop();
                int u = this->get_leader(glb::U[start_eid]);
                int v = this->get_leader(glb::V[start_eid]);
                std::map<int, std::map<int, bool>> used_edge1;
                for (pii nei1 : this->graph[u])
                {
                    int s = nei1.first;

                    if (merged[s] or in_group[s] or used_edge1[u][s] or used_edge1[s][u])
                        continue;
                    used_edge1[u][s] = true;
                    std::map<int, std::map<int, bool>> used_edge2;
                    for (pii nei2 : this->graph[s])
                    {
                        int t = nei2.first;

                        if (merged[t] or used_edge2[s][t] or used_edge2[t][s])
                            continue;
                        used_edge2[s][t] = true;

                        if (t == v)
                        {
                            int n_add_edge = this->edges[u][s].size() + this->edges[s][t].size();
                            if (n_add_edge + (int)new_node.edges.size() <= glb::D)
                            {
                                // updates for new node
                                for (int eid : this->edges[u][s])
                                    new_node.edges.push_back(eid);
                                for (int eid : this->edges[s][t])
                                    new_node.edges.push_back(eid);

                                new_node.children.push_back(s);
                                in_group[s] = true;
                                q.push(this->edges[u][s][0]);
                                q.push(this->edges[s][v][0]);
                            }
                        }
                    }
                }
            } // while

            this->node_map[city_cnt] = new_node;
            for (int nid : new_node.children)
            {
                merged[nid] = true;
                this->node_map[nid].parent = new_node.id;
            }
            this->city_cnt++;
        } // triangle loop

        // update new nodes
        std::set<int> new_node_ids;
        for (int i = 0; i < glb::N; i++)
        {
            int nid = this->get_leader(i);
            new_node_ids.insert(nid);
        }
        Vec tmp;
        for (auto v : new_node_ids)
            tmp.push_back(v);
        std::sort(tmp.begin(), tmp.end());
        this->nodes = tmp;

        // update new edges and graph
        std::map<int, std::map<int, Vec>> new_edges;
        std::map<int, std::vector<pii>> new_graph;

        for (int i = 0; i < glb::M; i++)
        {
            int uid = this->get_leader(glb::U[i]);
            int vid = this->get_leader(glb::V[i]);
            if (uid == vid)
                continue;
            new_edges[uid][vid].push_back(i);
            new_edges[vid][uid].push_back(i);
            new_graph[uid].push_back(pii{vid, i});
            new_graph[vid].push_back(pii{uid, i});
        }

        this->edges = new_edges;
        this->graph = new_graph;
    }

    // another algorithm of merge cities
    void merge_cities2(const std::vector<Triangle> &triangles)
    {
        std::map<int, bool> merged;
        std::map<int, std::map<int, bool>> used_edge;
        std::vector<Node> new_nodes;
        std::map<int, bool> added_edge;
        std::map<int, bool> added_node;
        for (Triangle tri : triangles)
        {
            int s, t, u;
            Vec est, etu, eus;
            s = tri.s;
            t = tri.t;
            u = tri.u;
            est = tri.est;
            etu = tri.etu;
            eus = tri.eus;
            if (merged[s] or merged[t] or merged[u])
                continue;
            Node new_node(this->city_cnt);
            std::map<int, bool> same_group;
            new_node.children.push_back(s);
            if (added_node[s])
            {
                std::cerr << "This node is alread Added! (" << s << ")" << std::endl;
            }
            added_node[s] = true;
            new_node.children.push_back(t);
            if (added_node[t])
            {
                std::cerr << "This node is alread Added! (" << t << ")" << std::endl;
            }
            added_node[t] = true;
            new_node.children.push_back(u);
            if (added_node[u])
            {
                std::cerr << "this node is alread added! (" << u << ")" << std::endl;
            }
            added_node[u] = true;
            for (int eid : est)
            {
                if (added_edge[eid])
                {
                    std::cerr << "This edge is alread Added!" << std::endl;
                    std::cerr << "From for(int eid: est) section" << std::endl;
                }
                added_edge[eid] = true;
                new_node.edges.push_back(eid);
            }
            for (int eid : etu)
            {
                if (added_edge[eid])
                {
                    std::cerr << "This edge is alread Added!" << std::endl;
                    std::cerr << "From for(int eid: etu) section" << std::endl;
                }
                added_edge[eid] = true;
                new_node.edges.push_back(eid);
            }
            for (int eid : eus)
            {
                if (added_edge[eid])
                {
                    std::cerr << "This edge is alread Added!" << std::endl;
                    std::cerr << "From for(int eid: eus) section" << std::endl;
                }
                added_edge[eid] = true;
                new_node.edges.push_back(eid);
            }

            if ((int)new_node.edges.size() > glb::D)
                continue;

            used_edge[s][t] = used_edge[t][s] = true;
            used_edge[t][u] = used_edge[u][t] = true;
            used_edge[u][s] = used_edge[s][u] = true;
            same_group[s] = true;
            same_group[t] = true;
            same_group[u] = true;
            std::queue<int> q;
            q.push(est[0]);
            q.push(etu[0]);
            q.push(eus[0]);
            while (!q.empty())
            {
                int start_eid = q.front();
                q.pop();
                int u = this->get_leader(glb::U[start_eid]);
                int v = this->get_leader(glb::V[start_eid]);

                for (pii nei1 : this->graph[u])
                {
                    int s = nei1.first;

                    if (merged[s] or (used_edge[u][s] and used_edge[v][s]))
                        continue;

                    int nus = this->edges[u][s].size();
                    int nvs = this->edges[v][s].size();
                    int n_cur_edge = new_node.edges.size();

                    if (nus > 0 and nvs > 0)
                    {
                        if (used_edge[u][s])
                            nus = 0;
                        if (used_edge[v][s])
                            nvs = 0;

                        int total = n_cur_edge + nus + nvs;
                        // check around s node. this process is needed for collect merge
                        for (auto p : this->graph[s])
                        {
                            int t = p.first;
                            int nst = this->edges[s][t].size();
                            if (not used_edge[s][t] and same_group[t])
                                total += nst;
                        }
                        if (total <= glb::D)
                        {
                            if (not used_edge[u][s] and not used_edge[v][s])
                            {
                                new_node.children.push_back(s);
                                same_group[s] = true;
                                if (added_node[s])
                                {
                                    std::cerr << "this node is alread added! (" << s << "). [nuvs section]" << std::endl;
                                }
                                added_node[s] = true;
                            }

                            if (not used_edge[u][s])
                            {
                                for (int eid : this->edges[u][s])
                                {
                                    if (added_edge[eid])
                                    {
                                        std::cerr << "This edge is alread Added! (" << eid << ")" << std::endl;
                                        std::cerr << "From for(int eid: this->edges[u][s]) section" << std::endl;
                                    }
                                    added_edge[eid] = true;
                                    new_node.edges.push_back(eid);
                                }
                                q.push(this->edges[u][s][0]);
                                used_edge[u][s] = used_edge[s][u] = true;
                            }

                            if (not used_edge[v][s])
                            {
                                for (int eid : this->edges[v][s])
                                {
                                    if (added_edge[eid])
                                    {
                                        std::cerr << "This edge is alread Added! (" << eid << ")" << std::endl;
                                        std::cerr << "From for(int eid: this->edges[v][s]) section" << std::endl;
                                    }
                                    added_edge[eid] = true;
                                    new_node.edges.push_back(eid);
                                }
                                q.push(this->edges[v][s][0]);
                                used_edge[v][s] = used_edge[s][v] = true;
                            }

                            // collect eIDs around s
                            for (auto p : this->graph[s])
                            {
                                int t = p.first;
                                if (same_group[t] and not used_edge[s][t])
                                {
                                    for (int eid : this->edges[s][t])
                                    {
                                        if (added_edge[eid])
                                        {
                                            std::cerr << "This edge is alread Added! (" << eid << ")" << std::endl;
                                            std::cerr << "From collect section" << std::endl;
                                        }
                                        added_edge[eid] = true;
                                        new_node.edges.push_back(eid);
                                    }
                                    used_edge[s][t] = used_edge[t][s] = true;
                                    q.push(this->edges[s][t][0]);
                                }
                            }
                        }
                    }
                    // else if (nus >= 2 and not used_edge[u][s])
                    // {
                    //     int total = n_cur_edge + nus;
                    //     // check around s node. this process is needed for collect merge
                    //     for (auto p : this->graph[s])
                    //     {
                    //         int t = p.first;
                    //         int nst = this->edges[s][t].size();
                    //         if (not used_edge[s][t] and same_group[t])
                    //             total += nst;
                    //     }
                    //     if (total <= glb::D)
                    //     {
                    //         new_node.children.push_back(s);
                    //         same_group[s] = true;
                    //         if (added_node[s])
                    //         {
                    //             std::cerr << "this node is alread added! (" << s << "). [nus section]" << std::endl;
                    //         }
                    //         added_node[s] = true;

                    //         // collect eids around s
                    //         for (auto p : this->graph[s])
                    //         {
                    //             int t = p.first;
                    //             int eid = p.second;
                    //             if (same_group[t] and not used_edge[s][t])
                    //             {
                    //                 for (int eid : this->edges[s][t])
                    //                 {
                    //                     if (added_edge[eid])
                    //                     {
                    //                         std::cerr << "This edge is alread Added! (" << eid << ")" << std::endl;
                    //                         std::cerr << "From for(int eid: this->edges[s][t]) section" << std::endl;
                    //                     }
                    //                     added_edge[eid] = true;
                    //                     new_node.edges.push_back(eid);
                    //                 }
                    //                 used_edge[s][t] = used_edge[t][s] = true;
                    //                 q.push(this->edges[s][t][0]);
                    //             }
                    //         }

                    //         for (int eid : this->edges[u][s])
                    //         {
                    //             if (added_edge[eid])
                    //             {
                    //                 std::cerr << "This edge is alread Added! (" << eid << ")" << std::endl;
                    //                 std::cerr << "From for(int eid: this->edges[u][s]) section" << std::endl;
                    //             }
                    //             added_edge[eid] = true;
                    //             new_node.edges.push_back(eid);
                    //         }
                    //         q.push(this->edges[u][s][0]);
                    //         used_edge[u][s] = used_edge[s][u] = true;
                    //     }
                    // }
                }

                // edges around v
                // for (pii nei1 : this->graph[v])
                // {
                //     int s = nei1.first;
                //     int eid1 = nei1.second;

                //     if (merged[s] or (used_edge[s][u] and used_edge[s][v]))
                //         continue;

                //     int nvs = this->edges[v][s].size();
                //     int n_cur_edge = new_node.edges.size();

                //     if (nvs >= 2 and not used_edge[v][s])
                //     {
                //         int total = n_cur_edge + nvs;
                //         // check around s node. this process is needed for collect merge
                //         for (auto p : this->graph[s])
                //         {
                //             int t = p.first;
                //             int nst = this->edges[s][t].size();
                //             if (not used_edge[s][t] and same_group[t])
                //                 total += nst;
                //         }
                //         if (total <= glb::D)
                //         {
                //             new_node.children.push_back(s);
                //             same_group[s] = true;
                //             if (added_node[s])
                //             {
                //                 std::cerr << "this node is alread added! (" << s << "). [nvs section]" << std::endl;
                //             }
                //             added_node[s] = true;

                //             // collect eids around s
                //             for (auto p : this->graph[s])
                //             {
                //                 int t = p.first;
                //                 int eid = p.second;
                //                 if (same_group[t] and not used_edge[s][t])
                //                 {
                //                     for (int eid : this->edges[s][t])
                //                     {
                //                         if (added_edge[eid])
                //                         {
                //                             std::cerr << "This edge is alread Added! (" << eid << ")" << std::endl;
                //                             std::cerr << "From for(int eid: this->edges[s][t]) section" << std::endl;
                //                         }
                //                         added_edge[eid] = true;

                //                         new_node.edges.push_back(eid);
                //                     }
                //                     used_edge[s][t] = used_edge[t][s] = true;
                //                     q.push(this->edges[s][t][0]);
                //                 }
                //             }

                //             for (int eid : this->edges[v][s])
                //                 new_node.edges.push_back(eid);
                //             q.push(this->edges[v][s][0]);
                //             used_edge[v][s] = used_edge[s][v] = true;
                //         }
                //     }
                // }
            } // while

            for (int nid : new_node.children)
                merged[nid] = true;
            new_nodes.push_back(new_node);

            this->city_cnt++;
        } // triangle loop

        // update new nodes
        for (Node new_node : new_nodes)
        {
            this->node_map[new_node.id] = new_node;
            for (int chid : new_node.children)
                this->node_map[chid].parent = new_node.id;
        }

        std::set<int> new_node_ids;
        for (int i = 0; i < glb::N; i++)
        {
            int nid = this->get_leader(i);
            new_node_ids.insert(nid);
        }
        Vec tmp;
        for (auto v : new_node_ids)
            tmp.push_back(v);
        std::sort(tmp.begin(), tmp.end());
        this->nodes = tmp;

        // update new edges and graph
        std::map<int, std::map<int, Vec>> new_edges;
        std::map<int, std::vector<pii>> new_graph;

        for (int i = 0; i < glb::M; i++)
        {
            int uid = this->get_leader(glb::U[i]);
            int vid = this->get_leader(glb::V[i]);
            if (uid == vid)
                continue;
            new_edges[uid][vid].push_back(i);
            new_edges[vid][uid].push_back(i);
            new_graph[uid].push_back(pii{vid, i});
            new_graph[vid].push_back(pii{uid, i});
        }

        this->edges = new_edges;
        this->graph = new_graph;
    }
};

struct Result
{
    std::vector<int> assign;
    double eval_score;
    Result()
    {
        this->eval_score = 1e18;
    }
    Result(std::vector<int> _assign, long long _eval_score)
    {
        this->assign = _assign;
        this->eval_score = _eval_score;
    }
};

class Solver
{
public:
    Parameters params;
    Network network;
    Result best_result;

    Solver()
    {
        std::cerr << "## ENTERING INITIALIZATION " << std::endl;
        std::cin >> glb::N >> glb::M >> glb::D >> glb::K;
        glb::G.resize(glb::N);
        glb::dist.resize(glb::N);
        for (int i = 0; i < glb::M; i++)
        {
            int u, v, w;
            cin >> u >> v >> w;
            u--;
            v--;
            glb::U[i] = u;
            glb::V[i] = v;
            glb::W[i] = w;
            glb::G[u].push_back({v, w});
            glb::G[v].push_back({u, w});
            glb::edge_id_matrix[u][v] = glb::edge_id_matrix[v][u] = i;
        }
        for (int i = 0; i < glb::N; i++)
        {
            std::cin >> glb::X[i] >> glb::Y[i];
        }

        // calculating shortest length. too heavy. Need reducing complexity.
        for (int u = 0; u < glb::N; u++)
        {
            // glb::dist[u] = dijkstra(u);
        }

        this->network = Network();

        std::cerr << "## FINISHED INIT ##" << std::endl;
    }

    // find small triangle and assign 1, 2, 3
    void solve2()
    {
        std::vector<Triangle> triangles = network.find_triangle(1);
        std::sort(triangles.begin(), triangles.end());
        std::vector<int> tmp_assign(glb::M, -1);
        int cnt = 0;
        for (int i = 0; i < (int)triangles.size() and cnt < glb::K * 2 / 3; i++)
        {

            Triangle tri = triangles[i];
            int s = tri.s;
            int t = tri.t;
            int u = tri.u;
            if (tmp_assign[s] == -1 and tmp_assign[t] == -1 and tmp_assign[u] == -1)
            {
                tmp_assign[s] = 1;
                tmp_assign[t] = 2;
                tmp_assign[u] = 3;
                cnt++;
            }
        }
        std::vector<int> cur_assigns = this->assign_at_random2(tmp_assign);
        double cur_score = this->calc_score(cur_assigns);
        this->best_result.assign = cur_assigns;
        this->best_result.eval_score = cur_score;
        std::cerr << cur_score << std::endl;
        while (get_time() < 20.0)
        {
            cur_assigns = this->assign_at_random2(tmp_assign);
            cur_score = this->calc_score(cur_assigns);
            std::cerr << cur_score << std::endl;
            if (cur_score < this->best_result.eval_score)
            {
                this->best_result.assign = cur_assigns;
                this->best_result.eval_score = cur_score;
            }
        }
    }

    void solve()
    {

        network.step(1);
        int tmp = -1;
        std::vector<Node> nodes;
        for (int i = 0; i < (int)this->network.nodes.size(); i++)
        {
            int nid = this->network.nodes[i];
            if (tmp == -1 and nid >= glb::N)
            {
                tmp = i;
            }
            if (nid >= glb::N)
            {
                nodes.push_back(this->network.node_map[nid]);
            }
        }
        std::vector<int> buf;
        std::sample(
            this->network.nodes.begin(),
            this->network.nodes.begin() + tmp,
            std::back_inserter(buf),
            std::min(200, tmp),
            mt);
        for (int nid : buf)
        {
            nodes.push_back(this->network.node_map[nid]);
        }

        std::cerr << "Number of reduced Node: " << nodes.size() << std::endl;

        // define representative
        std::vector<int> repr;
        for (Node node : nodes)
        {
            int nid = node.id;
            if (node.children.size() == 0)
            {
                repr.push_back(nid);
            }
            else
            {
                // find node which is cloest to nodes' center.
                double x = 0., y = 0.;
                double r2max = 1e18;
                int candidate = -1;
                for (int chid : node.children)
                {
                    x += (double)glb::X[chid];
                    y += (double)glb::Y[chid];
                }
                x /= (double)node.children.size();
                y /= (double)node.children.size();
                for (int chid : node.children)
                {
                    double r2 = pow2(glb::X[chid]) + pow2(glb::Y[chid]);
                    if (r2 < r2max)
                    {
                        candidate = chid;
                        r2max = r2;
                    }
                }
                repr.push_back(candidate);
            }
        }

        for (int nid : repr)
            glb::dist[nid] = std::get<0>(dijkstra(nid));

        std::vector<int> cur_assign = this->assign_at_random(nodes);
        long long cur_score = this->eval_score(cur_assign, repr, nodes);
        Result cur_result(cur_assign, cur_score);
        this->best_result = cur_result;

        std::cerr << "INITIAL EVAL SCORE: " << cur_score << std::endl;

        while (get_time() < 5.0)
        {
            cur_assign = this->assign_at_random(nodes);
            cur_score = this->eval_score(cur_assign, repr, nodes);
            std::cerr << "CURRENT SCORE: " << cur_score << std::endl;
            if (cur_score < best_result.eval_score)
            {
                best_result.assign = cur_assign;
                best_result.eval_score = cur_score;
                std::cerr << "HIGHEST UPDATED" << std::endl;
            }
        }
    }

    // consider edges pair. Do not have to consider distances many times.
    // Instead, evaluate edges pair combination score.
    // calc (e1, e2) pair scores, then use these for SA.
    void solve3()
    {
        std::vector<std::vector<long long>> edge_comb_score(glb::M, std::vector<long long>(glb::M, 0LL));
        std::vector<std::vector<int>> edge_comb_list(glb::M);
        std::vector<int> tmp_assign(glb::M, -1);

        for (int e1 = 0; e1 < glb::M; e1++)
        {
            int s = glb::U[e1];
            int t = glb::V[e1];
            std::vector<int> ss, ts;
            ss.push_back(s);
            ts.push_back(t);
            for (pii p : glb::G[s])
            {
                int nid = p.first;
                if (nid == t)
                    continue;
                ss.push_back(nid);
            }
            for (pii p : glb::G[t])
            {
                int nid = p.first;
                if (nid == s)
                    continue;
                ts.push_back(nid);
            }
            for (int s : ss)
            {
                for (int t : ts)
                {
                    if (s == t)
                        continue;
                    std::tuple<int, std::vector<int>> shortest2 = dijkstra(s, t, e1);
                    int dist2 = std::get<0>(shortest2);
                    std::vector<int> pre2 = std::get<1>(shortest2);
                    std::vector<int> edges_on_2nd;
                    int v = t;
                    while (pre2[v] != -1)
                    {
                        int eid = glb::edge_id_matrix[v][pre2[v]];
                        edges_on_2nd.push_back(eid);
                        v = pre2[v];
                    }

                    for (int e2 : edges_on_2nd)
                    {
                        std::tuple<int, std::vector<int>> shortest3 = dijkstra(s, t, e1, e2);
                        int dist3 = std::get<0>(shortest3);
                        // std::vector<int> pre3 = std::get<1>(shortest3);
                        if (edge_comb_score[e1][e2] == 0)
                            edge_comb_list[e1].push_back(e2);
                        edge_comb_score[e1][e2] += dist3 - dist2; //
                    }
                }
            }
        }

        std::vector<int> points = this->get_points_on_peri();
        for (int i = 0; i < (int)points.size(); i++)
        {
            for (int j = i + 1; j < (int)points.size(); j++)
            {
                std::tuple<int, std::vector<int>> shortest1 = dijkstra(i, j);
                int dist1 = std::get<0>(shortest1);
                std::vector<int> pre1 = std::get<1>(shortest1);
                std::vector<int> edges_on_1st;
                int v = j;
                while (pre1[v] != -1)
                {
                    int eid = glb::edge_id_matrix[v][pre1[v]];
                    edges_on_1st.push_back(eid);
                    v = pre1[v];
                }
                for (int eid1 : edges_on_1st)
                {
                    std::tuple<int, std::vector<int>> shortest2 = dijkstra(i, j, eid1);
                    int dist2 = std::get<0>(shortest2);
                    std::vector<int> pre2 = std::get<1>(shortest2);
                    std::vector<int> edges_on_2nd;
                    int vv = j;
                    while (pre2[vv] != -1)
                    {
                        int eid = glb::edge_id_matrix[vv][pre2[vv]];
                        edges_on_2nd.push_back(eid);
                        vv = pre2[vv];
                    }
                    for (int eid2 : edges_on_2nd)
                    {
                        std::tuple<int, std::vector<int>> shortest3 = dijkstra(i, j, eid1, eid2);
                        int dist3 = std::get<0>(shortest3);
                        if (edge_comb_score[eid1][eid2] == 0)
                            edge_comb_list[eid1].push_back(eid2);
                        edge_comb_score[eid1][eid2] += dist3 - dist2;
                    }
                }
            }
        }

        std::vector<int> init_assign = this->assign_at_random();
        long long init_score = this->eval_edge_comb_score(init_assign, edge_comb_score, edge_comb_list);
        this->best_result = {init_assign, (double)init_score};
        std::cerr << "INITIAL SCORE: " << init_score << "!!!!" << std::endl;
        std::vector<int> cur_assign = init_assign;
        long long cur_score = init_score;
        while (get_time() < params.TIME_LIMIT_95)
        {
            if (mt() % 100 < 5)
            {
                cur_assign = this->revise_assign_at_random(cur_assign);
                cur_score = this->eval_edge_comb_score(cur_assign, edge_comb_score, edge_comb_list);
            }
            else
            {
                cur_assign = this->revise_greedy_assign(cur_assign, edge_comb_score, edge_comb_list);
                cur_score = this->eval_edge_comb_score(cur_assign, edge_comb_score, edge_comb_list);
            }

            if ((double)cur_score < this->best_result.eval_score)
            {
                std::cerr << "HIGHEST UPDATED!!!" << std::endl;
                std::cerr << this->best_result.eval_score << " ->" << cur_score << std::endl;
                this->best_result = {cur_assign, cur_score};
            }
        }
    }

    // too heavy, not good.
    void solve4()
    {
        this->dist_init();
        int center_id = -1;
        long long r2min = 1e18;
        for (int i = 0; i < glb::N; i++)
        {
            int x = glb::X[i];
            int y = glb::Y[i];
            long long r2 = pow2(500 - x) + pow2(500 - y);
            if (r2 < r2min)
            {
                r2min = r2;
                center_id = i;
            }
        }
        int start = center_id;
        std::vector<bool> checked(glb::M, false);
        Vec node_list;
        while (true)
        {
            node_list.push_back(start);
            std::tuple<Vec, Vec> shortest = dijkstra(start);
            Vec pre = std::get<1>(shortest);
            for (int i = 0; i < glb::N; i++)
            {
                if (pre[i] == -1 or pre[i] == i)
                    continue;
                int u = i;
                int v = pre[i];
                int eid = glb::edge_id_matrix[u][v];
                checked[eid] = true;
            }

            bool flag = true;
            for (int i = 0; i < glb::M; i++)
            {
                if (!checked[i])
                {
                    start = glb::U[i];
                    flag = false;
                    break;
                }
            }

            if (flag)
                break;
        }
        std::cerr << "CHOOSED NODE NUMBER: " << (int)node_list.size() << std::endl;
        Vec cur_assign = this->assign_at_random();
        long long cur_score = this->calc_score(cur_assign, node_list);
        this->best_result = {cur_assign, (double)cur_score};
        std::cerr << "INITIAL SCORE: " << cur_score << std::endl;
        while (get_time() < this->params.TIME_LIMIT_95)
        {
            Vec next_assign;
            if (mt() % 100 < 95)
            {
                next_assign = this->revise_greedy_assign(cur_assign, node_list);
            }
            else
            {
                next_assign = this->revise_assign_at_random(cur_assign);
            }
            long long next_score = this->calc_score(next_assign, node_list);
            if ((double)next_score < this->best_result.eval_score)
            {
                std::cerr << "HIGHEST UPDATED: " << this->best_result.eval_score << "->" << next_score << std::endl;
                this->best_result = {next_assign, next_score};
            }
            else
            {
                if (mt() % 10 < 9)
                {
                    cur_assign = next_assign;
                }
            }
        }
    }

    // similart to solution 3.
    // In addition, considering (e1, e2, e3) pair contribution
    void solve5()
    {
        std::vector<std::vector<long long>> edge_comb_score(glb::M, std::vector<long long>(glb::M, 0LL));
        std::map<int, std::map<int, std::map<int, long long>>> edge_tri_score;

        std::vector<std::vector<int>> edge_comb_list(glb::M);
        std::map<int, std::map<int, Vec>> edge_tri_list;

        for (int e1 = 0; e1 < glb::M; e1++)
        {
            int s = glb::U[e1];
            int t = glb::V[e1];
            Vec ss, ts;
            ss.push_back(s);
            ts.push_back(t);
            for (pii p : glb::G[s])
            {
                int nid = p.first;
                if (nid == t)
                    continue;
                ss.push_back(nid);
            }
            for (pii p : glb::G[t])
            {
                int nid = p.first;
                if (nid == s)
                    continue;
                ts.push_back(nid);
            }
            for (int s : ss)
            {
                for (int t : ts)
                {
                    if (s == t)
                        continue;
                    std::tuple<int, std::vector<int>> shortest2 = dijkstra(s, t, e1);
                    int dist2 = std::get<0>(shortest2);
                    std::vector<int> pre2 = std::get<1>(shortest2);
                    std::vector<int> edges_on_2nd;
                    int v = t;
                    while (pre2[v] != -1)
                    {
                        int eid = glb::edge_id_matrix[v][pre2[v]];
                        edges_on_2nd.push_back(eid);
                        v = pre2[v];
                    }

                    for (int e2 : edges_on_2nd)
                    {
                        std::tuple<int, std::vector<int>> shortest3 = dijkstra(s, t, e1, e2);
                        int dist3 = std::get<0>(shortest3);
                        std::vector<int> pre3 = std::get<1>(shortest3);
                        if (edge_comb_score[e1][e2] == 0)
                            edge_comb_list[e1].push_back(e2);
                        edge_comb_score[e1][e2] += dist3 - dist2;

                        Vec edges_on_3rd;
                        int vv = t;
                        while (pre3[vv] != -1)
                        {
                            int eid = glb::edge_id_matrix[vv][pre3[vv]];
                            edges_on_3rd.push_back(eid);
                            vv = pre3[vv];
                        }

                        for (int e3 : edges_on_3rd)
                        {
                            std::tuple<int, std::vector<int>> shortest4 = dijkstra(s, t, e1, e2, e3);
                            int dist4 = std::get<0>(shortest4);
                            // std::vector<int> pre4 = std::get<1>(shortest4);
                            if (edge_tri_score[e1][e2][e3] == 0)
                                edge_tri_list[e1][e2].push_back(e3);
                            edge_tri_score[e1][e2][e3] += dist4 - dist3;
                            cerr << dist2 << " " << dist3 << " " << dist4 << endl;
                        }
                    }
                }
            }

        }

        // calculate (e1, e2) pair score near peri
        std::vector<int> points = this->get_points_on_peri();
        for (int i = 0; i < (int)points.size(); i++)
        {
            for (int j = i + 1; j < (int)points.size(); j++)
            {
                std::tuple<int, std::vector<int>> shortest1 = dijkstra(i, j);
                int dist1 = std::get<0>(shortest1);
                std::vector<int> pre1 = std::get<1>(shortest1);
                std::vector<int> edges_on_1st;
                int v = j;
                while (pre1[v] != -1)
                {
                    int eid = glb::edge_id_matrix[v][pre1[v]];
                    edges_on_1st.push_back(eid);
                    v = pre1[v];
                }
                for (int eid1 : edges_on_1st)
                {
                    std::tuple<int, std::vector<int>> shortest2 = dijkstra(i, j, eid1);
                    int dist2 = std::get<0>(shortest2);
                    std::vector<int> pre2 = std::get<1>(shortest2);
                    std::vector<int> edges_on_2nd;
                    int vv = j;
                    while (pre2[vv] != -1)
                    {
                        int eid = glb::edge_id_matrix[vv][pre2[vv]];
                        edges_on_2nd.push_back(eid);
                        vv = pre2[vv];
                    }
                    for (int eid2 : edges_on_2nd)
                    {
                        std::tuple<int, std::vector<int>> shortest3 = dijkstra(i, j, eid1, eid2);
                        int dist3 = std::get<0>(shortest3);
                        Vec pre3 = std::get<1>(shortest3);
                        if (edge_comb_score[eid1][eid2] == 0)
                            edge_comb_list[eid1].push_back(eid2);
                        edge_comb_score[eid1][eid2] += dist3 - dist2;

                        Vec edges_on_3rd;
                        int vv = j;
                        while (pre3[v] != -1)
                        {
                            int eid = glb::edge_id_matrix[v][pre3[v]];
                            edges_on_3rd.push_back(eid);
                            v = pre3[v];
                        }

                        for (int eid3 : edges_on_3rd)
                        {
                            std::tuple<int, std::vector<int>> shortest4 = dijkstra(i, j, eid1, eid2, eid3);
                            int dist4 = std::get<0>(shortest4);
                            // std::vector<int> pre4 = std::get<1>(shortest4);
                            if (edge_tri_score[eid1][eid2][eid3] == 0)
                                edge_tri_list[eid1][eid2].push_back(eid3);
                            edge_tri_score[eid1][eid2][eid3] += dist4 - dist3;
                        }
                    }
                }
            }
        }

        // assign section and acending
        std::vector<int> init_assign = this->assign_by_edge_score(edge_comb_score, edge_comb_list);
        long long init_score = this->eval_edge_comb_score(init_assign, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
        for (int i = 0; i < 30; i++)
        {
            Vec tmp_assign = this->assign_at_random();
            long long tmp_score = this->eval_edge_comb_score(tmp_assign, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
            if (tmp_score < init_score)
            {
                init_assign = tmp_assign;
                init_score = tmp_score;
            }
        }
        this->best_result = {init_assign, (double)init_score};
        std::cerr << "INITIAL SCORE: " << init_score << "!!!!" << std::endl;
        std::vector<int> cur_assign = init_assign;
        long long cur_score = init_score;
        while (get_time() < params.TIME_LIMIT_95)
        {
            // for(int a: cur_assign) cerr << a << " ";
            // cerr << endl;
            if (mt() % 100 < 5)
            {
                cur_assign = this->revise_assign_at_random(cur_assign);
                cur_score = this->eval_edge_comb_score(cur_assign, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
            }
            else
            {
                cur_assign = this->revise_greedy_assign(cur_assign, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
                cur_score = this->eval_edge_comb_score(cur_assign, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
            }

            if ((double)cur_score < this->best_result.eval_score)
            {
                std::cerr << "HIGHEST UPDATED!!!" << std::endl;
                std::cerr << this->best_result.eval_score << " ->" << cur_score << std::endl;
                this->best_result = {cur_assign, cur_score};
            }
        }
    }

    void dist_init()
    {
        for (int i = 0; i < glb::N; i++)
        {
            glb::dist[i] = std::get<0>(dijkstra(i));
        }
    }

    std::vector<int> get_points_on_peri()
    {
        int center_nid = -1;
        double r2max = 1e18;
        std::vector<int> ret(8);
        std::vector<double> r2max_list(8, 1e18);
        for (int nid = 0; nid < glb::N; nid++)
        {
            int x = glb::X[nid];
            int y = glb::Y[nid];
            double r2 = pow2(500 - x) + pow2(500 - y);
            if (r2 < r2max)
            {
                r2max = r2;
                center_nid = nid;
            }
            for (int j = 0; j < 8; j++)
            {
                double xx = glb::R * std::cos(2 * j / 8 * PI);
                double yy = glb::R * std::sin(2 * j / 8 * PI);
                double r2_peri = pow2(xx - x) + pow2(yy - y);
                if (r2_peri < r2max_list[j])
                {
                    r2max_list[j] = r2_peri;
                    ret[j] = nid;
                }
            }
        }
        ret.push_back(center_nid);
        return ret;
    }

    Vec assign_by_edge_score(const std::vector<std::vector<long long>> &edge_comb_score, const std::vector<Vec> &edge_comb_list)
    {
        Vec ret(glb::M, -1), color(glb::M, -1);
        Vec remain(glb::D, glb::K);

        for (int eid1 = 0; eid1 < glb::M; eid1++)
        {
            long long score = 0;
            int target_eid;
            for (int eid2 : edge_comb_list[eid1])
            {
                if (edge_comb_score[eid1][eid2] > score)
                {
                    target_eid = eid2;
                    score = edge_comb_score[eid1][eid2];
                }
            }

            int maxIter = 20;
            int day = mt() % glb::D + 1;
            bool flag = false;
            while (maxIter--)
            {
                day = mt() % glb::D + 1;
                if (ret[target_eid] == day or remain[day - 1] == 0)
                    continue;
                ret[eid1] = day;
                remain[day - 1]--;
                flag = true;
                break;
            }
            if (!flag)
            {
                Vec cum(glb::D + 1, 0);
                for (int i = 0; i < glb::D; i++)
                    cum[i + 1] = cum[i] + remain[i];
                int v = mt() % cum[glb::D];
                int next_day = std::upper_bound(cum.begin(), cum.end(), v) - cum.begin();
                ret[eid1] = next_day;
                remain[next_day - 1]--;
            }
        }
        return ret;
    }

    std::vector<int> revise_assign_at_random(const std::vector<int> &assign)
    {
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> ret = assign;

        for (int day : assign)
            remain[day - 1]--;
        if (mt() % 10 < 5)
        {
            // change e1 day to another
            int eid = mt() % glb::M;
            int next_day = mt() % glb::D + 1;
            while (next_day == assign[eid] or remain[next_day - 1] == 0)
            {
                next_day = mt() % glb::D + 1;
            }
            ret[eid] = next_day;
        }
        else
        {
            // swap (e1, e2) assign day.
            int maxIter = 100;
            while (maxIter--)
            {
                int e1 = mt() % glb::M;
                int e2 = mt() % glb::M;
                if (e1 == e2 or assign[e1] == assign[e2])
                    continue;
                std::swap(ret[e1], ret[e2]);
                break;
            }
        }
        return ret;
    }

    std::vector<int> revise_greedy_assign(const std::vector<int> &assign, const std::vector<std::vector<long long>> &edge_comb_score, const std::vector<std::vector<int>> &edge_comb_list)
    {
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> ret = assign;
        for (int day : assign)
            remain[day - 1]--;

        int target_eid = mt() % glb::M;
        int best_day = ret[target_eid];
        long long best_score = LINF;

        for (int d = 1; d <= glb::D; d++)
        {
            ret[target_eid] = d;
            long long cur_score = this->eval_edge_comb_score(ret, edge_comb_score, edge_comb_list);
            if (cur_score < best_score and remain[d - 1] > 0)
            {
                best_day = d;
                best_score = cur_score;
            }
        }
        ret[target_eid] = best_day;
        return ret;
    }

    Vec revise_greedy_assign(const std::vector<int> &assign, const std::vector<std::vector<long long>> &edge_comb_score, const std::vector<std::vector<int>> &edge_comb_list, std::map<int, std::map<int, std::map<int, long long>>> &edge_tri_score, std::map<int, std::map<int, Vec>> &edge_tri_list)
    {
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> ret = assign;
        for (int day : assign)
            remain[day - 1]--;

        // one edge
        int target_eid = mt() % glb::M;

        if ((mt() % 10) < 9)
        {

            remain[ret[target_eid] - 1]++;

            int best_day = ret[target_eid];
            long long best_score = LINF;

            for (int d = 1; d <= glb::D; d++)
            {
                ret[target_eid] = d;
                long long cur_score = this->eval_edge_comb_score(ret, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
                if (cur_score < best_score and remain[d - 1] > 0)
                {
                    best_day = d;
                    best_score = cur_score;
                }
            }
            ret[target_eid] = best_day;
        }
        else
        {
            // two edge
            Vec edge_list;
            int u = glb::U[target_eid];
            int v = glb::V[target_eid];
            for (pii p : glb::G[u])
            {
                edge_list.push_back(p.first);
            }
            for (pii p : glb::G[v])
            {
                edge_list.push_back(p.first);
            }
            int idx = mt() % (unsigned int)edge_list.size();
            int another_edge_id = edge_list[idx];

            // add remain for reassign
            remain[ret[target_eid] - 1]++;
            remain[ret[another_edge_id] - 1]++;

            long long best_score = LINF;
            pii best_pair;
            for (int d1 = 1; d1 <= glb::D; d1++)
            {
                for (int d2 = 1; d2 <= glb::D; d2++)
                {
                    ret[target_eid] = d1;
                    ret[another_edge_id] = d2;
                    long long cur_score = this->eval_edge_comb_score(ret, edge_comb_score, edge_comb_list, edge_tri_score, edge_tri_list);
                    if (cur_score < best_score and remain[d1 - 1] > 0 and remain[d2 - 1] > 0)
                    {
                        if (d1 == d2 and remain[d1 - 1] <= 1)
                            continue;
                        best_pair = {d1, d2};
                        best_score = cur_score;
                    }
                }
            }
            ret[target_eid] = best_pair.first;
            ret[another_edge_id] = best_pair.second;
        }
        return ret;
    }

    // functions for solve4
    std::vector<int> revise_greedy_assign(const Vec &assign, const Vec &node_list)
    {
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> ret = assign;
        for (int day : assign)
            remain[day - 1]--;

        int target_eid = mt() % glb::M;
        int best_day = ret[target_eid];
        long long best_score = LINF;

        for (int d = 1; d <= glb::D; d++)
        {
            ret[target_eid] = d;
            long long cur_score = this->calc_score(ret, node_list);
            if (cur_score < best_score and remain[d - 1] > 0)
            {
                best_day = d;
                best_score = cur_score;
            }
        }
        ret[target_eid] = best_day;
        return ret;
    }

    std::vector<int> assign_at_random()
    {
        std::vector<int> ret(glb::M, -1);
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> cum(glb::D + 1, 0);
        for (int i = 0; i < glb::M; i++)
        {
            if (ret[i] == -1)
            {
                for (int j = 0; j < glb::D; j++)
                    cum[j + 1] = cum[j] + remain[j];
                int M = cum[glb::D];
                int v = mt() % M;
                int day = std::upper_bound(cum.begin(), cum.end(), v) - cum.begin();
                ret[i] = day;
                remain[day - 1]--;
            }
        }
        return ret;
    }

    std::vector<int> assign_at_random(const std::vector<Node> &nodes)
    {
        std::vector<int> ret(glb::M, -1);
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> cum(glb::D + 1, 0);

        // 1. define edges inside big node.
        for (Node node : nodes)
        {
            std::vector<int> days;
            for (int i = 1; i <= glb::D; i++)
                days.push_back(i);
            std::shuffle(days.begin(), days.end(), mt);
            for (int eid : node.edges)
            {
                int day = days.back();
                days.pop_back();
                ret[eid] = day;
                remain[day - 1]--;
            }
        }

        // 2. define remain edges.
        for (int i = 0; i < glb::M; i++)
        {
            if (ret[i] == -1)
            {
                for (int j = 0; j < glb::D; j++)
                    cum[j + 1] = cum[j] + remain[j];
                int M = cum[glb::D];
                int v = mt() % M;
                int day = std::upper_bound(cum.begin(), cum.end(), v) - cum.begin();
                ret[i] = day;
                remain[day - 1]--;
            }
        }

        return ret;
    }

    std::vector<int> assign_at_random2(const std::vector<int> &tmp_assign)
    {
        std::vector<int> ret = tmp_assign;
        std::vector<int> remain(glb::D, glb::K);
        std::vector<int> cum(glb::D + 1, 0);
        for (int tassign : tmp_assign)
        {
            if (tassign != -1)
            {
                remain[tassign - 1]--;
            }
        }
        for (int i = 0; i < glb::M; i++)
        {
            if (ret[i] == -1)
            {
                for (int j = 0; j < glb::D; j++)
                    cum[j + 1] = cum[j] + remain[j];
                int M = cum[glb::D];
                int v = mt() % M;
                int day = std::upper_bound(cum.begin(), cum.end(), v) - cum.begin();
                ret[i] = day;
                remain[day - 1]--;
            }
        }
        return ret;
    }

    ll calc_score(const std::vector<int> &assign)
    {
        ll ret = 0;
        for (int d = 1; d <= glb::D; d++)
        {
            std::vector<std::vector<bool>> deleted(glb::N, std::vector<bool>(glb::N, false));
            std::vector<std::vector<int>> dist(glb::N);
            for (int eid = 0; eid < glb::M; eid++)
            {
                if (assign[eid] == d)
                {
                    int u = glb::U[eid];
                    int v = glb::V[eid];
                    deleted[u][v] = deleted[v][u] = true;
                }
            }
            ll f = 0;
            for (int nid = 0; nid < glb::N; nid++)
            {
                dist[nid] = dijkstra(nid, deleted);
                for (int nx = 0; nx < glb::N; nx++)
                {
                    ret += dist[nid][nx] - glb::dist[nid][nx];
                }
            }
        }
        ret = std::round((double)ret / glb::N / (glb::N - 1) / glb::D * 1e3);

        return ret;
    }

    ll calc_score(const Vec &assign, const Vec &node_list)
    {
        ll ret = 0;
        for (int d = 1; d <= glb::D; d++)
        {
            std::vector<std::vector<bool>> deleted(glb::N, std::vector<bool>(glb::N, false));
            std::vector<std::vector<int>> dist(glb::N);
            for (int eid = 0; eid < glb::M; eid++)
            {
                if (assign[eid] == d)
                {
                    int u = glb::U[eid];
                    int v = glb::V[eid];
                    deleted[u][v] = deleted[v][u] = true;
                }
            }
            ll f = 0;
            for (int nid : node_list)
            {
                dist[nid] = dijkstra(nid, deleted);
                for (int nx = 0; nx < glb::N; nx++)
                {
                    ret += dist[nid][nx] - glb::dist[nid][nx];
                }
            }
        }
        ret = std::round((double)ret / glb::N / (glb::N - 1) / glb::D * 1e3);

        return ret;
    }

    double eval_score(const std::vector<int> &assign, const std::vector<int> &repr, const std::vector<Node> &nodes)
    {
        long long ret = 0.0;
        int total_node = 0;
        for (Node node : nodes)
        {
            if (node.children.size() > 0)
                total_node += (int)node.children.size();
            else
                total_node++;
        }

        // day 1 to D score
        for (int d = 1; d <= glb::D; d++)
        {
            std::vector<std::vector<bool>> deleted(glb::N, std::vector<bool>(glb::N, false));
            for (int idx = 0; idx < glb::N; idx++)
            {
                if (assign[idx] == d)
                {
                    int u = glb::U[idx];
                    int v = glb::V[idx];
                    deleted[u][v] = deleted[v][u] = true;
                }
            }
            for (int r = 0; r < (int)repr.size(); r++)
            {
                int rid = repr[r];
                std::vector<int> dist = dijkstra(rid, deleted);
                double f = 0.0;
                for (int nx = 0; nx < glb::N; nx++)
                {
                    f += dist[nx] - glb::dist[rid][nx];
                }
                if (nodes[r].children.size() > 0)
                    f *= (int)nodes[r].children.size();

                ret += f;
            }
        }

        ret /= total_node * (total_node - 1) * glb::D;
        ret *= 1000.;
        return ret;
    }

    long long eval_edge_comb_score(const std::vector<int> &assign, const std::vector<std::vector<long long>> &edge_comb_score, const std::vector<std::vector<int>> &edge_comb_list)
    {
        long long ret = 0;
        for (int e1 = 0; e1 < glb::M; e1++)
        {
            for (int e2 : edge_comb_list[e1])
            {
                if (e1 == e2 or assign[e1] != assign[e2])
                    continue;
                ret += edge_comb_score[e1][e2];
            }
        }
        return ret;
    }

    long long eval_edge_comb_score(const Vec &assign, const std::vector<std::vector<long long>> &edge_comb_score, const std::vector<Vec> &edge_comb_list, std::map<int, std::map<int, std::map<int, long long>>> &edge_tri_score, std::map<int, std::map<int, Vec>> &edge_tri_list)
    {
        long long ret = 0;
        for (int e1 = 0; e1 < glb::M; e1++)
        {
            for (int e2 : edge_comb_list[e1])
            {
                if (e1 == e2 or assign[e1] != assign[e2])
                    continue;
                ret += edge_comb_score[e1][e2];

                for (int e3 : edge_tri_list[e1][e2])
                {
                    if (assign[e1] == assign[e2] == assign[e3])
                    {
                        ret += edge_tri_score[e1][e2][e3];
                    }
                }
            }
        }
        return ret;
    }

    void summary()
    {
        std::cerr << "###################### SUMMARY ######################" << std::endl;
        std::cerr << "ELAPSED TIME: " << ((double)clock() - start_time) / CLOCKS_PER_SEC << " s" << std::endl;
        std::cerr << "BRST EVAL SCORE: " << this->best_result.eval_score << std::endl;
        std::cerr << "#####################################################" << std::endl;
    }

    void init()
    {
    }

    void debug()
    {
        std::map<int, int> edge_cnt;
        int cnt = 0;
        int num_child = 0;
        int total_node = 0;
        std::vector<bool> found_edge(glb::M, false);

        for (auto p : this->network.node_map)
        {
            Node node = p.second;
            num_child += node.children.size();
            for (int eid : node.edges)
            {
                edge_cnt[eid]++;
                found_edge[eid] = true;
                cnt++;
                if (edge_cnt[eid] > 1)
                {
                    std::cerr << "Invalid number of edges." << std::endl;
                    return;
                }
            }
        }
        for (auto p : this->network.edges)
        {
            int u = p.first;
            for (auto p2 : p.second)
            {
                int v = p2.first;
                if (u > v)
                    continue;
                for (int eid : p2.second)
                {
                    edge_cnt[eid]++;
                    found_edge[eid] = true;
                    cnt++;
                    if (edge_cnt[eid] > 1)
                    {
                        std::cerr << "Invalid number of edges." << std::endl;
                    }
                }
            }
        }

        auto dfs = [&](auto dfs, int nid, int &c) -> void
        {
            if (this->network.node_map[nid].id < glb::N)
            {
                c++;
                return;
            }
            for (int chid : this->network.node_map[nid].children)
            {
                dfs(dfs, chid, c);
            }
        };
        for (int nid : this->network.nodes)
        {
            dfs(dfs, nid, total_node);
        }

        std::cerr << "M: " << glb::M << ", cnt: " << cnt << std::endl;
        std::cerr << "N: " << glb::N << ", total node: " << total_node << std::endl;
        std::cerr << "Number of child: " << num_child << std::endl;
        for (int i = 0; i < glb::M; i++)
        {
            if (found_edge[i] == false)
            {
                int u = glb::U[i];
                int v = glb::V[i];
                int pu = this->network.get_leader(u);
                int pv = this->network.get_leader(v);
                std::cerr << "Edge " << i << " is missing." << std::endl;
                std::cerr << "(U,V):" << u << ", " << v << ", parent: (" << pu << "," << pv << ")" << std::endl;
            }
        }
    }

    void show_best()
    {
        for (int i = 0; i < glb::M; i++)
        {
            std::cout << this->best_result.assign[i] << std::endl;
        }
    }

    void show_info()
    {
        printf("%d %d %d %d\n", glb::N, glb::M, glb::D, glb::K);
        for (int nid : network.nodes)
        {
            std::cout << "NID: " << nid << ", Parent: " << this->network.node_map[nid].parent << std::endl;
            for (int eid : this->network.node_map[nid].edges)
            {
                std::cout << eid << " ";
            }
            std::cout << std::endl;
            std::cout << "Edges:" << std::endl;
            for (auto p : this->network.edges[nid])
            {
                std::cout << p.first << " ";
            }
            std::cout << std::endl;
        }
    }
};

int main(int argv, char **argc)
{
    start_time = clock();
    Solver solver;
    solver.solve5();

#ifdef _TEST
    solver.summary();
#endif

#ifdef _DEBUG
    solver.debug();
#endif

#ifdef _PYTHON
    solver.solve();
    solver.network.show_nodes_pos();
    solver.network.show_edge_conn();
#else
    solver.show_best();
#endif
    return 0;
}