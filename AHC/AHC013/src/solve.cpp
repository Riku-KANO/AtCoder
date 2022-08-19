#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<queue>
#include<algorithm>
#include<set>
#include<cmath>
#include<chrono>
#include<random>
#include<iomanip>
#include<fstream>
#include<cassert>
#include<iterator>
#include<functional>
#include<utility>
#include<atcoder/dsu>
// #define DEBUG
// #define TEST

#define rep(i, n) for(int i = 0; i < (n); i++)
#define pii pair<int,int>
#define size_dir pair<int,int> // group sizeとdirのペア
#define leaderNum map<int, size_dir> //union findのleaderをkeyとしてグループサイズとdirectionのpairをvalueとして返す。
using namespace std;
using namespace atcoder;
using ll = long long;

struct PC {
    int i, j;
    char label;
    int id=-1;
    int adj[4]={0,0,0,0}; // 1: PCと隣接している。 0: PCと隣接していない。
    bool operator==(const PC&rhs){
        if(rhs.i == i and rhs.j == j)return true;
        else return false;
    };
};


struct Move {
    int ci, cj, ni, nj;
    Move(int _ci, int _cj, int _ni, int _nj) : ci{_ci}, cj{_cj}, ni{_ni}, nj{_n}{};
};

struct Connect {
    int ci, cj, ni, nj;
    Connect(int _ci, int _cj, int _ni, int _nj) : ci{_ci}, cj{_cj}, ni{_ni}, nj{_nj}{};
};

struct UnionFind {
    map<pair<int,int>, pair<int, int>> parent;
    UnionFind() :parent() {}

    pair<int, int> find(pair<int, int> x){
        if (parent.find(x) == parent.end()) {
            parent[x] = x;
            return x;
        } else if (parent[x] == x) {
            return x;
        } else {
            parent[x] = find(parent[x]);
            return parent[x];
        }
    }

    void unite(pair<int, int> x, pair<int, int> y){
        x = find(x);
        y = find(y);
        if (x != y) {
            parent[x] = y;
        }
    }

    bool same(pair<int, int> x, pair<int, int> y) {
        if(parent[x] == parent[y]) return true;
        else return false;
    }
};

const int dx[4] = {1, 0, -1, 0};
const int dy[4] = {0, 1, 0, -1};
const int dx8[8]={-1, -1, -1, 0, 1, 1, 1, 0};
const int dy8[8]={-1, 0, 1, -1, -1, 0, 1, 1};
const int INF = 1<<30;


// global variable//////
double TL;
int N, K;
string S_ini[50];
string S[50];
double PROPAGATION;

vector<vector<PC>> pc_pos(50, vector<PC>(50));
vector<Move> move_list;
vector<Connect> connect_list;
vector<PC> movable_pcs;

PC* pcid[5][100];
////////////////////////

class Config {
    public:
        double tl=3.0;
        double propagation = 0.6;
        string file_name="../output/config.txt";
        Config() {}
        void init() {
            TL = this->tl;
            PROPAGATION = this->propagation;
            return ;
        }

        void save() {
            ofstream ofs;
            ofs.open(this->file_name, ios::out);
            cout << "TL: " << this->tl << "\n";
            cout << "PROPACATION: " << this->propagation << "\n";
            ofs.close();
            return ;
        }
};

random_device seed_gen;
mt19937 mt(seed_gen());

bool out_bound(int ni, int nj, int n) {
    if(ni < 0 || nj < 0 || ni >= n || nj >=n) return true;
    return false;
}

int calc_score(string field[], vector<Connect> &c_list){
    for (auto r : move_list) {
        // assert(field[r.ci][r.cj] != '0');
        // assert(field[r.nj][r.ni] == '0');
        swap(field[r.ci][r.cj], field[r.ni][r.nj]);
    }

    UnionFind uf;
    for (auto r : c_list) {
        pair<int, int> p1(r.ci, r.cj), p2(r.ni, r.nj);
        uf.unite(p1, p2);
    }

    vector<pair<int, int>> computers;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (field[i][j] != '0') {
                computers.emplace_back(i, j);
            }
        }
    }

    int score = 0;
    for (int i = 0; i < (int)computers.size(); i++) {
        for (int j = i+1; j < (int)computers.size(); j++) {
            auto c1 = computers[i];
            auto c2 = computers[j];
            if (uf.find(c1) == uf.find(c2)) {
                score += (field[c1.first][c1.second] == field[c2.first][c2.second]) ? 1 : -1;
            }
        }
    }
    return max(score, 0);
}

void debug() {
    cerr << "### DEBUG ###" << "\n";
    return ;
}

void test() {
    return ;
}

void input() {
    cin >> N >> K;
    rep(i, N) cin >> S_ini[i];
}

void init() {
    cerr << "\n####### INIT START ######" << endl;
    rep(i, N) S[i] = S_ini[i];
    map<char,int> pc_cnt;
    cerr << "init: " << N << " " << K << endl;
    rep(i, N) {
        rep(j, N) {
            int label = S[i][j]-'0'; label--;
            int id = pc_cnt[S[i][j]];
            if(S[i][j]!='0') {
                pc_pos[i][j] = {i, j, S[i][j], id};
                pcid[label][id] = &pc_pos[i][j];
            }
            else pc_pos[i][j] = {i, j, S[i][j], -1};
            pc_cnt[S[i][j]]++;
            for(int dir = 0; dir < 4; ++dir) {
                int ni = i + dy[dir];
                int nj = j + dx[dir];
                if(out_bound(ni, nj, N)) {
                    pc_pos[i][j].adj[dir]=1;
                    continue;
                }
                if(S[ni][nj]!='0') pc_pos[i][j].adj[dir]=1;
            }
        }
    }
    Config config;
    config.init();
    #ifdef TEST
        cerr << "Parameter saving..." << endl;
        config.save();
    #endif
    cerr << "\n#####################" << endl;
    cerr << "N: " << N << ", K: " << K << endl;
    cerr << "TL: " << TL << endl;
    cerr << "PROPAGATION: " << PROPAGATION << endl;
    cerr << "#####################\n\n";
    return;
}


void update_pc_pos(PC& cur, PC& next, vector<PC> &mvl, vector<Move> &move_list) {
    // cur(PCあり) と next(空)の位置のPCの位置を交換(移動)する。
    vector<PC> to_be_erased;
    vector<PC> to_be_added;
    int ci = cur.i; int cj = cur.j;
    int next_i = next.i, next_j = next.j;
    move_list.push_back({ci, cj, next_i, next_j});
    // 最初の位置のPCについて
    for(int dir = 0; dir < 4; ++dir) {
        int ni = ci + dy[dir];
        int nj = cj + dx[dir];
        if(out_bound(ni, nj, N)) continue;
        if(S[ni][nj]!='0')pc_pos[ni][nj].adj[(dir+2)%4]=0;
        int sum = 0;
        for(int d = 0; d < 4; ++d) sum += pc_pos[ni][nj].adj[d];
        if(sum == 3) to_be_added.push_back(pc_pos[ni][nj]);
    }
    // 移動後のPCについて
    for(int dir = 0; dir < 4; ++dir) {
        int ni = next_i + dy[dir];
        int nj = next_j + dx[dir];
        if(out_bound(ni, nj, N)) continue;
        if(S[ni][nj]!='0')pc_pos[ni][nj].adj[(dir+2)%4]=1;
        int sum = 0;
        for(int d = 0; d < 4; ++d)sum += pc_pos[ni][nj].adj[d];
        if(sum == 4) to_be_erased.push_back(pc_pos[ni][nj]);
    }
    #ifdef DEBUG
        cerr << "FROM(i,j): " << ci << " " << cj << ", TO(i,j): " << next_i << " " << next_j << endl;
        cerr << pc_pos[ci][cj].i << " " << pc_pos[ci][cj].j << " " << pc_pos[ci][cj].id << " " << pc_pos[ci][cj].label << endl;
    #endif

    pcid[cur.label-'1'][cur.id] = &pc_pos[next_i][next_j];
    swap(S[ci][cj], S[next_i][next_j]);
    swap(pc_pos[ci][cj], pc_pos[next_i][next_j]);
    pc_pos[ci][cj].i = ci; pc_pos[ci][cj].j = cj;
    pc_pos[next_i][next_j].i = next_i;pc_pos[next_i][next_j].j=next_j;
    cerr << pc_pos[ci][cj].i << " " << pc_pos[ci][cj].j << " " << pc_pos[ci][cj].id <<" " << pc_pos[ci][cj].label << endl;
    // 動かした後のPCがない箇所の位置について
    for(int dir = 0; dir < 4; ++ dir){
        int ni = ci+dy[dir];
        int nj = cj+dx[dir];
        if(out_bound(ni, nj, N)) {
            pc_pos[ci][cj].adj[dir]=0;
            continue;
        }
        if(S[ni][nj]=='0') pc_pos[ci][cj].adj[dir]=1;
        else pc_pos[ci][cj].adj[dir]=0;
    }
    //動かした後のPCがある箇所の位置について
    for(int dir = 0; dir < 4; ++ dir){
        int ni = next_i+dy[dir];
        int nj = next_j+dx[dir];
        if(out_bound(ni, nj, N)) {
            pc_pos[next_i][next_j].adj[dir]=0;
            continue;
        }
        if(S[ni][nj]=='0') pc_pos[next_i][next_j].adj[dir]=1;
        else pc_pos[next_i][next_j].adj[dir]=0;
    }
    for(PC pc: to_be_erased) {
        for(int i = 0; i < mvl.size(); ++i) {
            if(mvl[i].i == pc.i and mvl[i].j == pc.j) {
                // cerr << mvl[i].i << " " << mvl[i].j << endl;
                // cerr << pc.i << " " << pc.j << endl;
                mvl.erase(mvl.begin()+i);
                break;
            }
        }
    }
    for(PC pc: to_be_added) mvl.push_back(pc);
    
    return ;
}

void pc_pos_debug() {
    for(int l = 0; l < K; ++l) {
        for(int i = 0; i < N; ++i) {
            // cerr << pcid[l][i]->i << " " << pcid[l][i]->j << " " << pcid[l][i]->id << " " << pcid[l][i]->label << endl;
            assert(pcid[l][i]->id == i);
            assert(pcid[l][i]->label-'0' == l+1);
        }
    }

    cerr << "FOR position:" <<endl;
    rep(i, N) {
        rep(j, N) {
            assert(pc_pos[i][j].i == i);
            assert(pc_pos[i][j].j == j);
            assert(pc_pos[i][j].label == S[i][j]);
        }
    }
    cerr << "POS DEBUG DONE\n\n";
    return ;
}

// メモ
// 操作の回数に上限があるのでやりたい放題できない。
// AHC011のように最適な状態を前もって計算しておいてそれに近づけた方がよさそうか。
// 100*Kなので移動回数は少ない方がよさそう。
// 同じ種類のパソコンをなるべく多く繋げられればスコアがとても良くなる。
// 同じ種類のパソコンをつなげることができたら100*99/2=4950点獲得できる。
// 細分化(2*2or3*3)して大きなスコアを作ってみるか。
// BFS + dequeでコストを調べようとしたがK=5の時の計算(特にN<30)がきつい。K＝4も少し遅いが許容できそうなレベル。K=2,K=3は余裕。
// 各地点でのポテンシャル的なものを計算してシミュレーション的なことをしてみる(動的計画法?)。
// 最初は１手でつながりそうなものがあれば繋げておいた方が後の操作が容易になる。

void solve1() {
    clock_t start, now;
    start = clock();
    


    #ifdef DEBUG
    debug();
    #elif TEST
    test();
    #endif
}

void align(){
    return ;
}

void trial_bfs() {
    vector<vector<PC>> pcs(K);
    rep(i, N) {
        rep(j, N) {
            if(S[i][j]=='0') continue;
            PC tmp = {i, j, S[i][j], 0};
            pcs[S[i][j]-'1'].push_back(tmp);
        }
    }
    
    for(int l = 1; l <= K; ++l) {
        vector<vector<int>> bfs_state(N, vector<int>(N, -1));
        vector<vector<int>> groups(N, vector<int>(N,-1));
        ll best_score = INF;
        char label = l+'0';
        cerr <<"LABEL: " << label << endl;
        for(int j = 0; j < N; ++j) {
            vector<vector<int>> bfs_tmp_state(N, vector<int>(N, INF));
            vector<vector<int>> groups_tmp(N, vector<int>(N,-1));

            PC start = pcs[l-1][j];
            ll score=0;
            deque<pair<PC, int>> dq;
            dq.push_back({start, 0});
            int gid=0;
            while(!dq.empty()) {
                auto p = dq.front(); dq.pop_front();
                PC now = p.first;
                int d = p.second;
                bfs_tmp_state[now.i][now.j]=d;
                groups[now.i][now.j]=gid;
                for(int dir = 0; dir < 4; dir++) {
                    int k = 0;
                    while(1) {
                        k++;
                        int ni = now.i + dy[dir]*k;
                        int nj = now.j + dx[dir]*k;
                        if(ni < 0 || nj < 0 || ni >= N || nj >= N) break;
                        if(S[ni][nj]=='0')continue;
                        if(S[ni][nj]==label && d<bfs_tmp_state[ni][nj]){
                            dq.push_front((pair<PC,int>){{ni,nj, S[ni][nj]}, d});
                            if(S[ni][nj]==S[now.i][now.j]) groups_tmp[ni][nj] = groups_tmp[now.i][now.j];
                            else groups_tmp[ni][nj]=++gid;
                            break;
                        }
                        else if(S[ni][nj]!=label && d+1<bfs_tmp_state[ni][nj]) {
                            dq.push_back((pair<PC,int>){{ni, nj, S[ni][nj]}, d+1});
                            break;
                        }
                    }
                }
            }
            rep(y,N)rep(x,N) {
                if(bfs_tmp_state[y][x]==INF) continue;
                score += bfs_tmp_state[y][x];
            }
            if(score<best_score) {
                rep(y,N)rep(x,N)bfs_state[y][x]=bfs_tmp_state[y][x];
                best_score=score;
                groups=groups_tmp;
            }
        }
        cout << "## LABEL: " << label  << ", SCORE: " << best_score << endl;
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                if (bfs_state[i][j]==0 and S[i][j]==label) printf("%2s ", "S");
                else if(S[i][j]=='0') printf("%2s ", ".");
                else if(S[i][j]!=label) printf("%2s ", "+");
                else if(bfs_state[i][j]==INF) printf("%2s ", 'X');
                else printf("%2d ", bfs_state[i][j]);
            }
            cout << endl;
        }
        cout << endl;
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                if(S[i][j]=='0') printf("%2s ", ".");
                else if(S[i][j]!=label) printf("%2s ", "+");
                else printf("%2d ", groups[i][j]);
            }
            cout << endl;
        }
    }
}

bool search_and_move(const pair<int, int> dest, const char target, dsu& uf, int source_id) {
    bool halt=false;
    if(move_list.size()>100*K-140){
        halt=true;
        return false;
    }
    cerr << "DESTINATION(i, j): " << dest.first << " " << dest.second << endl; 
    int tol = 10;
    vector< vector< pair<int,pii> > > traj(N, vector<pair<int, pii>>(N, {INF, {-1, -1}}));

    // bfs part
    deque<pii> dq;
    traj[dest.first][dest.second] = {0, {dest}};
    dq.push_front(dest);
    bool flag = true;
    bool found = false;
    pii pos_emp; //何もない場所の座標
    while(!dq.empty() and flag) {
        pii p = dq.front(); dq.pop_front();
        int ci = p.first; 
        int cj = p.second;
        cerr << "SEAECHING...(i,j): " << ci << " " << cj << ", label: " << S[ci][cj] << ", pre is: " << traj[ci][cj].second.first << " "<< traj[ci][cj].second.second <<endl;
        if(traj[ci][cj].first >= tol) continue;
        for(int dir = 0; dir < 4 and flag; ++ dir) {
            int ni = ci + dy[dir];
            int nj = cj + dx[dir];
            if(out_bound(ni, nj, N)) continue;
            if(S[ni][nj] == target) continue;
            if(traj[ni][nj].first > traj[ci][cj].first + 1) {
                traj[ni][nj] = {traj[ci][cj].first + 1, pii{ci, cj}};
                if(S[ni][nj]=='0') {
                    cerr << "FOUND!: " << ni << " " << nj << endl;
                    found = true;
                    flag = false;
                    pos_emp = {ni,nj};
                    break;
                }
                dq.push_back({ni, nj});
            }
        }
    }
    if(!found) return false;
    // traj detail
    pii cur = pos_emp;
    while(1){
        cerr << "(" << cur.first << "," << cur.second << "), ";
        if(cur == dest) break;
        cur = traj[cur.first][cur.second].second;
    }
    cerr << endl;
    cerr << "\n### EMPTY SPACE MOVE TURN ####\n";
    while(1) {
        const pii pc_pre = traj[pos_emp.first][pos_emp.second].second;
        if(pos_emp == dest) break;
        // cerr << "emp before: " << pos_emp.first << " " << pos_emp.second << endl;
        // cerr << "pc before: " << pc_pre.first << " " << pc_pre.second << endl;
        PC pc_cur = pc_pos[pc_pre.first][pc_pre.second];
        PC pc_next = pc_pos[pos_emp.first][pos_emp.second];
        // cerr << "cur label: " << pc_cur.label << ", pc next label: " << pc_next.label << endl;
        if(pc_cur.label == '0') {
            pos_emp = pc_pre;
            continue;
        } else {   
            assert(abs(pc_cur.i - pc_next.i) + abs(pc_cur.j - pc_next.j) == 1);
            update_pc_pos(pc_cur, pc_next, movable_pcs, move_list);
            //cerr << "emp after: " << pos_emp.first << " " << pos_emp.second << endl;
            //cerr << "pc before: " << pc_pre.first << " " << pc_pre.second << endl;
            pos_emp = pc_pre;
        }
        if(move_list.size()>100*K-100) {
            cerr << "MOVE SIZE: " << move_list.size() << endl;
            halt=true;
            break;
        }
    }
    cerr << "\n### EMPTY SPACE MOVE DONE ####\n"; 

    if(halt) return false;
    // dfs part: 0のPCを特定の方向に移動させる。
    vector<vector<bool>> visited(N, vector<bool>(N, false));
    auto dfs = [&](auto dfs, const pii start, const pii pre)->void {
        bool flag = false;
        pii next;
        int ci = start.first;
        int cj = start.second;
        visited[ci][cj]=true;
        for(int dir = 0; dir < 4 and !flag; ++dir) {
            int ni = start.first + dy[dir];
            int nj = start.second + dx[dir];
            if(out_bound(ni, nj, N)) continue;
            if(ni == pre.first and nj == pre.second) continue;
            if(S[ni][nj]==target and !visited[ni][nj]) {
                PC pc_cur = pc_pos[ni][nj]; // PCの順番に注意
                PC pc_next = pc_pos[ci][cj];
                update_pc_pos(pc_cur, pc_next, movable_pcs, move_list);
                next = {ni, nj};
                dfs(dfs, next, start);
                return;
            }
        }
    };
    dfs(dfs, dest, {-1, -1});
    return true;
}

// Nが小さい時に有効
void find_and_merge(const char target) {
    cerr << "FIND AND MERGE" << endl;
    dsu uf(100);
    vector<pii> vp;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            // cerr << target << " " << S[i][j] << endl;
            if(S[i][j] != target) continue;
            PC pc_cur = pc_pos[i][j];
            vp.push_back({i, j});
            for(int dir = 0; dir < 2; ++dir) {
                int ni = i + dy[dir];
                int nj = j + dx[dir];                
                PC pc_next = pc_pos[ni][nj];
                if(out_bound(ni, nj, N)) continue;
                if(S[i][j] == S[ni][nj]) uf.merge(pc_cur.id, pc_next.id);
            }
        }
    }
    for(auto v: uf.groups()) {
        for(auto vv: v) {
            cerr << vv << " ";
        }
        cerr << endl;
    }
    // ランダムにBFSで繋げられそうなものを見つけていく。
    int NUM_LOOP=50;
    int depth = 10;
    for(int loop = 0; loop < NUM_LOOP; ++loop) {
        cerr << "BFS LOOP: " << loop + 1 << endl;
        uniform_int_distribution<> randM(0, uf.groups().size()-1);
        vector<int> start_group = uf.groups()[randM(mt)];
        int start = start_group[0];
        cerr << "START group: " << start << endl;
        int boss = uf.leader(start);
        pair<int,int> goal_pos; // S[ci][cj] == targetとなる地点
        vector<vector<pair<int, pii>>> traj(N, vector<pair<int, pii>>(N, {INF, pii{-1, -1}})); // pair<int,pii> はスタートからの距離と1つ前にいた箇所の座標。
        PC pc_start = *pcid[target-'1'][start];
        traj[pc_start.i][pc_start.j] = {0, {pc_start.i, pc_start.j}};
        deque<pii> dq;
        dq.push_back({pc_start.i, pc_start.j});
        bool flag = true;
        while(!dq.empty() && flag) {
            auto p = dq.front(); dq.pop_front();
            int ci = p.first;
            int cj = p.second;
            if(traj[ci][cj].first >= depth) continue;
            for(int dir = 0; dir < 4; ++ dir) {
                int ni = ci + dy[dir];
                int nj = cj + dx[dir];
                if(out_bound(ni, nj, N)) continue;
                if(S[ni][nj] == target and traj[ni][nj].first > traj[ci][cj].first) {
                    traj[ni][nj] = {traj[ci][cj].first, pii{ci, cj}};
                    dq.push_front({ni, nj});
                } else if(S[ni][nj] != target and traj[ni][nj].first > traj[ci][cj].first + 1) {
                    traj[ni][nj] = {traj[ci][cj].first + 1, pii{ci, cj}};
                    dq.push_back({ni, nj});
                }
                if(S[ni][nj]==target) {
                    if(uf.leader(pc_pos[ni][nj].id) != boss) {
                        goal_pos = {ni, nj};
                        flag = false;
                        break;
                    }
                }
            }
        }
        if(flag) continue; //範囲内で見つからなければ次に移る。
        cerr << "START PC(i, j): "  << pc_start.i << " " << pc_start.j << endl;
        cerr << "GOAL POS(i, j): " << goal_pos.first << " " << goal_pos.second << endl;
        int target_id = pc_pos[goal_pos.first][goal_pos.second].id;
        int source_id = boss;
        pair<int,int> dest = goal_pos;
        while(1) {
            dest = traj[dest.first][dest.second].second;
            cerr << "TRAJ (i, j): " << dest.first << " " << dest.second << endl;
            if(dest.first == pc_start.i and dest.second == pc_start.j) {
                uf.merge(target_id, source_id);
                break;
            }
            if(search_and_move(dest, target, uf, boss)==false) break;
            // pc_pos_debug();
        }
    }
}

void greedy_init_move(){
    int NUM_LOOP=10;
    for(int loop = 0; loop < NUM_LOOP; ++loop) {
        dsu uf(N * N + 1); // i * N+jがID
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                if(S[i][j]=='0') continue;
                int cid = i * N + j;
                for(int dir = 0; dir < 2; ++dir) {
                    int ni = i + dy[dir];
                    int nj = j + dx[dir];
                    int nid = ni * N + nj;
                    if(out_bound(ni,nj,N)) continue;
                    if(S[i][j]==S[ni][nj])uf.merge(cid, nid);
                }
            }
        }
        rep(i, N) {
            rep(j, N) {
                if(S[i][j]=='0') {
                    map<char, leaderNum> mp;
                    for(int dir = 0; dir < 4; ++dir) {
                        int ni = i + dy[dir];
                        int nj = j + dx[dir];
                        int nid = ni * N + nj;
                        if(out_bound(ni,nj,N)) continue;
                        if(S[ni][nj]=='0') continue;
                        mp[S[ni][nj]][uf.leader(nid)] = {uf.size(nid), dir};
                        
                    }
                    for(pair<char, leaderNum> p: mp) {
                        if(p.second.size()>=2) {
                            vector<pii> cand;
                            for(auto pp: p.second) cand.push_back(pp.second);
                            sort(cand.begin(), cand.end());
                            int dir = cand[0].second;
                            int ni = i + dy[dir];
                            int nj = j + dx[dir];

                            PC pc_next = pc_pos[i][j];
                            PC pc_cur = pc_pos[ni][nj];
                            update_pc_pos(pc_cur, pc_next, movable_pcs, move_list);
                            i++; j++;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void state_score() {
    // ポテンシャルを計算してシミュレーション・焼きなまし。
    uniform_int_distribution<> randN(0, N-1);
    uniform_real_distribution<> rand01(0, 1);
    map<tuple<char,int,int,int,int>, bool> move_memo;

    rep(i, N) rep(j, N) {
        bool movable=false;
        if(S[i][j]!='0') {
            for(int dir = 0; dir < 4 and !movable; ++dir) {
                int ni =i+dy[dir];
                int nj =j+dx[dir];
                if(out_bound(ni,nj,N)) continue;
                if (S[ni][nj] == '0') movable=true;
            }
        }
        if(movable) movable_pcs.push_back({i,j,S[i][j]});
    }

    // double strain[N][N]={};
    int mlis = move_list.size(); // move_list_init_size
    for(int iter = 0; iter < 100*K-100-mlis; ++iter) {
        double state[K][N][N] = {};
        for(int l = 1; l <= K; ++l) {
            char label = l + '0';
            rep(i, N) {
                rep(j, N) {
                    if(S[i][j]!=label)continue;
                    int m = 1;
                    // for(int dir = 0; dir < 4; ++dir) {
                    //     int ai=i+dy[dir];
                    //     int aj=j+dx[dir];
                    //     if(ai<0 or aj<0 or ai >= N or aj >= N) continue;
                    //     if(S[ai][aj]==label)m++;
                    // }
                    state[l-1][i][j]+=(double)m/sqrt(0.5/N);
                    int cnt=0;
                    for(int ii=1;ii+i<N; ++ii) {
                        state[l-1][ii+i][j]+=(double)m/sqrt((double)ii/N)*pow(PROPAGATION, cnt);
                        if(S[ii+i][j]!='0' and S[ii+i][j]!=label) cnt++;
                    }
                    cnt=0;
                    for(int ii=1;i-ii>=0; ++ii) {
                        state[l-1][i-ii][j]+=(double)m/sqrt((double)ii/N)*pow(PROPAGATION, cnt);
                        if(S[i-ii][j]!='0' and S[i-ii][j]!=label) cnt++;
                    }
                    cnt=0;
                    for(int jj=1;jj+j<N; ++jj) {
                        state[l-1][i][j+jj]+=(double)m/sqrt((double)jj/N)*pow(PROPAGATION, cnt);
                        if(S[i][j+jj]!='0' and S[i][j+jj]!=label) cnt++;
                    }
                    cnt=0;
                    for(int jj=1;j-jj>=0; ++jj) {
                        state[l-1][i][j-jj]+=(double)m/sqrt((double)j/N)*pow(PROPAGATION, cnt);
                        if(S[i][j-jj]!='0' and S[i][j-jj]!=label) cnt++;
                    }
                }
            }
        }

        #ifdef DEBUG
        rep(i, N) {
            rep(j, N) {
                printf("%c:%6.2lf, ", S[i][j], state[i][j]);
            }
            cout << endl;
        }
        #endif
        bool ok=true;
        int mi = 0;
        while(ok) {
            mi++;
            if(mi >= 100) break;
            int ci, cj;
            uniform_int_distribution<> randM(0, movable_pcs.size()-1);
            int pc_id=randM(mt);
            PC pc_cur=movable_pcs[randM(mt)];
            ci = pc_cur.i; cj = pc_cur.j;
            double cur_score=0;
            if(S[ci][cj]=='1') {
                bool flag=true;
                rep(ll, K) (ll==0) ? cur_score += (K-1)*state[ll][ci][cj] : cur_score-=state[ll][ci][cj];
                for(int dir = 0; dir < 4 && flag;++dir) {
                    int ni=ci+dy[dir];
                    int nj=cj+dx[dir];
                    tuple<char,int,int,int,int> m={S[ci][cj], ni,nj,ci,cj};
                    if(move_memo[m])continue;
                    double nx_score=0;
                    rep(ll, K)(ll==0) ? nx_score+=(K-1)*state[ll][ni][nj] : nx_score-=state[ll][ni][nj];
                    double p = exp(nx_score-cur_score);
                    if(out_bound(ni, nj, N)) continue;
                    if(p>rand01(mt) && S[ni][nj]=='0' and rand01(mt)<0.9) {
                        PC pc_next=pc_pos[ni][nj];
                        update_pc_pos(pc_cur, pc_next, movable_pcs, move_list);
                        flag=false;
                        ok=false;
                        move_memo[m]=true;
                    } 
                }
            } else if(S[ci][cj]!='1' and S[ci][cj]!='0'){
                bool flag=true;
                rep(ll, K) (ll+1==S[ci][cj]-'0') ? cur_score+=(K-1)*state[ll][ci][cj] : cur_score-=state[ll][ci][cj];
                for(int dir = 0; dir < 4 && flag;++dir) {
                    int ni=ci+dy[dir];
                    int nj=cj+dx[dir];
                    tuple<char,int,int,int,int> m={S[ci][cj],ni,nj,ci,cj};
                    if(move_memo[m])continue;
                    double nx_score=0;
                    rep(ll, K)(ll+1==S[ci][cj]-'0') ? nx_score+=(K-1)*state[ll][ni][nj] : nx_score-=state[ll][ni][nj];
                    double p = exp(nx_score-cur_score);
                    if(out_bound(ni, nj, N)) continue;
                    if(p*rand01(mt)>3 && S[ni][nj]=='0' and rand01(mt)<0.9) {
                        PC pc_next = pc_pos[ni][nj];
                        update_pc_pos(pc_cur, pc_next, movable_pcs, move_list);
                        flag=false;
                        ok=false;
                        move_memo[m]=true;
                    } 
                }
            }
        }
    }
}

void connect(char target, int remain, vector<Connect> & c_result, vector<vector<bool>>& wire) {
    // cerr << "CONNETCT TURN" << endl;
    vector<Connect> c_list, best_c_list;
    dsu uf(100);
    rep(i, N) {
        rep(j, N) {
            if(S[i][j] != target) continue;
            int cid = pc_pos[i][j].id;
            for(int dir = 0; dir < 2; ++dir) {
                int ni = i + dy[dir];
                int nj = j + dx[dir];
                int nid = pc_pos[ni][nj].id;
                if(out_bound(ni,nj,N)) continue;
                if(S[ni][nj] == S[i][j]) {
                    if(uf.same(cid, nid)) continue;
                    uf.merge(cid, nid);
                    c_result.push_back({i, j, ni, nj});
                    for(int ii = i; ii <= ni; ++ii) {
                        for (int jj = j; jj<= nj; ++jj) {
                            wire[ii][jj]=true;
                        }
                    }
                }
            } 
        }
    }
    // cerr << "merge done\n";
    vector<vector<bool>> wire_copy = wire;

    int best_score = 0;
    int TRIAL = 10;
    int DEPTH = 2;
    uniform_int_distribution<> rand100(0, 99);
    for(int t = 0; t < TRIAL; ++t) {
        c_list.clear();
        wire_copy=wire;
        int start_id = rand100(mt);
        int r = remain-c_result.size();
    
        bool flag = true;
        PC start_pc = *pcid[target-'1'][start_id];
        vector<vector<pair<int, pii>>> traj(N, vector<pair<int, pii> >(N, {INF, {-1,-1}}));
        traj[start_pc.i][start_pc.j] = {0, {start_pc.i, start_pc.j}};
        deque<pii> dq;
        dq.push_front({start_pc.i, start_pc.j});
        cerr << start_pc.i << " " << start_pc.j << " start position" << endl;
        cerr << "BEFORE: dquw\n";
        while(!dq.empty() and flag) {
            if(r==0) break;
            auto p = dq.front(); 
            int ci = p.first; int cj = p.second;
            // cerr << "state: " << ci << " " << cj << endl;
            dq.pop_front();
            for(int dir = 0; dir < 4; ++dir) {
                int ni = ci + dy[dir];
                int nj = cj + dx[dir];
                // cerr << ni << " " << nj << " ninj" << endl;
                if(out_bound(ni, nj, N)) continue;
                // if(wire_copy[ni][nj]) continue;
                if(S[ni][nj] == target and traj[ni][nj].first > traj[ci][cj].first) {
                    // c_list.push_back({ci, cj, ni, nj}); r--;
                    traj[ni][nj] = {traj[ci][cj].first, {ci, cj}};
                    dq.push_front({ni,nj});
                    if(r == 0) {
                        flag=false;
                        break;
                    } 
                } else if(S[ni][nj] != target and S[ni][nj] != '0' and traj[ni][nj].first > traj[ci][cj].first + 1) {
                    traj[ni][nj] = {traj[ci][cj].first + 1, {ci, cj}};
                    dq.push_back({ni, nj});
                } else if(S[ni][nj] == '0') {
                    int k = 2;
                    while(1) {
                        ni = ci + k * dy[dir];
                        nj = ci + k * dx[dir];
                        if(out_bound(ni,nj,N)) break;
                        if(wire_copy[ni][nj])break;
                        if(S[ni][nj] == '0') {
                            k++;
                            continue;
                        }
                        if(S[ni][nj] == target) {
                            dq.push_front({ni, nj});
                            traj[ni][nj] = {traj[ci][cj].first, {ci,cj}};
                            int nni = ni; int nnj = nj;
                            while(1) {
                                auto bp = traj[nni][nnj].second;
                                int bi = bp.first;
                                int bj = bp.second;
                                c_list.push_back({bi, bj, nni, nnj}); r--;
                                if(r==0) {flag=false; break;}
                                bool ok = false;
                                if(wire_copy[bi][bj]) ok=true;
                                for(int y = min(bi, nni); y <= max(bi, nni); ++y) {
                                    for(int x = min(bj, nnj); x <= max(bj, nnj); ++x) {
                                        wire_copy[y][x] = true;
                                    }
                                }
                                if(ok)break;
                                nni = bi;
                                nnj = bj;
                            }
                            break;
                        } else if(S[ni][nj] != target and traj[ni][nj].first == INF){
                            dq.push_back({ni, nj});
                            traj[ni][nj] = {traj[ci][cj].first + 1, {ci, cj}};
                            break;
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        int score = calc_score(S_ini, c_list);
        if(score > best_score) {
            best_score = score;
            best_c_list = c_list;
        }
        for(auto co: best_c_list) {
            if(c_result.size() + move_list.size() == 100 * 2) break;
            c_result.push_back(co);
        }
    }
    return ;
}

void show_result(vector<Connect>& c_result) {
    assert(move_list.size() + c_result.size() <= 100 * K);
    cout << move_list.size() << endl;
    for(auto m: move_list) printf("%d %d %d %d\n", m.ci, m.cj, m.ni, m.nj);
    cout << c_result.size() << endl;
    for(auto c: c_result) printf("%d %d %d %d\n", c.ci, c.cj, c.ni, c.nj);
}


int main() {
    cerr << "*******************************************START**************************************" << endl;
    // init part
    input();
    cerr << N << " " << K << endl;
    init();

    greedy_init_move();
    // pc_pos_debug();
    find_and_merge('1');
    // align();
    // solve1();
    // pc_pos_debug();
    #ifdef DEBUG
        trial_bfs();
    #endif
    state_score();
    int remain_num_operation = 100*5-move_list.size();
    vector<Connect> c_result;
    vector<vector<bool>> wire(N, vector<bool>(N, false));
    rep(i, N) {
        rep(j, N) {
            cerr << pc_pos[i][j].id << " ";
        }
        cerr << endl;
    }
    connect('1', remain_num_operation, c_result, wire);
    show_result(c_result);
    // trial_bfs();
    #ifdef DEBUG
        trial_bfs();
    #endif
    cerr << "SCORE: " << calc_score(S_ini, c_result) << endl;
    #ifdef TEST
        calc_score();
    #endif 
    return 0;
}