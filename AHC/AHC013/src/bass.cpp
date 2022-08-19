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
// #define DEBUG
// #define TEST

#define rep(i, n) for(int i = 0; i < (n); i++)
#define pii pair<int,int>
using namespace std;
using ll = long long;

// global variable//////
double TL;
int N, K;
string S[50];
double PROPAGATION=0.90;
////////////////////////

class Config {
    public:
        double tl=3.0;

        string file_name="../output/config.txt";
        Config() {}
        void init() {
            TL=this->tl;
            return ;
        }

        void save() {
            ofstream ofs;
            ofs.open(this->file_name, ios::out);
            cout << "TL: " << this->tl << "\n";
            ofs.close();
            return ;
        }
};

random_device seed_gen;
mt19937 mt(seed_gen());

const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, -1, 0, 1};
const int dx8[8]={-1, -1, -1, 0, 1, 1, 1, 0};
const int dy8[8]={-1, 0, 1, -1, -1, 0, 1, 1};
const int INF = 1<<30;

struct PC {
    int i, j;
    char laebl;
    int gid=-1;
};

struct Move {
    int ci, cj, ni, nj;
};

bool out_bound(int ni, int nj, int n) {
    if(ni < 0 || nj < 0 || ni >= n || nj >=n) return true;
    return false;
}

ll calc_score() {
    return 0;
}

ll eval_score() {
    return 0;
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
    rep(i, N) cin >> S[i];
}

void init() {
    rep(i, N) {
        rep(j, N) {
            pc_pos[i][j] = {i, j, S[i][j], -1};
        }
    }
    return;
}

// IDEA
// 操作の回数に上限があるのでやりたい放題できない。
// AHC011のように最適な状態を前もって計算しておいてそれに近づけた方がよさそうか。
// 100*Kなので移動回数は少ない方がよさそう。
// 同じ種類のパソコンをなるべく多く繋げられればスコアがとても良くなる。
// 同じ種類のパソコンをつなげることができたら100*99/2=4950点獲得できる。
// 細分化(2*2or3*3)して大きなスコアを作ってみるか。
// BFS + dequeでコストを調べようとしたがK=5の時の計算(特にN<30)がきつい。K＝4も少し遅いが許容できそうなレベル。K=2,K=3は余裕。
// solve1はとりあえずラベル1のPCに注目して最適化してみる。
vector<vector<PC>> pc_pos;

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
            PC tmp = {i, j, S[i][j]};
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
                            dq.push_front((pair<PC,int>){{ni,nj,S[ni][nj]}, d});
                            if(S[ni][nj]==S[now.i][now.j]) groups_tmp[ni][nj]=groups_tmp[now.i][now.j];
                            else groups_tmp[ni][nj]=++gid;
                            break;
                        }
                        else if(S[ni][nj]!=label && d+1<bfs_tmp_state[ni][nj]) {
                            dq.push_back((pair<PC,int>){{ni,nj,S[ni][nj]}, d+1});
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

void state_score() {
    vector<Move> move_list;
    vector<PC> movable_pcs;
    uniform_int_distribution<> randN(0, N-1);
    uniform_real_distribution<> rand01(0, 1);
    
    rep(i, N) rep(j, N) {
        bool movable=true;
        if(S[i][j]!='0') {
            for(int dir = 0; dir < 4 and movable; ++dir) {
                int ni =i+dy[dir];
                int nj =j+dx[dir];
                if(out_bound(ni,nj,N)) continue;
                if (S[ni][nj] != '0') movable=false;
            }
        }
        if(movable) movable_pcs.push_back({i,j,S[i][j]})
    }

    double strain[N][N]={};
    for(int iter = 0; iter < 100*K/2; ++iter) {
        double state[K][N][N] = {};
        for(int l = 1; l <= K; ++l) {
            char label = l + '0';
            rep(i, N) {
                rep(j, N) {
                    if(S[i][j]!=label)continue;
                    int m = 1;
                    for(int dir = 0; dir < 4; ++dir) {
                        int ai=i+dy[dir];
                        int aj=j+dx[dir];
                        if(ai<0 or aj<0 or ai >= N or aj >= N) continue;
                        if(S[ai][aj]==label)m++;
                    }
                    int cnt=0;
                    for(int ii=1;ii+i<N; ++ii) {
                        state[l-1][ii+i][j]+=(double)m/sqrt(ii)*pow(PROPAGATION, cnt);
                        if(S[ii+i][j]!='0' and S[ii+i][j]!=label) cnt++;
                    }
                    cnt=0;
                    for(int ii=1;i-ii>=0; ++ii) {
                        state[l-1][i-ii][j]+=(double)m/sqrt(ii)*pow(PROPAGATION, cnt);
                        if(S[i-ii][j]!='0' and S[i-ii][j]!=label) cnt++;
                    }
                    cnt=0;
                    for(int jj=1;jj+j<N; ++jj) {
                        state[l-1][i][j+jj]+=(double)m/sqrt(jj)*pow(PROPAGATION, cnt);
                        if(S[i][j+jj]!='0' and S[i][jj+j]!=label) cnt++;
                    }
                    cnt=0;
                    for(int jj=1;j-jj>=0; ++jj) {
                        state[l-1][i][j-jj]+=(double)m/sqrt(jj)*pow(PROPAGATION, cnt);
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
            if(mi >= 80) break;
            int ci, cj;
            ci=randN(mt); cj=randN(mt);
            double cur_score=0;
            rep(ll, K) (ll+1==1) ? cur_score+=state[ll][ci][cj] : cur_score-=state[ll][ci][cj];
            if(S[ci][cj]=='1') {
                bool flag=true;
                for(int dir = 0; dir < 4 && flag;++dir) {
                    int ni=ci+dy[dir];
                    int nj=cj+dx[dir];
                    double nx_score=0;
                    rep(ll, K) (ll+1==1) ? nx_score+=state[ll][ni][nj] : nx_score-=state[ll][ni][nj];
                    if(out_bound(ni, nj, N)) continue;
                    if(nx_score > cur_score && S[ni][nj]=='0' && rand01(mt)<0.5) {
                        move_list.push_back({ci, cj, ni, nj});
                        swap(S[ci][cj], S[ni][nj]);
                        flag=false;
                        ok=false;
                    } 
                }
            } else if(S[ci][cj]!='1' and S[ci][cj]!='0'){
                bool flag=true;
                for(int dir = 0; dir < 4 && flag;++dir) {
                    int ni=ci+dy[dir];
                    int nj=cj+dx[dir];
                    double nx_score=0;
                    rep(ll, K) (ll+1==S[ci][cj]-'0') ? nx_score+=state[ll][ni][nj] : nx_score-=state[ll][ni][nj];
                    if(out_bound(ni, nj, N)) continue;
                    if(nx_score < cur_score && S[ni][nj]=='0' && rand01(mt)<0.5) {
                        move_list.push_back({ci, cj, ni, nj});
                        swap(S[ci][cj], S[ni][nj]);
                        flag=false;
                        ok=false;
                    } 
                }
            }
        }
    }
    cout << move_list.size() << endl;
    for(auto move: move_list) printf("%d %d %d %d\n", move.ci, move.cj, move.ni, move.nj);
    
}

int main() {
    input();
    init();
    // align();
    // solve1();
    trial_bfs();
    state_score();
    trial_bfs();
    return 0;
}