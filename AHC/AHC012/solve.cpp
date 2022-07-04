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

#define rep(i,n) for(int i = 0; i < n; i++)

using namespace std;
using ll = long long;

random_device seed_gen;
mt19937 mt(seed_gen());

const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, -1, 0, 1};
const int dx8[8]={-1, -1, -1, 0, 1, 1, 1, 0};
const int dy8[8]={-1, 0, 1, -1, -1, 0, 1, 1};
const double TL=2.8;
const int INIT_POS=1000;

int N;
int K;

struct Cut {
    ll px, py, qx, qy;
};

struct Strawberry {
    ll x, y;
    string label;
    bool online=false;
};

map<ll, vector<string>> num_to_class;
map<string, ll> label_num;

int A[10];

ll calc_score(map<ll, vector<string>> &num_to_class) {
    ll ret=0;
    ll bunshi=0;
    ll bunbo=0;
    rep(i, 10) {
        bunbo+=A[i];
        int b = num_to_class[i].size();
        bunshi+=min(A[i], b);
    }
    ret=round(1e6*(double)bunshi/bunbo);
    return ret;
}

ll eval_score(map<ll, vector<string>> &num_to_class) {
    ll ret=0;
    ll bunshi=0;
    ll bunbo=0;
    for(auto m: num_to_class) {
        int num = m.first;
        int b=m.second.size();
        if(num<=10) {
            bunbo+=A[num-1];
            int b = num_to_class[num-1].size();
            bunshi+=min(A[num-1], b);
            // regularizer
            ret-=abs(A[num-1]-b)*10000;
        } else {
            ret -= num*b*10000;
        }
    }
    ret+=round(1e6*(double)bunshi/bunbo);
    return ret;
}

// IDEA
// やきなまし
void solve() {
    clock_t start, now;
    start = clock();
    cin >> N >> K;

    Cut c[K];
    Strawberry straw[N], straw_tmp[N];
    rep(i, 10) cin >> A[i];
    rep(i, N) {
        cin >> straw[i].x >> straw[i].y;
        straw[i].label = string(K, '.');
        straw_tmp[i].x=straw[i].x;
        straw_tmp[i].y=straw[i].y;
        straw_tmp[i].label=string(K, '.');
    }
    const int diff=700;
    // cut init
    rep(i, K/2) {
        c[i].px=i*diff-diff*25;
        c[i].qx=i*diff-diff*25;
        c[i].py=INIT_POS;
        c[i].qy=-INIT_POS;
    }
    rep(i, K/2) {
        c[i+50].py=i*diff-diff*25;
        c[i+50].qy=i*diff-diff*25;
        c[i+50].px=INIT_POS;
        c[i+50].qx=-INIT_POS;
        
    }

    // label init
    for(int k = 0; k < K; ++k) {
        for(int i = 0; i < N; ++i) {
            if(k<50) {
                if(straw[i].x < c[k].px) straw[i].label[k]='-';
                else if(straw[i].x > c[k].px) straw[i].label[k]='+';
                else straw[i].label[k]='.';
            } else {
                if(straw[i].y < c[k].py) straw[i].label[k]='-';
                else if(straw[i].y > c[k].py) straw[i].label[k]='+';
                else straw[i].label[k]='.';
            }
        }
    }
    // label count
    for(auto st: straw) label_num[st.label]++;
    for(auto ln: label_num) num_to_class[ln.second].push_back(ln.first);
  
    // for(auto ntc: num_to_class) {
    //     cerr << ntc.first << " " << ntc.second[0] << endl;
    // }
    now=clock();

    ll cur_score=eval_score(num_to_class);
    double init_temp=3000.0;
    const int SCALE_FACTOR=1000;
    uniform_int_distribution<> randK(0, K-1);
    uniform_int_distribution<> rand_scale(0, SCALE_FACTOR);
    uniform_real_distribution<double> rand01(0, 1);


    cerr << "INITIAL SCORE:" << cur_score << endl;

    int iteration = 0;
    double temp=init_temp;
    const double DECAY=0.999;
    while(((double)now-start)/CLOCKS_PER_SEC < TL) {
        // cerr << "ITERATION: " << ++iteration << endl;
        if(rand01(mt)<0.5) {
            int target_cut_id = randK(mt);
            int scale=rand_scale(mt);
            Cut before=c[target_cut_id];
            c[target_cut_id].px += rand01(mt)*scale;
            c[target_cut_id].py += rand01(mt)*scale;
            c[target_cut_id].qx += rand01(mt)*scale;
            c[target_cut_id].qy += rand01(mt)*scale;
            
            for(int k = 0; k < K; ++k) {
                for(int i = 0; i < N; ++i) {
                    if(straw_tmp[i].online) continue;
                    int l = (c[k].px-c[k].qx)*(straw[i].y-c[k].py);
                    int r = (c[k].py-c[k].qy)*(straw[i].x-c[k].px);
                    if(l>r) straw_tmp[i].label[k]='+';
                    else if(r>l) straw_tmp[i].label[k]='-';
                    else {
                        straw_tmp[i].label[k]='.';
                        straw_tmp[i].online=true;
                    }
                }
            }
           
            // label count
            map<string, ll> label_num_tmp;
            map<ll, vector<string>> num_to_class_tmp;
            for(auto st: straw_tmp) {
                if(st.online) st.online=false;
                else label_num_tmp[st.label]++;
            }
            for(auto ln: label_num_tmp) num_to_class_tmp[ln.second].push_back(ln.first);
        
            // for(auto ntc: num_to_class_tmp) {
            //     cerr << ntc.first << " " << ntc.second.size() << endl;
            // }
            
            ll nx_score=eval_score(num_to_class_tmp);
            // cerr << "NEXT SCORE: " << nx_score << endl;
            if(nx_score>cur_score) {
                cur_score=nx_score;
                num_to_class=num_to_class_tmp;
            } else {
                ll dE=cur_score-nx_score;
                double p = exp(-dE/temp);
                // cerr << "dE: " << dE << ", TEMP: " << temp << ", P: " << p << endl;
                if(p>rand01(mt)) {
                    cur_score=nx_score;
                    num_to_class=num_to_class_tmp;
                } else {
                    c[target_cut_id]=before;
                }
            }
        } else {
            int target_cut_id1 = randK(mt);
            int target_cut_id2 = randK(mt);
            int scale=rand_scale(mt);
            Cut before1=c[target_cut_id1];
            c[target_cut_id1].px += rand01(mt)*scale;
            c[target_cut_id1].py += rand01(mt)*scale;
            c[target_cut_id1].qx += rand01(mt)*scale;
            c[target_cut_id1].qy += rand01(mt)*scale;
            Cut before2=c[target_cut_id2];
            c[target_cut_id2].px += rand01(mt)*scale;
            c[target_cut_id2].py += rand01(mt)*scale;
            c[target_cut_id2].qx += rand01(mt)*scale;
            c[target_cut_id2].qy += rand01(mt)*scale;
            for(int k = 0; k < K; ++k) {
                for(int i = 0; i < N; ++i) {
                    if(straw_tmp[i].online) continue;
                    int l = (c[k].px-c[k].qx)*(straw[i].y-c[k].py);
                    int r = (c[k].py-c[k].qy)*(straw[i].x-c[k].px);
                    if(l>r) straw_tmp[i].label[k]='+';
                    else if(r>l) straw_tmp[i].label[k]='-';
                    else {
                        straw_tmp[i].label[k]='.';
                        straw_tmp[i].online=true;
                    }
                }
            }

            // label count
            map<string, ll> label_num_tmp;
            map<ll, vector<string>> num_to_class_tmp;
            for(auto st: straw_tmp) {
                if(st.online) st.online=false;
                else label_num_tmp[st.label]++;
            }
            for(auto ln: label_num_tmp) num_to_class_tmp[ln.second].push_back(ln.first);
        
            // for(auto ntc: num_to_class_tmp) {
            //     cerr << ntc.first << " " << ntc.second.size() << endl;
            // }
            
            ll nx_score=eval_score(num_to_class_tmp);
            // cerr << "NEXT SCORE: " << nx_score << endl;
            if(nx_score>cur_score) {
                cur_score=nx_score;
                num_to_class=num_to_class_tmp;
            } else {
                ll dE=cur_score-nx_score;
                double p = exp(-dE/temp);
                // cerr << "dE: " << dE << ", TEMP: " << temp << ", P: " << p << endl;
                if(p>rand01(mt)) {
                    cur_score=nx_score;
                    num_to_class=num_to_class_tmp;
                } else {
                    c[target_cut_id1]=before1;
                    c[target_cut_id2]=before2;
                }
            }

        }
        temp*=DECAY;
        now=clock();
    }

    //output
    cout << K << endl;
    for(auto cc: c) printf("%lld %lld %lld %lld\n", cc.px, cc.py, cc.qx, cc.qy);
    cerr << "FINAL SCORE: " << cur_score << endl; 
    // for(auto aa: num_to_class) {
    //     if(aa.first > 10) break;
    //     cerr << "N: " << aa.first << endl;
    //     for(auto s: aa.second) cerr << s <<" ";
    //     cerr << endl;
    // }
}

int main() {
    solve();
    return 0;
}
