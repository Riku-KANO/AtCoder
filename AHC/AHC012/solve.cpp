// 78M
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
#include<tuple>

#define rep(i,n) for(int i = 0; i < n; i++)

using namespace std;
using ll = long long;

random_device seed_gen;
mt19937 mt(seed_gen());

const int dx[4] = {-1, 0, 1, 0};
const int dy[4] = {0, -1, 0, 1};
const int dx8[8]={-1, -1, -1, 0, 1, 1, 1, 0};
const int dy8[8]={-1, 0, 1, -1, -1, 0, 1, 1};
constexpr double TL=2.9;
constexpr int INIT_POS=1e9;
const int NUM_HORIZON_CUT = 8;
constexpr int RADIUS = 10000;
int N;
int K;

struct Point {
  int x, y;
};

struct Cut {
  Point p1, p2;
};

struct Strawberry {
  ll x, y;
  string label;
};

map<ll, vector<string>> num_to_class;
map<string, ll> label_num;

int A[100];

template<typename T> T pow2(T a) {return a * a;}

ll calc_score(map<ll, vector<string>> &num_to_class) {
  ll ret=0;
  ll bunshi=0;
  ll bunbo=0;
  rep(i, 10) {
    bunbo+=A[i];
    int b = num_to_class[i+1].size();
    bunshi+=min(A[i], b);
  }
  ret=round(1e6*(double)bunshi/bunbo);
  return ret;
}

ll eval_score(map<ll, vector<string>> &num_to_class) {
  ll ret=0;
  ll bunshi=0;
  ll bunbo=0;
  for(auto& p: num_to_class) {
    int num = p.first;
    int b = p.second.size();
    bunbo+=A[num-1];
    bunshi+=pow2(A[num-1]-b);
  }
  ret+=round(1e6*(double)bunshi/bunbo);
  return -ret;
}

ll eval_score2(map<ll, vector<string>> &num_to_class) {
  ll ret=0;
  ll bunshi=0;
  ll bunbo=0;
  for(auto& p: num_to_class) {
    int num = p.first;
    int b = p.second.size();
    bunbo+=A[num-1];
    bunshi+= num * min(A[num-1],b);
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

  Cut cutH[NUM_HORIZON_CUT];   // horizontal cut
  Cut cutV[K-NUM_HORIZON_CUT]; // vertical cut
  Cut nxCutH[NUM_HORIZON_CUT];
  Cut nxCutV[K-NUM_HORIZON_CUT];
  Cut bestCutH[NUM_HORIZON_CUT];
  Cut bestCutV[K-NUM_HORIZON_CUT];

  Strawberry berry[N]; 
  rep(i, 10) cin >> A[i];
  for(int i = 10; i < 100; i++) A[i] = 0;
  rep(i, N) {
      cin >> berry[i].x >> berry[i].y;
      berry[i].label = string(K, '.');
  }



  // cut init
  const int dH = RADIUS * 2 / (NUM_HORIZON_CUT + 1);
  const int dW = RADIUS * 2 / (K - NUM_HORIZON_CUT + 1);

  for(int i = 0; i < NUM_HORIZON_CUT; i++) {
    int px = -1e9, qx = 1e9;
    int py = -RADIUS + dH * (i + 1);
    int qy = -RADIUS + dH * (i + 1) + 1;
    cutH[i] = nxCutH[i] = bestCutH[i] = {Point{px, py}, Point{qx, qy}};
  }

  for(int i = 0; i < K - NUM_HORIZON_CUT; i++) {
    int py = -1e9, qy = 1e9;
    int px = - RADIUS + dW * (i + 1);
    int qx = - RADIUS + dW * (i + 1) + 1;
    cutV[i] = nxCutV[i] = bestCutV[i] = {Point{px, py}, Point{qx, qy}};
  }

  // label init
  for(int i = 0; i < NUM_HORIZON_CUT; i++) {
    for(int j = 0; j < N; j++) {
      if(berry[j].y >= cutH[i].p2.y) berry[j].label[i] = '+';
      else                           berry[j].label[i] = '-';
    }
  }

  for(int i = NUM_HORIZON_CUT; i < K; i++) {
    for(int j = 0; j < N; j++) {
      if(berry[j].x >= cutV[i-NUM_HORIZON_CUT].p2.x) berry[j].label[i] = '+';
      else                                           berry[j].label[i] = '-';
    }
  }

  // label count
  for(auto st: berry) label_num[st.label]++;
  for(auto ln: label_num) num_to_class[ln.second].push_back(ln.first);

  now=clock();

  ll cur_score=eval_score2(num_to_class);
  ll best_score = cur_score;
  
  double init_temp=4000000.0;
  uniform_int_distribution<> randK(0, K-1);
  uniform_int_distribution<> randH(-dH, dH);
  uniform_int_distribution<> randW(-dW, dW);
  uniform_int_distribution<> randCut(1, 10);
  uniform_real_distribution<double> rand01(0, 1);


  cerr << "INITIAL SCORE:" << cur_score << endl;

  int iteration = 0;
  double temp=init_temp;
  const double DECAY=0.995;

  // annealing
  while(((double)now-start) / CLOCKS_PER_SEC < TL) {

    // horizontal cut change
    if(rand01(mt) < 0.1) {
      for(int i = 0; i < NUM_HORIZON_CUT; i++) {
        int diff = randH(mt) / 5;
        bool ok = true;
        if(i>0) if(nxCutH[i-1].p1.y >= nxCutH[i].p1.y + diff) ok = false;
        if(i < NUM_HORIZON_CUT-1) if(nxCutH[i+1].p1.y <= nxCutH[i].p1.y + diff) ok = false;
        if(ok) {
          nxCutH[i].p1.y += diff;
          nxCutH[i].p2.y += diff;
        }
      }
    }

    // vertical cut change
    if(rand01(mt) < 1) {
      for(int i = 0; i < K - NUM_HORIZON_CUT; i++) {
        if(mt()%5>0) continue;
        int diff = randW(mt) / 5;
        bool ok = true;
        if(i>0) if(nxCutV[i-1].p1.x >= nxCutV[i].p1.x + diff) ok = false;
        if(i < K-NUM_HORIZON_CUT-1) if(nxCutV[i+1].p1.x <= nxCutV[i].p1.x + diff) ok = false;
        if(ok) {
          nxCutV[i].p1.x += diff;
          nxCutV[i].p2.x += diff;
        }
      }
    }

    for(int i = 0; i < NUM_HORIZON_CUT; i++) {
      for(int j = 0; j < N; j++) {
        if(berry[j].y >= nxCutH[i].p2.y) berry[j].label[i] = '+';
        else                             berry[j].label[i] = '-';
      }
    }

    for(int i = NUM_HORIZON_CUT; i < K; i++) {
      for(int j = 0; j < N; j++) {
        if(berry[j].x >= nxCutV[i-NUM_HORIZON_CUT].p2.x) berry[j].label[i] = '+';
        else                                             berry[j].label[i] = '-';
      }
    }

    // label count
    map<string, ll> nx_label_num;
    map<ll, vector<string>> nx_num_to_class;
    for(auto b: berry) nx_label_num[b.label]++;
    for(auto ln: nx_label_num) nx_num_to_class[ln.second].push_back(ln.first);

    ll nx_score=eval_score2(nx_num_to_class);
        
    cerr << nx_score << " " << best_score << endl;

    if(nx_score > best_score) {
      best_score = nx_score;
      num_to_class = nx_num_to_class;
      for(int i = 0; i < K; i++) {
        if(i < NUM_HORIZON_CUT) bestCutH[i] = nxCutH[i];
        else bestCutV[i-NUM_HORIZON_CUT] = nxCutV[i-NUM_HORIZON_CUT];
      }

      #ifdef LOCAL
      cout << K << endl;
      for(int i= 0; i < K; i++) {
        Cut c;
        if(i < NUM_HORIZON_CUT) c = bestCutH[i];
        else c = bestCutV[i-NUM_HORIZON_CUT];
        cout << c.p1.x << " " << c.p1.y << " " << c.p2.x << " " << c.p2.y << "\n";
      }
      #endif
    }

    if(nx_score>cur_score) {
      cur_score=nx_score;
      num_to_class=nx_num_to_class;
      for(int i = 0; i < K; i++) {
        if(i < NUM_HORIZON_CUT) cutH[i] = nxCutH[i];
        else cutV[i-NUM_HORIZON_CUT] = nxCutV[i-NUM_HORIZON_CUT];
      }
    } else {
      double dE=cur_score-nx_score;
      double p = exp(-dE/temp);
      if(p>rand01(mt)) {
        num_to_class=nx_num_to_class;
        cur_score = nx_score;
        for(int i = 0; i < K; i++) {
          if(i < NUM_HORIZON_CUT) cutH[i] = nxCutH[i];
          else cutV[i-NUM_HORIZON_CUT] = nxCutV[i-NUM_HORIZON_CUT];
        }
      } else {
        for(int i = 0; i < K; i++) {
          if(i < NUM_HORIZON_CUT) nxCutH[i] = cutH[i];
          else nxCutV[i-NUM_HORIZON_CUT] = cutV[i-NUM_HORIZON_CUT];
        }
      }
    }

    
    temp*=DECAY;
    now=clock();
  }


  //output
  cout << K << endl;
  for(int i = 0; i < K; i++) {
    Cut c;
    if(i < NUM_HORIZON_CUT) c = bestCutH[i];
    else c = bestCutV[i-NUM_HORIZON_CUT];
    printf("%d %d %d %d\n", c.p1.x, c.p1.y, c.p2.x, c.p2.y);
  }
  cerr << "FINAL SCORE: " << calc_score(num_to_class) <<", EVAL:" << cur_score << endl; 
    
}

int main() {
    solve();
    return 0;
}
