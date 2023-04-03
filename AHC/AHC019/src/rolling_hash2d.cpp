#include <iostream>
#include <vector>

using namespace std;

const int mod = 1000000009;
const int p = 13;
const int q = 19;

vector<int> p_table;
vector<int> q_table;

void prepare(int L) {
    p_table.resize(L+1); q_table.resize(L+1);
    p_table[0] = q_table[0] = 1;
    for (int i = 0; i < L; i++) {
        p_table[i+1] = (1LL * p_table[i] * p) % mod;
        q_table[i+1] = (1LL * q_table[i] * q) % mod;
    }
}

vector<vector<int>> rolling_hash(vector<vector<int>> S, int W, int H) {
    vector<vector<int>> D(H+1, vector<int>(W+1, 0));
    for (int i = 0; i < H; i++) {
        int su = 0;
        std::vector<int>& dp = D[i];
        std::vector<int>& di = D[i+1];
        std::vector<int>& si = S[i];
        for (int j = 0; j < W; j++) {
            int v = si[j];
            su = (1LL * su * p % mod + v) % mod;
            di[j+1] = (su + 1LL * dp[j+1] * q % mod) % mod;
        }
    }
    return D;
}

int get(vector<vector<int>> S, int x0, int y0, int x1, int y1) {
    int P = p_table[x1 - x0], Q = q_table[y1 - y0];
    return (S[y1][x1] - 1LL * S[y1][x0] * P % mod - 1LL * S[y0][x1] * Q % mod + 1LL * S[y0][x0] * (P * Q % mod) % mod + mod) % mod;
}

int main() {
    vector<vector<int>> data1 = {
        {1, 0, 1, 0, 1},
        {1, 1, 1, 0, 1},
        {0, 1, 1, 1, 1},
    };
    vector<vector<int>> data2 = {
        {1, 0, 1},
        {1, 0, 1},
    };
    prepare(5);
    auto rh1 = rolling_hash(data1, 5, 3);
    auto rh2 = rolling_hash(data2, 3, 2);
    int hash1 = get(rh1, 2, 0, 5, 2);
    int hash2 = get(rh2, 0, 0, 3, 2);
    cout << "hash1: " << hash1 << endl;
    cout << "hash2: " << hash2 << endl;
    cout << (hash1 == hash2) << endl; // output: 1
    return 0;
}
