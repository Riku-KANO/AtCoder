// N<=10^11
// N=sqrt(N)あたりでb進数が2桁になることに気づくことがポイント
#include <iostream>
#include <map>
using namespace std;

const long long LINF = 1LL << 60;
map<long long, long long> memo;

long long rec(long long b, long long n)
{
    if (memo[b] != 0)
        return memo[b];
    if (b < 2)
        return LINF;
    long long ret = 0;
    while (n != 0)
    {
        ret += n % b;
        n /= b;
    }
    return memo[b] = ret;
}

void solve()
{
    long long n, s;
    cin >> n >> s;
    if (n == s)
    {
        cout << n + 1 << endl;
        return;
    }
    else if (n < s)
    {
        cout << -1 << endl;
        return;
    }
    for (long long b = 2; b * b <= n; b++)
    {
        long long f = rec(b, n);
        if (f == s)
        {
            cout << b << endl;
            return;
        }
    }

    long long ans = LINF;
    for (long long p = 1; p * p < n; p++)
    {
        long long b = (n - s) / p + 1;
        long long f = rec(b, n);
        if (f == s)
            ans = min(ans, b);
    }
    if (ans == LINF)
        cout << -1 << endl;
    else
    {
        cout << ans << endl;
    }
}

int main()
{
    solve();
    return 0;
}