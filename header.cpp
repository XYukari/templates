#pragma GCC optimize(3)
#include <bits/stdc++.h>
using namespace std;
#define __LOCAL_DEBUG__

#ifdef __LOCAL_DEBUG__
#define _debug(fmt, ...) fprintf(stderr, "\033[91m[%s␣%3d]:␣" fmt "\n\033[0m", __func__, __LINE__, ##__VA_ARGS__)
#else
#define _debug(...) (void(0))
#endif
#define MP(a, b) make_pair(a, b)
#define GC() getchar()
#define PB(x) push_back(x)
#define EB(x) emplace_back(x)
#define rep(i, a, b) for (int i = a, _ = b; i <= _; i++)
#define dep(i, a, b) for (int i = a, _ = b; i >= _; i++)
#define all(x) x.begin(), x.end()
#define fi first
#define se second
#define untie                        \
    do {                             \
        ios::sync_with_stdio(false); \
        cin.tie(nullptr);            \
        cout.tie(nullptr);           \
    } while (0)
typedef long long ll;
typedef long long LL;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ll> vl;
typedef long double db;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
const int INF = 0x3f3f3f3f;
const LL INF_LL = 0x3f3f3f3f3f3f3f3fLL;
const LL mod = 998244353;
mt19937 gen(time(0));

int rnd(int x) { return mrand48() % x; }
ll gcd(ll a, ll b) { return b ? gcd(b, a % b) : a; }

ll qpow(ll a, ll n) {
    ll s = 1;
    for (; n; n >>= 1) {
        if (n & 1) s = s * a % mod;
        a = a * a % mod;
    }
    return s;
}

ll read() {
    ll x = 0, f = 1;
    char c = GC();
    for (; !isdigit(c); c = GC())
        if (c == '-') f = 0;
    for (; isdigit(c); c = GC())
        x = (x << 1) + (x << 3) + (c ^ 48);
    return f ? x : -x;
}
/********** header ************/

void solve() {
    return;
}

signed main() {
    int t = 1;
    // t = read();
    while (t--) solve();
    return 0;
}