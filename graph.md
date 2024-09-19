# 网络流

### 最大流、最小割
```cpp
template<class T>
struct MaxFlow {
    struct _Edge {
        int to;
        T cap;
        _Edge(int to, T cap) : to(to), cap(cap) {}
    };

    int n;
    std::vector<_Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<int> cur, h;

    MaxFlow() {}
    MaxFlow(int n) {
        init(n);
    }

    void init(int n) {
        this->n = n;
        e.clear();
        g.assign(n, {});
        cur.resize(n);
        h.resize(n);
    }

    bool bfs(int s, int t) {
        h.assign(n, -1);
        std::queue<int> que;
        h[s] = 0;
        que.push(s);
        while (!que.empty()) {
            const int u = que.front();
            que.pop();
            for (int i : g[u]) {
                auto [v, c] = e[i];
                if (c > 0 && h[v] == -1) {
                    h[v] = h[u] + 1;
                    if (v == t) {
                        return true;
                    }
                    que.push(v);
                }
            }
        }
        return false;
    }

    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        auto r = f;
        for (int &i = cur[u]; i < int(g[u].size()); ++i) {
            const int j = g[u][i];
            auto [v, c] = e[j];
            if (c > 0 && h[v] == h[u] + 1) {
                auto a = dfs(v, t, std::min(r, c));
                e[j].cap -= a;
                e[j ^ 1].cap += a;
                r -= a;
                if (r == 0) {
                    return f;
                }
            }
        }
        return f - r;
    }
    void addEdge(int u, int v, T c) {
        g[u].push_back(e.size());
        e.emplace_back(v, c);
        g[v].push_back(e.size());
        e.emplace_back(u, 0);
    }
    T flow(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, std::numeric_limits<T>::max());
        }
        return ans;
    }

    std::vector<bool> minCut() {
        std::vector<bool> c(n);
        for (int i = 0; i < n; i++) {
            c[i] = (h[i] != -1);
        }
        return c;
    }

    struct Edge {
        int from;
        int to;
        T cap;
        T flow;
    };
    std::vector<Edge> edges() {
        std::vector<Edge> a;
        for (int i = 0; i < e.size(); i += 2) {
            Edge x;
            x.from = e[i + 1].to;
            x.to = e[i].to;
            x.cap = e[i].cap + e[i + 1].cap;
            x.flow = e[i + 1].cap;
            a.push_back(x);
        }
        return a;
    }
};
```
### 最小费用最大流
```cpp
template<class T>
struct MaxFlow {
    struct _Edge {
        int to;
        T cap;
        _Edge(int to, T cap) : to(to), cap(cap) {}
    };
    
    int n;
    std::vector<_Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<int> cur, h;
    
    MaxFlow() {}
    MaxFlow(int n) {
        init(n);
    }
    
    void init(int n) {
        this->n = n;
        e.clear();
        g.assign(n, {});
        cur.resize(n);
        h.resize(n);
    }
    
    bool bfs(int s, int t) {
        h.assign(n, -1);
        std::queue<int> que;
        h[s] = 0;
        que.push(s);
        while (!que.empty()) {
            const int u = que.front();
            que.pop();
            for (int i : g[u]) {
                auto [v, c] = e[i];
                if (c > 0 && h[v] == -1) {
                    h[v] = h[u] + 1;
                    if (v == t) {
                        return true;
                    }
                    que.push(v);
                }
            }
        }
        return false;
    }
    
    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        auto r = f;
        for (int &i = cur[u]; i < int(g[u].size()); ++i) {
            const int j = g[u][i];
            auto [v, c] = e[j];
            if (c > 0 && h[v] == h[u] + 1) {
                auto a = dfs(v, t, std::min(r, c));
                e[j].cap -= a;
                e[j ^ 1].cap += a;
                r -= a;
                if (r == 0) {
                    return f;
                }
            }
        }
        return f - r;
    }
    void addEdge(int u, int v, T c) {
        g[u].push_back(e.size());
        e.emplace_back(v, c);
        g[v].push_back(e.size());
        e.emplace_back(u, 0);
    }
    T flow(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, std::numeric_limits<T>::max());
        }
        return ans;
    }
    
    std::vector<bool> minCut() {
        std::vector<bool> c(n);
        for (int i = 0; i < n; i++) {
            c[i] = (h[i] != -1);
        }
        return c;
    }
    
    struct Edge {
        int from;
        int to;
        T cap;
        T flow;
    };
    std::vector<Edge> edges() {
        std::vector<Edge> a;
        for (int i = 0; i < e.size(); i += 2) {
            Edge x;
            x.from = e[i + 1].to;
            x.to = e[i].to;
            x.cap = e[i].cap + e[i + 1].cap;
            x.flow = e[i + 1].cap;
            a.push_back(x);
        }
        return a;
    }
};
```

### tarjan ecc 缩点，求桥
```cpp
struct edge
{
    int v, ne;
} e[M];
int h[N], idx = 1;
int dfn[N], low[N], tot;
stack<int> stk;
int dcc[N], cnt;
int bri[M], d[N];

void add(int a, int b)
{
    e[++idx].v = b;
    e[idx].ne = h[a];
    h[a] = idx;
}
void tarjan(int x, int in_edg)
{
    dfn[x] = low[x] = ++tot;
    stk.push(x);
    for (int i = h[x]; i; i = e[i].ne)
    {
        int y = e[i].v;
        if (!dfn[y])
        {
            tarjan(y, i);
            low[x] = min(low[x], low[y]);
            if (low[y] > dfn[x])
                bri[i] = bri[i ^ 1] = true;
        }
        else if (i != (in_edg ^ 1))
            low[x] = min(low[x], dfn[y]);
    }
    if (dfn[x] == low[x])
    {
        ++cnt;
        while (1)
        {
            int y = stk.top();
            stk.pop();
            dcc[y] = cnt;
            if (y == x)
                break;
        }
    }
}
```
### tarjan vcc 缩点，求割点
```cpp
const int MAXN = 1e5 + 10;
int n, m;
vector<int> edge[MAXN];
vector<int> vcc[MAXN];
vector<int> vcc_edge[MAXN];
vector<int> id(MAXN);
int dfn[MAXN], low[MAXN];
bool cut[MAXN];
int tot = 0,cnt = 0;
stack<int> stk;
void tarjan(int u){
    dfn[u] = low[u] = ++tot;
    stk.emplace(u);
    if (edge[u].empty()){
        vcc[++cnt].emplace_back(u);
        return ;
    }
    for (auto &it : edge[u]){
        if (!dfn[it]){
            tarjan(it);
            low[u] = min(low[u],low[it]);
            if (low[it] >= dfn[u]){
                cut[u] = true;
            }
            cnt++;
            int y;
            do
            {
                y = stk.top();
                vcc[cnt].emplace_back(y);
                stk.pop();
            } while (y == it);
            
            vcc[cnt].emplace_back(u);
        }
        else{
            low[u] = min(low[u],dfn[it]);
        }
    }
}
```

### tarjan scc 缩点，求强连通分量
```cpp
vector<int> edge[N];
int dfn[N], low[N], top = 0, stk[N];
bool instk[N];
int timer = 0, scc_cnt = 0;
//无论对于一个有向图还是无向图来说，缩点后的结果一定是一个森林
void tarjan(int u)
{
    low[u] = dfn[u] = ++timer;
    stk[++top] = u;
    instk[u] = true;
    for (auto v : edge[u])
    {
        if (!dfn[v])
        {
            tarjan(v);
            low[u] = min(low[u], low[v]);
        }
        else if (instk[v])
        {
            low[u] = min(low[u], low[v]);
        }
    }

    if (dfn[u] == low[u])
    {
        scc_cnt ++;
        scc[u] = scc_cnt;
        instk[u] = false;
        while (stk[top] != u) {
            scc[stk[top]] = scc_cnt;
            instk[stk[top--]] = false;
        }
        top --;
    }
}
```