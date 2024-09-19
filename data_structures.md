## 数据结构

### 并查集

```cpp
vector<int> fa(n + 1); //扩展域并查集注意开n*3+1
iota(fa.begin(), fa.end(), 0);
// 带权并查集则同时更新d[x],siz[x]
function<int(int)> find = [&](int x) { return x == fa[x] ? x : fa[x] = find(fa[x]); };
auto unite = [&](int x, int y) { fa[find(x)] = find(y); };
```

### 线段树

```cpp
void refresh(int p) { sum[p] = sum[ls[p]] + sum[rs[p]]; }
//动态开点权值线段树
int add(int p, int l, int r, int x) {
    if (!p) p = ++tot;
    if (l == r) { sum[p]++; return p; }
    int mid = (l + r) >> 1;
    x <= mid ? ls[p] = add(ls[p], l, mid, x) : rs[p] = add(rs[p], mid + 1, r, x);
    return refresh(p), p;
}
//线段树合并
int merge(int p, int q, int l, int r) {
    if (!p) return q;
    if (!q) return p;
    if (l == r) { sum[p] += sum[q]; return p; }
    int mid = (l + r) >> 1;
    if (l <= mid) ls[p] = merge(ls[p], ls[q], l, mid);
    if (r > mid) rs[p] = merge(rs[p], rs[q], mid + 1, r);
    return refresh(p), p;
}
//线段树上二分
int ask(int p, int l, int r, int k) {
    if (l == r) return l;
    int mid = (l + r) >> 1;
    if (sum[ls[p]] >= k) return ask(ls[p], l, mid, k);
    else return ask(rs[p], mid + 1, r, k - sum[ls[p]]);
}
```

### 吉老师线段树

```cpp
struct seg {
	int l, r;
	int vmax, vmaxh, vsec, mxcnt, sum;
	int tmax, telse, tmaxh, telseh;
} t[MAXN << 2];
#define ls p << 1
#define rs p << 1 | 1
#define mid ((t[p].l + t[p].r) >> 1)

void refresh(int p) {
	t[p].vmax = max(t[ls].vmax, t[rs].vmax);
	t[p].vmaxh = max(t[ls].vmaxh, t[rs].vmaxh);
	t[p].sum = t[ls].sum + t[rs].sum;
	if (t[ls].vmax == t[rs].vmax) {
		t[p].vsec = max(t[ls].vsec, t[rs].vsec);
		t[p].mxcnt = t[ls].mxcnt + t[rs].mxcnt;
	} else if (t[ls].vmax > t[rs].vmax) {
		t[p].vsec = max(t[ls].vsec, t[rs].vmax);
		t[p].mxcnt = t[ls].mxcnt;
	} else {
		t[p].vsec = max(t[ls].vmax, t[rs].vsec);
		t[p].mxcnt = t[rs].mxcnt;
	}
}
void build(int p, int l, int r) {
	t[p].l = l, t[p].r = r;
	if (l == r) {
		cin >> t[p].vmaxh;
		t[p].sum = t[p].vmax = t[p].vmaxh;
		t[p].vsec = -INF;
		t[p].mxcnt = 1;
		return;
	}
	build(ls, l, mid), build(rs, mid + 1, r);
	refresh(p);
}
void addtag(int p, int tmax, int tmaxh, int telse, int telseh) {
	t[p].sum += t[p].mxcnt * tmax + (t[p].r - t[p].l + 1 - t[p].mxcnt) * telse;
	t[p].vmaxh = max(t[p].vmaxh, t[p].vmax + tmaxh);
	t[p].tmaxh = max(t[p].tmaxh, t[p].tmax + tmaxh);
	t[p].telseh = max(t[p].telseh, t[p].telse + telseh);
	t[p].vmax += tmax;
	t[p].tmax += tmax;
	t[p].telse += telse;
	if (t[p].vsec != -INF) t[p].vsec += telse;
}
void pushdown(int p) {
	int tmp = max(t[ls].vmax, t[rs].vmax);
	if (t[ls].vmax == tmp) addtag(ls, t[p].tmax, t[p].tmaxh, t[p].telse, t[p].telseh);
	else addtag(ls, t[p].telse, t[p].telseh, t[p].telse, t[p].telseh);
	if (t[rs].vmax == tmp) addtag(rs, t[p].tmax, t[p].tmaxh, t[p].telse, t[p].telseh);
	else add_tag(rs, t[p].telse, t[p].telseh, t[p].telse, t[p].telseh);
	t[p].tmax = t[p].tmaxh = t[p].telse = t[p].telseh = 0;
}
void cadd(int p, int l, int r, int val) {
	if (l <= t[p].l && t[p].r <= r) {
		add_tag(p, val, val, val, val);
		return;
	}
	pushdown(p);
	if (l <= mid) cadd(ls, l, r, val);
	if (r > mid) cadd(rs, l, r, val);
	refresh(p);
}
void cmin(int p, int l, int r, int val) {
	if (val >= t[p].vmax) return;
	if (l <= t[p].l && t[p].r <= r && val > t[p].vsec) {
		add_tag(p, val - t[p].vmax, val - t[p].vmax, 0, 0);
		return;
	}
	pushdown(p);
	if (l <= mid) cmin(ls, l, r, val);
	if (r > mid) cmin(rs, l, r, val);
	refresh(p);
}
int qsum(int p, int l, int r) {
	if (l <= t[p].l && t[p].r <= r) return t[p].sum;
	pushdown(p);
	int ans = 0;
	if (l <= mid) ans += qsum(ls, l, r);
	if (r > mid) ans += qsum(rs, l, r);
	return ans;
}
int qmax(int p, int l, int r) {
	if (l <= t[p].l && t[p].r <= r) return t[p].vmax;
	pushdown(p);
	int ans = -INF;
	if (l <= mid) ans = max(ans, qmax(ls, l, r));
	if (r > mid) ans = max(ans, qmax(rs, l, r));
	return ans;
}
int qmaxh(int p, int l, int r) {
	if (l <= t[p].l && t[p].r <= r) return t[p].vmaxh;
	pushdown(p);
	int ans = -INF;
	if (l <= mid) ans = max(ans, qmaxh(ls, l, r));
	if (r > mid) ans = max(ans, qmaxh(rs, l, r));
	return ans;
}
```

### 李超树

```cpp
struct seg {
    double k, b;
    void init(int x0, int y0, int x1, int y1) {
        if (x0 == x1) k = 0, b = max(y0, y1);
        else k = 1.0 * (y1 - y0) / (x1 - x0), b = y0 - k * x0;
    }
} a[MAXN];
struct node {
    int id;
    double y;
    bool operator<(const node &a) const {
        if (fabs(y - a.y) < eps) return id > a.id;
        return y < a.y;
    }
};
//第id条直线在x处取值
double f(int id, int x) { return a[id].k * x + a[id].b; }
struct segtree {
    int l, r, id;
} t[MAXN << 2];
#define ls p << 1
#define rs p << 1 | 1
#define mid ((t[p].l + t[p].r) >> 1)
//建树build(1,1,n)
void build(int p, int l, int r) {
    t[p].l = l, t[p].r = r;
    if (l == r) return;
    build(ls, l, mid), build(rs, mid + 1, r);
}
void change(int p, int l, int r, int id) {
    if (l <= t[p].l && t[p].r <= r) {
        if (!t[p].id) { t[p].id = id; return; }
        if (t[p].l == t[p].r) {
            if (f(id, t[p].l) > f(t[p].id, t[p].l)) t[p].id = id;
            return;
        }
        if (fabs(a[id].k - a[t[p].id].k) < eps) {
            if (f(id, mid) > f(t[p].id, mid)) t[p].id = id;
        } else if (a[id].k > a[t[p].id].k) {
            if (f(id, mid) > f(t[p].id, mid)) change(ls, l, r, t[p].id), t[p].id = id;
            else change(rs, l, r, id);
        } else {
            if (f(id, mid) > f(t[p].id, mid)) change(rs, l, r, t[p].id), t[p].id = id;
            else change(ls, l, r, id);
        }
        return;
    }
    if (l <= mid) change(ls, l, r, id);
    if (r > mid) change(rs, l, r, id);
}
node query(int p, int x) {
    node ans = {t[p].id, f(t[p].id, x)};
    if (t[p].l == t[p].r) return ans;
    return max(x <= mid ? query(ls, x) : query(rs, x), ans);
}
//查询x=k处编号最大的线段（可修改node operator）
query(1, x).id;
//插入线段
if (x0 > x1) swap(x0, x1), swap(y0, y1);
a[++cnt].init(x0, y0, x1, y1), change(1, x0, x1, cnt);
```

### 主席树 (静态区间第k小)

```cpp
struct segtree { int ls, rs, sum; } t[N << 5];
// 初始化rt[0]=build(1,n)为离散化后的值域
int build(int l, int r) {
	int p = ++tot;
	if (l == r) return p;
	int mid = l + r >> 1;
	t[p].ls = build(l, mid), t[p].rs = build(mid + 1, r);
	return root;
}
//插入rt[i]=update(x,1,n,rt[i-1])
int update(int x, int l, int r, int o) {
	int p = ++tot;
	t[p] = t[o], t[p].sum++;
	if (l == r) return dir;
	int mid = l + r >> 1;
	if (x <= mid) t[p].ls = update(x, l, mid, t[p].ls);
	else t[p].rs = update(x, mid + 1, r, t[p].rs);
	return p;
}
//查询query(rt[l-1],rt[r],1,n,k)
int query(int p, int q, int l, int r, int k) {
	if (l == r) return l;
	int mid = l + r >> 1, x = t[t[q].ls].sum - t[t[p].ls].sum;
	if (k <= x) return query(t[p].ls, t[q].ls, l, mid, k);
	else return query(t[p].rs, t[q].rs, mid + 1, r, k - x);
}
```



### FHQ平衡树

```cpp
struct node {
    int siz, ls, rs, prio, val;
} t[N];
int n, m, tot, root;
int newnode(int k) {
    t[++tot] = (node){1, 0, 0, rand(), k};
    return tot;
}
void refresh(int p) { t[p].siz = t[t[p].ls].siz + t[t[p].rs].siz + 1; }
int merge(int p, int q) {
    if (!p || !q) return p + q;
    if (t[p].prio < t[q].prio) {
        t[p].rs = merge(t[p].rs, q);
        refresh(p);
        return p;
    }
    t[q].ls = merge(p, t[q].ls);
    refresh(q);
    return q;
}
void split(int now, int k, int &p, int &q) {
    if (!now) { p = q = 0; return; }
    if (t[now].val <= k) p = now, split(t[now].rs, k, t[now].rs, q);
    else q = now, split(t[now].ls, k, p, t[now].ls);
    refresh(now);
}
int kth(int p, int k) {
    while (1) {
        if (k <= t[t[p].ls].siz) p = t[p].ls;
        else if (k > t[t[p].ls].siz + 1) k -= t[t[p].ls].siz + 1, p = t[p].rs;
        else return p;
    }
}
//插入整数k
split(root, k, x, y);
root = merge(merge(x, newnode(k)), y);
//删除整数k
split(root, k, x, z);
split(x, k - 1, x, y);
y = merge(t[y].ls, t[y].rs);
root = merge(merge(x, y), z);
//查询k的排名
split(root, k - 1, x, y);
ans = t[x].siz + 1;
root = merge(x, y);
//查询排名为k的数（如果不存在，则认为是排名小于k的最大数）
ans = t[kth(root, k)].val;
//k的前驱
split(root, k - 1, x, y);
ans = t[kth(x, t[x].siz)].val;
root = merge(x, y);
//k的后继
split(root, k, x, y);
ans = t[kth(y, 1)].val;
root = merge(x, y);
```

### 文艺平衡树 (splay)

```cpp
struct Splay {
    int val[N], ch[N][2], fa[N], tag[N], siz[N], root;
    Splay() { root = 0; }
    bool chk(int x) { return ch[fa[x]][1] == x; }
    void refresh(int x) { siz[x] = siz[ch[x][0]] + siz[ch[x][1]] + 1; }
    void pushdown(int x) {
        if (tag[x]) {
            if (ch[x][0]) {
                tag[ch[x][0]] ^= 1;
                swap(ch[ch[x][0]][0], ch[ch[x][0]][1]);
            }
            if (ch[x][1]) {
                tag[ch[x][1]] ^= 1;
                swap(ch[ch[x][1]][0], ch[ch[x][1]][1]);
            }
            tag[x] = 0;
        }
    }
    void rotate(int now) {
        int f = fa[now], gf = fa[f], k = chk(now), w = ch[now][k ^ 1];
        fa[w] = f, ch[f][k] = w;
        fa[now] = gf, ch[gf][chk(f)] = now;
        fa[f] = now, ch[now][k ^ 1] = f;
        refresh(f), refresh(now);
    }
    void splay(int now, int goal = 0) {
        while (fa[now] != goal) {
            int f = fa[now], gf = fa[f];
            if (gf != goal) {
                if (chk(f) == chk(now)) rotate(f);
                else rotate(now);
            }
            rotate(now);
        }
        if (!goal) root = now;
    }
    int kth(int k) {
        int now = root;
        while (1) {
            pushdown(now);
            if (ch[now][0] && siz[ch[now][0]] >= k) now = ch[now][0];
            else if (siz[ch[now][0]] + 1 < k) {
                k -= siz[ch[now][0]] + 1;
                now = ch[now][1];
            } else return now;
        }
    }
    int split(int l, int r) {
        l = kth(l), r = kth(r + 2);
        splay(l), splay(r, l);
        return ch[r][0];
    }
    void reverse(int l, int r) {
        int now = split(l, r);
        tag[now] ^= 1;
        swap(ch[now][0], ch[now][1]);
    }
    int query(int x) { return val[split(x, x)]; }
    int build(int l, int r, int f) {
        if (l > r) return 0;
        int mid = (l + r) >> 1;
        val[mid] = a[mid]， fa[mid] = f;
        ch[mid][0] = build(l, mid - 1, mid);
        ch[mid][1] = build(mid + 1, r, mid);
        refresh(mid);
        return mid;
    }
} T;
int main() {
    cin >> n >> m;
    a[1] = a[n + 2] = -INF;
    for (int i = 2; i <= n + 1; ++i) a[i] = i - 1;
    T.root = T.build(1, n + 2, 0);
    while (m--) {
        int l = read(), r = read();
        T.reverse(l, r);
    }
    for (int i = 1; i <= n; ++i) pintf("%d ", T.query(i));
}
```

### 左偏树

```cpp
struct heap { int val, fa, d, ch[2]; }t[N];
//dist大的视为左儿子，dist小的视为右儿子
int& rs(int x) { return t[x].ch[t[t[x].ch[1]].d < t[t[x].ch[0]].d]; }
//插入时直接merge一个新节点,O(log n)复杂度
int merge(int x, int y) {
  if (!x || !y) return x | y;
  if (t[x].val < t[y].val) swap(x, y);
  int& rs_ref = rs(x);
  rs_ref = merge(rs_ref, y);
  t[rs_ref].fa = x;
  t[x].d = t[rs(x)].d + 1;
  return x;
}
//配合并查集维护是否在一个堆里（旧根指向新根），用数组记录是否被删除
//合并代码
x = find(x), y = find(y);
if (x != y) rt[x] = rt[y] = merge(x, y);
//删除根
x = find(x);
rt[x] = rt[t[x].ch[0]] = rt[t[x].ch[1]] = merge(t[x].ch[0], t[x].ch[1]);
// 由于堆中的点会 find 到 x，所以 rt[x] 也要修改
```

### ST表

```cpp
void pre() {  //预处理log
  Logn[1] = 0, Logn[2] = 1;
  for (int i = 3; i < maxn; i++) Logn[i] = Logn[i / 2] + 1;
}
// 求静态区间[l,r]最大值/最小值/gcd
for (int i = 1; i <= n; i++) f[i][0] = read();
pre();
//根据maxn确定logn=21,22..
for (int j = 1; j <= logn; j++)
for (int i = 1; i + (1 << j) - 1 <= n; i++)
  f[i][j] = max(f[i][j - 1], f[i + (1 << (j - 1))][j - 1]);
for (int i = 1; i <= m; i++) {
    int x = read(), y = read();
    int s = Logn[y - x + 1];
    printf("%d\n", max(f[x][s], f[y - (1 << s) + 1][s]));
}
```

### 树状数组

```cpp
#define lowbit(x) (x)&(-x)
void add(int x, int k) { for (x <= n; x += lowbit(x)) c[x] += k; } //单点加
int ask(int x) {  // a[1]..a[x]的和，用前缀和实现区间查询
  int ans = 0;
  for (x > 0; x -= lowbit(x)) ans += c[x];
  return ans;
}
```

### 莫队

```cpp
//如果能把当前区间[l,r]的答案转移到与其相邻的区间，则可以实现O(n\sqrt(n))的复杂度
void add(int x);
void del(int x);

void solve() {
  BLOCK_SIZE = int(ceil(pow(n, 0.5)));
  //l为第一关键字，r为第二关键字从小到大排序
  sort(querys, querys + m);
  for (int i = 0; i < m; ++i) {
    const query &q = querys[i];
    while (l > a[i].l) add(c[--l]);
    while (r < a[i].r) add(c[++r]);
    while (l < a[i].l) del(c[l++]);
    while (r > a[i].r) del(c[r--]);
    ans[q.id] = res;
  }
}
```

### 树链剖分

```cpp
void dfs1(int u) {
    dep[u] = dep[fa[u]] + 1;
    siz[u] = 1;
    int mxson = -1;
    for (int i = h[u]; i; i = e[i].nxt) {
        int v = e[i].v;
        if (v == fa[u]) continue;
        fa[v] = u;
        dfs1(v);
        siz[u] += siz[v]; //子树大小
        if (!son[u] || siz[v] > mxson) {
            son[u] = v; //重儿子
            mxson = siz[v];
        }
    }
}
void dfs2(int u, int topf) {
    b[id[u] = ++dfn] = a[u];
    top[u] = topf; //链顶
    if (!son[u]) return;
    dfs2(son[u], topf);
    for (int i = h[u]; i; i = e[i].nxt) {
        int v = e[i].v;
        if (v != fa[u] && v != son[u]) dfs2(v, v);
    }
}
//树上路径+v
void tadd(int x, int y, int v) {
    while (top[x] != top[y]) {
        if (dep[top[x]] < dep[top[y]]) swap(x, y);
        add(1, id[top[x]], id[x], v);
        x = fa[top[x]];
    }
    if (dep[x] < dep[y]) swap(x, y);
    add(1, id[y], id[x], v);
}
//树上路径求和
int tquery(int x, int y) {
    int ans = 0;
    while (top[x] != top[y]) {
        if (dep[top[x]] < dep[top[y]]) swap(x, y);
        ans += query(1, id[top[x]], id[x]);
        x = fa[top[x]];
    }
    if (dep[x] < dep[y]) swap(x, y);
    ans += query(1, id[y], id[x]);
    return ans;
}
//子树内+v
add(1, id[x], id[x] + siz[x] - 1, k);
//子树内求和
query(1, id[x], id[x] + siz[x] - 1)；
```

### 原地堆化

```java
void heapify(int[] h) { // 从下往上遍历非叶子节点
    for (int i = h.length / 2 - 1; i >= 0; i--)  sink(h, i);
}
void sink(int[] h, int i) {
    int n = h.length;
    while (2 * i + 1 < n) {
        int j = 2 * i + 1; // 左儿子（因为起始下标为0所以+1）
        if (j + 1 < n && h[j + 1] > h[j]) j++; // 跟左右儿子中较大的交换
        if (h[i] >= h[j])  break; // 无法交换，递归结束
        swap(h, i, j), i = j; // 交换后继续向下递归
    }
}
```