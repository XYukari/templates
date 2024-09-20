## 计算几何

### 必有的
```cpp
const double eps = 1e-10;
const long double PI = acos(-1.0);
```

### Point、Vector类
```cpp
struct Point {
    int x, y;
    Point(int x = 0, int y = 0) : x(x), y(y) {}
};
ostream& operator<<(ostream &os, const Point &p) {return os << p.x << " " << p.y;}
typedef Point Vector;

Vector operator+(const Vector &p1, const Vector &p2) {return Vector(p1.x + p2.x, p1.y + p2.y);}
Vector operator-(const Vector &p1, const Vector &p2) {return Vector(p1.x - p2.x, p1.y - p2.y);}
bool operator==(const Point &p1, const Point &p2) {return p1.x == p2.x && p1.y == p2.y;}
bool operator<(const Point &p1, const Point &p2) {return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);}
int Cross(const Vector &a, const Vector &b) {return a.x * b.y - a.y * b.x;}
double Dist(const Point &a, const Point &b) {
    return sqrt(1.0 * (a.x - b.x) * (a.x - b.x) + 1.0 * (a.y - b.y) * (a.y - b.y));
}
Vector Rotate(const Vector &a, const double rad) {
    return Vector(a.x * cos(rad) - a.y * sin(rad), a.x * sin(rad) + a.y * cos(rad));
}
int Dot(const Vector &a, const Vector &b) { return a.x * b.x + a.y * b.y;}
double Length(const Vector& a) {return sqrt(1.0 * Dot(a, a)); }
double Angle(const Vector &a, const Vector &b) {
    return asin(Cross(a, b) / Length(a) / Length(b));
}
int dcmp(long double a) {if (fabs(a) <= eps) return 0; return a > 0 ? 1 : -1;}
int Area(const Point &a, const Point &b, const Point &c) {
    Vector u = a - b, v = a - c;
    return abs(Cross(u, v));
}
```
### 求凸包
```cpp
vector<Point> ConvexHell(vector<Point> &p) {
    sort(p.begin(), p.end());
    p.erase(unique(p.begin(), p.end()), p.end());
    int n = p.size(), m = 0;
    vector<Point> stk(n + 1);
    for (int i = 0;i < n;i++) {
        while (m > 1 && Cross(stk[m - 1] - stk[m - 2], p[i] - stk[m - 2]) <= 0) m--;
        stk[m++] = p[i];
    }
    int k = m;
    for (int i = n - 2;i >= 0;i--) {
        while (m > k && Cross(stk[m - 1] - stk[m - 2], p[i] - stk[m - 2]) <= 0) m--;
        stk[m++] = p[i];
    }
    if (n > 1) m--;
    stk.resize(m);
    return stk;
}
```

### 旋转卡壳求凸包的直径
```cpp
int diameter2(vector<Point> &points) {
    auto p = ConvexHell(points);
    int n = p.size();
    if (n == 1) return 0;
    else if (n == 2) return Dist2(p[0], p[1]);
    p.push_back(p[0]);
    int ans = 0;
    for (int u = 0, v = 1;u < n;u++) {
        for (;;) {
            int diff = Cross(p[u + 1] - p[u], p[v + 1] - p[v]);
            if (diff <= 0) {
                ans = max(ans, Dist2(p[u], p[v]));
                if (diff == 0) ans = max(ans, Dist2(p[u], p[v + 1]));
                break;
            }
            v = (v + 1) % n;
        }
    }
    return ans;
}
```