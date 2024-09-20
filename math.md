## 数论

### 扩展中国剩余定理
```cpp
template<class T> 
T exgcd(T a, T b, T &x, T &y) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    T d = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}

template<class T> 
T mul(T a, T b, T mod) {
	a = (a % mod + mod) % mod;
	b = (b % mod + mod) % mod;
	T ans = 0;
	while (b != 0) {
		if ((b & 1) != 0) {
			ans = (ans + a) % mod;
		}
		a = (a + a) % mod;
		b >>= 1;
	}
	return ans;
}

template<class T> 
T excrt(vector<T> &m, vector<T> &r) {
	T tail = 0, lcm = 1, tmp, b, c, x0, x, y, d;
	for (int i = 0;i < m.size();i++) {
		b = m[i]; c = ((r[i] - tail) % b + b) % b;
		d = exgcd(lcm, b, x, y);
		if (c % d != 0) return -1;
		x0 = mul(x, c / d, b / d);
		tmp = lcm * (b / d);
		tail = (tail + mul(x0, lcm, tmp)) % tmp;
		lcm = tmp;
	}
	return tail;
}
```
### exgcd
```cpp
template <class T> T sign(const T &a) {
    return a == 0 ? 0 : (a < 0 ? -1 : 1);
}
template <class T> T ceil(const T &a, const T &b) {
    T A = abs(a), B = abs(b);
    assert(b != 0);
    return sign(a) * sign(b) > 0 ? (A + B - 1) / B : -A / B;
}
template<class T> 
T exgcd(T a, T b, T &x, T &y) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    T d = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}
template<class T>
void exgcd(T a, T b, T c) {
	T x, y, d = exgcd(a, b, x, y);
	if (c % d != 0) {
		cout << "Impossible" << endl;
		return ;
	}
	x *= c / d, y *= c / d;
	T p = b / d, q = a / d, k;
	if (x < 0) {
		k = ceil(1 - x, p);
		x += p * k;
		y -= q * k;
	}
	else if (x >= 0) { //将x提高到最小正整数
		k = (x - 1) / p;
		x -= p * k; //将x降低到最小正整数
		y += q * k;
	}
	if (y > 0) { //有正整数解
		cout << (y - 1) / q + 1 << endl; //将y减到1的方案数即为解的个数
		cout << x << endl; //当前的x即为最小正整数x
		cout << (y - 1) % q + 1 << endl; //将y取到最小正整数
		cout << x + (y - 1) / q * p << endl; //将x提升到最大
		cout << y << endl; //特解即为y最大值      
	} else { //无整数解
		cout << x << endl; //当前的x即为最小的正整数x
		cout << y + q * ceil(1 - y, q) << endl; //将y提高到正整数
	}
}
```

### exlucas
```cpp
template<class T> 
T qpow(T a, T b, T mod) {
	T ans = 1;
	for (;b;b >>= 1, a = a * a % mod)
		if (b & 1) ans = ans * a % mod;
	return ans;
} 

template<class T> 
T exgcd(T a, T b, T &x, T &y) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    T d = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}

template<class T> 
T mul(T a, T b, T mod) {
	a = (a % mod + mod) % mod;
	b = (b % mod + mod) % mod;
	T ans = 0;
	while (b != 0) {
		if ((b & 1) != 0) {
			ans = (ans + a) % mod;
		}
		a = (a + a) % mod;
		b >>= 1;
	}
	return ans;
}

template<class T> 
T excrt(vector<T> &m, vector<T> &r) {
	T tail = 0, lcm = 1, tmp, b, c, x0, x, y, d;
	for (int i = 0;i < m.size();i++) {
		b = m[i]; c = ((r[i] - tail) % b + b) % b;
		d = exgcd(lcm, b, x, y);
		if (c % d != 0) return -1;
		x0 = mul(x, c / d, b / d);
		tmp = lcm * (b / d);
		tail = (tail + mul(x0, lcm, tmp)) % tmp;
		lcm = tmp;
	}
	return tail;
}

const int MAXN = 4e4 + 10;
i64 fact[MAXN];

template<class T>
T inv(T a, T mod) {
	return qpow(a, mod - 2, mod);
}

template<class T> 
T C(T n, T m, T mod) {
	if (n < m) return 0;
	return fact[n] * inv(fact[m], mod) % mod * inv(fact[n - m], mod) % mod;
}

template<class T>
T lucas(T n, T m, T mod) {
	if (n < m) return 0;
	if (!n) return 1;
	return lucas(n / mod, m / mod, mod) * C(n % mod, m % mod, mod) % mod;
}
```
### 高斯消元
```cpp
const double eps = 1e-7;
inline int dcmp(double a, double b) {
    if (fabs(a - b) < eps) return 0;
    return a > b ? 1 : -1;
}
//最普通的高斯消元
void gauss(vector<vector<double>> &matrix) {
    int n = matrix.size();
    for (int i = 0;i < n;i++) {
        int r = i;
        for (int j = i + 1;j < n;j ++) 
            if (dcmp(fabs(matrix[r][i]), fabs(matrix[j][i])) == -1) r = j; 
        if (dcmp(matrix[r][i], 0.0) == 0) {
            printf("No Solution\n");
            return ;
        } 
        if (i != r) swap(matrix[i], matrix[r]);
        double div = matrix[i][i];
        for (int j = i;j <= n;j++) matrix[i][j] /= div;
        for (int j = i + 1;j < n;j++) {
            div = matrix[j][i];
            for (int k = i;k <= n;k++) matrix[j][k] -= matrix[i][k] * div;
        }
    }
    vector<double> ans(n);
    ans[n - 1] = matrix[n - 1][n];
    for (int i = n - 2;i >= 0;i--) {
        ans[i] = matrix[i][n];
        for (int j = i + 1;j < n;j++) ans[i] -= matrix[i][j] * ans[j];
    }
    for (int i = 0;i < n;i++) printf("%.2lf\n", ans[i]);
}
//区分无解，多解和唯一解
void gauss2(vector<vector<double>> &matrix) {   
    int n = matrix.size();
    for (int i = 0;i < n;i++) {
        int mx = i;
        for (int j = 0;j < n;j++) {
            if (j < i && fabs(matrix[j][j]) >= eps) continue;
            if (dcmp(fabs(matrix[j][i]), fabs(matrix[mx][i])) == 1) mx = j;
        }
        if (mx != i) swap(matrix[mx], matrix[i]);
        if (fabs(matrix[i][i]) >= eps) {
            double div = matrix[i][i];
            for (int j = i;j <= n;j++) matrix[i][j] /= div; 
            for (int j = 0;j < n;j++) {
                if (i == j) continue;
                div = matrix[j][i] / matrix[i][i];
                for (int k = i;k <= n;k++) matrix[j][k] -= div * matrix[i][k];
            }
        }
    }
    bool flag = true;
    for (int i = 0;i < n;i++) {
        if ((fabs(matrix[i][i]) < eps) && (fabs(matrix[i][n]) >= eps)) {
            printf("-1\n");
            return ;
        }
        if (fabs(matrix[i][i]) < eps) flag = false;
    }
    if (!flag) printf("0\n");
    else {
        for (int i = 0;i < n;i++) printf("x%d=%.2lf\n", i + 1, matrix[i][n]);
    }
}
// 异或方程组
void Gauss3(vector<vector<int>> &matrix) {
    int n = matrix.size();
    for (int i = 0;i < n;i++) {
        int is = i;
        for (int j = 0;j < n;j++) {
            if (j < i && matrix[j][j] == 1) continue;
            if (matrix[j][i] == 1) {
                swap(matrix[i], matrix[j]);
                break;
            }
        }
        if (matrix[i][i] != 1) continue;
        for (int j = 0;j < n;j++) {
            if (i == j) continue;
            if (matrix[j][i] == 1) {
                for (int k = i;k <= n;k++) matrix[j][k] ^= matrix[i][k];
            }
        }
    }
}
//同余方程组（要求mod相同）
const int mod = 1e9 + 7;
void Gauss4(vector<vector<int>> &matrix) {
    int n = matrix.size(), m = matrix[0].size();
    for (int i = 1;i <= n;i++) {
        for (int j = 1;j <= n;j++) {
            if (j < i && matrix[j][j] != 0) continue;
            if (matrix[j][i] != 0) {
                swap(matrix[j], matrix[i]);
                break;
            }
        }
        if (matrix[i][i] == 0) break;
        for (int j = 1;j <= n;j++) {
            if (i == j || matrix[j][i] == 0) continue;
            int gcdNum = __gcd(matrix[i][i], matrix[j][i]);
            int a = matrix[i][i] / gcdNum, b = matrix[j][i] / gcdNum;
            if (j < i && matrix[j][j] != 0) {
                for (int k = j;k < i;k++) matrix[j][k] = (matrix[j][k] * a) % mod;
            }
            for (int k = i;k <= n + 1;k++) 
                matrix[j][k] = (matrix[j][k] * a % mod - matrix[i][k] * b % mod + mod) % mod;
        }
    } 
}
```

### 欧拉函数
```cpp
constexpr int N = 1E7;
constexpr int P = 1000003;

bool isprime[N + 1];
int phi[N + 1];
std::vector<int> primes;

std::fill(isprime + 2, isprime + N + 1, true);
phi[1] = 1;
for (int i = 2; i <= N; i++) {
    if (isprime[i]) {
        primes.push_back(i);
        phi[i] = i - 1;
    }
    for (auto p : primes) {
        if (i * p > N) {
            break;
        }
        isprime[i * p] = false;
        if (i % p == 0) {
            phi[i * p] = phi[i] * p;
            break;
        }
        phi[i * p] = phi[i] * (p - 1);
    }
}


// 求单个数的欧拉函数
int phi(int n) {
    int res = n;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0) {
                n /= i;
            }
            res = res / i * (i - 1);
        }
    }
    if (n > 1) {
        res = res / n * (n - 1);
    }
    return res;
}
```