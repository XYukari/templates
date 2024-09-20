## 杂项


### 离散化

```cpp
for (int i = 1; i <= n; i++) 
    lsh[i] = a[i] = read();
sort(lsh + 1, lsh + n + 1);
int m = unique(lsh + 1, lsh + n + 1) - lsh - 1;
for (int i = 1; i <= n; i++) 
    a[i] = lower_bound(lsh + 1, lsh + m + 1, a[i]) - lsh;

```

### 公式

```cpp
C(n,m)%2 = (n&m)==m
约瑟夫问题：
F(n,m) = 有 n 个人 (0,1,2,..,n−1)，每次杀掉编号为(x + m) % n的人，最终的幸存者。
F(n,m) = (F(n − 1, m) + m) % n
```