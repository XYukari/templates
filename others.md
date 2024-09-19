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