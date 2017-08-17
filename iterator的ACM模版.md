# Iterator Code Lib

## 实用数据结构

### 加权并查集

>解决集合问题中，集合内元素有关系并且关系具有传递性的问题。</br>
从集合中删除节点的方法：消除该点对集合的影响(如集合中的点个数、和、最值)，然后给它分配一个新的编号(原来的编号不管)

#### 头文件&宏&全局变量

```c++
#define MAXN 100000//最大点数
int p[MAXN];//父节点
int v[MAXN];//到父节点边的权值(加权解决集合中点的相互关系的问题)
```

#### 初始化

```c++
void union_init(int minn,int maxn)//编号最小值到最大值
{
    for(int i=minn;i<=maxn;i++)
    {
        p[i]=i;
        v[i]=0;
    }
}
```

#### 查找

>执行后p[a]为a所在集合的根节点,v[a]为a到其集合的根节点的权值

```c++
void union_find(int a)
{
    if(p[a]==a)return;
    union_find(p[a]);
    v[a]+=v[p[a]];
    p[a]=p[p[a]];
}
```

#### 合并

>合并a,b所在集合,vab为 a到b的权值

```c++
void union_merge(int a,int b,int vab)
{
    union_find(a);
    union_find(b);
    v[p[a]]=vab-v[a]+v[b];
    p[p[a]]=p[b];
}
```

### 树状数组

>要求所有数的和不能超出范围,也可修改为记录最值

#### 精准覆盖

```c++
struct DLX {
    int n, m, si; //n行数m列数si目前有的节点数
    //十字链表组成部分
    int U[MNN], D[MNN], L[MNN], R[MNN], Row[MNN], Col[MNN];
    //第i个结点的U向上指针D下L左R右，所在位置Row行Col列
    int H[MN], S[MM]; //记录行的选择情况和列的覆盖情况
    int ansd, ans[MN];
    void init(int _n, int _m) //初始化空表
    {
        n = _n;
        m = _m;
        for (int i = 0; i <= m; i++) //初始化第一横行（表头）
        {
            S[i] = 0;
            U[i] = D[i] = i; //目前纵向的链是空的
            L[i] = i - 1;
            R[i] = i + 1; //横向的连起来
        }
        R[m] = 0;
        L[0] = m;
        si = m; //目前用了前0~m个结点
        for (int i = 1; i <= n; i++)
            H[i] = -1;
    }
    void link(int r, int c) //插入点(r,c)
    {
        ++S[Col[++si] = c]; //si++;Col[si]=c;S[c]++;
        Row[si] = r;
        D[si] = D[c];
        U[D[c]] = si;
        U[si] = c;
        D[c] = si;
        if (H[r] < 0)
            H[r] = L[si] = R[si] = si;
        else {
            R[si] = R[H[r]];
            L[R[H[r]]] = si;
            L[si] = H[r];
            R[H[r]] = si;
        }
    }
    void remove(int c) //列表中删掉c列
    {
        L[R[c]] = L[c]; //表头操作
        R[L[c]] = R[c];
        for (int i = D[c]; i != c; i = D[i])
            for (int j = R[i]; j != i; j = R[j]) {
                U[D[j]] = U[j];
                D[U[j]] = D[j];
                --S[Col[j]];
            }
    }
    void resume(int c) //恢复c列
    {
        for (int i = U[c]; i != c; i = U[i])
            for (int j = L[i]; j != i; j = L[j])
                ++S[Col[U[D[j]] = D[U[j]] = j]];
        L[R[c]] = R[L[c]] = c;
    }

    // HUST-1017
    //仅仅判断有无解
    bool dance(int d) //选取了d行
    {
        if (R[0] == 0) //全部覆盖了
        {
            //全覆盖了之后的操作
            ansd = d;
            return 1;
        }
        int c = R[0];
        for (int i = R[0]; i != 0; i = R[i])
            if (S[i] < S[c])
                c = i;
        remove(c);
        for (int i = D[c]; i != c; i = D[i]) {
            ans[d] = Row[i];
            for (int j = R[i]; j != i; j = R[j])
                remove(Col[j]);
            if (dance(d + 1))
                return 1;
            for (int j = L[i]; j != i; j = L[j])
                resume(Col[j]);
        }
        resume(c);
        return 0;
    }

    //ZOJ 3209
    //求最小解
    void danceLeast(int d) //选取了d行
    {
        if (ansd != -1 && ansd <= d) return;
        if (R[0] == 0) //全部覆盖了
        {
            //全覆盖了之后的操作
            if (ansd == -1)
                ansd = d;
            else if (d < ansd)
                ansd = d;
            return;
        }
        int c = R[0];
        for (int i = R[0]; i != 0; i = R[i])
            if (S[i] < S[c])
                c = i;
        remove(c);
        for (int i = D[c]; i != c; i = D[i]) {
            ans[d] = Row[i];
            for (int j = R[i]; j != i; j = R[j])
                remove(Col[j]);
            danceLeast(d + 1);
            for (int j = L[i]; j != i; j = L[j])
                resume(Col[j]);
        }
        resume(c);
    }
} dlx;
```

#### 重复覆盖

```c++
struct DLX {//成员变量，init(),link()同上
    void remove(int c) //列表中删掉c列
    {
        for(int i = D[c];i != c;i = D[i])
             L[R[i]] = L[i], R[L[i]] = R[i];
    }
    void resume(int c) //恢复c列
    {
        for(int i = U[c];i != c;i = U[i])
             L[R[i]]=R[L[i]]=i;
    }

    bool v[MNN];
    int f()
    {
        int ret = 0;
        for(int c = R[0];c != 0;c = R[c])v[c] = true;
        for(int c = R[0];c != 0;c = R[c])
            if(v[c])
            {
                ret++;
                v[c] = false;
                for(int i = D[c];i != c;i = D[i])
                    for(int j = R[i];j != i;j = R[j])
                        v[Col[j]] = false;
            }
        return ret;

    }

    //HDU 2295
    //是否有解
    bool dance(int d)
    {
        if(d + f() > K)return false; //此处K为题目要求最多能选择的数量
        if(R[0] == 0)return d <= K;
        int c = R[0];
        for(int i = R[0];i != 0;i = R[i])
            if(S[i] < S[c])
                c = i;
        for(int i = D[c];i != c;i = D[i])
        {
            remove(i);
            for(int j = R[i];j != i;j = R[j])remove(j);
            if(Dance(d+1))return true;
            for(int j = L[i];j != i;j = L[j])resume(j);
            resume(i);
        }
        return false;
    }

    //FZU 1686
    //求最小解
    void danceLeast(int d) {
        if (d + f() >= ansd) return;
        if (R[0] == 0) //全部覆盖了
        {
            //全覆盖了之后的操作
            if (d < ansd)
                ansd = d;
            return;
        }
        int c = R[0];
        for (int i = R[0]; i != 0; i = R[i])
            if (S[i] < S[c])
                c = i;
        for (int i = D[c]; i != c; i = D[i]) {
            remove(i);
            for (int j = R[i]; j != i; j = R[j])
                remove(j);
            danceLeast(d + 1);
            for (int j = L[i]; j != i; j = L[j])
                resume(j);
            resume(i);
        }
    }
} dlx;
```

## 数论

### 扩展欧几里得

#### 定义

> 对于不完全为 0 的非负整数 a，b，gcd（a，b）表示 a，b 的最大公约数，必然存在整数对 x，y ，使得 gcd（a，b）=ax+by。

#### 代码

```c++
int exgcd(int a,int b,int &x,int &y){
    if (b==0){
        x=1,y=0;
        return a;
    }
    int q=exgcd(b,a%b,y,x);
    y-=a/b*x;
    return q;
}
```

#### 求逆元

>求a对b的逆元，即(a^(-1))mod b </br>
int x,y;</br>
exgcd(a,b,x,y);</br>
x即为a对b的逆元

### 中国剩余定理

#### 定义&通式

>给出了以下的一元线性同余方程组：</br>
$$
\left ( S \right ) :
\left\{
\begin{matrix}
x \equiv a_1 \left ( mod\ m_1 \right )\\
x \equiv a_2 \left ( mod\ m_2 \right )\\
\vdots \\
x \equiv a_n \left ( mod\ m_n \right )
\end{matrix}
\right.
$$ </br>
有解的判定条件，并用构造法给出了在有解情况下解的具体形式。</br>
中国剩余定理说明：假设整数$m_1,m_2, \cdots ,m_n$两两互质，则对任意的整数：$a1,a2, \cdots ,an$，方程组 有解，并且通解可以用如下方式构造得到：</br>
设</br>
$$ M = m_1 \times m_2 \times m_3 \times \cdots \times m_n = \prod_{i=1}^n m_i $$ </br>
是整数$m_1,m_2, \cdots ,m_n$的乘积，并设</br>
$$ M_i = M \div m_i \ , \forall i \in $$