# Iterator Code Lib

<!-- vim-markdown-toc GFM -->
* [实用数据结构](#实用数据结构)
	- [加权并查集](#加权并查集)
		+ [头文件&宏&全局变量](#头文件宏全局变量)
		+ [初始化](#初始化)
		+ [查找](#查找)
		+ [合并](#合并)
	- [树状数组](#树状数组)
		+ [精准覆盖](#精准覆盖)
		+ [重复覆盖](#重复覆盖)
	- [splay](#splay)
		+ [头文件&宏&全局变量&结构体](#头文件宏全局变量结构体)
		+ [辅助函数](#辅助函数)
		+ [初始化](#初始化-1)
		+ [可使用函数](#可使用函数)
		+ [使用方法](#使用方法)
	- [矩阵快速幂](#矩阵快速幂)
		+ [代码](#代码)
* [数论](#数论)
	- [扩展欧几里得](#扩展欧几里得)
		+ [定义](#定义)
		+ [代码](#代码-1)
		+ [求逆元](#求逆元)
	- [中国剩余定理](#中国剩余定理)
		+ [定义&通式](#定义通式)
		+ [代码](#代码-2)
	- [miller-rabin素性判断](#miller-rabin素性判断)
		+ [代码](#代码-3)
* [STL](#stl)
	- [求合并,交集,并集，差集](#求合并交集并集差集)
	- [二分查找](#二分查找)
	- [字符串操作](#字符串操作)
	- [读入优化](#读入优化)
* [Java](#java)
	- [a+b problem](#ab-problem)
	- [BigInteger](#biginteger)
		+ [构造函数](#构造函数)
		+ [方法](#方法)

<!-- vim-markdown-toc -->

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

### splay

#### 头文件&宏&全局变量&结构体

```c++
#define MAXN 600010//MAXN是插入操作数量

struct node{
    int num;
    node* p;
    node* son[2];
    int lazy;
    bool lazyr;
    int size;
    int maxnum;
}tree[MAXN],*nil,*root;
int top=0;
//用于伪内存分配
```

#### 辅助函数

```c++
node* new_node()//伪内存分配
{
    return &tree[top++];
}

void push_down(node* nown)//下移懒惰标记
{
    if(nown->lazyr)
    {
        swap(nown->son[0],nown->son[1]);
        for(int i=0;i<2;i++)nown->son[i]->lazyr=!nown->son[i]->lazyr;
        nown->lazyr=false;
    }
    for(int i=0;i<2;i++)nown->son[i]->lazy+=nown->lazy;
    nown->num+=nown->lazy;
    nown->maxnum+=nown->lazy;
    nown->lazy=0;
}

void push_up(node* nown)//上移计算
{
    node *left=nown->son[0],*right=nown->son[1];
    nown->maxnum=nown->num;
    if(left!=nil)
    {
        push_down(left);
        nown->maxnum=max(left->maxnum,nown->maxnum);
    }
    if(right!=nil)
    {
        push_down(right);
        nown->maxnum=max(right->maxnum,nown->maxnum);
    }
    nown->size=left->size+right->size+1;
}

void push_up_parents(node* nown)
{
    while(nown!=root)
    {
        nown=nown->p;
        push_up(nown);
    }
}

void plant(node *nown,node *p,int i)
{
    nown->p=p;
    p->son[i]=nown;
}

void rotate(node *nown)//旋转操作
{
    node *p=nown->p;
    if(p==root)root=nown;
    else plant(nown,p->p,p==p->p->son[0]?0:1);
    int i=(nown==p->son[0])?0:1;
    plant(nown->son[i^1],p,i);
    plant(p,nown,i^1);
    push_up(p);
    push_up(nown);
}

void splay(node *nown,node* &r)//splay操作，把nown伸展到根r
{
    while(nown!=r)
    {
        if(nown->p==r)rotate(nown);
        else
        {
            int i=(nown->p==nown->p->p->son[0])?0:1;
            int j=(nown==nown->p->son[0])?0:1;
            if(i^j)rotate(nown);
            else rotate(nown->p);
            rotate(nown);
        }
    }
}

void insert(node *nown,int k)//把树nown插入到位置k
{
    if(root==nil)
    {
        root=nown;
        return;
    }
    node* p=root;
    while(1)
    {
        push_down(p);
        int i=(k<=p->son[0]->size)?0:1;
        if(p->son[i]==nil)
        {
            plant(nown,p,i);
            break;
        }
        if(i)k-=p->son[0]->size+1;
        p=p->son[i];
    }
    push_up_parents(nown);
}

node *node_find(int k)//返回第k个数的节点(不要单独使用)
{
    node *nown=root;
    while(nown!=nil)
    {
        push_down(nown);
        if(nown->son[0]->size==k)return nown;
        else
        {
            int i=(k<nown->son[0]->size)?0:1;
            if(i)k-=nown->son[0]->size+1;
            nown=nown->son[i];
        }
    }
    return nil;
}

node *interval_find(int l,int r)//返回一棵splay树，包含区间[l,r]的节点
{
    if(l==0&&r==root->size-1)return root;
    else if(l==0)
    {
        splay(node_find(r+1),root);
        return root->son[0];
    }
    else if(r==root->size-1)
    {
        splay(node_find(l-1),root);
        return root->son[1];
    }
    splay(node_find(l-1),root);
    splay(node_find(r+1),root->son[1]);
    return root->son[1]->son[0];
}
```

#### 初始化

```c++
void init_tree()//初始化
{
    top=0;
    nil=new_node();
    nil->size=0;
    root=nil;
}
```

#### 可使用函数

```c++
void insert(int x,int k)//把x插入到位置k
{
    node *nown=new_node();
    nown->num=x;
    nown->size=1;
    nown->lazy=0;
    nown->lazyr=false;
    nown->maxnum=x;
    nown->son[0]=nown->son[1]=nil;
    insert(nown,k);
    splay(nown,root);
}

int tree_find(int k)//返回第k个数的值
{
    node *nown=node_find(k);
    splay(nown,root);
    return nown->num;
}

node* erase(int l,int r)//删除区间[l,r],并返回被删除的树的根
{
    node *nown=interval_find(l,r);
    if(nown==root)root=nil;
    else
    {
        if(nown==nown->p->son[0])nown->p->son[0]=nil;
        else nown->p->son[1]=nil;
        push_up_parents(nown);
    }
    return nown;
}

void flip(int l,int r)//把区间[l,r]的数翻转
{
    node *nown=interval_find(l,r);
    nown->lazyr=!nown->lazyr;
}

void add_num(int l,int r,int x)//把区间[l,r]的数都加x
{
    node *nown=interval_find(l,r);
    nown->lazy+=x;
    push_up_parents(nown);
}

void cut(int l,int r,int k)//把区间[l,r]的数裁剪下来，并放到剩下树的第k个位置
{
    node *nown=erase(l,r);
    insert(nown,k);
}

int max_num(int l,int r)
{
    node *nown=interval_find(l,r);
    push_down(nown);
    return nown->maxnum;
}
```

#### 使用方法

```c++
insert(x,k);//把x插入到第k个位置，建树只用循环插入就行了O(n)!
erase(l,r);//删除区间[l,r],如果想删第k个数，只用erase(k,k);
tree_find(k);//返回第k个数的值
flip(l,r);//翻转区间[l,r]
add_num(l,r,x);//把区间[l,r]都加上x
/*若想要把区间都赋成某个值或者都乘上一个数，添加对应懒惰标记,修改push_down,insert函数并创建对应的修改函数*/
max_num(l,r);//查找区间[l,r]最大值
/*若想维护区间和或最小值等，添加对应成员变量,修改push_down,push_up,insert函数并创建对应的修改函数*/
```

### 矩阵快速幂
#### 代码
```cpp
#define ll long long
const ll MOD = 1000000007;
struct Matrix{
    ll a[N][N];
    int r, c;
}ori, res;

void init(){
    memset(res.a, 0, sizeof(res.a));
    res.r = 1; res.c = 2;
    res.a[1][1] = p;
    res.a[1][2] = 2;
    ori.r = 2; ori.c = 2;//构造矩阵
    ori.a[1][1] = p;
    ori.a[1][2] = 1;
    ori.a[2][1] = -q;
    ori.a[2][2] = 0;
}

Matrix multi(Matrix x, Matrix y)//矩阵乘法
{
    Matrix z;
    memset(z.a, 0, sizeof(z.a));
    z.r = x.r, z.c = y.c;
    for(int i = 1; i <= x.r; i++){
        for(int k = 1; k <= x.c; k++)//加速优化
        {
            if(x.a[i][k] == 0) continue;
            for(int j = 1; j<= y.c; j++)
                z.a[i][j] = (z.a[i][j] + (x.a[i][k] * y.a[k][j]) % MOD) % MOD;
        }
    }
    return z;
}

void Matrix_pow(int n)//矩阵快速幂
{
    while(n){
        if(n & 1)
            res = multi(res, ori);
        ori = multi(ori, ori);
        n >>= 1;
    }
    printf("%llu\n", res.a[1][1] % MOD);
}
```

<!--TODO:-->

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

>求a对b的逆元，即(a^(-1))mod b</br>
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
$$ 有解的判定条件，并用构造法给出了在有解情况下解的具体形式。</br>
中国剩余定理说明：假设整数$m_1,m_2, \cdots ,m_n$两两互质，则对任意的整数：$a1,a2, \cdots ,an$，方程组 有解，并且通解可以用如下方式构造得到：</br>
设</br>
$$ M = m_1 \times m_2 \times m_3 \times \cdots \times m_n = \prod_{i=1}^n m_i $$ 是整数$m_1,m_2, \cdots ,m_n$的乘积，并设</br>
$$ M_i = M \div m_i \ , \forall i \in \left \{ 1, 2, \cdots, n \right \} $$ 是除了mi以外的n- 1个整数的乘积。</br>
设$t_i=M_i^{-1}$为$M_i$模$m_i$的数论倒数($t_i$为$M_i$意义下的逆元) </br>
$$ M_it_i \equiv 1 \left ( mod \ m_i \right ), \forall i \in \left \{ 1,2,\cdots,n \right \} $$ 方程组$\left ( S \right )$的通解形式为</br>
$$
\begin{aligned}
x &= a_1t_1M_1 + a_2t_2M_2 + \cdots + a_nt_nM_n + kM \\
&= kM + \sum_{i=1}^na_it_iM_i, \ k \in \mathbb{Z}
\end{aligned}
$$ 在模$M_i$的意义下，方程组$\left ( S \right )$只有一个解:</br>
$$ x \equiv \left ( a_1t_1M_1 + a_2t_2M_2 + \cdots + a_nt_nM_n \right ) \ mod \ M $$

#### 代码

```c++
#include <iostream>
#include <cstdio>
#include <cstring>

using namespace std;

void extend_Euclid(int a, int b, int &x, int &y)
{
    if(b==0){
        x = 1;
        y = 0;
        return;
    }
    extend_Euclid(b, a%b, x, y);
    int tmp = x;
    x = y;
    y = tmp-(a/b)*y;
}

int CRT(int *a, int *m, int n)
{
    int M = 1, ans = 0;
    for(int i = 1; i <= n; i++) M*=m[i];
    for(int i = 1; i <= n; i++){
        int x, y;
        int Mi = M/m[i];
        extend_Euclid(Mi, m[i], x, y);
        ans = (ans+Mi*x*a[i])%M;
    }
    if(ans < 0) ans+=M;
    return ans;
}

int main()
{
    int a[5], m[5] = {0, 23, 28, 33}, d, kase = 0;
    while(cin>>a[1]>>a[2]>>a[3]>>d)
    {
        if(a[1]==-1&&a[2]==-1&&a[3]==-1&&d==-1)break;
        int ans = CRT(a, m, 3);
        if(ans <= d) ans+=21252;
        printf("Case %d: the next triple peak occurs in %d days.\n", ++kase, ans-d);
    }
    return 0;
}
```

### miller-rabin素性判断

#### 代码

```cpp
const int TIMES = 10;//随机次数

//返回[0, n]之间的一个随机数
ll random0(ll n){
    return ((double)rand() / RAND_MAX*n + 0.5);
}

//快速乘a*b%mod
ll quick_mul(ll a, ll b, ll mod){
    ll ans = 0;
    while(b){
        if(b&1){
            b--;
            ans = (ans+a)%mod;
        }
        b >>= 1;
        a = (a+a) % mod;
    }
    return ans;
}

//快速幂a^b%mod
ll quick_pow(ll a, ll n, ll mod){
    ll ans = 1;
    while(n){
        if(n&1)ans = quick_mul(ans, a, mod);
        a = quick_mul(a, a, mod);
        n >>= 1;
    }
    return ans;
}

bool witness(ll a, ll n){
    ll tmp = n-1;
    int i = 0;
    while(tmp % 2 == 0){
        tmp >>= 1;
        i++;
    }
    ll x = quick_pow(a, tmp, n);
    if(x == 1 || x == n-1)return true;
    while(i--){
        x = quick_mul(x, x, n);
        if(x == n-1)return true;
    }
    return false;
}

bool miller_rabin(ll n){
    if(n == 2)return true;
    if(n < 2 || n % 2 == 0)return false;
    for(int i  = 1; i <= TIMES; i++){
        ll a = random0(n-2)+1;
        if(!witness(a, n))
              return false;
    }
    return true;
}

```

## STL

### 求合并,交集,并集，差集

```cpp
template<class _InIt1,class _InIt2,class _OutIt>
inline_OutIt set_intersection(       //参数格式
  _InIt1 _First1, _InIt1 _Last1,
   _InIt2 _First2, _InIt2 _Last2,
  _OutIt _Dest)
//传进去的两个容器必须是有序的

merge()        //合并
set_intersection()        //交集        A∩B
set_union()                 //并集         A∪B
set_difference()           //差集          A-B
set_symmetric_difference() //并集减去交集  (A-B)∪(B-A)=A∪B - A∩B

用法:
merge(a.begin(),a.end(),b.begin(),b.end(),inserter(c,c.begin()));
```

### 二分查找

```cpp
lower_bound()     //第一个大于等于
upper_bound()    //第一个大于
用法:
lower_bound(a.begin(),a.end(),x); //返回一个迭代器
lower_bound(a,a+n,x) //返回找到元素的指针
```

### 字符串操作

```cpp
strstr(a,b)//在a中找b
```

### 读入优化

```cpp
#include <cctype>

template<class TN>
inline void kread(TN &x)
{
    x=0;
    char c;
    while(!isdigit(c=getchar()));
    do{
        x=x*10+c-'0';
    }while(isdigit(c=getchar()));
}

template<class TN,class... ARGS>
inline void kread(TN &first,ARGS& ... args)
{
    kread(first);
    kread(args...);
}
```

## Java

### a+b problem

```java
import java.util.Scanner;
public class Main{
    public static void main(String args[]){
        Scanner cin = new Scanner(System.in);
        int a, b;
        while (cin.hasNext()){
            a = cin.nextInt(); b = cin.nextInt();
            System.out.println(a + b);
        }
    }
}
```

### BigInteger

#### 构造函数

```java
BigInteger(String val, int radix)
Translates the String representation of a BigInteger in the specified radix into a BigInteger.
```
#### 方法


