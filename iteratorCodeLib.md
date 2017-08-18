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
	- [可持久化线段树](#可持久化线段树)
		+ [全局变量](#全局变量)
		+ [离散化](#离散化)
		+ [线段树相关](#线段树相关)
		+ [实例](#实例)
		+ [用法](#用法)
* [算法](#算法)
	- [矩阵快速幂](#矩阵快速幂)
		+ [代码](#代码)
	- [生成树计数](#生成树计数)
		+ [定理](#定理)
		+ [代码](#代码-1)
	- [次小生成树](#次小生成树)
		+ [全局变量&结构体](#全局变量结构体)
		+ [算法](#算法-1)
	- [最小树形图](#最小树形图)
		+ [宏&常量&结构体&变量](#宏常量结构体变量)
		+ [算法](#算法-2)
* [数论](#数论)
	- [扩展欧几里得](#扩展欧几里得)
		+ [定义](#定义)
		+ [代码](#代码-2)
		+ [求逆元](#求逆元)
	- [中国剩余定理](#中国剩余定理)
		+ [定义&通式](#定义通式)
		+ [代码](#代码-3)
	- [miller-rabin素性判断](#miller-rabin素性判断)
		+ [代码](#代码-4)
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
	- [BigDecimal](#bigdecimal)
		+ [舍入方式](#舍入方式)
		+ [方法](#方法-1)

<!-- vim-markdown-toc -->

## 实用数据结构

### 加权并查集

>解决集合问题中，集合内元素有关系并且关系具有传递性的问题。
>从集合中删除节点的方法：消除该点对集合的影响(如集合中的点个数、和、最值)，然后给它分配一个新的编号(原来的编号不管)

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
### 可持久化线段树

#### 全局变量

```c++
const int maxn = 1e5 + 10;
const int M = maxn * 30;

int n, q, m, tot;
// tot:节点总个数
int a[maxn], t[maxn];
// a:原数组元素 t:将原数组元素按大小去重映射到t数组上
int T[maxn], lson[M], rson[M], c[M];
// T:树的根节点 lson:节点左孩子指针 rson:节点右孩子指针 c:树节点上的权值
```

#### 离散化

```c++
//将原数组a元素按大小去重映射到t数组上
void init_hash() {
	for (int i = 1; i <= n; i++)
	t[i] = a[i];
	sort(t + 1, t + n + 1);
	m = unique(t + 1, t + n + 1) - t - 1;
}

//在t数组上二分查找x，返回位置
int hash(int x) {
	return lower_bound(t + 1, t + 1 + m, x) - t;
}
```

#### 线段树相关

```c++
//建空树
int build(int l, int r) {
	int rt = tot++;
	c[rt] = 0;
	if (l != r) {
		int mid = (l + r) >> 1;
		lson[rt] = build(l, mid);
		rson[rt] = build(mid + 1, r);
	}
	return rt;
}

//更新节点，建立新树
int update(int rt, int pos, int val) {
	int newrt = tot++, tmp = newrt;
	c[newrt] = c[rt] + val;
	int l = 1, r = m;
	while (l < r) {
		int mid = (l + r) >> 1;
		if (pos <= mid) {
			//若更新节点在左子树，则新建左子树的根节点，
			//右子树的根节点利用已有的节点，同时新树与原树同时移向左子树
			lson[newrt] = tot++;
			rson[newrt] = rson[rt];
			newrt = lson[newrt];
			rt = lson[rt];
			r = mid;
		} else { //同理
			rson[newrt] = tot++;
			lson[newrt] = lson[rt];
			newrt = rson[newrt];
			rt = rson[rt];
			l = mid + 1;
		}
		c[newrt] = c[rt] + val; //更新新建节点
	}
	return tmp; //返回新树的根节点
}

//询问[lrt,rrt]区间第k大
int query(int lrt, int rrt, int k) {
	int l = 1, r = m;
	//二分查找t[i]，使得[lrt,rrt]中小于等于t[i]的数的个数为k个
	//则t[i]为区间第k大
	while (l < r) {
		int mid = (l + r) >> 1;
		if (c[lson[lrt]] - c[lson[rrt]] >= k) {
			r = mid;
			lrt = lson[lrt];
			rrt = lson[rrt];
		} else {
			l = mid + 1;
			k -= c[lson[lrt]] - c[lson[rrt]];
			lrt = rson[lrt];
			rrt = rson[rrt];
		}
	}
	return l;
}
```

#### 实例

```c++
//poj 2104
int main(){
	while (~scanf("%d%d", &n, &q)) {

		for (int i = 1; i <= n; i++)
			scanf("%d", a + i);
		init_hash();
		tot = 0;
		T[n + 1] = build(1, m); //建空树

		//向空树中插入原数组中的元素
		for (int i = n; i; i--) {
			int pos = hash(a[i]);
			T[i] = update(T[i + 1], pos, 1);
		}

		//处理询问
		while (q--) {
			int l, r, k;
			scanf("%d%d%d", &l, &r, &k);
			printf("%d\n", t[query(T[l], T[r + 1], k)]);
		}
	}
	return 0;
}
```

#### 用法

```c++
int build(int l, int r)//在[l,r]上建立空树；返回空树的根
int update(int rt, int pos, int val)//建立新树更新以rt为根节点的树上，pos节点，权值+val；返回新树的根
int query(int lrt, int rrt, int k)//返回区间[lrt,rrt]上的第k大
```

## 算法

### 矩阵快速幂

#### 代码

```c++
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

### 生成树计数

#### 定理

>算法引入：
>给定一个无向图G，求它生成树的个数t(G);
>
>算法思想：
>(1)G的度数矩阵D[G]是一个n*n的矩阵,并且满足:当i≠j时,dij=0;当i=j时,dij等于vi的度数;
>(2)G的邻接矩阵A[G]是一个n*n的矩阵,并且满足:如果vi,vj之间有边直接相连,则aij=1,否则为0;
>定义图G的Kirchhoff矩阵C[G]为C[G]=D[G]-A[G];
>Matrix-Tree定理:G的所有不同的生成树的个数等于其Kirchhoff矩阵C[G]任何一个n-1阶主子式的行列式的绝对值；
>所谓n-1阶主子式,就是对于r(1≤r≤n),将C[G]的第r行,第r列同时去掉后得到的新矩阵,用Cr[G]表示;
>
>Kirchhoff矩阵的特殊性质：
>(1)对于任何一个图G,它的Kirchhoff矩阵C的行列式总是0,这是因为C每行每列所有元素的和均为0;
>(2)如果G是不连通的,则它的Kirchhoff矩阵C的任一个主子式的行列式均为0;
>(3)如果G是一颗树,那么它的Kirchhoff矩阵C的任一个n-1阶主子式的行列式均为1;
>
>算法举例： sd:
>SPOJ104(Highways)
>
>题目地址：
><http://www.spoj.com/problems/HIGH/>
>
>题目大意：
>一个有n座城市的组成国家,城市1至n编号,其中一些城市之间可以修建高速公路;
>需要有选择的修建一些高速公路,从而组成一个交通网络;
>计算有多少种方案,使得任意两座城市之间恰好只有一条路径;

#### 代码

```c++
const int N=15;

typedef long long LL;

int degree[N];
LL C[N][N];

LL det(LL a[][N],int n)//生成树计数:Matrix-Tree定理
{
	LL ret=1;
	for(int i=1; i<n; i++)
	{
		for(int j=i+1; j<n; j++)
			while(a[j][i])
			{
				LL t=a[i][i]/a[j][i];
				for(int k=i; k<n; k++)
					a[i][k]=(a[i][k]-a[j][k]*t);
				for(int k=i; k<n; k++)
					swap(a[i][k],a[j][k]);
				ret=-ret;
			}
		if(a[i][i]==0)
			return 0;
		ret=ret*a[i][i];
	}
	if(ret<0)
		ret=-ret;
	return ret;
}

int main()
{
	int tcase;
	scanf("%d",&tcase);
	while(tcase--)
	{
		memset(degree,0,sizeof(degree));
		memset(C,0,sizeof(C));
		int n,m;
		scanf("%d%d",&n,&m);
		int u,v;
		while(m--)
		{
			scanf("%d%d",&u,&v);
			u--;
			v--;
			C[u][v]=C[v][u]=-1;
			degree[u]++;
			degree[v]++;
		}
		for(int i=0; i<n; ++i)
			C[i][i]=degree[i];
		printf("%lld\n",det(C,n));
	}
	return 0;
}
```

### 次小生成树

#### 全局变量&结构体

```c++
const int maxn = 1003;
const double DIS_INF = 999999;

int t, n, x[maxn], y[maxn], p[maxn], cas, book[maxn] = {0}, St[maxn], topSt, used[maxn][maxn] = {0};
//cas:样例数   book:标记点是否在生成树内     St、topSt:储存已经在生成树内的点    used:标记生成树内部的边
double dis[maxn][maxn], tot, maxDis[maxn][maxn], ans, low[maxn];
//dis:权值    tot:生成树总权值  maxDis[i][j]:生成树上i-j路径上最大权

struct Edge {
	int f, t;
	Edge(int _f = 0, int _t = 0) : f(_f), t(_t) {}
	bool operator<(const Edge &Right) const {
		return dis[f][t] > dis[Right.f][Right.t];
	}
};

priority_queue<Edge> pq;
```

#### 算法

```c++
void Prim() {
	book[0] = cas;
	St[topSt++] = 0;
	for (int i = 0; i < n; i++) {
		low[i] = dis[0][i];
		pq.push(Edge(0, i));
	}
	while (!pq.empty()) {
		Edge Front = pq.top();
		pq.pop();
		int f = Front.f, t = Front.t;
		if (book[t] == cas) continue;
		book[t] = cas;
		for (int i = 0; i < topSt; i++) {
			int u = St[i];
			maxDis[u][t] = maxDis[t][u] = max(dis[f][t], maxDis[u][f]); //dp求每一条路径上的最大边
		}
		St[topSt++] = t;
		tot += dis[f][t];
		used[f][t] = used[t][f] = cas;
		for (int i = 0; i < n; i++) {
			if (book[i] != cas && dis[t][i] < low[i]) {
				low[i] = dis[t][i];
				pq.push(Edge(t, i));
			}
		}
	}
}

int main() {
	scanf("%d", &t);
	for (cas = 1; cas <= t; cas++) {
		scanf("%d", &n);
		topSt = 0;
		tot = 0;
		ans = -1;
		for (int i = 0; i < n; i++) {
			scanf("%d%d%d", x + i, y + i, p + i);
			for (int j = 0; j <= i; j++) {
				dis[i][j] = dis[j][i] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
			}
		}
		Prim();
		for(int i=0; i<n; i++){
			for(int j=0; j<i; j++){
				double B;
				if(used[i][j] == cas){
					B = 1.0*(p[i]+p[j])/(tot - dis[i][j]);
				}else{
					B = 1.0*(p[i]+p[j])/(tot - maxDis[i][j]);
				}
				ans = max(B, ans);
			}
		}
		printf("%.2f\n", ans);
	}
	return 0;
}
```

### 最小树形图

#### 宏&常量&结构体&变量

```c++
#define M 109
#define type int

const double eps = 1e-10;
const type inf = (1) << 30;

struct point {
	double x, y;
} p[M];

struct Node {
	int u, v;
	type cost;
} E[M * M + 5];

int pre[M], ID[M], vis[M];
type In[M];
int n, m;
```

#### 算法

```c++
type Directed_MST(int root, int NV, int NE) {
	type ret = 0;
	while (true) {
		//1.找最小入边
		for (int i = 0; i < NV; i++)
			In[i] = inf;
		for (int i = 0; i < NE; i++) {
			int u = E[i].u;
			int v = E[i].v;
			if (E[i].cost < In[v] && u != v) {
				pre[v] = u;
				In[v] = E[i].cost;
			}
		}
		for (int i = 0; i < NV; i++) {
			if (i == root) continue;
			if (In[i] == inf) return -1; //除了跟以外有点没有入边,则根无法到达它
		}
		//2.找环
		int cntnode = 0;
		memset(ID, -1, sizeof(ID));
		memset(vis, -1, sizeof(vis));
		In[root] = 0;
		for (int i = 0; i < NV; i++) { //标记每个环
			ret += In[i];
			int v = i;
			while (vis[v] != i && ID[v] == -1 && v != root) {
				vis[v] = i;
				v = pre[v];
			}
			if (v != root && ID[v] == -1) {
				for (int u = pre[v]; u != v; u = pre[u]) {
					ID[u] = cntnode;
				}
				ID[v] = cntnode++;
			}
		}
		if (cntnode == 0) break; //无环
		for (int i = 0; i < NV; i++)
			if (ID[i] == -1) {
				ID[i] = cntnode++;
			}
		//3.缩点,重新标记
		for (int i = 0; i < NE; i++) {
			int v = E[i].v;
			E[i].u = ID[E[i].u];
			E[i].v = ID[E[i].v];
			if (E[i].u != E[i].v) {
				E[i].cost -= In[v];
			}
		}
		NV = cntnode;
		root = ID[root];
	}
	return ret;
}

int main() {
	while (scanf("%d%d", &n, &m), n + m) {
		for (int i = 0; i < m; i++) {
			scanf("%d%d%d", &E[i].u, &E[i].v, &E[i].cost);
			E[i].u--;
			E[i].v--;
		}
		type ans = Directed_MST(0, n, m);
		if (ans == -1)
			printf("impossible\n");
		else
			printf("%d\n", ans);
	}
	return 0;
}
```

## 数论

### 扩展欧几里得

#### 定义

>对于不完全为0的非负整数ab,gcd(a,b)表示a,b的最大公约数,必然存在整数对x,y,使得gcd(a,b)=ax+by。

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

>求a对b的逆元，即(a^(-1))mod b
>int x,y;
>exgcd(a,b,x,y);
>x即为a对b的逆元

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
$$ 
有解的判定条件，并用构造法给出了在有解情况下解的具体形式。</br>
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

### 欧拉函数

#### 定义&通式

>欧拉函数是小于等于 $n$ 的正整数中与 $n$ 互质的数的数目（$\varphi \left ( 1 \right )=1$）。</br>
通式：$\varphi \left ( x \right ) = x\left ( 1 - \frac{1}{p_1} \right )\left ( 1 - \frac{1}{p_2} \right )\left ( 1 - \frac{1}{p_3} \right )\cdots\left ( 1 - \frac{1}{p_n} \right )$ </br>
应用：欧拉降幂公式</br>$$
a^b \equiv a^{b \  \% \  \varphi \left( n\right) + \varphi \left( n \right)} (mod\ n)\ (b > \varphi (n))
$$

#### 代码

```c++
/*线性筛O(n)时间复杂度内筛出maxn内欧拉函数值*/
int m[maxn],phi[maxn],p[maxn],pt;//m[i]是i的最小素因数，p是素数，pt是素数个数
int make()
{
    phi[1]=1;
    int N=maxn;
    int k;
    for(int i=2;i<N;i++)
    {
        if(!m[i])//i是素数
            p[pt++]=m[i]=i,phi[i]=i-1;
        for(int j=0;j<pt&&(k=p[j]*i)<N;j++)
        {
            m[k]=p[j];
            if(m[i]==p[j])//为了保证以后的数不被再筛，要break
            {
                phi[k]=phi[i]*p[j];
/*这里的phi[k]与phi[i]后面的∏(p[i]-1)/p[i]都一样（m[i]==p[j]）只差一个p[j]，就可以保证∏(p[i]-1)/p[i]前面也一样了*/
                break;    
            }
            else
                phi[k]=phi[i]*(p[j]-1);//积性函数性质，f(i*k)=f(i)*f(k)
        }
    }
}


/*直接求解欧拉函数*/
int phi(int n)
{
    int ans=1;
    for(int i=2;i*i<=n;i++)
    {
        if(n%i==0)
        {
            ans*=i-1;
            n/=i;
            while(n%i==0)
            {
                ans*=i;
                n/=i;
            }
        }
    }
    if(n>1)ans*=(n-1);
    return ans;
}
```

### 素数筛法

#### 线形筛

```c++
int top;
int prime[MAXN];
bool notPrime[MAXN];
 
void calcPrime()
{
    for(int i=2;i<MAXN;i++)
    {
        if(!notPrime[i])prime[top++]=i;
        for(int j=0;j<top&&i*prime[j]<MAXN;j++)
        {
            notPrime[i*prime[j]]=true;
            if(i%prime[j]==0)break;
        }
    }
}
```

#### 复杂度 $O(n^{\frac{3}{4}})$

```c++
#include <bits/stdc++.h>  
#define ll long long

using namespace std;  

ll f[340000],g[340000],n; 

void init(){  
    ll i,j,m;  
    for(m=1;m*m<=n;++m)f[m]=n/m-1;  
    for(i=1;i<=m;++i)g[i]=i-1;  
    for(i=2;i<=m;++i){  
        if(g[i]==g[i-1])continue;  
        for(j=1;j<=min(m-1,n/i/i);++j){  
            if(i*j<m)f[j]-=f[i*j]-g[i-1];  
            else f[j]-=g[n/i/j]-g[i-1];  
        }  
        for(j=m;j>=i*i;--j)g[j]-=g[j/i]-g[i-1];  
    }  
}  

int main(){  
    while(scanf("%I64d",&n)!=EOF){  
        init();  
        cout<<f[1]<<endl;  
    }  
    return 0;  
}
```
#### 复杂度 $O(n^{\frac{2}{3}})$

```c++
#include<cstdio>  
#include<cmath>  
using namespace std;  
#define LL long long  

const int N = 5e6 + 2;  
bool np[N];  
int prime[N], pi[N];  

int getprime()  
{  
    int cnt = 0;  
    np[0] = np[1] = true;  
    pi[0] = pi[1] = 0;  
    for(int i = 2; i < N; ++i)  
    {  
        if(!np[i]) prime[++cnt] = i;  
        pi[i] = cnt;  
        for(int j = 1; j <= cnt && i * prime[j] < N; ++j)  
        {  
            np[i * prime[j]] = true;  
            if(i % prime[j] == 0)   break;  
        }  
    }  
    return cnt;  
}  

const int M = 7;  
const int PM = 2 * 3 * 5 * 7 * 11 * 13 * 17;  
int phi[PM + 1][M + 1], sz[M + 1]; 

void init()  
{  
    getprime();  
    sz[0] = 1;  
    for(int i = 0; i <= PM; ++i)  phi[i][0] = i;  
    for(int i = 1; i <= M; ++i)  
    {  
        sz[i] = prime[i] * sz[i - 1];  
        for(int j = 1; j <= PM; ++j) phi[j][i] = phi[j][i - 1] - phi[j / prime[i]][i - 1];  
    }  
}

int sqrt2(LL x)  
{  
    LL r = (LL)sqrt(x - 0.1);  
    while(r * r <= x)   ++r;  
    return int(r - 1);  
}

int sqrt3(LL x)  
{  
    LL r = (LL)cbrt(x - 0.1);  
    while(r * r * r <= x)   ++r;  
    return int(r - 1);  
}

LL getphi(LL x, int s)  
{  
    if(s == 0)  return x;  
    if(s <= M)  return phi[x % sz[s]][s] + (x / sz[s]) * phi[sz[s]][s];  
    if(x <= prime[s]*prime[s])   return pi[x] - s + 1;  
    if(x <= prime[s]*prime[s]*prime[s] && x < N)  
    {  
        int s2x = pi[sqrt2(x)];  
        LL ans = pi[x] - (s2x + s - 2) * (s2x - s + 1) / 2;  
        for(int i = s + 1; i <= s2x; ++i) ans += pi[x / prime[i]];  
        return ans;  
    }  
    return getphi(x, s - 1) - getphi(x / prime[s], s - 1);  
}

LL getpi(LL x)  
{  
    if(x < N)   return pi[x];  
    LL ans = getphi(x, pi[sqrt3(x)]) + pi[sqrt3(x)] - 1;  
    for(int i = pi[sqrt3(x)] + 1, ed = pi[sqrt2(x)]; i <= ed; ++i) ans -= getpi(x / prime[i]) - i + 1;  
    return ans;  
}

LL lehmer_pi(LL x)  
{  
    if(x < N)   return pi[x];  
    int a = (int)lehmer_pi(sqrt2(sqrt2(x)));  
    int b = (int)lehmer_pi(sqrt2(x));  
    int c = (int)lehmer_pi(sqrt3(x));  
    LL sum = getphi(x, a) +(LL)(b + a - 2) * (b - a + 1) / 2;  
    for (int i = a + 1; i <= b; i++)  
    {  
        LL w = x / prime[i];  
        sum -= lehmer_pi(w);  
        if (i > c) continue;  
        LL lim = lehmer_pi(sqrt2(w));  
        for (int j = i; j <= lim; j++) sum -= lehmer_pi(w / prime[j]) - (j - 1);  
    }  
    return sum;  
}

int main()  
{  
    init();  
    LL n;  
    while(~scanf("%lld",&n))  
    {  
        printf("%lld\n",lehmer_pi(n));  
    }  
    return 0;  
}
```

### miller-rabin素性判断

#### 代码

```c++
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
### 莫比乌斯函数

#### 定义

> $$ \mu = \begin{cases} 1 & n=1 \\ (-1)^k & n = p_1p_2\cdots p_k \\ 0 & other \end{cases}$$

#### 莫比乌斯反演

> $$f(n) = \sum_{d,n}g(d)=\sum_{d,n} g(\frac{n}{d})$$ 
> $$ g(n) = \sum_{d,n} \mu(d) f(\frac{n}{d}) = \sum_{d,n} \mu(\frac{n}{d})f(d) $$
>倍数形式只用把$\frac{n}{d}$变为$\frac{d}{n}$

#### 技巧
>若$g(d)=[\frac n d]*[\frac m d]$之类的阶梯状函数</br>
记录$\mu$的前缀和

```c++
nt d=1;
int ans=0;
while(d<=min(n,m))
{
    int last=min(n/(n/d),m/(m/d));
    ans+=(sum[last]-sum[d-1])*(n/d)*(m/d);
    d=last+1;
}
```

## STL

### 求合并,交集,并集，差集

```c++
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

```c++
lower_bound()     //第一个大于等于
upper_bound()    //第一个大于
用法:
lower_bound(a.begin(),a.end(),x); //返回一个迭代器
lower_bound(a,a+n,x) //返回找到元素的指针
```

### 字符串操作

```c++
strstr(a,b)//在a中找b
```

### 读入优化

```c++
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

| 返回值            | 函数                                      | 简介                                                                                       |
|:------------------|:------------------------------------------|:-------------------------------------------------------------------------------------------|
| BigInteger        | abs()                                     | Returns a BigInteger whose value is the absolute value of this BigInteger.                 |
| BigInteger        | add(BigInteger val)                       | Returns a BigInteger whose value is (this + val).                                          |
| BigInteger        | and(BigInteger val)                       | Returns a BigInteger whose value is (this & val).                                          |
| BigInteger        | andNot(BigInteger val)                    | Returns a BigInteger whose value is (this & ~val).                                         |
| int               | compareTo(BigInteger val)                 | Compares this BigInteger with the specified BigInteger.                                    |
| BigInteger        | divide(BigInteger val)                    | Returns a BigInteger whose value is (this / val).                                          |
| BigInteger[]      | divideAndRemainder(BigInteger val)        | Returns an array of two BigIntegers containing (this / val) followed by (this % val).      |
| double            | doubleValue()                             | Converts this BigInteger to a double.                                                      |
| boolean           | equals(Object x)                          | Compares this BigInteger with the specified Object for equality.                           |
| BigInteger        | gcd(BigInteger val)                       | Returns a BigInteger whose value is the greatest common divisor of abs(this) and abs(val). |
| BigInteger        | max(BigInteger val)                       | Returns the maximum of this BigInteger and val.                                            |
| BigInteger        | min(BigInteger val)                       | Returns the minimum of this BigInteger and val.                                            |
| BigInteger        | mod(BigInteger m)                         | Returns a BigInteger whose value is (this mod m).                                          |
| BigInteger        | modInverse(BigInteger m)                  | Returns a BigInteger whose value is (this ^ -1 mod m).                                     |
| BigInteger        | modPow(BigInteger exponent, BigInteger m) | Returns a BigInteger whose value is (this ^ exponent mod m).                               |
| BigInteger        | multiply(BigInteger val)                  | Returns a BigInteger whose value is (this * val).                                          |
| BigInteger        | negate()                                  | Returns a BigInteger whose value is (-this).                                               |
| BigInteger        | or(BigInteger val)                        | Returns a BigInteger whose value is (this &#124; val).                                     |
| BigInteger        | pow(int exponent)                         | Returns a BigInteger whose value is (this ^ exponent).                                     |
| BigInteger        | remainder(BigInteger val)                 | Returns a BigInteger whose value is (this % val).                                          |
| BigInteger        | shiftLeft(int n)                          | Returns a BigInteger whose value is (this << n).                                           |
| BigInteger        | shiftRight(int n)                         | Returns a BigInteger whose value is (this >> n).                                           |
| BigInteger        | subtract(BigInteger val)                  | Returns a BigInteger whose value is (this - val).                                          |
| String            | toString()                                | Returns the decimal String representation of this BigInteger.                              |
| String            | toString(int radix)                       | Returns the String representation of this BigInteger in the given radix.                   |
| static BigInteger | valueOf(long val)                         | Returns a BigInteger whose value is equal to that of the specified long.                   |
| BigInteger        | xor(BigInteger val)                       | Returns a BigInteger whose value is (this ^ val).                                          |

### BigDecimal

#### 舍入方式

>以下在roundingMode参数填入
>ROUND_CEILING向正无穷方向舍入 
>ROUND_DOWN向零方向舍入
>ROUND_FLOOR向负无穷方向舍入 
>ROUND_HALF_DOWN 
>向（距离）最近的一边舍入，除非两边（的距离）是相等,如果是这样，向下舍入, 例如1.55 保留一位小数结果为1.5 
>
>ROUND_HALF_EVEN 
>向（距离）最近的一边舍入，除非两边（的距离）是相等,如果是这样，如果保留位数是奇数，使用ROUND_HALF_UP ，如果是偶数，使用ROUND_HALF_DOWN 
>
>ROUND_HALF_UP 
>向（距离）最近的一边舍入，除非两边（的距离）是相等,如果是这样，向上舍入, 1.55保留一位小数结果为1.6 
>
>ROUND_UNNECESSARY 计算结果是精确的，不需要舍入模式

#### 方法

| 返回值     | 函数                                                    |
|:-----------|:--------------------------------------------------------|
| BigDecimal | divide(BigDecimal divisor, int roundingMode)            |
| BigDecimal | divide(BigDecimal divisor, int scale, int roundingMode) |
| BigDecimal | setScale(int newScale)                                  |
| BigDecimal | setScale(int newScale, int roundingMode)                |



#### 方法