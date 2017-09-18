# 实用数据结构

## 加权并查集

解决集合问题中，集合内元素有关系并且关系具有传递性的问题  
从集合中删除节点的方法：消除该点对集合的影响(如集合中的点个数、和、最值)，然后给它分配一个新的编号(原来的编号不管)

### 头文件&宏&全局变量

```c++
#define MAXN 100000//最大点数
int p[MAXN];//父节点
int v[MAXN];//到父节点边的权值(加权解决集合中点的相互关系的问题)
```

### 初始化

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

### 查找

执行后p[a]为a所在集合的根节点,v[a]为a到其集合的根节点的权值

```c++
void union_find(int a)
{
    if(p[a]==a)return;
    union_find(p[a]);
    v[a]+=v[p[a]];
    p[a]=p[p[a]];
}
```

### 合并

合并a,b所在集合,vab为 a到b的权值

```c++
void union_merge(int a,int b,int vab)
{
    union_find(a);
    union_find(b);
    v[p[a]]=vab-v[a]+v[b];
    p[p[a]]=p[b];
}
```

## 树状数组

要求所有数的和不能超出范围,也可修改为记录最值  
数组下标应从1开始

### 头文件&宏&全局变量

```c++
#define MAXN 10000//数组大小 
int tree[MAXN];//树状数组 
```

### 辅助函数

```c++
int lowbit(int a)
{
    return a&-a;
}
```

### 初始化

```c++
memset(tree,0,sizeof(tree));
```

### 单点修改

把a处的值增加b（如果是修改，需要记录原始数组，转化为增加就行了）

```c++
void tree_add(int a,int b)
{
    for(int i=a;i<MAXN;i+=lowbit(i))tree[i]+=b;
}
```

### 区间查询

查询1到a之间的值的和

```c++
int tree_find(int a)
{
    int ans=0;
    for(int i=a;i;i-=lowbit(i))ans+=tree[i];
    return ans;
}
```

## 线段树

### 宏&结构体&全局变量

```c++
#define MAXN 10000

struct Tree{
    int v;//此区间存的值 
    int lazy_inc;//整个区间被增加的值 
    bool lazy;//区间是否被整体修改过 
    int lazy_chg;//区间被整体修改后的值(lazy==true时有效) 
}tree[MAXN*4];

int kkke[MAXN];//用于初始化的数组
```

### 辅助函数

```c++
int lson(int k){return k<<1;}
int rson(int k){return (k<<1)|1;}
void tree_update(int k)//更新此区间存的值 
{
    tree[k].v=max(tree[lson(k)].v,tree[rson(k)].v)+tree[k].lazy_inc;
    /**可更改(max/min/sum)**/
}
void tree_chg(int left,int right,int k,int v)//将编号k的区间全变为v 
{
    tree[k].v=v;/**可更改(求和则为v*(right-left+1))**/
    tree[k].lazy=true;
    tree[k].lazy_chg=v;
    tree[k].lazy_inc=0;
}
void tree_pushdown(int left,int right,int k)//将整体修改的信息向下传递 
{
    if(tree[k].lazy)
    {
        tree[k].lazy=false;
        int mid=(left+right)>>1;
        tree_chg(left,mid,lson(k),tree[k].lazy_chg);
        tree_chg(mid+1,right,rson(k),tree[k].lazy_chg);
    }
}
```

### 初始化

```c++
void tree_init()
{
    memset(tree,0,sizeof(tree));//清0 
}
void tree_build(int left,int right,int k)//初始化维护某个数组 
{
    tree[k].lazy=false;
    tree[k].lazy_inc=0;
    if(left==right)
    {
        tree[k].v=kkke[left];/**需对应为原数组的名称**/
    }
    else
    {
        int mid=(left+right)>>1;
        tree_build(left,mid,lson(k));
        tree_build(mid+1,right,rson(k));
        tree_update(k);
    }
}
```

### 主要使用函数

```c++
//区间增加
//把l到r间的值都增加v 
void tree_add(int left,int right,int k,int l,int r,int v)
{
    if(l<=left&&right<=r)
    {
        tree[k].v+=v;
        tree[k].lazy_inc+=v;
    }
    else
    {
        tree_pushdown(left,right,k);
        int mid=(left+right)>>1;
        if(l<=mid)tree_add(left,mid,lson(k),l,r,v);
        if(r>mid)tree_add(mid+1,right,rson(k),l,r,v);
        tree_update(k);
    }
}
//区间修改 
//把l到r间的值都修改为v 
void tree_change(int left,int right,int k,int l,int r,int v)
{
    if(l<=left&&right<=r)tree_chg(left,right,k,v);
    else
    {
        tree_pushdown(left,right,k);
        int mid=(left+right)>>1;
        if(l<=mid)tree_change(left,mid,lson(k),l,r,v);
        if(r>mid)tree_change(mid+1,right,rson(k),l,r,v);
        tree_update(k);
    }
}
//区间查询
//查询区间[l,r]维护的值
int tree_find(int left,int right,int k,int l,int r,int v)
{
    if(l<=left&&right<=r)return tree[k].v;
    else
    {
        tree_pushdown(left,right,k);
        int mid=(left+right)>>1;
        if(l<=mid&&r>mid)
            return max(tree_find(left,mid,lson(k),l,r,v)
                      ,tree_find(mid+1,right,rson(k),l,r,v));
        /**可更改(max/min/sum)**/
        if(l<=mid)return tree_find(left,mid,lson(k),l,r,v);
        return tree_find(mid+1,right,rson(k),l,r,v);
    }
}
```

## 树链剖分

### 宏&全局变量&结构体

```c++
int n;//点数

int size[MAXN];//子树大小
int dep[MAXN];//节点深度
int pa[MAXN];//直系父节点
int PA[MAXN];//重链开始处
int id[MAXN];//编号
int ID;

struct Edge{
    int to;
    int v;
    int next;
}edge[MAXN*2];//大小至少为点的二倍
int head[MAXN],top;
int kkke[MAXN];//储存点权
```

### 初始化&加边函数

```c++
void init()
{
    kkke[0]=0;
    ID=0;
    top=0;
    memset(head,-1,sizeof(head));
}

void addEdge(int a,int b,int v)
{
    edge[top].to=a;
    edge[top].v=v;
    edge[top].next=head[b];
    head[b]=top++;

    edge[top].to=b;
    edge[top].v=v;
    edge[top].next=head[a];
    head[a]=top++;
}
```

### 主要函数

```c++
void calcSize(int nown=1,int p=-1,int DEP=0)//计算每棵子树大小
{
    pa[nown]=p;
    size[nown]=1;
    dep[nown]=DEP;
    for(int i=head[nown];i!=-1;i=edge[i].next)
    {
        if(edge[i].to==p)continue;
        calcSize(edge[i].to,nown,DEP+1);
        size[nown]+=size[edge[i].to];
    }
}

//树链剖分，id[PA[nown]]~id[nown]间都属于这条链
void dfs(int nown=1,int p=1)
{
    id[nown]=ID++;
    PA[nown]=p;
    int maxi=-1;
    for(int i=head[nown];i!=-1;i=edge[i].next)
    {
        if(edge[i].to==pa[nown])continue;
        if(maxi==-1||size[edge[i].to]>size[maxi])
            maxi=edge[i].to;
    }
    if(maxi==-1)return;
    dfs(maxi,p);
    for(int i=head[nown];i!=-1;i=edge[i].next)
    {
        if(edge[i].to==pa[nown]||edge[i].to==maxi)continue;
        dfs(edge[i].to,edge[i].to);
    }
    for(int i=head[nown];i!=-1;i=edge[i].next)
    {
        if(edge[i].to==pa[nown])continue;
        kkke[id[edge[i].to]]=edge[i].v;//初始化点权
    }
}
```

### 参考查找函数

```c++
//线段树函数省略
//a->b路径长度
int findDist(int a,int b)
{
    int ans=0;
    while(PA[a]!=PA[b])
    {
        if(dep[PA[a]]<dep[PA[b]])swap(a,b);
        ans+=treeFind(0,n-1,1,id[PA[a]],id[a]);
        a=pa[PA[a]];
    }
    if(dep[a]<dep[b])ans+=treeFind(0,n-1,1,id[a]+1,id[b]);
    else if(dep[a]>dep[b])ans+=treeFind(0,n-1,1,id[b]+1,id[a]);
    return ans;
}
```

### 使用方法

```c++
init()
for(each a->b)addEdge(a,b,v);
calcSize();
dfs();
```

## splay

### 头文件&宏&全局变量&结构体

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

### 辅助函数

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
    push_down(nown);
    splay(nown,root);
}

node *kth_node(int k)//返回第k个数的节点,并旋转至根
{
    node *nown=root;
    while(nown!=nil)
    {
        push_down(nown);
        if(nown->son[0]->size==k)
        {
            splay(nown,root);
            return nown;
        }
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
        splay(kth_node(r+1),root);
        return root->son[0];
    }
    else if(r==root->size-1)
    {
        splay(kth_node(l-1),root);
        return root->son[1];
    }
    splay(kth_node(l-1),root);
    splay(kth_node(r+1),root->son[1]);
    return root->son[1]->son[0];
}
```

### 初始化

```c++
void init_tree()//初始化
{
    top=0;
    nil=new_node();
    nil->size=0;
    root=nil;
}
```

### 可使用函数

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
}

int kth(int k)//返回第k个数的值
{
    return kth_node(k)->num;
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

### 使用方法

```c++
insert(x,k);//把x插入到第k个位置，建树只用循环插入就行了O(n)!
erase(l,r);//删除区间[l,r],如果想删第k个数，只用erase(k,k);
kth(k);//返回第k个数的值
flip(l,r);//翻转区间[l,r]
add_num(l,r,x);//把区间[l,r]都加上x
/*若想要把区间都赋成某个值或者都乘上一个数，添加对应懒惰标记,
修改push_down,insert函数并创建对应的修改函数*/
max_num(l,r);//查找区间[l,r]最大值
/*若想维护区间和或最小值等，添加对应成员变量,修改push_down,
push_up,insert函数并创建对应的修改函数*/
```
## 可持久化线段树

### 全局变量

```c++
#define MAXM 6666666//插入次数*log(表示范围)

using namespace std;

struct Tree{
    int num;
    int lson;
    int rson;
}tree[MAXM];//线段树
int top;
```

### 主要代码

```c++
void treeInit()
{
    tree[0].num=tree[0].lson=tree[0].rson=0;//用0节点表示NULL,便于处理
    top=1;
}

int treeAdd(int ori,int left,int right,int x,int a)
{//在ori[left,right]树上x位置加a,并返回新的根
    int nown=top++;
    tree[nown]=tree[ori];
    tree[nown].num+=a;
    if(left<right)
    {
        int mid=(left+right)>>1;
        if(x<=mid)tree[nown].lson=treeAdd(tree[nown].lson,left,mid,x,a);
        else tree[nown].rson=treeAdd(tree[nown].rson,mid+1,right,x,a);
    }
    return nown;
}

int treeFind(int nown,int left,int right,int l,int r)//查询区间[l,r]
{
    if(nown==0)return 0;
    if(l<=left&&right<=r)return tree[nown].num;
    int mid=(left+right)>>1;
    int ans=0;
    if(l<=mid)ans+=treeFind(tree[nown].lson,left,mid,l,r);
    if(r>mid)ans+=treeFind(tree[nown].rson,mid+1,right,l,r);
    return ans;
}
```

### 用法

```c++
//使用前先treeInit()初始化
int treeAdd(int ori,int left,int right,int x,int a)//x位置加a并返回新的根
int treeFind(int nown,int left,int right,int l,int r)//查询区间[l,r]
//需要保证传入x,l,r∈[left,right],并且每次left与right应相同
```

## 舞蹈链

行列下标皆从1开始 

### 头文件&宏&全局变量

```c++
const int MN = 1005; //最大行数
const int MM = 1005; //最大列数
const int MNN = MN*MM; //最大点数
```

### 精准覆盖

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

### 重复覆盖

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

# 字符串

## 最小循环表示

### 代码

```c++
//最小循环表示
//input: str[2*len] 原串扩展了一倍的串。如原串为“abc”， str为“abcabc”。
//         len 原串的长度。
//output：ptr 最小循环表示法的起始下标。
int min_representation(int len){
    int i = 0, j = 1, k = 0;
    while(i < len && j < len && k < len){
        if(str[i+k] == str[j+k])k++;
        else{
            if(str[i+k] < str[j+k]) j += k+1;
            else i += k+1;
            k = 0;
            if(i == j)j++;
        }
    }
    return min(i, j);
}
```

## 最长回文Manacher

### 代码

```c++
const int MAXN = 111234;

char orign[MAXN], str[2 * MAXN]; //字符串
int radius[2 * MAXN];            //对称轴为i的最长回文半径

//orign:初始字符串、str:插入间隔符的字符串（长度为orign的两倍加一）
//raidus:对称轴为i的最长回文半径、mark:间隔符
int Manacher(char *orign, char *str, int *radius, char mark) {
    //------------插入间隔符------------
    int len = strlen(orign);
    for (int i = 0; i < len; i++) {
        str[2 * i + 1] = orign[i];
        str[2 * i] = mark;
    }
    str[2 * len] = mark;
    str[2 * len + 1] = '\0';
    len = 2 * len + 1;
    //------------插入间隔符------------

    int ans = 2; //答案至少为2，即至少为 #a# 形式
    int max_right = 0, pos = 0;
    // max_right表示当前已知所有回文串右端点最大值，pos表示该回文的对称轴
    for (int i = 0; i < len; i++) {
        if (i < max_right)
            radius[i] = min(radius[2 * pos - i], max_right - i);
        else
            radius[i] = 1;
        //判断边界、对应字符是否相等
        while(i-radius[i] >= 0 && i+radius[i] < len 
                && str[i-radius[i]] == str[i+radius[i]]){
            radius[i]++;
        }
        ans = max(ans, radius[i]);
        //更新max_right、pos
        if (i + radius[i] - 1 > max_right) {
            max_right = i + radius[i] - 1;
            pos = i;
        }
    }
    return ans - 1;
}
```

## KMP

### 宏&全局变量

```c++
#define MAXN 6666666

int nextn[MAXN];
```

### 核心代码

```c++
void initNext(const char *pattern)
{
    nextn[0]=-1;
    int i=0,j=-1;
    while(pattern[i])
    {
        while(j!=-1&&pattern[i]!=pattern[j])j=nextn[j];
        i++;
        j++;
        nextn[i]=j;
    }
}

int kmp(const char *s,const char *pattern,bool flag=true)
{
    if(flag)initNext(pattern);
    int i=0,j=0,cnt=0;
    while(s[i])
    {
        while(j!=-1&&s[i]!=pattern[j])j=nextn[j];
        i++;
        j++;
        if(!pattern[j])cnt++;
    } 
    return cnt;
}
```

## AC自动机

### 头文件&宏&全局变量

```c++
#include <queue>

#define MAXN 666666
#define MAXK 26//字符数量

struct Node{
    Node *son[MAXK];
    Node *fail;
    int num;//以此节点为末尾的模式串数量
    bool flag;//去重用,可选
}node[MAXN],*root,*top;
queue<Node*>q;//建立自动机时使用
```

### 辅助函数

```c++
int mapToK(char c)//把字符隐射到0~MAXK-1
{
    return c-'a';
}

Node *newNode()
{
    memset(top->son,0,sizeof(top->son));
    top->num=0;
    return top++;
}

void initNode()//初始化节点分配
{
    top=node;
    root=newNode();
}
```

### 主要函数

```c++
void addPattern(char *s)//添加模式串
{
    Node *nown=root;
    while(*s)
    {
        int k=mapToK(*s);
        if(!nown->son[k])nown->son[k]=newNode();
        nown=nown->son[k];
        s++;
    }
    nown->num++;
}

void buildACAutoMaton()//计算fail
{
    root->fail=nullptr;
    q.push(root);
    while(!q.empty())
    {
        Node *nown=q.front();q.pop();
        for(int i=0;i<MAXK;i++)
        {
            if(nown->son[i])
            {
                Node *p=nown->fail;
                while(p&&!p->son[i])p=p->fail;
                nown->son[i]->fail=p?p->son[i]:root;
                q.push(nown->son[i]);
            }
        }
    }
}
```

### 可选参考函数

```c++
void initFlag()//初始化去重标记
{
    root->flag=false;
    for(Node *i=root+1;i<top;i++)
        i->flag=true;
}

int match(char *s)//返回匹配次数
{
    initFlag();
    Node *nown=root;
    int ans=0;
    while(*s)
    {
        int k=mapToK(*s);
        while(nown&&!nown->son[k])nown=nown->fail;
        nown=nown?nown->son[k]:root;
        for(Node *i=nown;i->flag;i=i->fail)
        {
            ans+=i->num;
            i->flag=false;
        }
        s++;
    }
    return ans;
}
```

### 用法

```c++
//先initNode();初始化
//然后addPattern添加模式串
//最后buildACAutoMaton
```

## 后缀数组

### 宏&全局变量

```c++
#define MAXN 666666//大于字符串长度二倍

int krank[MAXN];//第i个元素是第几大 1~n
int SA[MAXN];//第i大的元素在原数组中位置 1~n
int height[MAXN];
int *const tmp=height;//一开始height没用,使用它当tmp
int cnt[MAXN];//用于基数排序,统计
int st[MAXN][30];//st表
int LOG[MAXN];//log表
```

### 辅助函数

```c++
void initHeight(char *s,int n)//计算height数组
{
    int j,k=0;
    for(int i=1;i<=n;height[krank[i++]]=k)
        for(k=max(k-1,0),j=SA[krank[i]-1];krank[i]>1&&s[i+k-1]==s[j+k-1];k++)
            ;
}

void initLOG()
{
    if(LOG[0]==-1)return;
    LOG[0]=-1;
    for(int i=1;i<MAXN;i++)
        LOG[i]=(i&(i-1))?LOG[i-1]:LOG[i-1]+1;
}

void initSt(int n)
{
    initLOG();
    for(int i=0;i<n;i++)st[i][0]=height[i+1];
    for(int j=1;(1<<j)<=n;j++)
        for(int i=0;i+(1<<j)<=n;i++)
            st[i][j]=min(st[i][j-1],st[i+(1<<(j-1))][j-1]);
}

bool comp(int n,int a,int b,int w)
{
    //判断a和b,a+w和b+w的第一关键字是否对应相等
    if(tmp[a]==tmp[b])
    {
        if(a+w>n||b+w>n)
        {
            if(a+w>n&&b+w>n)return true;
            return false;
        }
        if(tmp[a+w]==tmp[b+w])return true;
    }
    return false;
}

bool rSort(int n,int &m,int w)
{
    //krank当作第一关键字，tmp相当于第二关键字的SA
    //此时第二关键字已有序,顺序是tmp
    memset(cnt+1,0,m*sizeof(cnt[0]));
    for(int i=1;i<=n;i++)cnt[krank[i]]++;//统计
    for(int i=2;i<=m;i++)cnt[i]+=cnt[i-1];
    for(int i=n;i;i--)//比其第一关键字小的数量就是其新位置
        SA[cnt[krank[tmp[i]]]--]=tmp[i];

    //用tmp的空间暂存rank
    memcpy(tmp+1,krank+1,n*sizeof(krank[0]));
    krank[SA[1]]=m=1;
    for(int i=2;i<=n;i++)//生成新的rank
        krank[SA[i]]=comp(n,SA[i],SA[i-1],w)?m:++m;
    return m>=n;//分为n类,排序完成

}
```

### 主要函数

```c++
void initSA(char *s,int n)//初始化后缀数组
{
    int m=0;
    for(int i=1;i<=n;i++)
    {
        krank[i]=s[i-1];
        m=max(m,krank[i]);
        tmp[i]=i;
    }
    int w=0;
    while(!rSort(n,m,w))
    {
        if(w)w<<=1;
        else w=1;
        //重新计算tmp
        int top=0;
        for(int i=n-w+1;i<=n;i++)tmp[++top]=i;//越界的最小
        for(int i=1;i<=n;i++)
            if(SA[i]>w)//不越界的从小到大排
                tmp[++top]=SA[i]-w;
    }
    initHeight(s,n);
    initSt(n);
}

int calcLCP(int l,int r)//后缀l到后缀r的最长公共前缀
{
    l=krank[l];r=krank[r];
    if(l>r)swap(l,r);
    int k=LOG[r-l];
    return min(st[l][k],st[r-(1<<k)][k]);
}
```

### 用法

```c++
//调用initSA后height,SA,krank数组都计算好了
//下标都从1开始
//调用calcLCP计算LCP,不需要可以去掉LOG表和st表
```

# 图论

## 最短路

### dijkstra配对堆优化

复杂度$\Theta \left ( m \right )$

#### 头文件&宏&全局变量

```c++
#include <ext/pb_ds/priority_queue.hpp>

#define MAXN 666
#define MAXM 6666666
#define INF 0x3f3f3f3f

struct Edge{
    int to;
    int v;
    int next;
}edge[MAXM];
int head[MAXN];
int top;

typedef __gnu_pbds::priority_queue<pair<int,int>,greater<pair<int,int>>
,__gnu_pbds::pairing_heap_tag>
Heap;

Heap heap;
Heap::point_iterator pit[MAXN];
```

#### 初始化&加边

```c++
void initEdge()
{
    memset(head,-1,sizeof(head));
    top=0;
}

void addEdge(int a,int b,int v)
{
    edge[top].to=b;
    edge[top].v=v;
    edge[top].next=head[a];
    head[a]=top++;
}
```

#### 核心代码

```c++
void dijkstra(int n,int S,int dist[])//点标号从0开始
{
    for(int i=0;i<n;i++)
        dist[i]=INF;
    dist[S]=0;
    for(int i=0;i<n;i++)
        pit[i]=heap.push(make_pair(dist[i],i));
    while(!heap.empty())
    {
        int nown=heap.top().second;heap.pop();
        for(int i=head[nown];i!=-1;i=edge[i].next)
            if(dist[edge[i].to]>dist[nown]+edge[i].v)
                heap.modify(pit[edge[i].to],
                                make_pair(dist[edge[i].to]=dist[nown]+edge[i].v,
                                          edge[i].to));
    }
}
```

## 差分约束

* 第一：  
    * 感觉难点在于建图  

* 第二：  
    * ①：对于差分不等式，a - b <= c ，建一条 b 到 a 的权值为 c 的边，求的是最短路，得到的是最大值  
    * ②：对于不等式 a - b >= c ，建一条 b 到 a 的权值为 c 的边，求的是最长路，得到的是最小值  
    * ③：存在负环的话是无解  
    * ④：求不出最短路（dist[ ]没有得到更新）的话是任意解  

* 第三：  

    * 一种建图方法：  
    设x[i]是第i位置（或时刻）的值（跟所求值的属性一样），那么把x[i]看成数列，前n项和为s[n]，则x[i] = s[i] - s[i-1]；  
    那么这样就可以最起码建立起类似这样的一个关系：0 <= s[i] - s[i-1] <= 1;  

## 最大权匹配Kuhn-Munkres

复杂度$O(n^{2}m)$

### 头文件&宏&结构体&全局变量

```c++
#define MAXN 666
#define MAXM 666666
#define INF 0x3f3f3f3f

using namespace std;

struct Edge{
    int to;
    int v;
    int next;
}edge[MAXM];
int head[MAXN],top;

int from[MAXN];
int X[MAXN];
int Y[MAXN];
bool visX[MAXN];
bool visY[MAXN];
```

### 初始化&加边
```c++
void initEdge()
{
    memset(head,-1,sizeof(head));
    top=0;
}

void addEdge(int a,int b,int v)//第一个集合的a连向第二个集合的b
{
    edge[top].to=b;
    edge[top].v=v;
    edge[top].next=head[a];
    head[a]=top++;
}
```

### 辅助函数

```c++
bool dfs(int nown)//匈牙利找增广路
{
    visX[nown]=true;
    for(int i=head[nown];i!=-1;i=edge[i].next)
    {
        if(visY[edge[i].to])continue;
        if(X[nown]+Y[edge[i].to]!=edge[i].v)continue;
        visY[edge[i].to]=true;
        if(from[edge[i].to]==-1||dfs(from[edge[i].to]))
        {
            from[edge[i].to]=nown;
            return true;
        }
    }
    return false;
}
```

### 核心代码

```c++
int KM(int n)//n为点数,需保证有完备匹配,标号从0开始
{
    int ans=0;
    for(int i=0;i<n;i++)
    {
        from[i]=-1;
        Y[i]=0;
        X[i]=-INF;
        for(int j=head[i];j!=-1;j=edge[j].next)
            X[i]=max(X[i],edge[j].v);
        ans+=X[i];
    }
    for(int k=0;k<n;)
    {
        memset(visX,0,sizeof(visX));
        memset(visY,0,sizeof(visY));
        if(dfs(k))k++;
        else
        {
            int d=INF;
            for(int i=0;i<n;i++)
                if(visX[i])
                    for(int j=head[i];j!=-1;j=edge[j].next)
                        if(!visY[edge[j].to])
                            d=min(d,X[i]+Y[edge[j].to]-edge[j].v);
            ans-=d;
            for(int i=0;i<n;i++)
            {
                if(visX[i])X[i]-=d;
                if(visY[i])Y[i]+=d;
            }
        }
    }
    return ans;
}
```

### 用法

```c++
//需保证是二分图且有完备匹配(两边点数相同)
//两边点的标号都从0开始
initEdge();//初始化
for(each a->b)addEdge(a,b,v);//加边
int ans=KM();//求解
```

## 全局最小割SW

复杂度 $O(nm)$

### 头文件&宏&全局变量

```c++
#include <algorithm>
#include <ext/pb_ds/priority_queue.hpp>

#define MAXN 3333
#define MAXM 444444//最好是边数的两倍

const long long INF=0x3f3f3f3f;

int n,m;

struct Edge{
    long long v;
    int to;
    int next;
    int re;
}edge[MAXM];//边
int head[MAXN],top;//邻接链表
int dist[MAXN];

typedef __gnu_pbds::priority_queue<pair<int,int>,
    less<pair<int,int>>,__gnu_pbds::pairing_heap_tag
    > Heap;
Heap pq;
Heap::point_iterator pqIterator[MAXN];
```

### 建图

```c++
void init()//初始化链表
{
    top=0;
    memset(head,-1,sizeof(head));
}

void addEdge(int a,int b,long long v)//a->b,容量为v的边
{
    if(v==0)return;
    edge[top].v=v;
    edge[top].to=b;
    edge[top].re=top+1;
    edge[top].next=head[a];
    head[a]=top++;

    edge[top].v=v;
    edge[top].to=a;
    edge[top].re=top-1;
    edge[top].next=head[b];
    head[b]=top++;
}
```

### 核心代码

```c++
int findST(int &s,int &t)//找到某一s点和t点间最小割
{
    for(int i=1;i<=n;i++)
    {
        if(head[i]!=-2)
        {
            dist[i]=0;
            pqIterator[i]=pq.push(make_pair(dist[i],i));
        }
    }//初始化
    while(pq.size()>1)
    {
        s=pq.top().second;pq.pop();
        pqIterator[s]=pq.end();
        for(int j=head[s];j!=-1;j=edge[j].next)
            if(pqIterator[edge[j].to]!=pq.end())
                pq.modify(pqIterator[edge[j].to]
                         ,make_pair(dist[edge[j].to]+=edge[j].v
                                    ,edge[j].to));
    }
    t=pq.top().second;pq.pop();
    return dist[t];//dist[t]为s-t最小割
}

void merge(int s,int t)//合并s和t点
{
    int i=head[t],next;
    while(i!=-1)
    {
        next=edge[i].next;
        edge[i].next=head[s];
        head[s]=i;
        edge[edge[i].re].to=s;
        i=next;
    }
    head[t]=-2;//标记t被合并
}

int StoerWagner()
{
    int mincut=INF,s,t;
    for(int i=1;i<n;i++)//最多合并n-1次
    {
        mincut=min(mincut,findST(s,t));
        if(mincut==0)return 0;//达到下限
        merge(s,t);
    }
    return mincut;
}
```

### 用法

```c++
int work()
{
    init();
    int a,b;
    long long v;
    for(int i=0;i<m;i++)
    {
        kread(a,b,v);
        addEdge(a,b,v);
    }
    return StoerWagner();
}

int main()
{
    while(~scanf("%d%d",&n,&m))printf("%d\n",work());
    return 0;
}
```

## 网络流Dinic

### 头文件&全局变量&宏

```c++
#include <algorithm>
#include <cstring>
#include <queue>

using namespace std;

const int MAXN = 6666;
const int MAXM=66666;
const int INF = 0x3f3f3f3f;

int S,T;
int n;

int head[MAXN*2],top;
int cur[MAXN*2];
int level[MAXN*2];
queue<int>q;
struct the_edge{
    int next;
    int to;
    int v;
    int re;
}edge[MAXM];
int va[MAXN];
int vb[MAXN];
int m;
```

### 建图

```c++
void init_edge()
{
    memset(head,-1,sizeof(head));
    top=0;
}

void add_edge(int a,int b,int v)
{
    edge[top].to=b;
    edge[top].v=v;
    edge[top].next=head[a];
    head[a]=top++;

    edge[top].to=a;
    edge[top].v=0;
    edge[top].next=head[b];
    head[b]=top++;
}
```

### 辅助函数

```c++
bool bfs()
{
    memset(level,-1,sizeof(level));
    level[S]=0;
    q.push(S);
    while(!q.empty())
    {
        int nown=q.front();q.pop();
        for(int i=head[nown];i!=-1;i=edge[i].next)
        {
            if(!edge[i].v||level[edge[i].to]!=-1)continue;
            level[edge[i].to]=level[nown]+1;
            q.push(edge[i].to);
        }
    }
    return level[T]!=-1;
}

int dfs(int nown,int maxf)
{
    if(nown==T)return maxf;
    int nowf=0,flow;
    for(int &i=cur[nown];i!=-1;i=edge[i].next)
    {
        if(!edge[i].v||
            level[edge[i].to]!=level[nown]+1)continue;
        if((flow=dfs(edge[i].to,min(maxf-nowf,edge[i].v)))!=0)
        {
            nowf+=flow;
            edge[i].v-=flow;
            edge[i^1].v+=flow;
            if(nowf==maxf)return maxf;
        }
    }
    return nowf;
}
```

### 核心代码

``` c++
int dinic()
{
    int ans=0;
    while(bfs())
    {
        memcpy(cur,head,sizeof(cur));
        ans+=dfs(S,INF);
    }
    return ans;
}
```

### 用法

```c++
int main(){
    kread(n,m);
    S=1;T=n;
    init_edge();
    for(int i=2;i<n;i++)kread(va[i]);
    for(int i=2;i<n;i++)kread(vb[i]);
    int x,y,a,b,c;
    for(int i=0;i<m;i++)
    {
        kread(x,y,a,b,c);
        add_edge(x,y,c);
        add_edge(y,x,c);
    }
    printf("%d\n",-dinic());
    return 0;
}
```

## 最小费用流

### 头文件&宏&全局变量

```c++
#include <queue>

#define MAXN 2222
#define MAXM MAXN*MAXN
#define INF 0x3f3f3f3f 

using namespace std;

int S,T; //源点 汇点   
struct Edge  
{  
    int from,to,flow,worth,next; //结点，流量，费用，链表   
    Edge(){}  
    Edge(int fr,int ro,int fl,int wo,int ne)  
    {  
        from=fr,to=ro,flow=fl,worth=wo,next=ne;  
    }  
}edge[MAXM];  
int head[MAXN]; // 建立链表  
int top;  //边数
bool visque[MAXN]; //查看是否入队  
int dis[MAXN]; //最小距离  
int pre[MAXN],prx[MAXN]; //记录路线用于更新残量图   
queue<int>q;  
```

### 建图

```c++
void init() //初始化  
{  
    memset(head,0,sizeof(head));  
    top=2;//必须是2
}

void addEdge(int from,int to,int flow,int worth)  //建图   
{  
    edge[top]=Edge(from,to,flow,worth,head[from]);  
    head[from]=top++;  
    edge[top]=Edge(to,from,0,-worth,head[to]);    //反向弧   
    head[to]=top++;  
}  
```

### 辅助函数

```c++
int bfs() //寻找最短路  
{  

    while(!q.empty()) q.pop(); //初始化队列  
    for(int i=0;i<=MAXN;i++) dis[i]=INF; //初始化距离   
    q.push(S); //源点入队  
    dis[S]=0;  
    visque[S]=true;  
    while(!q.empty())  
    {  
        int u=q.front();  
        q.pop();  
        for(int i=head[u];i;i=edge[i].next)  
        {  
            if(edge[i].flow>0
                &&dis[u]+edge[i].worth<dis[edge[i].to]) 
                //更新最短路   
            {  
                dis[edge[i].to]=dis[u]+edge[i].worth;  
                pre[edge[i].to]=u;  
                prx[edge[i].to]=i;  
                if(!visque[edge[i].to])  
                {  
                    visque[edge[i].to]=true;  
                    q.push(edge[i].to);  
                }  
            }  
        }  
        visque[u]=false; //前面已经让u出队了所以这里要写一下   
    }   
    return dis[T]!=INF; //判断是否可以到达汇点   
}   
int dfs()  
{  
    int u=T;  
    int ans=INF;  
    while(u!=S) //找当前路中的最小流量   
    {  
        if(edge[prx[u]].flow<ans) ans=edge[prx[u]].flow;  
        u=pre[u];  
    }  
    u=T;  
    while(u!=S) //更新残量图   
    {  
        edge[prx[u]].flow-=ans;  
        edge[prx[u]^1].flow+=ans;  
        u=pre[u];  
    }  
    return ans*dis[T];  
}
```

### 主要函数

```c++
int solve()  
{  
    int ans=0;  
    while(bfs())  
    {  
        ans+=dfs();  
    }   
    return ans;  
}
```

## 生成树计数

### 定理

算法引入：  
给定一个无向图G，求它生成树的个数t(G);
 
算法思想：  
(1)G的度数矩阵D[G]是一个n*n的矩阵,并且满足:当i≠j时,dij=0;当i=j时,dij等于vi的度数;  
(2)G的邻接矩阵A[G]是一个n*n的矩阵,并且满足:如果vi,vj之间有边直接相连,则aij=1,否则为0;  
定义图G的Kirchhoff矩阵C[G]为C[G]=D[G]-A[G];  
Matrix-Tree定理:G的所有不同的生成树的个数等于其Kirchhoff矩阵C[G]任何一个n-1阶主子式的行列式的绝对值；  
所谓n-1阶主子式,就是对于r(1≤r≤n),将C[G]的第r行,第r列同时去掉后得到的新矩阵,用Cr[G]表示;
 
Kirchhoff矩阵的特殊性质：  
(1)对于任何一个图G,它的Kirchhoff矩阵C的行列式总是0,这是因为C每行每列所有元素的和均为0;  
(2)如果G是不连通的,则它的Kirchhoff矩阵C的任一个主子式的行列式均为0;  
(3)如果G是一颗树,那么它的Kirchhoff矩阵C的任一个n-1阶主子式的行列式均为1;
 
算法举例： sd:  
SPOJ104(Highways)  
 
题目地址：  
<http://www.spoj.com/problems/HIGH/>
 
题目大意：  
一个有n座城市的组成国家,城市1至n编号,其中一些城市之间可以修建高速公路;  
需要有选择的修建一些高速公路,从而组成一个交通网络;  
计算有多少种方案,使得任意两座城市之间恰好只有一条路径;

### 代码

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

## 次小生成树

### 全局变量&结构体

```c++
const int maxn = 1003;
const double DIS_INF = 999999;

int t, n, x[maxn], y[maxn], p[maxn], cas, book[maxn] = {0}
//cas:样例数 book:标记点是否在生成树内   
int St[maxn], topSt, used[maxn][maxn] = {0};
//St、topSt:储存已经在生成树内的点 used:标记生成树内部的边
double dis[maxn][maxn], tot, maxDis[maxn][maxn], ans, low[maxn];
//dis:权值 tot:生成树总权值 maxDis[i][j]:生成树上i-j路径上最大权

struct Edge {
    int f, t;
    Edge(int _f = 0, int _t = 0) : f(_f), t(_t) {}
    bool operator<(const Edge &Right) const {
        return dis[f][t] > dis[Right.f][Right.t];
    }
};

priority_queue<Edge> pq;
```

### 算法

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
            maxDis[u][t] = maxDis[t][u]
                         = max(dis[f][t], maxDis[u][f]); 
                         //dp求每一条路径上的最大边
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
                dis[i][j] = dis[j][i] 
                          = sqrt((x[i] - x[j]) *(x[i] - x[j])
                              + (y[i] - y[j])* (y[i] - y[j]));
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

## 最小树形图

### 宏&常量&结构体&变量

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

### 算法

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
            if (In[i] == inf) return -1;
            //除了跟以外有点没有入边,则根无法到达它
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

## 强连通分量

### 变量&宏

```c++
#include <cstring>
#include <vector>
const int N = 10010;
vector<int> G[N];//邻接表存图
vector<int> rG[N];//存反向图
vector<int> vs;//后序遍历顺序的顶点列表
bool vis[N];
int cmp[N];//所属强连通分量的拓扑序
```

### 主要函数

```c++
void add_edge(int u, int v){
    G[u].push_back(v);
    rG[v].push_back(u);
}

//input: u 顶点
//output: vs 后序遍历顺序的顶点列表
void dfs(int u){
    vis[u] = true;
    for(int i = 0; i < G[u].size(); i++){
        int v = G[u][i];
        if(!vis[v])
              dfs(v);
    }
    vs.push_back(u);
}

//input: u 顶点编号; k 拓扑序号
//output: cmp[] 强连通分量拓扑序
void rdfs(int u, int k){
    vis[u] = true;
    cmp[u] = k;
    for(int i = 0; i < rG[u].size(); i++){
        int v = rG[u][i];
        if(!vis[v])
              rdfs(v, k);
    }
}

//Strongly Connected Component 强连通分量
//input: n 顶点个数
//output: k 强连通分量数;
int scc(int n){
    memset(vis, 0, sizeof(vis));
    vs.clear();
    for(int u = 0; u < n; u++)
        if(!vis[u])
              dfs(u);
    int k = 0;
    memset(vis, 0, sizeof(vis));
    for(int i = vs.size()-1; i >= 0; i--)
          if(!vis[vs[i]])
              rdfs(vs[i], k++);
    return k;
}
```

## 2-SAT问题

【2-SAT问题】  
​    现有一个由N个布尔值组成的序列A，给出一些限制关系，比如A[x] AND A[y]=0、A[x] OR A[y] OR A[z]=1等，要确定A[0..N-1]的值，使得其满足所有限制关系。这个称为SAT问题，特别的，若每种限制关系中最多只对两个元素进行限制，则称为2-SAT问题。  

由于在2-SAT问题中，最多只对两个元素进行限制，所以可能的限制关系共有11种：  
A[x]  
NOT A[x]  
A[x] AND A[y]  
A[x] AND NOT A[y]  
A[x] OR A[y]  
A[x] OR NOT A[y]  
NOT (A[x] AND A[y])  
NOT (A[x] OR A[y])  
A[x] XOR A[y]  
NOT (A[x] XOR A[y])  
A[x] XOR NOT A[y]  
进一步，A[x] AND A[y]相当于(A[x]) AND (A[y])（也就是可以拆分成A[x]与A[y]两个限制关系），NOT(A[x] OR A[y])相当于NOT A[x] AND NOT A[y]（也就是可以拆分成NOT A[x]与NOT A[y]两个限制关系）。因此，可能的限制关系最多只有9种。  

在实际问题中，2-SAT问题在大多数时候表现成以下形式：有N对物品，每对物品中必须选取一个，也只能选取一个，并且它们之间存在某些限制关系（如某两个物品不能都选，某两个物品不能都不选，某两个物品必须且只能选一个，某个物品必选）等，这时，可以将每对物品当成一个布尔值（选取第一个物品相当于0，选取第二个相当于1），如果所有的限制关系最多只对两个物品进行限制，则它们都可以转化成9种基本限制关系，从而转化为2-SAT模型。  

【建模】  
其实2-SAT问题的建模是和实际问题非常相似的。  
建立一个2N阶的有向图，其中的点分为N对，每对点表示布尔序列A的一个元素的0、1取值（以下将代表A[i]的0取值的点称为i，代表A[i]的1取值的点称为i'）。显然每对点必须且只能选取一个。然后，图中的边具有特定含义。若图中存在边<i, j>，则表示若选了i必须选j。可以发现，上面的9种限制关系中，后7种二元限制关系都可以用连边实现，比如NOT(A[x] AND A[y])需要连两条边<x, y'>和<y, x'>，A[x] OR A[y]需要连两条边<x', y>和<y', x>。而前两种一元关系，对于A[x]（即x必选），可以通过连边<x', x>来实现，而NOT A[x]（即x不能选），可以通过连边<x, x'>来实现。  

【O(NM)算法：求字典序最小的解】  
根据2-SAT建成的图中边的定义可以发现，若图中i到j有路径，则若i选，则j也要选；或者说，若j不选，则i也不能选；  
因此得到一个很直观的算法：  
（1）给每个点设置一个状态V，V=0表示未确定，V=1表示确定选取，V=2表示确定不选取。称一个点是已确定的当且仅当其V值非0。设立两个队列Q1和Q2，分别存放本次尝试选取的点的编号和尝试不选的点的编号。  
（2）若图中所有的点均已确定，则找到一组解，结束，否则，将Q1、Q2清空，并任选一个未确定的点i，将i加入队列Q1，将i'加入队列Q2；  
（3）找到i的所有后继。对于后继j，若j未确定，则将j加入队列Q1；若j'（这里的j'是指与j在同一对的另一个点）未确定，则将j'加入队列Q2；  
（4）遍历Q2中的每个点，找到该点的所有前趋（这里需要先建一个补图），若该前趋未确定，则将其加入队列Q2；  
（5）在（3）（4）步操作中，出现以下情况之一，则本次尝试失败，否则本次尝试成功：  
<1>某个已被加入队列Q1的点被加入队列Q2；  
<2>某个已被加入队列Q2的点被加入队列Q1;  
<3>某个j的状态为2；  
<4>某个i'或j'的状态为1或某个i'或j'的前趋的状态为1；  
（6）若本次尝试成功，则将Q1中的所有点的状态改为1，将Q2中所有点的状态改为2，转（2），否则尝试点i'，若仍失败则问题无解。  
该算法的时间复杂度为O(NM)（最坏情况下要尝试所有的点，每次尝试要遍历所有的边），但是在多数情况下，远远达不到这个上界。  
具体实现时，可以用一个数组vst来表示队列Q1和Q2。设立两个标志变量i1和i2（要求对于不同的i，i1和i2均不同，这样可以避免每次尝试都要初始化一次，节省时间），若vst[i]=i1则表示i已被加入Q1，若vst[i]=i2则表示i已被加入Q2。不过Q1和Q2仍然是要设立的，因为遍历（BFS）的时候需要队列，为了防止重复遍历，加入Q1（或Q2）中的点的vst值必然不等于i1（或i2）。中间一旦发生矛盾，立即中止尝试，宣告失败。  

该算法虽然在多数情况下时间复杂度到不了O(NM)，但是综合性能仍然不如下面的O(M)算法。不过，该算法有一个很重要的用处：求字典序最小的解！
如果原图中的同一对点编号都是连续的（01、23、45……）则可以依次尝试第0对、第1对……点，每对点中先尝试编号小的，若失败再尝试编号大的。这样一定能求出字典序最小的解（如果有解的话），因为一个点一旦被确定，则不可更改。  
如果原图中的同一对点编号不连续（比如03、25、14……）则按照该对点中编号小的点的编号递增顺序将每对点排序，然后依次扫描排序后的每对点，先尝试其编号小的点，若成功则将这个点选上，否则尝试编号大的点，若成功则选上，否则（都失败）无解。  

```c++
bool check(int limit){
    //1.建图
    //2.强连通分量缩点
    scc(2n)
    //3.判断
    for(int i = 0; i < n; i++)
        if(cmp[i] == cmp[i+n])
            return false;
    return true;
}
```

## tarjan缩点/找割边/找割点

> 复杂度$\Theta\left(m\right)$

### 环境

```c++
#define MAXN 100001//点数
#define MAXM 200002//边数

//此题是无向图找割边,然后把环缩点
int n,m;
int DFN;
int dfn[MAXN];//dfs序
int low[MAXN];//不经过父节点能访问的最早的点的dfn
int stk[MAXN],tp;//栈
int fa[MAXN];//并查集,fa[i]为缩点后i所在的点

struct Edge{
    int to;
    int next;
}edge[MAXM];//无向图
int head[MAXN],top;
```

### 初始化&加边

```c++
inline void initTarjan(int n)//点编号从1开始
{
    tp=DFN=top=0;
    for(int i=1;i<=n;i++)
    {
        fa[i]=i;
        head[i]=dfn[i]=-1;
    }
}

inline void addEdge(int a,int b)
{
    edge[top].to=b;
    edge[top].next=head[a];
    head[a]=top++;

    edge[top].to=a;
    edge[top].next=head[b];
    head[b]=top++;
}
```

### 辅助函数(并查集,栈)

```c++
inline void stkPush(int x)
{
    stk[tp++]=x;
}

inline int stkPop()
{
    return stk[--tp];
}

inline int unionFind(int x)
{
    while(fa[fa[x]]!=fa[x])fa[x]=fa[fa[x]];
    return fa[x];
}

inline int reAddEdge(int e,int a)//边e交给a
{
    if(unionFind(edge[e].to)==a)return edge[e].next;
    int tmp=edge[e].next;
    edge[e].next=head[a];
    head[a]=e;
    edge[e^1].to=a;
    return tmp;
}
```

### 核心代码

```c++
void tarjan(int nown,int p)//缩点
{
    if(dfn[nown]!=-1)return;
    dfn[nown]=low[nown]=DFN++;
    stkPush(nown);
    for(int i=head[nown];i!=-1;i=edge[i].next)
    {
        if(edge[i].to==p)continue;
        tarjan(edge[i].to,nown);
        low[nown]=min(low[nown],low[edge[i].to]);
    }
    if(dfn[nown]==low[nown])
    {
        //此时p为割点
        //把栈顶到nown的点合并到nown上
        int i;
        for(i=tp-1;stk[i]!=nown;i--)
            fa[stk[i]]=nown;
        //重连边
        i=head[nown];head[nown]=-1;
        while(i!=-1)i=reAddEdge(i,nown);
        while((i=stkPop())!=nown)
            for(int j=head[i];j!=-1;j=reAddEdge(j,nown))
                ;
    }
}
```

### 用法

```c++
//要求无重边,无自环!!!
initTarjan();
for(e:u->v)addEdge(u,v);
tarjan(1,0);
//之后图变为以1为根节点的树
```

## 二分图最大匹配

### 建图&变量

```c++
const int N = 1010;
int head[N], tot;
struct Edge{
    int to, next;
}edge[N*N];

void init(){
    tot = 0;
    memset(head, -1, sizeof(head));
}

void add_edge(int u, int v){
    edge[tot].to = v;
    edge[tot].next = head[u];
    head[u] = tot++;

    edge[tot].to = u;
    edge[tot].next = head[v];
    head[v] = tot++;
}
```

### 匈牙利算法代码

复杂度： $O(nm)$

```c++
int matching[N];
int check[N];

bool dfs(int u){
    for(int i =  head[u]; i != -1; i = edge[i].next){
        int v = edge[i].to;
        if(!check[v]){//要求不在交替路
            check[v] = 1;//放入交替路
            if(matching[v] == -1 || dfs(matching[v])){
                //如果是未匹配点，说明交替路为增广路，则交换路径，并返回成功
                matching[u] = v;
                matching[v] = u;
                return true;
            }
        }
    }
    return false;//不存在增广路
}

//hungarian: 二分图最大匹配匈牙利算法
//input: null
//output: ans 最大匹配数
int hungarian(){
    int ans = 0;
    memset(matching, -1, sizeof(matching));
    for(int u = 0; u < p; u++){
        if(matching[u] == -1){
            memset(check, 0, sizeof(check));
            if(dfs(u))
              ans++;
        }
    }
    return ans;
}
```

### Hopcroft_Karp算法 变量&建图

```c++
const int N = 100000;
const int M = 20000000;
const int INF = 0x3f3f3f3f;
int head[N], tot;
struct Edge{
    int to, next;
}edge[M];

void add_edge(int u, int v){
    edge[tot].to = v;
    edge[tot].next = head[u];
    head[u] = tot++;
}
```

### Hopcroft_Karp算法代码

复杂度：$O(\sqrt{V}E)$

```c++
//xlink[i]表示左集合顶点i匹配的右集合的点，ylink[i]表示右集合顶点i匹配的左集合的点
int xlink[N], ylink[N];
//xlevel[i]表示左集合顶点i的所在层数，ylevel[i]表示右集合顶点i的所在层数
int xlevel[N], ylevel[N];
bool vis[N];
struct Hopcroft_Karp{
    int dis, xn, yn;//xn表示左集合顶点个数，yn表示右集合顶点个数
    void init(int _xn, int _yn){
        tot = 0;
        xn = _xn;
        yn = _yn;
        memset(head, -1, sizeof(head));
        memset(xlink, -1, sizeof(xlink));
        memset(ylink, -1, sizeof(ylink));
    }
    bool bfs(){
        queue<int> que;
        dis = INF;
        memset(xlevel, -1, sizeof(xlevel));
        memset(ylevel, -1, sizeof(ylevel));
        for(int i = 0; i < xn; i++)
            if(xlink[i] == -1){
                que.push(i);
                xlevel[i] = 0;
            }
        while(!que.empty()){
            int u = que.front();
            que.pop();
            if(xlevel[u] > dis)break;
            for(int i = head[u]; i != -1; i = edge[i].next){
                int v = edge[i].to;
                if(ylevel[v] == -1){
                    ylevel[v] = xlevel[u] + 1;
                    if(ylink[v] == -1)
                          dis = ylevel[v];
                    else{
                        xlevel[ylink[v]] = ylevel[v]+1;
                        que.push(ylink[v]);
                    }
                }
            }
        }
        return dis != INF;
    }
    int dfs(int u){
        for(int i = head[u]; i != -1; i = edge[i].next){
            int v = edge[i].to;
            if(!vis[v] && ylevel[v] == xlevel[u]+1){
                vis[v] = 1;
                if(ylink[v] != -1 && ylevel[v] == dis)
                      continue;
                if(ylink[v] == -1 || dfs(ylink[v])){
                    xlink[u] = v;
                    ylink[v] = u;
                    return 1;
                }
            }
        }
        return 0;
    }
    //二分图最大匹配
    //input：建好的二分图
    //output：ans 最大匹配数
    int max_match(){
        int ans = 0;
        while(bfs()){
            memset(vis, 0, sizeof(vis));
            for(int i = 0; i < xn; i++)
                  if(xlink[i] == -1)
                      ans += dfs(i);
        }
        return ans;
    }
}hk_match;

```

### Hopcroft_Karp算法 用法

```c++
hk_match.init(n, m);//n为左集合大小，m为右集合大小
for(e in E)
    add_edge(u, v);
hk_match.max_match();
```

# 数论

## 二分等比数列求和

### 代码

```c++
ll sum(ll a, ll p, ll n) {
    if (n <= 1) return a;
    if (n & 1)
        return (sum(a, p, n / 2) * (1 + quickPow(p, n / 2)) % MOD
                + quickPow(p, n - 1)) % MOD;
    else
        return (sum(a, p, n / 2) * (1 + quickPow(p, n / 2))) % MOD;
}
```

### 用法

调用函数 $sum(a,p,n)$ ,求 $\sum_{i=0}^{n-1}ap^i$

## 扩展欧几里得

### 定义

对于不完全为0的非负整数ab,gcd(a,b)表示a,b的最大公约数,必然存在整数对x,y,使得gcd(a,b)=ax+by。

### 代码

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

### 求逆元

求a对b的逆元，即(a^(-1))mod b  
int x,y;  
exgcd(a,b,x,y);  
x即为a对b的逆元

## 矩阵快速幂

### 代码

```c++
#define MAXN 111
typedef long long ll;

const ll MOD = 1000000007;
struct Matrix{
    ll a[MAXN][MAXN];
    int r, c;
};

Matrix multi(const Matrix &x, const Matrix &y)//矩阵乘法
{
    Matrix z;
    memset(z.a, 0, sizeof(z.a));
    z.r = x.r, z.c = y.c;
    for(int i = 0; i < x.r; i++){
        for(int k = 0; k < x.c; k++)//加速优化
        {
            if(x.a[i][k] == 0) continue;
            for(int j = 0; j< y.c; j++)
                z.a[i][j] =(
                    z.a[i][j] + (x.a[i][k] * y.a[k][j]) % MOD
                    ) % MOD;
        }
    }
    return z;
}

Matrix kpow(Matrix a,Matrix b,int n)//a*b^n
{
    while(n){
        if(n & 1)
            a = multi(a, b);
        b = multi(b, b);
        n >>= 1;
    }
    return a;
}
```

## 中国剩余定理

### 定义&通式

给出了以下的一元线性同余方程组：

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

有解的判定条件，并用构造法给出了在有解情况下解的具体形式。  
中国剩余定理说明：假设整数$m_1,m_2, \cdots ,m_n$两两互质，则对任意的整数：$a1,a2, \cdots ,an$，方程组 有解，并且通解可以用如下方式构造得到：  
设

$$
M = m_1 \times m_2 \times m_3 \times \cdots \times m_n = \prod_{i=1}^n m_i 
$$

是整数$m_1,m_2, \cdots ,m_n$的乘积，并设

$$
M_i = M \div m_i \ , \forall i \in \left \{ 1, 2, \cdots, n \right \} 
$$

是除了$m_i$以外的$n-1$个整数的乘积。  
设$t_i=M_i^{-1}$为$M_i$模$m_i$的数论倒数($t_i$为$M_i$意义下的逆元) 

$$
M_it_i \equiv 1 \left ( mod \ m_i \right ), \forall i \in \left \{ 1,2,\cdots,n \right \}
$$

方程组$\left ( S \right )$的通解形式为

$$
\begin{aligned}
x &= a_1t_1M_1 + a_2t_2M_2 + \cdots + a_nt_nM_n + kM \\
&= kM + \sum_{i=1}^na_it_iM_i, \ k \in \mathbb{Z}
\end{aligned}
$$

在模$M_i$的意义下，方程组$\left ( S \right )$只有一个解:

$$
x \equiv \left ( a_1t_1M_1 + a_2t_2M_2 + \cdots + a_nt_nM_n \right ) \ mod \ M
$$

### 代码

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
        printf("Case %d: the next triple peak occurs in %d days.\n",
        ++kase, ans-d);
    }
    return 0;
}
```

## 欧拉函数

### 定义&通式

欧拉函数是小于等于 $n$ 的正整数中与 $n$ 互质的数的数目（$\varphi \left ( 1 \right )=1$）。  
通式：$\varphi \left ( x \right ) = x\left ( 1 - \frac{1}{p_1} \right )\left ( 1 - \frac{1}{p_2} \right )\left ( 1 - \frac{1}{p_3} \right )\cdots\left ( 1 - \frac{1}{p_n} \right )$  
性质: 若a与b**互质**,则$\varphi \left( ab \right) = \varphi \left(a \right)\varphi \left(b \right)$  
应用：欧拉降幂公式  
$a^b \equiv a^{b \  \% \  \varphi \left( n\right) + \varphi \left( n \right)} (mod\ n)\ (b > \varphi (n))$

### 代码

```c++
/*线性筛O(n)时间复杂度内筛出maxn内欧拉函数值*/
int m[maxn],phi[maxn],p[maxn],pt;
//m[i]是i的最小素因数，p是素数，pt是素数个数
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
    /*这里的phi[k]与phi[i]后面的∏(p[i]-1)/p[i]都一样
    （m[i]==p[j]）只差一个p[j]，就可以保证∏(p[i]-1)/p[i]前面也一样了*/
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

## 素数筛法

### 线形筛

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

### 复杂度 $O(n^{\frac{3}{4}})$

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
### 复杂度 $O(n^{\frac{2}{3}})$

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
        for(int j = 1; j <= PM; ++j) 
            phi[j][i] = phi[j][i - 1] - phi[j / prime[i]][i - 1];
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
    for(int i = pi[sqrt3(x)] + 1, ed = pi[sqrt2(x)]; i <= ed; ++i) 
        ans -= getpi(x / prime[i]) - i + 1;
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
        for (int j = i; j <= lim; j++) 
            sum -= lehmer_pi(w / prime[j]) - (j - 1);
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

## miller-rabin素性判断

### 代码

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
## 莫比乌斯函数

### 定义

$$ \mu = \begin{cases} 1 & n=1 \\ (-1)^k & n = p_1p_2\cdots p_k \\ 0 & other \end{cases}$$

### 莫比乌斯反演

$$f(n) = \sum_{d,n}g(d)=\sum_{d,n} g(\frac{n}{d})$$  
$$ g(n) = \sum_{d,n} \mu(d) f(\frac{n}{d}) = \sum_{d,n} \mu(\frac{n}{d})f(d) $$  
倍数形式只用把$\frac{n}{d}$变为$\frac{d}{n}$  

### 技巧
若$g(d)=[\frac n d]*[\frac m d]$之类的阶梯状函数  
记录$\mu$的前缀和

```c++
int d=1;
int ans=0;
while(d<=min(n,m))
{
    int last=min(n/(n/d),m/(m/d));
    ans+=(sum[last]-sum[d-1])*(n/d)*(m/d);
    d=last+1;
}
```

## 求原根

### 定义
给定一个数$n$，若存在一个与 $n$互素的 $a$,使得 $a^i(i=0,1,\cdots,\varphi(n))$在模$n$ 下两两不同,那么称$a$是$n$的一个原根。

### 性质
$1,2,4,p^n,2p^n$有原根，其中$p$是奇素数  
一个数$n$如果有原根，原根个数为 $\varphi(\varphi(n))$  
一个数$n$的全体原根的乘积模 $n$余1  
一个数$n$的全体原根的总和模 $n$余 $\mu(n-1)$(莫比乌斯函数)

### 头文件&全局变量
 
```c++
#include <algorithm>

const int maxn = 1e6 + 6;

long long m[maxn], phi[maxn], p[maxn], pt; 
//m[i]是i的最小素因数，p是素数，pt是素数个数
int n, T;
long long sum[maxn];
int prime[maxn], ptop;
bool book[maxn]={0};
int pr[maxn];
```
### 辅助函数

```c++
int gcd(int a, int b){
    return (b>0)?gcd(b,a%b):a;
}

void getPrime(int n) {
    ptop = 0;
    for (int i = 0; p[i] <= n; i++) {
        if (n % p[i] == 0) {
            prime[ptop++] = p[i];
            while (n % p[i] == 0) {
                n /= p[i];
            }
        }
    }
}

void make() {
    phi[1] = 1;
    int N = maxn;
    int k;
    for (int i = 2; i < N; i++) {
        if (!m[i]) //i是素数
            p[pt++] = m[i] = i, phi[i] = i - 1, book[i] = 1;
        for (int j = 0; j < pt && (k = p[j] * i) < N; j++) {
            m[k] = p[j];
            if (m[i] == p[j]) //为了保证以后的数不被再筛，要break
            {
                phi[k] = phi[i] * p[j];
                /*这里的phi[k]与phi[i]后面的∏(p[i]-1)/p[i]都一样
                （m[i]==p[j]）只差一个p[j]，
                就可以保证∏(p[i]-1)/p[i]前面也一样了*/
                break;
            } else
                phi[k] = phi[i] * (p[j] - 1); 
                //积性函数性质，f(i*k)=f(i)*f(k)
        }
    }
}

long long quickPowMod(long long a, int k, int mod) {
    long long ans = 1;
    while (k) {
        if (k & 1) {
            ans *= a;
            ans %= mod;
        }
        a *= a;
        a %= mod;
        k /= 2;
    }
    return ans;
}
```

### 核心代码

```c++
//判断是否有原根
bool hasPrimitiveRoot(int n){
    if(n == 2 || n == 4) return true;
    if(book[n]) return true;
    if(n % 2 == 0) n /= 2;
    for(int i=1; p[i] <= n; i++){
        if(n % p[i] == 0){
            while(n % p[i] == 0){
                n /= p[i];
            }
            if(n==1) return true;
            else return false;
        }
    }
    return false;
}

int cntPrimitiveRoot(int n) {
    if(!hasPrimitiveRoot(n)){
        printf("-1\n");
        return -1;
    }
    int cnt = 0;
    int phi_n = phi[n];
    getPrime(phi_n);
    //枚举
    if(n == 2){
        cnt = 1;
        printf("1\n");
        return 1;
    }
    for (int a = 2; a < n; a++) {
        //判断a是否为n的原根
        bool flag = true;
        if (quickPowMod(a, phi_n, n) != 1) continue;
        for (int i = 0; i < ptop; i++) {
            int k = phi_n / prime[i];
            if (quickPowMod(a, k, n) == 1) {
                flag = false;
                break;
            }
        }
        if (flag) {
            pr[cnt++] = a;
            break;
        }
    }
    
    for(int i=2; i<phi_n; i++){
        if(gcd(i, phi_n) == 1) pr[cnt++] = quickPowMod(pr[0], i, n);
    }

    sort(pr, pr+cnt);
    int prt = unique(pr, pr+cnt) - pr;
    for(int i=0; i<prt; i++){
        printf("%d", pr[i]);
        if(i!=prt-1)printf(" ");
        else printf("\n"); 
    }
    //printf("%d %lld\n", cnt, phi[phi[n]]);
    return phi[phi[n]];
}
```
### 用法

```c++
int main() {
    make();
    while (~scanf("%d", &n)) {
        cntPrimitiveRoot(n);
    }
    return 0;
}
```

## 高斯消元

### 异或-高斯消元

```c++
template <class TN>                    //TN为bitset
int xor_gauss(TN bits[], int n, int m) //n行m列
{
    int i = 0, k = m - 1; //先消高位
    //i枚举行,k枚举列
    while (i < n && k >= 0) {
        int nown = i;
        while (nown < n && !bits[nown].test(k))
            nown++;
        if (nown < n) {
            for (int j = nown + 1; j < n; j++)
                if (bits[j].test(k))
                    bits[j] ^= bits[nown];
            swap(bits[nown], bits[i]);
            i++;
        }
        k--;
    }
    //返回秩
    return i;
}
```

## 二次剩余

### 引用自
http://blog.csdn.net/xf_zhen/article/details/52097988

### 代码

```c++
#include <iostream>  
#include <cstdio>  
using namespace std;  

#define LL long long  

LL quick_mod(LL a,LL b,LL p)//快速幂  
{  
    LL ans = 1;  
    a%=p;  
    while(b)  
    {  
        if(b&1)  
        {  
            ans = ans * a%p;  
            b--;  
        }  
        b>>=1;  
        a = a*a%p;  
    }  
    return (ans+p)%p;  
}  
LL Legendre(LL a,LL p)//求勒让得符号（-1,0,1）这里-1返回p-1  
{  
    return quick_mod(a,(p-1)>>1,p);  
}  
struct T //二次域  
{  
    LL p,d;  
};  
LL w; //二次域第二个单位参数  

LL mod(LL t, LL p)  
{  
    t %=p;  
    if(t<0) t+=p;  
    return t;  
}  
T multi_er(T a,T b, LL p)//二次域乘法  
{  
    T ans;  
    ans.p = (a.p*b.p%p+ a.d*b.d%p*w%p)%p;  
    ans.d = (a.p*b.d%p +a.d*b.p%p)%p;  
    return ans;  
}  
T pow_er(T a,LL b,LL p)//二次域上的快速幂  
{  
     T ans;  
     ans.p=1;  
     ans.d =0;  

    while(b)  
    {  
        if(b&1)  
        {  
           ans = multi_er(ans, a, p);  
            b--;  
        }  
        b>>=1;  
        a = multi_er(a, a, p);  
    }  
    return ans;  
}  
int solve(int n,int p)  
{  
    if(p==2) return 1;  
    if(Legendre(n,p)+1 == p) return -1;  
   // printf("solve ...\n");  
    LL a = -1,t;  
    while(true)  
    {  
        a = rand()%p;  
        t = a*a -n;  
        w = mod(t,p);  
        if(Legendre(w,p)+1==p) break;  
    }  
    T tmp;  
    tmp.p = a;  
    tmp.d = 1;  
    T ans = pow_er(tmp,(p+1)>>1,p);  
    return ans.p;  
}  
int main()  
{  
   int cas;  
   scanf("%d",&cas);  
   while(cas--)  
   {  
       int n,p;  
       scanf("%d %d",&n,&p);  
       n %=p;  
       int a = solve(n,p);  
       if(a == -1)  
       {  
           puts("No root");  
           continue;  
       }  
       else  
       {  
           int b = p - a;  
           if(a==b)  
           {  
               printf("%d\n",a);  
           }  
           else  
           {  
               if(a>b) swap(a,b);  
               printf("%d %d\n",a,b);  
           }  
       }  
   }  
    return 0;  
}
```
## 离散对数

### 问题描述
给定B、N、P，求一个整数L满足 $B^L \equiv N \ (mod \ P)$

### 全局变量 & 宏

```c++
#include <stdio.h>
#include <string.h>
#include <math.h>
#define MOD 76543
int hs[MOD],head[MOD],next[MOD],id[MOD],top;
```

### 辅助函数

```c++
void insert(int x,int y)
{
    int k = x%MOD;
    hs[top] = x, id[top] = y, next[top] = head[k], head[k] = top++;
}
int find(int x)
{
    int k = x%MOD;
    for(int i = head[k]; i != -1; i = next[i])
        if(hs[i] == x)
            return id[i];
    return -1;
}
```

### BSGS算法

```c++
//a^ans = b (mod n)
int BSGS(int a,int b,int n)
{
    memset(head,-1,sizeof(head));
    top = 1;
    if(b == 1)return 0;
    int m = sqrt(n*1.0), j;
    long long x = 1, p = 1;
    for(int i = 0; i < m; ++i, p = p*a%n)insert(p*b%n,i);
    for(long long i = m; ;i += m)
    {
        if( (j = find(x = x*p%n)) != -1 )return i-j;
        if(i > n)break;
    }
    return -1;
}
```

### 用法

```c++
int main()
{
    int P,B,N;
    while(scanf("%d%d%d",&P,&B,&N) == 3)
    {
        int ans = BSGS(B,N,P);
        if(ans == -1)printf("no solution\n");
        else printf("%d\n",ans);
    }
    return 0;
}
```

## 扩展离散对数

### 问题描述
给定B、N、P，求一个整数L满足 $B^L \equiv N \ (mod \ P)$

### 全局变量&结构体&宏

```c++
#include <string.h>
#include <math.h>
 
#define MAXN 65536

typedef long long ll;
 
struct LINK{
    ll data;
    ll j;
    ll next;
}HASH_LINK[1000000];
ll ad, head[MAXN];
```

### 辅助函数

```c++
ll Gcd(ll a, ll b){
return b ? Gcd(b, a % b) : a;
}
 
ll Ext_Gcd(ll a, ll b, ll &x, ll &y){
    if(!b){
       x = 1; y = 0;
       return a;
    }
    ll r = Ext_Gcd(b, a % b, x, y);
    ll t = x; x = y; y = t - a / b * y;
    return r;
}
 
ll POWER(ll a, ll b, ll c){
    ll ans = 1;
    while(b){
       if(b & 1) ans = ans * a % c;
       a = a * a % c;
       b >>= 1;
    }
    return ans;
}
 
void clear(){
    memset(head, -1, sizeof(head));
    ad = 0;
}
 
ll hash(ll a){
    return a % MAXN;
}
 
void INSERT_HASH(ll i, ll buf){
    ll hs = hash(buf), tail;
    for(tail = head[hs]; ~tail; tail = HASH_LINK[tail]. next)
       if(buf == HASH_LINK[tail]. data) return;
    HASH_LINK[ad]. data = buf;
    HASH_LINK[ad]. j    = i;
    HASH_LINK[ad]. next = head[hs];
    head[hs] = ad ++;
}
```

### 扩展BSGS算法

```c++
// a^ans = b(mod c)
ll bady_step_giant_step(ll a, ll b, ll c){
    ll i, buf, m, temp, g, D, x, y, n = 0;
    for(i = 0, buf = 1; i < 100; i ++, buf = buf * a % c)
       if(buf == b) return i;
    D = 1;
    while((g = Gcd(a, c)) != 1){
       if(b % g) return -1; // g | b 不满足，则说明无解
       b /= g;
       c /= g;
       D = D * a / g % c;
       ++ n;
    }
    clear();
    m = ceil(sqrt((long double) c));
    for(i = 0, buf = 1; i <= m; buf = buf * a % c, i ++) INSERT_HASH(i, buf);
    for(i = 0, temp = POWER(a, m, c), buf = D; i <= m; i ++, buf = temp * buf % c){
       Ext_Gcd(buf, c, x, y);
       x = ((x * b) % c + c) % c;
       for(ll tail = head[hash(x)]; ~tail; tail = HASH_LINK[tail].next)
           if(HASH_LINK[tail]. data == x) return HASH_LINK[tail].j + n + i * m;
    }
    return -1;
}
```

### 用法

```c++
//k^ans = n(mod p) 
int main(){
    ll k, p, n, ans;
    while(~scanf("%lld %lld %lld", &k, &p, &n)){
       if(n >= p){ printf("Orz,I can’t find D!\n"); continue; }
       ans = bady_step_giant_step(k, n, p);
       ans == -1 ? printf("Orz,I can’t find D!\n") : printf("%lld\n", ans);
    }
    return 0;
}
```

## Pollard_rho 质因子分解

### 环境

```c++
#include<cstdlib>
#include<ctime>

typedef long long LL;
#define MAXN 10000

LL factor[MAXN];
int tot;
const int S=20;
```

### 辅助函数

```c++
LL muti_mod(LL a,LL b,LL c){    //返回(a*b) mod c,a,b,c<2^63
    a%=c;
    b%=c;
    LL ret=0;
    while (b){
        if (b&1){
            ret+=a;
            if (ret>=c) ret-=c;
        }
        a<<=1;
        if (a>=c) a-=c;
        b>>=1;
    }
    return ret;
}

LL pow_mod(LL x,LL n,LL mod){  //返回x^n mod c ,非递归版
    if (n==1) return x%mod;
    int bit[90],k=0;
    while (n){
        bit[k++]=n&1;
        n>>=1;
    }
    LL ret=1;
    for (k=k-1;k>=0;k--){
        ret=muti_mod(ret,ret,mod);
        if (bit[k]==1) ret=muti_mod(ret,x,mod);
    }
    return ret;
}

bool check(LL a,LL n,LL x,LL t){   //以a为基，n-1=x*2^t，检验n是不是合数
    LL ret=pow_mod(a,x,n),last=ret;
    for (int i=1;i<=t;i++){
        ret=muti_mod(ret,ret,n);
        if (ret==1 && last!=1 && last!=n-1) return 1;
        last=ret;
    }
    if (ret!=1) return 1;
    return 0;
}

bool Miller_Rabin(LL n){
    LL x=n-1,t=0;
    while ((x&1)==0) x>>=1,t++;
    bool flag=1;
    if (t>=1 && (x&1)==1){
        for (int k=0;k<S;k++){
            LL a=rand()%(n-1)+1;
            if (check(a,n,x,t)) {flag=1;break;}
            flag=0;
        }
    }
    if (!flag || n==2) return 0;
    return 1;
}

LL gcd(LL a,LL b){
    if (a==0) return 1;
    if (a<0) return gcd(-a,b);
    while (b){
        LL t=a%b; a=b; b=t;
    }
    return a;
}
```

### 核心代码

```c++
//找出任意质因数
LL Pollard_rho(LL x,LL c){ 
    LL i=1,x0=rand()%x,y=x0,k=2;
    while (1){
        i++;
        x0=(muti_mod(x0,x0,x)+c)%x;
        LL d=gcd(y-x0,x);
        if (d!=1 && d!=x){
            return d;
        }
        if (y==x0) return x;
        if (i==k){
            y=x0;
            k+=k;
        }
    }
}

//递归进行质因数分解N
void findfac(LL n){           
    if (!Miller_Rabin(n)){
        factor[tot++] = n;
        return;
    }
    LL p=n;
    while (p>=n) p=Pollard_rho(p,rand() % (n-1) +1);
    findfac(p);
    findfac(n/p);
}
```

### 用法

初始化tot=0，调用findfac(n)，factor[0~tot-1]中为n的全部质因子

## FFT

### 定义&公式

时间复杂度$\Theta \left ( nlog \left ( n \right ) \right )$  
FFT（Fast Fourier Transformation）是离散傅氏变换（DFT）的快速算法。即为快速傅氏变换。它是根据离散傅氏变换的奇、偶、虚、实等特性，对离散傅立叶变换的算法进行改进获得的。

令 $W_N=e^{-i \frac{2\pi }{N} }$  
DFT : $X \left( k \right) = \sum_{n=0}^{N-1}x \left( n \right)W_N^{nk} ,k=0,1,\cdots,N-1$  
IDFT: $x \left( n \right) = \frac{1}{N} \sum_{k=0}^{N-1}X \left( k \right)W_N^{-nk} ,n=0,1,\cdots,N-1$

### 推导

若N为偶数,且$k<\frac{n}{2}$,则

$$
\begin{aligned}
X\left ( k \right ) 
 &=\sum_{n=0}^{N-1}x \left( n \right ) e^{-i2\pi kn/N} \\ 
 &=\sum_{m=0}^{\frac{N}{2}-1}x \left( 2m \right ) e^{-i2\pi km/\frac{N}{2}}
   +
   e^{-i2\pi k/N}\sum_{m=0}^{\frac{N}{2}-1}x \left( 2m+1 \right ) e^{-i2\pi km/\frac{N}{2}} \\ 
 &=X_1\left( k \right ) + W_N^{k}X_2\left( k \right )
\end{aligned}
$$

其中对偶数项做FFT得到$X_1\left( k \right )$,对奇数项做FFT得到$X_2\left( k \right )$  
其中$X_1\left( k \right )$为第k个偶数项的DFT值  
由于$e^{i\pi 2n}=1,n\in \mathbb{Z}$且$e^{i\pi \left(2n+1\right)}=-1,n\in \mathbb{Z}$

$$
\begin{aligned}
X\left ( k+\frac{N}{2} \right ) 
 &=e^{-i\pi 2m}\sum_{m=0}^{\frac{N}{2}-1}x \left( 2m \right ) e^{-i2\pi km/\frac{N}{2}}
   +
   e^{-i\pi \left(2m+1\right)}\cdot e^{-i2\pi k/N}\sum_{m=0}^{\frac{N}{2}-1}x \left( 2m+1 \right ) e^{-i2\pi km/\frac{N}{2}} \\ 
 &=X_1\left( k \right ) - W_N^{k}X_2\left( k \right )
\end{aligned}
$$

这样$X$每一项都可通过$X_1$和$X_2$求出

### 头文件&全局变量

```c++
#include <complex>
#include <cmath>

typedef complex<double> Complex;
const double PI=acos(-1.0);
```

### 辅助函数

```c++
int lowbit(int x)
{
    return x&-x;
}
```

### 可选函数

```c++
int getN(int n)//返回大于等于n的最小的2^t的值
{
    int nn=1;
    while(nn<n)nn<<=1;
    return nn;
}

//规则化数组,使长度n变为2^t的形式,并在多余部分填充0
void normalize(Complex x[],int &n)
{
    int nn=getN(n);
    fill(x+n,x+nn,0);
    n=nn;
}
```

### 核心代码

```c++
//n必须为2^t,若不是,用normalize转换
//on==1为DFT,on==-1为IDFT
void fft(Complex x[],int n,int on=1)
{
    //循环模拟递归向下数组最终样子,位置i和(i二进制对称的数)互换就行了
    for(int i=1;i<n;i++)
    {
        int j=0;
        for(int k=i;k;k^=lowbit(k))
            j|=n/lowbit(k);
        j>>=1;
        if(i<j)swap(x[i],x[j]);
    }

    //循环模拟递归回溯
    double ww=-2.0*PI*on;//中间变量,字面上意思
    for(int i=1,j=2;j<=n;i<<=1,j<<=1)
    {
        Complex wn(cos(ww/j),sin(ww/j));
        for(int k=0;k<n;k+=j)
        {
            Complex w(1.0);
            for(int t=k;t<k+i;t++)
            {
                Complex a(x[t]);
                Complex b(x[t+i]*w);
                x[t]=a+b;
                x[t+i]=a-b;
                w*=wn;
            }
        }
    }

    if(on==-1)
        for(int i=0;i<n;i++)
            x[i]/=n;
}
```

### 使用例子

```c++
//求两个大数相乘a[0~n-1]*b[0~m-1]
int nn=getN(max(n,m))<<1;//大于等于位数的两倍
fill(a+n,a+nn,0);
fill(b+m,b+nn,0);//末尾填充0
n=nn;
fft(a,n);
fft(b,n);
for(int i=0;i<n;i++)
    a[i]*=b[i];
fft(a,n,-1);
for(int i=0;i<n;i++)
    ans[i]=floor(a[i].real()+0.5);
for(int i=0;i<n;i++)
{
    ans[i+1]+=ans[i]/10;
    ans[i]%=10;
}
```

# 组合数学

## 第一类斯特林数

### 定义

把n个元素分成k个环排列的方法数。  
边界值：$S_1(n,0)=0$ $S_1(n,n)=1$  
递推式：$S_1(n+1, k) = S_1(n, k-1)+nS_1(n,k)$

### 递推初始化

```c++
void init(){
    for(ll i=1; i<MAXN; i++){
        S[i][0] = 0;
        S[i][i] = 1;
        for(ll j=1; j<i; j++){
            S[i][j] = (((i-1)*S[i-1][j])%MOD+S[i-1][j-1])%MOD;
        }
    }
}
```

## 第二类斯特林数

### 定义

把n个元素分为k个非空子集的方案数  
边界值：$S_2(n,0)=0$ $S_2(n,n)=1$  
递推式：$S_2(n+1, k) = S_2(n, k-1)+kS_2(n,k)$

### 递推初始化

```c++
void init() {
    for (int i = 1; i <= 100; i++) {
        S2[i][0] = 0, S2[i][i] = 1;
        for (int j = 1; j < i; j++) {
            S2[i][j] = (S2[i - 1][j - 1] + (j * S2[i - 1][j]) % MOD) % MOD;
        }
    }
}
```

### 常见用法

把n个小球放入k个盒子中
*  盒子不能为空，不可区分： $S_2(n, k)$
*  盒子能为空，不可区分： $\sum_{i=1}^k{S_2(n, i)}$
*  盒子不能为空，可区分： $k!S_2(n, k)$

## 分拆数

### 定义

把正整数n拆分为至多k个正整数之和的方案数。  
把正整数n拆分为不超过k的若干个正整数之和的方案数。

### DP递推

复杂度： $O(n^2)$  
初始值：$dp[0][0]=1\ \ \ dp[n][0]=0$  
递推：

$$
\begin{aligned}
& dp[n][m]=dp[n][m-1]+dp[n-m][m]\ \ \ (n \geqslant  m) \\
& dp[n][m] = dp[n][n] \ \ \ (n < m)
\end{aligned}
$$

### 母函数

复杂度： $O(n\sqrt{n})$ 

#### 推导

正整数n的划分可以用母函数表示为

$$
\sum_{i=0}^{\infty}p(i)x^i = (1+x^1+x^2+\cdots)(1+x^2+x^4+\cdots)(1+x^3+x^6+\cdots)\cdots
$$

对其进行等比数列求和，可得

$$
\sum_{i=0}^{\infty}p(i)x^i = \prod_{k=1}^{\infty}(\frac{1}{1-x^k})
$$

将$\prod^{\infty}_{k=1}(\frac{1}{1-x^k})$展开可得：

$$
(1-x)(1-x^2)(1-x^3)\cdots = 1 -x - x^2 + x^5 + x^7 - x^{12} - x^{15} + x^{22} + x^{26} + \cdots
$$

观察可得：  
指数即扩展五边形数（下一列为对应下标）：

$$
\begin{matrix}
& 0 & 1 & 2 & 5 & 7 & 12 & 15 \\
& 0 & 1 & -1 & 2 & -2 & 3 & -3 
\end{matrix}
$$

*扩展五边形数：$Five(x)=\frac{3x^2-x}{2}$

令$Q(x)=\prod^{\infty}_{k=1}{(1-x^k)}$，所以$Q(x)P(x) = 1$，即：

$$
(p(0)x^0+p(1)x^1+p(2)x^2+\cdots)(1-x-x^2+x^5+x^7-x^22\cdots) = 1
$$

所以$p(k) = p(k-1) + p(k-2) - p(k-5) \cdots$

#### 应用

要求拆分的数中每个数出现的次数不能大于等于k此，则

$$
\begin{aligned}
P_k(x) \ 
& = (1+x+x^2+ \cdots +x^{k-1})(1+x^2+x^4+\cdots + c^{2(k-1)})\cdots \\
& = \prod^{\infty}_{i=1}{\frac {1-x^{ki}} {1-x^{i}}} \\
& = \frac {Q(x^k)} {Q(x)} \\
& = Q(x^k)P(x)
\end{aligned}
$$

例如：当$n=8, k=4$时：

$$
Q(x^4)P(x) = (1-x^4-x^8+x^{12}+\cdots)P(x)
$$

满足指数为8的项的和：$1\times 22x^8 - x^4 \times 5 x^4 - x^8 \times 1 = 16 x^8$

*  示例代码：

```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<string>
#include<cmath>
#include<algorithm>

using namespace std;

typedef __int64 LL;
const int Maxn=100010;
const LL MOD=1000000007;
LL Q[Maxn],P[Maxn];
LL GetQ(LL x)
{
    LL ans=(LL)x*x*3-x;
    return (ans/2)%MOD;
}
void _init()
{
    Q[0]=0;
    for(int i=1;i<Maxn;i++)
    {
        if(i&1) Q[i]=GetQ(i/2+1);
        else Q[i]=GetQ(i/2*(-1));
    }
    P[0]=P[1]=1;
    for(int i=2;i<Maxn;i++)
    {
        for(int j=1;;j++)
        {
            if(Q[j]>i) break;
            int t=j;
            if(t&1) t=t/2+1;
            else t=t/2;
            if(t&1)
                P[i]=(P[i]+P[i-Q[j]]);
            else
                P[i]=(P[i]-P[i-Q[j]]);
            if(P[i]>=MOD) P[i]%=MOD;
            if(P[i]<0) P[i]+=MOD;
        }
    }
}
void solved(LL n,LL k)
{
    LL ans=0;
    for(int i=0;;i++)
    {
        if(Q[i]*k>n) break;
        int t=i;
        if(t&1) t=t/2+1;
        else t=t/2;
        if(t&1) ans=(ans-P[n-Q[i]*k]);
        else ans=(ans+P[n-Q[i]*k]);
        if(ans>=MOD) ans%=MOD;
        if(ans<0) ans+=MOD;
    }
    printf("%I64d\n",ans);
}
int main()
{
    _init();
    int T;
    LL n,k;
    scanf("%d",&T);
    while(T--)
    {
        scanf("%I64d%I64d",&n,&k);
        solved(n,k);
    }
    return 0;
}
```

### 性质
  
限制拆分结论：  
*  “将一个正整数n拆为若干个两两不同的正整数之和” 与 “将一个正整数n拆为若干个奇数之和”的方案数相同
*  “把正整数n拆分为至多k个正整数之和” 与 “把正整数n拆分为不超过k的若干个正整数之和”的方案数相同
*  “把正整数n拆分为若干个相等的正整数之和”的方案数为n的因子个数

### 代码

#### 环境

```c++
typedef long long LL;

const int MAXN = 100010;
const LL MOD = 1000000007;

LL dp[MAXN]={0};
```

#### 辅助函数

```c++
LL Five(LL x) { //计算正五边形数
    LL ans = 3 * x * x - x;
    return ans / 2;
}

bool change(int i, int j) {
    LL k = Five(j);
    if (k > i) return false;

    if (j % 2 == 0) //判断奇偶从而判断系数为正还是为负
        dp[i] = (dp[i] - dp[i - k]);
    else
        dp[i] = (dp[i] + dp[i - k]);
    (dp[i] += MOD) %= MOD;

    return true;
}
```

#### 递推母函数初始化

```c++
void init() {
    dp[0] = 1;
    for (int i = 1; i < MAXN; i++)
        for (int j = 1;change(i,j)&&change(i,-j); j++);
}
```

#### 用法

调用init()后，dp数组内dp[i]表示将正整数划分为若干个正整数之和的方案数。

## 错排

### 定义

将n个不同的元素重新排列后每个元素都不放在自己原来位置上的方法数

### 递推

$$F(i) = (i-1)\times (\ F(i-1)+F(i-2)\ )$$

### 代码

```c++
//初始化错排函数，循环节为2*mod
void initF(LL n, LL mod) {
    F[0] = 1, F[1] = 0;
    for (int i = 2; i < 2 * mod; i++)
        F[i] = LL(i - 1) * (F[i - 1] + F[i - 2]) % mod;
}
```

## 扩展Lucas定理

快速求组合数取模

### 预处理加速

空间复杂度O(MOD)

#### 环境

```c++
typedef long long LL;

const int MAXF = 200003;
const int MAXMOD = 100003;

LL num[MAXMOD], factor[MAXMOD], power[MAXMOD], tot, pri[MAXMOD];
```

#### 辅助函数

```c++
//对模数分解质因子，记录每个质因子factor的指数power和对应的值num
void div(LL mod) {
    tot = 0;
    LL tmp = mod;
    for (int i = 2; i * i <= mod; i++) {
        if (tmp % i != 0) continue;
        factor[++tot] = i, num[tot] = 1, power[tot] = 0;
        while (tmp % i == 0)
            num[tot] *= i, tmp /= i, power[tot]++;
    }
    if (tmp != 1) factor[++tot] = tmp, num[tot] = tmp, power[tot] = 1;
}

//扩展欧几里得求逆元
void ext_gcd(LL a, LL b, LL &x, LL &y) {
    if (b == 0) {
        x = 1, y = 0;
        return;
    }
    ext_gcd(b, a % b, x, y);
    LL tp = x;
    x = y;
    y = tp - a / b * y;
}

LL quickPow(LL a, LL b, LL mod) {
    LL res = 1;
    for (; b > 0; b /= 2) {
        if (b % 2 == 1) res = res * a % mod;
        a = a * a % mod;
    }
    return res;
}

LL calc(LL n, LL p, LL _p) {
    if (n < _p) return pri[n];
    return (LL(pri[n % p]) * calc(n / _p, p, _p) % p
        * quickPow(pri[p - 1], n / p, p) % p + p) % p;
}

LL count(LL n, LL p) {
    if (n < p) return 0;
    return n / p + count(n / p, p);
}

LL C(LL n, LL k, LL p, LL _p, LL t) {
    pri[0] = 1;
    for (int i = 1; i <= p; i++)
        if (i % _p != 0)
            pri[i] = LL(pri[i - 1]) * i % p;
        else
            pri[i] = pri[i - 1];
    LL t3 = count(n, _p) - count(k, _p) - count(n - k, _p);
    if (t3 >= t) return 0;
    LL t1 = calc(n, p, _p), t2 = calc(k, p, _p) * calc(n - k, p, _p) % p;
    LL x, y;
    ext_gcd(t2, p, x, y);
    return (t1 * x % p * quickPow(_p, t3, p) + p) % p;
}
```

#### 核心代码

```c++
//求n个元素中任意取m个的方法数模mod
LL exLucas(LL n, LL m, LL mod) {
    LL res = 0, tp;
    LL x, y;
    for (int i = 1; i <= tot; i++) {
        ext_gcd(mod / num[i], num[i], x, y);
        tp = C(n, m, num[i], factor[i], power[i]);
        res += x * tp * (mod / num[i]) % mod, res %= mod;
    }
    return (res + mod) % mod;
}
```
#### 用法

调用$div(mod)$对模数质因数分解，随后调用$exLucas(n, m, mod)$求n个元素中任意取m个的方法数模mod

## 二项式反演

### 定义

若

$$
f(n) = \sum^n_{j=0} C^j_n \ g(j)
$$

可由二项式反演得

$$
g(n) = \sum ^n _{i=0} (-1)^i C^i_n \ f(n-i)
$$

## 卡特兰数

### 定义

*  n个元素进出栈方案数
*  n个左括号与n个右括号的匹配方案数
*  一个正n多边形用n-3条不相交的对角线划分成n-2个三角形的方案数
*  一棵体积为n的有根二叉树有多少种形态
*  $\cdots \cdots$

$$
Catalon(n) = \frac{(2n)!}{n!(n+1)!} \ \ (0! = 1)
$$

### 代码

#### 环境

```c++
typedef long long ll;

const int MAXN = 2000010;

int prime[MAXN], size[MAXN]={0}, p, n, len, num[MAXN];
bool isnot[MAXN];
```

#### 辅助函数

```c++
void init() { //素数筛
    for (int i = 2; i <= (n << 1); i++) {
        if (!isnot[i]) prime[++len] = i, num[i] = len;
        for (int j = 1; prime[j] * i <= (n << 1); j++) {
            isnot[prime[j] * i] = 1, num[prime[j] * i] = j;
            if (i % prime[j] == 0) break;
        }
    }
}

void div(int x, int s) {
    while (x != 1) {
        size[num[x]] += s;
        x /= prime[num[x]];
    }
}

inline void quickPow(long long &ans, int x, int y) {
    while (y) {
        if (y & 1) (ans *= x) %= p;
        y >>= 1;
        (x *= x) %= p;
    }
}
```

#### 核心代码

```c++
ll catalon(ll n) {
    for (int i = n + 2; i <= (n << 1); i++)
        div(i, 1);
    for (int i = 2; i <= n; i++)
        div(i, -1);
    long long ans = 1;
    for (int i = 1; i <= len; i++)
        if (size[i])
            quickPow(ans, prime[i], size[i]);
    return ans;
}
```

# 计算几何

> 来自[ACM计算几何模板__HIT_jerrybond](https://wenku.baidu.com/view/194472899e314332396893b0.html?qq-pf-to=pcqq.c2c)(**Author:jerrybond**)  
> 手动转码

## 几何公式

### 三角形

1.  半周长 $P=\frac{a+b+c}{2}$

2.  面积 $S=\frac{aH_a}{2}=\frac{ab\cdot sin(C)}{2}=\sqrt{P(P-a)(P-b)(P-c)}$

3.  中线 $M_a=\frac{\sqrt{2(b^2+c^2)-a^2}}{2}=\frac{\sqrt{b^2+c^2+2bc\cdot cos(A)}}{2}$

4.  角平分线 $T_a=\frac{\sqrt{bc((b+c)^2-a^2)}}{b+c}=\frac{2bc\cdot cos(A/2)}{b+c}$

5.  高线 $H_a=b\cdot sin(C)=c\cdot sin(B)=\sqrt{b^2-(\frac{a^2+b^2-c^2}{2a}))^2}$

6.  内切圆半径

$$
\begin{aligned}
r=\frac{S}{P}&=a\cdot sin(\frac{B}{2})sin(\frac{C}{2})/sin(\frac{B+C}{2})\\
             &=4R\cdot sin(\frac{A}{2})sin(\frac{B}{2})sin(\frac{C}{2})\\
              &=\sqrt{(P-a)(P-b)(P-c)/P} \\
             &=P\cdot tan(\frac{A}{2})tan(\frac{B}{2})tan(\frac{C}{2})
\end{aligned}
$$

7.  外接圆半径 $R=\frac{abc}{4S}=\frac{a}{2sin(A)}=\frac{b}{2sin(B)}=\frac{c}{2sin(C)}$

### 四边形

> $D_1$,$D_2$ 为对角线,$M$ 对角线中点连线,$\theta$ 为对角线夹角

1.  $a^2+b^2+c^2+d^2=D_1^2+D_2^2+4M^2$

2.  $S=\frac{D_1D_2sin(\theta)}{2}$

> (以下对圆的内接四边形)

1.  $ac+bd=D_1D_2$

2.  $S=\sqrt{(P-a)(P-b)(P-c)(P-d)}$,P 为半周长

### 正 n 边形

> $R$ 为外接圆半径,$r$ 为内切圆半径

1.  中心角 $\theta=2\pi/n$

2.  内角 $C=(n-2)\pi/n$

3.  边长 $a=2\sqrt{R^2-r^2}=2R\cdot sin(\frac{\theta}{2})=2r\cdot tan(\frac{\theta}{2})$

4.  面积 $S=\frac{nar}{2}=nr^2tan(\frac{\theta}{2})=\frac{nR^2sin(\theta)}{2}=\frac{na^2}{4tan(\theta/2)}$

### 圆

1.  弧长 $l=r\theta$

2.  弦长 $a=2\sqrt{2hr-h^2}=2r\cdot sin(\frac{\theta}{2})$

3.  弓形高 $h=r-\sqrt{r^2-\frac{a^2}{4}}=r(1-cos(\frac{\theta}{2}))=\frac{a\cdot tan(\frac{\theta}{4})}{2}$

4.  扇形面积 $S_1=\frac{rl}{2}=\frac{r^2\theta}{2}$

5.  弓形面积 $S_2=\frac{rl-a(r-h)}{2}=\frac{r^2(\theta-sin(\theta))}{2}$

### 棱柱

1.  体积 $V=Ah$,$A$ 为底面积,$h$ 为高

2.  侧面积 $S=lp$,$l$ 为棱长,$p$ 为直截面周长

3.  全面积 $T=S+2A$

### 棱锥

1.  体积 $V=\frac{Ah}{3}$,$A$ 为底面积,$h$ 为高(以下对正棱锥)

2.  侧面积 $S=\frac{lp}{2}$,$l$ 为斜高,$p$ 为底面周长

3.  全面积 $T=S+A$

### 棱台

1.  体积 $V=\frac{\left(A_1+A_2+\sqrt{A_1A_2}\right)h}{3}$,$A_1$.$A_2$ 为上下底面积,$h$ 为高(以下为正棱台)

2.  侧面积 $S=\frac{\left(p_1+p_2\right)l}{2}$,$p_1$.$p_2$ 为上下底面周长,$l$ 为斜高

3.  全面积 $T=S+A_1+A_2$

### 圆柱

1.  侧面积 $S=2\pi rh$

2.  全面积 $T=2\pi r(h+r)$

3.  体积 $V=\pi r^2h$

### 圆锥

1.  母线 $l=\sqrt{h^2+r^2}$

2.  侧面积 $S=\pi rl$

3.  全面积 $T=\pi r(l+r)$

4.  体积 $V=\frac{\pi r^2h}{3}$

### 圆台

1.  母线 $l=\sqrt{h^2+(r_1-r_2)^2}$

2.  侧面积 $S=\pi(r_1+r_2)l$

3.  全面积 $T=\pi r_1(l+r_1)+\pi r_2(l+r2)$

4.  体积 $V=\frac{\pi\left(r_1^2+r_2^2+r_1r_2\right)h}{3}$

### 球

1.  全面积 $T=4\pi r^2$

2.  体积 $V=4\pi r^3/3$

### 球台

1.  侧面积 $S=2\pi rh$

2.  全面积 $T=\pi (2rh+r_1^2+r_2^2)$

3.  体积 $V=\frac{\pi h\left(3\left(r_1^2+r_2^2\right)+h^2\right)}{6}$

### 球扇形

1.  全面积 $T=\pi r(2h+r_0)$,$h$ 为球冠高,$r_0$ 为球冠底面半径

2.  体积 $V=\frac{2\pi r^2h}{3}$

## 直线与线段

### 预备函数

```c++
//结构定义与宏定义
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>

const double eps=1e-8

bool zero(int x)
{
    return fabs(x)<eps;
}

struct point{
    double x,y;
};

struct line{
    point a,b;
};

//计算 cross product (P1-P0)x(P2-P0)
double xmult(point p1,point p2,point p0)
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}

double xmult(double x1,double y1,double x2,double y2,double x0,double y0)
{
    return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
}

//计算 dot product (P1-P0).(P2-P0)
double dmult(point p1,point p2,point p0)
{
    return (p1.x-p0.x)*(p2.x-p0.x)+(p1.y-p0.y)*(p2.y-p0.y);
}

double dmult(double x1,double y1,double x2,double y2,double x0,double y0)
{
    return (x1-x0)*(x2-x0)+(y1-y0)*(y2-y0);
}

//两点距离
double distance(point p1,point p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}

double distance(double x1,double y1,double x2,double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
```

## 判三点是否共线

```c++
int dots_inline(point p1,point p2,point p3)
{
    return zero(xmult(p1,p2,p3));
}
```

### 判点是否在线段上

```c++
//判点是否在线段上,包括端点（下面为两种接口模式）
int dot_online_in(point p,line l)
{
    return zero(xmult(p,l.a,l.b))
            &&(l.a.x-p.x)*(l.b.x-p.x)<eps
            &&(l.a.y-p.y)*(l.b.y-p.y)<eps;
}

int dot_online_in(point p,point l1,point l2)
{
    return zero(xmult(p,l1,l2))
        &&(l1.x-p.x)*(l2.x-p.x)<eps
        &&(l1.y-p.y)*(l2.y-p.y)<eps;
}

//判点是否在线段上,不包括端点
int dot_online_ex(point p,line l)
{
    return dot_online_in(p,l)&&(!zero(p.x-l.a.x)||!zero(p.y-l.a.y))
        &&(!zero(p.x-l.b.x)||!zero(p.y-l.b.y));
}
```

### 判断两点在线段的同一侧

```c++
//判两点在线段同侧,点在线段上返回 0
int same_side(point p1,point p2,line l)
{
    return xmult(l.a,p1,l.b)*xmult(l.a,p2,l.b)>eps;
}

int same_side(point p1,point p2,point l1,point l2)
{
    return xmult(l1,p1,l2)*xmult(l1,p2,l2)>eps;
}
```

### 判断两点是否在线段的异侧

```c++
//判两点在线段异侧,点在线段上返回 0
int opposite_side(point p1,point p2,line l)
{
    return xmult(l.a,p1,l.b)*xmult(l.a,p2,l.b)<-eps;
}

int opposite_side(point p1,point p2,point l1,point l2)
{
    return xmult(l1,p1,l2)*xmult(l1,p2,l2)<-eps;
}
```

### 求点关于直线的对称点

-   **点关于直线的对称点 // by lyt**

-   **缺点：用了斜率**

-   **也可以利用"点到直线上的最近点"来做，避免使用斜率。**

```c++
point symmetric_point(point p1, point l1, point l2)
{
    point ret;
    if (l1.x > l2.x - eps && l1.x < l2.x + eps)
    {
        ret.x = (2 * l1.x - p1.x);
        ret.y = p1.y;
    }
    else
    {
        double k = (l1.y - l2.y ) / (l1.x - l2.x);
        ret.x = (2*k*k*l1.x + 2*k*p1.y - 2*k*l1.y - k*k*p1.x + p1.x)
            / (1 + k*k); ret.y = p1.y - (ret.x - p1.x ) / k;
    }
    return ret;
}
```

### 判断两线段是否相交

#### 常用版

```c++
//定义点
struct Point
{
    double x;
    double y;
};

typedef struct Point point;

//叉积
double multi(point p0, point p1, point p2)
{
    return ( p1.x - p0.x )*( p2.y - p0.y)
          -( p2.x - p0.x )*( p1.y - p0.y);
}

//相交返回 true,否则为 false,接口为两线段的端点
bool isIntersected(point s1,point e1, point s2,point e2)
{
    return (max(s1.x,e1.x) >= min(s2.x,e2.x)) &&
        (max(s2.x,e2.x) >= min(s1.x,e1.x)) &&
        (max(s1.y,e1.y) >= min(s2.y,e2.y)) &&
        (max(s2.y,e2.y) >= min(s1.y,e1.y)) &&
        (multi(s1,s2,e1)*multi(s1,e1,e2)>0) &&
        (multi(s2,s1,e2)*multi(s2,e2,e1)>0);
}
```

#### 不常用版

```c++
//判两线段相交,包括端点和部分重合
int intersect_in(line u,line v)
{
    if (!dots_inline(u.a,u.b,v.a)||!dots_inline(u.a,u.b,v.b))
        return !same_side(u.a,u.b,v)&&!same_side(v.a,v.b,u);
    return dot_online_in(u.a,v)
        ||dot_online_in(u.b,v)
        ||dot_online_in(v.a,u)
        ||dot_online_in(v.b,u);
}

int intersect_in(point u1,point u2,point v1,point v2)
{
    if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
        return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2);
    return dot_online_in(u1,v1,v2)
        ||dot_online_in(u2,v1,v2)
        ||dot_online_in(v1,u1,u2)
        ||dot_online_in(v2,u1,u2);
}

//判两线段相交,不包括端点和部分重合
int intersect_ex(line u,line v)
{
    return opposite_side(u.a,u.b,v)&&opposite_side(v.a,v.b,u);
}

int intersect_ex(point u1,point u2,point v1,point v2)
{
    return opposite_side(u1,u2,v1,v2)&&opposite_side(v1,v2,u1,u2);
}
```

### 求两条直线的交点

```c++
//计算两直线交点,注意事先判断直线是否平行!
//线段交点请另外判线段相交(同时还是要判断是否平行!)
point intersection(point u1,point u2,point v1,point v2)
{
    point ret=u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
        /((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x+=(u2.x-u1.x)*t;
    ret.y+=(u2.y-u1.y)*t;
    return ret;
}
```

### 点到直线的最近距离

```c++
point ptoline(point p,point l1,point l2)
{
    point t=p;
    t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
    return intersection(p,t,l1,l2);
}
```

### 点到线段的最近距离

```c++
point ptoseg(point p,point l1,point l2)
{
    point t=p;
    t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
    if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
        return distance(p,l1)<distance(p,l2)?l1:l2;
    return intersection(p,t,l1,l2);
}
```

## 多边形

### 预备浮点函数

```c++
#include <stdlib.h>
#include<stdio.h>
#include<string.h>
#include <math.h>

#define MAXN 1000
//offset 为多变形坐标的最大绝对值
#define offset 10000
#define eps 1e-8
//浮点数判 0
#define zero(x) (((x)>0?(x):-(x))<eps)
//浮点数判断符
#define _sign(x) ((x)>eps?1:((x)<-eps?2:0))

//定义点
struct point
{
    double x,y;
}pt[MAXN ];

//定义线段
struct line
{
    point a,b;
};

//叉积
double xmult(point p1,point p2,point p0)
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
```

### 判定是否是凸多边形

```c++
//判定凸多边形,顶点按顺时针或逆时针给出,允许相邻边共线,是凸多边形返回1，否则返回0
int is_convex(int n,point* p)
{
    int i,s[3]={1,1,1};
    for (i=0;i<n&&s[1]|s[2];i++)
        s[_sign(xmult(p[(i+1)%n],p[(i+2)%n],p[i]))]=0;
    return s[1]|s[2];
}

//判凸边行，顶点按顺时针或逆时针给出,不允许相邻边共线,是凸多边形返回1，否则返回 0
int is_convex_v2(int n,point* p)
{
    int i,s[3]={1,1,1};
    for (i=0;i<n&&s[0]&&s[1]|s[2];i++)
        s[_sign(xmult(p[(i+1)%n],p[(i+2)%n],p[i]))]=0;
    return s[0]&&s[1]|s[2];
}
```

### 判定点是否在多边形内

```c++
//判点在凸多边形内或多边形边上时返回 1，严格在凸多边形外返回0
int inside_convex(point q,int n,point* p) {
    int i,s[3]={1,1,1};
    for (i=0;i<n&&s[1]|s[2];i++)
        s[_sign(xmult(p[(i+1)%n],q,p[i]))]=0;
    return s[1]|s[2];
}

//判点严格在凸多边形内返回 1,在边上或者严格在外返回0
int inside_convex_v2(point q,int n,point* p) {
    int i,s[3]={1,1,1};
    for (i=0;i<n&&s[0]&&s[1]|s[2];i++)
        s[_sign(xmult(p[(i+1)%n],q,p[i]))]=0;
    return s[0]&&s[1]|s[2];
}

//判点在任意多边形内,顶点按顺时针或逆时针给出
//on_edge 表示点在多边形边上时的返回值
//offset为多边形坐标上限,严格在内返回1，严格在外返回0
int inside_polygon(point q,int n,point* p,int on_edge=2)
{
    point q2;
    int i=0,count;
    while (i<n)
        for (count=i=0,q2.x=rand()+offset,q2.y=rand()+offset;i<n;i++)
        {
            if(zero(xmult(q,p[i],p[(i+1)%n]))
                &&(p[i].x-q.x)*(p[(i+1)%n].x-q.x)<eps
                &&(p[i].y-q.y)*(p[(i+1)%n].y-q.y)<eps)
                    return on_edge;
            else if (zero(xmult(q,q2,p[i])))
                break;
            else if (xmult(q,p[i],q2)*xmult(q,p[(i+1)%n],q2)<-eps&&
                     xmult(p[i],q,p[(i+1)%n])*xmult(p[i],q2,p[(i+1)%n])<-eps)
                count++;
        }
    return count&1;
}
```

### 判定一条线段是否在一个任意多边形内

```c++
//预备函数
inline int opposite_side(point p1,point p2,point l1,point l2)
{
    return xmult(l1,p1,l2)*xmult(l1,p2,l2)<-eps;
}

inline int dot_online_in(point p,point l1,point l2)
{
    return zero(xmult(p,l1,l2))
            &&(l1.x-p.x)*(l2.x-p.x)<eps
            &&(l1.y-p.y)*(l2.y-p.y)<eps;
}

//判线段在任意多边形内,顶点按顺时针或逆时针给出,与边界相交返回 1
int inside_polygon(point l1,point l2,int n,point* p) {
    point t[MAXN],tt;
    int i,j,k=0;
    if (!inside_polygon(l1,n,p)||!inside_polygon(l2,n,p))
        return 0;
    for (i=0;i<n;i++)
    {
        if(opposite_side(l1,l2,p[i],p[(i+1)%n])
            &&opposite_side(p[i],p[(i+1)%n],l1,l2))
            return 0;
        else if (dot_online_in(l1,p[i],p[(i+1)%n]))
            t[k++]=l1;
        else if (dot_online_in(l2,p[i],p[(i+1)%n]))
            t[k++]=l2;
        else if (dot_online_in(p[i],l1,l2))
            t[k++]=p[i];
    }
    for (i=0;i<k;i++)
        for (j=i+1;j<k;j++)
        {
            tt.x=(t[i].x+t[j].x)/2;
            tt.y=(t[i].y+t[j].y)/2;
            if (!inside_polygon(tt,n,p))
                return 0;
        }
    return 1;
}
```

## 三角形

### 预备函数

```c++
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include<stdio.h>

//定义点
struct point
{
    double x,y;
};

//定义直线
struct line
{
    point a,b;
};

//两点距离
double distance(point p1,point p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}

//两直线求交点
point intersection(line u,line v)
{
    point ret=u.a;
    double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
        /((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
    ret.x+=(u.b.x-u.a.x)*t;
    ret.y+=(u.b.y-u.a.y)*t;
    return ret;
}
```

### 求三角形的外心

```c++
point circumcenter(point a,point b,point c)
{
    line u,v;
    u.a.x=(a.x+b.x)/2;
    u.a.y=(a.y+b.y)/2;
    u.b.x=u.a.x-a.y+b.y;
    u.b.y=u.a.y+a.x-b.x;
    v.a.x=(a.x+c.x)/2;
    v.a.y=(a.y+c.y)/2;
    v.b.x=v.a.x-a.y+c.y;
    v.b.y=v.a.y+a.x-c.x;
    return intersection(u,v);
}
```

### 求三角形内心

```c++
point incenter(point a,point b,point c)
{
    line u,v;
    double m,n;
    u.a=a;
    m=atan2(b.y-a.y,b.x-a.x);
    n=atan2(c.y-a.y,c.x-a.x);
    u.b.x=u.a.x+cos((m+n)/2);
    u.b.y=u.a.y+sin((m+n)/2);
    v.a=b;
    m=atan2(a.y-b.y,a.x-b.x);
    n=atan2(c.y-b.y,c.x-b.x);
    v.b.x=v.a.x+cos((m+n)/2);
    v.b.y=v.a.y+sin((m+n)/2);
    return intersection(u,v);
}
```

### 求三角形垂心

```c++
point perpencenter(point a,point b,point c)
{
    line u,v;
    u.a=c;
    u.b.x=u.a.x-a.y+b.y;
    u.b.y=u.a.y+a.x-b.x;
    v.a=b;
    v.b.x=v.a.x-a.y+c.y;
    v.b.y=v.a.y+a.x-c.x;
    return intersection(u,v);
}
```

## 圆

### 预备函数

```c++
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define eps 1e-8

struct point
{
    double x,y;
};

double xmult(point p1,point p2,point p0)
{
    return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}

double distance(point p1,point p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}

//点到直线的距离
double disptoline(point p,point l1,point l2)
{
    return fabs(xmult(p,l1,l2))/distance(l1,l2);
}

//求两直线交点
point intersection(point u1,point u2,point v1,point v2)
{
    point ret=u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
        /((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x+=(u2.x-u1.x)*t;
    ret.y+=(u2.y-u1.y)*t;
    return ret;
}
```

### 判定直线是否与圆相交

```c++
//判直线和圆相交,包括相切
int intersect_line_circle(point c,double r,point l1,point l2)
{
    return disptoline(c,l1,l2)<r+eps;
}
```

### 判定线段与圆相交

```c++
int intersect_seg_circle(point c,double r, point l1,point l2)
{
    double t1=distance(c,l1)-r,t2=distance(c,l2)-r;
    point t=c;
    if (t1<eps||t2<eps)
        return t1>-eps||t2>-eps;
    t.x+=l1.y-l2.y;
    t.y+=l2.x-l1.x;
    return xmult(l1,c,t)*xmult(l2,c,t)<eps&&disptoline(c,l1,l2)-r<eps;
}
```

### 判圆和圆相交

```c++
int intersect_circle_circle(point c1,double r1,point c2,double r2)
{
    return distance(c1,c2)<r1+r2+eps&&distance(c1,c2)>fabs(r1-r2)-eps;
}

```

### 计算圆上到点 p 最近点

```c++
//当 p 为圆心时，返回圆心本身
point dot_to_circle(point c,double r,point p)
{
    point u,v;
    if (distance(p,c)<eps)
        return p;
    u.x=c.x+r*fabs(c.x-p.x)/distance(c,p);
    u.y=c.y+r*fabs(c.y-p.y)/distance(c,p)*((c.x-p.x)*(c.y-p.y)<0?-1:1);
    v.x=c.x-r*fabs(c.x-p.x)/distance(c,p);
    v.y=c.y-r*fabs(c.y-p.y)/distance(c,p)*((c.x-p.x)*(c.y-p.y)<0?-1:1);
    return distance(u,p)<distance(v,p)?u:v;
}
```

### 计算直线与圆的交点

```c++
//计算直线与圆的交点,保证直线与圆有交点
//计算线段与圆的交点可用这个函数后判点是否在线段上
void intersection_line_circle(point c,double r,point l1
                              ,point l2,point& p1,point& p2)
{
    point p=c;
    double t;
    p.x+=l1.y-l2.y;
    p.y+=l2.x-l1.x;
    p=intersection(p,c,l1,l2);
    t=sqrt(r*r-distance(p,c)*distance(p,c))/distance(l1,l2);
    p1.x=p.x+(l2.x-l1.x)*t;
    p1.y=p.y+(l2.y-l1.y)*t;
    p2.x=p.x-(l2.x-l1.x)*t;
    p2.y=p.y-(l2.y-l1.y)*t;
}
```

### 计算两个圆的交点

```c++
//计算圆与圆的交点,保证圆与圆有交点,圆心不重合
void intersection_circle_circle(point c1,double r1,point c2
                                ,double r2,point& p1,point& p2)
{
    point u,v;
    double t;
    t=(1+(r1*r1-r2*r2)/distance(c1,c2)/distance(c1,c2))/2;
    u.x=c1.x+(c2.x-c1.x)*t;
    u.y=c1.y+(c2.y-c1.y)*t;
    v.x=u.x+c1.y-c2.y;
    v.y=u.y-c1.x+c2.x;
    intersection_line_circle(c1,r1,u,v,p1,p2);
}
```

## 球面

### 给出地球经度纬度，计算圆心角

```c++
#include <math.h>

const double pi=acos(-1);

//计算圆心角 lat 表示纬度,-90<=w<=90,lng 表示经度
//返回两点所在大圆劣弧对应圆心角,0<=angle<=pi
double angle(double lng1,double lat1,double lng2,double lat2)
{
    double dlng=fabs(lng1-lng2)*pi/180;
    while (dlng>=pi+pi)
        dlng-=pi+pi;
    if (dlng>pi)
        dlng=pi+pi-dlng;
    lat1*=pi/180,lat2*=pi/180;
    return acos(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2));
}
```

### 已知经纬度，计算地球上两点直线距离

```c++
//计算距离,r 为球半径
double line_dist(double r,double lng1,double lat1,double lng2,double lat2)
{
    double dlng=fabs(lng1-lng2)*pi/180;
    while (dlng>=pi+pi)
        dlng-=pi+pi;
    if (dlng>pi)
        dlng=pi+pi-dlng;
    lat1*=pi/180,lat2*=pi/180;
    return
        r*sqrt(2-2*(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2)));
}
```

### 已知经纬度，计算地球上两点球面距离

```c++
//计算球面距离,r 为球半径
inline double sphere_dist(double r,double lng1
                          ,double lat1,double lng2,double lat2)
{
    return r*angle(lng1,lat1,lng2,lat2);
}
```

## 三维几何的若干模板

### 预备函数

```c++
//三维几何函数库
#include <math.h>

#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)

struct point3{double x,y,z;};
struct line3{point3 a,b;};
struct plane3{point3 a,b,c;};

//计算 cross product U x V
point3 xmult(point3 u,point3 v){
    point3 ret;
    ret.x=u.y*v.z-v.y*u.z;
    ret.y=u.z*v.x-u.x*v.z;
    ret.z=u.x*v.y-u.y*v.x;
    return ret;
}

//计算 dot product U . V
double dmult(point3 u,point3 v){
    return u.x*v.x+u.y*v.y+u.z*v.z;
}

//矢量差 U - V
point3 subt(point3 u,point3 v){
    point3 ret;
    ret.x=u.x-v.x;
    ret.y=u.y-v.y;
    ret.z=u.z-v.z;
    return ret;
}

//取平面法向量
point3 pvec(plane3 s){
    return xmult(subt(s.a,s.b),subt(s.b,s.c));
}

point3 pvec(point3 s1,point3 s2,point3 s3){
    return xmult(subt(s1,s2),subt(s2,s3));
}

//两点距离,单参数取向量大小
double distance(point3 p1,point3 p2){
    return
        sqrt((p1.x-p2.x)*(p1.x-p2.x)+
             (p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
}

//向量大小
double vlen(point3 p){
    return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}
```

### 判定三点是否共线

```c++
//判三点共线
int dots_inline(point3 p1,point3 p2,point3 p3)
{
    return vlen(xmult(subt(p1,p2),subt(p2,p3)))<eps;
}
```

### 判定四点是否共面

```c++
//判四点共面
int dots_onplane(point3 a,point3 b,point3 c,point3 d)
{
    return zero(dmult(pvec(a,b,c),subt(d,a)));
}
```

### 判定点是否在线段上

```c++

//判点是否在线段上,包括端点和共线
int dot_online_in(point3 p,line3 l){
    return zero(vlen(xmult(subt(p,l.a),subt(p,l.b))))
            &&(l.a.x-p.x)*(l.b.x-p.x)<eps
            &&(l.a.y-p.y)*(l.b.y-p.y)<eps
            &&(l.a.z-p.z)*(l.b.z-p.z)<eps;
}

int dot_online_in(point3 p,point3 l1,point3 l2){
    return zero(vlen(xmult(subt(p,l1),subt(p,l2))))
            &&(l1.x-p.x)*(l2.x-p.x)<eps
            &&(l1.y-p.y)*(l2.y-p.y)<eps
            &&(l1.z-p.z)*(l2.z-p.z)<eps;
}

//判点是否在线段上,不包括端点
int dot_online_ex(point3 p,line3 l){
    return
        dot_online_in(p,l)
        &&(!zero(p.x-l.a.x)||!zero(p.y-l.a.y)||!zero(p.z-l.a.z))
        &&(!zero(p.x-l.b.x)||!zero(p.y-l.b.y)||!zero(p.z-l.b.z));
}

int dot_online_ex(point3 p,point3 l1,point3 l2){
    return
        dot_online_in(p,l1,l2)
        &&(!zero(p.x-l1.x)||!zero(p.y-l1.y)||!zero(p.z-l1.z))
        &&(!zero(p.x-l2.x)||!zero(p.y-l2.y)||!zero(p.z-l2.z));
}
```

### 判断点是否在空间三角形上

```c++
//判点是否在空间三角形上,包括边界,三点共线无意义
int dot_inplane_in(point3 p,plane3 s){
    return zero(vlen(xmult(subt(s.a,s.b),subt(s.a,s.c)))
             -vlen(xmult(subt(p,s.a),subt(p,s.b)))
             -vlen(xmult(subt(p,s.b),subt(p,s.c)))
             -vlen(xmult(subt(p,s.c),subt(p,s.a))));
}
int dot_inplane_in(point3 p,point3 s1,point3 s2,point3 s3){
    return zero(vlen(xmult(subt(s1,s2),subt(s1,s3)))
             -vlen(xmult(subt(p,s1),subt(p,s2)))
             -vlen(xmult(subt(p,s2),subt(p,s3)))
             -vlen(xmult(subt(p,s3),subt(p,s1))));
}
//判点是否在空间三角形上,不包括边界,三点共线无意义
int dot_inplane_ex(point3 p,plane3 s){
    return dot_inplane_in(p,s)
        &&vlen(xmult(subt(p,s.a),subt(p,s.b)))>eps
        &&vlen(xmult(subt(p,s.b),subt(p,s.c)))>eps
        &&vlen(xmult(subt(p,s.c),subt(p,s.a)))>eps;
}
int dot_inplane_ex(point3 p,point3 s1,point3 s2,point3 s3){
    return dot_inplane_in(p,s1,s2,s3)
        &&vlen(xmult(subt(p,s1),subt(p,s2)))>eps
        &&vlen(xmult(subt(p,s2),subt(p,s3)))>eps
        &&vlen(xmult(subt(p,s3),subt(p,s1)))>eps;
}
```

### 判断两点是否在线段同侧

```c++
int same_side(point3 p1,point3 p2,line3 l){
    return dmult(xmult(subt(l.a,l.b),subt(p1,l.b))
                ,xmult(subt(l.a,l.b),subt(p2,l.b)))>eps;
}

int same_side(point3 p1,point3 p2,point3 l1,point3 l2){
    return dmult(xmult(subt(l1,l2),subt(p1,l2))
                ,xmult(subt(l1,l2),subt(p2,l2)))>eps;
}
```

### 判断两点是否在线段异侧

```c++
//判两点在线段异侧,点在线段上返回 0,不共面无意义
int opposite_side(point3 p1,point3 p2,line3 l){
    return dmult(xmult(subt(l.a,l.b),subt(p1,l.b))
                ,xmult(subt(l.a,l.b),subt(p2,l.b)))<-eps;
}

int opposite_side(point3 p1,point3 p2,point3 l1,point3 l2){
    return dmult(xmult(subt(l1,l2),subt(p1,l2))
                ,xmult(subt(l1,l2),subt(p2,l2)))<-eps;
}
```

### 判断两点是否在平面同侧

```c++
//判两点在平面同侧,点在平面上返回 0
int same_side(point3 p1,point3 p2,plane3 s){
    return dmult(pvec(s),subt(p1,s.a))*dmult(pvec(s),subt(p2,s.a))>eps;
}

int same_side(point3 p1,point3 p2,point3 s1,point3 s2,point3 s3){
    return dmult(pvec(s1,s2,s3),subt(p1,s1))
            *dmult(pvec(s1,s2,s3),subt(p2,s1))>eps;
}
```

### 判断两点是否在平面异侧

```c++
//判两点在平面异侧,点在平面上返回 0
int opposite_side(point3 p1,point3 p2,plane3 s){
    return dmult(pvec(s),subt(p1,s.a))*dmult(pvec(s),subt(p2,s.a))<-eps;
}

int opposite_side(point3 p1,point3 p2,point3 s1,point3 s2,point3 s3){
    return dmult(pvec(s1,s2,s3),subt(p1,s1))
            *dmult(pvec(s1,s2,s3),subt(p2,s1))<-eps;
}
```

### 判断两空间直线是否平行

```c++
//判两直线平行
int parallel(line3 u,line3 v){
    return vlen(xmult(subt(u.a,u.b),subt(v.a,v.b)))<eps;
}

int parallel(point3 u1,point3 u2,point3 v1,point3 v2){
    return vlen(xmult(subt(u1,u2),subt(v1,v2)))<eps;
}
```

### 判断两平面是否平行

```c++
//判两平面平行
int parallel(plane3 u,plane3 v){
    return vlen(xmult(pvec(u),pvec(v)))<eps;
}

int parallel(point3 u1,point3 u2,point3 u3
             ,point3 v1,point3 v2,point3 v3)
{ 
    return vlen(xmult(pvec(u1,u2,u3),pvec(v1,v2,v3)))<eps;
}
```

### 判断直线是否与平面平行

```c++
//判直线与平面平行
int parallel(line3 l,plane3 s){
    return zero(dmult(subt(l.a,l.b),pvec(s)));
}

int parallel(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
    return zero(dmult(subt(l1,l2),pvec(s1,s2,s3)));
}
```

### 判断两直线是否垂直

```c++
//判两直线垂直
int perpendicular(line3 u,line3 v){
    return zero(dmult(subt(u.a,u.b),subt(v.a,v.b)));
}

int perpendicular(point3 u1,point3 u2,point3 v1,point3 v2){
    return zero(dmult(subt(u1,u2),subt(v1,v2)));
}
```

### 判断两平面是否垂直

```c++
//判两平面垂直
int perpendicular(plane3 u,plane3 v){
    return zero(dmult(pvec(u),pvec(v)));
}

int perpendicular(point3 u1,point3 u2,point3 u3
                  ,point3 v1,point3 v2,point3 v3)
{
    return zero(dmult(pvec(u1,u2,u3),pvec(v1,v2,v3)));
}
```

### 判断两条空间线段是否相交

```c++
//判两线段相交,包括端点和部分重合
int intersect_in(line3 u,line3 v){
    if (!dots_onplane(u.a,u.b,v.a,v.b))
        return 0;
    if (!dots_inline(u.a,u.b,v.a)||!dots_inline(u.a,u.b,v.b)) 
        return !same_side(u.a,u.b,v)&&!same_side(v.a,v.b,u);
    return dot_online_in(u.a,v)
        ||dot_online_in(u.b,v)
        ||dot_online_in(v.a,u)
        ||dot_online_in(v.b,u);
}

int intersect_in(point3 u1,point3 u2,point3 v1,point3 v2){ 
    if(!dots_onplane(u1,u2,v1,v2))
        return 0;
    if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
        return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2); 
    return dot_online_in(u1,v1,v2)
        ||dot_online_in(u2,v1,v2)
        ||dot_online_in(v1,u1,u2)
        ||dot_online_in(v2,u1,u2);
}

//判两线段相交,不包括端点和部分重合
int intersect_ex(line3 u,line3 v){
    return dots_onplane(u.a,u.b,v.a,v.b)
        &&opposite_side(u.a,u.b,v)
        &&opposite_side(v.a,v.b,u);
}

int intersect_ex(point3 u1,point3 u2,point3 v1,point3 v2){
    return dots_onplane(u1,u2,v1,v2)
        &&opposite_side(u1,u2,v1,v2)
        &&opposite_side(v1,v2,u1,u2);
}
```

### 判断线段是否与空间三角形相交

```c++
//判线段与空间三角形相交,包括交于边界和(部分)包含
int intersect_in(line3 l,plane3 s){
    return !same_side(l.a,l.b,s)&&!same_side(s.a,s.b,l.a,l.b,s.c)&&
        !same_side(s.b,s.c,l.a,l.b,s.a)&&!same_side(s.c,s.a,l.a,l.b,s.b);
}

int intersect_in(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
    return !same_side(l1,l2,s1,s2,s3)&&!same_side(s1,s2,l1,l2,s3)&&
        !same_side(s2,s3,l1,l2,s1)&&!same_side(s3,s1,l1,l2,s2);
}

//判线段与空间三角形相交,不包括交于边界和(部分)包含
int intersect_ex(line3 l,plane3 s){
    return
        opposite_side(l.a,l.b,s)&&opposite_side(s.a,s.b,l.a,l.b,s.c)&&
        opposite_side(s.b,s.c,l.a,l.b,s.a)&&opposite_side(s.c,s.a,l.a,l.b,s.b);
}

int intersect_ex(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
    return
        opposite_side(l1,l2,s1,s2,s3)&&opposite_side(s1,s2,l1,l2,s3)&&
        opposite_side(s2,s3,l1,l2,s1)&&opposite_side(s3,s1,l1,l2,s2);
}
```

### 计算两条直线的交点

```c++
//计算两直线交点,注意事先判断直线是否共面和平行 !
//线段交点请另外判线段相交(同时还是要判断是否平行!)
point3 intersection(line3 u,line3 v){
    point3 ret=u.a;
    double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
        /((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
    ret.x+=(u.b.x-u.a.x)*t;
    ret.y+=(u.b.y-u.a.y)*t;
    ret.z+=(u.b.z-u.a.z)*t;
    return ret;
}

point3 intersection(point3 u1,point3 u2,point3 v1,point3 v2){
    point3 ret=u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
        /((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x+=(u2.x-u1.x)*t;
    ret.y+=(u2.y-u1.y)*t;
    ret.z+=(u2.z-u1.z)*t;
    return ret;
}
```

### 计算直线与平面的交点

```c++
//计算直线与平面交点,注意事先判断是否平行,并保证三点不共线!
//线段和空间三角形交点请另外判断
point3 intersection(line3 l,plane3 s){
    point3 ret=pvec(s);
    double t=(ret.x*(s.a.x-l.a.x)+ret.y*(s.a.y-l.a.y)+ret.z*(s.a.z-l.a.z))/
        (ret.x*(l.b.x-l.a.x)+ret.y*(l.b.y-l.a.y)+ret.z*(l.b.z-l.a.z));
    ret.x=l.a.x+(l.b.x-l.a.x)*t;
    ret.y=l.a.y+(l.b.y-l.a.y)*t;
    ret.z=l.a.z+(l.b.z-l.a.z)*t;
    return ret;
}

point3 intersection(point3 l1,point3 l2
                    ,point3 s1,point3 s2,point3 s3)
{
    point3 ret=pvec(s1,s2,s3);
    double t=(ret.x*(s1.x-l1.x)+ret.y*(s1.y-l1.y)+ret.z*(s1.z-l1.z))/
        (ret.x*(l2.x-l1.x)+ret.y*(l2.y-l1.y)+ret.z*(l2.z-l1.z));
    ret.x=l1.x+(l2.x-l1.x)*t;
    ret.y=l1.y+(l2.y-l1.y)*t;
    ret.z=l1.z+(l2.z-l1.z)*t;
    return ret;
}
```

### 计算两平面的交线

```c++
//计算两平面交线,注意事先判断是否平行,并保证三点不共线!
line3 intersection(plane3 u,plane3 v)
{
    line3 ret;
    ret.a=parallel(v.a,v.b,u.a,u.b,u.c)
        ?intersection(v.b,v.c,u.a,u.b,u.c)
        :intersection(v.a,v.b,u.a,u.b,u.c);
    ret.b=parallel(v.c,v.a,u.a,u.b,u.c)
        ?intersection(v.b,v.c,u.a,u.b,u.c)
        :intersection(v.c,v.a,u.a,u.b,u.c);
    return ret;
}

line3 intersection(point3 u1,point3 u2,point3 u3
                   ,point3 v1,point3 v2,point3 v3)
{
    line3 ret;
    ret.a=parallel(v1,v2,u1,u2,u3)
        ?intersection(v2,v3,u1,u2,u3)
        :intersection(v1,v2,u1,u2,u3);
    ret.b=parallel(v3,v1,u1,u2,u3)
        ?intersection(v2,v3,u1,u2,u3)
        :intersection(v3,v1,u1,u2,u3);
    return ret;
}
```

### 点到直线的距离

```c++
//点到直线距离
double ptoline(point3 p,line3 l){
    return vlen(xmult(subt(p,l.a),subt(l.b,l.a)))/distance(l.a,l.b);
}

double ptoline(point3 p,point3 l1,point3 l2){
    return vlen(xmult(subt(p,l1),subt(l2,l1)))/distance(l1,l2);
}
```

### 计算点到平面的距离

```c++
//点到平面距离
double ptoplane(point3 p,plane3 s){
    return fabs(dmult(pvec(s),subt(p,s.a)))/vlen(pvec(s));
}

double ptoplane(point3 p,point3 s1,point3 s2,point3 s3){
    return fabs(dmult(pvec(s1,s2,s3),subt(p,s1)))/vlen(pvec(s1,s2,s3));
}
```

### 计算直线到直线的距离

```c++
//直线到直线距离
double linetoline(line3 u,line3 v){
    point3 n=xmult(subt(u.a,u.b),subt(v.a,v.b));
    return fabs(dmult(subt(u.a,v.a),n))/vlen(n);
}

double linetoline(point3 u1,point3 u2,point3 v1,point3 v2){ 
    point3 n=xmult(subt(u1,u2),subt(v1,v2));
    return fabs(dmult(subt(u1,v1),n))/vlen(n);
}
```

### 空间两直线夹角的 cos 值

```c++
//两直线夹角 cos 值
double angle_cos(line3 u,line3 v){
    return dmult(subt(u.a,u.b),subt(v.a,v.b))
            /vlen(subt(u.a,u.b))
             /vlen(subt(v.a,v.b));
}

double angle_cos(point3 u1,point3 u2,point3 v1,point3 v2){
    return dmult(subt(u1,u2),subt(v1,v2))
            /vlen(subt(u1,u2))
            /vlen(subt(v1,v2));
}
```

### 两平面夹角的 cos 值

```c++
//两平面夹角 cos 值
double angle_cos(plane3 u,plane3 v){
    return dmult(pvec(u),pvec(v))/vlen(pvec(u))/vlen(pvec(v));
}

double angle_cos(point3 u1,point3 u2,point3 u3
                 ,point3 v1,point3 v2,point3 v3)
{
    return dmult(pvec(u1,u2,u3),pvec(v1,v2,v3))
        /vlen(pvec(u1,u2,u3))
        /vlen(pvec(v1,v2,v3));
}
```

### 直线与平面夹角 sin 值

```c++
//直线平面夹角 sin 值
double angle_sin(line3 l,plane3 s){
    return dmult(subt(l.a,l.b),pvec(s))/vlen(subt(l.a,l.b))/vlen(pvec(s));
}

double angle_sin(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
    return dmult(subt(l1,l2),pvec(s1,s2,s3))
            /vlen(subt(l1,l2))
            /vlen(pvec(s1,s2,s3));
}
```

## 最远曼哈顿距离

> 计算5维坐标系中曼哈顿距离最远的两个点的

```c++
#include <stdio.h>

#define INF 9999999999999.0

struct Point{
    double x[5];
}pt[100005];

double dis[32][100005], coe[5], minx[32], maxx[32];

//去掉绝对值后有 2^D 种可能
void GetD(int N, int D)
{
    int s, i, j, tot=(1<<D);
    for (s=0;s<tot;s++)
    {
        for (i=0;i<D;i++)
            if (s&(1<<i))
                coe[i]=-1.0;
            else coe[i]=1.0;
        for (i=0;i<N;i++)
        {
            dis[s][i]=0.0;
            for (j=0;j<D;j++)
                dis[s][i]=dis[s][i]+coe[j]*pt[i].x[j];
        }
    }
}

//取每种可能中的最大差距
void Solve(int N, int D)
{
    int s, i, tot=(1<<D);
    double tmp, ans;
    for (s=0;s<tot;s++)
    {
        minx[s]=INF;
        maxx[s]=-INF;
        for (i=0; i<N; i++)
        {
            if (minx[s]>dis[s][i]) minx[s]=dis[s][i];
            if (maxx[s]<dis[s][i]) maxx[s]=dis[s][i];
        }
    }
    ans=0.0;
    for (s=0; s<tot; s++)
    {
        tmp=maxx[s]-minx[s];
        if (tmp>ans) ans=tmp;
    }
    printf("%.2lf\n", ans);
}

int main (void)
{
    int n, i;
    while (scanf("%d",&n)==1)
    {
        for (i=0;i<n;i++)
            scanf("%lf%lf%lf%lf%lf",&pt[i].x[0]
                  ,&pt[i].x[1],&pt[i].x[2]
                  ,&pt[i].x[3],&pt[i].x[4]);
        GetD(n, 5);
        Solve(n, 5);
    }
    return 0;
}
```

## 最近点对

```c++
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define Max(x,y) (x)>(y)?(x):(y)

struct Q{
    double x, y;
}q[100001], sl[10], sr[10];

int cntl, cntr, lm, rm;
double ans;

int cmp(const void*p1, const void*p2)
{
    struct Q*a1=(struct Q*)p1;
    struct Q*a2=(struct Q*)p2;
    if (a1->x<a2->x)return -1;
    else if (a1->x==a2->x)return 0;
    else return 1;
}

double CalDis(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

void MinDis(int l, int r)
{
    if (l==r) return;
    double dis;
    if (l+1==r)
    {
        dis=CalDis(q[l].x,q[l].y,q[r].x,q[r].y);
        if (ans>dis) ans=dis;
        return;
    }
    int mid=(l+r)>>1, i, j;
    MinDis(l,mid);
    MinDis(mid+1,r);
    lm=mid+1-5;
    if (lm<l) lm=l;
    rm=mid+5;
    if (rm>r) rm=r;
    cntl=cntr=0;
    for (i=mid;i>=lm;i--)
    {
        if (q[mid+1].x-q[i].x>=ans)break;
        sl[++cntl]=q[i];
    }
    for (i=mid+1;i<=rm;i++)
    {
        if (q[i].x-q[mid].x>=ans)break;
        sr[++cntr]=q[i];
    }
    for (i=1;i<=cntl;i++)
        for (j=1;j<=cntr;j++)
        {
            dis=CalDis(sl[i].x,sl[i].y,sr[j].x,sr[j].y);
            if (dis<ans) ans=dis;
        }
}

int main (void)
{
    int n, i;
    while (scanf("%d",&n)==1&&n)
    {
        for (i=1;i<=n;i++)
            scanf("%lf %lf", &q[i].x,&q[i].y);
        qsort(q+1,n,sizeof(struct Q),cmp);
        ans=CalDis(q[1].x,q[1].y,q[2].x,q[2].y);
        MinDis(1,n);
        printf("%.2lf\n",ans/2.0);
    }
    return 0;
}
```

## 最小包围圆

```c++
#include<stdio.h>
#include<string.h>
#include<math.h>

struct Point{
    double x;
    double y;
}pt[1005];

struct Traingle{
    struct Point p[3];
};

struct Circle{
    struct Point center;
    double r;
}ans;

//计算两点距离
double Dis(struct Point p, struct Point q)
{
    double dx=p.x-q.x;
    double dy=p.y-q.y;
    return sqrt(dx*dx+dy*dy);
}

//计算三角形面积
double Area(struct Traingle ct)
{
    return fabs((ct.p[1].x-ct.p[0].x)*(ct.p[2].y-ct.p[0].y)
                -(ct.p[2].x-ct.p[0].x)*(ct.p[1].y-ct.p[0].y))/2.0;
}

//求三角形的外接圆，返回圆心和半径 (存在结构体"圆"中)
struct Circle CircumCircle(struct Traingle t)
{
    struct Circle tmp;
    double a, b, c, c1, c2;
    double xA, yA, xB, yB, xC, yC;
    a = Dis(t.p[0], t.p[1]);
    b = Dis(t.p[1], t.p[2]);
    c = Dis(t.p[2], t.p[0]);
    //根据 S = a * b * c / R / 4;求半径 R
    tmp.r = (a*b*c)/(Area(t)*4.0);
    xA = t.p[0].x;
    yA = t.p[0].y;
    xB = t.p[1].x;
    yB = t.p[1].y;
    xC = t.p[2].x;
    yC = t.p[2].y;
    c1 = (xA*xA+yA*yA - xB*xB-yB*yB) / 2;
    c2 = (xA*xA+yA*yA - xC*xC-yC*yC) / 2;
    tmp.center.x = (c1*(yA - yC)-c2*(yA - yB))
        / ((xA - xB)*(yA - yC)-(xA - xC)*(yA - yB));
    tmp.center.y = (c1*(xA - xC)-c2*(xA - xB))
        / ((yA - yB)*(xA - xC)-(yA - yC)*(xA - xB));
    return tmp;
}

//确定最小包围圆
struct Circle MinCircle(int num, struct Traingle ct)
{
    struct Circle ret;
    if (num==0) ret.r = 0.0;
    else if (num==1)
    {
        ret.center = ct.p[0];
        ret.r = 0.0;
    }
    else if (num==2)
    {
        ret.center.x = (ct.p[0].x+ct.p[1].x)/2.0;
        ret.center.y = (ct.p[0].y+ct.p[1].y)/2.0;
        ret.r = Dis(ct.p[0], ct.p[1])/2.0;
    }
    else if(num==3) ret = CircumCircle(ct);
    return ret;
}

//递归实现增量算法
void Dfs(int x, int num, struct Traingle ct)
{
    int i, j;
    struct Point tmp;
    ans = MinCircle(num, ct);
    if (num==3) return;
    for (i=1; i<=x; i++)
        if (Dis(pt[i], ans.center)>ans.r)
        {
            ct.p[num]=pt[i];
            Dfs(i-1, num+1, ct);
            tmp=pt[i];
            for (j=i;j>=2;j--)
                pt[j]=pt[j-1];
            pt[1]=tmp;
        }
}

void Solve(int n)
{
    struct Traingle ct;
    Dfs(n, 0, ct);
}

int main (void)
{
    int n, i;
    while (scanf("%d", &n)!=EOF && n)
    {
        for (i=1;i<=n;i++)
            scanf("%lf %lf", &pt[i].x, &pt[i].y);
        Solve(n);
        printf("%.2lf %.2lf %.2lf\n", ans.center.x, ans.center.y, ans.r);
    }
    return 0;
}
```

## 求两个圆的交点

```c++
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>

const double eps = 1e-8;
const double PI = acos(-1.0);

struct Point{
    double x;
    double y;
};

struct Line{
    double s, t;
};

struct Circle{
    Point center;
    double r;
    Line line[505];
    int cnt;
    bool covered;
}circle[105];

double distance(point p1, point p2)
{
    double dx = p1.x-p2.x;
    double dy = p1.y-p2.y;
    return sqrt(dx*dx + dy*dy);
}

point intersection(point u1,point u2, point v1,point v2)
{
    point ret = u1;
    double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
        / ((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
    ret.x += (u2.x-u1.x)*t;
    ret.y += (u2.y-u1.y)*t;
    return ret;
}

void intersection_line_circle(point c,double r
                              ,point l1,point l2
                              ,point& p1,point& p2)
{
    point p=c;
    double t;
    p.x+=l1.y-l2.y;
    p.y+=l2.x-l1.x;
    p=intersection(p,c,l1,l2);
    t=sqrt(r*r-distance(p,c)*distance(p,c))/distance(l1,l2);
    p1.x=p.x+(l2.x-l1.x)*t;
    p1.y=p.y+(l2.y-l1.y)*t;
    p2.x=p.x-(l2.x-l1.x)*t;
    p2.y=p.y-(l2.y-l1.y)*t;
}

//计算圆与圆的交点,保证圆与圆有交点,圆心不重合
void intersection_circle_circle(point c1,double r1
                                ,point c2,double r2
                                ,point& p1,point& p2)
{
    point u,v;
    double t;
    t=(1+(r1*r1-r2*r2)/distance(c1,c2)/distance(c1,c2))/2;
    u.x=c1.x+(c2.x-c1.x)*t;
    u.y=c1.y+(c2.y-c1.y)*t;
    v.x=u.x+c1.y-c2.y;
    v.y=u.y-c1.x+c2.x;
    intersection_line_circle(c1,r1,u,v,p1,p2);
}
```

## 求三角形外接圆圆心

```c++
struct Point{
    double x;
    double y;
}pt[1005];

struct Traingle{
    struct Point p[3];
};

struct Circle{
    struct Point center;
    double r;
}ans;

//计算两点距离
double Dis(struct Point p, struct Point q)
{
    double dx=p.x-q.x;
    double dy=p.y-q.y;
    return sqrt(dx*dx+dy*dy);
}

//计算三角形面积
double Area(struct Traingle ct)
{
    return fabs((ct.p[1].x-ct.p[0].x)*(ct.p[2].y-ct.p[0].y)
                -(ct.p[2].x-ct.p[0].x)*(ct.p[1].y-ct.p[0].y))/2.0;
}

//求三角形的外接圆，返回圆心和半径 (存在结构体"圆"中)
struct Circle CircumCircle(struct Traingle t)
{
    struct Circle tmp;
    double a, b, c, c1, c2;
    double xA, yA, xB, yB, xC, yC;
    a = Dis(t.p[0], t.p[1]);
    b = Dis(t.p[1], t.p[2]);
    c = Dis(t.p[2], t.p[0]);
    //根据 S = a * b * c / R / 4;求半径 R
    tmp.r = (a*b*c)/(Area(t)*4.0);
    xA = t.p[0].x;
    yA = t.p[0].y;
    xB = t.p[1].x;
    yB = t.p[1].y;
    xC = t.p[2].x;
    yC = t.p[2].y;
    c1 = (xA*xA+yA*yA - xB*xB-yB*yB) / 2;
    c2 = (xA*xA+yA*yA - xC*xC-yC*yC) / 2;
    tmp.center.x = (c1*(yA - yC)-c2*(yA - yB))
        / ((xA - xB)*(yA - yC)-(xA - xC)*(yA - yB));
    tmp.center.y = (c1*(xA - xC)-c2*(xA - xB))
        / ((yA - yB)*(xA - xC)-(yA - yC)*(xA - xB));
    return tmp;
}
```

## 求凸包

```c++
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define INF 999999999.9
#define PI acos(-1.0)

struct Point{
    double x, y, dis;
}pt[1005], stack[1005], p0;

int top, tot;

//计算几何距离
double Dis(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//极角比较， 返回-1: p0p1 在 p0p2 的右侧，返回 0:p0,p1,p2 共线
int Cmp_PolarAngel(struct Point p1, struct Point p2, struct Point pb)
{
    double delta=(p1.x-pb.x)*(p2.y-pb.y)-(p2.x-pb.x)*(p1.y-pb.y);
    if(delta<0.0) return 1;
    else if (delta==0.0) return 0;
    else return -1;
}

//判断向量 p2p3 是否对 p1p2 构成左旋
bool Is_LeftTurn(struct Point p3, struct Point p2, struct Point p1)
{
    int type=Cmp_PolarAngel(p3, p1, p2);
    if (type<0) return true;
    return false;
}

//先按极角排，再按距离由小到大排
int Cmp(const void*p1, const void*p2)
{
    struct Point*a1=(struct Point*)p1;
    struct Point*a2=(struct Point*)p2;
    int type=Cmp_PolarAngel(*a1, *a2, p0);
    if (type<0) return -1;
    else if (type==0)
    {
        if (a1->dis<a2->dis) return -1;
        else if (a1->dis==a2->dis) return 0;
        else return 1;
    }
    else return 1;
}

//求凸包
void Solve(int n)
{
    int i, k;
    p0.x=p0.y=INF;
    for (i=0;i<n;i++)
    {
        scanf("%lf %lf",&pt[i].x, &pt[i].y);
        if (pt[i].y < p0.y)
        {
            p0.y=pt[i].y;
            p0.x=pt[i].x;
            k=i;
        }
        else if (pt[i].y==p0.y)
        {
            if (pt[i].x<p0.x)
            {
                p0.x=pt[i].x;
                k=i;
            }
        }
    }
    pt[k]=pt[0];
    pt[0]=p0;
    for (i=1;i<n;i++)
        pt[i].dis=Dis(pt[i].x,pt[i].y, p0.x,p0.y);
    qsort(pt+1, n-1, sizeof(struct Point), Cmp);
    //去掉极角相同的点
    tot=1;
    for (i=2;i<n;i++)
        if (Cmp_PolarAngel(pt[i], pt[i-1], p0))
            pt[tot++]=pt[i-1];
    pt[tot++]=pt[n-1];
    //求凸包
    top=1;
    stack[0]=pt[0];
    stack[1]=pt[1];
    for (i=2;i<tot;i++)
    {
        while (top>=1
               &&
               Is_LeftTurn(pt[i], stack[top],stack[top-1])==false)
            top--;
        stack[++top]=pt[i];
    }
}

int main (void)
{
    int n;
    while (scanf("%d",&n)==2)
    {
        Solve(n);
    }
    return 0;
}
```

## 凸包卡壳旋转求出所有对踵点、最远点对

```c++
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define INF 999999999.9
#define PI acos(-1.0)

struct Point{
    double x, y, dis;
}pt[6005], stack[6005], p0;

int top, tot;

//计算几何距离
double Dis(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//极角比较， 返回-1: p0p1 在 p0p2 的右侧，返回 0:p0,p1,p2 共线
int Cmp_PolarAngel(struct Point p1, struct Point p2, struct Point pb)
{
    double delta=(p1.x-pb.x)*(p2.y-pb.y)-(p2.x-pb.x)*(p1.y-pb.y);
    if(delta<0.0) return 1;
    else if (delta==0.0) return 0;
    else return -1;
}

//判断向量 p2p3 是否对 p1p2 构成左旋
bool Is_LeftTurn(struct Point p3, struct Point p2, struct Point p1)
{
    int type=Cmp_PolarAngel(p3, p1, p2);
    if (type<0) return true;
    return false;
}

//先按极角排，再按距离由小到大排
int Cmp(const void*p1, const void*p2)
{
    struct Point*a1=(struct Point*)p1;
    struct Point*a2=(struct Point*)p2;
    int type=Cmp_PolarAngel(*a1, *a2, p0);
    if (type<0) return -1;
    else if (type==0)
    {
        if (a1->dis<a2->dis) return -1;
        else if (a1->dis==a2->dis) return 0;
        else return 1;
    }
    else return 1;
}

//求凸包
void Hull(int n)
{
    int i, k;
    p0.x=p0.y=INF;
    for (i=0;i<n;i++)
    {
        scanf("%lf %lf",&pt[i].x, &pt[i].y);
        if (pt[i].y < p0.y)
        {
            p0.y=pt[i].y;
            p0.x=pt[i].x;
            k=i;
        }
        else if (pt[i].y==p0.y)
        {
            if (pt[i].x<p0.x)
            {
                p0.x=pt[i].x;
                k=i;
            }
        }
    }
    pt[k]=pt[0];
    pt[0]=p0;
    for (i=1;i<n;i++)
        pt[i].dis=Dis(pt[i].x,pt[i].y, p0.x,p0.y);
    qsort(pt+1, n-1, sizeof(struct Point), Cmp);
    //去掉极角相同的点
    tot=1;
    for (i=2;i<n;i++)
        if (Cmp_PolarAngel(pt[i], pt[i-1], p0))
            pt[tot++]=pt[i-1];
    pt[tot++]=pt[n-1];
    //求凸包
    top=1;
    stack[0]=pt[0];
    stack[1]=pt[1];
    for (i=2;i<tot;i++)
    {
        while (top>=1 && Is_LeftTurn(pt[i], stack[top],
                                     stack[top-1])==false)
            top--;
        stack[++top]=pt[i];
    }
}

//计算叉积
double CrossProduct(struct Point p1, struct Point p2, struct Point p3)
{
    return (p1.x-p3.x)*(p2.y-p3.y)-(p2.x-p3.x)*(p1.y-p3.y);
}

//卡壳旋转，求出凸多边形所有对踵点
void Rotate(struct Point*ch, int n)
{
    int i, p=1;
    double t1, t2, ans=0.0, dif;
    ch[n]=ch[0];
    for (i=0;i<n;i++)
    {
        //如果下一个点与当前边构成的三角形的面积更大，则说明此时不构成对踵点
        while (fabs(CrossProduct(ch[i],ch[i+1],ch[p+1])) >
               fabs(CrossProduct(ch[i],ch[i+1],ch[p])))
            p=(p+1)%n;
        dif=fabs(CrossProduct(ch[i],ch[i+1],ch[p+1])) -
            fabs(CrossProduct(ch[i],ch[i+1],ch[p]));
        //如果当前点和下一个点分别构成的三角形面积相等
        //则说明两条边即为平行线，对角线两端都可能是对踵点
        if (dif==0.0)
        {
            t1=Dis(ch[p].x, ch[p].y, ch[i].x, ch[i].y);
            t2=Dis(ch[p+1].x, ch[p+1].y, ch[i+1].x, ch[i+1].y);
            if (t1>ans)ans=t1;
            if (t2>ans)ans=t2;
        }
        //说明 p，i 是对踵点
        else if (dif<0.0)
        {
            t1=Dis(ch[p].x, ch[p].y, ch[i].x, ch[i].y);
            if (t1>ans)ans=t1;
        }
    }
    printf("%.2lf\n",ans);
}

int main (void)
{
    int n;
    while (scanf("%d",&n)==1)
    {
        Hull(n);
        Rotate(stack, top+1);
    }
    return 0;
}
```

## 凸包+旋转卡壳求平面面积最大三角

```c++
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define INF 99999999999.9
#define PI acos(-1.0)

struct Point{
    double x, y, dis;
}pt[50005], stack[50005], p0;

int top, tot;

double Dis(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

int Cmp_PolarAngel(struct Point p1, struct Point p2, struct Point pb)
{
    double delta=(p1.x-pb.x)*(p2.y-pb.y)-(p2.x-pb.x)*(p1.y-pb.y);
    if(delta<0.0) return 1;
    else if (delta==0.0) return 0;
    else return -1;
}

bool Is_LeftTurn(struct Point p3, struct Point p2, struct Point p1)
{
    int type=Cmp_PolarAngel(p3, p1, p2);
    if (type<0) return true;
    return false;
}

int Cmp(const void*p1, const void*p2)
{
    struct Point*a1=(struct Point*)p1;
    struct Point*a2=(struct Point*)p2;
    int type=Cmp_PolarAngel(*a1, *a2, p0);
    if (type<0) return -1;
    else if (type==0)
    {
        if (a1->dis<a2->dis) return -1;
        else if (a1->dis==a2->dis) return 0;
        else return 1;
    }
    else return 1;
}

void Hull(int n)
{
    int i, k;
    p0.x=p0.y=INF;
    for (i=0;i<n;i++)
    {
        scanf("%lf %lf",&pt[i].x, &pt[i].y);
        if (pt[i].y < p0.y)
        {
            p0.y=pt[i].y;
            p0.x=pt[i].x;
            k=i;
        }
        else if (pt[i].y==p0.y)
        {
            if (pt[i].x<p0.x)
            {
                p0.x=pt[i].x;
                k=i;
            }
        }
    }
    pt[k]=pt[0];
    pt[0]=p0;
    for (i=1;i<n;i++)
        pt[i].dis=Dis(pt[i].x,pt[i].y, p0.x,p0.y);
    qsort(pt+1, n-1, sizeof(struct Point), Cmp);
    tot=1;
    for (i=2;i<n;i++)
        if (Cmp_PolarAngel(pt[i], pt[i-1], p0))
            pt[tot++]=pt[i-1];
    pt[tot++]=pt[n-1];
    top=1;
    stack[0]=pt[0];
    stack[1]=pt[1];
    for (i=2;i<tot;i++)
    {
        while (top>=1 && Is_LeftTurn(pt[i], stack[top],
                                     stack[top-1])==false)
            top--;
        stack[++top]=pt[i];
    }
}

double TArea(struct Point p1, struct Point p2, struct Point p3)
{
    return fabs((p1.x-p3.x)*(p2.y-p3.y)-(p2.x-p3.x)*(p1.y-p3.y));
}

void Rotate(struct Point*ch, int n)
{
    if (n<3)
    {
        printf("0.00\n");
        return;
    }
    int i, j, k;
    double ans=0.0, tmp;
    ch[n]=ch[0];
    for (i=0;i<n;i++)
    {
        j=(i+1)%n;
        k=(j+1)%n;
        while ((j!=k) && (k!=i))
        {
            while
                (TArea(ch[i],ch[j],ch[k+1])>TArea(ch[i],ch[j],ch[k]))
                    k=(k+1)%n;
            tmp=TArea(ch[i],ch[j], ch[k]);
            if (tmp>ans) ans=tmp;
            j=(j+1)%n;
        }
    }
    printf("%.2lf\n",ans/2.0);
}

int main (void)
{
    int n;
    while (scanf("%d",&n)==1)
    {
        if (n==-1)break;
        Hull(n);
        Rotate(stack, top+1);
    }
    return 0;
}
```

## Pick 定理

**Pick 定理求整点多边形内部整点数目**

1. 给定顶点坐标均是整点（或正方形格点）的简单多边形，皮克定理说明了其面积$A$ 和内部格点数目 $i$、边上格点数目 $b$ 的关系：$A = i + \frac{b}{2} - 1$

2. 在两点$(x_1,y_1)$，$(x_2,y_2)$连线之间的整点个数（包含一个端点）为$gcd\left(\left|x_1－x_2\right|,\left|y_1－y_2\right|\right)$

3. 求三角形面积用叉乘

```c++
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

long long x[3], y[3], area, b;

long long My_Abs(long long t)
{
    if (t<0) return -t;
    return t;
}

long long Gcd(long long x, long long y)
{
    if (y==0) return x;
    long long mod=x%y;
    while (mod)
    {
        x=y;
        y=mod;
        mod=x%y;
    }
    return y;
}

int main (void)
{
    int i;
    while (1)
    {
        for (i = 0;i < 3;i ++)
            scanf("%lld %lld", &x[i], &y[i]);
        if(x[0]==0&&y[0]==0&&x[1]==0&&y[1]==0&&x[2]==0&&y[2]==0)
            break;
        area = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
        area = My_Abs(area);
        b=0;
        b=Gcd(My_Abs(x[1]-x[0]), My_Abs(y[1]-y[0])) +
            Gcd(My_Abs(x[2]-x[0]), My_Abs(y[2]-y[0])) +
            Gcd(My_Abs(x[1]-x[2]), My_Abs(y[1]-y[2]));
        printf("%lld\n", (area-b+2)/2);
    }
    return 0;
}
```

## 求多边形面积和重心

```c++
#include <stdio.h>
#include <math.h>

int x[1000003], y[1000003];
double A, tx, ty, tmp;

int main (void)
{
    int cases, n, i;
    scanf ("%d", &cases);
    while (cases --)
    {
        scanf ("%d", &n);
        A = 0.0;
        x[0] = y[0] = 0;
        for (i = 1; i <= n; i ++)
        {
            scanf ("%d %d", &x[i], &y[i]);
            A += (x[i-1]*y[i] - x[i]*y[i-1]);
        }
        A += x[n]*y[1] - x[1]*y[n];
        A = A / 2.0;
        tx = ty = 0.0;
        for (i = 1; i < n; i ++)
        {
            tmp = x[i]*y[i+1] - x[i+1]*y[i];
            tx += (x[i]+x[i+1]) * tmp;
            ty += (y[i]+y[i+1]) * tmp;
        }
        tmp = x[n]*y[1] - x[1]*y[n];
        tx += (x[n]+x[1])*tmp;
        ty += (y[n]+y[1])*tmp;
        printf ("%.2lf %.2lf\n", tx/(6.0*A), ty/(6.0*A));
    }
    return 0;
}
```

## 判断一个简单多边形是否有核

```c++
#include <stdio.h>
#include <string.h>

const int INF = (1<<30);

struct Point{
    int x, y;
}pt[150];

bool turn_right[150];

int det(Point s1, Point t1, Point s2, Point t2)
{
    int d1x = t1.x-s1.x;
    int d1y = t1.y-s1.y;
    int d2x = t2.x-s2.x;
    int d2y = t2.y-s2.y;
    return d1x*d2y - d2x*d1y;
}

void Swap(int &a, int &b)
{
    if (a>b)
    {
        int t=a;
        a=b;
        b=t;
    }
}

int main (void)
{
    int n, i, cross, maxx, minx, maxy, miny, maxn, minn, countn=0;
    while(scanf("%d", &n)==1&&n)
    {
        maxx=maxy=-INF;
        minx=miny=INF;
        //点按顺时针给出
        for (i=1; i<=n; i++)
        {
            scanf("%d %d", &pt[i].x, &pt[i].y);
            if (maxx<pt[i].x) maxx=pt[i].x;
            if (maxy<pt[i].y) maxy=pt[i].y;
            if (minx>pt[i].x) minx=pt[i].x;
            if (miny>pt[i].y) miny=pt[i].y;
        }
        pt[n+1]=pt[1];
        pt[n+2]=pt[2];
        pt[n+3]=pt[3];
        pt[n+4]=pt[4];
        //求每条线段的转向
        for (i=1; i<=n+1; i ++)
        {
            cross = det(pt[i],pt[i+1], pt[i+1], pt[i+2]);
            if (cross<0)turn_right[i+1]=true;
            else turn_right[i+1]=false;
        }
        for (i=2; i<= n+1; i++)
            if (turn_right[i] && turn_right[i+1])
            {
                if (pt[i].x==pt[i+1].x)
                {
                    minn=pt[i].y;
                    maxn=pt[i+1].y;
                    Swap(minn, maxn);
                    if (minn>miny) miny=minn;
                    if (maxn<maxy) maxy=maxn;
                }
                else
                {
                    minn=pt[i].x;
                    maxn=pt[i+1].x;
                    Swap(minn, maxn);
                    if (minn>minx) minx=minn;
                    if (maxn<maxx) maxx=maxn;
                }
            }
        if (minx<=maxx && miny<=maxy)
            printf("Floor #%d\nSurveillance is possible.\n\n", ++countn);
        else printf("Floor #%d\nSurveillance is impossible.\n\n", ++countn);
    }
    return 0;
}
```

## 模拟退火

```c++
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Lim 0.999999
#define EPS 1e-2
#define PI acos(-1.0)

double Temp, maxx, minx, maxy, miny, lx, ly, dif;
int nt, ns, nc;

struct Target{
    double x, y;
}T[105];

struct Solution{
    double x, y;
    double f;
}S[25], P, A;

double Dis(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

void Seed(void)
{
    int i, j;
    for (i=0;i<ns;i++)
    {
        S[i].x=minx+((double)(rand()%1000+1)/1000.0)*lx;
        S[i].y=miny+((double)(rand()%1000+1)/1000.0)*ly;
        S[i].f=0.0;
        for (j=0;j<nt;j++)
            S[i].f=S[i].f+Dis(S[i].x,S[i].y, T[j].x, T[j].y);
    }
}

void Trans(void)
{
    int i, j, k;
    double theta;
    for (i=0;i<ns;i++)
    {
        P=S[i];
        for (j=0;j<nc;j++)
        {
            theta=(((double)(rand()%1000+1))/1000.0)*2.0*PI;
            A.x=P.x+Temp*cos(theta);
            A.y=P.y+Temp*sin(theta);
            if (A.x<minx||A.x>maxx||A.y<miny||A.y>maxy)
                continue;
            A.f=0.0;
            for (k=0;k<nt;k++)
                A.f=A.f+Dis(A.x,A.y,T[k].x,T[k].y);
            dif=A.f-S[i].f;
            if (dif<0.0)S[i]=A;
            else
            {
                dif=exp(-dif/Temp);
                if (dif>Lim) S[i]=A;
            }
        }
    }
}

int main (void)
{
    int i, k;
    while (scanf("%d",&nt)==1&&nt)
    {
        maxx=maxy=0;
        minx=miny=(1<<20);
        for (i=0;i<nt;i++)
        {
            scanf("%lf %lf",&T[i].x,&T[i].y);
            if (maxx<T[i].x)maxx=T[i].x;
            if (minx>T[i].x)minx=T[i].x;
            if (maxy<T[i].y)maxy=T[i].y;
            if (miny>T[i].y)miny=T[i].y;
        }
        lx=maxx-minx;
        ly=maxy-miny;
        Temp=sqrt(lx*lx+ly*ly)/3.0;
        ns=5, nc=10;
        Seed();
        while (Temp>EPS)
        {
            Trans();
            Temp=Temp*0.40;
        }
        k=0;
        for (i=1;i<ns;i++)
            if (S[k].f>S[i].f)
                k=i;
        printf ("%.0lf\n", S[k].f);
    }
    return 0;
}
```

## 六边形坐标系

```c++
//第一种六边形坐标系
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

double Dis(double x1, double y1, double x2, double y2)
{
    double dx=x1-x2;
    double dy=y1-y2;
    return sqrt(dx*dx+dy*dy);
}

void Get_KL(double L, double x, double y, int &k, int &l, double &cd)
{
    k=floor((2.0*x)/(3.0*L));
    l=floor((2.0*y)/(sqrt(3.0)*L));
    double d1, d2, x1, y1, x2, y2;
    if ((k+l)&1)
    {
        x1=k*L*1.5;
        y1=(l+1.0)*L*sqrt(3.0)*0.5;
        x2=(k+1.0)*L*1.5;
        y2=l*L*sqrt(3.0)*0.5;
        d1=Dis(x1,y1, x,y);
        d2=Dis(x2,y2, x,y);
        if (d1>d2)
        {
            k++;
            cd=d2;
        }
        else
        {
            l++;
            cd=d1;
        }
    }
    else
    {
        x1=k*L*1.5;
        y1=l*L*sqrt(3.0)*0.5;
        x2=(k+1.0)*L*1.5;
        y2=(l+1.0)*L*sqrt(3.0)*0.5;
        d1=Dis(x1,y1, x,y);
        d2=Dis(x2,y2, x,y);
        if (d1>d2)
        {
            k++,l++;
            cd=d2;
        }
        else cd=d1;
    }
}

int My_Abs(int x)
{
    if (x<0) return -x;
    return x;
}

int main (void)
{
    double L, x1, y1, x2, y2, ans, cd1, cd2;
    int k1, l1, k2, l2;
    while (scanf("%lf %lf %lf %lf %lf",&L,&x1,&y1,&x2,&y2)==5)
    {
        if (L==0.0&&x1==0.0&&y1==0.0&&x2==0.0&&y2==0.0) break;
        Get_KL(L, x1, y1, k1, l1, cd1);
        Get_KL(L, x2, y2, k2, l2, cd2);
        if (k1==k2&&l1==l2) printf("%.3lf\n", Dis(x1,y1, x2,y2));
        else
        {
            ans=cd1+cd2;
            if (My_Abs(k1-k2) > My_Abs(l1-l2))
                ans=ans+sqrt(3.0)*L*My_Abs(k1-k2);
            else
                ans=ans+sqrt(3.0)*L*My_Abs(k1-k2)
                    +sqrt(3.0)*L*(double)(My_Abs(l1-l2)-My_Abs(k1-k2))/2.0
                    ;
            printf("%.3lf\n", ans);
        }
    }
    return 0;
}

//第二种六边形坐标系
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

struct A{
    int x, y, num;
}a[10001];

const int dec[6][2] = {{-1,1},{-1,0},{0,-1},{1,-1},{1,0},{0,1}};

bool adj(int x1, int y1, int x2, int y2) {
    if (x1 == x2 && abs(y1-y2) == 1) return true;
    if (y1 == y2 && abs(x1-x2) == 1) return true;
    if (x1 == x2 + 1 && y1 == y2 -1) return true;
    if (x1 == x2 - 1 && y1 == y2 +1) return true;
    return false;
}

bool flag[10001];

int main ()
{
    int i, j, k, x, u, v, cut, minn, cnt[6];
    memset(cnt, 0, sizeof(cnt));
    a[1].num = 1, cnt[1] = 1;
    a[1].x = a[1].y = 0;
    for (i = 2; i < 10001; i ++)
    {
        k = (int)((3.0+sqrt(12.0*i - 3.0))/6.0+0.0000001);
        if (i == 3*(k-1)*(k-1)+3*(k-1)+1) k --;
        j = i - (3*(k-1)*(k-1)+3*(k-1)+1);
        //当前的六边形是第 k 层的第 j 个六边形
        if (j == 1) a[i].x = a[i-1].x, a[i].y = a[i-1].y + 1;
        else
        {
            x = (j-1) / k;
            a[i].x = a[i-1].x + dec[x][0];
            a[i].y = a[i-1].y + dec[x][1];
        }
        memset(flag, false, sizeof(flag));
        x = 12*k-6, cut = 0;
        for (u = i-1, v = 0; u>=1&&v<x; u --, v ++)
            if (adj(a[u].x, a[u].y, a[i].x, a[i].y))
            {
                cut ++;
                flag[a[u].num] = true;
                if (cut == 3) break;
            }
        minn = 10001;
        for (u = 1; u < 6; u ++)
            if ((!flag[u])&&minn > cnt[u])
            {
                minn = cnt[u];
                x = u;
            }
        a[i].num = x;
        cnt[x] ++;
    }
    scanf ("%d", &x);
    while (x --)
    {
        scanf ("%d", &i);
        printf ("%d\n", a[i].num);
    }
    return 0;
}
```

## 用一个给定半径的圆覆盖最多的点

```c++
//同半径圆的圆弧表示
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define PI acos(-1.0)

struct Point{
    double x, y;
}pt[2005];

double dis[2005][2005];

struct List{
    double a;
    bool flag;
    int id;
}list[8005];

int cnt;

double Dis(int i, int j)
{
    double dx=pt[i].x-pt[j].x;
    double dy=pt[i].y-pt[j].y;
    return sqrt(dx*dx+dy*dy);
}

int Cmp(const void*p1, const void*p2)
{
    struct List*a1=(struct List*)p1;
    struct List*a2=(struct List*)p2;
    if (a1->a<a2->a)return -1;
    else if (a1->a==a2->a) return a1->id-a2->id;
    else return 1;
}

int main (void)
{
    int n, i, j, ans, num;
    double r, theta, delta, a1, a2;
    while (scanf("%d %lf",&n,&r)==2)
    {
        if (n==0&&r==0.0) break;
        r=r+0.001;
        r=r*2.0;
        for (i=1;i<=n;i++)
            scanf("%lf %lf", &pt[i].x, &pt[i].y);
        for (i=1;i<n;i++)
            for (j=i+1;j<=n;j++)
            {
                dis[i][j]=Dis(i, j);
                dis[j][i]=dis[i][j];
            }
        ans=0;
        for (i=1;i<=n;i++)
        {
            cnt=0;
            for (j=1;j<=n;j++)
                if ((j!=i)&&(dis[i][j]<=r))
                {
                    theta=atan2(pt[j].y-pt[i].y, pt[j].x-pt[i].x);
                    if (theta<0.0) theta=theta+2.0*PI;
                    delta=acos(dis[i][j]/r);
                    a1=theta-delta;
                    a2=theta+delta;
                    list[++cnt].a=a1;
                    list[cnt].flag=true;
                    list[cnt].id=cnt;
                    list[++cnt].a=a2;
                    list[cnt].flag=false;
                    list[cnt].id=cnt;
                }
            qsort(list+1,cnt,sizeof(struct List),Cmp);
            num=0;
            for (j=1;j<=cnt;j++)
                if (list[j].flag)
                {
                    num++;
                    if (num>ans) ans=num;
                }
                else num--;
        }
        printf("It is possible to cover %d points.\n", ans+1);
    }
    return 0;
}
```

## 最近圆对

```c++
#include<iostream>
#include<stdlib.h>
#include<string.h>
#include<set>
#include<math.h>

using namespace std;
set <int>tree;
set <int>::iterator iter;

struct Point{
    double x;
    int id, flag;
}p1[100001], p2[100001];

int tot1, tot2;

struct Q{
    double x,y, r;
}q[50001];

int cmp(const void*p1, const void*p2)
{
    struct Point*a1=(struct Point*)p1;
    struct Point*a2=(struct Point*)p2;
    if (a1->x<a2->x) return -1;
    else if (a1->x==a2->x) return a2->flag-a1->flag;
    else return 1;
}

int cmp1(const void*p1, const void*p2)
{
    struct Q*a1=(struct Q*)p1;
    struct Q*a2=(struct Q*)p2;
    if (a1->y<a2->y)return -1;
    else if (a1->y==a2->y)return 0;
    else return 1;
}

double dis(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

bool judge(int i, int j, double d)
{
    if (dis(q[i].x, q[i].y, q[j].x,
            q[j].y)<=q[i].r+q[j].r+2.0*d)
        return true;
    return false;
}

bool insert(int v,double d)
{
    iter = tree.insert(v).first;
    if (iter != tree.begin())
    {
        if (judge(v, *--iter,d))
        {
            return true;
        }
        ++iter;
    }
    if (++iter != tree.end())
    {
        if (judge(v, *iter,d))
        {
            return true;
        }
    }
    return false;
}

bool remove(int v,double d)
{
    iter = tree.find(v);
    if (iter != tree.begin() && iter != --tree.end())
    {
        int a = *--iter;
        ++iter;
        int b = *++iter;
        if (judge(a, b,d))
        {
            return true;
        }
    }
    tree.erase(v);
    return false;
}

bool check(double d)
{
    int i=1, j=1;
    while (i<=tot1&&j<=tot2)
    {
        if (p1[i].x-d<=p2[j].x+d)
        {
            if (insert(p1[i++].id, d))
                return true;
        }
        else
        {
            if (remove(p2[j++].id, d))
                return true;
        }
    }
    while (i<=tot1)
    {
        if (insert(p1[i++].id, d))
            return true;
    }
    while (j<=tot2)
    {
        if (remove(p2[j++].id, d))
            return true;
    }
    return false;
}

int main ()
{
    int cases, n, i;
    scanf("%d",&cases);
    while (cases--)
    {
        scanf("%d",&n);
        tot1=tot2=0;
        for (i=1;i<=n;i++)
            scanf("%lf %lf %lf",&q[i].x,&q[i].y, &q[i].r);
        qsort(q+1,n,sizeof(struct Q),cmp1);
        for (i=1;i<=n;i++)
        {
            tot1++;
            p1[tot1].x=q[i].x-q[i].r;
            p1[tot1].id=i;
            p1[tot1].flag=1;
            tot2++;
            p2[tot2].x=q[i].x+q[i].r;
            p2[tot2].id=i;
            p2[tot2].flag=-1;
        }
        qsort(p1+1,tot1,sizeof(struct Point),cmp);
        qsort(p2+1,tot2,sizeof(struct Point),cmp);
        double head=0.0, tail=dis(q[1].x,q[1].y,q[2].x,q[2].y)+1.0,
        mid;
        while (tail-head>1e-8)
        {
            tree.clear();
            mid=(head+tail)/2.0;
            if (check(mid))
            {
                tail=mid;
            }
            else head=mid;
        }
        printf ("%.6lf\n",2.0*head);
    }
    return 0;
}
```

## 求两个圆的面积交

```c++
double area_of_overlap(point c1, double r1, point c2, double r2)
{
    double a = distance(c1, c2), b = r1, c = r2;
    double cta1 = acos((a * a + b * b - c * c) / 2 / (a * b));
    double cta2 = acos((a * a + c * c - b * b) / 2 / (a * c));
    double s1 = r1*r1*cta1 - r1*r1*sin(cta1)*(a * a + b * b - c * c) / 2 / (a * b);
    double s2 = r2*r2*cta2 - r2*r2*sin(cta2)*(a * a + c * c - b * b) / 2 / (a * c);
    return s1 + s2;
}
```

# DP

## 插头DP

### 示例代码

```c++
inline long long work()
{
    memset(dp,0,sizeof(dp));
    dp[0][0][0]=1;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            if(mp[i][j])
            {
                for(int k=0;k<=mask;k++)
                {
                    if(dp[i][j][k]==0)continue;
                    int tmp=k&(3<<j);
                    if(0<tmp&&tmp<(3<<j))
                        dp[i][j+1][k]+=dp[i][j][k];
                    dp[i][j+1][k^(3<<j)]+=dp[i][j][k];
                }
            }
            else
            {
                for(int k=0;k<=mask;k++)
                {
                    if(dp[i][j][k]==0)continue;
                    if(k&(3<<j))continue;
                    dp[i][j+1][k]+=dp[i][j][k];
                }
            }
        }
        for(int k=mask>>1;k>=0;k--)dp[i+1][0][k<<1]=dp[i][m][k];//处理换行
    }
    return dp[n][0][0];
}
```

# 一些套路题

## 区间第k小

### 静态不修改(可持久化线段树+树上二分)

#### 题目

[HDU 2665 Kth number](http://acm.hdu.edu.cn/showproblem.php?pid=2665)  
数组下标从1开始,k从1开始  
复杂度$\Theta \left ( nlog \left ( n \right ) \right )$

#### 代码

```c++
#define MAXN 111111//数组大小
#define MAXM 6666666//MAXN*log(MAXN)

struct Tree{
    int num;
    int lson;
    int rson;
}tree[MAXM];//线段树
int top;

void treeInit()
{
    tree[0].num=tree[0].lson=tree[0].rson=0;//用0节点表示NULL,便于处理
    top=1;
}

int treeAdd(int ori,int left,int right,int x,int a)
{//在ori[left,right]树上x位置加a,并返回新的根
    int nown=top++;
    tree[nown]=tree[ori];
    tree[nown].num+=a;
    if(left<right)
    {
        int mid=(left+right)>>1;
        if(x<=mid)tree[nown].lson=treeAdd(tree[nown].lson,left,mid,x,a);
        else tree[nown].rson=treeAdd(tree[nown].rson,mid+1,right,x,a);
    }
    return nown;
}

int treeFind(int nown,int left,int right,int l,int r)//查询区间[l,r]
{
    if(nown==0)return 0;
    if(l<=left&&right<=r)return tree[nown].num;
    int mid=(left+right)>>1;
    int ans=0;
    if(l<=mid)ans+=treeFind(tree[nown].lson,left,mid,l,r);
    if(r>mid)ans+=treeFind(tree[nown].rson,mid+1,right,l,r);
    return ans;
}

int root[MAXN];//第i棵线段树的根,下标从1开始
int hehe[MAXN],nn;

inline int khash(int x)
{
    return lower_bound(hehe,hehe+nn,x)-hehe;
}

void treeBuild(int arr[],int n)
{
    treeInit();
    root[0]=0;
    nn=0;
    for(int i=1;i<=n;i++)
        hehe[nn++]=arr[i];
    sort(hehe,hehe+nn);
    nn=unique(hehe,hehe+n)-hehe;
    for(int i=1;i<=n;i++)
        root[i]=treeAdd(root[i-1],0,nn,khash(arr[i]),1);
}

int findKth(int lr,int rr,int left,int right,int k)
{//左边树为lr,右边树为rr,区间[left,right]中的第k大,k从1开始
    if(left==right)return left;
    int mid=(left+right)>>1;
    int tmp=tree[tree[rr].lson].num-tree[tree[lr].lson].num;
    if(tmp>=k)return findKth(tree[lr].lson,tree[rr].lson,left,mid,k);
    return findKth(tree[lr].rson,tree[rr].rson,mid+1,right,k-tmp);
}

int n,m;
int kkke[MAXN];

int main()
{
    int t,l,r,k;
    kread(t);
    while(t--)
    {
        kread(n,m);
        for(int i=1;i<=n;i++)
            kread(kkke[i]);
        treeBuild(kkke,n);
        while(m--)
        {
            kread(l,r,k);
            printf("%d\n",hehe[findKth(root[l-1],root[r],0,nn,k)]);
        }
    }
    return 0;
}
```

### 动态可修改(树状数组套线段树)

#### 题目

[ZOJ 2112 Dynamic Rankings](http://acm.zju.edu.cn/onlinejudge/showProblem.do?problemId=1112)  
数组长度n=50000,询问次数m=10000  
此题空间卡得很死,如果直接全部插入空间复杂度为$\Theta \left ( \left ( m+n \right ) log^{2}\left ( n \right ) \right )$  
所以最开始建树使用线段树合并$\Theta \left ( nlog\left ( n \right ) \right )$,然后修改时插入$\Theta \left ( mlog^{2}\left ( n \right ) \right )$  
查询时间复杂度$\Theta \left ( nlog^{3}\left ( n \right ) \right )$,若树上二分可少一个$log\left ( n \right )$

#### 代码(省略头文件和读入优化)

```c++
#define MAXN 50010//数组长度
#define MAXM 10010//询问次数
#define MAXMM 2666666//MAXM*log(MAXX)*log(MAXX)
#define MINX 0//最小数
#define MAXX nn//最大数

struct Tree{
    int num;
    int lson;
    int rson;
}tree[MAXMM];//线段树
int top;

int root[MAXN];//树状数组(存线段树的根)

int n;
int arr[MAXN];//原数组
int nn;//离散化后数的数量
int kkke[MAXN+MAXM];//离散化用

inline int lowbit(int x)
{
    return x&-x;
}

inline void treeInit()//初始化
{
    top=0;
    memset(root,-1,sizeof(root));
}

inline void treeAdd2(int& nown,int left,int right,int x,int c)
{
    if(nown==-1)
    {
        tree[top].lson=tree[top].rson=-1;
        tree[top].num=c;
    }
    else
    {
        tree[top].lson=tree[nown].lson;
        tree[top].rson=tree[nown].rson;
        tree[top].num=tree[nown].num+c;
    }
    nown=top++;
    if(left<right)
    {
        int mid=(left+right)>>1;
        if(x<=mid)treeAdd2(tree[nown].lson,left,mid,x,c);
        else treeAdd2(tree[nown].rson,mid+1,right,x,c);
    }
}

inline void treeAdd(int x,int y,int c)//时空复杂度log(MAXX)*log(MAXX)
{
    for(int i=x;i<=n;i+=lowbit(i))
        treeAdd2(root[i],MINX,MAXX,y,c);
}

inline int treeFind2(int nown,int left,int right,int l,int r)
{
    if(nown==-1)return 0;
    if(l<=left&&right<=r)return tree[nown].num;
    int mid=(left+right)>>1;
    int ans=0;
    if(l<=mid)ans+=treeFind2(tree[nown].lson,left,mid,l,r);
    if(r>mid)ans+=treeFind2(tree[nown].rson,mid+1,right,l,r);
    return ans;
}

inline int treeFind(int x,int y1,int y2)
{
    int ans=0;
    for(int i=x;i;i^=lowbit(i))
        ans+=treeFind2(root[i],MINX,MAXX,y1,y2);
    return ans;
}

inline int treeMerge(int a,int b)//合并a,b两线段树
{
    if(a==-1)return b;
    if(b==-1)return a;
    int nown=top++;
    tree[nown].num=tree[a].num+tree[b].num;
    tree[nown].lson=treeMerge(tree[a].lson,tree[b].lson);
    tree[nown].rson=treeMerge(tree[a].rson,tree[b].rson);
    return nown;
}

inline void treeBuild()//通过共用节点合并建树n*log(n)
{
    for(int i=1;i<=n;i++)treeAdd2(root[i],MINX,MAXX,arr[i],1);
    for(int i=1;i<=n;i++)
    {
        int j=i+lowbit(i);
        if(j<=n)root[j]=treeMerge(root[i],root[j]);
    }
}

inline int findKth(int l,int r,int k)//找[l,r]第k小,1<=k<=l+r
{
    int left=MINX,right=MAXX,mid;
    while(left<right)
    {
        mid=(left+right)>>1;
        int ans=treeFind(r,MINX,mid)-treeFind(l-1,MINX,mid);
        if(ans<k)left=mid+1;
        else right=mid;
    }
    return left;
}

inline void arrChange(int a,int b)//a位置变为b
{
    treeAdd(a,arr[a],-1);
    treeAdd(a,b,1);
    arr[a]=b;
}

char s[MAXM][2];
int a[MAXM];
int b[MAXM];
int c[MAXM];

int main()
{
    int t;
    int m;
    kread(t);
    while(t--)
    {
        treeInit();
        nn=0;
        kread(n,m);
        for(int i=1;i<=n;i++)
        {
            kread(arr[i]);
            kkke[nn++]=arr[i];
        }
        for(int i=0;i<m;i++)
        {
            scanf("%s",s[i]);
            if(s[i][0]=='Q')
            {
                kread(a[i],b[i],c[i]);
                if(a[i]>b[i])swap(a[i],b[i]);
            }
            else
            {
                kread(a[i],b[i]);
                kkke[nn++]=b[i];
            }
        }
        sort(kkke,kkke+nn);
        nn=unique(kkke,kkke+nn)-kkke;
        for(int i=1;i<=n;i++)arr[i]=lower_bound(kkke,kkke+nn,arr[i])-kkke;
        treeBuild();
        for(int i=0;i<m;i++)
        {
            if(s[i][0]=='Q')
            {
                printf("%d\n",kkke[findKth(a[i],b[i],c[i])]);
            }
            else
            {
                arrChange(a[i],lower_bound(kkke,kkke+nn,b[i])-kkke);
            }
        }
    }
    return 0;
}
```

### 动态可修改(整体二分)

题目同上
复杂度$\Theta \left( nlog^2 \left( n \right) \right)$

#### 代码

```c++
#define MAXN 66666

int tree[MAXN];//树状数组 

int lowbit(int a)
{
    return a&-a;
}

void tree_add(int a,int b)
{
    for(int i=a;i<MAXN;i+=lowbit(i))tree[i]+=b;
}

int tree_find(int a)
{
    int ans=0;
    for(int i=a;i;i-=lowbit(i))ans+=tree[i];
    return ans;
}

int kkke[MAXN],nk;//离散化用

int khash(int x)
{
    return lower_bound(kkke,kkke+nk,x)-kkke;
}

int n,m;
int arr[MAXN];//原数组
struct Query{
    int a,b,k,id;
    //id>=0 : 查询[a,b]中第k小,是第id个询问
    //id==-1 : a位置+1 位置上变为b
    //id==-2 : a位置-1 位置上原是b
}query[MAXN],q1[MAXN],q2[MAXN];
int ID;//询问数
int mm;//修改拆为删除和添加后的总操作数
int ans[MAXN];//储存答案

char S[22];

void work(int head,int tail,int left,int right)
{
    if(head>tail)return;
    if(left>=right)
    {
        for(int i=head;i<=tail;i++)
            if(query[i].id>=0)
                ans[query[i].id]=left;
        return;
    }
    int mid=(left+right)>>1;
    int n1=0,n2=0;
    for(int i=head;i<=tail;i++)
    {
        if(query[i].id<0)
        {
            if(query[i].b<=mid)
            {
                tree_add(query[i].a,query[i].id==-1?1:-1);
                q1[n1++]=query[i];
            }
            else q2[n2++]=query[i];
        }
        else
        {
            int tmp=tree_find(query[i].b)-tree_find(query[i].a-1);
            if(tmp>=query[i].k)q1[n1++]=query[i];
            else
            {
                query[i].k-=tmp;
                q2[n2++]=query[i];
            }
        }
    }

    for(int i=head;i<=tail;i++)//还原
        if(query[i].id<0&&query[i].b<=mid)
            tree_add(query[i].a,query[i].id==-1?-1:1);

    int nn=head;
    for(int i=0;i<n1;i++)
        query[nn++]=q1[i];
    for(int i=0;i<n2;i++)
        query[nn++]=q2[i];
    work(head,head+n1-1,left,mid);
    work(head+n1,tail,mid+1,right);
}

int main()
{
    int t;
    kread(t);
    while(t--)
    {
        memset(tree,0,sizeof(tree));
        kread(n,m);
        ID=mm=nk=0;
        for(int i=1;i<=n;i++)
        {
            kread(arr[i]);

            query[mm].a=i;
            query[mm].b=arr[i];
            query[mm++].id=-1;

            kkke[nk++]=arr[i];
        }
        int a,b;
        for(int i=0;i<m;i++)
        {
            scanf("%s",S);
            if(S[0]=='Q')
            {
                kread(query[mm].a,query[mm].b,query[mm].k);
                query[mm++].id=ID++;
            }
            else
            {
                kread(a,b);
                query[mm].a=a;
                query[mm].b=b;
                query[mm++].id=-1;

                query[mm].a=a;
                query[mm].b=arr[a];
                query[mm++].id=-2;

                kkke[nk++]=arr[a]=b;
            }
        }
        sort(kkke,kkke+nk);
        nk=unique(kkke,kkke+nk)-kkke;
        for(int i=0;i<mm;i++)
            if(query[i].id<0)
                query[i].b=khash(query[i].b);
        work(0,mm-1,0,nk-1);
        for(int i=0;i<ID;i++)
            printf("%d\n",kkke[ans[i]]);
    }
    return 0;
}
```

# STL

## 求合并,交集,并集，差集

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

## 二分查找

```c++
lower_bound()     //第一个大于等于
upper_bound()    //第一个大于
用法:
lower_bound(a.begin(),a.end(),x); //返回一个迭代器
lower_bound(a,a+n,x) //返回找到元素的指针
```

## 随机排列

```c++
template<class RandomAccessIterator>
   void random_shuffle(
      RandomAccessIterator _First, //指向序列首元素的迭代器
      RandomAccessIterator _Last  //指向序列最后一个元素的下一个位置的迭代器
   );
```

## 字符串操作

```c++
strstr(a,b)//在a中找b
```

## 读入优化

```c++
#include <cctype>

template<class TN>
inline void kread(TN &x)
{
    x=0;
    char c;
    bool flag=false;
    while(!isdigit(c=getchar()))
        if(c=='-')
            flag=true;
    do{
        x=x*10+c-48;
    }while(isdigit(c=getchar()));
    if(flag)x=-x;
}

template<class TN,class... ARGS>
inline void kread(TN &first,ARGS& ... args)
{
    kread(first);
    kread(args...);
}
```

## 哈希表unordered_set & unordered_map

### 声明

```c++
template < class Key,  
    class Hash = hash<Key>,  
    class Pred = equal_to<Key>,  
    class Alloc = allocator<Key>  
>class unordered_set;  
```

### 特例化hash类

```c++
template<class TN>
inline void hash_combine(size_t& seed, const TN &v)
{
    seed ^= hash<TN>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

//以pair<int,int>为例
namespace std{
        template<>
            struct hash<pair<int,int>>{
                typedef pair<int,int> type;
                explicit hash(){};
                size_t operator()(const type &p)const
                {
                    size_t seed=0;
                    hash_combine(seed,p.first);
                    hash_combine(seed,p.second);
                    return seed;
                }
            };
}
```

## pbds

### 优先队列

```c++
#include <ext/pb_ds/priority_queue.hpp>
typedef __gnu_pbds::priority_queue<int ,less<int>
                        ,__gnu_pbds::pairing_heap_tag>
                        Heap;
//thin_heap_tag 斐波那契堆
//pairing_heap_tag 配对堆
```

## boost

### dynamic_bitset

```c++
#include <boost/dynamic_bitset.hpp>
using namespace boost;
//可以像vector一样用push_back,resize改变长度
```

### 平衡树

```c++
#include <ext/pb_ds/assoc_container.hpp>
typedef __gnu_pbds::tree<int,__gnu_pbds::null_type, less<int>,
                    __gnu_pbds::rb_tree_tag
                    , __gnu_pbds::tree_order_statistics_node_update>
                    Tree;
//rb_tree_tag 红黑树
//splay_tag splay树
```

# Java

## a+b problem

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

## BigInteger

### 构造函数

BigInteger(String val, int radix)  
Translates the String representation of a BigInteger in the specified radix into a BigInteger.

### 方法

| 返回值            | 函数                                      | 简介                                                                                       |
|:------------------|:------------------------------------------|:-------------------------------------------------------------------------------------------|
| BigInteger        | abs()                                     | Returns a BigInteger whose value is the absolute value of this BigInteger.                 |
| BigInteger        | add(BigInteger val)                       | Returns a BigInteger whose value is (this + val).                                          |
| BigInteger        | and(BigInteger val)                       | Returns a BigInteger whose value is (this & val).                                          |
| BigInteger        | andNot(BigInteger val)                    | Returns a BigInteger whose value is (this & ~val).                                         |
| int               | compareTo(BigInteger val)                 | Compares this BigInteger with the specified BigInteger.                                    |
| BigInteger        | divide(BigInteger val)                    | Returns a BigInteger whose value is (this / val).                                          |
| BigInteger[]      | divideAndRemainder     (BigInteger val)   | Returns an array of two BigIntegers containing (this / val) followed by (this % val).      |
| double            | doubleValue()                             | Converts this BigInteger to a double.                                                      |
| boolean           | equals(Object x)                          | Compares this BigInteger with the specified Object for equality.                           |
| BigInteger        | gcd(BigInteger val)                       | Returns a BigInteger whose value is the greatest common divisor of abs(this) and abs(val). |
| BigInteger        | max(BigInteger val)                       | Returns the maximum of this BigInteger and val.                                            |
| BigInteger        | min(BigInteger val)                       | Returns the minimum of this BigInteger and val.                                            |
| BigInteger        | mod(BigInteger m)                         | Returns a BigInteger whose value is (this mod m).                                          |
| BigInteger        | modInverse(BigInteger m)                  | Returns a BigInteger whose value is (this ^ -1 mod m).                                     |
| BigInteger        | modPow(BigInteger exponent, BigInteger m) | Returns a BigInteger whose value is (this ^ exponent mod m).                               |
| BigInteger        | multiply(BigInteger val)                  | Returns a BigInteger whose value is (this * val).                                          |
| BigInteger        | negate()                                  | Returns a BigInteger whose value is (-this).                                               |
| BigInteger        | or(BigInteger val)                        | Returns a BigInteger whose value is (this &#124; val).                                     |
| BigInteger        | pow(int exponent)                         | Returns a BigInteger whose value is (this ^ exponent).                                     |
| BigInteger        | remainder(BigInteger val)                 | Returns a BigInteger whose value is (this % val).                                          |
| BigInteger        | shiftLeft(int n)                          | Returns a BigInteger whose value is (this << n).                                           |
| BigInteger        | shiftRight(int n)                         | Returns a BigInteger whose value is (this >> n).                                           |
| BigInteger        | subtract(BigInteger val)                  | Returns a BigInteger whose value is (this - val).                                          |
| String            | toString()                                | Returns the decimal String representation of this BigInteger.                              |
| String            | toString(int radix)                       | Returns the String representation of this BigInteger in the given radix.                   |
| static BigInteger | valueOf(long val)                         | Returns a BigInteger whose value is equal to that of the specified long.                   |
| BigInteger        | xor(BigInteger val)                       | Returns a BigInteger whose value is (this ^ val).                                          |

## BigDecimal

### 舍入方式

以下在roundingMode参数填入  
ROUND_CEILING向正无穷方向舍入  
ROUND_DOWN向零方向舍入  
ROUND_FLOOR向负无穷方向舍入  
ROUND_HALF_DOWN  
向（距离）最近的一边舍入，除非两边（的距离）是相等,如果是这样，向下舍入, 例如1.55 保留一位小数结果为1.5  

ROUND_HALF_EVEN  
向（距离）最近的一边舍入，除非两边（的距离）是相等,如果是这样，如果保留位数是奇数，使用ROUND_HALF_UP ，如果是偶数，使用ROUND_HALF_DOWN  

ROUND_HALF_UP  
向（距离）最近的一边舍入，除非两边（的距离）是相等,如果是这样，向上舍入, 1.55保留一位小数结果为1.6  

ROUND_UNNECESSARY 计算结果是精确的，不需要舍入模式

### 方法

| 返回值     | 函数                                                    |
|:-----------|:--------------------------------------------------------|
| BigDecimal | divide(BigDecimal divisor, int roundingMode)            |
| BigDecimal | divide(BigDecimal divisor, int scale, int roundingMode) |
| BigDecimal | setScale(int newScale)                                  |
| BigDecimal | setScale(int newScale, int roundingMode)                |

