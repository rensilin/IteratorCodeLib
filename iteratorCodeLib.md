# 实用数据结构

## 加权并查集

>解决集合问题中，集合内元素有关系并且关系具有传递性的问题。

>从集合中删除节点的方法：消除该点对集合的影响(如集合中的点个数、和、最值)，然后给它分配一个新的编号(原来的编号不管)

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

### 合并

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

## 树状数组

>要求所有数的和不能超出范围,也可修改为记录最值

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
const int maxn = 1e5 + 10;
const int M = maxn * 30;

int n, q, m, tot;
// tot:节点总个数
int a[maxn], t[maxn];
// a:原数组元素 t:将原数组元素按大小去重映射到t数组上
int T[maxn], lson[M], rson[M], c[M];
// T:树的根节点 lson:节点左孩子指针 rson:节点右孩子指针 c:树节点上的权值
```

### 离散化

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

### 线段树相关

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

### 实例

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

### 用法

```c++
int build(int l, int r)//在[l,r]上建立空树；返回空树的根
int update(int rt, int pos, int val)
//建立新树更新以rt为根节点的树上，pos节点，权值+val；返回新树的根
int query(int lrt, int rrt, int k)//返回区间[lrt,rrt]上的第k大
```

# 字符串

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

## 最大权匹配KM

### 头文件&宏&结构体&全局变量

```c++
#include <vector>

#define MAXN 66666
#define INF 0x3f3f3f3f

struct Edge{
	int to;
	int v;
	Edge(int to,int v):to(to),v(v)
	{
	}
};

int n;
vector<Edge> edge[MAXN];
int from[MAXN];
int X[MAXN];
int Y[MAXN];
bool L[MAXN];
bool R[MAXN];
int ans;
```

### 初始化
```c++
void init()
{
	ans=0;
	for(int i=0;i<n;i++)
	{
		from[i]=-1;
		Y[i]=0;
		X[i]=-INF;
		for(auto j=edge[i].begin();j!=edge[i].end();j++)
		{
			X[i]=max(X[i],j->v);
		}
		ans+=X[i];
	}
}
```

### 辅助函数

```c++
bool dfs(int nown)
{
	L[nown]=true;
	for(auto i=edge[nown].begin();i!=edge[nown].end();i++)
	{
		if(X[nown]+Y[i->to]!=i->v)continue;

		if(R[i->to])continue;
		R[i->to]=true;
		if(from[i->to]==-1||dfs(from[i->to]))
		{
			from[i->to]=nown;
			return true;
		}
	}
	return false;
}
```

### 核心代码

```c++
int KM() //自己建图edge
{
	init();
	for(int k=0;k<n;k++)
	{
		for(int i=0;i<n;i++)R[i]=L[i]=false;
		if(dfs(k))continue;
		int d=INF;
		for(int i=0;i<n;i++)
			if(L[i])
				for(auto j=edge[i].begin();j!=edge[i].end();j++)
					if(!R[j->to])
						d=min(d,X[i]+Y[j->to]-j->v);
		ans-=d;
		for(int i=0;i<n;i++)
			if(L[i])
				X[i]-=d;
		for(int i=0;i<n;i++)
			if(R[i])
				Y[i]+=d;
		k--;
	}
	return ans;
}
```

## 全局最小割SW

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
	freopen("./divide/divide5.in","r",stdin);
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
	top=0;  
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
	int cur=0;  
	int ans=INF;  
	while(bfs())  
	{  
		cur+=dfs();  
		if(cur<ans) ans=cur;  
	}   
	return ans;  
}
```

## 生成树计数

### 定理

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

### 复杂度：V×E

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

### 复杂度：sqrt(V)×E

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

## 扩展欧几里得

### 定义

>对于不完全为0的非负整数ab,gcd(a,b)表示a,b的最大公约数,必然存在整数对x,y,使得gcd(a,b)=ax+by。

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

>求a对b的逆元，即(a^(-1))mod b
>int x,y;
>exgcd(a,b,x,y);
>x即为a对b的逆元

## 矩阵快速幂

### 代码

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
                z.a[i][j] =(
						z.a[i][j] + (x.a[i][k] * y.a[k][j]) % MOD
					) % MOD;
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

## 中国剩余定理

### 定义&通式

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

> 欧拉函数是小于等于 $n$ 的正整数中与 $n$ 互质的数的数目（$\varphi \left ( 1 \right )=1$）。

> 通式：$\varphi \left ( x \right ) = x\left ( 1 - \frac{1}{p_1} \right )\left ( 1 - \frac{1}{p_2} \right )\left ( 1 - \frac{1}{p_3} \right )\cdots\left ( 1 - \frac{1}{p_n} \right )$

> 应用：欧拉降幂公式

> $a^b \equiv a^{b \  \% \  \varphi \left( n\right) + \varphi \left( n \right)} (mod\ n)\ (b > \varphi (n))$

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
## 莫比乌斯函数

### 定义

> $$ \mu = \begin{cases} 1 & n=1 \\ (-1)^k & n = p_1p_2\cdots p_k \\ 0 & other \end{cases}$$

### 莫比乌斯反演

> $$f(n) = \sum_{d,n}g(d)=\sum_{d,n} g(\frac{n}{d})$$
> $$ g(n) = \sum_{d,n} \mu(d) f(\frac{n}{d}) = \sum_{d,n} \mu(\frac{n}{d})f(d) $$
>倍数形式只用把$\frac{n}{d}$变为$\frac{d}{n}$

### 技巧
>若$g(d)=[\frac n d]*[\frac m d]$之类的阶梯状函数</br>
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
>给定一个数$n$，若存在一个与 $n$互素的 $a$,使得 $a^i(i=0,1,\cdots,\varphi(n))$在模$n$ 下两两不同,那么称$a$是$n$的一个原根。

### 性质
>$1,2,4,p^n,2p^n$有原根，其中$p$是奇素数
>一个数$n$如果有原根，原根个数为 $\varphi(\varphi(n))$
>一个数$n$的全体原根的乘积模 $n$余1
>一个数$n$的全体原根的总和模 $n$余 $\mu(n-1)$(莫比乌斯函数)

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
template<class TN>//TN为bitset
void xor_gauss(TN bits[],int n,int m)//n行m列
{
	int i=0,k=m-1;//先消高位
	//i枚举行,k枚举列
	while(i<n&&k>=0)
	{
		int nown=i;
		while(nown<n&&!bits[nown].test(k))nown++;
		if(nown<n)
		{
			for(int j=nown+1;j<n;j++)
				if(bits[j].test(k))
					bits[j]^=bits[nown];
			swap(bits[nown],bits[i]);
			i++;
		}
		k--;
	}
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
>给定B、N、P，求一个整数L满足 $B^L \equiv N \ (mod \ P)$

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
>给定B、N、P，求一个整数L满足 $B^L \equiv N \ (mod \ P)$

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

# DP

## 插头DP

### 示例

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

## 哈希表unordered_set & unordered_map

### 声明

```c++
template < class Key,  
    class Hash = hash<Key>,  
    class Pred = equal_to<Key>,  
    class Alloc = allocator<Key>  
> class unordered_set;  
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

>BigInteger(String val, int radix)
>Translates the String representation of a BigInteger in the specified radix into a BigInteger.

### 方法

| 返回值            | 函数                                      | 简介                                                                                       |
|:------------------|:------------------------------------------|:-------------------------------------------------------------------------------------------|
| BigInteger        | abs()                                     | Returns a BigInteger whose value is the absolute value of this BigInteger.                 |
| BigInteger        | add(BigInteger val)                       | Returns a BigInteger whose value is (this + val).                                          |
| BigInteger        | and(BigInteger val)                       | Returns a BigInteger whose value is (this & val).                                          |
| BigInteger        | andNot(BigInteger val)                    | Returns a BigInteger whose value is (this & ~val).                                         |
| int               | compareTo(BigInteger val)                 | Compares this BigInteger with the specified BigInteger.                                    |
| BigInteger        | divide(BigInteger val)                    | Returns a BigInteger whose value is (this / val).                                          |
| BigInteger[]      | divideAndRemainder     (BigInteger val)   | Returns an array of two BigIntegers containing (this / val) followed by (this % val).      |
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

## BigDecimal

### 舍入方式

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

### 方法

| 返回值     | 函数                                                    |
|:-----------|:--------------------------------------------------------|
| BigDecimal | divide(BigDecimal divisor, int roundingMode)            |
| BigDecimal | divide(BigDecimal divisor, int scale, int roundingMode) |
| BigDecimal | setScale(int newScale)                                  |
| BigDecimal | setScale(int newScale, int roundingMode)                |
