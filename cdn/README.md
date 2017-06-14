# Huawei_RT
华为软挑

# 说明
除了sdk以外，主要有Net,Flow,SmartPtr,flow_algori等文件。  

## ***net.h/net.cpp***  
包含有ConsumNode,NetNode,Net类。Net对象持有大量的CosumeNode和NetNode，以及其他的计算所需要的变量，以及依赖于Net计算的函数方法。  

```cpp
struct ConsumNode
{
    int linkNode;	//相邻节点号
    int demand;	//消费需求
    int feasibleServiceNodeNum;
    int feasibleServiceNode[MAX_NUM_NET_NODE]={0};
    ConsumNode();
};

struct Node
{
    int  linkedNum=0;				//与之有直接连接的节点数
    int  state[MAX_NUM_LINK_NUM]={0};						//capacity de yuan shi zhuang tai
    int  linkedNode[MAX_NUM_LINK_NUM]={0};			//存储与之有直接连接的节点
    int  cost[MAX_NUM_LINK_NUM]={0};						//对应于 与之有直接连接的节点 之间 链路的 单位网络租用费
    int  capacity[MAX_NUM_LINK_NUM]={0};				//对应于 与之有直接连接的节点 之间 链路的 总带宽
    int  reachedConNodeNum=0;	//记录能达到的消费节点个数
    int	 reachedConNode[MAX_NUM_CONSUM_NODE]={0};	//记录所能到达的消费节点
    Node();
};

struct Net{
    //网络节点
    ConsumNode conNode[500];//消费节点,索引号就是节点号
    Node node[1000];//网络节点,索引号就是节点号
    int  indexOfconNode[1000]={0};//根据相邻的 节点号 索引消费节点
    int nodeNum=0,consumNodeNum=0,linkNum=0;//网络节点,消费节点,链路数量
    int serverPrice=0;//服务器成本
    SmartPtr<Flow> flow;//流配置
    /********funciton********/
    Net();
    Net(SmartPtr<Flow>& flow);
    //初始化
    void init(char** , int);
    //读入数据,构建构建
    void build(char **, int);
    //公共节点搜索函数:搜索潜在服务器节点,并标记能到达的消费节点
    void initFeasibleNode();
    //插入函数:将给定的数组,做为超级汇点输出的第一级
    void insertSuperT(void);
    //删除函数:将超级源点的输出全部删除
    void freeSuperSource(void);
    //插入函数:将给定的数组,做为超级源点输出的第一级
    void insertSuperSource(int* a,int len);
    /********搜索,排序等需要用到的数据********/
    bool ifClogged[1000];					//用于dfs函数中剪枝
    bool ifInTheTable[500][1000];	//用于标记是否能够到达 普通节点
    int  commonFlag[1000]={0};					//用于记录 可能公共节点 被几个消费节点共有
    int  countNum=0;								//用于dfs 搜索深度的计数
    int	tmpCost=0;//暂时记费用
    int serviceTableNum=0;	//每次插入的服务器个数
    int serviceTable[500];
private:
    //深度优先搜索函数:ifInTheTable中标记消费节点能到达的普通节点
    void dfs(int topNode,int demand,int n);
};

std::ostream& operator<<(std::ostream& out,Net& net);
```

## ***flow.h/flow.cpp***  
与费用流计算相关的变量都封装在Flow类内。

```cpp
struct Flow{
    Net* net;
    int MaxFlow;
    int minCost;//记录每次最小费用流的总费用
    int miniCostUntilNow;
    int	FeasiblePath[MAX_PATH_NUM][1000];	//潜在最优路径
    int FeasiblePathNum;//路径计数
    std::stack<int>	bestPath[MAX_PATH_NUM];//目前为止最优路径
    int bestPathNum;//路径计数
    Flow();
    void init(Net& net);
};
```  

## flow_algori.h/flow_algori.cpp  
所有计算的算法都被放在这个文件内，有些比较大的算法也可以写在独立的文件内。但是注意，所有的算法函数的前两个入参是`Net&`和`Flow&`类型，表示所要处理的网络和流对象。在Net和Flow对象上连续地运行这些算法，那么每个步骤都将影响其最终的结果。当开发新的算法时，只要不断地编写形如`void algorithm(Net&,Flow&)`的函数。实现算法和具体数据(Net、Flow对象)的独立，使得算法可重用并且互不影响。

```cpp
struct GeneChain						//一个个体
{
    float		rate;							//所处位置的累计比率
    long		cost;							//本条基因的费用
    int			gene[500];				//定义基因链:即染色体
    int			childGene[500];		//子代基因链
};
void geneticAlgorithm(Net&,Flow&,int geneLength);
void initFirstGeneration(Net&,Flow&,int length);//初代种群初始化
void fitnessFunction(Net&,Flow&);//定义适合度函数
void Xover(GeneChain& a,GeneChain& b,int& geneLength);
void swapArray(GeneChain& a,GeneChain& b,int length);

int minCostMaxFlow(Net&,Flow&,int s,int t);
bool spfa(Net&,Flow&,int s,int t);
```

## Timer.h/Timer.cpp
定时器模块.

```cpp
#创建对象
Timer* timer = new Timer();
//do something
timer->pass();//得到从创建对象开始时,到目前所经过的时间.单位秒数,类型double
//设置事件
timer->set_event(e,time);//e是Event类型的指针,time是一个时间值,单位秒,类型double
timer->watch_event();//当内部计时超过设定的时间time,则会出发事件
```

## SmartPtr  
智能指针类，主要是为了自动管理对象的内存释放问题。
