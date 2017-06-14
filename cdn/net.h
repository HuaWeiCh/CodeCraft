#include "SmartPtr.h"
#include <iostream>

using namespace ffun;

struct Flow;

#define MAX_NUM_NET_NODE 1005
#define MAX_NUM_CONSUM_NODE 505
#define MAX_NUM_LINK_NUM 1005

#ifndef NET_H
#define NET_H

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

//输出网络,主要用于Net.init()是否正常
std::ostream& operator<<(std::ostream& out,Net& net);

#endif

