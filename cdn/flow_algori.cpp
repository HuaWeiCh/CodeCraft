//
// Created by fang.junpeng on 2017/3/30.
//
#include "flow_algori.h"
#include <iostream>
#include <queue>
#include <deque>
#include <cstring>
#include <unordered_map>
#include <queue>
using namespace std;

#define groupNum 3					//定义种群大小
#define NG 1000000							//种群的最大繁殖代数

float good_point[1002];
float compair_base[1002];
int		dis[1002];				//spfa函数中到各点最短距离的记录数组
int 	pre[1002];				//spfa函数中前向节点的记录数组
bool	ifVisited[1002];	//spfa函数中记录节点是否访问过

int conNode_status[MAX_NUM_CONSUM_NODE];
queue<int> uselessPath;

extern ffun::SmartPtr<Timer> timer;

unordered_map<int,pair<int,int>> cache;
unordered_map<int,float> compair[100];
int compair_count=0;

GeneChain geneChain[groupNum];	//定义种群
GeneChain solution_geneChain[groupNum];//保存目前最优的6个解
int wrost_solution_cost=INT32_MAX;
bool force=false;
void geneticAlgorithm(Net& net,Flow& flow,int geneLength){
    srand((unsigned int) time(NULL)); 				//定义随机数生成的种子
    int 		generation=0; 		//当前繁殖的最大代数
    float  	pc=0.0; 					//定义交叉的概率
    float  	base_pm=0.2; 					//定义变异的概率
    float   pm=base_pm;
    float		pro;

    int xxx= (int) clock();
    cout<<"begin()"<<endl;
    point_evaluation(net,flow);//统计每个备选节点的 flow/cost ，用于影响动态变异率
    cout<<(clock()-xxx)<<endl;

    for(int i=0;i<1002;i++) compair_base[i]=good_point[i];

    net.insertSuperT();//超级汇点现在才加 因为point_evalution要不停改变汇点

    initFirstGeneration(net,flow,geneLength);			//初代的初始化
    //当停止准则不满足 即繁殖代数没到最大代数 ,继续繁殖
    while(generation<=NG)
    {
        timer->watch_event();
        if(timer->pass()>40) {
            pm=0.1;
        }
        //选择双亲
        for(int i=0;i<groupNum;i++)//8选6?;随机,看落点
        {
            if(i==0)		//保留最优		实现位置待定
                for(int j=0;j<geneLength;j++)	geneChain[0].childGene[j]=geneChain[0].gene[j];
            else
            {
                pro= (float) ((float)rand() / (RAND_MAX + 1.0));
                if(pro<solution_geneChain[0].rate)
                    for(int j=0;j<geneLength;j++)	geneChain[i].childGene[j]=solution_geneChain[0].gene[j];
                else if(pro>=solution_geneChain[0].rate&&pro<solution_geneChain[1].rate)
                    for(int j=0;j<geneLength;j++)	geneChain[i].childGene[j]=solution_geneChain[1].gene[j];
                else if(pro>=solution_geneChain[1].rate&&pro<solution_geneChain[2].rate)
                    for(int j=0;j<geneLength;j++)	geneChain[i].childGene[j]=solution_geneChain[2].gene[j];
                else if(pro>=solution_geneChain[2].rate&&pro<solution_geneChain[3].rate)
                    for(int j=0;j<geneLength;j++)	geneChain[i].childGene[j]=solution_geneChain[3].gene[j];
                else if(pro>=solution_geneChain[4].rate&&pro<solution_geneChain[5].rate)
                    for(int j=0;j<geneLength;j++)	geneChain[i].childGene[j]=solution_geneChain[4].gene[j];
                else
                    for(int j=0;j<geneLength;j++)	geneChain[i].childGene[j]=solution_geneChain[5].gene[j];
            }
        }
        //***************************杂交算子***************************
        int r=0;
        int z=0;
        for(int j=0;j<groupNum;j++)
        {
            pro=rand()/(RAND_MAX+1.0);
            if(pro<pc)
            {
                ++z;
                if(z%2==0)	Xover(geneChain[r],geneChain[j],geneLength);
                else				r=j;
            }
        }
        //***************************变异算子***************************
        //int pos;
        for(int i=0;i<groupNum;i++)
        {
            float dynamic_mutation[1002];
            for(int j=0;j<net.nodeNum;++j) dynamic_mutation[j]=good_point[j];
            for(int j=0;j<geneLength;j++)
            {
                pro= (float) (rand() / (RAND_MAX + 1.0));//在[0,1]区间产生随机数
                if(pro<pm)
                {	//在基因链的j位进行变异
                    float total=0;
                    for(int feasible_node=0;feasible_node<net.conNode[j].feasibleServiceNodeNum;feasible_node++)
                        total+=dynamic_mutation[net.conNode[j].feasibleServiceNode[feasible_node]];
                    float x= (rand() / (float) (RAND_MAX + 1.0));
                    float cur=0;
                    for(int feasible_node=0;feasible_node<net.conNode[j].feasibleServiceNodeNum;feasible_node++){
                        cur+=(dynamic_mutation[net.conNode[j].feasibleServiceNode[feasible_node]]/total);
                        if(x<cur) {
                            geneChain[i].childGene[j]=net.conNode[j].feasibleServiceNode[feasible_node];
                            dynamic_mutation[net.conNode[j].feasibleServiceNode[feasible_node]]/=2;
                            break;
                        }
                    }
                }
            }
        }
        for(int i=0;i<groupNum;i++)
        {
            for(int j=0;j<geneLength;j++)
            {
                geneChain[i].gene[j]=geneChain[i].childGene[j];
            }
        }
        fitnessFunction(net,flow);
        generation++;
    }
}

/*********************************************************************
******初始种群初始化函数:length:基因链长度;保证每个基因链都可行*******
*********************************************************************/
void initFirstGeneration(Net& net,Flow& flow,int length)	//初代种群初始化,
{
/*    for(int i=0;i<groupNum;i++)
    {
        if(0==i) {
            for (int j = 0; j < length; j++) {
                geneChain[0].childGene[j] = net.conNode[j].linkNode;
            }
            geneChain[0].cost = flow.miniCostUntilNow;
        }
        else
        {
            for(int j=0;j<length;j++)
            {
                geneChain[i].childGene[j]=net.conNode[j].feasibleServiceNode[(int)\
					((float)net.conNode[j].feasibleServiceNodeNum*rand()/(RAND_MAX+1.0))];
            }
        }
    }*/
    for(int i=0;i<groupNum;i++) {
        for(int j=0;j<net.consumNodeNum;j++) {
            solution_geneChain[i].gene[j] = net.conNode[j].linkNode;
            geneChain[0].gene[j]=net.conNode[j].linkNode;
        }
        solution_geneChain[i].cost=flow.miniCostUntilNow;
    }
    solution_geneChain[0].rate=1;
    wrost_solution_cost=flow.miniCostUntilNow;
}

void fitnessFunction(Net& net,Flow& flow)//定义适合度函数
{
    long tmpFlow;
    double	sum;
    long tmpCost=(long)(1.01*wrost_solution_cost);

    for(int i=1;i<groupNum;i++)	//对后8个基因链依次处理
    {
        //**************根据每个个体基因链初始化服务器列表**************
        /*/--------test-------------
        geneChain[1].gene[0]=1;
        geneChain[1].gene[1]=26;
        geneChain[1].gene[2]=28;
        geneChain[1].gene[3]=22;
        geneChain[1].gene[4]=48;
        geneChain[1].gene[5]=15;
        geneChain[1].gene[6]=37;
        geneChain[1].gene[7]=1;
        //--------test-------------*/
        net.serviceTableNum=0;
        for(int j=0,k;j<net.consumNodeNum;j++)
        {
            if(j!=0)
            {
                k=0;
                while(k<net.serviceTableNum)
                {
                    if(net.serviceTable[k]==geneChain[i].gene[j])	goto	ignore;
                    k++;
                }
            }
            net.serviceTable[net.serviceTableNum]=geneChain[i].gene[j];
            net.serviceTableNum++;
            ignore:	;
        }

        //**********************最小费用最大流处理**********************
        flow.minCost=0;//操作之前应完成服务器节点的更新
        for(int j=0;j<net.nodeNum;j++)//重新初始化带宽capacity
        {
            for(int k=0;k<net.node[j].linkedNum;k++)
                net.node[j].capacity[k]=net.node[j].state[k];
        }
        for(int j=0;j<net.consumNodeNum;j++) conNode_status[j]=0;

        net.freeSuperSource();
        net.insertSuperSource(net.serviceTable,net.serviceTableNum);


        tmpFlow=minCostMaxFlow(net,flow,net.nodeNum+1,net.nodeNum);
        while(tmpFlow==-1)
            tmpFlow=minCostMaxFlow(net,flow,net.nodeNum+1,net.nodeNum);
        //****************************检查更新*************************
        if(flow.MaxFlow==tmpFlow)
        {
            geneChain[i].cost=flow.minCost;//minCost全局变量,记录最小费用流的费用
            if(flow.minCost<flow.miniCostUntilNow) {
                flow.miniCostUntilNow = flow.minCost;//更新最优费用
                cout << "miniCostUntilNow:" << flow.miniCostUntilNow <<"----"<<net.serviceTableNum<< endl;

                //更新最优路径
                for (int k = 0; k < flow.bestPathNum; k++)//清空原来的
                    while (!flow.bestPath[k].empty()) flow.bestPath[k].pop();

                for (int k = 0, v; k < flow.FeasiblePathNum; k++)//更新
                {
                    v = 0;
                    while (flow.FeasiblePath[k][v] != -1) {
                        flow.bestPath[k].push(flow.FeasiblePath[k][v]);
                        v++;
                    }
                }
                flow.bestPathNum = flow.FeasiblePathNum;

                //更新最优基因链
                swapArray(geneChain[0], geneChain[i], net.consumNodeNum);
             //   swap(geneChain[0].cost, geneChain[i].cost);
                geneChain[0].cost=geneChain[i].cost;

                for(auto it=cache.begin();it!=cache.end();++it) {
                    compair[compair_count][it->first]=it->second.first/(float)it->second.second;
                    if(it->second.second!=net.serverPrice) {
                        float v=it->second.first/(float)it->second.second;
                        if(v>good_point[it->first])
                            good_point[it->first]=v;
                    }
                }
                compair_count++;
            }
            if(flow.minCost<wrost_solution_cost) {//如果结果优于 保存结果中最差的一个解，就替换
                for(int k=0;k<groupNum;k++) {
                    if(solution_geneChain[k].cost==wrost_solution_cost) {
                        for(int l=0;l<net.consumNodeNum;l++)
                            solution_geneChain[k].gene[l]=geneChain[i].gene[l];
                        solution_geneChain[k].cost=flow.minCost;
                        break;
                    }
                }
                int max_cost=INT32_MIN;
                for(int k=0;k<groupNum;k++) {
                    if(solution_geneChain[k].cost>max_cost){
                        max_cost=solution_geneChain[k].cost;
                        break;
                    }
                }
                max_cost=INT32_MIN;
                for(int k=0;k<groupNum;k++) {
                    if(solution_geneChain[k].cost>max_cost)
                        max_cost=solution_geneChain[k].cost;
                }
                wrost_solution_cost=max_cost;
            }
        }
        cache.clear();
    }
    //************************计算适应度*****************************

    sum=0.0;
    //因为总费用越小越优，这里做了下处理
    for(int i=0;i<groupNum;i++) sum+=tmpCost-solution_geneChain[i].cost;
    for(int i=0;i<groupNum;i++)
    {
        int extra = (int) (tmpCost - solution_geneChain[i].cost);
        if(0==i)
            solution_geneChain[i].rate= (float) (extra / sum);
        else
            solution_geneChain[i].rate= (float) (extra / sum + solution_geneChain[i - 1].rate);
    }
}

/*********************************************************************
***********基因杂交函数:根据选定的双亲，完成基因的交叉操作************
*********************************************************************/
void Xover(GeneChain& a,GeneChain& b,int& geneLength)
{
    int pos; 				//随机生成杂交点 即第几个分量进行相互交换
    int tmp;				//用于交换的暂存变量
    pos=(int)(geneLength*rand()/(RAND_MAX+1.0)); //在n个分量中，随机确定第pos个分量
    for(int i=0;i<pos;i++)
    {
        tmp=a.childGene[i];
        a.childGene[i]=b.childGene[i];
        b.childGene[i]=tmp;
    }
}



void point_evaluation(Net& net,Flow& flow) {
    for(int i=0;i<net.nodeNum;i++) good_point[i]=0;//初始化

    for(int i_con=0;i_con<net.consumNodeNum;i_con++) {//所有消费节点
        for(int j_node=0;j_node<net.conNode[i_con].feasibleServiceNodeNum;j_node++) {//消费节点i_con的备选网络网络节点
            int n = net.conNode[i_con].feasibleServiceNode[j_node];
            if(good_point[n]==0) {
                int n_flow=0;
                int n_cost=net.serverPrice/100;

                for(int k_con=0;k_con<net.node[n].reachedConNodeNum;k_con++) {//n所能reach的消费节点的 cost和flow量
                    int user = net.node[n].reachedConNode[k_con];//消费节点ID
                    int tmp=net.conNode[user].linkNode;//网络节点ID
                    net.node[tmp].linkedNode[net.node[tmp].linkedNum]=net.nodeNum;
                    net.node[tmp].capacity[net.node[tmp].linkedNum]=net.conNode[user].demand;
                    net.node[tmp].state[net.node[tmp].linkedNum]=net.conNode[user].demand;
                    net.node[tmp].cost[net.node[tmp].linkedNum]=0;
                    net.node[tmp].linkedNum++;
                }

                for(int i=0;i<net.nodeNum;i++)//重新初始化带宽capacity
                {
                    for(int j=0;j<net.node[i].linkedNum;j++)
                        net.node[i].capacity[j]=net.node[i].state[j];
                }

                while(spfa(net,flow,n,net.nodeNum)) {
                    int minflow = INT32_MAX;
                    for(int i = pre[net.nodeNum],j = net.nodeNum;j!=n;i = pre[i])
                    {
                        for(int k=0;k<net.node[i].linkedNum;k++)
                        {
                            if((net.node[i].linkedNode[k]==j) && (net.node[i].capacity[k]<minflow))
                                minflow = net.node[i].capacity[k];
                        }
                        j=i;
                    }
                    n_flow+=minflow;
                    n_cost+=dis[net.nodeNum]*minflow;
                    for(int i = pre[net.nodeNum],j = net.nodeNum;j!=n;i = pre[i])
                    {
                        for(int k=0;k<net.node[i].linkedNum;k++)
                        {
                            if(net.node[i].linkedNode[k]==j)
                                net.node[i].capacity[k] -= minflow;
                        }
                        j=i;
                    }
                }

                good_point[n]=(float)n_flow/n_cost;
                if(net.conNode[i_con].linkNode==n) good_point[n]*=5;
                for(int k_con=0;k_con<net.node[n].reachedConNodeNum;k_con++) {//去掉超级汇点
                    int usr=net.node[n].reachedConNode[k_con];
                    int tmp=net.conNode[usr].linkNode;
                    net.node[tmp].linkedNum--;
                }

            }
        }
    }
}

void swapArray(GeneChain& a,GeneChain& b,int length)
{
    int tmp;
    for(int i=0;i<length;i++)
    {
    //    tmp=a.gene[i];
        a.gene[i]=b.gene[i];
    //    b.gene[i]=tmp;
    }
}



int minCostMaxFlow(Net& net,Flow& flow,int s,int t){
    int	count;
    int max_flow = 0; // 总流量
    int minflow;
    if(net.serviceTableNum==net.consumNodeNum) return (int) (flow.miniCostUntilNow * 1.5);

    if(flow.FeasiblePathNum>0)
    {
        for(int i=0,j;i<flow.FeasiblePathNum;i++)//重新初始化FeasiblePath及计数
        {
            j=0;
            while(flow.FeasiblePath[i][j] != -1 && j<1000){
                flow.FeasiblePath[i][j]=-1;
                j++;
            }
        }
        flow.FeasiblePathNum=0;
    }

    //spfa初始化:每次费用流的首次，以后重写数据将有所不同
    //for(int i=0;i<nodeNum+2;i++)  pre[i]=-1;
    //for(int i=0;i<nodeNum+2;i++)  dis[i]=INF;
    //for(int i=0;i<nodeNum+2;i++)  ifVisited[i]=false;

    //spfa(s,t);
    std::unordered_map<int,int> m;
    while(spfa(net,flow,s,t))
    {
        minflow = INF + 1;
        for(int i = pre[t],j = t;i != s;i = pre[i])
        {
            for(int k=0;k<net.node[i].linkedNum;k++)
            {
                if(net.node[i].linkedNode[k]==j)
                {
                    if(net.node[i].capacity[k] < minflow)
                    {
                        minflow = net.node[i].capacity[k];
                        //reWriteStandard = dis[j];
                        //cout<<reWriteStandard<<" ";
                    }
                }
            }
            j=i;
        }

        //以下为记录路径等信息
        flow.FeasiblePath[flow.FeasiblePathNum][0]=minflow;//加入本条路的流量
        flow.FeasiblePath[flow.FeasiblePathNum][1]=net.indexOfconNode[pre[t]];//加入消费节点
        count=2;
        int Server = 999;
        for(int i = pre[t],j = t;i != s;i = pre[i])//倒推路径求出这次用了那个服务器，便于下面求出真正使用的服务器数量//看着好冗余，可优化
        {
            for(int k=0;k<net.node[i].linkedNum;k++)
            {
                if(net.node[i].linkedNode[k]==j)
                {
                    net.node[i].capacity[k] -= minflow;
                    flow.FeasiblePath[flow.FeasiblePathNum][count]=i;
                    Server=i;
                    count++;
                }
            }
            j=i;
        }
        int conNode_ID=flow.FeasiblePath[flow.FeasiblePathNum][1];
        conNode_status[conNode_ID]+=dis[t]*minflow;
        //------------------------------------
        if(conNode_status[conNode_ID]>net.serverPrice && net.serviceTableNum!=net.consumNodeNum) {
            flow.minCost=0;//操作之前应完成服务器节点的更新
            for(int j=0;j<net.nodeNum;j++)//重新初始化带宽capacity
            {
                for(int k=0;k<net.node[j].linkedNum;k++)
                    net.node[j].capacity[k]=net.node[j].state[k];
            }
            for(int j=0;j<net.consumNodeNum;j++) conNode_status[j]=0;
            net.node[net.nodeNum+1].linkedNode[net.node[net.nodeNum+1].linkedNum]=net.conNode[conNode_ID].linkNode;
            net.node[net.nodeNum+1].capacity[net.node[net.nodeNum+1].linkedNum]=INF;
            net.node[net.nodeNum+1].state[net.node[net.nodeNum+1].linkedNum]=INF;
            net.node[net.nodeNum+1].cost[net.node[net.nodeNum+1].linkedNum]=0;
            net.node[net.nodeNum+1].linkedNum++;

            return -1;
        }
        //------------------------------------
        flow.FeasiblePathNum++;
        m[Server]++;
        if(m[Server]==1){
            flow.minCost+=net.serverPrice;
        }
        max_flow += minflow;
        flow.minCost += dis[t] * minflow;
        //----------------------
        if(m[Server]==1) {
            cache[Server].first=0;
            cache[Server].second=net.serverPrice/100;
        }
        cache[Server].first+=minflow;
        cache[Server].second+=dis[t]*minflow;
        //----------------------
        if(flow.minCost>wrost_solution_cost){
            return max_flow;
        }
    }
    int breakpoint=0;
    return	max_flow;// 最大流
}

bool spfa(Net& net,Flow& flow,int s,int t)
{
    int u,v;
    deque<int>q;

    //spfa³õÊŒ»¯
    for(int i=0;i<net.nodeNum+2;i++)  pre[i]=-1;
    for(int i=0;i<net.nodeNum+2;i++)  dis[i]=INF;
    for(int i=0;i<net.nodeNum+2;i++)  ifVisited[i]=false;

    dis[s] = 0;
    ifVisited[s] = true;
    q.push_back(s);          //superSourceÈë¶Ó

    //spfa³õÊŒ»¯£¬»òÖØÐŽÊý×éÄÚÒ»Ð©ÊýŸÝ£¬Ÿù·ÅÔÚspfaº¯ÊýÍâ
    while(!q.empty())
    {
        u = q.front();
        q.pop_front();
        ifVisited[u] = false;//release from queue q
        if(u != net.nodeNum){   //³¬Œ¶»ãµãÃ»ÓÐÏòÍâ·Ö·¢¹ŠÄÜ
            for(int i = 0; i < net.node[u].linkedNum; i++) {
                v = net.node[u].linkedNode[i];
                if((net.node[u].capacity[i]>0) && (dis[v] >= dis[u]+net.node[u].cost[i])) {
                    dis[v] = dis[u] + net.node[u].cost[i];
                    pre[v] = u;				//ÓÃÓÚÇ°ÏòËÑË÷
                    if(!ifVisited[v]) {
                        ifVisited[v] = true;
                        if(!q.empty()) {
                            if(dis[v]>dis[q.front()]) {
                                q.push_back(v);
                            } else {
                                q.push_front(v);
                            }
                        } else {
                            q.push_back(v);
                        }
                    }
                }
            }
        }
        ifVisited[u]=false;
    }
    if(dis[t] == INF)	return false;
    else               return true;
}
