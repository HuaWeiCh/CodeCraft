//
// Created by fang.junpeng on 2017/3/30.
//
#include "net.h"
#include "flow.h"
#include <cstdio>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>

using namespace std;

ConsumNode::ConsumNode() {

}

Node::Node() {

}


Net::Net() { }

Net::Net(SmartPtr<Flow> &flow) {
    this->flow = flow;
}

void Net::init(char **topo, int edge_num) {
    build(topo,edge_num);
    //一些数组的初始化
    for(int i=0;i<nodeNum;i++){
        ifClogged[i]=false;
    }
    for(int i=0;i<consumNodeNum;i++){
        for(int j=0;j<nodeNum;j++){
            ifInTheTable[i][j]=false;
        }
    }
    initFeasibleNode();
}

void Net::build(char **topo, int edge_num) {
    //*************************相关数组的初始化**********************
    sscanf(topo[0],"%d %d %d",&nodeNum,&linkNum,&consumNodeNum);//第一行
    //************************服务器单位费用************************
    sscanf(topo[2],"%d",&serverPrice);											//第三行
    //************************普通节点信息**************************
    int start,end,capacity,cost;
    for(int i=4;i<linkNum+4;i++)//存储节点、链路信息;链路包含在节点内
    {														//认为编号从0开始，且连续
        sscanf(topo[i],"%d %d %d %d",&start,&end,&capacity,&cost);
        node[start].linkedNode[node[start].linkedNum]=end;
        node[start].capacity[node[start].linkedNum]=capacity;
        node[start].state[node[start].linkedNum]=capacity;
        node[start].cost[node[start].linkedNum]=cost;
        node[start].linkedNum++;

        node[end].linkedNode[node[end].linkedNum]=start;
        node[end].capacity[node[end].linkedNum]=capacity;
        node[end].state[node[end].linkedNum]=capacity;
        node[end].cost[node[end].linkedNum]=cost;
        node[end].linkedNum++;
        //printf("%d %d %d\n",start,end,capacity);
    }
    //************************消费节点信息**************************
    for(int i=linkNum+5;i<edge_num;i++)
    {
        sscanf(topo[i],"%d %d %d",&start,&end,&capacity);
        conNode[start].linkNode=end;
        conNode[start].demand=capacity;
        //printf("%d %d %d %d\n",start,end,capacity,cost);
        //****************************************
        indexOfconNode[end]=start; //设置反向索引
        //****************************************
        flow->MaxFlow+=capacity;         //设置总的消费流量，用做最大流的检测
    }
}

void Net::initFeasibleNode() {
    for(int i=0;i<consumNodeNum;i++)//标记每个消费节点能到达的节点
    {
        ifInTheTable[i][conNode[i].linkNode]=true;

        countNum=0;
        dfs(conNode[i].linkNode,conNode[i].demand,i);
    }

    for(int i=0;i<consumNodeNum;i++)//查找可能的公共节点
    {												//commonFlag中保存节点所能到达的消费节点数
        for(int j=0;j<nodeNum;j++)
        {
            if(ifInTheTable[i][j]){
                commonFlag[j]++;//累计消费点个数
                node[j].reachedConNode[node[j].reachedConNodeNum]=i;
                node[j].reachedConNodeNum++;
            }
        }
    }
    //***保证与消费节点相邻的节点 在 潜在服务器节点表*****
//    for(int i=0;i<consumNodeNum;i++)	commonFlag[conNode[i].linkNode] +=2;

    //***************创建潜在服务器节点表*****************
    for(int i=0;i<consumNodeNum;i++)//记录每个消费节点的 潜在服务器节点
    {
        for(int j=0;j<nodeNum;j++)
        {
            if(ifInTheTable[i][j]&&((commonFlag[j]>1) || (j==conNode[i].linkNode)))//新条件代替上面的+=2
            {
                conNode[i].feasibleServiceNode[conNode[i].feasibleServiceNodeNum]=j;
                conNode[i].feasibleServiceNodeNum++;
            }
        }
    }
}

void Net::dfs(int topNode, int demand, int n) {
    if(countNum>=nodeNum)	return;
    for(int i=0;i<node[topNode].linkedNum;i++)
    {
        if(countNum>=nodeNum)	return;
        tmpCost+=node[topNode].cost[i]*demand;
        if(serverPrice>=tmpCost)
        {
            if(!ifClogged[node[topNode].linkedNode[i]])
            {
                countNum++;
                ifInTheTable[n][node[topNode].linkedNode[i]]=true;

                ifClogged[node[topNode].linkedNode[i]]=true;
                dfs(node[topNode].linkedNode[i],demand,n);
                ifClogged[node[topNode].linkedNode[i]]=false;
            }
        }
        tmpCost-=node[topNode].cost[i]*demand;
    }
}

void Net::insertSuperT() {
    for(int i=0;i<consumNodeNum;i++)
    {
        node[conNode[i].linkNode].linkedNode[node[conNode[i].linkNode].linkedNum]=nodeNum;
        node[conNode[i].linkNode].capacity[node[conNode[i].linkNode].linkedNum]=conNode[i].demand;
        node[conNode[i].linkNode].state[node[conNode[i].linkNode].linkedNum]=conNode[i].demand;
        node[conNode[i].linkNode].cost[node[conNode[i].linkNode].linkedNum]=0;
        node[conNode[i].linkNode].linkedNum++;
    }
}

void Net::freeSuperSource() {
    node[nodeNum+1].linkedNum=0;
}

void Net::insertSuperSource(int *a, int len) {
    for(int i=0;i<len;i++)
    {
        node[nodeNum+1].linkedNode[node[nodeNum+1].linkedNum]=a[i];
        node[nodeNum+1].capacity[node[nodeNum+1].linkedNum]=INF;
        node[nodeNum+1].state[node[nodeNum+1].linkedNum]=INF;
        node[nodeNum+1].cost[node[nodeNum+1].linkedNum]=0;
        node[nodeNum+1].linkedNum++;
    }
}

std::ostream& operator<<(std::ostream& out,Net& net){
    if( 0 == net.consumNodeNum || 0 == net.nodeNum){
        out<<"The net is empty"<<endl;
        return out;
    }
    //打印基本信息
    out<<net.nodeNum<<" "<<net.consumNodeNum<<" "<<net.linkNum<<endl;
    out<<endl;
    //打印服务器价格
    out<<net.serverPrice<<endl<<endl;
    //打印网络节点
    for(int i = 0;i<net.nodeNum;i++){
        Node& node = net.node[i];
        for(int j = 0;j<node.linkedNum;j++){
            out<<i<<" "<<node.linkedNode[j]<<" "<<node.capacity[j]<<" "<<node.cost[j]<<endl;
        }
    }
    out<<endl;
    //打印消费节点
    for(int i=0;i<net.consumNodeNum;i++){
        ConsumNode& node = net.conNode[i];
        out<<i<<" "<<node.linkNode<<" "<<node.demand<<endl;
    }
    return out;
}