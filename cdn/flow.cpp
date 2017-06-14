//
// Created by fang.junpeng on 2017/3/30.
//

#include "flow.h"
#include "net.h"
#include "printResult.h"
#include <iostream>

extern char topo_file[1000*MAX_PATH_NUM];

using namespace std;

Flow::Flow() { }

void Flow::init(Net &net) {
    for(int i=0;i<net.consumNodeNum;i++) {
        for(int j=0;j<net.nodeNum;j++)	net.ifInTheTable[i][j]=false;
    }
    //路径暂存数组的初始化
    for(int i=0;i<MAX_PATH_NUM;i++) {
        for(int j=0;j<net.nodeNum;j++)
            FeasiblePath[i][j]=-1;
    }
    bestPathNum=net.consumNodeNum;
    miniCostUntilNow=net.serverPrice*net.consumNodeNum;//设置每个消费节点一个服务器的 初始状态
    for(int i=0;i<bestPathNum;i++)
    {
        bestPath[i].push(net.conNode[i].demand);
        bestPath[i].push(i);
        bestPath[i].push(net.conNode[i].linkNode);
    }
}

FlowEvent::FlowEvent(Flow *flow, char *filename) {
    this->flow = flow;
    this->filename = filename;
}

void FlowEvent::do_event() {
    write_result_exit(*flow,topo_file,filename);
    //cout<<"write file OK"<<endl;
}
