#include "deploy.h"
#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <iostream>
#include "SmartPtr.h"
#include "net.h"
#include "flow.h"
#include "flow_algori.h"
#include "printResult.h"
#include "Timer.h"

using namespace ffun;

using namespace std;

char topo_file[1000*MAX_PATH_NUM];	//输出的字符数组
extern ffun::SmartPtr<Timer> timer;

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num,char * filename)
{
    SmartPtr<Flow> flow = new Flow();
    SmartPtr<Net> net = new Net(flow);
    SmartPtr<FlowEvent> event = new FlowEvent(flow.ptr(),filename);
    //设置事件
    timer->set_event(event.ptr(),85);
    //监听事件
    timer->watch_event();

    net->init(topo, line_num);
    flow->init(*net);

    cout<<"until init timepass:"<<timer->pass()<<endl;

    //cout<<*net<<endl;

    geneticAlgorithm(*net,*flow,net->consumNodeNum);//遗传算法

    cout<<"until init geneticAlgorithm timepass:"<<timer->pass()<<endl;

    cout<<"miniCostUntilNow"<<flow->miniCostUntilNow<<endl;
    record_result(*flow,topo_file);
    write_result(topo_file,filename);
}
