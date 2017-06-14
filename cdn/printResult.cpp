//
// Created by fang.junpeng on 2017/3/29.
//
#include "printResult.h"
#include "lib/lib_io.h"
#include <stdio.h>

extern char topo_file[1000*MAX_PATH_NUM];

/*********************************************************************
**********结果输转化函数:完成路径等信息，按示例格式 写入字符串********
	//ps请注意格式的正确性，如果有解:
	//第一行只有一个数据;
	//第二行为空;
	//第三行开始才是具体的数据，数据之间用一个空格分隔开.
*********************************************************************/
//这里忽略了：“NA”情况的出现
void record_result(Flow& flow,char* str){
    size_t sizeStr=sizeof(str);
    size_t len = 0;

    if (len > (sizeStr - 10))	return;//是否超出

    len += sprintf(topo_file + len,"%d",flow.bestPathNum);
    len += sprintf(topo_file + len,"\n");

    for(int i=0;i<flow.bestPathNum;i++)
    {
        len += sprintf(topo_file + len,"\n");
        //---------格式奇怪是为了满足如例程的输出形式----------
        len += sprintf(topo_file + len,"%d",flow.bestPath[i].top());
        flow.bestPath[i].pop();
        while(!flow.bestPath[i].empty())
        {
            len += sprintf(topo_file + len," ");
            len += sprintf(topo_file + len,"%d",flow.bestPath[i].top());
            flow.bestPath[i].pop();
        }
    }
}

/*********************************************************************
***********打印结果函数:用于调试是输出，和record_result有冲突***********
*********************************************************************/
void printResult(Flow& flow){
    printf("Print path:%d\n",flow.bestPathNum);
    printf("MinCost:%d\n",flow.miniCostUntilNow);
    for(int i=0;i<flow.bestPathNum;i++)
    {
        while(!flow.bestPath[i].empty())
        {
            printf("%4d ",flow.bestPath[i].top());
            flow.bestPath[i].pop();
        }
        printf("\n");
    }
}

void write_result_exit(Flow& flow,char* str,char* filename){
    record_result(flow,str);
    write_result(str,filename);
    exit(0);
}