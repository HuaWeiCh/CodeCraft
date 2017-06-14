//
// Created by fang.junpeng on 2017/3/30.
//
#include "net.h"
#include "flow.h"
#include <unordered_map>
#ifndef FJP_HWRT_FLOW_ALGORI_H
#define FJP_HWRT_FLOW_ALGORI_H

//遗传算法
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
void point_evaluation(Net&,Flow&);
void swapArray(GeneChain& a,GeneChain& b,int length);

int minCostMaxFlow(Net&,Flow&,int s,int t);
bool spfa(Net&,Flow&,int s,int t);

#endif //FJP_HWRT_FLOW_ALGORI_H
