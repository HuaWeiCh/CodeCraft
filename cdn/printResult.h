//
// Created by fang.junpeng on 2017/3/29.
//
#include "flow.h"
#include "net.h"

#ifndef HUAWEI_RT_PRINTRESULT_H
#define HUAWEI_RT_PRINTRESULT_H

void record_result(Flow&,char* str);

void printResult(Flow&);

//写文件,并且退出进程
void write_result_exit(Flow&,char* str,char* filename);

#endif //HUAWEI_RT_PRINTRESULT_H
