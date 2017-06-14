#include "deploy.h"
#include "lib/lib_io.h"
#include "lib/lib_time.h"
#include "stdio.h"
#include "Timer.h"
#include "SmartPtr.h"

ffun::SmartPtr<Timer> timer;

int main(int argc, char *argv[])
{
    timer = new Timer();

    print_time("Begin");
    char *topo[MAX_EDGE_NUM];
    int line_num;

    char *topo_file = argv[1];
//    char *topo_file = "/Users/fang/Desktop/HUAWEI_Code_Craft_2017/case_example/b1/case0.txt";

    line_num = read_file(topo, MAX_EDGE_NUM, topo_file);

    printf("line num is :%d \n", line_num);
    if (line_num == 0)
    {
        printf("Please input valid topo file.\n");
        return -1;
    }

    char *result_file = argv[2];
    //char *result_file = "/Users/fang/Desktop/HUAWEI_Code_Craft_2017/case_example/result/fjp_r0.txt";

    deploy_server(topo, line_num, result_file);

    release_buff(topo, line_num);

    print_time("End");
    printf("all time pass:%lf",timer->pass());
	return 0;
}

