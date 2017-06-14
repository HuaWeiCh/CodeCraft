//
// Created by fang.junpeng on 2017/3/30.
//

#include "Timer.h"
#include <time.h>
#include <sys/timeb.h>

#define HUAWEI_TIME

Timer::Timer() {
    reset();
    set_event(nullptr,0.0);
}

Timer::Timer(TimerEvent* e, double time) {
    reset();
    set_event(e,time);
}

void Timer::set_event(TimerEvent* e, double time) {
    _event = e;
    _event_time = time;
}

double Timer::pass() {
#ifdef C_TIME
    double time = (double)(clock()-_start)/CLOCKS_PER_SEC;
    return time;
#endif

#ifdef HUAWEI_TIME
    //这个版本的代码是和华为一样的计时方式
    struct timeb rawtime;
    struct tm * timeinfo;
    ftime(&rawtime);
    double out_ms = (double)(rawtime.millitm - _ms);
    int out_s = rawtime.time - _s;

    if (out_ms < 0)
    {
        out_ms += 1000;
        out_s -= 1;
    }
    double time = out_s+out_ms/1000;
    return time;
#endif
}

void Timer::reset() {
#ifdef C_TIME
    _start = (double)(clock());
#endif
#ifdef HUAWEI_TIME
    struct timeb rawtime;
    struct tm * timeinfo;
    ftime(&rawtime);
    _ms = rawtime.millitm;
    _s = rawtime.time;
#endif
}

void Timer::watch_event() {
    if(nullptr != _event){
        if(pass() > _event_time){
            _event->do_event();
        }
    }
}