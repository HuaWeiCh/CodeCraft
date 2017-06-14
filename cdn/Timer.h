//
// Created by fang.junpeng on 2017/3/30.
//

#ifndef FJP_HWRT_TIMER_H
#define FJP_HWRT_TIMER_H

struct TimerEvent{
    virtual void do_event()=0;
};

class Timer{
private:
    double _start;//从创建对象时的时间.单位,硬件滴答
    int _ms;
    unsigned long _s;
    double _event_time;
    TimerEvent* _event;
public:
    Timer();
    Timer(TimerEvent* e,double time);
    double pass();//从创建对象开始,到目前所经过的时间.单位:秒
    void reset();//重新开始计时,原有计时清零
    void set_event(TimerEvent* e, double time);
    void watch_event();
};

#endif //FJP_HWRT_TIMER_H
