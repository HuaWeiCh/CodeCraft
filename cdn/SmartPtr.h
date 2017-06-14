//
// Created by fang.junpeng on 2017/3/8.
//

#ifndef CODETRAIN_SMARTPTR_H
#define CODETRAIN_SMARTPTR_H

namespace ffun{
    template <typename T>
    class SmartPtr;

    template <typename T>
    class RefPtr{
    private:
        T* _ptr;//对象指针
        int _use;//引用计数
        inline void ref_inc(void){_use++;}
        inline void ref_dec(void){_use--;}
        inline int ref_cnt(void){ return _use;}
        inline T* ptr(){ return _ptr;}
    public:
        RefPtr(T* ptr):_ptr(ptr),_use(1){}
        ~RefPtr(){ delete _ptr;}
        friend class SmartPtr<T>;
    };

    /**
     * 1.声明智能指针后,应当尽快使用对应类型的指针对其进行初始化.
     * 2.类型指针和智能指针是一对一关系.即:不能使用指向同一对象的指针初始化2个以上的智能指针
     * 3.使用类型指针初始化智能指针后,由智能指针接管该对象的内存管理.即:不要在程序运行过程中主动delete该指针
     * 4.智能指针主要应用于对象类型.(基本类型也可以用,但不推荐)
     */
    template<typename T>
    class SmartPtr {
    private:
        RefPtr<T>* _rp;
        void release(void){
            _rp->ref_dec();//引用计数减1
            if(0 == _rp->ref_cnt()){delete(_rp);}//如果引用计数为0,删除对象
        }
    public:
        SmartPtr():_rp(nullptr){}
        SmartPtr(T* obj):_rp(new RefPtr<T>(obj)){}
        SmartPtr(const SmartPtr& sp):_rp(sp._rp){ _rp->ref_inc(); }
        ~SmartPtr(){ release(); }

        inline T &operator*() { return *_rp->ptr();}

        inline T *operator->() { return _rp->ptr(); }

        SmartPtr& operator=(const SmartPtr& sp){
            if(_rp == sp._rp) { return *this; }
            else{
                if(nullptr != _rp){ release();}//若已经持有_rp,则释放rp
            }
            _rp = sp._rp;
            _rp->ref_inc();//引用计数加1
            return *this;
        }
        SmartPtr& operator=(T* obj){
            if(nullptr != _rp && _rp->ptr() == obj){ return *this;}//先判断_rp是否存在,&&有短路功能
            if(nullptr != _rp){release();}//释放已持有对象
            _rp = new RefPtr<T>(obj);
            return *this;
        }

        bool isNone(){ return _rp == nullptr; }

        T* ptr(){ return  _rp == nullptr? nullptr:_rp->ptr();}

        bool operator==(const SmartPtr& sp){ return _rp==sp._rp;}
        bool operator==(const T* obj){ return  _rp != nullptr && _rp->ptr() == obj; }
    };
}
#endif //CODETRAIN_SMARTPTR_H
