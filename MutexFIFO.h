#include <boost/thread.hpp>
#include <csignal>
#include <sstream>

// Define MutexFIFO Class
template<class T> class MutexFIFO {

// Public Members
public:
    // Destructor
    ~MutexFIFO(void);

    // Queue Operations
    void push(T);
    bool pop(T*);
    int size();

    // Mutex Operations
    void lock();
    void unlock();

// Private Members
private:
    boost::recursive_mutex mtx;         //允许同一个线程对互斥量多次上锁,即即递归上锁
    std::list<T> queue;
};


// Class Destructor
template <class T> MutexFIFO<T>::~MutexFIFO(void) {
    queue.clear();
    //mtx.destroy();
}

// Push an entity onto the FIFO queue
template <class T> void MutexFIFO<T>::push(T entry) {
    boost::lock_guard<boost::recursive_mutex> scoped_lock(mtx);
    queue.push_back(entry);            //在list的末尾添加一个元素
}

// Return the size of the FIFO queue
template <class T> int MutexFIFO<T>::size(void) {
    mtx.lock();
    int qsize = queue.size();
    mtx.unlock();
    return qsize;
}

// Pop the top element from the FIFO queue
template <class T> bool MutexFIFO<T>::pop(T* entry) {
    mtx.lock();
    bool is_empty = queue.empty();
    if (not is_empty) {
        *entry = queue.front();       //返回第一个元素
        queue.pop_front();            //删除第一个元素 
    }
    mtx.unlock();
    if (is_empty) {
        return false;
    }
    else {
        return true;
    }
}

// Lock the FIFO with a mutex
template <class T> void MutexFIFO<T>::lock(void) {
    mtx.lock();
}

// Unlock the FIFO with a mutex
template <class T> void MutexFIFO<T>::unlock(void) {
    mtx.unlock();
}