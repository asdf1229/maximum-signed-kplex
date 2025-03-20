#ifndef _TIMER_H_
#define _TIMER_H_

#include <cstdlib>
#include <sys/time.h>

class Timer {
public:
    Timer() { m_start = timestamp(); }
    Timer(string s) : log(s) { m_start = timestamp(); }
    void restart() { m_start = timestamp(); }
    long long elapsed() { return timestamp() - m_start; }
    // void print_time() { printf("%s : time cost = %d\n", log, integer_to_string(elapsed())); }

private:
    long long m_start;
    string log;

    // Returns a timestamp ('now') in microseconds
    long long timestamp()
    {
        struct timeval tp;
        gettimeofday(&tp, nullptr);
        return ((long long)(tp.tv_sec)) * 1000000 + tp.tv_usec;
    }
};

#endif /* _TIMER_H_ */