#pragma once

namespace proj {

    class StopWatch {
    
        public:
            StopWatch() : running(false) {}
            ~StopWatch() {}

            void   start();
            double stop();
            double peek();
            
        private:
            bool    running;
            clock_t started;
            clock_t stopped;
    };

    inline void StopWatch::start() {
        assert(!running);
        started = clock();
        running = true;
    }

    inline double StopWatch::stop() {
        assert(running);
        stopped = clock();
        double seconds = (double)(stopped - started)/CLOCKS_PER_SEC;
        running = false;
        started = stopped = 0L;
        return seconds;
    }

    inline double StopWatch::peek() {
        assert(running);
        clock_t peeked = clock();
        double seconds = (double)(peeked - started)/CLOCKS_PER_SEC;
        return seconds;
    }

}
