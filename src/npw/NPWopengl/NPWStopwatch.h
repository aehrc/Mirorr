/*
 * NPWStopwatch.h
 *
 *  Created on: 25/05/2010
 *      Author: bro86j
 */

#ifndef NPW_STOPWATCH_H_
#define NPW_STOPWATCH_H_


// UNIX implementation
#if defined (__unix__)

    #include <sys/time.h>

    class NPWStopwatch
    {
        private:
            struct timeval start, end;
            struct timezone tz;

        public:
            /**
             * constructor, starts stopwatch
             */
            NPWStopwatch()
            {
                gettimeofday( &start, &tz );
                end = start;
            }

            /**
             * destructor
             */
            virtual ~NPWStopwatch() {}

            /**
             * starts the stopwatch, returns the current time
             */
            float reset()
            {
                float elapsedTime = lap();

                gettimeofday( &start, &tz );

                return elapsedTime;
            }

            /**
             * returns the elapsed time
             */
            float lap()
            {
                gettimeofday( &end, &tz );

                double t1, t2;
                t1 =  (double) start.tv_sec + (double) start.tv_usec * 0.000001;
                t2 =  (double) end.tv_sec + (double) end.tv_usec * 0.000001;

                return (float)( t2 - t1 );
            }
    };

#elif defined(WIN32)

    #include <windows.h>

    class NPWStopwatch
    {
        private:
            LARGE_INTEGER start;
            LARGE_INTEGER end;

        public:
            /**
             * constructor
             */
            NPWStopwatch()
            {
                start.HighPart = 0;
                start.LowPart = 0;
                start.QuadPart = 0;

                end.HighPart = 0;
                end.LowPart = 0;
                end.QuadPart = 0;

                QueryPerformanceCounter(&start);
            }

            /**
             * destructor
             */
            virtual ~NPWStopwatch() {}

            /**
             * resets the stopwatch
             */
            void reset()
            {
                float elapsedTime = lap();

                start = end;

                lap();
            }

            /**
             * returns the elapsed time
             */
            float lap()
            {
                QueryPerformanceCounter(&end);

                LARGE_INTEGER time;
                time.QuadPart = end.QuadPart - start.QuadPart;

                LARGE_INTEGER frequency;
                QueryPerformanceFrequency(&frequency);

                return (float) ((double)time.QuadPart / (double)frequency.QuadPart);
            }
    };

#endif // WINDOWS


#endif /* NPW_STOPWATCH_H_ */
