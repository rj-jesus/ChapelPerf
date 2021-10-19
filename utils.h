#ifndef UTILS_H
#define UTILS_H

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while(0)


double elapsed_time()
{
    static _Thread_local struct timespec previous_time;
    struct timespec current_time;
    double tdiff;

    if(clock_gettime(CLOCK_MONOTONIC, &current_time) != 0) handle_error("clock_gettime");
    tdiff = (double)current_time.tv_sec-previous_time.tv_sec+1e-9*((double)current_time.tv_nsec-previous_time.tv_nsec);
    previous_time = current_time;

    return tdiff;
}

#endif
