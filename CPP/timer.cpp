// ****************************************************************************
// File: timer.cpp
//
// Purpose:
//   Basic timer method
//
// Programmer:  Kyle Spafford
// Creation:    November 19, 2010
//
// ****************************************************************************

#include <stddef.h>
#include <sys/time.h>
double curr_second (void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}
