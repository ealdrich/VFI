//////////////////////////////////////////////////////////////////////////////
///
/// @file timer.cpp
///
/// @brief File containing timer function.
///
//////////////////////////////////////////////////////////////////////////////

#include <stddef.h>
#include <sys/time.h>
#include "global.h"

//////////////////////////////////////////////////////////////////////////////
///
/// @fn REAL curr_second (void)
///
/// @brief Basic timer method.
///
/// @author Kyle Spafford \n
///
/// @version 1.0
///
/// @date 19 November 2010
///
/// @copyright Copyright Kyle Spafford 2012 \n
///            Distributed under the Boost Software License, Version 1.0
///            (See accompanying file LICENSE_1_0.txt or copy at \n
///            http://www.boost.org/LICENSE_1_0.txt)
///
//////////////////////////////////////////////////////////////////////////////
REAL curr_second (void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (REAL)tv.tv_sec + (REAL)tv.tv_usec / 1000000.0;
}
