/*!
        @file    bridge_defs.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#ifndef DEFS_INCLUDED
#define DEFS_INCLUDED

#include <string>

// for debug
//#define LOG printf(">>> %s\n", __PRETTY_FUNCTION__)
#define LOG

// direction label
enum Direction
{
  XDIR = 0,
  YDIR = 1,
  ZDIR = 2,
  TDIR = 3,
  WDIR = 4
};

enum ForwardBackward
{
  Forward  =  1,  // +mu
  Backward = -1   // -mu
};

namespace Element_type
{
  enum type
  {
    REAL = 1, COMPLEX = 2
  };
}

#endif /* DEFS_INCLUDED */
