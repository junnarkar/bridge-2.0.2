/*!
        @file    bridgeIO.h

        @brief

        @author  Satoru Ueda (maintained by I.Kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
 */

#ifndef BRIDGEIO_INCLUDED
#define BRIDGEIO_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <cstdarg>
#include "Communicator/communicator.h"

//! BridgeIO for output under parallel environment with verbose level control.

/**
   BridgeIO provides output under parallel environment with verbose level control.
   There are 4 verbose levels as follows from strong to weak:

    CRUCIAL,
    GENERAL,
    DETAILED,
    PARANOIAC

   Multi-threaded.  [3 Mar 2015 Y.Namekawa]
 */

#define ENABLE_ILDG_TAG


namespace Bridge {
  enum VerboseLevel
  {
    CRUCIAL,
    GENERAL,
    DETAILED,
    PARANOIAC
  };

  class BridgeIO {
   public:
    static const std::string class_name;

   public:
    // Constructor
    BridgeIO(const std::string& filename = "stdout");
    virtual ~BridgeIO();

    // set output to file or stream
    void init(const std::string& filename);
    void init(const std::ostream& os);

    // pop out current output stream
    void unset();

    // verbose output for c style
    // default verbose level, node 0
    void crucial(const char *format, ...);
    void general(const char *format, ...);
    void detailed(const char *format, ...);
    void paranoiac(const char *format, ...);

    // input verbose level, node 0
    void crucial(VerboseLevel vl, const char *format, ...);
    void general(VerboseLevel vl, const char *format, ...);
    void detailed(VerboseLevel vl, const char *format, ...);
    void paranoiac(VerboseLevel vl, const char *format, ...);

    // input verbose level, input node
    void crucial(VerboseLevel vl, int node, const char *format, ...);
    void general(VerboseLevel vl, int node, const char *format, ...);
    void detailed(VerboseLevel vl, int node, const char *format, ...);
    void paranoiac(VerboseLevel vl, int node, const char *format, ...);

    void increase_indent() { ++m_indent_level; }
    void decrease_indent() { --m_indent_level; }
    int indent_level() { return m_indent_level; }
    void set_indent(const int level) { m_indent_level = level; }

    bool isOpen();

    std::ostream& getStream();

#ifdef ENABLE_ILDG_TAG
    void ildg_init(const std::ostream& os);
    void ildg_init(const std::string& filename);

    void ildg(const char *format, ...);

    std::ostream& getILDGStream();
#endif

    // convert between VerboseLevel and string expression
    static VerboseLevel set_verbose_level(const std::string& str);
    static std::string get_verbose_level(const VerboseLevel vl);

   private:

    // Hide copy constructor and assignment.
    BridgeIO(const BridgeIO&);
    BridgeIO& operator=(const BridgeIO&);

    // main method for verbose output for c style
    inline void print(VerboseLevel level, VerboseLevel write_level,
                      int node, const char *format, va_list& arg);

    // internal methods
    void rewind_();
    void tidyup_();

   private:

    // current output target
    std::ostream *os_;

    // previous output targets
    std::stack<std::ostream *> stack_;

    // workarea
    char buff_[1024];

    static int m_indent_level;

#ifdef ENABLE_ILDG_TAG
    std::ostream *ildg_os_;
#endif
  };

  extern BridgeIO vout;
}
#endif //BRIDGE_IO_INCLUDED
