/*!
        @file    gaugeConfig.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/


#ifndef GAUGECONFIG_INCLUDED
#define GAUGECONFIG_INCLUDED

#include <string>
using std::string;

#include "Field/field_G.h"

#include "fieldIO_Text.h"
#include "fieldIO_Text_4x4x4x8.h"
#include "fieldIO_Binary.h"
#include "fieldIO_Binary_Parallel.h"
#include "fieldIO_Binary_Distributed.h"
#include "fieldIO_Fortran.h"
#include "fieldIO_LIME.h"
#include "fieldIO_LIME_Parallel.h"
#include "fieldIO_NERSC.h"
#include "fieldIO_None.h"

#include "io_format.h"
// #include "io_format_gauge.h"

#include "bridgeIO.h"
using Bridge::vout;

//! GaugeConfig class for file I/O of gauge configuration.

/*!
   This class family is used to setup and output the gauge
   configuration.
   This is the base class, which implements common functions,
   and read/write methods are implemented in subclasses.
   At present, cutting off the gauge field for each node
   and deliver it to the node is implemented in this class.
   It may be better to separate that to other class for
   general usage for other field objects.
                                [28 Dec 2011 H.Matsufuru]

   GaugeConfig class provides file I/O of gauge configuration.
   It provides an interface to underlying FieldIO class family;

   The file format is specified by a string argument type to
   the constructor (in a somewhat similar manner as factory).
   Data layout is ILDG layout for most of the cases except
   for Fortran_JLQCD in which JLQCD layout is applied.

   "NO_OUTPUT" is added to GaugeConfig::write_file()
                                [22 Feb 2015 Y.Namekawa]
   unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
   Config types "Unit" and "Random" are introduced to generate
   unit and random gauge configurations for cold and hot start,
   respectively.
   The methods are renamed to read/write. read_file/write_file
   are left for compatibility. these methods now take args of
   Field_G* type instead of Field*.
                                [24 June 2016 T.Aoyama]
   "Null" config type is introduced that suppresses output.
   Specifying "NO_OUTPUT" for filename is also valid.
                                [8 July 2016 T.Aoyama]
   "Null" renamed by "None" in compliance with
   yaml specifications. (kept for backward compatilibity)
                                [21 November 2018 T.Aoyama]
 */

class GaugeConfig
{
 public:
  static const std::string class_name;

 public:
  GaugeConfig(const string& type);
  virtual ~GaugeConfig();

 private:
  // non-copyable
  GaugeConfig(const GaugeConfig&);
  GaugeConfig& operator=(const GaugeConfig&);

 public:
  void read(Field_G& U, const string& filename = string());

  void write(Field_G& U, const string& filename = string());


  void read_file(Field_G& U, const string& filename)
  { return read(U, filename); }

  void write_file(Field_G& U, const string& filename)
  { return write(U, filename); }

 protected:
  Bridge::VerboseLevel m_vl;
  FieldIO *m_fieldio;
  IO_Format::Format *m_format;

 private:
  class DataSource;
  class DataSource_Unit;
  class DataSource_Random;

  DataSource *m_datasource;
};
#endif
