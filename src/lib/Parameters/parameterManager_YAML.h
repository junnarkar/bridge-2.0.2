/*!
        @file    parameterManager_YAML.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
 */

#ifndef PARAMETERMANAGER_YAML_INCLUDED
#define PARAMETERMANAGER_YAML_INCLUDED

#include "parameterManager.h"
#include <string>

//! Parameter manager with YAML parser.

/*!
   This is a simple parser to read parameters from a file
   prepared with YAML format.
   Only simple cases were checked.
                                      [17 Jul 2012 H.Matsufuru]
 */
class ParameterManager_YAML : public ParameterManager
{
 public:
  static const std::string class_name;

  ParameterManager_YAML() {}

  //! read parameters from file.
  void read_params(const std::string& params_file, Parameters& params);

  //! read parameters from input stream.
  void read_params(std::istream&, Parameters& params);
};
#endif
