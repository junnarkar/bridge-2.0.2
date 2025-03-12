/*!
        @file    force_F.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef FORCE_F_INCLUDED
#define FORCE_F_INCLUDED

#include "Field/field_G.h"
#include "Parameters/parameters.h"

//! Base class of fermion force calculation.

/*!
    This class defines the interface of fermion force calculation.
    force_udiv() and force_udiv1() are used recursively
    to determine the smeared fermion force.
                                       [28 Dec 2011 H.Matsufuru]
    set_mode() is added. This is for the cases when the force
    calculation is nonhermitian.       [18 Jan 2012 H.Matsufuru]
*/

class Force
{
 public:
  Force() : m_U(0) {}

  virtual ~Force() {}

 private:
  // non-copyable
  Force(const Force&);
  Force& operator=(const Force&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  virtual void set_config(Field *) = 0;

  //! in Force, setting the mode is optional when H is nonhermitian.
  virtual void set_mode(const std::string& mode)
  {
    // do nothing if not defined in a subclass.
  }

  virtual void force_core(Field&, const Field&);
  virtual void force_core1(Field&, const Field&, const Field&);

  virtual void force_udiv(Field&, const Field&) {}
  virtual void force_udiv1(Field&, const Field&, const Field&) {}

 private:
  void mult_generator(Field_G&);

 protected:
  Field_G *m_U;
};
#endif
