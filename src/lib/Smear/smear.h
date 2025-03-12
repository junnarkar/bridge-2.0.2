/*!
        @file    smear.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef SMEAR_INCLUDED
#define SMEAR_INCLUDED

#include "projection.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif


//! base class for smearing of link variables.

/*!
                            [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                            [21 Mar 2015 Y.Namekawa]
 */

class Smear
{
 public:
  Smear() {}
  virtual ~Smear() {}

 private:
  Smear(const Smear&);
  Smear& operator=(const Smear&);

 public:
  virtual void smear(Field_G&, const Field_G&) = 0;

  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

#ifdef USE_FACTORY
 public:
  typedef Smear *(*ProductCreator)(Projection *);
  typedef Smear *(*ProductCreator_params)(Projection *, const Parameters&);

  typedef FactoryTemplate<Smear, ProductCreator>          Factory;
  typedef FactoryTemplate<Smear, ProductCreator_params>   Factory_params;

  static Smear *New(const IdentifierType& subtype, Projection *proj)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)(proj) : 0;
  }

  static Smear *New(const IdentifierType& subtype, Projection *proj, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(proj, params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
#endif
