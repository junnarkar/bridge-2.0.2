/*!
        @file    bridge_init_factory_alt_QXS.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-02-28 16:09:41 #$

        @version $LastChangedRevision: 2492 $
*/

#ifndef BRIDGE_INIT_FACTORY_ALT_QXS_INCLUDED
#define BRIDGE_INIT_FACTORY_ALT_QXS_INCLUDED


#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else

bool bridge_init_factory_alt_QXS();

#endif /* USE_FACTORY_AUTOREGISTER */

#ifdef DEBUG
void bridge_report_factory_alt_QXS();

#endif

bool bridge_fin_factory_alt_QXS();

#endif /* USE_FACTORY */


#endif
