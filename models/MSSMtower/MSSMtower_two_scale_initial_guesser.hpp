// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Tue 5 Sep 2017 09:47:55

#ifndef MSSMtower_STANDARD_MODEL_TWO_SCALE_INITIAL_GUESSER_H
#define MSSMtower_STANDARD_MODEL_TWO_SCALE_INITIAL_GUESSER_H

#include "MSSMtower_initial_guesser.hpp"
#include "MSSMtower_two_scale_susy_scale_constraint.hpp"
#include "standard_model_two_scale_low_scale_constraint.hpp"
#include "two_scale_initial_guesser.hpp"
#include "lowe.h"

#include <sstream>

namespace flexiblesusy {

class Two_scale;

template <class T>
class MSSMtower;

template <class T>
class StandardModel;

template <class T>
class MSSMtower_standard_model_initial_guesser;

/**
 * @class MSSMtower_standard_model_initial_guesser<Two_scale>
 * @brief initial guesser for the MSSMtower tower
 */

template<>
class MSSMtower_standard_model_initial_guesser<Two_scale> : public Initial_guesser<Two_scale> {
public:
   MSSMtower_standard_model_initial_guesser(MSSMtower<Two_scale>*,
                               standard_model::StandardModel<Two_scale>*,
                               const softsusy::QedQcd&,
                               const standard_model::Standard_model_low_scale_constraint<Two_scale>&,
                               const MSSMtower_susy_scale_constraint<Two_scale>&);
   virtual ~MSSMtower_standard_model_initial_guesser();
   virtual void guess(); ///< initial guess

   void set_running_precision(double p) { running_precision = p; }

private:
   MSSMtower<Two_scale>* model; ///< pointer to model class
   standard_model::StandardModel<Two_scale>* eft; ///< pointer to effective model class
   softsusy::QedQcd qedqcd; ///< Standard Model low-energy data
   double mu_guess; ///< guessed DR-bar mass of up-quark
   double mc_guess; ///< guessed DR-bar mass of charm-quark
   double mt_guess; ///< guessed DR-bar mass of top-quark
   double md_guess; ///< guessed DR-bar mass of down-quark
   double ms_guess; ///< guessed DR-bar mass of strange-quark
   double mb_guess; ///< guessed DR-bar mass of bottom-quark
   double me_guess; ///< guessed DR-bar mass of electron
   double mm_guess; ///< guessed DR-bar mass of muon
   double mtau_guess; ///< guessed DR-bar mass of tau
   double running_precision; ///< Runge-Kutta RG running precision
   standard_model::Standard_model_low_scale_constraint<Two_scale> low_constraint;
   MSSMtower_susy_scale_constraint<Two_scale> susy_constraint;

   void guess_eft_parameters();
   void guess_model_parameters();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
};

} // namespace flexiblesusy

#endif
