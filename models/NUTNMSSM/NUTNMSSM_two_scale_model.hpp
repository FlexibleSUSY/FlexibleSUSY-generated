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

// File generated at Sun 28 Aug 2016 15:17:50

/**
 * @file NUTNMSSM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Sun 28 Aug 2016 15:17:50 with FlexibleSUSY
 * 1.6.1 (git commit: de78c399bea838d002fe54cbc91eb6f88fa62626) and SARAH 4.9.1 .
 */

#ifndef NUTNMSSM_TWO_SCALE_H
#define NUTNMSSM_TWO_SCALE_H

#include "NUTNMSSM_model.hpp"
#include "NUTNMSSM_mass_eigenstates.hpp"
#include "two_scale_model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class NUTNMSSM<Two_scale>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
template<>
class NUTNMSSM<Two_scale> : public Two_scale_model, public NUTNMSSM_mass_eigenstates {
public:
   explicit NUTNMSSM(const NUTNMSSM_input_parameters& input_ = NUTNMSSM_input_parameters());
   virtual ~NUTNMSSM();

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream& out = std::cout) const;
   virtual void set_precision(double);
};

std::ostream& operator<<(std::ostream&, const NUTNMSSM<Two_scale>&);

} // namespace flexiblesusy

#endif
