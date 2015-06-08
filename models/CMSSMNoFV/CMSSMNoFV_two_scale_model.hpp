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

// File generated at Mon 8 Jun 2015 18:05:01

/**
 * @file CMSSMNoFV_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 8 Jun 2015 18:05:01 with FlexibleSUSY
 * 1.1.1 (git commit: v1.1.1) and SARAH 4.5.6 .
 */

#ifndef CMSSMNoFV_TWO_SCALE_H
#define CMSSMNoFV_TWO_SCALE_H

#include "CMSSMNoFV_model.hpp"
#include "CMSSMNoFV_mass_eigenstates.hpp"
#include "two_scale_model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class CMSSMNoFV<Two_scale>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
template<>
class CMSSMNoFV<Two_scale> : public Two_scale_model, public CMSSMNoFV_mass_eigenstates {
public:
   explicit CMSSMNoFV(const CMSSMNoFV_input_parameters& input_ = CMSSMNoFV_input_parameters());
   virtual ~CMSSMNoFV();

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);

private:
};

} // namespace flexiblesusy

#endif
