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

// File generated at Mon 27 Feb 2017 13:24:06

/**
 * @file HTHDMIIMSSMBC_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 27 Feb 2017 13:24:06 with FlexibleSUSY
 * 1.7.3 (git commit: 622a80d5da461a0a259a094325cd734ff8e79c61) and SARAH 4.9.3 .
 */

#ifndef HTHDMIIMSSMBC_TWO_SCALE_H
#define HTHDMIIMSSMBC_TWO_SCALE_H

#include "HTHDMIIMSSMBC_model.hpp"
#include "HTHDMIIMSSMBC_mass_eigenstates.hpp"
#include "two_scale_model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class HTHDMIIMSSMBC<Two_scale>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
template<>
class HTHDMIIMSSMBC<Two_scale> : public Two_scale_model, public HTHDMIIMSSMBC_mass_eigenstates {
public:
   explicit HTHDMIIMSSMBC(const HTHDMIIMSSMBC_input_parameters& input_ = HTHDMIIMSSMBC_input_parameters());
   virtual ~HTHDMIIMSSMBC();

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream& out = std::cout) const;
   virtual void set_precision(double);
};

std::ostream& operator<<(std::ostream&, const HTHDMIIMSSMBC<Two_scale>&);

} // namespace flexiblesusy

#endif
