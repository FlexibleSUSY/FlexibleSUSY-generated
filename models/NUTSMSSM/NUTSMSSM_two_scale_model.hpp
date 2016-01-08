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

// File generated at Fri 8 Jan 2016 15:17:30

/**
 * @file NUTSMSSM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Fri 8 Jan 2016 15:17:30 with FlexibleSUSY
 * 1.3.1 (git commit: v1.3.1) and SARAH 4.6.0 .
 */

#ifndef NUTSMSSM_TWO_SCALE_H
#define NUTSMSSM_TWO_SCALE_H

#include "NUTSMSSM_model.hpp"
#include "NUTSMSSM_mass_eigenstates.hpp"
#include "two_scale_model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class NUTSMSSM<Two_scale>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
template<>
class NUTSMSSM<Two_scale> : public Two_scale_model, public NUTSMSSM_mass_eigenstates {
public:
   explicit NUTSMSSM(const NUTSMSSM_input_parameters& input_ = NUTSMSSM_input_parameters());
   virtual ~NUTSMSSM();

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream&) const;
   virtual void set_precision(double);
};

} // namespace flexiblesusy

#endif