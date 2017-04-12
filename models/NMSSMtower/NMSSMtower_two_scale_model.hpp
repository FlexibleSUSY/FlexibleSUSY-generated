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

// File generated at Wed 12 Apr 2017 10:47:06

/**
 * @file NMSSMtower_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Wed 12 Apr 2017 10:47:06 with FlexibleSUSY
 * 1.7.4 (git commit: bf9e92a2ddb43c203483621f6150a96a16f51536) and SARAH 4.11.0 .
 */

#ifndef NMSSMtower_TWO_SCALE_H
#define NMSSMtower_TWO_SCALE_H

#include "NMSSMtower_model.hpp"
#include "NMSSMtower_mass_eigenstates.hpp"
#include "two_scale_model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class NMSSMtower<Two_scale>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
template<>
class NMSSMtower<Two_scale> : public Two_scale_model, public NMSSMtower_mass_eigenstates {
public:
   explicit NMSSMtower(const NMSSMtower_input_parameters& input_ = NMSSMtower_input_parameters());
   virtual ~NMSSMtower();

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream& out = std::cout) const;
   virtual void set_precision(double);
};

std::ostream& operator<<(std::ostream&, const NMSSMtower<Two_scale>&);

} // namespace flexiblesusy

#endif
