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


/**
 * @file MRSSM2_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#ifndef MRSSM2_TWO_SCALE_H
#define MRSSM2_TWO_SCALE_H

#include "MRSSM2_model.hpp"
#include "MRSSM2_model_slha.hpp"
#include "MRSSM2_input_parameters.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class MRSSM2<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class MRSSM2<Two_scale> : public Model, public MRSSM2_slha {
public:
   explicit MRSSM2(const MRSSM2_input_parameters& input_ = MRSSM2_input_parameters(), bool do_convert_masses_to_slha = true);
   explicit MRSSM2(const MRSSM2_slha&, bool do_convert_masses_to_slha = true);
   MRSSM2(const MRSSM2&) = default;
   MRSSM2(MRSSM2&&) = default;
   virtual ~MRSSM2() = default;
   MRSSM2& operator=(const MRSSM2&) = default;
   MRSSM2& operator=(MRSSM2&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream&) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const MRSSM2<Two_scale>&);

} // namespace flexiblesusy

#endif
