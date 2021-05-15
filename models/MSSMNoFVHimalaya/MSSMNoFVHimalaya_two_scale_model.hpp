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
 * @file MSSMNoFVHimalaya_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.4 .
 */

#ifndef MSSMNoFVHimalaya_TWO_SCALE_H
#define MSSMNoFVHimalaya_TWO_SCALE_H

#include "MSSMNoFVHimalaya_model.hpp"
#include "MSSMNoFVHimalaya_model_slha.hpp"
#include "MSSMNoFVHimalaya_input_parameters.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class MSSMNoFVHimalaya<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class MSSMNoFVHimalaya<Two_scale> : public Model, public MSSMNoFVHimalaya_slha {
public:
   explicit MSSMNoFVHimalaya(const MSSMNoFVHimalaya_input_parameters& input_ = MSSMNoFVHimalaya_input_parameters(), bool do_convert_masses_to_slha = true);
   explicit MSSMNoFVHimalaya(const MSSMNoFVHimalaya_slha&, bool do_convert_masses_to_slha = true);
   MSSMNoFVHimalaya(const MSSMNoFVHimalaya&) = default;
   MSSMNoFVHimalaya(MSSMNoFVHimalaya&&) = default;
   virtual ~MSSMNoFVHimalaya() = default;
   MSSMNoFVHimalaya& operator=(const MSSMNoFVHimalaya&) = default;
   MSSMNoFVHimalaya& operator=(MSSMNoFVHimalaya&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream&) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const MSSMNoFVHimalaya<Two_scale>&);

} // namespace flexiblesusy

#endif
