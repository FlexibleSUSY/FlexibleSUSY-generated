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
 * @file TMSSM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#ifndef TMSSM_TWO_SCALE_H
#define TMSSM_TWO_SCALE_H

#include "TMSSM_model.hpp"
#include "TMSSM_model_slha.hpp"
#include "TMSSM_input_parameters.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class TMSSM<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class TMSSM<Two_scale> : public Model, public TMSSM_slha {
public:
   explicit TMSSM(const TMSSM_input_parameters& input_ = TMSSM_input_parameters(), bool do_convert_masses_to_slha = true);
   explicit TMSSM(const TMSSM_slha&, bool do_convert_masses_to_slha = true);
   TMSSM(const TMSSM&) = default;
   TMSSM(TMSSM&&) = default;
   virtual ~TMSSM() = default;
   TMSSM& operator=(const TMSSM&) = default;
   TMSSM& operator=(TMSSM&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream&) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const TMSSM<Two_scale>&);

} // namespace flexiblesusy

#endif
