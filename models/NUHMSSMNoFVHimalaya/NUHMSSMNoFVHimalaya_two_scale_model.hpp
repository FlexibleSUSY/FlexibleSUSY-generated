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

// File generated at Mon 5 Mar 2018 15:00:08

/**
 * @file NUHMSSMNoFVHimalaya_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 5 Mar 2018 15:00:08 with FlexibleSUSY
 * 2.1.0 (git commit: 8f20f6c9c42c159c1588fbc0bb3e15ce5ab6ace3) and SARAH 4.12.3 .
 */

#ifndef NUHMSSMNoFVHimalaya_TWO_SCALE_H
#define NUHMSSMNoFVHimalaya_TWO_SCALE_H

#include "NUHMSSMNoFVHimalaya_model.hpp"
#include "NUHMSSMNoFVHimalaya_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class NUHMSSMNoFVHimalaya<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class NUHMSSMNoFVHimalaya<Two_scale> : public Model, public NUHMSSMNoFVHimalaya_mass_eigenstates {
public:
   explicit NUHMSSMNoFVHimalaya(const NUHMSSMNoFVHimalaya_input_parameters& input_ = NUHMSSMNoFVHimalaya_input_parameters());
   NUHMSSMNoFVHimalaya(const NUHMSSMNoFVHimalaya&) = default;
   NUHMSSMNoFVHimalaya(NUHMSSMNoFVHimalaya&&) = default;
   virtual ~NUHMSSMNoFVHimalaya() = default;
   NUHMSSMNoFVHimalaya& operator=(const NUHMSSMNoFVHimalaya&) = default;
   NUHMSSMNoFVHimalaya& operator=(NUHMSSMNoFVHimalaya&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const NUHMSSMNoFVHimalaya<Two_scale>&);

} // namespace flexiblesusy

#endif
