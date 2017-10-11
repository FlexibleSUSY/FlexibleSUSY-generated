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

// File generated at Tue 10 Oct 2017 22:34:37

/**
 * @file E6SSMEFTHiggs_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Tue 10 Oct 2017 22:34:37 with FlexibleSUSY
 * 2.0.0 (git commit: e7cd01524dc37f9ba34ce6090bb584b8c724259f) and SARAH 4.12.0 .
 */

#ifndef E6SSMEFTHiggs_TWO_SCALE_H
#define E6SSMEFTHiggs_TWO_SCALE_H

#include "E6SSMEFTHiggs_model.hpp"
#include "E6SSMEFTHiggs_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class E6SSMEFTHiggs<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class E6SSMEFTHiggs<Two_scale> : public Model, public E6SSMEFTHiggs_mass_eigenstates {
public:
   explicit E6SSMEFTHiggs(const E6SSMEFTHiggs_input_parameters& input_ = E6SSMEFTHiggs_input_parameters());
   E6SSMEFTHiggs(const E6SSMEFTHiggs&) = default;
   E6SSMEFTHiggs(E6SSMEFTHiggs&&) = default;
   virtual ~E6SSMEFTHiggs() = default;
   E6SSMEFTHiggs& operator=(const E6SSMEFTHiggs&) = default;
   E6SSMEFTHiggs& operator=(E6SSMEFTHiggs&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const E6SSMEFTHiggs<Two_scale>&);

} // namespace flexiblesusy

#endif