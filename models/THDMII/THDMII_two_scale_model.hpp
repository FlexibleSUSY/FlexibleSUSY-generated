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

// File generated at Sun 4 Aug 2019 19:21:21

/**
 * @file THDMII_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Sun 4 Aug 2019 19:21:21 with FlexibleSUSY
 * 2.4.0 (git commit: 544c83a2e6b5f23da8d0b6ccdb06f1c91f75d6eb) and SARAH 4.14.2 .
 */

#ifndef THDMII_TWO_SCALE_H
#define THDMII_TWO_SCALE_H

#include "THDMII_model.hpp"
#include "THDMII_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class THDMII<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class THDMII<Two_scale> : public Model, public THDMII_mass_eigenstates {
public:
   explicit THDMII(const THDMII_input_parameters& input_ = THDMII_input_parameters());
   THDMII(const THDMII&) = default;
   THDMII(THDMII&&) = default;
   virtual ~THDMII() = default;
   THDMII& operator=(const THDMII&) = default;
   THDMII& operator=(THDMII&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const THDMII<Two_scale>&);

} // namespace flexiblesusy

#endif
