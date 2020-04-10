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

// File generated at Fri 10 Apr 2020 20:46:18

/**
 * @file MSSM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Fri 10 Apr 2020 20:46:18 with FlexibleSUSY
 * 2.4.2 (git commit: a94199e5620b8684f5d30d0eece5757a5a72c4a4) and SARAH 4.14.3 .
 */

#ifndef MSSM_TWO_SCALE_H
#define MSSM_TWO_SCALE_H

#include "MSSM_model.hpp"
#include "MSSM_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class MSSM<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class MSSM<Two_scale> : public Model, public MSSM_mass_eigenstates {
public:
   explicit MSSM(const MSSM_input_parameters& input_ = MSSM_input_parameters());
   MSSM(const MSSM&) = default;
   MSSM(MSSM&&) = default;
   virtual ~MSSM() = default;
   MSSM& operator=(const MSSM&) = default;
   MSSM& operator=(MSSM&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const MSSM<Two_scale>&);

} // namespace flexiblesusy

#endif
