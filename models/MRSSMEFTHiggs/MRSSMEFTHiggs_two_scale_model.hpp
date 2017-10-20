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

// File generated at Fri 20 Oct 2017 08:40:02

/**
 * @file MRSSMEFTHiggs_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Fri 20 Oct 2017 08:40:02 with FlexibleSUSY
 * 2.0.1 (git commit: 5296739235bd0ef7020eda218da9c069270c3f45) and SARAH 4.12.0 .
 */

#ifndef MRSSMEFTHiggs_TWO_SCALE_H
#define MRSSMEFTHiggs_TWO_SCALE_H

#include "MRSSMEFTHiggs_model.hpp"
#include "MRSSMEFTHiggs_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class MRSSMEFTHiggs<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class MRSSMEFTHiggs<Two_scale> : public Model, public MRSSMEFTHiggs_mass_eigenstates {
public:
   explicit MRSSMEFTHiggs(const MRSSMEFTHiggs_input_parameters& input_ = MRSSMEFTHiggs_input_parameters());
   MRSSMEFTHiggs(const MRSSMEFTHiggs&) = default;
   MRSSMEFTHiggs(MRSSMEFTHiggs&&) = default;
   virtual ~MRSSMEFTHiggs() = default;
   MRSSMEFTHiggs& operator=(const MRSSMEFTHiggs&) = default;
   MRSSMEFTHiggs& operator=(MRSSMEFTHiggs&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const MRSSMEFTHiggs<Two_scale>&);

} // namespace flexiblesusy

#endif
