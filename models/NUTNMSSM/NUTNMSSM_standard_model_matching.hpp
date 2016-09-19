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

// File generated at Mon 19 Sep 2016 10:05:21

#ifndef NUTNMSSM_STANDRD_MODEL_TWO_SCALE_MATCHING_H
#define NUTNMSSM_STANDRD_MODEL_TWO_SCALE_MATCHING_H

#include "NUTNMSSM_mass_eigenstates.hpp"
#include "standard_model.hpp"

namespace flexiblesusy {

/**
 * @class standard_model_two_scale_matching
 * @brief provides matching conditions from the NUTNMSSM to the Standard Model and reverse
 */
class NUTNMSSM_standard_model_matching
{
public:
   static void match_high_to_low_scale_model_tree_level(standard_model::Standard_model&, NUTNMSSM_mass_eigenstates&, unsigned);
   static void match_high_to_low_scale_model(standard_model::Standard_model&, NUTNMSSM_mass_eigenstates&, unsigned, unsigned);
   static void match_low_to_high_scale_model_tree_level(NUTNMSSM_mass_eigenstates&, standard_model::Standard_model&);
   static void match_low_to_high_scale_model(NUTNMSSM_mass_eigenstates&, standard_model::Standard_model&, unsigned);
protected:
   static void match_high_to_low_scale_model(standard_model::Standard_model&, NUTNMSSM_mass_eigenstates&, unsigned);
   static void match_low_to_high_scale_model(NUTNMSSM_mass_eigenstates&, standard_model::Standard_model&);
};

} // namespace flexiblesusy

#endif