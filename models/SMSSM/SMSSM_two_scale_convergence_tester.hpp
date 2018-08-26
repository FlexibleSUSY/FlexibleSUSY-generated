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

// File generated at Sun 26 Aug 2018 14:34:04

#ifndef SMSSM_TWO_SCALE_CONVERGENCE_TESTER_H
#define SMSSM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "SMSSM_convergence_tester.hpp"
#include "SMSSM_two_scale_model.hpp"

#include "convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class SMSSM_convergence_tester<Two_scale> : public Convergence_tester_DRbar<SMSSM<Two_scale> > {
public:
   using Scale_getter = Convergence_tester_DRbar<SMSSM<Two_scale>>::Scale_getter;

   SMSSM_convergence_tester(SMSSM<Two_scale>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~SMSSM_convergence_tester() = default;

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
