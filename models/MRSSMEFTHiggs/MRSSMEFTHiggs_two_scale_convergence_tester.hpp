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

// File generated at Tue 10 Oct 2017 21:40:34

#ifndef MRSSMEFTHiggs_TWO_SCALE_CONVERGENCE_TESTER_H
#define MRSSMEFTHiggs_TWO_SCALE_CONVERGENCE_TESTER_H

#include "MRSSMEFTHiggs_convergence_tester.hpp"
#include "MRSSMEFTHiggs_two_scale_model.hpp"

#include "convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Two_scale;

template<>
class MRSSMEFTHiggs_convergence_tester<Two_scale> : public Convergence_tester_DRbar<MRSSMEFTHiggs<Two_scale> > {
public:
   using Scale_getter = Convergence_tester_DRbar<MRSSMEFTHiggs<Two_scale>>::Scale_getter;

   MRSSMEFTHiggs_convergence_tester(MRSSMEFTHiggs<Two_scale>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~MRSSMEFTHiggs_convergence_tester() = default;

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
