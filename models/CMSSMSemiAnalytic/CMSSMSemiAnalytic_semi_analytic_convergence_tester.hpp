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


#ifndef CMSSMSemiAnalytic_SEMI_ANALYTIC_CONVERGENCE_TESTER_H
#define CMSSMSemiAnalytic_SEMI_ANALYTIC_CONVERGENCE_TESTER_H

#include "CMSSMSemiAnalytic_convergence_tester.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_model.hpp"

#include "convergence_tester_drbar.hpp"

namespace flexiblesusy {

class Semi_analytic;

template<>
class CMSSMSemiAnalytic_convergence_tester<Semi_analytic> : public Convergence_tester_DRbar<CMSSMSemiAnalytic<Semi_analytic> > {
public:
   using Scale_getter = Convergence_tester_DRbar<CMSSMSemiAnalytic<Semi_analytic>>::Scale_getter;

   CMSSMSemiAnalytic_convergence_tester(CMSSMSemiAnalytic<Semi_analytic>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~CMSSMSemiAnalytic_convergence_tester() = default;

protected:
   virtual double max_rel_diff() const;
};

} // namespace flexiblesusy

#endif
