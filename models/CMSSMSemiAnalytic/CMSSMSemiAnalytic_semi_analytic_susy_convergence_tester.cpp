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


#include "CMSSMSemiAnalytic_semi_analytic_susy_convergence_tester.hpp"

#include "wrappers.hpp"

#include <array>
#include <cmath>
#include <algorithm>

namespace flexiblesusy {

#define OLD(p) ol.get_##p()
#define NEW(p) ne.get_##p()

#define OLD1(p,i) ol.get_##p()(i)
#define NEW1(p,i) ne.get_##p()(i)

#define OLD2(p,i,j) ol.get_##p(i,j)
#define NEW2(p,i,j) ne.get_##p(i,j)

#define OLD3(p,i,j,k) ol.get_##p(i,j,k)
#define NEW3(p,i,j,k) ne.get_##p(i,j,k)

#define OLD4(p,i,j,k,l) ol.get_##p(i,j,k,l)
#define NEW4(p,i,j,k,l) ne.get_##p(i,j,k,l)

CMSSMSemiAnalytic_susy_convergence_tester<Semi_analytic>::CMSSMSemiAnalytic_susy_convergence_tester(
   CMSSMSemiAnalytic<Semi_analytic>* model, double accuracy_goal, const Scale_getter& sg)
   : Convergence_tester_DRbar<CMSSMSemiAnalytic<Semi_analytic> >(model, accuracy_goal, sg)
{
}

double CMSSMSemiAnalytic_susy_convergence_tester<Semi_analytic>::max_rel_diff() const
{
   const CMSSMSemiAnalytic<Semi_analytic>& ol = get_last_iteration_model();
   const CMSSMSemiAnalytic<Semi_analytic>& ne = get_current_iteration_model();

   std::array<double, 33> diff{};

   diff[0] = MaxRelDiff(OLD(g1),NEW(g1));
   diff[1] = MaxRelDiff(OLD(g2),NEW(g2));
   diff[2] = MaxRelDiff(OLD(g3),NEW(g3));
   diff[3] = MaxRelDiff(OLD(vd),NEW(vd));
   diff[4] = MaxRelDiff(OLD(vu),NEW(vu));
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 5] = MaxRelDiff(OLD2(Yd,i,j),NEW2(Yd,i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 14] = MaxRelDiff(OLD2(Ye,i,j),NEW2(Ye,i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 23] = MaxRelDiff(OLD2(Yu,i,j),NEW2(Yu,i,j));
      }
   }
   diff[32] = MaxRelDiff(OLD(Mu),NEW(Mu));

   return *std::max_element(diff.cbegin(), diff.cend());

}

} // namespace flexiblesusy
