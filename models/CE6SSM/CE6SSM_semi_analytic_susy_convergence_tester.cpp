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


#include "CE6SSM_semi_analytic_susy_convergence_tester.hpp"

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

CE6SSM_susy_convergence_tester<Semi_analytic>::CE6SSM_susy_convergence_tester(
   CE6SSM<Semi_analytic>* model, double accuracy_goal, const Scale_getter& sg)
   : Convergence_tester_DRbar<CE6SSM<Semi_analytic> >(model, accuracy_goal, sg)
{
}

double CE6SSM_susy_convergence_tester<Semi_analytic>::max_rel_diff() const
{
   const CE6SSM<Semi_analytic>& ol = get_last_iteration_model();
   const CE6SSM<Semi_analytic>& ne = get_current_iteration_model();

   std::array<double, 49> diff{};

   diff[0] = MaxRelDiff(OLD(g1),NEW(g1));
   diff[1] = MaxRelDiff(OLD(g2),NEW(g2));
   diff[2] = MaxRelDiff(OLD(g3),NEW(g3));
   diff[3] = MaxRelDiff(OLD(gN),NEW(gN));
   diff[4] = MaxRelDiff(OLD(vd),NEW(vd));
   diff[5] = MaxRelDiff(OLD(vs),NEW(vs));
   diff[6] = MaxRelDiff(OLD(vu),NEW(vu));
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 7] = MaxRelDiff(OLD2(Yd,i,j),NEW2(Yd,i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 16] = MaxRelDiff(OLD2(Ye,i,j),NEW2(Ye,i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 25] = MaxRelDiff(OLD2(Yu,i,j),NEW2(Yu,i,j));
      }
   }
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         diff[i + 3*j + 34] = MaxRelDiff(OLD2(Kappa,i,j),NEW2(Kappa,i,j));
      }
   }
   diff[43] = MaxRelDiff(OLD(Lambdax),NEW(Lambdax));
   for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
         diff[i + 2*j + 44] = MaxRelDiff(OLD2(Lambda12,i,j),NEW2(Lambda12,i,j));
      }
   }
   diff[48] = MaxRelDiff(OLD(MuPr),NEW(MuPr));

   return *std::max_element(diff.cbegin(), diff.cend());

}

} // namespace flexiblesusy
