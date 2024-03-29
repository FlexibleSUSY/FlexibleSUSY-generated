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


#include "E6SSM_two_scale_convergence_tester.hpp"
#include <array>
#include <cmath>
#include <algorithm>
#include "wrappers.hpp"

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

E6SSM_convergence_tester<Two_scale>::E6SSM_convergence_tester(
   E6SSM<Two_scale>* model, double accuracy_goal, const Scale_getter& sg)
   : Convergence_tester_DRbar<E6SSM<Two_scale> >(model, accuracy_goal, sg)
{
}

double E6SSM_convergence_tester<Two_scale>::max_rel_diff() const
{
   const E6SSM<Two_scale>& ol = get_last_iteration_model();
   const E6SSM<Two_scale>& ne = get_current_iteration_model();

   std::array<double, 73> diff{};

   diff[0] = MaxRelDiff(OLD(MGlu),NEW(MGlu));
   diff[1] = MaxRelDiff(OLD(MChaP),NEW(MChaP));
   diff[2] = MaxRelDiff(OLD(MVZp),NEW(MVZp));
   for (int i = 0; i < 6; ++i) {
      diff[i + 3] = MaxRelDiff(OLD1(MSd,i),NEW1(MSd,i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 9] = MaxRelDiff(OLD1(MSv,i),NEW1(MSv,i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 12] = MaxRelDiff(OLD1(MSu,i),NEW1(MSu,i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 18] = MaxRelDiff(OLD1(MSe,i),NEW1(MSe,i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 24] = MaxRelDiff(OLD1(MSDX,i),NEW1(MSDX,i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 30] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (int i = 2; i < 3; ++i) {
      diff[i + 33] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 36] = MaxRelDiff(OLD1(MHpm,i),NEW1(MHpm,i));
   }
   for (int i = 0; i < 6; ++i) {
      diff[i + 38] = MaxRelDiff(OLD1(MChi,i),NEW1(MChi,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 44] = MaxRelDiff(OLD1(MCha,i),NEW1(MCha,i));
   }
   for (int i = 0; i < 3; ++i) {
      diff[i + 46] = MaxRelDiff(OLD1(MFDX,i),NEW1(MFDX,i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 49] = MaxRelDiff(OLD1(MSHI0,i),NEW1(MSHI0,i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 53] = MaxRelDiff(OLD1(MSHIp,i),NEW1(MSHIp,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 57] = MaxRelDiff(OLD1(MChaI,i),NEW1(MChaI,i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 59] = MaxRelDiff(OLD1(MChiI,i),NEW1(MChiI,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 63] = MaxRelDiff(OLD1(MSSI0,i),NEW1(MSSI0,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 65] = MaxRelDiff(OLD1(MFSI,i),NEW1(MFSI,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 67] = MaxRelDiff(OLD1(MSHp0,i),NEW1(MSHp0,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 69] = MaxRelDiff(OLD1(MSHpp,i),NEW1(MSHpp,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 71] = MaxRelDiff(OLD1(MChiP,i),NEW1(MChiP,i));
   }

   return *std::max_element(diff.cbegin(), diff.cend());

}

} // namespace flexiblesusy
