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

// File generated at Fri 10 Apr 2020 17:35:16

#include "MSSMNoFVHimalaya_two_scale_convergence_tester.hpp"
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

MSSMNoFVHimalaya_convergence_tester<Two_scale>::MSSMNoFVHimalaya_convergence_tester(
   MSSMNoFVHimalaya<Two_scale>* model, double accuracy_goal, const Scale_getter& sg)
   : Convergence_tester_DRbar<MSSMNoFVHimalaya<Two_scale> >(model, accuracy_goal, sg)
{
}

double MSSMNoFVHimalaya_convergence_tester<Two_scale>::max_rel_diff() const
{
   const MSSMNoFVHimalaya<Two_scale>& ol = get_last_iteration_model();
   const MSSMNoFVHimalaya<Two_scale>& ne = get_current_iteration_model();

   std::array<double, 34> diff{};

   diff[0] = MaxRelDiff(OLD(MGlu),NEW(MGlu));
   diff[1] = MaxRelDiff(OLD(MSveL),NEW(MSveL));
   diff[2] = MaxRelDiff(OLD(MSvmL),NEW(MSvmL));
   diff[3] = MaxRelDiff(OLD(MSvtL),NEW(MSvtL));
   for (int i = 0; i < 2; ++i) {
      diff[i + 4] = MaxRelDiff(OLD1(MSd,i),NEW1(MSd,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 6] = MaxRelDiff(OLD1(MSu,i),NEW1(MSu,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 8] = MaxRelDiff(OLD1(MSe,i),NEW1(MSe,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 10] = MaxRelDiff(OLD1(MSm,i),NEW1(MSm,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 12] = MaxRelDiff(OLD1(MStau,i),NEW1(MStau,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 14] = MaxRelDiff(OLD1(MSs,i),NEW1(MSs,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 16] = MaxRelDiff(OLD1(MSc,i),NEW1(MSc,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 18] = MaxRelDiff(OLD1(MSb,i),NEW1(MSb,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 20] = MaxRelDiff(OLD1(MSt,i),NEW1(MSt,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 22] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 24] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (int i = 1; i < 2; ++i) {
      diff[i + 26] = MaxRelDiff(OLD1(MHpm,i),NEW1(MHpm,i));
   }
   for (int i = 0; i < 4; ++i) {
      diff[i + 28] = MaxRelDiff(OLD1(MChi,i),NEW1(MChi,i));
   }
   for (int i = 0; i < 2; ++i) {
      diff[i + 32] = MaxRelDiff(OLD1(MCha,i),NEW1(MCha,i));
   }

   return *std::max_element(diff.cbegin(), diff.cend());

}

} // namespace flexiblesusy
