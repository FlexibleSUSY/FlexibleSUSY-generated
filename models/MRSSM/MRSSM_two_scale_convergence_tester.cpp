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

// File generated at Tue 24 Feb 2015 17:31:32

#include "MRSSM_two_scale_convergence_tester.hpp"
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

MRSSM_convergence_tester<Two_scale>::MRSSM_convergence_tester(MRSSM<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<MRSSM<Two_scale> >(model, accuracy_goal)
{
}

MRSSM_convergence_tester<Two_scale>::~MRSSM_convergence_tester()
{
}

double MRSSM_convergence_tester<Two_scale>::max_rel_diff() const
{
   const MRSSM<Two_scale>& ol = get_last_iteration_model();
   const MRSSM<Two_scale>& ne = get_model();

   double diff[47] = { 0 };

   diff[0] = MaxRelDiff(OLD(MGlu),NEW(MGlu));
   diff[1] = MaxRelDiff(OLD(MSOc),NEW(MSOc));
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 2] = MaxRelDiff(OLD1(MSd,i),NEW1(MSd,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 8] = MaxRelDiff(OLD1(MSv,i),NEW1(MSv,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 11] = MaxRelDiff(OLD1(MSu,i),NEW1(MSu,i));
   }
   for (unsigned i = 0; i < 6; i++) {
      diff[i + 17] = MaxRelDiff(OLD1(MSe,i),NEW1(MSe,i));
   }
   for (unsigned i = 0; i < 4; i++) {
      diff[i + 23] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (unsigned i = 1; i < 4; i++) {
      diff[i + 27] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 31] = MaxRelDiff(OLD1(MRh,i),NEW1(MRh,i));
   }
   for (unsigned i = 1; i < 4; i++) {
      diff[i + 33] = MaxRelDiff(OLD1(MHpm,i),NEW1(MHpm,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 37] = MaxRelDiff(OLD1(MRpm,i),NEW1(MRpm,i));
   }
   for (unsigned i = 0; i < 4; i++) {
      diff[i + 39] = MaxRelDiff(OLD1(MChi,i),NEW1(MChi,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 43] = MaxRelDiff(OLD1(MCha1,i),NEW1(MCha1,i));
   }
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 45] = MaxRelDiff(OLD1(MCha2,i),NEW1(MCha2,i));
   }

   return *std::max_element(diff, diff + 47);

}

} // namespace flexiblesusy
