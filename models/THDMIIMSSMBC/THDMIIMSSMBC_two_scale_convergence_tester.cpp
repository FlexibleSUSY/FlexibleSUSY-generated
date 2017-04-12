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

// File generated at Wed 12 Apr 2017 11:15:54

#include "THDMIIMSSMBC_two_scale_convergence_tester.hpp"
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

THDMIIMSSMBC_convergence_tester<Two_scale>::THDMIIMSSMBC_convergence_tester(THDMIIMSSMBC<Two_scale>* model, double accuracy_goal)
   : Convergence_tester_DRbar<THDMIIMSSMBC<Two_scale> >(model, accuracy_goal)
{
}

THDMIIMSSMBC_convergence_tester<Two_scale>::~THDMIIMSSMBC_convergence_tester()
{
}

double THDMIIMSSMBC_convergence_tester<Two_scale>::max_rel_diff() const
{
   const THDMIIMSSMBC<Two_scale>& ol = get_last_iteration_model();
   const THDMIIMSSMBC<Two_scale>& ne = get_model();

   double diff[17] = { 0 };

   diff[0] = MaxRelDiff(OLD(MVZ),NEW(MVZ));
   for (unsigned i = 0; i < 2; i++) {
      diff[i + 1] = MaxRelDiff(OLD1(Mhh,i),NEW1(Mhh,i));
   }
   for (unsigned i = 1; i < 2; i++) {
      diff[i + 3] = MaxRelDiff(OLD1(MAh,i),NEW1(MAh,i));
   }
   for (unsigned i = 1; i < 2; i++) {
      diff[i + 5] = MaxRelDiff(OLD1(MHm,i),NEW1(MHm,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 7] = MaxRelDiff(OLD1(MFd,i),NEW1(MFd,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 10] = MaxRelDiff(OLD1(MFu,i),NEW1(MFu,i));
   }
   for (unsigned i = 0; i < 3; i++) {
      diff[i + 13] = MaxRelDiff(OLD1(MFe,i),NEW1(MFe,i));
   }
   diff[16] = MaxRelDiff(OLD(MVWm),NEW(MVWm));

   return *std::max_element(diff, diff + 17);

}

} // namespace flexiblesusy
