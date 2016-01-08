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

// File generated at Fri 8 Jan 2016 15:31:35

#include "MSSMNoFV_observables.hpp"
#include "MSSMNoFV_mass_eigenstates.hpp"
#include "gm2calc_interface.hpp"
#include "eigen_utils.hpp"
#include "numerics2.hpp"
#include "lowe.h"

#define MODEL model
#define AMUGM2CALC a_muon_gm2calc
#define AMUGM2CALCUNCERTAINTY a_muon_gm2calc_uncertainty

#define ALPHA_S_MZ qedqcd.displayAlpha(softsusy::ALPHAS)
#define MWPole qedqcd.displayPoleMW()
#define MZPole qedqcd.displayPoleMZ()
#define MTPole qedqcd.displayPoleMt()
#define MBMB qedqcd.displayMbMb()
#define MTauPole qedqcd.displayPoleMtau()
#define MMPole qedqcd.displayPoleMmuon()

namespace flexiblesusy {

Observables calculate_observables(const MSSMNoFV_mass_eigenstates& model,
                                  const softsusy::QedQcd& qedqcd)
{
   Observables observables;

   GM2Calc_data gm2calc_data;
   gm2calc_data.alpha_s_MZ = ALPHA_S_MZ;
   gm2calc_data.MZ    = MZPole;
   if (!is_zero(MODEL.get_physical().MVWm))
      gm2calc_data.MW = MODEL.get_physical().MVWm;
   else if (!is_zero(MWPole))
      gm2calc_data.MW = MWPole;
   gm2calc_data.mb_mb = MBMB;
   gm2calc_data.MT    = MTPole;
   gm2calc_data.MTau  = MTauPole;
   gm2calc_data.MM    = MMPole;
   gm2calc_data.MA0   = MODEL.get_physical().MAh(1);
   gm2calc_data.MSvm  = MODEL.get_physical().MSvmL;
   gm2calc_data.MSm   = MODEL.get_physical().MSm;
   gm2calc_data.MCha  = MODEL.get_physical().MCha;
   gm2calc_data.MChi  = MODEL.get_physical().MChi;
   gm2calc_data.scale = MODEL.get_scale();
   gm2calc_data.TB    = MODEL.get_vu() / MODEL.get_vd();
   gm2calc_data.Mu    = MODEL.get_Mu();
   gm2calc_data.M1    = MODEL.get_MassB();
   gm2calc_data.M2    = MODEL.get_MassWB();
   gm2calc_data.M3    = MODEL.get_MassG();
   gm2calc_data.mq2   = MODEL.get_mq2();
   gm2calc_data.mu2   = MODEL.get_mu2();
   gm2calc_data.md2   = MODEL.get_md2();
   gm2calc_data.ml2   = MODEL.get_ml2();
   gm2calc_data.me2   = MODEL.get_me2();
   gm2calc_data.Au    = div_save(MODEL.get_TYu(), MODEL.get_Yu());
   gm2calc_data.Ad    = div_save(MODEL.get_TYd(), MODEL.get_Yd());
   gm2calc_data.Ae    = div_save(MODEL.get_TYe(), MODEL.get_Ye());

   observables.AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);
   observables.AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);

   return observables;
}

} // namespace flexiblesusy
