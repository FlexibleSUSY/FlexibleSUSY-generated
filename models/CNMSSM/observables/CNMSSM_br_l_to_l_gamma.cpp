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


/**
 * @file CNMSSM_br_l_to_l_gamma.cpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include <valarray>
#include <complex>

#include "CNMSSM_br_l_to_l_gamma.hpp"
#include "CNMSSM_mass_eigenstates.hpp"

#include "cxx_qft/CNMSSM_qft.hpp"
#include "CNMSSM_FFV_form_factors.hpp"

#include "lowe.h"
#include "physical_input.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

namespace {

double get_MSUSY(const CNMSSM_mass_eigenstates_interface& model)
{
   return Min(model.get_MSd().tail<6>().minCoeff(), model.get_MSu().tail<6>().
      minCoeff(), model.get_MSe().tail<6>().minCoeff(), model.get_MHpm().tail<1>()
      .minCoeff(), model.get_MCha().tail<2>().minCoeff());
}

void run_to_MSUSY(CNMSSM_mass_eigenstates& model)
{
   const double precision_goal = model.get_precision();
   const int Nmax = static_cast<int>(-std::log10(precision_goal) * 10);
   int it = 0;
   double precision = 0.;
   double MSUSY_old = 0., MSUSY_new = 0.;

   VERBOSE_MSG("CNMSSM_l_to_lgamma: iterate to run to MSUSY ...");

   do {
      MSUSY_new = get_MSUSY(model);

      if (MSUSY_new <= 0.)
         return;

      VERBOSE_MSG("CNMSSM_l_to_lgamma:    running to new MSUSY = " << MSUSY_new);

      model.run_to(MSUSY_new);
      model.calculate_DRbar_masses();

      precision = MaxRelDiff(MSUSY_old, MSUSY_new);
      MSUSY_old = MSUSY_new;
   } while (precision > precision_goal && ++it < Nmax);

   VERBOSE_MSG("CNMSSM_l_to_lgamma: iteration finished. MSUSY = " << MSUSY_new);

   if (precision > precision_goal) {
      WARNING("CNMSSM_l_to_lgamma: run_to_MSUSY(): "
              "Iteration did not converge after " << Nmax
              << " iterations (precision goal: " << precision_goal << ")");
   }
}

/**
 * @param[in] nI     Generation index of incoming lepton.
 * @return Total decay width in GeV, according to PDG data.
 */
double get_total_width(int nI) {
   // https://pdg.lbl.gov/2020/reviews/rpp2020-rev-phys-constants.pdf
   constexpr double hbar = 6.582119569e-25, // [GeV*s]
   // https://pdg.lbl.gov/2020/listings/rpp2020-list-muon.pdf
                    muon = 2.1969811e-6,    // [s]
   // https://pdg.lbl.gov/2020/listings/rpp2020-list-tau.pdf
                    tau  = 290.3e-15;       // [s]
   switch (nI) {
      case 1: return hbar/muon;
      case 2: return hbar/tau;
      default: throw std::invalid_argument("Unrecognized lepton");
   }
}

} // anonymous namespace

using namespace CNMSSM_cxx_diagrams;
using namespace CNMSSM_cxx_diagrams::fields;
using namespace CNMSSM_FFV_form_factors;

namespace CNMSSM_br_l_to_l_gamma {

template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const& indices1, T2 const& indices2,
      const CNMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd) {
   return  get_total_width(indices1[0]);
}



}
} // namespace flexiblesusy
