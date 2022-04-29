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
 * @file MSSMNoFV_b_to_s_gamma.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#include <array>
#include <complex>
#include <fstream>
#include <iomanip>

#include "MSSMNoFV_b_to_s_gamma.hpp"
#include "MSSMNoFV_mass_eigenstates.hpp"
#include "cxx_qft/MSSMNoFV_qft.hpp"
#include "MSSMNoFV_FFV_form_factors.hpp"

#include "lowe.h"
#include "wrappers.hpp"
#include "json.hpp"


#define MODELPARAMETER(p) context.model.get_##p()

namespace flexiblesusy {

using namespace MSSMNoFV_cxx_diagrams;
using namespace MSSMNoFV_cxx_diagrams::fields;
using namespace MSSMNoFV_FFV_form_factors;

namespace MSSMNoFV_b_to_s_gamma {

/**
 * @brief Writes the Wilson coefficients into a json file to make it usable
 *        for the wilson and flavio package
 */
void write_wilsoncoeffs(const std::complex<double>& C7NP_bs, const std::complex<double>& C7pNP_bs,
   const std::complex<double>& C8NP_bs, const std::complex<double>& C8pNP_bs,
   const double& matching_scale)
{
   nlohmann::json j;
   j["eft"] = "WET";
   j["basis"] = "flavio";
   j["scale"] = matching_scale;
   j["values"] = {
      {"C7_bs",  {{"Re", Re(C7NP_bs)},  {"Im", Im(C7NP_bs)}}},
      {"C7p_bs", {{"Re", Re(C7pNP_bs)}, {"Im", Im(C7pNP_bs)}}},
      {"C8_bs",  {{"Re", Re(C8NP_bs)},  {"Im", Im(C8NP_bs)}}},
      {"C8p_bs", {{"Re", Re(C8pNP_bs)}, {"Im", Im(C8pNP_bs)}}}
   };
   std::ofstream wc_json("WC_MSSMNoFV.json");
   wc_json << std::setw(4) << j << std::endl;
   wc_json.close();
}

/**
 * @brief Calculates the Wilson coefficients C7, C7p, C8 and C8p
 *        which are required to compute the branching ratio
 *        BR(B->Xs gamma)
 */
std::array<std::complex<double>, 4> calculate_b_to_s_gamma(
   const MSSMNoFV_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   context_base context { model };

   const auto complex_CKM = qedqcd.get_complex_ckm();
   const std::complex<double> Vtb (complex_CKM(2, 2));
   std::complex<double> Vts;
   // avoid division by zero
   if (is_zero(std::abs(complex_CKM(2, 1)))) {
       Vts = std::complex<double>(-0.0404063888996, -0.00072107415577);
   } else {
       Vts = complex_CKM(2, 1);
   }
   const auto g2 = MODELPARAMETER(g2);
   const auto g3 = MODELPARAMETER(g3);

   return {0, 0, 0, 0};
}

} // namespace MSSMNoFV_b_to_s_gamma
} // namespace flexiblesusy
