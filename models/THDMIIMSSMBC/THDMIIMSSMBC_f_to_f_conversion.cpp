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
 * @file THDMIIMSSMBC_f_to_f_conversion.cpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#include <valarray>
#include <complex>

#include "THDMIIMSSMBC_mass_eigenstates.hpp"
#include "cxx_qft/THDMIIMSSMBC_qft.hpp"

#include "THDMIIMSSMBC_f_to_f_conversion.hpp"
#include "THDMIIMSSMBC_FFV_form_factors.hpp"

#include "lowe.h"
#include "wrappers.hpp"



using namespace flexiblesusy;
using namespace THDMIIMSSMBC_cxx_diagrams;
using namespace THDMIIMSSMBC_FFV_form_factors;

namespace flexiblesusy {
using namespace THDMIIMSSMBC_cxx_diagrams::fields;

namespace THDMIIMSSMBC_cxx_diagrams
{
namespace npointfunctions
{

} // namespace npointfunctions
} // namespace THDMIIMSSMBC_cxx_diagrams
namespace THDMIIMSSMBC_f_to_f_conversion {

struct overlap_integrals {
   double D;
   double Vp;
   double Vn;
   double Sp;
   double Sn;
};

overlap_integrals get_overlap_integrals(const Nucleus N, const softsusy::QedQcd& qedqcd) {
   overlap_integrals res;

   // get muon pole mass from input slha file
   const auto muon_pole_mass_5o2 = pow(qedqcd.displayMass(softsusy::mMuon), 5./2.);

   // Tab. 2 of hep-ph/0203110
   switch (N) {
      case Nucleus::Au:
         res.D  = 0.1670 * muon_pole_mass_5o2;
         res.Vp = 0.0859 * muon_pole_mass_5o2;
         res.Vn = 0.1080 * muon_pole_mass_5o2;
         res.Sp = 0.0523 * muon_pole_mass_5o2;
         res.Sn = 0.0610 * muon_pole_mass_5o2;
         break;
      case Nucleus::Al:
         res.D  = 0.0357 * muon_pole_mass_5o2;
         res.Vp = 0.0159 * muon_pole_mass_5o2;
         res.Vn = 0.0169 * muon_pole_mass_5o2;
         res.Sp = 0.0153 * muon_pole_mass_5o2;
         res.Sn = 0.0163 * muon_pole_mass_5o2;
         break;
      case Nucleus::Ti:
         res.D  = 0. * muon_pole_mass_5o2;
         res.Vp = 0. * muon_pole_mass_5o2;
         res.Vn = 0. * muon_pole_mass_5o2;
         res.Sp = 0. * muon_pole_mass_5o2;
         res.Sn = 0. * muon_pole_mass_5o2;
         break;
      default:
         throw std::invalid_argument("Unknown nucleus");
   }
   
   return res;
}

double get_capture_rate (const Nucleus N) {
   switch (N) {
      case Nucleus::Au:
         return 8.84868e-18;
      case Nucleus::Al:
         return 4.64079e-19;
      case Nucleus::Ti:
	 return 0.;
      default:
         throw std::invalid_argument("Unknown nucleus");
   }
}

template <class FOutBar, class FIn, class VBar>
std::complex<double> vectorCurrent(const THDMIIMSSMBC_mass_eigenstates& model) {
    context_base context {model};
    using vertex = Vertex<FOutBar, FIn, VBar>;
    std::array<int, 2> indices {0, 0};
    const auto value =  vertex::evaluate(indices, context);
    return 0.5*(value.left() + value.right());
}


} // namespace THDMIIMSSMBC_f_to_f_conversion
} // namespace flexiblesusy
