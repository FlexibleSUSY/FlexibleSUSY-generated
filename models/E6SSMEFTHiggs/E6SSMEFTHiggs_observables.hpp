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

// File generated at Tue 10 Oct 2017 22:34:48

#ifndef E6SSMEFTHiggs_OBSERVABLES_H
#define E6SSMEFTHiggs_OBSERVABLES_H

#include <string>
#include <vector>
#include <Eigen/Core>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class E6SSMEFTHiggs_mass_eigenstates;
class Physical_input;

struct E6SSMEFTHiggs_observables {
   static const int NUMBER_OF_OBSERVABLES = 17;

   E6SSMEFTHiggs_observables();
   Eigen::ArrayXd get() const; ///< returns vector of all observables
   static std::vector<std::string> get_names(); ///< returns vector of all observable names
   void clear(); ///< sets all observables to zero
   void set(const Eigen::ArrayXd&); ///< sets all observables from given vector

   double a_muon; ///< a_muon = (g-2)/2 of the muon (calculated with FlexibleSUSY)
   Eigen::Array<std::complex<double>,3,1> eff_cp_higgs_photon_photon; ///< effective H-Photon-Photon coupling
   Eigen::Array<std::complex<double>,3,1> eff_cp_higgs_gluon_gluon; ///< effective H-Gluon-Gluon coupling
   std::complex<double> eff_cp_pseudoscalar_photon_photon; ///< effective A-Photon-Photon coupling
   std::complex<double> eff_cp_pseudoscalar_gluon_gluon; ///< effective A-Gluon-Gluon coupling

};

E6SSMEFTHiggs_observables calculate_observables(
   const E6SSMEFTHiggs_mass_eigenstates&, const softsusy::QedQcd&,
   const Physical_input&);

E6SSMEFTHiggs_observables calculate_observables(
   const E6SSMEFTHiggs_mass_eigenstates&, const softsusy::QedQcd&,
   const Physical_input&, double scale);

} // namespace flexiblesusy

#endif