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

// File generated at Mon 19 Sep 2016 09:42:23

#ifndef MRSSMtower_EFFECTIVE_COUPLINGS_H
#define MRSSMtower_EFFECTIVE_COUPLINGS_H

#include "MRSSMtower_mass_eigenstates.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

class MRSSMtower_effective_couplings {
public:
   MRSSMtower_effective_couplings(const MRSSMtower_mass_eigenstates&,
                                   const softsusy::QedQcd&,
                                   const Physical_input&);

   void do_run_couplings(bool flag) { rg_improve = flag; }
   bool do_run_couplings() const { return rg_improve; }
   void do_include_qcd_corrections(bool flag) { include_qcd_corrections = flag; }
   bool do_include_qcd_corrections() const { return include_qcd_corrections; }
   void set_physical_inputs(const Physical_input& inputs_) { physical_input = inputs_; }
   void set_low_energy_data(const softsusy::QedQcd& qedqcd_) { qedqcd = qedqcd_; }
   void set_model(const MRSSMtower_mass_eigenstates& model_);

   double get_hhVPVP_partial_width(unsigned gO1) const;
   double get_hhVGVG_partial_width(unsigned gO1) const;
   double get_AhVPVP_partial_width(unsigned gO1) const;
   double get_AhVGVG_partial_width(unsigned gO1) const;
   std::complex<double> get_eff_CphhVPVP(unsigned gO1) const { return eff_CphhVPVP(gO1); }
   std::complex<double> get_eff_CphhVGVG(unsigned gO1) const { return eff_CphhVGVG(gO1); }
   std::complex<double> get_eff_CpAhVPVP(unsigned gO1) const { return eff_CpAhVPVP(gO1); }
   std::complex<double> get_eff_CpAhVGVG(unsigned gO1) const { return eff_CpAhVGVG(gO1); }

   void calculate_effective_couplings();

   std::complex<double> CphhSRdpconjSRdp(unsigned gt1) const;
   std::complex<double> CphhSRumconjSRum(unsigned gt1) const;
   std::complex<double> CphhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CphhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CphhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CphhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpCha1hhbarCha1PL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpCha2hhbarCha2PL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpFehhbarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpFdhhbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpFuhhbarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CphhVWmconjVWm(unsigned gt1) const;
   std::complex<double> CpAhSRdpconjSRdp(unsigned gt1) const;
   std::complex<double> CpAhSRumconjSRum(unsigned gt1) const;
   std::complex<double> CpAhSdconjSd(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhSuconjSu(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhSeconjSe(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhHpmconjHpm(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhCha1barCha1PL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhCha2barCha2PL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhFebarFePL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhFdbarFdPL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhFubarFuPL(unsigned gt1, unsigned gt2, unsigned gt3) const;
   void calculate_eff_CphhVPVP(unsigned gO1);
   void calculate_eff_CphhVGVG(unsigned gO1);
   void calculate_eff_CpAhVPVP(unsigned gO1);
   void calculate_eff_CpAhVGVG(unsigned gO1);

private:
   MRSSMtower_mass_eigenstates model;
   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   bool rg_improve;
   bool include_qcd_corrections;

   void copy_mixing_matrices_from_model();

   void run_SM_strong_coupling_to(double m);

   // higher order corrections to the amplitudes for
   // effective coupling to photons
   std::complex<double> scalar_scalar_qcd_factor(double, double) const;
   std::complex<double> scalar_fermion_qcd_factor(double, double) const;
   std::complex<double> pseudoscalar_fermion_qcd_factor(double, double) const;

   // higher order corrections to the leading order
   // effective couplings to gluons
   double number_of_active_flavours(double) const;
   double scalar_scaling_factor(double) const;
   double pseudoscalar_scaling_factor(double) const;

   Eigen::Matrix<double,6,6> ZD;
   Eigen::Matrix<double,3,3> ZV;
   Eigen::Matrix<double,6,6> ZU;
   Eigen::Matrix<double,6,6> ZE;
   Eigen::Matrix<double,4,4> ZH;
   Eigen::Matrix<double,4,4> ZA;
   Eigen::Matrix<double,2,2> ZHR;
   Eigen::Matrix<double,4,4> ZP;
   Eigen::Matrix<std::complex<double>,4,4> ZN1;
   Eigen::Matrix<std::complex<double>,4,4> ZN2;
   Eigen::Matrix<std::complex<double>,2,2> UM1;
   Eigen::Matrix<std::complex<double>,2,2> UP1;
   Eigen::Matrix<std::complex<double>,2,2> UM2;
   Eigen::Matrix<std::complex<double>,2,2> UP2;
   Eigen::Matrix<std::complex<double>,3,3> ZEL;
   Eigen::Matrix<std::complex<double>,3,3> ZER;
   Eigen::Matrix<std::complex<double>,3,3> ZDL;
   Eigen::Matrix<std::complex<double>,3,3> ZDR;
   Eigen::Matrix<std::complex<double>,3,3> ZUL;
   Eigen::Matrix<std::complex<double>,3,3> ZUR;
   Eigen::Matrix<double,2,2> ZZ;

   Eigen::Array<std::complex<double>,4,1> eff_CphhVPVP;
   Eigen::Array<std::complex<double>,4,1> eff_CphhVGVG;
   Eigen::Array<std::complex<double>,4,1> eff_CpAhVPVP;
   Eigen::Array<std::complex<double>,4,1> eff_CpAhVGVG;

};

} // namespace flexiblesusy

#endif