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
 * @file MRSSM_two_scale_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

// File generated at Tue 24 Feb 2015 17:34:04

#ifndef MRSSM_TWO_SCALE_SLHA_H
#define MRSSM_TWO_SCALE_SLHA_H

#include "MRSSM_two_scale_model.hpp"
#include "MRSSM_physical.hpp"
#include "MRSSM_model_slha.hpp"

namespace flexiblesusy {

class Two_scale;

/**
 * @class MRSSM_slha<Two_scale>
 * @brief model class wrapper in SLHA convention
 */

template<>
class MRSSM_slha<Two_scale> : public MRSSM<Two_scale> {
public:
   explicit MRSSM_slha(const MRSSM_input_parameters& input_ = MRSSM_input_parameters());
   explicit MRSSM_slha(const MRSSM<Two_scale>&);
   virtual ~MRSSM_slha();

   virtual void clear();
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   const MRSSM_physical& get_physical_slha() const; ///< returns pole masses to SLHA convention
   MRSSM_physical& get_physical_slha(); ///< returns pole masses to SLHA convention

   // interface functions
   virtual void calculate_spectrum();
   virtual void print(std::ostream&) const;

   double get_MVG_pole_slha() const { return physical_slha.MVG; }
   double get_MGlu_pole_slha() const { return physical_slha.MGlu; }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return physical_slha.MFv; }
   double get_MFv_pole_slha(int i) const { return physical_slha.MFv(i); }
   double get_MSOc_pole_slha() const { return physical_slha.MSOc; }
   double get_MVP_pole_slha() const { return physical_slha.MVP; }
   double get_MVZ_pole_slha() const { return physical_slha.MVZ; }
   const Eigen::Array<double,6,1>& get_MSd_pole_slha() const { return physical_slha.MSd; }
   double get_MSd_pole_slha(int i) const { return physical_slha.MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv_pole_slha() const { return physical_slha.MSv; }
   double get_MSv_pole_slha(int i) const { return physical_slha.MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu_pole_slha() const { return physical_slha.MSu; }
   double get_MSu_pole_slha(int i) const { return physical_slha.MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe_pole_slha() const { return physical_slha.MSe; }
   double get_MSe_pole_slha(int i) const { return physical_slha.MSe(i); }
   const Eigen::Array<double,4,1>& get_Mhh_pole_slha() const { return physical_slha.Mhh; }
   double get_Mhh_pole_slha(int i) const { return physical_slha.Mhh(i); }
   const Eigen::Array<double,4,1>& get_MAh_pole_slha() const { return physical_slha.MAh; }
   double get_MAh_pole_slha(int i) const { return physical_slha.MAh(i); }
   const Eigen::Array<double,2,1>& get_MRh_pole_slha() const { return physical_slha.MRh; }
   double get_MRh_pole_slha(int i) const { return physical_slha.MRh(i); }
   const Eigen::Array<double,4,1>& get_MHpm_pole_slha() const { return physical_slha.MHpm; }
   double get_MHpm_pole_slha(int i) const { return physical_slha.MHpm(i); }
   const Eigen::Array<double,2,1>& get_MRpm_pole_slha() const { return physical_slha.MRpm; }
   double get_MRpm_pole_slha(int i) const { return physical_slha.MRpm(i); }
   const Eigen::Array<double,4,1>& get_MChi_pole_slha() const { return physical_slha.MChi; }
   double get_MChi_pole_slha(int i) const { return physical_slha.MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha1_pole_slha() const { return physical_slha.MCha1; }
   double get_MCha1_pole_slha(int i) const { return physical_slha.MCha1(i); }
   const Eigen::Array<double,2,1>& get_MCha2_pole_slha() const { return physical_slha.MCha2; }
   double get_MCha2_pole_slha(int i) const { return physical_slha.MCha2(i); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return physical_slha.MFe; }
   double get_MFe_pole_slha(int i) const { return physical_slha.MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return physical_slha.MFd; }
   double get_MFd_pole_slha(int i) const { return physical_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return physical_slha.MFu; }
   double get_MFu_pole_slha(int i) const { return physical_slha.MFu(i); }
   double get_MVWm_pole_slha() const { return physical_slha.MVWm; }

   const Eigen::Matrix<double,6,6>& get_ZD_pole_slha() const { return physical_slha.ZD; }
   double get_ZD_pole_slha(int i, int k) const { return physical_slha.ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV_pole_slha() const { return physical_slha.ZV; }
   double get_ZV_pole_slha(int i, int k) const { return physical_slha.ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU_pole_slha() const { return physical_slha.ZU; }
   double get_ZU_pole_slha(int i, int k) const { return physical_slha.ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE_pole_slha() const { return physical_slha.ZE; }
   double get_ZE_pole_slha(int i, int k) const { return physical_slha.ZE(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZH_pole_slha() const { return physical_slha.ZH; }
   double get_ZH_pole_slha(int i, int k) const { return physical_slha.ZH(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZA_pole_slha() const { return physical_slha.ZA; }
   double get_ZA_pole_slha(int i, int k) const { return physical_slha.ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZHR_pole_slha() const { return physical_slha.ZHR; }
   double get_ZHR_pole_slha(int i, int k) const { return physical_slha.ZHR(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZP_pole_slha() const { return physical_slha.ZP; }
   double get_ZP_pole_slha(int i, int k) const { return physical_slha.ZP(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZRP_pole_slha() const { return physical_slha.ZRP; }
   double get_ZRP_pole_slha(int i, int k) const { return physical_slha.ZRP(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN1_pole_slha() const { return physical_slha.ZN1; }
   const std::complex<double>& get_ZN1_pole_slha(int i, int k) const { return physical_slha.ZN1(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN2_pole_slha() const { return physical_slha.ZN2; }
   const std::complex<double>& get_ZN2_pole_slha(int i, int k) const { return physical_slha.ZN2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM1_pole_slha() const { return physical_slha.UM1; }
   const std::complex<double>& get_UM1_pole_slha(int i, int k) const { return physical_slha.UM1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP1_pole_slha() const { return physical_slha.UP1; }
   const std::complex<double>& get_UP1_pole_slha(int i, int k) const { return physical_slha.UP1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM2_pole_slha() const { return physical_slha.UM2; }
   const std::complex<double>& get_UM2_pole_slha(int i, int k) const { return physical_slha.UM2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP2_pole_slha() const { return physical_slha.UP2; }
   const std::complex<double>& get_UP2_pole_slha(int i, int k) const { return physical_slha.UP2(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL_pole_slha() const { return physical_slha.ZEL; }
   const std::complex<double>& get_ZEL_pole_slha(int i, int k) const { return physical_slha.ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER_pole_slha() const { return physical_slha.ZER; }
   const std::complex<double>& get_ZER_pole_slha(int i, int k) const { return physical_slha.ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL_pole_slha() const { return physical_slha.ZDL; }
   const std::complex<double>& get_ZDL_pole_slha(int i, int k) const { return physical_slha.ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR_pole_slha() const { return physical_slha.ZDR; }
   const std::complex<double>& get_ZDR_pole_slha(int i, int k) const { return physical_slha.ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL_pole_slha() const { return physical_slha.ZUL; }
   const std::complex<double>& get_ZUL_pole_slha(int i, int k) const { return physical_slha.ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR_pole_slha() const { return physical_slha.ZUR; }
   const std::complex<double>& get_ZUR_pole_slha(int i, int k) const { return physical_slha.ZUR(i,k); }


private:
   MRSSM_physical physical_slha; ///< contains the pole masses and mixings in slha convention
};

} // namespace flexiblesusy

#endif
