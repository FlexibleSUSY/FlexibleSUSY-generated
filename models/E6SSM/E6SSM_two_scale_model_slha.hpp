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
 * @file E6SSM_two_scale_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

// File generated at Mon 23 Feb 2015 13:30:30

#ifndef E6SSM_TWO_SCALE_SLHA_H
#define E6SSM_TWO_SCALE_SLHA_H

#include "E6SSM_two_scale_model.hpp"
#include "E6SSM_physical.hpp"
#include "E6SSM_model_slha.hpp"

namespace flexiblesusy {

class Two_scale;

/**
 * @class E6SSM_slha<Two_scale>
 * @brief model class wrapper in SLHA convention
 */

template<>
class E6SSM_slha<Two_scale> : public E6SSM<Two_scale> {
public:
   explicit E6SSM_slha(const E6SSM_input_parameters& input_ = E6SSM_input_parameters());
   explicit E6SSM_slha(const E6SSM<Two_scale>&);
   virtual ~E6SSM_slha();

   virtual void clear();
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   const E6SSM_physical& get_physical_slha() const; ///< returns pole masses to SLHA convention
   E6SSM_physical& get_physical_slha(); ///< returns pole masses to SLHA convention

   // interface functions
   virtual void calculate_spectrum();
   virtual void print(std::ostream&) const;

   double get_MGlu_pole_slha() const { return physical_slha.MGlu; }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return physical_slha.MFv; }
   double get_MFv_pole_slha(int i) const { return physical_slha.MFv(i); }
   double get_MChaP_pole_slha() const { return physical_slha.MChaP; }
   double get_MVZ_pole_slha() const { return physical_slha.MVZ; }
   double get_MVZp_pole_slha() const { return physical_slha.MVZp; }
   const Eigen::Array<double,6,1>& get_MSd_pole_slha() const { return physical_slha.MSd; }
   double get_MSd_pole_slha(int i) const { return physical_slha.MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv_pole_slha() const { return physical_slha.MSv; }
   double get_MSv_pole_slha(int i) const { return physical_slha.MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu_pole_slha() const { return physical_slha.MSu; }
   double get_MSu_pole_slha(int i) const { return physical_slha.MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe_pole_slha() const { return physical_slha.MSe; }
   double get_MSe_pole_slha(int i) const { return physical_slha.MSe(i); }
   const Eigen::Array<double,6,1>& get_MSDX_pole_slha() const { return physical_slha.MSDX; }
   double get_MSDX_pole_slha(int i) const { return physical_slha.MSDX(i); }
   const Eigen::Array<double,3,1>& get_Mhh_pole_slha() const { return physical_slha.Mhh; }
   double get_Mhh_pole_slha(int i) const { return physical_slha.Mhh(i); }
   const Eigen::Array<double,3,1>& get_MAh_pole_slha() const { return physical_slha.MAh; }
   double get_MAh_pole_slha(int i) const { return physical_slha.MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm_pole_slha() const { return physical_slha.MHpm; }
   double get_MHpm_pole_slha(int i) const { return physical_slha.MHpm(i); }
   const Eigen::Array<double,6,1>& get_MChi_pole_slha() const { return physical_slha.MChi; }
   double get_MChi_pole_slha(int i) const { return physical_slha.MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha_pole_slha() const { return physical_slha.MCha; }
   double get_MCha_pole_slha(int i) const { return physical_slha.MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return physical_slha.MFe; }
   double get_MFe_pole_slha(int i) const { return physical_slha.MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return physical_slha.MFd; }
   double get_MFd_pole_slha(int i) const { return physical_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return physical_slha.MFu; }
   double get_MFu_pole_slha(int i) const { return physical_slha.MFu(i); }
   const Eigen::Array<double,3,1>& get_MFDX_pole_slha() const { return physical_slha.MFDX; }
   double get_MFDX_pole_slha(int i) const { return physical_slha.MFDX(i); }
   const Eigen::Array<double,4,1>& get_MSHI0_pole_slha() const { return physical_slha.MSHI0; }
   double get_MSHI0_pole_slha(int i) const { return physical_slha.MSHI0(i); }
   const Eigen::Array<double,4,1>& get_MSHIp_pole_slha() const { return physical_slha.MSHIp; }
   double get_MSHIp_pole_slha(int i) const { return physical_slha.MSHIp(i); }
   const Eigen::Array<double,2,1>& get_MChaI_pole_slha() const { return physical_slha.MChaI; }
   double get_MChaI_pole_slha(int i) const { return physical_slha.MChaI(i); }
   const Eigen::Array<double,4,1>& get_MChiI_pole_slha() const { return physical_slha.MChiI; }
   double get_MChiI_pole_slha(int i) const { return physical_slha.MChiI(i); }
   const Eigen::Array<double,2,1>& get_MSSI0_pole_slha() const { return physical_slha.MSSI0; }
   double get_MSSI0_pole_slha(int i) const { return physical_slha.MSSI0(i); }
   const Eigen::Array<double,2,1>& get_MFSI_pole_slha() const { return physical_slha.MFSI; }
   double get_MFSI_pole_slha(int i) const { return physical_slha.MFSI(i); }
   const Eigen::Array<double,2,1>& get_MSHp0_pole_slha() const { return physical_slha.MSHp0; }
   double get_MSHp0_pole_slha(int i) const { return physical_slha.MSHp0(i); }
   const Eigen::Array<double,2,1>& get_MSHpp_pole_slha() const { return physical_slha.MSHpp; }
   double get_MSHpp_pole_slha(int i) const { return physical_slha.MSHpp(i); }
   const Eigen::Array<double,2,1>& get_MChiP_pole_slha() const { return physical_slha.MChiP; }
   double get_MChiP_pole_slha(int i) const { return physical_slha.MChiP(i); }
   double get_MVG_pole_slha() const { return physical_slha.MVG; }
   double get_MVP_pole_slha() const { return physical_slha.MVP; }
   double get_MVWm_pole_slha() const { return physical_slha.MVWm; }

   const Eigen::Matrix<double,6,6>& get_ZD_pole_slha() const { return physical_slha.ZD; }
   double get_ZD_pole_slha(int i, int k) const { return physical_slha.ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV_pole_slha() const { return physical_slha.ZV; }
   double get_ZV_pole_slha(int i, int k) const { return physical_slha.ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU_pole_slha() const { return physical_slha.ZU; }
   double get_ZU_pole_slha(int i, int k) const { return physical_slha.ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE_pole_slha() const { return physical_slha.ZE; }
   double get_ZE_pole_slha(int i, int k) const { return physical_slha.ZE(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZDX_pole_slha() const { return physical_slha.ZDX; }
   double get_ZDX_pole_slha(int i, int k) const { return physical_slha.ZDX(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZH_pole_slha() const { return physical_slha.ZH; }
   double get_ZH_pole_slha(int i, int k) const { return physical_slha.ZH(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZA_pole_slha() const { return physical_slha.ZA; }
   double get_ZA_pole_slha(int i, int k) const { return physical_slha.ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP_pole_slha() const { return physical_slha.ZP; }
   double get_ZP_pole_slha(int i, int k) const { return physical_slha.ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,6,6>& get_ZN_pole_slha() const { return physical_slha.ZN; }
   const std::complex<double>& get_ZN_pole_slha(int i, int k) const { return physical_slha.ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM_pole_slha() const { return physical_slha.UM; }
   const std::complex<double>& get_UM_pole_slha(int i, int k) const { return physical_slha.UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP_pole_slha() const { return physical_slha.UP; }
   const std::complex<double>& get_UP_pole_slha(int i, int k) const { return physical_slha.UP(i,k); }
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
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXL_pole_slha() const { return physical_slha.ZDXL; }
   const std::complex<double>& get_ZDXL_pole_slha(int i, int k) const { return physical_slha.ZDXL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXR_pole_slha() const { return physical_slha.ZDXR; }
   const std::complex<double>& get_ZDXR_pole_slha(int i, int k) const { return physical_slha.ZDXR(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHI0_pole_slha() const { return physical_slha.UHI0; }
   double get_UHI0_pole_slha(int i, int k) const { return physical_slha.UHI0(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHIp_pole_slha() const { return physical_slha.UHIp; }
   double get_UHIp_pole_slha(int i, int k) const { return physical_slha.UHIp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZMI_pole_slha() const { return physical_slha.ZMI; }
   const std::complex<double>& get_ZMI_pole_slha(int i, int k) const { return physical_slha.ZMI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZPI_pole_slha() const { return physical_slha.ZPI; }
   const std::complex<double>& get_ZPI_pole_slha(int i, int k) const { return physical_slha.ZPI(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZNI_pole_slha() const { return physical_slha.ZNI; }
   const std::complex<double>& get_ZNI_pole_slha(int i, int k) const { return physical_slha.ZNI(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZSSI_pole_slha() const { return physical_slha.ZSSI; }
   double get_ZSSI_pole_slha(int i, int k) const { return physical_slha.ZSSI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZFSI_pole_slha() const { return physical_slha.ZFSI; }
   const std::complex<double>& get_ZFSI_pole_slha(int i, int k) const { return physical_slha.ZFSI(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHp0_pole_slha() const { return physical_slha.UHp0; }
   double get_UHp0_pole_slha(int i, int k) const { return physical_slha.UHp0(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHpp_pole_slha() const { return physical_slha.UHpp; }
   double get_UHpp_pole_slha(int i, int k) const { return physical_slha.UHpp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZNp_pole_slha() const { return physical_slha.ZNp; }
   const std::complex<double>& get_ZNp_pole_slha(int i, int k) const { return physical_slha.ZNp(i,k); }


private:
   E6SSM_physical physical_slha; ///< contains the pole masses and mixings in slha convention
};

} // namespace flexiblesusy

#endif
