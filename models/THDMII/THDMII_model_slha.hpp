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
 * @file THDMII_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */


#ifndef THDMII_SLHA_H
#define THDMII_SLHA_H

#include "THDMII_input_parameters.hpp"
#include "THDMII_mass_eigenstates.hpp"
#include "THDMII_physical.hpp"
#include "wrappers.hpp"

#define LOCALPHYSICAL(p) physical.p
#define MODELPARAMETER(p) this->p
#define PHYSICAL_SLHA(p) physical_slha.p
#define PHYSICAL_SLHA_REAL(p) Re(physical_slha.p)

namespace flexiblesusy {

/**
 * @class THDMII_slha<T>
 * @brief model class wrapper in SLHA convention
 *
 * @tparam Model model class to wrap
 */

class THDMII_slha : public THDMII_mass_eigenstates {
public:
   explicit THDMII_slha(const THDMII_input_parameters& input_ = THDMII_input_parameters(), bool do_convert_masses_to_slha = true);
   explicit THDMII_slha(const THDMII_mass_eigenstates&, bool do_convert_masses_to_slha = true);
   THDMII_slha(const THDMII_slha&) = default;
   THDMII_slha(THDMII_slha&&) = default;
   virtual ~THDMII_slha() = default;
   THDMII_slha& operator=(const THDMII_slha&) = default;
   THDMII_slha& operator=(THDMII_slha&&) = default;

   virtual void clear() override;
   void convert_to_slha(); ///< converts pole masses and couplings to SLHA convention
   const Eigen::Matrix<std::complex<double>,3,3>& get_ckm_matrix() const { return ckm; }
   const Eigen::Matrix<std::complex<double>,3,3>& get_pmns_matrix() const { return pmns; }
   const THDMII_physical& get_physical_slha() const; ///< returns pole masses to SLHA convention
   THDMII_physical& get_physical_slha(); ///< returns pole masses to SLHA convention
   void set_convert_masses_to_slha(bool); ///< allow/disallow for negative majoran fermion masses

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void print(std::ostream&) const override;

   double get_MVG_pole_slha() const { return PHYSICAL_SLHA(MVG); }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return PHYSICAL_SLHA(MFv); }
   double get_MFv_pole_slha(int i) const { return PHYSICAL_SLHA(MFv(i)); }
   const Eigen::Array<double,2,1>& get_Mhh_pole_slha() const { return PHYSICAL_SLHA(Mhh); }
   double get_Mhh_pole_slha(int i) const { return PHYSICAL_SLHA(Mhh(i)); }
   const Eigen::Array<double,2,1>& get_MAh_pole_slha() const { return PHYSICAL_SLHA(MAh); }
   double get_MAh_pole_slha(int i) const { return PHYSICAL_SLHA(MAh(i)); }
   const Eigen::Array<double,2,1>& get_MHm_pole_slha() const { return PHYSICAL_SLHA(MHm); }
   double get_MHm_pole_slha(int i) const { return PHYSICAL_SLHA(MHm(i)); }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return PHYSICAL_SLHA(MFd); }
   double get_MFd_pole_slha(int i) const { return PHYSICAL_SLHA(MFd(i)); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return PHYSICAL_SLHA(MFu); }
   double get_MFu_pole_slha(int i) const { return PHYSICAL_SLHA(MFu(i)); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return PHYSICAL_SLHA(MFe); }
   double get_MFe_pole_slha(int i) const { return PHYSICAL_SLHA(MFe(i)); }
   double get_MVWm_pole_slha() const { return PHYSICAL_SLHA(MVWm); }
   double get_MVP_pole_slha() const { return PHYSICAL_SLHA(MVP); }
   double get_MVZ_pole_slha() const { return PHYSICAL_SLHA(MVZ); }

   const Eigen::Matrix<double,2,2>& get_ZH_pole_slha() const { return PHYSICAL_SLHA(ZH); }
   double get_ZH_pole_slha(int i, int k) const { return PHYSICAL_SLHA(ZH(i,k)); }
   const Eigen::Matrix<double,2,2>& get_ZA_pole_slha() const { return PHYSICAL_SLHA(ZA); }
   double get_ZA_pole_slha(int i, int k) const { return PHYSICAL_SLHA(ZA(i,k)); }
   const Eigen::Matrix<double,2,2>& get_ZP_pole_slha() const { return PHYSICAL_SLHA(ZP); }
   double get_ZP_pole_slha(int i, int k) const { return PHYSICAL_SLHA(ZP(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd_pole_slha() const { return PHYSICAL_SLHA(Vd); }
   std::complex<double> get_Vd_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Vd(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud_pole_slha() const { return PHYSICAL_SLHA(Ud); }
   std::complex<double> get_Ud_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Ud(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu_pole_slha() const { return PHYSICAL_SLHA(Vu); }
   std::complex<double> get_Vu_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Vu(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu_pole_slha() const { return PHYSICAL_SLHA(Uu); }
   std::complex<double> get_Uu_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Uu(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve_pole_slha() const { return PHYSICAL_SLHA(Ve); }
   std::complex<double> get_Ve_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Ve(i,k)); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue_pole_slha() const { return PHYSICAL_SLHA(Ue); }
   std::complex<double> get_Ue_pole_slha(int i, int k) const { return PHYSICAL_SLHA(Ue(i,k)); }
   const Eigen::Matrix<double,2,2>& get_ZZ_pole_slha() const { return PHYSICAL_SLHA(ZZ); }
   double get_ZZ_pole_slha(int i, int k) const { return PHYSICAL_SLHA(ZZ(i,k)); }

   const Eigen::Array<double,3,1>& get_Yu_slha() const { return Yu_slha; }
   double get_Yu_slha(int i) const { return Yu_slha(i); }
   const Eigen::Array<double,3,1>& get_Yd_slha() const { return Yd_slha; }
   double get_Yd_slha(int i) const { return Yd_slha(i); }
   const Eigen::Array<double,3,1>& get_Ye_slha() const { return Ye_slha; }
   double get_Ye_slha(int i) const { return Ye_slha(i); }



   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd_slha() const { return Vd_slha; }
   std::complex<double> get_Vd_slha(int i, int k) const { return Vd_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu_slha() const { return Vu_slha; }
   std::complex<double> get_Vu_slha(int i, int k) const { return Vu_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud_slha() const { return Ud_slha; }
   std::complex<double> get_Ud_slha(int i, int k) const { return Ud_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu_slha() const { return Uu_slha; }
   std::complex<double> get_Uu_slha(int i, int k) const { return Uu_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve_slha() const { return Ve_slha; }
   std::complex<double> get_Ve_slha(int i, int k) const { return Ve_slha(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue_slha() const { return Ue_slha; }
   std::complex<double> get_Ue_slha(int i, int k) const { return Ue_slha(i,k); }


private:
   THDMII_physical physical_slha{}; ///< contains the pole masses and mixings in slha convention
   Eigen::Matrix<std::complex<double>,3,3> ckm{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   Eigen::Matrix<std::complex<double>,3,3> pmns{Eigen::Matrix<std::complex<double>,3,3>::Identity()};
   bool convert_masses_to_slha{true};    ///< allow/disallow for negative majoran fermion masses
   Eigen::Array<double,3,1> Yu_slha{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> Yd_slha{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> Ye_slha{Eigen::Array<double,3,1>::Zero()};

   Eigen::Matrix<std::complex<double>,3,3> Vd_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue_slha{Eigen::Matrix<std::complex<double>,3,3>::Zero()};




   void calculate_ckm_matrix();
   void calculate_pmns_matrix();
   void convert_yukawa_couplings_to_slha();
   void convert_trilinear_couplings_to_slha();
   void convert_soft_squared_masses_to_slha();
};

} // namespace flexiblesusy

#undef LOCALPHYSICAL
#undef MODELPARAMETER
#undef PHYSICAL_SLHA
#undef PHYSICAL_SLHA_REAL

#endif
