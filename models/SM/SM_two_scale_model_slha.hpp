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
 * @file SM_two_scale_model_slha.hpp
 * @brief contains wrapper class for model class in SLHA convention
 */

// File generated at Mon 23 Feb 2015 12:25:45

#ifndef SM_TWO_SCALE_SLHA_H
#define SM_TWO_SCALE_SLHA_H

#include "SM_two_scale_model.hpp"
#include "SM_physical.hpp"
#include "SM_model_slha.hpp"

namespace flexiblesusy {

class Two_scale;

/**
 * @class SM_slha<Two_scale>
 * @brief model class wrapper in SLHA convention
 */

template<>
class SM_slha<Two_scale> : public SM<Two_scale> {
public:
   explicit SM_slha(const SM_input_parameters& input_ = SM_input_parameters());
   explicit SM_slha(const SM<Two_scale>&);
   virtual ~SM_slha();

   virtual void clear();
   void convert_to_slha(); ///< converts pole masses to SLHA convention
   const SM_physical& get_physical_slha() const; ///< returns pole masses to SLHA convention
   SM_physical& get_physical_slha(); ///< returns pole masses to SLHA convention

   // interface functions
   virtual void calculate_spectrum();
   virtual void print(std::ostream&) const;

   double get_MHp_pole_slha() const { return physical_slha.MHp; }
   const Eigen::Array<double,3,1>& get_MFv_pole_slha() const { return physical_slha.MFv; }
   double get_MFv_pole_slha(int i) const { return physical_slha.MFv(i); }
   double get_MAh_pole_slha() const { return physical_slha.MAh; }
   double get_Mhh_pole_slha() const { return physical_slha.Mhh; }
   double get_MVZ_pole_slha() const { return physical_slha.MVZ; }
   const Eigen::Array<double,3,1>& get_MFd_pole_slha() const { return physical_slha.MFd; }
   double get_MFd_pole_slha(int i) const { return physical_slha.MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu_pole_slha() const { return physical_slha.MFu; }
   double get_MFu_pole_slha(int i) const { return physical_slha.MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe_pole_slha() const { return physical_slha.MFe; }
   double get_MFe_pole_slha(int i) const { return physical_slha.MFe(i); }
   double get_MVG_pole_slha() const { return physical_slha.MVG; }
   double get_MVP_pole_slha() const { return physical_slha.MVP; }
   double get_MVWp_pole_slha() const { return physical_slha.MVWp; }

   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd_pole_slha() const { return physical_slha.Vd; }
   const std::complex<double>& get_Vd_pole_slha(int i, int k) const { return physical_slha.Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud_pole_slha() const { return physical_slha.Ud; }
   const std::complex<double>& get_Ud_pole_slha(int i, int k) const { return physical_slha.Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu_pole_slha() const { return physical_slha.Vu; }
   const std::complex<double>& get_Vu_pole_slha(int i, int k) const { return physical_slha.Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu_pole_slha() const { return physical_slha.Uu; }
   const std::complex<double>& get_Uu_pole_slha(int i, int k) const { return physical_slha.Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve_pole_slha() const { return physical_slha.Ve; }
   const std::complex<double>& get_Ve_pole_slha(int i, int k) const { return physical_slha.Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue_pole_slha() const { return physical_slha.Ue; }
   const std::complex<double>& get_Ue_pole_slha(int i, int k) const { return physical_slha.Ue(i,k); }


private:
   SM_physical physical_slha; ///< contains the pole masses and mixings in slha convention
};

} // namespace flexiblesusy

#endif
