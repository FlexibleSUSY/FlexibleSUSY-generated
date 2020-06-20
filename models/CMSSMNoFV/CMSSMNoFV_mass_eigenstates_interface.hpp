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
 * @file CMSSMNoFV_mass_eigenstates_interface.hpp
 *
 * @brief Contains the mass eigenstates interface class definition
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#ifndef CMSSMNoFV_MASS_EIGENSTATES_INTERFACE_H
#define CMSSMNoFV_MASS_EIGENSTATES_INTERFACE_H

#include "CMSSMNoFV_physical.hpp"

#include <memory>

#include <Eigen/Core>

namespace flexiblesusy {

class Problems;
struct CMSSMNoFV_input_parameters;

/**
 * @class CMSSMNoFV_mass_eigenstates_interface
 * @brief Interface definition for model parameters, masses and mixings
 */
class CMSSMNoFV_mass_eigenstates_interface {
public:
   virtual ~CMSSMNoFV_mass_eigenstates_interface() {}

   virtual std::unique_ptr<CMSSMNoFV_mass_eigenstates_interface> clone() const = 0;

   virtual void calculate_tree_level_mass_spectrum() = 0;
   virtual void calculate_pole_mass_spectrum() = 0;
   virtual void calculate_mass_spectrum() = 0;

   virtual int solve_ewsb_equations_tree_level() = 0;
   virtual int solve_ewsb_equations() = 0;

   virtual Eigen::ArrayXd get_tree_level_masses() const = 0;
   virtual Eigen::ArrayXd get_tree_level_masses_and_mixings() const = 0;
   virtual const CMSSMNoFV_input_parameters& get_input_parameters() const = 0;
   virtual CMSSMNoFV_input_parameters& get_input_parameters() = 0;
   virtual Eigen::ArrayXd get_extra_parameters() const = 0;
   virtual const CMSSMNoFV_physical& get_physical() const = 0;
   virtual CMSSMNoFV_physical& get_physical() = 0;
   virtual const Problems& get_problems() const = 0;
   virtual Problems& get_problems() = 0;
   virtual void set_tree_level_masses(const Eigen::ArrayXd&) = 0;
   virtual void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) = 0;
   virtual void set_extra_parameters(const Eigen::ArrayXd&) = 0;
   virtual void set_physical(const CMSSMNoFV_physical&) = 0;
   virtual void clear_problems() = 0;

   virtual const Eigen::Matrix<double,3,3>& get_Yd() const = 0;
   virtual double get_Yd(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Ye() const = 0;
   virtual double get_Ye(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yu() const = 0;
   virtual double get_Yu(int i, int k) const = 0;
   virtual double get_Mu() const = 0;
   virtual double get_g1() const = 0;
   virtual double get_g2() const = 0;
   virtual double get_g3() const = 0;
   virtual double get_vd() const = 0;
   virtual double get_vu() const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_TYd() const = 0;
   virtual double get_TYd(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_TYe() const = 0;
   virtual double get_TYe(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_TYu() const = 0;
   virtual double get_TYu(int i, int k) const = 0;
   virtual double get_BMu() const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_mq2() const = 0;
   virtual double get_mq2(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_ml2() const = 0;
   virtual double get_ml2(int i, int k) const = 0;
   virtual double get_mHd2() const = 0;
   virtual double get_mHu2() const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_md2() const = 0;
   virtual double get_md2(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_mu2() const = 0;
   virtual double get_mu2(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_me2() const = 0;
   virtual double get_me2(int i, int k) const = 0;
   virtual double get_MassB() const = 0;
   virtual double get_MassWB() const = 0;
   virtual double get_MassG() const = 0;
   virtual void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) = 0;
   virtual void set_Yd(int i, int k, const double& value) = 0;
   virtual void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) = 0;
   virtual void set_Ye(int i, int k, const double& value) = 0;
   virtual void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) = 0;
   virtual void set_Yu(int i, int k, const double& value) = 0;
   virtual void set_Mu(double Mu_) = 0;
   virtual void set_g1(double g1_) = 0;
   virtual void set_g2(double g2_) = 0;
   virtual void set_g3(double g3_) = 0;
   virtual void set_vd(double vd_) = 0;
   virtual void set_vu(double vu_) = 0;
   virtual void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) = 0;
   virtual void set_TYd(int i, int k, const double& value) = 0;
   virtual void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) = 0;
   virtual void set_TYe(int i, int k, const double& value) = 0;
   virtual void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) = 0;
   virtual void set_TYu(int i, int k, const double& value) = 0;
   virtual void set_BMu(double BMu_) = 0;
   virtual void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) = 0;
   virtual void set_mq2(int i, int k, const double& value) = 0;
   virtual void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) = 0;
   virtual void set_ml2(int i, int k, const double& value) = 0;
   virtual void set_mHd2(double mHd2_) = 0;
   virtual void set_mHu2(double mHu2_) = 0;
   virtual void set_md2(const Eigen::Matrix<double,3,3>& md2_) = 0;
   virtual void set_md2(int i, int k, const double& value) = 0;
   virtual void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) = 0;
   virtual void set_mu2(int i, int k, const double& value) = 0;
   virtual void set_me2(const Eigen::Matrix<double,3,3>& me2_) = 0;
   virtual void set_me2(int i, int k, const double& value) = 0;
   virtual void set_MassB(double MassB_) = 0;
   virtual void set_MassWB(double MassWB_) = 0;
   virtual void set_MassG(double MassG_) = 0;
   virtual double get_MVG() const = 0;
   virtual double get_MGlu() const = 0;
   virtual double get_MFd() const = 0;
   virtual double get_MFs() const = 0;
   virtual double get_MFb() const = 0;
   virtual double get_MFu() const = 0;
   virtual double get_MFc() const = 0;
   virtual double get_MFt() const = 0;
   virtual double get_MFve() const = 0;
   virtual double get_MFvm() const = 0;
   virtual double get_MFvt() const = 0;
   virtual double get_MFe() const = 0;
   virtual double get_MFm() const = 0;
   virtual double get_MFtau() const = 0;
   virtual double get_MSveL() const = 0;
   virtual double get_MSvmL() const = 0;
   virtual double get_MSvtL() const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSd() const = 0;
   virtual double get_MSd(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSu() const = 0;
   virtual double get_MSu(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSe() const = 0;
   virtual double get_MSe(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSm() const = 0;
   virtual double get_MSm(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MStau() const = 0;
   virtual double get_MStau(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSs() const = 0;
   virtual double get_MSs(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSc() const = 0;
   virtual double get_MSc(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSb() const = 0;
   virtual double get_MSb(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MSt() const = 0;
   virtual double get_MSt(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_Mhh() const = 0;
   virtual double get_Mhh(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MAh() const = 0;
   virtual double get_MAh(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MHpm() const = 0;
   virtual double get_MHpm(int i) const = 0;
   virtual const Eigen::Array<double,4,1>& get_MChi() const = 0;
   virtual double get_MChi(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MCha() const = 0;
   virtual double get_MCha(int i) const = 0;
   virtual double get_MVWm() const = 0;
   virtual double get_MVP() const = 0;
   virtual double get_MVZ() const = 0;
   
   virtual Eigen::Array<double,1,1> get_MChargedHiggs() const = 0;
   virtual Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZD() const = 0;
   virtual double get_ZD(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZU() const = 0;
   virtual double get_ZU(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZE() const = 0;
   virtual double get_ZE(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZM() const = 0;
   virtual double get_ZM(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZTau() const = 0;
   virtual double get_ZTau(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZS() const = 0;
   virtual double get_ZS(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZC() const = 0;
   virtual double get_ZC(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZB() const = 0;
   virtual double get_ZB(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZT() const = 0;
   virtual double get_ZT(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZH() const = 0;
   virtual double get_ZH(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZA() const = 0;
   virtual double get_ZA(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZP() const = 0;
   virtual double get_ZP(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,4,4>& get_ZN() const = 0;
   virtual std::complex<double> get_ZN(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,2,2>& get_UM() const = 0;
   virtual std::complex<double> get_UM(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,2,2>& get_UP() const = 0;
   virtual std::complex<double> get_UP(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZZ() const = 0;
   virtual double get_ZZ(int i, int k) const = 0;
   virtual void set_PhaseGlu(std::complex<double> PhaseGlu_) = 0;
   virtual std::complex<double> get_PhaseGlu() const = 0;


   virtual double get_mass_matrix_VG() const = 0;
   virtual void calculate_MVG() = 0;
   virtual double get_mass_matrix_Glu() const = 0;
   virtual void calculate_MGlu() = 0;
   virtual double get_mass_matrix_Fd() const = 0;
   virtual void calculate_MFd() = 0;
   virtual double get_mass_matrix_Fs() const = 0;
   virtual void calculate_MFs() = 0;
   virtual double get_mass_matrix_Fb() const = 0;
   virtual void calculate_MFb() = 0;
   virtual double get_mass_matrix_Fu() const = 0;
   virtual void calculate_MFu() = 0;
   virtual double get_mass_matrix_Fc() const = 0;
   virtual void calculate_MFc() = 0;
   virtual double get_mass_matrix_Ft() const = 0;
   virtual void calculate_MFt() = 0;
   virtual double get_mass_matrix_Fve() const = 0;
   virtual void calculate_MFve() = 0;
   virtual double get_mass_matrix_Fvm() const = 0;
   virtual void calculate_MFvm() = 0;
   virtual double get_mass_matrix_Fvt() const = 0;
   virtual void calculate_MFvt() = 0;
   virtual double get_mass_matrix_Fe() const = 0;
   virtual void calculate_MFe() = 0;
   virtual double get_mass_matrix_Fm() const = 0;
   virtual void calculate_MFm() = 0;
   virtual double get_mass_matrix_Ftau() const = 0;
   virtual void calculate_MFtau() = 0;
   virtual double get_mass_matrix_SveL() const = 0;
   virtual void calculate_MSveL() = 0;
   virtual double get_mass_matrix_SvmL() const = 0;
   virtual void calculate_MSvmL() = 0;
   virtual double get_mass_matrix_SvtL() const = 0;
   virtual void calculate_MSvtL() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Sd() const = 0;
   virtual void calculate_MSd() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Su() const = 0;
   virtual void calculate_MSu() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Se() const = 0;
   virtual void calculate_MSe() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Sm() const = 0;
   virtual void calculate_MSm() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Stau() const = 0;
   virtual void calculate_MStau() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Ss() const = 0;
   virtual void calculate_MSs() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Sc() const = 0;
   virtual void calculate_MSc() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Sb() const = 0;
   virtual void calculate_MSb() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_St() const = 0;
   virtual void calculate_MSt() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_hh() const = 0;
   virtual void calculate_Mhh() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const = 0;
   virtual void calculate_MAh() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Hpm() const = 0;
   virtual void calculate_MHpm() = 0;
   virtual Eigen::Matrix<double,4,4> get_mass_matrix_Chi() const = 0;
   virtual void calculate_MChi() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Cha() const = 0;
   virtual void calculate_MCha() = 0;
   virtual double get_mass_matrix_VWm() const = 0;
   virtual void calculate_MVWm() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const = 0;
   virtual void calculate_MVPVZ() = 0;
   virtual double get_ewsb_eq_hh_1() const = 0;
   virtual double get_ewsb_eq_hh_2() const = 0;
   virtual double v() const = 0;
   virtual double Betax() const = 0;
   virtual double Alpha() const = 0;
   virtual double ThetaW() const = 0;
   virtual double VEV() const = 0;
};

} // namespace flexiblesusy

#endif
