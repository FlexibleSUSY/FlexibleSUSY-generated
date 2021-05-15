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
 * @file SplitMSSM_mass_eigenstates_interface.hpp
 *
 * @brief Contains the mass eigenstates interface class definition
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.4 .
 */

#ifndef SplitMSSM_MASS_EIGENSTATES_INTERFACE_H
#define SplitMSSM_MASS_EIGENSTATES_INTERFACE_H

#include "SplitMSSM_physical.hpp"

#include <memory>

#include <Eigen/Core>

namespace flexiblesusy {

class Problems;
struct SplitMSSM_input_parameters;

/**
 * @class SplitMSSM_mass_eigenstates_interface
 * @brief Interface definition for model parameters, masses and mixings
 */
class SplitMSSM_mass_eigenstates_interface {
public:
   virtual ~SplitMSSM_mass_eigenstates_interface() {}

   virtual std::unique_ptr<SplitMSSM_mass_eigenstates_interface> clone() const = 0;

   virtual void calculate_tree_level_mass_spectrum() = 0;
   virtual void calculate_pole_mass_spectrum() = 0;
   virtual void calculate_mass_spectrum() = 0;

   virtual int solve_ewsb_equations_tree_level() = 0;
   virtual int solve_ewsb_equations() = 0;

   virtual Eigen::ArrayXd get_tree_level_masses() const = 0;
   virtual Eigen::ArrayXd get_tree_level_masses_and_mixings() const = 0;
   virtual const SplitMSSM_input_parameters& get_input_parameters() const = 0;
   virtual SplitMSSM_input_parameters& get_input_parameters() = 0;
   virtual Eigen::ArrayXd get_extra_parameters() const = 0;
   virtual const SplitMSSM_physical& get_physical() const = 0;
   virtual SplitMSSM_physical& get_physical() = 0;
   virtual const Problems& get_problems() const = 0;
   virtual Problems& get_problems() = 0;
   virtual void set_tree_level_masses(const Eigen::ArrayXd&) = 0;
   virtual void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) = 0;
   virtual void set_extra_parameters(const Eigen::ArrayXd&) = 0;
   virtual void set_physical(const SplitMSSM_physical&) = 0;
   virtual void clear_problems() = 0;

   virtual double get_g1() const = 0;
   virtual double get_g2() const = 0;
   virtual double get_g3() const = 0;
   virtual double get_Lambdax() const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yu() const = 0;
   virtual double get_Yu(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yd() const = 0;
   virtual double get_Yd(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Ye() const = 0;
   virtual double get_Ye(int i, int k) const = 0;
   virtual double get_gYd() const = 0;
   virtual double get_g2d() const = 0;
   virtual double get_gYu() const = 0;
   virtual double get_g2u() const = 0;
   virtual double get_MassB() const = 0;
   virtual double get_MassG() const = 0;
   virtual double get_MassWB() const = 0;
   virtual double get_Mu() const = 0;
   virtual double get_mu2() const = 0;
   virtual double get_v() const = 0;
   virtual void set_g1(double g1_) = 0;
   virtual void set_g2(double g2_) = 0;
   virtual void set_g3(double g3_) = 0;
   virtual void set_Lambdax(double Lambdax_) = 0;
   virtual void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) = 0;
   virtual void set_Yu(int i, int k, const double& value) = 0;
   virtual void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) = 0;
   virtual void set_Yd(int i, int k, const double& value) = 0;
   virtual void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) = 0;
   virtual void set_Ye(int i, int k, const double& value) = 0;
   virtual void set_gYd(double gYd_) = 0;
   virtual void set_g2d(double g2d_) = 0;
   virtual void set_gYu(double gYu_) = 0;
   virtual void set_g2u(double g2u_) = 0;
   virtual void set_MassB(double MassB_) = 0;
   virtual void set_MassG(double MassG_) = 0;
   virtual void set_MassWB(double MassWB_) = 0;
   virtual void set_Mu(double Mu_) = 0;
   virtual void set_mu2(double mu2_) = 0;
   virtual void set_v(double v_) = 0;
   virtual double get_MVG() const = 0;
   virtual double get_MHp() const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFv() const = 0;
   virtual double get_MFv(int i) const = 0;
   virtual double get_MGlu() const = 0;
   virtual double get_MAh() const = 0;
   virtual double get_Mhh() const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFd() const = 0;
   virtual double get_MFd(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFu() const = 0;
   virtual double get_MFu(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFe() const = 0;
   virtual double get_MFe(int i) const = 0;
   virtual const Eigen::Array<double,4,1>& get_MChi() const = 0;
   virtual double get_MChi(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MCha() const = 0;
   virtual double get_MCha(int i) const = 0;
   virtual double get_MVWp() const = 0;
   virtual double get_MVP() const = 0;
   virtual double get_MVZ() const = 0;

   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const = 0;
   virtual std::complex<double> get_Vd(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const = 0;
   virtual std::complex<double> get_Ud(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const = 0;
   virtual std::complex<double> get_Vu(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const = 0;
   virtual std::complex<double> get_Uu(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const = 0;
   virtual std::complex<double> get_Ve(int i, int k) const = 0;
   virtual const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const = 0;
   virtual std::complex<double> get_Ue(int i, int k) const = 0;
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
   virtual double get_mass_matrix_Hp() const = 0;
   virtual void calculate_MHp() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const = 0;
   virtual void calculate_MFv() = 0;
   virtual double get_mass_matrix_Glu() const = 0;
   virtual void calculate_MGlu() = 0;
   virtual double get_mass_matrix_Ah() const = 0;
   virtual void calculate_MAh() = 0;
   virtual double get_mass_matrix_hh() const = 0;
   virtual void calculate_Mhh() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const = 0;
   virtual void calculate_MFd() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const = 0;
   virtual void calculate_MFu() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const = 0;
   virtual void calculate_MFe() = 0;
   virtual Eigen::Matrix<double,4,4> get_mass_matrix_Chi() const = 0;
   virtual void calculate_MChi() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Cha() const = 0;
   virtual void calculate_MCha() = 0;
   virtual double get_mass_matrix_VWp() const = 0;
   virtual void calculate_MVWp() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const = 0;
   virtual void calculate_MVPVZ() = 0;
   virtual double get_ewsb_eq_hh_1() const = 0;
   virtual double ThetaW() const = 0;
   virtual double VEV() const = 0;
};

} // namespace flexiblesusy

#endif
