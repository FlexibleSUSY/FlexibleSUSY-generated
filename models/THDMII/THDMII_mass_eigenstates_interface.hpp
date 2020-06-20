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
 * @file THDMII_mass_eigenstates_interface.hpp
 *
 * @brief Contains the mass eigenstates interface class definition
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.3 .
 */

#ifndef THDMII_MASS_EIGENSTATES_INTERFACE_H
#define THDMII_MASS_EIGENSTATES_INTERFACE_H

#include "THDMII_physical.hpp"

#include <memory>

#include <Eigen/Core>

namespace flexiblesusy {

class Problems;
struct THDMII_input_parameters;

/**
 * @class THDMII_mass_eigenstates_interface
 * @brief Interface definition for model parameters, masses and mixings
 */
class THDMII_mass_eigenstates_interface {
public:
   virtual ~THDMII_mass_eigenstates_interface() {}

   virtual std::unique_ptr<THDMII_mass_eigenstates_interface> clone() const = 0;

   virtual void calculate_tree_level_mass_spectrum() = 0;
   virtual void calculate_pole_mass_spectrum() = 0;
   virtual void calculate_mass_spectrum() = 0;

   virtual int solve_ewsb_equations_tree_level() = 0;
   virtual int solve_ewsb_equations() = 0;

   virtual Eigen::ArrayXd get_tree_level_masses() const = 0;
   virtual Eigen::ArrayXd get_tree_level_masses_and_mixings() const = 0;
   virtual const THDMII_input_parameters& get_input_parameters() const = 0;
   virtual THDMII_input_parameters& get_input_parameters() = 0;
   virtual Eigen::ArrayXd get_extra_parameters() const = 0;
   virtual const THDMII_physical& get_physical() const = 0;
   virtual THDMII_physical& get_physical() = 0;
   virtual const Problems& get_problems() const = 0;
   virtual Problems& get_problems() = 0;
   virtual void set_tree_level_masses(const Eigen::ArrayXd&) = 0;
   virtual void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) = 0;
   virtual void set_extra_parameters(const Eigen::ArrayXd&) = 0;
   virtual void set_physical(const THDMII_physical&) = 0;
   virtual void clear_problems() = 0;

   virtual double get_g1() const = 0;
   virtual double get_g2() const = 0;
   virtual double get_g3() const = 0;
   virtual double get_Lambda6() const = 0;
   virtual double get_Lambda5() const = 0;
   virtual double get_Lambda7() const = 0;
   virtual double get_Lambda1() const = 0;
   virtual double get_Lambda4() const = 0;
   virtual double get_Lambda3() const = 0;
   virtual double get_Lambda2() const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yu() const = 0;
   virtual double get_Yu(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Yd() const = 0;
   virtual double get_Yd(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,3,3>& get_Ye() const = 0;
   virtual double get_Ye(int i, int k) const = 0;
   virtual double get_M122() const = 0;
   virtual double get_M112() const = 0;
   virtual double get_M222() const = 0;
   virtual double get_v1() const = 0;
   virtual double get_v2() const = 0;
   virtual void set_g1(double g1_) = 0;
   virtual void set_g2(double g2_) = 0;
   virtual void set_g3(double g3_) = 0;
   virtual void set_Lambda6(double Lambda6_) = 0;
   virtual void set_Lambda5(double Lambda5_) = 0;
   virtual void set_Lambda7(double Lambda7_) = 0;
   virtual void set_Lambda1(double Lambda1_) = 0;
   virtual void set_Lambda4(double Lambda4_) = 0;
   virtual void set_Lambda3(double Lambda3_) = 0;
   virtual void set_Lambda2(double Lambda2_) = 0;
   virtual void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) = 0;
   virtual void set_Yu(int i, int k, const double& value) = 0;
   virtual void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) = 0;
   virtual void set_Yd(int i, int k, const double& value) = 0;
   virtual void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) = 0;
   virtual void set_Ye(int i, int k, const double& value) = 0;
   virtual void set_M122(double M122_) = 0;
   virtual void set_M112(double M112_) = 0;
   virtual void set_M222(double M222_) = 0;
   virtual void set_v1(double v1_) = 0;
   virtual void set_v2(double v2_) = 0;
   virtual double get_MVG() const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFv() const = 0;
   virtual double get_MFv(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_Mhh() const = 0;
   virtual double get_Mhh(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MAh() const = 0;
   virtual double get_MAh(int i) const = 0;
   virtual const Eigen::Array<double,2,1>& get_MHm() const = 0;
   virtual double get_MHm(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFd() const = 0;
   virtual double get_MFd(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFu() const = 0;
   virtual double get_MFu(int i) const = 0;
   virtual const Eigen::Array<double,3,1>& get_MFe() const = 0;
   virtual double get_MFe(int i) const = 0;
   virtual double get_MVWm() const = 0;
   virtual double get_MVP() const = 0;
   virtual double get_MVZ() const = 0;
   
   virtual Eigen::Array<double,1,1> get_MChargedHiggs() const = 0;
   virtual Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZH() const = 0;
   virtual double get_ZH(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZA() const = 0;
   virtual double get_ZA(int i, int k) const = 0;
   virtual const Eigen::Matrix<double,2,2>& get_ZP() const = 0;
   virtual double get_ZP(int i, int k) const = 0;
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
   virtual const Eigen::Matrix<double,2,2>& get_ZZ() const = 0;
   virtual double get_ZZ(int i, int k) const = 0;



   virtual double get_mass_matrix_VG() const = 0;
   virtual void calculate_MVG() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const = 0;
   virtual void calculate_MFv() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_hh() const = 0;
   virtual void calculate_Mhh() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const = 0;
   virtual void calculate_MAh() = 0;
   virtual Eigen::Matrix<double,2,2> get_mass_matrix_Hm() const = 0;
   virtual void calculate_MHm() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const = 0;
   virtual void calculate_MFd() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const = 0;
   virtual void calculate_MFu() = 0;
   virtual Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const = 0;
   virtual void calculate_MFe() = 0;
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
