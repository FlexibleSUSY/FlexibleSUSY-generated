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
 * @file HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme.hpp
 *
 * @brief Defines model class for Stöckinger/Kotlarski decoupling scheme.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.4 .
 */

#ifndef HGTHDMIIMSSMBC_MASS_EIGENSTATES_DECOUPLING_SCHEME_H
#define HGTHDMIIMSSMBC_MASS_EIGENSTATES_DECOUPLING_SCHEME_H

#include "HGTHDMIIMSSMBC_info.hpp"
#include "HGTHDMIIMSSMBC_physical.hpp"
#include "HGTHDMIIMSSMBC_soft_parameters.hpp"
#include "HGTHDMIIMSSMBC_mass_eigenstates_interface.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>

#include <Eigen/Core>

#define SUPER(p) HGTHDMIIMSSMBC_soft_parameters::p

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

struct HGTHDMIIMSSMBC_input_parameters;
class HGTHDMIIMSSMBC_mass_eigenstates;

/**
 * @class HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme
 *
 * @brief model class with routines for determing masses and mixings
 * and EWSB in the Stöckinger/Kotlarski decoupling scheme
 */
class HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme
   : private HGTHDMIIMSSMBC_soft_parameters
   , public HGTHDMIIMSSMBC_mass_eigenstates_interface
{
public:
   explicit HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme(const HGTHDMIIMSSMBC_input_parameters& input_ = HGTHDMIIMSSMBC_input_parameters());
   explicit HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme(const HGTHDMIIMSSMBC_mass_eigenstates&);
   HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme(const HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme&) = default;
   HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme(HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme&&) = default;
   virtual ~HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme() = default;
   HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme& operator=(const HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme&) = default;
   HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme& operator=(HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme&&) = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   std::unique_ptr<HGTHDMIIMSSMBC_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 2;

   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void fill_from(const standard_model::Standard_model&);
   void fill_from(const HGTHDMIIMSSMBC_mass_eigenstates&);
   void reorder_tree_level_masses();
   void reorder_pole_masses();
   void print(std::ostream&) const override;
   void set_precision(double);
   double get_precision() const;
   void clear() override;

   // mass_eigenstates_interface functions

   void calculate_tree_level_mass_spectrum() override;
   void calculate_pole_mass_spectrum() override;
   void calculate_mass_spectrum() override;
   int solve_ewsb_equations_tree_level() override;
   int solve_ewsb_equations() override;
   Eigen::ArrayXd get_tree_level_masses() const override;
   Eigen::ArrayXd get_tree_level_masses_and_mixings() const override;
   const HGTHDMIIMSSMBC_input_parameters& get_input_parameters() const override;
   HGTHDMIIMSSMBC_input_parameters& get_input_parameters() override;
   Eigen::ArrayXd get_extra_parameters() const override;
   const HGTHDMIIMSSMBC_physical& get_physical() const override;
   HGTHDMIIMSSMBC_physical& get_physical() override;
   const Problems& get_problems() const override;
   Problems& get_problems() override;
   void set_tree_level_masses(const Eigen::ArrayXd&) override;
   void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) override;
   void set_extra_parameters(const Eigen::ArrayXd&) override;
   void set_physical(const HGTHDMIIMSSMBC_physical&) override;
   void clear_problems() override;

   double get_g1() const override { return SUPER(g1); }
   double get_g2() const override { return SUPER(g2); }
   double get_g3() const override { return SUPER(g3); }
   double get_Lambda6() const override { return SUPER(Lambda6); }
   double get_Lambda5() const override { return SUPER(Lambda5); }
   double get_Lambda7() const override { return SUPER(Lambda7); }
   double get_Lambda1() const override { return SUPER(Lambda1); }
   double get_Lambda4() const override { return SUPER(Lambda4); }
   double get_Lambda3() const override { return SUPER(Lambda3); }
   double get_Lambda2() const override { return SUPER(Lambda2); }
   const Eigen::Matrix<double,3,3>& get_Yu() const override { return SUPER(Yu); }
   double get_Yu(int i, int k) const override { return SUPER(Yu(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Yd() const override { return SUPER(Yd); }
   double get_Yd(int i, int k) const override { return SUPER(Yd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Ye() const override { return SUPER(Ye); }
   double get_Ye(int i, int k) const override { return SUPER(Ye(i,k)); }
   double get_g1dp() const override { return SUPER(g1dp); }
   double get_g1d() const override { return SUPER(g1d); }
   double get_g2up() const override { return SUPER(g2up); }
   double get_g2u() const override { return SUPER(g2u); }
   double get_MassB() const override { return SUPER(MassB); }
   double get_MassG() const override { return SUPER(MassG); }
   double get_MassWB() const override { return SUPER(MassWB); }
   double get_Mu() const override { return SUPER(Mu); }
   double get_M122() const override { return SUPER(M122); }
   double get_M112() const override { return SUPER(M112); }
   double get_M222() const override { return SUPER(M222); }
   double get_v1() const override { return SUPER(v1); }
   double get_v2() const override { return SUPER(v2); }

   void set_g1(double g1_) override { g1 = g1_; }
   void set_g2(double g2_) override { g2 = g2_; }
   void set_g3(double g3_) override { g3 = g3_; }
   void set_Lambda6(double Lambda6_) override { Lambda6 = Lambda6_; }
   void set_Lambda5(double Lambda5_) override { Lambda5 = Lambda5_; }
   void set_Lambda7(double Lambda7_) override { Lambda7 = Lambda7_; }
   void set_Lambda1(double Lambda1_) override { Lambda1 = Lambda1_; }
   void set_Lambda4(double Lambda4_) override { Lambda4 = Lambda4_; }
   void set_Lambda3(double Lambda3_) override { Lambda3 = Lambda3_; }
   void set_Lambda2(double Lambda2_) override { Lambda2 = Lambda2_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) override { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) override { Yu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) override { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) override { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) override { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) override { Ye(i,k) = value; }
   void set_g1dp(double g1dp_) override { g1dp = g1dp_; }
   void set_g1d(double g1d_) override { g1d = g1d_; }
   void set_g2up(double g2up_) override { g2up = g2up_; }
   void set_g2u(double g2u_) override { g2u = g2u_; }
   void set_MassB(double MassB_) override { MassB = MassB_; }
   void set_MassG(double MassG_) override { MassG = MassG_; }
   void set_MassWB(double MassWB_) override { MassWB = MassWB_; }
   void set_Mu(double Mu_) override { Mu = Mu_; }
   void set_M122(double M122_) override { M122 = M122_; }
   void set_M112(double M112_) override { M112 = M112_; }
   void set_M222(double M222_) override { M222 = M222_; }
   void set_v1(double v1_) override { v1 = v1_; }
   void set_v2(double v2_) override { v2 = v2_; }

   double get_MVG() const override { return MVG; }
   const Eigen::Array<double,3,1>& get_MFv() const override { return MFv; }
   double get_MFv(int i) const override { return MFv(i); }
   double get_MGlu() const override { return MGlu; }
   const Eigen::Array<double,2,1>& get_Mhh() const override { return Mhh; }
   double get_Mhh(int i) const override { return Mhh(i); }
   const Eigen::Array<double,2,1>& get_MAh() const override { return MAh; }
   double get_MAh(int i) const override { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHm() const override { return MHm; }
   double get_MHm(int i) const override { return MHm(i); }
   const Eigen::Array<double,3,1>& get_MFd() const override { return MFd; }
   double get_MFd(int i) const override { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const override { return MFu; }
   double get_MFu(int i) const override { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const override { return MFe; }
   double get_MFe(int i) const override { return MFe(i); }
   const Eigen::Array<double,4,1>& get_MChi() const override { return MChi; }
   double get_MChi(int i) const override { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha() const override { return MCha; }
   double get_MCha(int i) const override { return MCha(i); }
   double get_MVWm() const override { return MVWm; }
   double get_MVP() const override { return MVP; }
   double get_MVZ() const override { return MVZ; }

   
   Eigen::Array<double,1,1> get_MChargedHiggs() const override;

   Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const override;

   const Eigen::Matrix<double,2,2>& get_ZH() const override { return ZH; }
   double get_ZH(int i, int k) const override { return ZH(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZA() const override { return ZA; }
   double get_ZA(int i, int k) const override { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const override { return ZP; }
   double get_ZP(int i, int k) const override { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const override { return Vd; }
   std::complex<double> get_Vd(int i, int k) const override { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const override { return Ud; }
   std::complex<double> get_Ud(int i, int k) const override { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const override { return Vu; }
   std::complex<double> get_Vu(int i, int k) const override { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const override { return Uu; }
   std::complex<double> get_Uu(int i, int k) const override { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const override { return Ve; }
   std::complex<double> get_Ve(int i, int k) const override { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const override { return Ue; }
   std::complex<double> get_Ue(int i, int k) const override { return Ue(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN() const override { return ZN; }
   std::complex<double> get_ZN(int i, int k) const override { return ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM() const override { return UM; }
   std::complex<double> get_UM(int i, int k) const override { return UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP() const override { return UP; }
   std::complex<double> get_UP(int i, int k) const override { return UP(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const override { return ZZ; }
   double get_ZZ(int i, int k) const override { return ZZ(i,k); }

   void set_PhaseGlu(std::complex<double> PhaseGlu_) override { PhaseGlu = PhaseGlu_; }
   std::complex<double> get_PhaseGlu() const override { return PhaseGlu; }



   double get_mass_matrix_VG() const override;
   void calculate_MVG() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const override;
   void calculate_MFv() override;
   double get_mass_matrix_Glu() const override;
   void calculate_MGlu() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_hh() const override;
   void calculate_Mhh() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const override;
   void calculate_MAh() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Hm() const override;
   void calculate_MHm() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const override;
   void calculate_MFd() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const override;
   void calculate_MFu() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const override;
   void calculate_MFe() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_Chi() const override;
   void calculate_MChi() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha() const override;
   void calculate_MCha() override;
   double get_mass_matrix_VWm() const override;
   void calculate_MVWm() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const override;
   void calculate_MVPVZ() override;

   double get_ewsb_eq_hh_1() const override;
   double get_ewsb_eq_hh_2() const override;

   double v() const override;
   double Betax() const override;
   double Alpha() const override;
   double ThetaW() const override;
   double VEV() const override;


private:
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< mass eigenstate precision
   HGTHDMIIMSSMBC_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{HGTHDMIIMSSMBC_info::model_name,
                     &HGTHDMIIMSSMBC_info::particle_names_getter,
                     &HGTHDMIIMSSMBC_info::parameter_names_getter}; ///< problems

   void clear_tree_level_parameters();
   void copy_tree_level_masses_to_pole_masses();

   // DR-bar masses
   double MVG{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double MGlu{};
   Eigen::Array<double,2,1> Mhh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MAh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MHm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,4,1> MChi{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MCha{Eigen::Array<double,2,1>::Zero()};
   double MVWm{};
   double MVP{};
   double MVZ{};

   // DR-bar mixing matrices
   Eigen::Matrix<double,2,2> ZH{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZA{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,4,4> ZN{Eigen::Matrix<std::complex<double>,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};

   // phases
   std::complex<double> PhaseGlu{1.,0.};

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const HGTHDMIIMSSMBC_mass_eigenstates_decoupling_scheme&);

} // namespace flexiblesusy

#undef SUPER

#endif
