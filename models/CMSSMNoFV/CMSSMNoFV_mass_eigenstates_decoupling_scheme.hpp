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
 * @file CMSSMNoFV_mass_eigenstates_decoupling_scheme.hpp
 *
 * @brief Defines model class for Stöckinger/Kotlarski decoupling scheme.
 *
 * This file was generated with FlexibleSUSY 2.5.0 and SARAH 4.14.4 .
 */

#ifndef CMSSMNoFV_MASS_EIGENSTATES_DECOUPLING_SCHEME_H
#define CMSSMNoFV_MASS_EIGENSTATES_DECOUPLING_SCHEME_H

#include "CMSSMNoFV_info.hpp"
#include "CMSSMNoFV_physical.hpp"
#include "CMSSMNoFV_soft_parameters.hpp"
#include "CMSSMNoFV_mass_eigenstates_interface.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>

#include <Eigen/Core>

#define SUPER(p) CMSSMNoFV_soft_parameters::p

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

struct CMSSMNoFV_input_parameters;
class CMSSMNoFV_mass_eigenstates;

/**
 * @class CMSSMNoFV_mass_eigenstates_decoupling_scheme
 *
 * @brief model class with routines for determing masses and mixings
 * and EWSB in the Stöckinger/Kotlarski decoupling scheme
 */
class CMSSMNoFV_mass_eigenstates_decoupling_scheme
   : private CMSSMNoFV_soft_parameters
   , public CMSSMNoFV_mass_eigenstates_interface
{
public:
   explicit CMSSMNoFV_mass_eigenstates_decoupling_scheme(const CMSSMNoFV_input_parameters& input_ = CMSSMNoFV_input_parameters());
   explicit CMSSMNoFV_mass_eigenstates_decoupling_scheme(const CMSSMNoFV_mass_eigenstates&);
   CMSSMNoFV_mass_eigenstates_decoupling_scheme(const CMSSMNoFV_mass_eigenstates_decoupling_scheme&) = default;
   CMSSMNoFV_mass_eigenstates_decoupling_scheme(CMSSMNoFV_mass_eigenstates_decoupling_scheme&&) = default;
   virtual ~CMSSMNoFV_mass_eigenstates_decoupling_scheme() = default;
   CMSSMNoFV_mass_eigenstates_decoupling_scheme& operator=(const CMSSMNoFV_mass_eigenstates_decoupling_scheme&) = default;
   CMSSMNoFV_mass_eigenstates_decoupling_scheme& operator=(CMSSMNoFV_mass_eigenstates_decoupling_scheme&&) = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   std::unique_ptr<CMSSMNoFV_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 2;

   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void fill_from(const standard_model::Standard_model&);
   void fill_from(const CMSSMNoFV_mass_eigenstates&);
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
   const CMSSMNoFV_input_parameters& get_input_parameters() const override;
   CMSSMNoFV_input_parameters& get_input_parameters() override;
   Eigen::ArrayXd get_extra_parameters() const override;
   const CMSSMNoFV_physical& get_physical() const override;
   CMSSMNoFV_physical& get_physical() override;
   const Problems& get_problems() const override;
   Problems& get_problems() override;
   void set_tree_level_masses(const Eigen::ArrayXd&) override;
   void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) override;
   void set_extra_parameters(const Eigen::ArrayXd&) override;
   void set_physical(const CMSSMNoFV_physical&) override;
   void clear_problems() override;

   const Eigen::Matrix<double,3,3>& get_Yd() const override { return SUPER(Yd); }
   double get_Yd(int i, int k) const override { return SUPER(Yd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Ye() const override { return SUPER(Ye); }
   double get_Ye(int i, int k) const override { return SUPER(Ye(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Yu() const override { return SUPER(Yu); }
   double get_Yu(int i, int k) const override { return SUPER(Yu(i,k)); }
   double get_Mu() const override { return SUPER(Mu); }
   double get_g1() const override { return SUPER(g1); }
   double get_g2() const override { return SUPER(g2); }
   double get_g3() const override { return SUPER(g3); }
   double get_vd() const override { return SUPER(vd); }
   double get_vu() const override { return SUPER(vu); }
   const Eigen::Matrix<double,3,3>& get_TYd() const override { return SUPER(TYd); }
   double get_TYd(int i, int k) const override { return SUPER(TYd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_TYe() const override { return SUPER(TYe); }
   double get_TYe(int i, int k) const override { return SUPER(TYe(i,k)); }
   const Eigen::Matrix<double,3,3>& get_TYu() const override { return SUPER(TYu); }
   double get_TYu(int i, int k) const override { return SUPER(TYu(i,k)); }
   double get_BMu() const override { return SUPER(BMu); }
   const Eigen::Matrix<double,3,3>& get_mq2() const override { return SUPER(mq2); }
   double get_mq2(int i, int k) const override { return SUPER(mq2(i,k)); }
   const Eigen::Matrix<double,3,3>& get_ml2() const override { return SUPER(ml2); }
   double get_ml2(int i, int k) const override { return SUPER(ml2(i,k)); }
   double get_mHd2() const override { return SUPER(mHd2); }
   double get_mHu2() const override { return SUPER(mHu2); }
   const Eigen::Matrix<double,3,3>& get_md2() const override { return SUPER(md2); }
   double get_md2(int i, int k) const override { return SUPER(md2(i,k)); }
   const Eigen::Matrix<double,3,3>& get_mu2() const override { return SUPER(mu2); }
   double get_mu2(int i, int k) const override { return SUPER(mu2(i,k)); }
   const Eigen::Matrix<double,3,3>& get_me2() const override { return SUPER(me2); }
   double get_me2(int i, int k) const override { return SUPER(me2(i,k)); }
   double get_MassB() const override { return SUPER(MassB); }
   double get_MassWB() const override { return SUPER(MassWB); }
   double get_MassG() const override { return SUPER(MassG); }

   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) override { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) override { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) override { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) override { Ye(i,k) = value; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) override { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) override { Yu(i,k) = value; }
   void set_Mu(double Mu_) override { Mu = Mu_; }
   void set_g1(double g1_) override { g1 = g1_; }
   void set_g2(double g2_) override { g2 = g2_; }
   void set_g3(double g3_) override { g3 = g3_; }
   void set_vd(double vd_) override { vd = vd_; }
   void set_vu(double vu_) override { vu = vu_; }
   void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) override { TYd = TYd_; }
   void set_TYd(int i, int k, const double& value) override { TYd(i,k) = value; }
   void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) override { TYe = TYe_; }
   void set_TYe(int i, int k, const double& value) override { TYe(i,k) = value; }
   void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) override { TYu = TYu_; }
   void set_TYu(int i, int k, const double& value) override { TYu(i,k) = value; }
   void set_BMu(double BMu_) override { BMu = BMu_; }
   void set_mq2(const Eigen::Matrix<double,3,3>& mq2_) override { mq2 = mq2_; }
   void set_mq2(int i, int k, const double& value) override { mq2(i,k) = value; }
   void set_ml2(const Eigen::Matrix<double,3,3>& ml2_) override { ml2 = ml2_; }
   void set_ml2(int i, int k, const double& value) override { ml2(i,k) = value; }
   void set_mHd2(double mHd2_) override { mHd2 = mHd2_; }
   void set_mHu2(double mHu2_) override { mHu2 = mHu2_; }
   void set_md2(const Eigen::Matrix<double,3,3>& md2_) override { md2 = md2_; }
   void set_md2(int i, int k, const double& value) override { md2(i,k) = value; }
   void set_mu2(const Eigen::Matrix<double,3,3>& mu2_) override { mu2 = mu2_; }
   void set_mu2(int i, int k, const double& value) override { mu2(i,k) = value; }
   void set_me2(const Eigen::Matrix<double,3,3>& me2_) override { me2 = me2_; }
   void set_me2(int i, int k, const double& value) override { me2(i,k) = value; }
   void set_MassB(double MassB_) override { MassB = MassB_; }
   void set_MassWB(double MassWB_) override { MassWB = MassWB_; }
   void set_MassG(double MassG_) override { MassG = MassG_; }

   double get_MVG() const override { return MVG; }
   double get_MGlu() const override { return MGlu; }
   double get_MFd() const override { return MFd; }
   double get_MFs() const override { return MFs; }
   double get_MFb() const override { return MFb; }
   double get_MFu() const override { return MFu; }
   double get_MFc() const override { return MFc; }
   double get_MFt() const override { return MFt; }
   double get_MFve() const override { return MFve; }
   double get_MFvm() const override { return MFvm; }
   double get_MFvt() const override { return MFvt; }
   double get_MFe() const override { return MFe; }
   double get_MFm() const override { return MFm; }
   double get_MFtau() const override { return MFtau; }
   double get_MSveL() const override { return MSveL; }
   double get_MSvmL() const override { return MSvmL; }
   double get_MSvtL() const override { return MSvtL; }
   const Eigen::Array<double,2,1>& get_MSd() const override { return MSd; }
   double get_MSd(int i) const override { return MSd(i); }
   const Eigen::Array<double,2,1>& get_MSu() const override { return MSu; }
   double get_MSu(int i) const override { return MSu(i); }
   const Eigen::Array<double,2,1>& get_MSe() const override { return MSe; }
   double get_MSe(int i) const override { return MSe(i); }
   const Eigen::Array<double,2,1>& get_MSm() const override { return MSm; }
   double get_MSm(int i) const override { return MSm(i); }
   const Eigen::Array<double,2,1>& get_MStau() const override { return MStau; }
   double get_MStau(int i) const override { return MStau(i); }
   const Eigen::Array<double,2,1>& get_MSs() const override { return MSs; }
   double get_MSs(int i) const override { return MSs(i); }
   const Eigen::Array<double,2,1>& get_MSc() const override { return MSc; }
   double get_MSc(int i) const override { return MSc(i); }
   const Eigen::Array<double,2,1>& get_MSb() const override { return MSb; }
   double get_MSb(int i) const override { return MSb(i); }
   const Eigen::Array<double,2,1>& get_MSt() const override { return MSt; }
   double get_MSt(int i) const override { return MSt(i); }
   const Eigen::Array<double,2,1>& get_Mhh() const override { return Mhh; }
   double get_Mhh(int i) const override { return Mhh(i); }
   const Eigen::Array<double,2,1>& get_MAh() const override { return MAh; }
   double get_MAh(int i) const override { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm() const override { return MHpm; }
   double get_MHpm(int i) const override { return MHpm(i); }
   const Eigen::Array<double,4,1>& get_MChi() const override { return MChi; }
   double get_MChi(int i) const override { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha() const override { return MCha; }
   double get_MCha(int i) const override { return MCha(i); }
   double get_MVWm() const override { return MVWm; }
   double get_MVP() const override { return MVP; }
   double get_MVZ() const override { return MVZ; }

   
   Eigen::Array<double,1,1> get_MChargedHiggs() const override;

   Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const override;

   const Eigen::Matrix<double,2,2>& get_ZD() const override { return ZD; }
   double get_ZD(int i, int k) const override { return ZD(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZU() const override { return ZU; }
   double get_ZU(int i, int k) const override { return ZU(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZE() const override { return ZE; }
   double get_ZE(int i, int k) const override { return ZE(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZM() const override { return ZM; }
   double get_ZM(int i, int k) const override { return ZM(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZTau() const override { return ZTau; }
   double get_ZTau(int i, int k) const override { return ZTau(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZS() const override { return ZS; }
   double get_ZS(int i, int k) const override { return ZS(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZC() const override { return ZC; }
   double get_ZC(int i, int k) const override { return ZC(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZB() const override { return ZB; }
   double get_ZB(int i, int k) const override { return ZB(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZT() const override { return ZT; }
   double get_ZT(int i, int k) const override { return ZT(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZH() const override { return ZH; }
   double get_ZH(int i, int k) const override { return ZH(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZA() const override { return ZA; }
   double get_ZA(int i, int k) const override { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const override { return ZP; }
   double get_ZP(int i, int k) const override { return ZP(i,k); }
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
   double get_mass_matrix_Glu() const override;
   void calculate_MGlu() override;
   double get_mass_matrix_Fd() const override;
   void calculate_MFd() override;
   double get_mass_matrix_Fs() const override;
   void calculate_MFs() override;
   double get_mass_matrix_Fb() const override;
   void calculate_MFb() override;
   double get_mass_matrix_Fu() const override;
   void calculate_MFu() override;
   double get_mass_matrix_Fc() const override;
   void calculate_MFc() override;
   double get_mass_matrix_Ft() const override;
   void calculate_MFt() override;
   double get_mass_matrix_Fve() const override;
   void calculate_MFve() override;
   double get_mass_matrix_Fvm() const override;
   void calculate_MFvm() override;
   double get_mass_matrix_Fvt() const override;
   void calculate_MFvt() override;
   double get_mass_matrix_Fe() const override;
   void calculate_MFe() override;
   double get_mass_matrix_Fm() const override;
   void calculate_MFm() override;
   double get_mass_matrix_Ftau() const override;
   void calculate_MFtau() override;
   double get_mass_matrix_SveL() const override;
   void calculate_MSveL() override;
   double get_mass_matrix_SvmL() const override;
   void calculate_MSvmL() override;
   double get_mass_matrix_SvtL() const override;
   void calculate_MSvtL() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Sd() const override;
   void calculate_MSd() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Su() const override;
   void calculate_MSu() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Se() const override;
   void calculate_MSe() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Sm() const override;
   void calculate_MSm() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Stau() const override;
   void calculate_MStau() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Ss() const override;
   void calculate_MSs() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Sc() const override;
   void calculate_MSc() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Sb() const override;
   void calculate_MSb() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_St() const override;
   void calculate_MSt() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_hh() const override;
   void calculate_Mhh() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const override;
   void calculate_MAh() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Hpm() const override;
   void calculate_MHpm() override;
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
   CMSSMNoFV_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{CMSSMNoFV_info::model_name,
                     &CMSSMNoFV_info::particle_names_getter,
                     &CMSSMNoFV_info::parameter_names_getter}; ///< problems

   void clear_tree_level_parameters();
   void copy_tree_level_masses_to_pole_masses();

   // DR-bar masses
   double MVG{};
   double MGlu{};
   double MFd{};
   double MFs{};
   double MFb{};
   double MFu{};
   double MFc{};
   double MFt{};
   double MFve{};
   double MFvm{};
   double MFvt{};
   double MFe{};
   double MFm{};
   double MFtau{};
   double MSveL{};
   double MSvmL{};
   double MSvtL{};
   Eigen::Array<double,2,1> MSd{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSu{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSe{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MStau{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSs{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSc{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSb{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSt{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> Mhh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MAh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MHpm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,4,1> MChi{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MCha{Eigen::Array<double,2,1>::Zero()};
   double MVWm{};
   double MVP{};
   double MVZ{};

   // DR-bar mixing matrices
   Eigen::Matrix<double,2,2> ZD{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZU{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZE{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZM{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZTau{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZS{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZC{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZB{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZT{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZH{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZA{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,4,4> ZN{Eigen::Matrix<std::complex<double>,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};

   // phases
   std::complex<double> PhaseGlu{1.,0.};

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const CMSSMNoFV_mass_eigenstates_decoupling_scheme&);

} // namespace flexiblesusy

#undef SUPER

#endif
