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
 * @file MRSSM2_mass_eigenstates_decoupling_scheme.hpp
 *
 * @brief Defines model class for Stöckinger/Kotlarski decoupling scheme.
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#ifndef MRSSM2_MASS_EIGENSTATES_DECOUPLING_SCHEME_H
#define MRSSM2_MASS_EIGENSTATES_DECOUPLING_SCHEME_H

#include "MRSSM2_info.hpp"
#include "MRSSM2_physical.hpp"
#include "MRSSM2_soft_parameters.hpp"
#include "MRSSM2_mass_eigenstates_interface.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>

#include <Eigen/Core>

#define SUPER(p) MRSSM2_soft_parameters::p

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

struct MRSSM2_input_parameters;
class MRSSM2_mass_eigenstates;

/**
 * @class MRSSM2_mass_eigenstates_decoupling_scheme
 *
 * @brief model class with routines for determing masses and mixings
 * and EWSB in the Stöckinger/Kotlarski decoupling scheme
 */
class MRSSM2_mass_eigenstates_decoupling_scheme
   : private MRSSM2_soft_parameters
   , public MRSSM2_mass_eigenstates_interface
{
public:
   explicit MRSSM2_mass_eigenstates_decoupling_scheme(const MRSSM2_input_parameters& input_ = MRSSM2_input_parameters());
   explicit MRSSM2_mass_eigenstates_decoupling_scheme(const MRSSM2_mass_eigenstates&);
   MRSSM2_mass_eigenstates_decoupling_scheme(const MRSSM2_mass_eigenstates_decoupling_scheme&) = default;
   MRSSM2_mass_eigenstates_decoupling_scheme(MRSSM2_mass_eigenstates_decoupling_scheme&&) = default;
   virtual ~MRSSM2_mass_eigenstates_decoupling_scheme() = default;
   MRSSM2_mass_eigenstates_decoupling_scheme& operator=(const MRSSM2_mass_eigenstates_decoupling_scheme&) = default;
   MRSSM2_mass_eigenstates_decoupling_scheme& operator=(MRSSM2_mass_eigenstates_decoupling_scheme&&) = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   std::unique_ptr<MRSSM2_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 4;

   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void fill_from(const standard_model::Standard_model&);
   void fill_from(const MRSSM2_mass_eigenstates&);
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
   const MRSSM2_input_parameters& get_input_parameters() const override;
   MRSSM2_input_parameters& get_input_parameters() override;
   Eigen::ArrayXd get_extra_parameters() const override;
   const MRSSM2_physical& get_physical() const override;
   MRSSM2_physical& get_physical() override;
   const Problems& get_problems() const override;
   Problems& get_problems() override;
   void set_tree_level_masses(const Eigen::ArrayXd&) override;
   void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) override;
   void set_extra_parameters(const Eigen::ArrayXd&) override;
   void set_physical(const MRSSM2_physical&) override;
   void clear_problems() override;

   const Eigen::Matrix<double,3,3>& get_Yd() const override { return SUPER(Yd); }
   double get_Yd(int i, int k) const override { return SUPER(Yd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Ye() const override { return SUPER(Ye); }
   double get_Ye(int i, int k) const override { return SUPER(Ye(i,k)); }
   double get_LamTD() const override { return SUPER(LamTD); }
   double get_LamTU() const override { return SUPER(LamTU); }
   double get_LamSD() const override { return SUPER(LamSD); }
   double get_LamSU() const override { return SUPER(LamSU); }
   const Eigen::Matrix<double,3,3>& get_Yu() const override { return SUPER(Yu); }
   double get_Yu(int i, int k) const override { return SUPER(Yu(i,k)); }
   double get_Mu() const override { return SUPER(Mu); }
   double get_MuD() const override { return SUPER(MuD); }
   double get_MuU() const override { return SUPER(MuU); }
   double get_g1() const override { return SUPER(g1); }
   double get_g2() const override { return SUPER(g2); }
   double get_g3() const override { return SUPER(g3); }
   double get_vd() const override { return SUPER(vd); }
   double get_vu() const override { return SUPER(vu); }
   double get_vT() const override { return SUPER(vT); }
   double get_vS() const override { return SUPER(vS); }
   double get_BMu() const override { return SUPER(BMu); }
   double get_BMuD() const override { return SUPER(BMuD); }
   double get_BMuU() const override { return SUPER(BMuU); }
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
   double get_mS2() const override { return SUPER(mS2); }
   double get_mT2() const override { return SUPER(mT2); }
   double get_moc2() const override { return SUPER(moc2); }
   double get_mRd2() const override { return SUPER(mRd2); }
   double get_mRu2() const override { return SUPER(mRu2); }
   double get_MDBS() const override { return SUPER(MDBS); }
   double get_MDWBT() const override { return SUPER(MDWBT); }
   double get_MDGoc() const override { return SUPER(MDGoc); }

   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) override { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) override { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) override { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) override { Ye(i,k) = value; }
   void set_LamTD(double LamTD_) override { LamTD = LamTD_; }
   void set_LamTU(double LamTU_) override { LamTU = LamTU_; }
   void set_LamSD(double LamSD_) override { LamSD = LamSD_; }
   void set_LamSU(double LamSU_) override { LamSU = LamSU_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) override { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) override { Yu(i,k) = value; }
   void set_Mu(double Mu_) override { Mu = Mu_; }
   void set_MuD(double MuD_) override { MuD = MuD_; }
   void set_MuU(double MuU_) override { MuU = MuU_; }
   void set_g1(double g1_) override { g1 = g1_; }
   void set_g2(double g2_) override { g2 = g2_; }
   void set_g3(double g3_) override { g3 = g3_; }
   void set_vd(double vd_) override { vd = vd_; }
   void set_vu(double vu_) override { vu = vu_; }
   void set_vT(double vT_) override { vT = vT_; }
   void set_vS(double vS_) override { vS = vS_; }
   void set_BMu(double BMu_) override { BMu = BMu_; }
   void set_BMuD(double BMuD_) override { BMuD = BMuD_; }
   void set_BMuU(double BMuU_) override { BMuU = BMuU_; }
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
   void set_mS2(double mS2_) override { mS2 = mS2_; }
   void set_mT2(double mT2_) override { mT2 = mT2_; }
   void set_moc2(double moc2_) override { moc2 = moc2_; }
   void set_mRd2(double mRd2_) override { mRd2 = mRd2_; }
   void set_mRu2(double mRu2_) override { mRu2 = mRu2_; }
   void set_MDBS(double MDBS_) override { MDBS = MDBS_; }
   void set_MDWBT(double MDWBT_) override { MDWBT = MDWBT_; }
   void set_MDGoc(double MDGoc_) override { MDGoc = MDGoc_; }

   double get_MVG() const override { return MVG; }
   double get_MGlu() const override { return MGlu; }
   const Eigen::Array<double,3,1>& get_MFv() const override { return MFv; }
   double get_MFv(int i) const override { return MFv(i); }
   double get_MSRdp() const override { return MSRdp; }
   double get_MSRum() const override { return MSRum; }
   double get_MsigmaO() const override { return MsigmaO; }
   double get_MphiO() const override { return MphiO; }
   const Eigen::Array<double,6,1>& get_MSd() const override { return MSd; }
   double get_MSd(int i) const override { return MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv() const override { return MSv; }
   double get_MSv(int i) const override { return MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu() const override { return MSu; }
   double get_MSu(int i) const override { return MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe() const override { return MSe; }
   double get_MSe(int i) const override { return MSe(i); }
   const Eigen::Array<double,4,1>& get_Mhh() const override { return Mhh; }
   double get_Mhh(int i) const override { return Mhh(i); }
   const Eigen::Array<double,4,1>& get_MAh() const override { return MAh; }
   double get_MAh(int i) const override { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MRh() const override { return MRh; }
   double get_MRh(int i) const override { return MRh(i); }
   const Eigen::Array<double,4,1>& get_MHpm() const override { return MHpm; }
   double get_MHpm(int i) const override { return MHpm(i); }
   const Eigen::Array<double,4,1>& get_MChi() const override { return MChi; }
   double get_MChi(int i) const override { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha1() const override { return MCha1; }
   double get_MCha1(int i) const override { return MCha1(i); }
   const Eigen::Array<double,2,1>& get_MCha2() const override { return MCha2; }
   double get_MCha2(int i) const override { return MCha2(i); }
   const Eigen::Array<double,3,1>& get_MFe() const override { return MFe; }
   double get_MFe(int i) const override { return MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd() const override { return MFd; }
   double get_MFd(int i) const override { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const override { return MFu; }
   double get_MFu(int i) const override { return MFu(i); }
   double get_MVWm() const override { return MVWm; }
   double get_MVP() const override { return MVP; }
   double get_MVZ() const override { return MVZ; }

   
   Eigen::Array<double,3,1> get_MChargedHiggs() const override;

   Eigen::Array<double,3,1> get_MPseudoscalarHiggs() const override;

   const Eigen::Matrix<double,6,6>& get_ZD() const override { return ZD; }
   double get_ZD(int i, int k) const override { return ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV() const override { return ZV; }
   double get_ZV(int i, int k) const override { return ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU() const override { return ZU; }
   double get_ZU(int i, int k) const override { return ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE() const override { return ZE; }
   double get_ZE(int i, int k) const override { return ZE(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZH() const override { return ZH; }
   double get_ZH(int i, int k) const override { return ZH(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZA() const override { return ZA; }
   double get_ZA(int i, int k) const override { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZHR() const override { return ZHR; }
   double get_ZHR(int i, int k) const override { return ZHR(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZP() const override { return ZP; }
   double get_ZP(int i, int k) const override { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN1() const override { return ZN1; }
   std::complex<double> get_ZN1(int i, int k) const override { return ZN1(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN2() const override { return ZN2; }
   std::complex<double> get_ZN2(int i, int k) const override { return ZN2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM1() const override { return UM1; }
   std::complex<double> get_UM1(int i, int k) const override { return UM1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP1() const override { return UP1; }
   std::complex<double> get_UP1(int i, int k) const override { return UP1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM2() const override { return UM2; }
   std::complex<double> get_UM2(int i, int k) const override { return UM2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP2() const override { return UP2; }
   std::complex<double> get_UP2(int i, int k) const override { return UP2(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL() const override { return ZEL; }
   std::complex<double> get_ZEL(int i, int k) const override { return ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER() const override { return ZER; }
   std::complex<double> get_ZER(int i, int k) const override { return ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL() const override { return ZDL; }
   std::complex<double> get_ZDL(int i, int k) const override { return ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR() const override { return ZDR; }
   std::complex<double> get_ZDR(int i, int k) const override { return ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL() const override { return ZUL; }
   std::complex<double> get_ZUL(int i, int k) const override { return ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR() const override { return ZUR; }
   std::complex<double> get_ZUR(int i, int k) const override { return ZUR(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const override { return ZZ; }
   double get_ZZ(int i, int k) const override { return ZZ(i,k); }




   double get_mass_matrix_VG() const override;
   void calculate_MVG() override;
   double get_mass_matrix_Glu() const override;
   void calculate_MGlu() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const override;
   void calculate_MFv() override;
   double get_mass_matrix_SRdp() const override;
   void calculate_MSRdp() override;
   double get_mass_matrix_SRum() const override;
   void calculate_MSRum() override;
   double get_mass_matrix_sigmaO() const override;
   void calculate_MsigmaO() override;
   double get_mass_matrix_phiO() const override;
   void calculate_MphiO() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Sd() const override;
   void calculate_MSd() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Sv() const override;
   void calculate_MSv() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Su() const override;
   void calculate_MSu() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Se() const override;
   void calculate_MSe() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_hh() const override;
   void calculate_Mhh() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_Ah() const override;
   void calculate_MAh() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Rh() const override;
   void calculate_MRh() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_Hpm() const override;
   void calculate_MHpm() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_Chi() const override;
   void calculate_MChi() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha1() const override;
   void calculate_MCha1() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha2() const override;
   void calculate_MCha2() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const override;
   void calculate_MFe() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const override;
   void calculate_MFd() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const override;
   void calculate_MFu() override;
   double get_mass_matrix_VWm() const override;
   void calculate_MVWm() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const override;
   void calculate_MVPVZ() override;

   double get_ewsb_eq_hh_1() const override;
   double get_ewsb_eq_hh_2() const override;
   double get_ewsb_eq_hh_3() const override;
   double get_ewsb_eq_hh_4() const override;

   double v() const override;
   double Betax() const override;
   double ThetaW() const override;
   double VEV() const override;


private:
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< mass eigenstate precision
   MRSSM2_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{MRSSM2_info::model_name,
                     &MRSSM2_info::particle_names_getter,
                     &MRSSM2_info::parameter_names_getter}; ///< problems

   void clear_tree_level_parameters();
   void copy_tree_level_masses_to_pole_masses();

   // DR-bar masses
   double MVG{};
   double MGlu{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double MSRdp{};
   double MSRum{};
   double MsigmaO{};
   double MphiO{};
   Eigen::Array<double,6,1> MSd{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,3,1> MSv{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,6,1> MSu{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,6,1> MSe{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,4,1> Mhh{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,4,1> MAh{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MRh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,4,1> MHpm{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,4,1> MChi{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MCha1{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MCha2{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   double MVWm{};
   double MVP{};
   double MVZ{};

   // DR-bar mixing matrices
   Eigen::Matrix<double,6,6> ZD{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,3,3> ZV{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,6,6> ZU{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,6,6> ZE{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,4,4> ZH{Eigen::Matrix<double,4,4>::Zero()};
   Eigen::Matrix<double,4,4> ZA{Eigen::Matrix<double,4,4>::Zero()};
   Eigen::Matrix<double,2,2> ZHR{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,4,4> ZP{Eigen::Matrix<double,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,4,4> ZN1{Eigen::Matrix<std::complex<double>,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,4,4> ZN2{Eigen::Matrix<std::complex<double>,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM1{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP1{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM2{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP2{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZEL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZER{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZUL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZUR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};

   // phases

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const MRSSM2_mass_eigenstates_decoupling_scheme&);

} // namespace flexiblesusy

#undef SUPER

#endif
