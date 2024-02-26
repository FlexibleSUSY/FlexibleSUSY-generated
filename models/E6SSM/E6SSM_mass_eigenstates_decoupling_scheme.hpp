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
 * @file E6SSM_mass_eigenstates_decoupling_scheme.hpp
 *
 * @brief Defines model class for Stöckinger/Kotlarski decoupling scheme.
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef E6SSM_MASS_EIGENSTATES_DECOUPLING_SCHEME_H
#define E6SSM_MASS_EIGENSTATES_DECOUPLING_SCHEME_H

#include "E6SSM_info.hpp"
#include "E6SSM_physical.hpp"
#include "E6SSM_soft_parameters.hpp"
#include "E6SSM_mass_eigenstates_interface.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>

#include <Eigen/Core>

#define SUPER(p) E6SSM_soft_parameters::p

namespace flexiblesusy {

namespace standard_model {
class Standard_model;
}

struct E6SSM_input_parameters;
class E6SSM_mass_eigenstates;

/**
 * @class E6SSM_mass_eigenstates_decoupling_scheme
 *
 * @brief model class with routines for determing masses and mixings
 * and EWSB in the Stöckinger/Kotlarski decoupling scheme
 */
class E6SSM_mass_eigenstates_decoupling_scheme
   : private E6SSM_soft_parameters
   , public E6SSM_mass_eigenstates_interface
{
public:
   explicit E6SSM_mass_eigenstates_decoupling_scheme(const E6SSM_input_parameters& input_ = E6SSM_input_parameters());
   explicit E6SSM_mass_eigenstates_decoupling_scheme(const E6SSM_mass_eigenstates&);
   E6SSM_mass_eigenstates_decoupling_scheme(const E6SSM_mass_eigenstates_decoupling_scheme&) = default;
   E6SSM_mass_eigenstates_decoupling_scheme(E6SSM_mass_eigenstates_decoupling_scheme&&) = default;
   virtual ~E6SSM_mass_eigenstates_decoupling_scheme() = default;
   E6SSM_mass_eigenstates_decoupling_scheme& operator=(const E6SSM_mass_eigenstates_decoupling_scheme&) = default;
   E6SSM_mass_eigenstates_decoupling_scheme& operator=(E6SSM_mass_eigenstates_decoupling_scheme&&) = default;

   std::unique_ptr<E6SSM_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static constexpr int number_of_ewsb_equations = 3;

   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void fill_from(const standard_model::Standard_model&);
   void fill_from(const E6SSM_mass_eigenstates&);
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
   const E6SSM_input_parameters& get_input_parameters() const override;
   E6SSM_input_parameters& get_input_parameters() override;
   Eigen::ArrayXd get_extra_parameters() const override;
   const E6SSM_physical& get_physical() const override;
   E6SSM_physical& get_physical() override;
   const Problems& get_problems() const override;
   Problems& get_problems() override;
   void set_tree_level_masses(const Eigen::ArrayXd&) override;
   void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) override;
   void set_extra_parameters(const Eigen::ArrayXd&) override;
   void set_physical(const E6SSM_physical&) override;
   void clear_problems() override;

   const Eigen::Matrix<double,3,3>& get_Yd() const override { return SUPER(Yd); }
   double get_Yd(int i, int k) const override { return SUPER(Yd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Ye() const override { return SUPER(Ye); }
   double get_Ye(int i, int k) const override { return SUPER(Ye(i,k)); }
   const Eigen::Matrix<double,3,3>& get_Kappa() const override { return SUPER(Kappa); }
   double get_Kappa(int i, int k) const override { return SUPER(Kappa(i,k)); }
   const Eigen::Matrix<double,2,2>& get_Lambda12() const override { return SUPER(Lambda12); }
   double get_Lambda12(int i, int k) const override { return SUPER(Lambda12(i,k)); }
   double get_Lambdax() const override { return SUPER(Lambdax); }
   const Eigen::Matrix<double,3,3>& get_Yu() const override { return SUPER(Yu); }
   double get_Yu(int i, int k) const override { return SUPER(Yu(i,k)); }
   double get_MuPr() const override { return SUPER(MuPr); }
   double get_g1() const override { return SUPER(g1); }
   double get_g2() const override { return SUPER(g2); }
   double get_g3() const override { return SUPER(g3); }
   double get_gN() const override { return SUPER(gN); }
   double get_vd() const override { return SUPER(vd); }
   double get_vu() const override { return SUPER(vu); }
   double get_vs() const override { return SUPER(vs); }
   const Eigen::Matrix<double,3,3>& get_TYd() const override { return SUPER(TYd); }
   double get_TYd(int i, int k) const override { return SUPER(TYd(i,k)); }
   const Eigen::Matrix<double,3,3>& get_TYe() const override { return SUPER(TYe); }
   double get_TYe(int i, int k) const override { return SUPER(TYe(i,k)); }
   const Eigen::Matrix<double,3,3>& get_TKappa() const override { return SUPER(TKappa); }
   double get_TKappa(int i, int k) const override { return SUPER(TKappa(i,k)); }
   const Eigen::Matrix<double,2,2>& get_TLambda12() const override { return SUPER(TLambda12); }
   double get_TLambda12(int i, int k) const override { return SUPER(TLambda12(i,k)); }
   double get_TLambdax() const override { return SUPER(TLambdax); }
   const Eigen::Matrix<double,3,3>& get_TYu() const override { return SUPER(TYu); }
   double get_TYu(int i, int k) const override { return SUPER(TYu(i,k)); }
   double get_BMuPr() const override { return SUPER(BMuPr); }
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
   double get_ms2() const override { return SUPER(ms2); }
   const Eigen::Matrix<double,2,2>& get_mH1I2() const override { return SUPER(mH1I2); }
   double get_mH1I2(int i, int k) const override { return SUPER(mH1I2(i,k)); }
   const Eigen::Matrix<double,2,2>& get_mH2I2() const override { return SUPER(mH2I2); }
   double get_mH2I2(int i, int k) const override { return SUPER(mH2I2(i,k)); }
   const Eigen::Matrix<double,2,2>& get_msI2() const override { return SUPER(msI2); }
   double get_msI2(int i, int k) const override { return SUPER(msI2(i,k)); }
   const Eigen::Matrix<double,3,3>& get_mDx2() const override { return SUPER(mDx2); }
   double get_mDx2(int i, int k) const override { return SUPER(mDx2(i,k)); }
   const Eigen::Matrix<double,3,3>& get_mDxbar2() const override { return SUPER(mDxbar2); }
   double get_mDxbar2(int i, int k) const override { return SUPER(mDxbar2(i,k)); }
   double get_mHp2() const override { return SUPER(mHp2); }
   double get_mHpbar2() const override { return SUPER(mHpbar2); }
   double get_MassB() const override { return SUPER(MassB); }
   double get_MassWB() const override { return SUPER(MassWB); }
   double get_MassG() const override { return SUPER(MassG); }
   double get_MassBp() const override { return SUPER(MassBp); }

   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) override { Yd = Yd_; }
   void set_Yd(int i, int k, const double& value) override { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) override { Ye = Ye_; }
   void set_Ye(int i, int k, const double& value) override { Ye(i,k) = value; }
   void set_Kappa(const Eigen::Matrix<double,3,3>& Kappa_) override { Kappa = Kappa_; }
   void set_Kappa(int i, int k, const double& value) override { Kappa(i,k) = value; }
   void set_Lambda12(const Eigen::Matrix<double,2,2>& Lambda12_) override { Lambda12 = Lambda12_; }
   void set_Lambda12(int i, int k, const double& value) override { Lambda12(i,k) = value; }
   void set_Lambdax(double Lambdax_) override { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) override { Yu = Yu_; }
   void set_Yu(int i, int k, const double& value) override { Yu(i,k) = value; }
   void set_MuPr(double MuPr_) override { MuPr = MuPr_; }
   void set_g1(double g1_) override { g1 = g1_; }
   void set_g2(double g2_) override { g2 = g2_; }
   void set_g3(double g3_) override { g3 = g3_; }
   void set_gN(double gN_) override { gN = gN_; }
   void set_vd(double vd_) override { vd = vd_; }
   void set_vu(double vu_) override { vu = vu_; }
   void set_vs(double vs_) override { vs = vs_; }
   void set_TYd(const Eigen::Matrix<double,3,3>& TYd_) override { TYd = TYd_; }
   void set_TYd(int i, int k, const double& value) override { TYd(i,k) = value; }
   void set_TYe(const Eigen::Matrix<double,3,3>& TYe_) override { TYe = TYe_; }
   void set_TYe(int i, int k, const double& value) override { TYe(i,k) = value; }
   void set_TKappa(const Eigen::Matrix<double,3,3>& TKappa_) override { TKappa = TKappa_; }
   void set_TKappa(int i, int k, const double& value) override { TKappa(i,k) = value; }
   void set_TLambda12(const Eigen::Matrix<double,2,2>& TLambda12_) override { TLambda12 = TLambda12_; }
   void set_TLambda12(int i, int k, const double& value) override { TLambda12(i,k) = value; }
   void set_TLambdax(double TLambdax_) override { TLambdax = TLambdax_; }
   void set_TYu(const Eigen::Matrix<double,3,3>& TYu_) override { TYu = TYu_; }
   void set_TYu(int i, int k, const double& value) override { TYu(i,k) = value; }
   void set_BMuPr(double BMuPr_) override { BMuPr = BMuPr_; }
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
   void set_ms2(double ms2_) override { ms2 = ms2_; }
   void set_mH1I2(const Eigen::Matrix<double,2,2>& mH1I2_) override { mH1I2 = mH1I2_; }
   void set_mH1I2(int i, int k, const double& value) override { mH1I2(i,k) = value; }
   void set_mH2I2(const Eigen::Matrix<double,2,2>& mH2I2_) override { mH2I2 = mH2I2_; }
   void set_mH2I2(int i, int k, const double& value) override { mH2I2(i,k) = value; }
   void set_msI2(const Eigen::Matrix<double,2,2>& msI2_) override { msI2 = msI2_; }
   void set_msI2(int i, int k, const double& value) override { msI2(i,k) = value; }
   void set_mDx2(const Eigen::Matrix<double,3,3>& mDx2_) override { mDx2 = mDx2_; }
   void set_mDx2(int i, int k, const double& value) override { mDx2(i,k) = value; }
   void set_mDxbar2(const Eigen::Matrix<double,3,3>& mDxbar2_) override { mDxbar2 = mDxbar2_; }
   void set_mDxbar2(int i, int k, const double& value) override { mDxbar2(i,k) = value; }
   void set_mHp2(double mHp2_) override { mHp2 = mHp2_; }
   void set_mHpbar2(double mHpbar2_) override { mHpbar2 = mHpbar2_; }
   void set_MassB(double MassB_) override { MassB = MassB_; }
   void set_MassWB(double MassWB_) override { MassWB = MassWB_; }
   void set_MassG(double MassG_) override { MassG = MassG_; }
   void set_MassBp(double MassBp_) override { MassBp = MassBp_; }

   double get_MVG() const override { return MVG; }
   double get_MGlu() const override { return MGlu; }
   const Eigen::Array<double,3,1>& get_MFv() const override { return MFv; }
   double get_MFv(int i) const override { return MFv(i); }
   double get_MChaP() const override { return MChaP; }
   const Eigen::Array<double,6,1>& get_MSd() const override { return MSd; }
   double get_MSd(int i) const override { return MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv() const override { return MSv; }
   double get_MSv(int i) const override { return MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu() const override { return MSu; }
   double get_MSu(int i) const override { return MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe() const override { return MSe; }
   double get_MSe(int i) const override { return MSe(i); }
   const Eigen::Array<double,6,1>& get_MSDX() const override { return MSDX; }
   double get_MSDX(int i) const override { return MSDX(i); }
   const Eigen::Array<double,3,1>& get_Mhh() const override { return Mhh; }
   double get_Mhh(int i) const override { return Mhh(i); }
   const Eigen::Array<double,3,1>& get_MAh() const override { return MAh; }
   double get_MAh(int i) const override { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm() const override { return MHpm; }
   double get_MHpm(int i) const override { return MHpm(i); }
   const Eigen::Array<double,6,1>& get_MChi() const override { return MChi; }
   double get_MChi(int i) const override { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha() const override { return MCha; }
   double get_MCha(int i) const override { return MCha(i); }
   const Eigen::Array<double,3,1>& get_MFe() const override { return MFe; }
   double get_MFe(int i) const override { return MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd() const override { return MFd; }
   double get_MFd(int i) const override { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const override { return MFu; }
   double get_MFu(int i) const override { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFDX() const override { return MFDX; }
   double get_MFDX(int i) const override { return MFDX(i); }
   const Eigen::Array<double,4,1>& get_MSHI0() const override { return MSHI0; }
   double get_MSHI0(int i) const override { return MSHI0(i); }
   const Eigen::Array<double,4,1>& get_MSHIp() const override { return MSHIp; }
   double get_MSHIp(int i) const override { return MSHIp(i); }
   const Eigen::Array<double,2,1>& get_MChaI() const override { return MChaI; }
   double get_MChaI(int i) const override { return MChaI(i); }
   const Eigen::Array<double,4,1>& get_MChiI() const override { return MChiI; }
   double get_MChiI(int i) const override { return MChiI(i); }
   const Eigen::Array<double,2,1>& get_MSSI0() const override { return MSSI0; }
   double get_MSSI0(int i) const override { return MSSI0(i); }
   const Eigen::Array<double,2,1>& get_MFSI() const override { return MFSI; }
   double get_MFSI(int i) const override { return MFSI(i); }
   const Eigen::Array<double,2,1>& get_MSHp0() const override { return MSHp0; }
   double get_MSHp0(int i) const override { return MSHp0(i); }
   const Eigen::Array<double,2,1>& get_MSHpp() const override { return MSHpp; }
   double get_MSHpp(int i) const override { return MSHpp(i); }
   const Eigen::Array<double,2,1>& get_MChiP() const override { return MChiP; }
   double get_MChiP(int i) const override { return MChiP(i); }
   double get_MVWm() const override { return MVWm; }
   double get_MVP() const override { return MVP; }
   double get_MVZ() const override { return MVZ; }
   double get_MVZp() const override { return MVZp; }

   
   Eigen::Array<double,1,1> get_MChargedHiggs() const override;

   Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const override;

   const Eigen::Matrix<double,6,6>& get_ZD() const override { return ZD; }
   double get_ZD(int i, int k) const override { return ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV() const override { return ZV; }
   double get_ZV(int i, int k) const override { return ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU() const override { return ZU; }
   double get_ZU(int i, int k) const override { return ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE() const override { return ZE; }
   double get_ZE(int i, int k) const override { return ZE(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZDX() const override { return ZDX; }
   double get_ZDX(int i, int k) const override { return ZDX(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZH() const override { return ZH; }
   double get_ZH(int i, int k) const override { return ZH(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZA() const override { return ZA; }
   double get_ZA(int i, int k) const override { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const override { return ZP; }
   double get_ZP(int i, int k) const override { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,6,6>& get_ZN() const override { return ZN; }
   std::complex<double> get_ZN(int i, int k) const override { return ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM() const override { return UM; }
   std::complex<double> get_UM(int i, int k) const override { return UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP() const override { return UP; }
   std::complex<double> get_UP(int i, int k) const override { return UP(i,k); }
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
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXL() const override { return ZDXL; }
   std::complex<double> get_ZDXL(int i, int k) const override { return ZDXL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDXR() const override { return ZDXR; }
   std::complex<double> get_ZDXR(int i, int k) const override { return ZDXR(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHI0() const override { return UHI0; }
   double get_UHI0(int i, int k) const override { return UHI0(i,k); }
   const Eigen::Matrix<double,4,4>& get_UHIp() const override { return UHIp; }
   double get_UHIp(int i, int k) const override { return UHIp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZMI() const override { return ZMI; }
   std::complex<double> get_ZMI(int i, int k) const override { return ZMI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZPI() const override { return ZPI; }
   std::complex<double> get_ZPI(int i, int k) const override { return ZPI(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZNI() const override { return ZNI; }
   std::complex<double> get_ZNI(int i, int k) const override { return ZNI(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZSSI() const override { return ZSSI; }
   double get_ZSSI(int i, int k) const override { return ZSSI(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZFSI() const override { return ZFSI; }
   std::complex<double> get_ZFSI(int i, int k) const override { return ZFSI(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHp0() const override { return UHp0; }
   double get_UHp0(int i, int k) const override { return UHp0(i,k); }
   const Eigen::Matrix<double,2,2>& get_UHpp() const override { return UHpp; }
   double get_UHpp(int i, int k) const override { return UHpp(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_ZNp() const override { return ZNp; }
   std::complex<double> get_ZNp(int i, int k) const override { return ZNp(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZZ() const override { return ZZ; }
   double get_ZZ(int i, int k) const override { return ZZ(i,k); }

   void set_PhaseGlu(std::complex<double> PhaseGlu_) override { PhaseGlu = PhaseGlu_; }
   std::complex<double> get_PhaseGlu() const override { return PhaseGlu; }
   void set_PhaseFHpup(std::complex<double> PhaseFHpup_) override { PhaseFHpup = PhaseFHpup_; }
   std::complex<double> get_PhaseFHpup() const override { return PhaseFHpup; }



   double get_mass_matrix_VG() const override;
   void calculate_MVG() override;
   double get_mass_matrix_Glu() const override;
   void calculate_MGlu() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const override;
   void calculate_MFv() override;
   double get_mass_matrix_ChaP() const override;
   void calculate_MChaP() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Sd() const override;
   void calculate_MSd() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Sv() const override;
   void calculate_MSv() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Su() const override;
   void calculate_MSu() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Se() const override;
   void calculate_MSe() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_SDX() const override;
   void calculate_MSDX() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_hh() const override;
   void calculate_Mhh() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Ah() const override;
   void calculate_MAh() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Hpm() const override;
   void calculate_MHpm() override;
   Eigen::Matrix<double,6,6> get_mass_matrix_Chi() const override;
   void calculate_MChi() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha() const override;
   void calculate_MCha() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const override;
   void calculate_MFe() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const override;
   void calculate_MFd() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const override;
   void calculate_MFu() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_FDX() const override;
   void calculate_MFDX() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_SHI0() const override;
   void calculate_MSHI0() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_SHIp() const override;
   void calculate_MSHIp() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_ChaI() const override;
   void calculate_MChaI() override;
   Eigen::Matrix<double,4,4> get_mass_matrix_ChiI() const override;
   void calculate_MChiI() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_SSI0() const override;
   void calculate_MSSI0() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_FSI() const override;
   void calculate_MFSI() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_SHp0() const override;
   void calculate_MSHp0() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_SHpp() const override;
   void calculate_MSHpp() override;
   Eigen::Matrix<double,2,2> get_mass_matrix_ChiP() const override;
   void calculate_MChiP() override;
   double get_mass_matrix_VWm() const override;
   void calculate_MVWm() override;
   Eigen::Matrix<double,3,3> get_mass_matrix_VPVZVZp() const override;
   void calculate_MVPVZVZp() override;

   double get_ewsb_eq_hh_1() const override;
   double get_ewsb_eq_hh_2() const override;
   double get_ewsb_eq_hh_3() const override;

   double v() const override;
   double Betax() const override;
   double ThetaW() const override;
   double ThetaWp() const override;
   double VEV() const override;


private:
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< mass eigenstate precision
   E6SSM_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{E6SSM_info::model_name,
                     &E6SSM_info::particle_names_getter,
                     &E6SSM_info::parameter_names_getter}; ///< problems

   void clear_tree_level_parameters();
   void copy_tree_level_masses_to_pole_masses();

   // DR-bar masses
   double MVG{};
   double MGlu{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double MChaP{};
   Eigen::Array<double,6,1> MSd{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,3,1> MSv{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,6,1> MSu{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,6,1> MSe{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,6,1> MSDX{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,3,1> Mhh{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MAh{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,2,1> MHpm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,6,1> MChi{Eigen::Array<double,6,1>::Zero()};
   Eigen::Array<double,2,1> MCha{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFDX{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,4,1> MSHI0{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,4,1> MSHIp{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MChaI{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,4,1> MChiI{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MSSI0{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MFSI{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSHp0{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSHpp{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MChiP{Eigen::Array<double,2,1>::Zero()};
   double MVWm{};
   double MVP{};
   double MVZ{};
   double MVZp{};

   // DR-bar mixing matrices
   Eigen::Matrix<double,6,6> ZD{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,3,3> ZV{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,6,6> ZU{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,6,6> ZE{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,6,6> ZDX{Eigen::Matrix<double,6,6>::Zero()};
   Eigen::Matrix<double,3,3> ZH{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> ZA{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,6,6> ZN{Eigen::Matrix<std::complex<double>,6,6>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZEL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZER{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZUL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZUR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDXL{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> ZDXR{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,4,4> UHI0{Eigen::Matrix<double,4,4>::Zero()};
   Eigen::Matrix<double,4,4> UHIp{Eigen::Matrix<double,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> ZMI{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> ZPI{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,4,4> ZNI{Eigen::Matrix<std::complex<double>,4,4>::Zero()};
   Eigen::Matrix<double,2,2> ZSSI{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> ZFSI{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<double,2,2> UHp0{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> UHpp{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> ZNp{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<double,3,3> ZZ{Eigen::Matrix<double,3,3>::Zero()};

   // phases
   std::complex<double> PhaseGlu{1.,0.};
   std::complex<double> PhaseFHpup{1.,0.};

   // extra parameters

};

std::ostream& operator<<(std::ostream&, const E6SSM_mass_eigenstates_decoupling_scheme&);

} // namespace flexiblesusy

#undef SUPER

#endif
