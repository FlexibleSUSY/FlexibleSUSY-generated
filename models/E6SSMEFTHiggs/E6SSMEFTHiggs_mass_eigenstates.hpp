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
 * @file E6SSMEFTHiggs_mass_eigenstates.hpp
 *
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef E6SSMEFTHiggs_MASS_EIGENSTATES_H
#define E6SSMEFTHiggs_MASS_EIGENSTATES_H

#include "E6SSMEFTHiggs_info.hpp"
#include "E6SSMEFTHiggs_physical.hpp"
#include "E6SSMEFTHiggs_soft_parameters.hpp"
#include "E6SSMEFTHiggs_mass_eigenstates_interface.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "problems.hpp"

#include <iosfwd>
#include <memory>
#include <string>

#include <Eigen/Core>

#define SUPER(p) E6SSMEFTHiggs_soft_parameters::p

namespace flexiblesusy {

class E6SSMEFTHiggs_ewsb_solver_interface;
/**
 * @class E6SSMEFTHiggs_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class E6SSMEFTHiggs_mass_eigenstates
   : public E6SSMEFTHiggs_soft_parameters
   , public E6SSMEFTHiggs_mass_eigenstates_interface
{
public:
   explicit E6SSMEFTHiggs_mass_eigenstates(const E6SSMEFTHiggs_input_parameters& input_ = E6SSMEFTHiggs_input_parameters());
   E6SSMEFTHiggs_mass_eigenstates(const E6SSMEFTHiggs_mass_eigenstates&) = default;
   E6SSMEFTHiggs_mass_eigenstates(E6SSMEFTHiggs_mass_eigenstates&&) = default;
   virtual ~E6SSMEFTHiggs_mass_eigenstates() = default;
   E6SSMEFTHiggs_mass_eigenstates& operator=(const E6SSMEFTHiggs_mass_eigenstates&) = default;
   E6SSMEFTHiggs_mass_eigenstates& operator=(E6SSMEFTHiggs_mass_eigenstates&&) = default;
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   std::unique_ptr<E6SSMEFTHiggs_mass_eigenstates_interface> clone() const override;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 3;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   virtual void clear() override;
   void clear_DRbar_parameters();
   Eigen::ArrayXd get_DRbar_masses() const;
   Eigen::ArrayXd get_DRbar_masses_and_mixings() const;
   void do_calculate_sm_pole_masses(bool);
   bool do_calculate_sm_pole_masses() const;
   void do_calculate_bsm_pole_masses(bool);
   bool do_calculate_bsm_pole_masses() const;
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_DRbar_masses();
   void reorder_pole_masses();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(int);
   void set_loop_corrections(const Loop_corrections&);
   const Loop_corrections& get_loop_corrections() const;
   void set_threshold_corrections(const Threshold_corrections&);
   const Threshold_corrections& get_threshold_corrections() const;
   void set_DRbar_masses(const Eigen::ArrayXd&);
   void set_DRbar_masses_and_mixings(const Eigen::ArrayXd&);
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   void set_ewsb_solver(const std::shared_ptr<E6SSMEFTHiggs_ewsb_solver_interface>&);
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level
   int solve_ewsb_tree_level_custom();
   
   virtual void calculate_spectrum();
   std::string name() const;
   void run_to(double scale, double eps = -1.0) override;
   void print(std::ostream&) const override;
   void set_precision(double);
   double get_precision() const;

   // mass_eigenstates_interface functions

   void calculate_tree_level_mass_spectrum() override;
   void calculate_pole_mass_spectrum() override;
   void calculate_mass_spectrum() override;
   int solve_ewsb_equations_tree_level() override;
   int solve_ewsb_equations() override;
   Eigen::ArrayXd get_tree_level_masses() const override;
   Eigen::ArrayXd get_tree_level_masses_and_mixings() const override;
   const E6SSMEFTHiggs_input_parameters& get_input_parameters() const override;
   E6SSMEFTHiggs_input_parameters& get_input_parameters() override;
   Eigen::ArrayXd get_extra_parameters() const override;
   const E6SSMEFTHiggs_physical& get_physical() const override;
   E6SSMEFTHiggs_physical& get_physical() override;
   const Problems& get_problems() const override;
   Problems& get_problems() override;
   void set_tree_level_masses(const Eigen::ArrayXd&) override;
   void set_tree_level_masses_and_mixings(const Eigen::ArrayXd&) override;
   void set_extra_parameters(const Eigen::ArrayXd&) override;
   void set_physical(const E6SSMEFTHiggs_physical&) override;
   void clear_problems() override;

   double get_lsp(E6SSMEFTHiggs_info::Particles&) const;

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

   std::complex<double> CpUSdconjUSdVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSdconjUSdVZpVZp(int gO1, int gO2) const;
   double CpUSdconjUSdconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSdconjHpmconjUSd(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSdSHp0conjUSdconjSHp0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSdSHppconjUSdconjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSdSSI0conjUSdconjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhAhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSdSvconjUSdconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpChaFuconjUSdPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChaFuconjUSdPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpChiFdconjUSdPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFdconjUSdPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpUSdSHI0conjUSdconjSHI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSdSHIpconjUSdconjSHIp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSdconjUSdconjSdSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSdconjUSdconjSDXSDX(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSdconjUSdconjSuSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSdSeconjUSdconjSe(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpHpmSuconjUSd(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhSdconjUSd(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSdconjUSd(int gI2, int gI1, int gO2) const;
   std::complex<double> CpGluFdconjUSdPR(int gI2, int gO2) const;
   std::complex<double> CpGluFdconjUSdPL(int gI2, int gO1) const;
   std::complex<double> CpSdconjUSdVG(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSdVP(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSdVZ(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSdVZp(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSdVWm(int gI2, int gO2) const;
   std::complex<double> CpUSvconjUSvVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSvconjUSvVZpVZp(int gO1, int gO2) const;
   double CpUSvconjUSvconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSvconjHpmconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHp0USvconjSHp0conjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHppUSvconjSHppconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   double CpSSI0USvconjSSI0conjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpbarChaFeconjUSvPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaFeconjUSvPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSeconjHpmconjUSv(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhAhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpSvUSvconjSvconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CphhSvconjUSv(int gI2, int gI1, int gO2) const;
   double CpChiFvconjUSvPR(int , int , int ) const;
   std::complex<double> CpChiFvconjUSvPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpSHI0USvconjSHI0conjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHIpUSvconjSHIpconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSdUSvconjSdconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSvconjSDXconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSvconjSeconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSuUSvconjSuconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSvconjUSvVZ(int gI2, int gO2) const;
   std::complex<double> CpSvconjUSvVZp(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSvconjVWm(int gI2, int gO2) const;
   std::complex<double> CpUSuconjUSuVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSuconjUSuVZpVZp(int gO1, int gO2) const;
   double CpUSuconjUSuconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSuconjHpmconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHp0USuconjSHp0conjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHppUSuconjSHppconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   double CpSSI0USuconjSSI0conjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpbarChaFdconjUSuPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaFdconjUSuPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSdconjHpmconjUSu(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhAhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSuSvconjUSuconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiFuconjUSuPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFuconjUSuPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpSHI0USuconjSHI0conjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHIpUSuconjSHIpconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSuconjSeconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSuconjUSuconjSdSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSuconjUSuconjSDXSDX(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSuconjUSuconjSuSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpAhSuconjUSu(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSuconjUSu(int gI2, int gI1, int gO2) const;
   std::complex<double> CpGluFuconjUSuPR(int gI2, int gO2) const;
   std::complex<double> CpGluFuconjUSuPL(int gI2, int gO1) const;
   std::complex<double> CpSdconjUSuconjVWm(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVG(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVP(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVZ(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVZp(int gI2, int gO2) const;
   std::complex<double> CpUSeconjUSeVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSeconjUSeVZpVZp(int gO1, int gO2) const;
   double CpUSeconjUSeconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSeconjHpmconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSeSHp0conjUSeconjSHp0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSeSHppconjUSeconjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSeSSI0conjUSeconjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhAhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSeSvconjUSeconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpHpmSvconjUSe(int gI2, int gI1, int gO2) const;
   double CpChaFvconjUSePR(int , int , int ) const;
   std::complex<double> CpChaFvconjUSePL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpChiFeconjUSePR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFeconjUSePL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpUSeSHI0conjUSeconjSHI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSeSHIpconjUSeconjSHIp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSdUSeconjSdconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSeconjSDXconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSeconjSeconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSeSuconjUSeconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhSeconjUSe(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSeconjUSe(int gI2, int gI1, int gO2) const;
   std::complex<double> CpSvconjUSeVWm(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSeVP(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSeVZ(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSeVZp(int gI2, int gO2) const;
   std::complex<double> CpUSDXconjUSDXVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSDXconjUSDXVZpVZp(int gO1, int gO2) const;
   std::complex<double> CpHpmUSDXconjHpmconjUSDX(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSDXSHp0conjUSDXconjSHp0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSDXSHppconjUSDXconjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSDXSSI0conjUSDXconjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhAhUSDXconjUSDX(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSDXconjUSDX(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSDXSvconjUSDXconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiFDXconjUSDXPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFDXconjUSDXPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpUSDXSHI0conjUSDXconjSHI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSDXSHIpconjUSDXconjSHIp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSdSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSDXSDX(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSDXconjUSDXconjSuSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSDXSeconjUSDXconjSe(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhSDXconjUSDX(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSDXconjUSDX(int gI2, int gI1, int gO2) const;
   std::complex<double> CpGluFDXconjUSDXPR(int gI2, int gO2) const;
   std::complex<double> CpGluFDXconjUSDXPL(int gI2, int gO1) const;
   std::complex<double> CpSDXconjUSDXVG(int gI2, int gO2) const;
   std::complex<double> CpSDXconjUSDXVP(int gI2, int gO2) const;
   std::complex<double> CpSDXconjUSDXVZ(int gI2, int gO2) const;
   std::complex<double> CpSDXconjUSDXVZp(int gI2, int gO2) const;
   std::complex<double> CpbargWmgWmUhh(int gO1) const;
   std::complex<double> CpbargWmCgWmCUhh(int gO1) const;
   std::complex<double> CpbargZgZUhh(int gO1) const;
   std::complex<double> CpbargZpgZUhh(int gO1) const;
   std::complex<double> CpbargZpgZpUhh(int gO1) const;
   std::complex<double> CpUhhVZVZ(int gO2) const;
   std::complex<double> CpUhhVZVZp(int gO2) const;
   std::complex<double> CpUhhVZpVZp(int gO2) const;
   std::complex<double> CpUhhconjVWmVWm(int gO2) const;
   std::complex<double> CpUhhUhhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhVZpVZp(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSHp0conjSHp0(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSHppconjSHpp(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSSI0conjSSI0(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhHpmconjHpm(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSHp0conjSHp0(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSHppconjSHpp(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSSI0conjSSI0(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarChaChaUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaChaUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarChaIChaIUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaIChaIUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUhhUhhSvconjSv(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpAhAhUhh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhhhUhh(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhhhUhh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpUhhSvconjSv(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFDXFDXUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFDXFDXUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUhhUhhSHI0conjSHI0(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSHIpconjSHIp(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhSHI0conjSHI0(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSHIpconjSHIp(int gO2, int gI2, int gI1) const;
   std::complex<double> CpChiIChiIUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiIChiIUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUhhUhhSdconjSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSDXconjSDX(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSeconjSe(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSuconjSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhSdconjSd(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSDXconjSDX(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSeconjSe(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSuconjSu(int gO2, int gI2, int gI1) const;
   std::complex<double> CpChiChiUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiChiUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUhhHpmconjVWm(int gO2, int gI2) const;
   std::complex<double> CpAhUhhVZ(int gI2, int gO2) const;
   std::complex<double> CpAhUhhVZp(int gI2, int gO2) const;
   std::complex<double> CpbargWmgWmUAh(int gO1) const;
   std::complex<double> CpbargWmCgWmCUAh(int gO1) const;
   std::complex<double> CpUAhUAhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhVZpVZp(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSHp0conjSHp0(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSHppconjSHpp(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSSI0conjSSI0(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhHpmconjHpm(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarChaChaUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaChaUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarChaIChaIUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaIChaIUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSvconjSv(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpAhAhUAh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhUAhhh(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUAhhhhh(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFDXFDXUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFDXFDXUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUAhUAhSHI0conjSHI0(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSHIpconjSHIp(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhSHI0conjSHI0(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhSHIpconjSHIp(int gO2, int gI2, int gI1) const;
   std::complex<double> CpChiIChiIUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiIChiIUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUAhUAhSdconjSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSDXconjSDX(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSeconjSe(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSuconjSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhSdconjSd(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhSDXconjSDX(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhSeconjSe(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhSuconjSu(int gO2, int gI2, int gI1) const;
   std::complex<double> CpChiChiUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiChiUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUAhHpmconjVWm(int gO2, int gI2) const;
   std::complex<double> CpUAhhhVZ(int gO2, int gI2) const;
   std::complex<double> CpUAhhhVZp(int gO2, int gI2) const;
   std::complex<double> CpbargWmgZUHpm(int gO2) const;
   std::complex<double> CpbargZgWmconjUHpm(int gO1) const;
   std::complex<double> CpbargWmCgZconjUHpm(int gO1) const;
   std::complex<double> CpbargZgWmCUHpm(int gO2) const;
   std::complex<double> CpbargWmgZpUHpm(int gO2) const;
   std::complex<double> CpbargZpgWmconjUHpm(int gO1) const;
   std::complex<double> CpbargWmCgZpconjUHpm(int gO1) const;
   std::complex<double> CpbargZpgWmCUHpm(int gO2) const;
   std::complex<double> CpconjUHpmVPVWm(int gO2) const;
   std::complex<double> CpconjUHpmVWmVZ(int gO2) const;
   std::complex<double> CpconjUHpmVWmVZp(int gO2) const;
   std::complex<double> CpUHpmconjUHpmVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUHpmconjUHpmVZpVZp(int gO1, int gO2) const;
   std::complex<double> CpUHpmconjUHpmconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUHpmconjHpmconjUHpm(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUHpmSHp0conjUHpmconjSHp0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSHppconjUHpmconjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSSI0conjUHpmconjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHppconjUHpmconjSHp0(int gI2, int gO2, int gI1) const;
   std::complex<double> CpAhHpmconjUHpm(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhHpmconjUHpm(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhAhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUHpmSvconjUHpmconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarFuFdconjUHpmPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFdconjUHpmPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFvFeconjUHpmPR(int gI1, int gI2, int gO2) const;
   double CpbarFvFeconjUHpmPL(int , int , int ) const;
   std::complex<double> CpSeconjUHpmconjSv(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUHpmSHI0conjUHpmconjSHI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSHIpconjUHpmconjSHIp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHIpconjUHpmconjSHI0(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUHpmSdconjUHpmconjSd(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSDXconjUHpmconjSDX(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSeconjUHpmconjSe(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSuconjUHpmconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiChaconjUHpmPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiChaconjUHpmPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSdconjUHpmconjSu(int gI2, int gO2, int gI1) const;
   std::complex<double> CpHpmconjUHpmVP(int gI2, int gO2) const;
   std::complex<double> CpHpmconjUHpmVZ(int gI2, int gO2) const;
   std::complex<double> CpHpmconjUHpmVZp(int gI2, int gO2) const;
   std::complex<double> CpAhconjUHpmVWm(int gI2, int gO2) const;
   std::complex<double> CphhconjUHpmVWm(int gI2, int gO2) const;
   std::complex<double> CpUSHI0conjUSHI0VZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSHI0conjUSHI0VZpVZp(int gO1, int gO2) const;
   double CpUSHI0conjUSHI0conjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSHI0conjHpmconjUSHI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHI0SHp0conjUSHI0conjSHp0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSHI0SHppconjUSHI0conjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSHI0SSI0conjUSHI0conjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChaChaIconjUSHI0PR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChaChaIconjUSHI0PL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSHIpconjHpmconjUSHI0(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhAhUSHI0conjUSHI0(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSHI0conjUSHI0(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSHI0SvconjUSHI0conjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHI0USHI0conjSHI0conjUSHI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHI0SHIpconjUSHI0conjSHIp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhSHI0conjUSHI0(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSHI0conjUSHI0(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiChiIconjUSHI0PR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiChiIconjUSHI0PL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpSdUSHI0conjSdconjUSHI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSHI0conjSDXconjUSHI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSHI0conjSeconjUSHI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHI0SuconjUSHI0conjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHI0conjUSHI0VZ(int gI2, int gO2) const;
   std::complex<double> CpSHI0conjUSHI0VZp(int gI2, int gO2) const;
   std::complex<double> CpSHIpconjUSHI0conjVWm(int gI2, int gO2) const;
   std::complex<double> CpUSHIpconjUSHIpVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSHIpconjUSHIpVZpVZp(int gO1, int gO2) const;
   double CpUSHIpconjUSHIpconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSHIpconjHpmconjUSHIp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHIpSHp0conjUSHIpconjSHp0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSHIpSHppconjUSHIpconjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSHIpSSI0conjUSHIpconjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpAhAhUSHIpconjUSHIp(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSHIpconjUSHIp(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSHIpSvconjUSHIpconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHI0USHIpconjSHI0conjUSHIp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHIpUSHIpconjSHIpconjUSHIp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpHpmSHI0conjUSHIp(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiIChaconjUSHIpPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiIChaconjUSHIpPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhSHIpconjUSHIp(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSHIpconjUSHIp(int gI2, int gI1, int gO2) const;
   std::complex<double> CpSdUSHIpconjSdconjUSHIp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSHIpconjSDXconjUSHIp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSHIpconjSeconjUSHIp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHIpSuconjUSHIpconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiChaIconjUSHIpPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiChaIconjUSHIpPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSHI0conjUSHIpVWm(int gI2, int gO2) const;
   std::complex<double> CpSHIpconjUSHIpVP(int gI2, int gO2) const;
   std::complex<double> CpSHIpconjUSHIpVZ(int gI2, int gO2) const;
   std::complex<double> CpSHIpconjUSHIpVZp(int gI2, int gO2) const;
   double CpUSSI0conjUSSI0VZVZ(int gO1, int gO2) const;
   double CpUSSI0conjUSSI0VZpVZp(int gO1, int gO2) const;
   std::complex<double> CpHpmUSSI0conjHpmconjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHp0USSI0conjSHp0conjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHppUSSI0conjSHppconjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSSI0USSI0conjSSI0conjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CphhSSI0conjUSSI0(int gI2, int gI1, int gO2) const;
   double CpChiFSIconjUSSI0PR(int , int , int ) const;
   std::complex<double> CpChiFSIconjUSSI0PL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhAhUSSI0conjUSSI0(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSSI0conjUSSI0(int gI1, int gI2, int gO1, int gO2) const;
   double CpUSSI0SvconjUSSI0conjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHI0USSI0conjSHI0conjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHIpUSSI0conjSHIpconjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSdUSSI0conjSdconjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSSI0conjSDXconjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSSI0conjSeconjUSSI0(int gI1, int gO1, int gI2, int gO2) const;
   double CpUSSI0SuconjUSSI0conjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSSI0conjUSSI0VZ(int gI2, int gO2) const;
   std::complex<double> CpSSI0conjUSSI0VZp(int gI2, int gO2) const;
   std::complex<double> CpUSHp0conjUSHp0VZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSHp0conjUSHp0VZpVZp(int gO1, int gO2) const;
   std::complex<double> CpUSHp0conjUSHp0conjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSHp0conjHpmconjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHp0USHp0conjSHp0conjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHp0SHppconjUSHp0conjSHpp(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUSHp0SSI0conjUSHp0conjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHppconjHpmconjUSHp0(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSHp0conjUSHp0(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiChiPconjUSHp0PR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiChiPconjUSHp0PL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpbarChaChaPconjUSHp0PR(int gI1, int gO2) const;
   std::complex<double> CpbarChaChaPconjUSHp0PL(int gI1, int gO1) const;
   std::complex<double> CpAhAhUSHp0conjUSHp0(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSHp0conjUSHp0(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSHp0SvconjUSHp0conjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHI0USHp0conjSHI0conjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHIpUSHp0conjSHIpconjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSdUSHp0conjSdconjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSHp0conjSDXconjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSHp0conjSeconjUSHp0(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHp0SuconjUSHp0conjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHp0conjUSHp0VZ(int gI2, int gO2) const;
   std::complex<double> CpSHp0conjUSHp0VZp(int gI2, int gO2) const;
   std::complex<double> CpSHppconjUSHp0conjVWm(int gI2, int gO2) const;
   std::complex<double> CpUSHppconjUSHppVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUSHppconjUSHppVZpVZp(int gO1, int gO2) const;
   std::complex<double> CpUSHppconjUSHppconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpHpmUSHppconjHpmconjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHp0USHppconjSHp0conjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHppUSHppconjSHppconjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHppSSI0conjUSHppconjSSI0(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpHpmSHp0conjUSHpp(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiPChaconjUSHppPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiPChaconjUSHppPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CphhSHppconjUSHpp(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhAhUSHppconjUSHpp(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSHppconjUSHpp(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUSHppSvconjUSHppconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSHI0USHppconjSHI0conjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSHIpUSHppconjSHIpconjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSdUSHppconjSdconjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSDXUSHppconjSDXconjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSHppconjSeconjUSHpp(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSHppSuconjUSHppconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiChaPconjUSHppPR(int gI1, int gO2) const;
   std::complex<double> CpChiChaPconjUSHppPL(int gI1, int gO1) const;
   std::complex<double> CpSHp0conjUSHppVWm(int gI2, int gO2) const;
   std::complex<double> CpSHppconjUSHppVP(int gI2, int gO2) const;
   std::complex<double> CpSHppconjUSHppVZ(int gI2, int gO2) const;
   std::complex<double> CpSHppconjUSHppVZp(int gI2, int gO2) const;
   std::complex<double> CpVGVGVG() const;
   std::complex<double> CpbargGgGVG() const;
   double CpbarFdFdVGPL(int gI1, int gI2) const;
   double CpbarFdFdVGPR(int gI1, int gI2) const;
   double CpbarFDXFDXVGPL(int gI1, int gI2) const;
   double CpbarFDXFDXVGPR(int gI1, int gI2) const;
   double CpbarFuFuVGPL(int gI1, int gI2) const;
   double CpbarFuFuVGPR(int gI1, int gI2) const;
   double CpSdconjSdVGVG(int gI1, int gI2) const;
   double CpSDXconjSDXVGVG(int gI1, int gI2) const;
   double CpSuconjSuVGVG(int gI1, int gI2) const;
   double CpSdconjSdVG(int gI2, int gI1) const;
   double CpSDXconjSDXVG(int gI2, int gI1) const;
   double CpSuconjSuVG(int gI2, int gI1) const;
   std::complex<double> CpGluGluVGPL() const;
   std::complex<double> CpGluGluVGPR() const;
   double CpVGVGVGVG1() const;
   double CpVGVGVGVG2() const;
   double CpVGVGVGVG3() const;
   double CpbargWmgWmVP() const;
   double CpbargWmCgWmCVP() const;
   double CpconjVWmVPVWm() const;
   double CpbarChaPChaPVPPL() const;
   double CpbarChaPChaPVPPR() const;
   std::complex<double> CpHpmconjHpmVPVP(int gI1, int gI2) const;
   std::complex<double> CpSHppconjSHppVPVP(int gI1, int gI2) const;
   double CpHpmconjHpmVP(int gI2, int gI1) const;
   double CpSHppconjSHppVP(int gI2, int gI1) const;
   std::complex<double> CpbarChaChaVPPL(int gI1, int gI2) const;
   std::complex<double> CpbarChaChaVPPR(int gI1, int gI2) const;
   double CpbarChaIChaIVPPL(int gI1, int gI2) const;
   double CpbarChaIChaIVPPR(int gI1, int gI2) const;
   double CpbarFdFdVPPL(int gI1, int gI2) const;
   double CpbarFdFdVPPR(int gI1, int gI2) const;
   double CpbarFDXFDXVPPL(int gI1, int gI2) const;
   double CpbarFDXFDXVPPR(int gI1, int gI2) const;
   double CpbarFeFeVPPL(int gI1, int gI2) const;
   double CpbarFeFeVPPR(int gI1, int gI2) const;
   double CpbarFuFuVPPL(int gI1, int gI2) const;
   double CpbarFuFuVPPR(int gI1, int gI2) const;
   std::complex<double> CpSHIpconjSHIpVPVP(int gI1, int gI2) const;
   double CpSHIpconjSHIpVP(int gI2, int gI1) const;
   std::complex<double> CpSdconjSdVPVP(int gI1, int gI2) const;
   double CpSDXconjSDXVPVP(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeVPVP(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuVPVP(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVP(int gI2, int gI1) const;
   double CpSDXconjSDXVP(int gI2, int gI1) const;
   std::complex<double> CpSeconjSeVP(int gI2, int gI1) const;
   std::complex<double> CpSuconjSuVP(int gI2, int gI1) const;
   std::complex<double> CpHpmconjVWmVP(int gI2) const;
   double CpconjVWmVPVPVWm1() const;
   double CpconjVWmVPVPVWm2() const;
   double CpconjVWmVPVPVWm3() const;
   double CpbargWmgWmVZ() const;
   double CpbargWmCgWmCVZ() const;
   double CpconjVWmVWmVZ() const;
   double CpbarChaPChaPVZPL() const;
   double CpbarChaPChaPVZPR() const;
   std::complex<double> CpHpmconjHpmVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSHp0conjSHp0VZVZ(int gI1, int gI2) const;
   std::complex<double> CpSHppconjSHppVZVZ(int gI1, int gI2) const;
   double CpSSI0conjSSI0VZVZ(int gI1, int gI2) const;
   std::complex<double> CpHpmconjHpmVZ(int gI2, int gI1) const;
   double CpSHp0conjSHp0VZ(int gI2, int gI1) const;
   double CpSHppconjSHppVZ(int gI2, int gI1) const;
   double CpSSI0conjSSI0VZ(int gI2, int gI1) const;
   std::complex<double> CpbarChaChaVZPL(int gI1, int gI2) const;
   std::complex<double> CpbarChaChaVZPR(int gI1, int gI2) const;
   double CpbarChaIChaIVZPL(int gI1, int gI2) const;
   double CpbarChaIChaIVZPR(int gI1, int gI2) const;
   std::complex<double> CpChiPChiPVZPL(int gI1, int gI2) const;
   std::complex<double> CpChiPChiPVZPR(int gI1, int gI2) const;
   double CpFSIFSIVZPL(int gI1, int gI2) const;
   double CpFSIFSIVZPR(int gI1, int gI2) const;
   std::complex<double> CpAhAhVZVZ(int gI1, int gI2) const;
   std::complex<double> CphhhhVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSvconjSvVZVZ(int gI1, int gI2) const;
   std::complex<double> CpAhhhVZ(int gI2, int gI1) const;
   double CpSvconjSvVZ(int gI2, int gI1) const;
   double CpbarFdFdVZPL(int gI1, int gI2) const;
   double CpbarFdFdVZPR(int gI1, int gI2) const;
   double CpbarFDXFDXVZPL(int gI1, int gI2) const;
   double CpbarFDXFDXVZPR(int gI1, int gI2) const;
   double CpbarFeFeVZPL(int gI1, int gI2) const;
   double CpbarFeFeVZPR(int gI1, int gI2) const;
   double CpbarFuFuVZPL(int gI1, int gI2) const;
   double CpbarFuFuVZPR(int gI1, int gI2) const;
   double CpbarFvFvVZPL(int gI1, int gI2) const;
   double CpbarFvFvVZPR(int , int ) const;
   std::complex<double> CpSHI0conjSHI0VZVZ(int gI1, int gI2) const;
   std::complex<double> CpSHIpconjSHIpVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSHI0conjSHI0VZ(int gI2, int gI1) const;
   std::complex<double> CpSHIpconjSHIpVZ(int gI2, int gI1) const;
   std::complex<double> CpChiIChiIVZPL(int gI1, int gI2) const;
   std::complex<double> CpChiIChiIVZPR(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSDXconjSDXVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVZ(int gI2, int gI1) const;
   std::complex<double> CpSDXconjSDXVZ(int gI2, int gI1) const;
   std::complex<double> CpSeconjSeVZ(int gI2, int gI1) const;
   std::complex<double> CpSuconjSuVZ(int gI2, int gI1) const;
   std::complex<double> CpChiChiVZPL(int gI1, int gI2) const;
   std::complex<double> CpChiChiVZPR(int gI1, int gI2) const;
   std::complex<double> CpHpmconjVWmVZ(int gI2) const;
   std::complex<double> CphhVZVZ(int gI2) const;
   std::complex<double> CphhVZVZp(int gI2) const;
   double CpconjVWmVWmVZVZ1() const;
   double CpconjVWmVWmVZVZ2() const;
   double CpconjVWmVWmVZVZ3() const;
   double CpbargWmgWmVZp() const;
   double CpbargWmCgWmCVZp() const;
   double CpconjVWmVWmVZp() const;
   double CpbarChaPChaPVZpPL() const;
   double CpbarChaPChaPVZpPR() const;
   std::complex<double> CpHpmconjHpmVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSHp0conjSHp0VZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSHppconjSHppVZpVZp(int gI1, int gI2) const;
   double CpSSI0conjSSI0VZpVZp(int gI1, int gI2) const;
   std::complex<double> CpHpmconjHpmVZp(int gI2, int gI1) const;
   double CpSHp0conjSHp0VZp(int gI2, int gI1) const;
   double CpSHppconjSHppVZp(int gI2, int gI1) const;
   double CpSSI0conjSSI0VZp(int gI2, int gI1) const;
   std::complex<double> CpbarChaChaVZpPL(int gI1, int gI2) const;
   std::complex<double> CpbarChaChaVZpPR(int gI1, int gI2) const;
   double CpbarChaIChaIVZpPL(int gI1, int gI2) const;
   double CpbarChaIChaIVZpPR(int gI1, int gI2) const;
   std::complex<double> CpChiPChiPVZpPL(int gI1, int gI2) const;
   std::complex<double> CpChiPChiPVZpPR(int gI1, int gI2) const;
   double CpFSIFSIVZpPL(int gI1, int gI2) const;
   double CpFSIFSIVZpPR(int gI1, int gI2) const;
   std::complex<double> CpAhAhVZpVZp(int gI1, int gI2) const;
   std::complex<double> CphhhhVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSvconjSvVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpAhhhVZp(int gI2, int gI1) const;
   double CpSvconjSvVZp(int gI2, int gI1) const;
   double CpbarFdFdVZpPL(int gI1, int gI2) const;
   double CpbarFdFdVZpPR(int gI1, int gI2) const;
   double CpbarFDXFDXVZpPL(int gI1, int gI2) const;
   double CpbarFDXFDXVZpPR(int gI1, int gI2) const;
   double CpbarFeFeVZpPL(int gI1, int gI2) const;
   double CpbarFeFeVZpPR(int gI1, int gI2) const;
   double CpbarFuFuVZpPL(int gI1, int gI2) const;
   double CpbarFuFuVZpPR(int gI1, int gI2) const;
   double CpbarFvFvVZpPL(int gI1, int gI2) const;
   double CpbarFvFvVZpPR(int , int ) const;
   std::complex<double> CpSHI0conjSHI0VZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSHIpconjSHIpVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSHI0conjSHI0VZp(int gI2, int gI1) const;
   std::complex<double> CpSHIpconjSHIpVZp(int gI2, int gI1) const;
   std::complex<double> CpChiIChiIVZpPL(int gI1, int gI2) const;
   std::complex<double> CpChiIChiIVZpPR(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSDXconjSDXVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuVZpVZp(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVZp(int gI2, int gI1) const;
   std::complex<double> CpSDXconjSDXVZp(int gI2, int gI1) const;
   std::complex<double> CpSeconjSeVZp(int gI2, int gI1) const;
   std::complex<double> CpSuconjSuVZp(int gI2, int gI1) const;
   std::complex<double> CpChiChiVZpPL(int gI1, int gI2) const;
   std::complex<double> CpChiChiVZpPR(int gI1, int gI2) const;
   std::complex<double> CpHpmconjVWmVZp(int gI2) const;
   std::complex<double> CphhVZpVZp(int gI2) const;
   double CpconjVWmVWmVZpVZp1() const;
   double CpconjVWmVWmVZpVZp2() const;
   double CpconjVWmVWmVZpVZp3() const;
   double CpbargPgWmconjVWm() const;
   double CpbargWmCgPconjVWm() const;
   double CpbargWmCgZconjVWm() const;
   double CpbargWmCgZpconjVWm() const;
   double CpbargZgWmconjVWm() const;
   double CpbargZpgWmconjVWm() const;
   std::complex<double> CpHpmconjHpmconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSHp0conjSHp0conjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSHppconjSHppconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSHppconjSHp0conjVWm(int gI2, int gI1) const;
   std::complex<double> CpAhHpmconjVWm(int gI2, int gI1) const;
   std::complex<double> CphhHpmconjVWm(int gI2, int gI1) const;
   std::complex<double> CpChiPChaPconjVWmPL(int gI1) const;
   std::complex<double> CpChiPChaPconjVWmPR(int gI1) const;
   std::complex<double> CpAhAhconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CphhhhconjVWmVWm(int gI1, int gI2) const;
   double CpSvconjSvconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpbarFuFdconjVWmPL(int gI1, int gI2) const;
   double CpbarFuFdconjVWmPR(int , int ) const;
   std::complex<double> CpbarFvFeconjVWmPL(int gI1, int gI2) const;
   double CpbarFvFeconjVWmPR(int , int ) const;
   std::complex<double> CpSeconjSvconjVWm(int gI2, int gI1) const;
   double CpSHI0conjSHI0conjVWmVWm(int gI1, int gI2) const;
   double CpSHIpconjSHIpconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpChiIChaIconjVWmPL(int gI1, int gI2) const;
   std::complex<double> CpChiIChaIconjVWmPR(int gI1, int gI2) const;
   std::complex<double> CpSHIpconjSHI0conjVWm(int gI2, int gI1) const;
   std::complex<double> CpSdconjSdconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpChiChaconjVWmPL(int gI1, int gI2) const;
   std::complex<double> CpChiChaconjVWmPR(int gI1, int gI2) const;
   std::complex<double> CpSdconjSuconjVWm(int gI2, int gI1) const;
   std::complex<double> CphhconjVWmVWm(int gI2) const;
   double CpconjVWmconjVWmVWmVWm1() const;
   double CpconjVWmconjVWmVWmVWm2() const;
   double CpconjVWmconjVWmVWmVWm3() const;
   std::complex<double> CpbarChaUChiHpmPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChaUChiHpmPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpUChiChaconjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiChaconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiChiPconjSHp0PL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiChiPconjSHp0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiChiPSHp0PL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiChiPSHp0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiFSIconjSSI0PL(int gO2, int gI2, int gI1) const;
   double CpUChiFSIconjSSI0PR(int , int , int ) const;
   double CpUChiFSISSI0PL(int , int , int ) const;
   std::complex<double> CpUChiFSISSI0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChaIUChiSHIpPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChaIUChiSHIpPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarChaUChiVWmPL(int gI1, int gO2) const;
   std::complex<double> CpbarChaUChiVWmPR(int gI1, int gO1) const;
   std::complex<double> CpUChiChaPconjSHppPL(int gO2, int gI1) const;
   std::complex<double> CpUChiChaPconjSHppPR(int gO1, int gI1) const;
   double CpbarFvUChiSvPL(int , int , int ) const;
   std::complex<double> CpbarFvUChiSvPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpUChiFvconjSvPL(int gO2, int gI2, int gI1) const;
   double CpUChiFvconjSvPR(int , int , int ) const;
   std::complex<double> CpbarFdUChiSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarFdUChiSdPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFDXUChiSDXPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarFDXUChiSDXPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFeUChiSePL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarFeUChiSePR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFuUChiSuPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarFuUChiSuPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpChiUChihhPL(int gI2, int gO2, int gI1) const;
   std::complex<double> CpChiUChihhPR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpUChiChaIconjSHIpPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiChaIconjSHIpPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiChiIconjSHI0PL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiChiIconjSHI0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiChiISHI0PL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiChiISHI0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpChiUChiAhPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiUChiAhPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpUChiFdconjSdPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiFdconjSdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiFDXconjSDXPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiFDXconjSDXPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiFeconjSePL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiFeconjSePR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpUChiFuconjSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiFuconjSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChaPUChiSHppPL(int gO2, int gI2) const;
   std::complex<double> CpbarChaPUChiSHppPR(int gO1, int gI2) const;
   std::complex<double> CpUChiChaconjVWmPR(int gO2, int gI2) const;
   std::complex<double> CpUChiChaconjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpChiUChiVZPL(int gI2, int gO2) const;
   std::complex<double> CpChiUChiVZPR(int gI2, int gO1) const;
   std::complex<double> CpChiUChiVZpPL(int gI2, int gO2) const;
   std::complex<double> CpChiUChiVZpPR(int gI2, int gO1) const;
   std::complex<double> CpbarUChaChiPSHppPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChiPSHppPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChaAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChaChaAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChaChiHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChiHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChaPconjSHp0PL(int gO2, int gI1) const;
   std::complex<double> CpbarUChaChaPconjSHp0PR(int gO1, int gI1) const;
   std::complex<double> CpbarUChaChahhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChahhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaFeconjSvPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaFeconjSvPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChabarFuSdPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChabarFuSdPR(int gO1, int gI1, int gI2) const;
   double CpbarUChabarFvSePL(int , int , int ) const;
   std::complex<double> CpbarUChabarFvSePR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChaChaIconjSHI0PL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChaIconjSHI0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChiISHIpPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChiISHIpPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaFdconjSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaFdconjSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaChaVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaChaVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaChaVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaChaVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaChaVZpPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaChaVZpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaChiVWmPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaChiVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFvHpmPL(int gO2, int gI2, int gI1) const;
   double CpbarUFeFvHpmPR(int , int , int ) const;
   std::complex<double> CpbarUFeChaSvPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFeChaSvPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFehhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFeChiSePL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFeChiSePR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFeVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZpPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZpPL(int gO1, int gI2) const;
   double CpbarUFeFvVWmPR(int , int ) const;
   double CpbarUFeFvVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFuHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFuHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdChaSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdChaSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdChiSdPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdChiSdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdGluSdPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFdGluSdPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFdFdVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZpPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZpPL(int gO1, int gI2) const;
   double CpbarUFdFuVWmPR(int , int ) const;
   std::complex<double> CpbarUFdFuVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFdconjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFdconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChabarUFuSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChabarUFuSdPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFuChiSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuChiSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFuGluSuPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFuGluSuPR(int gO1, int gI1) const;
   double CpbarUFuFdconjVWmPR(int , int ) const;
   std::complex<double> CpbarUFuFdconjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZpPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFDXFDXAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFDXFDXAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFDXFDXhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFDXFDXhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFDXChiSDXPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFDXChiSDXPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFDXGluSDXPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFDXGluSDXPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFDXFDXVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFDXFDXVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFDXFDXVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFDXFDXVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFDXFDXVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFDXFDXVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFDXFDXVZpPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFDXFDXVZpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaIChaIAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChaIChaIAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChaIChaIhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaIChaIhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaIChaSHI0PL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaIChaSHI0PR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaIChiSHIpPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChaIChiSHIpPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChaIChaIVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaIChaIVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaIChaIVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaIChaIVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaIChaIVZpPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaIChaIVZpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChaIChiIVWmPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChaIChiIVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarChaUChiISHIpPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChaUChiISHIpPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarChaIUChiIVWmPL(int gI1, int gO2) const;
   std::complex<double> CpbarChaIUChiIVWmPR(int gI1, int gO1) const;
   std::complex<double> CpChiIUChiIhhPL(int gI2, int gO2, int gI1) const;
   std::complex<double> CpChiIUChiIhhPR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpUChiIChaconjSHIpPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiIChaconjSHIpPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpChiIUChiIAhPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpChiIUChiIAhPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpChiUChiIconjSHI0PL(int gI2, int gO2, int gI1) const;
   std::complex<double> CpChiUChiIconjSHI0PR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpChiUChiISHI0PL(int gI2, int gO2, int gI1) const;
   std::complex<double> CpChiUChiISHI0PR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpUChiIChaIconjVWmPR(int gO2, int gI2) const;
   std::complex<double> CpUChiIChaIconjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpChiIUChiIVZPL(int gI2, int gO2) const;
   std::complex<double> CpChiIUChiIVZPR(int gI2, int gO1) const;
   std::complex<double> CpChiIUChiIVZpPL(int gI2, int gO2) const;
   std::complex<double> CpChiIUChiIVZpPR(int gI2, int gO1) const;
   std::complex<double> CpChiUFSIconjSSI0PL(int gI2, int gO2, int gI1) const;
   double CpChiUFSIconjSSI0PR(int , int , int ) const;
   double CpChiUFSISSI0PL(int , int , int ) const;
   std::complex<double> CpChiUFSISSI0PR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpFSIUFSIVZPL(int gI2, int gO2) const;
   std::complex<double> CpFSIUFSIVZPR(int gI2, int gO1) const;
   std::complex<double> CpFSIUFSIVZpPL(int gI2, int gO2) const;
   std::complex<double> CpFSIUFSIVZpPR(int gI2, int gO1) const;
   std::complex<double> CpbarChaUChiPSHppPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChaUChiPSHppPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpUChiPChaconjSHppPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUChiPChaconjSHppPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpChiUChiPconjSHp0PL(int gI2, int gO2, int gI1) const;
   std::complex<double> CpChiUChiPconjSHp0PR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpChiUChiPSHp0PL(int gI2, int gO2, int gI1) const;
   std::complex<double> CpChiUChiPSHp0PR(int gI2, int gO1, int gI1) const;
   std::complex<double> CpChiPUChiPVZPL(int gI2, int gO2) const;
   std::complex<double> CpChiPUChiPVZPR(int gI2, int gO1) const;
   std::complex<double> CpChiPUChiPVZpPL(int gI2, int gO2) const;
   std::complex<double> CpChiPUChiPVZpPR(int gI2, int gO1) const;
   double CpbarChaPUChiPVWmPL(int gO2) const;
   double CpbarChaPUChiPVWmPR(int gO1) const;
   double CpUChiPChaPconjVWmPR(int gO2) const;
   double CpUChiPChaPconjVWmPL(int gO1) const;
   std::complex<double> CpbarFdGluSdPL(int gI1, int gI2) const;
   std::complex<double> CpbarFdGluSdPR(int gI1, int gI2) const;
   std::complex<double> CpbarFDXGluSDXPL(int gI1, int gI2) const;
   std::complex<double> CpbarFDXGluSDXPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuGluSuPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuGluSuPR(int gI1, int gI2) const;
   std::complex<double> CpGluFdconjSdPL(int gI2, int gI1) const;
   std::complex<double> CpGluFdconjSdPR(int gI2, int gI1) const;
   std::complex<double> CpGluFDXconjSDXPL(int gI2, int gI1) const;
   std::complex<double> CpGluFDXconjSDXPR(int gI2, int gI1) const;
   std::complex<double> CpGluFuconjSuPL(int gI2, int gI1) const;
   std::complex<double> CpGluFuconjSuPR(int gI2, int gI1) const;
   double CpbarFvFeconjHpmPL(int , int , int ) const;
   std::complex<double> CpbarFvFeconjHpmPR(int gO1, int gI2, int gI1) const;
   double CpbarChabarFvSePL(int , int , int ) const;
   std::complex<double> CpbarChabarFvSePR(int gI1, int gO1, int gI2) const;
   double CpbarFvChiSvPL(int , int , int ) const;
   std::complex<double> CpbarFvChiSvPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChaPChaSHp0PL(int gI2, int gI1) const;
   std::complex<double> CpbarChaPChaSHp0PR(int gI2, int gI1) const;
   std::complex<double> CpbarChaPChiSHppPL(int gI2, int gI1) const;
   std::complex<double> CpbarChaPChiSHppPR(int gI2, int gI1) const;
   std::complex<double> CpbarChaPChiPVWmPR(int gI2) const;
   std::complex<double> CpbarChaPChiPVWmPL(int gI2) const;
   std::complex<double> CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const;
   double CpbarFeFvHpmPR(int , int , int ) const;
   std::complex<double> CpbarFeChaSvPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeChaSvPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeFehhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFeChiSePL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeChiSePR(int gO1, int gI2, int gI1) const;
   double CpbarFeFvVWmPR(int , int ) const;
   std::complex<double> CpbarFeFvVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarFdFuHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFuHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdChaSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdChaSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdChiSdPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdChiSdPR(int gO1, int gI2, int gI1) const;
   double CpbarFdFuVWmPR(int , int ) const;
   std::complex<double> CpbarFdFuVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarFuFdconjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFdconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChabarFuSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChabarFuSdPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFuhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFuChiSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuChiSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> self_energy_Sd_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Sd_1loop(double p) const;
   std::complex<double> self_energy_Sv_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Sv_1loop(double p) const;
   std::complex<double> self_energy_Su_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Su_1loop(double p) const;
   std::complex<double> self_energy_Se_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Se_1loop(double p) const;
   std::complex<double> self_energy_SDX_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_SDX_1loop(double p) const;
   std::complex<double> self_energy_hh_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_hh_1loop(double p) const;
   std::complex<double> self_energy_Ah_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Ah_1loop(double p) const;
   std::complex<double> self_energy_Hpm_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Hpm_1loop(double p) const;
   std::complex<double> self_energy_SHI0_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_SHI0_1loop(double p) const;
   std::complex<double> self_energy_SHIp_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_SHIp_1loop(double p) const;
   std::complex<double> self_energy_SSI0_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_SSI0_1loop(double p) const;
   std::complex<double> self_energy_SHp0_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_SHp0_1loop(double p) const;
   std::complex<double> self_energy_SHpp_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_SHpp_1loop(double p) const;
   std::complex<double> self_energy_VG_1loop(double p ) const;
   std::complex<double> self_energy_VP_1loop(double p ) const;
   std::complex<double> self_energy_VZ_1loop(double p ) const;
   std::complex<double> self_energy_VZp_1loop(double p ) const;
   std::complex<double> self_energy_VWm_1loop(double p ) const;
   std::complex<double> self_energy_Chi_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Chi_1loop_1(double p) const;
   std::complex<double> self_energy_Chi_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Chi_1loop_PR(double p) const;
   std::complex<double> self_energy_Chi_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Chi_1loop_PL(double p) const;
   std::complex<double> self_energy_Cha_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha_1loop_1(double p) const;
   std::complex<double> self_energy_Cha_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha_1loop_PR(double p) const;
   std::complex<double> self_energy_Cha_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha_1loop_PL(double p) const;
   std::complex<double> self_energy_Fe_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1(double p) const;
   std::complex<double> self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR(double p) const;
   std::complex<double> self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL(double p) const;
   std::complex<double> self_energy_Fd_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1(double p) const;
   std::complex<double> self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR(double p) const;
   std::complex<double> self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL(double p) const;
   std::complex<double> self_energy_Fu_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL(double p) const;
   std::complex<double> self_energy_FDX_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_FDX_1loop_1(double p) const;
   std::complex<double> self_energy_FDX_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_FDX_1loop_PR(double p) const;
   std::complex<double> self_energy_FDX_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_FDX_1loop_PL(double p) const;
   std::complex<double> self_energy_ChaI_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_ChaI_1loop_1(double p) const;
   std::complex<double> self_energy_ChaI_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_ChaI_1loop_PR(double p) const;
   std::complex<double> self_energy_ChaI_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_ChaI_1loop_PL(double p) const;
   std::complex<double> self_energy_ChiI_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_ChiI_1loop_1(double p) const;
   std::complex<double> self_energy_ChiI_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_ChiI_1loop_PR(double p) const;
   std::complex<double> self_energy_ChiI_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_ChiI_1loop_PL(double p) const;
   std::complex<double> self_energy_FSI_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_FSI_1loop_1(double p) const;
   std::complex<double> self_energy_FSI_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_FSI_1loop_PR(double p) const;
   std::complex<double> self_energy_FSI_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_FSI_1loop_PL(double p) const;
   std::complex<double> self_energy_ChiP_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_ChiP_1loop_1(double p) const;
   std::complex<double> self_energy_ChiP_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_ChiP_1loop_PR(double p) const;
   std::complex<double> self_energy_ChiP_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_ChiP_1loop_PL(double p) const;
   std::complex<double> self_energy_Glu_1loop_1(double p ) const;
   std::complex<double> self_energy_Glu_1loop_PR(double p ) const;
   std::complex<double> self_energy_Glu_1loop_PL(double p ) const;
   std::complex<double> self_energy_Fv_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_1(double p) const;
   std::complex<double> self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PR(double p) const;
   std::complex<double> self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PL(double p) const;
   std::complex<double> self_energy_ChaP_1loop_1(double p ) const;
   std::complex<double> self_energy_ChaP_1loop_PR(double p ) const;
   std::complex<double> self_energy_ChaP_1loop_PL(double p ) const;
   std::complex<double> self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy_rotated(double p) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy(double p) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy(double p) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy(double p) const;
   std::complex<double> tadpole_hh_1loop(int gO1) const;


   /// calculates the tadpoles at current loop order
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations() const;
   /// calculates the tadpoles divided by VEVs at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations_over_vevs() const;

   void calculate_MSu_2nd_generation(double&, double&, double&) const;
   void calculate_MSd_2nd_generation(double&, double&, double&) const;
   void calculate_MSv_2nd_generation(double&, double&, double&) const;
   void calculate_MSe_2nd_generation(double&, double&, double&) const;

   void calculate_MSu_3rd_generation(double&, double&, double&) const;
   void calculate_MSd_3rd_generation(double&, double&, double&) const;
   void calculate_MSv_3rd_generation(double&, double&, double&) const;
   void calculate_MSe_3rd_generation(double&, double&, double&) const;

   Eigen::Matrix<double,3,3> self_energy_hh_2loop() const;
   Eigen::Matrix<double,3,3> self_energy_Ah_2loop() const;

   Eigen::Matrix<double,3,1> tadpole_hh_2loop() const;



   void calculate_MVG_pole();
   void calculate_MGlu_pole();
   void calculate_MFv_pole();
   void calculate_MChaP_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MVZp_pole();
   void calculate_MSd_pole();
   void calculate_MSv_pole();
   void calculate_MSu_pole();
   void calculate_MSe_pole();
   void calculate_MSDX_pole();
   void calculate_Mhh_pole();
   void calculate_MAh_pole();
   void calculate_MHpm_pole();
   void calculate_MChi_pole();
   void calculate_MCha_pole();
   void calculate_MFe_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFDX_pole();
   void calculate_MSHI0_pole();
   void calculate_MSHIp_pole();
   void calculate_MChaI_pole();
   void calculate_MChiI_pole();
   void calculate_MSSI0_pole();
   void calculate_MFSI_pole();
   void calculate_MSHp0_pole();
   void calculate_MSHpp_pole();
   void calculate_MChiP_pole();
   void calculate_MVWm_pole();
   double calculate_MVWm_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MVP_DRbar(double) const;
   double calculate_MVZ_DRbar(double) const;
   double calculate_MVWm_DRbar(double) const;

   double v() const override;
   double Betax() const override;
   double ThetaW() const override;
   double ThetaWp() const override;
   double VEV() const override;


private:
   int ewsb_loop_order{4};           ///< loop order for EWSB
   int pole_mass_loop_order{4};      ///< loop order for pole masses
   bool calculate_sm_pole_masses{false};  ///< switch to calculate the pole masses of the Standard Model particles
   bool calculate_bsm_pole_masses{true};  ///< switch to calculate the pole masses of the BSM particles
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-4};               ///< RG running precision
   double ewsb_iteration_precision{1.e-5};///< precision goal of EWSB solution
   E6SSMEFTHiggs_physical physical{}; ///< contains the pole masses and mixings
   mutable Problems problems{E6SSMEFTHiggs_info::model_name,
                             &E6SSMEFTHiggs_info::particle_names_getter,
                             &E6SSMEFTHiggs_info::parameter_names_getter}; ///< problems
   Loop_corrections loop_corrections{}; ///< used pole mass corrections
   std::shared_ptr<E6SSMEFTHiggs_ewsb_solver_interface> ewsb_solver{};
   Threshold_corrections threshold_corrections{}; ///< used threshold corrections

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   void copy_DRbar_masses_to_pole_masses();

   // Passarino-Veltman loop functions
   double A0(double) const noexcept;
   double B0(double, double, double) const noexcept;
   double B1(double, double, double) const noexcept;
   double B00(double, double, double) const noexcept;
   double B22(double, double, double) const noexcept;
   double H0(double, double, double) const noexcept;
   double F0(double, double, double) const noexcept;
   double G0(double, double, double) const noexcept;

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

std::ostream& operator<<(std::ostream&, const E6SSMEFTHiggs_mass_eigenstates&);

} // namespace flexiblesusy

#undef SUPER

#endif
