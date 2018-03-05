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

// File generated at Mon 5 Mar 2018 17:51:22

/**
 * @file MRSSM_mass_eigenstates.hpp
 *
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 5 Mar 2018 17:51:22 with FlexibleSUSY
 * 2.1.0 (git commit: 8f20f6c9c42c159c1588fbc0bb3e15ce5ab6ace3) and SARAH 4.12.3 .
 */

#ifndef MRSSM_MASS_EIGENSTATES_H
#define MRSSM_MASS_EIGENSTATES_H

#include "MRSSM_info.hpp"
#include "MRSSM_physical.hpp"
#include "MRSSM_soft_parameters.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "error.hpp"
#include "problems.hpp"
#include "config.h"

#include <iosfwd>
#include <memory>
#include <string>

#include <Eigen/Core>

namespace flexiblesusy {

class MRSSM_ewsb_solver_interface;
/**
 * @class MRSSM_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class MRSSM_mass_eigenstates : public MRSSM_soft_parameters {
public:
   explicit MRSSM_mass_eigenstates(const MRSSM_input_parameters& input_ = MRSSM_input_parameters());
   MRSSM_mass_eigenstates(const MRSSM_mass_eigenstates&) = default;
   MRSSM_mass_eigenstates(MRSSM_mass_eigenstates&&) = default;
   virtual ~MRSSM_mass_eigenstates() = default;
   MRSSM_mass_eigenstates& operator=(const MRSSM_mass_eigenstates&) = default;
   MRSSM_mass_eigenstates& operator=(MRSSM_mass_eigenstates&&) = default;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 4;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   virtual void clear() override;
   void clear_DRbar_parameters();
   Eigen::ArrayXd get_DRbar_masses() const;
   Eigen::ArrayXd get_DRbar_masses_and_mixings() const;
   Eigen::ArrayXd get_extra_parameters() const;
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
   void set_extra_parameters(const Eigen::ArrayXd&);
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   void set_physical(const MRSSM_physical&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   const MRSSM_physical& get_physical() const;
   MRSSM_physical& get_physical();
   const Problems& get_problems() const;
   Problems& get_problems();
   void set_ewsb_solver(const std::shared_ptr<MRSSM_ewsb_solver_interface>&);
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   void calculate_spectrum();
   void clear_problems();
   std::string name() const;
   void run_to(double scale, double eps = -1.0) override;
   void print(std::ostream& out = std::cerr) const override;
   void set_precision(double);
   double get_precision() const;

   double get_lsp(MRSSM_info::Particles&) const;

   double get_MVG() const { return MVG; }
   double get_MGlu() const { return MGlu; }
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   double get_MSRdp() const { return MSRdp; }
   double get_MSRum() const { return MSRum; }
   double get_MsigmaO() const { return MsigmaO; }
   double get_MphiO() const { return MphiO; }
   const Eigen::Array<double,6,1>& get_MSd() const { return MSd; }
   double get_MSd(int i) const { return MSd(i); }
   const Eigen::Array<double,3,1>& get_MSv() const { return MSv; }
   double get_MSv(int i) const { return MSv(i); }
   const Eigen::Array<double,6,1>& get_MSu() const { return MSu; }
   double get_MSu(int i) const { return MSu(i); }
   const Eigen::Array<double,6,1>& get_MSe() const { return MSe; }
   double get_MSe(int i) const { return MSe(i); }
   const Eigen::Array<double,4,1>& get_Mhh() const { return Mhh; }
   double get_Mhh(int i) const { return Mhh(i); }
   const Eigen::Array<double,4,1>& get_MAh() const { return MAh; }
   double get_MAh(int i) const { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MRh() const { return MRh; }
   double get_MRh(int i) const { return MRh(i); }
   const Eigen::Array<double,4,1>& get_MHpm() const { return MHpm; }
   double get_MHpm(int i) const { return MHpm(i); }
   const Eigen::Array<double,4,1>& get_MChi() const { return MChi; }
   double get_MChi(int i) const { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha1() const { return MCha1; }
   double get_MCha1(int i) const { return MCha1(i); }
   const Eigen::Array<double,2,1>& get_MCha2() const { return MCha2; }
   double get_MCha2(int i) const { return MCha2(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   double get_MVWm() const { return MVWm; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }

   
   Eigen::Array<double,3,1> get_MChargedHiggs() const;

   Eigen::Array<double,3,1> get_MPseudoscalarHiggs() const;

   const Eigen::Matrix<double,6,6>& get_ZD() const { return ZD; }
   double get_ZD(int i, int k) const { return ZD(i,k); }
   const Eigen::Matrix<double,3,3>& get_ZV() const { return ZV; }
   double get_ZV(int i, int k) const { return ZV(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZU() const { return ZU; }
   double get_ZU(int i, int k) const { return ZU(i,k); }
   const Eigen::Matrix<double,6,6>& get_ZE() const { return ZE; }
   double get_ZE(int i, int k) const { return ZE(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZH() const { return ZH; }
   double get_ZH(int i, int k) const { return ZH(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZA() const { return ZA; }
   double get_ZA(int i, int k) const { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZHR() const { return ZHR; }
   double get_ZHR(int i, int k) const { return ZHR(i,k); }
   const Eigen::Matrix<double,4,4>& get_ZP() const { return ZP; }
   double get_ZP(int i, int k) const { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN1() const { return ZN1; }
   std::complex<double> get_ZN1(int i, int k) const { return ZN1(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN2() const { return ZN2; }
   std::complex<double> get_ZN2(int i, int k) const { return ZN2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM1() const { return UM1; }
   std::complex<double> get_UM1(int i, int k) const { return UM1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP1() const { return UP1; }
   std::complex<double> get_UP1(int i, int k) const { return UP1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM2() const { return UM2; }
   std::complex<double> get_UM2(int i, int k) const { return UM2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP2() const { return UP2; }
   std::complex<double> get_UP2(int i, int k) const { return UP2(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL() const { return ZEL; }
   std::complex<double> get_ZEL(int i, int k) const { return ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER() const { return ZER; }
   std::complex<double> get_ZER(int i, int k) const { return ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL() const { return ZDL; }
   std::complex<double> get_ZDL(int i, int k) const { return ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR() const { return ZDR; }
   std::complex<double> get_ZDR(int i, int k) const { return ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL() const { return ZUL; }
   std::complex<double> get_ZUL(int i, int k) const { return ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR() const { return ZUR; }
   std::complex<double> get_ZUR(int i, int k) const { return ZUR(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const { return ZZ; }
   double get_ZZ(int i, int k) const { return ZZ(i,k); }




   double get_mass_matrix_VG() const;
   void calculate_MVG();
   double get_mass_matrix_Glu() const;
   void calculate_MGlu();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   double get_mass_matrix_SRdp() const;
   void calculate_MSRdp();
   double get_mass_matrix_SRum() const;
   void calculate_MSRum();
   double get_mass_matrix_sigmaO() const;
   void calculate_MsigmaO();
   double get_mass_matrix_phiO() const;
   void calculate_MphiO();
   Eigen::Matrix<double,6,6> get_mass_matrix_Sd() const;
   void calculate_MSd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Sv() const;
   void calculate_MSv();
   Eigen::Matrix<double,6,6> get_mass_matrix_Su() const;
   void calculate_MSu();
   Eigen::Matrix<double,6,6> get_mass_matrix_Se() const;
   void calculate_MSe();
   Eigen::Matrix<double,4,4> get_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<double,4,4> get_mass_matrix_Ah() const;
   void calculate_MAh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Rh() const;
   void calculate_MRh();
   Eigen::Matrix<double,4,4> get_mass_matrix_Hpm() const;
   void calculate_MHpm();
   Eigen::Matrix<double,4,4> get_mass_matrix_Chi() const;
   void calculate_MChi();
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha1() const;
   void calculate_MCha1();
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha2() const;
   void calculate_MCha2();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   double get_mass_matrix_VWm() const;
   void calculate_MVWm();
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const;
   void calculate_MVPVZ();

   double get_ewsb_eq_hh_1() const;
   double get_ewsb_eq_hh_2() const;
   double get_ewsb_eq_hh_3() const;
   double get_ewsb_eq_hh_4() const;

   std::complex<double> CpSRdpUSdconjSRdpconjUSd(int gO1, int gO2) const;
   std::complex<double> CpSRumUSdconjSRumconjUSd(int gO1, int gO2) const;
   std::complex<double> CpUSdconjUSdVZVZ(int gO1, int gO2) const;
   double CpUSdconjUSdconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpRhUSdconjRhconjUSd(int gI1, int gO1, int gI2, int gO2) const;
   double CpbarCha1FuconjUSdPR(int , int , int ) const;
   std::complex<double> CpbarCha1FuconjUSdPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUSdSvconjUSdconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpCha2FuconjUSdPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpCha2FuconjUSdPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpChiFdconjUSdPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFdconjUSdPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhAhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpHpmUSdconjHpmconjUSd(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpbarChiFdconjUSdPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChiFdconjUSdPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpphiOSdconjUSd(int gI1, int gO2) const;
   std::complex<double> CpSRumSuconjUSd(int gI1, int gO2) const;
   std::complex<double> CpUSdconjUSdconjSdSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSdconjUSdconjSuSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSdSeconjUSdconjSe(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpRhSdconjUSd(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhSdconjUSd(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSdconjUSd(int gI2, int gI1, int gO2) const;
   std::complex<double> CpHpmSuconjUSd(int gI2, int gI1, int gO2) const;
   double CpGluFdconjUSdPR(int , int ) const;
   std::complex<double> CpGluFdconjUSdPL(int gI2, int gO1) const;
   std::complex<double> CpbarGluFdconjUSdPR(int gI2, int gO2) const;
   double CpbarGluFdconjUSdPL(int , int ) const;
   std::complex<double> CpsigmaOSdconjUSd(int gI2, int gO2) const;
   std::complex<double> CpSuconjSRdpconjUSd(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSdVG(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSdVP(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSdVZ(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSdVWm(int gI2, int gO2) const;
   std::complex<double> CpSRdpUSvconjSRdpconjUSv(int gO1, int gO2) const;
   std::complex<double> CpSRumUSvconjSRumconjUSv(int gO1, int gO2) const;
   std::complex<double> CpUSvconjUSvVZVZ(int gO1, int gO2) const;
   double CpUSvconjUSvconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpRhUSvconjRhconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSvUSvconjSvconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpCha1FeconjUSvPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpCha1FeconjUSvPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhSvconjUSv(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSvconjUSv(int gI2, int gI1, int gO2) const;
   double CpChiFvconjUSvPR(int , int , int ) const;
   std::complex<double> CpChiFvconjUSvPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhAhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpHpmUSvconjHpmconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeconjHpmconjUSv(int gI2, int gI1, int gO2) const;
   std::complex<double> CpSdUSvconjSdconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSvconjSeconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSuUSvconjSuconjUSv(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSvconjUSvVZ(int gI2, int gO2) const;
   std::complex<double> CpSRdpSeconjUSv(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSvconjVWm(int gI2, int gO2) const;
   std::complex<double> CpSRdpUSuconjSRdpconjUSu(int gO1, int gO2) const;
   std::complex<double> CpSRumUSuconjSRumconjUSu(int gO1, int gO2) const;
   std::complex<double> CpUSuconjUSuVZVZ(int gO1, int gO2) const;
   double CpUSuconjUSuconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpRhUSuconjRhconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   double CpbarCha2FdconjUSuPR(int , int , int ) const;
   std::complex<double> CpbarCha2FdconjUSuPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUSuSvconjUSuconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpCha1FdconjUSuPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpCha1FdconjUSuPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpChiFuconjUSuPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFuconjUSuPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhAhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpHpmUSuconjHpmconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpbarChiFuconjUSuPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChiFuconjUSuPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSdconjHpmconjUSu(int gI2, int gI1, int gO2) const;
   std::complex<double> CpphiOSuconjUSu(int gI1, int gO2) const;
   std::complex<double> CpsigmaOSuconjUSu(int gI1, int gO2) const;
   std::complex<double> CpSeUSuconjSeconjUSu(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSuconjUSuconjSdSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUSuconjUSuconjSuSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpRhSuconjUSu(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhSuconjUSu(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSuconjUSu(int gI2, int gI1, int gO2) const;
   double CpGluFuconjUSuPR(int , int ) const;
   std::complex<double> CpGluFuconjUSuPL(int gI2, int gO1) const;
   std::complex<double> CpbarGluFuconjUSuPR(int gI2, int gO2) const;
   double CpbarGluFuconjUSuPL(int , int ) const;
   std::complex<double> CpSRdpSdconjUSu(int gI2, int gO2) const;
   std::complex<double> CpSdconjSRumconjUSu(int gI2, int gO2) const;
   std::complex<double> CpSdconjUSuconjVWm(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVG(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVP(int gI2, int gO2) const;
   std::complex<double> CpSuconjUSuVZ(int gI2, int gO2) const;
   std::complex<double> CpSRdpUSeconjSRdpconjUSe(int gO1, int gO2) const;
   std::complex<double> CpSRumUSeconjSRumconjUSe(int gO1, int gO2) const;
   std::complex<double> CpUSeconjUSeVZVZ(int gO1, int gO2) const;
   double CpUSeconjUSeconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpRhUSeconjRhconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   double CpbarCha1FvconjUSePR(int , int , int ) const;
   std::complex<double> CpbarCha1FvconjUSePL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUSeSvconjUSeconjSv(int gO1, int gI1, int gO2, int gI2) const;
   double CpCha2FvconjUSePR(int , int , int ) const;
   std::complex<double> CpCha2FvconjUSePL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpHpmSvconjUSe(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFeconjUSePR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiFeconjUSePL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhAhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpHpmUSeconjHpmconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpbarChiFeconjUSePR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChiFeconjUSePL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpSdUSeconjSdconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpSeUSeconjSeconjUSe(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpUSeSuconjUSeconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpRhSeconjUSe(int gI2, int gI1, int gO2) const;
   std::complex<double> CpAhSeconjUSe(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhSeconjUSe(int gI2, int gI1, int gO2) const;
   std::complex<double> CpSvconjSRdpconjUSe(int gI2, int gO2) const;
   std::complex<double> CpSvconjUSeVWm(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSeVP(int gI2, int gO2) const;
   std::complex<double> CpSeconjUSeVZ(int gI2, int gO2) const;
   std::complex<double> CpSRdpUhhconjSRdp(int gO2) const;
   std::complex<double> CpSRumUhhconjSRum(int gO2) const;
   std::complex<double> CpbargWmgWmUhh(int gO1) const;
   std::complex<double> CpbargWmCgWmCUhh(int gO1) const;
   std::complex<double> CpbargZgZUhh(int gO1) const;
   std::complex<double> CpUhhVZVZ(int gO2) const;
   std::complex<double> CpUhhconjVWmVWm(int gO2) const;
   std::complex<double> CpSRdpUhhUhhconjSRdp(int gO1, int gO2) const;
   std::complex<double> CpSRumUhhUhhconjSRum(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUhhUhhRhconjRh(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhRhconjRh(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarCha1Cha1UhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarCha1Cha1UhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarCha2Cha2UhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarCha2Cha2UhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhUhhconjRh(int gI2, int gO2, int gI1) const;
   std::complex<double> CphhUhhconjRh(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUhhUhhSvconjSv(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhSvconjSv(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUhhUhhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpAhAhUhh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhhhUhh(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhhhUhh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpUhhHpmconjHpm(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarChiChiUhhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChiChiUhhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUhhUhhSdconjSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSeconjSe(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhUhhSuconjSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUhhSdconjSd(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSeconjSe(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUhhSuconjSu(int gO2, int gI2, int gI1) const;
   std::complex<double> CpSRdpUhhHpm(int gO2, int gI2) const;
   std::complex<double> CpUhhHpmconjSRum(int gO2, int gI2) const;
   std::complex<double> CpAhUhhVZ(int gI2, int gO2) const;
   std::complex<double> CpUhhHpmconjVWm(int gO2, int gI2) const;
   std::complex<double> CpSRdpUAhconjSRdp(int gO2) const;
   std::complex<double> CpSRumUAhconjSRum(int gO2) const;
   std::complex<double> CpbargWmgWmUAh(int gO1) const;
   std::complex<double> CpbargWmCgWmCUAh(int gO1) const;
   std::complex<double> CpSRdpUAhUAhconjSRdp(int gO1, int gO2) const;
   std::complex<double> CpSRumUAhUAhconjSRum(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpUAhUAhRhconjRh(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhRhconjRh(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarCha1Cha1UAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarCha1Cha1UAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarCha2Cha2UAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarCha2Cha2UAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhUAhconjRh(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUAhhhconjRh(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhUAhSvconjSv(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhSvconjSv(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpAhAhUAh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhUAhhh(int gI2, int gO2, int gI1) const;
   std::complex<double> CpUAhhhhh(int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhHpmconjHpm(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarChiChiUAhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChiChiUAhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUAhUAhSdconjSd(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSeconjSe(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhUAhSuconjSu(int gO1, int gO2, int gI1, int gI2) const;
   std::complex<double> CpUAhSdconjSd(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhSeconjSe(int gO2, int gI2, int gI1) const;
   std::complex<double> CpUAhSuconjSu(int gO2, int gI2, int gI1) const;
   std::complex<double> CpSRdpUAhHpm(int gO2, int gI2) const;
   std::complex<double> CpUAhHpmconjSRum(int gO2, int gI2) const;
   std::complex<double> CpUAhhhVZ(int gO2, int gI2) const;
   std::complex<double> CpUAhHpmconjVWm(int gO2, int gI2) const;
   std::complex<double> CpSRdpURhconjSRdpconjURh(int gO1, int gO2) const;
   std::complex<double> CpSRumURhconjSRumconjURh(int gO1, int gO2) const;
   std::complex<double> CpURhconjURhVZVZ(int gO1, int gO2) const;
   std::complex<double> CpURhconjURhconjVWmVWm(int gO1, int gO2) const;
   double CpSRdpconjURhVWm(int gO2) const;
   double CpSRumconjURhconjVWm(int gO2) const;
   std::complex<double> CpRhURhconjRhconjURh(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpCha1Cha2conjURhPR(int gI2, int gI1, int gO2) const;
   std::complex<double> CpCha1Cha2conjURhPL(int gI2, int gI1, int gO1) const;
   std::complex<double> CpAhRhconjURh(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhRhconjURh(int gI2, int gI1, int gO2) const;
   std::complex<double> CpURhSvconjURhconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSRumconjHpmconjURh(int gI1, int gO2) const;
   std::complex<double> CpAhAhURhconjURh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhURhconjURh(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpHpmURhconjHpmconjURh(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpAhAhconjURh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpAhhhconjURh(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhhhconjURh(int gI1, int gI2, int gO2) const;
   std::complex<double> CpHpmconjHpmconjURh(int gI2, int gI1, int gO2) const;
   std::complex<double> CpChiChiconjURhPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpChiChiconjURhPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpURhSdconjURhconjSd(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpURhSeconjURhconjSe(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpURhSuconjURhconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSdconjURhconjSd(int gI2, int gO2, int gI1) const;
   std::complex<double> CpSeconjURhconjSe(int gI2, int gO2, int gI1) const;
   std::complex<double> CpSuconjURhconjSu(int gI2, int gO2, int gI1) const;
   std::complex<double> CpRhconjURhVZ(int gI2, int gO2) const;
   std::complex<double> CpSRdpHpmconjURh(int gI2, int gO2) const;
   std::complex<double> CpbargWmgZUHpm(int gO2) const;
   std::complex<double> CpbargZgWmconjUHpm(int gO1) const;
   std::complex<double> CpbargWmCgZconjUHpm(int gO1) const;
   std::complex<double> CpbargZgWmCUHpm(int gO2) const;
   std::complex<double> CpconjUHpmVPVWm(int gO2) const;
   std::complex<double> CpconjUHpmVWmVZ(int gO2) const;
   std::complex<double> CpSRdpUHpmconjSRdpconjUHpm(int gO1, int gO2) const;
   std::complex<double> CpSRumUHpmconjSRumconjUHpm(int gO1, int gO2) const;
   std::complex<double> CpUHpmconjUHpmVZVZ(int gO1, int gO2) const;
   std::complex<double> CpUHpmconjUHpmconjVWmVWm(int gO1, int gO2) const;
   std::complex<double> CpSRumconjUHpmconjRh(int gO2, int gI1) const;
   std::complex<double> CpUHpmRhconjUHpmconjRh(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpHpmconjUHpmconjRh(int gI2, int gO2, int gI1) const;
   std::complex<double> CpbarCha1ChiconjUHpmPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarCha1ChiconjUHpmPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpUHpmSvconjUHpmconjSv(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarFuFdconjUHpmPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarFuFdconjUHpmPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpbarFvFeconjUHpmPR(int gI1, int gI2, int gO2) const;
   double CpbarFvFeconjUHpmPL(int , int , int ) const;
   std::complex<double> CpSeconjUHpmconjSv(int gI2, int gO2, int gI1) const;
   std::complex<double> CpAhAhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CphhhhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const;
   std::complex<double> CpHpmUHpmconjHpmconjUHpm(int gI1, int gO1, int gI2, int gO2) const;
   std::complex<double> CpbarChiCha2conjUHpmPR(int gI1, int gI2, int gO2) const;
   std::complex<double> CpbarChiCha2conjUHpmPL(int gI1, int gI2, int gO1) const;
   std::complex<double> CpAhHpmconjUHpm(int gI2, int gI1, int gO2) const;
   std::complex<double> CphhHpmconjUHpm(int gI2, int gI1, int gO2) const;
   std::complex<double> CpUHpmSdconjUHpmconjSd(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSeconjUHpmconjSe(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpUHpmSuconjUHpmconjSu(int gO1, int gI1, int gO2, int gI2) const;
   std::complex<double> CpSdconjUHpmconjSu(int gI2, int gO2, int gI1) const;
   std::complex<double> CpRhconjSRdpconjUHpm(int gI2, int gO2) const;
   std::complex<double> CpSRumAhconjUHpm(int gI2, int gO2) const;
   std::complex<double> CpSRumhhconjUHpm(int gI2, int gO2) const;
   std::complex<double> CpAhconjSRdpconjUHpm(int gI2, int gO2) const;
   std::complex<double> CphhconjSRdpconjUHpm(int gI2, int gO2) const;
   std::complex<double> CpAhconjUHpmVWm(int gI2, int gO2) const;
   std::complex<double> CphhconjUHpmVWm(int gI2, int gO2) const;
   std::complex<double> CpHpmconjUHpmVP(int gI2, int gO2) const;
   std::complex<double> CpHpmconjUHpmVZ(int gI2, int gO2) const;
   double CpSRdpSRdpconjSRdpconjSRdp() const;
   double CpSRdpSRumconjSRdpconjSRum() const;
   std::complex<double> CpSRdpconjSRdpVZVZ() const;
   double CpSRdpconjSRdpconjVWmVWm() const;
   double CpSRdpconjSRdpVP() const;
   double CpSRdpconjSRdpVZ() const;
   std::complex<double> CpSRdpRhconjSRdpconjRh(int gI1, int gI2) const;
   double CpSRdpSvconjSRdpconjSv(int gI1, int gI2) const;
   std::complex<double> CpSRdpAhAhconjSRdp(int gI1, int gI2) const;
   std::complex<double> CpSRdphhhhconjSRdp(int gI1, int gI2) const;
   std::complex<double> CpSRdpHpmconjSRdpconjHpm(int gI1, int gI2) const;
   std::complex<double> CpRhconjSRdpconjHpm(int gI2, int gI1) const;
   std::complex<double> CpCha1ChiconjSRdpPR(int gI2, int gI1) const;
   std::complex<double> CpCha1ChiconjSRdpPL(int gI2, int gI1) const;
   std::complex<double> CpAhconjSRdpconjHpm(int gI2, int gI1) const;
   std::complex<double> CphhconjSRdpconjHpm(int gI2, int gI1) const;
   std::complex<double> CpSRdpSdconjSRdpconjSd(int gI1, int gI2) const;
   std::complex<double> CpSRdpSeconjSRdpconjSe(int gI1, int gI2) const;
   std::complex<double> CpSRdpSuconjSRdpconjSu(int gI1, int gI2) const;
   std::complex<double> CpSvconjSRdpconjSe(int gI2, int gI1) const;
   std::complex<double> CpSuconjSRdpconjSd(int gI2, int gI1) const;
   std::complex<double> CpRhconjSRdpconjVWm(int gI2) const;
   std::complex<double> CpSRdpAhconjSRdp(int gI2) const;
   std::complex<double> CpSRdphhconjSRdp(int gI2) const;
   double CpSRumSRumconjSRumconjSRum() const;
   std::complex<double> CpSRumconjSRumVZVZ() const;
   double CpSRumconjSRumconjVWmVWm() const;
   double CpSRumconjSRumVP() const;
   double CpSRumconjSRumVZ() const;
   std::complex<double> CpSRumRhconjSRumconjRh(int gI1, int gI2) const;
   std::complex<double> CpHpmRhconjSRum(int gI2, int gI1) const;
   double CpSRumSvconjSRumconjSv(int gI1, int gI2) const;
   std::complex<double> CpSRumAhAhconjSRum(int gI1, int gI2) const;
   std::complex<double> CpSRumhhhhconjSRum(int gI1, int gI2) const;
   std::complex<double> CpSRumHpmconjSRumconjHpm(int gI1, int gI2) const;
   std::complex<double> CpCha2ChiconjSRumPR(int gI2, int gI1) const;
   std::complex<double> CpCha2ChiconjSRumPL(int gI2, int gI1) const;
   std::complex<double> CpAhHpmconjSRum(int gI2, int gI1) const;
   std::complex<double> CphhHpmconjSRum(int gI2, int gI1) const;
   std::complex<double> CpSRumSdconjSRumconjSd(int gI1, int gI2) const;
   std::complex<double> CpSRumSeconjSRumconjSe(int gI1, int gI2) const;
   std::complex<double> CpSRumSuconjSRumconjSu(int gI1, int gI2) const;
   std::complex<double> CpSdconjSRumconjSu(int gI2, int gI1) const;
   std::complex<double> CpRhconjSRumVWm(int gI2) const;
   std::complex<double> CpSRumAhconjSRum(int gI2) const;
   std::complex<double> CpSRumhhconjSRum(int gI2) const;
   double CpsigmaOsigmaOphiOphiO() const;
   std::complex<double> CpsigmaOsigmaOVG() const;
   std::complex<double> CpsigmaOSdconjSd(int gI2, int gI1) const;
   std::complex<double> CpsigmaOSuconjSu(int gI2, int gI1) const;
   double CpbarGluGlusigmaOPR() const;
   double CpbarGluGlusigmaOPL() const;
   double CpphiOphiOsigmaOsigmaO() const;
   std::complex<double> CpphiOphiOVG() const;
   std::complex<double> CpphiOSdconjSd(int gI2, int gI1) const;
   std::complex<double> CpphiOSuconjSu(int gI2, int gI1) const;
   std::complex<double> CpbarGluGluphiOPR() const;
   std::complex<double> CpbarGluGluphiOPL() const;
   std::complex<double> CpVGVGVG() const;
   std::complex<double> CpbargGgGVG() const;
   double CpphiOphiOVGVG() const;
   double CpsigmaOsigmaOVGVG() const;
   double CpbarFdFdVGPL(int gI1, int gI2) const;
   double CpbarFdFdVGPR(int gI1, int gI2) const;
   double CpbarFuFuVGPL(int gI1, int gI2) const;
   double CpbarFuFuVGPR(int gI1, int gI2) const;
   double CpSdconjSdVGVG(int gI1, int gI2) const;
   double CpSuconjSuVGVG(int gI1, int gI2) const;
   double CpSdconjSdVG(int gI2, int gI1) const;
   double CpSuconjSuVG(int gI2, int gI1) const;
   std::complex<double> CpbarGluGluVGPL() const;
   std::complex<double> CpbarGluGluVGPR() const;
   double CpVGVGVGVG1() const;
   double CpVGVGVGVG2() const;
   double CpVGVGVGVG3() const;
   double CpbargWmgWmVP() const;
   double CpbargWmCgWmCVP() const;
   std::complex<double> CpSRdpconjSRdpVPVP() const;
   std::complex<double> CpSRumconjSRumVPVP() const;
   double CpconjVWmVPVWm() const;
   std::complex<double> CpbarCha1Cha1VPPL(int gI1, int gI2) const;
   std::complex<double> CpbarCha1Cha1VPPR(int gI1, int gI2) const;
   std::complex<double> CpbarCha2Cha2VPPL(int gI1, int gI2) const;
   std::complex<double> CpbarCha2Cha2VPPR(int gI1, int gI2) const;
   double CpbarFdFdVPPL(int gI1, int gI2) const;
   double CpbarFdFdVPPR(int gI1, int gI2) const;
   double CpbarFeFeVPPL(int gI1, int gI2) const;
   double CpbarFeFeVPPR(int gI1, int gI2) const;
   double CpbarFuFuVPPL(int gI1, int gI2) const;
   double CpbarFuFuVPPR(int gI1, int gI2) const;
   std::complex<double> CpHpmconjHpmVPVP(int gI1, int gI2) const;
   std::complex<double> CpHpmconjHpmVP(int gI2, int gI1) const;
   std::complex<double> CpSdconjSdVPVP(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeVPVP(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuVPVP(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVP(int gI2, int gI1) const;
   std::complex<double> CpSeconjSeVP(int gI2, int gI1) const;
   std::complex<double> CpSuconjSuVP(int gI2, int gI1) const;
   std::complex<double> CpHpmconjVWmVP(int gI2) const;
   double CpconjVWmVPVPVWm1() const;
   double CpconjVWmVPVPVWm2() const;
   double CpconjVWmVPVPVWm3() const;
   double CpbargWmgWmVZ() const;
   double CpbargWmCgWmCVZ() const;
   double CpconjVWmVWmVZ() const;
   std::complex<double> CpRhconjRhVZVZ(int gI1, int gI2) const;
   std::complex<double> CpRhconjRhVZ(int gI2, int gI1) const;
   std::complex<double> CpbarCha1Cha1VZPL(int gI1, int gI2) const;
   std::complex<double> CpbarCha1Cha1VZPR(int gI1, int gI2) const;
   std::complex<double> CpbarCha2Cha2VZPL(int gI1, int gI2) const;
   std::complex<double> CpbarCha2Cha2VZPR(int gI1, int gI2) const;
   std::complex<double> CpSvconjSvVZVZ(int gI1, int gI2) const;
   double CpSvconjSvVZ(int gI2, int gI1) const;
   double CpbarFdFdVZPL(int gI1, int gI2) const;
   double CpbarFdFdVZPR(int gI1, int gI2) const;
   double CpbarFeFeVZPL(int gI1, int gI2) const;
   double CpbarFeFeVZPR(int gI1, int gI2) const;
   double CpbarFuFuVZPL(int gI1, int gI2) const;
   double CpbarFuFuVZPR(int gI1, int gI2) const;
   double CpbarFvFvVZPL(int gI1, int gI2) const;
   double CpbarFvFvVZPR(int , int ) const;
   std::complex<double> CpAhAhVZVZ(int gI1, int gI2) const;
   std::complex<double> CphhhhVZVZ(int gI1, int gI2) const;
   std::complex<double> CpHpmconjHpmVZVZ(int gI1, int gI2) const;
   std::complex<double> CpAhhhVZ(int gI2, int gI1) const;
   std::complex<double> CpHpmconjHpmVZ(int gI2, int gI1) const;
   std::complex<double> CpbarChiChiVZPL(int gI1, int gI2) const;
   std::complex<double> CpbarChiChiVZPR(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuVZVZ(int gI1, int gI2) const;
   std::complex<double> CpSdconjSdVZ(int gI2, int gI1) const;
   std::complex<double> CpSeconjSeVZ(int gI2, int gI1) const;
   std::complex<double> CpSuconjSuVZ(int gI2, int gI1) const;
   std::complex<double> CphhVZVZ(int gI2) const;
   std::complex<double> CpHpmconjVWmVZ(int gI2) const;
   double CpconjVWmVWmVZVZ1() const;
   double CpconjVWmVWmVZVZ2() const;
   double CpconjVWmVWmVZVZ3() const;
   double CpbargPgWmconjVWm() const;
   double CpbargWmCgPconjVWm() const;
   double CpbargWmCgZconjVWm() const;
   double CpbargZgWmconjVWm() const;
   std::complex<double> CpSRumconjRhconjVWm(int gI1) const;
   std::complex<double> CpRhconjRhconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpbarCha1ChiconjVWmPL(int gI1, int gI2) const;
   std::complex<double> CpbarCha1ChiconjVWmPR(int gI1, int gI2) const;
   double CpSvconjSvconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpbarFuFdconjVWmPL(int gI1, int gI2) const;
   double CpbarFuFdconjVWmPR(int , int ) const;
   std::complex<double> CpbarFvFeconjVWmPL(int gI1, int gI2) const;
   double CpbarFvFeconjVWmPR(int , int ) const;
   std::complex<double> CpSeconjSvconjVWm(int gI2, int gI1) const;
   std::complex<double> CpAhAhconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CphhhhconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpHpmconjHpmconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpbarChiCha2conjVWmPL(int gI1, int gI2) const;
   std::complex<double> CpbarChiCha2conjVWmPR(int gI1, int gI2) const;
   std::complex<double> CpAhHpmconjVWm(int gI2, int gI1) const;
   std::complex<double> CphhHpmconjVWm(int gI2, int gI1) const;
   std::complex<double> CpSdconjSdconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSeconjSeconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSuconjSuconjVWmVWm(int gI1, int gI2) const;
   std::complex<double> CpSdconjSuconjVWm(int gI2, int gI1) const;
   std::complex<double> CphhconjVWmVWm(int gI2) const;
   double CpconjVWmconjVWmVWmVWm1() const;
   double CpconjVWmconjVWmVWmVWm2() const;
   double CpconjVWmconjVWmVWmVWm3() const;
   std::complex<double> CpbarCha1barUChiSRdpPL(int gI1, int gO2) const;
   std::complex<double> CpbarCha1barUChiSRdpPR(int gI1, int gO1) const;
   std::complex<double> CpbarCha2barUChiSRumPL(int gI1, int gO2) const;
   std::complex<double> CpbarCha2barUChiSRumPR(int gI1, int gO1) const;
   double CpbarUChibarFvSvPL(int , int , int ) const;
   std::complex<double> CpbarUChibarFvSvPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChibarFdSdPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChibarFdSdPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChibarFeSePL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChibarFeSePR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChibarFuSuPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChibarFuSuPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarChibarUChiRhPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarUChiRhPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUChiCha1HpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChiCha1HpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChiCha2conjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChiCha2conjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChiChiAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUChiChiAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUChiChihhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChiChihhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChiFdconjSdPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChiFdconjSdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChiFeconjSePL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChiFeconjSePR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChiFuconjSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUChiFuconjSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUChiCha1VWmPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChiCha1VWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChiCha2conjVWmPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChiCha2conjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUChiChiVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUChiChiVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUCha1barCha2RhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1barCha2RhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1Cha1AhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1Cha1AhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1barFeSvPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1barFeSvPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1barFdSuPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1barFdSuPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUCha1Cha1hhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUCha1Cha1hhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUCha1ChiconjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUCha1ChiconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUCha1barChiSRdpPL(int gO2, int gI1) const;
   std::complex<double> CpbarUCha1barChiSRdpPR(int gO1, int gI1) const;
   std::complex<double> CpbarUCha1FuconjSdPL(int gO2, int gI2, int gI1) const;
   double CpbarUCha1FuconjSdPR(int , int , int ) const;
   std::complex<double> CpbarUCha1FvconjSePL(int gO2, int gI2, int gI1) const;
   double CpbarUCha1FvconjSePR(int , int , int ) const;
   std::complex<double> CpbarUCha1Cha1VPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUCha1Cha1VPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUCha1Cha1VZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUCha1Cha1VZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUCha1ChiconjVWmPR(int gO2, int gI2) const;
   std::complex<double> CpbarUCha1ChiconjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarCha1barUCha2RhPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha1barUCha2RhPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUCha2Cha2AhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUCha2Cha2AhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUCha2barFuSdPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUCha2barFuSdPR(int gO1, int gI1, int gI2) const;
   double CpbarUCha2barFvSePL(int , int , int ) const;
   std::complex<double> CpbarUCha2barFvSePR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUCha2Cha2hhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUCha2Cha2hhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUCha2ChiHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUCha2ChiHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUCha2barChiSRumPL(int gO2, int gI1) const;
   std::complex<double> CpbarUCha2barChiSRumPR(int gO1, int gI1) const;
   std::complex<double> CpbarUCha2FdconjSuPL(int gO2, int gI2, int gI1) const;
   double CpbarUCha2FdconjSuPR(int , int , int ) const;
   std::complex<double> CpbarUCha2Cha2VPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUCha2Cha2VPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUCha2Cha2VZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUCha2Cha2VZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUCha2ChiVWmPR(int gO2, int gI2) const;
   std::complex<double> CpbarUCha2ChiVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarCha1barUFeSvPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha1barUFeSvPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFehhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFvHpmPL(int gO2, int gI2, int gI1) const;
   double CpbarUFeFvHpmPR(int , int , int ) const;
   std::complex<double> CpbarChibarUFeSePL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarUFeSePR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUFeChiSePL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFeChiSePR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFeFeVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZPL(int gO1, int gI2) const;
   double CpbarUFeFvVWmPR(int , int ) const;
   double CpbarUFeFvVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarCha1barUFdSuPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha1barUFdSuPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFuHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdFuHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChibarUFdSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarUFdSdPR(int gI1, int gO1, int gI2) const;
   double CpbarUFdCha2SuPL(int , int , int ) const;
   std::complex<double> CpbarUFdCha2SuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdChiSdPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFdChiSdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFdGluSdPL(int gO2, int gI1) const;
   double CpbarUFdGluSdPR(int , int ) const;
   std::complex<double> CpbarUFdFdVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZPL(int gO1, int gI2) const;
   double CpbarUFdFuVWmPR(int , int ) const;
   std::complex<double> CpbarUFdFuVWmPL(int gO1, int gI2) const;
   double CpbarGlubarUFdSdPL(int , int ) const;
   std::complex<double> CpbarGlubarUFdSdPR(int gO1, int gI2) const;
   std::complex<double> CpbarCha2barUFuSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha2barUFuSdPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarUFuFdconjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFdconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChibarUFuSuPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarUFuSuPR(int gI1, int gO1, int gI2) const;
   double CpbarUFuCha1SdPL(int , int , int ) const;
   std::complex<double> CpbarUFuCha1SdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFuChiSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarUFuChiSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarUFuGluSuPL(int gO2, int gI1) const;
   double CpbarUFuGluSuPR(int , int ) const;
   double CpbarUFuFdconjVWmPR(int , int ) const;
   std::complex<double> CpbarUFuFdconjVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZPL(int gO1, int gI2) const;
   double CpbarGlubarUFuSuPL(int , int ) const;
   std::complex<double> CpbarGlubarUFuSuPR(int gO1, int gI2) const;
   double CpbarGlubarFdSdPL(int , int ) const;
   std::complex<double> CpbarGlubarFdSdPR(int gI1, int gI2) const;
   double CpbarGlubarFuSuPL(int , int ) const;
   std::complex<double> CpbarGlubarFuSuPR(int gI1, int gI2) const;
   double CpbarGluFdconjSdPL(int , int ) const;
   std::complex<double> CpbarGluFdconjSdPR(int gI2, int gI1) const;
   double CpbarGluFuconjSuPL(int , int ) const;
   std::complex<double> CpbarGluFuconjSuPR(int gI2, int gI1) const;
   double CpbarCha2barFvSePL(int , int , int ) const;
   std::complex<double> CpbarCha2barFvSePR(int gI1, int gO1, int gI2) const;
   double CpbarChibarFvSvPL(int , int , int ) const;
   std::complex<double> CpbarChibarFvSvPR(int gI1, int gO1, int gI2) const;
   double CpbarFvFeconjHpmPL(int , int , int ) const;
   std::complex<double> CpbarFvFeconjHpmPR(int gO1, int gI2, int gI1) const;
   double CpbarFvCha1SePL(int , int , int ) const;
   std::complex<double> CpbarFvCha1SePR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarCha1barFeSvPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha1barFeSvPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFeFeAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeFehhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const;
   double CpbarFeFvHpmPR(int , int , int ) const;
   std::complex<double> CpbarChibarFeSePL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarFeSePR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFeChiSePL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFeChiSePR(int gO1, int gI2, int gI1) const;
   double CpbarFeFvVWmPR(int , int ) const;
   std::complex<double> CpbarFeFvVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarCha1barFdSuPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha1barFdSuPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFdFdAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFdhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdFuHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdFuHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChibarFdSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarFdSdPR(int gI1, int gO1, int gI2) const;
   double CpbarFdCha2SuPL(int , int , int ) const;
   std::complex<double> CpbarFdCha2SuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdChiSdPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFdChiSdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFdGluSdPL(int gO2, int gI1) const;
   double CpbarFdGluSdPR(int , int ) const;
   double CpbarFdFuVWmPR(int , int ) const;
   std::complex<double> CpbarFdFuVWmPL(int gO1, int gI2) const;
   std::complex<double> CpbarCha2barFuSdPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarCha2barFuSdPR(int gI1, int gO1, int gI2) const;
   std::complex<double> CpbarFuFuAhPL(int gO2, int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPR(int gO1, int gI1, int gI2) const;
   std::complex<double> CpbarFuFdconjHpmPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFdconjHpmPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFuFuhhPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuFuhhPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarChibarFuSuPL(int gI1, int gO2, int gI2) const;
   std::complex<double> CpbarChibarFuSuPR(int gI1, int gO1, int gI2) const;
   double CpbarFuCha1SdPL(int , int , int ) const;
   std::complex<double> CpbarFuCha1SdPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFuChiSuPL(int gO2, int gI2, int gI1) const;
   std::complex<double> CpbarFuChiSuPR(int gO1, int gI2, int gI1) const;
   std::complex<double> CpbarFuGluSuPL(int gO2, int gI1) const;
   double CpbarFuGluSuPR(int , int ) const;
   std::complex<double> self_energy_Sd_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Sd_1loop(double p) const;
   std::complex<double> self_energy_Sv_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Sv_1loop(double p) const;
   std::complex<double> self_energy_Su_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Su_1loop(double p) const;
   std::complex<double> self_energy_Se_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,6,6> self_energy_Se_1loop(double p) const;
   std::complex<double> self_energy_hh_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_hh_1loop(double p) const;
   std::complex<double> self_energy_Ah_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_Ah_1loop(double p) const;
   std::complex<double> self_energy_Rh_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Rh_1loop(double p) const;
   std::complex<double> self_energy_Hpm_1loop(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_Hpm_1loop(double p) const;
   std::complex<double> self_energy_SRdp_1loop(double p ) const;
   std::complex<double> self_energy_SRum_1loop(double p ) const;
   std::complex<double> self_energy_sigmaO_1loop(double p ) const;
   std::complex<double> self_energy_phiO_1loop(double p ) const;
   std::complex<double> self_energy_VG_1loop(double p ) const;
   std::complex<double> self_energy_VP_1loop(double p ) const;
   std::complex<double> self_energy_VZ_1loop(double p ) const;
   std::complex<double> self_energy_VWm_1loop(double p ) const;
   std::complex<double> self_energy_Chi_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_Chi_1loop_1(double p) const;
   std::complex<double> self_energy_Chi_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_Chi_1loop_PR(double p) const;
   std::complex<double> self_energy_Chi_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,4,4> self_energy_Chi_1loop_PL(double p) const;
   std::complex<double> self_energy_Cha1_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha1_1loop_1(double p) const;
   std::complex<double> self_energy_Cha1_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha1_1loop_PR(double p) const;
   std::complex<double> self_energy_Cha1_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha1_1loop_PL(double p) const;
   std::complex<double> self_energy_Cha2_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha2_1loop_1(double p) const;
   std::complex<double> self_energy_Cha2_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha2_1loop_PR(double p) const;
   std::complex<double> self_energy_Cha2_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,2,2> self_energy_Cha2_1loop_PL(double p) const;
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
   std::complex<double> self_energy_Glu_1loop_1(double p ) const;
   std::complex<double> self_energy_Glu_1loop_PR(double p ) const;
   std::complex<double> self_energy_Glu_1loop_PL(double p ) const;
   std::complex<double> self_energy_Fv_1loop_1(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_1(double p) const;
   std::complex<double> self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PR(double p) const;
   std::complex<double> self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PL(double p) const;
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
   std::complex<double> tadpole_phiO_1loop() const;


   /// calculates the tadpoles at current loop order
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations() const;
   /// calculates the tadpoles divided by VEVs at current loop order
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole_equations_over_vevs() const;







   void calculate_MVG_pole();
   void calculate_MGlu_pole();
   void calculate_MFv_pole();
   void calculate_MSRdp_pole();
   void calculate_MSRum_pole();
   void calculate_MsigmaO_pole();
   void calculate_MphiO_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MSd_pole();
   void calculate_MSv_pole();
   void calculate_MSu_pole();
   void calculate_MSe_pole();
   void calculate_Mhh_pole();
   void calculate_MAh_pole();
   void calculate_MRh_pole();
   void calculate_MHpm_pole();
   void calculate_MChi_pole();
   void calculate_MCha1_pole();
   void calculate_MCha2_pole();
   void calculate_MFe_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MVWm_pole();
   double calculate_MVWm_pole(double);
   double calculate_MVZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_MVP_DRbar(double);
   double calculate_MVZ_DRbar(double);
   double calculate_MVWm_DRbar(double);

   double v() const;
   double Betax() const;
   double ThetaW() const;


private:
   int ewsb_loop_order{2};           ///< loop order for EWSB
   int pole_mass_loop_order{2};      ///< loop order for pole masses
   bool calculate_sm_pole_masses{false};  ///< switch to calculate the pole masses of the Standard Model particles
   bool calculate_bsm_pole_masses{true};  ///< switch to calculate the pole masses of the BSM particles
   bool force_output{false};              ///< switch to force output of pole masses
   double precision{1.e-3};               ///< RG running precision
   double ewsb_iteration_precision{1.e-5};///< precision goal of EWSB solution
   MRSSM_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{MRSSM_info::model_name,
                     &MRSSM_info::particle_names_getter,
                     &MRSSM_info::parameter_names_getter}; ///< problems
   Loop_corrections loop_corrections{}; ///< used pole mass corrections
   std::shared_ptr<MRSSM_ewsb_solver_interface> ewsb_solver{};
   Threshold_corrections threshold_corrections{}; ///< used threshold corrections

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   int solve_ewsb_tree_level_custom();
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

std::ostream& operator<<(std::ostream&, const MRSSM_mass_eigenstates&);

} // namespace flexiblesusy

#endif
