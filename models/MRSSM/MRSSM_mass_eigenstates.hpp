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

// File generated at Mon 19 Sep 2016 09:58:14

/**
 * @file MRSSM_mass_eigenstates.hpp
 *
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Mon 19 Sep 2016 09:58:14 with FlexibleSUSY
 * 1.7.0 (git commit: 5938cc5e9320fd7a22b1a853dc2285c56e40a49f) and SARAH 4.9.1 .
 */

#ifndef MRSSM_MASS_EIGENSTATES_H
#define MRSSM_MASS_EIGENSTATES_H

#include "MRSSM_two_scale_soft_parameters.hpp"
#include "MRSSM_physical.hpp"
#include "MRSSM_info.hpp"
#include "two_loop_corrections.hpp"
#include "error.hpp"
#include "problems.hpp"
#include "config.h"

#include <iosfwd>
#include <string>

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

class EWSB_solver;
/**
 * @class MRSSM_mass_eigenstates
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
class MRSSM_mass_eigenstates : public MRSSM_soft_parameters {
public:
   explicit MRSSM_mass_eigenstates(const MRSSM_input_parameters& input_ = MRSSM_input_parameters());
   virtual ~MRSSM_mass_eigenstates();

   /// number of EWSB equations
   static const std::size_t number_of_ewsb_equations = 4;

   void calculate_DRbar_masses();
   void calculate_DRbar_parameters();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   virtual void clear();
   void clear_DRbar_parameters();
   Eigen::ArrayXd get_DRbar_masses() const;
   void do_calculate_sm_pole_masses(bool);
   bool do_calculate_sm_pole_masses() const;
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_DRbar_masses();
   void reorder_pole_masses();
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(unsigned);
   void set_two_loop_corrections(const Two_loop_corrections&);
   const Two_loop_corrections& get_two_loop_corrections() const;
   void set_DRbar_masses(const Eigen::ArrayXd&);
   void set_number_of_ewsb_iterations(std::size_t);
   void set_number_of_mass_iterations(std::size_t);
   std::size_t get_number_of_ewsb_iterations() const;
   std::size_t get_number_of_mass_iterations() const;
   void set_pole_mass_loop_order(unsigned);
   unsigned get_pole_mass_loop_order() const;
   void set_physical(const MRSSM_physical&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   const MRSSM_physical& get_physical() const;
   MRSSM_physical& get_physical();
   const Problems<MRSSM_info::NUMBER_OF_PARTICLES>& get_problems() const;
   Problems<MRSSM_info::NUMBER_OF_PARTICLES>& get_problems();
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   void calculate_spectrum();
   void clear_problems();
   std::string name() const;
   void run_to(double scale, double eps = -1.0);
   void print(std::ostream& out = std::cout) const;
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
   const std::complex<double>& get_ZN1(int i, int k) const { return ZN1(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN2() const { return ZN2; }
   const std::complex<double>& get_ZN2(int i, int k) const { return ZN2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM1() const { return UM1; }
   const std::complex<double>& get_UM1(int i, int k) const { return UM1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP1() const { return UP1; }
   const std::complex<double>& get_UP1(int i, int k) const { return UP1(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM2() const { return UM2; }
   const std::complex<double>& get_UM2(int i, int k) const { return UM2(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP2() const { return UP2; }
   const std::complex<double>& get_UP2(int i, int k) const { return UP2(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZEL() const { return ZEL; }
   const std::complex<double>& get_ZEL(int i, int k) const { return ZEL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZER() const { return ZER; }
   const std::complex<double>& get_ZER(int i, int k) const { return ZER(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDL() const { return ZDL; }
   const std::complex<double>& get_ZDL(int i, int k) const { return ZDL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZDR() const { return ZDR; }
   const std::complex<double>& get_ZDR(int i, int k) const { return ZDR(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUL() const { return ZUL; }
   const std::complex<double>& get_ZUL(int i, int k) const { return ZUL(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_ZUR() const { return ZUR; }
   const std::complex<double>& get_ZUR(int i, int k) const { return ZUR(i,k); }
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

   std::complex<double> CpUSdconjUSdVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSdconjUSdconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSdconjUSdconjSRumSRum(unsigned gO1, unsigned gO2) const;
   double CpUSdconjUSdconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSdconjUSdconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSdbarCha1FuPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSdbarCha1FuPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFuCha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFuCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdbarChiFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdbarChiFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSdphiO(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpconjUSdSuSRum(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSdRh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSdbarGluFdPR(unsigned gO2, unsigned , unsigned gI2) const;
   double CpconjUSdbarGluFdPL(unsigned , unsigned , unsigned ) const;
   double CpconjUSdGluFdPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSdGluFdPL(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSdsigmaOSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdconjSRdpSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVGSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVPSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVZSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSdVWmSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSvconjUSvconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSvconjUSvconjSRumSRum(unsigned gO1, unsigned gO2) const;
   double CpUSvconjUSvconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSvconjUSvconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvFeCha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvFeCha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvSvAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSvFvChiPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSvVZSv(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSvSRdpSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSvconjVWmSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSuconjUSuconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSuconjUSuconjSRumSRum(unsigned gO1, unsigned gO2) const;
   double CpUSuconjUSuconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSuconjUSuconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSubarCha2FdPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSubarCha2FdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuFdCha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuFdCha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSubarChiFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSubarChiFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuSuphiO(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpconjUSuSusigmaO(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuSuRh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSubarGluFuPR(unsigned gO2, unsigned , unsigned gI2) const;
   double CpconjUSubarGluFuPL(unsigned , unsigned , unsigned ) const;
   double CpconjUSuGluFuPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSuGluFuPL(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpconjUSuSRdpSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuconjSRumSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuconjVWmSd(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVGSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVPSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSuVZSu(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSeconjUSeconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSeconjUSeconjSRumSRum(unsigned gO1, unsigned gO2) const;
   double CpUSeconjUSeconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUSeconjUSeconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSebarCha1FvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSebarCha1FvPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUSeFvCha2PR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUSeFvCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSebarChiFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSebarChiFePL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUSeconjUSeconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSeRh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUSeconjSRdpSv(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSeVWmSv(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSeVPSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUSeVZSe(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhVZVZ(unsigned gO2) const;
   std::complex<double> CpUhhconjSRdpSRdp(unsigned gO2) const;
   std::complex<double> CpUhhconjSRumSRum(unsigned gO2) const;
   std::complex<double> CpUhhconjVWmVWm(unsigned gO2) const;
   std::complex<double> CpUhhbargWmgWm(unsigned gO1) const;
   std::complex<double> CpUhhbargWmCgWmC(unsigned gO1) const;
   std::complex<double> CpUhhbargZgZ(unsigned gO1) const;
   std::complex<double> CpUhhUhhVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhconjSRumSRum(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUhhUhhconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjRhRh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarCha1Cha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarCha1Cha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarCha2Cha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarCha2Cha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjRhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjRhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhbarChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUhhSRdpHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhconjSRumHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhVZAh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUhhconjVWmHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhconjSRdpSRdp(unsigned gO2) const;
   std::complex<double> CpUAhconjSRumSRum(unsigned gO2) const;
   std::complex<double> CpUAhbargWmgWm(unsigned gO1) const;
   std::complex<double> CpUAhbargWmCgWmC(unsigned gO1) const;
   std::complex<double> CpUAhUAhVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhconjSRumSRum(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUAhUAhconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjRhRh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarCha1Cha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarCha1Cha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarCha2Cha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarCha2Cha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjRhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjRhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhbarChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUAhSRdpHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhconjSRumHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhVZhh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpUAhconjVWmHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpURhconjURhVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpURhconjURhconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpURhconjURhconjSRumSRum(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpURhconjURhconjVWmVWm(unsigned gO1, unsigned gO2) const;
   double CpconjURhVWmSRdp(unsigned gO2) const;
   double CpconjURhconjVWmSRum(unsigned gO2) const;
   std::complex<double> CpURhconjURhconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhCha2Cha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhCha2Cha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhRhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhRhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpURhconjURhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhconjHpmSRum(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpURhconjURhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpURhconjURhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpURhconjURhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpURhconjURhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpURhconjURhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpURhconjURhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjURhVZRh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjURhSRdpHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVWmVP(unsigned gO2) const;
   std::complex<double> CpconjUHpmVZVWm(unsigned gO2) const;
   std::complex<double> CpconjUHpmbargWmCgZ(unsigned gO1) const;
   std::complex<double> CpUHpmgWmCbargZ(unsigned gO2) const;
   std::complex<double> CpconjUHpmbargZgWm(unsigned gO1) const;
   std::complex<double> CpUHpmgZbargWm(unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmVZVZ(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmconjSRdpSRdp(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmconjSRumSRum(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpUHpmconjUHpmconjVWmVWm(unsigned gO1, unsigned gO2) const;
   std::complex<double> CpconjUHpmconjRhSRum(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpUHpmconjUHpmconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjRhHpm(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarCha1ChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarCha1ChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpconjUHpmbarFvFePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpconjUHpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarChiCha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmbarChiCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSRdpRh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmSRumAh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmSRumhh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSRdpAh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmconjSRdphh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVPHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpconjUHpmVZHpm(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdpVZVZ() const;
   double CpSRdpconjSRdpconjSRdpSRdp() const;
   double CpSRdpconjSRdpconjSRumSRum() const;
   double CpSRdpconjSRdpconjVWmVWm() const;
   double CpconjSRdpVPSRdp() const;
   double CpconjSRdpVZSRdp() const;
   std::complex<double> CpSRdpconjSRdpconjRhRh(unsigned gI1, unsigned gI2) const;
   double CpSRdpconjSRdpconjSvSv(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdpAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdpconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdphhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpconjHpmRh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpChiCha1PR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpChiCha1PL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpconjHpmAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpconjHpmhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdpconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdpconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRdpconjSRdpconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpconjSeSv(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpconjSdSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRdpconjVWmRh(unsigned gI2) const;
   std::complex<double> CpconjSRdpSRdpAh(unsigned gI2) const;
   std::complex<double> CpconjSRdpSRdphh(unsigned gI2) const;
   std::complex<double> CpSRumconjSRumVZVZ() const;
   double CpSRumconjSRumconjSRdpSRdp() const;
   double CpSRumconjSRumconjSRumSRum() const;
   double CpSRumconjSRumconjVWmVWm() const;
   double CpconjSRumVPSRum() const;
   double CpconjSRumVZSRum() const;
   std::complex<double> CpSRumconjSRumconjRhRh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumRhHpm(unsigned gI1, unsigned gI2) const;
   double CpSRumconjSRumconjSvSv(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRumconjSRumAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRumconjSRumconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRumconjSRumhhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumChiCha2PR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumChiCha2PL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumHpmAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumHpmhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRumconjSRumconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRumconjSRumconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpSRumconjSRumconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumconjSuSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjSRumVWmRh(unsigned gI2) const;
   std::complex<double> CpconjSRumSRumAh(unsigned gI2) const;
   std::complex<double> CpconjSRumSRumhh(unsigned gI2) const;
   double CpsigmaOsigmaOphiOphiO() const;
   std::complex<double> CpsigmaOVGsigmaO() const;
   std::complex<double> CpsigmaOconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpsigmaOconjSuSu(unsigned gI1, unsigned gI2) const;
   double CpsigmaObarGluGluPR(unsigned , unsigned ) const;
   double CpsigmaObarGluGluPL(unsigned , unsigned ) const;
   double CpphiOphiOsigmaOsigmaO() const;
   std::complex<double> CpphiOVGphiO() const;
   std::complex<double> CpphiOconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpphiOconjSuSu(unsigned gI1, unsigned gI2) const;
   double CpphiObarGluGluPR(unsigned , unsigned ) const;
   double CpphiObarGluGluPL(unsigned , unsigned ) const;
   double CpVZbargWmgWm() const;
   double CpVZbargWmCgWmC() const;
   double CpVZconjSRdpSRdp() const;
   double CpVZconjSRumSRum() const;
   std::complex<double> CpVZVZconjSRdpSRdp() const;
   std::complex<double> CpVZVZconjSRumSRum() const;
   double CpVZconjVWmVWm() const;
   std::complex<double> CpVZVZconjRhRh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjRhRh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarCha1Cha1PL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarCha1Cha1PR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarCha2Cha2PL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarCha2Cha2PR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const;
   double CpVZconjSvSv(unsigned gI1, unsigned gI2) const;
   double CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFeFePL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFeFePR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const;
   double CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const;
   double CpVZbarFvFvPR(unsigned , unsigned ) const;
   std::complex<double> CpVZVZAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZhhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZhhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarChiChiPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZbarChiChiPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVZVZhh(unsigned gI2) const;
   std::complex<double> CpVZconjVWmHpm(unsigned gI2) const;
   double CpVZVZconjVWmVWm1() const;
   double CpVZVZconjVWmVWm2() const;
   double CpVZVZconjVWmVWm3() const;
   double CpconjVWmbargPgWm() const;
   double CpconjVWmbargWmCgP() const;
   double CpconjVWmbargWmCgZ() const;
   double CpconjVWmbargZgWm() const;
   double CpVWmconjVWmconjSRdpSRdp() const;
   double CpVWmconjVWmconjSRumSRum() const;
   double CpconjVWmVWmVP() const;
   double CpconjVWmVZVWm() const;
   std::complex<double> CpconjVWmconjRhSRum(unsigned gI1) const;
   std::complex<double> CpVWmconjVWmconjRhRh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmbarCha1ChiPL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmbarCha1ChiPR(unsigned gI1, unsigned gI2) const;
   double CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const;
   double CpconjVWmbarFuFdPR(unsigned , unsigned ) const;
   std::complex<double> CpconjVWmbarFvFePL(unsigned gI1, unsigned gI2) const;
   double CpconjVWmbarFvFePR(unsigned , unsigned ) const;
   std::complex<double> CpconjVWmconjSvSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmbarChiCha2PL(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmbarChiCha2PR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpconjVWmconjSRdpRh(unsigned gI2) const;
   std::complex<double> CpconjVWmVPHpm(unsigned gI2) const;
   std::complex<double> CpconjVWmVWmhh(unsigned gI2) const;
   std::complex<double> CpconjVWmVZHpm(unsigned gI2) const;
   double CpVWmconjVWmVPVP1() const;
   double CpVWmconjVWmVPVP2() const;
   double CpVWmconjVWmVPVP3() const;
   double CpVWmconjVWmVZVZ1() const;
   double CpVWmconjVWmVZVZ2() const;
   double CpVWmconjVWmVZVZ3() const;
   double CpVWmconjVWmconjVWmVWm1() const;
   double CpVWmconjVWmconjVWmVWm2() const;
   double CpVWmconjVWmconjVWmVWm3() const;
   std::complex<double> CpbarUChibarCha1SRdpPL(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpbarUChibarCha1SRdpPR(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpbarUChibarCha2SRumPL(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpbarUChibarCha2SRumPR(unsigned gO1, unsigned gI1) const;
   double CpbarUChibarFvSvPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUChibarFvSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarFdSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarFdSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarFeSePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarFeSePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarFuSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarFuSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarChiRhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChibarChiRhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjHpmCha2PL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjHpmCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiHpmCha1PL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiHpmCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUChiVWmCha1PR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChiVWmCha1PL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChiconjVWmCha2PR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChiconjVWmCha2PL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUChiVZChiPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUChiVZChiPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barCha2RhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barCha2RhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1Cha1AhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1Cha1AhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barFeSvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barFeSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barFdSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barFdSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1hhCha1PL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1hhCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1conjHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1conjHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha1barChiSRdpPL(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpbarUCha1barChiSRdpPR(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpbarUCha1conjSdFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUCha1conjSdFuPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUCha1conjSeFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUCha1conjSeFvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUCha1VPCha1PR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUCha1VPCha1PL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUCha1VZCha1PR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUCha1VZCha1PL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUCha1conjVWmChiPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUCha1conjVWmChiPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUCha2barCha1RhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2barCha1RhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2Cha2AhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2Cha2AhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2barFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2barFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUCha2barFvSePL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUCha2barFvSePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2hhCha2PL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2hhCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2HpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2HpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUCha2barChiSRumPL(unsigned gO2, unsigned gI1) const;
   std::complex<double> CpbarUCha2barChiSRumPR(unsigned gO1, unsigned gI1) const;
   std::complex<double> CpbarUCha2conjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUCha2conjSuFdPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUCha2VPCha2PR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUCha2VPCha2PL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUCha2VZCha2PR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUCha2VZCha2PL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUCha2VWmChiPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUCha2VWmChiPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFebarCha1SvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFebarCha1SvPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarUFeHpmFvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFebarChiSePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFebarChiSePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const;
   double CpbarUFeVWmFvPR(unsigned , unsigned ) const;
   double CpbarUFeVWmFvPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFeVZFePR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFdbarCha1SuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdbarCha1SuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdbarChiSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdbarChiSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUFdSuCha2PL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFdSuCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   double CpbarUFdSdGluPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFdVWmFuPR(unsigned , unsigned ) const;
   std::complex<double> CpbarUFdVWmFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFdbarGluSdPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFdbarGluSdPR(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpbarUFubarCha2SdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFubarCha2SdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFubarChiSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFubarChiSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarUFuSdCha1PL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFuSdCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarUFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   double CpbarUFuSuGluPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const;
   std::complex<double> CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFuconjVWmFdPR(unsigned , unsigned ) const;
   std::complex<double> CpbarUFuconjVWmFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarUFubarGluSuPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarUFubarGluSuPR(unsigned gO1, unsigned , unsigned gI2) const;
   double CpbarGlubarFdSdPL(unsigned , unsigned ) const;
   std::complex<double> CpbarGlubarFdSdPR(unsigned gI1, unsigned gI2) const;
   double CpbarGlubarFuSuPL(unsigned , unsigned ) const;
   std::complex<double> CpbarGlubarFuSuPR(unsigned gI1, unsigned gI2) const;
   double CpbarGluconjSdFdPL(unsigned , unsigned ) const;
   std::complex<double> CpbarGluconjSdFdPR(unsigned gI1, unsigned gI2) const;
   double CpbarGluconjSuFuPL(unsigned , unsigned ) const;
   std::complex<double> CpbarGluconjSuFuPR(unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarGluphiOGluPL() const;
   std::complex<double> CpbarGluphiOGluPR() const;
   double CpbarGlusigmaOGluPL() const;
   double CpbarGlusigmaOGluPR() const;
   std::complex<double> CpbarGluVGGluPR() const;
   std::complex<double> CpbarGluVGGluPL() const;
   std::complex<double> CpbarFebarCha1SvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFebarCha1SvPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   double CpbarFeHpmFvPR(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFebarChiSePL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFebarChiSePR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarFeVWmFvPR(unsigned , unsigned ) const;
   std::complex<double> CpbarFeVWmFvPL(unsigned gO1, unsigned gI2) const;
   double CpbarFeVZFePR(unsigned gO2, unsigned gI2) const;
   double CpbarFeVZFePL(unsigned gO1, unsigned gI2) const;
   std::complex<double> CpbarFdbarCha1SuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdbarCha1SuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdbarChiSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdbarChiSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarFdSuCha2PL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFdSuCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   double CpbarFdSdGluPR(unsigned , unsigned , unsigned ) const;
   double CpbarFdVWmFuPR(unsigned , unsigned ) const;
   std::complex<double> CpbarFdVWmFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const;
   double CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarFdbarGluSdPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFdbarGluSdPR(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> CpbarFubarCha2SdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFubarCha2SdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFubarChiSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFubarChiSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   double CpbarFuSdCha1PL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFuSdCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const;
   std::complex<double> CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const;
   double CpbarFuSuGluPR(unsigned , unsigned , unsigned ) const;
   double CpbarFuVPFuPR(unsigned gO2, unsigned gI2) const;
   double CpbarFuVPFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFuVZFuPR(unsigned gO2, unsigned gI2) const;
   double CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const;
   double CpbarFuconjVWmFdPR(unsigned , unsigned ) const;
   std::complex<double> CpbarFuconjVWmFdPL(unsigned gO1, unsigned gI2) const;
   double CpbarFubarGluSuPL(unsigned , unsigned , unsigned ) const;
   std::complex<double> CpbarFubarGluSuPR(unsigned gO1, unsigned , unsigned gI2) const;
   std::complex<double> self_energy_Sd(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Sv(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Su(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Se(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_hh(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Ah(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Rh(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Hpm(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_SRdp(double p ) const;
   std::complex<double> self_energy_SRum(double p ) const;
   std::complex<double> self_energy_sigmaO(double p ) const;
   std::complex<double> self_energy_phiO(double p ) const;
   std::complex<double> self_energy_VZ(double p ) const;
   std::complex<double> self_energy_VWm(double p ) const;
   std::complex<double> self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha1_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha1_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha1_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha2_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha2_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Cha2_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Glu_1(double p ) const;
   std::complex<double> self_energy_Glu_PR(double p ) const;
   std::complex<double> self_energy_Glu_PL(double p ) const;
   std::complex<double> self_energy_VZ_heavy(double p ) const;
   std::complex<double> self_energy_VWm_heavy(double p ) const;
   std::complex<double> self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const;
   std::complex<double> tadpole_hh(unsigned gO1) const;
   std::complex<double> tadpole_phiO() const;


   /// calculates the tadpoles at current loop order
   void tadpole_equations(double[number_of_ewsb_equations]) const;
   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole_equations() const;





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
   struct EWSB_args {
      MRSSM_mass_eigenstates* model;
      unsigned ewsb_loop_order;
   };

   class EEWSBStepFailed : public Error {
   public:
      virtual ~EEWSBStepFailed() {}
      virtual std::string what() const { return "Could not perform EWSB step."; }
   };

   std::size_t number_of_ewsb_iterations;
   std::size_t number_of_mass_iterations;
   unsigned ewsb_loop_order;
   unsigned pole_mass_loop_order;
   bool calculate_sm_pole_masses; ///< switch to calculate the pole masses of the Standard Model particles
   bool force_output;             ///< switch to force output of pole masses
   double precision;              ///< RG running precision
   double ewsb_iteration_precision;
   MRSSM_physical physical; ///< contains the pole masses and mixings
   Problems<MRSSM_info::NUMBER_OF_PARTICLES> problems;
   Two_loop_corrections two_loop_corrections; ///< used 2-loop corrections

   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(unsigned);
   int solve_ewsb_iteratively_with(EWSB_solver*, const Eigen::Matrix<double, number_of_ewsb_equations, 1>&);
   int solve_ewsb_tree_level_custom();
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_initial_guess();
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_step() const;
   static int ewsb_step(const gsl_vector*, void*, gsl_vector*);
   static int tadpole_equations(const gsl_vector*, void*, gsl_vector*);
   void copy_DRbar_masses_to_pole_masses();

   // Passarino-Veltman loop functions
   double A0(double) const;
   double B0(double, double, double) const;
   double B1(double, double, double) const;
   double B00(double, double, double) const;
   double B22(double, double, double) const;
   double H0(double, double, double) const;
   double F0(double, double, double) const;
   double G0(double, double, double) const;

   // DR-bar masses
   double MVG;
   double MGlu;
   Eigen::Array<double,3,1> MFv;
   double MSRdp;
   double MSRum;
   double MsigmaO;
   double MphiO;
   Eigen::Array<double,6,1> MSd;
   Eigen::Array<double,3,1> MSv;
   Eigen::Array<double,6,1> MSu;
   Eigen::Array<double,6,1> MSe;
   Eigen::Array<double,4,1> Mhh;
   Eigen::Array<double,4,1> MAh;
   Eigen::Array<double,2,1> MRh;
   Eigen::Array<double,4,1> MHpm;
   Eigen::Array<double,4,1> MChi;
   Eigen::Array<double,2,1> MCha1;
   Eigen::Array<double,2,1> MCha2;
   Eigen::Array<double,3,1> MFe;
   Eigen::Array<double,3,1> MFd;
   Eigen::Array<double,3,1> MFu;
   double MVWm;
   double MVP;
   double MVZ;

   // DR-bar mixing matrices
   Eigen::Matrix<double,6,6> ZD;
   Eigen::Matrix<double,3,3> ZV;
   Eigen::Matrix<double,6,6> ZU;
   Eigen::Matrix<double,6,6> ZE;
   Eigen::Matrix<double,4,4> ZH;
   Eigen::Matrix<double,4,4> ZA;
   Eigen::Matrix<double,2,2> ZHR;
   Eigen::Matrix<double,4,4> ZP;
   Eigen::Matrix<std::complex<double>,4,4> ZN1;
   Eigen::Matrix<std::complex<double>,4,4> ZN2;
   Eigen::Matrix<std::complex<double>,2,2> UM1;
   Eigen::Matrix<std::complex<double>,2,2> UP1;
   Eigen::Matrix<std::complex<double>,2,2> UM2;
   Eigen::Matrix<std::complex<double>,2,2> UP2;
   Eigen::Matrix<std::complex<double>,3,3> ZEL;
   Eigen::Matrix<std::complex<double>,3,3> ZER;
   Eigen::Matrix<std::complex<double>,3,3> ZDL;
   Eigen::Matrix<std::complex<double>,3,3> ZDR;
   Eigen::Matrix<std::complex<double>,3,3> ZUL;
   Eigen::Matrix<std::complex<double>,3,3> ZUR;
   Eigen::Matrix<double,2,2> ZZ;

   // phases

};

std::ostream& operator<<(std::ostream&, const MRSSM_mass_eigenstates&);

} // namespace flexiblesusy

#endif
