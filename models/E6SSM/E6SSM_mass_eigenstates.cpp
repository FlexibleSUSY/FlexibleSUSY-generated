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
 * @file E6SSM_mass_eigenstates.cpp
 * @brief implementation of the E6SSM model class
 *
 * Contains the definition of the E6SSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "E6SSM_mass_eigenstates.hpp"
#include "E6SSM_ewsb_solver_interface.hpp"
#include "config.h"
#include "eigen_utils.hpp"
#include "ewsb_solver.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "loop_libraries/loop_library.hpp"
#include "raii.hpp"

#ifdef ENABLE_THREADS
#include "thread_pool.hpp"
#endif

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "E6SSM_two_scale_ewsb_solver.hpp"
#endif

#include "sfermions.hpp"
#include "mssm_twoloophiggs.hpp"
#include "nmssm_twoloophiggs.hpp"





#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>

namespace flexiblesusy {

#define STRINGIFY(s) XSTRINGIFY(s)
#define XSTRINGIFY(s) #s
#define CLASSNAME E6SSM_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS       loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS       loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT       loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU   loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION            loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS    loop_corrections.higgs_at_as_as
#define HIGGS_3LOOP_CORRECTION_AB_AS_AS    loop_corrections.higgs_ab_as_as
#define HIGGS_3LOOP_SCHEME                 loop_corrections.higgs_3L_scheme
#define HIGGS_3LOOP_CORRECTION_AT_AT_AS    loop_corrections.higgs_at_at_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AT    loop_corrections.higgs_at_at_at
#define HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS loop_corrections.higgs_at_as_as_as

CLASSNAME::CLASSNAME(const E6SSM_input_parameters& input_)
   : E6SSM_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new E6SSM_ewsb_solver<Two_scale>())
#endif
{
}

std::unique_ptr<E6SSM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<E6SSM_mass_eigenstates>(*this);
}

void CLASSNAME::do_calculate_sm_pole_masses(bool flag)
{
   calculate_sm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_sm_pole_masses() const
{
   return calculate_sm_pole_masses;
}

void CLASSNAME::do_calculate_bsm_pole_masses(bool flag)
{
   calculate_bsm_pole_masses = flag;
}

bool CLASSNAME::do_calculate_bsm_pole_masses() const
{
   return calculate_bsm_pole_masses;
}

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

void CLASSNAME::set_ewsb_loop_order(int loop_order)
{
   ewsb_loop_order = loop_order;
   if (ewsb_solver) {
      ewsb_solver->set_loop_order(ewsb_loop_order);
   }
}

void CLASSNAME::set_loop_corrections(const Loop_corrections& loop_corrections_)
{
   loop_corrections = loop_corrections_;
}

const Loop_corrections& CLASSNAME::get_loop_corrections() const
{
   return loop_corrections;
}

void CLASSNAME::set_threshold_corrections(const Threshold_corrections& tc)
{
   threshold_corrections = tc;
}

const Threshold_corrections& CLASSNAME::get_threshold_corrections() const
{
   return threshold_corrections;
}

int CLASSNAME::get_number_of_ewsb_iterations() const
{
   return static_cast<int>(std::abs(-log10(ewsb_iteration_precision) * 10));
}

int CLASSNAME::get_number_of_mass_iterations() const
{
   return static_cast<int>(std::abs(-log10(precision) * 10));
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
   if (ewsb_solver) {
      ewsb_solver->set_precision(precision_);
   }
}

void CLASSNAME::set_pole_mass_loop_order(int loop_order)
{
   pole_mass_loop_order = loop_order;
}

int CLASSNAME::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
   if (ewsb_solver) {
      ewsb_solver->set_precision(precision);
   }
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double CLASSNAME::get_precision() const
{
   return precision;
}

double CLASSNAME::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const E6SSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

E6SSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const E6SSM_physical& physical_)
{
   physical = physical_;
}

const Problems& CLASSNAME::get_problems() const
{
   return problems;
}

Problems& CLASSNAME::get_problems()
{
   return problems;
}

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<E6SSM_ewsb_solver_interface>& solver)
{
   ewsb_solver = solver;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   const auto tadpole_(tadpole_equations());
   std::copy(tadpole_.data(), tadpole_.data() + number_of_ewsb_equations, tadpole);
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double,CLASSNAME::number_of_ewsb_equations,1> CLASSNAME::tadpole_equations() const
{
   Eigen::Matrix<double,number_of_ewsb_equations,1> tadpole(
      Eigen::Matrix<double,number_of_ewsb_equations,1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();
   tadpole[2] = get_ewsb_eq_hh_3();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop(0));
      tadpole[1] -= Re(tadpole_hh_1loop(1));
      tadpole[2] -= Re(tadpole_hh_1loop(2));

      if (ewsb_loop_order > 1) {
         const auto tadpole_2l(tadpole_hh_2loop());
         tadpole[0] -= tadpole_2l(0);
         tadpole[1] -= tadpole_2l(1);
         tadpole[2] -= tadpole_2l(2);

      }
   }

   return tadpole;
}

/**
 * This function returns the vector of tadpoles, each divided by the
 * corresponding VEV.  Thus, the returned tadpoles have the dimension
 * GeV^2 each.
 *
 * @return vector of tadpoles
 */
Eigen::Matrix<double,CLASSNAME::number_of_ewsb_equations,1> CLASSNAME::tadpole_equations_over_vevs() const
{
   auto tadpole = tadpole_equations();

   tadpole[0] /= vd;
   tadpole[1] /= vu;
   tadpole[2] /= vs;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;

   
   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_ms2 = ms2;

   mHd2 = Re((0.0125*(28.284271247461902*vs*vu*Conj(TLambdax) - 6*Cube(vd)*Sqr(g1)
      - 10*Cube(vd)*Sqr(g2) - 9*Cube(vd)*Sqr(gN) - 40*vd*AbsSqr(Lambdax)*Sqr(vs) +
      15*vd*Sqr(gN)*Sqr(vs) - 40*vd*AbsSqr(Lambdax)*Sqr(vu) + 6*vd*Sqr(g1)*Sqr(vu)
      + 10*vd*Sqr(g2)*Sqr(vu) - 6*vd*Sqr(gN)*Sqr(vu) + 28.284271247461902*vs*vu*
      TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vs*Conj(TLambdax) - 3*Cube(vu)*Sqr(g1)
      - 5*Cube(vu)*Sqr(g2) - 2*Cube(vu)*Sqr(gN) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) +
      3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 3*vu*Sqr(gN)*Sqr(vd) - 20*vu*
      AbsSqr(Lambdax)*Sqr(vs) + 5*vu*Sqr(gN)*Sqr(vs) + 14.142135623730951*vd*vs*
      TLambdax))/vu);
   ms2 = Re((0.0625*(5.656854249492381*vd*vu*Conj(TLambdax) - 5*Cube(vs)*Sqr(gN) -
      8*vs*AbsSqr(Lambdax)*Sqr(vd) + 3*vs*Sqr(gN)*Sqr(vd) - 8*vs*AbsSqr(Lambdax)*
      Sqr(vu) + 2*vs*Sqr(gN)*Sqr(vu) + 5.656854249492381*vd*vu*TLambdax))/vs);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      ms2 = old_ms2;
      error = EWSB_solver::FAIL;
   }

   return error;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   if (!ewsb_solver) {
      throw SetupError(STRINGIFY(CLASSNAME) "::solve_ewsb_tree_level: "
                       "no EWSB solver set");
   }

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(0);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

int CLASSNAME::solve_ewsb_one_loop()
{
   if (!ewsb_solver) {
      throw SetupError(STRINGIFY(CLASSNAME) "::solve_ewsb_one_loop: "
                       "no EWSB solver set");
   }

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(1);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

int CLASSNAME::solve_ewsb()
{
   if (!ewsb_solver) {
      throw SetupError(STRINGIFY(CLASSNAME) "::solve_ewsb: "
                       "no EWSB solver set");
   }

   VERBOSE_MSG("\t\tSolving E6SSM EWSB at " << ewsb_loop_order << "-loop order");

   const int old_loop_order = ewsb_solver->get_loop_order();
   const auto save_loop_order = make_raii_guard(
      [this, old_loop_order] () {
         this->ewsb_solver->set_loop_order(old_loop_order);
      });

   const int old_iterations = ewsb_solver->get_number_of_iterations();
   const auto save_iterations = make_raii_guard(
      [this, old_iterations] () {
         this->ewsb_solver->set_number_of_iterations(old_iterations);
      });

   const double old_precision = ewsb_solver->get_precision();
   const auto save_precision = make_raii_guard(
      [this, old_precision] () {
         this->ewsb_solver->set_precision(old_precision);
      });

   ewsb_solver->set_loop_order(ewsb_loop_order);
   ewsb_solver->set_number_of_iterations(get_number_of_ewsb_iterations());
   ewsb_solver->set_precision(ewsb_iteration_precision);

   return ewsb_solver->solve(*this);
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "E6SSM\n"
           "========================================\n";
   E6SSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MChaP = " << MChaP << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSDX = " << MSDX.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFDX = " << MFDX.transpose() << '\n';
   ostr << "MSHI0 = " << MSHI0.transpose() << '\n';
   ostr << "MSHIp = " << MSHIp.transpose() << '\n';
   ostr << "MChaI = " << MChaI.transpose() << '\n';
   ostr << "MChiI = " << MChiI.transpose() << '\n';
   ostr << "MSSI0 = " << MSSI0.transpose() << '\n';
   ostr << "MFSI = " << MFSI.transpose() << '\n';
   ostr << "MSHp0 = " << MSHp0.transpose() << '\n';
   ostr << "MSHpp = " << MSHpp.transpose() << '\n';
   ostr << "MChiP = " << MChiP.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MVZp = " << MVZp << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZDX = " << ZDX << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
   ostr << "ZDXL = " << ZDXL << '\n';
   ostr << "ZDXR = " << ZDXR << '\n';
   ostr << "UHI0 = " << UHI0 << '\n';
   ostr << "UHIp = " << UHIp << '\n';
   ostr << "ZMI = " << ZMI << '\n';
   ostr << "ZPI = " << ZPI << '\n';
   ostr << "ZNI = " << ZNI << '\n';
   ostr << "ZSSI = " << ZSSI << '\n';
   ostr << "ZFSI = " << ZFSI << '\n';
   ostr << "UHp0 = " << UHp0 << '\n';
   ostr << "UHpp = " << UHpp << '\n';
   ostr << "ZNp = " << ZNp << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 * @note: They take squared arguments!
 */

double CLASSNAME::A0(double m) const noexcept
{
   return Loop_library::get().A0(m, Sqr(get_scale())).real();
}

double CLASSNAME::B0(double p, double m1, double m2) const noexcept
{
   return Loop_library::get().B0(p, m1, m2, Sqr(get_scale())).real();
}

double CLASSNAME::B1(double p, double m1, double m2) const noexcept
{
   return Loop_library::get().B1(p, m1, m2, Sqr(get_scale())).real();
}

double CLASSNAME::B00(double p, double m1, double m2) const noexcept
{
   return Loop_library::get().B00(p, m1, m2, Sqr(get_scale())).real();
}

double CLASSNAME::B22(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return (
      Loop_library::get().B00(p, m1, m2, scl2) - Loop_library::get().A0(m1, scl2)/4.0 -
      Loop_library::get().A0(m2, scl2)/4.0
   ).real();
}

double CLASSNAME::H0(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return 4.0*Loop_library::get().B00(p, m1, m2, scl2).real() + CLASSNAME::G0(p, m1, m2);
}

double CLASSNAME::F0(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return (
      Loop_library::get().A0(m1, scl2) - 2.0*Loop_library::get().A0(m2, scl2)
      - (2.0*p + 2.0*m1 - m2) * Loop_library::get().B0(p, m1, m2, scl2)
   ).real();
}

double CLASSNAME::G0(double p, double m1, double m2) const noexcept
{
   const double scl2 = Sqr(get_scale());
   return (
      (p - m1 - m2) * Loop_library::get().B0(p, m1, m2, scl2)
      - Loop_library::get().A0(m1, scl2) - Loop_library::get().A0(m2, scl2)
   ).real();
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto save_mHd2_raii = make_raii_save(mHd2);
   const auto save_mHu2_raii = make_raii_save(mHu2);
   const auto save_ms2_raii = make_raii_save(ms2);

   solve_ewsb_tree_level_custom();

   calculate_MVPVZVZp();
   calculate_MVWm();
   calculate_MChiP();
   calculate_MSHpp();
   calculate_MSHp0();
   calculate_MFSI();
   calculate_MSSI0();
   calculate_MChiI();
   calculate_MChaI();
   calculate_MSHIp();
   calculate_MSHI0();
   calculate_MFDX();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSDX();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
   calculate_MChaP();
   calculate_MFv();
   calculate_MGlu();
   calculate_MVG();

}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 31u));

   if (calculate_bsm_pole_masses) {
      tp.run_task([this] () { calculate_MAh_pole(); });
      tp.run_task([this] () { calculate_MCha_pole(); });
      tp.run_task([this] () { calculate_MChaI_pole(); });
      tp.run_task([this] () { calculate_MChaP_pole(); });
      tp.run_task([this] () { calculate_MChi_pole(); });
      tp.run_task([this] () { calculate_MChiI_pole(); });
      tp.run_task([this] () { calculate_MChiP_pole(); });
      tp.run_task([this] () { calculate_MFDX_pole(); });
      tp.run_task([this] () { calculate_MFSI_pole(); });
      tp.run_task([this] () { calculate_MGlu_pole(); });
      tp.run_task([this] () { calculate_Mhh_pole(); });
      tp.run_task([this] () { calculate_MHpm_pole(); });
      tp.run_task([this] () { calculate_MSd_pole(); });
      tp.run_task([this] () { calculate_MSDX_pole(); });
      tp.run_task([this] () { calculate_MSe_pole(); });
      tp.run_task([this] () { calculate_MSHI0_pole(); });
      tp.run_task([this] () { calculate_MSHIp_pole(); });
      tp.run_task([this] () { calculate_MSHp0_pole(); });
      tp.run_task([this] () { calculate_MSHpp_pole(); });
      tp.run_task([this] () { calculate_MSSI0_pole(); });
      tp.run_task([this] () { calculate_MSu_pole(); });
      tp.run_task([this] () { calculate_MSv_pole(); });
      tp.run_task([this] () { calculate_MVZp_pole(); });
   }

   if (calculate_sm_pole_masses) {
      tp.run_task([this] () { calculate_MVG_pole(); });
      tp.run_task([this] () { calculate_MFv_pole(); });
      tp.run_task([this] () { calculate_MVP_pole(); });
      tp.run_task([this] () { calculate_MVZ_pole(); });
      tp.run_task([this] () { calculate_MFe_pole(); });
      tp.run_task([this] () { calculate_MFd_pole(); });
      tp.run_task([this] () { calculate_MFu_pole(); });
      tp.run_task([this] () { calculate_MVWm_pole(); });
   }

#else
   if (calculate_bsm_pole_masses) {
      calculate_MAh_pole();
      calculate_MCha_pole();
      calculate_MChaI_pole();
      calculate_MChaP_pole();
      calculate_MChi_pole();
      calculate_MChiI_pole();
      calculate_MChiP_pole();
      calculate_MFDX_pole();
      calculate_MFSI_pole();
      calculate_MGlu_pole();
      calculate_Mhh_pole();
      calculate_MHpm_pole();
      calculate_MSd_pole();
      calculate_MSDX_pole();
      calculate_MSe_pole();
      calculate_MSHI0_pole();
      calculate_MSHIp_pole();
      calculate_MSHp0_pole();
      calculate_MSHpp_pole();
      calculate_MSSI0_pole();
      calculate_MSu_pole();
      calculate_MSv_pole();
      calculate_MVZp_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFe_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MVWm_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MChaP) = MChaP;
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSv) = MSv;
   PHYSICAL(ZV) = ZV;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(MSDX) = MSDX;
   PHYSICAL(ZDX) = ZDX;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MCha) = MCha;
   PHYSICAL(UM) = UM;
   PHYSICAL(UP) = UP;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(ZEL) = ZEL;
   PHYSICAL(ZER) = ZER;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(ZDL) = ZDL;
   PHYSICAL(ZDR) = ZDR;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(ZUL) = ZUL;
   PHYSICAL(ZUR) = ZUR;
   PHYSICAL(MFDX) = MFDX;
   PHYSICAL(ZDXL) = ZDXL;
   PHYSICAL(ZDXR) = ZDXR;
   PHYSICAL(MSHI0) = MSHI0;
   PHYSICAL(UHI0) = UHI0;
   PHYSICAL(MSHIp) = MSHIp;
   PHYSICAL(UHIp) = UHIp;
   PHYSICAL(MChaI) = MChaI;
   PHYSICAL(ZMI) = ZMI;
   PHYSICAL(ZPI) = ZPI;
   PHYSICAL(MChiI) = MChiI;
   PHYSICAL(ZNI) = ZNI;
   PHYSICAL(MSSI0) = MSSI0;
   PHYSICAL(ZSSI) = ZSSI;
   PHYSICAL(MFSI) = MFSI;
   PHYSICAL(ZFSI) = ZFSI;
   PHYSICAL(MSHp0) = MSHp0;
   PHYSICAL(UHp0) = UHp0;
   PHYSICAL(MSHpp) = MSHpp;
   PHYSICAL(UHpp) = UHpp;
   PHYSICAL(MChiP) = MChiP;
   PHYSICAL(ZNp) = ZNp;
   PHYSICAL(MVWm) = MVWm;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(MVZp) = MVZp;
   PHYSICAL(ZZ) = ZZ;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(1, MVZp, MAh, ZA);
   move_goldstone_to(0, MVWm, MHpm, ZP);

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{
   move_goldstone_to(0, MVZ, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(1, MVZp, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(0, MVWm, PHYSICAL(MHpm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::Sd); }
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::Se); }
   if (PHYSICAL(MSDX).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::SDX); }
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::Hpm); }
   if (PHYSICAL(MSHI0).tail<4>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::SHI0); }
   if (PHYSICAL(MSHIp).tail<4>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::SHIp); }
   if (PHYSICAL(MSSI0).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::SSI0); }
   if (PHYSICAL(MSHp0).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::SHp0); }
   if (PHYSICAL(MSHpp).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(E6SSM_info::SHpp); }
}

/**
 * calculates spectrum for model once the DRbar parameters at
 * at low energies are known
 */
void CLASSNAME::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   // move goldstone bosons to the front
   reorder_DRbar_masses();
   if (pole_mass_loop_order == 0)
      copy_DRbar_masses_to_pole_masses();
   else
      reorder_pole_masses();

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void CLASSNAME::clear_DRbar_parameters()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MChaP = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   MSDX = Eigen::Matrix<double,6,1>::Zero();
   ZDX = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,6,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,6,6>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFDX = Eigen::Matrix<double,3,1>::Zero();
   ZDXL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDXR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MSHI0 = Eigen::Matrix<double,4,1>::Zero();
   UHI0 = Eigen::Matrix<double,4,4>::Zero();
   MSHIp = Eigen::Matrix<double,4,1>::Zero();
   UHIp = Eigen::Matrix<double,4,4>::Zero();
   MChaI = Eigen::Matrix<double,2,1>::Zero();
   ZMI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   ZPI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MChiI = Eigen::Matrix<double,4,1>::Zero();
   ZNI = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MSSI0 = Eigen::Matrix<double,2,1>::Zero();
   ZSSI = Eigen::Matrix<double,2,2>::Zero();
   MFSI = Eigen::Matrix<double,2,1>::Zero();
   ZFSI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MSHp0 = Eigen::Matrix<double,2,1>::Zero();
   UHp0 = Eigen::Matrix<double,2,2>::Zero();
   MSHpp = Eigen::Matrix<double,2,1>::Zero();
   UHpp = Eigen::Matrix<double,2,2>::Zero();
   MChiP = Eigen::Matrix<double,2,1>::Zero();
   ZNp = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;

   PhaseGlu = std::complex<double>(1.,0.);
   PhaseFHpup = std::complex<double>(1.,0.);


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   E6SSM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MChaP = pars(5);
   MSd(0) = pars(6);
   MSd(1) = pars(7);
   MSd(2) = pars(8);
   MSd(3) = pars(9);
   MSd(4) = pars(10);
   MSd(5) = pars(11);
   MSv(0) = pars(12);
   MSv(1) = pars(13);
   MSv(2) = pars(14);
   MSu(0) = pars(15);
   MSu(1) = pars(16);
   MSu(2) = pars(17);
   MSu(3) = pars(18);
   MSu(4) = pars(19);
   MSu(5) = pars(20);
   MSe(0) = pars(21);
   MSe(1) = pars(22);
   MSe(2) = pars(23);
   MSe(3) = pars(24);
   MSe(4) = pars(25);
   MSe(5) = pars(26);
   MSDX(0) = pars(27);
   MSDX(1) = pars(28);
   MSDX(2) = pars(29);
   MSDX(3) = pars(30);
   MSDX(4) = pars(31);
   MSDX(5) = pars(32);
   Mhh(0) = pars(33);
   Mhh(1) = pars(34);
   Mhh(2) = pars(35);
   MAh(0) = pars(36);
   MAh(1) = pars(37);
   MAh(2) = pars(38);
   MHpm(0) = pars(39);
   MHpm(1) = pars(40);
   MChi(0) = pars(41);
   MChi(1) = pars(42);
   MChi(2) = pars(43);
   MChi(3) = pars(44);
   MChi(4) = pars(45);
   MChi(5) = pars(46);
   MCha(0) = pars(47);
   MCha(1) = pars(48);
   MFe(0) = pars(49);
   MFe(1) = pars(50);
   MFe(2) = pars(51);
   MFd(0) = pars(52);
   MFd(1) = pars(53);
   MFd(2) = pars(54);
   MFu(0) = pars(55);
   MFu(1) = pars(56);
   MFu(2) = pars(57);
   MFDX(0) = pars(58);
   MFDX(1) = pars(59);
   MFDX(2) = pars(60);
   MSHI0(0) = pars(61);
   MSHI0(1) = pars(62);
   MSHI0(2) = pars(63);
   MSHI0(3) = pars(64);
   MSHIp(0) = pars(65);
   MSHIp(1) = pars(66);
   MSHIp(2) = pars(67);
   MSHIp(3) = pars(68);
   MChaI(0) = pars(69);
   MChaI(1) = pars(70);
   MChiI(0) = pars(71);
   MChiI(1) = pars(72);
   MChiI(2) = pars(73);
   MChiI(3) = pars(74);
   MSSI0(0) = pars(75);
   MSSI0(1) = pars(76);
   MFSI(0) = pars(77);
   MFSI(1) = pars(78);
   MSHp0(0) = pars(79);
   MSHp0(1) = pars(80);
   MSHpp(0) = pars(81);
   MSHpp(1) = pars(82);
   MChiP(0) = pars(83);
   MChiP(1) = pars(84);
   MVWm = pars(85);
   MVP = pars(86);
   MVZ = pars(87);
   MVZp = pars(88);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(89);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MChaP;
   pars(6) = MSd(0);
   pars(7) = MSd(1);
   pars(8) = MSd(2);
   pars(9) = MSd(3);
   pars(10) = MSd(4);
   pars(11) = MSd(5);
   pars(12) = MSv(0);
   pars(13) = MSv(1);
   pars(14) = MSv(2);
   pars(15) = MSu(0);
   pars(16) = MSu(1);
   pars(17) = MSu(2);
   pars(18) = MSu(3);
   pars(19) = MSu(4);
   pars(20) = MSu(5);
   pars(21) = MSe(0);
   pars(22) = MSe(1);
   pars(23) = MSe(2);
   pars(24) = MSe(3);
   pars(25) = MSe(4);
   pars(26) = MSe(5);
   pars(27) = MSDX(0);
   pars(28) = MSDX(1);
   pars(29) = MSDX(2);
   pars(30) = MSDX(3);
   pars(31) = MSDX(4);
   pars(32) = MSDX(5);
   pars(33) = Mhh(0);
   pars(34) = Mhh(1);
   pars(35) = Mhh(2);
   pars(36) = MAh(0);
   pars(37) = MAh(1);
   pars(38) = MAh(2);
   pars(39) = MHpm(0);
   pars(40) = MHpm(1);
   pars(41) = MChi(0);
   pars(42) = MChi(1);
   pars(43) = MChi(2);
   pars(44) = MChi(3);
   pars(45) = MChi(4);
   pars(46) = MChi(5);
   pars(47) = MCha(0);
   pars(48) = MCha(1);
   pars(49) = MFe(0);
   pars(50) = MFe(1);
   pars(51) = MFe(2);
   pars(52) = MFd(0);
   pars(53) = MFd(1);
   pars(54) = MFd(2);
   pars(55) = MFu(0);
   pars(56) = MFu(1);
   pars(57) = MFu(2);
   pars(58) = MFDX(0);
   pars(59) = MFDX(1);
   pars(60) = MFDX(2);
   pars(61) = MSHI0(0);
   pars(62) = MSHI0(1);
   pars(63) = MSHI0(2);
   pars(64) = MSHI0(3);
   pars(65) = MSHIp(0);
   pars(66) = MSHIp(1);
   pars(67) = MSHIp(2);
   pars(68) = MSHIp(3);
   pars(69) = MChaI(0);
   pars(70) = MChaI(1);
   pars(71) = MChiI(0);
   pars(72) = MChiI(1);
   pars(73) = MChiI(2);
   pars(74) = MChiI(3);
   pars(75) = MSSI0(0);
   pars(76) = MSSI0(1);
   pars(77) = MFSI(0);
   pars(78) = MFSI(1);
   pars(79) = MSHp0(0);
   pars(80) = MSHp0(1);
   pars(81) = MSHpp(0);
   pars(82) = MSHpp(1);
   pars(83) = MChiP(0);
   pars(84) = MChiP(1);
   pars(85) = MVWm;
   pars(86) = MVP;
   pars(87) = MVZ;
   pars(88) = MVZp;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   ZD(0,0) = pars(89);
   ZD(0,1) = pars(90);
   ZD(0,2) = pars(91);
   ZD(0,3) = pars(92);
   ZD(0,4) = pars(93);
   ZD(0,5) = pars(94);
   ZD(1,0) = pars(95);
   ZD(1,1) = pars(96);
   ZD(1,2) = pars(97);
   ZD(1,3) = pars(98);
   ZD(1,4) = pars(99);
   ZD(1,5) = pars(100);
   ZD(2,0) = pars(101);
   ZD(2,1) = pars(102);
   ZD(2,2) = pars(103);
   ZD(2,3) = pars(104);
   ZD(2,4) = pars(105);
   ZD(2,5) = pars(106);
   ZD(3,0) = pars(107);
   ZD(3,1) = pars(108);
   ZD(3,2) = pars(109);
   ZD(3,3) = pars(110);
   ZD(3,4) = pars(111);
   ZD(3,5) = pars(112);
   ZD(4,0) = pars(113);
   ZD(4,1) = pars(114);
   ZD(4,2) = pars(115);
   ZD(4,3) = pars(116);
   ZD(4,4) = pars(117);
   ZD(4,5) = pars(118);
   ZD(5,0) = pars(119);
   ZD(5,1) = pars(120);
   ZD(5,2) = pars(121);
   ZD(5,3) = pars(122);
   ZD(5,4) = pars(123);
   ZD(5,5) = pars(124);
   ZV(0,0) = pars(125);
   ZV(0,1) = pars(126);
   ZV(0,2) = pars(127);
   ZV(1,0) = pars(128);
   ZV(1,1) = pars(129);
   ZV(1,2) = pars(130);
   ZV(2,0) = pars(131);
   ZV(2,1) = pars(132);
   ZV(2,2) = pars(133);
   ZU(0,0) = pars(134);
   ZU(0,1) = pars(135);
   ZU(0,2) = pars(136);
   ZU(0,3) = pars(137);
   ZU(0,4) = pars(138);
   ZU(0,5) = pars(139);
   ZU(1,0) = pars(140);
   ZU(1,1) = pars(141);
   ZU(1,2) = pars(142);
   ZU(1,3) = pars(143);
   ZU(1,4) = pars(144);
   ZU(1,5) = pars(145);
   ZU(2,0) = pars(146);
   ZU(2,1) = pars(147);
   ZU(2,2) = pars(148);
   ZU(2,3) = pars(149);
   ZU(2,4) = pars(150);
   ZU(2,5) = pars(151);
   ZU(3,0) = pars(152);
   ZU(3,1) = pars(153);
   ZU(3,2) = pars(154);
   ZU(3,3) = pars(155);
   ZU(3,4) = pars(156);
   ZU(3,5) = pars(157);
   ZU(4,0) = pars(158);
   ZU(4,1) = pars(159);
   ZU(4,2) = pars(160);
   ZU(4,3) = pars(161);
   ZU(4,4) = pars(162);
   ZU(4,5) = pars(163);
   ZU(5,0) = pars(164);
   ZU(5,1) = pars(165);
   ZU(5,2) = pars(166);
   ZU(5,3) = pars(167);
   ZU(5,4) = pars(168);
   ZU(5,5) = pars(169);
   ZE(0,0) = pars(170);
   ZE(0,1) = pars(171);
   ZE(0,2) = pars(172);
   ZE(0,3) = pars(173);
   ZE(0,4) = pars(174);
   ZE(0,5) = pars(175);
   ZE(1,0) = pars(176);
   ZE(1,1) = pars(177);
   ZE(1,2) = pars(178);
   ZE(1,3) = pars(179);
   ZE(1,4) = pars(180);
   ZE(1,5) = pars(181);
   ZE(2,0) = pars(182);
   ZE(2,1) = pars(183);
   ZE(2,2) = pars(184);
   ZE(2,3) = pars(185);
   ZE(2,4) = pars(186);
   ZE(2,5) = pars(187);
   ZE(3,0) = pars(188);
   ZE(3,1) = pars(189);
   ZE(3,2) = pars(190);
   ZE(3,3) = pars(191);
   ZE(3,4) = pars(192);
   ZE(3,5) = pars(193);
   ZE(4,0) = pars(194);
   ZE(4,1) = pars(195);
   ZE(4,2) = pars(196);
   ZE(4,3) = pars(197);
   ZE(4,4) = pars(198);
   ZE(4,5) = pars(199);
   ZE(5,0) = pars(200);
   ZE(5,1) = pars(201);
   ZE(5,2) = pars(202);
   ZE(5,3) = pars(203);
   ZE(5,4) = pars(204);
   ZE(5,5) = pars(205);
   ZDX(0,0) = pars(206);
   ZDX(0,1) = pars(207);
   ZDX(0,2) = pars(208);
   ZDX(0,3) = pars(209);
   ZDX(0,4) = pars(210);
   ZDX(0,5) = pars(211);
   ZDX(1,0) = pars(212);
   ZDX(1,1) = pars(213);
   ZDX(1,2) = pars(214);
   ZDX(1,3) = pars(215);
   ZDX(1,4) = pars(216);
   ZDX(1,5) = pars(217);
   ZDX(2,0) = pars(218);
   ZDX(2,1) = pars(219);
   ZDX(2,2) = pars(220);
   ZDX(2,3) = pars(221);
   ZDX(2,4) = pars(222);
   ZDX(2,5) = pars(223);
   ZDX(3,0) = pars(224);
   ZDX(3,1) = pars(225);
   ZDX(3,2) = pars(226);
   ZDX(3,3) = pars(227);
   ZDX(3,4) = pars(228);
   ZDX(3,5) = pars(229);
   ZDX(4,0) = pars(230);
   ZDX(4,1) = pars(231);
   ZDX(4,2) = pars(232);
   ZDX(4,3) = pars(233);
   ZDX(4,4) = pars(234);
   ZDX(4,5) = pars(235);
   ZDX(5,0) = pars(236);
   ZDX(5,1) = pars(237);
   ZDX(5,2) = pars(238);
   ZDX(5,3) = pars(239);
   ZDX(5,4) = pars(240);
   ZDX(5,5) = pars(241);
   ZH(0,0) = pars(242);
   ZH(0,1) = pars(243);
   ZH(0,2) = pars(244);
   ZH(1,0) = pars(245);
   ZH(1,1) = pars(246);
   ZH(1,2) = pars(247);
   ZH(2,0) = pars(248);
   ZH(2,1) = pars(249);
   ZH(2,2) = pars(250);
   ZA(0,0) = pars(251);
   ZA(0,1) = pars(252);
   ZA(0,2) = pars(253);
   ZA(1,0) = pars(254);
   ZA(1,1) = pars(255);
   ZA(1,2) = pars(256);
   ZA(2,0) = pars(257);
   ZA(2,1) = pars(258);
   ZA(2,2) = pars(259);
   ZP(0,0) = pars(260);
   ZP(0,1) = pars(261);
   ZP(1,0) = pars(262);
   ZP(1,1) = pars(263);
   ZN(0,0) = std::complex<double>(pars(264), pars(265));
   ZN(0,1) = std::complex<double>(pars(266), pars(267));
   ZN(0,2) = std::complex<double>(pars(268), pars(269));
   ZN(0,3) = std::complex<double>(pars(270), pars(271));
   ZN(0,4) = std::complex<double>(pars(272), pars(273));
   ZN(0,5) = std::complex<double>(pars(274), pars(275));
   ZN(1,0) = std::complex<double>(pars(276), pars(277));
   ZN(1,1) = std::complex<double>(pars(278), pars(279));
   ZN(1,2) = std::complex<double>(pars(280), pars(281));
   ZN(1,3) = std::complex<double>(pars(282), pars(283));
   ZN(1,4) = std::complex<double>(pars(284), pars(285));
   ZN(1,5) = std::complex<double>(pars(286), pars(287));
   ZN(2,0) = std::complex<double>(pars(288), pars(289));
   ZN(2,1) = std::complex<double>(pars(290), pars(291));
   ZN(2,2) = std::complex<double>(pars(292), pars(293));
   ZN(2,3) = std::complex<double>(pars(294), pars(295));
   ZN(2,4) = std::complex<double>(pars(296), pars(297));
   ZN(2,5) = std::complex<double>(pars(298), pars(299));
   ZN(3,0) = std::complex<double>(pars(300), pars(301));
   ZN(3,1) = std::complex<double>(pars(302), pars(303));
   ZN(3,2) = std::complex<double>(pars(304), pars(305));
   ZN(3,3) = std::complex<double>(pars(306), pars(307));
   ZN(3,4) = std::complex<double>(pars(308), pars(309));
   ZN(3,5) = std::complex<double>(pars(310), pars(311));
   ZN(4,0) = std::complex<double>(pars(312), pars(313));
   ZN(4,1) = std::complex<double>(pars(314), pars(315));
   ZN(4,2) = std::complex<double>(pars(316), pars(317));
   ZN(4,3) = std::complex<double>(pars(318), pars(319));
   ZN(4,4) = std::complex<double>(pars(320), pars(321));
   ZN(4,5) = std::complex<double>(pars(322), pars(323));
   ZN(5,0) = std::complex<double>(pars(324), pars(325));
   ZN(5,1) = std::complex<double>(pars(326), pars(327));
   ZN(5,2) = std::complex<double>(pars(328), pars(329));
   ZN(5,3) = std::complex<double>(pars(330), pars(331));
   ZN(5,4) = std::complex<double>(pars(332), pars(333));
   ZN(5,5) = std::complex<double>(pars(334), pars(335));
   UM(0,0) = std::complex<double>(pars(336), pars(337));
   UM(0,1) = std::complex<double>(pars(338), pars(339));
   UM(1,0) = std::complex<double>(pars(340), pars(341));
   UM(1,1) = std::complex<double>(pars(342), pars(343));
   UP(0,0) = std::complex<double>(pars(344), pars(345));
   UP(0,1) = std::complex<double>(pars(346), pars(347));
   UP(1,0) = std::complex<double>(pars(348), pars(349));
   UP(1,1) = std::complex<double>(pars(350), pars(351));
   ZEL(0,0) = std::complex<double>(pars(352), pars(353));
   ZEL(0,1) = std::complex<double>(pars(354), pars(355));
   ZEL(0,2) = std::complex<double>(pars(356), pars(357));
   ZEL(1,0) = std::complex<double>(pars(358), pars(359));
   ZEL(1,1) = std::complex<double>(pars(360), pars(361));
   ZEL(1,2) = std::complex<double>(pars(362), pars(363));
   ZEL(2,0) = std::complex<double>(pars(364), pars(365));
   ZEL(2,1) = std::complex<double>(pars(366), pars(367));
   ZEL(2,2) = std::complex<double>(pars(368), pars(369));
   ZER(0,0) = std::complex<double>(pars(370), pars(371));
   ZER(0,1) = std::complex<double>(pars(372), pars(373));
   ZER(0,2) = std::complex<double>(pars(374), pars(375));
   ZER(1,0) = std::complex<double>(pars(376), pars(377));
   ZER(1,1) = std::complex<double>(pars(378), pars(379));
   ZER(1,2) = std::complex<double>(pars(380), pars(381));
   ZER(2,0) = std::complex<double>(pars(382), pars(383));
   ZER(2,1) = std::complex<double>(pars(384), pars(385));
   ZER(2,2) = std::complex<double>(pars(386), pars(387));
   ZDL(0,0) = std::complex<double>(pars(388), pars(389));
   ZDL(0,1) = std::complex<double>(pars(390), pars(391));
   ZDL(0,2) = std::complex<double>(pars(392), pars(393));
   ZDL(1,0) = std::complex<double>(pars(394), pars(395));
   ZDL(1,1) = std::complex<double>(pars(396), pars(397));
   ZDL(1,2) = std::complex<double>(pars(398), pars(399));
   ZDL(2,0) = std::complex<double>(pars(400), pars(401));
   ZDL(2,1) = std::complex<double>(pars(402), pars(403));
   ZDL(2,2) = std::complex<double>(pars(404), pars(405));
   ZDR(0,0) = std::complex<double>(pars(406), pars(407));
   ZDR(0,1) = std::complex<double>(pars(408), pars(409));
   ZDR(0,2) = std::complex<double>(pars(410), pars(411));
   ZDR(1,0) = std::complex<double>(pars(412), pars(413));
   ZDR(1,1) = std::complex<double>(pars(414), pars(415));
   ZDR(1,2) = std::complex<double>(pars(416), pars(417));
   ZDR(2,0) = std::complex<double>(pars(418), pars(419));
   ZDR(2,1) = std::complex<double>(pars(420), pars(421));
   ZDR(2,2) = std::complex<double>(pars(422), pars(423));
   ZUL(0,0) = std::complex<double>(pars(424), pars(425));
   ZUL(0,1) = std::complex<double>(pars(426), pars(427));
   ZUL(0,2) = std::complex<double>(pars(428), pars(429));
   ZUL(1,0) = std::complex<double>(pars(430), pars(431));
   ZUL(1,1) = std::complex<double>(pars(432), pars(433));
   ZUL(1,2) = std::complex<double>(pars(434), pars(435));
   ZUL(2,0) = std::complex<double>(pars(436), pars(437));
   ZUL(2,1) = std::complex<double>(pars(438), pars(439));
   ZUL(2,2) = std::complex<double>(pars(440), pars(441));
   ZUR(0,0) = std::complex<double>(pars(442), pars(443));
   ZUR(0,1) = std::complex<double>(pars(444), pars(445));
   ZUR(0,2) = std::complex<double>(pars(446), pars(447));
   ZUR(1,0) = std::complex<double>(pars(448), pars(449));
   ZUR(1,1) = std::complex<double>(pars(450), pars(451));
   ZUR(1,2) = std::complex<double>(pars(452), pars(453));
   ZUR(2,0) = std::complex<double>(pars(454), pars(455));
   ZUR(2,1) = std::complex<double>(pars(456), pars(457));
   ZUR(2,2) = std::complex<double>(pars(458), pars(459));
   ZDXL(0,0) = std::complex<double>(pars(460), pars(461));
   ZDXL(0,1) = std::complex<double>(pars(462), pars(463));
   ZDXL(0,2) = std::complex<double>(pars(464), pars(465));
   ZDXL(1,0) = std::complex<double>(pars(466), pars(467));
   ZDXL(1,1) = std::complex<double>(pars(468), pars(469));
   ZDXL(1,2) = std::complex<double>(pars(470), pars(471));
   ZDXL(2,0) = std::complex<double>(pars(472), pars(473));
   ZDXL(2,1) = std::complex<double>(pars(474), pars(475));
   ZDXL(2,2) = std::complex<double>(pars(476), pars(477));
   ZDXR(0,0) = std::complex<double>(pars(478), pars(479));
   ZDXR(0,1) = std::complex<double>(pars(480), pars(481));
   ZDXR(0,2) = std::complex<double>(pars(482), pars(483));
   ZDXR(1,0) = std::complex<double>(pars(484), pars(485));
   ZDXR(1,1) = std::complex<double>(pars(486), pars(487));
   ZDXR(1,2) = std::complex<double>(pars(488), pars(489));
   ZDXR(2,0) = std::complex<double>(pars(490), pars(491));
   ZDXR(2,1) = std::complex<double>(pars(492), pars(493));
   ZDXR(2,2) = std::complex<double>(pars(494), pars(495));
   UHI0(0,0) = pars(496);
   UHI0(0,1) = pars(497);
   UHI0(0,2) = pars(498);
   UHI0(0,3) = pars(499);
   UHI0(1,0) = pars(500);
   UHI0(1,1) = pars(501);
   UHI0(1,2) = pars(502);
   UHI0(1,3) = pars(503);
   UHI0(2,0) = pars(504);
   UHI0(2,1) = pars(505);
   UHI0(2,2) = pars(506);
   UHI0(2,3) = pars(507);
   UHI0(3,0) = pars(508);
   UHI0(3,1) = pars(509);
   UHI0(3,2) = pars(510);
   UHI0(3,3) = pars(511);
   UHIp(0,0) = pars(512);
   UHIp(0,1) = pars(513);
   UHIp(0,2) = pars(514);
   UHIp(0,3) = pars(515);
   UHIp(1,0) = pars(516);
   UHIp(1,1) = pars(517);
   UHIp(1,2) = pars(518);
   UHIp(1,3) = pars(519);
   UHIp(2,0) = pars(520);
   UHIp(2,1) = pars(521);
   UHIp(2,2) = pars(522);
   UHIp(2,3) = pars(523);
   UHIp(3,0) = pars(524);
   UHIp(3,1) = pars(525);
   UHIp(3,2) = pars(526);
   UHIp(3,3) = pars(527);
   ZMI(0,0) = std::complex<double>(pars(528), pars(529));
   ZMI(0,1) = std::complex<double>(pars(530), pars(531));
   ZMI(1,0) = std::complex<double>(pars(532), pars(533));
   ZMI(1,1) = std::complex<double>(pars(534), pars(535));
   ZPI(0,0) = std::complex<double>(pars(536), pars(537));
   ZPI(0,1) = std::complex<double>(pars(538), pars(539));
   ZPI(1,0) = std::complex<double>(pars(540), pars(541));
   ZPI(1,1) = std::complex<double>(pars(542), pars(543));
   ZNI(0,0) = std::complex<double>(pars(544), pars(545));
   ZNI(0,1) = std::complex<double>(pars(546), pars(547));
   ZNI(0,2) = std::complex<double>(pars(548), pars(549));
   ZNI(0,3) = std::complex<double>(pars(550), pars(551));
   ZNI(1,0) = std::complex<double>(pars(552), pars(553));
   ZNI(1,1) = std::complex<double>(pars(554), pars(555));
   ZNI(1,2) = std::complex<double>(pars(556), pars(557));
   ZNI(1,3) = std::complex<double>(pars(558), pars(559));
   ZNI(2,0) = std::complex<double>(pars(560), pars(561));
   ZNI(2,1) = std::complex<double>(pars(562), pars(563));
   ZNI(2,2) = std::complex<double>(pars(564), pars(565));
   ZNI(2,3) = std::complex<double>(pars(566), pars(567));
   ZNI(3,0) = std::complex<double>(pars(568), pars(569));
   ZNI(3,1) = std::complex<double>(pars(570), pars(571));
   ZNI(3,2) = std::complex<double>(pars(572), pars(573));
   ZNI(3,3) = std::complex<double>(pars(574), pars(575));
   ZSSI(0,0) = pars(576);
   ZSSI(0,1) = pars(577);
   ZSSI(1,0) = pars(578);
   ZSSI(1,1) = pars(579);
   ZFSI(0,0) = std::complex<double>(pars(580), pars(581));
   ZFSI(0,1) = std::complex<double>(pars(582), pars(583));
   ZFSI(1,0) = std::complex<double>(pars(584), pars(585));
   ZFSI(1,1) = std::complex<double>(pars(586), pars(587));
   UHp0(0,0) = pars(588);
   UHp0(0,1) = pars(589);
   UHp0(1,0) = pars(590);
   UHp0(1,1) = pars(591);
   UHpp(0,0) = pars(592);
   UHpp(0,1) = pars(593);
   UHpp(1,0) = pars(594);
   UHpp(1,1) = pars(595);
   ZNp(0,0) = std::complex<double>(pars(596), pars(597));
   ZNp(0,1) = std::complex<double>(pars(598), pars(599));
   ZNp(1,0) = std::complex<double>(pars(600), pars(601));
   ZNp(1,1) = std::complex<double>(pars(602), pars(603));
   ZZ(0,0) = pars(604);
   ZZ(0,1) = pars(605);
   ZZ(0,2) = pars(606);
   ZZ(1,0) = pars(607);
   ZZ(1,1) = pars(608);
   ZZ(1,2) = pars(609);
   ZZ(2,0) = pars(610);
   ZZ(2,1) = pars(611);
   ZZ(2,2) = pars(612);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(613);

   pars(89) = ZD(0,0);
   pars(90) = ZD(0,1);
   pars(91) = ZD(0,2);
   pars(92) = ZD(0,3);
   pars(93) = ZD(0,4);
   pars(94) = ZD(0,5);
   pars(95) = ZD(1,0);
   pars(96) = ZD(1,1);
   pars(97) = ZD(1,2);
   pars(98) = ZD(1,3);
   pars(99) = ZD(1,4);
   pars(100) = ZD(1,5);
   pars(101) = ZD(2,0);
   pars(102) = ZD(2,1);
   pars(103) = ZD(2,2);
   pars(104) = ZD(2,3);
   pars(105) = ZD(2,4);
   pars(106) = ZD(2,5);
   pars(107) = ZD(3,0);
   pars(108) = ZD(3,1);
   pars(109) = ZD(3,2);
   pars(110) = ZD(3,3);
   pars(111) = ZD(3,4);
   pars(112) = ZD(3,5);
   pars(113) = ZD(4,0);
   pars(114) = ZD(4,1);
   pars(115) = ZD(4,2);
   pars(116) = ZD(4,3);
   pars(117) = ZD(4,4);
   pars(118) = ZD(4,5);
   pars(119) = ZD(5,0);
   pars(120) = ZD(5,1);
   pars(121) = ZD(5,2);
   pars(122) = ZD(5,3);
   pars(123) = ZD(5,4);
   pars(124) = ZD(5,5);
   pars(125) = ZV(0,0);
   pars(126) = ZV(0,1);
   pars(127) = ZV(0,2);
   pars(128) = ZV(1,0);
   pars(129) = ZV(1,1);
   pars(130) = ZV(1,2);
   pars(131) = ZV(2,0);
   pars(132) = ZV(2,1);
   pars(133) = ZV(2,2);
   pars(134) = ZU(0,0);
   pars(135) = ZU(0,1);
   pars(136) = ZU(0,2);
   pars(137) = ZU(0,3);
   pars(138) = ZU(0,4);
   pars(139) = ZU(0,5);
   pars(140) = ZU(1,0);
   pars(141) = ZU(1,1);
   pars(142) = ZU(1,2);
   pars(143) = ZU(1,3);
   pars(144) = ZU(1,4);
   pars(145) = ZU(1,5);
   pars(146) = ZU(2,0);
   pars(147) = ZU(2,1);
   pars(148) = ZU(2,2);
   pars(149) = ZU(2,3);
   pars(150) = ZU(2,4);
   pars(151) = ZU(2,5);
   pars(152) = ZU(3,0);
   pars(153) = ZU(3,1);
   pars(154) = ZU(3,2);
   pars(155) = ZU(3,3);
   pars(156) = ZU(3,4);
   pars(157) = ZU(3,5);
   pars(158) = ZU(4,0);
   pars(159) = ZU(4,1);
   pars(160) = ZU(4,2);
   pars(161) = ZU(4,3);
   pars(162) = ZU(4,4);
   pars(163) = ZU(4,5);
   pars(164) = ZU(5,0);
   pars(165) = ZU(5,1);
   pars(166) = ZU(5,2);
   pars(167) = ZU(5,3);
   pars(168) = ZU(5,4);
   pars(169) = ZU(5,5);
   pars(170) = ZE(0,0);
   pars(171) = ZE(0,1);
   pars(172) = ZE(0,2);
   pars(173) = ZE(0,3);
   pars(174) = ZE(0,4);
   pars(175) = ZE(0,5);
   pars(176) = ZE(1,0);
   pars(177) = ZE(1,1);
   pars(178) = ZE(1,2);
   pars(179) = ZE(1,3);
   pars(180) = ZE(1,4);
   pars(181) = ZE(1,5);
   pars(182) = ZE(2,0);
   pars(183) = ZE(2,1);
   pars(184) = ZE(2,2);
   pars(185) = ZE(2,3);
   pars(186) = ZE(2,4);
   pars(187) = ZE(2,5);
   pars(188) = ZE(3,0);
   pars(189) = ZE(3,1);
   pars(190) = ZE(3,2);
   pars(191) = ZE(3,3);
   pars(192) = ZE(3,4);
   pars(193) = ZE(3,5);
   pars(194) = ZE(4,0);
   pars(195) = ZE(4,1);
   pars(196) = ZE(4,2);
   pars(197) = ZE(4,3);
   pars(198) = ZE(4,4);
   pars(199) = ZE(4,5);
   pars(200) = ZE(5,0);
   pars(201) = ZE(5,1);
   pars(202) = ZE(5,2);
   pars(203) = ZE(5,3);
   pars(204) = ZE(5,4);
   pars(205) = ZE(5,5);
   pars(206) = ZDX(0,0);
   pars(207) = ZDX(0,1);
   pars(208) = ZDX(0,2);
   pars(209) = ZDX(0,3);
   pars(210) = ZDX(0,4);
   pars(211) = ZDX(0,5);
   pars(212) = ZDX(1,0);
   pars(213) = ZDX(1,1);
   pars(214) = ZDX(1,2);
   pars(215) = ZDX(1,3);
   pars(216) = ZDX(1,4);
   pars(217) = ZDX(1,5);
   pars(218) = ZDX(2,0);
   pars(219) = ZDX(2,1);
   pars(220) = ZDX(2,2);
   pars(221) = ZDX(2,3);
   pars(222) = ZDX(2,4);
   pars(223) = ZDX(2,5);
   pars(224) = ZDX(3,0);
   pars(225) = ZDX(3,1);
   pars(226) = ZDX(3,2);
   pars(227) = ZDX(3,3);
   pars(228) = ZDX(3,4);
   pars(229) = ZDX(3,5);
   pars(230) = ZDX(4,0);
   pars(231) = ZDX(4,1);
   pars(232) = ZDX(4,2);
   pars(233) = ZDX(4,3);
   pars(234) = ZDX(4,4);
   pars(235) = ZDX(4,5);
   pars(236) = ZDX(5,0);
   pars(237) = ZDX(5,1);
   pars(238) = ZDX(5,2);
   pars(239) = ZDX(5,3);
   pars(240) = ZDX(5,4);
   pars(241) = ZDX(5,5);
   pars(242) = ZH(0,0);
   pars(243) = ZH(0,1);
   pars(244) = ZH(0,2);
   pars(245) = ZH(1,0);
   pars(246) = ZH(1,1);
   pars(247) = ZH(1,2);
   pars(248) = ZH(2,0);
   pars(249) = ZH(2,1);
   pars(250) = ZH(2,2);
   pars(251) = ZA(0,0);
   pars(252) = ZA(0,1);
   pars(253) = ZA(0,2);
   pars(254) = ZA(1,0);
   pars(255) = ZA(1,1);
   pars(256) = ZA(1,2);
   pars(257) = ZA(2,0);
   pars(258) = ZA(2,1);
   pars(259) = ZA(2,2);
   pars(260) = ZP(0,0);
   pars(261) = ZP(0,1);
   pars(262) = ZP(1,0);
   pars(263) = ZP(1,1);
   pars(264) = Re(ZN(0,0));
   pars(265) = Im(ZN(0,0));
   pars(266) = Re(ZN(0,1));
   pars(267) = Im(ZN(0,1));
   pars(268) = Re(ZN(0,2));
   pars(269) = Im(ZN(0,2));
   pars(270) = Re(ZN(0,3));
   pars(271) = Im(ZN(0,3));
   pars(272) = Re(ZN(0,4));
   pars(273) = Im(ZN(0,4));
   pars(274) = Re(ZN(0,5));
   pars(275) = Im(ZN(0,5));
   pars(276) = Re(ZN(1,0));
   pars(277) = Im(ZN(1,0));
   pars(278) = Re(ZN(1,1));
   pars(279) = Im(ZN(1,1));
   pars(280) = Re(ZN(1,2));
   pars(281) = Im(ZN(1,2));
   pars(282) = Re(ZN(1,3));
   pars(283) = Im(ZN(1,3));
   pars(284) = Re(ZN(1,4));
   pars(285) = Im(ZN(1,4));
   pars(286) = Re(ZN(1,5));
   pars(287) = Im(ZN(1,5));
   pars(288) = Re(ZN(2,0));
   pars(289) = Im(ZN(2,0));
   pars(290) = Re(ZN(2,1));
   pars(291) = Im(ZN(2,1));
   pars(292) = Re(ZN(2,2));
   pars(293) = Im(ZN(2,2));
   pars(294) = Re(ZN(2,3));
   pars(295) = Im(ZN(2,3));
   pars(296) = Re(ZN(2,4));
   pars(297) = Im(ZN(2,4));
   pars(298) = Re(ZN(2,5));
   pars(299) = Im(ZN(2,5));
   pars(300) = Re(ZN(3,0));
   pars(301) = Im(ZN(3,0));
   pars(302) = Re(ZN(3,1));
   pars(303) = Im(ZN(3,1));
   pars(304) = Re(ZN(3,2));
   pars(305) = Im(ZN(3,2));
   pars(306) = Re(ZN(3,3));
   pars(307) = Im(ZN(3,3));
   pars(308) = Re(ZN(3,4));
   pars(309) = Im(ZN(3,4));
   pars(310) = Re(ZN(3,5));
   pars(311) = Im(ZN(3,5));
   pars(312) = Re(ZN(4,0));
   pars(313) = Im(ZN(4,0));
   pars(314) = Re(ZN(4,1));
   pars(315) = Im(ZN(4,1));
   pars(316) = Re(ZN(4,2));
   pars(317) = Im(ZN(4,2));
   pars(318) = Re(ZN(4,3));
   pars(319) = Im(ZN(4,3));
   pars(320) = Re(ZN(4,4));
   pars(321) = Im(ZN(4,4));
   pars(322) = Re(ZN(4,5));
   pars(323) = Im(ZN(4,5));
   pars(324) = Re(ZN(5,0));
   pars(325) = Im(ZN(5,0));
   pars(326) = Re(ZN(5,1));
   pars(327) = Im(ZN(5,1));
   pars(328) = Re(ZN(5,2));
   pars(329) = Im(ZN(5,2));
   pars(330) = Re(ZN(5,3));
   pars(331) = Im(ZN(5,3));
   pars(332) = Re(ZN(5,4));
   pars(333) = Im(ZN(5,4));
   pars(334) = Re(ZN(5,5));
   pars(335) = Im(ZN(5,5));
   pars(336) = Re(UM(0,0));
   pars(337) = Im(UM(0,0));
   pars(338) = Re(UM(0,1));
   pars(339) = Im(UM(0,1));
   pars(340) = Re(UM(1,0));
   pars(341) = Im(UM(1,0));
   pars(342) = Re(UM(1,1));
   pars(343) = Im(UM(1,1));
   pars(344) = Re(UP(0,0));
   pars(345) = Im(UP(0,0));
   pars(346) = Re(UP(0,1));
   pars(347) = Im(UP(0,1));
   pars(348) = Re(UP(1,0));
   pars(349) = Im(UP(1,0));
   pars(350) = Re(UP(1,1));
   pars(351) = Im(UP(1,1));
   pars(352) = Re(ZEL(0,0));
   pars(353) = Im(ZEL(0,0));
   pars(354) = Re(ZEL(0,1));
   pars(355) = Im(ZEL(0,1));
   pars(356) = Re(ZEL(0,2));
   pars(357) = Im(ZEL(0,2));
   pars(358) = Re(ZEL(1,0));
   pars(359) = Im(ZEL(1,0));
   pars(360) = Re(ZEL(1,1));
   pars(361) = Im(ZEL(1,1));
   pars(362) = Re(ZEL(1,2));
   pars(363) = Im(ZEL(1,2));
   pars(364) = Re(ZEL(2,0));
   pars(365) = Im(ZEL(2,0));
   pars(366) = Re(ZEL(2,1));
   pars(367) = Im(ZEL(2,1));
   pars(368) = Re(ZEL(2,2));
   pars(369) = Im(ZEL(2,2));
   pars(370) = Re(ZER(0,0));
   pars(371) = Im(ZER(0,0));
   pars(372) = Re(ZER(0,1));
   pars(373) = Im(ZER(0,1));
   pars(374) = Re(ZER(0,2));
   pars(375) = Im(ZER(0,2));
   pars(376) = Re(ZER(1,0));
   pars(377) = Im(ZER(1,0));
   pars(378) = Re(ZER(1,1));
   pars(379) = Im(ZER(1,1));
   pars(380) = Re(ZER(1,2));
   pars(381) = Im(ZER(1,2));
   pars(382) = Re(ZER(2,0));
   pars(383) = Im(ZER(2,0));
   pars(384) = Re(ZER(2,1));
   pars(385) = Im(ZER(2,1));
   pars(386) = Re(ZER(2,2));
   pars(387) = Im(ZER(2,2));
   pars(388) = Re(ZDL(0,0));
   pars(389) = Im(ZDL(0,0));
   pars(390) = Re(ZDL(0,1));
   pars(391) = Im(ZDL(0,1));
   pars(392) = Re(ZDL(0,2));
   pars(393) = Im(ZDL(0,2));
   pars(394) = Re(ZDL(1,0));
   pars(395) = Im(ZDL(1,0));
   pars(396) = Re(ZDL(1,1));
   pars(397) = Im(ZDL(1,1));
   pars(398) = Re(ZDL(1,2));
   pars(399) = Im(ZDL(1,2));
   pars(400) = Re(ZDL(2,0));
   pars(401) = Im(ZDL(2,0));
   pars(402) = Re(ZDL(2,1));
   pars(403) = Im(ZDL(2,1));
   pars(404) = Re(ZDL(2,2));
   pars(405) = Im(ZDL(2,2));
   pars(406) = Re(ZDR(0,0));
   pars(407) = Im(ZDR(0,0));
   pars(408) = Re(ZDR(0,1));
   pars(409) = Im(ZDR(0,1));
   pars(410) = Re(ZDR(0,2));
   pars(411) = Im(ZDR(0,2));
   pars(412) = Re(ZDR(1,0));
   pars(413) = Im(ZDR(1,0));
   pars(414) = Re(ZDR(1,1));
   pars(415) = Im(ZDR(1,1));
   pars(416) = Re(ZDR(1,2));
   pars(417) = Im(ZDR(1,2));
   pars(418) = Re(ZDR(2,0));
   pars(419) = Im(ZDR(2,0));
   pars(420) = Re(ZDR(2,1));
   pars(421) = Im(ZDR(2,1));
   pars(422) = Re(ZDR(2,2));
   pars(423) = Im(ZDR(2,2));
   pars(424) = Re(ZUL(0,0));
   pars(425) = Im(ZUL(0,0));
   pars(426) = Re(ZUL(0,1));
   pars(427) = Im(ZUL(0,1));
   pars(428) = Re(ZUL(0,2));
   pars(429) = Im(ZUL(0,2));
   pars(430) = Re(ZUL(1,0));
   pars(431) = Im(ZUL(1,0));
   pars(432) = Re(ZUL(1,1));
   pars(433) = Im(ZUL(1,1));
   pars(434) = Re(ZUL(1,2));
   pars(435) = Im(ZUL(1,2));
   pars(436) = Re(ZUL(2,0));
   pars(437) = Im(ZUL(2,0));
   pars(438) = Re(ZUL(2,1));
   pars(439) = Im(ZUL(2,1));
   pars(440) = Re(ZUL(2,2));
   pars(441) = Im(ZUL(2,2));
   pars(442) = Re(ZUR(0,0));
   pars(443) = Im(ZUR(0,0));
   pars(444) = Re(ZUR(0,1));
   pars(445) = Im(ZUR(0,1));
   pars(446) = Re(ZUR(0,2));
   pars(447) = Im(ZUR(0,2));
   pars(448) = Re(ZUR(1,0));
   pars(449) = Im(ZUR(1,0));
   pars(450) = Re(ZUR(1,1));
   pars(451) = Im(ZUR(1,1));
   pars(452) = Re(ZUR(1,2));
   pars(453) = Im(ZUR(1,2));
   pars(454) = Re(ZUR(2,0));
   pars(455) = Im(ZUR(2,0));
   pars(456) = Re(ZUR(2,1));
   pars(457) = Im(ZUR(2,1));
   pars(458) = Re(ZUR(2,2));
   pars(459) = Im(ZUR(2,2));
   pars(460) = Re(ZDXL(0,0));
   pars(461) = Im(ZDXL(0,0));
   pars(462) = Re(ZDXL(0,1));
   pars(463) = Im(ZDXL(0,1));
   pars(464) = Re(ZDXL(0,2));
   pars(465) = Im(ZDXL(0,2));
   pars(466) = Re(ZDXL(1,0));
   pars(467) = Im(ZDXL(1,0));
   pars(468) = Re(ZDXL(1,1));
   pars(469) = Im(ZDXL(1,1));
   pars(470) = Re(ZDXL(1,2));
   pars(471) = Im(ZDXL(1,2));
   pars(472) = Re(ZDXL(2,0));
   pars(473) = Im(ZDXL(2,0));
   pars(474) = Re(ZDXL(2,1));
   pars(475) = Im(ZDXL(2,1));
   pars(476) = Re(ZDXL(2,2));
   pars(477) = Im(ZDXL(2,2));
   pars(478) = Re(ZDXR(0,0));
   pars(479) = Im(ZDXR(0,0));
   pars(480) = Re(ZDXR(0,1));
   pars(481) = Im(ZDXR(0,1));
   pars(482) = Re(ZDXR(0,2));
   pars(483) = Im(ZDXR(0,2));
   pars(484) = Re(ZDXR(1,0));
   pars(485) = Im(ZDXR(1,0));
   pars(486) = Re(ZDXR(1,1));
   pars(487) = Im(ZDXR(1,1));
   pars(488) = Re(ZDXR(1,2));
   pars(489) = Im(ZDXR(1,2));
   pars(490) = Re(ZDXR(2,0));
   pars(491) = Im(ZDXR(2,0));
   pars(492) = Re(ZDXR(2,1));
   pars(493) = Im(ZDXR(2,1));
   pars(494) = Re(ZDXR(2,2));
   pars(495) = Im(ZDXR(2,2));
   pars(496) = UHI0(0,0);
   pars(497) = UHI0(0,1);
   pars(498) = UHI0(0,2);
   pars(499) = UHI0(0,3);
   pars(500) = UHI0(1,0);
   pars(501) = UHI0(1,1);
   pars(502) = UHI0(1,2);
   pars(503) = UHI0(1,3);
   pars(504) = UHI0(2,0);
   pars(505) = UHI0(2,1);
   pars(506) = UHI0(2,2);
   pars(507) = UHI0(2,3);
   pars(508) = UHI0(3,0);
   pars(509) = UHI0(3,1);
   pars(510) = UHI0(3,2);
   pars(511) = UHI0(3,3);
   pars(512) = UHIp(0,0);
   pars(513) = UHIp(0,1);
   pars(514) = UHIp(0,2);
   pars(515) = UHIp(0,3);
   pars(516) = UHIp(1,0);
   pars(517) = UHIp(1,1);
   pars(518) = UHIp(1,2);
   pars(519) = UHIp(1,3);
   pars(520) = UHIp(2,0);
   pars(521) = UHIp(2,1);
   pars(522) = UHIp(2,2);
   pars(523) = UHIp(2,3);
   pars(524) = UHIp(3,0);
   pars(525) = UHIp(3,1);
   pars(526) = UHIp(3,2);
   pars(527) = UHIp(3,3);
   pars(528) = Re(ZMI(0,0));
   pars(529) = Im(ZMI(0,0));
   pars(530) = Re(ZMI(0,1));
   pars(531) = Im(ZMI(0,1));
   pars(532) = Re(ZMI(1,0));
   pars(533) = Im(ZMI(1,0));
   pars(534) = Re(ZMI(1,1));
   pars(535) = Im(ZMI(1,1));
   pars(536) = Re(ZPI(0,0));
   pars(537) = Im(ZPI(0,0));
   pars(538) = Re(ZPI(0,1));
   pars(539) = Im(ZPI(0,1));
   pars(540) = Re(ZPI(1,0));
   pars(541) = Im(ZPI(1,0));
   pars(542) = Re(ZPI(1,1));
   pars(543) = Im(ZPI(1,1));
   pars(544) = Re(ZNI(0,0));
   pars(545) = Im(ZNI(0,0));
   pars(546) = Re(ZNI(0,1));
   pars(547) = Im(ZNI(0,1));
   pars(548) = Re(ZNI(0,2));
   pars(549) = Im(ZNI(0,2));
   pars(550) = Re(ZNI(0,3));
   pars(551) = Im(ZNI(0,3));
   pars(552) = Re(ZNI(1,0));
   pars(553) = Im(ZNI(1,0));
   pars(554) = Re(ZNI(1,1));
   pars(555) = Im(ZNI(1,1));
   pars(556) = Re(ZNI(1,2));
   pars(557) = Im(ZNI(1,2));
   pars(558) = Re(ZNI(1,3));
   pars(559) = Im(ZNI(1,3));
   pars(560) = Re(ZNI(2,0));
   pars(561) = Im(ZNI(2,0));
   pars(562) = Re(ZNI(2,1));
   pars(563) = Im(ZNI(2,1));
   pars(564) = Re(ZNI(2,2));
   pars(565) = Im(ZNI(2,2));
   pars(566) = Re(ZNI(2,3));
   pars(567) = Im(ZNI(2,3));
   pars(568) = Re(ZNI(3,0));
   pars(569) = Im(ZNI(3,0));
   pars(570) = Re(ZNI(3,1));
   pars(571) = Im(ZNI(3,1));
   pars(572) = Re(ZNI(3,2));
   pars(573) = Im(ZNI(3,2));
   pars(574) = Re(ZNI(3,3));
   pars(575) = Im(ZNI(3,3));
   pars(576) = ZSSI(0,0);
   pars(577) = ZSSI(0,1);
   pars(578) = ZSSI(1,0);
   pars(579) = ZSSI(1,1);
   pars(580) = Re(ZFSI(0,0));
   pars(581) = Im(ZFSI(0,0));
   pars(582) = Re(ZFSI(0,1));
   pars(583) = Im(ZFSI(0,1));
   pars(584) = Re(ZFSI(1,0));
   pars(585) = Im(ZFSI(1,0));
   pars(586) = Re(ZFSI(1,1));
   pars(587) = Im(ZFSI(1,1));
   pars(588) = UHp0(0,0);
   pars(589) = UHp0(0,1);
   pars(590) = UHp0(1,0);
   pars(591) = UHp0(1,1);
   pars(592) = UHpp(0,0);
   pars(593) = UHpp(0,1);
   pars(594) = UHpp(1,0);
   pars(595) = UHpp(1,1);
   pars(596) = Re(ZNp(0,0));
   pars(597) = Im(ZNp(0,0));
   pars(598) = Re(ZNp(0,1));
   pars(599) = Im(ZNp(0,1));
   pars(600) = Re(ZNp(1,0));
   pars(601) = Im(ZNp(1,0));
   pars(602) = Re(ZNp(1,1));
   pars(603) = Im(ZNp(1,1));
   pars(604) = ZZ(0,0);
   pars(605) = ZZ(0,1);
   pars(606) = ZZ(0,2);
   pars(607) = ZZ(1,0);
   pars(608) = ZZ(1,1);
   pars(609) = ZZ(1,2);
   pars(610) = ZZ(2,0);
   pars(611) = ZZ(2,1);
   pars(612) = ZZ(2,2);


   return pars;
}

void CLASSNAME::set_extra_parameters(const Eigen::ArrayXd& pars)
{

}

Eigen::ArrayXd CLASSNAME::get_extra_parameters() const
{
   return Eigen::ArrayXd();

}

std::string CLASSNAME::name() const
{
   return "E6SSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   E6SSM_soft_parameters::run_to(scale, eps);
}

void CLASSNAME::calculate_tree_level_mass_spectrum()
{
   calculate_DRbar_masses();
}

void CLASSNAME::calculate_pole_mass_spectrum()
{
   calculate_pole_masses();
}

void CLASSNAME::calculate_mass_spectrum()
{
   calculate_spectrum();
}

int CLASSNAME::solve_ewsb_equations_tree_level()
{
   return solve_ewsb_tree_level();
}

int CLASSNAME::solve_ewsb_equations()
{
   return solve_ewsb();
}

const E6SSM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

E6SSM_input_parameters& CLASSNAME::get_input_parameters()
{
   return get_input();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses() const
{
   return get_DRbar_masses();
}

Eigen::ArrayXd CLASSNAME::get_tree_level_masses_and_mixings() const
{
   return get_DRbar_masses_and_mixings();
}

void CLASSNAME::set_tree_level_masses(const Eigen::ArrayXd& m)
{
   set_DRbar_masses(m);
}

void CLASSNAME::set_tree_level_masses_and_mixings(const Eigen::ArrayXd& m)
{
   set_DRbar_masses_and_mixings(m);
}


Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_goldstone;
   MHpm_goldstone(0) = MVWm;

   return remove_if_equal(MHpm, MHpm_goldstone);
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,2,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;
   MAh_goldstone(1) = MVZp;

   return remove_if_equal(MAh, MAh_goldstone);
}


/**
 * @brief finds the LSP and returns it's mass
 *
 * This function finds the lightest supersymmetric particle (LSP) and
 * returns it's mass.  The corresponding particle type is retured in
 * the reference parameter.  The list of potential LSPs is set in the
 * model file varible PotentialLSPParticles.  For this model it is set
 * to:
 * {Chi, Sv, Su, Sd, Se, Cha, Glu}
 *
 * @param particle_type particle type
 * @return mass of LSP
 */
double CLASSNAME::get_lsp(E6SSM_info::Particles& particle_type) const
{
   double lsp_mass = std::numeric_limits<double>::max();
   double tmp_mass;
   particle_type = E6SSM_info::NUMBER_OF_PARTICLES;

   tmp_mass = Abs(PHYSICAL(MChi(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Chi;
   }

   tmp_mass = Abs(PHYSICAL(MSv(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Sv;
   }

   tmp_mass = Abs(PHYSICAL(MSu(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Su;
   }

   tmp_mass = Abs(PHYSICAL(MSd(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Sd;
   }

   tmp_mass = Abs(PHYSICAL(MSe(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Se;
   }

   tmp_mass = Abs(PHYSICAL(MCha(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Cha;
   }

   tmp_mass = Abs(PHYSICAL(MGlu));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = E6SSM_info::Glu;
   }

   return lsp_mass;
}


double CLASSNAME::get_mass_matrix_VG() const
{

   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{

   const auto mass_matrix_VG = get_mass_matrix_VG();
   MVG = mass_matrix_VG;
}

double CLASSNAME::get_mass_matrix_Glu() const
{

   const double mass_matrix_Glu = Re(MassG);

   return mass_matrix_Glu;
}

void CLASSNAME::calculate_MGlu()
{

   const auto mass_matrix_Glu = get_mass_matrix_Glu();
   MGlu = calculate_majorana_singlet_mass(mass_matrix_Glu, PhaseGlu);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,2) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{

   MFv.setConstant(0);
}

double CLASSNAME::get_mass_matrix_ChaP() const
{

   const double mass_matrix_ChaP = Re(MuPr);

   return mass_matrix_ChaP;
}

void CLASSNAME::calculate_MChaP()
{

   const auto mass_matrix_ChaP = get_mass_matrix_ChaP();
   MChaP = calculate_dirac_singlet_mass(mass_matrix_ChaP, PhaseFHpup);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*AbsSqr(Yd(0,0))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN
      )*Sqr(vs) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN)
      *Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1);
   mass_matrix_Sd(0,2) = mq2(0,2);
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) - 0.5*vs*vu*Conj(
      Yd(0,0))*Lambdax;
   mass_matrix_Sd(0,4) = 0;
   mass_matrix_Sd(0,5) = 0;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*AbsSqr(Yd(1,1))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN
      )*Sqr(vs) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN)
      *Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2);
   mass_matrix_Sd(1,3) = 0;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) - 0.5*vs*vu*Conj(
      Yd(1,1))*Lambdax;
   mass_matrix_Sd(1,5) = 0;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*AbsSqr(Yd(2,2))*Sqr(vd) - 0.025*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN
      )*Sqr(vs) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN)
      *Sqr(vu);
   mass_matrix_Sd(2,3) = 0;
   mass_matrix_Sd(2,4) = 0;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vs*vu*Conj(
      Yd(2,2))*Lambdax;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*AbsSqr(Yd(0,0))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) + 0.05*Sqr(g1)*
      Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1);
   mass_matrix_Sd(3,5) = md2(0,2);
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*AbsSqr(Yd(1,1))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) + 0.05*Sqr(g1)*
      Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2);
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*AbsSqr(Yd(2,2))*Sqr(vd) - 0.05*Sqr(g1)*
      Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) + 0.05*Sqr(g1)*
      Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Sd, eigenvalue_error > precision * Abs(
      MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Sv, eigenvalue_error > precision * Abs(
      MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(0,0
      ))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN
      )*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1);
   mass_matrix_Su(0,2) = mq2(0,2);
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) - 0.5*vd*vs*Conj(
      Yu(0,0))*Lambdax;
   mass_matrix_Su(0,4) = 0;
   mass_matrix_Su(0,5) = 0;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(1,1
      ))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN
      )*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2);
   mass_matrix_Su(1,3) = 0;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vs*Conj(
      Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(2,2
      ))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.025*Sqr(gN
      )*Sqr(vu);
   mass_matrix_Su(2,3) = 0;
   mass_matrix_Su(2,4) = 0;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vs*Conj(
      Yu(2,2))*Lambdax;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd
      ) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(0,0))*Sqr(vu) - 0.1*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1);
   mass_matrix_Su(3,5) = mu2(0,2);
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd
      ) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(1,1))*Sqr(vu) - 0.1*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2);
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd
      ) + 0.0625*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Yu(2,2))*Sqr(vu) - 0.1*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Su, eigenvalue_error > precision * Abs(
      MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*AbsSqr(Ye(0,0))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1);
   mass_matrix_Se(0,2) = ml2(0,2);
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) - 0.5*vs*vu*Conj(
      Ye(0,0))*Lambdax;
   mass_matrix_Se(0,4) = 0;
   mass_matrix_Se(0,5) = 0;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*AbsSqr(Ye(1,1))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2);
   mass_matrix_Se(1,3) = 0;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) - 0.5*vs*vu*Conj(
      Ye(1,1))*Lambdax;
   mass_matrix_Se(1,5) = 0;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*AbsSqr(Ye(2,2))*Sqr(vd) + 0.075*Sqr(g1)
      *Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Se(2,3) = 0;
   mass_matrix_Se(2,4) = 0;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vs*vu*Conj(
      Ye(2,2))*Lambdax;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*AbsSqr(Ye(0,0))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.15*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1);
   mass_matrix_Se(3,5) = me2(0,2);
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*AbsSqr(Ye(1,1))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.15*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2);
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*AbsSqr(Ye(2,2))*Sqr(vd) - 0.15*Sqr(g1)*
      Sqr(vd) - 0.0375*Sqr(gN)*Sqr(vd) + 0.0625*Sqr(gN)*Sqr(vs) + 0.15*Sqr(g1)*
      Sqr(vu) - 0.025*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Se, eigenvalue_error > precision * Abs(
      MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_SDX() const
{

   Eigen::Matrix<double,6,6> mass_matrix_SDX;

   mass_matrix_SDX(0,0) = mDx2(0,0) + 0.05*Sqr(g1)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(
      vd) + 0.5*AbsSqr(Kappa(0,0))*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) - 0.05*Sqr(
      g1)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(0,1) = 0;
   mass_matrix_SDX(0,2) = 0;
   mass_matrix_SDX(0,3) = 0.7071067811865475*vs*Conj(TKappa(0,0)) - 0.5*vd*vu*
      Conj(Kappa(0,0))*Lambdax;
   mass_matrix_SDX(0,4) = 0;
   mass_matrix_SDX(0,5) = 0;
   mass_matrix_SDX(1,1) = mDx2(1,1) + 0.05*Sqr(g1)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(
      vd) + 0.5*AbsSqr(Kappa(1,1))*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) - 0.05*Sqr(
      g1)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(1,2) = 0;
   mass_matrix_SDX(1,3) = 0;
   mass_matrix_SDX(1,4) = 0.7071067811865475*vs*Conj(TKappa(1,1)) - 0.5*vd*vu*
      Conj(Kappa(1,1))*Lambdax;
   mass_matrix_SDX(1,5) = 0;
   mass_matrix_SDX(2,2) = mDx2(2,2) + 0.05*Sqr(g1)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(
      vd) + 0.5*AbsSqr(Kappa(2,2))*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) - 0.05*Sqr(
      g1)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(2,3) = 0;
   mass_matrix_SDX(2,4) = 0;
   mass_matrix_SDX(2,5) = 0.7071067811865475*vs*Conj(TKappa(2,2)) - 0.5*vd*vu*
      Conj(Kappa(2,2))*Lambdax;
   mass_matrix_SDX(3,3) = mDxbar2(0,0) - 0.05*Sqr(g1)*Sqr(vd) + 0.1125*Sqr(gN)*
      Sqr(vd) + 0.5*AbsSqr(Kappa(0,0))*Sqr(vs) - 0.1875*Sqr(gN)*Sqr(vs) + 0.05*
      Sqr(g1)*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(3,4) = 0;
   mass_matrix_SDX(3,5) = 0;
   mass_matrix_SDX(4,4) = mDxbar2(1,1) - 0.05*Sqr(g1)*Sqr(vd) + 0.1125*Sqr(gN)*
      Sqr(vd) + 0.5*AbsSqr(Kappa(1,1))*Sqr(vs) - 0.1875*Sqr(gN)*Sqr(vs) + 0.05*
      Sqr(g1)*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SDX(4,5) = 0;
   mass_matrix_SDX(5,5) = mDxbar2(2,2) - 0.05*Sqr(g1)*Sqr(vd) + 0.1125*Sqr(gN)*
      Sqr(vd) + 0.5*AbsSqr(Kappa(2,2))*Sqr(vs) - 0.1875*Sqr(gN)*Sqr(vs) + 0.05*
      Sqr(g1)*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SDX);

   return mass_matrix_SDX;
}

void CLASSNAME::calculate_MSDX()
{
   const auto mass_matrix_SDX(get_mass_matrix_SDX());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SDX, MSDX, ZDX, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::SDX, eigenvalue_error > precision * Abs(
      MSDX(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SDX, MSDX, ZDX);
#endif
   normalize_to_interval(ZDX);


   if (MSDX.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::SDX);
   }

   MSDX = AbsSqrt(MSDX);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      0.3375*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.1875*Sqr(gN)*Sqr
      (vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2
      )*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vs*Conj(
      TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) + 0.15*vd*vu*Sqr(gN)
      - 0.35355339059327373*vs*TLambdax;
   mass_matrix_hh(0,2) = vd*vs*AbsSqr(Lambdax) - 0.35355339059327373*vu*Conj(
      TLambdax) - 0.375*vd*vs*Sqr(gN) - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax
      )*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) + 0.225*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)
      *Sqr(vu) + 0.15*Sqr(gN)*Sqr(vu);
   mass_matrix_hh(1,2) = vs*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*Conj(
      TLambdax) - 0.25*vs*vu*Sqr(gN) - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.1875*Sqr(gN)*Sqr
      (vd) + 0.9375*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.125*Sqr(
      gN)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::hh, eigenvalue_error > precision * Abs(
      Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) +
      0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.1875*Sqr(gN)*Sqr
      (vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2
      )*Sqr(vu) + 0.075*Sqr(gN)*Sqr(vu) + 0.3872983346207417*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())*Sqr(vd)*Sqr(Cos(ThetaWp())) + 0.225*Sqr(gN)*Sqr(vd)*Sqr(
      Cos(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp
      ())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vd)*Sqr(Sin(
      ThetaWp())) + 0.225*Sqr(gN)*Sqr(vd)*Sqr(Sin(ThetaWp())) + 0.25*Sqr(g2)*
      Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr
      (Sin(ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vs*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW())*Sqr(Cos(
      ThetaWp())) + 0.15*vd*vu*Sqr(gN)*Sqr(Cos(ThetaWp())) - 0.25*vd*vu*Sqr(g2)
      *Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 0.15*vd*vu*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) - 0.3872983346207417*g1*g2*vd*vu*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 0.15*vd*vu*Sqr(gN)*Sqr(Sin(
      ThetaWp())) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) -
      0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vs*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) - 0.375*vd*vs*
      Sqr(gN)*Sqr(Cos(ThetaWp())) - 0.375*vd*vs*Sqr(gN)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax
      )*Sqr(vs) - 0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)
      *Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu) + 0.3872983346207417*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())*Sqr(vu)*Sqr(Cos(ThetaWp())) + 0.1*Sqr(gN)*Sqr(vu)*Sqr(Cos(
      ThetaWp())) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu)*Sqr(Sin(
      ThetaWp())) + 0.1*Sqr(gN)*Sqr(vu)*Sqr(Sin(ThetaWp())) + 0.25*Sqr(g2)*Sqr(
      vu)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin
      (ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) - 0.25*vs*vu*Sqr
      (gN)*Sqr(Cos(ThetaWp())) - 0.25*vs*vu*Sqr(gN)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vd*TLambdax;
   mass_matrix_Ah(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.1875*Sqr(gN)*Sqr
      (vd) + 0.3125*Sqr(gN)*Sqr(vs) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.125*Sqr(
      gN)*Sqr(vu) + 0.625*Sqr(gN)*Sqr(vs)*Sqr(Cos(ThetaWp())) + 0.625*Sqr(gN)*
      Sqr(vs)*Sqr(Sin(ThetaWp()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Ah, eigenvalue_error > precision * Abs(
      MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd)
      + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.1875*Sqr(gN)*
      Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.075*Sqr(gN)*
      Sqr(vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475*vs*
      Conj(TLambdax);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd)
      + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vs) - 0.125*Sqr(gN)*Sqr
      (vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(
      vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Hpm, eigenvalue_error > precision * Abs(
      MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Chi() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(0,4) = 0;
   mass_matrix_Chi(0,5) = 0;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(1,4) = 0;
   mass_matrix_Chi(1,5) = 0;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -0.7071067811865475*vs*Lambdax;
   mass_matrix_Chi(2,4) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Chi(2,5) = -0.4743416490252569*gN*vd;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.7071067811865475*vd*Lambdax;
   mass_matrix_Chi(3,5) = -0.31622776601683794*gN*vu;
   mass_matrix_Chi(4,4) = 0;
   mass_matrix_Chi(4,5) = 0.7905694150420949*gN*vs;
   mass_matrix_Chi(5,5) = MassBp;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Chi, eigenvalue_error > precision * Abs(
      MChi(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
   normalize_to_interval(ZN);

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = 0.7071067811865475*vs*Lambdax;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Cha, eigenvalue_error > precision * Abs(
      MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*vd*Ye(0,0);
   mass_matrix_Fe(0,1) = 0;
   mass_matrix_Fe(0,2) = 0;
   mass_matrix_Fe(1,1) = 0.7071067811865475*vd*Ye(1,1);
   mass_matrix_Fe(1,2) = 0;
   mass_matrix_Fe(2,2) = 0.7071067811865475*vd*Ye(2,2);

   Symmetrize(mass_matrix_Fe);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Fe, eigenvalue_error > precision * Abs(
      MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*vd*Yd(0,0);
   mass_matrix_Fd(0,1) = 0;
   mass_matrix_Fd(0,2) = 0;
   mass_matrix_Fd(1,1) = 0.7071067811865475*vd*Yd(1,1);
   mass_matrix_Fd(1,2) = 0;
   mass_matrix_Fd(2,2) = 0.7071067811865475*vd*Yd(2,2);

   Symmetrize(mass_matrix_Fd);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Fd, eigenvalue_error > precision * Abs(
      MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*vu*Yu(0,0);
   mass_matrix_Fu(0,1) = 0;
   mass_matrix_Fu(0,2) = 0;
   mass_matrix_Fu(1,1) = 0.7071067811865475*vu*Yu(1,1);
   mass_matrix_Fu(1,2) = 0;
   mass_matrix_Fu(2,2) = 0.7071067811865475*vu*Yu(2,2);

   Symmetrize(mass_matrix_Fu);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::Fu, eigenvalue_error > precision * Abs(
      MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_FDX() const
{

   Eigen::Matrix<double,3,3> mass_matrix_FDX;

   mass_matrix_FDX(0,0) = 0.7071067811865475*vs*Kappa(0,0);
   mass_matrix_FDX(0,1) = 0;
   mass_matrix_FDX(0,2) = 0;
   mass_matrix_FDX(1,1) = 0.7071067811865475*vs*Kappa(1,1);
   mass_matrix_FDX(1,2) = 0;
   mass_matrix_FDX(2,2) = 0.7071067811865475*vs*Kappa(2,2);

   Symmetrize(mass_matrix_FDX);

   return mass_matrix_FDX;
}

void CLASSNAME::calculate_MFDX()
{
   const auto mass_matrix_FDX(get_mass_matrix_FDX());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_FDX, MFDX, ZDXL, ZDXR, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::FDX, eigenvalue_error > precision * Abs(
      MFDX(0)));
#else
   fs_svd(mass_matrix_FDX, MFDX, ZDXL, ZDXR);
#endif

}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_SHI0() const
{

   Eigen::Matrix<double,4,4> mass_matrix_SHI0;

   mass_matrix_SHI0(0,0) = mH1I2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHI0(0,1) = 0;
   mass_matrix_SHI0(0,2) = -0.7071067811865475*vs*Conj(TLambda12(0,0)) + 0.5*vd
      *vu*Conj(Lambda12(0,0))*Lambdax;
   mass_matrix_SHI0(0,3) = 0;
   mass_matrix_SHI0(1,1) = mH1I2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHI0(1,2) = 0;
   mass_matrix_SHI0(1,3) = -0.7071067811865475*vs*Conj(TLambda12(1,1)) + 0.5*vd
      *vu*Conj(Lambda12(1,1))*Lambdax;
   mass_matrix_SHI0(2,2) = mH2I2(0,0) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHI0(2,3) = 0;
   mass_matrix_SHI0(3,3) = mH2I2(1,1) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHI0);

   return mass_matrix_SHI0;
}

void CLASSNAME::calculate_MSHI0()
{
   const auto mass_matrix_SHI0(get_mass_matrix_SHI0());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHI0, MSHI0, UHI0, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::SHI0, eigenvalue_error > precision * Abs(
      MSHI0(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHI0, MSHI0, UHI0);
#endif
   normalize_to_interval(UHI0);


   if (MSHI0.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::SHI0);
   }

   MSHI0 = AbsSqrt(MSHI0);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_SHIp() const
{

   Eigen::Matrix<double,4,4> mass_matrix_SHIp;

   mass_matrix_SHIp(0,0) = mH1I2(0,0) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHIp(0,1) = 0;
   mass_matrix_SHIp(0,2) = 0.7071067811865475*vs*Conj(TLambda12(0,0)) - 0.5*vd*
      vu*Conj(Lambda12(0,0))*Lambdax;
   mass_matrix_SHIp(0,3) = 0;
   mass_matrix_SHIp(1,1) = mH1I2(1,1) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*
      Sqr(vd) + 0.1125*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.1875*Sqr(gN)*Sqr(vs) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) +
      0.075*Sqr(gN)*Sqr(vu);
   mass_matrix_SHIp(1,2) = 0;
   mass_matrix_SHIp(1,3) = 0.7071067811865475*vs*Conj(TLambda12(1,1)) - 0.5*vd*
      vu*Conj(Lambda12(1,1))*Lambdax;
   mass_matrix_SHIp(2,2) = mH2I2(0,0) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(0,0))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHIp(2,3) = 0;
   mass_matrix_SHIp(3,3) = mH2I2(1,1) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) + 0.5*AbsSqr(Lambda12(1,1))*Sqr(vs) -
      0.125*Sqr(gN)*Sqr(vs) + 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) +
      0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHIp);

   return mass_matrix_SHIp;
}

void CLASSNAME::calculate_MSHIp()
{
   const auto mass_matrix_SHIp(get_mass_matrix_SHIp());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHIp, MSHIp, UHIp, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::SHIp, eigenvalue_error > precision * Abs(
      MSHIp(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHIp, MSHIp, UHIp);
#endif
   normalize_to_interval(UHIp);


   if (MSHIp.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::SHIp);
   }

   MSHIp = AbsSqrt(MSHIp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_ChaI() const
{

   Eigen::Matrix<double,2,2> mass_matrix_ChaI;

   mass_matrix_ChaI(0,0) = 0.7071067811865475*vs*Lambda12(0,0);
   mass_matrix_ChaI(0,1) = 0;
   mass_matrix_ChaI(1,1) = 0.7071067811865475*vs*Lambda12(1,1);

   Symmetrize(mass_matrix_ChaI);

   return mass_matrix_ChaI;
}

void CLASSNAME::calculate_MChaI()
{
   const auto mass_matrix_ChaI(get_mass_matrix_ChaI());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_ChaI, MChaI, ZMI, ZPI, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::ChaI, eigenvalue_error > precision * Abs(
      MChaI(0)));
#else
   fs_svd(mass_matrix_ChaI, MChaI, ZMI, ZPI);
#endif

}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_ChiI() const
{

   Eigen::Matrix<double,4,4> mass_matrix_ChiI;

   mass_matrix_ChiI(0,0) = 0;
   mass_matrix_ChiI(0,1) = 0;
   mass_matrix_ChiI(0,2) = -0.7071067811865475*vs*Lambda12(0,0);
   mass_matrix_ChiI(0,3) = 0;
   mass_matrix_ChiI(1,1) = 0;
   mass_matrix_ChiI(1,2) = 0;
   mass_matrix_ChiI(1,3) = -0.7071067811865475*vs*Lambda12(1,1);
   mass_matrix_ChiI(2,2) = 0;
   mass_matrix_ChiI(2,3) = 0;
   mass_matrix_ChiI(3,3) = 0;

   Symmetrize(mass_matrix_ChiI);

   return mass_matrix_ChiI;
}

void CLASSNAME::calculate_MChiI()
{
   const auto mass_matrix_ChiI(get_mass_matrix_ChiI());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_ChiI, MChiI, ZNI, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::ChiI, eigenvalue_error > precision * Abs(
      MChiI(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_ChiI, MChiI, ZNI);
#endif
   normalize_to_interval(ZNI);

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_SSI0() const
{

   Eigen::Matrix<double,2,2> mass_matrix_SSI0;

   mass_matrix_SSI0(0,0) = msI2(0,0) - 0.1875*Sqr(gN)*Sqr(vd) + 0.3125*Sqr(gN)*
      Sqr(vs) - 0.125*Sqr(gN)*Sqr(vu);
   mass_matrix_SSI0(0,1) = 0;
   mass_matrix_SSI0(1,1) = msI2(1,1) - 0.1875*Sqr(gN)*Sqr(vd) + 0.3125*Sqr(gN)*
      Sqr(vs) - 0.125*Sqr(gN)*Sqr(vu);

   Symmetrize(mass_matrix_SSI0);

   return mass_matrix_SSI0;
}

void CLASSNAME::calculate_MSSI0()
{
   const auto mass_matrix_SSI0(get_mass_matrix_SSI0());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SSI0, MSSI0, ZSSI, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::SSI0, eigenvalue_error > precision * Abs(
      MSSI0(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SSI0, MSSI0, ZSSI);
#endif
   normalize_to_interval(ZSSI);


   if (MSSI0.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::SSI0);
   }

   MSSI0 = AbsSqrt(MSSI0);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_FSI() const
{

   Eigen::Matrix<double,2,2> mass_matrix_FSI;

   mass_matrix_FSI(0,0) = 0;
   mass_matrix_FSI(0,1) = 0;
   mass_matrix_FSI(1,1) = 0;

   Symmetrize(mass_matrix_FSI);

   return mass_matrix_FSI;
}

void CLASSNAME::calculate_MFSI()
{
   const auto mass_matrix_FSI(get_mass_matrix_FSI());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_FSI, MFSI, ZFSI, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::FSI, eigenvalue_error > precision * Abs(
      MFSI(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_FSI, MFSI, ZFSI);
#endif
   normalize_to_interval(ZFSI);

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_SHp0() const
{

   Eigen::Matrix<double,2,2> mass_matrix_SHp0;

   mass_matrix_SHp0(0,0) = mHp2 + AbsSqr(MuPr) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*
      Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*
      Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHp0(0,1) = -Conj(BMuPr);
   mass_matrix_SHp0(1,1) = mHpbar2 + AbsSqr(MuPr) - 0.075*Sqr(g1)*Sqr(vd) -
      0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) - 0.125*Sqr(gN)*Sqr(vs) +
      0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHp0);

   return mass_matrix_SHp0;
}

void CLASSNAME::calculate_MSHp0()
{
   const auto mass_matrix_SHp0(get_mass_matrix_SHp0());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHp0, MSHp0, UHp0, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::SHp0, eigenvalue_error > precision * Abs(
      MSHp0(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHp0, MSHp0, UHp0);
#endif
   normalize_to_interval(UHp0);


   if (MSHp0.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::SHp0);
   }

   MSHp0 = AbsSqrt(MSHp0);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_SHpp() const
{

   Eigen::Matrix<double,2,2> mass_matrix_SHpp;

   mass_matrix_SHpp(0,0) = mHp2 + AbsSqr(MuPr) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*
      Sqr(g2)*Sqr(vd) - 0.075*Sqr(gN)*Sqr(vd) + 0.125*Sqr(gN)*Sqr(vs) - 0.075*
      Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) - 0.05*Sqr(gN)*Sqr(vu);
   mass_matrix_SHpp(0,1) = Conj(BMuPr);
   mass_matrix_SHpp(1,1) = mHpbar2 + AbsSqr(MuPr) - 0.075*Sqr(g1)*Sqr(vd) +
      0.125*Sqr(g2)*Sqr(vd) + 0.075*Sqr(gN)*Sqr(vd) - 0.125*Sqr(gN)*Sqr(vs) +
      0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.05*Sqr(gN)*Sqr(vu);

   Hermitianize(mass_matrix_SHpp);

   return mass_matrix_SHpp;
}

void CLASSNAME::calculate_MSHpp()
{
   const auto mass_matrix_SHpp(get_mass_matrix_SHpp());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_SHpp, MSHpp, UHpp, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::SHpp, eigenvalue_error > precision * Abs(
      MSHpp(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_SHpp, MSHpp, UHpp);
#endif
   normalize_to_interval(UHpp);


   if (MSHpp.minCoeff() < 0.) {
      problems.flag_running_tachyon(E6SSM_info::SHpp);
   }

   MSHpp = AbsSqrt(MSHpp);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_ChiP() const
{

   Eigen::Matrix<double,2,2> mass_matrix_ChiP;

   mass_matrix_ChiP(0,0) = 0;
   mass_matrix_ChiP(0,1) = -MuPr;
   mass_matrix_ChiP(1,1) = 0;

   Symmetrize(mass_matrix_ChiP);

   return mass_matrix_ChiP;
}

void CLASSNAME::calculate_MChiP()
{
   const auto mass_matrix_ChiP(get_mass_matrix_ChiP());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_ChiP, MChiP, ZNp, eigenvalue_error);
   problems.flag_bad_mass(E6SSM_info::ChiP, eigenvalue_error > precision * Abs(
      MChiP(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_ChiP, MChiP, ZNp);
#endif
   normalize_to_interval(ZNp);

}

double CLASSNAME::get_mass_matrix_VWm() const
{

   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(vd) + Sqr(vu)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{

   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(E6SSM_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_VPVZVZp() const
{

   Eigen::Matrix<double,3,3> mass_matrix_VPVZVZp;

   mass_matrix_VPVZVZp(0,0) = 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_VPVZVZp(0,1) = -0.19364916731037085*g1*g2*Sqr(vd) -
      0.19364916731037085*g1*g2*Sqr(vu);
   mass_matrix_VPVZVZp(0,2) = 0.18371173070873834*g1*gN*Sqr(vd) -
      0.1224744871391589*g1*gN*Sqr(vu);
   mass_matrix_VPVZVZp(1,1) = 0.25*Sqr(g2)*Sqr(vd) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_VPVZVZp(1,2) = -0.23717082451262844*g2*gN*Sqr(vd) +
      0.15811388300841897*g2*gN*Sqr(vu);
   mass_matrix_VPVZVZp(2,2) = 0.225*Sqr(gN)*Sqr(vd) + 0.625*Sqr(gN)*Sqr(vs) +
      0.1*Sqr(gN)*Sqr(vu);

   Symmetrize(mass_matrix_VPVZVZp);

   return mass_matrix_VPVZVZp;
}

void CLASSNAME::calculate_MVPVZVZp()
{
   const auto mass_matrix_VPVZVZp(get_mass_matrix_VPVZVZp());
   Eigen::Array<double,3,1> MVPVZVZp;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZVZp, MVPVZVZp, ZZ, eigenvalue_error
      );
#else

   fs_diagonalize_hermitian(mass_matrix_VPVZVZp, MVPVZVZp, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);


   MVPVZVZp = AbsSqrt(MVPVZVZp);

   MVP = 0.;
   MVZ = MVPVZVZp(1);
   MVZp = MVPVZVZp(2);
}



double CLASSNAME::get_ewsb_eq_hh_1() const
{
   
   double result = Re(mHd2*vd - 0.35355339059327373*vs*vu*Conj(TLambdax) + 0.075*
      Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.1125*Cube(vd)*Sqr(gN) + 0.5*vd
      *AbsSqr(Lambdax)*Sqr(vs) - 0.1875*vd*Sqr(gN)*Sqr(vs) + 0.5*vd*AbsSqr(Lambdax
      )*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu) + 0.075*vd*
      Sqr(gN)*Sqr(vu) - 0.35355339059327373*vs*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(mHu2*vu - 0.35355339059327373*vd*vs*Conj(TLambdax) + 0.075*
      Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) + 0.05*Cube(vu)*Sqr(gN) + 0.5*vu*
      AbsSqr(Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd
      ) + 0.075*vu*Sqr(gN)*Sqr(vd) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vs) - 0.125*vu*Sqr
      (gN)*Sqr(vs) - 0.35355339059327373*vd*vs*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = Re(ms2*vs - 0.35355339059327373*vd*vu*Conj(TLambdax) + 0.3125*
      Cube(vs)*Sqr(gN) + 0.5*vs*AbsSqr(Lambdax)*Sqr(vd) - 0.1875*vs*Sqr(gN)*Sqr(vd
      ) + 0.5*vs*AbsSqr(Lambdax)*Sqr(vu) - 0.125*vs*Sqr(gN)*Sqr(vu) -
      0.35355339059327373*vd*vu*TLambdax);

   return result;
}



std::complex<double> CLASSNAME::CpUSdconjUSdVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.31622776601683794*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 3,-0.08164965809277262*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,0.2581988897471611*g1*g2*Cos(ThetaW
      ())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(gO1 <
      3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()))
      ,0) + IF(gO1 < 3,0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos
      (ThetaWp()))*Sqr(Sin(ThetaW())),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)
      *Sqr(gN)*Sqr(Sin(ThetaWp())),0) - 0.32659863237109044*g1*gN*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 0.13333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp())
      )*Sqr(Sin(ThetaW()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1)) + 0.2*Sqr(gN)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.31622776601683794*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 3,0.08164965809277262*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,0.2581988897471611*g1*g2*Cos(ThetaW())
      *KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0
      ) + IF(gO1 < 3,0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())),0) + 0.32659863237109044*g1*gN*Cos(ThetaWp())
      *Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 0.2*Sqr(gN)*Sqr(Cos(ThetaWp()))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.13333333333333333
      *Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSdconjUSdconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSdconjHpmconjUSd(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-(AbsSqr(Yu(gO1,gO1))*KroneckerDelta(
      gO1,gO2)*ZP(gI1,1)*ZP(gI2,1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr
      (g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*ZP(gI1,1)*ZP(gI2,1),0) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)
      *KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) + 0.15*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) -
      SUM(j2,0,2,AbsSqr(Yd(j2,j2))*KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
       + j2))*ZP(gI1,0)*ZP(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) + 0.1*Sqr(gN)*SUM(j1,0,2
      ,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSdSHp0conjUSdconjSHp0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*Conj(UHp0(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,0),0) + IF(gO1 < 3,0.25*Conj(UHp0(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,0),0) + IF(gO1 < 3,-0.05*
      Conj(UHp0(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,0),0) + IF(gO1 <
      3,-0.05*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,1),0) +
      IF(gO1 < 3,-0.25*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,
      1),0) + IF(gO1 < 3,0.05*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      UHp0(gI2,1),0) + 0.1*Conj(UHp0(gI1,0))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,0) - 0.1*Conj(UHp0(gI1,0))*Sqr
      (gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(
      gI2,0) - 0.1*Conj(UHp0(gI1,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)
      *KroneckerDelta(gO2,3 + j1))*UHp0(gI2,1) + 0.1*Conj(UHp0(gI1,1))*Sqr(gN)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSdSHppconjUSdconjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*Conj(UHpp(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,0),0) + IF(gO1 < 3,-0.25*Conj(UHpp(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,0),0) + IF(gO1 < 3,-0.05*
      Conj(UHpp(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,0),0) + IF(gO1 <
      3,-0.05*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,1),0) +
      IF(gO1 < 3,0.25*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,1
      ),0) + IF(gO1 < 3,0.05*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      UHpp(gI2,1),0) + 0.1*Conj(UHpp(gI1,0))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,0) - 0.1*Conj(UHpp(gI1,0))*Sqr
      (gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(
      gI2,0) - 0.1*Conj(UHpp(gI1,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)
      *KroneckerDelta(gO2,3 + j1))*UHpp(gI2,1) + 0.1*Conj(UHpp(gI1,1))*Sqr(gN)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSdSSI0conjUSdconjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.125*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN),0) - 0.25*KroneckerDelta(gI1,gI2)*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(AbsSqr(Yd(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZA(gI1,0)*ZA(gI2,0)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.5*
      Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Yd(gO1,gO1)*ZA(gI1,2)*ZA(gI2,1),0)
      + IF(gO1 < 3,-0.5*Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Yd(gO1,gO1)*ZA(
      gI1,1)*ZA(gI2,2),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(
      gI1,2)*ZA(gI2,2),0) + IF(gO2 < 3,-0.5*Conj(Yd(gO2,gO2))*KroneckerDelta(gO1,3
       + gO2)*Lambdax*ZA(gI1,2)*ZA(gI2,1),0) + IF(gO2 < 3,-0.5*Conj(Yd(gO2,gO2))*
      KroneckerDelta(gO1,3 + gO2)*Lambdax*ZA(gI1,1)*ZA(gI2,2),0) + 0.1*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(
      gI2,0) + 0.15*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) - SUM(j2,0,2,AbsSqr(Yd(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZA(gI1,0)*ZA(gI2,0) -
      0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*ZA(gI1,1)*ZA(gI2,1) + 0.1*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1) - 0.25*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,2)*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(AbsSqr(Yd(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZH(gI1,0)*ZH(gI2,0)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.5*Conj
      (Lambdax)*KroneckerDelta(3 + gO1,gO2)*Yd(gO1,gO1)*ZH(gI1,2)*ZH(gI2,1),0) +
      IF(gO1 < 3,0.5*Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Yd(gO1,gO1)*ZH(gI1,
      1)*ZH(gI2,2),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,2
      )*ZH(gI2,2),0) + IF(gO2 < 3,0.5*Conj(Yd(gO2,gO2))*KroneckerDelta(gO1,3 + gO2
      )*Lambdax*ZH(gI1,2)*ZH(gI2,1),0) + IF(gO2 < 3,0.5*Conj(Yd(gO2,gO2))*
      KroneckerDelta(gO1,3 + gO2)*Lambdax*ZH(gI1,1)*ZH(gI2,2),0) + 0.1*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(
      gI2,0) + 0.15*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) - SUM(j2,0,2,AbsSqr(Yd(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZH(gI1,0)*ZH(gI2,0) -
      0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*ZH(gI1,1)*ZH(gI2,1) + 0.1*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1) - 0.25*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,2)*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUSdSvconjUSdconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gI1,gI2)
      *KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gI1,
      gI2)*KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.1*KroneckerDelta(gI1,gI2)*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.1*
      KroneckerDelta(gI1,gI2)*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChaFuconjUSdPR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(Yu(gO2,gO2))*UP(gI2,1)*ZUR(
      gI1,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpChaFuconjUSdPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UM(gI2,0))*Conj(ZUL(
      gI1,gO1))),0) + Conj(UM(gI2,1))*SUM(j1,0,2,Conj(ZUL(gI1,j1))*KroneckerDelta(
      gO1,3 + j1)*Yd(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(Conj(Yd(gO2,gO2))*ZDR(gI1,gO2)
      *ZN(gI2,2)),0) - 0.14907119849998596*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*
      ZDR(gI1,j1))*(2.449489742783178*g1*ZN(gI2,0) + 3*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZDL
      (gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZDL(
      gI1,gO1))*Conj(ZN(gI2,1)),0) + IF(gO1 < 3,-0.22360679774997896*gN*Conj(ZDL(
      gI1,gO1))*Conj(ZN(gI2,5)),0) - Conj(ZN(gI2,2))*SUM(j1,0,2,Conj(ZDL(gI1,j1))*
      KroneckerDelta(gO1,3 + j1)*Yd(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSdSHI0conjUSdconjSHI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),
      0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(
      gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,-0.125
      *KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,
      2 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),
      0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(
      gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) +
      0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(
      gI1,2 + j2))*UHI0(gI2,2 + j2)) + 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(
      gI2,2 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.075*Sqr(gN
      )*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,
      3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,
      2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.05*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))
      *UHI0(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,
      3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSdSHIpconjUSdconjSHIp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),
      0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(
      gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2
       + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),
      0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(
      gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) +
      0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(
      gI1,2 + j2))*UHIp(gI2,2 + j2)) + 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(
      gI2,2 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.075*Sqr(gN
      )*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,
      3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,
      2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.05*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))
      *UHIp(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,
      3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666
      *Conj(ZD(gI2,gO2))*Sqr(g1)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*
      Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-
      1.3333333333333333*Conj(ZD(gI2,gO2))*Sqr(g3)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,
      IF(gO2 < 3,-0.025*Conj(ZD(gI2,gO2))*Sqr(gN)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,
      IF(gO2 < 3,-(Conj(Yd(gO2,gO2))*Conj(ZD(gI2,3 + gO2))*Yd(gO1,gO1)*ZD(gI1,3 +
      gO1)),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.375*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.0375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 +
      j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) +
      IF(gO1 < 3,-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI2,j2)
      )*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2
      ,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,
      gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + IF(gO1 < 3
      ,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(
      gI1,3 + j2)),0) + IF(gO1 < 3,-3*KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,2,Conj(
      Yd(j2,j2))*Conj(ZD(gI2,3 + j2))*ZD(gI1,j2))*Yd(gO1,gO1),0) + IF(gO1 < 3,-
      0.016666666666666666*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j1,0
      ,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1
      < 3,-0.025*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1
      ))*ZD(gI1,gO1),0) + IF(gO1 < 3,-0.016666666666666666*Sqr(g1)*SUM(j2,0,2,Conj
      (ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO1 < 3,
      0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(
      gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO1 < 3,-0.025*Sqr(gN)*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO2 < 3,-
      0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZD(gI2,gO2)
      )*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2
      < 3,-0.025*Conj(ZD(gI2,gO2))*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      ZD(gI1,3 + j1)),0) + IF(gO2 < 3,-3*Conj(Yd(gO2,gO2))*KroneckerDelta(gO1,3 +
      gO2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*Yd(j1,j1)*ZD(gI1,3 + j1)),0) + IF(gO2 < 3,-
      0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1
      ,3 + j2)*ZD(gI1,3 + j2)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZD(gI2,gO2)
      )*Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)),0) + IF(gO2
      < 3,-0.025*Conj(ZD(gI2,gO2))*Sqr(gN)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      ZD(gI1,3 + j2)),0) - 0.03333333333333333*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(
      gO2,3 + j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 +
      j2)) - 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1))*
      SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*
      SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2
      )*KroneckerDelta(gO2,3 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD
      (gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))
      - 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.15*Sqr(gN)*SUM(j1
      ,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2))
      - 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) - SUM(j1,0,2,Conj(ZD(gI2,j1))*
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*ZD(gI1,j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*ZD(gI1,3 + j2)) - 0.15*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 +
      j2)) - 0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 +
      j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 +
      j2)) - 0.05*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSDXSDX(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)),0)
      + 0.025*(40*IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj
      (ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)),0) + 40*IF(gO1 < 3,0.1125*KroneckerDelta(
      gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)),0) + 40*
      IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI2,j2))
      *ZDX(gI1,j2)),0) + 40*IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(
      j2,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2)),0) + 40*IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 +
      j2)),0) + 40*IF(gO1 < 3,0.1125*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,
      Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)),0) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(
      gI2,j1))*ZDX(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) + 6*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1))*SUM(j2,0,
      2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 4*Sqr(g1)*SUM(j1,
      0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) + 9*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))
      *ZDX(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
       + j2)) + 4*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2
      ,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2)) + 6*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(
      gI2,j2))*ZDX(gI1,j2)) - 4*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)
      ) + 9*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZU(gI2,gO2)
      )*Sqr(g2)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-(Conj(Yu(gO2,gO2))*Conj
      (ZU(gI2,3 + gO2))*Yu(gO1,gO1)*ZU(gI1,3 + gO1)),0),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) +
      IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI2,j1))
      *ZU(gI1,j1)),0) + IF(gO1 < 3,-0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,
      0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-
      0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(
      gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0
      ,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2
      )*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,-0.0375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2
      ))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,-0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(
      gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) +
      0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.075*Sqr(gN)*SUM(
      j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,
      j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) - SUM(j1,0,2,Conj(ZU(
      gI2,j1))*KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*SUM(j2,0,2,Conj(Yd(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*ZU(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*ZU(gI1,3 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3
      + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSdSeconjUSdconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1)
      )*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0
      ,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-0.0125*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-0.0125*KroneckerDelta(gO1,gO2)*Sqr(gN)
      *SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(
      KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gI1,3 + j2))*
      ZE(gI2,j2))*Yd(gO1,gO1)),0) + IF(gO2 < 3,-(Conj(Yd(gO2,gO2))*KroneckerDelta(
      gO1,3 + gO2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*Ye(j1,j1)*ZE(gI2,3 + j1))),0) +
      0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(gN)*SUM(j1
      ,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE
      (gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.025*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,
      0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.05*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2
      ,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) -
      0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.025*Sqr(gN)*SUM(j1,0,2
      ,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpHpmSuconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*vd*AbsSqr(Yd(
      gO2,gO2))*Conj(ZU(gI1,gO2))*ZP(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*vs*
      Conj(Yu(gO2,gO2))*Conj(ZU(gI1,3 + gO2))*Lambdax*ZP(gI2,0),0) + IF(gO2 < 3,-
      0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,
      0.7071067811865475*vu*AbsSqr(Yu(gO2,gO2))*Conj(ZU(gI1,gO2))*ZP(gI2,1),0) +
      IF(gO2 < 3,Conj(ZU(gI1,3 + gO2))*Conj(TYu(gO2,gO2))*ZP(gI2,1),0) + IF(gO2 <
      3,-0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1),0) + SUM(j1,0,
      2,Conj(ZU(gI1,j1))*KroneckerDelta(gO2,3 + j1)*TYd(j1,j1))*ZP(gI2,0) +
      0.7071067811865475*vu*SUM(j2,0,2,Conj(Yu(j2,j2))*Conj(ZU(gI1,3 + j2))*
      KroneckerDelta(gO2,3 + j2)*Yd(j2,j2))*ZP(gI2,0) + 0.7071067811865475*vs*Conj
      (Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*
      ZP(gI2,1) + 0.7071067811865475*vd*SUM(j2,0,2,Conj(Yu(j2,j2))*Conj(ZU(gI1,3 +
      j2))*KroneckerDelta(gO2,3 + j2)*Yd(j2,j2))*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhSdconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZD(gI1,3 + gO2))*Conj(TYd(gO2,gO2))*ZA(gI2,0),0) +
      IF(gO2 < 3,std::complex<double>(0,0.5)*vs*Conj(Yd(gO2,gO2))*Conj(ZD(gI1,3 +
      gO2))*Lambdax*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(0,0.5)*vu*Conj(
      Yd(gO2,gO2))*Conj(ZD(gI1,3 + gO2))*Lambdax*ZA(gI2,2),0) - std::complex<
      double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(ZD(gI1,j1))*KroneckerDelta(
      gO2,3 + j1)*TYd(j1,j1))*ZA(gI2,0) - std::complex<double>(0,0.5)*vs*Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZD(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*
      ZA(gI2,1) - std::complex<double>(0,0.5)*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(
      gI1,j1))*KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSdconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(vd*AbsSqr(Yd(gO2,gO2))*Conj(ZD
      (gI1,gO2))*ZH(gI2,0)),0) + IF(gO2 < 3,-0.7071067811865475*Conj(ZD(gI1,3 +
      gO2))*Conj(TYd(gO2,gO2))*ZH(gI2,0),0) + IF(gO2 < 3,0.05*vd*Conj(ZD(gI1,gO2))
      *Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 3,0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,
      0),0) + IF(gO2 < 3,0.075*vd*Conj(ZD(gI1,gO2))*Sqr(gN)*ZH(gI2,0),0) + IF(gO2
      < 3,0.5*vs*Conj(Yd(gO2,gO2))*Conj(ZD(gI1,3 + gO2))*Lambdax*ZH(gI2,1),0) + IF
      (gO2 < 3,-0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gO2 < 3,-0.25*
      vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF(gO2 < 3,0.05*vu*Conj(ZD(gI1,
      gO2))*Sqr(gN)*ZH(gI2,1),0) + IF(gO2 < 3,0.5*vu*Conj(Yd(gO2,gO2))*Conj(ZD(gI1
      ,3 + gO2))*Lambdax*ZH(gI2,2),0) + IF(gO2 < 3,-0.125*vs*Conj(ZD(gI1,gO2))*Sqr
      (gN)*ZH(gI2,2),0) + 0.1*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) + 0.15*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZD(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) - 0.7071067811865475*SUM(
      j1,0,2,Conj(ZD(gI1,j1))*KroneckerDelta(gO2,3 + j1)*TYd(j1,j1))*ZH(gI2,0) -
      vd*SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(gI1,3 + j2))*KroneckerDelta(gO2,3 +
      j2))*ZH(gI2,0) - 0.1*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) + 0.1*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZD(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) + 0.5*vs*Conj(Lambdax)*
      SUM(j1,0,2,Conj(ZD(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*ZH(gI2,1)
      - 0.25*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)
      )*ZH(gI2,2) + 0.5*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j1))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpGluFdconjUSdPR(int gI2, int gO2) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*Conj(PhaseGlu)*SUM(j1
      ,0,2,KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFdconjUSdPL(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*PhaseGlu*
      Conj(ZDL(gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVG(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 6,g3*Conj(ZD(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.12909944487358055*g1*Conj(ZD(
      gI2,gO2))*Cos(ThetaW()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZD(gI2,gO2))*Sin(ThetaW
      ()),0) - 0.2581988897471611*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))
      *KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Conj(ZD(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZD(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,0.15811388300841897*gN*
      Conj(ZD(gI2,gO2))*Sin(ThetaWp()),0) + 0.2581988897471611*g1*Cos(ThetaWp())*
      Sin(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) -
      0.31622776601683794*gN*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.15811388300841897*gN*Conj(ZD(
      gI2,gO2))*Cos(ThetaWp()),0) + IF(gO2 < 3,0.5*g2*Conj(ZD(gI2,gO2))*Cos(ThetaW
      ())*Sin(ThetaWp()),0) + IF(gO2 < 3,0.12909944487358055*g1*Conj(ZD(gI2,gO2))*
      Sin(ThetaW())*Sin(ThetaWp()),0) - 0.31622776601683794*gN*Cos(ThetaWp())*SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) - 0.2581988897471611
      *g1*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSdVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZU(
      gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(gO1,gO2)*(
      3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*
      g1*Sin(ThetaW())*Sin(2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(gO1,gO2)*(-
      3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 2*Sqr(gN)*Sqr(Cos(
      ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*(-
      2.449489742783178*gN*Sin(ThetaW())*Sin(2*ThetaWp()) + 3.872983346207417*g2*
      Sin(2*ThetaW())*Sqr(Sin(ThetaWp())) + 3*g1*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSvconjHpmconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(AbsSqr(Ye(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZP(gI1,0)*ZP(gI2,0)),0) + 0.05*KroneckerDelta(gO1,
      gO2)*((-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0) + (3*Sqr(g1)
      - 5*Sqr(g2) + 2*Sqr(gN))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0USvconjSHp0conjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) +
      5*Sqr(g2) + 2*Sqr(gN))*(Conj(UHp0(gI1,0))*UHp0(gI2,0) - Conj(UHp0(gI1,1))*
      UHp0(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppUSvconjSHppconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) -
      5*Sqr(g2) + 2*Sqr(gN))*(Conj(UHpp(gI1,0))*UHpp(gI2,0) - Conj(UHpp(gI1,1))*
      UHpp(gI2,1));

   return result;
}

double CLASSNAME::CpSSI0USvconjSSI0conjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const double result = -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr
      (gN);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFeconjUSvPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(Ye(gO2,gO2))*UM(gI1,1)*ZER(
      gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFeconjUSvPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UP(gI1,0))*Conj(ZEL(
      gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpSeconjHpmconjUSv(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*vd*AbsSqr(Ye(
      gO2,gO2))*Conj(ZE(gI2,gO2))*ZP(gI1,0),0) + IF(gO2 < 3,Conj(ZE(gI2,3 + gO2))*
      Conj(TYe(gO2,gO2))*ZP(gI1,0),0) + IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZE
      (gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vs*Conj(Ye(
      gO2,gO2))*Conj(ZE(gI2,3 + gO2))*Lambdax*ZP(gI1,1),0) + IF(gO2 < 3,-
      0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((-3*Sqr(g1) -
      5*Sqr(g2) + 3*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) + (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(
      gN))*ZA(gI1,1)*ZA(gI2,1) - 5*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((-3*Sqr(g1) -
      5*Sqr(g2) + 3*Sqr(gN))*ZH(gI1,0)*ZH(gI2,0) + (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(
      gN))*ZH(gI1,1)*ZH(gI2,1) - 5*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpSvUSvconjSvconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,IF(gI2 < 3,-0.15*Conj(ZV(gI1,gO2
      ))*Sqr(g1)*ZV(gI2,gO1),0),0) + IF(gI1 < 3,IF(gI2 < 3,-0.25*Conj(ZV(gI1,gO2))
      *Sqr(g2)*ZV(gI2,gO1),0),0) + IF(gI1 < 3,IF(gI2 < 3,-0.1*Conj(ZV(gI1,gO2))*
      Sqr(gN)*ZV(gI2,gO1),0),0) - 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,
      gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CphhSvconjUSv(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-0.15*vd*Conj(ZV(gI1,gO2))*Sqr(
      g1)*ZH(gI2,0),0) + IF(gI1 < 3,-0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,0),0
      ) + IF(gI1 < 3,0.15*vd*Conj(ZV(gI1,gO2))*Sqr(gN)*ZH(gI2,0),0) + IF(gI1 < 3,
      0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gI1 < 3,0.25*vu*Conj(ZV(
      gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF(gI1 < 3,0.1*vu*Conj(ZV(gI1,gO2))*Sqr(gN)
      *ZH(gI2,1),0) + IF(gI1 < 3,-0.25*vs*Conj(ZV(gI1,gO2))*Sqr(gN)*ZH(gI2,2),0);

   return result;
}

double CLASSNAME::CpChiFvconjUSvPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpChiFvconjUSvPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gI1 < 3,0.5477225575051661*g1*Conj(ZN(
      gI2,0))*KroneckerDelta(gI1,gO1),0) + IF(gI1 < 3,-0.7071067811865475*g2*Conj(
      ZN(gI2,1))*KroneckerDelta(gI1,gO1),0) + IF(gI1 < 3,-0.4472135954999579*gN*
      Conj(ZN(gI2,5))*KroneckerDelta(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USvconjSHI0conjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((-3*Sqr(g1) -
      5*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) + (3*Sqr(
      g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 +
      j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpUSvconjSHIpconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((-3*Sqr(g1) +
      5*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) + (3*Sqr(
      g1) - 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 +
      j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSvconjSdconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((Sqr(g1) + 5*
      Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) -
      Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSvconjSDXconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*(-2*(Sqr(g1) -
      Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) + (2*Sqr(g1) + 3*Sqr(gN))
      *SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSvconjSeconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZE(gI1,gO2)
      )*Sqr(g2)*ZE(gI2,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-(Conj(Ye(gO2,gO2))*Conj
      (ZE(gI1,3 + gO2))*Ye(gO1,gO1)*ZE(gI2,3 + gO1)),0),0) - 0.05*KroneckerDelta(
      gO1,gO2)*((3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE
      (gI2,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3
      + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuUSvconjSuconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((Sqr(g1) - 5*
      Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - (4*Sqr(g1) +
      Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSvVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZV(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,0.3872983346207417*g1*Conj(ZV(gI2,
      gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,0.31622776601683794*gN*
      Conj(ZV(gI2,gO2))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSvVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.31622776601683794*gN*Conj(ZV(
      gI2,gO2))*Cos(ThetaWp()),0) + IF(gI2 < 3,-0.5*g2*Conj(ZV(gI2,gO2))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(ZV(gI2,
      gO2))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSvconjVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZE(
      gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.31622776601683794*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 3,-0.08164965809277262*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,-0.2581988897471611*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(
      gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())),0) + IF(gO1 < 3,0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(Sin(ThetaWp())),0) + 0.32659863237109044
      *g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta
      (gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.5333333333333333*Sqr(g1)*Sqr(
      Cos(ThetaWp()))*Sqr(Sin(ThetaW()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 0.05*Sqr(gN)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.31622776601683794*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 3,0.08164965809277262*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,-0.2581988897471611*g1*g2*Cos(ThetaW()
      )*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0
      ) + IF(gO1 < 3,0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())),0) - 0.32659863237109044*g1*gN*Cos(ThetaWp())
      *Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 0.05*Sqr(gN)*Sqr(Cos(ThetaWp()))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.5333333333333333*
      Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSuconjUSuconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSuconjHpmconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(AbsSqr(Yd(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZP(gI1,0)*ZP(gI2,0)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZP(gI1,1)*ZP(gI2,1),0) - 0.2*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,
      0) + 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,
      3 + j1))*ZP(gI1,0)*ZP(gI2,0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) + 0.05*Sqr(gN)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1)
      - SUM(j2,0,2,AbsSqr(Yu(j2,j2))*KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2
      ,3 + j2))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSHp0USuconjSHp0conjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*Conj(UHp0(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,0),0) + IF(gO1 < 3,-0.25*Conj(UHp0(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,0),0) + IF(gO1 < 3,-0.05*
      Conj(UHp0(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,0),0) + IF(gO1 <
      3,-0.05*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,1),0) +
      IF(gO1 < 3,0.25*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,1
      ),0) + IF(gO1 < 3,0.05*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      UHp0(gI2,1),0) - 0.2*Conj(UHp0(gI1,0))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,0) - 0.05*Conj(UHp0(gI1,0))*
      Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*
      UHp0(gI2,0) + 0.2*Conj(UHp0(gI1,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,1) + 0.05*Conj(UHp0(gI1,1))*Sqr(
      gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(
      gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSHppUSuconjSHppconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*Conj(UHpp(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,0),0) + IF(gO1 < 3,0.25*Conj(UHpp(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,0),0) + IF(gO1 < 3,-0.05*
      Conj(UHpp(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,0),0) + IF(gO1 <
      3,-0.05*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,1),0) +
      IF(gO1 < 3,-0.25*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,
      1),0) + IF(gO1 < 3,0.05*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      UHpp(gI2,1),0) - 0.2*Conj(UHpp(gI1,0))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,0) - 0.05*Conj(UHpp(gI1,0))*
      Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*
      UHpp(gI2,0) + 0.2*Conj(UHpp(gI1,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,1) + 0.05*Conj(UHpp(gI1,1))*Sqr(
      gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(
      gI2,1);

   return result;
}

double CLASSNAME::CpSSI0USuconjSSI0conjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const double result = -0.125*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*
      Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFdconjUSuPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(Yd(gO2,gO2))*UM(gI1,1)*ZDR(
      gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFdconjUSuPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UP(gI1,0))*Conj(ZDL(
      gI2,gO1))),0) + Conj(UP(gI1,1))*SUM(j1,0,2,Conj(ZDL(gI2,j1))*KroneckerDelta(
      gO1,3 + j1)*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjHpmconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*vd*AbsSqr(Yd(
      gO2,gO2))*Conj(ZD(gI2,gO2))*ZP(gI1,0),0) + IF(gO2 < 3,Conj(ZD(gI2,3 + gO2))*
      Conj(TYd(gO2,gO2))*ZP(gI1,0),0) + IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZD
      (gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vu*AbsSqr(Yu
      (gO2,gO2))*Conj(ZD(gI2,gO2))*ZP(gI1,1),0) + IF(gO2 < 3,0.7071067811865475*vs
      *Conj(Yd(gO2,gO2))*Conj(ZD(gI2,3 + gO2))*Lambdax*ZP(gI1,1),0) + IF(gO2 < 3,-
      0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) +
      0.7071067811865475*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gI2,j1))*
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j1))*ZP(gI1,0) + 0.7071067811865475*vu*SUM(
      j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)*Yu(j2
      ,j2))*ZP(gI1,0) + SUM(j1,0,2,Conj(ZD(gI2,j1))*KroneckerDelta(gO2,3 + j1)*TYu
      (j1,j1))*ZP(gI1,1) + 0.7071067811865475*vd*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(
      ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)*Yu(j2,j2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.5*Conj(Lambdax)*KroneckerDelta(3 +
      gO1,gO2)*Yu(gO1,gO1)*ZA(gI1,2)*ZA(gI2,0),0) + IF(gO1 < 3,-(AbsSqr(Yu(gO1,gO1
      ))*KroneckerDelta(gO1,gO2)*ZA(gI1,1)*ZA(gI2,1)),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.5*
      Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Yu(gO1,gO1)*ZA(gI1,0)*ZA(gI2,2),0)
      + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2),0) +
      IF(gO2 < 3,-0.5*Conj(Yu(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*Lambdax*ZA(gI1
      ,2)*ZA(gI2,0),0) + IF(gO2 < 3,-0.5*Conj(Yu(gO2,gO2))*KroneckerDelta(gO1,3 +
      gO2)*Lambdax*ZA(gI1,0)*ZA(gI2,2),0) - 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) + 0.075*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*
      ZA(gI2,0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1) + 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1) - SUM(j2,0,2,
      AbsSqr(Yu(j2,j2))*KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZA(
      gI1,1)*ZA(gI2,1) - 0.125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZA(gI1,2)*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*KroneckerDelta(3 +
      gO1,gO2)*Yu(gO1,gO1)*ZH(gI1,2)*ZH(gI2,0),0) + IF(gO1 < 3,-(AbsSqr(Yu(gO1,gO1
      ))*KroneckerDelta(gO1,gO2)*ZH(gI1,1)*ZH(gI2,1)),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.5*Conj
      (Lambdax)*KroneckerDelta(3 + gO1,gO2)*Yu(gO1,gO1)*ZH(gI1,0)*ZH(gI2,2),0) +
      IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2),0) +
      IF(gO2 < 3,0.5*Conj(Yu(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*Lambdax*ZH(gI1,
      2)*ZH(gI2,0),0) + IF(gO2 < 3,0.5*Conj(Yu(gO2,gO2))*KroneckerDelta(gO1,3 +
      gO2)*Lambdax*ZH(gI1,0)*ZH(gI2,2),0) - 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) + 0.075*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*
      ZH(gI2,0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1) + 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1) - SUM(j2,0,2,
      AbsSqr(Yu(j2,j2))*KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZH(
      gI1,1)*ZH(gI2,1) - 0.125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZH(gI1,2)*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUSuSvconjUSuconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gI1,gI2
      )*KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gI1,
      gI2)*KroneckerDelta(gO1,gO2)*Sqr(gN),0) - 0.2*KroneckerDelta(gI1,gI2)*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.05*
      KroneckerDelta(gI1,gI2)*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(Conj(Yu(gO2,gO2))*ZN(gI2,3)*
      ZUR(gI1,gO2)),0) + 0.07453559924999298*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)
      *ZUR(gI1,j1))*(9.797958971132712*g1*ZN(gI2,0) - 3*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZN(
      gI2,0))*Conj(ZUL(gI1,gO1)),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZN(
      gI2,1))*Conj(ZUL(gI1,gO1)),0) + IF(gO1 < 3,-0.22360679774997896*gN*Conj(ZN(
      gI2,5))*Conj(ZUL(gI1,gO1)),0) - Conj(ZN(gI2,3))*SUM(j1,0,2,Conj(ZUL(gI1,j1))
      *KroneckerDelta(gO1,3 + j1)*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USuconjSHI0conjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),
      0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(
      gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2
       + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),
      0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(
      gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) +
      0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(
      gI1,2 + j2))*UHI0(gI2,2 + j2)) + 0.025*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(
      gI2,2 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.0375*Sqr(
      gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(
      gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.025*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1)
      )*UHI0(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2
      ,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpUSuconjSHIpconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),
      0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(
      gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,-0.125
      *KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,
      2 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),
      0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(
      gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) +
      0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(
      gI1,2 + j2))*UHIp(gI2,2 + j2)) + 0.025*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(
      gI2,2 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.0375*Sqr(
      gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(
      gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.025*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1)
      )*UHIp(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2
      ,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSuconjSeconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1)
      )*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0
      ,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-0.0125*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-0.0125*KroneckerDelta(gO1,gO2)*Sqr(gN)
      *SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) - 0.1*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.025*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(
      gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) +
      0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.0125*Sqr(gN)*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,
      j2)) - 0.025*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) + 0.2*Sqr(g1)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3
       + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZD(gI2,gO2)
      )*Sqr(g2)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-(Conj(Yd(gO2,gO2))*Conj
      (ZD(gI2,3 + gO2))*Yd(gO1,gO1)*ZD(gI1,3 + gO1)),0),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) +
      IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))
      *ZD(gI1,j1)),0) + IF(gO1 < 3,-0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,
      0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,
      gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3
      ,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(
      gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0
      ,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2
      )*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.0375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,3 +
      j2))*ZD(gI1,3 + j2)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + 0.1*Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.0375*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(
      gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) +
      0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.075*Sqr(gN)*SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,
      j2)) - 0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) - SUM(j1,0,2,Conj(ZD(
      gI2,j1))*KroneckerDelta(gO2,3 + j1)*Yu(j1,j1))*SUM(j2,0,2,Conj(Yu(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*ZD(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*ZD(gI1,3 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3
      + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSDXSDX(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)),0)
      + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3
       + j1))*ZDX(gI1,3 + j1)),0) + IF(gO1 < 3,0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2)),0)
      + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI2,
      j2))*ZDX(gI1,j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)),0) + IF(gO1 < 3,0.1125*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 +
      j2)),0) - 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.075*Sqr(gN)*SUM(
      j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*
      ZDX(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
      + j2)) + 0.1125*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.2*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2
      ,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2)) + 0.075*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(
      gI2,j2))*ZDX(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)
      ) + 0.1125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,
      3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666
      *Conj(ZU(gI2,gO2))*Sqr(g1)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*
      Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-
      1.3333333333333333*Conj(ZU(gI2,gO2))*Sqr(g3)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,
      IF(gO2 < 3,-0.025*Conj(ZU(gI2,gO2))*Sqr(gN)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,
      IF(gO2 < 3,-(Conj(Yu(gO2,gO2))*Conj(ZU(gI2,3 + gO2))*Yu(gO1,gO1)*ZU(gI1,3 +
      gO1)),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,
      Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,-0.375*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,-0.0375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1
      ))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) +
      IF(gO1 < 3,-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI2,j2)
      )*ZU(gI1,j2)),0) + IF(gO1 < 3,-0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2
      ,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,-
      0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(
      gI1,3 + j2)),0) + IF(gO1 < 3,-3*KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,2,Conj(
      Yu(j2,j2))*Conj(ZU(gI2,3 + j2))*ZU(gI1,j2))*Yu(gO1,gO1),0) + IF(gO1 < 3,
      0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j1,0
      ,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1
      < 3,-0.0125*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 +
      j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,0.03333333333333333*Sqr(g1)*SUM(j2,0,2,Conj
      (ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO1 < 3,
      0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*KroneckerDelta(
      gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO1 < 3,-0.0125*Sqr(gN)*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO2 < 3,
      0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*ZU(gI1,3 + j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZU(gI2,gO2))
      *Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1)),0) + IF(gO2 <
      3,-0.0125*Conj(ZU(gI2,gO2))*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU
      (gI1,3 + j1)),0) + IF(gO2 < 3,-3*Conj(Yu(gO2,gO2))*KroneckerDelta(gO1,3 +
      gO2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*Yu(j1,j1)*ZU(gI1,3 + j1)),0) + IF(gO2 < 3,
      0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1,
      3 + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZU(gI2,gO2))
      *Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 <
      3,-0.0125*Conj(ZU(gI2,gO2))*Sqr(gN)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU
      (gI1,3 + j2)),0) - 0.13333333333333333*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*ZU(gI1,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,
      3 + j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      ZU(gI1,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))
      - 0.0125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1))*SUM(
      j2,0,2,Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1
      ,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.0375*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(
      gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) -
      0.4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.0375*Sqr(gN)*SUM(
      j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,
      j2)) - 0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) - SUM(j1,0,2,Conj(ZU(
      gI2,j1))*KroneckerDelta(gO2,3 + j1)*Yu(j1,j1))*SUM(j2,0,2,Conj(Yu(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*ZU(gI1,j2)) - 0.4*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*ZU(gI1,3 + j2)) - 0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3
       + j2)) - 0.13333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 +
      j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 +
      j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpAhSuconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,0.5)*vs*
      Conj(Yu(gO2,gO2))*Conj(ZU(gI1,3 + gO2))*Lambdax*ZA(gI2,0),0) + IF(gO2 < 3,
      std::complex<double>(0.,0.7071067811865475)*Conj(ZU(gI1,3 + gO2))*Conj(TYu(
      gO2,gO2))*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(0,0.5)*vd*Conj(Yu(
      gO2,gO2))*Conj(ZU(gI1,3 + gO2))*Lambdax*ZA(gI2,2),0) - std::complex<double>(
      0,0.5)*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*KroneckerDelta(gO2,3 +
      j1)*Yu(j1,j1))*ZA(gI2,0) - std::complex<double>(0.,0.7071067811865475)*SUM(
      j1,0,2,Conj(ZU(gI1,j1))*KroneckerDelta(gO2,3 + j1)*TYu(j1,j1))*ZA(gI2,1) -
      std::complex<double>(0,0.5)*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSuconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.5*vs*Conj(Yu(gO2,gO2))*Conj(ZU
      (gI1,3 + gO2))*Lambdax*ZH(gI2,0),0) + IF(gO2 < 3,0.05*vd*Conj(ZU(gI1,gO2))*
      Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 3,-0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,
      0),0) + IF(gO2 < 3,0.075*vd*Conj(ZU(gI1,gO2))*Sqr(gN)*ZH(gI2,0),0) + IF(gO2
      < 3,-(vu*AbsSqr(Yu(gO2,gO2))*Conj(ZU(gI1,gO2))*ZH(gI2,1)),0) + IF(gO2 < 3,-
      0.7071067811865475*Conj(ZU(gI1,3 + gO2))*Conj(TYu(gO2,gO2))*ZH(gI2,1),0) +
      IF(gO2 < 3,-0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gO2 < 3,0.25
      *vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF(gO2 < 3,0.05*vu*Conj(ZU(gI1,
      gO2))*Sqr(gN)*ZH(gI2,1),0) + IF(gO2 < 3,0.5*vd*Conj(Yu(gO2,gO2))*Conj(ZU(gI1
      ,3 + gO2))*Lambdax*ZH(gI2,2),0) + IF(gO2 < 3,-0.125*vs*Conj(ZU(gI1,gO2))*Sqr
      (gN)*ZH(gI2,2),0) - 0.2*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) + 0.075*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZU(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) + 0.5*vs*Conj(Lambdax)*
      SUM(j1,0,2,Conj(ZU(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Yu(j1,j1))*ZH(gI2,0)
      + 0.2*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))
      *ZH(gI2,1) + 0.05*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZH(gI2,1) - 0.7071067811865475*SUM(j1,0,2,Conj(ZU(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*TYu(j1,j1))*ZH(gI2,1) - vu*SUM(j2,0,2,AbsSqr(Yu(
      j2,j2))*Conj(ZU(gI1,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZH(gI2,1) - 0.125*
      vs*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(
      gI2,2) + 0.5*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*KroneckerDelta(gO2
      ,3 + j1)*Yu(j1,j1))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpGluFuconjUSuPR(int gI2, int gO2) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*Conj(PhaseGlu)*SUM(j1
      ,0,2,KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFuconjUSuPL(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*PhaseGlu*
      Conj(ZUL(gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSuconjVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZD(
      gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVG(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 6,g3*Conj(ZU(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.12909944487358055*g1*Conj(ZU(
      gI2,gO2))*Cos(ThetaW()),0) + IF(gO2 < 3,0.5*g2*Conj(ZU(gI2,gO2))*Sin(ThetaW(
      )),0) + 0.5163977794943222*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.5*g2*Conj(ZU(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZU(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,0.15811388300841897*gN*
      Conj(ZU(gI2,gO2))*Sin(ThetaWp()),0) - 0.5163977794943222*g1*Cos(ThetaWp())*
      Sin(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) -
      0.15811388300841897*gN*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.15811388300841897*gN*Conj(ZU(
      gI2,gO2))*Cos(ThetaWp()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZU(gI2,gO2))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO2 < 3,0.12909944487358055*g1*Conj(ZU(gI2,
      gO2))*Sin(ThetaW())*Sin(ThetaWp()),0) - 0.15811388300841897*gN*Cos(ThetaWp()
      )*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) +
      0.5163977794943222*g1*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZU(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.6324555320336759*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 3,0.4898979485566356*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,-0.7745966692414834*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(
      gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())),0) + IF(gO1 < 3,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())),0) + IF(gO1 < 3,0.2*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*Sqr(Sin(ThetaWp())),0) - 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin
      (ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 1.2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW
      ()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) +
      0.05*Sqr(gN)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.4898979485566356*g1*gN*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp()),0) + IF(gO1
      < 3,0.31622776601683794*g2*gN*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*Sin(2*
      ThetaWp()),0) + IF(gO1 < 3,0.2*KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(Cos(
      ThetaWp())),0) + IF(gO1 < 3,-0.7745966692414834*g1*g2*Cos(ThetaW())*
      KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0
      ) + IF(gO1 < 3,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp())),0) + 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1)) + 0.05*Sqr(gN)*Sqr(Cos(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) + 1.2*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      );

   return result;
}

double CLASSNAME::CpUSeconjUSeconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSeconjHpmconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr
      (gN)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*
      ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZP(
      gI1,1)*ZP(gI2,1),0) + 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) + 0.075*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) -
      SUM(j2,0,2,AbsSqr(Ye(j2,j2))*KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
       + j2))*ZP(gI1,0)*ZP(gI2,0) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) + 0.05*Sqr(gN)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUSeSHp0conjUSeconjSHp0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.15*Conj(UHp0(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,0),0) + IF(gO1 < 3,0.25*Conj(UHp0(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,0),0) + IF(gO1 < 3,-0.1*
      Conj(UHp0(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,0),0) + IF(gO1 <
      3,0.15*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,1),0) + IF
      (gO1 < 3,-0.25*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,1)
      ,0) + IF(gO1 < 3,0.1*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(
      gI2,1),0) + 0.3*Conj(UHp0(gI1,0))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,0) - 0.05*Conj(UHp0(gI1,0))*Sqr(gN)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,
      0) - 0.3*Conj(UHp0(gI1,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*UHp0(gI2,1) + 0.05*Conj(UHp0(gI1,1))*Sqr(gN)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHp0(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSeSHppconjUSeconjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.15*Conj(UHpp(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,0),0) + IF(gO1 < 3,-0.25*Conj(UHpp(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,0),0) + IF(gO1 < 3,-0.1*
      Conj(UHpp(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,0),0) + IF(gO1 <
      3,0.15*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,1),0) + IF
      (gO1 < 3,0.25*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,1),
      0) + IF(gO1 < 3,0.1*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(
      gI2,1),0) + 0.3*Conj(UHpp(gI1,0))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,0) - 0.05*Conj(UHpp(gI1,0))*Sqr(gN)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,
      0) - 0.3*Conj(UHpp(gI1,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*UHpp(gI2,1) + 0.05*Conj(UHpp(gI1,1))*Sqr(gN)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*UHpp(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSeSSI0conjUSeconjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.25*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN),0) - 0.125*KroneckerDelta(gI1,gI2)*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(AbsSqr(Ye(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZA(gI1,0)*ZA(gI2,0)),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,0.1*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.5*
      Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Ye(gO1,gO1)*ZA(gI1,2)*ZA(gI2,1),0)
      + IF(gO1 < 3,-0.5*Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Ye(gO1,gO1)*ZA(
      gI1,1)*ZA(gI2,2),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(
      gI1,2)*ZA(gI2,2),0) + IF(gO2 < 3,-0.5*Conj(Ye(gO2,gO2))*KroneckerDelta(gO1,3
       + gO2)*Lambdax*ZA(gI1,2)*ZA(gI2,1),0) + IF(gO2 < 3,-0.5*Conj(Ye(gO2,gO2))*
      KroneckerDelta(gO1,3 + gO2)*Lambdax*ZA(gI1,1)*ZA(gI2,2),0) + 0.3*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(
      gI2,0) + 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) - SUM(j2,0,2,AbsSqr(Ye(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZA(gI1,0)*ZA(gI2,0) -
      0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*ZA(gI1,1)*ZA(gI2,1) + 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1) - 0.125*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,2)*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(AbsSqr(Ye(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZH(gI1,0)*ZH(gI2,0)),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.1*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.5*Conj
      (Lambdax)*KroneckerDelta(3 + gO1,gO2)*Ye(gO1,gO1)*ZH(gI1,2)*ZH(gI2,1),0) +
      IF(gO1 < 3,0.5*Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)*Ye(gO1,gO1)*ZH(gI1,
      1)*ZH(gI2,2),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,2)
      *ZH(gI2,2),0) + IF(gO2 < 3,0.5*Conj(Ye(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)
      *Lambdax*ZH(gI1,2)*ZH(gI2,1),0) + IF(gO2 < 3,0.5*Conj(Ye(gO2,gO2))*
      KroneckerDelta(gO1,3 + gO2)*Lambdax*ZH(gI1,1)*ZH(gI2,2),0) + 0.3*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(
      gI2,0) + 0.075*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) - SUM(j2,0,2,AbsSqr(Ye(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZH(gI1,0)*ZH(gI2,0) -
      0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*ZH(gI1,1)*ZH(gI2,1) + 0.05*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1) - 0.125*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,2)*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUSeSvconjUSeconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZV(gI1,gO2)
      )*Sqr(g2)*ZV(gI2,gO1),0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gI1,gI2)
      *KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gI1,gI2
      )*KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.05*
      KroneckerDelta(gI1,gI2)*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) - SUM(j1,0,2,Conj(ZV(gI1,j1))*KroneckerDelta(gO2
      ,3 + j1)*Ye(j1,j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*KroneckerDelta(gO1,3 + j2)*ZV
      (gI2,j2));

   return result;
}

std::complex<double> CLASSNAME::CpHpmSvconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*vd*AbsSqr(Ye(
      gO2,gO2))*Conj(ZV(gI1,gO2))*ZP(gI2,0),0) + IF(gO2 < 3,-0.35355339059327373*
      vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,-0.35355339059327373*
      vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1),0) + SUM(j1,0,2,Conj(ZV(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*TYe(j1,j1))*ZP(gI2,0) + 0.7071067811865475*vs*
      Conj(Lambdax)*SUM(j1,0,2,Conj(ZV(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Ye(j1,
      j1))*ZP(gI2,1);

   return result;
}

double CLASSNAME::CpChaFvconjUSePR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpChaFvconjUSePL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-(g2*Conj(UM(gI2,0))*
      KroneckerDelta(gI1,gO1)),0) + IF(gI1 < 3,Conj(UM(gI2,1))*KroneckerDelta(3 +
      gI1,gO1)*Ye(gI1,gI1),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(Conj(Ye(gO2,gO2))*ZER(gI1,gO2)
      *ZN(gI2,2)),0) - 0.22360679774997896*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*
      ZER(gI1,j1))*(4.898979485566356*g1*ZN(gI2,0) + gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZEL(
      gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZEL(gI1
      ,gO1))*Conj(ZN(gI2,1)),0) + IF(gO1 < 3,-0.4472135954999579*gN*Conj(ZEL(gI1,
      gO1))*Conj(ZN(gI2,5)),0) - Conj(ZN(gI2,2))*SUM(j1,0,2,Conj(ZEL(gI1,j1))*
      KroneckerDelta(gO1,3 + j1)*Ye(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSHI0conjUSeconjSHI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),
      0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(
      gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2
       + j1)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj
      (UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(
      gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3
      ,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(
      gI2,j2)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,
      Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,
      gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + IF(gO1
       < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2)
      )*UHI0(gI2,2 + j2)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM
      (j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + 0.15*Sqr(g1)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(
      UHI0(gI1,j2))*UHI0(gI2,j2)) + 0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2
      )) - 0.15*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
       + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)) + 0.025*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1
      ,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)) + 0.15*Sqr(g1)*SUM(j1,0,1,Conj(
      UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.0375*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*
      UHI0(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.15*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.025*Sqr(gN
      )*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSHIpconjUSeconjSHIp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),
      0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(
      gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2
       + j1)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj
      (UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(
      gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3
      ,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(
      gI2,j2)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,
      Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,
      gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + IF(gO1
       < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))
      *UHIp(gI2,2 + j2)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(
      j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + 0.15*Sqr(g1)*SUM(j1,0,2
      ,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp
      (gI1,j2))*UHIp(gI2,j2)) + 0.0375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2))
      - 0.15*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)) + 0.025*Sqr(gN)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,
      Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)) + 0.15*Sqr(g1)*SUM(j1,0,1,Conj(UHIp
      (gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta
      (gO2,3 + j2)) + 0.0375*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.15*Sqr
      (g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.025*Sqr(gN)*SUM(
      j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSeconjSdconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,j1)
      )*ZD(gI2,j1)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,
      2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(
      gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)),0) + IF(gO1
       < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(
      gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 +
      j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(
      ZD(gI1,3 + j2))*ZD(gI2,3 + j2)),0) + IF(gO1 < 3,-(KroneckerDelta(3 + gO1,gO2
      )*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gI1,3 + j2))*ZD(gI2,j2))*Ye(gO1,gO1)),0
      ) + IF(gO2 < 3,-(Conj(Ye(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*SUM(j1,0,2,
      Conj(ZD(gI1,j1))*Yd(j1,j1)*ZD(gI2,3 + j1))),0) - 0.05*Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) -
      0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.025*Sqr(gN)*SUM(
      j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,
      j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) - 0.1*Sqr(g1)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI1,3 + j2))*ZD(gI2,3 + j2)) - 0.025*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3
      + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSeconjSDXconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)),0)
      + 0.0125*(80*IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj
      (ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)),0) + 80*IF(gO1 < 3,0.075*KroneckerDelta(
      gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)),0) + 80*
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI1,j2)
      )*ZDX(gI2,j2)),0) + 80*IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(
      j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)),0) + 80*IF(gO1 < 3,0.05*KroneckerDelta
      (gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 + j2)),0) + 80*
      IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI1,3 +
      j2))*ZDX(gI2,3 + j2)),0) + 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1
      ))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 2*Sqr
      (gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1
      ,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,3 +
      j1))*ZDX(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 8*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0
      ,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) + 2*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)
      ) - 8*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX
      (gI1,3 + j2))*ZDX(gI2,3 + j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSeconjSeconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.15*Conj(ZE(gI1,gO2
      ))*Sqr(g1)*ZE(gI2,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZE(gI1,gO2))
      *Sqr(g2)*ZE(gI2,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.1*Conj(ZE(gI1,gO2))*
      Sqr(gN)*ZE(gI2,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-(Conj(Ye(gO2,gO2))*Conj(
      ZE(gI1,3 + gO2))*Ye(gO1,gO1)*ZE(gI2,3 + gO1)),0),0) + IF(gO1 < 3,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1)
      )*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0
      ,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-
      0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2
      ,3 + j1)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(
      KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,2,Conj(Ye(j2,j2))*Conj(ZE(gI1,3 + j2))*
      ZE(gI2,j2))*Ye(gO1,gO1)),0) + IF(gO1 < 3,0.15*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1
      ,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZE(gI2,gO1),0) + IF(gO1 < 3,-0.025*Sqr
      (gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZE(gI2,gO1)
      ,0) + IF(gO1 < 3,0.15*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*KroneckerDelta
      (gO2,3 + j2))*ZE(gI2,gO1),0) + IF(gO1 < 3,-0.025*Sqr(gN)*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZE(gI2,gO1),0) + IF(gO2 < 3,0.15*
      Conj(ZE(gI1,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 +
      j1)),0) + IF(gO2 < 3,-0.025*Conj(ZE(gI1,gO2))*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1)),0) + IF(gO2 < 3,-(Conj(Ye(gO2,gO2
      ))*KroneckerDelta(gO1,3 + gO2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*Ye(j1,j1)*ZE(gI2,
      3 + j1))),0) + IF(gO2 < 3,0.15*Conj(ZE(gI1,gO2))*Sqr(g1)*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)),0) + IF(gO2 < 3,-0.025*Conj(ZE(
      gI1,gO2))*Sqr(gN)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)),0) -
      0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1))*SUM(j2,0,2
      ,Conj(ZE(gI1,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))
      *KroneckerDelta(gO2,3 + j2)) + 0.15*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(
      gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) -
      0.025*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.3*Sqr(g1)*SUM(j1,
      0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 +
      j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) + 0.15*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 0.025*
      Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - SUM(j1,0,2,Conj(ZE(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j1))*SUM(j2,0,2,Conj(Ye(j2,j2))*
      KroneckerDelta(gO1,3 + j2)*ZE(gI2,j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3
       + j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)) - 0.0125*Sqr(gN)*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSuconjUSeconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1)),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,
      2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta
      (gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(
      gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU
      (gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,-0.1*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 +
      j2)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(
      ZU(gI1,3 + j2))*ZU(gI2,3 + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,
      0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 +
      j1))*ZU(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)) - 0.0125
      *Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*
      SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI1,3 + j2))*ZU(gI2,3 + j2)) - 0.0125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3
       + j2));

   return result;
}

std::complex<double> CLASSNAME::CpAhSeconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZE(gI1,3 + gO2))*Conj(TYe(gO2,gO2))*ZA(gI2,0),0) +
      IF(gO2 < 3,std::complex<double>(0,0.5)*vs*Conj(Ye(gO2,gO2))*Conj(ZE(gI1,3 +
      gO2))*Lambdax*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(0,0.5)*vu*Conj(
      Ye(gO2,gO2))*Conj(ZE(gI1,3 + gO2))*Lambdax*ZA(gI2,2),0) - std::complex<
      double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(ZE(gI1,j1))*KroneckerDelta(
      gO2,3 + j1)*TYe(j1,j1))*ZA(gI2,0) - std::complex<double>(0,0.5)*vs*Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZE(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Ye(j1,j1))*
      ZA(gI2,1) - std::complex<double>(0,0.5)*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(
      gI1,j1))*KroneckerDelta(gO2,3 + j1)*Ye(j1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSeconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(vd*AbsSqr(Ye(gO2,gO2))*Conj(ZE
      (gI1,gO2))*ZH(gI2,0)),0) + IF(gO2 < 3,-0.7071067811865475*Conj(ZE(gI1,3 +
      gO2))*Conj(TYe(gO2,gO2))*ZH(gI2,0),0) + IF(gO2 < 3,-0.15*vd*Conj(ZE(gI1,gO2)
      )*Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 3,0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2
      ,0),0) + IF(gO2 < 3,0.15*vd*Conj(ZE(gI1,gO2))*Sqr(gN)*ZH(gI2,0),0) + IF(gO2
      < 3,0.5*vs*Conj(Ye(gO2,gO2))*Conj(ZE(gI1,3 + gO2))*Lambdax*ZH(gI2,1),0) + IF
      (gO2 < 3,0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gO2 < 3,-0.25*
      vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF(gO2 < 3,0.1*vu*Conj(ZE(gI1,
      gO2))*Sqr(gN)*ZH(gI2,1),0) + IF(gO2 < 3,0.5*vu*Conj(Ye(gO2,gO2))*Conj(ZE(gI1
      ,3 + gO2))*Lambdax*ZH(gI2,2),0) + IF(gO2 < 3,-0.25*vs*Conj(ZE(gI1,gO2))*Sqr(
      gN)*ZH(gI2,2),0) + 0.3*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) + 0.075*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZE(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) - 0.7071067811865475*SUM(
      j1,0,2,Conj(ZE(gI1,j1))*KroneckerDelta(gO2,3 + j1)*TYe(j1,j1))*ZH(gI2,0) -
      vd*SUM(j2,0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(gI1,3 + j2))*KroneckerDelta(gO2,3 +
      j2))*ZH(gI2,0) - 0.3*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) + 0.05*vu*Sqr(gN)*SUM(j1,0,2,Conj(ZE(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) + 0.5*vs*Conj(Lambdax)*
      SUM(j1,0,2,Conj(ZE(gI1,j1))*KroneckerDelta(gO2,3 + j1)*Ye(j1,j1))*ZH(gI2,1)
      - 0.125*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1
      ))*ZH(gI2,2) + 0.5*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j1))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSeVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZV(
      gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSeVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.3872983346207417*g1*Conj(ZE(
      gI2,gO2))*Cos(ThetaW()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZE(gI2,gO2))*Sin(ThetaW
      ()),0) - 0.7745966692414834*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))
      *KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSeVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Conj(ZE(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,0.3872983346207417*g1*Conj(ZE(gI2,
      gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,0.31622776601683794*gN*
      Conj(ZE(gI2,gO2))*Sin(ThetaWp()),0) + 0.7745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) -
      0.15811388300841897*gN*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSeVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.31622776601683794*gN*Conj(ZE(
      gI2,gO2))*Cos(ThetaWp()),0) + IF(gO2 < 3,0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW
      ())*Sin(ThetaWp()),0) + IF(gO2 < 3,-0.3872983346207417*g1*Conj(ZE(gI2,gO2))*
      Sin(ThetaW())*Sin(ThetaWp()),0) - 0.15811388300841897*gN*Cos(ThetaWp())*SUM(
      j1,0,2,Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) - 0.7745966692414834
      *g1*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXconjUSDXVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.32659863237109044*g1*gN*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp()),0) + IF(gO1
      < 3,0.13333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*
      Sqr(Sin(ThetaW())),0) + IF(gO1 < 3,0.2*KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(
      Sin(ThetaWp())),0) + 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1)) + 0.13333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.45*Sqr
      (gN)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXconjUSDXVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.32659863237109044*g1*gN*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp()),0) + IF(gO1
      < 3,0.2*KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,
      0.13333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp())),0) - 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1)) + 0.45*Sqr(gN)*Sqr(Cos(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) + 0.13333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()
      ))*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSDXconjHpmconjUSDX(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZP(
      gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(Conj(Lambdax)*KroneckerDelta(3 + gO1,gO2)
      *ZP(gI1,1)*ZP(gI2,0)*Kappa(gO1,gO1)),0) + IF(gO2 < 3,-(Conj(Kappa(gO2,gO2))*
      KroneckerDelta(gO1,3 + gO2)*Lambdax*ZP(gI1,0)*ZP(gI2,1)),0) + 0.1*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*
      ZP(gI2,0) - 0.225*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) -
      0.15*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1
      ))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSHp0conjUSDXconjSHp0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.1*Conj(UHp0(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,0),0) + IF(gO1 < 3,0.1*Conj(UHp0(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,0),0) + IF(gO1 < 3,0.1*Conj
      (UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,1),0) + IF(gO1 < 3,-
      0.1*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,1),0) + 0.05*
      (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*(Conj(UHp0(gI1,0))*UHp0(gI2,0) - Conj(UHp0(gI1,1))*UHp0(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSHppconjUSDXconjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.1*Conj(UHpp(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,0),0) + IF(gO1 < 3,0.1*Conj(UHpp(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,0),0) + IF(gO1 < 3,0.1*Conj
      (UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,1),0) + IF(gO1 < 3,-
      0.1*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,1),0) + 0.05*
      (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*(Conj(UHpp(gI1,0))*UHpp(gI2,0) - Conj(UHpp(gI1,1))*UHpp(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSSI0conjUSDXconjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.25*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.375*KroneckerDelta(gI1,gI2)*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSDXconjUSDX(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(
      gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-(AbsSqr(Kappa(gO1,gO1))*KroneckerDelta(gO1
      ,gO2)*ZA(gI1,2)*ZA(gI2,2)),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZA(gI1,2)*ZA(gI2,2),0) + IF(gO1 < 3,-0.5*Conj(Lambdax)*KroneckerDelta(3
      + gO1,gO2)*ZA(gI1,1)*ZA(gI2,0)*Kappa(gO1,gO1),0) + IF(gO1 < 3,-0.5*Conj(
      Lambdax)*KroneckerDelta(3 + gO1,gO2)*ZA(gI1,0)*ZA(gI2,1)*Kappa(gO1,gO1),0) +
      IF(gO2 < 3,-0.5*Conj(Kappa(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*Lambdax*ZA(
      gI1,1)*ZA(gI2,0),0) + IF(gO2 < 3,-0.5*Conj(Kappa(gO2,gO2))*KroneckerDelta(
      gO1,3 + gO2)*Lambdax*ZA(gI1,0)*ZA(gI2,1),0) + 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) -
      0.225*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*ZA(gI1,0)*ZA(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)
      *KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1) - 0.15*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1) +
      0.375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*ZA(gI1,2)*ZA(gI2,2) - SUM(j2,0,2,AbsSqr(Kappa(j2,j2))*KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2))*ZA(gI1,2)*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSDXconjUSDX(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(
      gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-(AbsSqr(Kappa(gO1,gO1))*KroneckerDelta(gO1
      ,gO2)*ZH(gI1,2)*ZH(gI2,2)),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZH(gI1,2)*ZH(gI2,2),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*KroneckerDelta(3 +
      gO1,gO2)*ZH(gI1,1)*ZH(gI2,0)*Kappa(gO1,gO1),0) + IF(gO1 < 3,0.5*Conj(Lambdax
      )*KroneckerDelta(3 + gO1,gO2)*ZH(gI1,0)*ZH(gI2,1)*Kappa(gO1,gO1),0) + IF(gO2
       < 3,0.5*Conj(Kappa(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*Lambdax*ZH(gI1,1)*
      ZH(gI2,0),0) + IF(gO2 < 3,0.5*Conj(Kappa(gO2,gO2))*KroneckerDelta(gO1,3 +
      gO2)*Lambdax*ZH(gI1,0)*ZH(gI2,1),0) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) - 0.225*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*
      ZH(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1) - 0.15*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1) + 0.375*Sqr(gN)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,2)*
      ZH(gI2,2) - SUM(j2,0,2,AbsSqr(Kappa(j2,j2))*KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2))*ZH(gI1,2)*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSvconjUSDXconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.1*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,0.1*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.05*KroneckerDelta(gI1,gI2)*(2*Sqr(g1)
      + 3*Sqr(gN))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpChiFDXconjUSDXPR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(Conj(Kappa(gO2,gO2))*ZDXR(gI1,
      gO2)*ZN(gI2,4)),0) + 0.07453559924999298*SUM(j1,0,2,KroneckerDelta(gO2,3 +
      j1)*ZDXR(gI1,j1))*(-4.898979485566356*g1*ZN(gI2,0) + 9*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiFDXconjUSDXPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.3651483716701107*g1*Conj(ZDXL(
      gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,0.4472135954999579*gN*Conj(ZDXL(
      gI1,gO1))*Conj(ZN(gI2,5)),0) - Conj(ZN(gI2,4))*SUM(j1,0,2,Conj(ZDXL(gI1,j1))
      *KroneckerDelta(gO1,3 + j1)*Kappa(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSHI0conjUSDXconjSHI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 3,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),
      0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(
      gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)
      *Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 3
      ,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(
      gI2,j2)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,
      Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + IF(gO1 <
      3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*
      UHI0(gI2,2 + j2)),0) + IF(gO1 < 3,KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,1,
      Conj(UHI0(gI1,2 + j2))*Conj(Lambda12(j2,j2))*UHI0(gI2,j2))*Kappa(gO1,gO1),0)
      + IF(gO2 < 3,Conj(Kappa(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*SUM(j1,0,1,
      Conj(UHI0(gI1,j1))*UHI0(gI2,2 + j1)*Lambda12(j1,j1)),0) + 0.05*Sqr(g1)*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,
      Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) - 0.1125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(
      gI2,j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)) - 0.075*
      Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*
      SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,
      1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.1125*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*
      UHI0(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.05*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.075*Sqr(gN
      )*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSHIpconjUSDXconjSHIp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 3,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),
      0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(
      gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)
      *Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 3
      ,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(
      gI2,j2)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,
      Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + IF(gO1 <
      3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*
      UHIp(gI2,2 + j2)),0) + IF(gO1 < 3,-(KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,1,
      Conj(UHIp(gI1,2 + j2))*Conj(Lambda12(j2,j2))*UHIp(gI2,j2))*Kappa(gO1,gO1)),0
      ) + IF(gO2 < 3,-(Conj(Kappa(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*SUM(j1,0,1
      ,Conj(UHIp(gI1,j1))*UHIp(gI2,2 + j1)*Lambda12(j1,j1))),0) + 0.05*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,
      Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) - 0.1125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(
      gI2,j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)) - 0.075*
      Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*
      SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,
      1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.1125*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*
      UHIp(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.05*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(
      j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.075*Sqr(gN
      )*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXconjUSDXconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1
      ))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) +
      IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI2,j2))
      *ZD(gI1,j2)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,
      Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1
      ,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) - 0.05*Sqr(
      g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1125*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,j1
      ))*ZD(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,
      2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.225*Sqr(gN)*SUM
      (j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
       + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,
      j2)) + 0.1125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) - 0.1*Sqr(g1)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*ZD(gI1,3 + j2)) + 0.225*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3
      + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXconjUSDXconjSDXSDX(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.06666666666666667*
      Conj(ZDX(gI2,gO2))*Sqr(g1)*ZDX(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-
      1.3333333333333333*Conj(ZDX(gI2,gO2))*Sqr(g3)*ZDX(gI1,gO1),0),0) + IF(gO1 <
      3,IF(gO2 < 3,-0.1*Conj(ZDX(gI2,gO2))*Sqr(gN)*ZDX(gI1,gO1),0),0) + IF(gO1 < 3
      ,IF(gO2 < 3,-(Conj(ZDX(gI2,3 + gO2))*Conj(Kappa(gO2,gO2))*ZDX(gI1,3 + gO1)*
      Kappa(gO1,gO1)),0),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)),0) + IF(gO1 < 3,-0.15*KroneckerDelta(
      gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)),0) + IF(gO1 < 3,
      0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1
      ,3 + j1)),0) + IF(gO1 < 3,-0.225*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,
      Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)),0) + IF(gO1 < 3,-0.1*KroneckerDelta(
      gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2)),0) + IF(gO1 < 3,-
      0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2
      )),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(
      gI2,3 + j2))*ZDX(gI1,3 + j2)),0) + IF(gO1 < 3,-0.225*KroneckerDelta(gO1,gO2)
      *Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)),0) + IF(gO1 < 3,
      0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZDX(gI1,gO1),0) + IF(gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j1,
      0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZDX(gI1,gO1),0) + IF(
      gO1 < 3,-0.075*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(gO2,3
       + j1))*ZDX(gI1,gO1),0) + IF(gO1 < 3,0.03333333333333333*Sqr(g1)*SUM(j2,0,2,
      Conj(ZDX(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZDX(gI1,gO1),0) + IF(gO1 <
      3,0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*KroneckerDelta
      (gO2,3 + j2))*ZDX(gI1,gO1),0) + IF(gO1 < 3,-0.075*Sqr(gN)*SUM(j2,0,2,Conj(
      ZDX(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZDX(gI1,gO1),0) + IF(gO1 < 3,-3
      *KroneckerDelta(3 + gO1,gO2)*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*Conj(Kappa(j2,
      j2))*ZDX(gI1,j2))*Kappa(gO1,gO1),0) + IF(gO2 < 3,0.03333333333333333*Conj(
      ZDX(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZDX(gI1,3 + j1))
      ,0) + IF(gO2 < 3,0.6666666666666666*Conj(ZDX(gI2,gO2))*Sqr(g3)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZDX(gI1,3 + j1)),0) + IF(gO2 < 3,-0.075*Conj(ZDX(
      gI2,gO2))*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZDX(gI1,3 + j1)),0)
      + IF(gO2 < 3,-3*Conj(Kappa(gO2,gO2))*KroneckerDelta(gO1,3 + gO2)*SUM(j1,0,2,
      Conj(ZDX(gI2,j1))*ZDX(gI1,3 + j1)*Kappa(j1,j1)),0) + IF(gO2 < 3,
      0.03333333333333333*Conj(ZDX(gI2,gO2))*Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1
      ,3 + j2)*ZDX(gI1,3 + j2)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZDX(gI2,
      gO2))*Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZDX(gI1,3 + j2)),0) + IF
      (gO2 < 3,-0.075*Conj(ZDX(gI2,gO2))*Sqr(gN)*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*ZDX(gI1,3 + j2)),0) - 0.03333333333333333*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZDX(gI1,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))
      *KroneckerDelta(gO2,3 + j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZDX(gI1,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))
      *KroneckerDelta(gO2,3 + j2)) - 0.1125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*ZDX(gI1,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*KroneckerDelta(gO2
      ,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1))*SUM(j2,0,2
      ,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.225*Sqr(gN)*SUM(
      j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*
      ZDX(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
      + j2)) - 0.3375*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2
      ,0,2,Conj(ZDX(gI2,j2))*ZDX(gI1,j2)) - 0.225*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(
      gI2,j2))*ZDX(gI1,j2)) - SUM(j1,0,2,Conj(ZDX(gI2,j1))*KroneckerDelta(gO2,3 +
      j1)*Kappa(j1,j1))*SUM(j2,0,2,Conj(Kappa(j2,j2))*KroneckerDelta(gO1,3 + j2)*
      ZDX(gI1,j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)
      ) - 0.3375*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,
      3 + j1))*SUM(j2,0,2,Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2)) -
      0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZDX(gI1,3 + j2)) -
      0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZDX(gI1,3 + j2)) - 0.1125
      *Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*SUM(j2
      ,0,2,KroneckerDelta(gO1,3 + j2)*ZDX(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXconjUSDXconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) +
      IF(gO1 < 3,-0.2*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 +
      j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) +
      IF(gO1 < 3,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI2,j2))
      *ZU(gI1,j2)),0) + IF(gO1 < 3,-0.2*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2
      ,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,0.075*KroneckerDelta(
      gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) - 0.05*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1125*Sqr(gN)*SUM(j1,0,2,Conj(ZU(
      gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1125*
      Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1
      ,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(
      ZU(gI2,j2))*ZU(gI1,j2)) + 0.1125*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) +
      0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)) + 0.1125*Sqr(gN)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*ZU(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSDXSeconjUSDXconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      0.0125*(80*IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE
      (gI1,3 + j1))*ZE(gI2,3 + j1)),0) + 80*IF(gO1 < 3,0.025*KroneckerDelta(gO1,
      gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + 80*IF(gO1
      < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2
      ,j2)),0) + 80*IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + 80*IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + 80*IF(gO1 < 3,
      0.025*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2
      ,3 + j2)),0) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 6*Sqr(gN)*SUM(j1,0,
      2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 8*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(
      gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2
      )) + 3*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 4*Sqr(g1)*SUM(j1,0,
      2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*ZE(gI2,j2)) + 6*Sqr(gN)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 8*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2
      ,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)) + 3*Sqr(gN)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)));

   return result;
}

std::complex<double> CLASSNAME::CpAhSDXconjUSDX(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,0.5)*vu*
      Conj(ZDX(gI1,3 + gO2))*Conj(Kappa(gO2,gO2))*Lambdax*ZA(gI2,0),0) + IF(gO2 <
      3,std::complex<double>(0,0.5)*vd*Conj(ZDX(gI1,3 + gO2))*Conj(Kappa(gO2,gO2))
      *Lambdax*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZDX(gI1,3 + gO2))*Conj(TKappa(gO2,gO2))*ZA(gI2,2),0
      ) - std::complex<double>(0,0.5)*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gI1,j1)
      )*KroneckerDelta(gO2,3 + j1)*Kappa(j1,j1))*ZA(gI2,0) - std::complex<double>(
      0,0.5)*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*KroneckerDelta(gO2,3 +
      j1)*Kappa(j1,j1))*ZA(gI2,1) - std::complex<double>(0.,0.7071067811865475)*
      SUM(j1,0,2,Conj(ZDX(gI1,j1))*KroneckerDelta(gO2,3 + j1)*TKappa(j1,j1))*ZA(
      gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSDXconjUSDX(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.5*vu*Conj(ZDX(gI1,3 + gO2))*
      Conj(Kappa(gO2,gO2))*Lambdax*ZH(gI2,0),0) + IF(gO2 < 3,-0.1*vd*Conj(ZDX(gI1,
      gO2))*Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 3,-0.15*vd*Conj(ZDX(gI1,gO2))*Sqr(gN)*
      ZH(gI2,0),0) + IF(gO2 < 3,0.5*vd*Conj(ZDX(gI1,3 + gO2))*Conj(Kappa(gO2,gO2))
      *Lambdax*ZH(gI2,1),0) + IF(gO2 < 3,0.1*vu*Conj(ZDX(gI1,gO2))*Sqr(g1)*ZH(gI2,
      1),0) + IF(gO2 < 3,-0.1*vu*Conj(ZDX(gI1,gO2))*Sqr(gN)*ZH(gI2,1),0) + IF(gO2
      < 3,-(vs*AbsSqr(Kappa(gO2,gO2))*Conj(ZDX(gI1,gO2))*ZH(gI2,2)),0) + IF(gO2 <
      3,-0.7071067811865475*Conj(ZDX(gI1,3 + gO2))*Conj(TKappa(gO2,gO2))*ZH(gI2,2)
      ,0) + IF(gO2 < 3,0.25*vs*Conj(ZDX(gI1,gO2))*Sqr(gN)*ZH(gI2,2),0) + 0.1*vd*
      Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,
      0) - 0.225*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*KroneckerDelta(gO2,3
      + j1))*ZH(gI2,0) + 0.5*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      KroneckerDelta(gO2,3 + j1)*Kappa(j1,j1))*ZH(gI2,0) - 0.1*vu*Sqr(g1)*SUM(j1,0
      ,2,Conj(ZDX(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) - 0.15*vu*Sqr
      (gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1)
      + 0.5*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*KroneckerDelta(gO2,3 +
      j1)*Kappa(j1,j1))*ZH(gI2,1) + 0.375*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 +
      j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,2) - 0.7071067811865475*SUM(j1,0,2,
      Conj(ZDX(gI1,j1))*KroneckerDelta(gO2,3 + j1)*TKappa(j1,j1))*ZH(gI2,2) - vs*
      SUM(j2,0,2,AbsSqr(Kappa(j2,j2))*Conj(ZDX(gI1,3 + j2))*KroneckerDelta(gO2,3 +
      j2))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpGluFDXconjUSDXPR(int gI2, int gO2) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*Conj(PhaseGlu)*SUM(j1
      ,0,2,KroneckerDelta(gO2,3 + j1)*ZDXR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFDXconjUSDXPL(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*PhaseGlu*
      Conj(ZDXL(gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjUSDXVG(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 6,g3*Conj(ZDX(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjUSDXVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 6,-0.2581988897471611*g1*Conj(ZDX(
      gI2,gO2))*Cos(ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjUSDXVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.2581988897471611*g1*Conj(ZDX(
      gI2,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,-0.31622776601683794*
      gN*Conj(ZDX(gI2,gO2))*Sin(ThetaWp()),0) + 0.03726779962499649*(
      6.928203230275509*g1*Cos(ThetaWp())*Sin(ThetaW()) + 12.727922061357857*gN*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjUSDXVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.31622776601683794*gN*Conj(ZDX
      (gI2,gO2))*Cos(ThetaWp()),0) + IF(gO2 < 3,-0.2581988897471611*g1*Conj(ZDX(
      gI2,gO2))*Sin(ThetaW())*Sin(ThetaWp()),0) + 0.03726779962499649*(
      12.727922061357857*gN*Cos(ThetaWp()) - 6.928203230275509*g1*Sin(ThetaW())*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgWmUhh(int gO1) const
{
   
   const std::complex<double> result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgWmCUhh(int gO1) const
{
   
   const std::complex<double> result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargZgZUhh(int gO1) const
{
   
   const std::complex<double> result = 0.025*(-25*vs*KroneckerDelta(2,gO1)*Sqr(gN)
      *Sqr(Sin(ThetaWp())) - vd*KroneckerDelta(0,gO1)*(-14.696938456699067*g1*gN*
      Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*
      Sqr(Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp()
      )*Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp())
      )) + 6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(
      ThetaWp()))) - 2*vu*KroneckerDelta(1,gO1)*(3.1622776601683795*g2*gN*Cos(
      ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(
      2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgZUhh(int gO1) const
{
   
   const std::complex<double> result = 0.025*(-25*vs*Cos(ThetaWp())*KroneckerDelta
      (2,gO1)*Sin(ThetaWp())*Sqr(gN) + vd*KroneckerDelta(0,gO1)*(9.486832980505138
      *g2*gN*Cos(ThetaW())*Cos(2*ThetaWp()) - 9*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(
      gN) + 5*Sin(2*ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) + 7.348469228349534*g1*
      gN*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + g1*(3.872983346207417*g2*Sin(2*ThetaW
      ())*Sin(2*ThetaWp()) + 3*g1*Sin(2*ThetaWp())*Sqr(Sin(ThetaW())) -
      7.348469228349534*gN*Sin(ThetaW())*Sqr(Sin(ThetaWp())))) + vu*KroneckerDelta
      (1,gO1)*(-6.324555320336759*g2*gN*Cos(ThetaW())*Cos(2*ThetaWp()) - 4*Cos(
      ThetaWp())*Sin(ThetaWp())*Sqr(gN) + 5*Sin(2*ThetaWp())*Sqr(g2)*Sqr(Cos(
      ThetaW())) - 4.898979485566356*g1*gN*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + g1*
      (3.872983346207417*g2*Sin(2*ThetaW())*Sin(2*ThetaWp()) + 3*g1*Sin(2*ThetaWp(
      ))*Sqr(Sin(ThetaW())) + 4.898979485566356*gN*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      )))));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgZpUhh(int gO1) const
{
   
   const std::complex<double> result = 0.025*(-25*vs*KroneckerDelta(2,gO1)*Sqr(gN)
      *Sqr(Cos(ThetaWp())) - vd*KroneckerDelta(0,gO1)*(9*Sqr(gN)*Sqr(Cos(ThetaWp()
      )) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*
      (2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp(
      )))) + Cos(ThetaW())*(9.486832980505138*g2*gN*Sin(2*ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Sin(ThetaWp())))) - 2*vu*
      KroneckerDelta(1,gO1)*(-3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp
      ()) + 2*Sqr(gN)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp())) + g1*(-2.449489742783178*gN*Sin(ThetaW())*Sin(2*ThetaWp()) +
      3.872983346207417*g2*Sin(2*ThetaW())*Sqr(Sin(ThetaWp())) + 3*g1*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*vs*KroneckerDelta(2,gO2)*Sqr(gN)*
      Sqr(Sin(ThetaWp())) + vd*KroneckerDelta(0,gO2)*(-14.696938456699067*g1*gN*
      Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*
      Sqr(Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp()
      )*Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp())
      )) + 6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(
      ThetaWp()))) + 2*vu*KroneckerDelta(1,gO2)*(3.1622776601683795*g2*gN*Cos(
      ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(
      2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZp(int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*vs*Cos(ThetaWp())*KroneckerDelta(2
      ,gO2)*Sin(ThetaWp())*Sqr(gN) - vd*KroneckerDelta(0,gO2)*(9.486832980505138*
      g2*gN*Cos(ThetaW())*Cos(2*ThetaWp()) - 9*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(
      gN) + 5*Sin(2*ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) + 7.348469228349534*g1*
      gN*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + g1*(3.872983346207417*g2*Sin(2*ThetaW
      ())*Sin(2*ThetaWp()) + 3*g1*Sin(2*ThetaWp())*Sqr(Sin(ThetaW())) -
      7.348469228349534*gN*Sin(ThetaW())*Sqr(Sin(ThetaWp())))) + vu*KroneckerDelta
      (1,gO2)*(6.324555320336759*g2*gN*Cos(ThetaW())*Cos(2*ThetaWp()) + 4*Cos(
      ThetaWp())*Sin(ThetaWp())*Sqr(gN) - 5*Sin(2*ThetaWp())*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 4.898979485566356*g1*gN*Sin(ThetaW())*Sqr(Cos(ThetaWp())) - g1*
      (3.872983346207417*g2*Sin(2*ThetaW())*Sin(2*ThetaWp()) + 3*g1*Sin(2*ThetaWp(
      ))*Sqr(Sin(ThetaW())) + 4.898979485566356*gN*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      )))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZpVZp(int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*vs*KroneckerDelta(2,gO2)*Sqr(gN)*
      Sqr(Cos(ThetaWp())) + vd*KroneckerDelta(0,gO2)*(9*Sqr(gN)*Sqr(Cos(ThetaWp())
      ) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(
      2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      ))) + Cos(ThetaW())*(9.486832980505138*g2*gN*Sin(2*ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Sin(ThetaWp())))) + 2*vu*
      KroneckerDelta(1,gO2)*(-3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp
      ()) + 2*Sqr(gN)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp())) + g1*(-2.449489742783178*gN*Sin(ThetaW())*Sin(2*ThetaWp()) +
      3.872983346207417*g2*Sin(2*ThetaW())*Sqr(Sin(ThetaWp())) + 3*g1*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(int gO2) const
{
   
   const std::complex<double> result = 0.5*(vd*KroneckerDelta(0,gO2) + vu*
      KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gN)*Sqr(Sin(ThetaWp())) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(-14.696938456699067*g1*gN*Cos(ThetaWp())*Sin(ThetaW()
      )*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + Cos(
      ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp())*Sin(ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 6*Sqr(g1)*Sqr(
      Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(ThetaWp()))) + 2*
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3.1622776601683795*g2*gN*Cos(
      ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(
      2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gN)*Sqr(Cos(ThetaWp())) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(9*Sqr(gN)*Sqr(Cos(ThetaWp())) + 10*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(2.449489742783178*gN*
      Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))) + Cos(ThetaW())*(
      9.486832980505138*g2*gN*Sin(2*ThetaWp()) + 15.491933384829668*g1*g2*Sin(
      ThetaW())*Sqr(Sin(ThetaWp())))) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(-3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 2*Sqr(gN)*
      Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*
      (-2.449489742783178*gN*Sin(ThetaW())*Sin(2*ThetaWp()) + 3.872983346207417*g2
      *Sin(2*ThetaW())*Sqr(Sin(ThetaWp())) + 3*g1*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0)
      + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP(gI1,1)*ZP(gI2,1)) + 2*KroneckerDelta(1
      ,gO1)*(-5*KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP
      (gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(
      g2) - 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*
      ZP(gI1,1)*ZP(gI2,1))) - KroneckerDelta(0,gO1)*(10*KroneckerDelta(1,gO2)*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gI1,0)*ZP(gI2
      ,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSHp0conjSHp0(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN)) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gN) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*
      (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN)))*(Conj(UHp0(gI1,0))*UHp0(gI2,0) - Conj(
      UHp0(gI1,1))*UHp0(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSHppconjSHpp(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN)) + 5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gN))*(Conj(UHpp(gI1,0))*UHpp
      (gI2,0) - Conj(UHpp(gI1,1))*UHpp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSSI0conjSSI0(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.125*(3*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) - 5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2))*KroneckerDelta(gI1,gI2)*Sqr(gN)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(5*KroneckerDelta(2,gO2)*(ZP(gI1,0)*(
      vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gI2,0) - 5.656854249492381*Conj(
      TLambdax)*ZP(gI2,1)) + 2*ZP(gI1,1)*(-2.8284271247461903*TLambdax*ZP(gI2,0) +
      vs*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP(gI2,1))) + 2*KroneckerDelta(1,gO2)*(ZP(
      gI1,0)*(vu*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZP(gI2,0) - 5*vd*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gI2,1)) - ZP(gI1,1)*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(
      g2))*ZP(gI2,0) + vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gI2,1))) -
      KroneckerDelta(0,gO2)*(ZP(gI1,0)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP
      (gI2,0) + 10*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,1)) + 2*ZP(gI1,1)*(5*
      vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 3
      *Sqr(gN))*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSHp0conjSHp0(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.05*(vd*KroneckerDelta(0,gO2)*(3*Sqr(g1)
      + 5*Sqr(g2) - 3*Sqr(gN)) + 5*vs*KroneckerDelta(2,gO2)*Sqr(gN) - vu*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN)))*(Conj(UHp0(gI2,0)
      )*UHp0(gI1,0) - Conj(UHp0(gI2,1))*UHp0(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSHppconjSHpp(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.05*(vd*KroneckerDelta(0,gO2)*(3*Sqr(g1)
      - 5*Sqr(g2) - 3*Sqr(gN)) + vu*KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)
      - 2*Sqr(gN)) + 5*vs*KroneckerDelta(2,gO2)*Sqr(gN))*(Conj(UHpp(gI2,0))*UHpp(
      gI1,0) - Conj(UHpp(gI2,1))*UHpp(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSSI0conjSSI0(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.125*(3*vd*KroneckerDelta(0,gO2) + 2*vu*
      KroneckerDelta(1,gO2) - 5*vs*KroneckerDelta(2,gO2))*KroneckerDelta(gI1,gI2)*
      Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*(g2*KroneckerDelta(0,
      gO2)*UM(gI1,1)*UP(gI2,0) + (g2*KroneckerDelta(1,gO2)*UM(gI1,0) + Conj(
      Lambdax)*KroneckerDelta(2,gO2)*UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*
      Conj(UP(gI1,1))*KroneckerDelta(1,gO1) + Conj(UM(gI2,1))*(g2*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) + Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIChaIUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(2,gO2)*
      SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZMI(gI1,j1)*ZPI(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIChaIUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(2,gO1)*
      SUM(j1,0,1,Conj(ZMI(gI2,j1))*Conj(ZPI(gI1,j1))*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZA(gI1,0)*ZA(gI2
      ,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZA(gI1,1)*
      ZA(gI2,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZA(gI1,2)*ZA(gI2,2))) + 5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN)
      )*ZA(gI1,0)*ZA(gI2,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gI1,1)*ZA(gI2,1)
      - 5*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2)) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gI1,0)*ZA
      (gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZA(gI1,1)*ZA(gI2,1) + 5*(-4*
      AbsSqr(Lambdax) + Sqr(gN))*ZA(gI1,2)*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*(-2*
      KroneckerDelta(1,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(
      gN))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) - 5*KroneckerDelta(2,gO2)*(
      -8*AbsSqr(Lambdax) + 3*Sqr(gN))*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2))
      + KroneckerDelta(0,gO2)*(3*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gI1,0)*ZH
      (gI2,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gI1,
      1)*ZH(gI2,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gI1,2)*ZH(gI2,2)))) + 5*
      KroneckerDelta(2,gO1)*(KroneckerDelta(0,gO2)*(-8*AbsSqr(Lambdax) + 3*Sqr(gN)
      )*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + 2*KroneckerDelta(1,gO2)*(-4*
      AbsSqr(Lambdax) + Sqr(gN))*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) +
      KroneckerDelta(2,gO2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gI1,0)*ZH(gI2,0)
      + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gI1,1)*ZH(gI2,1) - 15*Sqr(gN)*ZH(gI1,2
      )*ZH(gI2,2))) + 2*KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,
      0)*ZH(gI2,1)) + 5*KroneckerDelta(2,gO2)*(-4*AbsSqr(Lambdax) + Sqr(gN))*(ZH(
      gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*((-20*AbsSqr
      (Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gI1,0)*ZH(gI2,0) - 3*(3*
      Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gI1,1)*ZH(gI2,1) + 5*(-4*AbsSqr(Lambdax)
      + Sqr(gN))*ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gI1,gI2)*(-5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gN) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN)) + KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(-3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN)));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.025*(5*KroneckerDelta(2,gO2)*(-
      2.8284271247461903*Conj(TLambdax)*ZA(gI1,1)*ZA(gI2,0) - 2.8284271247461903*
      TLambdax*ZA(gI1,1)*ZA(gI2,0) - 8*vs*AbsSqr(Lambdax)*ZA(gI1,1)*ZA(gI2,1) + 2*
      vs*Sqr(gN)*ZA(gI1,1)*ZA(gI2,1) + ZA(gI1,0)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(
      gN))*ZA(gI2,0) - 2.8284271247461903*(Conj(TLambdax) + TLambdax)*ZA(gI2,1)) -
      5*vs*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2)) + 2*KroneckerDelta(1,gO2)*(-
      7.0710678118654755*Conj(TLambdax)*ZA(gI1,2)*ZA(gI2,0) - 7.0710678118654755*
      TLambdax*ZA(gI1,2)*ZA(gI2,0) - 3*vu*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1) - 5*vu*Sqr(
      g2)*ZA(gI1,1)*ZA(gI2,1) - 2*vu*Sqr(gN)*ZA(gI1,1)*ZA(gI2,1) - 20*vu*AbsSqr(
      Lambdax)*ZA(gI1,2)*ZA(gI2,2) + 5*vu*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2) + ZA(gI1,0)*
      (vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gI2,0) -
      7.0710678118654755*(Conj(TLambdax) + TLambdax)*ZA(gI2,2))) - KroneckerDelta(
      0,gO2)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) + 5*ZA(
      gI1,2)*(2.8284271247461903*Conj(TLambdax)*ZA(gI2,1) + 2.8284271247461903*
      TLambdax*ZA(gI2,1) + vd*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZA(gI2,2)) + 2*ZA(
      gI1,1)*(vd*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZA(gI2,1
      ) + 7.0710678118654755*(Conj(TLambdax) + TLambdax)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhUhh(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(Conj(TLambdax) - TLambdax)*(KroneckerDelta(2,gO2)*(ZA(
      gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + KroneckerDelta(1,gO2)*(ZA(gI2,2)*
      ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2)) + KroneckerDelta(0,gO2)*(ZA(gI2,2)*ZH(gI1,1
      ) + ZA(gI2,1)*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.025*(5*KroneckerDelta(2,gO2)*(
      2.8284271247461903*Conj(TLambdax)*ZH(gI1,1)*ZH(gI2,0) + 2.8284271247461903*
      TLambdax*ZH(gI1,1)*ZH(gI2,0) - 8*vd*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,0) + 3*
      vd*Sqr(gN)*ZH(gI1,2)*ZH(gI2,0) - 8*vs*AbsSqr(Lambdax)*ZH(gI1,1)*ZH(gI2,1) +
      2*vs*Sqr(gN)*ZH(gI1,1)*ZH(gI2,1) - 8*vu*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,1)
      + 2*vu*Sqr(gN)*ZH(gI1,2)*ZH(gI2,1) - 8*vu*AbsSqr(Lambdax)*ZH(gI1,1)*ZH(gI2,2
      ) + 2*vu*Sqr(gN)*ZH(gI1,1)*ZH(gI2,2) - 15*vs*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2) +
      ZH(gI1,0)*(vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gI2,0) +
      2.8284271247461903*Conj(TLambdax)*ZH(gI2,1) + 2.8284271247461903*TLambdax*ZH
      (gI2,1) - 8*vd*AbsSqr(Lambdax)*ZH(gI2,2) + 3*vd*Sqr(gN)*ZH(gI2,2))) +
      KroneckerDelta(0,gO2)*(-(ZH(gI1,0)*(3*vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN)
      )*ZH(gI2,0) + 2*vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*
      ZH(gI2,1) + 5*vs*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gI2,2))) + 5*ZH(gI1,2)*(
      vs*(-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZH(gI2,0) + 2.8284271247461903*Conj(
      TLambdax)*ZH(gI2,1) + 2.8284271247461903*TLambdax*ZH(gI2,1) - 8*vd*AbsSqr(
      Lambdax)*ZH(gI2,2) + 3*vd*Sqr(gN)*ZH(gI2,2)) + 2*ZH(gI1,1)*(vu*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gI2,0) + vd*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gI2,1) + 7.0710678118654755
      *(Conj(TLambdax) + TLambdax)*ZH(gI2,2))) + 2*KroneckerDelta(1,gO2)*(ZH(gI1,1
      )*(vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gI2,0) -
      3*vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gI2,1) + 5*vs*(-4*AbsSqr(Lambdax
      ) + Sqr(gN))*ZH(gI2,2)) + ZH(gI1,0)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5
      *Sqr(g2) - 3*Sqr(gN))*ZH(gI2,0) + vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*
      Sqr(g2) - 3*Sqr(gN))*ZH(gI2,1) + 7.0710678118654755*(Conj(TLambdax) +
      TLambdax)*ZH(gI2,2)) + 5*ZH(gI1,2)*(1.4142135623730951*Conj(TLambdax)*ZH(gI2
      ,0) + 1.4142135623730951*TLambdax*ZH(gI2,0) + (-4*AbsSqr(Lambdax) + Sqr(gN))
      *(vs*ZH(gI2,1) + vu*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSvconjSv(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gI1,gI2)*(-5*vs*
      KroneckerDelta(2,gO2)*Sqr(gN) + vu*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(
      g2) + 2*Sqr(gN)) + vd*KroneckerDelta(0,gO2)*(-3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(
      gN)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO2)*
      SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gI1,j1)*ZDR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO1)*
      SUM(j1,0,2,Conj(ZDL(gI2,j1))*Conj(ZDR(gI1,j1))*Yd(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXFDXUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(2,gO2)*
      SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDXL(gI1,j1)*ZDXR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXFDXUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(2,gO1)*
      SUM(j1,0,2,Conj(ZDXL(gI2,j1))*Conj(ZDXR(gI1,j1))*Kappa(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO2)*
      SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gI1,j1)*ZER(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO1)*
      SUM(j1,0,2,Conj(ZEL(gI2,j1))*Conj(ZER(gI1,j1))*Ye(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO2)*
      SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gI1,j1)*ZUR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO1)*
      SUM(j1,0,2,Conj(ZUL(gI2,j1))*Conj(ZUR(gI1,j1))*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSHI0conjSHI0(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(
      UHI0(gI1,j1))*UHI0(gI2,j1)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1
      ,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))) - 10*KroneckerDelta(0,gO2)*(
      Lambdax*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gI2,j1)
      ) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,2 + j1)*Lambda12(j1
      ,j1)))) - KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(
      g2) + 9*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) - 2*(3*Sqr(g1)
      + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))
      ) + 20*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj
      (Lambda12(j1,j1))*UHI0(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gI1,j1)
      )*UHI0(gI2,2 + j1)*Lambda12(j1,j1)))) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(3*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))
      + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)) - 8*(SUM(j2,
      0,1,AbsSqr(Lambda12(j2,j2))*Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) + SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSHIpconjSHIp(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(
      UHIp(gI1,j1))*UHIp(gI2,j1)) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,
      1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))) + 10*KroneckerDelta(0,gO2)*(
      Lambdax*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gI2,j1)
      ) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,2 + j1)*Lambda12(j1
      ,j1)))) + KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*((6*Sqr(g1) - 10*
      Sqr(g2) + 9*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) + 2*(-3*Sqr
      (g1) + 5*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 +
      j1)))) + 20*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))
      *Conj(Lambda12(j1,j1))*UHIp(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(
      gI1,j1))*UHIp(gI2,2 + j1)*Lambda12(j1,j1)))) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(3*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))
      + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)) - 8*(SUM(j2,
      0,1,AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) + SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSHI0conjSHI0(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO2)*(vu*(3*Sqr(
      g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(gI1,j1)) -
      10*vd*Lambdax*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(
      gI1,j1)) - 3*vu*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*UHI0(gI1,2 + j1))
      - 5*vu*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*UHI0(gI1,2 + j1)) - 2*vu*
      Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*UHI0(gI1,2 + j1)) - 10*vd*Conj(
      Lambdax)*SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(gI1,2 + j1)*Lambda12(j1,j1))) -
      KroneckerDelta(0,gO2)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*SUM(j1,0,1,
      Conj(UHI0(gI2,j1))*UHI0(gI1,j1)) + 20*vu*Lambdax*SUM(j1,0,1,Conj(UHI0(gI2,2
      + j1))*Conj(Lambda12(j1,j1))*UHI0(gI1,j1)) - 6*vd*Sqr(g1)*SUM(j1,0,1,Conj(
      UHI0(gI2,2 + j1))*UHI0(gI1,2 + j1)) - 10*vd*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI2
      ,2 + j1))*UHI0(gI1,2 + j1)) + 6*vd*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))
      *UHI0(gI1,2 + j1)) + 20*vu*Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(
      gI1,2 + j1)*Lambda12(j1,j1))) + 5*KroneckerDelta(2,gO2)*(3*vs*Sqr(gN)*SUM(j1
      ,0,1,Conj(UHI0(gI2,j1))*UHI0(gI1,j1)) + 5.656854249492381*SUM(j1,0,1,Conj(
      UHI0(gI2,2 + j1))*Conj(TLambda12(j1,j1))*UHI0(gI1,j1)) + 2*vs*Sqr(gN)*SUM(j1
      ,0,1,Conj(UHI0(gI2,2 + j1))*UHI0(gI1,2 + j1)) + 5.656854249492381*SUM(j1,0,1
      ,Conj(UHI0(gI2,j1))*UHI0(gI1,2 + j1)*TLambda12(j1,j1)) - 8*vs*SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHI0(gI2,j2))*UHI0(gI1,j2)) - 8*vs*SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHI0(gI2,2 + j2))*UHI0(gI1,2 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSHIpconjSHIp(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO2)*(vu*(3*Sqr(
      g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(gI1,j1)) +
      10*vd*Lambdax*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(
      gI1,j1)) - 3*vu*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*UHIp(gI1,2 + j1))
      + 5*vu*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*UHIp(gI1,2 + j1)) - 2*vu*
      Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*UHIp(gI1,2 + j1)) + 10*vd*Conj(
      Lambdax)*SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(gI1,2 + j1)*Lambda12(j1,j1))) +
      KroneckerDelta(0,gO2)*(vd*(-6*Sqr(g1) + 10*Sqr(g2) - 9*Sqr(gN))*SUM(j1,0,1,
      Conj(UHIp(gI2,j1))*UHIp(gI1,j1)) + 20*vu*Lambdax*SUM(j1,0,1,Conj(UHIp(gI2,2
      + j1))*Conj(Lambda12(j1,j1))*UHIp(gI1,j1)) + 6*vd*Sqr(g1)*SUM(j1,0,1,Conj(
      UHIp(gI2,2 + j1))*UHIp(gI1,2 + j1)) - 10*vd*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI2
      ,2 + j1))*UHIp(gI1,2 + j1)) - 6*vd*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))
      *UHIp(gI1,2 + j1)) + 20*vu*Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(
      gI1,2 + j1)*Lambda12(j1,j1))) + 5*KroneckerDelta(2,gO2)*(3*vs*Sqr(gN)*SUM(j1
      ,0,1,Conj(UHIp(gI2,j1))*UHIp(gI1,j1)) - 5.656854249492381*SUM(j1,0,1,Conj(
      UHIp(gI2,2 + j1))*Conj(TLambda12(j1,j1))*UHIp(gI1,j1)) + 2*vs*Sqr(gN)*SUM(j1
      ,0,1,Conj(UHIp(gI2,2 + j1))*UHIp(gI1,2 + j1)) - 5.656854249492381*SUM(j1,0,1
      ,Conj(UHIp(gI2,j1))*UHIp(gI1,2 + j1)*TLambda12(j1,j1)) - 8*vs*SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gI2,j2))*UHIp(gI1,j2)) - 8*vs*SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gI2,2 + j2))*UHIp(gI1,2 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(2,gO2)*(
      SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gI1,2 + j1)*ZNI(gI2,j1)) + SUM(j1,0,1,
      Conj(Lambda12(j1,j1))*ZNI(gI1,j1)*ZNI(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.7071067811865475*KroneckerDelta(2,gO1)*(
      SUM(j1,0,1,Conj(ZNI(gI1,2 + j1))*Conj(ZNI(gI2,j1))*Lambda12(j1,j1)) + SUM(j1
      ,0,1,Conj(ZNI(gI1,j1))*Conj(ZNI(gI2,2 + j1))*Lambda12(j1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*Sqr(gN)*(SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) - 4*KroneckerDelta(1,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI1,3 + j1))*ZD(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZD(gI1,j1))*Yd(j1,j1)*ZD(gI2,3 + j1)))) +
      KroneckerDelta(1,gO1)*(-2*KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2) - Sqr(
      gN))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) - Sqr(gN))*SUM(j1,
      0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) + 20*KroneckerDelta(2,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI1,3 + j1))*ZD(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZD(gI1,j1))*Yd(j1,j1)*ZD(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr
      (gN))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + (4*Sqr(g1) + 6*Sqr(gN))*SUM(
      j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) - 40*(SUM(j2,0,2,AbsSqr(Yd(j2,j2
      ))*Conj(ZD(gI1,j2))*ZD(gI2,j2)) + SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(gI1,3
       + j2))*ZD(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSDXconjSDX(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*(4*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      ZDX(gI2,j1)) - 2*(2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*
      ZDX(gI2,3 + j1))) + 20*KroneckerDelta(0,gO2)*(Lambdax*SUM(j1,0,2,Conj(ZDX(
      gI1,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,2,Conj
      (ZDX(gI1,j1))*ZDX(gI2,3 + j1)*Kappa(j1,j1)))) + KroneckerDelta(0,gO1)*(-(
      KroneckerDelta(0,gO2)*((4*Sqr(g1) + 6*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      ZDX(gI2,j1)) + (-4*Sqr(g1) + 9*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX
      (gI2,3 + j1)))) + 20*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,2,Conj(ZDX(gI1,
      3 + j1))*Conj(Kappa(j1,j1))*ZDX(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX
      (gI1,j1))*ZDX(gI2,3 + j1)*Kappa(j1,j1)))) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(2*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) +
      3*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)) - 8*(SUM(j2,0,2,
      AbsSqr(Kappa(j2,j2))*Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) + SUM(j2,0,2,AbsSqr(
      Kappa(j2,j2))*Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*Sqr(gN)*(2*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) +
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) - 4*KroneckerDelta(1,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI1,3 + j1))*ZE(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZE(gI1,j1))*Ye(j1,j1)*ZE(gI2,3 + j1)))) + 2*
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr
      (gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) + 10*KroneckerDelta(2,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI1,3 + j1))*ZE(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZE(gI1,j1))*Ye(j1,j1)*ZE(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-6*Sqr(g1) + 10*Sqr(g2) + 6*
      Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 3*(4*Sqr(g1) + Sqr(gN))*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) - 40*(SUM(j2,0,2,AbsSqr(Ye(
      j2,j2))*Conj(ZE(gI1,j2))*ZE(gI2,j2)) + SUM(j2,0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*Sqr(gN)*(SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + SUM
      (j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) - 4*KroneckerDelta(0,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI1,3 + j1))*ZU(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*Yu(j1,j1)*ZU(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((2*Sqr(g1) - 10*Sqr(g2) + 3*
      Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (-8*Sqr(g1) + 3*Sqr(gN))*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) + 20*KroneckerDelta(2,gO2)*
      (Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI1,3 + j1))*ZU(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*Yu(j1,j1)*ZU(gI2,3 + j1)))) + 2*
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-Sqr(g1) + 5*Sqr(g2) + Sqr(gN)
      )*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2
      ,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) - 20*(SUM(j2,0,2,AbsSqr(Yu(j2,j2))*
      Conj(ZU(gI1,j2))*ZU(gI2,j2)) + SUM(j2,0,2,AbsSqr(Yu(j2,j2))*Conj(ZU(gI1,3 +
      j2))*ZU(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSdconjSd(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO2)*(vs*Sqr(gN)
      *SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) - 4*vu*Lambdax*SUM(j1,0,2,Conj(Yd(
      j1,j1))*Conj(ZD(gI2,3 + j1))*ZD(gI1,j1)) + 2*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZD(
      gI2,3 + j1))*ZD(gI1,3 + j1)) - 4*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gI2,j1)
      )*Yd(j1,j1)*ZD(gI1,3 + j1))) + KroneckerDelta(1,gO2)*(-2*vu*(Sqr(g1) + 5*Sqr
      (g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) + 4*(5*vs*Lambdax*
      SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI2,3 + j1))*ZD(gI1,j1)) + vu*(-Sqr(g1) +
      Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)) + 5*vs*Conj(Lambdax
      )*SUM(j1,0,2,Conj(ZD(gI2,j1))*Yd(j1,j1)*ZD(gI1,3 + j1)))) + KroneckerDelta(0
      ,gO2)*(vd*(2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI2,j1))*
      ZD(gI1,j1)) - 28.284271247461902*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(TYd(j1
      ,j1))*ZD(gI1,j1)) + 4*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 +
      j1)) + 6*vd*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)) -
      28.284271247461902*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,3 + j1)*TYd(j1,j1)) -
      40*vd*SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(gI2,j2))*ZD(gI1,j2)) - 40*vd*SUM(
      j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSDXconjSDX(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(KroneckerDelta(1,gO2)*(4*vu*(Sqr(g1)
      - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)) + 20*vd*Lambdax*SUM(j1,
      0,2,Conj(ZDX(gI2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gI1,j1)) - 4*vu*Sqr(g1)*SUM
      (j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)) - 6*vu*Sqr(gN)*SUM(j1,0,2,
      Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)) + 20*vd*Conj(Lambdax)*SUM(j1,0,2,Conj
      (ZDX(gI2,j1))*ZDX(gI1,3 + j1)*Kappa(j1,j1))) + KroneckerDelta(0,gO2)*(-2*vd*
      (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)) + 20*vu*
      Lambdax*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gI1,j1)) + 4
      *vd*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)) - 9*vd*Sqr(gN)
      *SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)) + 20*vu*Conj(Lambdax)*SUM
      (j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,3 + j1)*Kappa(j1,j1))) + 5*KroneckerDelta(
      2,gO2)*(2*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)) -
      5.656854249492381*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(TKappa(j1,j1))*ZDX(
      gI1,j1)) + 3*vs*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*ZDX(gI1,3 + j1)) -
      5.656854249492381*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,3 + j1)*TKappa(j1,j1)
      ) - 8*vs*SUM(j2,0,2,AbsSqr(Kappa(j2,j2))*Conj(ZDX(gI2,j2))*ZDX(gI1,j2)) - 8*
      vs*SUM(j2,0,2,AbsSqr(Kappa(j2,j2))*Conj(ZDX(gI2,3 + j2))*ZDX(gI1,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSeconjSe(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO2)*(vu*(3*Sqr(
      g1) - 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) + 10*vs
      *Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI2,3 + j1))*ZE(gI1,j1)) - 6*vu*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) + vu*Sqr(gN)*SUM(j1,
      0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) + 10*vs*Conj(Lambdax)*SUM(j1,0,2,
      Conj(ZE(gI2,j1))*Ye(j1,j1)*ZE(gI1,3 + j1))) - 5*KroneckerDelta(2,gO2)*(2*vs*
      Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) - 4*vu*Lambdax*SUM(j1,0,2,
      Conj(Ye(j1,j1))*Conj(ZE(gI2,3 + j1))*ZE(gI1,j1)) + vs*Sqr(gN)*SUM(j1,0,2,
      Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) - 4*vu*Conj(Lambdax)*SUM(j1,0,2,Conj(ZE
      (gI2,j1))*Ye(j1,j1)*ZE(gI1,3 + j1))) + KroneckerDelta(0,gO2)*(2*vd*(-3*Sqr(
      g1) + 5*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) -
      28.284271247461902*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j1))*ZE(gI1,
      j1)) + 12*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) + 3*vd*
      Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) - 28.284271247461902
      *SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,3 + j1)*TYe(j1,j1)) - 40*vd*SUM(j2,0,2,
      AbsSqr(Ye(j2,j2))*Conj(ZE(gI2,j2))*ZE(gI1,j2)) - 40*vd*SUM(j2,0,2,AbsSqr(Ye(
      j2,j2))*Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSuconjSu(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO2)*(vs*Sqr(gN)
      *SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*vd*Lambdax*SUM(j1,0,2,Conj(Yu(
      j1,j1))*Conj(ZU(gI2,3 + j1))*ZU(gI1,j1)) + vs*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI2
      ,3 + j1))*ZU(gI1,3 + j1)) - 4*vd*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI2,j1))*
      Yu(j1,j1)*ZU(gI1,3 + j1))) + KroneckerDelta(0,gO2)*(vd*(2*Sqr(g1) - 10*Sqr(
      g2) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) + 20*vs*Lambdax*SUM
      (j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI2,3 + j1))*ZU(gI1,j1)) - 8*vd*Sqr(g1)*SUM(
      j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + 3*vd*Sqr(gN)*SUM(j1,0,2,Conj(
      ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + 20*vs*Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI2
      ,j1))*Yu(j1,j1)*ZU(gI1,3 + j1))) - 2*KroneckerDelta(1,gO2)*(vu*(Sqr(g1) - 5*
      Sqr(g2) - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) +
      14.142135623730951*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j1))*ZU(gI1,
      j1)) - 4*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) - vu*Sqr
      (gN)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + 14.142135623730951*
      SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,3 + j1)*TYu(j1,j1)) + 20*vu*SUM(j2,0,2,
      AbsSqr(Yu(j2,j2))*Conj(ZU(gI2,j2))*ZU(gI1,j2)) + 20*vu*SUM(j2,0,2,AbsSqr(Yu(
      j2,j2))*Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(2*KroneckerDelta(1,gO2)*((-
      3.872983346207417*g1*ZN(gI1,0) + 5*g2*ZN(gI1,1) + 3.1622776601683795*gN*ZN(
      gI1,5))*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2,2) +
      ZN(gI1,2)*ZN(gI2,4)) + ZN(gI1,3)*(-3.872983346207417*g1*ZN(gI2,0) + 5*g2*ZN(
      gI2,1) + 3.1622776601683795*gN*ZN(gI2,5))) + KroneckerDelta(0,gO2)*(
      7.745966692414834*g1*ZN(gI1,0)*ZN(gI2,2) - 10*g2*ZN(gI1,1)*ZN(gI2,2) +
      9.486832980505138*gN*ZN(gI1,5)*ZN(gI2,2) + 14.142135623730951*Conj(Lambdax)*
      ZN(gI1,4)*ZN(gI2,3) + 14.142135623730951*Conj(Lambdax)*ZN(gI1,3)*ZN(gI2,4) +
      ZN(gI1,2)*(7.745966692414834*g1*ZN(gI2,0) - 10*g2*ZN(gI2,1) +
      9.486832980505138*gN*ZN(gI2,5))) + 7.0710678118654755*KroneckerDelta(2,gO2)*
      (2*Conj(Lambdax)*(ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,3)) -
      2.23606797749979*gN*(ZN(gI1,5)*ZN(gI2,4) + ZN(gI1,4)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.05*(-10*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2)
      )*KroneckerDelta(0,gO1) + 9.486832980505138*gN*Conj(ZN(gI1,5))*Conj(ZN(gI2,2
      ))*KroneckerDelta(0,gO1) - 7.745966692414834*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,
      0))*KroneckerDelta(1,gO1) + 10*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*
      KroneckerDelta(1,gO1) + 10*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta
      (1,gO1) + 6.324555320336759*gN*Conj(ZN(gI1,5))*Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1) + 6.324555320336759*gN*Conj(ZN(gI1,3))*Conj(ZN(gI2,5))
      *KroneckerDelta(1,gO1) + 7.745966692414834*g1*Conj(ZN(gI1,0))*(Conj(ZN(gI2,2
      ))*KroneckerDelta(0,gO1) - Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)) -
      15.811388300841898*gN*Conj(ZN(gI1,5))*Conj(ZN(gI2,4))*KroneckerDelta(2,gO1)
      - 15.811388300841898*gN*Conj(ZN(gI1,4))*Conj(ZN(gI2,5))*KroneckerDelta(2,gO1
      ) + 14.142135623730951*Conj(ZN(gI1,4))*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1)
      *Lambdax + 14.142135623730951*Conj(ZN(gI1,3))*Conj(ZN(gI2,4))*KroneckerDelta
      (0,gO1)*Lambdax + 14.142135623730951*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*
      KroneckerDelta(1,gO1)*Lambdax + 14.142135623730951*Conj(ZN(gI1,3))*Conj(ZN(
      gI2,2))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(7.745966692414834*
      g1*Conj(ZN(gI2,0))*KroneckerDelta(0,gO1) - 10*g2*Conj(ZN(gI2,1))*
      KroneckerDelta(0,gO1) + 9.486832980505138*gN*Conj(ZN(gI2,5))*KroneckerDelta(
      0,gO1) + 14.142135623730951*(Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) + Conj(ZN
      (gI2,3))*KroneckerDelta(2,gO1))*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjVWm(int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.5*(-(g2*KroneckerDelta(0,gO2)*ZP(gI2,0))
      + g2*KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhUhhVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      KroneckerDelta(0,gO2)*(10*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*Sin
      (ThetaWp()))*ZA(gI2,0) - 2*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      3.1622776601683795*gN*Sin(ThetaWp()))*ZA(gI2,1) + 15.811388300841898*gN*
      KroneckerDelta(2,gO2)*Sin(ThetaWp())*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpAhUhhVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.05)*(
      KroneckerDelta(0,gO2)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZA(gI2,0) +
      2*KroneckerDelta(1,gO2)*(3.1622776601683795*gN*Cos(ThetaWp()) - (5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZA(gI2,1) -
      15.811388300841898*gN*Cos(ThetaWp())*KroneckerDelta(2,gO2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgWmUAh(int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.25)*(vd*
      KroneckerDelta(0,gO1) - vu*KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgWmCUAh(int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(vd*
      KroneckerDelta(0,gO1) - vu*KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gN)*Sqr(Sin(ThetaWp())) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(-14.696938456699067*g1*gN*Cos(ThetaWp())*Sin(ThetaW()
      )*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + Cos(
      ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp())*Sin(ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 6*Sqr(g1)*Sqr(
      Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(ThetaWp()))) + 2*
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3.1622776601683795*g2*gN*Cos(
      ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(
      2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(25*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gN)*Sqr(Cos(ThetaWp())) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(9*Sqr(gN)*Sqr(Cos(ThetaWp())) + 10*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(2.449489742783178*gN*
      Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))) + Cos(ThetaW())*(
      9.486832980505138*g2*gN*Sin(2*ThetaWp()) + 15.491933384829668*g1*g2*Sin(
      ThetaW())*Sqr(Sin(ThetaWp())))) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(-3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 2*Sqr(gN)*
      Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*
      (-2.449489742783178*gN*Sin(ThetaW())*Sin(2*ThetaWp()) + 3.872983346207417*g2
      *Sin(2*ThetaW())*Sqr(Sin(ThetaWp())) + 3*g1*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0)
      + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP(gI1,1)*ZP(gI2,1)) + 2*KroneckerDelta(1
      ,gO1)*(5*KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(
      gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2
      ) - 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(
      gI1,1)*ZP(gI2,1))) - KroneckerDelta(0,gO1)*(-10*KroneckerDelta(1,gO2)*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gI1,0)*ZP(gI2
      ,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSHp0conjSHp0(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN)) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gN) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*
      (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN)))*(Conj(UHp0(gI1,0))*UHp0(gI2,0) - Conj(
      UHp0(gI1,1))*UHp0(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSHppconjSHpp(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN)) + 5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gN))*(Conj(UHpp(gI1,0))*UHpp
      (gI2,0) - Conj(UHpp(gI1,1))*UHpp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSSI0conjSSI0(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.125*(3*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) - 5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2))*KroneckerDelta(gI1,gI2)*Sqr(gN)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUAhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(vu*
      KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) -
      ZP(gI1,0)*ZP(gI2,1)) + vd*KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2
      ))*(ZP(gI1,1)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,1)) + 2.8284271247461903*
      KroneckerDelta(2,gO2)*(-(TLambdax*ZP(gI1,1)*ZP(gI2,0)) + Conj(TLambdax)*ZP(
      gI1,0)*ZP(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*(g2*KroneckerDelta(0,gO2)*UM(gI1,1)*UP(gI2,0) + (g2*KroneckerDelta(1,gO2)*
      UM(gI1,0) - Conj(Lambdax)*KroneckerDelta(2,gO2)*UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *(g2*Conj(UM(gI2,0))*Conj(UP(gI1,1))*KroneckerDelta(1,gO1) + Conj(UM(gI2,1))
      *(g2*Conj(UP(gI1,0))*KroneckerDelta(0,gO1) - Conj(UP(gI1,1))*KroneckerDelta(
      2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIChaIUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(2,gO2)*SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZMI(gI1,j1)*ZPI(gI2,
      j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIChaIUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(2,gO1)*SUM(j1,0,1,Conj(ZMI(gI2,j1))*Conj(ZPI(gI1,j1))*
      Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*(-2*
      KroneckerDelta(1,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(
      gN))*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) - 5*KroneckerDelta(2,gO2)*(
      -8*AbsSqr(Lambdax) + 3*Sqr(gN))*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2))
      + KroneckerDelta(0,gO2)*(3*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZA(gI1,0)*ZA
      (gI2,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZA(gI1,
      1)*ZA(gI2,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZA(gI1,2)*ZA(gI2,2)))) + 5*
      KroneckerDelta(2,gO1)*(KroneckerDelta(0,gO2)*(-8*AbsSqr(Lambdax) + 3*Sqr(gN)
      )*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) + 2*KroneckerDelta(1,gO2)*(-4*
      AbsSqr(Lambdax) + Sqr(gN))*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)) +
      KroneckerDelta(2,gO2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0)
      + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gI1,1)*ZA(gI2,1) - 15*Sqr(gN)*ZA(gI1,2
      )*ZA(gI2,2))) + 2*KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,
      0)*ZA(gI2,1)) + 5*KroneckerDelta(2,gO2)*(-4*AbsSqr(Lambdax) + Sqr(gN))*(ZA(
      gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)) + KroneckerDelta(1,gO2)*((-20*AbsSqr
      (Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) - 3*(3*
      Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZA(gI1,1)*ZA(gI2,1) + 5*(-4*AbsSqr(Lambdax)
      + Sqr(gN))*ZA(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH(gI1,0)*ZH(gI2
      ,0) + 2*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(gI1,1)*
      ZH(gI2,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gI1,2)*ZH(gI2,2))) + 5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*((-8*AbsSqr(Lambdax) + 3*Sqr(gN)
      )*ZH(gI1,0)*ZH(gI2,0) + 2*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gI1,1)*ZH(gI2,1)
      - 5*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2)) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZH(gI1,0)*ZH
      (gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gI1,1)*ZH(gI2,1) + 5*(-4*
      AbsSqr(Lambdax) + Sqr(gN))*ZH(gI1,2)*ZH(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gI1,gI2)*(-5*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gN) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN)) + KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(-3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN)));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373
      )*(Conj(TLambdax) - TLambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,1)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0
      )*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2
      ,2)));

   return result;
}

std::complex<double> CLASSNAME::CpAhUAhhh(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO2)*(
      2.8284271247461903*Conj(TLambdax)*(ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)
      ) + 2.8284271247461903*TLambdax*(ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1))
      + ZA(gI2,2)*((8*vd*AbsSqr(Lambdax) - 3*vd*Sqr(gN))*ZH(gI1,0) - 2*vu*(-4*
      AbsSqr(Lambdax) + Sqr(gN))*ZH(gI1,1) + 5*vs*Sqr(gN)*ZH(gI1,2))) + 2*
      KroneckerDelta(1,gO2)*(ZA(gI2,1)*(vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*
      Sqr(g2) - 3*Sqr(gN))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(
      gI1,1) + 5*vs*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(gI1,2)) - 7.0710678118654755
      *(Conj(TLambdax) + TLambdax)*(ZA(gI2,2)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2))) -
      KroneckerDelta(0,gO2)*(ZA(gI2,0)*(vd*(6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZH
      (gI1,0) + 2*vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZH(
      gI1,1) + 5*vs*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gI1,2)) +
      14.142135623730951*(Conj(TLambdax) + TLambdax)*(ZA(gI2,2)*ZH(gI1,1) + ZA(gI2
      ,1)*ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(Conj(TLambdax) - TLambdax)*(KroneckerDelta(2,gO2)*(ZH(
      gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(1,gO2)*(ZH(gI1,2)*
      ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1
      ) + ZH(gI1,1)*ZH(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gI1,j1)*ZDR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(0,gO1)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*Conj(ZDR(gI1,j1))*Yd(j1
      ,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXFDXUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(Kappa(j1,j1))*ZDXL(gI1,j1)*ZDXR(gI2,
      j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXFDXUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(ZDXL(gI2,j1))*Conj(ZDXR(gI1,j1))*
      Kappa(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gI1,j1)*ZER(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(0,gO1)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*Conj(ZER(gI1,j1))*Ye(j1
      ,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gI1,j1)*ZUR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*Conj(ZUR(gI1,j1))*Yu(j1
      ,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSHI0conjSHI0(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(
      UHI0(gI1,j1))*UHI0(gI2,j1)) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1
      ,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))) + 10*KroneckerDelta(0,gO2)*(
      Lambdax*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gI2,j1)
      ) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,2 + j1)*Lambda12(j1
      ,j1)))) + KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*
      Sqr(g2) + 9*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) - 2*(3*Sqr(
      g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 +
      j1)))) + 20*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))
      *Conj(Lambda12(j1,j1))*UHI0(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(
      gI1,j1))*UHI0(gI2,2 + j1)*Lambda12(j1,j1)))) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(3*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))
      + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)) - 8*(SUM(j2,
      0,1,AbsSqr(Lambda12(j2,j2))*Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) + SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSHIpconjSHIp(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(
      UHIp(gI1,j1))*UHIp(gI2,j1)) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,
      1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))) - 10*KroneckerDelta(0,gO2)*(
      Lambdax*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gI2,j1)
      ) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,2 + j1)*Lambda12(j1
      ,j1)))) - KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((6*Sqr(g1) - 10*Sqr(
      g2) + 9*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) + 2*(-3*Sqr(g1)
      + 5*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))
      ) + 20*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*Conj
      (Lambda12(j1,j1))*UHIp(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gI1,j1)
      )*UHIp(gI2,2 + j1)*Lambda12(j1,j1)))) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(3*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))
      + 2*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)) - 8*(SUM(j2,
      0,1,AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) + SUM(j2,0,1,
      AbsSqr(Lambda12(j2,j2))*Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSHI0conjSHI0(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(vu*
      KroneckerDelta(0,gO2)*(Lambdax*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*Conj(
      Lambda12(j1,j1))*UHI0(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gI2,j1))
      *UHI0(gI1,2 + j1)*Lambda12(j1,j1))) + vd*KroneckerDelta(1,gO2)*(Lambdax*SUM(
      j1,0,1,Conj(UHI0(gI2,2 + j1))*Conj(Lambda12(j1,j1))*UHI0(gI1,j1)) - Conj(
      Lambdax)*SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(gI1,2 + j1)*Lambda12(j1,j1))) +
      1.4142135623730951*KroneckerDelta(2,gO2)*(SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*
      Conj(TLambda12(j1,j1))*UHI0(gI1,j1)) - SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(
      gI1,2 + j1)*TLambda12(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSHIpconjSHIp(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(vu*
      KroneckerDelta(0,gO2)*(Lambdax*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*Conj(
      Lambda12(j1,j1))*UHIp(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gI2,j1))
      *UHIp(gI1,2 + j1)*Lambda12(j1,j1))) + vd*KroneckerDelta(1,gO2)*(Lambdax*SUM(
      j1,0,1,Conj(UHIp(gI2,2 + j1))*Conj(Lambda12(j1,j1))*UHIp(gI1,j1)) - Conj(
      Lambdax)*SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(gI1,2 + j1)*Lambda12(j1,j1))) +
      1.4142135623730951*KroneckerDelta(2,gO2)*(SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*
      Conj(TLambda12(j1,j1))*UHIp(gI1,j1)) - SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(
      gI1,2 + j1)*TLambda12(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(2,gO2)*(SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gI1,2 + j1)*
      ZNI(gI2,j1)) + SUM(j1,0,1,Conj(Lambda12(j1,j1))*ZNI(gI1,j1)*ZNI(gI2,2 + j1))
      );

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(2,gO1)*(SUM(j1,0,1,Conj(ZNI(gI1,2 + j1))*Conj(ZNI(gI2,j1))*
      Lambda12(j1,j1)) + SUM(j1,0,1,Conj(ZNI(gI1,j1))*Conj(ZNI(gI2,2 + j1))*
      Lambda12(j1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*Sqr(gN)*(SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) + 4*KroneckerDelta(1,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI1,3 + j1))*ZD(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZD(gI1,j1))*Yd(j1,j1)*ZD(gI2,3 + j1)))) - 2*
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2) - Sqr(gN)
      )*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2
      ,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) + 10*KroneckerDelta(2,gO2)*(Lambdax*
      SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI1,3 + j1))*ZD(gI2,j1)) + Conj(Lambdax)*
      SUM(j1,0,2,Conj(ZD(gI1,j1))*Yd(j1,j1)*ZD(gI2,3 + j1)))) + KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,2,
      Conj(ZD(gI1,j1))*ZD(gI2,j1)) + (4*Sqr(g1) + 6*Sqr(gN))*SUM(j1,0,2,Conj(ZD(
      gI1,3 + j1))*ZD(gI2,3 + j1)) - 40*(SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(gI1,
      j2))*ZD(gI2,j2)) + SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZD(gI1,3 + j2))*ZD(gI2,
      3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSDXconjSDX(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*(4*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      ZDX(gI2,j1)) - 2*(2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*
      ZDX(gI2,3 + j1))) - 20*KroneckerDelta(0,gO2)*(Lambdax*SUM(j1,0,2,Conj(ZDX(
      gI1,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,2,Conj
      (ZDX(gI1,j1))*ZDX(gI2,3 + j1)*Kappa(j1,j1)))) - KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*((4*Sqr(g1) + 6*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      ZDX(gI2,j1)) + (-4*Sqr(g1) + 9*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX
      (gI2,3 + j1))) + 20*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,2,Conj(ZDX(gI1,3
       + j1))*Conj(Kappa(j1,j1))*ZDX(gI2,j1)) + Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(
      gI1,j1))*ZDX(gI2,3 + j1)*Kappa(j1,j1)))) + 5*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(2*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) +
      3*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)) - 8*(SUM(j2,0,2,
      AbsSqr(Kappa(j2,j2))*Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) + SUM(j2,0,2,AbsSqr(
      Kappa(j2,j2))*Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*Sqr(gN)*(2*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) +
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) + 4*KroneckerDelta(1,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI1,3 + j1))*ZE(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZE(gI1,j1))*Ye(j1,j1)*ZE(gI2,3 + j1)))) + 2*
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr
      (gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) - 10*KroneckerDelta(2,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI1,3 + j1))*ZE(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZE(gI1,j1))*Ye(j1,j1)*ZE(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-6*Sqr(g1) + 10*Sqr(g2) + 6*
      Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 3*(4*Sqr(g1) + Sqr(gN))*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) - 40*(SUM(j2,0,2,AbsSqr(Ye(
      j2,j2))*Conj(ZE(gI1,j2))*ZE(gI2,j2)) + SUM(j2,0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-5*KroneckerDelta(2,gO1)*(
      KroneckerDelta(2,gO2)*Sqr(gN)*(SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + SUM
      (j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) + 4*KroneckerDelta(0,gO2)*(
      Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI1,3 + j1))*ZU(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*Yu(j1,j1)*ZU(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((2*Sqr(g1) - 10*Sqr(g2) + 3*
      Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (-8*Sqr(g1) + 3*Sqr(gN))*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) - 20*KroneckerDelta(2,gO2)*
      (Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI1,3 + j1))*ZU(gI2,j1)) + Conj(
      Lambdax)*SUM(j1,0,2,Conj(ZU(gI1,j1))*Yu(j1,j1)*ZU(gI2,3 + j1)))) + 2*
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-Sqr(g1) + 5*Sqr(g2) + Sqr(gN)
      )*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2
      ,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) - 20*(SUM(j2,0,2,AbsSqr(Yu(j2,j2))*
      Conj(ZU(gI1,j2))*ZU(gI2,j2)) + SUM(j2,0,2,AbsSqr(Yu(j2,j2))*Conj(ZU(gI1,3 +
      j2))*ZU(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSdconjSd(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(vs*
      KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI2,3 + j1
      ))*ZD(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gI2,j1))*Yd(j1,j1)*ZD(gI1,
      3 + j1))) + vu*KroneckerDelta(2,gO2)*(Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*
      Conj(ZD(gI2,3 + j1))*ZD(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZD(gI2,j1))
      *Yd(j1,j1)*ZD(gI1,3 + j1))) + 1.4142135623730951*KroneckerDelta(0,gO2)*(SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j1))*ZD(gI1,j1)) - SUM(j1,0,2,Conj(
      ZD(gI2,j1))*ZD(gI1,3 + j1)*TYd(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSDXconjSDX(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(vu*
      KroneckerDelta(0,gO2)*(Lambdax*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(Kappa(
      j1,j1))*ZDX(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,3
      + j1)*Kappa(j1,j1))) + vd*KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,2,Conj(ZDX
      (gI2,3 + j1))*Conj(Kappa(j1,j1))*ZDX(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,
      Conj(ZDX(gI2,j1))*ZDX(gI1,3 + j1)*Kappa(j1,j1))) + 1.4142135623730951*
      KroneckerDelta(2,gO2)*(SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(TKappa(j1,j1))*
      ZDX(gI1,j1)) - SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,3 + j1)*TKappa(j1,j1))))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUAhSeconjSe(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(vs*
      KroneckerDelta(1,gO2)*(Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI2,3 + j1
      ))*ZE(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gI2,j1))*Ye(j1,j1)*ZE(gI1,
      3 + j1))) + vu*KroneckerDelta(2,gO2)*(Lambdax*SUM(j1,0,2,Conj(Ye(j1,j1))*
      Conj(ZE(gI2,3 + j1))*ZE(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZE(gI2,j1))
      *Ye(j1,j1)*ZE(gI1,3 + j1))) + 1.4142135623730951*KroneckerDelta(0,gO2)*(SUM(
      j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j1))*ZE(gI1,j1)) - SUM(j1,0,2,Conj(
      ZE(gI2,j1))*ZE(gI1,3 + j1)*TYe(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSuconjSu(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(vs*
      KroneckerDelta(0,gO2)*(Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI2,3 + j1
      ))*ZU(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI2,j1))*Yu(j1,j1)*ZU(gI1,
      3 + j1))) + vd*KroneckerDelta(2,gO2)*(Lambdax*SUM(j1,0,2,Conj(Yu(j1,j1))*
      Conj(ZU(gI2,3 + j1))*ZU(gI1,j1)) - Conj(Lambdax)*SUM(j1,0,2,Conj(ZU(gI2,j1))
      *Yu(j1,j1)*ZU(gI1,3 + j1))) + 1.4142135623730951*KroneckerDelta(1,gO2)*(SUM(
      j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j1))*ZU(gI1,j1)) - SUM(j1,0,2,Conj(
      ZU(gI2,j1))*ZU(gI1,3 + j1)*TYu(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.05)*(
      KroneckerDelta(0,gO2)*(-7.745966692414834*g1*ZN(gI1,0)*ZN(gI2,2) + 10*g2*ZN(
      gI1,1)*ZN(gI2,2) - 9.486832980505138*gN*ZN(gI1,5)*ZN(gI2,2) +
      14.142135623730951*Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,3) + 14.142135623730951*
      Conj(Lambdax)*ZN(gI1,3)*ZN(gI2,4) + ZN(gI1,2)*(-7.745966692414834*g1*ZN(gI2,
      0) + 10*g2*ZN(gI2,1) - 9.486832980505138*gN*ZN(gI2,5))) + 2*KroneckerDelta(1
      ,gO2)*((3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1) - 3.1622776601683795
      *gN*ZN(gI1,5))*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*(ZN(gI1,4)*ZN(
      gI2,2) + ZN(gI1,2)*ZN(gI2,4)) + ZN(gI1,3)*(3.872983346207417*g1*ZN(gI2,0) -
      5*g2*ZN(gI2,1) - 3.1622776601683795*gN*ZN(gI2,5))) + 7.0710678118654755*
      KroneckerDelta(2,gO2)*(2*Conj(Lambdax)*(ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(
      gI2,3)) + 2.23606797749979*gN*(ZN(gI1,5)*ZN(gI2,4) + ZN(gI1,4)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.05)*(10*g2*Conj(ZN
      (gI1,1))*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) - 9.486832980505138*gN*Conj(
      ZN(gI1,5))*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) + 7.745966692414834*g1*Conj
      (ZN(gI1,3))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) - 10*g2*Conj(ZN(gI1,3))*
      Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) - 10*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3)
      )*KroneckerDelta(1,gO1) - 6.324555320336759*gN*Conj(ZN(gI1,5))*Conj(ZN(gI2,3
      ))*KroneckerDelta(1,gO1) - 6.324555320336759*gN*Conj(ZN(gI1,3))*Conj(ZN(gI2,
      5))*KroneckerDelta(1,gO1) + 7.745966692414834*g1*Conj(ZN(gI1,0))*(-(Conj(ZN(
      gI2,2))*KroneckerDelta(0,gO1)) + Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)) +
      15.811388300841898*gN*Conj(ZN(gI1,5))*Conj(ZN(gI2,4))*KroneckerDelta(2,gO1)
      + 15.811388300841898*gN*Conj(ZN(gI1,4))*Conj(ZN(gI2,5))*KroneckerDelta(2,gO1
      ) + 14.142135623730951*Conj(ZN(gI1,4))*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1)
      *Lambdax + 14.142135623730951*Conj(ZN(gI1,3))*Conj(ZN(gI2,4))*KroneckerDelta
      (0,gO1)*Lambdax + 14.142135623730951*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*
      KroneckerDelta(1,gO1)*Lambdax + 14.142135623730951*Conj(ZN(gI1,3))*Conj(ZN(
      gI2,2))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(-7.745966692414834*
      g1*Conj(ZN(gI2,0))*KroneckerDelta(0,gO1) + 10*g2*Conj(ZN(gI2,1))*
      KroneckerDelta(0,gO1) - 9.486832980505138*gN*Conj(ZN(gI2,5))*KroneckerDelta(
      0,gO1) + 14.142135623730951*(Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) + Conj(ZN
      (gI2,3))*KroneckerDelta(2,gO1))*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHpmconjVWm(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZP(gI2,0) + KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhVZ(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      KroneckerDelta(0,gO2)*(10*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*Sin
      (ThetaWp()))*ZH(gI2,0) - 2*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      3.1622776601683795*gN*Sin(ThetaWp()))*ZH(gI2,1) + 15.811388300841898*gN*
      KroneckerDelta(2,gO2)*Sin(ThetaWp())*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhVZp(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.05)*(
      KroneckerDelta(0,gO2)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZH(gI2,0) +
      2*KroneckerDelta(1,gO2)*(3.1622776601683795*gN*Cos(ThetaWp()) - (5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZH(gI2,1) -
      15.811388300841898*gN*Cos(ThetaWp())*KroneckerDelta(2,gO2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgZUHpm(int gO2) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO2)*(-5*g2
      *Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp())) + vd*KroneckerDelta(0,gO2)
      *(10*g2*Cos(ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 9.486832980505138*gN*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmconjUHpm(int gO1) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO1)*(5*g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp())) + vd*KroneckerDelta(0,gO1)
      *(-10*g2*Cos(ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 9.486832980505138*gN*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgZconjUHpm(int gO1) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO1)*(-5*g2
      *Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp())) + vd*KroneckerDelta(0,gO1)
      *(10*g2*Cos(ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 9.486832980505138*gN*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmCUHpm(int gO2) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO2)*(5*g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp())) + vd*KroneckerDelta(0,gO2)
      *(-10*g2*Cos(ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 9.486832980505138*gN*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgZpUHpm(int gO2) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO2)*(
      3.1622776601683795*gN*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) + vd*KroneckerDelta(0,
      gO2)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgWmconjUHpm(int gO1) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO1)*(
      3.1622776601683795*gN*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) + vd*KroneckerDelta(0,
      gO1)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgZpconjUHpm(int gO1) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO1)*(
      3.1622776601683795*gN*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) + vd*KroneckerDelta(0,
      gO1)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgWmCUHpm(int gO2) const
{
   
   const std::complex<double> result = 0.025*g2*(2*vu*KroneckerDelta(1,gO2)*(
      3.1622776601683795*gN*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) + vd*KroneckerDelta(0,
      gO2)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVPVWm(int gO2) const
{
   
   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd
      *KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVZ(int gO2) const
{
   
   const std::complex<double> result = 0.11180339887498948*g2*(vd*KroneckerDelta(0
      ,gO2)*(3.4641016151377544*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      4.242640687119286*gN*Sin(ThetaWp())) - 2*vu*KroneckerDelta(1,gO2)*(
      1.7320508075688772*g1*Cos(ThetaWp())*Sin(ThetaW()) + 1.4142135623730951*gN*
      Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVZp(int gO2) const
{
   
   const std::complex<double> result = -0.11180339887498948*g2*(2*vu*
      KroneckerDelta(1,gO2)*(1.4142135623730951*gN*Cos(ThetaWp()) -
      1.7320508075688772*g1*Sin(ThetaW())*Sin(ThetaWp())) + vd*KroneckerDelta(0,
      gO2)*(4.242640687119286*gN*Cos(ThetaWp()) + 3.4641016151377544*g1*Sin(ThetaW
      ())*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(2*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(2.449489742783178*g1*gN*Sin(ThetaW())*Sin(2*ThetaWp()
      ) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 2*Cos(ThetaW())*(
      3.1622776601683795*g2*gN*Cos(ThetaWp())*Sin(ThetaWp()) + 3.872983346207417*
      g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr
      (Sin(ThetaW())) + 2*Sqr(gN)*Sqr(Sin(ThetaWp()))) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(-14.696938456699067*g1*gN*Cos(ThetaWp())*Sin(ThetaW()
      )*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + Cos(
      ThetaW())*(18.973665961010276*g2*gN*Cos(ThetaWp())*Sin(ThetaWp()) -
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 6*Sqr(g1)*Sqr(
      Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)*(-9.486832980505138*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 9*Sqr(gN)*
      Sqr(Cos(ThetaWp())) - 15.491933384829668*g1*g2*Cos(ThetaW())*Sin(ThetaW())*
      Sqr(Sin(ThetaWp())) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*
      g1*Sin(ThetaW())*(2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())
      *Sqr(Sin(ThetaWp())))) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-
      4.898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 2*Sqr(
      gN)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) +
      3*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) + Cos(ThetaW())*(
      3.1622776601683795*g2*gN*Sin(2*ThetaWp()) - 7.745966692414834*g1*g2*Sin(
      ThetaW())*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUHpmconjHpmconjUHpm(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*(KroneckerDelta
      (0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*ZP(gI1,0)*
      ZP(gI2,1) + KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(
      g2) - 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0) - 2*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))
      *ZP(gI1,1)*ZP(gI2,1))) - KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*(20*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZP(gI1,1)*ZP(gI2,0) +
      KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gI1,0)*ZP(gI2
      ,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 3*Sqr(gN))*ZP(gI1,1)*ZP(
      gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSHp0conjUHpmconjSHp0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*(3*Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN)))*(Conj(UHp0(gI1,0
      ))*UHp0(gI2,0) - Conj(UHp0(gI1,1))*UHp0(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSHppconjUHpmconjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN)) - KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN)))*(Conj(UHpp(gI1,0)
      )*UHpp(gI2,0) - Conj(UHpp(gI1,1))*UHpp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSSI0conjUHpmconjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.125*(3*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + 2*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjUHpmconjSHp0(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = -0.35355339059327373*(vd*KroneckerDelta(0,
      gO2) + vu*KroneckerDelta(1,gO2))*Sqr(g2)*(Conj(UHpp(gI2,0))*UHp0(gI1,0) +
      Conj(UHpp(gI2,1))*UHp0(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(
      KroneckerDelta(1,gO2)*(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI2,0) + vd*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZA(gI2,1) - 2.8284271247461903*TLambdax*ZA(gI2,2)
      )*ZP(gI1,0) + KroneckerDelta(0,gO2)*(-(vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(
      gI2,0)) - vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI2,1) + 2.8284271247461903*
      Conj(TLambdax)*ZA(gI2,2))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO2)*(5*ZH(gI2,2)
      *((8*vs*AbsSqr(Lambdax) - 3*vs*Sqr(gN))*ZP(gI1,0) + 5.656854249492381*Conj(
      TLambdax)*ZP(gI1,1)) + ZH(gI2,1)*(2*vu*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*
      ZP(gI1,0) + 10*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) + ZH(gI2,0)*(vd*
      (6*Sqr(g1) + 10*Sqr(g2) + 9*Sqr(gN))*ZP(gI1,0) + 10*vu*(-2*AbsSqr(Lambdax) +
      Sqr(g2))*ZP(gI1,1)))) - 2*KroneckerDelta(1,gO2)*(5*ZH(gI2,2)*(
      2.8284271247461903*TLambdax*ZP(gI1,0) - vs*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZP
      (gI1,1)) + ZH(gI2,1)*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) + vu*(3*
      Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZP(gI1,1)) + ZH(gI2,0)*(5*vu*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gI1,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZP(
      gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*(-10*
      KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*
      Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZA(gI1
      ,1)*ZA(gI2,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZA(gI1,2)*ZA(gI2,2)))) + 2
      *KroneckerDelta(1,gO1)*(5*KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2
      ))*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*((3*
      Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2
      ) + 2*Sqr(gN))*ZA(gI1,1)*ZA(gI2,1) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZA(gI1
      ,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*(10*
      KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) +
      ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(g2) + 9*
      Sqr(gN))*ZH(gI1,0)*ZH(gI2,0) + 2*(-3*Sqr(g1) + 5*Sqr(g2) + 3*Sqr(gN))*ZH(gI1
      ,1)*ZH(gI2,1) + 5*(8*AbsSqr(Lambdax) - 3*Sqr(gN))*ZH(gI1,2)*ZH(gI2,2)))) + 2
      *KroneckerDelta(1,gO1)*(-5*KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(
      g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(1,gO2)*((3
      *Sqr(g1) - 5*Sqr(g2) - 3*Sqr(gN))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(
      g2) + 2*Sqr(gN))*ZH(gI1,1)*ZH(gI2,1) + 5*(-4*AbsSqr(Lambdax) + Sqr(gN))*ZH(
      gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSvconjUHpmconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(
      1,gO2)*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN)) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(KroneckerDelta(gI1,gI2)*(-3*Sqr
      (g1) + 5*Sqr(g2) + 3*Sqr(gN)) - 20*SUM(j2,0,2,AbsSqr(Ye(j2,j2))*Conj(ZV(gI1,
      j2))*ZV(gI2,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHpmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(Yd(j1
      ,j1))*ZDR(gI2,j1)*ZUL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHpmPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZDL(
      gI2,j1))*Conj(ZUR(gI1,j1))*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjUHpmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,Conj(Ye(gI1,gI1))*KroneckerDelta
      (0,gO2)*ZER(gI2,gI1),0);

   return result;
}

double CLASSNAME::CpbarFvFeconjUHpmPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUHpmconjSv(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(1.4142135623730951*KroneckerDelta(1,
      gO2)*(-(vu*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZV(gI1,j1))) + 2*vs*Lambdax*
      SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI2,3 + j1))*ZV(gI1,j1))) +
      KroneckerDelta(0,gO2)*(-1.4142135623730951*vd*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2
      ,j1))*ZV(gI1,j1)) + 4*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j1))*ZV(
      gI1,j1)) + 2.8284271247461903*vd*SUM(j2,0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(gI2,j2
      ))*ZV(gI1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSHI0conjUHpmconjSHI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.025*(KroneckerDelta(0,gO1)*(40*
      KroneckerDelta(1,gO2)*Lambdax*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj(
      Lambda12(j1,j1))*UHI0(gI2,j1)) - KroneckerDelta(0,gO2)*((6*Sqr(g1) - 10*Sqr(
      g2) + 9*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) + 2*(-3*Sqr(g1)
      + 5*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))
      )) + 2*KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2)
      - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) + (-3*Sqr(g1) + 5*
      Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))) +
      20*Conj(Lambdax)*KroneckerDelta(0,gO2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(
      gI2,2 + j1)*Lambda12(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSHIpconjUHpmconjSHIp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*(40*
      KroneckerDelta(1,gO2)*Lambdax*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*Conj(
      Lambda12(j1,j1))*UHIp(gI2,j1)) + KroneckerDelta(0,gO2)*((6*Sqr(g1) + 10*Sqr(
      g2) + 9*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) - 2*(3*Sqr(g1)
      + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))
      ))) + 2*KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2)
      - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) - (3*Sqr(g1) + 5*
      Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))) -
      20*Conj(Lambdax)*KroneckerDelta(0,gO2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(
      gI2,2 + j1)*Lambda12(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjUHpmconjSHI0(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = -0.35355339059327373*(vd*KroneckerDelta(0,
      gO2) + vu*KroneckerDelta(1,gO2))*Sqr(g2)*SUM(j1,0,3,Conj(UHIp(gI2,j1))*UHI0(
      gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSdconjUHpmconjSd(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((-Sqr(g1) + 5*Sqr(g2) + Sqr(gN))*SUM(j1,0,2,Conj(ZD(
      gI1,j1))*ZD(gI2,j1)) - 2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))
      *ZD(gI2,3 + j1)) - 20*SUM(j2,0,2,AbsSqr(Yu(j2,j2))*Conj(ZD(gI1,j2))*ZD(gI2,
      j2))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((2*Sqr(g1) - 10*Sqr(g2)
      + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + (4*Sqr(g1) + 6*Sqr(gN
      ))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) - 40*SUM(j2,0,2,AbsSqr(Yd
      (j2,j2))*Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSDXconjUHpmconjSDX(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.025*(-(KroneckerDelta(0,gO1)*(40*
      KroneckerDelta(1,gO2)*Lambdax*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*Conj(Kappa(j1
      ,j1))*ZDX(gI2,j1)) + KroneckerDelta(0,gO2)*((4*Sqr(g1) + 6*Sqr(gN))*SUM(j1,0
      ,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) + (-4*Sqr(g1) + 9*Sqr(gN))*SUM(j1,0,2,Conj
      (ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1))))) + 2*KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*(2*(Sqr(g1) - Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      ZDX(gI2,j1)) - (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(
      gI2,3 + j1))) - 20*Conj(Lambdax)*KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(ZDX(
      gI1,j1))*ZDX(gI2,3 + j1)*Kappa(j1,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSeconjUHpmconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.025*(2*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(
      ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,3 +
      j1))*ZE(gI2,3 + j1))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-2*(3*
      Sqr(g1) + 5*Sqr(g2) - 3*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 3
      *(4*Sqr(g1) + Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) - 40*
      SUM(j2,0,2,AbsSqr(Ye(j2,j2))*Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSuconjUHpmconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.025*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*((2*Sqr(g1) + 10*Sqr(g2) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU
      (gI2,j1)) + (-8*Sqr(g1) + 3*Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,
      3 + j1)) - 40*SUM(j2,0,2,AbsSqr(Yd(j2,j2))*Conj(ZU(gI1,j2))*ZU(gI2,j2))) + 2
      *KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-Sqr(g1) - 5*Sqr(g2) + Sqr(gN
      ))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (4*Sqr(g1) + Sqr(gN))*SUM(j1,0,
      2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) - 20*SUM(j2,0,2,AbsSqr(Yu(j2,j2))*
      Conj(ZU(gI1,3 + j2))*ZU(gI2,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjUHpmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(Conj(Lambdax)*KroneckerDelta(0,gO2)*UP(
      gI2,1)*ZN(gI1,4)) - 0.1*KroneckerDelta(1,gO2)*(10*g2*UP(gI2,0)*ZN(gI1,3) +
      UP(gI2,1)*(5.477225575051661*g1*ZN(gI1,0) + 7.0710678118654755*g2*ZN(gI1,1)
      - 4.47213595499958*gN*ZN(gI1,5)));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjUHpmPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Conj(ZN(gI1,2))*
      KroneckerDelta(0,gO1)) + 0.1*Conj(UM(gI2,1))*(5.477225575051661*g1*Conj(ZN(
      gI1,0))*KroneckerDelta(0,gO1) + 7.0710678118654755*g2*Conj(ZN(gI1,1))*
      KroneckerDelta(0,gO1) + 6.708203932499369*gN*Conj(ZN(gI1,5))*KroneckerDelta(
      0,gO1) - 10*Conj(ZN(gI1,4))*KroneckerDelta(1,gO1)*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUHpmconjSu(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(KroneckerDelta(1,gO2)*(-
      1.4142135623730951*vu*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 2*(
      1.4142135623730951*vs*Lambdax*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI2,3 + j1)
      )*ZU(gI1,j1)) + 2*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,3 + j1)*TYu(j1,j1)) +
      1.4142135623730951*vu*SUM(j2,0,2,AbsSqr(Yu(j2,j2))*Conj(ZD(gI2,j2))*ZU(gI1,
      j2)) + 1.4142135623730951*vd*SUM(j2,0,2,Conj(Yd(j2,j2))*Conj(ZD(gI2,3 + j2))
      *Yu(j2,j2)*ZU(gI1,3 + j2)))) + KroneckerDelta(0,gO2)*(-1.4142135623730951*vd
      *Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 4*SUM(j1,0,2,Conj(ZD(gI2,
      3 + j1))*Conj(TYd(j1,j1))*ZU(gI1,j1)) + 2.8284271247461903*(vs*Conj(Lambdax)
      *SUM(j1,0,2,Conj(ZD(gI2,j1))*Yu(j1,j1)*ZU(gI1,3 + j1)) + vd*SUM(j2,0,2,
      AbsSqr(Yd(j2,j2))*Conj(ZD(gI2,j2))*ZU(gI1,j2)) + vu*SUM(j2,0,2,Conj(Yd(j2,j2
      ))*Conj(ZD(gI2,3 + j2))*Yu(j2,j2)*ZU(gI1,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.3872983346207417*g1*Cos(
      ThetaW())*ZP(gI2,gO2),0) + IF(gI2 < 2,-0.5*g2*Sin(ThetaW())*ZP(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO2)*(-10*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())
      - 9.486832980505138*gN*Sin(ThetaWp()))*ZP(gI2,0) + 2*KroneckerDelta(1,gO2)*(
      -5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin
      (ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO2)*(-
      9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZP(gI2,0) + 2*
      KroneckerDelta(1,gO2)*(3.1622776601683795*gN*Cos(ThetaWp()) + (5*g2*Cos(
      ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjUHpmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZA(gI2,0) + KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhconjUHpmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) -
      KroneckerDelta(1,gO2)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0conjUSHI0VZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.9486832980505138*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 2,-0.7348469228349533*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 2,0.7745966692414834*g1*g2*Cos(ThetaW
      ())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(gO1 <
      2,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()))
      ,0) + IF(gO1 < 2,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr
      (Sin(ThetaW())),0) + IF(gO1 < 2,0.45*KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(Sin
      (ThetaWp())),0) + 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))
      + 0.31622776601683794*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp())*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.7745966692414834*
      g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Cos(ThetaWp()))*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.5*Sqr(g2)*Sqr(Cos
      (ThetaW()))*Sqr(Cos(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1)) + 0.3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW
      ()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.2
      *Sqr(gN)*Sqr(Sin(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0conjUSHI0VZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.9486832980505138*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 2,0.7348469228349533*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 2,0.45*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 2,0.7745966692414834*g1*g2*Cos(ThetaW())
      *KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 2,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0
      ) + IF(gO1 < 2,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp())),0) - 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1)) - 0.31622776601683794*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp())*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.2*Sqr(gN)*Sqr(Cos
      (ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      )) + 0.7745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())
      )*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.5*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1
      ,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.3*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1));

   return result;
}

double CLASSNAME::CpUSHI0conjUSHI0conjVWmVWm(int gO1, int gO2) const
{
   
   const double result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSHI0conjHpmconjUSHI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 2,0.25*KroneckerDelta(gO1,gO2)*Sqr
      (g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 2,-0.225*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*
      ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZP
      (gI1,1)*ZP(gI2,1),0) + IF(gO1 < 2,Conj(Lambdax)*KroneckerDelta(2 + gO1,gO2)*
      ZP(gI1,1)*ZP(gI2,0)*Lambda12(gO1,gO1),0) + IF(gO2 < 2,Conj(Lambda12(gO2,gO2)
      )*KroneckerDelta(gO1,2 + gO2)*Lambdax*ZP(gI1,0)*ZP(gI2,1),0) + 0.15*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZP(gI1,0)*
      ZP(gI2,0) - 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZP(gI1,0)*ZP(gI2,0) - 0.15*Sqr(gN)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZP(gI1,0)*ZP(gI2,0) -
      0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*ZP(gI1,1)*ZP(gI2,1) + 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZP(gI1,1)*ZP(gI2,1) - 0.1*Sqr(gN)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0SHp0conjUSHI0conjSHp0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*Conj(UHp0(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,0),0) + IF(gO1 < 2,-0.25*Conj(UHp0(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,0),0) + IF(gO1 < 2,0.15*
      Conj(UHp0(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,0),0) + IF(gO1 <
      2,0.15*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,1),0) + IF
      (gO1 < 2,0.25*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,1),
      0) + IF(gO1 < 2,-0.15*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0
      (gI2,1),0) + 0.15*Conj(UHp0(gI1,0))*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2
      + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,0) + 0.25*Conj(UHp0(gI1,0))*Sqr(
      g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(
      gI2,0) + 0.1*Conj(UHp0(gI1,0))*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*UHp0(gI2,0) - 0.15*Conj(UHp0(gI1,1))*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,1
      ) - 0.25*Conj(UHp0(gI1,1))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*UHp0(gI2,1) - 0.1*Conj(UHp0(gI1,1))*Sqr(gN)*SUM(
      j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0SHppconjUSHI0conjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*Conj(UHpp(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,0),0) + IF(gO1 < 2,0.25*Conj(UHpp(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,0),0) + IF(gO1 < 2,0.15*
      Conj(UHpp(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,0),0) + IF(gO1 <
      2,0.15*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,1),0) + IF
      (gO1 < 2,-0.25*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,1)
      ,0) + IF(gO1 < 2,-0.15*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      UHpp(gI2,1),0) + 0.15*Conj(UHpp(gI1,0))*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(
      gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,0) - 0.25*Conj(UHpp(gI1,0))
      *Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      UHpp(gI2,0) + 0.1*Conj(UHpp(gI1,0))*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2
      + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,0) - 0.15*Conj(UHpp(gI1,1))*Sqr(
      g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(
      gI2,1) + 0.25*Conj(UHpp(gI1,1))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1
      )*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,1) - 0.1*Conj(UHpp(gI1,1))*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,1
      );

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0SSI0conjUSHI0conjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.375*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.25*KroneckerDelta(gI1,gI2)*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaIconjUSHI0PR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,1,KroneckerDelta(gO2,2 + j1)*
      ZPI(gI2,j1))*UM(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaIconjUSHI0PL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-(g2*Conj(UP(gI1,0))*Conj(ZMI(
      gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjHpmconjUSHI0(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 4,-0.35355339059327373*vd*Conj(
      UHIp(gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 4,-0.35355339059327373*vu*
      Conj(UHIp(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSHI0conjUSHI0(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 2,-0.225*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 2,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 2,-(AbsSqr(Lambda12(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZA(gI1,2)*ZA(gI2,2)),0) + IF(gO1 < 2,0.375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2),0) + IF(gO1 < 2,0.5*Conj
      (Lambdax)*KroneckerDelta(2 + gO1,gO2)*ZA(gI1,1)*ZA(gI2,0)*Lambda12(gO1,gO1),
      0) + IF(gO1 < 2,0.5*Conj(Lambdax)*KroneckerDelta(2 + gO1,gO2)*ZA(gI1,0)*ZA(
      gI2,1)*Lambda12(gO1,gO1),0) + IF(gO2 < 2,0.5*Conj(Lambda12(gO2,gO2))*
      KroneckerDelta(gO1,2 + gO2)*Lambdax*ZA(gI1,1)*ZA(gI2,0),0) + IF(gO2 < 2,0.5*
      Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*Lambdax*ZA(gI1,0)*ZA(gI2
      ,1),0) + 0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(
      gO2,2 + j1))*ZA(gI1,0)*ZA(gI2,0) + 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(
      gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZA(gI1,0)*ZA(gI2,0) - 0.15*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZA(gI1,0)*
      ZA(gI2,0) - 0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZA(gI1,1)*ZA(gI2,1) - 0.25*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZA(gI1,1)*ZA(gI2,1) -
      0.1*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)
      )*ZA(gI1,1)*ZA(gI2,1) + 0.25*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZA(gI1,2)*ZA(gI2,2) - SUM(j2,0,1,AbsSqr(Lambda12
      (j2,j2))*KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2))*ZA(gI1,2)*ZA
      (gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSHI0conjUSHI0(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 2,-0.225*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 2,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 2,-(AbsSqr(Lambda12(gO1,gO1))*
      KroneckerDelta(gO1,gO2)*ZH(gI1,2)*ZH(gI2,2)),0) + IF(gO1 < 2,0.375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2),0) + IF(gO1 < 2,-0.5*
      Conj(Lambdax)*KroneckerDelta(2 + gO1,gO2)*ZH(gI1,1)*ZH(gI2,0)*Lambda12(gO1,
      gO1),0) + IF(gO1 < 2,-0.5*Conj(Lambdax)*KroneckerDelta(2 + gO1,gO2)*ZH(gI1,0
      )*ZH(gI2,1)*Lambda12(gO1,gO1),0) + IF(gO2 < 2,-0.5*Conj(Lambda12(gO2,gO2))*
      KroneckerDelta(gO1,2 + gO2)*Lambdax*ZH(gI1,1)*ZH(gI2,0),0) + IF(gO2 < 2,-0.5
      *Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*Lambdax*ZH(gI1,0)*ZH(
      gI2,1),0) + 0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZH(gI1,0)*ZH(gI2,0) + 0.25*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZH(gI1,0)*ZH(gI2,0) -
      0.15*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*ZH(gI1,0)*ZH(gI2,0) - 0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZH(gI1,1)*ZH(gI2,1) - 0.25*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZH(gI1,1)*ZH(gI2,1) -
      0.1*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)
      )*ZH(gI1,1)*ZH(gI2,1) + 0.25*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZH(gI1,2)*ZH(gI2,2) - SUM(j2,0,1,AbsSqr(Lambda12
      (j2,j2))*KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2))*ZH(gI1,2)*ZH
      (gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0SvconjUSHI0conjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gI1,gI2
      )*KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 2,0.15*KroneckerDelta(gI1,
      gI2)*KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.15*KroneckerDelta(gI1,gI2)*Sqr(
      g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.25
      *KroneckerDelta(gI1,gI2)*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1)) + 0.1*KroneckerDelta(gI1,gI2)*Sqr(gN)*SUM(j1,0,1
      ,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USHI0conjSHI0conjUSHI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,IF(gO2 < 2,-0.15*Conj(UHI0(gI1,
      gO2))*Sqr(g1)*UHI0(gI2,gO1),0),0) + IF(gO1 < 2,IF(gO2 < 2,-0.25*Conj(UHI0(
      gI1,gO2))*Sqr(g2)*UHI0(gI2,gO1),0),0) + IF(gO1 < 2,IF(gO2 < 2,-0.225*Conj(
      UHI0(gI1,gO2))*Sqr(gN)*UHI0(gI2,gO1),0),0) + IF(gO1 < 2,IF(gO2 < 2,-(Conj(
      UHI0(gI1,2 + gO2))*Conj(Lambda12(gO2,gO2))*UHI0(gI2,2 + gO1)*Lambda12(gO1,
      gO1)),0),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,1,
      Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,
      gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 2,-
      0.1125*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(
      gI2,j1)),0) + IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,1,
      Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2
       + j1)),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),
      0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(
      gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 2,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + IF(gO1 < 2,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + IF(gO1 < 2,0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*UHI0(gI2,gO1),0) + IF(gO1 < 2,0.125*Sqr(g2)*SUM(
      j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*UHI0(gI2,gO1),0) +
      IF(gO1 < 2,-0.075*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(
      gO2,2 + j1))*UHI0(gI2,gO1),0) + IF(gO1 < 2,0.075*Sqr(g1)*SUM(j2,0,1,Conj(
      UHI0(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*UHI0(gI2,gO1),0) + IF(gO1 < 2,
      0.125*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*
      UHI0(gI2,gO1),0) + IF(gO1 < 2,-0.075*Sqr(gN)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2
      ))*KroneckerDelta(gO2,2 + j2))*UHI0(gI2,gO1),0) + IF(gO1 < 2,-(
      KroneckerDelta(2 + gO1,gO2)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*Conj(Lambda12(
      j2,j2))*UHI0(gI2,j2))*Lambda12(gO1,gO1)),0) + IF(gO2 < 2,0.075*Conj(UHI0(gI1
      ,gO2))*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHI0(gI2,2 + j1)),0) +
      IF(gO2 < 2,0.125*Conj(UHI0(gI1,gO2))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2
       + j1)*UHI0(gI2,2 + j1)),0) + IF(gO2 < 2,-0.075*Conj(UHI0(gI1,gO2))*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHI0(gI2,2 + j1)),0) + IF(gO2 < 2,-(
      Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*SUM(j1,0,1,Conj(UHI0(gI1
      ,j1))*UHI0(gI2,2 + j1)*Lambda12(j1,j1))),0) + IF(gO2 < 2,0.075*Conj(UHI0(gI1
      ,gO2))*Sqr(g1)*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHI0(gI2,2 + j2)),0) +
      IF(gO2 < 2,0.125*Conj(UHI0(gI1,gO2))*Sqr(g2)*SUM(j2,0,1,KroneckerDelta(gO1,2
       + j2)*UHI0(gI2,2 + j2)),0) + IF(gO2 < 2,-0.075*Conj(UHI0(gI1,gO2))*Sqr(gN)*
      SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHI0(gI2,2 + j2)),0) - 0.075*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHI0(gI2,2 + j1))*SUM(j2,0,1,Conj(UHI0
      (gI1,2 + j2))*KroneckerDelta(gO2,2 + j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*UHI0(gI2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2
      ))*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,
      2 + j1)*UHI0(gI2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*KroneckerDelta(
      gO2,2 + j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM
      (j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.125*Sqr(
      g2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,1,Conj(UHI0
      (gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta
      (gO2,2 + j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 +
      j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) -
      0.125*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,1
      ,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(
      j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(
      gI1,j2))*UHI0(gI2,j2)) + 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) -
      0.075*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) - SUM(j1,0,1,Conj(UHI0(gI1,
      j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(Lambda12(j2
      ,j2))*KroneckerDelta(gO1,2 + j2)*UHI0(gI2,j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(
      gI1,2 + j2))*UHI0(gI2,2 + j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1
      ,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(
      gI2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 +
      j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2
      + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHI0(gI2,2 + j2)) - 0.125*Sqr(
      g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1
      ,KroneckerDelta(gO1,2 + j2)*UHI0(gI2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,Conj
      (UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,KroneckerDelta(gO1
      ,2 + j2)*UHI0(gI2,2 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0SHIpconjUSHI0conjSHIp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,IF(gO2 < 2,-0.5*Conj(UHIp(gI1,
      gO2))*Sqr(g2)*UHIp(gI2,gO1),0),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),
      0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHIp
      (gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 2,-0.125
      *KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,
      2 + j1)),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),
      0) + IF(gO1 < 2,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(
      gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 2,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + IF(gO1 < 2,-0.25*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*UHIp(gI2,gO1),0) + IF(gO1 < 2,-0.25*Sqr(g2)*SUM(
      j2,0,1,Conj(UHIp(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*UHIp(gI2,gO1),0) +
      IF(gO1 < 2,KroneckerDelta(2 + gO1,gO2)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*
      Conj(Lambda12(j2,j2))*UHIp(gI2,j2))*Lambda12(gO1,gO1),0) + IF(gO2 < 2,-0.25*
      Conj(UHIp(gI1,gO2))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHIp(gI2,2
       + j1)),0) + IF(gO2 < 2,Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*
      SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,2 + j1)*Lambda12(j1,j1)),0) + IF(gO2
      < 2,-0.25*Conj(UHIp(gI1,gO2))*Sqr(g2)*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*
      UHIp(gI2,2 + j2)),0) - 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      UHIp(gI2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*KroneckerDelta(gO2,2 +
      j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.125*Sqr(g2)*SUM(
      j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2
      )*KroneckerDelta(gO2,2 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*
      UHIp(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM
      (j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.125*Sqr(
      g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1
      ,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,
      2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(
      gI1,j2))*UHIp(gI2,j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) -
      0.075*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(
      gI1,2 + j2))*UHIp(gI2,2 + j2)) + 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1
      ,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(
      gI2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 +
      j2)) - 0.25*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHIp(gI2,2 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpAhSHI0conjUSHI0(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,std::complex<double>(0,-0.5)*vu*
      Conj(UHI0(gI1,2 + gO2))*Conj(Lambda12(gO2,gO2))*Lambdax*ZA(gI2,0),0) + IF(
      gO2 < 2,std::complex<double>(0,-0.5)*vd*Conj(UHI0(gI1,2 + gO2))*Conj(
      Lambda12(gO2,gO2))*Lambdax*ZA(gI2,1),0) + IF(gO2 < 2,std::complex<double>(0.
      ,-0.7071067811865475)*Conj(UHI0(gI1,2 + gO2))*Conj(TLambda12(gO2,gO2))*ZA(
      gI2,2),0) + std::complex<double>(0,0.5)*vu*Conj(Lambdax)*SUM(j1,0,1,Conj(
      UHI0(gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZA(gI2,0) + std::
      complex<double>(0,0.5)*vd*Conj(Lambdax)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*
      KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZA(gI2,1) + std::complex<double>
      (0.,0.7071067811865475)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*KroneckerDelta(gO2,2 +
      j1)*TLambda12(j1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSHI0conjUSHI0(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.5*vu*Conj(UHI0(gI1,2 + gO2))*
      Conj(Lambda12(gO2,gO2))*Lambdax*ZH(gI2,0),0) + IF(gO2 < 2,-0.15*vd*Conj(UHI0
      (gI1,gO2))*Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 2,-0.25*vd*Conj(UHI0(gI1,gO2))*
      Sqr(g2)*ZH(gI2,0),0) + IF(gO2 < 2,-0.225*vd*Conj(UHI0(gI1,gO2))*Sqr(gN)*ZH(
      gI2,0),0) + IF(gO2 < 2,-0.5*vd*Conj(UHI0(gI1,2 + gO2))*Conj(Lambda12(gO2,gO2
      ))*Lambdax*ZH(gI2,1),0) + IF(gO2 < 2,0.15*vu*Conj(UHI0(gI1,gO2))*Sqr(g1)*ZH(
      gI2,1),0) + IF(gO2 < 2,0.25*vu*Conj(UHI0(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF
      (gO2 < 2,-0.15*vu*Conj(UHI0(gI1,gO2))*Sqr(gN)*ZH(gI2,1),0) + IF(gO2 < 2,-(vs
      *AbsSqr(Lambda12(gO2,gO2))*Conj(UHI0(gI1,gO2))*ZH(gI2,2)),0) + IF(gO2 < 2,
      0.7071067811865475*Conj(UHI0(gI1,2 + gO2))*Conj(TLambda12(gO2,gO2))*ZH(gI2,2
      ),0) + IF(gO2 < 2,0.375*vs*Conj(UHI0(gI1,gO2))*Sqr(gN)*ZH(gI2,2),0) + 0.15*
      vd*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*ZH(
      gI2,0) + 0.25*vd*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(
      gO2,2 + j1))*ZH(gI2,0) - 0.15*vd*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*ZH(gI2,0) - 0.5*vu*Conj(Lambdax)*SUM(j1,0,1,Conj
      (UHI0(gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZH(gI2,0) - 0.15*
      vu*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*ZH(
      gI2,1) - 0.25*vu*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(
      gO2,2 + j1))*ZH(gI2,1) - 0.1*vu*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*ZH(gI2,1) - 0.5*vd*Conj(Lambdax)*SUM(j1,0,1,Conj
      (UHI0(gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZH(gI2,1) + 0.25*
      vs*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*ZH(
      gI2,2) + 0.7071067811865475*SUM(j1,0,1,Conj(UHI0(gI1,j1))*KroneckerDelta(gO2
      ,2 + j1)*TLambda12(j1,j1))*ZH(gI2,2) - vs*SUM(j2,0,1,AbsSqr(Lambda12(j2,j2))
      *Conj(UHI0(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpChiChiIconjUSHI0PR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,Conj(Lambda12(gO2,gO2))*ZN(gI2,4
      )*ZNI(gI1,2 + gO2),0) + 0.1*SUM(j1,0,1,KroneckerDelta(gO2,2 + j1)*ZNI(gI1,2
      + j1))*(-5.477225575051661*g1*ZN(gI2,0) + 7.0710678118654755*g2*ZN(gI2,1) +
      4.47213595499958*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiIconjUSHI0PL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.5477225575051661*g1*Conj(ZN(
      gI2,0))*Conj(ZNI(gI1,gO1)),0) + IF(gO1 < 2,-0.7071067811865475*g2*Conj(ZN(
      gI2,1))*Conj(ZNI(gI1,gO1)),0) + IF(gO1 < 2,0.6708203932499369*gN*Conj(ZN(gI2
      ,5))*Conj(ZNI(gI1,gO1)),0) + Conj(ZN(gI2,4))*SUM(j1,0,1,Conj(ZNI(gI1,j1))*
      KroneckerDelta(gO1,2 + j1)*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSHI0conjSdconjUSHI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) +
      0.025*(40*IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(
      ZD(gI1,j1))*ZD(gI2,j1)),0) + 40*IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)),0) + 40*IF(gO1 < 2,0.075
      *KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 +
      j1)),0) + 40*IF(gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + 40*IF(gO1 < 2,0.125*KroneckerDelta(gO1,gO2
      )*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + 40*IF(gO1 < 2,0.0375*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) +
      40*IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI1,3
      + j2))*ZD(gI2,3 + j2)),0) + 40*IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2)),0) - Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*
      KroneckerDelta(gO2,2 + j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,
      j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + Sqr
      (gN)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2
       + j2)*KroneckerDelta(gO2,2 + j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1
      ))*ZD(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,
      2 + j2)) + 2*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))*SUM(j2,
      0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - Sqr(g1)*SUM(j1,
      0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(
      ZD(gI1,j2))*ZD(gI2,j2)) - 5*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) + Sqr(gN
      )*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0
      ,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) - 2*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2
      + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 +
      j2)) + 2*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2
      + j1))*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSHI0conjSDXconjUSHI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)),0)
      + IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,3
      + j1))*ZDX(gI2,3 + j1)),0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)),0) + IF(gO1 < 2,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)),0)
      + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI1,
      j2))*ZDX(gI2,j2)),0) + IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 + j2)),0) + IF(gO1 < 2,-0.1125*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 +
      j2)),0) + IF(gO1 < 2,KroneckerDelta(2 + gO1,gO2)*SUM(j2,0,2,Conj(ZDX(gI1,3 +
      j2))*Conj(Kappa(j2,j2))*ZDX(gI2,j2))*Lambda12(gO1,gO1),0) + IF(gO2 < 2,Conj(
      Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*
      ZDX(gI2,3 + j1)*Kappa(j1,j1)),0) + 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,j1))
      *ZDX(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) - 0.05*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(g1)*SUM(j1
      ,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2
      + j2)*KroneckerDelta(gO2,2 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,3
      + j1))*ZDX(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta
      (gO2,2 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) - 0.05
      *Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) - 0.05*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZDX(
      gI1,3 + j2))*ZDX(gI2,3 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,
      2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2
      ,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSHI0conjSeconjUSHI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      0.025*(40*IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(
      ZE(gI1,j1))*ZE(gI2,j1)),0) + 40*IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + 40*IF(gO1 < 2,
      0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(
      gI2,3 + j1)),0) + 40*IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + 40*IF(gO1 < 2,0.125*KroneckerDelta(
      gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + 40*IF(gO1 < 2,
      0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)
      ),0) + 40*IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE
      (gI1,3 + j2))*ZE(gI2,3 + j2)),0) + 40*IF(gO1 < 2,0.0375*KroneckerDelta(gO1,
      gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + 3*Sqr(g1)*
      SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2
      )*KroneckerDelta(gO2,2 + j2)) - 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2
      ,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 2*
      Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3
       + j1))*ZE(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta
      (gO2,2 + j2)) + Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(
      j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 3*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2
      ,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 5*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 +
      j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) + 2*
      Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 6*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)) + Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2))
      );

   return result;
}

std::complex<double> CLASSNAME::CpUSHI0SuconjUSHI0conjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 2,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1)),0) + IF(gO1 < 2,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,
      2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 2,0.0375*KroneckerDelta
      (gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(
      gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU
      (gI2,j2)),0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 2,-0.1*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 +
      j2)),0) + IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(
      ZU(gI1,3 + j2))*ZU(gI2,3 + j2)),0) - 0.025*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1
      ))*ZU(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) + 0.125*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.025*Sqr(gN)*SUM(
      j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*
      KroneckerDelta(gO2,2 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU
      (gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) + 0.025*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,
      0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.025*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2
      ,Conj(ZU(gI1,j2))*ZU(gI2,j2)) + 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,
      2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2))
      + 0.025*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)) + 0.1*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZU(
      gI1,3 + j2))*ZU(gI2,3 + j2)) + 0.025*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2
       + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3
      + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjUSHI0VZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.5*g2*Conj(UHI0(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 2,0.3872983346207417*g1*Conj(UHI0(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 2,-0.4743416490252569*gN*
      Conj(UHI0(gI2,gO2))*Sin(ThetaWp()),0) + 0.5*g2*Cos(ThetaW())*Cos(ThetaWp())*
      SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1)) +
      0.3872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*SUM(j1,0,1,Conj(UHI0(gI2,
      2 + j1))*KroneckerDelta(gO2,2 + j1)) + 0.31622776601683794*gN*Sin(ThetaWp())
      *SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjUSHI0VZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.4743416490252569*gN*Conj(UHI0
      (gI2,gO2))*Cos(ThetaWp()),0) + IF(gO2 < 2,-0.5*g2*Conj(UHI0(gI2,gO2))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO2 < 2,-0.3872983346207417*g1*Conj(UHI0(
      gI2,gO2))*Sin(ThetaW())*Sin(ThetaWp()),0) + 0.31622776601683794*gN*Cos(
      ThetaWp())*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1)) -
      0.5*g2*Cos(ThetaW())*Sin(ThetaWp())*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*
      KroneckerDelta(gO2,2 + j1)) - 0.3872983346207417*g1*Sin(ThetaW())*Sin(
      ThetaWp())*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjUSHI0conjVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.7071067811865475*g2*Conj(UHIp(
      gI2,gO2)),0) - 0.7071067811865475*g2*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*
      KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpconjUSHIpVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.9486832980505138*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 2,-0.7348469228349533*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 2,-0.7745966692414834*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(
      gO1 < 2,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())),0) + IF(gO1 < 2,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())),0) + IF(gO1 < 2,0.45*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*Sqr(Sin(ThetaWp())),0) - 0.6324555320336759*g2*gN*Cos(ThetaW())*Cos(
      ThetaWp())*Sin(ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1)) + 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(
      ThetaW())*Sin(ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1)) - 0.7745966692414834*g1*g2*Cos(ThetaW())*Sin(
      ThetaW())*Sqr(Cos(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1)) + 0.5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp
      ()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.3
      *Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))*SUM(j1,0,1,KroneckerDelta(
      gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.2*Sqr(gN)*Sqr(Sin(ThetaWp()))*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpconjUSHIpVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.9486832980505138*g2*gN*Cos(
      ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1
      < 2,0.7348469228349533*g1*gN*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 2,0.45*KroneckerDelta(gO1,gO2)*Sqr(gN
      )*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 2,-0.7745966692414834*g1*g2*Cos(ThetaW()
      )*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 2,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0
      ) + IF(gO1 < 2,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp())),0) - 0.4898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1)) + 0.31622776601683794*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp())*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.2*Sqr(gN)*Sqr(Cos
      (ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      )) - 0.7745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())
      )*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.5*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1
      ,2 + j1)*KroneckerDelta(gO2,2 + j1)) + 0.3*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp()))*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1));

   return result;
}

double CLASSNAME::CpUSHIpconjUSHIpconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSHIpconjHpmconjUSHIp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 2,-0.225*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 2,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 2,-(Conj(Lambdax)*KroneckerDelta(2 + gO1,
      gO2)*ZP(gI1,1)*ZP(gI2,0)*Lambda12(gO1,gO1)),0) + IF(gO2 < 2,-(Conj(Lambda12(
      gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*Lambdax*ZP(gI1,0)*ZP(gI2,1)),0) + 0.15
      *Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      ZP(gI1,0)*ZP(gI2,0) + 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZP(gI1,0)*ZP(gI2,0) - 0.15*Sqr(gN)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZP(gI1,0)*ZP(gI2,0) -
      0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*ZP(gI1,1)*ZP(gI2,1) - 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZP(gI1,1)*ZP(gI2,1) - 0.1*Sqr(gN)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpSHp0conjUSHIpconjSHp0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*Conj(UHp0(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,0),0) + IF(gO1 < 2,0.25*Conj(UHp0(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,0),0) + IF(gO1 < 2,0.15*
      Conj(UHp0(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHp0(gI2,0),0) + IF(gO1 <
      2,0.15*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHp0(gI2,1),0) + IF
      (gO1 < 2,-0.25*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHp0(gI2,1)
      ,0) + IF(gO1 < 2,-0.15*Conj(UHp0(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*
      UHp0(gI2,1),0) + 0.15*Conj(UHp0(gI1,0))*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(
      gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,0) - 0.25*Conj(UHp0(gI1,0))
      *Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      UHp0(gI2,0) + 0.1*Conj(UHp0(gI1,0))*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2
      + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,0) - 0.15*Conj(UHp0(gI1,1))*Sqr(
      g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(
      gI2,1) + 0.25*Conj(UHp0(gI1,1))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1
      )*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,1) - 0.1*Conj(UHp0(gI1,1))*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHp0(gI2,1
      );

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpSHppconjUSHIpconjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*Conj(UHpp(gI1,0))*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,0),0) + IF(gO1 < 2,-0.25*Conj(UHpp(
      gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,0),0) + IF(gO1 < 2,0.15*
      Conj(UHpp(gI1,0))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp(gI2,0),0) + IF(gO1 <
      2,0.15*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*UHpp(gI2,1),0) + IF
      (gO1 < 2,0.25*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*UHpp(gI2,1),
      0) + IF(gO1 < 2,-0.15*Conj(UHpp(gI1,1))*KroneckerDelta(gO1,gO2)*Sqr(gN)*UHpp
      (gI2,1),0) + 0.15*Conj(UHpp(gI1,0))*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2
      + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,0) + 0.25*Conj(UHpp(gI1,0))*Sqr(
      g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(
      gI2,0) + 0.1*Conj(UHpp(gI1,0))*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*UHpp(gI2,0) - 0.15*Conj(UHpp(gI1,1))*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,1
      ) - 0.25*Conj(UHpp(gI1,1))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*UHpp(gI2,1) - 0.1*Conj(UHpp(gI1,1))*Sqr(gN)*SUM(
      j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*UHpp(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpSSI0conjUSHIpconjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.375*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.25*KroneckerDelta(gI1,gI2)*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSHIpconjUSHIp(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 2,0.25*KroneckerDelta(gO1,gO2)*Sqr
      (g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 2,-0.225*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*
      ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZA
      (gI1,1)*ZA(gI2,1),0) + IF(gO1 < 2,-(AbsSqr(Lambda12(gO1,gO1))*KroneckerDelta
      (gO1,gO2)*ZA(gI1,2)*ZA(gI2,2)),0) + IF(gO1 < 2,0.375*KroneckerDelta(gO1,gO2)
      *Sqr(gN)*ZA(gI1,2)*ZA(gI2,2),0) + IF(gO1 < 2,-0.5*Conj(Lambdax)*
      KroneckerDelta(2 + gO1,gO2)*ZA(gI1,1)*ZA(gI2,0)*Lambda12(gO1,gO1),0) + IF(
      gO1 < 2,-0.5*Conj(Lambdax)*KroneckerDelta(2 + gO1,gO2)*ZA(gI1,0)*ZA(gI2,1)*
      Lambda12(gO1,gO1),0) + IF(gO2 < 2,-0.5*Conj(Lambda12(gO2,gO2))*
      KroneckerDelta(gO1,2 + gO2)*Lambdax*ZA(gI1,1)*ZA(gI2,0),0) + IF(gO2 < 2,-0.5
      *Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*Lambdax*ZA(gI1,0)*ZA(
      gI2,1),0) + 0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZA(gI1,0)*ZA(gI2,0) - 0.25*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZA(gI1,0)*ZA(gI2,0) -
      0.15*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*ZA(gI1,0)*ZA(gI2,0) - 0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZA(gI1,1)*ZA(gI2,1) + 0.25*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZA(gI1,1)*ZA(gI2,1) -
      0.1*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)
      )*ZA(gI1,1)*ZA(gI2,1) + 0.25*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZA(gI1,2)*ZA(gI2,2) - SUM(j2,0,1,AbsSqr(Lambda12
      (j2,j2))*KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2))*ZA(gI1,2)*ZA
      (gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSHIpconjUSHIp(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 2,0.25*KroneckerDelta(gO1,gO2)*Sqr
      (g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 2,-0.225*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 2,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*
      ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 2,-0.15*KroneckerDelta(gO1,gO2)*Sqr(gN)*ZH
      (gI1,1)*ZH(gI2,1),0) + IF(gO1 < 2,-(AbsSqr(Lambda12(gO1,gO1))*KroneckerDelta
      (gO1,gO2)*ZH(gI1,2)*ZH(gI2,2)),0) + IF(gO1 < 2,0.375*KroneckerDelta(gO1,gO2)
      *Sqr(gN)*ZH(gI1,2)*ZH(gI2,2),0) + IF(gO1 < 2,0.5*Conj(Lambdax)*
      KroneckerDelta(2 + gO1,gO2)*ZH(gI1,1)*ZH(gI2,0)*Lambda12(gO1,gO1),0) + IF(
      gO1 < 2,0.5*Conj(Lambdax)*KroneckerDelta(2 + gO1,gO2)*ZH(gI1,0)*ZH(gI2,1)*
      Lambda12(gO1,gO1),0) + IF(gO2 < 2,0.5*Conj(Lambda12(gO2,gO2))*KroneckerDelta
      (gO1,2 + gO2)*Lambdax*ZH(gI1,1)*ZH(gI2,0),0) + IF(gO2 < 2,0.5*Conj(Lambda12(
      gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*Lambdax*ZH(gI1,0)*ZH(gI2,1),0) + 0.15*
      Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZH
      (gI1,0)*ZH(gI2,0) - 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZH(gI1,0)*ZH(gI2,0) - 0.15*Sqr(gN)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZH(gI1,0)*ZH(gI2,0) -
      0.15*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*ZH(gI1,1)*ZH(gI2,1) + 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*ZH(gI1,1)*ZH(gI2,1) - 0.1*Sqr(gN)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*ZH(gI1,1)*ZH(gI2,1) +
      0.25*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*ZH(gI1,2)*ZH(gI2,2) - SUM(j2,0,1,AbsSqr(Lambda12(j2,j2))*KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2))*ZH(gI1,2)*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpSvconjUSHIpconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.15*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 2,0.25*KroneckerDelta(gI1,gI2)
      *KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 2,0.15*KroneckerDelta(gI1,gI2
      )*KroneckerDelta(gO1,gO2)*Sqr(gN),0) + 0.15*KroneckerDelta(gI1,gI2)*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1)) - 0.25*
      KroneckerDelta(gI1,gI2)*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1)) + 0.1*KroneckerDelta(gI1,gI2)*Sqr(gN)*SUM(j1,0,1
      ,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USHIpconjSHI0conjUSHIp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,IF(gO2 < 2,-0.5*Conj(UHI0(gI1,
      gO2))*Sqr(g2)*UHI0(gI2,gO1),0),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)),
      0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHI0
      (gI1,j1))*UHI0(gI2,j1)),0) + IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 2,-0.125
      *KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,
      2 + j1)),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),
      0) + IF(gO1 < 2,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHI0(
      gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)),0) + IF(gO1 < 2,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 + j2)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2
       + j2)),0) + IF(gO1 < 2,-0.25*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*UHI0(gI2,gO1),0) + IF(gO1 < 2,-0.25*Sqr(g2)*SUM(
      j2,0,1,Conj(UHI0(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*UHI0(gI2,gO1),0) +
      IF(gO1 < 2,KroneckerDelta(2 + gO1,gO2)*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*
      Conj(Lambda12(j2,j2))*UHI0(gI2,j2))*Lambda12(gO1,gO1),0) + IF(gO2 < 2,-0.25*
      Conj(UHI0(gI1,gO2))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHI0(gI2,2
       + j1)),0) + IF(gO2 < 2,Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*
      SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,2 + j1)*Lambda12(j1,j1)),0) + IF(gO2
      < 2,-0.25*Conj(UHI0(gI1,gO2))*Sqr(g2)*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*
      UHI0(gI2,2 + j2)),0) - 0.25*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      UHI0(gI2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*KroneckerDelta(gO2,2 +
      j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.125*Sqr(g2)*SUM(
      j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2
      )*KroneckerDelta(gO2,2 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*
      UHI0(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM
      (j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.125*Sqr(
      g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1
      ,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,
      2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(
      gI1,j2))*UHI0(gI2,j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) -
      0.075*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,1,Conj(UHI0(gI1,j2))*UHI0(gI2,j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(
      gI1,2 + j2))*UHI0(gI2,2 + j2)) + 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1
      ,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(
      gI2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHI0(gI1,2 + j2))*UHI0(gI2,2 +
      j2)) - 0.25*Sqr(g2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHI0(gI2,2 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpUSHIpconjSHIpconjUSHIp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,IF(gO2 < 2,-0.15*Conj(UHIp(gI1,
      gO2))*Sqr(g1)*UHIp(gI2,gO1),0),0) + IF(gO1 < 2,IF(gO2 < 2,-0.25*Conj(UHIp(
      gI1,gO2))*Sqr(g2)*UHIp(gI2,gO1),0),0) + IF(gO1 < 2,IF(gO2 < 2,-0.225*Conj(
      UHIp(gI1,gO2))*Sqr(gN)*UHIp(gI2,gO1),0),0) + IF(gO1 < 2,IF(gO2 < 2,-(Conj(
      UHIp(gI1,2 + gO2))*Conj(Lambda12(gO2,gO2))*UHIp(gI2,2 + gO1)*Lambda12(gO1,
      gO1)),0),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,1,
      Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,
      gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)),0) + IF(gO1 < 2,-
      0.1125*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(
      gI2,j1)),0) + IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,1,
      Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2
       + j1)),0) + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,1,
      Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),
      0) + IF(gO1 < 2,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(
      gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)),0) + IF(gO1 < 2,0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + IF(gO1 < 2,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,1,
      Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 + j2)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2
       + j2)),0) + IF(gO1 < 2,0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*UHIp(gI2,gO1),0) + IF(gO1 < 2,0.125*Sqr(g2)*SUM(
      j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*UHIp(gI2,gO1),0) +
      IF(gO1 < 2,-0.075*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(
      gO2,2 + j1))*UHIp(gI2,gO1),0) + IF(gO1 < 2,0.075*Sqr(g1)*SUM(j2,0,1,Conj(
      UHIp(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*UHIp(gI2,gO1),0) + IF(gO1 < 2,
      0.125*Sqr(g2)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*
      UHIp(gI2,gO1),0) + IF(gO1 < 2,-0.075*Sqr(gN)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2
      ))*KroneckerDelta(gO2,2 + j2))*UHIp(gI2,gO1),0) + IF(gO1 < 2,-(
      KroneckerDelta(2 + gO1,gO2)*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*Conj(Lambda12(
      j2,j2))*UHIp(gI2,j2))*Lambda12(gO1,gO1)),0) + IF(gO2 < 2,0.075*Conj(UHIp(gI1
      ,gO2))*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHIp(gI2,2 + j1)),0) +
      IF(gO2 < 2,0.125*Conj(UHIp(gI1,gO2))*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2
       + j1)*UHIp(gI2,2 + j1)),0) + IF(gO2 < 2,-0.075*Conj(UHIp(gI1,gO2))*Sqr(gN)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHIp(gI2,2 + j1)),0) + IF(gO2 < 2,-(
      Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*SUM(j1,0,1,Conj(UHIp(gI1
      ,j1))*UHIp(gI2,2 + j1)*Lambda12(j1,j1))),0) + IF(gO2 < 2,0.075*Conj(UHIp(gI1
      ,gO2))*Sqr(g1)*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHIp(gI2,2 + j2)),0) +
      IF(gO2 < 2,0.125*Conj(UHIp(gI1,gO2))*Sqr(g2)*SUM(j2,0,1,KroneckerDelta(gO1,2
       + j2)*UHIp(gI2,2 + j2)),0) + IF(gO2 < 2,-0.075*Conj(UHIp(gI1,gO2))*Sqr(gN)*
      SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHIp(gI2,2 + j2)),0) - 0.075*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*UHIp(gI2,2 + j1))*SUM(j2,0,1,Conj(UHIp
      (gI1,2 + j2))*KroneckerDelta(gO2,2 + j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*UHIp(gI2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2
      ))*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,
      2 + j1)*UHIp(gI2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*KroneckerDelta(
      gO2,2 + j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM
      (j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.125*Sqr(
      g2)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,1,Conj(UHIp
      (gI1,j1))*UHIp(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta
      (gO2,2 + j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 +
      j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) -
      0.125*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,1
      ,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(
      j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(
      gI1,j2))*UHIp(gI2,j2)) + 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) -
      0.075*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,1,Conj(UHIp(gI1,j2))*UHIp(gI2,j2)) - SUM(j1,0,1,Conj(UHIp(gI1,
      j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*SUM(j2,0,1,Conj(Lambda12(j2
      ,j2))*KroneckerDelta(gO1,2 + j2)*UHIp(gI2,j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(
      gI1,2 + j2))*UHIp(gI2,2 + j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1
      ,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(
      gI2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,Conj(UHIp(gI1,2 + j2))*UHIp(gI2,2 +
      j2)) - 0.075*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2
      + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*UHIp(gI2,2 + j2)) - 0.125*Sqr(
      g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1
      ,KroneckerDelta(gO1,2 + j2)*UHIp(gI2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,1,Conj
      (UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,1,KroneckerDelta(gO1
      ,2 + j2)*UHIp(gI2,2 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpHpmSHI0conjUSHIp(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 4,-0.35355339059327373*vd*Conj(
      UHI0(gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 4,-0.35355339059327373*vu*
      Conj(UHI0(gI1,gO2))*Sqr(g2)*ZP(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiIChaconjUSHIpPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,1,KroneckerDelta(gO2,2 + j1)*
      ZNI(gI1,2 + j1))*UP(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChaconjUSHIpPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-(g2*Conj(UM(gI2,0))*Conj(ZNI(
      gI1,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhSHIpconjUSHIp(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,std::complex<double>(0,0.5)*vu*
      Conj(UHIp(gI1,2 + gO2))*Conj(Lambda12(gO2,gO2))*Lambdax*ZA(gI2,0),0) + IF(
      gO2 < 2,std::complex<double>(0,0.5)*vd*Conj(UHIp(gI1,2 + gO2))*Conj(Lambda12
      (gO2,gO2))*Lambdax*ZA(gI2,1),0) + IF(gO2 < 2,std::complex<double>(0.,
      0.7071067811865475)*Conj(UHIp(gI1,2 + gO2))*Conj(TLambda12(gO2,gO2))*ZA(gI2,
      2),0) - std::complex<double>(0,0.5)*vu*Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(
      gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZA(gI2,0) - std::
      complex<double>(0,0.5)*vd*Conj(Lambdax)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*
      KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZA(gI2,1) - std::complex<double>
      (0.,0.7071067811865475)*SUM(j1,0,1,Conj(UHIp(gI1,j1))*KroneckerDelta(gO2,2 +
      j1)*TLambda12(j1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSHIpconjUSHIp(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.5*vu*Conj(UHIp(gI1,2 + gO2))*
      Conj(Lambda12(gO2,gO2))*Lambdax*ZH(gI2,0),0) + IF(gO2 < 2,-0.15*vd*Conj(UHIp
      (gI1,gO2))*Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 2,0.25*vd*Conj(UHIp(gI1,gO2))*Sqr
      (g2)*ZH(gI2,0),0) + IF(gO2 < 2,-0.225*vd*Conj(UHIp(gI1,gO2))*Sqr(gN)*ZH(gI2,
      0),0) + IF(gO2 < 2,0.5*vd*Conj(UHIp(gI1,2 + gO2))*Conj(Lambda12(gO2,gO2))*
      Lambdax*ZH(gI2,1),0) + IF(gO2 < 2,0.15*vu*Conj(UHIp(gI1,gO2))*Sqr(g1)*ZH(gI2
      ,1),0) + IF(gO2 < 2,-0.25*vu*Conj(UHIp(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF(
      gO2 < 2,-0.15*vu*Conj(UHIp(gI1,gO2))*Sqr(gN)*ZH(gI2,1),0) + IF(gO2 < 2,-(vs*
      AbsSqr(Lambda12(gO2,gO2))*Conj(UHIp(gI1,gO2))*ZH(gI2,2)),0) + IF(gO2 < 2,-
      0.7071067811865475*Conj(UHIp(gI1,2 + gO2))*Conj(TLambda12(gO2,gO2))*ZH(gI2,2
      ),0) + IF(gO2 < 2,0.375*vs*Conj(UHIp(gI1,gO2))*Sqr(gN)*ZH(gI2,2),0) + 0.15*
      vd*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*ZH(
      gI2,0) - 0.25*vd*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(
      gO2,2 + j1))*ZH(gI2,0) - 0.15*vd*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*ZH(gI2,0) + 0.5*vu*Conj(Lambdax)*SUM(j1,0,1,Conj
      (UHIp(gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZH(gI2,0) - 0.15*
      vu*Sqr(g1)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*ZH(
      gI2,1) + 0.25*vu*Sqr(g2)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(
      gO2,2 + j1))*ZH(gI2,1) - 0.1*vu*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1))*ZH(gI2,1) + 0.5*vd*Conj(Lambdax)*SUM(j1,0,1,Conj
      (UHIp(gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZH(gI2,1) + 0.25*
      vs*Sqr(gN)*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1))*ZH(
      gI2,2) - 0.7071067811865475*SUM(j1,0,1,Conj(UHIp(gI1,j1))*KroneckerDelta(gO2
      ,2 + j1)*TLambda12(j1,j1))*ZH(gI2,2) - vs*SUM(j2,0,1,AbsSqr(Lambda12(j2,j2))
      *Conj(UHIp(gI1,2 + j2))*KroneckerDelta(gO2,2 + j2))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpSdUSHIpconjSdconjUSHIp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 2,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) +
      0.025*(40*IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(
      ZD(gI1,j1))*ZD(gI2,j1)),0) + 40*IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)),0) + 40*IF(gO1 < 2,0.075
      *KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 +
      j1)),0) + 40*IF(gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + 40*IF(gO1 < 2,-0.125*KroneckerDelta(gO1,
      gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + 40*IF(gO1 < 2,
      0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2
      )),0) + 40*IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZD(gI1,3 + j2))*ZD(gI2,3 + j2)),0) + 40*IF(gO1 < 2,0.075*KroneckerDelta(gO1,
      gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2)),0) - Sqr(g1)*
      SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2
      )*KroneckerDelta(gO2,2 + j2)) + 5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2
      ,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) +
      Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3
       + j1))*ZD(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta
      (gO2,2 + j2)) + 2*Sqr(gN)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))*
      SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2
      ,Conj(ZD(gI1,j2))*ZD(gI2,j2)) + 5*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 +
      j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) +
      Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) - 2*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZD(
      gI1,3 + j2))*ZD(gI2,3 + j2)) + 2*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 +
      j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 +
      j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSHIpconjSDXconjUSHIp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)),0) + IF(gO1 < 2,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)),0)
      + IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(gI1,3
      + j1))*ZDX(gI2,3 + j1)),0) + IF(gO1 < 2,-0.1125*KroneckerDelta(gO1,gO2)*Sqr(
      gN)*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)),0) + IF(gO1 < 2,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)),0)
      + IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI1,
      j2))*ZDX(gI2,j2)),0) + IF(gO1 < 2,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 + j2)),0) + IF(gO1 < 2,-0.1125*
      KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2,3 +
      j2)),0) + IF(gO1 < 2,-(KroneckerDelta(2 + gO1,gO2)*SUM(j2,0,2,Conj(ZDX(gI1,3
       + j2))*Conj(Kappa(j2,j2))*ZDX(gI2,j2))*Lambda12(gO1,gO1)),0) + IF(gO2 < 2,-
      (Conj(Lambda12(gO2,gO2))*KroneckerDelta(gO1,2 + gO2)*SUM(j1,0,2,Conj(ZDX(gI1
      ,j1))*ZDX(gI2,3 + j1)*Kappa(j1,j1))),0) + 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZDX(
      gI1,j1))*ZDX(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(
      gO2,2 + j2)) - 0.05*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1))*SUM(j2
      ,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.05*Sqr(g1)*
      SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,2,Conj(ZDX(
      gI1,3 + j1))*ZDX(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*
      KroneckerDelta(gO2,2 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,1,KroneckerDelta(gO1,2 +
      j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) -
      0.05*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1
      ))*SUM(j2,0,2,Conj(ZDX(gI1,j2))*ZDX(gI2,j2)) - 0.05*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZDX(
      gI1,3 + j2))*ZDX(gI2,3 + j2)) - 0.075*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,
      2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZDX(gI1,3 + j2))*ZDX(gI2
      ,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSHIpconjSeconjUSHIp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 2,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      0.025*(40*IF(gO1 < 2,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(
      ZE(gI1,j1))*ZE(gI2,j1)),0) + 40*IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + 40*IF(gO1 < 2,
      0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(
      gI2,3 + j1)),0) + 40*IF(gO1 < 2,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + 40*IF(gO1 < 2,-0.125*KroneckerDelta
      (gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + 40*IF(gO1 < 2
      ,0.075*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2
      )),0) + 40*IF(gO1 < 2,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + 40*IF(gO1 < 2,0.0375*KroneckerDelta(gO1
      ,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + 3*Sqr(g1)
      *SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 +
      j2)*KroneckerDelta(gO2,2 + j2)) + 5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(
      gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) +
      2*Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,1,KroneckerDelta(
      gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3
       + j1))*ZE(gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta
      (gO2,2 + j2)) + Sqr(gN)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(
      j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 3*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2
      ,Conj(ZE(gI1,j2))*ZE(gI2,j2)) + 5*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,2 +
      j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) + 2*
      Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 6*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)) + Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)
      *KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2))
      );

   return result;
}

std::complex<double> CLASSNAME::CpUSHIpSuconjUSHIpconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 2,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1)),0) + IF(gO1 < 2,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,
      2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 2,0.0375*KroneckerDelta
      (gO1,gO2)*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(
      gO1 < 2,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU
      (gI2,j2)),0) + IF(gO1 < 2,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*
      Sqr(gN)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 2,-0.1*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 +
      j2)),0) + IF(gO1 < 2,0.0375*KroneckerDelta(gO1,gO2)*Sqr(gN)*SUM(j2,0,2,Conj(
      ZU(gI1,3 + j2))*ZU(gI2,3 + j2)),0) - 0.025*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1
      ))*ZU(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) - 0.125*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,1,
      KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) + 0.025*Sqr(gN)*SUM(
      j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*
      KroneckerDelta(gO2,2 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU
      (gI2,3 + j1))*SUM(j2,0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 +
      j2)) + 0.025*Sqr(gN)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,
      0,1,KroneckerDelta(gO1,2 + j2)*KroneckerDelta(gO2,2 + j2)) - 0.025*Sqr(g1)*
      SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2
      ,Conj(ZU(gI1,j2))*ZU(gI2,j2)) - 0.125*Sqr(g2)*SUM(j1,0,1,KroneckerDelta(gO1,
      2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2))
      + 0.025*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 +
      j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)) + 0.1*Sqr(g1)*SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZU(
      gI1,3 + j2))*ZU(gI2,3 + j2)) + 0.025*Sqr(gN)*SUM(j1,0,1,KroneckerDelta(gO1,2
       + j1)*KroneckerDelta(gO2,2 + j1))*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3
      + j2));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaIconjUSHIpPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-(Conj(Lambda12(gO2,gO2))*ZN(gI1
      ,4)*ZPI(gI2,gO2)),0) - 0.1*SUM(j1,0,1,KroneckerDelta(gO2,2 + j1)*ZPI(gI2,j1)
      )*(5.477225575051661*g1*ZN(gI1,0) + 7.0710678118654755*g2*ZN(gI1,1) -
      4.47213595499958*gN*ZN(gI1,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaIconjUSHIpPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.5477225575051661*g1*Conj(ZMI(
      gI2,gO1))*Conj(ZN(gI1,0)),0) + IF(gO1 < 2,0.7071067811865475*g2*Conj(ZMI(gI2
      ,gO1))*Conj(ZN(gI1,1)),0) + IF(gO1 < 2,0.6708203932499369*gN*Conj(ZMI(gI2,
      gO1))*Conj(ZN(gI1,5)),0) - Conj(ZN(gI1,4))*SUM(j1,0,1,Conj(ZMI(gI2,j1))*
      KroneckerDelta(gO1,2 + j1)*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjUSHIpVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.7071067811865475*g2*Conj(UHI0(
      gI2,gO2)),0) - 0.7071067811865475*g2*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*
      KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjUSHIpVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 4,-0.3872983346207417*g1*Conj(UHIp
      (gI2,gO2))*Cos(ThetaW()),0) + IF(gI2 < 4,-0.5*g2*Conj(UHIp(gI2,gO2))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjUSHIpVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.5*g2*Conj(UHIp(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 2,0.3872983346207417*g1*Conj(UHIp(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 2,-0.4743416490252569*gN*
      Conj(UHIp(gI2,gO2))*Sin(ThetaWp()),0) - 0.5*g2*Cos(ThetaW())*Cos(ThetaWp())*
      SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1)) +
      0.3872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*SUM(j1,0,1,Conj(UHIp(gI2,
      2 + j1))*KroneckerDelta(gO2,2 + j1)) + 0.31622776601683794*gN*Sin(ThetaWp())
      *SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjUSHIpVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.4743416490252569*gN*Conj(UHIp
      (gI2,gO2))*Cos(ThetaWp()),0) + IF(gO2 < 2,0.5*g2*Conj(UHIp(gI2,gO2))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO2 < 2,-0.3872983346207417*g1*Conj(UHIp(
      gI2,gO2))*Sin(ThetaW())*Sin(ThetaWp()),0) + 0.31622776601683794*gN*Cos(
      ThetaWp())*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1)) +
      0.5*g2*Cos(ThetaW())*Sin(ThetaWp())*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*
      KroneckerDelta(gO2,2 + j1)) - 0.3872983346207417*g1*Sin(ThetaW())*Sin(
      ThetaWp())*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*KroneckerDelta(gO2,2 + j1));

   return result;
}

double CLASSNAME::CpUSSI0conjUSSI0VZVZ(int gO1, int gO2) const
{
   
   const double result = 1.25*KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpUSSI0conjUSSI0VZpVZp(int gO1, int gO2) const
{
   
   const double result = 1.25*KroneckerDelta(gO1,gO2)*Sqr(gN)*Sqr(Cos(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSSI0conjHpmconjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(3*ZP
      (gI1,0)*ZP(gI2,0) + 2*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0USSI0conjSHp0conjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.25*KroneckerDelta(gO1,gO2)*Sqr(gN)*(-(
      Conj(UHp0(gI1,0))*UHp0(gI2,0)) + Conj(UHp0(gI1,1))*UHp0(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppUSSI0conjSHppconjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.25*KroneckerDelta(gO1,gO2)*Sqr(gN)*(-(
      Conj(UHpp(gI1,0))*UHpp(gI2,0)) + Conj(UHpp(gI1,1))*UHpp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSSI0USSI0conjSSI0conjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 2,IF(gI2 < 2,-0.625*Conj(ZSSI(gI1,
      gO2))*Sqr(gN)*ZSSI(gI2,gO1),0),0) - 0.625*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CphhSSI0conjUSSI0(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 2,0.375*vd*Conj(ZSSI(gI1,gO2))*Sqr
      (gN)*ZH(gI2,0),0) + IF(gI1 < 2,0.25*vu*Conj(ZSSI(gI1,gO2))*Sqr(gN)*ZH(gI2,1)
      ,0) + IF(gI1 < 2,-0.625*vs*Conj(ZSSI(gI1,gO2))*Sqr(gN)*ZH(gI2,2),0);

   return result;
}

double CLASSNAME::CpChiFSIconjUSSI0PR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpChiFSIconjUSSI0PL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-1.118033988749895*gN*Conj(ZFSI(
      gI1,gO1))*Conj(ZN(gI2,5)),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSSI0conjUSSI0(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(3*ZA
      (gI1,0)*ZA(gI2,0) + 2*ZA(gI1,1)*ZA(gI2,1) - 5*ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSSI0conjUSSI0(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(3*ZH
      (gI1,0)*ZH(gI2,0) + 2*ZH(gI1,1)*ZH(gI2,1) - 5*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpUSSI0SvconjUSSI0conjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const double result = -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr
      (gN);

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USSI0conjSHI0conjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(3*
      SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) + 2*SUM(j1,0,1,Conj(UHI0(gI1,2 +
      j1))*UHI0(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpUSSI0conjSHIpconjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(3*
      SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) + 2*SUM(j1,0,1,Conj(UHIp(gI1,2 +
      j1))*UHIp(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSSI0conjSdconjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(SUM
      (j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(
      gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSSI0conjSDXconjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(2*
      SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) + 3*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1
      ))*ZDX(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSSI0conjSeconjUSSI0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.125*KroneckerDelta(gO1,gO2)*Sqr(gN)*(2*
      SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE
      (gI2,3 + j1)));

   return result;
}

double CLASSNAME::CpUSSI0SuconjUSSI0conjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const double result = -0.125*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*
      Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CpSSI0conjUSSI0VZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.7905694150420949*gN*Conj(ZSSI(
      gI2,gO2))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSSI0conjUSSI0VZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.7905694150420949*gN*Conj(ZSSI(
      gI2,gO2))*Cos(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0conjUSHp0VZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.0125*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(-
      9.797958971132712*g1*gN*Cos(ThetaW() - 2*ThetaWp()) + 9.797958971132712*g1*
      gN*Cos(ThetaW() + 2*ThetaWp()) - 15.491933384829668*g1*g2*Sin(2*ThetaW()) +
      12.649110640673518*g2*gN*Sin(ThetaW() - 2*ThetaWp()) - 7.745966692414834*g1*
      g2*Sin(2*(ThetaW() - ThetaWp())) - 7.745966692414834*g1*g2*Sin(2*(ThetaW() +
      ThetaWp())) - 12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) - 6*Sqr(
      g1) + 3*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) +
      3*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g1) + 2*Cos(2*ThetaW())*(3*Sqr(g1) - 5*
      Sqr(g2)) - 10*Sqr(g2) - 5*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g2) - 10*Cos(2*
      ThetaWp())*Sqr(g2) - 5*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g2) - 8*Sqr(gN) + 8
      *Cos(2*ThetaWp())*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0conjUSHp0VZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.0125*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(-
      9.797958971132712*g1*gN*Cos(ThetaW() - 2*ThetaWp()) + 9.797958971132712*g1*
      gN*Cos(ThetaW() + 2*ThetaWp()) + 15.491933384829668*g1*g2*Sin(2*ThetaW()) +
      12.649110640673518*g2*gN*Sin(ThetaW() - 2*ThetaWp()) - 7.745966692414834*g1*
      g2*Sin(2*(ThetaW() - ThetaWp())) - 7.745966692414834*g1*g2*Sin(2*(ThetaW() +
      ThetaWp())) - 12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) + 6*Sqr(
      g1) + 3*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) +
      3*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g1) + 10*Sqr(g2) - 5*Cos(2*(ThetaW() -
      ThetaWp()))*Sqr(g2) - 10*Cos(2*ThetaWp())*Sqr(g2) - 5*Cos(2*(ThetaW() +
      ThetaWp()))*Sqr(g2) + Cos(2*ThetaW())*(-6*Sqr(g1) + 10*Sqr(g2)) + 8*Sqr(gN)
      + 8*Cos(2*ThetaWp())*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0conjUSHp0conjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSHp0conjHpmconjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) - 5*Sqr(
      g2) - 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*
      ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0USHp0conjSHp0conjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*(
      -(Conj(UHp0(gI1,1))*(-2*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*UHp0(gI2
      ,1) + KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*UHp0(gI2,0) +
      KroneckerDelta(0,gO2)*UHp0(gI2,1)))) + Conj(UHp0(gI1,0))*(2*KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*UHp0(gI2,0) - KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*UHp0(gI2,0) + KroneckerDelta(0,gO2)*UHp0(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0SHppconjUSHp0conjSHpp(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-(Conj(UHpp(gI1,0))*(KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*UHpp(gI2,0)
      + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) - 2*
      Sqr(gN))*UHpp(gI2,0) + 10*KroneckerDelta(0,gO2)*Sqr(g2)*UHpp(gI2,1)))) -
      Conj(UHpp(gI1,1))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) +
      5*Sqr(g2) + 2*Sqr(gN))*UHpp(gI2,1) + KroneckerDelta(0,gO1)*(10*
      KroneckerDelta(1,gO2)*Sqr(g2)*UHpp(gI2,0) + KroneckerDelta(0,gO2)*(-3*Sqr(g1
      ) + 5*Sqr(g2) - 2*Sqr(gN))*UHpp(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0SSI0conjUSHp0conjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.25*(-(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjHpmconjUSHp0(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.35355339059327373*(Conj(UHpp(gI2,0))*
      KroneckerDelta(0,gO2) + Conj(UHpp(gI2,1))*KroneckerDelta(1,gO2))*Sqr(g2)*(vd
      *ZP(gI1,0) + vu*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhSHp0conjUSHp0(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.05*(Conj(UHp0(gI1,0))*KroneckerDelta(0,
      gO2) - Conj(UHp0(gI1,1))*KroneckerDelta(1,gO2))*(vd*(3*Sqr(g1) + 5*Sqr(g2) -
      3*Sqr(gN))*ZH(gI2,0) - vu*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*ZH(gI2,1) + 5*
      vs*Sqr(gN)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiPconjUSHp0PR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(1,gO2)*(-
      5.477225575051661*g1*ZN(gI2,0) + 7.0710678118654755*g2*ZN(gI2,1) +
      4.47213595499958*gN*ZN(gI2,5))*ZNp(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpChiChiPconjUSHp0PL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = 0.1*(5.477225575051661*g1*Conj(ZN(gI2,0)) -
      7.0710678118654755*g2*Conj(ZN(gI2,1)) - 4.47213595499958*gN*Conj(ZN(gI2,5)))
      *Conj(ZNp(gI1,0))*KroneckerDelta(0,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaPconjUSHp0PR(int gI1, int gO2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(1,gO2)*UM(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaPconjUSHp0PL(int gI1, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(0,gO1))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSHp0conjUSHp0(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) + 5*Sqr(
      g2) - 3*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*
      ZA(gI1,1)*ZA(gI2,1) + 5*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSHp0conjUSHp0(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) + 5*Sqr(
      g2) - 3*Sqr(gN))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*
      ZH(gI1,1)*ZH(gI2,1) + 5*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0SvconjUSHp0conjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,
      gI2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USHp0conjSHI0conjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) + 5*Sqr(
      g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) - (3*Sqr(g1) +
      5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpUSHp0conjSHIpconjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) - 5*Sqr(
      g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) + (-3*Sqr(g1) +
      5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSHp0conjSdconjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((Sqr(g1) + 5*Sqr(g2)
      - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) - Sqr(gN))*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSHp0conjSDXconjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(2*(Sqr(g1) - Sqr(gN)
      )*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) - (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1
      ,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSHp0conjSeconjUSHp0(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) - 5*Sqr(
      g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + Sqr
      (gN))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUSHp0SuconjUSHp0conjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((Sqr(g1) - 5*Sqr(g2)
      - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - (4*Sqr(g1) + Sqr(gN))*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0conjUSHp0VZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.5*g2*Conj(UHp0(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 2,0.3872983346207417*g1*Conj(UHp0(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 2,0.31622776601683794*gN*
      Conj(UHp0(gI2,gO2))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHp0conjUSHp0VZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.31622776601683794*gN*Conj(UHp0
      (gI2,gO2))*Cos(ThetaWp()),0) + IF(gI2 < 2,-0.5*g2*Conj(UHp0(gI2,gO2))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 2,-0.3872983346207417*g1*Conj(UHp0(
      gI2,gO2))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjUSHp0conjVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*(Conj(UHpp(gI2,0))*
      KroneckerDelta(0,gO2) - Conj(UHpp(gI2,1))*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpUSHppconjUSHppVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.0125*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(-
      9.797958971132712*g1*gN*Cos(ThetaW() - 2*ThetaWp()) + 9.797958971132712*g1*
      gN*Cos(ThetaW() + 2*ThetaWp()) + 15.491933384829668*g1*g2*Sin(2*ThetaW()) -
      12.649110640673518*g2*gN*Sin(ThetaW() - 2*ThetaWp()) + 7.745966692414834*g1*
      g2*Sin(2*(ThetaW() - ThetaWp())) + 7.745966692414834*g1*g2*Sin(2*(ThetaW() +
      ThetaWp())) + 12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) - 6*Sqr(
      g1) + 3*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) +
      3*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g1) + 2*Cos(2*ThetaW())*(3*Sqr(g1) - 5*
      Sqr(g2)) - 10*Sqr(g2) - 5*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g2) - 10*Cos(2*
      ThetaWp())*Sqr(g2) - 5*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g2) - 8*Sqr(gN) + 8
      *Cos(2*ThetaWp())*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CpUSHppconjUSHppVZpVZp(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.0125*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(-
      9.797958971132712*g1*gN*Cos(ThetaW() - 2*ThetaWp()) + 9.797958971132712*g1*
      gN*Cos(ThetaW() + 2*ThetaWp()) - 15.491933384829668*g1*g2*Sin(2*ThetaW()) -
      12.649110640673518*g2*gN*Sin(ThetaW() - 2*ThetaWp()) + 7.745966692414834*g1*
      g2*Sin(2*(ThetaW() - ThetaWp())) + 7.745966692414834*g1*g2*Sin(2*(ThetaW() +
      ThetaWp())) + 12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) + 6*Sqr(
      g1) + 3*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) +
      3*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g1) + 10*Sqr(g2) - 5*Cos(2*(ThetaW() -
      ThetaWp()))*Sqr(g2) - 10*Cos(2*ThetaWp())*Sqr(g2) - 5*Cos(2*(ThetaW() +
      ThetaWp()))*Sqr(g2) + Cos(2*ThetaW())*(-6*Sqr(g1) + 10*Sqr(g2)) + 8*Sqr(gN)
      + 8*Cos(2*ThetaWp())*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CpUSHppconjUSHppconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSHppconjHpmconjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) + 5*Sqr(
      g2) - 3*Sqr(gN))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*
      ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0USHppconjSHp0conjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-(Conj(UHp0(gI1,0))*(KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*UHp0(gI2,0)
      + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) - 2*
      Sqr(gN))*UHp0(gI2,0) + 10*KroneckerDelta(0,gO2)*Sqr(g2)*UHp0(gI2,1)))) -
      Conj(UHp0(gI1,1))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) +
      5*Sqr(g2) + 2*Sqr(gN))*UHp0(gI2,1) + KroneckerDelta(0,gO1)*(10*
      KroneckerDelta(1,gO2)*Sqr(g2)*UHp0(gI2,0) + KroneckerDelta(0,gO2)*(-3*Sqr(g1
      ) + 5*Sqr(g2) - 2*Sqr(gN))*UHp0(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpSHppUSHppconjSHppconjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(3*Sqr(g1) + 5*Sqr(g2) + 2*Sqr(gN))*(
      -(Conj(UHpp(gI1,1))*(-2*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*UHpp(gI2
      ,1) + KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*UHpp(gI2,0) +
      KroneckerDelta(0,gO2)*UHpp(gI2,1)))) + Conj(UHpp(gI1,0))*(2*KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*UHpp(gI2,0) - KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*UHpp(gI2,0) + KroneckerDelta(0,gO2)*UHpp(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUSHppSSI0conjUSHppconjSSI0(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.25*(-(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*Sqr(gN);

   return result;
}

std::complex<double> CLASSNAME::CpHpmSHp0conjUSHpp(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.35355339059327373*(Conj(UHp0(gI1,0))*
      KroneckerDelta(0,gO2) + Conj(UHp0(gI1,1))*KroneckerDelta(1,gO2))*Sqr(g2)*(vd
      *ZP(gI2,0) + vu*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChaconjUSHppPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)*ZNp(
      gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChaconjUSHppPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Conj(ZNp(gI1,0))*
      KroneckerDelta(0,gO1));

   return result;
}

std::complex<double> CLASSNAME::CphhSHppconjUSHpp(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.05*(Conj(UHpp(gI1,0))*KroneckerDelta(0,
      gO2) - Conj(UHpp(gI1,1))*KroneckerDelta(1,gO2))*(vd*(3*Sqr(g1) - 5*Sqr(g2) -
      3*Sqr(gN))*ZH(gI2,0) + vu*(-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*ZH(gI2,1) + 5
      *vs*Sqr(gN)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSHppconjUSHpp(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) - 5*Sqr(
      g2) - 3*Sqr(gN))*ZA(gI1,0)*ZA(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*
      ZA(gI1,1)*ZA(gI2,1) + 5*Sqr(gN)*ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSHppconjUSHpp(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) - 5*Sqr(
      g2) - 3*Sqr(gN))*ZH(gI1,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) - 2*Sqr(gN))*
      ZH(gI1,1)*ZH(gI2,1) + 5*Sqr(gN)*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUSHppSvconjUSHppconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,
      gI2)*(3*Sqr(g1) - 5*Sqr(g2) + 2*Sqr(gN));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0USHppconjSHI0conjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) - 5*Sqr(
      g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) + (-3*Sqr(g1) +
      5*Sqr(g2) - 2*Sqr(gN))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpUSHppconjSHIpconjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) + 5*Sqr(
      g2) - 3*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) - (3*Sqr(g1) +
      5*Sqr(g2) + 2*Sqr(gN))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSHppconjSdconjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((Sqr(g1) - 5*Sqr(g2)
      - Sqr(gN))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) - Sqr(gN))*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXUSHppconjSDXconjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(2*(Sqr(g1) - Sqr(gN)
      )*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) - (2*Sqr(g1) + 3*Sqr(gN))*SUM(j1
      ,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSHppconjSeconjUSHpp(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr(g1) + 5*Sqr(
      g2) + 2*Sqr(gN))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + Sqr
      (gN))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUSHppSuconjUSHppconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((Sqr(g1) + 5*Sqr(g2)
      - Sqr(gN))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - (4*Sqr(g1) + Sqr(gN))*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaPconjUSHppPR(int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.1*KroneckerDelta(1,gO2)*(
      5.477225575051661*g1*ZN(gI1,0) + 7.0710678118654755*g2*ZN(gI1,1) -
      4.47213595499958*gN*ZN(gI1,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaPconjUSHppPL(int gI1, int gO1) const
{
   
   const std::complex<double> result = 0.1*(5.477225575051661*g1*Conj(ZN(gI1,0)) +
      7.0710678118654755*g2*Conj(ZN(gI1,1)) - 4.47213595499958*gN*Conj(ZN(gI1,5)))
      *KroneckerDelta(0,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpSHp0conjUSHppVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*(Conj(UHp0(gI2,0))*
      KroneckerDelta(0,gO2) - Conj(UHp0(gI2,1))*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjUSHppVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.3872983346207417*g1*Conj(UHpp
      (gI2,gO2))*Cos(ThetaW()),0) + IF(gI2 < 2,-0.5*g2*Conj(UHpp(gI2,gO2))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjUSHppVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.5*g2*Conj(UHpp(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 2,0.3872983346207417*g1*Conj(UHpp(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 2,0.31622776601683794*gN*
      Conj(UHpp(gI2,gO2))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjUSHppVZp(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.31622776601683794*gN*Conj(UHpp
      (gI2,gO2))*Cos(ThetaWp()),0) + IF(gI2 < 2,0.5*g2*Conj(UHpp(gI2,gO2))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 2,-0.3872983346207417*g1*Conj(UHpp(
      gI2,gO2))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpVGVGVG() const
{
   
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbargGgGVG() const
{
   
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpbarFdFdVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFdFdVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFDXFDXVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFDXFDXVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFuFuVGPL(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpbarFuFuVGPR(int gI1, int gI2) const
{
   
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpSdconjSdVGVG(int gI1, int gI2) const
{
   
   const double result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpSDXconjSDXVGVG(int gI1, int gI2) const
{
   
   const double result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpSuconjSuVGVG(int gI1, int gI2) const
{
   
   const double result = 6*KroneckerDelta(gI1,gI2)*Sqr(g3);

   return result;
}

double CLASSNAME::CpSdconjSdVG(int gI2, int gI1) const
{
   
   const double result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpSDXconjSDXVG(int gI2, int gI1) const
{
   
   const double result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpSuconjSuVG(int gI2, int gI1) const
{
   
   const double result = g3*KroneckerDelta(gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpGluGluVGPL() const
{
   
   const std::complex<double> result = std::complex<double>(0,1)*g3*AbsSqr(
      PhaseGlu);

   return result;
}

std::complex<double> CLASSNAME::CpGluGluVGPR() const
{
   
   const std::complex<double> result = std::complex<double>(0,1)*g3*AbsSqr(
      PhaseGlu);

   return result;
}

double CLASSNAME::CpVGVGVGVG1() const
{
   
   const double result = -16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGVGVG2() const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG3() const
{
   
   const double result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpbargWmgWmVP() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmCgWmCVP() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmVPVWm() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarChaPChaPVPPL() const
{
   
   const double result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW(
      )));

   return result;
}

double CLASSNAME::CpbarChaPChaPVPPR() const
{
   
   const double result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW(
      )));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + Cos(2*ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) + 5*Sqr(g2))*(ZP(
      gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjSHppVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + Cos(2*ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) + 5*Sqr(g2))*(Conj(
      UHpp(gI1,0))*UHpp(gI2,0) + Conj(UHpp(gI1,1))*UHpp(gI2,1));

   return result;
}

double CLASSNAME::CpHpmconjHpmVP(int gI2, int gI1) const
{
   
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpSHppconjSHppVP(int gI2, int gI1) const
{
   
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVPPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*Sin(ThetaW())*UM(gI1,0)
      + 0.5*Conj(UM(gI2,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()
      ))*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVPPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UP(gI1,0))*Sin(ThetaW())*UP(gI2,0)
      + 0.5*Conj(UP(gI1,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()
      ))*UP(gI2,1);

   return result;
}

double CLASSNAME::CpbarChaIChaIVPPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarChaIChaIVPPR(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVPPL(int gI1, int gI2) const
{
   
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVPPR(int gI1, int gI2) const
{
   
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFDXFDXVPPL(int gI1, int gI2) const
{
   
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFDXFDXVPPR(int gI1, int gI2) const
{
   
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFeFeVPPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVPPR(int gI1, int gI2) const
{
   
   const double result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpbarFuFuVPPL(int gI1, int gI2) const
{
   
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVPPR(int gI1, int gI2) const
{
   
   const double result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjSHIpVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(gI1,gI2)*(g2*Sin(ThetaW(
      ))*(7.745966692414834*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr
      (Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpSHIpconjSHIpVP(int gI2, int gI1) const
{
   
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.03333333333333333*((-7.745966692414834*g1
      *g2*Cos(ThetaW())*Sin(ThetaW()) + Sqr(g1)*Sqr(Cos(ThetaW())) + 15*Sqr(g2)*
      Sqr(Sin(ThetaW())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*Sqr(g1)*Sqr(
      Cos(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

double CLASSNAME::CpSDXconjSDXVPVP(int gI1, int gI2) const
{
   
   const double result = 0.13333333333333333*KroneckerDelta(gI1,gI2)*Sqr(g1)*Sqr(
      Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*((g2*Sin(ThetaW())*(7.745966692414834*
      g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())))*SUM(
      j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 12*Sqr(g1)*Sqr(Cos(ThetaW()))*SUM(j1,0
      ,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.03333333333333333*((g2*Sin(ThetaW())*(
      7.745966692414834*g1*Cos(ThetaW()) + 15*g2*Sin(ThetaW())) + Sqr(g1)*Sqr(Cos(
      ThetaW())))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + 16*Sqr(g1)*Sqr(Cos(
      ThetaW()))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVP(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.16666666666666666*(0.7745966692414834*g1*
      Cos(ThetaW()) - 3*g2*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))
      - 0.2581988897471611*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1
      ,3 + j1));

   return result;
}

double CLASSNAME::CpSDXconjSDXVP(int gI2, int gI1) const
{
   
   const double result = -0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,
      gI2);

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVP(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5*(0.7745966692414834*g1*Cos(ThetaW()) +
      g2*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) -
      0.7745966692414834*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3
       + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVP(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.16666666666666666*(0.7745966692414834*g1*
      Cos(ThetaW()) + 3*g2*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))
      + 0.5163977794943222*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1
      ,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjVWmVP(int gI2) const
{
   
   const std::complex<double> result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd
      *ZP(gI2,0) - vu*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpconjVWmVPVPVWm1() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVPVPVWm2() const
{
   
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVPVPVWm3() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmgWmVZ() const
{
   
   const double result = -(g2*Cos(ThetaW())*Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpbargWmCgWmCVZ() const
{
   
   const double result = g2*Cos(ThetaW())*Cos(ThetaWp());

   return result;
}

double CLASSNAME::CpconjVWmVWmVZ() const
{
   
   const double result = -(g2*Cos(ThetaW())*Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarChaPChaPVZPL() const
{
   
   const double result = 0.1*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) -
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 3.1622776601683795*gN*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarChaPChaPVZPR() const
{
   
   const double result = 0.1*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) -
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 3.1622776601683795*gN*
      Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + Cos(ThetaW())*(18.973665961010276*g2*gN*Cos(ThetaWp())*Sin
      (ThetaWp()) - 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) +
      6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(ThetaWp
      ())))*ZP(gI1,0)*ZP(gI2,0) + 2*(2.449489742783178*g1*gN*Sin(ThetaW())*Sin(2*
      ThetaWp()) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 2*Cos(ThetaW
      ())*(3.1622776601683795*g2*gN*Cos(ThetaWp())*Sin(ThetaWp()) +
      3.872983346207417*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 3*Sqr(g1)*Sqr(
      Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 2*Sqr(gN)*Sqr(Sin(ThetaWp())))*ZP(gI1,1
      )*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0conjSHp0VZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.0125*(-9.797958971132712*g1*gN*Cos(
      ThetaW() - 2*ThetaWp()) + 9.797958971132712*g1*gN*Cos(ThetaW() + 2*ThetaWp()
      ) - 15.491933384829668*g1*g2*Sin(2*ThetaW()) + 12.649110640673518*g2*gN*Sin(
      ThetaW() - 2*ThetaWp()) - 7.745966692414834*g1*g2*Sin(2*(ThetaW() - ThetaWp(
      ))) - 7.745966692414834*g1*g2*Sin(2*(ThetaW() + ThetaWp())) -
      12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) - 6*Sqr(g1) + 3*Cos(2*(
      ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) + 3*Cos(2*(
      ThetaW() + ThetaWp()))*Sqr(g1) + 2*Cos(2*ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) -
      10*Sqr(g2) - 5*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g2) - 10*Cos(2*ThetaWp())*
      Sqr(g2) - 5*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g2) - 8*Sqr(gN) + 8*Cos(2*
      ThetaWp())*Sqr(gN))*(Conj(UHp0(gI1,0))*UHp0(gI2,0) + Conj(UHp0(gI1,1))*UHp0(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjSHppVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.0125*(-9.797958971132712*g1*gN*Cos(
      ThetaW() - 2*ThetaWp()) + 9.797958971132712*g1*gN*Cos(ThetaW() + 2*ThetaWp()
      ) + 15.491933384829668*g1*g2*Sin(2*ThetaW()) - 12.649110640673518*g2*gN*Sin(
      ThetaW() - 2*ThetaWp()) + 7.745966692414834*g1*g2*Sin(2*(ThetaW() - ThetaWp(
      ))) + 7.745966692414834*g1*g2*Sin(2*(ThetaW() + ThetaWp())) +
      12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) - 6*Sqr(g1) + 3*Cos(2*(
      ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) + 3*Cos(2*(
      ThetaW() + ThetaWp()))*Sqr(g1) + 2*Cos(2*ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) -
      10*Sqr(g2) - 5*Cos(2*(ThetaW() - ThetaWp()))*Sqr(g2) - 10*Cos(2*ThetaWp())*
      Sqr(g2) - 5*Cos(2*(ThetaW() + ThetaWp()))*Sqr(g2) - 8*Sqr(gN) + 8*Cos(2*
      ThetaWp())*Sqr(gN))*(Conj(UHpp(gI1,0))*UHpp(gI2,0) + Conj(UHpp(gI1,1))*UHpp(
      gI2,1));

   return result;
}

double CLASSNAME::CpSSI0conjSSI0VZVZ(int gI1, int gI2) const
{
   
   const double result = 1.25*KroneckerDelta(gI1,gI2)*Sqr(gN)*Sqr(Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*((-10*g2*Cos(ThetaW())*Cos(ThetaWp())
      + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*
      Sin(ThetaWp()))*ZP(gI1,0)*ZP(gI2,0) + 2*(-5*g2*Cos(ThetaW())*Cos(ThetaWp())
      + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*
      Sin(ThetaWp()))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpSHp0conjSHp0VZ(int gI2, int gI1) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpSHppconjSHppVZ(int gI2, int gI1) const
{
   
   const double result = -0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpSSI0conjSSI0VZ(int gI2, int gI1) const
{
   
   const double result = 0.7905694150420949*gN*KroneckerDelta(gI1,gI2)*Sin(ThetaWp
      ());

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*Cos(
      ThetaWp())*UM(gI1,0) + 0.05*Conj(UM(gI2,1))*(10*g2*Cos(ThetaW())*Cos(ThetaWp
      ()) - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 9.486832980505138*
      gN*Sin(ThetaWp()))*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UP(gI1,0))*Cos(ThetaW())*Cos(
      ThetaWp())*UP(gI2,0) + 0.1*Conj(UP(gI1,1))*(5*g2*Cos(ThetaW())*Cos(ThetaWp()
      ) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 3.1622776601683795*
      gN*Sin(ThetaWp()))*UP(gI2,1);

   return result;
}

double CLASSNAME::CpbarChaIChaIVZPL(int gI1, int gI2) const
{
   
   const double result = 0.05*KroneckerDelta(gI1,gI2)*(10*g2*Cos(ThetaW())*Cos(
      ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      9.486832980505138*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarChaIChaIVZPR(int gI1, int gI2) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChiPVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.1*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*
      Sin(ThetaWp()))*(Conj(ZNp(gI2,0))*ZNp(gI1,0) - Conj(ZNp(gI2,1))*ZNp(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChiPVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*
      Sin(ThetaWp()))*(Conj(ZNp(gI1,0))*ZNp(gI2,0) - Conj(ZNp(gI1,1))*ZNp(gI2,1));

   return result;
}

double CLASSNAME::CpFSIFSIVZPL(int gI1, int gI2) const
{
   
   const double result = -0.7905694150420949*gN*KroneckerDelta(gI1,gI2)*Sin(
      ThetaWp());

   return result;
}

double CLASSNAME::CpFSIFSIVZPR(int gI1, int gI2) const
{
   
   const double result = 0.7905694150420949*gN*KroneckerDelta(gI1,gI2)*Sin(ThetaWp
      ());

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp())*
      Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp())))
      + 6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(
      ThetaWp())))*ZA(gI1,0)*ZA(gI2,0) + 2*(3.1622776601683795*g2*gN*Cos(ThetaW())
      *Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3
      *g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(2*ThetaWp()) +
      2*gN*Sqr(Sin(ThetaWp()))))*ZA(gI1,1)*ZA(gI2,1) + 25*Sqr(gN)*Sqr(Sin(ThetaWp(
      )))*ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp())*
      Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp())))
      + 6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(
      ThetaWp())))*ZH(gI1,0)*ZH(gI2,0) + 2*(3.1622776601683795*g2*gN*Cos(ThetaW())
      *Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3
      *g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(2*ThetaWp()) +
      2*gN*Sqr(Sin(ThetaWp()))))*ZH(gI1,1)*ZH(gI2,1) + 25*Sqr(gN)*Sqr(Sin(ThetaWp(
      )))*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(gI1,gI2)*(
      3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*
      g1*Sin(ThetaW())*Sin(2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.05)*((10*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())
      - 9.486832980505138*gN*Sin(ThetaWp()))*ZA(gI2,0)*ZH(gI1,0) - 2*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 3.1622776601683795*gN*Sin(ThetaWp()))*ZA(gI2,1)*ZH(gI1,1) +
      15.811388300841898*gN*Sin(ThetaWp())*ZA(gI2,2)*ZH(gI1,2));

   return result;
}

double CLASSNAME::CpSvconjSvVZ(int gI2, int gI1) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdFdVZPL(int gI1, int gI2) const
{
   
   const double result = 0.016666666666666666*KroneckerDelta(gI1,gI2)*(30*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())
      - 9.486832980505138*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdFdVZPR(int gI1, int gI2) const
{
   
   const double result = KroneckerDelta(gI1,gI2)*(-0.2581988897471611*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 0.31622776601683794*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFDXFDXVZPL(int gI1, int gI2) const
{
   
   const double result = KroneckerDelta(gI1,gI2)*(-0.2581988897471611*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 0.31622776601683794*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFDXFDXVZPR(int gI1, int gI2) const
{
   
   const double result = -0.016666666666666666*KroneckerDelta(gI1,gI2)*(
      15.491933384829668*g1*Cos(ThetaWp())*Sin(ThetaW()) + 28.460498941515414*gN*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPL(int gI1, int gI2) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPR(int gI1, int gI2) const
{
   
   const double result = -0.05*KroneckerDelta(gI1,gI2)*(15.491933384829668*g1*Cos(
      ThetaWp())*Sin(ThetaW()) - 3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPL(int gI1, int gI2) const
{
   
   const double result = -0.016666666666666666*KroneckerDelta(gI1,gI2)*(30*g2*Cos(
      ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 9.486832980505138*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPR(int gI1, int gI2) const
{
   
   const double result = 0.016666666666666666*KroneckerDelta(gI1,gI2)*(
      30.983866769659336*g1*Cos(ThetaWp())*Sin(ThetaW()) + 9.486832980505138*gN*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZPL(int gI1, int gI2) const
{
   
   const double result = -0.1*KroneckerDelta(gI1,gI2)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      3.1622776601683795*gN*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjSHI0VZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp())*
      Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp())))
      + 6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(
      ThetaWp())))*SUM(j1,0,1,Conj(UHI0(gI1,j1))*UHI0(gI2,j1)) + 2*(
      3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + gN*(2.449489742783178*
      g1*Sin(ThetaW())*Sin(2*ThetaWp()) + 2*gN*Sqr(Sin(ThetaWp()))))*SUM(j1,0,1,
      Conj(UHI0(gI1,2 + j1))*UHI0(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjSHIpVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-14.696938456699067*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + Cos(ThetaW())*(18.973665961010276*g2*gN*Cos(ThetaWp())*Sin
      (ThetaWp()) - 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) +
      6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(ThetaWp
      ())))*SUM(j1,0,1,Conj(UHIp(gI1,j1))*UHIp(gI2,j1)) + 2*(2.449489742783178*g1*
      gN*Sin(ThetaW())*Sin(2*ThetaWp()) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) - 2*Cos(ThetaW())*(3.1622776601683795*g2*gN*Cos(ThetaWp())*Sin(
      ThetaWp()) + 3.872983346207417*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 3*
      Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 2*Sqr(gN)*Sqr(Sin(ThetaWp()
      )))*SUM(j1,0,1,Conj(UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjSHI0VZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*Sin
      (ThetaWp()))*SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(gI1,j1)) + 2*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 3.1622776601683795*gN*Sin(ThetaWp()))*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*
      UHI0(gI1,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjSHIpVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*((-10*g2*Cos(ThetaW())*Cos(ThetaWp())
      + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*
      Sin(ThetaWp()))*SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(gI1,j1)) + 2*(-5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 3.1622776601683795*gN*Sin(ThetaWp()))*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*
      UHIp(gI1,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-10*g2*Cos(ThetaW())*Cos(ThetaWp())
      - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 9.486832980505138*gN*
      Sin(ThetaWp()))*SUM(j1,0,1,Conj(ZNI(gI2,j1))*ZNI(gI1,j1)) + 2*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 3.1622776601683795*gN*Sin(ThetaWp()))*SUM(j1,0,1,Conj(ZNI(gI2,2 + j1))*ZNI
      (gI1,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((10*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*Sin
      (ThetaWp()))*SUM(j1,0,1,Conj(ZNI(gI1,j1))*ZNI(gI2,j1)) - 2*(5*g2*Cos(ThetaW(
      ))*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      3.1622776601683795*gN*Sin(ThetaWp()))*SUM(j1,0,1,Conj(ZNI(gI1,2 + j1))*ZNI(
      gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.016666666666666666*((-4.898979485566356*
      g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 30*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*
      Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(
      Cos(ThetaWp()))) + 2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 3*Sqr(
      gN)*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*(-
      4.898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 2*Sqr(
      g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 3*Sqr(gN)*Sqr(Sin(ThetaWp())))*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjSDXVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.016666666666666666*(4*(-4.898979485566356
      *g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 2*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) + 3*Sqr(gN)*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,
      Conj(ZDX(gI1,j1))*ZDX(gI2,j1)) + (8*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(
      ThetaW())) + 3*gN*(4.898979485566356*g1*Sin(ThetaW())*Sin(2*ThetaWp()) + 9*
      gN*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(2*(2.449489742783178*g1*gN*Sin(ThetaW
      ())*Sin(2*ThetaWp()) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 2*
      Cos(ThetaW())*(3.1622776601683795*g2*gN*Cos(ThetaWp())*Sin(ThetaWp()) +
      3.872983346207417*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) + 3*Sqr(g1)*Sqr(
      Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 2*Sqr(gN)*Sqr(Sin(ThetaWp())))*SUM(j1,0
      ,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-9.797958971132712*g1*gN*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + 24*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW
      ())) + Sqr(gN)*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3
       + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.016666666666666666*((-4.898979485566356*
      g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 30*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + Cos(ThetaW())*(18.973665961010276*g2*gN*Cos
      (ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(
      ThetaWp()))) + 2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 3*Sqr(gN)*
      Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (
      19.595917942265423*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 32*
      Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 3*Sqr(gN)*Sqr(Sin(ThetaWp()
      )))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.016666666666666666*(-((30*g2*Cos(ThetaW()
      )*Cos(ThetaWp()) + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      9.486832980505138*gN*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))
      ) + 2*(7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 9.486832980505138
      *gN*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjSDXVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.016666666666666666*((15.491933384829668*
      g1*Cos(ThetaWp())*Sin(ThetaW()) - 18.973665961010276*gN*Sin(ThetaWp()))*SUM(
      j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)) + (15.491933384829668*g1*Cos(ThetaWp()
      )*Sin(ThetaW()) + 28.460498941515414*gN*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZDX(
      gI2,3 + j1))*ZDX(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(2*(-5*g2*Cos(ThetaW())*Cos(ThetaWp())
      + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) + 2.23606797749979*(
      6.928203230275509*g1*Cos(ThetaWp())*Sin(ThetaW()) - 1.4142135623730951*gN*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.016666666666666666*((30*g2*Cos(ThetaW())*
      Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      9.486832980505138*gN*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))
      - (30.983866769659336*g1*Cos(ThetaWp())*Sin(ThetaW()) + 9.486832980505138*gN
      *Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(Conj(ZN(gI2,2))*(-10*g2*Cos(ThetaW())
      *Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) +
      9.486832980505138*gN*Sin(ThetaWp()))*ZN(gI1,2) + 2*Conj(ZN(gI2,3))*(5*g2*Cos
      (ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()
      ) + 3.1622776601683795*gN*Sin(ThetaWp()))*ZN(gI1,3) - 15.811388300841898*gN*
      Conj(ZN(gI2,4))*Sin(ThetaWp())*ZN(gI1,4));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(Conj(ZN(gI1,2))*(10*g2*Cos(ThetaW())*
      Cos(ThetaWp()) + 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      9.486832980505138*gN*Sin(ThetaWp()))*ZN(gI2,2) - 2*Conj(ZN(gI1,3))*(5*g2*Cos
      (ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()
      ) + 3.1622776601683795*gN*Sin(ThetaWp()))*ZN(gI2,3) + 15.811388300841898*gN*
      Conj(ZN(gI1,4))*Sin(ThetaWp())*ZN(gI2,4));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjVWmVZ(int gI2) const
{
   
   const std::complex<double> result = 0.05*g2*(vd*(7.745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) - 9.486832980505138*gN*Sin(ThetaWp()))*ZP(gI2,0) -
      2*vu*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795
      *gN*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhVZVZ(int gI2) const
{
   
   const std::complex<double> result = 0.05*(vd*(-14.696938456699067*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + Cos(ThetaW())*(-18.973665961010276*g2*gN*Cos(ThetaWp())*
      Sin(ThetaWp()) + 15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Cos(ThetaWp())))
      + 6*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 9*Sqr(gN)*Sqr(Sin(
      ThetaWp())))*ZH(gI2,0) + 2*vu*(3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + gN*(2.449489742783178*g1*Sin(ThetaW())*Sin(2*ThetaWp()) + 2*gN
      *Sqr(Sin(ThetaWp()))))*ZH(gI2,1) + 25*vs*Sqr(gN)*Sqr(Sin(ThetaWp()))*ZH(gI2,
      2));

   return result;
}

std::complex<double> CLASSNAME::CphhVZVZp(int gI2) const
{
   
   const std::complex<double> result = 0.05*(-(vd*(9.486832980505138*g2*gN*Cos(
      ThetaW())*Cos(2*ThetaWp()) - 9*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(gN) + 5*Sin
      (2*ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) + 7.348469228349534*g1*gN*Sin(
      ThetaW())*Sqr(Cos(ThetaWp())) + g1*(3.872983346207417*g2*Sin(2*ThetaW())*Sin
      (2*ThetaWp()) + 3*g1*Sin(2*ThetaWp())*Sqr(Sin(ThetaW())) - 7.348469228349534
      *gN*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*ZH(gI2,0)) + vu*(6.324555320336759*
      g2*gN*Cos(ThetaW())*Cos(2*ThetaWp()) + 4*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(
      gN) - 5*Sin(2*ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) + 4.898979485566356*g1*
      gN*Sin(ThetaW())*Sqr(Cos(ThetaWp())) - g1*(3.872983346207417*g2*Sin(2*ThetaW
      ())*Sin(2*ThetaWp()) + 3*g1*Sin(2*ThetaWp())*Sqr(Sin(ThetaW())) +
      4.898979485566356*gN*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*ZH(gI2,1) + 25*vs*
      Cos(ThetaWp())*Sin(ThetaWp())*Sqr(gN)*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ1() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ2() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ3() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpbargWmgWmVZp() const
{
   
   const double result = g2*Cos(ThetaW())*Sin(ThetaWp());

   return result;
}

double CLASSNAME::CpbargWmCgWmCVZp() const
{
   
   const double result = -(g2*Cos(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZp() const
{
   
   const double result = g2*Cos(ThetaW())*Sin(ThetaWp());

   return result;
}

double CLASSNAME::CpbarChaPChaPVZpPL() const
{
   
   const double result = 0.1*(-3.1622776601683795*gN*Cos(ThetaWp()) - 5*g2*Cos(
      ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())
      );

   return result;
}

double CLASSNAME::CpbarChaPChaPVZpPR() const
{
   
   const double result = 0.1*(-3.1622776601683795*gN*Cos(ThetaWp()) - 5*g2*Cos(
      ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())
      );

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-9.486832980505138*g2*gN*Cos(ThetaW(
      ))*Sin(2*ThetaWp()) + 9*Sqr(gN)*Sqr(Cos(ThetaWp())) - 15.491933384829668*g1*
      g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 10*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(2.449489742783178*gN*
      Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*ZP(gI1,0)*ZP(gI2
      ,0) + 2*(-4.898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()
      ) + 2*Sqr(gN)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) + Cos(ThetaW(
      ))*(3.1622776601683795*g2*gN*Sin(2*ThetaWp()) - 7.745966692414834*g1*g2*Sin(
      ThetaW())*Sqr(Sin(ThetaWp()))))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0conjSHp0VZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.0125*(-9.797958971132712*g1*gN*Cos(ThetaW
      () - 2*ThetaWp()) + 9.797958971132712*g1*gN*Cos(ThetaW() + 2*ThetaWp()) +
      15.491933384829668*g1*g2*Sin(2*ThetaW()) + 12.649110640673518*g2*gN*Sin(
      ThetaW() - 2*ThetaWp()) - 7.745966692414834*g1*g2*Sin(2*(ThetaW() - ThetaWp(
      ))) - 7.745966692414834*g1*g2*Sin(2*(ThetaW() + ThetaWp())) -
      12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) + 6*Sqr(g1) + 3*Cos(2*(
      ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) + 3*Cos(2*(
      ThetaW() + ThetaWp()))*Sqr(g1) + 10*Sqr(g2) - 5*Cos(2*(ThetaW() - ThetaWp())
      )*Sqr(g2) - 10*Cos(2*ThetaWp())*Sqr(g2) - 5*Cos(2*(ThetaW() + ThetaWp()))*
      Sqr(g2) + Cos(2*ThetaW())*(-6*Sqr(g1) + 10*Sqr(g2)) + 8*Sqr(gN) + 8*Cos(2*
      ThetaWp())*Sqr(gN))*(Conj(UHp0(gI1,0))*UHp0(gI2,0) + Conj(UHp0(gI1,1))*UHp0(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjSHppVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.0125*(-9.797958971132712*g1*gN*Cos(ThetaW
      () - 2*ThetaWp()) + 9.797958971132712*g1*gN*Cos(ThetaW() + 2*ThetaWp()) -
      15.491933384829668*g1*g2*Sin(2*ThetaW()) - 12.649110640673518*g2*gN*Sin(
      ThetaW() - 2*ThetaWp()) + 7.745966692414834*g1*g2*Sin(2*(ThetaW() - ThetaWp(
      ))) + 7.745966692414834*g1*g2*Sin(2*(ThetaW() + ThetaWp())) +
      12.649110640673518*g2*gN*Sin(ThetaW() + 2*ThetaWp()) + 6*Sqr(g1) + 3*Cos(2*(
      ThetaW() - ThetaWp()))*Sqr(g1) - 6*Cos(2*ThetaWp())*Sqr(g1) + 3*Cos(2*(
      ThetaW() + ThetaWp()))*Sqr(g1) + 10*Sqr(g2) - 5*Cos(2*(ThetaW() - ThetaWp())
      )*Sqr(g2) - 10*Cos(2*ThetaWp())*Sqr(g2) - 5*Cos(2*(ThetaW() + ThetaWp()))*
      Sqr(g2) + Cos(2*ThetaW())*(-6*Sqr(g1) + 10*Sqr(g2)) + 8*Sqr(gN) + 8*Cos(2*
      ThetaWp())*Sqr(gN))*(Conj(UHpp(gI1,0))*UHpp(gI2,0) + Conj(UHpp(gI1,1))*UHpp(
      gI2,1));

   return result;
}

double CLASSNAME::CpSSI0conjSSI0VZpVZp(int gI1, int gI2) const
{
   
   const double result = 1.25*KroneckerDelta(gI1,gI2)*Sqr(gN)*Sqr(Cos(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*((-9.486832980505138*gN*Cos(ThetaWp())
      + 2*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())
      )*ZP(gI1,0)*ZP(gI2,0) + 2*(3.1622776601683795*gN*Cos(ThetaWp()) + (5*g2*Cos(
      ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZP(gI1,1)*ZP
      (gI2,1));

   return result;
}

double CLASSNAME::CpSHp0conjSHp0VZp(int gI2, int gI1) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpSHppconjSHppVZp(int gI2, int gI1) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpSSI0conjSSI0VZp(int gI2, int gI1) const
{
   
   const double result = 0.7905694150420949*gN*Cos(ThetaWp())*KroneckerDelta(gI1,
      gI2);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Cos(ThetaW())*Sin(
      ThetaWp())*UM(gI1,0)) + 0.05*Conj(UM(gI2,1))*(9.486832980505138*gN*Cos(
      ThetaWp()) + 2*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*
      Sin(ThetaWp()))*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZpPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*Cos(ThetaW())*Sin(
      ThetaWp())*UP(gI2,0)) - 0.1*Conj(UP(gI1,1))*(3.1622776601683795*gN*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*UP(gI2,1);

   return result;
}

double CLASSNAME::CpbarChaIChaIVZpPL(int gI1, int gI2) const
{
   
   const double result = 0.05*KroneckerDelta(gI1,gI2)*(9.486832980505138*gN*Cos(
      ThetaWp()) + 2*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarChaIChaIVZpPR(int gI1, int gI2) const
{
   
   const double result = -0.1*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChiPVZpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.1*(3.1622776601683795*gN*Cos(ThetaWp())
      - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*
      (Conj(ZNp(gI2,0))*ZNp(gI1,0) - Conj(ZNp(gI2,1))*ZNp(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChiPVZpPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*(3.1622776601683795*gN*Cos(ThetaWp()) -
      (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*(
      Conj(ZNp(gI1,0))*ZNp(gI2,0) - Conj(ZNp(gI1,1))*ZNp(gI2,1));

   return result;
}

double CLASSNAME::CpFSIFSIVZpPL(int gI1, int gI2) const
{
   
   const double result = -0.7905694150420949*gN*Cos(ThetaWp())*KroneckerDelta(gI1,
      gI2);

   return result;
}

double CLASSNAME::CpFSIFSIVZpPR(int gI1, int gI2) const
{
   
   const double result = 0.7905694150420949*gN*Cos(ThetaWp())*KroneckerDelta(gI1,
      gI2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp())) + 10*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(
      2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      ))) + Cos(ThetaW())*(9.486832980505138*g2*gN*Sin(2*ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*ZA(gI1,0)*ZA(
      gI2,0) + 2*(-4.898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) - 3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 2*Sqr
      (gN)*Sqr(Cos(ThetaWp())) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW
      ()) + 3*g1*Sin(ThetaW()))*Sqr(Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))
      *Sqr(Sin(ThetaWp())))*ZA(gI1,1)*ZA(gI2,1) + 25*Sqr(gN)*Sqr(Cos(ThetaWp()))*
      ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp())) + 10*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(
      2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      ))) + Cos(ThetaW())*(9.486832980505138*g2*gN*Sin(2*ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*ZH(gI1,0)*ZH(
      gI2,0) + 2*(-4.898979485566356*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) - 3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 2*Sqr
      (gN)*Sqr(Cos(ThetaWp())) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW
      ()) + 3*g1*Sin(ThetaW()))*Sqr(Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))
      *Sqr(Sin(ThetaWp())))*ZH(gI1,1)*ZH(gI2,1) + 25*Sqr(gN)*Sqr(Cos(ThetaWp()))*
      ZH(gI1,2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(gI1,gI2)*(-
      3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 2*Sqr(gN)*Sqr(Cos(
      ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*(-
      2.449489742783178*gN*Sin(ThetaW())*Sin(2*ThetaWp()) + 3.872983346207417*g2*
      Sin(2*ThetaW())*Sqr(Sin(ThetaWp())) + 3*g1*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.05)*((
      9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZA(gI2,0)*ZH(gI1,0) + 2*
      (3.1622776601683795*gN*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZA(gI2,1)*ZH(gI1,1) -
      15.811388300841898*gN*Cos(ThetaWp())*ZA(gI2,2)*ZH(gI1,2));

   return result;
}

double CLASSNAME::CpSvconjSvVZp(int gI2, int gI1) const
{
   
   const double result = 0.1*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdFdVZpPL(int gI1, int gI2) const
{
   
   const double result = -0.016666666666666666*KroneckerDelta(gI1,gI2)*(
      9.486832980505138*gN*Cos(ThetaWp()) + 2*(15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdFdVZpPR(int gI1, int gI2) const
{
   
   const double result = KroneckerDelta(gI1,gI2)*(0.31622776601683794*gN*Cos(
      ThetaWp()) + 0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFDXFDXVZpPL(int gI1, int gI2) const
{
   
   const double result = KroneckerDelta(gI1,gI2)*(0.31622776601683794*gN*Cos(
      ThetaWp()) + 0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFDXFDXVZpPR(int gI1, int gI2) const
{
   
   const double result = -0.016666666666666666*KroneckerDelta(gI1,gI2)*(
      28.460498941515414*gN*Cos(ThetaWp()) - 15.491933384829668*g1*Sin(ThetaW())*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZpPL(int gI1, int gI2) const
{
   
   const double result = -0.1*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZpPR(int gI1, int gI2) const
{
   
   const double result = 0.05*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) + 15.491933384829668*g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZpPL(int gI1, int gI2) const
{
   
   const double result = -0.016666666666666666*KroneckerDelta(gI1,gI2)*(
      9.486832980505138*gN*Cos(ThetaWp()) + 2*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZpPR(int gI1, int gI2) const
{
   
   const double result = 0.016666666666666666*KroneckerDelta(gI1,gI2)*(
      9.486832980505138*gN*Cos(ThetaWp()) - 30.983866769659336*g1*Sin(ThetaW())*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZpPL(int gI1, int gI2) const
{
   
   const double result = -0.1*KroneckerDelta(gI1,gI2)*(3.1622776601683795*gN*Cos(
      ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZpPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjSHI0VZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((9*Sqr(gN)*Sqr(Cos(ThetaWp())) + 10*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(
      2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      ))) + Cos(ThetaW())*(9.486832980505138*g2*gN*Sin(2*ThetaWp()) +
      15.491933384829668*g1*g2*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,1,Conj
      (UHI0(gI1,j1))*UHI0(gI2,j1)) + 2*(-4.898979485566356*g1*gN*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) - 3.1622776601683795*g2*gN*Cos(ThetaW())*Sin(2*
      ThetaWp()) + 2*Sqr(gN)*Sqr(Cos(ThetaWp())) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Sin(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())))*SUM(j1,0,1,Conj(UHI0(gI1
      ,2 + j1))*UHI0(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjSHIpVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((-9.486832980505138*g2*gN*Cos(ThetaW(
      ))*Sin(2*ThetaWp()) + 9*Sqr(gN)*Sqr(Cos(ThetaWp())) - 15.491933384829668*g1*
      g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 10*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Sin(ThetaWp())) + 3*g1*Sin(ThetaW())*(2.449489742783178*gN*
      Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,1,Conj(
      UHIp(gI1,j1))*UHIp(gI2,j1)) + 2*(-4.898979485566356*g1*gN*Cos(ThetaWp())*Sin
      (ThetaW())*Sin(ThetaWp()) + 2*Sqr(gN)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(
      Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp())) + Cos(ThetaW())*(3.1622776601683795*g2*gN*Sin(2*ThetaWp()) -
      7.745966692414834*g1*g2*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,1,Conj(
      UHIp(gI1,2 + j1))*UHIp(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHI0conjSHI0VZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(-((9.486832980505138*gN*Cos(ThetaWp()
      ) + 2*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp(
      )))*SUM(j1,0,1,Conj(UHI0(gI2,j1))*UHI0(gI1,j1))) + 2*(3.1622776601683795*gN*
      Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*
      Sin(ThetaWp()))*SUM(j1,0,1,Conj(UHI0(gI2,2 + j1))*UHI0(gI1,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjSHIpVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*((-9.486832980505138*gN*Cos(ThetaWp())
      + 2*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())
      )*SUM(j1,0,1,Conj(UHIp(gI2,j1))*UHIp(gI1,j1)) + 2*(3.1622776601683795*gN*Cos
      (ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*UHIp(gI1,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIVZpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*((9.486832980505138*gN*Cos(ThetaWp())
      + 2*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())
      )*SUM(j1,0,1,Conj(ZNI(gI2,j1))*ZNI(gI1,j1)) + 2*(3.1622776601683795*gN*Cos(
      ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*SUM(j1,0,1,Conj(ZNI(gI2,2 + j1))*ZNI(gI1,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChiIVZpPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-((9.486832980505138*gN*Cos(ThetaWp()
      ) + 2*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp(
      )))*SUM(j1,0,1,Conj(ZNI(gI1,j1))*ZNI(gI2,j1))) + 2*(-3.1622776601683795*gN*
      Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*
      Sin(ThetaWp()))*SUM(j1,0,1,Conj(ZNI(gI1,2 + j1))*ZNI(gI2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.016666666666666666*((18.973665961010276*
      g2*gN*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 3*Sqr(gN)*Sqr(Cos(
      ThetaWp())) + 15.491933384829668*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(
      ThetaWp())) + 30*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*Sin(
      ThetaW())*(2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(
      Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*(3*Sqr(gN)*Sqr
      (Cos(ThetaWp())) + g1*Sin(ThetaW())*(2.449489742783178*gN*Sin(2*ThetaWp()) +
      2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(
      gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjSDXVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.016666666666666666*(4*(3*Sqr(gN)*Sqr(Cos(
      ThetaWp())) + g1*Sin(ThetaW())*(2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1
      *Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZDX(gI1,j1))*ZDX(gI2,j1
      )) + (27*Sqr(gN)*Sqr(Cos(ThetaWp())) + 2*g1*Sin(ThetaW())*(-
      7.348469228349534*gN*Sin(2*ThetaWp()) + 4*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      ))))*SUM(j1,0,2,Conj(ZDX(gI1,3 + j1))*ZDX(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(2*(-4.898979485566356*g1*gN*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 3.1622776601683795*g2*gN*Cos(
      ThetaW())*Sin(2*ThetaWp()) + 2*Sqr(gN)*Sqr(Cos(ThetaWp())) -
      7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 5*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())
      )*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (
      9.797958971132712*g1*gN*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + Sqr(gN
      )*Sqr(Cos(ThetaWp())) + 24*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())))*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZpVZp(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.016666666666666666*((-9.486832980505138*
      g2*gN*Cos(ThetaW())*Sin(2*ThetaWp()) + 3*Sqr(gN)*Sqr(Cos(ThetaWp())) -
      15.491933384829668*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) +
      30*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*Sin(ThetaW())*(
      2.449489742783178*gN*Sin(2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()
      ))))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + (-19.595917942265423*g1*gN*
      Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 3*Sqr(gN)*Sqr(Cos(ThetaWp()))
      + 32*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZU(gI1,
      3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.016666666666666666*((9.486832980505138*gN
      *Cos(ThetaWp()) + 2*(15*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()
      ))*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) - 2*(
      9.486832980505138*gN*Cos(ThetaWp()) + 7.745966692414834*g1*Sin(ThetaW())*Sin
      (ThetaWp()))*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSDXconjSDXVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.016666666666666666*(-2*(9.486832980505138
      *gN*Cos(ThetaWp()) + 7.745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM(
      j1,0,2,Conj(ZDX(gI2,j1))*ZDX(gI1,j1)) + (28.460498941515414*gN*Cos(ThetaWp()
      ) - 15.491933384829668*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZDX(
      gI2,3 + j1))*ZDX(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(2*(3.1622776601683795*gN*Cos(ThetaWp(
      )) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()
      ))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) - (3.1622776601683795*gN*Cos(
      ThetaWp()) + 15.491933384829668*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM(j1,0,2,
      Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZp(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.016666666666666666*((9.486832980505138*gN
      *Cos(ThetaWp()) + 2*(-15*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW(
      )))*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) + (-
      9.486832980505138*gN*Cos(ThetaWp()) + 30.983866769659336*g1*Sin(ThetaW())*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZpPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(Conj(ZN(gI2,2))*(9.486832980505138*gN
      *Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW())
      )*Sin(ThetaWp()))*ZN(gI1,2) + 2*Conj(ZN(gI2,3))*(3.1622776601683795*gN*Cos(
      ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*ZN(gI1,3) - 15.811388300841898*gN*Conj(ZN(gI2,4))*Cos(ThetaWp())
      *ZN(gI1,4));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZpPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-(Conj(ZN(gI1,2))*(9.486832980505138*
      gN*Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW(
      )))*Sin(ThetaWp()))*ZN(gI2,2)) + 2*Conj(ZN(gI1,3))*(-3.1622776601683795*gN*
      Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*
      Sin(ThetaW())*Sin(ThetaWp()))*ZN(gI2,3) + 15.811388300841898*gN*Conj(ZN(gI1,
      4))*Cos(ThetaWp())*ZN(gI2,4));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjVWmVZp(int gI2) const
{
   
   const std::complex<double> result = -0.05*g2*(vd*(9.486832980505138*gN*Cos(
      ThetaWp()) + 7.745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI2,0) +
      2*vu*(3.1622776601683795*gN*Cos(ThetaWp()) - 3.872983346207417*g1*Sin(ThetaW
      ())*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhVZpVZp(int gI2) const
{
   
   const std::complex<double> result = 0.05*(vd*(9.486832980505138*gN*(g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 9*Sqr(gN
      )*Sqr(Cos(ThetaWp())) + 10*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())))*ZH(gI2,0) + 2*vu*(-3.1622776601683795*gN*(g2
      *Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 2*
      Sqr(gN)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1
      *Sin(ThetaW()))*Sqr(Sin(ThetaWp())))*ZH(gI2,1) + 25*vs*Sqr(gN)*Sqr(Cos(
      ThetaWp()))*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZpVZp1() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZpVZp2() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZpVZp3() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbargPgWmconjVWm() const
{
   
   const double result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbargWmCgPconjVWm() const
{
   
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmCgZconjVWm() const
{
   
   const double result = -(g2*Cos(ThetaW())*Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpbargWmCgZpconjVWm() const
{
   
   const double result = g2*Cos(ThetaW())*Sin(ThetaWp());

   return result;
}

double CLASSNAME::CpbargZgWmconjVWm() const
{
   
   const double result = g2*Cos(ThetaW())*Cos(ThetaWp());

   return result;
}

double CLASSNAME::CpbargZpgWmconjVWm() const
{
   
   const double result = -(g2*Cos(ThetaW())*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1
      )*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHp0conjSHp0conjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(Conj(UHp0(gI1,0))*UHp0(gI2,0)
      + Conj(UHp0(gI1,1))*UHp0(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjSHppconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(Conj(UHpp(gI1,0))*UHpp(gI2,0)
      + Conj(UHpp(gI1,1))*UHpp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSHppconjSHp0conjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*(Conj(UHpp(gI2,0))*
      UHp0(gI1,0) - Conj(UHpp(gI2,1))*UHp0(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(ZA(gI2,0)*
      ZP(gI1,0) + ZA(gI2,1)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP
      (gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPChaPconjVWmPL(int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*ZNp(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiPChaPconjVWmPR(int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*Conj(ZNp(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1
      )*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1
      )*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpSvconjSvconjVWmVWm(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZDL(
      gI2,j1))*ZUL(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFuFdconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-0.7071067811865475*g2*Conj(ZEL(
      gI2,gI1)),0);

   return result;
}

double CLASSNAME::CpbarFvFeconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSvconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZE(
      gI2,j1))*ZV(gI1,j1));

   return result;
}

double CLASSNAME::CpSHI0conjSHI0conjVWmVWm(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

double CLASSNAME::CpSHIpconjSHIpconjVWmVWm(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpChiIChaIconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,1,Conj(ZMI(
      gI2,j1))*ZNI(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiIChaIconjVWmPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,1,Conj(ZNI(
      gI1,2 + j1))*ZPI(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSHIpconjSHI0conjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*(SUM(j1,0,1,Conj(UHIp
      (gI2,j1))*UHI0(gI1,j1)) - SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*UHI0(gI1,2 + j1)
      ));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(
      gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(
      gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,1) +
      1.4142135623730951*Conj(UM(gI2,1))*ZN(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjVWmPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(ZN(gI1,1))*UP(gI2,0)) +
      0.7071067811865475*g2*Conj(ZN(gI1,3))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSuconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*ZU(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CphhconjVWmVWm(int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpconjVWmconjVWmVWmVWm1() const
{
   
   const double result = 2*Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWmconjVWmVWmVWm2() const
{
   
   const double result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpconjVWmconjVWmVWmVWm3() const
{
   
   const double result = -Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiHpmPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(3,gO2)*
      ZP(gI2,1)) - 0.1*Conj(UP(gI1,1))*(10*KroneckerDelta(4,gO2)*Lambdax*ZP(gI2,0)
      + (5.477225575051661*g1*KroneckerDelta(0,gO2) + 7.0710678118654755*g2*
      KroneckerDelta(1,gO2) - 4.47213595499958*gN*KroneckerDelta(5,gO2))*ZP(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiHpmPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(2,gO1)*UM(gI1,0)*ZP(gI2
      ,0)) + 0.1*UM(gI1,1)*(5.477225575051661*g1*KroneckerDelta(0,gO1)*ZP(gI2,0) +
      7.0710678118654755*g2*KroneckerDelta(1,gO1)*ZP(gI2,0) + 6.708203932499369*gN
      *KroneckerDelta(5,gO1)*ZP(gI2,0) - 10*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZP
      (gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*KroneckerDelta(2,gO2)*
      ZP(gI1,0)) + 0.1*Conj(UM(gI2,1))*(5.477225575051661*g1*KroneckerDelta(0,gO2)
      *ZP(gI1,0) + 7.0710678118654755*g2*KroneckerDelta(1,gO2)*ZP(gI1,0) +
      6.708203932499369*gN*KroneckerDelta(5,gO2)*ZP(gI1,0) - 10*KroneckerDelta(4,
      gO2)*Lambdax*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(Conj(Lambdax)*KroneckerDelta(4,gO1)*UP(
      gI2,1)*ZP(gI1,0)) - 0.1*(10*g2*KroneckerDelta(3,gO1)*UP(gI2,0) + (
      5.477225575051661*g1*KroneckerDelta(0,gO1) + 7.0710678118654755*g2*
      KroneckerDelta(1,gO1) - 4.47213595499958*gN*KroneckerDelta(5,gO1))*UP(gI2,1)
      )*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiPconjSHp0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*Conj(ZNp(gI2,0))*(5.477225575051661*g1*
      KroneckerDelta(0,gO2) - 7.0710678118654755*g2*KroneckerDelta(1,gO2) -
      4.47213595499958*gN*KroneckerDelta(5,gO2))*UHp0(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiPconjSHp0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*(-5.477225575051661*g1*KroneckerDelta(0
      ,gO1) + 7.0710678118654755*g2*KroneckerDelta(1,gO1) + 4.47213595499958*gN*
      KroneckerDelta(5,gO1))*UHp0(gI1,1)*ZNp(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiPSHp0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*Conj(UHp0(gI1,1))*Conj(ZNp(gI2,1))*(-
      5.477225575051661*g1*KroneckerDelta(0,gO2) + 7.0710678118654755*g2*
      KroneckerDelta(1,gO2) + 4.47213595499958*gN*KroneckerDelta(5,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiPSHp0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*Conj(UHp0(gI1,0))*(5.477225575051661*g1
      *KroneckerDelta(0,gO1) - 7.0710678118654755*g2*KroneckerDelta(1,gO1) -
      4.47213595499958*gN*KroneckerDelta(5,gO1))*ZNp(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiFSIconjSSI0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.118033988749895*gN*KroneckerDelta(5,gO2)
      *SUM(j1,0,1,Conj(ZFSI(gI2,j1))*ZSSI(gI1,j1));

   return result;
}

double CLASSNAME::CpUChiFSIconjSSI0PR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpUChiFSISSI0PL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpUChiFSISSI0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.118033988749895*gN*KroneckerDelta(5,gO1)
      *SUM(j1,0,1,Conj(ZSSI(gI1,j1))*ZFSI(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIUChiSHIpPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.5477225575051661*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*Conj(ZPI(gI1,j1))) - 0.7071067811865475*
      g2*KroneckerDelta(1,gO2)*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*Conj(ZPI(gI1,j1))
      ) + 0.4472135954999579*gN*KroneckerDelta(5,gO2)*SUM(j1,0,1,Conj(UHIp(gI2,2 +
      j1))*Conj(ZPI(gI1,j1))) - KroneckerDelta(4,gO2)*SUM(j1,0,1,Conj(UHIp(gI2,j1)
      )*Conj(ZPI(gI1,j1))*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIUChiSHIpPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO1)
      *SUM(j1,0,1,Conj(UHIp(gI2,j1))*ZMI(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,1,Conj(UHIp(gI2,j1))*ZMI(gI1,j1)) +
      0.6708203932499369*gN*KroneckerDelta(5,gO1)*SUM(j1,0,1,Conj(UHIp(gI2,j1))*
      ZMI(gI1,j1)) - KroneckerDelta(4,gO1)*SUM(j1,0,1,Conj(UHIp(gI2,2 + j1))*Conj(
      Lambda12(j1,j1))*ZMI(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiVWmPL(int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*KroneckerDelta(1,gO2)*UM(gI1,0)
      + 1.4142135623730951*KroneckerDelta(2,gO2)*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiVWmPR(int gI1, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(1,gO1))
      + 0.7071067811865475*g2*Conj(UP(gI1,1))*KroneckerDelta(3,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaPconjSHppPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.1*(5.477225575051661*g1*KroneckerDelta(0,
      gO2) + 7.0710678118654755*g2*KroneckerDelta(1,gO2) - 4.47213595499958*gN*
      KroneckerDelta(5,gO2))*UHpp(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaPconjSHppPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = -0.1*(5.477225575051661*g1*KroneckerDelta(0
      ,gO1) + 7.0710678118654755*g2*KroneckerDelta(1,gO1) - 4.47213595499958*gN*
      KroneckerDelta(5,gO1))*UHpp(gI1,1);

   return result;
}

double CLASSNAME::CpbarFvUChiSvPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvUChiSvPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,0.5477225575051661*g1*Conj(ZV(
      gI2,gI1))*KroneckerDelta(0,gO1),0) + IF(gI1 < 3,-0.7071067811865475*g2*Conj(
      ZV(gI2,gI1))*KroneckerDelta(1,gO1),0) + IF(gI1 < 3,-0.4472135954999579*gN*
      Conj(ZV(gI2,gI1))*KroneckerDelta(5,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiFvconjSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5477225575051661*g1*
      KroneckerDelta(0,gO2)*ZV(gI1,gI2),0) + IF(gI2 < 3,-0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZV(gI1,gI2),0) + IF(gI2 < 3,-0.4472135954999579*gN*
      KroneckerDelta(5,gO2)*ZV(gI1,gI2),0);

   return result;
}

double CLASSNAME::CpUChiFvconjSvPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdUChiSdPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.3651483716701107*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1))) - 0.4472135954999579*gN
      *KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1))) -
      KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*Conj(ZDR(gI1,j1))*Yd(j1,j1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdUChiSdPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1)) -
      0.22360679774997896*gN*KroneckerDelta(5,gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL
      (gI1,j1)) - KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI2,3 +
      j1))*ZDL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXUChiSDXPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.3651483716701107*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(ZDXR(gI1,j1))) + 0.6708203932499369*
      gN*KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(ZDXR(gI1,j1))
      ) - KroneckerDelta(4,gO2)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*Conj(ZDXR(gI1,j1))*
      Kappa(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXUChiSDXPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = 0.3651483716701107*g1*KroneckerDelta(0,gO1)
      *SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDXL(gI1,j1)) + 0.4472135954999579*gN*
      KroneckerDelta(5,gO1)*SUM(j1,0,2,Conj(ZDX(gI2,j1))*ZDXL(gI1,j1)) -
      KroneckerDelta(4,gO1)*SUM(j1,0,2,Conj(ZDX(gI2,3 + j1))*Conj(Kappa(j1,j1))*
      ZDXL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUChiSePL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -1.0954451150103321*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(ZER(gI1,j1))) - 0.22360679774997896*
      gN*KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(ZER(gI1,j1)))
      - KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*Conj(ZER(gI1,j1))*Ye(j1,
      j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUChiSePR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO1)
      *SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(gI1,j1)) -
      0.4472135954999579*gN*KroneckerDelta(5,gO1)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(
      gI1,j1)) - KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI2,3 +
      j1))*ZEL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuUChiSuPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.7302967433402214*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1))) - 0.22360679774997896*gN
      *KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1))) -
      KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*Conj(ZUR(gI1,j1))*Yu(j1,j1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuUChiSuPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1)) - 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1)) -
      0.22360679774997896*gN*KroneckerDelta(5,gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL
      (gI1,j1)) - KroneckerDelta(3,gO1)*SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI2,3 +
      j1))*ZUL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChihhPL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(-10*g2*Conj(ZN(gI2,1))*KroneckerDelta
      (2,gO2)*ZH(gI1,0) + 9.486832980505138*gN*Conj(ZN(gI2,5))*KroneckerDelta(2,
      gO2)*ZH(gI1,0) + 14.142135623730951*Conj(ZN(gI2,4))*KroneckerDelta(3,gO2)*
      Lambdax*ZH(gI1,0) + 14.142135623730951*Conj(ZN(gI2,3))*KroneckerDelta(4,gO2)
      *Lambdax*ZH(gI1,0) - 7.745966692414834*g1*Conj(ZN(gI2,3))*KroneckerDelta(0,
      gO2)*ZH(gI1,1) + 10*g2*Conj(ZN(gI2,3))*KroneckerDelta(1,gO2)*ZH(gI1,1) + 10*
      g2*Conj(ZN(gI2,1))*KroneckerDelta(3,gO2)*ZH(gI1,1) + 6.324555320336759*gN*
      Conj(ZN(gI2,5))*KroneckerDelta(3,gO2)*ZH(gI1,1) + 6.324555320336759*gN*Conj(
      ZN(gI2,3))*KroneckerDelta(5,gO2)*ZH(gI1,1) + 14.142135623730951*Conj(ZN(gI2,
      4))*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1) + 7.745966692414834*g1*Conj(ZN(
      gI2,0))*(KroneckerDelta(2,gO2)*ZH(gI1,0) - KroneckerDelta(3,gO2)*ZH(gI1,1))
      - 15.811388300841898*gN*Conj(ZN(gI2,5))*KroneckerDelta(4,gO2)*ZH(gI1,2) -
      15.811388300841898*gN*Conj(ZN(gI2,4))*KroneckerDelta(5,gO2)*ZH(gI1,2) +
      14.142135623730951*Conj(ZN(gI2,3))*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,2) +
      Conj(ZN(gI2,2))*(7.745966692414834*g1*KroneckerDelta(0,gO2)*ZH(gI1,0) - 10*
      g2*KroneckerDelta(1,gO2)*ZH(gI1,0) + 9.486832980505138*gN*KroneckerDelta(5,
      gO2)*ZH(gI1,0) + 14.142135623730951*Lambdax*(KroneckerDelta(4,gO2)*ZH(gI1,1)
      + KroneckerDelta(3,gO2)*ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChihhPR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*KroneckerDelta(0
      ,gO1)*ZH(gI1,0)*ZN(gI2,2) - 10*g2*KroneckerDelta(1,gO1)*ZH(gI1,0)*ZN(gI2,2)
      + 9.486832980505138*gN*KroneckerDelta(5,gO1)*ZH(gI1,0)*ZN(gI2,2) +
      14.142135623730951*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZH(gI1,1)*ZN(gI2,2) +
      14.142135623730951*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZH(gI1,0)*ZN(gI2,3) -
      7.745966692414834*g1*KroneckerDelta(0,gO1)*ZH(gI1,1)*ZN(gI2,3) + 10*g2*
      KroneckerDelta(1,gO1)*ZH(gI1,1)*ZN(gI2,3) + 6.324555320336759*gN*
      KroneckerDelta(5,gO1)*ZH(gI1,1)*ZN(gI2,3) - 15.811388300841898*gN*
      KroneckerDelta(5,gO1)*ZH(gI1,2)*ZN(gI2,4) - 15.811388300841898*gN*
      KroneckerDelta(4,gO1)*ZH(gI1,2)*ZN(gI2,5) + 2*KroneckerDelta(3,gO1)*(
      7.0710678118654755*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,2) + ZH(gI1,0)*ZN(gI2,4))
      + ZH(gI1,1)*(-3.872983346207417*g1*ZN(gI2,0) + 5*g2*ZN(gI2,1) +
      3.1622776601683795*gN*ZN(gI2,5))) + KroneckerDelta(2,gO1)*(
      14.142135623730951*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,3) + ZH(gI1,1)*ZN(gI2,4))
      + ZH(gI1,0)*(7.745966692414834*g1*ZN(gI2,0) - 10*g2*ZN(gI2,1) +
      9.486832980505138*gN*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaIconjSHIpPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,1,Conj(ZMI(gI2,j1))*UHIp(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,1,Conj(ZMI(gI2,j1))*UHIp(gI1,j1)) +
      0.6708203932499369*gN*KroneckerDelta(5,gO2)*SUM(j1,0,1,Conj(ZMI(gI2,j1))*
      UHIp(gI1,j1)) - KroneckerDelta(4,gO2)*SUM(j1,0,1,Conj(ZMI(gI2,j1))*UHIp(gI1,
      2 + j1)*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaIconjSHIpPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(KroneckerDelta(4,gO1)*SUM(j1,0,1,Conj(
      Lambda12(j1,j1))*UHIp(gI1,j1)*ZPI(gI2,j1))) - 0.1*(5.477225575051661*g1*
      KroneckerDelta(0,gO1) + 7.0710678118654755*g2*KroneckerDelta(1,gO1) -
      4.47213595499958*gN*KroneckerDelta(5,gO1))*SUM(j1,0,1,UHIp(gI1,2 + j1)*ZPI(
      gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiIconjSHI0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,1,Conj(ZNI(gI2,j1))*UHI0(gI1,j1)) - 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,1,Conj(ZNI(gI2,j1))*UHI0(gI1,j1)) +
      0.6708203932499369*gN*KroneckerDelta(5,gO2)*SUM(j1,0,1,Conj(ZNI(gI2,j1))*
      UHI0(gI1,j1)) + KroneckerDelta(4,gO2)*SUM(j1,0,1,Conj(ZNI(gI2,j1))*UHI0(gI1,
      2 + j1)*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiIconjSHI0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = KroneckerDelta(4,gO1)*SUM(j1,0,1,Conj(
      Lambda12(j1,j1))*UHI0(gI1,j1)*ZNI(gI2,2 + j1)) + 0.1*(-5.477225575051661*g1*
      KroneckerDelta(0,gO1) + 7.0710678118654755*g2*KroneckerDelta(1,gO1) +
      4.47213595499958*gN*KroneckerDelta(5,gO1))*SUM(j1,0,1,UHI0(gI1,2 + j1)*ZNI(
      gI2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiISHI0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5477225575051661*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj(ZNI(gI2,2 + j1))) +
      0.7071067811865475*g2*KroneckerDelta(1,gO2)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1)
      )*Conj(ZNI(gI2,2 + j1))) + 0.4472135954999579*gN*KroneckerDelta(5,gO2)*SUM(
      j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj(ZNI(gI2,2 + j1))) + KroneckerDelta(4,gO2)
      *SUM(j1,0,1,Conj(UHI0(gI1,j1))*Conj(ZNI(gI2,2 + j1))*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiISHI0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO1)
      *SUM(j1,0,1,Conj(UHI0(gI1,j1))*ZNI(gI2,j1)) - 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*ZNI(gI2,j1)) +
      0.6708203932499369*gN*KroneckerDelta(5,gO1)*SUM(j1,0,1,Conj(UHI0(gI1,j1))*
      ZNI(gI2,j1)) + KroneckerDelta(4,gO1)*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*Conj(
      Lambda12(j1,j1))*ZNI(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiAhPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.05)*(10*g2*Conj(ZN
      (gI1,1))*KroneckerDelta(2,gO2)*ZA(gI2,0) - 9.486832980505138*gN*Conj(ZN(gI1,
      5))*KroneckerDelta(2,gO2)*ZA(gI2,0) + 14.142135623730951*Conj(ZN(gI1,4))*
      KroneckerDelta(3,gO2)*Lambdax*ZA(gI2,0) + 14.142135623730951*Conj(ZN(gI1,3))
      *KroneckerDelta(4,gO2)*Lambdax*ZA(gI2,0) + 7.745966692414834*g1*Conj(ZN(gI1,
      3))*KroneckerDelta(0,gO2)*ZA(gI2,1) - 10*g2*Conj(ZN(gI1,3))*KroneckerDelta(1
      ,gO2)*ZA(gI2,1) - 10*g2*Conj(ZN(gI1,1))*KroneckerDelta(3,gO2)*ZA(gI2,1) -
      6.324555320336759*gN*Conj(ZN(gI1,5))*KroneckerDelta(3,gO2)*ZA(gI2,1) -
      6.324555320336759*gN*Conj(ZN(gI1,3))*KroneckerDelta(5,gO2)*ZA(gI2,1) +
      14.142135623730951*Conj(ZN(gI1,4))*KroneckerDelta(2,gO2)*Lambdax*ZA(gI2,1) +
      7.745966692414834*g1*Conj(ZN(gI1,0))*(-(KroneckerDelta(2,gO2)*ZA(gI2,0)) +
      KroneckerDelta(3,gO2)*ZA(gI2,1)) + 15.811388300841898*gN*Conj(ZN(gI1,5))*
      KroneckerDelta(4,gO2)*ZA(gI2,2) + 15.811388300841898*gN*Conj(ZN(gI1,4))*
      KroneckerDelta(5,gO2)*ZA(gI2,2) + 14.142135623730951*Conj(ZN(gI1,3))*
      KroneckerDelta(2,gO2)*Lambdax*ZA(gI2,2) + Conj(ZN(gI1,2))*(-
      7.745966692414834*g1*KroneckerDelta(0,gO2)*ZA(gI2,0) + 10*g2*KroneckerDelta(
      1,gO2)*ZA(gI2,0) - 9.486832980505138*gN*KroneckerDelta(5,gO2)*ZA(gI2,0) +
      14.142135623730951*Lambdax*(KroneckerDelta(4,gO2)*ZA(gI2,1) + KroneckerDelta
      (3,gO2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiAhPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      7.745966692414834*g1*KroneckerDelta(0,gO1)*ZA(gI2,0)*ZN(gI1,2) - 10*g2*
      KroneckerDelta(1,gO1)*ZA(gI2,0)*ZN(gI1,2) + 9.486832980505138*gN*
      KroneckerDelta(5,gO1)*ZA(gI2,0)*ZN(gI1,2) - 14.142135623730951*Conj(Lambdax)
      *KroneckerDelta(4,gO1)*ZA(gI2,1)*ZN(gI1,2) - 14.142135623730951*Conj(Lambdax
      )*KroneckerDelta(4,gO1)*ZA(gI2,0)*ZN(gI1,3) - 7.745966692414834*g1*
      KroneckerDelta(0,gO1)*ZA(gI2,1)*ZN(gI1,3) + 10*g2*KroneckerDelta(1,gO1)*ZA(
      gI2,1)*ZN(gI1,3) + 6.324555320336759*gN*KroneckerDelta(5,gO1)*ZA(gI2,1)*ZN(
      gI1,3) - 15.811388300841898*gN*KroneckerDelta(5,gO1)*ZA(gI2,2)*ZN(gI1,4) -
      15.811388300841898*gN*KroneckerDelta(4,gO1)*ZA(gI2,2)*ZN(gI1,5) - 2*
      KroneckerDelta(3,gO1)*(7.0710678118654755*Conj(Lambdax)*(ZA(gI2,2)*ZN(gI1,2)
      + ZA(gI2,0)*ZN(gI1,4)) + ZA(gI2,1)*(3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN
      (gI1,1) - 3.1622776601683795*gN*ZN(gI1,5))) + KroneckerDelta(2,gO1)*(-
      14.142135623730951*Conj(Lambdax)*(ZA(gI2,2)*ZN(gI1,3) + ZA(gI2,1)*ZN(gI1,4))
      + ZA(gI2,0)*(7.745966692414834*g1*ZN(gI1,0) - 10*g2*ZN(gI1,1) +
      9.486832980505138*gN*ZN(gI1,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFdconjSdPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(gI1,j1)) -
      0.22360679774997896*gN*KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD
      (gI1,j1)) - KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*Yd(j1,j1)*ZD(
      gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFdconjSdPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(Yd(
      j1,j1))*ZD(gI1,j1)*ZDR(gI2,j1))) - 0.14907119849998596*(2.449489742783178*g1
      *KroneckerDelta(0,gO1) + 3*gN*KroneckerDelta(5,gO1))*SUM(j1,0,2,ZD(gI1,3 +
      j1)*ZDR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFDXconjSDXPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.3651483716701107*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,2,Conj(ZDXL(gI2,j1))*ZDX(gI1,j1)) + 0.4472135954999579*gN*
      KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZDXL(gI2,j1))*ZDX(gI1,j1)) -
      KroneckerDelta(4,gO2)*SUM(j1,0,2,Conj(ZDXL(gI2,j1))*ZDX(gI1,3 + j1)*Kappa(j1
      ,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFDXconjSDXPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(KroneckerDelta(4,gO1)*SUM(j1,0,2,Conj(
      Kappa(j1,j1))*ZDX(gI1,j1)*ZDXR(gI2,j1))) + 0.07453559924999298*(-
      4.898979485566356*g1*KroneckerDelta(0,gO1) + 9*gN*KroneckerDelta(5,gO1))*SUM
      (j1,0,2,ZDX(gI1,3 + j1)*ZDXR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFeconjSePL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(gI1,j1)) -
      0.4472135954999579*gN*KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(
      gI1,j1)) - KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*Ye(j1,j1)*ZE(
      gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFeconjSePR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(Ye(
      j1,j1))*ZE(gI1,j1)*ZER(gI2,j1))) - 0.22360679774997896*(4.898979485566356*g1
      *KroneckerDelta(0,gO1) + gN*KroneckerDelta(5,gO1))*SUM(j1,0,2,ZE(gI1,3 + j1)
      *ZER(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFuconjSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(gI1,j1)) - 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(gI1,j1)) -
      0.22360679774997896*gN*KroneckerDelta(5,gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU
      (gI1,j1)) - KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*Yu(j1,j1)*ZU(
      gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFuconjSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(KroneckerDelta(3,gO1)*SUM(j1,0,2,Conj(Yu(
      j1,j1))*ZU(gI1,j1)*ZUR(gI2,j1))) + 0.07453559924999298*(9.797958971132712*g1
      *KroneckerDelta(0,gO1) - 3*gN*KroneckerDelta(5,gO1))*SUM(j1,0,2,ZU(gI1,3 +
      j1)*ZUR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPUChiSHppPL(int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.1*Conj(UHpp(gI2,1))*(5.477225575051661*
      g1*KroneckerDelta(0,gO2) + 7.0710678118654755*g2*KroneckerDelta(1,gO2) -
      4.47213595499958*gN*KroneckerDelta(5,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPUChiSHppPR(int gO1, int gI2) const
{
   
   const std::complex<double> result = 0.1*Conj(UHpp(gI2,0))*(5.477225575051661*g1
      *KroneckerDelta(0,gO1) + 7.0710678118654755*g2*KroneckerDelta(1,gO1) -
      4.47213595499958*gN*KroneckerDelta(5,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)) +
      0.7071067811865475*g2*KroneckerDelta(3,gO2)*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(1
      ,gO1) + 1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(2,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(2,gO2)*(-10*g2*Cos(
      ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 9.486832980505138*gN*Sin(ThetaWp()))*ZN(gI2,2) + 2*KroneckerDelta(3,gO2)*(
      5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp()))*ZN(gI2,3) -
      15.811388300841898*gN*KroneckerDelta(4,gO2)*Sin(ThetaWp())*ZN(gI2,4));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.05*(15.811388300841898*gN*Conj(ZN(gI2,4))
      *KroneckerDelta(4,gO1)*Sin(ThetaWp()) + Conj(ZN(gI2,2))*KroneckerDelta(2,gO1
      )*(10*g2*Cos(ThetaW())*Cos(ThetaWp()) + 7.745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) - 9.486832980505138*gN*Sin(ThetaWp())) - 2*Conj(ZN(gI2,3))*
      KroneckerDelta(3,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417
      *g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZpPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(2,gO2)*(
      9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZN(gI2,2) + 2*
      KroneckerDelta(3,gO2)*(3.1622776601683795*gN*Cos(ThetaWp()) - (5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZN(gI2,3) -
      15.811388300841898*gN*Cos(ThetaWp())*KroneckerDelta(4,gO2)*ZN(gI2,4));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZpPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.05*(15.811388300841898*gN*Conj(ZN(gI2,4))
      *Cos(ThetaWp())*KroneckerDelta(4,gO1) + 2*Conj(ZN(gI2,3))*KroneckerDelta(3,
      gO1)*(-3.1622776601683795*gN*Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(ThetaWp
      ()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())) - Conj(ZN(gI2,2))*
      KroneckerDelta(2,gO1)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiPSHppPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*Conj(UHpp(gI1,1))*Conj(ZNp(gI2,1))*
      KroneckerDelta(0,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiPSHppPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*Conj(UHpp(gI1,0))*KroneckerDelta(0,gO1
      )*ZNp(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *(g2*Conj(UM(gI1,0))*KroneckerDelta(1,gO2)*ZA(gI2,1) + Conj(UM(gI1,1))*(g2*
      KroneckerDelta(0,gO2)*ZA(gI2,0) - KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*(g2*KroneckerDelta(0,gO1)*UP(gI1,1)*ZA(gI2,1) + KroneckerDelta(1,gO1)*(g2*
      UP(gI1,0)*ZA(gI2,0) - Conj(Lambdax)*UP(gI1,1)*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(Conj(ZN(gI2,4))*KroneckerDelta(1,gO2)*
      Lambdax*ZP(gI1,0)) - 0.1*(10*g2*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2) + (
      5.477225575051661*g1*Conj(ZN(gI2,0)) + 7.0710678118654755*g2*Conj(ZN(gI2,1))
      - 4.47213595499958*gN*Conj(ZN(gI2,5)))*KroneckerDelta(1,gO2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*ZN(gI2,2)*ZP(gI1
      ,0)) + 0.1*KroneckerDelta(1,gO1)*(5.477225575051661*g1*ZN(gI2,0)*ZP(gI1,0) +
      7.0710678118654755*g2*ZN(gI2,1)*ZP(gI1,0) + 6.708203932499369*gN*ZN(gI2,5)*
      ZP(gI1,0) - 10*Conj(Lambdax)*ZN(gI2,4)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaPconjSHp0PL(int gO2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*UHp0(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaPconjSHp0PR(int gO1, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*UHp0(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChahhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*
      KroneckerDelta(1,gO2)*ZH(gI1,1) + Conj(UM(gI2,1))*(g2*KroneckerDelta(0,gO2)*
      ZH(gI1,0) + KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChahhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*(g2*KroneckerDelta(0,
      gO1)*UP(gI2,1)*ZH(gI1,1) + KroneckerDelta(1,gO1)*(g2*UP(gI2,0)*ZH(gI1,0) +
      Conj(Lambdax)*UP(gI2,1)*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFeconjSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(
      ZEL(gI2,j1))*ZV(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFeconjSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(Ye(j1
      ,j1))*ZER(gI2,j1)*ZV(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*Conj(ZUR(gI1,j1))*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*SUM(j1,0,2,Conj(
      ZD(gI2,j1))*ZUL(gI1,j1))) + KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(Yd(j1,j1))
      *Conj(ZD(gI2,3 + j1))*ZUL(gI1,j1));

   return result;
}

double CLASSNAME::CpbarUChabarFvSePL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-(g2*Conj(ZE(gI2,gI1))*
      KroneckerDelta(0,gO1)),0) + IF(gI1 < 3,Conj(Ye(gI1,gI1))*Conj(ZE(gI2,3 + gI1
      ))*KroneckerDelta(1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaIconjSHI0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*SUM(j1,0,1,Conj(
      ZMI(gI2,j1))*UHI0(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaIconjSHI0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*SUM(j1,0,1,UHI0(
      gI1,2 + j1)*ZPI(gI2,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiISHIpPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*SUM(j1,0,1,Conj(
      UHIp(gI1,2 + j1))*Conj(ZNI(gI2,2 + j1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiISHIpPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*SUM(j1,0,1,Conj(
      UHIp(gI1,j1))*ZNI(gI2,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFdconjSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(
      ZDL(gI2,j1))*ZU(gI1,j1))) + KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1
      ))*Yu(j1,j1)*ZU(gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFdconjSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(Yd(j1
      ,j1))*ZDR(gI2,j1)*ZU(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*UP(
      gI2,0) + 0.1*KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*
      g2*Sin(ThetaW()))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*KroneckerDelta(0,gO1)*
      Sin(ThetaW()) + 0.1*Conj(UM(gI2,1))*KroneckerDelta(1,gO1)*(3.872983346207417
      *g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = g2*Cos(ThetaW())*Cos(ThetaWp())*
      KroneckerDelta(0,gO2)*UP(gI2,0) + 0.1*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW
      ())*Cos(ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) -
      3.1622776601683795*gN*Sin(ThetaWp()))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*Cos(
      ThetaWp())*KroneckerDelta(0,gO1) + 0.05*Conj(UM(gI2,1))*KroneckerDelta(1,gO1
      )*(10*g2*Cos(ThetaW())*Cos(ThetaWp()) - 7.745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 9.486832980505138*gN*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZpPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*
      Sin(ThetaWp())*UP(gI2,0)) - 0.1*KroneckerDelta(1,gO2)*(3.1622776601683795*gN
      *Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*
      Sin(ThetaWp()))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Cos(ThetaW())*
      KroneckerDelta(0,gO1)*Sin(ThetaWp())) + 0.05*Conj(UM(gI2,1))*KroneckerDelta(
      1,gO1)*(9.486832980505138*gN*Cos(ThetaWp()) + 2*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,1)) +
      0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZN(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(ZN(gI2,1))*KroneckerDelta(0
      ,gO1) + 1.4142135623730951*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gI2 < 3,KroneckerDelta(gI2,gO2)*Ye(gI2,
      gI2)*ZP(gI1,0),0);

   return result;
}

double CLASSNAME::CpbarUFeFvHpmPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChaSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UM(gI2,1))*Conj(ZV(gI1,gO2)
      )*Ye(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChaSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZV(gI1,gO1))*UP(gI2,0)
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZEL(gI1,gO2))*Ye(gO2,gO2)*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(Ye(gO1,gO1))*ZA(gI2,0)*ZER(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZEL(gI2
      ,gO2))*Ye(gO2,gO2)*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(Ye(gO1,
      gO1))*ZER(gI2,gO1)*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-1.0954451150103321*g1*Conj(ZE(
      gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,-0.22360679774997896*gN*Conj(
      ZE(gI1,3 + gO2))*Conj(ZN(gI2,5)),0) + IF(gO2 < 3,-(Conj(ZE(gI1,gO2))*Conj(ZN
      (gI2,2))*Ye(gO2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZE(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZE(gI1,gO1))*
      ZN(gI2,1),0) + IF(gO1 < 3,-(Conj(Ye(gO1,gO1))*Conj(ZE(gI1,3 + gO1))*ZN(gI2,2
      )),0) + IF(gO1 < 3,-0.4472135954999579*gN*Conj(ZE(gI1,gO1))*ZN(gI2,5),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.7745966692414834*g1*Cos(ThetaW
      ())*ZER(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.3872983346207417*g1*Conj(ZEL(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(ZEL(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZER(gI2,gO2),0) + IF(gI2 < 3,0.15811388300841897*gN
      *Sin(ThetaWp())*ZER(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZEL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(ZEL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-0.31622776601683794*gN*
      Conj(ZEL(gI2,gO1))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZpPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.15811388300841897*gN*Cos(
      ThetaWp())*ZER(gI2,gO2),0) + IF(gI2 < 3,0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp())*ZER(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.31622776601683794*gN*Conj(ZEL
      (gI2,gO1))*Cos(ThetaWp()),0) + IF(gI2 < 3,-0.5*g2*Conj(ZEL(gI2,gO1))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 3,0.3872983346207417*g1*Conj(ZEL(gI2,
      gO1))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

double CLASSNAME::CpbarUFeFvVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpbarUFeFvVWmPL(int gO1, int gI2) const
{
   
   const double result = IF(gI2 < 3,-0.7071067811865475*g2*KroneckerDelta(gI2,gO1)
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(ZUL(gI2,gO2))*Yd(gO2,gO2)*
      ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,Conj(Yu(gO1,gO1))*ZP(gI1,1)*ZUR(
      gI2,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZDL(gI1,gO2))*Yd(gO2,gO2)*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(Yd(gO1,gO1))*ZA(gI2,0)*ZDR(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZDL(gI2
      ,gO2))*Yd(gO2,gO2)*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(Yd(gO1,
      gO1))*ZDR(gI2,gO1)*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChaSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UM(gI2,1))*Conj(ZU(gI1,gO2)
      )*Yd(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChaSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZU(gI1,gO1))*UP(gI2,0)
      ),0) + IF(gO1 < 3,Conj(Yu(gO1,gO1))*Conj(ZU(gI1,3 + gO1))*UP(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.3651483716701107*g1*Conj(ZD(
      gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,-0.4472135954999579*gN*Conj(ZD
      (gI1,3 + gO2))*Conj(ZN(gI2,5)),0) + IF(gO2 < 3,-(Conj(ZD(gI1,gO2))*Conj(ZN(
      gI2,2))*Yd(gO2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZD(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZD(gI1,gO1))*
      ZN(gI2,1),0) + IF(gO1 < 3,-(Conj(Yd(gO1,gO1))*Conj(ZD(gI1,3 + gO1))*ZN(gI2,2
      )),0) + IF(gO1 < 3,-0.22360679774997896*gN*Conj(ZD(gI1,gO1))*ZN(gI2,5),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdGluSdPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,1.4142135623730951*g3*PhaseGlu*
      Conj(ZD(gI1,3 + gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdGluSdPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      PhaseGlu)*Conj(ZD(gI1,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*ZDR(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(ZDL(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(ThetaW
      ())*ZDR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(ZDL
      (gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(ZDL(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZDR(gI2,gO2),0) + IF(gI2 < 3,0.31622776601683794*gN
      *Sin(ThetaWp())*ZDR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZDL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZDL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-0.15811388300841897*gN*
      Conj(ZDL(gI2,gO1))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZpPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.31622776601683794*gN*Cos(
      ThetaWp())*ZDR(gI2,gO2),0) + IF(gI2 < 3,0.2581988897471611*g1*Sin(ThetaW())*
      Sin(ThetaWp())*ZDR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.15811388300841897*gN*Conj(ZDL
      (gI2,gO1))*Cos(ThetaWp()),0) + IF(gI2 < 3,-0.5*g2*Conj(ZDL(gI2,gO1))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 3,-0.12909944487358055*g1*Conj(ZDL(
      gI2,gO1))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

double CLASSNAME::CpbarUFdFuVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZUL(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(ZDL(gI2,gO2))*Yu(gO2,gO2)*
      ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,Conj(Yd(gO1,gO1))*ZDR(gI2,gO1)*
      ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarUFuSdPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UP(gI1,1))*Conj(ZD(gI2,gO2)
      )*Yu(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarUFuSdPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZD(gI2,gO1))*UM(gI1,0)
      ),0) + IF(gO1 < 3,Conj(Yd(gO1,gO1))*Conj(ZD(gI2,3 + gO1))*UM(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZUL(gI1,gO2))*Yu(gO2,gO2)*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(Yu(gO1,gO1))*ZA(gI2,1)*ZUR(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZUL(gI2
      ,gO2))*Yu(gO2,gO2)*ZH(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(Yu(gO1,
      gO1))*ZH(gI1,1)*ZUR(gI2,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7302967433402214*g1*Conj(ZN(
      gI2,0))*Conj(ZU(gI1,3 + gO2)),0) + IF(gO2 < 3,-0.22360679774997896*gN*Conj(
      ZN(gI2,5))*Conj(ZU(gI1,3 + gO2)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,3))*Conj(ZU(
      gI1,gO2))*Yu(gO2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZU(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZU(gI1,gO1))
      *ZN(gI2,1),0) + IF(gO1 < 3,-(Conj(Yu(gO1,gO1))*Conj(ZU(gI1,3 + gO1))*ZN(gI2,
      3)),0) + IF(gO1 < 3,-0.22360679774997896*gN*Conj(ZU(gI1,gO1))*ZN(gI2,5),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuGluSuPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,1.4142135623730951*g3*PhaseGlu*
      Conj(ZU(gI1,3 + gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuGluSuPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      PhaseGlu)*Conj(ZU(gI1,gO1)),0);

   return result;
}

double CLASSNAME::CpbarUFuFdconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZDL(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*ZUR(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(ZUL(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5163977794943222*g1*Cos(
      ThetaW())*ZUR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(ZUL
      (gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,-0.5*g2*Conj(ZUL(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZUR(gI2,gO2),0) + IF(gI2 < 3,0.15811388300841897*gN
      *Sin(ThetaWp())*ZUR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(ZUL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZUL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-0.15811388300841897*gN*
      Conj(ZUL(gI2,gO1))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZpPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.15811388300841897*gN*Cos(
      ThetaWp())*ZUR(gI2,gO2),0) + IF(gI2 < 3,-0.5163977794943222*g1*Sin(ThetaW())
      *Sin(ThetaWp())*ZUR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.15811388300841897*gN*Conj(ZUL
      (gI2,gO1))*Cos(ThetaWp()),0) + IF(gI2 < 3,0.5*g2*Conj(ZUL(gI2,gO1))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 3,-0.12909944487358055*g1*Conj(ZUL(
      gI2,gO1))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZDXL(gI1,gO2))*ZA(gI2,2)*Kappa(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(Kappa(gO1,gO1))*ZA(gI2,2)*ZDXR(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZDXL(
      gI2,gO2))*ZH(gI1,2)*Kappa(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(Kappa(
      gO1,gO1))*ZDXR(gI2,gO1)*ZH(gI1,2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXChiSDXPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.3651483716701107*g1*Conj(ZDX(
      gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,0.6708203932499369*gN*Conj(ZDX
      (gI1,3 + gO2))*Conj(ZN(gI2,5)),0) + IF(gO2 < 3,-(Conj(ZDX(gI1,gO2))*Conj(ZN(
      gI2,4))*Kappa(gO2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXChiSDXPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.3651483716701107*g1*Conj(ZDX(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-(Conj(ZDX(gI1,3 + gO1))*Conj(Kappa(gO1,
      gO1))*ZN(gI2,4)),0) + IF(gO1 < 3,0.4472135954999579*gN*Conj(ZDX(gI1,gO1))*ZN
      (gI2,5),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXGluSDXPL(int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,1.4142135623730951*g3*PhaseGlu*
      Conj(ZDX(gI1,3 + gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXGluSDXPR(int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      PhaseGlu)*Conj(ZDX(gI1,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVGPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*ZDXR(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVGPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(ZDXL(gI2,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(ThetaW
      ())*ZDXR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Conj(ZDXL(
      gI2,gO1))*Cos(ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZDXR(gI2,gO2),0) + IF(gI2 < 3,-0.4743416490252569*
      gN*Sin(ThetaWp())*ZDXR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Conj(ZDXL
      (gI2,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,0.31622776601683794*
      gN*Conj(ZDXL(gI2,gO1))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVZpPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.4743416490252569*gN*Cos(
      ThetaWp())*ZDXR(gI2,gO2),0) + IF(gI2 < 3,0.2581988897471611*g1*Sin(ThetaW())
      *Sin(ThetaWp())*ZDXR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFDXFDXVZpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.31622776601683794*gN*Conj(ZDXL
      (gI2,gO1))*Cos(ThetaWp()),0) + IF(gI2 < 3,0.2581988897471611*g1*Conj(ZDXL(
      gI2,gO1))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZMI(gI1,gO2))*ZA(gI2,2)*Lambda12(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,std::complex<double>(0.,
      0.7071067811865475)*Conj(Lambda12(gO1,gO1))*ZA(gI2,2)*ZPI(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.7071067811865475*Conj(ZMI(gI2
      ,gO2))*ZH(gI1,2)*Lambda12(gO2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.7071067811865475*Conj(
      Lambda12(gO1,gO1))*ZH(gI1,2)*ZPI(gI2,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaSHI0PL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-(g2*Conj(UHI0(gI1,2 + gO2))*
      Conj(UM(gI2,0))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaSHI0PR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-(g2*Conj(UHI0(gI1,gO1))*UP(gI2,
      0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChiSHIpPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.5477225575051661*g1*Conj(UHIp
      (gI1,2 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 2,-0.7071067811865475*g2*Conj(
      UHIp(gI1,2 + gO2))*Conj(ZN(gI2,1)),0) + IF(gO2 < 2,0.4472135954999579*gN*
      Conj(UHIp(gI1,2 + gO2))*Conj(ZN(gI2,5)),0) + IF(gO2 < 2,-(Conj(UHIp(gI1,gO2)
      )*Conj(ZN(gI2,4))*Lambda12(gO2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChiSHIpPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.5477225575051661*g1*Conj(UHIp(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 2,0.7071067811865475*g2*Conj(UHIp(gI1,gO1)
      )*ZN(gI2,1),0) + IF(gO1 < 2,-(Conj(UHIp(gI1,2 + gO1))*Conj(Lambda12(gO1,gO1)
      )*ZN(gI2,4)),0) + IF(gO1 < 2,0.6708203932499369*gN*Conj(UHIp(gI1,gO1))*ZN(
      gI2,5),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.3872983346207417*g1*Cos(ThetaW
      ())*ZPI(gI2,gO2),0) + IF(gI2 < 2,0.5*g2*Sin(ThetaW())*ZPI(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.3872983346207417*g1*Conj(ZMI(
      gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 2,0.5*g2*Conj(ZMI(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.5*g2*Cos(ThetaW())*Cos(ThetaWp
      ())*ZPI(gI2,gO2),0) + IF(gI2 < 2,-0.3872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW())*ZPI(gI2,gO2),0) + IF(gI2 < 2,-0.31622776601683794*gN*Sin(ThetaWp()
      )*ZPI(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.5*g2*Conj(ZMI(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 2,-0.3872983346207417*g1*Conj(ZMI(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 2,0.4743416490252569*gN*
      Conj(ZMI(gI2,gO1))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIVZpPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.31622776601683794*gN*Cos(
      ThetaWp())*ZPI(gI2,gO2),0) + IF(gI2 < 2,-0.5*g2*Cos(ThetaW())*Sin(ThetaWp())
      *ZPI(gI2,gO2),0) + IF(gI2 < 2,0.3872983346207417*g1*Sin(ThetaW())*Sin(
      ThetaWp())*ZPI(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChaIVZpPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.4743416490252569*gN*Conj(ZMI(
      gI2,gO1))*Cos(ThetaWp()),0) + IF(gI2 < 2,-0.5*g2*Conj(ZMI(gI2,gO1))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gI2 < 2,0.3872983346207417*g1*Conj(ZMI(gI2,
      gO1))*Sin(ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChiIVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.7071067811865475*g2*ZNI(gI2,2
      + gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaIChiIVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.7071067811865475*g2*Conj(ZNI(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiISHIpPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*SUM(j1,0,1,Conj(UHIp(
      gI2,2 + j1))*KroneckerDelta(gO2,2 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiISHIpPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-(g2*Conj(UHIp(gI2,gO1))*UM(gI1,
      0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIUChiIVWmPL(int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.7071067811865475*g2*ZMI(gI1,
      gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaIUChiIVWmPR(int gI1, int gO1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,1,Conj(ZPI(
      gI1,j1))*KroneckerDelta(gO1,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIhhPL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.7071067811865475*Conj(ZNI(gI2,
      2 + gO2))*ZH(gI1,2)*Lambda12(gO2,gO2),0) + 0.7071067811865475*SUM(j1,0,1,
      Conj(ZNI(gI2,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZH(gI1,2);

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIhhPR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.7071067811865475*Conj(Lambda12
      (gO1,gO1))*ZH(gI1,2)*ZNI(gI2,2 + gO1),0) + 0.7071067811865475*SUM(j1,0,1,
      Conj(Lambda12(j1,j1))*KroneckerDelta(gO1,2 + j1)*ZNI(gI2,j1))*ZH(gI1,2);

   return result;
}

std::complex<double> CLASSNAME::CpUChiIChaconjSHIpPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-(g2*Conj(UM(gI2,0))*UHIp(gI1,
      gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiIChaconjSHIpPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,1,KroneckerDelta(gO1,2 + j1)*
      UHIp(gI1,2 + j1))*UP(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIAhPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZNI(gI1,2 + gO2))*ZA(gI2,2)*Lambda12(gO2,gO2),0) +
      std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,1,Conj(ZNI(gI1,j1))*
      KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIAhPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,std::complex<double>(0.,-
      0.7071067811865475)*Conj(Lambda12(gO1,gO1))*ZA(gI2,2)*ZNI(gI1,2 + gO1),0) -
      std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,1,Conj(Lambda12(j1,j1))
      *KroneckerDelta(gO1,2 + j1)*ZNI(gI1,j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiIconjSHI0PL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.5477225575051661*g1*Conj(ZN(
      gI2,0))*UHI0(gI1,gO2),0) + IF(gO2 < 2,-0.7071067811865475*g2*Conj(ZN(gI2,1))
      *UHI0(gI1,gO2),0) + IF(gO2 < 2,0.6708203932499369*gN*Conj(ZN(gI2,5))*UHI0(
      gI1,gO2),0) + IF(gO2 < 2,Conj(ZN(gI2,4))*UHI0(gI1,2 + gO2)*Lambda12(gO2,gO2)
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiIconjSHI0PR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,1,Conj(Lambda12(j1,j1))*
      KroneckerDelta(gO1,2 + j1)*UHI0(gI1,j1))*ZN(gI2,4) + SUM(j1,0,1,
      KroneckerDelta(gO1,2 + j1)*UHI0(gI1,2 + j1))*(-0.5477225575051661*g1*ZN(gI2,
      0) + 0.7071067811865475*g2*ZN(gI2,1) + 0.4472135954999579*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiISHI0PL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = -0.5477225575051661*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1)) +
      0.7071067811865475*g2*Conj(ZN(gI2,1))*SUM(j1,0,1,Conj(UHI0(gI1,2 + j1))*
      KroneckerDelta(gO2,2 + j1)) + 0.4472135954999579*gN*Conj(ZN(gI2,5))*SUM(j1,0
      ,1,Conj(UHI0(gI1,2 + j1))*KroneckerDelta(gO2,2 + j1)) + Conj(ZN(gI2,4))*SUM(
      j1,0,1,Conj(UHI0(gI1,j1))*KroneckerDelta(gO2,2 + j1)*Lambda12(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiISHI0PR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.5477225575051661*g1*Conj(UHI0(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 2,-0.7071067811865475*g2*Conj(UHI0(gI1,gO1
      ))*ZN(gI2,1),0) + IF(gO1 < 2,Conj(UHI0(gI1,2 + gO1))*Conj(Lambda12(gO1,gO1))
      *ZN(gI2,4),0) + IF(gO1 < 2,0.6708203932499369*gN*Conj(UHI0(gI1,gO1))*ZN(gI2,
      5),0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiIChaIconjVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,1,
      KroneckerDelta(gO2,2 + j1)*ZPI(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiIChaIconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.7071067811865475*g2*Conj(ZMI(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIVZPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-0.5*g2*Cos(ThetaW())*Cos(
      ThetaWp())*ZNI(gI2,gO2),0) + IF(gO2 < 2,-0.3872983346207417*g1*Cos(ThetaWp()
      )*Sin(ThetaW())*ZNI(gI2,gO2),0) + IF(gO2 < 2,0.4743416490252569*gN*Sin(
      ThetaWp())*ZNI(gI2,gO2),0) + 0.5*g2*Cos(ThetaW())*Cos(ThetaWp())*SUM(j1,0,1,
      KroneckerDelta(gO2,2 + j1)*ZNI(gI2,2 + j1)) + 0.3872983346207417*g1*Cos(
      ThetaWp())*Sin(ThetaW())*SUM(j1,0,1,KroneckerDelta(gO2,2 + j1)*ZNI(gI2,2 +
      j1)) + 0.31622776601683794*gN*Sin(ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO2,2
       + j1)*ZNI(gI2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIVZPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,0.5*g2*Conj(ZNI(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO1 < 2,0.3872983346207417*g1*Conj(ZNI(gI2,
      gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO1 < 2,-0.4743416490252569*gN*
      Conj(ZNI(gI2,gO1))*Sin(ThetaWp()),0) - 0.5*g2*Cos(ThetaW())*Cos(ThetaWp())*
      SUM(j1,0,1,Conj(ZNI(gI2,2 + j1))*KroneckerDelta(gO1,2 + j1)) -
      0.3872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*SUM(j1,0,1,Conj(ZNI(gI2,2
       + j1))*KroneckerDelta(gO1,2 + j1)) - 0.31622776601683794*gN*Sin(ThetaWp())*
      SUM(j1,0,1,Conj(ZNI(gI2,2 + j1))*KroneckerDelta(gO1,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIVZpPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 2,0.4743416490252569*gN*Cos(
      ThetaWp())*ZNI(gI2,gO2),0) + IF(gO2 < 2,0.5*g2*Cos(ThetaW())*Sin(ThetaWp())*
      ZNI(gI2,gO2),0) + IF(gO2 < 2,0.3872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp
      ())*ZNI(gI2,gO2),0) + 0.31622776601683794*gN*Cos(ThetaWp())*SUM(j1,0,1,
      KroneckerDelta(gO2,2 + j1)*ZNI(gI2,2 + j1)) - 0.5*g2*Cos(ThetaW())*Sin(
      ThetaWp())*SUM(j1,0,1,KroneckerDelta(gO2,2 + j1)*ZNI(gI2,2 + j1)) -
      0.3872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,1,KroneckerDelta
      (gO2,2 + j1)*ZNI(gI2,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiIUChiIVZpPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-0.4743416490252569*gN*Conj(ZNI(
      gI2,gO1))*Cos(ThetaWp()),0) + IF(gO1 < 2,-0.5*g2*Conj(ZNI(gI2,gO1))*Cos(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 2,-0.3872983346207417*g1*Conj(ZNI(gI2
      ,gO1))*Sin(ThetaW())*Sin(ThetaWp()),0) - 0.31622776601683794*gN*Cos(ThetaWp(
      ))*SUM(j1,0,1,Conj(ZNI(gI2,2 + j1))*KroneckerDelta(gO1,2 + j1)) + 0.5*g2*Cos
      (ThetaW())*Sin(ThetaWp())*SUM(j1,0,1,Conj(ZNI(gI2,2 + j1))*KroneckerDelta(
      gO1,2 + j1)) + 0.3872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,1
      ,Conj(ZNI(gI2,2 + j1))*KroneckerDelta(gO1,2 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUFSIconjSSI0PL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 2,-1.118033988749895*gN*Conj(ZN(
      gI2,5))*ZSSI(gI1,gO2),0);

   return result;
}

double CLASSNAME::CpChiUFSIconjSSI0PR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

double CLASSNAME::CpChiUFSISSI0PL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpChiUFSISSI0PR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 2,-1.118033988749895*gN*Conj(ZSSI(
      gI1,gO1))*ZN(gI2,5),0);

   return result;
}

std::complex<double> CLASSNAME::CpFSIUFSIVZPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.7905694150420949*gN*Sin(
      ThetaWp())*ZFSI(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpFSIUFSIVZPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.7905694150420949*gN*Conj(ZFSI(
      gI2,gO1))*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpFSIUFSIVZpPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 2,-0.7905694150420949*gN*Cos(
      ThetaWp())*ZFSI(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpFSIUFSIVZpPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gI2 < 2,0.7905694150420949*gN*Conj(ZFSI(
      gI2,gO1))*Cos(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiPSHppPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UHpp(gI2,1))*Conj(UP(gI1,0))*
      KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiPSHppPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(UHpp(gI2,0))*KroneckerDelta(0,gO1
      )*UM(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpUChiPChaconjSHppPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*KroneckerDelta(0,gO2)*
      UHpp(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpUChiPChaconjSHppPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(1,gO1)*UHpp(gI1,1)*UP(
      gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiPconjSHp0PL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.1*(5.477225575051661*g1*Conj(ZN(gI2,0)) -
      7.0710678118654755*g2*Conj(ZN(gI2,1)) - 4.47213595499958*gN*Conj(ZN(gI2,5)))
      *KroneckerDelta(0,gO2)*UHp0(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiPconjSHp0PR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(1,gO1)*UHp0(gI1,1)*(-
      5.477225575051661*g1*ZN(gI2,0) + 7.0710678118654755*g2*ZN(gI2,1) +
      4.47213595499958*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiPSHp0PL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.1*Conj(UHp0(gI1,1))*(-5.477225575051661*
      g1*Conj(ZN(gI2,0)) + 7.0710678118654755*g2*Conj(ZN(gI2,1)) +
      4.47213595499958*gN*Conj(ZN(gI2,5)))*KroneckerDelta(1,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiPSHp0PR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = 0.1*Conj(UHp0(gI1,0))*KroneckerDelta(0,gO1)
      *(5.477225575051661*g1*ZN(gI2,0) - 7.0710678118654755*g2*ZN(gI2,1) -
      4.47213595499958*gN*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiPUChiPVZPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.1*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*
      Sin(ThetaWp()))*(KroneckerDelta(0,gO2)*ZNp(gI2,0) - KroneckerDelta(1,gO2)*
      ZNp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPUChiPVZPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.1*(Conj(ZNp(gI2,0))*KroneckerDelta(0,gO1)
      - Conj(ZNp(gI2,1))*KroneckerDelta(1,gO1))*(5*g2*Cos(ThetaW())*Cos(ThetaWp())
      + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3.1622776601683795*gN*
      Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpChiPUChiPVZpPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.1*(3.1622776601683795*gN*Cos(ThetaWp())
      - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*
      (KroneckerDelta(0,gO2)*ZNp(gI2,0) - KroneckerDelta(1,gO2)*ZNp(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpChiPUChiPVZpPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.1*(Conj(ZNp(gI2,0))*KroneckerDelta(0,gO1)
      - Conj(ZNp(gI2,1))*KroneckerDelta(1,gO1))*(3.1622776601683795*gN*Cos(ThetaWp
      ()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp(
      )));

   return result;
}

double CLASSNAME::CpbarChaPUChiPVWmPL(int gO2) const
{
   
   const double result = -0.7071067811865475*g2*KroneckerDelta(0,gO2);

   return result;
}

double CLASSNAME::CpbarChaPUChiPVWmPR(int gO1) const
{
   
   const double result = 0.7071067811865475*g2*KroneckerDelta(1,gO1);

   return result;
}

double CLASSNAME::CpUChiPChaPconjVWmPR(int gO2) const
{
   
   const double result = 0.7071067811865475*g2*KroneckerDelta(1,gO2);

   return result;
}

double CLASSNAME::CpUChiPChaPconjVWmPL(int gO1) const
{
   
   const double result = -0.7071067811865475*g2*KroneckerDelta(0,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdGluSdPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*PhaseGlu*SUM(j1,0,2,
      Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdGluSdPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -1.4142135623730951*g3*Conj(PhaseGlu)*SUM(
      j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXGluSDXPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*PhaseGlu*SUM(j1,0,2,
      Conj(ZDX(gI2,3 + j1))*Conj(ZDXR(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFDXGluSDXPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -1.4142135623730951*g3*Conj(PhaseGlu)*SUM(
      j1,0,2,Conj(ZDX(gI2,j1))*ZDXL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuGluSuPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*PhaseGlu*SUM(j1,0,2,
      Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuGluSuPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -1.4142135623730951*g3*Conj(PhaseGlu)*SUM(
      j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFdconjSdPL(int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.4142135623730951*g3*PhaseGlu*SUM(j1,0,2,
      Conj(ZDL(gI2,j1))*ZD(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFdconjSdPR(int gI2, int gI1) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*Conj(PhaseGlu)*SUM(j1
      ,0,2,ZD(gI1,3 + j1)*ZDR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFDXconjSDXPL(int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.4142135623730951*g3*PhaseGlu*SUM(j1,0,2,
      Conj(ZDXL(gI2,j1))*ZDX(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFDXconjSDXPR(int gI2, int gI1) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*Conj(PhaseGlu)*SUM(j1
      ,0,2,ZDX(gI1,3 + j1)*ZDXR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFuconjSuPL(int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.4142135623730951*g3*PhaseGlu*SUM(j1,0,2,
      Conj(ZUL(gI2,j1))*ZU(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpGluFuconjSuPR(int gI2, int gI1) const
{
   
   const std::complex<double> result = 1.4142135623730951*g3*Conj(PhaseGlu)*SUM(j1
      ,0,2,ZU(gI1,3 + j1)*ZUR(gI2,j1));

   return result;
}

double CLASSNAME::CpbarFvFeconjHpmPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,Conj(Ye(gO1,gO1))*ZER(gI2,gO1)*
      ZP(gI1,0),0);

   return result;
}

double CLASSNAME::CpbarChabarFvSePL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarFvSePR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZE(gI2,gO1))*UM(gI1,0)
      ),0) + IF(gO1 < 3,Conj(Ye(gO1,gO1))*Conj(ZE(gI2,3 + gO1))*UM(gI1,1),0);

   return result;
}

double CLASSNAME::CpbarFvChiSvPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvChiSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZV(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZV(gI1,gO1))
      *ZN(gI2,1),0) + IF(gO1 < 3,-0.4472135954999579*gN*Conj(ZV(gI1,gO1))*ZN(gI2,5
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPChaSHp0PL(int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*Conj(UHp0(gI1,1))*Conj(UM(gI2,0)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPChaSHp0PR(int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*Conj(UHp0(gI1,0))*UP(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPChiSHppPL(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.1*Conj(UHpp(gI1,1))*(5.477225575051661*
      g1*Conj(ZN(gI2,0)) + 7.0710678118654755*g2*Conj(ZN(gI2,1)) -
      4.47213595499958*gN*Conj(ZN(gI2,5)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPChiSHppPR(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*Conj(UHpp(gI1,0))*(5.477225575051661*g1
      *ZN(gI2,0) + 7.0710678118654755*g2*ZN(gI2,1) - 4.47213595499958*gN*ZN(gI2,5)
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPChiPVWmPR(int gI2) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*ZNp(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaPChiPVWmPL(int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*Conj(ZNp(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gI2 < 3,Conj(ZER(gO2,gI2))*Ye(gI2,gI2)*
      ZP(gI1,0),0);

   return result;
}

double CLASSNAME::CpbarFeFvHpmPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChaSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = Conj(UM(gI2,1))*SUM(j1,0,2,Conj(ZER(gO2,j1)
      )*Conj(ZV(gI1,j1))*Ye(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChaSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZEL(gO1,j1
      ))*UP(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZEL(gI1,j1))*Conj(ZER(gO2,j1))*Ye(j1,j1))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j1,0,2,Conj(Ye(j1,j1))*ZEL(gO1,j1)*ZER(gI1,j1))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZEL(gI2
      ,j1))*Conj(ZER(gO2,j1))*Ye(j1,j1))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(Ye(j1,
      j1))*ZEL(gO1,j1)*ZER(gI2,j1))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.0954451150103321*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1))) - 0.22360679774997896*gN*Conj
      (ZN(gI2,5))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1))) - Conj(ZN(gI2
      ,2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*Conj(ZER(gO2,j1))*Ye(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j1,0,2,Conj(Ye(j1,j1))*Conj(ZE(gI1,3
      + j1))*ZEL(gO1,j1))*ZN(gI2,2)) + SUM(j1,0,2,Conj(ZE(gI1,j1))*ZEL(gO1,j1))*(
      0.5477225575051661*g1*ZN(gI2,0) + 0.7071067811865475*g2*ZN(gI2,1) -
      0.4472135954999579*gN*ZN(gI2,5));

   return result;
}

double CLASSNAME::CpbarFeFvVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*ZEL(gO1,
      gI2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(ZDR(gO2,j1))*Conj(ZUL(gI2,
      j1))*Yd(j1,j1))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(Yu(j1,j1))*ZDL(gO1,j1)*ZUR(
      gI2,j1))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZDL(gI1,j1))*Conj(ZDR(gO2,j1))*Yd(j1,j1))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j1,0,2,Conj(Yd(j1,j1))*ZDL(gO1,j1)*ZDR(gI1,j1))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZDL(gI2
      ,j1))*Conj(ZDR(gO2,j1))*Yd(j1,j1))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(Yd(j1,
      j1))*ZDL(gO1,j1)*ZDR(gI2,j1))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChaSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = Conj(UM(gI2,1))*SUM(j1,0,2,Conj(ZDR(gO2,j1)
      )*Conj(ZU(gI1,j1))*Yd(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChaSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZDL(gO1,j1
      ))*UP(gI2,0)) + SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI1,3 + j1))*ZDL(gO1,j1))
      *UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.3651483716701107*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1))) - 0.4472135954999579*gN*Conj(
      ZN(gI2,5))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1))) - Conj(ZN(gI2,
      2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*Conj(ZDR(gO2,j1))*Yd(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI1,3
      + j1))*ZDL(gO1,j1))*ZN(gI2,2)) - 0.03333333333333333*SUM(j1,0,2,Conj(ZD(gI1,
      j1))*ZDL(gO1,j1))*(5.477225575051661*g1*ZN(gI2,0) - 21.213203435596427*g2*ZN
      (gI2,1) + 6.708203932499369*gN*ZN(gI2,5));

   return result;
}

double CLASSNAME::CpbarFdFuVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZUL(
      gI2,j1))*ZDL(gO1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(ZDL(gI2,j1))*Conj(ZUR(gO2,
      j1))*Yu(j1,j1))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(Yd(j1,j1))*ZDR(gI2,j1)*ZUL(
      gO1,j1))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarFuSdPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = Conj(UP(gI1,1))*SUM(j1,0,2,Conj(ZD(gI2,j1))
      *Conj(ZUR(gO2,j1))*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarFuSdPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZUL(gO1,j1
      ))*UM(gI1,0)) + SUM(j1,0,2,Conj(Yd(j1,j1))*Conj(ZD(gI2,3 + j1))*ZUL(gO1,j1))
      *UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j1,0,2,Conj(ZUL(gI1,j1))*Conj(ZUR(gO2,j1))*Yu(j1,j1))*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j1,0,2,Conj(Yu(j1,j1))*ZUL(gO1,j1)*ZUR(gI1,j1))*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(ZUL(gI2
      ,j1))*Conj(ZUR(gO2,j1))*Yu(j1,j1))*ZH(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j1,0,2,Conj(Yu(j1,
      j1))*ZUL(gO1,j1)*ZUR(gI2,j1))*ZH(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7302967433402214*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,2,Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1))) - 0.22360679774997896*gN*Conj
      (ZN(gI2,5))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1))) - Conj(ZN(gI2
      ,3))*SUM(j1,0,2,Conj(ZU(gI1,j1))*Conj(ZUR(gO2,j1))*Yu(j1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(SUM(j1,0,2,Conj(Yu(j1,j1))*Conj(ZU(gI1,3
      + j1))*ZUL(gO1,j1))*ZN(gI2,3)) - 0.03333333333333333*SUM(j1,0,2,Conj(ZU(gI1,
      j1))*ZUL(gO1,j1))*(5.477225575051661*g1*ZN(gI2,0) + 21.213203435596427*g2*ZN
      (gI2,1) + 6.708203932499369*gN*ZN(gI2,5));

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSdconjUSdVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSdconjUSdVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSdconjHpmconjUSd(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUSdSHp0conjUSdconjSHp0(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUSdSHppconjUSdconjSHpp(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSdSSI0conjUSdconjSSI0(gO1,gI1,gO2
      ,gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSdconjUSd(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSdconjUSd(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSdSvconjUSdconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,1,(Conj(CpChaFuconjUSdPL(gI2,gI1,gO2))*
      CpChaFuconjUSdPL(gI2,gI1,gO1) + Conj(CpChaFuconjUSdPR(gI2,gI1,gO2))*
      CpChaFuconjUSdPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MCha(
      gI2)))*(Conj(CpChaFuconjUSdPR(gI2,gI1,gO2))*CpChaFuconjUSdPL(gI2,gI1,gO1) +
      Conj(CpChaFuconjUSdPL(gI2,gI1,gO2))*CpChaFuconjUSdPR(gI2,gI1,gO1))*MCha(gI2)
      ));
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFdconjUSdPL(gI2,gI1,gO2))*
      CpChiFdconjUSdPL(gI2,gI1,gO1) + Conj(CpChiFdconjUSdPR(gI2,gI1,gO2))*
      CpChiFdconjUSdPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFdconjUSdPR(gI2,gI1,gO2))*CpChiFdconjUSdPL(gI2,gI1,gO1) +
      Conj(CpChiFdconjUSdPL(gI2,gI1,gO2))*CpChiFdconjUSdPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUSdSHI0conjUSdconjSHI0(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUSdSHIpconjUSdconjSHIp(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUSdconjUSdconjSDXSDX(gO1,gO2,gI1,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUSdSeconjUSdconjSe(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MHpm(gI2)))*Conj(
      CpHpmSuconjUSd(gI2,gI1,gO2))*CpHpmSuconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSdconjUSd(gI2,gI1,gO2))*CpAhSdconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSdconjUSd(gI2,gI1,gO2))*CphhSdconjUSd(gI2,gI1,gO1)));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpGluFdconjUSdPL(gI2,gO2))*
      CpGluFdconjUSdPL(gI2,gO1) + Conj(CpGluFdconjUSdPR(gI2,gO2))*CpGluFdconjUSdPR
      (gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFd(gI2))));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFd(gI2)
      ))*(Conj(CpGluFdconjUSdPR(gI2,gO2))*CpGluFdconjUSdPL(gI2,gO1) + Conj(
      CpGluFdconjUSdPL(gI2,gO2))*CpGluFdconjUSdPR(gI2,gO1))*MFd(gI2));
   result += 1.3333333333333333*SUM(gI2,0,5,Conj(CpSdconjUSdVG(gI2,gO2))*
      CpSdconjUSdVG(gI2,gO1)*F0(Sqr(p),Sqr(MSd(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSdconjUSdVP(gI2,gO2))*CpSdconjUSdVP(gI2,gO1)*F0(
      Sqr(p),Sqr(MSd(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSdconjUSdVZ(gI2,gO2))*CpSdconjUSdVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSd(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSdconjUSdVZp(gI2,gO2))*CpSdconjUSdVZp(gI2,gO1)*F0(
      Sqr(p),Sqr(MSd(gI2)),Sqr(MVZp)));
   result += SUM(gI2,0,5,Conj(CpSuconjUSdVWm(gI2,gO2))*CpSuconjUSdVWm(gI2,gO1)*F0(
      Sqr(p),Sqr(MSu(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Sd_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = i; k < 6; k++)
         self_energy(i, k) = self_energy_Sd_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Sv_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSvconjUSvVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSvconjUSvVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSvconjHpmconjUSv(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0USvconjSHp0conjUSv(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppUSvconjSHppconjUSv(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpSSI0USvconjSSI0conjUSv(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,2,(Conj(CpbarChaFeconjUSvPL(gI1,gI2,gO2))*
      CpbarChaFeconjUSvPL(gI1,gI2,gO1) + Conj(CpbarChaFeconjUSvPR(gI1,gI2,gO2))*
      CpbarChaFeconjUSvPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFe(gI2)))));
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFe
      (gI2)))*(Conj(CpbarChaFeconjUSvPR(gI1,gI2,gO2))*CpbarChaFeconjUSvPL(gI1,gI2,
      gO1) + Conj(CpbarChaFeconjUSvPL(gI1,gI2,gO2))*CpbarChaFeconjUSvPR(gI1,gI2,
      gO1))*MFe(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSe(gI2)))*Conj(
      CpSeconjHpmconjUSv(gI2,gI1,gO2))*CpSeconjHpmconjUSv(gI2,gI1,gO1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSvconjUSv(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSvconjUSv(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvUSvconjSvconjUSv(gI1,gO1,gI1,gO2))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSvconjUSv(gI2,gI1,gO2))*CphhSvconjUSv(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*
      CpChiFvconjUSvPL(gI2,gI1,gO1) + Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*
      CpChiFvconjUSvPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*CpChiFvconjUSvPL(gI2,gI1,gO1) +
      Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*CpChiFvconjUSvPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USvconjSHI0conjUSv(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpUSvconjSHIpconjUSv(gI1,gO1,gI1
      ,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSvconjSdconjUSv(gI1,gO1,gI1,gO2
      ));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSvconjSDXconjUSv(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSvconjSeconjUSv(gI1,gO1,gI1,gO2))
      ;
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuUSvconjSuconjUSv(gI1,gO1,gI1,gO2
      ));
   result += SUM(gI2,0,2,Conj(CpSvconjUSvVZ(gI2,gO2))*CpSvconjUSvVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,2,Conj(CpSvconjUSvVZp(gI2,gO2))*CpSvconjUSvVZp(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVZp)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSvconjVWm(gI2,gO2))*CpSeconjUSvconjVWm(gI2,
      gO1)*F0(Sqr(p),Sqr(MSe(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Sv_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = i; k < 3; k++)
         self_energy(i, k) = self_energy_Sv_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Su_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSuconjUSuVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSuconjUSuVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSuconjHpmconjUSu(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0USuconjSHp0conjUSu(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppUSuconjSHppconjUSu(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpSSI0USuconjSSI0conjUSu(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,2,(Conj(CpbarChaFdconjUSuPL(gI1,gI2,gO2))*
      CpbarChaFdconjUSuPL(gI1,gI2,gO1) + Conj(CpbarChaFdconjUSuPR(gI1,gI2,gO2))*
      CpbarChaFdconjUSuPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFd(gI2)))));
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFd
      (gI2)))*(Conj(CpbarChaFdconjUSuPR(gI1,gI2,gO2))*CpbarChaFdconjUSuPL(gI1,gI2,
      gO1) + Conj(CpbarChaFdconjUSuPL(gI1,gI2,gO2))*CpbarChaFdconjUSuPR(gI1,gI2,
      gO1))*MFd(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSd(gI2)))*Conj(
      CpSdconjHpmconjUSu(gI2,gI1,gO2))*CpSdconjHpmconjUSu(gI2,gI1,gO1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSuconjUSu(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSuconjUSu(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSuSvconjUSuconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*
      CpChiFuconjUSuPL(gI2,gI1,gO1) + Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*
      CpChiFuconjUSuPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*CpChiFuconjUSuPL(gI2,gI1,gO1) +
      Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*CpChiFuconjUSuPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USuconjSHI0conjUSu(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpUSuconjSHIpconjUSu(gI1,gO1,gI1
      ,gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSuconjSeconjUSu(gI1,gO1,gI1,gO2))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUSuconjUSuconjSDXSDX(gO1,gO2,gI1,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSuconjUSu(gI2,gI1,gO2))*CpAhSuconjUSu(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSuconjUSu(gI2,gI1,gO2))*CphhSuconjUSu(gI2,gI1,gO1)));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpGluFuconjUSuPL(gI2,gO2))*
      CpGluFuconjUSuPL(gI2,gO1) + Conj(CpGluFuconjUSuPR(gI2,gO2))*CpGluFuconjUSuPR
      (gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFu(gI2))));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFu(gI2)
      ))*(Conj(CpGluFuconjUSuPR(gI2,gO2))*CpGluFuconjUSuPL(gI2,gO1) + Conj(
      CpGluFuconjUSuPL(gI2,gO2))*CpGluFuconjUSuPR(gI2,gO1))*MFu(gI2));
   result += SUM(gI2,0,5,Conj(CpSdconjUSuconjVWm(gI2,gO2))*CpSdconjUSuconjVWm(gI2,
      gO1)*F0(Sqr(p),Sqr(MSd(gI2)),Sqr(MVWm)));
   result += 1.3333333333333333*SUM(gI2,0,5,Conj(CpSuconjUSuVG(gI2,gO2))*
      CpSuconjUSuVG(gI2,gO1)*F0(Sqr(p),Sqr(MSu(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSuconjUSuVP(gI2,gO2))*CpSuconjUSuVP(gI2,gO1)*F0(
      Sqr(p),Sqr(MSu(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSuconjUSuVZ(gI2,gO2))*CpSuconjUSuVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSu(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSuconjUSuVZp(gI2,gO2))*CpSuconjUSuVZp(gI2,gO1)*F0(
      Sqr(p),Sqr(MSu(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Su_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = i; k < 6; k++)
         self_energy(i, k) = self_energy_Su_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Se_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSeconjUSeVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSeconjUSeVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSeconjHpmconjUSe(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUSeSHp0conjUSeconjSHp0(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUSeSHppconjUSeconjSHpp(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSeSSI0conjUSeconjSSI0(gO1,gI1,gO2
      ,gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSeconjUSe(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSeconjUSe(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSeSvconjUSeconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MHpm(gI2)))*Conj(
      CpHpmSvconjUSe(gI2,gI1,gO2))*CpHpmSvconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,1,(Conj(CpChaFvconjUSePL(gI2,gI1,gO2))*
      CpChaFvconjUSePL(gI2,gI1,gO1) + Conj(CpChaFvconjUSePR(gI2,gI1,gO2))*
      CpChaFvconjUSePR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MCha(
      gI2)))*(Conj(CpChaFvconjUSePR(gI2,gI1,gO2))*CpChaFvconjUSePL(gI2,gI1,gO1) +
      Conj(CpChaFvconjUSePL(gI2,gI1,gO2))*CpChaFvconjUSePR(gI2,gI1,gO1))*MCha(gI2)
      ));
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFeconjUSePL(gI2,gI1,gO2))*
      CpChiFeconjUSePL(gI2,gI1,gO1) + Conj(CpChiFeconjUSePR(gI2,gI1,gO2))*
      CpChiFeconjUSePR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFeconjUSePR(gI2,gI1,gO2))*CpChiFeconjUSePL(gI2,gI1,gO1) +
      Conj(CpChiFeconjUSePL(gI2,gI1,gO2))*CpChiFeconjUSePR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUSeSHI0conjUSeconjSHI0(gO1,gI1,gO2
      ,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUSeSHIpconjUSeconjSHIp(gO1,gI1,gO2
      ,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSeconjSdconjUSe(gI1,gO1,gI1,gO2
      ));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSeconjSDXconjUSe(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSeconjSeconjUSe(gI1,gO1,gI1,gO2))
      ;
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSeSuconjUSeconjSu(gO1,gI1,gO2,gI1
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSeconjUSe(gI2,gI1,gO2))*CpAhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSeconjUSe(gI2,gI1,gO2))*CphhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI2,0,2,Conj(CpSvconjUSeVWm(gI2,gO2))*CpSvconjUSeVWm(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVP(gI2,gO2))*CpSeconjUSeVP(gI2,gO1)*F0(
      Sqr(p),Sqr(MSe(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVZ(gI2,gO2))*CpSeconjUSeVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSe(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVZp(gI2,gO2))*CpSeconjUSeVZp(gI2,gO1)*F0(
      Sqr(p),Sqr(MSe(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Se_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = i; k < 6; k++)
         self_energy(i, k) = self_energy_Se_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_SDX_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 2*A0(Sqr(MVZp))*CpUSDXconjUSDXVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSDXconjUSDXVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSDXconjHpmconjUSDX(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUSDXSHp0conjUSDXconjSHp0(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUSDXSHppconjUSDXconjSHpp(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSDXSSI0conjUSDXconjSSI0(gO1,gI1,
      gO2,gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSDXconjUSDX(gI1,gI1,gO1,gO2
      ));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSDXconjUSDX(gI1,gI1,gO1,gO2
      ));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSDXSvconjUSDXconjSv(gO1,gI1,gO2,gI1
      ));
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFDXconjUSDXPL(gI2,gI1,gO2))*
      CpChiFDXconjUSDXPL(gI2,gI1,gO1) + Conj(CpChiFDXconjUSDXPR(gI2,gI1,gO2))*
      CpChiFDXconjUSDXPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFDX(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiFDXconjUSDXPR(gI2,gI1,gO2))*CpChiFDXconjUSDXPL(gI2,
      gI1,gO1) + Conj(CpChiFDXconjUSDXPL(gI2,gI1,gO2))*CpChiFDXconjUSDXPR(gI2,gI1,
      gO1))*MChi(gI2)));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUSDXSHI0conjUSDXconjSHI0(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUSDXSHIpconjUSDXconjSHIp(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSDXconjUSDXconjSdSd(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUSDXconjUSDXconjSDXSDX(gO1,gO2,gI1,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSDXconjUSDXconjSuSu(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUSDXSeconjUSDXconjSe(gO1,gI1,gO2,gI1
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSDX(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSDXconjUSDX(gI2,gI1,gO2))*CpAhSDXconjUSDX(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSDX(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSDXconjUSDX(gI2,gI1,gO2))*CphhSDXconjUSDX(gI2,gI1,gO1)));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpGluFDXconjUSDXPL(gI2,gO2))*
      CpGluFDXconjUSDXPL(gI2,gO1) + Conj(CpGluFDXconjUSDXPR(gI2,gO2))*
      CpGluFDXconjUSDXPR(gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFDX(gI2))));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFDX(gI2
      )))*(Conj(CpGluFDXconjUSDXPR(gI2,gO2))*CpGluFDXconjUSDXPL(gI2,gO1) + Conj(
      CpGluFDXconjUSDXPL(gI2,gO2))*CpGluFDXconjUSDXPR(gI2,gO1))*MFDX(gI2));
   result += 1.3333333333333333*SUM(gI2,0,5,Conj(CpSDXconjUSDXVG(gI2,gO2))*
      CpSDXconjUSDXVG(gI2,gO1)*F0(Sqr(p),Sqr(MSDX(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSDXconjUSDXVP(gI2,gO2))*CpSDXconjUSDXVP(gI2,gO1)*
      F0(Sqr(p),Sqr(MSDX(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSDXconjUSDXVZ(gI2,gO2))*CpSDXconjUSDXVZ(gI2,gO1)*
      F0(Sqr(p),Sqr(MSDX(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSDXconjUSDXVZp(gI2,gO2))*CpSDXconjUSDXVZp(gI2,gO1)
      *F0(Sqr(p),Sqr(MSDX(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_SDX_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = i; k < 6; k++)
         self_energy(i, k) = self_energy_SDX_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_hh_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmCgWmCUhh(gO1)*
      CpbargWmCgWmCUhh(gO2));
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmgWmUhh(gO1)*CpbargWmgWmUhh(
      gO2));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*CpbargZgZUhh(gO1)*CpbargZgZUhh(gO2));
   result += -(B0(Sqr(p),Sqr(MVZp),Sqr(MVZp))*CpbargZpgZpUhh(gO1)*CpbargZpgZpUhh(
      gO2));
   result += -2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZp))*CpbargZpgZUhh(gO1)*CpbargZpgZUhh(gO2
      );
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*Conj(CpUhhconjVWmVWm(gO2))*
      CpUhhconjVWmVWm(gO1);
   result += 4*A0(Sqr(MVWm))*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUhhUhhVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUhhUhhVZVZ(gO1,gO2);
   result += 2*B0(Sqr(p),Sqr(MVZp),Sqr(MVZp))*Conj(CpUhhVZpVZp(gO2))*CpUhhVZpVZp(
      gO1);
   result += 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
   result += 4*B0(Sqr(p),Sqr(MVZ),Sqr(MVZp))*Conj(CpUhhVZVZp(gO2))*CpUhhVZVZp(gO1)
      ;
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpUhhUhhHpmconjHpm(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUhhUhhSHp0conjSHp0(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUhhUhhSHppconjSHpp(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUhhUhhSSI0conjSSI0(gO1,gO2,gI1,gI1
      ));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*Conj
      (CpUhhHpmconjHpm(gO2,gI2,gI1))*CpUhhHpmconjHpm(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSHp0(gI1)),Sqr(MSHp0(gI2)))*
      Conj(CpUhhSHp0conjSHp0(gO2,gI2,gI1))*CpUhhSHp0conjSHp0(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSHpp(gI1)),Sqr(MSHpp(gI2)))*
      Conj(CpUhhSHppconjSHpp(gO2,gI2,gI1))*CpUhhSHppconjSHpp(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSSI0(gI1)),Sqr(MSSI0(gI2)))*
      Conj(CpUhhSSI0conjSSI0(gO2,gI2,gI1))*CpUhhSSI0conjSSI0(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaChaUhhPL(gI1,gI2,gO2))*
      CpbarChaChaUhhPL(gI1,gI2,gO1) + Conj(CpbarChaChaUhhPR(gI1,gI2,gO2))*
      CpbarChaChaUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaIChaIUhhPL(gI1,gI2,gO2))*
      CpbarChaIChaIUhhPL(gI1,gI2,gO1) + Conj(CpbarChaIChaIUhhPR(gI1,gI2,gO2))*
      CpbarChaIChaIUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpbarChaChaUhhPR(gI1,gI2,gO2))*CpbarChaChaUhhPL(gI1,gI2,
      gO1) + Conj(CpbarChaChaUhhPL(gI1,gI2,gO2))*CpbarChaChaUhhPR(gI1,gI2,gO1))*
      MCha(gI2)));
   result += -2*SUM(gI1,0,1,MChaI(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(
      MChaI(gI2)))*(Conj(CpbarChaIChaIUhhPR(gI1,gI2,gO2))*CpbarChaIChaIUhhPL(gI1,
      gI2,gO1) + Conj(CpbarChaIChaIUhhPL(gI1,gI2,gO2))*CpbarChaIChaIUhhPR(gI1,gI2,
      gO1))*MChaI(gI2)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhhUhh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhhUhh(gI1,gI1,gO1,gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUhhUhhSvconjSv(gO1,gO2,gI1,gI1));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUhh(gI1,gI2,gO2))*CpAhAhUhh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhhhUhh(gI2,gI1,gO2))*CpAhhhUhh(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhhhUhh(gI1,gI2,gO2))*CphhhhUhh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSv(gI2)))*Conj(
      CpUhhSvconjSv(gO2,gI2,gI1))*CpUhhSvconjSv(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUhhPL(gI1,gI2,gO2))*
      CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*
      CpbarFdFdUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFDXFDXUhhPL(gI1,gI2,gO2))*
      CpbarFDXFDXUhhPL(gI1,gI2,gO1) + Conj(CpbarFDXFDXUhhPR(gI1,gI2,gO2))*
      CpbarFDXFDXUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUhhPL(gI1,gI2,gO2))*
      CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*
      CpbarFeFeUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUhhPL(gI1,gI2,gO2))*
      CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*
      CpbarFuFuUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUhhPL(gI1,gI2,gO2))*CpbarFdFdUhhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -6*SUM(gI1,0,2,MFDX(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(
      MFDX(gI2)))*(Conj(CpbarFDXFDXUhhPR(gI1,gI2,gO2))*CpbarFDXFDXUhhPL(gI1,gI2,
      gO1) + Conj(CpbarFDXFDXUhhPL(gI1,gI2,gO2))*CpbarFDXFDXUhhPR(gI1,gI2,gO1))*
      MFDX(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUhhPL(gI1,gI2,gO2))*CpbarFeFeUhhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUhhPL(gI1,gI2,gO2))*CpbarFuFuUhhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUhhUhhSHI0conjSHI0(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUhhUhhSHIpconjSHIp(gO1,gO2,gI1,gI1
      ));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSHI0(gI1)),Sqr(MSHI0(gI2)))*
      Conj(CpUhhSHI0conjSHI0(gO2,gI2,gI1))*CpUhhSHI0conjSHI0(gO1,gI2,gI1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSHIp(gI1)),Sqr(MSHIp(gI2)))*
      Conj(CpUhhSHIpconjSHIp(gO2,gI2,gI1))*CpUhhSHIpconjSHIp(gO1,gI2,gI1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,(Conj(CpChiIChiIUhhPL(gI1,gI2,gO2))*
      CpChiIChiIUhhPL(gI1,gI2,gO1) + Conj(CpChiIChiIUhhPR(gI1,gI2,gO2))*
      CpChiIChiIUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChiI(gI2)))));
   result += -SUM(gI1,0,3,MChiI(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(
      MChiI(gI2)))*(Conj(CpChiIChiIUhhPR(gI1,gI2,gO2))*CpChiIChiIUhhPL(gI1,gI2,gO1
      ) + Conj(CpChiIChiIUhhPL(gI1,gI2,gO2))*CpChiIChiIUhhPR(gI1,gI2,gO1))*MChiI(
      gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhUhhSdconjSd(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUhhUhhSDXconjSDX(gO1,gO2,gI1,gI1)
      );
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhUhhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhUhhSuconjSu(gO1,gO2,gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*Conj
      (CpUhhSdconjSd(gO2,gI2,gI1))*CpUhhSdconjSd(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSDX(gI1)),Sqr(MSDX(gI2)))*
      Conj(CpUhhSDXconjSDX(gO2,gI2,gI1))*CpUhhSDXconjSDX(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*Conj(
      CpUhhSeconjSe(gO2,gI2,gI1))*CpUhhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*Conj
      (CpUhhSuconjSu(gO2,gI2,gI1))*CpUhhSuconjSu(gO1,gI2,gI1)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(Conj(CpChiChiUhhPL(gI1,gI2,gO2))*
      CpChiChiUhhPL(gI1,gI2,gO1) + Conj(CpChiChiUhhPR(gI1,gI2,gO2))*CpChiChiUhhPR(
      gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiChiUhhPR(gI1,gI2,gO2))*CpChiChiUhhPL(gI1,gI2,gO1) + Conj(
      CpChiChiUhhPL(gI1,gI2,gO2))*CpChiChiUhhPR(gI1,gI2,gO1))*MChi(gI2)));
   result += 2*SUM(gI2,0,1,Conj(CpUhhHpmconjVWm(gO2,gI2))*CpUhhHpmconjVWm(gO1,gI2)
      *F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CpAhUhhVZ(gI2,gO2))*CpAhUhhVZ(gI2,gO1)*F0(Sqr(p),Sqr
      (MAh(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,2,Conj(CpAhUhhVZp(gI2,gO2))*CpAhUhhVZp(gI2,gO1)*F0(Sqr(p),
      Sqr(MAh(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_hh_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = i; k < 3; k++)
         self_energy(i, k) = self_energy_hh_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Ah_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmCgWmCUAh(gO1)*
      CpbargWmCgWmCUAh(gO2));
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmgWmUAh(gO1)*CpbargWmgWmUAh(
      gO2));
   result += 4*A0(Sqr(MVWm))*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUAhUAhVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUAhUAhVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpUAhUAhHpmconjHpm(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUAhUAhSHp0conjSHp0(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUAhUAhSHppconjSHpp(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUAhUAhSSI0conjSSI0(gO1,gO2,gI1,gI1
      ));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*Conj
      (CpUAhHpmconjHpm(gO2,gI2,gI1))*CpUAhHpmconjHpm(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaChaUAhPL(gI1,gI2,gO2))*
      CpbarChaChaUAhPL(gI1,gI2,gO1) + Conj(CpbarChaChaUAhPR(gI1,gI2,gO2))*
      CpbarChaChaUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaIChaIUAhPL(gI1,gI2,gO2))*
      CpbarChaIChaIUAhPL(gI1,gI2,gO1) + Conj(CpbarChaIChaIUAhPR(gI1,gI2,gO2))*
      CpbarChaIChaIUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpbarChaChaUAhPR(gI1,gI2,gO2))*CpbarChaChaUAhPL(gI1,gI2,
      gO1) + Conj(CpbarChaChaUAhPL(gI1,gI2,gO2))*CpbarChaChaUAhPR(gI1,gI2,gO1))*
      MCha(gI2)));
   result += -2*SUM(gI1,0,1,MChaI(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(
      MChaI(gI2)))*(Conj(CpbarChaIChaIUAhPR(gI1,gI2,gO2))*CpbarChaIChaIUAhPL(gI1,
      gI2,gO1) + Conj(CpbarChaIChaIUAhPL(gI1,gI2,gO2))*CpbarChaIChaIUAhPR(gI1,gI2,
      gO1))*MChaI(gI2)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUAhUAh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUAhUAhSvconjSv(gO1,gO2,gI1,gI1));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUAh(gI1,gI2,gO2))*CpAhAhUAh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhUAhhh(gI2,gO2,gI1))*CpAhUAhhh(gI2,gO1,gI1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CpUAhhhhh(gO2,gI1,gI2))*CpUAhhhhh(gO1,gI1,gI2)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUAhPL(gI1,gI2,gO2))*
      CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*
      CpbarFdFdUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFDXFDXUAhPL(gI1,gI2,gO2))*
      CpbarFDXFDXUAhPL(gI1,gI2,gO1) + Conj(CpbarFDXFDXUAhPR(gI1,gI2,gO2))*
      CpbarFDXFDXUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUAhPL(gI1,gI2,gO2))*
      CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*
      CpbarFeFeUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUAhPL(gI1,gI2,gO2))*
      CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*
      CpbarFuFuUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUAhPL(gI1,gI2,gO2))*CpbarFdFdUAhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -6*SUM(gI1,0,2,MFDX(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(
      MFDX(gI2)))*(Conj(CpbarFDXFDXUAhPR(gI1,gI2,gO2))*CpbarFDXFDXUAhPL(gI1,gI2,
      gO1) + Conj(CpbarFDXFDXUAhPL(gI1,gI2,gO2))*CpbarFDXFDXUAhPR(gI1,gI2,gO1))*
      MFDX(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUAhPL(gI1,gI2,gO2))*CpbarFeFeUAhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUAhPL(gI1,gI2,gO2))*CpbarFuFuUAhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUAhUAhSHI0conjSHI0(gO1,gO2,gI1,gI1
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUAhUAhSHIpconjSHIp(gO1,gO2,gI1,gI1
      ));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSHI0(gI1)),Sqr(MSHI0(gI2)))*
      Conj(CpUAhSHI0conjSHI0(gO2,gI2,gI1))*CpUAhSHI0conjSHI0(gO1,gI2,gI1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSHIp(gI1)),Sqr(MSHIp(gI2)))*
      Conj(CpUAhSHIpconjSHIp(gO2,gI2,gI1))*CpUAhSHIpconjSHIp(gO1,gI2,gI1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,(Conj(CpChiIChiIUAhPL(gI1,gI2,gO2))*
      CpChiIChiIUAhPL(gI1,gI2,gO1) + Conj(CpChiIChiIUAhPR(gI1,gI2,gO2))*
      CpChiIChiIUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChiI(gI2)))));
   result += -SUM(gI1,0,3,MChiI(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(
      MChiI(gI2)))*(Conj(CpChiIChiIUAhPR(gI1,gI2,gO2))*CpChiIChiIUAhPL(gI1,gI2,gO1
      ) + Conj(CpChiIChiIUAhPL(gI1,gI2,gO2))*CpChiIChiIUAhPR(gI1,gI2,gO1))*MChiI(
      gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUAhUAhSdconjSd(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUAhUAhSDXconjSDX(gO1,gO2,gI1,gI1)
      );
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUAhUAhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUAhUAhSuconjSu(gO1,gO2,gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*Conj
      (CpUAhSdconjSd(gO2,gI2,gI1))*CpUAhSdconjSd(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSDX(gI1)),Sqr(MSDX(gI2)))*
      Conj(CpUAhSDXconjSDX(gO2,gI2,gI1))*CpUAhSDXconjSDX(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*Conj(
      CpUAhSeconjSe(gO2,gI2,gI1))*CpUAhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*Conj
      (CpUAhSuconjSu(gO2,gI2,gI1))*CpUAhSuconjSu(gO1,gI2,gI1)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(Conj(CpChiChiUAhPL(gI1,gI2,gO2))*
      CpChiChiUAhPL(gI1,gI2,gO1) + Conj(CpChiChiUAhPR(gI1,gI2,gO2))*CpChiChiUAhPR(
      gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiChiUAhPR(gI1,gI2,gO2))*CpChiChiUAhPL(gI1,gI2,gO1) + Conj(
      CpChiChiUAhPL(gI1,gI2,gO2))*CpChiChiUAhPR(gI1,gI2,gO1))*MChi(gI2)));
   result += 2*SUM(gI2,0,1,Conj(CpUAhHpmconjVWm(gO2,gI2))*CpUAhHpmconjVWm(gO1,gI2)
      *F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CpUAhhhVZ(gO2,gI2))*CpUAhhhVZ(gO1,gI2)*F0(Sqr(p),Sqr
      (Mhh(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,2,Conj(CpUAhhhVZp(gO2,gI2))*CpUAhhhVZp(gO1,gI2)*F0(Sqr(p),
      Sqr(Mhh(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Ah_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = i; k < 3; k++)
         self_energy(i, k) = self_energy_Ah_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Hpm_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVZ))*CpbargWmgZUHpm(gO2)*
      CpbargZgWmconjUHpm(gO1));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVWm))*CpbargWmCgZconjUHpm(gO1)*
      CpbargZgWmCUHpm(gO2));
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVZp))*CpbargWmgZpUHpm(gO2)*
      CpbargZpgWmconjUHpm(gO1));
   result += -(B0(Sqr(p),Sqr(MVZp),Sqr(MVWm))*CpbargWmCgZpconjUHpm(gO1)*
      CpbargZpgWmCUHpm(gO2));
   result += 4*B0(Sqr(p),0,Sqr(MVWm))*Conj(CpconjUHpmVPVWm(gO2))*CpconjUHpmVPVWm(
      gO1);
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVZ))*Conj(CpconjUHpmVWmVZ(gO2))*
      CpconjUHpmVWmVZ(gO1);
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVZp))*Conj(CpconjUHpmVWmVZp(gO2))*
      CpconjUHpmVWmVZp(gO1);
   result += 4*A0(Sqr(MVWm))*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUHpmconjUHpmVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUHpmconjUHpmVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUHpmconjHpmconjUHpm(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUHpmSHp0conjUHpmconjSHp0(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUHpmSHppconjUHpmconjSHpp(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUHpmSSI0conjUHpmconjSSI0(gO1,gI1,
      gO2,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSHp0(gI1)),Sqr(MSHpp(gI2)))*
      Conj(CpSHppconjUHpmconjSHp0(gI2,gO2,gI1))*CpSHppconjUHpmconjSHp0(gI2,gO1,gI1
      )));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhHpmconjUHpm(gI2,gI1,gO2))*CpAhHpmconjUHpm(gI2,gI1,gO1)));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhHpmconjUHpm(gI2,gI1,gO2))*CphhHpmconjUHpm(gI2,gI1,gO1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUHpmconjUHpm(gI1,gI1,gO1,gO2
      ));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUHpmconjUHpm(gI1,gI1,gO1,gO2
      ));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUHpmSvconjUHpmconjSv(gO1,gI1,gO2,gI1
      ));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFdconjUHpmPL(gI1,gI2,gO2))*
      CpbarFuFdconjUHpmPL(gI1,gI2,gO1) + Conj(CpbarFuFdconjUHpmPR(gI1,gI2,gO2))*
      CpbarFuFdconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFvFeconjUHpmPL(gI1,gI2,gO2))*
      CpbarFvFeconjUHpmPL(gI1,gI2,gO1) + Conj(CpbarFvFeconjUHpmPR(gI1,gI2,gO2))*
      CpbarFvFeconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFuFdconjUHpmPR(gI1,gI2,gO2))*CpbarFuFdconjUHpmPL(gI1,gI2,
      gO1) + Conj(CpbarFuFdconjUHpmPL(gI1,gI2,gO2))*CpbarFuFdconjUHpmPR(gI1,gI2,
      gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFvFeconjUHpmPR(gI1,gI2,gO2))*CpbarFvFeconjUHpmPL(gI1,gI2,
      gO1) + Conj(CpbarFvFeconjUHpmPL(gI1,gI2,gO2))*CpbarFvFeconjUHpmPR(gI1,gI2,
      gO1))*MFe(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSe(gI2)))*Conj(
      CpSeconjUHpmconjSv(gI2,gO2,gI1))*CpSeconjUHpmconjSv(gI2,gO1,gI1)));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUHpmSHI0conjUHpmconjSHI0(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUHpmSHIpconjUHpmconjSHIp(gO1,gI1,
      gO2,gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSHI0(gI1)),Sqr(MSHIp(gI2)))*
      Conj(CpSHIpconjUHpmconjSHI0(gI2,gO2,gI1))*CpSHIpconjUHpmconjSHI0(gI2,gO1,gI1
      )));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUHpmSdconjUHpmconjSd(gO1,gI1,gO2,
      gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUHpmSDXconjUHpmconjSDX(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUHpmSeconjUHpmconjSe(gO1,gI1,gO2,gI1
      ));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUHpmSuconjUHpmconjSu(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,(Conj(CpChiChaconjUHpmPL(gI1,gI2,gO2))*
      CpChiChaconjUHpmPL(gI1,gI2,gO1) + Conj(CpChiChaconjUHpmPR(gI1,gI2,gO2))*
      CpChiChaconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpChiChaconjUHpmPR(gI1,gI2,gO2))*CpChiChaconjUHpmPL(gI1,
      gI2,gO1) + Conj(CpChiChaconjUHpmPL(gI1,gI2,gO2))*CpChiChaconjUHpmPR(gI1,gI2,
      gO1))*MCha(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSd(gI2)))*Conj
      (CpSdconjUHpmconjSu(gI2,gO2,gI1))*CpSdconjUHpmconjSu(gI2,gO1,gI1)));
   result += SUM(gI2,0,1,Conj(CpHpmconjUHpmVP(gI2,gO2))*CpHpmconjUHpmVP(gI2,gO1)*
      F0(Sqr(p),Sqr(MHpm(gI2)),0));
   result += SUM(gI2,0,1,Conj(CpHpmconjUHpmVZ(gI2,gO2))*CpHpmconjUHpmVZ(gI2,gO1)*
      F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,1,Conj(CpHpmconjUHpmVZp(gI2,gO2))*CpHpmconjUHpmVZp(gI2,gO1)
      *F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVZp)));
   result += SUM(gI2,0,2,Conj(CpAhconjUHpmVWm(gI2,gO2))*CpAhconjUHpmVWm(gI2,gO1)*
      F0(Sqr(p),Sqr(MAh(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CphhconjUHpmVWm(gI2,gO2))*CphhconjUHpmVWm(gI2,gO1)*
      F0(Sqr(p),Sqr(Mhh(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Hpm_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = i; k < 2; k++)
         self_energy(i, k) = self_energy_Hpm_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_SHI0_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSHI0conjUSHI0conjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSHI0conjUSHI0VZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSHI0conjUSHI0VZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSHI0conjHpmconjUSHI0(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUSHI0SHp0conjUSHI0conjSHp0(gO1,gI1
      ,gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUSHI0SHppconjUSHI0conjSHpp(gO1,gI1
      ,gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSHI0SSI0conjUSHI0conjSSI0(gO1,gI1
      ,gO2,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaChaIconjUSHI0PL(gI1,gI2,gO2))*
      CpbarChaChaIconjUSHI0PL(gI1,gI2,gO1) + Conj(CpbarChaChaIconjUSHI0PR(gI1,gI2,
      gO2))*CpbarChaChaIconjUSHI0PR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MChaI(gI2)))));
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MChaI(gI2)))*(Conj(CpbarChaChaIconjUSHI0PR(gI1,gI2,gO2))*
      CpbarChaChaIconjUSHI0PL(gI1,gI2,gO1) + Conj(CpbarChaChaIconjUSHI0PL(gI1,gI2,
      gO2))*CpbarChaChaIconjUSHI0PR(gI1,gI2,gO1))*MChaI(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSHIp(gI2)))*
      Conj(CpSHIpconjHpmconjUSHI0(gI2,gI1,gO2))*CpSHIpconjHpmconjUSHI0(gI2,gI1,gO1
      )));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSHI0conjUSHI0(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSHI0conjUSHI0(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSHI0SvconjUSHI0conjSv(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USHI0conjSHI0conjUSHI0(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUSHI0SHIpconjUSHI0conjSHIp(gO1,gI1
      ,gO2,gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSHI0(gI1)),Sqr(MAh(gI2)))*Conj
      (CpAhSHI0conjUSHI0(gI2,gI1,gO2))*CpAhSHI0conjUSHI0(gI2,gI1,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSHI0(gI1)),Sqr(Mhh(gI2)))*Conj
      (CphhSHI0conjUSHI0(gI2,gI1,gO2))*CphhSHI0conjUSHI0(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,5,(Conj(CpChiChiIconjUSHI0PL(gI2,gI1,gO2))*
      CpChiChiIconjUSHI0PL(gI2,gI1,gO1) + Conj(CpChiChiIconjUSHI0PR(gI2,gI1,gO2))*
      CpChiChiIconjUSHI0PR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChi(gI2)))
      ));
   result += -SUM(gI1,0,3,MChiI(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiChiIconjUSHI0PR(gI2,gI1,gO2))*CpChiChiIconjUSHI0PL(
      gI2,gI1,gO1) + Conj(CpChiChiIconjUSHI0PL(gI2,gI1,gO2))*CpChiChiIconjUSHI0PR(
      gI2,gI1,gO1))*MChi(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSHI0conjSdconjUSHI0(gI1,gO1,gI1
      ,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSHI0conjSDXconjUSHI0(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSHI0conjSeconjUSHI0(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSHI0SuconjUSHI0conjSu(gO1,gI1,gO2
      ,gI1));
   result += SUM(gI2,0,3,Conj(CpSHI0conjUSHI0VZ(gI2,gO2))*CpSHI0conjUSHI0VZ(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHI0(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,3,Conj(CpSHI0conjUSHI0VZp(gI2,gO2))*CpSHI0conjUSHI0VZp(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHI0(gI2)),Sqr(MVZp)));
   result += SUM(gI2,0,3,Conj(CpSHIpconjUSHI0conjVWm(gI2,gO2))*
      CpSHIpconjUSHI0conjVWm(gI2,gO1)*F0(Sqr(p),Sqr(MSHIp(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_SHI0_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = i; k < 4; k++)
         self_energy(i, k) = self_energy_SHI0_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_SHIp_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSHIpconjUSHIpconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSHIpconjUSHIpVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSHIpconjUSHIpVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSHIpconjHpmconjUSHIp(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUSHIpSHp0conjUSHIpconjSHp0(gO1,gI1
      ,gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUSHIpSHppconjUSHIpconjSHpp(gO1,gI1
      ,gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSHIpSSI0conjUSHIpconjSSI0(gO1,gI1
      ,gO2,gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSHIpconjUSHIp(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSHIpconjUSHIp(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSHIpSvconjUSHIpconjSv(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USHIpconjSHI0conjUSHIp(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpUSHIpconjSHIpconjUSHIp(gI1,gO1
      ,gI1,gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSHI0(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpHpmSHI0conjUSHIp(gI2,gI1,gO2))*CpHpmSHI0conjUSHIp(gI2,gI1,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,1,(Conj(CpChiIChaconjUSHIpPL(gI1,gI2,gO2))*
      CpChiIChaconjUSHIpPL(gI1,gI2,gO1) + Conj(CpChiIChaconjUSHIpPR(gI1,gI2,gO2))*
      CpChiIChaconjUSHIpPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MCha(gI2)))
      ));
   result += -2*SUM(gI1,0,3,MChiI(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpChiIChaconjUSHIpPR(gI1,gI2,gO2))*CpChiIChaconjUSHIpPL(
      gI1,gI2,gO1) + Conj(CpChiIChaconjUSHIpPL(gI1,gI2,gO2))*CpChiIChaconjUSHIpPR(
      gI1,gI2,gO1))*MCha(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSHIp(gI1)),Sqr(MAh(gI2)))*Conj
      (CpAhSHIpconjUSHIp(gI2,gI1,gO2))*CpAhSHIpconjUSHIp(gI2,gI1,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSHIp(gI1)),Sqr(Mhh(gI2)))*Conj
      (CphhSHIpconjUSHIp(gI2,gI1,gO2))*CphhSHIpconjUSHIp(gI2,gI1,gO1)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSHIpconjSdconjUSHIp(gI1,gO1,gI1
      ,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSHIpconjSDXconjUSHIp(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSHIpconjSeconjUSHIp(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSHIpSuconjUSHIpconjSu(gO1,gI1,gO2
      ,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,(Conj(CpChiChaIconjUSHIpPL(gI1,gI2,gO2))*
      CpChiChaIconjUSHIpPL(gI1,gI2,gO1) + Conj(CpChiChaIconjUSHIpPR(gI1,gI2,gO2))*
      CpChiChaIconjUSHIpPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChaI(gI2)))
      ));
   result += -2*SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MChaI(gI2)))*(Conj(CpChiChaIconjUSHIpPR(gI1,gI2,gO2))*CpChiChaIconjUSHIpPL(
      gI1,gI2,gO1) + Conj(CpChiChaIconjUSHIpPL(gI1,gI2,gO2))*CpChiChaIconjUSHIpPR(
      gI1,gI2,gO1))*MChaI(gI2)));
   result += SUM(gI2,0,3,Conj(CpSHI0conjUSHIpVWm(gI2,gO2))*CpSHI0conjUSHIpVWm(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHI0(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,3,Conj(CpSHIpconjUSHIpVP(gI2,gO2))*CpSHIpconjUSHIpVP(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHIp(gI2)),0));
   result += SUM(gI2,0,3,Conj(CpSHIpconjUSHIpVZ(gI2,gO2))*CpSHIpconjUSHIpVZ(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHIp(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,3,Conj(CpSHIpconjUSHIpVZp(gI2,gO2))*CpSHIpconjUSHIpVZp(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHIp(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_SHIp_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = i; k < 4; k++)
         self_energy(i, k) = self_energy_SHIp_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_SSI0_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 2*A0(Sqr(MVZp))*CpUSSI0conjUSSI0VZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSSI0conjUSSI0VZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSSI0conjHpmconjUSSI0(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0USSI0conjSHp0conjUSSI0(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppUSSI0conjSHppconjUSSI0(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpSSI0USSI0conjSSI0conjUSSI0(gI1,gO1
      ,gI1,gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSSI0(gI1)),Sqr(Mhh(gI2)))*Conj
      (CphhSSI0conjUSSI0(gI2,gI1,gO2))*CphhSSI0conjUSSI0(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,1,SUM(gI2,0,5,(Conj(CpChiFSIconjUSSI0PL(gI2,gI1,gO2))*
      CpChiFSIconjUSSI0PL(gI2,gI1,gO1) + Conj(CpChiFSIconjUSSI0PR(gI2,gI1,gO2))*
      CpChiFSIconjUSSI0PR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFSI(gI1)),Sqr(MChi(gI2)))))
      ;
   result += -SUM(gI1,0,1,MFSI(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFSI(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFSIconjUSSI0PR(gI2,gI1,gO2))*CpChiFSIconjUSSI0PL(gI2,gI1,
      gO1) + Conj(CpChiFSIconjUSSI0PL(gI2,gI1,gO2))*CpChiFSIconjUSSI0PR(gI2,gI1,
      gO1))*MChi(gI2)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSSI0conjUSSI0(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSSI0conjUSSI0(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSSI0SvconjUSSI0conjSv(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USSI0conjSHI0conjUSSI0(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpUSSI0conjSHIpconjUSSI0(gI1,gO1
      ,gI1,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSSI0conjSdconjUSSI0(gI1,gO1,gI1
      ,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSSI0conjSDXconjUSSI0(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSSI0conjSeconjUSSI0(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSSI0SuconjUSSI0conjSu(gO1,gI1,gO2
      ,gI1));
   result += SUM(gI2,0,1,Conj(CpSSI0conjUSSI0VZ(gI2,gO2))*CpSSI0conjUSSI0VZ(gI2,
      gO1)*F0(Sqr(p),Sqr(MSSI0(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,1,Conj(CpSSI0conjUSSI0VZp(gI2,gO2))*CpSSI0conjUSSI0VZp(gI2,
      gO1)*F0(Sqr(p),Sqr(MSSI0(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_SSI0_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = i; k < 2; k++)
         self_energy(i, k) = self_energy_SSI0_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_SHp0_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSHp0conjUSHp0conjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSHp0conjUSHp0VZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSHp0conjUSHp0VZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSHp0conjHpmconjUSHp0(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0USHp0conjSHp0conjUSHp0(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUSHp0SHppconjUSHp0conjSHpp(gO1,gI1
      ,gO2,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSHp0SSI0conjUSHp0conjSSI0(gO1,gI1
      ,gO2,gI1));
   result += SUM(gI1,0,1,(Conj(CpbarChaChaPconjUSHp0PL(gI1,gO2))*
      CpbarChaChaPconjUSHp0PL(gI1,gO1) + Conj(CpbarChaChaPconjUSHp0PR(gI1,gO2))*
      CpbarChaChaPconjUSHp0PR(gI1,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MChaP)));
   result += -2*MChaP*SUM(gI1,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MChaP))*(Conj(
      CpbarChaChaPconjUSHp0PR(gI1,gO2))*CpbarChaChaPconjUSHp0PL(gI1,gO1) + Conj(
      CpbarChaChaPconjUSHp0PL(gI1,gO2))*CpbarChaChaPconjUSHp0PR(gI1,gO1))*MCha(gI1
      ));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSHpp(gI2)))*
      Conj(CpSHppconjHpmconjUSHp0(gI2,gI1,gO2))*CpSHppconjHpmconjUSHp0(gI2,gI1,gO1
      )));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSHp0(gI1)),Sqr(Mhh(gI2)))*Conj
      (CphhSHp0conjUSHp0(gI2,gI1,gO2))*CphhSHp0conjUSHp0(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,1,SUM(gI2,0,5,(Conj(CpChiChiPconjUSHp0PL(gI2,gI1,gO2))*
      CpChiChiPconjUSHp0PL(gI2,gI1,gO1) + Conj(CpChiChiPconjUSHp0PR(gI2,gI1,gO2))*
      CpChiChiPconjUSHp0PR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MChi(gI2)))
      ));
   result += -SUM(gI1,0,1,MChiP(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChiP(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiChiPconjUSHp0PR(gI2,gI1,gO2))*CpChiChiPconjUSHp0PL(
      gI2,gI1,gO1) + Conj(CpChiChiPconjUSHp0PL(gI2,gI1,gO2))*CpChiChiPconjUSHp0PR(
      gI2,gI1,gO1))*MChi(gI2)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSHp0conjUSHp0(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSHp0conjUSHp0(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSHp0SvconjUSHp0conjSv(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USHp0conjSHI0conjUSHp0(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpUSHp0conjSHIpconjUSHp0(gI1,gO1
      ,gI1,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSHp0conjSdconjUSHp0(gI1,gO1,gI1
      ,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSHp0conjSDXconjUSHp0(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSHp0conjSeconjUSHp0(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSHp0SuconjUSHp0conjSu(gO1,gI1,gO2
      ,gI1));
   result += SUM(gI2,0,1,Conj(CpSHp0conjUSHp0VZ(gI2,gO2))*CpSHp0conjUSHp0VZ(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHp0(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,1,Conj(CpSHp0conjUSHp0VZp(gI2,gO2))*CpSHp0conjUSHp0VZp(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHp0(gI2)),Sqr(MVZp)));
   result += SUM(gI2,0,1,Conj(CpSHppconjUSHp0conjVWm(gI2,gO2))*
      CpSHppconjUSHp0conjVWm(gI2,gO1)*F0(Sqr(p),Sqr(MSHpp(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_SHp0_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = i; k < 2; k++)
         self_energy(i, k) = self_energy_SHp0_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_SHpp_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSHppconjUSHppconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZp))*CpUSHppconjUSHppVZpVZp(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSHppconjUSHppVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmUSHppconjHpmconjUSHpp(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0USHppconjSHp0conjUSHpp(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppUSHppconjSHppconjUSHpp(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUSHppSSI0conjUSHppconjSSI0(gO1,gI1
      ,gO2,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSHp0(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpHpmSHp0conjUSHpp(gI2,gI1,gO2))*CpHpmSHp0conjUSHpp(gI2,gI1,gO1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpChiPChaconjUSHppPL(gI1,gI2,gO2))*
      CpChiPChaconjUSHppPL(gI1,gI2,gO1) + Conj(CpChiPChaconjUSHppPR(gI1,gI2,gO2))*
      CpChiPChaconjUSHppPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MCha(gI2)))
      ));
   result += -2*SUM(gI1,0,1,MChiP(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpChiPChaconjUSHppPR(gI1,gI2,gO2))*CpChiPChaconjUSHppPL(
      gI1,gI2,gO1) + Conj(CpChiPChaconjUSHppPL(gI1,gI2,gO2))*CpChiPChaconjUSHppPR(
      gI1,gI2,gO1))*MCha(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSHpp(gI1)),Sqr(Mhh(gI2)))*Conj
      (CphhSHppconjUSHpp(gI2,gI1,gO2))*CphhSHppconjUSHpp(gI2,gI1,gO1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSHppconjUSHpp(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSHppconjUSHpp(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSHppSvconjUSHppconjSv(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0USHppconjSHI0conjUSHpp(gI1,gO1
      ,gI1,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpUSHppconjSHIpconjUSHpp(gI1,gO1
      ,gI1,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSHppconjSdconjUSHpp(gI1,gO1,gI1
      ,gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXUSHppconjSDXconjUSHpp(gI1,gO1,
      gI1,gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSHppconjSeconjUSHpp(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSHppSuconjUSHppconjSu(gO1,gI1,gO2
      ,gI1));
   result += SUM(gI1,0,5,(Conj(CpChiChaPconjUSHppPL(gI1,gO2))*CpChiChaPconjUSHppPL
      (gI1,gO1) + Conj(CpChiChaPconjUSHppPR(gI1,gO2))*CpChiChaPconjUSHppPR(gI1,gO1
      ))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChaP)));
   result += -2*MChaP*SUM(gI1,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChaP))*(Conj(
      CpChiChaPconjUSHppPR(gI1,gO2))*CpChiChaPconjUSHppPL(gI1,gO1) + Conj(
      CpChiChaPconjUSHppPL(gI1,gO2))*CpChiChaPconjUSHppPR(gI1,gO1))*MChi(gI1));
   result += SUM(gI2,0,1,Conj(CpSHp0conjUSHppVWm(gI2,gO2))*CpSHp0conjUSHppVWm(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHp0(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,1,Conj(CpSHppconjUSHppVP(gI2,gO2))*CpSHppconjUSHppVP(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHpp(gI2)),0));
   result += SUM(gI2,0,1,Conj(CpSHppconjUSHppVZ(gI2,gO2))*CpSHppconjUSHppVZ(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHpp(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,1,Conj(CpSHppconjUSHppVZp(gI2,gO2))*CpSHppconjUSHppVZp(gI2,
      gO1)*F0(Sqr(p),Sqr(MSHpp(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_SHpp_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = i; k < 2; k++)
         self_energy(i, k) = self_energy_SHpp_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_VG_1loop(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpbargGgGVG())*B00(Sqr(p),Sqr(MVG),Sqr(MVG));
   result += -3*AbsSqr(CpVGVGVG())*(5*B00(Sqr(p),0,0) + 2*B0(Sqr(p),0,0)*Sqr(p));
   result += 0;
   result += 1.5*(AbsSqr(CpGluGluVGPL()) + AbsSqr(CpGluGluVGPR()))*H0(Sqr(p),Sqr(
      MGlu),Sqr(MGlu)) + 6*B0(Sqr(p),Sqr(MGlu),Sqr(MGlu))*Re(Conj(CpGluGluVGPL())*
      CpGluGluVGPR())*Sqr(MGlu);
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVGPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVGPL(gI1,
      gI2))*CpbarFdFdVGPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFDXFDXVGPL(gI1,gI2)) +
      AbsSqr(CpbarFDXFDXVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2))) +
      4*B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2)))*MFDX(gI1)*MFDX(gI2)*Re(Conj(
      CpbarFDXFDXVGPL(gI1,gI2))*CpbarFDXFDXVGPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVGPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVGPL(gI1,
      gI2))*CpbarFuFuVGPR(gI1,gI2))));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVGVG(gI1,gI1));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXconjSDXVGVG(gI1,gI1));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVGVG(gI1,gI1));
   result += -2*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVG(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -2*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSDXconjSDXVG(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MSDX(gI1)),Sqr(MSDX(gI2)))));
   result += -2*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVG(gI2,gI1))*B00(Sqr(p),
      Sqr(MSu(gI1)),Sqr(MSu(gI2)))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VP_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWmCgWmCVP())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += AbsSqr(CpbargWmgWmVP())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += -(A0(Sqr(MVWm))*(CpconjVWmVPVPVWm1() + CpconjVWmVPVPVWm2() + 4*
      CpconjVWmVPVPVWm3()));
   result += (AbsSqr(CpbarChaPChaPVPPL()) + AbsSqr(CpbarChaPChaPVPPR()))*H0(Sqr(p)
      ,Sqr(MChaP),Sqr(MChaP));
   result += -2*AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm),
      Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += 4*B0(Sqr(p),Sqr(MChaP),Sqr(MChaP))*Re(Conj(CpbarChaPChaPVPPL())*
      CpbarChaPChaPVPPR())*Sqr(MChaP);
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVPVP(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppconjSHppVPVP(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHpmconjHpmVP(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSHppconjSHppVP(gI2,gI1))*B00(Sqr(
      p),Sqr(MSHpp(gI1)),Sqr(MSHpp(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaChaVPPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVPPL(gI1,gI2))*CpbarChaChaVPPR(gI1,gI2))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaIChaIVPPL(gI1,gI2)) + AbsSqr(
      CpbarChaIChaIVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2))) + 4*
      B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2)))*MChaI(gI1)*MChaI(gI2)*Re(Conj(
      CpbarChaIChaIVPPL(gI1,gI2))*CpbarChaIChaIVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVPPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVPPL(gI1,
      gI2))*CpbarFdFdVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFDXFDXVPPL(gI1,gI2)) + AbsSqr(
      CpbarFDXFDXVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2))) + 4*B0(
      Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2)))*MFDX(gI1)*MFDX(gI2)*Re(Conj(
      CpbarFDXFDXVPPL(gI1,gI2))*CpbarFDXFDXVPPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVPPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVPPL(gI1,
      gI2))*CpbarFeFeVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVPPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVPPL(gI1,
      gI2))*CpbarFuFuVPPR(gI1,gI2))));
   result += SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpconjSHIpVPVP(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpSHIpconjSHIpVP(gI2,gI1))*B00(Sqr(
      p),Sqr(MSHIp(gI1)),Sqr(MSHIp(gI2)))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVPVP(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXconjSDXVPVP(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVPVP(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVPVP(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSDXconjSDXVP(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSDX(gI1)),Sqr(MSDX(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += 2*SUM(gI2,0,1,AbsSqr(CpHpmconjVWmVP(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(
      MHpm(gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VZ_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWmCgWmCVZ())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += AbsSqr(CpbargWmgWmVZ())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += -(A0(Sqr(MVWm))*(4*CpconjVWmVWmVZVZ1() + CpconjVWmVWmVZVZ2() +
      CpconjVWmVWmVZVZ3()));
   result += (AbsSqr(CpbarChaPChaPVZPL()) + AbsSqr(CpbarChaPChaPVZPR()))*H0(Sqr(p)
      ,Sqr(MChaP),Sqr(MChaP));
   result += -2*AbsSqr(CpconjVWmVWmVZ())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm),
      Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += 4*B0(Sqr(p),Sqr(MChaP),Sqr(MChaP))*Re(Conj(CpbarChaPChaPVZPL())*
      CpbarChaPChaPVZPR())*Sqr(MChaP);
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVZVZ(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0conjSHp0VZVZ(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppconjSHppVZVZ(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpSSI0conjSSI0VZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHpmconjHpmVZ(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSHp0conjSHp0VZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSHp0(gI1)),Sqr(MSHp0(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSHppconjSHppVZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSHpp(gI1)),Sqr(MSHpp(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSSI0conjSSI0VZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSSI0(gI1)),Sqr(MSSI0(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaChaVZPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVZPL(gI1,gI2))*CpbarChaChaVZPR(gI1,gI2))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaIChaIVZPL(gI1,gI2)) + AbsSqr(
      CpbarChaIChaIVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2))) + 4*
      B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2)))*MChaI(gI1)*MChaI(gI2)*Re(Conj(
      CpbarChaIChaIVZPL(gI1,gI2))*CpbarChaIChaIVZPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpChiPChiPVZPL(gI1,gI2)) + AbsSqr
      (CpChiPChiPVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MChiP(gI2))) + 4*B0
      (Sqr(p),Sqr(MChiP(gI1)),Sqr(MChiP(gI2)))*MChiP(gI1)*MChiP(gI2)*Re(Conj(
      CpChiPChiPVZPL(gI1,gI2))*CpChiPChiPVZPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpFSIFSIVZPL(gI1,gI2)) + AbsSqr(
      CpFSIFSIVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFSI(gI1)),Sqr(MFSI(gI2))) + 4*B0(Sqr(
      p),Sqr(MFSI(gI1)),Sqr(MFSI(gI2)))*MFSI(gI1)*MFSI(gI2)*Re(Conj(CpFSIFSIVZPL(
      gI1,gI2))*CpFSIFSIVZPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhVZVZ(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhVZVZ(gI1,gI1));
   result += SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvconjSvVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpAhhhVZ(gI2,gI1))*B00(Sqr(p),Sqr(
      MAh(gI2)),Sqr(Mhh(gI1)))));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpSvconjSvVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSv(gI1)),Sqr(MSv(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZPL(gI1,
      gI2))*CpbarFdFdVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFDXFDXVZPL(gI1,gI2)) + AbsSqr(
      CpbarFDXFDXVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2))) + 4*B0(
      Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2)))*MFDX(gI1)*MFDX(gI2)*Re(Conj(
      CpbarFDXFDXVZPL(gI1,gI2))*CpbarFDXFDXVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVZPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVZPL(gI1,
      gI2))*CpbarFeFeVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVZPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVZPL(gI1,
      gI2))*CpbarFuFuVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFvVZPL(gI1,gI2)) + AbsSqr(
      CpbarFvFvVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFv(gI1)),Sqr(MFv(gI2)))*MFv(gI1)*MFv(gI2)*Re(Conj(CpbarFvFvVZPL(gI1,
      gI2))*CpbarFvFvVZPR(gI1,gI2))));
   result += SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0conjSHI0VZVZ(gI1,gI1));
   result += SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpconjSHIpVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpSHI0conjSHI0VZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSHI0(gI1)),Sqr(MSHI0(gI2)))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpSHIpconjSHIpVZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSHIp(gI1)),Sqr(MSHIp(gI2)))));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,(AbsSqr(CpChiIChiIVZPL(gI1,gI2)) + AbsSqr
      (CpChiIChiIVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChiI(gI2))) + 4*B0
      (Sqr(p),Sqr(MChiI(gI1)),Sqr(MChiI(gI2)))*MChiI(gI1)*MChiI(gI2)*Re(Conj(
      CpChiIChiIVZPL(gI1,gI2))*CpChiIChiIVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVZVZ(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXconjSDXVZVZ(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVZVZ(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVZVZ(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSDXconjSDXVZ(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSDX(gI1)),Sqr(MSDX(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(AbsSqr(CpChiChiVZPL(gI1,gI2)) + AbsSqr(
      CpChiChiVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2))) + 4*B0(Sqr(
      p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))*MChi(gI1)*MChi(gI2)*Re(Conj(CpChiChiVZPL(
      gI1,gI2))*CpChiChiVZPR(gI1,gI2))));
   result += 2*SUM(gI2,0,1,AbsSqr(CpHpmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(
      MHpm(gI2))));
   result += SUM(gI2,0,2,AbsSqr(CphhVZVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh(gI2))));
   result += SUM(gI2,0,2,AbsSqr(CphhVZVZp(gI2))*B0(Sqr(p),Sqr(MVZp),Sqr(Mhh(gI2)))
      );

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VZp_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWmCgWmCVZp())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += AbsSqr(CpbargWmgWmVZp())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += -(A0(Sqr(MVWm))*(4*CpconjVWmVWmVZpVZp1() + CpconjVWmVWmVZpVZp2() +
      CpconjVWmVWmVZpVZp3()));
   result += (AbsSqr(CpbarChaPChaPVZpPL()) + AbsSqr(CpbarChaPChaPVZpPR()))*H0(Sqr(
      p),Sqr(MChaP),Sqr(MChaP));
   result += -2*AbsSqr(CpconjVWmVWmVZp())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm),
      Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += 4*B0(Sqr(p),Sqr(MChaP),Sqr(MChaP))*Re(Conj(CpbarChaPChaPVZpPL())*
      CpbarChaPChaPVZpPR())*Sqr(MChaP);
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVZpVZp(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0conjSHp0VZpVZp(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppconjSHppVZpVZp(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpSSI0conjSSI0VZpVZp(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHpmconjHpmVZp(gI2,gI1))*B00(Sqr(p
      ),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSHp0conjSHp0VZp(gI2,gI1))*B00(Sqr
      (p),Sqr(MSHp0(gI1)),Sqr(MSHp0(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSHppconjSHppVZp(gI2,gI1))*B00(Sqr
      (p),Sqr(MSHpp(gI1)),Sqr(MSHpp(gI2)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSSI0conjSSI0VZp(gI2,gI1))*B00(Sqr
      (p),Sqr(MSSI0(gI1)),Sqr(MSSI0(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaChaVZpPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVZpPL(gI1,gI2))*CpbarChaChaVZpPR(gI1,gI2))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaIChaIVZpPL(gI1,gI2)) + AbsSqr
      (CpbarChaIChaIVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2))) +
      4*B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MChaI(gI2)))*MChaI(gI1)*MChaI(gI2)*Re(Conj(
      CpbarChaIChaIVZpPL(gI1,gI2))*CpbarChaIChaIVZpPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpChiPChiPVZpPL(gI1,gI2)) +
      AbsSqr(CpChiPChiPVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MChiP(gI2)))
      + 4*B0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MChiP(gI2)))*MChiP(gI1)*MChiP(gI2)*Re(Conj
      (CpChiPChiPVZpPL(gI1,gI2))*CpChiPChiPVZpPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpFSIFSIVZpPL(gI1,gI2)) + AbsSqr(
      CpFSIFSIVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFSI(gI1)),Sqr(MFSI(gI2))) + 4*B0(Sqr
      (p),Sqr(MFSI(gI1)),Sqr(MFSI(gI2)))*MFSI(gI1)*MFSI(gI2)*Re(Conj(CpFSIFSIVZpPL
      (gI1,gI2))*CpFSIFSIVZpPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhVZpVZp(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhVZpVZp(gI1,gI1));
   result += SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvconjSvVZpVZp(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpAhhhVZp(gI2,gI1))*B00(Sqr(p),Sqr(
      MAh(gI2)),Sqr(Mhh(gI1)))));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpSvconjSvVZp(gI2,gI1))*B00(Sqr(p),
      Sqr(MSv(gI1)),Sqr(MSv(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZpPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(
      p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZpPL(gI1
      ,gI2))*CpbarFdFdVZpPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFDXFDXVZpPL(gI1,gI2)) + AbsSqr
      (CpbarFDXFDXVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2))) + 4*B0
      (Sqr(p),Sqr(MFDX(gI1)),Sqr(MFDX(gI2)))*MFDX(gI1)*MFDX(gI2)*Re(Conj(
      CpbarFDXFDXVZpPL(gI1,gI2))*CpbarFDXFDXVZpPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVZpPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(
      p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVZpPL(gI1
      ,gI2))*CpbarFeFeVZpPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVZpPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVZpPL(gI1
      ,gI2))*CpbarFuFuVZpPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFvVZpPL(gI1,gI2)) + AbsSqr(
      CpbarFvFvVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2))) + 4*B0(Sqr(
      p),Sqr(MFv(gI1)),Sqr(MFv(gI2)))*MFv(gI1)*MFv(gI2)*Re(Conj(CpbarFvFvVZpPL(gI1
      ,gI2))*CpbarFvFvVZpPR(gI1,gI2))));
   result += SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0conjSHI0VZpVZp(gI1,gI1));
   result += SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpconjSHIpVZpVZp(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpSHI0conjSHI0VZp(gI2,gI1))*B00(Sqr
      (p),Sqr(MSHI0(gI1)),Sqr(MSHI0(gI2)))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpSHIpconjSHIpVZp(gI2,gI1))*B00(Sqr
      (p),Sqr(MSHIp(gI1)),Sqr(MSHIp(gI2)))));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,(AbsSqr(CpChiIChiIVZpPL(gI1,gI2)) +
      AbsSqr(CpChiIChiIVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChiI(gI2)))
      + 4*B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChiI(gI2)))*MChiI(gI1)*MChiI(gI2)*Re(Conj
      (CpChiIChiIVZpPL(gI1,gI2))*CpChiIChiIVZpPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVZpVZp(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpSDXconjSDXVZpVZp(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVZpVZp(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVZpVZp(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVZp(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSDXconjSDXVZp(gI2,gI1))*B00(Sqr(
      p),Sqr(MSDX(gI1)),Sqr(MSDX(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVZp(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVZp(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,(AbsSqr(CpChiChiVZpPL(gI1,gI2)) + AbsSqr(
      CpChiChiVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2))) + 4*B0(Sqr
      (p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))*MChi(gI1)*MChi(gI2)*Re(Conj(CpChiChiVZpPL
      (gI1,gI2))*CpChiChiVZpPR(gI1,gI2))));
   result += 2*SUM(gI2,0,1,AbsSqr(CpHpmconjVWmVZp(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(
      MHpm(gI2))));
   result += SUM(gI2,0,2,AbsSqr(CphhVZpVZp(gI2))*B0(Sqr(p),Sqr(MVZp),Sqr(Mhh(gI2))
      ));
   result += SUM(gI2,0,2,AbsSqr(CphhVZVZp(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh(gI2))))
      ;

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VWm_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargPgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVP));
   result += AbsSqr(CpbargWmCgPconjVWm())*B00(Sqr(p),Sqr(MVP),Sqr(MVWm));
   result += AbsSqr(CpbargWmCgZconjVWm())*B00(Sqr(p),Sqr(MVZ),Sqr(MVWm));
   result += AbsSqr(CpbargWmCgZpconjVWm())*B00(Sqr(p),Sqr(MVZp),Sqr(MVWm));
   result += AbsSqr(CpbargZgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVZ));
   result += AbsSqr(CpbargZpgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVZp));
   result += -(A0(Sqr(MVWm))*(CpconjVWmconjVWmVWmVWm1() + 4*
      CpconjVWmconjVWmVWmVWm2() + CpconjVWmconjVWmVWmVWm3()));
   result += 0;
   result += -0.5*A0(Sqr(MVZp))*(4*CpconjVWmVWmVZpVZp1() + CpconjVWmVWmVZpVZp2() +
      CpconjVWmVWmVZpVZp3());
   result += -0.5*A0(Sqr(MVZ))*(4*CpconjVWmVWmVZVZ1() + CpconjVWmVWmVZVZ2() +
      CpconjVWmVWmVZVZ3());
   result += -(AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 10*B00(Sqr(p),Sqr(MVWm),0
      ) + B0(Sqr(p),Sqr(MVWm),0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += AbsSqr(CpconjVWmVWmVZ())*(-A0(Sqr(MVWm)) - A0(Sqr(MVZ)) - 10*B00(Sqr(
      p),Sqr(MVZ),Sqr(MVWm)) - B0(Sqr(p),Sqr(MVZ),Sqr(MVWm))*(Sqr(MVWm) + Sqr(MVZ)
      + 4*Sqr(p)));
   result += AbsSqr(CpconjVWmVWmVZp())*(-A0(Sqr(MVWm)) - A0(Sqr(MVZp)) - 10*B00(
      Sqr(p),Sqr(MVZp),Sqr(MVWm)) - B0(Sqr(p),Sqr(MVZp),Sqr(MVWm))*(Sqr(MVWm) +
      Sqr(MVZp) + 4*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpSHp0conjSHp0conjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpSHppconjSHppconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,1,(AbsSqr(CpChiPChaPconjVWmPL(gI1)) + AbsSqr(
      CpChiPChaPconjVWmPR(gI1)))*H0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MChaP)) + 4*MChaP*
      B0(Sqr(p),Sqr(MChiP(gI1)),Sqr(MChaP))*MChiP(gI1)*Re(Conj(CpChiPChaPconjVWmPL
      (gI1))*CpChiPChaPconjVWmPR(gI1)));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpSHppconjSHp0conjVWm(gI2,gI1))*B00
      (Sqr(p),Sqr(MSHpp(gI2)),Sqr(MSHp0(gI1)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,2,AbsSqr(CpAhHpmconjVWm(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MAh(gI2)),Sqr(MHpm(gI1)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,2,AbsSqr(CphhHpmconjVWm(gI2,gI1))*B00(Sqr(p)
      ,Sqr(Mhh(gI2)),Sqr(MHpm(gI1)))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhconjVWmVWm(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvconjSvconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFdconjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarFuFdconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))
      + 4*B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))*MFd(gI2)*MFu(gI1)*Re(Conj(
      CpbarFuFdconjVWmPL(gI1,gI2))*CpbarFuFdconjVWmPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFeconjVWmPL(gI1,gI2)) + AbsSqr
      (CpbarFvFeconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2))) + 4*B0
      (Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))*MFe(gI2)*MFv(gI1)*Re(Conj(
      CpbarFvFeconjVWmPL(gI1,gI2))*CpbarFvFeconjVWmPR(gI1,gI2))));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpSeconjSvconjVWm(gI2,gI1))*B00(Sqr
      (p),Sqr(MSe(gI2)),Sqr(MSv(gI1)))));
   result += SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpSHI0conjSHI0conjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpSHIpconjSHIpconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,1,(AbsSqr(CpChiIChaIconjVWmPL(gI1,gI2)) +
      AbsSqr(CpChiIChaIconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChaI(
      gI2))) + 4*B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MChaI(gI2)))*MChaI(gI2)*MChiI(gI1)*
      Re(Conj(CpChiIChaIconjVWmPL(gI1,gI2))*CpChiIChaIconjVWmPR(gI1,gI2))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpSHIpconjSHI0conjVWm(gI2,gI1))*B00
      (Sqr(p),Sqr(MSHIp(gI2)),Sqr(MSHI0(gI1)))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,(AbsSqr(CpChiChaconjVWmPL(gI1,gI2)) + AbsSqr(
      CpChiChaconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2))) + 4*B0
      (Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2)))*MCha(gI2)*MChi(gI1)*Re(Conj(
      CpChiChaconjVWmPL(gI1,gI2))*CpChiChaconjVWmPR(gI1,gI2))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSuconjVWm(gI2,gI1))*B00(
      Sqr(p),Sqr(MSd(gI2)),Sqr(MSu(gI1)))));
   result += SUM(gI2,0,1,AbsSqr(CpHpmconjVWmVP(gI2))*B0(Sqr(p),0,Sqr(MHpm(gI2))));
   result += SUM(gI2,0,1,AbsSqr(CpHpmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(MHpm(
      gI2))));
   result += SUM(gI2,0,1,AbsSqr(CpHpmconjVWmVZp(gI2))*B0(Sqr(p),Sqr(MVZp),Sqr(MHpm
      (gI2))));
   result += SUM(gI2,0,2,AbsSqr(CphhconjVWmVWm(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(Mhh(
      gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += MChaP*SUM(gI1,0,1,B0(Sqr(p),Sqr(MChaP),Sqr(MSHpp(gI1)))*Conj(
      CpUChiChaPconjSHppPL(gO2,gI1))*CpUChiChaPconjSHppPR(gO1,gI1));
   result += -4*SUM(gI1,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPL(gI1,gO2))*CpbarChaUChiVWmPR(gI1,gO1)*MCha(gI1));
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(
      gI2)))*Conj(CpbarChaUChiHpmPL(gI1,gO2,gI2))*CpbarChaUChiHpmPR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))*Conj
      (CpUChiChaconjHpmPL(gO2,gI2,gI1))*CpUChiChaconjHpmPR(gO1,gI2,gI1)*MCha(gI2))
      );
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHp0(gI1)))*
      Conj(CpUChiChiPconjSHp0PL(gO2,gI2,gI1))*CpUChiChiPconjSHp0PR(gO1,gI2,gI1)*
      MChiP(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHp0(gI1)))*
      Conj(CpUChiChiPSHp0PL(gO2,gI2,gI1))*CpUChiChiPSHp0PR(gO1,gI2,gI1)*MChiP(gI2)
      ));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MFSI(gI2)),Sqr(MSSI0(gI1)))*
      Conj(CpUChiFSIconjSSI0PL(gO2,gI2,gI1))*CpUChiFSIconjSSI0PR(gO1,gI2,gI1)*MFSI
      (gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MFSI(gI2)),Sqr(MSSI0(gI1)))*
      Conj(CpUChiFSISSI0PL(gO2,gI2,gI1))*CpUChiFSISSI0PR(gO1,gI2,gI1)*MFSI(gI2)));
   result += SUM(gI1,0,1,MChaI(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(
      MSHIp(gI2)))*Conj(CpbarChaIUChiSHIpPL(gI1,gO2,gI2))*CpbarChaIUChiSHIpPR(gI1,
      gO1,gI2)));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)
      ))*Conj(CpbarFvUChiSvPL(gI1,gO2,gI2))*CpbarFvUChiSvPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*Conj(
      CpUChiFvconjSvPL(gO2,gI2,gI1))*CpUChiFvconjSvPR(gO1,gI2,gI1)*MFv(gI2)));
   result += 3*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarFdUChiSdPL(gI1,gO2,gI2))*CpbarFdUChiSdPR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,2,MFDX(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MSDX
      (gI2)))*Conj(CpbarFDXUChiSDXPL(gI1,gO2,gI2))*CpbarFDXUChiSDXPR(gI1,gO1,gI2))
      );
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarFeUChiSePL(gI1,gO2,gI2))*CpbarFeUChiSePR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(
      gI2)))*Conj(CpbarFuUChiSuPL(gI1,gO2,gI2))*CpbarFuUChiSuPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpChiUChihhPL(gI2,gO2,gI1))*CpChiUChihhPR(gI2,gO1,gI1)*MChi(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),Sqr(MSHIp(gI1)))*
      Conj(CpUChiChaIconjSHIpPL(gO2,gI2,gI1))*CpUChiChaIconjSHIpPR(gO1,gI2,gI1)*
      MChaI(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHI0(gI1)))*
      Conj(CpUChiChiIconjSHI0PL(gO2,gI2,gI1))*CpUChiChiIconjSHI0PR(gO1,gI2,gI1)*
      MChiI(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHI0(gI1)))*
      Conj(CpUChiChiISHI0PL(gO2,gI2,gI1))*CpUChiChiISHI0PR(gO1,gI2,gI1)*MChiI(gI2)
      ));
   result += SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(
      gI2)))*Conj(CpChiUChiAhPL(gI1,gO2,gI2))*CpChiUChiAhPR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*Conj
      (CpUChiFdconjSdPL(gO2,gI2,gI1))*CpUChiFdconjSdPR(gO1,gI2,gI1)*MFd(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),Sqr(MSDX(gI1)))*
      Conj(CpUChiFDXconjSDXPL(gO2,gI2,gI1))*CpUChiFDXconjSDXPR(gO1,gI2,gI1)*MFDX(
      gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*Conj(
      CpUChiFeconjSePL(gO2,gI2,gI1))*CpUChiFeconjSePR(gO1,gI2,gI1)*MFe(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*Conj
      (CpUChiFuconjSuPL(gO2,gI2,gI1))*CpUChiFuconjSuPR(gO1,gI2,gI1)*MFu(gI2)));
   result += MChaP*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaP),Sqr(MSHpp(gI2)))*Conj(
      CpbarChaPUChiSHppPL(gO2,gI2))*CpbarChaPUChiSHppPR(gO1,gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MVWm))*Conj(
      CpUChiChaconjVWmPR(gO2,gI2))*CpUChiChaconjVWmPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZp))*Conj(
      CpChiUChiVZpPL(gI2,gO2))*CpChiUChiVZpPR(gI2,gO1)*MChi(gI2));
   result += -4*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(CpChiUChiVZPL(
      gI2,gO2))*CpChiUChiVZPR(gI2,gO1)*MChi(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Chi_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = 0; k < 6; k++)
         self_energy(i, k) = self_energy_Chi_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -SUM(gI1,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPR(gI1,gO2))*CpbarChaUChiVWmPR(gI1,gO1));
   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MChaP),Sqr(MSHpp(gI1)))*Conj(
      CpUChiChaPconjSHppPR(gO2,gI1))*CpUChiChaPconjSHppPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(gI2)))
      *Conj(CpbarChaUChiHpmPR(gI1,gO2,gI2))*CpbarChaUChiHpmPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpUChiChaconjHpmPR(gO2,gI2,gI1))*CpUChiChaconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHp0(gI1)
      ))*Conj(CpUChiChiPconjSHp0PR(gO2,gI2,gI1))*CpUChiChiPconjSHp0PR(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHp0(gI1)
      ))*Conj(CpUChiChiPSHp0PR(gO2,gI2,gI1))*CpUChiChiPSHp0PR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpUChiFSIconjSSI0PR(gO2,gI2,gI1))*CpUChiFSIconjSSI0PR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpUChiFSISSI0PR(gO2,gI2,gI1))*CpUChiFSISSI0PR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChaI(gI1)),Sqr(MSHIp(gI2)
      ))*Conj(CpbarChaIUChiSHIpPR(gI1,gO2,gI2))*CpbarChaIUChiSHIpPR(gI1,gO1,gI2)))
      ;
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpbarFvUChiSvPR(gI1,gO2,gI2))*CpbarFvUChiSvPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*
      Conj(CpUChiFvconjSvPR(gO2,gI2,gI1))*CpUChiFvconjSvPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarFdUChiSdPR(gI1,gO2,gI2))*CpbarFdUChiSdPR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFDX(gI1)),Sqr(MSDX(gI2)))
      *Conj(CpbarFDXUChiSDXPR(gI1,gO2,gI2))*CpbarFDXUChiSDXPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarFeUChiSePR(gI1,gO2,gI2))*CpbarFeUChiSePR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpbarFuUChiSuPR(gI1,gO2,gI2))*CpbarFuUChiSuPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpChiUChihhPR(gI2,gO2,gI1))*CpChiUChihhPR(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MSHIp(gI1)
      ))*Conj(CpUChiChaIconjSHIpPR(gO2,gI2,gI1))*CpUChiChaIconjSHIpPR(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHI0(gI1)
      ))*Conj(CpUChiChiIconjSHI0PR(gO2,gI2,gI1))*CpUChiChiIconjSHI0PR(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHI0(gI1)
      ))*Conj(CpUChiChiISHI0PR(gO2,gI2,gI1))*CpUChiChiISHI0PR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)))*
      Conj(CpChiUChiAhPR(gI1,gO2,gI2))*CpChiUChiAhPR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpUChiFdconjSdPR(gO2,gI2,gI1))*CpUChiFdconjSdPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(MSDX(gI1)))
      *Conj(CpUChiFDXconjSDXPR(gO2,gI2,gI1))*CpUChiFDXconjSDXPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpUChiFeconjSePR(gO2,gI2,gI1))*CpUChiFeconjSePR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpUChiFuconjSuPR(gO2,gI2,gI1))*CpUChiFuconjSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaP),Sqr(MSHpp(gI2)))*Conj(
      CpbarChaPUChiSHppPR(gO2,gI2))*CpbarChaPUChiSHppPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVWm))*Conj(
      CpUChiChaconjVWmPL(gO2,gI2))*CpUChiChaconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZp))*Conj(CpChiUChiVZpPR(
      gI2,gO2))*CpChiUChiVZpPR(gI2,gO1));
   result += -SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(CpChiUChiVZPR(
      gI2,gO2))*CpChiUChiVZPR(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Chi_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = 0; k < 6; k++)
         self_energy(i, k) = self_energy_Chi_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -SUM(gI1,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPL(gI1,gO2))*CpbarChaUChiVWmPL(gI1,gO1));
   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MChaP),Sqr(MSHpp(gI1)))*Conj(
      CpUChiChaPconjSHppPL(gO2,gI1))*CpUChiChaPconjSHppPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(gI2)))
      *Conj(CpbarChaUChiHpmPL(gI1,gO2,gI2))*CpbarChaUChiHpmPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpUChiChaconjHpmPL(gO2,gI2,gI1))*CpUChiChaconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHp0(gI1)
      ))*Conj(CpUChiChiPconjSHp0PL(gO2,gI2,gI1))*CpUChiChiPconjSHp0PL(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHp0(gI1)
      ))*Conj(CpUChiChiPSHp0PL(gO2,gI2,gI1))*CpUChiChiPSHp0PL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpUChiFSIconjSSI0PL(gO2,gI2,gI1))*CpUChiFSIconjSSI0PL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpUChiFSISSI0PL(gO2,gI2,gI1))*CpUChiFSISSI0PL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChaI(gI1)),Sqr(MSHIp(gI2)
      ))*Conj(CpbarChaIUChiSHIpPL(gI1,gO2,gI2))*CpbarChaIUChiSHIpPL(gI1,gO1,gI2)))
      ;
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpbarFvUChiSvPL(gI1,gO2,gI2))*CpbarFvUChiSvPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*
      Conj(CpUChiFvconjSvPL(gO2,gI2,gI1))*CpUChiFvconjSvPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarFdUChiSdPL(gI1,gO2,gI2))*CpbarFdUChiSdPL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFDX(gI1)),Sqr(MSDX(gI2)))
      *Conj(CpbarFDXUChiSDXPL(gI1,gO2,gI2))*CpbarFDXUChiSDXPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarFeUChiSePL(gI1,gO2,gI2))*CpbarFeUChiSePL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpbarFuUChiSuPL(gI1,gO2,gI2))*CpbarFuUChiSuPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpChiUChihhPL(gI2,gO2,gI1))*CpChiUChihhPL(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MSHIp(gI1)
      ))*Conj(CpUChiChaIconjSHIpPL(gO2,gI2,gI1))*CpUChiChaIconjSHIpPL(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHI0(gI1)
      ))*Conj(CpUChiChiIconjSHI0PL(gO2,gI2,gI1))*CpUChiChiIconjSHI0PL(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHI0(gI1)
      ))*Conj(CpUChiChiISHI0PL(gO2,gI2,gI1))*CpUChiChiISHI0PL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)))*
      Conj(CpChiUChiAhPL(gI1,gO2,gI2))*CpChiUChiAhPL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpUChiFdconjSdPL(gO2,gI2,gI1))*CpUChiFdconjSdPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(MSDX(gI1)))
      *Conj(CpUChiFDXconjSDXPL(gO2,gI2,gI1))*CpUChiFDXconjSDXPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpUChiFeconjSePL(gO2,gI2,gI1))*CpUChiFeconjSePL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpUChiFuconjSuPL(gO2,gI2,gI1))*CpUChiFuconjSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaP),Sqr(MSHpp(gI2)))*Conj(
      CpbarChaPUChiSHppPL(gO2,gI2))*CpbarChaPUChiSHppPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVWm))*Conj(
      CpUChiChaconjVWmPR(gO2,gI2))*CpUChiChaconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(CpChiUChiVZPL(
      gI2,gO2))*CpChiUChiVZPL(gI2,gO1));
   result += -SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZp))*Conj(CpChiUChiVZpPL(
      gI2,gO2))*CpChiUChiVZpPL(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Chi_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = 0; k < 6; k++)
         self_energy(i, k) = self_energy_Chi_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += MChaP*SUM(gI1,0,1,B0(Sqr(p),Sqr(MChaP),Sqr(MSHp0(gI1)))*Conj(
      CpbarUChaChaPconjSHp0PL(gO2,gI1))*CpbarUChaChaPconjSHp0PR(gO1,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHpp(gI1)))*
      Conj(CpbarUChaChiPSHppPL(gO2,gI2,gI1))*CpbarUChaChiPSHppPR(gO1,gI2,gI1)*
      MChiP(gI2)));
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUChaChaAhPL(gO2,gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))*Conj
      (CpbarUChaChiHpmPL(gO2,gI2,gI1))*CpbarUChaChiHpmPR(gO1,gI2,gI1)*MChi(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUChaChahhPL(gO2,gI2,gI1))*CpbarUChaChahhPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUChaFeconjSvPL(gO2,gI2,gI1))*CpbarUChaFeconjSvPR(gO1,gI2,gI1)*MFe(gI2))
      );
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarUChabarFuSdPL(gO2,gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2)
      ));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarUChabarFvSePL(gO2,gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),Sqr(MSHI0(gI1)))*
      Conj(CpbarUChaChaIconjSHI0PL(gO2,gI2,gI1))*CpbarUChaChaIconjSHI0PR(gO1,gI2,
      gI1)*MChaI(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHIp(gI1)))*
      Conj(CpbarUChaChiISHIpPL(gO2,gI2,gI1))*CpbarUChaChiISHIpPR(gO1,gI2,gI1)*
      MChiI(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*Conj
      (CpbarUChaFdconjSuPL(gO2,gI2,gI1))*CpbarUChaFdconjSuPR(gO1,gI2,gI1)*MFd(gI2)
      ));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),0)*Conj(CpbarUChaChaVPPR(gO2,
      gI2))*CpbarUChaChaVPPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZ))*Conj(
      CpbarUChaChaVZPR(gO2,gI2))*CpbarUChaChaVZPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZp))*Conj(
      CpbarUChaChaVZpPR(gO2,gI2))*CpbarUChaChaVZpPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaChiVWmPR(gO2,gI2))*CpbarUChaChiVWmPL(gO1,gI2)*MChi(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MChaP),Sqr(MSHp0(gI1)))*Conj(
      CpbarUChaChaPconjSHp0PR(gO2,gI1))*CpbarUChaChaPconjSHp0PR(gO1,gI1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHpp(gI1)
      ))*Conj(CpbarUChaChiPSHppPR(gO2,gI2,gI1))*CpbarUChaChiPSHppPR(gO1,gI2,gI1)))
      ;
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUChaChaAhPR(gO2,gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpbarUChaChiHpmPR(gO2,gI2,gI1))*CpbarUChaChiHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChaChahhPR(gO2,gI2,gI1))*CpbarUChaChahhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUChaFeconjSvPR(gO2,gI2,gI1))*CpbarUChaFeconjSvPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarUChabarFuSdPR(gO2,gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarUChabarFvSePR(gO2,gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MSHI0(gI1)
      ))*Conj(CpbarUChaChaIconjSHI0PR(gO2,gI2,gI1))*CpbarUChaChaIconjSHI0PR(gO1,
      gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHIp(gI1)
      ))*Conj(CpbarUChaChiISHIpPR(gO2,gI2,gI1))*CpbarUChaChiISHIpPR(gO1,gI2,gI1)))
      ;
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChaFdconjSuPR(gO2,gI2,gI1))*CpbarUChaFdconjSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),0)*Conj(CpbarUChaChaVPPL(gO2,
      gI2))*CpbarUChaChaVPPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZ))*Conj(CpbarUChaChaVZPL
      (gO2,gI2))*CpbarUChaChaVZPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZp))*Conj(
      CpbarUChaChaVZpPL(gO2,gI2))*CpbarUChaChaVZpPL(gO1,gI2));
   result += -SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaChiVWmPL(gO2,gI2))*CpbarUChaChiVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MChaP),Sqr(MSHp0(gI1)))*Conj(
      CpbarUChaChaPconjSHp0PL(gO2,gI1))*CpbarUChaChaPconjSHp0PL(gO1,gI1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MSHpp(gI1)
      ))*Conj(CpbarUChaChiPSHppPL(gO2,gI2,gI1))*CpbarUChaChiPSHppPL(gO1,gI2,gI1)))
      ;
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUChaChaAhPL(gO2,gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpbarUChaChiHpmPL(gO2,gI2,gI1))*CpbarUChaChiHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChaChahhPL(gO2,gI2,gI1))*CpbarUChaChahhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUChaFeconjSvPL(gO2,gI2,gI1))*CpbarUChaFeconjSvPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarUChabarFuSdPL(gO2,gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarUChabarFvSePL(gO2,gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MSHI0(gI1)
      ))*Conj(CpbarUChaChaIconjSHI0PL(gO2,gI2,gI1))*CpbarUChaChaIconjSHI0PL(gO1,
      gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MSHIp(gI1)
      ))*Conj(CpbarUChaChiISHIpPL(gO2,gI2,gI1))*CpbarUChaChiISHIpPL(gO1,gI2,gI1)))
      ;
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChaFdconjSuPL(gO2,gI2,gI1))*CpbarUChaFdconjSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),0)*Conj(CpbarUChaChaVPPR(gO2,
      gI2))*CpbarUChaChaVPPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZp))*Conj(
      CpbarUChaChaVZpPR(gO2,gI2))*CpbarUChaChaVZpPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZ))*Conj(CpbarUChaChaVZPR
      (gO2,gI2))*CpbarUChaChaVZPR(gO1,gI2));
   result += -SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaChiVWmPR(gO2,gI2))*CpbarUChaChiVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFeFvHpmPL(gO2,gI2,gI1))*CpbarUFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUFeChaSvPL(gO2,gI2,gI1))*CpbarUFeChaSvPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*Conj(
      CpbarUFeChiSePL(gO2,gI2,gI1))*CpbarUFeChiSePR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPR(gO2,gI2
      ))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPR(
      gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZp))*Conj(
      CpbarUFeFeVZpPR(gO2,gI2))*CpbarUFeFeVZpPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(
      CpbarUFeFvVWmPR(gO2,gI2))*CpbarUFeFvVWmPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFeFvHpmPR(gO2,gI2,gI1))*CpbarUFeFvHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFeChaSvPR(gO2,gI2,gI1))*CpbarUFeChaSvPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPR(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPR(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUFeChiSePR(gO2,gI2,gI1))*CpbarUFeChiSePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPL(gO2,gI2))
      *CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPL(
      gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZp))*Conj(CpbarUFeFeVZpPL(
      gO2,gI2))*CpbarUFeFeVZpPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(CpbarUFeFvVWmPL(
      gO2,gI2))*CpbarUFeFvVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFeFvHpmPL(gO2,gI2,gI1))*CpbarUFeFvHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFeChaSvPL(gO2,gI2,gI1))*CpbarUFeChaSvPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUFeChiSePL(gO2,gI2,gI1))*CpbarUFeChiSePL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPR(gO2,gI2))
      *CpbarUFeFeVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZp))*Conj(CpbarUFeFeVZpPR(
      gO2,gI2))*CpbarUFeFeVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPR(
      gO2,gI2))*CpbarUFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(CpbarUFeFvVWmPR(
      gO2,gI2))*CpbarUFeFvVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFdFuHpmPL(gO2,gI2,gI1))*CpbarUFdFuHpmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1))
      )*Conj(CpbarUFdGluSdPL(gO2,gI1))*CpbarUFdGluSdPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarUFdChaSuPL(gO2,gI2,gI1))*CpbarUFdChaSuPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*Conj(
      CpbarUFdChiSdPL(gO2,gI2,gI1))*CpbarUFdChiSdPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPR(gO2,gI2
      ))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPR(
      gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZp))*Conj(
      CpbarUFdFdVZpPR(gO2,gI2))*CpbarUFdFdVZpPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarUFdFuVWmPR(gO2,gI2))*CpbarUFdFuVWmPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFdFuHpmPR(gO2,gI2,gI1))*CpbarUFdFuHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFdFdAhPR(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPR(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarUFdGluSdPR(gO2,gI1))*CpbarUFdGluSdPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFdChaSuPR(gO2,gI2,gI1))*CpbarUFdChaSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFdChiSdPR(gO2,gI2,gI1))*CpbarUFdChiSdPR(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPL(gO2,gI2))
      *CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPL(
      gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZp))*Conj(CpbarUFdFdVZpPL(
      gO2,gI2))*CpbarUFdFdVZpPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(CpbarUFdFuVWmPL(
      gO2,gI2))*CpbarUFdFuVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFdFuHpmPL(gO2,gI2,gI1))*CpbarUFdFuHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarUFdGluSdPL(gO2,gI1))*CpbarUFdGluSdPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFdChaSuPL(gO2,gI2,gI1))*CpbarUFdChaSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFdChiSdPL(gO2,gI2,gI1))*CpbarUFdChiSdPL(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPR(gO2,gI2))
      *CpbarUFdFdVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZp))*Conj(CpbarUFdFdVZpPR(
      gO2,gI2))*CpbarUFdFdVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPR(
      gO2,gI2))*CpbarUFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(CpbarUFdFuVWmPR(
      gO2,gI2))*CpbarUFdFuVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)*MFd(gI2))
      );
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)
      ));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1))
      )*Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2
      ))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(
      CpbarUFuFuVZpPR(gO2,gI2))*CpbarUFuFuVZpPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPR(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPR(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPR(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPR(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPL(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPL(gO2,gI2))
      *CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPL(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarUFuFuVZpPL(
      gO2,gI2))*CpbarUFuFuVZpPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2))
      *CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarUFuFuVZpPR(
      gO2,gI2))*CpbarUFuFuVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR(
      gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_FDX_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,MFDX(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUFDXFDXAhPL(gO2,gI1,gI2))*CpbarUFDXFDXAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFDXFDXhhPL(gO2,gI2,gI1))*CpbarUFDXFDXhhPR(gO1,gI2,gI1)*MFDX(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSDX(gI1)
      ))*Conj(CpbarUFDXGluSDXPL(gO2,gI1))*CpbarUFDXGluSDXPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSDX(gI1)))*Conj
      (CpbarUFDXChiSDXPL(gO2,gI2,gI1))*CpbarUFDXChiSDXPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),0)*Conj(
      CpbarUFDXFDXVGPR(gO2,gI2))*CpbarUFDXFDXVGPL(gO1,gI2)*MFDX(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),0)*Conj(CpbarUFDXFDXVPPR(gO2,
      gI2))*CpbarUFDXFDXVPPL(gO1,gI2)*MFDX(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),Sqr(MVZ))*Conj(
      CpbarUFDXFDXVZPR(gO2,gI2))*CpbarUFDXFDXVZPL(gO1,gI2)*MFDX(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),Sqr(MVZp))*Conj(
      CpbarUFDXFDXVZpPR(gO2,gI2))*CpbarUFDXFDXVZpPL(gO1,gI2)*MFDX(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_FDX_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_FDX_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_FDX_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFDXFDXAhPR(gO2,gI1,gI2))*CpbarUFDXFDXAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFDXFDXhhPR(gO2,gI2,gI1))*CpbarUFDXFDXhhPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSDX(gI1)))*
      Conj(CpbarUFDXGluSDXPR(gO2,gI1))*CpbarUFDXGluSDXPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSDX(gI1)))
      *Conj(CpbarUFDXChiSDXPR(gO2,gI2,gI1))*CpbarUFDXChiSDXPR(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),0)*Conj(
      CpbarUFDXFDXVGPL(gO2,gI2))*CpbarUFDXFDXVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),0)*Conj(CpbarUFDXFDXVPPL(gO2,
      gI2))*CpbarUFDXFDXVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(MVZ))*Conj(CpbarUFDXFDXVZPL
      (gO2,gI2))*CpbarUFDXFDXVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(MVZp))*Conj(
      CpbarUFDXFDXVZpPL(gO2,gI2))*CpbarUFDXFDXVZpPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_FDX_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_FDX_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_FDX_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFDXFDXAhPL(gO2,gI1,gI2))*CpbarUFDXFDXAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFDXFDXhhPL(gO2,gI2,gI1))*CpbarUFDXFDXhhPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSDX(gI1)))*
      Conj(CpbarUFDXGluSDXPL(gO2,gI1))*CpbarUFDXGluSDXPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSDX(gI1)))
      *Conj(CpbarUFDXChiSDXPL(gO2,gI2,gI1))*CpbarUFDXChiSDXPL(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),0)*Conj(
      CpbarUFDXFDXVGPR(gO2,gI2))*CpbarUFDXFDXVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),0)*Conj(CpbarUFDXFDXVPPR(gO2,
      gI2))*CpbarUFDXFDXVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(MVZp))*Conj(
      CpbarUFDXFDXVZpPR(gO2,gI2))*CpbarUFDXFDXVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFDX(gI2)),Sqr(MVZ))*Conj(CpbarUFDXFDXVZPR
      (gO2,gI2))*CpbarUFDXFDXVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_FDX_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_FDX_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChaI_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,MChaI(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUChaIChaIAhPL(gO2,gI1,gI2))*CpbarUChaIChaIAhPR(gO1,gI1,gI2)
      ));
   result += SUM(gI1,0,2,SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),Sqr(Mhh(gI1)))*Conj
      (CpbarUChaIChaIhhPL(gO2,gI2,gI1))*CpbarUChaIChaIhhPR(gO1,gI2,gI1)*MChaI(gI2)
      ));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHI0(gI1)))*
      Conj(CpbarUChaIChaSHI0PL(gO2,gI2,gI1))*CpbarUChaIChaSHI0PR(gO1,gI2,gI1)*MCha
      (gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHIp(gI1)))*
      Conj(CpbarUChaIChiSHIpPL(gO2,gI2,gI1))*CpbarUChaIChiSHIpPR(gO1,gI2,gI1)*MChi
      (gI2)));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),0)*Conj(CpbarUChaIChaIVPPR(
      gO2,gI2))*CpbarUChaIChaIVPPL(gO1,gI2)*MChaI(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVZ))*Conj(
      CpbarUChaIChaIVZPR(gO2,gI2))*CpbarUChaIChaIVZPL(gO1,gI2)*MChaI(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVZp))*Conj(
      CpbarUChaIChaIVZpPR(gO2,gI2))*CpbarUChaIChaIVZpPL(gO1,gI2)*MChaI(gI2));
   result += -4*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaIChiIVWmPR(gO2,gI2))*CpbarUChaIChiIVWmPL(gO1,gI2)*MChiI(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_ChaI_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_ChaI_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChaI_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChaI(gI1)),Sqr(MAh(gI2)))
      *Conj(CpbarUChaIChaIAhPR(gO2,gI1,gI2))*CpbarUChaIChaIAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(Mhh(gI1)))
      *Conj(CpbarUChaIChaIhhPR(gO2,gI2,gI1))*CpbarUChaIChaIhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHI0(gI1))
      )*Conj(CpbarUChaIChaSHI0PR(gO2,gI2,gI1))*CpbarUChaIChaSHI0PR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHIp(gI1))
      )*Conj(CpbarUChaIChiSHIpPR(gO2,gI2,gI1))*CpbarUChaIChiSHIpPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),0)*Conj(CpbarUChaIChaIVPPL(gO2
      ,gI2))*CpbarUChaIChaIVPPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVZ))*Conj(
      CpbarUChaIChaIVZPL(gO2,gI2))*CpbarUChaIChaIVZPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVZp))*Conj(
      CpbarUChaIChaIVZpPL(gO2,gI2))*CpbarUChaIChaIVZpPL(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaIChiIVWmPL(gO2,gI2))*CpbarUChaIChiIVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_ChaI_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_ChaI_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChaI_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChaI(gI1)),Sqr(MAh(gI2)))
      *Conj(CpbarUChaIChaIAhPL(gO2,gI1,gI2))*CpbarUChaIChaIAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(Mhh(gI1)))
      *Conj(CpbarUChaIChaIhhPL(gO2,gI2,gI1))*CpbarUChaIChaIhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHI0(gI1))
      )*Conj(CpbarUChaIChaSHI0PL(gO2,gI2,gI1))*CpbarUChaIChaSHI0PL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHIp(gI1))
      )*Conj(CpbarUChaIChiSHIpPL(gO2,gI2,gI1))*CpbarUChaIChiSHIpPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),0)*Conj(CpbarUChaIChaIVPPR(gO2
      ,gI2))*CpbarUChaIChaIVPPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVZp))*Conj(
      CpbarUChaIChaIVZpPR(gO2,gI2))*CpbarUChaIChaIVZpPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVZ))*Conj(
      CpbarUChaIChaIVZPR(gO2,gI2))*CpbarUChaIChaIVZPR(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaIChiIVWmPR(gO2,gI2))*CpbarUChaIChiIVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_ChaI_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_ChaI_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChiI_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -4*SUM(gI1,0,1,B0(Sqr(p),Sqr(MChaI(gI1)),Sqr(MVWm))*Conj(
      CpbarChaIUChiIVWmPL(gI1,gO2))*CpbarChaIUChiIVWmPR(gI1,gO1)*MChaI(gI1));
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSHIp(
      gI2)))*Conj(CpbarChaUChiISHIpPL(gI1,gO2,gI2))*CpbarChaUChiISHIpPR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(Mhh(gI1)))*Conj
      (CpChiIUChiIhhPL(gI2,gO2,gI1))*CpChiIUChiIhhPR(gI2,gO1,gI1)*MChiI(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHIp(gI1)))*
      Conj(CpUChiIChaconjSHIpPL(gO2,gI2,gI1))*CpUChiIChaconjSHIpPR(gO1,gI2,gI1)*
      MCha(gI2)));
   result += SUM(gI1,0,3,MChiI(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChiI(gI1)),Sqr(MAh(
      gI2)))*Conj(CpChiIUChiIAhPL(gI1,gO2,gI2))*CpChiIUChiIAhPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHI0(gI1)))*
      Conj(CpChiUChiIconjSHI0PL(gI2,gO2,gI1))*CpChiUChiIconjSHI0PR(gI2,gO1,gI1)*
      MChi(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHI0(gI1)))*
      Conj(CpChiUChiISHI0PL(gI2,gO2,gI1))*CpChiUChiISHI0PR(gI2,gO1,gI1)*MChi(gI2))
      );
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVWm))*Conj(
      CpUChiIChaIconjVWmPR(gO2,gI2))*CpUChiIChaIconjVWmPL(gO1,gI2)*MChaI(gI2));
   result += -4*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVZp))*Conj(
      CpChiIUChiIVZpPL(gI2,gO2))*CpChiIUChiIVZpPR(gI2,gO1)*MChiI(gI2));
   result += -4*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVZ))*Conj(
      CpChiIUChiIVZPL(gI2,gO2))*CpChiIUChiIVZPR(gI2,gO1)*MChiI(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_ChiI_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++)
         self_energy(i, k) = self_energy_ChiI_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChiI_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -SUM(gI1,0,1,B1(Sqr(p),Sqr(MChaI(gI1)),Sqr(MVWm))*Conj(
      CpbarChaIUChiIVWmPR(gI1,gO2))*CpbarChaIUChiIVWmPR(gI1,gO1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSHIp(gI2))
      )*Conj(CpbarChaUChiISHIpPR(gI1,gO2,gI2))*CpbarChaUChiISHIpPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(Mhh(gI1)))
      *Conj(CpChiIUChiIhhPR(gI2,gO2,gI1))*CpChiIUChiIhhPR(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHIp(gI1))
      )*Conj(CpUChiIChaconjSHIpPR(gO2,gI2,gI1))*CpUChiIChaconjSHIpPR(gO1,gI2,gI1))
      );
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChiI(gI1)),Sqr(MAh(gI2)))
      *Conj(CpChiIUChiIAhPR(gI1,gO2,gI2))*CpChiIUChiIAhPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHI0(gI1))
      )*Conj(CpChiUChiIconjSHI0PR(gI2,gO2,gI1))*CpChiUChiIconjSHI0PR(gI2,gO1,gI1))
      );
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHI0(gI1))
      )*Conj(CpChiUChiISHI0PR(gI2,gO2,gI1))*CpChiUChiISHI0PR(gI2,gO1,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVWm))*Conj(
      CpUChiIChaIconjVWmPL(gO2,gI2))*CpUChiIChaIconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVZp))*Conj(
      CpChiIUChiIVZpPR(gI2,gO2))*CpChiIUChiIVZpPR(gI2,gO1));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVZ))*Conj(CpChiIUChiIVZPR
      (gI2,gO2))*CpChiIUChiIVZPR(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_ChiI_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++)
         self_energy(i, k) = self_energy_ChiI_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChiI_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -SUM(gI1,0,1,B1(Sqr(p),Sqr(MChaI(gI1)),Sqr(MVWm))*Conj(
      CpbarChaIUChiIVWmPL(gI1,gO2))*CpbarChaIUChiIVWmPL(gI1,gO1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSHIp(gI2))
      )*Conj(CpbarChaUChiISHIpPL(gI1,gO2,gI2))*CpbarChaUChiISHIpPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(Mhh(gI1)))
      *Conj(CpChiIUChiIhhPL(gI2,gO2,gI1))*CpChiIUChiIhhPL(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHIp(gI1))
      )*Conj(CpUChiIChaconjSHIpPL(gO2,gI2,gI1))*CpUChiIChaconjSHIpPL(gO1,gI2,gI1))
      );
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChiI(gI1)),Sqr(MAh(gI2)))
      *Conj(CpChiIUChiIAhPL(gI1,gO2,gI2))*CpChiIUChiIAhPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHI0(gI1))
      )*Conj(CpChiUChiIconjSHI0PL(gI2,gO2,gI1))*CpChiUChiIconjSHI0PL(gI2,gO1,gI1))
      );
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHI0(gI1))
      )*Conj(CpChiUChiISHI0PL(gI2,gO2,gI1))*CpChiUChiISHI0PL(gI2,gO1,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChaI(gI2)),Sqr(MVWm))*Conj(
      CpUChiIChaIconjVWmPR(gO2,gI2))*CpUChiIChaIconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVZ))*Conj(CpChiIUChiIVZPL
      (gI2,gO2))*CpChiIUChiIVZPL(gI2,gO1));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChiI(gI2)),Sqr(MVZp))*Conj(
      CpChiIUChiIVZpPL(gI2,gO2))*CpChiIUChiIVZpPL(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_ChiI_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++)
         self_energy(i, k) = self_energy_ChiI_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_FSI_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSSI0(gI1)))*
      Conj(CpChiUFSIconjSSI0PL(gI2,gO2,gI1))*CpChiUFSIconjSSI0PR(gI2,gO1,gI1)*MChi
      (gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSSI0(gI1)))*
      Conj(CpChiUFSISSI0PL(gI2,gO2,gI1))*CpChiUFSISSI0PR(gI2,gO1,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFSI(gI2)),Sqr(MVZp))*Conj(
      CpFSIUFSIVZpPL(gI2,gO2))*CpFSIUFSIVZpPR(gI2,gO1)*MFSI(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFSI(gI2)),Sqr(MVZ))*Conj(CpFSIUFSIVZPL(
      gI2,gO2))*CpFSIUFSIVZPR(gI2,gO1)*MFSI(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_FSI_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_FSI_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_FSI_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpChiUFSIconjSSI0PR(gI2,gO2,gI1))*CpChiUFSIconjSSI0PR(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpChiUFSISSI0PR(gI2,gO2,gI1))*CpChiUFSISSI0PR(gI2,gO1,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MVZp))*Conj(CpFSIUFSIVZpPR(
      gI2,gO2))*CpFSIUFSIVZpPR(gI2,gO1));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MVZ))*Conj(CpFSIUFSIVZPR(
      gI2,gO2))*CpFSIUFSIVZPR(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_FSI_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_FSI_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_FSI_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpChiUFSIconjSSI0PL(gI2,gO2,gI1))*CpChiUFSIconjSSI0PL(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSSI0(gI1))
      )*Conj(CpChiUFSISSI0PL(gI2,gO2,gI1))*CpChiUFSISSI0PL(gI2,gO1,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MVZ))*Conj(CpFSIUFSIVZPL(
      gI2,gO2))*CpFSIUFSIVZPL(gI2,gO1));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MFSI(gI2)),Sqr(MVZp))*Conj(CpFSIUFSIVZpPL(
      gI2,gO2))*CpFSIUFSIVZpPL(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_FSI_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_FSI_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChiP_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -4*MChaP*B0(Sqr(p),Sqr(MChaP),Sqr(MVWm))*Conj(CpbarChaPUChiPVWmPL(gO2
      ))*CpbarChaPUChiPVWmPR(gO1);
   result += -4*MChaP*B0(Sqr(p),Sqr(MChaP),Sqr(MVWm))*Conj(CpUChiPChaPconjVWmPR(
      gO2))*CpUChiPChaPconjVWmPL(gO1);
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSHpp(
      gI2)))*Conj(CpbarChaUChiPSHppPL(gI1,gO2,gI2))*CpbarChaUChiPSHppPR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHpp(gI1)))*
      Conj(CpUChiPChaconjSHppPL(gO2,gI2,gI1))*CpUChiPChaconjSHppPR(gO1,gI2,gI1)*
      MCha(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHp0(gI1)))*
      Conj(CpChiUChiPconjSHp0PL(gI2,gO2,gI1))*CpChiUChiPconjSHp0PR(gI2,gO1,gI1)*
      MChi(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHp0(gI1)))*
      Conj(CpChiUChiPSHp0PL(gI2,gO2,gI1))*CpChiUChiPSHp0PR(gI2,gO1,gI1)*MChi(gI2))
      );
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVZp))*Conj(
      CpChiPUChiPVZpPL(gI2,gO2))*CpChiPUChiPVZpPR(gI2,gO1)*MChiP(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVZ))*Conj(
      CpChiPUChiPVZPL(gI2,gO2))*CpChiPUChiPVZPR(gI2,gO1)*MChiP(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_ChiP_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_ChiP_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChiP_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B1(Sqr(p),Sqr(MChaP),Sqr(MVWm))*Conj(CpbarChaPUChiPVWmPR(gO2))*
      CpbarChaPUChiPVWmPR(gO1));
   result += -(B1(Sqr(p),Sqr(MChaP),Sqr(MVWm))*Conj(CpUChiPChaPconjVWmPL(gO2))*
      CpUChiPChaPconjVWmPL(gO1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSHpp(gI2))
      )*Conj(CpbarChaUChiPSHppPR(gI1,gO2,gI2))*CpbarChaUChiPSHppPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHpp(gI1))
      )*Conj(CpUChiPChaconjSHppPR(gO2,gI2,gI1))*CpUChiPChaconjSHppPR(gO1,gI2,gI1))
      );
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHp0(gI1))
      )*Conj(CpChiUChiPconjSHp0PR(gI2,gO2,gI1))*CpChiUChiPconjSHp0PR(gI2,gO1,gI1))
      );
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHp0(gI1))
      )*Conj(CpChiUChiPSHp0PR(gI2,gO2,gI1))*CpChiUChiPSHp0PR(gI2,gO1,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVZp))*Conj(
      CpChiPUChiPVZpPR(gI2,gO2))*CpChiPUChiPVZpPR(gI2,gO1));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVZ))*Conj(CpChiPUChiPVZPR
      (gI2,gO2))*CpChiPUChiPVZPR(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_ChiP_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_ChiP_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChiP_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B1(Sqr(p),Sqr(MChaP),Sqr(MVWm))*Conj(CpbarChaPUChiPVWmPL(gO2))*
      CpbarChaPUChiPVWmPL(gO1));
   result += -(B1(Sqr(p),Sqr(MChaP),Sqr(MVWm))*Conj(CpUChiPChaPconjVWmPR(gO2))*
      CpUChiPChaPconjVWmPR(gO1));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSHpp(gI2))
      )*Conj(CpbarChaUChiPSHppPL(gI1,gO2,gI2))*CpbarChaUChiPSHppPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHpp(gI1))
      )*Conj(CpUChiPChaconjSHppPL(gO2,gI2,gI1))*CpUChiPChaconjSHppPL(gO1,gI2,gI1))
      );
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHp0(gI1))
      )*Conj(CpChiUChiPconjSHp0PL(gI2,gO2,gI1))*CpChiUChiPconjSHp0PL(gI2,gO1,gI1))
      );
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHp0(gI1))
      )*Conj(CpChiUChiPSHp0PL(gI2,gO2,gI1))*CpChiUChiPSHp0PL(gI2,gO1,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVZ))*Conj(CpChiPUChiPVZPL
      (gI2,gO2))*CpChiPUChiPVZPL(gI2,gO1));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVZp))*Conj(
      CpChiPUChiPVZpPL(gI2,gO2))*CpChiPUChiPVZpPL(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_ChiP_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_ChiP_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Glu_1loop_1(double p ) const
{
   std::complex<double> result;

   result += -12*MGlu*B0(Sqr(p),Sqr(MGlu),0)*Conj(CpGluGluVGPR())*CpGluGluVGPL();
   result += 0.5*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarFdGluSdPL(gI1,gI2))*CpbarFdGluSdPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,2,MFDX(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFDX(gI1)),Sqr(
      MSDX(gI2)))*Conj(CpbarFDXGluSDXPL(gI1,gI2))*CpbarFDXGluSDXPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(
      gI2)))*Conj(CpbarFuGluSuPL(gI1,gI2))*CpbarFuGluSuPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpGluFdconjSdPL(gI2,gI1))*CpGluFdconjSdPR(gI2,gI1)*MFd(gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFDX(gI2)),Sqr(MSDX(gI1)))*
      Conj(CpGluFDXconjSDXPL(gI2,gI1))*CpGluFDXconjSDXPR(gI2,gI1)*MFDX(gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpGluFuconjSuPL(gI2,gI1))*CpGluFuconjSuPR(gI2,gI1)*MFu(gI2)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Glu_1loop_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluGluVGPL())*B1(Sqr(p),Sqr(MGlu),0);
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFdGluSdPR(gI1,gI2))*B1(Sqr(
      p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))));
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFDXGluSDXPR(gI1,gI2))*B1(
      Sqr(p),Sqr(MFDX(gI1)),Sqr(MSDX(gI2)))));
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFuGluSuPR(gI1,gI2))*B1(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFdconjSdPR(gI2,gI1))*B1(Sqr
      (p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFDXconjSDXPR(gI2,gI1))*B1(
      Sqr(p),Sqr(MFDX(gI2)),Sqr(MSDX(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFuconjSuPR(gI2,gI1))*B1(Sqr
      (p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Glu_1loop_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluGluVGPR())*B1(Sqr(p),Sqr(MGlu),0);
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFdGluSdPL(gI1,gI2))*B1(Sqr(
      p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))));
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFDXGluSDXPL(gI1,gI2))*B1(
      Sqr(p),Sqr(MFDX(gI1)),Sqr(MSDX(gI2)))));
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFuGluSuPL(gI1,gI2))*B1(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFdconjSdPL(gI2,gI1))*B1(Sqr
      (p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFDXconjSDXPL(gI2,gI1))*B1(
      Sqr(p),Sqr(MFDX(gI2)),Sqr(MSDX(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFuconjSuPL(gI2,gI1))*B1(Sqr
      (p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFvFeconjHpmPL(gO2,gI2,gI1))*CpbarFvFeconjHpmPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(
      gI2)))*Conj(CpbarChabarFvSePL(gI1,gO2,gI2))*CpbarChabarFvSePR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarFvChiSvPL(gO2,gI2,gI1))*CpbarFvChiSvPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPR(
      gO2,gI2))*CpbarFvFvVZPL(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZp))*Conj(CpbarFvFvVZpPR
      (gO2,gI2))*CpbarFvFvVZpPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fv_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFvFeconjHpmPR(gO2,gI2,gI1))*CpbarFvFeconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarChabarFvSePR(gI1,gO2,gI2))*CpbarChabarFvSePR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFvChiSvPR(gO2,gI2,gI1))*CpbarFvChiSvPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPL(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPL(gO2
      ,gI2))*CpbarFvFvVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZp))*Conj(CpbarFvFvVZpPL(
      gO2,gI2))*CpbarFvFvVZpPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fv_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFvFeconjHpmPL(gO2,gI2,gI1))*CpbarFvFeconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarChabarFvSePL(gI1,gO2,gI2))*CpbarChabarFvSePL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFvChiSvPL(gO2,gI2,gI1))*CpbarFvChiSvPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZp))*Conj(CpbarFvFvVZpPR(
      gO2,gI2))*CpbarFvFvVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPR(gO2
      ,gI2))*CpbarFvFvVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fv_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_ChaP_1loop_1(double p ) const
{
   std::complex<double> result;

   result += -4*MChaP*B0(Sqr(p),Sqr(MChaP),0)*Conj(CpbarChaPChaPVPPR())*
      CpbarChaPChaPVPPL();
   result += -4*MChaP*B0(Sqr(p),Sqr(MChaP),Sqr(MVZ))*Conj(CpbarChaPChaPVZPR())*
      CpbarChaPChaPVZPL();
   result += -4*MChaP*B0(Sqr(p),Sqr(MChaP),Sqr(MVZp))*Conj(CpbarChaPChaPVZpPR())*
      CpbarChaPChaPVZpPL();
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSHp0(gI1)))*
      Conj(CpbarChaPChaSHp0PL(gI2,gI1))*CpbarChaPChaSHp0PR(gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSHpp(gI1)))*
      Conj(CpbarChaPChiSHppPL(gI2,gI1))*CpbarChaPChiSHppPR(gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChiP(gI2)),Sqr(MVWm))*Conj(
      CpbarChaPChiPVWmPR(gI2))*CpbarChaPChiPVWmPL(gI2)*MChiP(gI2));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_ChaP_1loop_PR(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpbarChaPChaPVPPL())*B1(Sqr(p),Sqr(MChaP),0));
   result += -(AbsSqr(CpbarChaPChaPVZPL())*B1(Sqr(p),Sqr(MChaP),Sqr(MVZ)));
   result += -(AbsSqr(CpbarChaPChaPVZpPL())*B1(Sqr(p),Sqr(MChaP),Sqr(MVZp)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpbarChaPChaSHp0PR(gI2,gI1))*B1(
      Sqr(p),Sqr(MCha(gI2)),Sqr(MSHp0(gI1)))));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,AbsSqr(CpbarChaPChiSHppPR(gI2,gI1))*B1(
      Sqr(p),Sqr(MChi(gI2)),Sqr(MSHpp(gI1)))));
   result += -SUM(gI2,0,1,AbsSqr(CpbarChaPChiPVWmPL(gI2))*B1(Sqr(p),Sqr(MChiP(gI2)
      ),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_ChaP_1loop_PL(double p ) const
{
   std::complex<double> result;

   result += -(AbsSqr(CpbarChaPChaPVPPR())*B1(Sqr(p),Sqr(MChaP),0));
   result += -(AbsSqr(CpbarChaPChaPVZpPR())*B1(Sqr(p),Sqr(MChaP),Sqr(MVZp)));
   result += -(AbsSqr(CpbarChaPChaPVZPR())*B1(Sqr(p),Sqr(MChaP),Sqr(MVZ)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpbarChaPChaSHp0PL(gI2,gI1))*B1(
      Sqr(p),Sqr(MCha(gI2)),Sqr(MSHp0(gI1)))));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,AbsSqr(CpbarChaPChiSHppPL(gI2,gI1))*B1(
      Sqr(p),Sqr(MChi(gI2)),Sqr(MSHpp(gI1)))));
   result += -SUM(gI2,0,1,AbsSqr(CpbarChaPChiPVWmPR(gI2))*B1(Sqr(p),Sqr(MChiP(gI2)
      ),Sqr(MVWm)));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarFeChaSvPL(gO2,gI2,gI1))*CpbarFeChaSvPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*Conj(
      CpbarFeChiSePL(gO2,gI2,gI1))*CpbarFeChiSePR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPR(
      gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZp))*Conj(CpbarFeFeVZpPR
      (gO2,gI2))*CpbarFeFeVZpPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(CpbarFeFvVWmPR
      (gO2,gI2))*CpbarFeFvVWmPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFeFvHpmPR(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFeChaSvPR(gO2,gI2,gI1))*CpbarFeChaSvPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPR(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPR(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarFeChiSePR(gO2,gI2,gI1))*CpbarFeChiSePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPL(gO2
      ,gI2))*CpbarFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZp))*Conj(CpbarFeFeVZpPL(
      gO2,gI2))*CpbarFeFeVZpPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(CpbarFeFvVWmPL(
      gO2,gI2))*CpbarFeFvVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFeChaSvPL(gO2,gI2,gI1))*CpbarFeChaSvPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarFeChiSePL(gO2,gI2,gI1))*CpbarFeChiSePL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZp))*Conj(CpbarFeFeVZpPR(
      gO2,gI2))*CpbarFeFeVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPR(gO2
      ,gI2))*CpbarFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(CpbarFeFvVWmPR(
      gO2,gI2))*CpbarFeFvVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fe_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFdFuHpmPL(gO2,gI2,gI1))*CpbarFdFuHpmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1))
      )*Conj(CpbarFdGluSdPL(gO2,gI1))*CpbarFdGluSdPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarFdChaSuPL(gO2,gI2,gI1))*CpbarFdChaSuPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*Conj(
      CpbarFdChiSdPL(gO2,gI2,gI1))*CpbarFdChiSdPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPR(
      gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZp))*Conj(CpbarFdFdVZpPR
      (gO2,gI2))*CpbarFdFdVZpPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(CpbarFdFuVWmPR
      (gO2,gI2))*CpbarFdFuVWmPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFdFuHpmPR(gO2,gI2,gI1))*CpbarFdFuHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFdFdAhPR(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPR(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarFdGluSdPR(gO2,gI1))*CpbarFdGluSdPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFdChaSuPR(gO2,gI2,gI1))*CpbarFdChaSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarFdChiSdPR(gO2,gI2,gI1))*CpbarFdChiSdPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPL(gO2
      ,gI2))*CpbarFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZp))*Conj(CpbarFdFdVZpPL(
      gO2,gI2))*CpbarFdFdVZpPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(CpbarFdFuVWmPL(
      gO2,gI2))*CpbarFdFuVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fd_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFdFuHpmPL(gO2,gI2,gI1))*CpbarFdFuHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarFdGluSdPL(gO2,gI1))*CpbarFdGluSdPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFdChaSuPL(gO2,gI2,gI1))*CpbarFdChaSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarFdChiSdPL(gO2,gI2,gI1))*CpbarFdChiSdPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZp))*Conj(CpbarFdFdVZpPR(
      gO2,gI2))*CpbarFdFdVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPR(gO2
      ,gI2))*CpbarFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(CpbarFdFuVWmPR(
      gO2,gI2))*CpbarFdFuVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fd_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFuFdconjHpmPL(gO2,gI2,gI1))*CpbarFuFdconjHpmPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarChabarFuSdPL(gI1,gO2,gI2))*CpbarChabarFuSdPR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1))
      )*Conj(CpbarFuGluSuPL(gO2,gI1))*CpbarFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarFuChiSuPL(gO2,gI2,gI1))*CpbarFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPR(gO2,gI2)
      )*CpbarFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPR(
      gO2,gI2))*CpbarFuFuVZPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarFuFuVZpPR
      (gO2,gI2))*CpbarFuFuVZpPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFuFdconjHpmPR(gO2,gI2,gI1))*CpbarFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarFuSdPR(gI1,gO2,gI2))*CpbarChabarFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFuFuAhPR(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPR(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarFuGluSuPR(gO2,gI1))*CpbarFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFuChiSuPR(gO2,gI2,gI1))*CpbarFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPL(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPL(gO2,gI2))*
      CpbarFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPL(gO2
      ,gI2))*CpbarFuFuVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarFuFuVZpPL(
      gO2,gI2))*CpbarFuFuVZpPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFuFdconjHpmPL(gO2,gI2,gI1))*CpbarFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarFuSdPL(gI1,gO2,gI2))*CpbarChabarFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarFuGluSuPL(gO2,gI1))*CpbarFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFuChiSuPL(gO2,gI2,gI1))*CpbarFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPR(gO2,gI2))*
      CpbarFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarFuFuVZpPR(
      gO2,gI2))*CpbarFuFuVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPR(gO2
      ,gI2))*CpbarFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL_heavy_rotated(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy_rotated(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_1_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)*MFd(gI2))
      );
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)
      ));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1))
      )*Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2
      ))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(
      CpbarUFuFuVZpPR(gO2,gI2))*CpbarUFuFuVZpPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_1_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PR_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPR(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPR(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPR(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPR(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPL(gO2,gI2))
      *CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPL(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarUFuFuVZpPL(
      gO2,gI2))*CpbarUFuFuVZpPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PR_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fu_1loop_PL_heavy(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2))
      *CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZp))*Conj(CpbarUFuFuVZpPR(
      gO2,gI2))*CpbarUFuFuVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR(
      gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Fu_1loop_PL_heavy(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::tadpole_hh_1loop(int gO1) const
{
   std::complex<double> result;

   result += A0(Sqr(MVWm))*CpbargWmCgWmCUhh(gO1);
   result += A0(Sqr(MVWm))*CpbargWmgWmUhh(gO1);
   result += A0(Sqr(MVZ))*CpbargZgZUhh(gO1);
   result += A0(Sqr(MVZp))*CpbargZpgZpUhh(gO1);
   result += 4*A0(Sqr(MVWm))*CpUhhconjVWmVWm(gO1);
   result += 2*A0(Sqr(MVZp))*CpUhhVZpVZp(gO1);
   result += 2*A0(Sqr(MVZ))*CpUhhVZVZ(gO1);
   result += -SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpUhhHpmconjHpm(gO1,gI1,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHp0(gI1)))*CpUhhSHp0conjSHp0(gO1,gI1,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSHpp(gI1)))*CpUhhSHppconjSHpp(gO1,gI1,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MSSI0(gI1)))*CpUhhSSI0conjSSI0(gO1,gI1,gI1));
   result += 2*SUM(gI1,0,1,A0(Sqr(MCha(gI1)))*(CpbarChaChaUhhPL(gI1,gI1,gO1) +
      CpbarChaChaUhhPR(gI1,gI1,gO1))*MCha(gI1));
   result += 2*SUM(gI1,0,1,A0(Sqr(MChaI(gI1)))*(CpbarChaIChaIUhhPL(gI1,gI1,gO1) +
      CpbarChaIChaIUhhPR(gI1,gI1,gO1))*MChaI(gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhh(gI1,gI1,gO1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhh(gI1,gI1,gO1));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUhhSvconjSv(gO1,gI1,gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdUhhPL(gI1,gI1,gO1) +
      CpbarFdFdUhhPR(gI1,gI1,gO1))*MFd(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFDX(gI1)))*(CpbarFDXFDXUhhPL(gI1,gI1,gO1) +
      CpbarFDXFDXUhhPR(gI1,gI1,gO1))*MFDX(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFeUhhPL(gI1,gI1,gO1) +
      CpbarFeFeUhhPR(gI1,gI1,gO1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuUhhPL(gI1,gI1,gO1) +
      CpbarFuFuUhhPR(gI1,gI1,gO1))*MFu(gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHI0(gI1)))*CpUhhSHI0conjSHI0(gO1,gI1,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MSHIp(gI1)))*CpUhhSHIpconjSHIp(gO1,gI1,gI1));
   result += SUM(gI1,0,3,A0(Sqr(MChiI(gI1)))*(CpChiIChiIUhhPL(gI1,gI1,gO1) +
      CpChiIChiIUhhPR(gI1,gI1,gO1))*MChiI(gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhSdconjSd(gO1,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSDX(gI1)))*CpUhhSDXconjSDX(gO1,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhSeconjSe(gO1,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhSuconjSu(gO1,gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MChi(gI1)))*(CpChiChiUhhPL(gI1,gI1,gO1) +
      CpChiChiUhhPR(gI1,gI1,gO1))*MChi(gI1));

   return result * oneOver16PiSqr;
}


void CLASSNAME::calculate_MSu_2nd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(mq2(1,1));
   sf_data.mr2 = Re(mu2(1,1));
   sf_data.yf  = Re(Yu(1,1));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYu(1,1));
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::up];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::up];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::up];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: stop tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSd_2nd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(mq2(1,1));
   sf_data.mr2 = Re(md2(1,1));
   sf_data.yf  = Re(Yd(1,1));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYd(1,1));
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::down];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::down];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::down];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: sbottom tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSv_2nd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(ml2(1,1));
   sf_data.mr2 = 0.;
   sf_data.yf  = 0.;
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = 0.;
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::neutrino];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::neutrino];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::neutrino];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: sneutrino tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSe_2nd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(ml2(1,1));
   sf_data.mr2 = Re(me2(1,1));
   sf_data.yf  = Re(Ye(1,1));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYe(1,1));
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::electron];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::electron];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::electron];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: selecton tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSu_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(mq2(2,2));
   sf_data.mr2 = Re(mu2(2,2));
   sf_data.yf  = Re(Yu(2,2));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYu(2,2));
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::up];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::up];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::up];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: stop tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSd_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(mq2(2,2));
   sf_data.mr2 = Re(md2(2,2));
   sf_data.yf  = Re(Yd(2,2));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYd(2,2));
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::down];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::down];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::down];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: sbottom tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSv_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(ml2(2,2));
   sf_data.mr2 = 0.;
   sf_data.yf  = 0.;
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = 0.;
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::neutrino];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::neutrino];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::neutrino];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: sneutrino tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}

void CLASSNAME::calculate_MSe_3rd_generation(double& msf1, double& msf2, double& theta) const {
   sfermions::Mass_data sf_data;
   sf_data.ml2 = Re(ml2(2,2));
   sf_data.mr2 = Re(me2(2,2));
   sf_data.yf  = Re(Ye(2,2));
   sf_data.vd  = Re(vd);
   sf_data.vu  = Re(vu);
   sf_data.gY  = 0.7745966692414834*g1;
   sf_data.g2  = g2;
   sf_data.Tyf = Re(TYe(2,2));
   sf_data.mu  = Re(0.7071067811865475*vs*Lambdax);
   sf_data.T3  = sfermions::Isospin[sfermions::electron];
   sf_data.Yl  = sfermions::Hypercharge_left[sfermions::electron];
   sf_data.Yr  = sfermions::Hypercharge_right[sfermions::electron];

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf1, msf2);

   if (msf1 < 0 || msf2 < 0) {
      VERBOSE_MSG("diagonalize_sfermions_2x2: selecton tachyon");
   }

   msf1 = AbsSqrt(msf1);
   msf2 = AbsSqrt(msf2);
}


Eigen::Matrix<double,3,3> CLASSNAME::self_energy_hh_2loop() const
{
   using namespace flexiblesusy::mssm_twoloophiggs;
   using namespace flexiblesusy::nmssm_twoloophiggs;

   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   const double msnusq = Sqr(msnu_2);
   const double sxt = Sin(theta_t), cxt = Cos(theta_t);
   const double sxb = Sin(theta_b), cxb = Cos(theta_b);
   const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
   const double gs = g3;
   const double as = Sqr(gs) / (4.0 * Pi);
   const double rmt = MFu(2);
   const double rmtsq = Sqr(rmt);
   const double scalesq = Sqr(get_scale());
   const double vev2 = Sqr(vd) + Sqr(vu);
   const double tanb = vu/vd;
   const double amu = Re(-0.7071067811865475*vs*Lambdax);
   const double mg = MGlu;
   const double mAsq = (0.7071067811865475*vs*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
   const double cotb = 1.0 / tanb;
   const double rmb = MFd(2);
   const double rmbsq = Sqr(rmb);
   const double rmtausq = Sqr(MFe(2));
   const double lam = Re(Lambdax);
   const double svev = Abs(amu / lam);

   Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,3>::Zero());

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy_2l += self_energy_higgs_2loop_at_as_nmssm(
         rmt, mg, mst1sq, mst2sq, sxt, cxt,
         scalesq, tanb, vev2, lam, svev, as, amu);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy_2l += self_energy_higgs_2loop_ab_as_nmssm(
         rmb, mg, msb1sq, msb2sq, sxb, cxb,
         scalesq, cotb, vev2, lam, svev, as, amu);
   }

   // Corrections as in MSSM, not corrected for NMSSM,
   // should be OK for MSSM states when S state is close to decoupled

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_at_at_mssm(
         rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq,
         msb2sq, sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_atau_atau_mssm(
         rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
         costau, scalesq, amu, tanb, vev2);
   }

   return self_energy_2l;
}

Eigen::Matrix<double,3,3> CLASSNAME::self_energy_Ah_2loop() const
{
   using namespace flexiblesusy::mssm_twoloophiggs;
   using namespace flexiblesusy::nmssm_twoloophiggs;

   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   const double msnusq = Sqr(msnu_2);
   const double sxt = Sin(theta_t), cxt = Cos(theta_t);
   const double sxb = Sin(theta_b), cxb = Cos(theta_b);
   const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
   const double gs = g3;
   const double as = Sqr(gs) / (4.0 * Pi);
   const double rmt = MFu(2);
   const double rmtsq = Sqr(rmt);
   const double scalesq = Sqr(get_scale());
   const double vev2 = Sqr(vd) + Sqr(vu);
   const double tanb = vu/vd;
   const double amu = Re(-0.7071067811865475*vs*Lambdax);
   const double mg = MGlu;
   const double mAsq = (0.7071067811865475*vs*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
   const double cotb = 1.0 / tanb;
   const double rmb = MFd(2);
   const double rmbsq = Sqr(rmb);
   const double rmtausq = Sqr(MFe(2));
   const double lam = Re(Lambdax);
   const double svev = Abs(amu / lam);

   Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,3>::Zero());

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy_2l += self_energy_pseudoscalar_2loop_at_as_nmssm(
         rmt, mg, mst1sq, mst2sq, sxt, cxt,
         scalesq, tanb, vev2, lam, svev, as, amu);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy_2l += self_energy_pseudoscalar_2loop_ab_as_nmssm(
         rmb, mg, msb1sq, msb2sq, sxb, cxb,
         scalesq, cotb, vev2, lam, svev, as, amu);
   }

   // Corrections as in MSSM, not corrected for NMSSM,
   // should be OK for MSSM states when S state is close to decoupled

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_at_at_mssm(
         rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
         sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_atau_atau_mssm(
         rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
         costau, scalesq, amu, tanb, vev2);
   }

   return self_energy_2l;
}



Eigen::Matrix<double,3,1> CLASSNAME::tadpole_hh_2loop() const
{
   using namespace flexiblesusy::mssm_twoloophiggs;
   using namespace flexiblesusy::nmssm_twoloophiggs;

   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   const double msnusq = Sqr(msnu_2);
   const double sxt = Sin(theta_t), cxt = Cos(theta_t);
   const double sxb = Sin(theta_b), cxb = Cos(theta_b);
   const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
   const double gs = g3;
   const double rmtsq = Sqr(MFu(2));
   const double scalesq = Sqr(get_scale());
   const double vev2 = Sqr(vd) + Sqr(vu);
   const double tanb = vu/vd;
   const double amu = Re(-0.7071067811865475*vs*Lambdax);
   const double mg = MGlu;
   const double mAsq = (0.7071067811865475*vs*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
   const double cotbeta = 1.0 / tanb;
   const double rmbsq = Sqr(MFd(2));
   const double rmtausq = Sqr(MFe(2));
   const double lam = Re(Lambdax);
   const double svev = Abs(amu / lam);

   Eigen::Matrix<double,3,1> tadpole_2l(Eigen::Matrix<double,3,1>::Zero());

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      tadpole_2l += tadpole_higgs_2loop_at_as_nmssm(
         rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq,
         amu, tanb, vev2, gs, svev);
   }

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      tadpole_2l.head<2>() += tadpole_higgs_2loop_at_at_mssm(
         rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
         sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      tadpole_2l += tadpole_higgs_2loop_ab_as_nmssm(
         rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
         amu, cotbeta, vev2, gs, svev);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      tadpole_2l.head<2>() += tadpole_higgs_2loop_atau_atau_mssm(
         rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
         costau, scalesq, amu, tanb, vev2);
   }

   tadpole_2l(0) *= vd;
   tadpole_2l(1) *= vu;
   tadpole_2l(2) *= vs;

   return tadpole_2l;
}






void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MGlu_pole()
{
   // diagonalization with medium precision
   const double M_tree(MGlu);
   const double p = MGlu;
   const double self_energy_1  = Re(self_energy_Glu_1loop_1(p));
   const double self_energy_PL = Re(self_energy_Glu_1loop_PL(p));
   const double self_energy_PR = Re(self_energy_Glu_1loop_PR(p));
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL +
      self_energy_PR);
   PHYSICAL(MGlu) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void CLASSNAME::calculate_MChaP_pole()
{
   // diagonalization with medium precision
   const double M_tree(MChaP);
   const double p = MChaP;
   const double self_energy_1  = Re(self_energy_ChaP_1loop_1(p));
   const double self_energy_PL = Re(self_energy_ChaP_1loop_PL(p));
   const double self_energy_PR = Re(self_energy_ChaP_1loop_PR(p));
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL +
      self_energy_PR);
   PHYSICAL(MChaP) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VZ);
   }

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MVZp_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::VZp))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZp));
   const double p = MVZp;
   const double self_energy = Re(self_energy_VZp_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VZp);
   }

   PHYSICAL(MVZp) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::Sd))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sd());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSd(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Sd_1loop(p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZD;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZD,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::Sd, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZD);
      #endif
         normalize_to_interval(mix_ZD);

      PHYSICAL(MSd(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZD) = mix_ZD;
   }
}

void CLASSNAME::calculate_MSv_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::Sv))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Sv());

   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MSv(es));
      Eigen::Matrix<double,3,3> self_energy = Re(self_energy_Sv_1loop(p));
      const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
      Eigen::Array<double,3,1> eigen_values;
      Eigen::Matrix<double,3,3> mix_ZV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZV,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::Sv, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZV);
      #endif
         normalize_to_interval(mix_ZV);

      PHYSICAL(MSv(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZV) = mix_ZV;
   }
}

void CLASSNAME::calculate_MSu_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::Su))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Su());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSu(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Su_1loop(p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZU;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZU,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::Su, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZU);
      #endif
         normalize_to_interval(mix_ZU);

      PHYSICAL(MSu(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZU) = mix_ZU;
   }
}

void CLASSNAME::calculate_MSe_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::Se))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Se());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSe(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Se_1loop(p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZE;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZE,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::Se, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZE);
      #endif
         normalize_to_interval(mix_ZE);

      PHYSICAL(MSe(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZE) = mix_ZE;
   }
}

void CLASSNAME::calculate_MSDX_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::SDX))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_SDX());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSDX(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_SDX_1loop(p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZDX;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZDX,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::SDX, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZDX);
      #endif
         normalize_to_interval(mix_ZDX);

      PHYSICAL(MSDX(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZDX) = mix_ZDX;
   }
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   // two-loop Higgs self-energy contributions
   Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,3>::Zero());

   if (pole_mass_loop_order > 1) {
      self_energy_2l = self_energy_hh_2loop();
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            if (!std::isfinite(self_energy_2l(i,k))) {
               self_energy_2l(i,k) = 0.;
               problems.flag_bad_mass(E6SSM_info::hh);
            }
         }
      }
   }

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_hh());

      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_Mhh(es));
         Eigen::Matrix<double,3,3> self_energy = Re(self_energy_hh_1loop(p));
         self_energy += self_energy_2l;
         const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZH,
               eigenvalue_error);
            problems.flag_bad_mass(E6SSM_info::hh, eigenvalue_error > precision
                * Abs(eigen_values(0)));
         #else

            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZH);
         #endif
            normalize_to_interval(mix_ZH);

         PHYSICAL(Mhh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZH) = mix_ZH;
      }

      new_Mhh = PHYSICAL(Mhh);
      diff = MaxRelDiff(new_Mhh, old_Mhh);
      old_Mhh = new_Mhh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(E6SSM_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(E6SSM_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::Ah))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   // two-loop Higgs self-energy contributions
   Eigen::Matrix<double,3,3> self_energy_2l(Eigen::Matrix<double,3,3>::Zero());

   if (pole_mass_loop_order > 1) {
      self_energy_2l = self_energy_Ah_2loop();
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            if (!std::isfinite(self_energy_2l(i,k))) {
               self_energy_2l(i,k) = 0.;
               problems.flag_bad_mass(E6SSM_info::Ah);
            }
         }
      }
   }

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Ah());

      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_MAh(es));
         Eigen::Matrix<double,3,3> self_energy = Re(self_energy_Ah_1loop(p));
         self_energy += self_energy_2l;
         const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA,
               eigenvalue_error);
            problems.flag_bad_mass(E6SSM_info::Ah, eigenvalue_error > precision
                * Abs(eigen_values(0)));
         #else

            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA);
         #endif
            normalize_to_interval(mix_ZA);

         PHYSICAL(MAh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 2)
            PHYSICAL(ZA) = mix_ZA;
      }

      new_MAh = PHYSICAL(MAh);
      diff = MaxRelDiff(new_MAh, old_MAh);
      old_MAh = new_MAh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(E6SSM_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(E6SSM_info::Ah);
}

void CLASSNAME::calculate_MHpm_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::Hpm))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MHpm) old_MHpm(MHpm), new_MHpm(MHpm);

   do {
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Hpm());

      for (int es = 0; es < 2; ++es) {
         const double p = Abs(old_MHpm(es));
         Eigen::Matrix<double,2,2> self_energy = Re(self_energy_Hpm_1loop(p));
         const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP,
               eigenvalue_error);
            problems.flag_bad_mass(E6SSM_info::Hpm, eigenvalue_error >
               precision * Abs(eigen_values(0)));
         #else

            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP);
         #endif
            normalize_to_interval(mix_ZP);

         PHYSICAL(MHpm(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 1)
            PHYSICAL(ZP) = mix_ZP;
      }

      new_MHpm = PHYSICAL(MHpm);
      diff = MaxRelDiff(new_MHpm, old_MHpm);
      old_MHpm = new_MHpm;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(E6SSM_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(E6SSM_info::Hpm);
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Chi());
   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MChi(es));
      const Eigen::Matrix<double,6,6> self_energy_1  = Re(
         self_energy_Chi_1loop_1(p));
      const Eigen::Matrix<double,6,6> self_energy_PL = Re(
         self_energy_Chi_1loop_PL(p));
      const Eigen::Matrix<double,6,6> self_energy_PR = Re(
         self_energy_Chi_1loop_PR(p));
      const Eigen::Matrix<double,6,6> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,6,6> M_loop(M_tree + 0.5 * (delta_M + delta_M.
         transpose()));
      Eigen::Array<double,6,1> eigen_values;
      decltype(ZN) mix_ZN;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZN,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::Chi, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZN);
      #endif
         normalize_to_interval(mix_ZN);
      if (es == 0)
         PHYSICAL(ZN) = mix_ZN;
      PHYSICAL(MChi(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MCha_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha());
   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MCha(es));
      const Eigen::Matrix<double,2,2> self_energy_1  = Re(
         self_energy_Cha_1loop_1(p));
      const Eigen::Matrix<double,2,2> self_energy_PL = Re(
         self_energy_Cha_1loop_PL(p));
      const Eigen::Matrix<double,2,2> self_energy_PR = Re(
         self_energy_Cha_1loop_PR(p));
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM) mix_UM;
      decltype(UP) mix_UP;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_UM, mix_UP, eigenvalue_error);
      problems.flag_bad_mass(E6SSM_info::Cha, eigenvalue_error > precision *
         Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_UM, mix_UP);
   #endif
      if (es == 0) {
         PHYSICAL(UM) = mix_UM;
         PHYSICAL(UP) = mix_UP;
      }
      PHYSICAL(MCha(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFe(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fe_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fe_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fe_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZEL) mix_ZEL;
      decltype(ZER) mix_ZER;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER, eigenvalue_error);
      problems.flag_bad_mass(E6SSM_info::Fe, eigenvalue_error > precision * Abs
         (eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER);
   #endif
      if (es == 0) {
         PHYSICAL(ZEL) = mix_ZEL;
         PHYSICAL(ZER) = mix_ZER;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFd(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fd_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fd_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fd_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZDL) mix_ZDL;
      decltype(ZDR) mix_ZDR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR, eigenvalue_error);
      problems.flag_bad_mass(E6SSM_info::Fd, eigenvalue_error > precision * Abs
         (eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR);
   #endif
      if (es == 0) {
         PHYSICAL(ZDL) = mix_ZDL;
         PHYSICAL(ZDR) = mix_ZDR;
      }
      PHYSICAL(MFd(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFu_pole()
{
   // diagonalization with medium precision
   double qcd_1l = 0.;

   {
      const double currentScale = get_scale();
      qcd_1l = 0.008443431970194815*(-5. + 3.*Log(Sqr(MFu(2))/Sqr(currentScale)
         ))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = 2.2278607323533713e-6*Quad(g3)*(-2330.129769909197 + 1476.*Log(
         Sqr(MFu(2))/Sqr(currentScale)) - 396.*Sqr(Log(Sqr(MFu(2))/Sqr(
         currentScale))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = 0;
   }

   double qcd_4l = 0.;

   if (pole_mass_loop_order > 3 && TOP_POLE_QCD_CORRECTION > 2) {
      const double currentScale = get_scale();
      qcd_4l = 0;
   }

   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      for (int i1 = 0; i1 < 3; ++i1) {
         for (int i2 = 0; i2 < 3; ++i2) {
            if (i1 == 2 && i2 == 2) {
               self_energy_1(i1,i2)  = Re(self_energy_Fu_1loop_1_heavy(p,i1,i2)
                  );
               self_energy_PL(i1,i2) = Re(self_energy_Fu_1loop_PL_heavy(p,i1,i2
                  ));
               self_energy_PR(i1,i2) = Re(self_energy_Fu_1loop_PR_heavy(p,i1,i2
                  ));
            } else {
               self_energy_1(i1,i2)  = Re(self_energy_Fu_1loop_1(p,i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fu_1loop_PL(p,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fu_1loop_PR(p,i1,i2));
            }
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree - M_tree *
         self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l + qcd_4l);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZUL) mix_ZUL;
      decltype(ZUR) mix_ZUR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR, eigenvalue_error);
      problems.flag_bad_mass(E6SSM_info::Fu, eigenvalue_error > precision * Abs
         (eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR);
   #endif
      if (es == 0) {
         PHYSICAL(ZUL) = mix_ZUL;
         PHYSICAL(ZUR) = mix_ZUR;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFDX_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_FDX());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFDX(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_FDX_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_FDX_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_FDX_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZDXL) mix_ZDXL;
      decltype(ZDXR) mix_ZDXR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZDXL, mix_ZDXR, eigenvalue_error);
      problems.flag_bad_mass(E6SSM_info::FDX, eigenvalue_error > precision *
         Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZDXL, mix_ZDXR);
   #endif
      if (es == 0) {
         PHYSICAL(ZDXL) = mix_ZDXL;
         PHYSICAL(ZDXR) = mix_ZDXR;
      }
      PHYSICAL(MFDX(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MSHI0_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::SHI0))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_SHI0());

   for (int es = 0; es < 4; ++es) {
      const double p = Abs(MSHI0(es));
      Eigen::Matrix<double,4,4> self_energy = Re(self_energy_SHI0_1loop(p));
      const Eigen::Matrix<double,4,4> M_loop(M_tree - self_energy);
      Eigen::Array<double,4,1> eigen_values;
      Eigen::Matrix<double,4,4> mix_UHI0;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHI0,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::SHI0, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHI0);
      #endif
         normalize_to_interval(mix_UHI0);

      PHYSICAL(MSHI0(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(UHI0) = mix_UHI0;
   }
}

void CLASSNAME::calculate_MSHIp_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::SHIp))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_SHIp());

   for (int es = 0; es < 4; ++es) {
      const double p = Abs(MSHIp(es));
      Eigen::Matrix<double,4,4> self_energy = Re(self_energy_SHIp_1loop(p));
      const Eigen::Matrix<double,4,4> M_loop(M_tree - self_energy);
      Eigen::Array<double,4,1> eigen_values;
      Eigen::Matrix<double,4,4> mix_UHIp;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHIp,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::SHIp, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHIp);
      #endif
         normalize_to_interval(mix_UHIp);

      PHYSICAL(MSHIp(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(UHIp) = mix_UHIp;
   }
}

void CLASSNAME::calculate_MChaI_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_ChaI());
   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MChaI(es));
      const Eigen::Matrix<double,2,2> self_energy_1  = Re(
         self_energy_ChaI_1loop_1(p));
      const Eigen::Matrix<double,2,2> self_energy_PL = Re(
         self_energy_ChaI_1loop_PL(p));
      const Eigen::Matrix<double,2,2> self_energy_PR = Re(
         self_energy_ChaI_1loop_PR(p));
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(ZMI) mix_ZMI;
      decltype(ZPI) mix_ZPI;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZMI, mix_ZPI, eigenvalue_error);
      problems.flag_bad_mass(E6SSM_info::ChaI, eigenvalue_error > precision *
         Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZMI, mix_ZPI);
   #endif
      if (es == 0) {
         PHYSICAL(ZMI) = mix_ZMI;
         PHYSICAL(ZPI) = mix_ZPI;
      }
      PHYSICAL(MChaI(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MChiI_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_ChiI());
   for (int es = 0; es < 4; ++es) {
      const double p = Abs(MChiI(es));
      const Eigen::Matrix<double,4,4> self_energy_1  = Re(
         self_energy_ChiI_1loop_1(p));
      const Eigen::Matrix<double,4,4> self_energy_PL = Re(
         self_energy_ChiI_1loop_PL(p));
      const Eigen::Matrix<double,4,4> self_energy_PR = Re(
         self_energy_ChiI_1loop_PR(p));
      const Eigen::Matrix<double,4,4> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,4,4> M_loop(M_tree + 0.5 * (delta_M + delta_M.
         transpose()));
      Eigen::Array<double,4,1> eigen_values;
      decltype(ZNI) mix_ZNI;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZNI,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::ChiI, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZNI);
      #endif
         normalize_to_interval(mix_ZNI);
      if (es == 0)
         PHYSICAL(ZNI) = mix_ZNI;
      PHYSICAL(MChiI(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MSSI0_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::SSI0))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_SSI0());

   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MSSI0(es));
      Eigen::Matrix<double,2,2> self_energy = Re(self_energy_SSI0_1loop(p));
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZSSI;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZSSI,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::SSI0, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZSSI);
      #endif
         normalize_to_interval(mix_ZSSI);

      PHYSICAL(MSSI0(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZSSI) = mix_ZSSI;
   }
}

void CLASSNAME::calculate_MFSI_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_FSI());
   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MFSI(es));
      const Eigen::Matrix<double,2,2> self_energy_1  = Re(
         self_energy_FSI_1loop_1(p));
      const Eigen::Matrix<double,2,2> self_energy_PL = Re(
         self_energy_FSI_1loop_PL(p));
      const Eigen::Matrix<double,2,2> self_energy_PR = Re(
         self_energy_FSI_1loop_PR(p));
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + 0.5 * (delta_M + delta_M.
         transpose()));
      Eigen::Array<double,2,1> eigen_values;
      decltype(ZFSI) mix_ZFSI;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZFSI,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::FSI, eigenvalue_error > precision *
            Abs(eigen_values(0)));
      #else

         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZFSI);
      #endif
         normalize_to_interval(mix_ZFSI);
      if (es == 0)
         PHYSICAL(ZFSI) = mix_ZFSI;
      PHYSICAL(MFSI(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MSHp0_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::SHp0))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_SHp0());

   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MSHp0(es));
      Eigen::Matrix<double,2,2> self_energy = Re(self_energy_SHp0_1loop(p));
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_UHp0;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHp0,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::SHp0, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHp0);
      #endif
         normalize_to_interval(mix_UHp0);

      PHYSICAL(MSHp0(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(UHp0) = mix_UHp0;
   }
}

void CLASSNAME::calculate_MSHpp_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::SHpp))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_SHpp());

   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MSHpp(es));
      Eigen::Matrix<double,2,2> self_energy = Re(self_energy_SHpp_1loop(p));
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_UHpp;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHpp,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::SHpp, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_hermitian(M_loop, eigen_values, mix_UHpp);
      #endif
         normalize_to_interval(mix_UHpp);

      PHYSICAL(MSHpp(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(UHpp) = mix_UHpp;
   }
}

void CLASSNAME::calculate_MChiP_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_ChiP());
   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MChiP(es));
      const Eigen::Matrix<double,2,2> self_energy_1  = Re(
         self_energy_ChiP_1loop_1(p));
      const Eigen::Matrix<double,2,2> self_energy_PL = Re(
         self_energy_ChiP_1loop_PL(p));
      const Eigen::Matrix<double,2,2> self_energy_PR = Re(
         self_energy_ChiP_1loop_PR(p));
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + 0.5 * (delta_M + delta_M.
         transpose()));
      Eigen::Array<double,2,1> eigen_values;
      decltype(ZNp) mix_ZNp;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZNp,
            eigenvalue_error);
         problems.flag_bad_mass(E6SSM_info::ChiP, eigenvalue_error > precision
            * Abs(eigen_values(0)));
      #else

         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZNp);
      #endif
         normalize_to_interval(mix_ZNp);
      if (es == 0)
         PHYSICAL(ZNp) = mix_ZNp;
      PHYSICAL(MChiP(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MVWm_pole()
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VWm);
   }

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VWm);
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(E6SSM_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VZ);
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1 - 0.0023747152416172916*(0.6*Sqr(g1) - Sqr
      (g2));
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mf_1loop = - self_energy_1/m_sm_drbar - self_energy_PL -
      self_energy_PR;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mf_1loop);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_1loop_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   double qcd_1l = 0., qcd_2l = 0., qcd_3l = 0., qcd_4l = 0.;
   double atas_S_2l = 0., atas_LR_2l = 0., atat_S_2l = 0., atat_LR_2l = 0.;

   qcd_1l = 0.008443431970194815*(-5. + 3.*Log(Sqr(MFu(idx))/Sqr(currentScale))
      )*Sqr(g3);

   if (get_thresholds() > 1 && threshold_corrections.mt > 1) {
      const double q_2l = 2.2278607323533713e-6*Quad(g3)*(2330.129769909197 -
         1476.*Log(Sqr(MFu(idx))/Sqr(currentScale)) + 396.*Sqr(Log(Sqr(MFu(idx)
         )/Sqr(currentScale))));

      qcd_2l = -q_2l + qcd_1l * qcd_1l;
   }

   const double m_susy_drbar = m_pole + self_energy_1 + atas_S_2l + atat_S_2l +
      m_pole * (self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l +
      qcd_4l + atas_LR_2l + atat_LR_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1loop_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1 + 0.0006860288475783287*Sqr(g1) +
      0.0023747152416172916*Sqr(g2) - 0.008443431970194815*Sqr(g3);
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mb_1loop = - self_energy_1/m_tree - self_energy_PL -
      self_energy_PR;
   double qcd_2l = 0.;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mb_1loop + qcd_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double) const
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole) const
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VZ);return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole) const
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(E6SSM_info::VWm);return m_pole;
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::v() const
{

   return Sqrt(Sqr(vd) + Sqr(vu));
}

double CLASSNAME::Betax() const
{

   return -ArcSin(ZA(0,1));
}

double CLASSNAME::ThetaW() const
{

   return ArcCos(Abs(ZZ(0,0)));
}

double CLASSNAME::ThetaWp() const
{

   return ArcCos(Abs(ZZ(2,2)));
}

double CLASSNAME::VEV() const
{

   return Sqrt(Sqr(vd) + Sqr(vu));
}



std::ostream& operator<<(std::ostream& ostr, const CLASSNAME& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
