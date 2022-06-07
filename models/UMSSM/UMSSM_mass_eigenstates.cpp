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
 * @file UMSSM_mass_eigenstates.cpp
 * @brief implementation of the UMSSM model class
 *
 * Contains the definition of the UMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "UMSSM_mass_eigenstates.hpp"
#include "UMSSM_ewsb_solver_interface.hpp"
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
#include "UMSSM_two_scale_ewsb_solver.hpp"
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
#define CLASSNAME UMSSM_mass_eigenstates

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

CLASSNAME::CLASSNAME(const UMSSM_input_parameters& input_)
   : UMSSM_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new UMSSM_ewsb_solver<Two_scale>())
#endif
{
}

std::unique_ptr<UMSSM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<UMSSM_mass_eigenstates>(*this);
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

const UMSSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

UMSSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const UMSSM_physical& physical_)
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

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<UMSSM_ewsb_solver_interface>& solver)
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
   tadpole[2] /= vS;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;

   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_ms2 = ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) - 3*Cube(vd)*Sqr(g1)
      - 5*Cube(vd)*Sqr(g2) - 20*Cube(vd)*Sqr(gp)*Sqr(QHd) - 20*vd*AbsSqr(Lambdax)*
      Sqr(vS) - 20*QHd*Qs*vd*Sqr(gp)*Sqr(vS) - 20*vd*AbsSqr(Lambdax)*Sqr(vu) + 3*
      vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) - 20*QHd*QHu*vd*Sqr(gp)*Sqr(vu) +
      14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) - 3*Cube(vu)*Sqr(g1)
      - 5*Cube(vu)*Sqr(g2) - 20*Cube(vu)*Sqr(gp)*Sqr(QHu) - 20*vu*AbsSqr(Lambdax)*
      Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*QHd*QHu*vu*Sqr(gp
      )*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) - 20*QHu*Qs*vu*Sqr(gp)*Sqr(vS) +
      14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(1.4142135623730951*vd*vu*Conj(TLambdax) - 2*Cube(vS)*Sqr(gp)*
      Sqr(Qs) - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 2*QHd*Qs*vS*Sqr(gp)*Sqr(vd) - 2*vS*
      AbsSqr(Lambdax)*Sqr(vu) - 2*QHu*Qs*vS*Sqr(gp)*Sqr(vu) + 1.4142135623730951*
      vd*vu*TLambdax))/vS);

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

   VERBOSE_MSG("\t\tSolving UMSSM EWSB at " << ewsb_loop_order << "-loop order");

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
           "UMSSM\n"
           "========================================\n";
   UMSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
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
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "ZVL = " << ZVL << '\n';
   ostr << "ZVR = " << ZVR << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
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
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MFv();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
   calculate_MGlu();
   calculate_MVG();

}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 19u));

   if (calculate_bsm_pole_masses) {
      tp.run_task([this] () { calculate_MAh_pole(); });
      tp.run_task([this] () { calculate_MCha_pole(); });
      tp.run_task([this] () { calculate_MChi_pole(); });
      tp.run_task([this] () { calculate_MGlu_pole(); });
      tp.run_task([this] () { calculate_Mhh_pole(); });
      tp.run_task([this] () { calculate_MHpm_pole(); });
      tp.run_task([this] () { calculate_MSd_pole(); });
      tp.run_task([this] () { calculate_MSe_pole(); });
      tp.run_task([this] () { calculate_MSu_pole(); });
      tp.run_task([this] () { calculate_MSv_pole(); });
      tp.run_task([this] () { calculate_MVZp_pole(); });
   }

   if (calculate_sm_pole_masses) {
      tp.run_task([this] () { calculate_MVG_pole(); });
      tp.run_task([this] () { calculate_MVP_pole(); });
      tp.run_task([this] () { calculate_MVZ_pole(); });
      tp.run_task([this] () { calculate_MFv_pole(); });
      tp.run_task([this] () { calculate_MFe_pole(); });
      tp.run_task([this] () { calculate_MFd_pole(); });
      tp.run_task([this] () { calculate_MFu_pole(); });
      tp.run_task([this] () { calculate_MVWm_pole(); });
   }

#else
   if (calculate_bsm_pole_masses) {
      calculate_MAh_pole();
      calculate_MCha_pole();
      calculate_MChi_pole();
      calculate_MGlu_pole();
      calculate_Mhh_pole();
      calculate_MHpm_pole();
      calculate_MSd_pole();
      calculate_MSe_pole();
      calculate_MSu_pole();
      calculate_MSv_pole();
      calculate_MVZp_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFv_pole();
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
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSv) = MSv;
   PHYSICAL(ZV) = ZV;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(ZVL) = ZVL;
   PHYSICAL(ZVR) = ZVR;
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
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Sd); }
   if (PHYSICAL(MSv).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Se); }
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::hh); }
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Ah); }
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) { problems.flag_pole_tachyon(UMSSM_info::Hpm); }
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
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,6,1>::Zero();
   ZV = Eigen::Matrix<double,6,6>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,6,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,6,6>::Zero();
   MFv = Eigen::Matrix<double,3,1>::Zero();
   ZVL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZVR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
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
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;

   PhaseGlu = std::complex<double>(1.,0.);


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   UMSSM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MSd(0) = pars(2);
   MSd(1) = pars(3);
   MSd(2) = pars(4);
   MSd(3) = pars(5);
   MSd(4) = pars(6);
   MSd(5) = pars(7);
   MSv(0) = pars(8);
   MSv(1) = pars(9);
   MSv(2) = pars(10);
   MSv(3) = pars(11);
   MSv(4) = pars(12);
   MSv(5) = pars(13);
   MSu(0) = pars(14);
   MSu(1) = pars(15);
   MSu(2) = pars(16);
   MSu(3) = pars(17);
   MSu(4) = pars(18);
   MSu(5) = pars(19);
   MSe(0) = pars(20);
   MSe(1) = pars(21);
   MSe(2) = pars(22);
   MSe(3) = pars(23);
   MSe(4) = pars(24);
   MSe(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   Mhh(2) = pars(28);
   MAh(0) = pars(29);
   MAh(1) = pars(30);
   MAh(2) = pars(31);
   MHpm(0) = pars(32);
   MHpm(1) = pars(33);
   MChi(0) = pars(34);
   MChi(1) = pars(35);
   MChi(2) = pars(36);
   MChi(3) = pars(37);
   MChi(4) = pars(38);
   MChi(5) = pars(39);
   MFv(0) = pars(40);
   MFv(1) = pars(41);
   MFv(2) = pars(42);
   MCha(0) = pars(43);
   MCha(1) = pars(44);
   MFe(0) = pars(45);
   MFe(1) = pars(46);
   MFe(2) = pars(47);
   MFd(0) = pars(48);
   MFd(1) = pars(49);
   MFd(2) = pars(50);
   MFu(0) = pars(51);
   MFu(1) = pars(52);
   MFu(2) = pars(53);
   MVWm = pars(54);
   MVP = pars(55);
   MVZ = pars(56);
   MVZp = pars(57);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(58);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MSd(0);
   pars(3) = MSd(1);
   pars(4) = MSd(2);
   pars(5) = MSd(3);
   pars(6) = MSd(4);
   pars(7) = MSd(5);
   pars(8) = MSv(0);
   pars(9) = MSv(1);
   pars(10) = MSv(2);
   pars(11) = MSv(3);
   pars(12) = MSv(4);
   pars(13) = MSv(5);
   pars(14) = MSu(0);
   pars(15) = MSu(1);
   pars(16) = MSu(2);
   pars(17) = MSu(3);
   pars(18) = MSu(4);
   pars(19) = MSu(5);
   pars(20) = MSe(0);
   pars(21) = MSe(1);
   pars(22) = MSe(2);
   pars(23) = MSe(3);
   pars(24) = MSe(4);
   pars(25) = MSe(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = Mhh(2);
   pars(29) = MAh(0);
   pars(30) = MAh(1);
   pars(31) = MAh(2);
   pars(32) = MHpm(0);
   pars(33) = MHpm(1);
   pars(34) = MChi(0);
   pars(35) = MChi(1);
   pars(36) = MChi(2);
   pars(37) = MChi(3);
   pars(38) = MChi(4);
   pars(39) = MChi(5);
   pars(40) = MFv(0);
   pars(41) = MFv(1);
   pars(42) = MFv(2);
   pars(43) = MCha(0);
   pars(44) = MCha(1);
   pars(45) = MFe(0);
   pars(46) = MFe(1);
   pars(47) = MFe(2);
   pars(48) = MFd(0);
   pars(49) = MFd(1);
   pars(50) = MFd(2);
   pars(51) = MFu(0);
   pars(52) = MFu(1);
   pars(53) = MFu(2);
   pars(54) = MVWm;
   pars(55) = MVP;
   pars(56) = MVZ;
   pars(57) = MVZp;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   ZD(0,0) = pars(58);
   ZD(0,1) = pars(59);
   ZD(0,2) = pars(60);
   ZD(0,3) = pars(61);
   ZD(0,4) = pars(62);
   ZD(0,5) = pars(63);
   ZD(1,0) = pars(64);
   ZD(1,1) = pars(65);
   ZD(1,2) = pars(66);
   ZD(1,3) = pars(67);
   ZD(1,4) = pars(68);
   ZD(1,5) = pars(69);
   ZD(2,0) = pars(70);
   ZD(2,1) = pars(71);
   ZD(2,2) = pars(72);
   ZD(2,3) = pars(73);
   ZD(2,4) = pars(74);
   ZD(2,5) = pars(75);
   ZD(3,0) = pars(76);
   ZD(3,1) = pars(77);
   ZD(3,2) = pars(78);
   ZD(3,3) = pars(79);
   ZD(3,4) = pars(80);
   ZD(3,5) = pars(81);
   ZD(4,0) = pars(82);
   ZD(4,1) = pars(83);
   ZD(4,2) = pars(84);
   ZD(4,3) = pars(85);
   ZD(4,4) = pars(86);
   ZD(4,5) = pars(87);
   ZD(5,0) = pars(88);
   ZD(5,1) = pars(89);
   ZD(5,2) = pars(90);
   ZD(5,3) = pars(91);
   ZD(5,4) = pars(92);
   ZD(5,5) = pars(93);
   ZV(0,0) = pars(94);
   ZV(0,1) = pars(95);
   ZV(0,2) = pars(96);
   ZV(0,3) = pars(97);
   ZV(0,4) = pars(98);
   ZV(0,5) = pars(99);
   ZV(1,0) = pars(100);
   ZV(1,1) = pars(101);
   ZV(1,2) = pars(102);
   ZV(1,3) = pars(103);
   ZV(1,4) = pars(104);
   ZV(1,5) = pars(105);
   ZV(2,0) = pars(106);
   ZV(2,1) = pars(107);
   ZV(2,2) = pars(108);
   ZV(2,3) = pars(109);
   ZV(2,4) = pars(110);
   ZV(2,5) = pars(111);
   ZV(3,0) = pars(112);
   ZV(3,1) = pars(113);
   ZV(3,2) = pars(114);
   ZV(3,3) = pars(115);
   ZV(3,4) = pars(116);
   ZV(3,5) = pars(117);
   ZV(4,0) = pars(118);
   ZV(4,1) = pars(119);
   ZV(4,2) = pars(120);
   ZV(4,3) = pars(121);
   ZV(4,4) = pars(122);
   ZV(4,5) = pars(123);
   ZV(5,0) = pars(124);
   ZV(5,1) = pars(125);
   ZV(5,2) = pars(126);
   ZV(5,3) = pars(127);
   ZV(5,4) = pars(128);
   ZV(5,5) = pars(129);
   ZU(0,0) = pars(130);
   ZU(0,1) = pars(131);
   ZU(0,2) = pars(132);
   ZU(0,3) = pars(133);
   ZU(0,4) = pars(134);
   ZU(0,5) = pars(135);
   ZU(1,0) = pars(136);
   ZU(1,1) = pars(137);
   ZU(1,2) = pars(138);
   ZU(1,3) = pars(139);
   ZU(1,4) = pars(140);
   ZU(1,5) = pars(141);
   ZU(2,0) = pars(142);
   ZU(2,1) = pars(143);
   ZU(2,2) = pars(144);
   ZU(2,3) = pars(145);
   ZU(2,4) = pars(146);
   ZU(2,5) = pars(147);
   ZU(3,0) = pars(148);
   ZU(3,1) = pars(149);
   ZU(3,2) = pars(150);
   ZU(3,3) = pars(151);
   ZU(3,4) = pars(152);
   ZU(3,5) = pars(153);
   ZU(4,0) = pars(154);
   ZU(4,1) = pars(155);
   ZU(4,2) = pars(156);
   ZU(4,3) = pars(157);
   ZU(4,4) = pars(158);
   ZU(4,5) = pars(159);
   ZU(5,0) = pars(160);
   ZU(5,1) = pars(161);
   ZU(5,2) = pars(162);
   ZU(5,3) = pars(163);
   ZU(5,4) = pars(164);
   ZU(5,5) = pars(165);
   ZE(0,0) = pars(166);
   ZE(0,1) = pars(167);
   ZE(0,2) = pars(168);
   ZE(0,3) = pars(169);
   ZE(0,4) = pars(170);
   ZE(0,5) = pars(171);
   ZE(1,0) = pars(172);
   ZE(1,1) = pars(173);
   ZE(1,2) = pars(174);
   ZE(1,3) = pars(175);
   ZE(1,4) = pars(176);
   ZE(1,5) = pars(177);
   ZE(2,0) = pars(178);
   ZE(2,1) = pars(179);
   ZE(2,2) = pars(180);
   ZE(2,3) = pars(181);
   ZE(2,4) = pars(182);
   ZE(2,5) = pars(183);
   ZE(3,0) = pars(184);
   ZE(3,1) = pars(185);
   ZE(3,2) = pars(186);
   ZE(3,3) = pars(187);
   ZE(3,4) = pars(188);
   ZE(3,5) = pars(189);
   ZE(4,0) = pars(190);
   ZE(4,1) = pars(191);
   ZE(4,2) = pars(192);
   ZE(4,3) = pars(193);
   ZE(4,4) = pars(194);
   ZE(4,5) = pars(195);
   ZE(5,0) = pars(196);
   ZE(5,1) = pars(197);
   ZE(5,2) = pars(198);
   ZE(5,3) = pars(199);
   ZE(5,4) = pars(200);
   ZE(5,5) = pars(201);
   ZH(0,0) = pars(202);
   ZH(0,1) = pars(203);
   ZH(0,2) = pars(204);
   ZH(1,0) = pars(205);
   ZH(1,1) = pars(206);
   ZH(1,2) = pars(207);
   ZH(2,0) = pars(208);
   ZH(2,1) = pars(209);
   ZH(2,2) = pars(210);
   ZA(0,0) = pars(211);
   ZA(0,1) = pars(212);
   ZA(0,2) = pars(213);
   ZA(1,0) = pars(214);
   ZA(1,1) = pars(215);
   ZA(1,2) = pars(216);
   ZA(2,0) = pars(217);
   ZA(2,1) = pars(218);
   ZA(2,2) = pars(219);
   ZP(0,0) = pars(220);
   ZP(0,1) = pars(221);
   ZP(1,0) = pars(222);
   ZP(1,1) = pars(223);
   ZN(0,0) = std::complex<double>(pars(224), pars(225));
   ZN(0,1) = std::complex<double>(pars(226), pars(227));
   ZN(0,2) = std::complex<double>(pars(228), pars(229));
   ZN(0,3) = std::complex<double>(pars(230), pars(231));
   ZN(0,4) = std::complex<double>(pars(232), pars(233));
   ZN(0,5) = std::complex<double>(pars(234), pars(235));
   ZN(1,0) = std::complex<double>(pars(236), pars(237));
   ZN(1,1) = std::complex<double>(pars(238), pars(239));
   ZN(1,2) = std::complex<double>(pars(240), pars(241));
   ZN(1,3) = std::complex<double>(pars(242), pars(243));
   ZN(1,4) = std::complex<double>(pars(244), pars(245));
   ZN(1,5) = std::complex<double>(pars(246), pars(247));
   ZN(2,0) = std::complex<double>(pars(248), pars(249));
   ZN(2,1) = std::complex<double>(pars(250), pars(251));
   ZN(2,2) = std::complex<double>(pars(252), pars(253));
   ZN(2,3) = std::complex<double>(pars(254), pars(255));
   ZN(2,4) = std::complex<double>(pars(256), pars(257));
   ZN(2,5) = std::complex<double>(pars(258), pars(259));
   ZN(3,0) = std::complex<double>(pars(260), pars(261));
   ZN(3,1) = std::complex<double>(pars(262), pars(263));
   ZN(3,2) = std::complex<double>(pars(264), pars(265));
   ZN(3,3) = std::complex<double>(pars(266), pars(267));
   ZN(3,4) = std::complex<double>(pars(268), pars(269));
   ZN(3,5) = std::complex<double>(pars(270), pars(271));
   ZN(4,0) = std::complex<double>(pars(272), pars(273));
   ZN(4,1) = std::complex<double>(pars(274), pars(275));
   ZN(4,2) = std::complex<double>(pars(276), pars(277));
   ZN(4,3) = std::complex<double>(pars(278), pars(279));
   ZN(4,4) = std::complex<double>(pars(280), pars(281));
   ZN(4,5) = std::complex<double>(pars(282), pars(283));
   ZN(5,0) = std::complex<double>(pars(284), pars(285));
   ZN(5,1) = std::complex<double>(pars(286), pars(287));
   ZN(5,2) = std::complex<double>(pars(288), pars(289));
   ZN(5,3) = std::complex<double>(pars(290), pars(291));
   ZN(5,4) = std::complex<double>(pars(292), pars(293));
   ZN(5,5) = std::complex<double>(pars(294), pars(295));
   ZVL(0,0) = std::complex<double>(pars(296), pars(297));
   ZVL(0,1) = std::complex<double>(pars(298), pars(299));
   ZVL(0,2) = std::complex<double>(pars(300), pars(301));
   ZVL(1,0) = std::complex<double>(pars(302), pars(303));
   ZVL(1,1) = std::complex<double>(pars(304), pars(305));
   ZVL(1,2) = std::complex<double>(pars(306), pars(307));
   ZVL(2,0) = std::complex<double>(pars(308), pars(309));
   ZVL(2,1) = std::complex<double>(pars(310), pars(311));
   ZVL(2,2) = std::complex<double>(pars(312), pars(313));
   ZVR(0,0) = std::complex<double>(pars(314), pars(315));
   ZVR(0,1) = std::complex<double>(pars(316), pars(317));
   ZVR(0,2) = std::complex<double>(pars(318), pars(319));
   ZVR(1,0) = std::complex<double>(pars(320), pars(321));
   ZVR(1,1) = std::complex<double>(pars(322), pars(323));
   ZVR(1,2) = std::complex<double>(pars(324), pars(325));
   ZVR(2,0) = std::complex<double>(pars(326), pars(327));
   ZVR(2,1) = std::complex<double>(pars(328), pars(329));
   ZVR(2,2) = std::complex<double>(pars(330), pars(331));
   UM(0,0) = std::complex<double>(pars(332), pars(333));
   UM(0,1) = std::complex<double>(pars(334), pars(335));
   UM(1,0) = std::complex<double>(pars(336), pars(337));
   UM(1,1) = std::complex<double>(pars(338), pars(339));
   UP(0,0) = std::complex<double>(pars(340), pars(341));
   UP(0,1) = std::complex<double>(pars(342), pars(343));
   UP(1,0) = std::complex<double>(pars(344), pars(345));
   UP(1,1) = std::complex<double>(pars(346), pars(347));
   ZEL(0,0) = std::complex<double>(pars(348), pars(349));
   ZEL(0,1) = std::complex<double>(pars(350), pars(351));
   ZEL(0,2) = std::complex<double>(pars(352), pars(353));
   ZEL(1,0) = std::complex<double>(pars(354), pars(355));
   ZEL(1,1) = std::complex<double>(pars(356), pars(357));
   ZEL(1,2) = std::complex<double>(pars(358), pars(359));
   ZEL(2,0) = std::complex<double>(pars(360), pars(361));
   ZEL(2,1) = std::complex<double>(pars(362), pars(363));
   ZEL(2,2) = std::complex<double>(pars(364), pars(365));
   ZER(0,0) = std::complex<double>(pars(366), pars(367));
   ZER(0,1) = std::complex<double>(pars(368), pars(369));
   ZER(0,2) = std::complex<double>(pars(370), pars(371));
   ZER(1,0) = std::complex<double>(pars(372), pars(373));
   ZER(1,1) = std::complex<double>(pars(374), pars(375));
   ZER(1,2) = std::complex<double>(pars(376), pars(377));
   ZER(2,0) = std::complex<double>(pars(378), pars(379));
   ZER(2,1) = std::complex<double>(pars(380), pars(381));
   ZER(2,2) = std::complex<double>(pars(382), pars(383));
   ZDL(0,0) = std::complex<double>(pars(384), pars(385));
   ZDL(0,1) = std::complex<double>(pars(386), pars(387));
   ZDL(0,2) = std::complex<double>(pars(388), pars(389));
   ZDL(1,0) = std::complex<double>(pars(390), pars(391));
   ZDL(1,1) = std::complex<double>(pars(392), pars(393));
   ZDL(1,2) = std::complex<double>(pars(394), pars(395));
   ZDL(2,0) = std::complex<double>(pars(396), pars(397));
   ZDL(2,1) = std::complex<double>(pars(398), pars(399));
   ZDL(2,2) = std::complex<double>(pars(400), pars(401));
   ZDR(0,0) = std::complex<double>(pars(402), pars(403));
   ZDR(0,1) = std::complex<double>(pars(404), pars(405));
   ZDR(0,2) = std::complex<double>(pars(406), pars(407));
   ZDR(1,0) = std::complex<double>(pars(408), pars(409));
   ZDR(1,1) = std::complex<double>(pars(410), pars(411));
   ZDR(1,2) = std::complex<double>(pars(412), pars(413));
   ZDR(2,0) = std::complex<double>(pars(414), pars(415));
   ZDR(2,1) = std::complex<double>(pars(416), pars(417));
   ZDR(2,2) = std::complex<double>(pars(418), pars(419));
   ZUL(0,0) = std::complex<double>(pars(420), pars(421));
   ZUL(0,1) = std::complex<double>(pars(422), pars(423));
   ZUL(0,2) = std::complex<double>(pars(424), pars(425));
   ZUL(1,0) = std::complex<double>(pars(426), pars(427));
   ZUL(1,1) = std::complex<double>(pars(428), pars(429));
   ZUL(1,2) = std::complex<double>(pars(430), pars(431));
   ZUL(2,0) = std::complex<double>(pars(432), pars(433));
   ZUL(2,1) = std::complex<double>(pars(434), pars(435));
   ZUL(2,2) = std::complex<double>(pars(436), pars(437));
   ZUR(0,0) = std::complex<double>(pars(438), pars(439));
   ZUR(0,1) = std::complex<double>(pars(440), pars(441));
   ZUR(0,2) = std::complex<double>(pars(442), pars(443));
   ZUR(1,0) = std::complex<double>(pars(444), pars(445));
   ZUR(1,1) = std::complex<double>(pars(446), pars(447));
   ZUR(1,2) = std::complex<double>(pars(448), pars(449));
   ZUR(2,0) = std::complex<double>(pars(450), pars(451));
   ZUR(2,1) = std::complex<double>(pars(452), pars(453));
   ZUR(2,2) = std::complex<double>(pars(454), pars(455));
   ZZ(0,0) = pars(456);
   ZZ(0,1) = pars(457);
   ZZ(0,2) = pars(458);
   ZZ(1,0) = pars(459);
   ZZ(1,1) = pars(460);
   ZZ(1,2) = pars(461);
   ZZ(2,0) = pars(462);
   ZZ(2,1) = pars(463);
   ZZ(2,2) = pars(464);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(465);

   pars(58) = ZD(0,0);
   pars(59) = ZD(0,1);
   pars(60) = ZD(0,2);
   pars(61) = ZD(0,3);
   pars(62) = ZD(0,4);
   pars(63) = ZD(0,5);
   pars(64) = ZD(1,0);
   pars(65) = ZD(1,1);
   pars(66) = ZD(1,2);
   pars(67) = ZD(1,3);
   pars(68) = ZD(1,4);
   pars(69) = ZD(1,5);
   pars(70) = ZD(2,0);
   pars(71) = ZD(2,1);
   pars(72) = ZD(2,2);
   pars(73) = ZD(2,3);
   pars(74) = ZD(2,4);
   pars(75) = ZD(2,5);
   pars(76) = ZD(3,0);
   pars(77) = ZD(3,1);
   pars(78) = ZD(3,2);
   pars(79) = ZD(3,3);
   pars(80) = ZD(3,4);
   pars(81) = ZD(3,5);
   pars(82) = ZD(4,0);
   pars(83) = ZD(4,1);
   pars(84) = ZD(4,2);
   pars(85) = ZD(4,3);
   pars(86) = ZD(4,4);
   pars(87) = ZD(4,5);
   pars(88) = ZD(5,0);
   pars(89) = ZD(5,1);
   pars(90) = ZD(5,2);
   pars(91) = ZD(5,3);
   pars(92) = ZD(5,4);
   pars(93) = ZD(5,5);
   pars(94) = ZV(0,0);
   pars(95) = ZV(0,1);
   pars(96) = ZV(0,2);
   pars(97) = ZV(0,3);
   pars(98) = ZV(0,4);
   pars(99) = ZV(0,5);
   pars(100) = ZV(1,0);
   pars(101) = ZV(1,1);
   pars(102) = ZV(1,2);
   pars(103) = ZV(1,3);
   pars(104) = ZV(1,4);
   pars(105) = ZV(1,5);
   pars(106) = ZV(2,0);
   pars(107) = ZV(2,1);
   pars(108) = ZV(2,2);
   pars(109) = ZV(2,3);
   pars(110) = ZV(2,4);
   pars(111) = ZV(2,5);
   pars(112) = ZV(3,0);
   pars(113) = ZV(3,1);
   pars(114) = ZV(3,2);
   pars(115) = ZV(3,3);
   pars(116) = ZV(3,4);
   pars(117) = ZV(3,5);
   pars(118) = ZV(4,0);
   pars(119) = ZV(4,1);
   pars(120) = ZV(4,2);
   pars(121) = ZV(4,3);
   pars(122) = ZV(4,4);
   pars(123) = ZV(4,5);
   pars(124) = ZV(5,0);
   pars(125) = ZV(5,1);
   pars(126) = ZV(5,2);
   pars(127) = ZV(5,3);
   pars(128) = ZV(5,4);
   pars(129) = ZV(5,5);
   pars(130) = ZU(0,0);
   pars(131) = ZU(0,1);
   pars(132) = ZU(0,2);
   pars(133) = ZU(0,3);
   pars(134) = ZU(0,4);
   pars(135) = ZU(0,5);
   pars(136) = ZU(1,0);
   pars(137) = ZU(1,1);
   pars(138) = ZU(1,2);
   pars(139) = ZU(1,3);
   pars(140) = ZU(1,4);
   pars(141) = ZU(1,5);
   pars(142) = ZU(2,0);
   pars(143) = ZU(2,1);
   pars(144) = ZU(2,2);
   pars(145) = ZU(2,3);
   pars(146) = ZU(2,4);
   pars(147) = ZU(2,5);
   pars(148) = ZU(3,0);
   pars(149) = ZU(3,1);
   pars(150) = ZU(3,2);
   pars(151) = ZU(3,3);
   pars(152) = ZU(3,4);
   pars(153) = ZU(3,5);
   pars(154) = ZU(4,0);
   pars(155) = ZU(4,1);
   pars(156) = ZU(4,2);
   pars(157) = ZU(4,3);
   pars(158) = ZU(4,4);
   pars(159) = ZU(4,5);
   pars(160) = ZU(5,0);
   pars(161) = ZU(5,1);
   pars(162) = ZU(5,2);
   pars(163) = ZU(5,3);
   pars(164) = ZU(5,4);
   pars(165) = ZU(5,5);
   pars(166) = ZE(0,0);
   pars(167) = ZE(0,1);
   pars(168) = ZE(0,2);
   pars(169) = ZE(0,3);
   pars(170) = ZE(0,4);
   pars(171) = ZE(0,5);
   pars(172) = ZE(1,0);
   pars(173) = ZE(1,1);
   pars(174) = ZE(1,2);
   pars(175) = ZE(1,3);
   pars(176) = ZE(1,4);
   pars(177) = ZE(1,5);
   pars(178) = ZE(2,0);
   pars(179) = ZE(2,1);
   pars(180) = ZE(2,2);
   pars(181) = ZE(2,3);
   pars(182) = ZE(2,4);
   pars(183) = ZE(2,5);
   pars(184) = ZE(3,0);
   pars(185) = ZE(3,1);
   pars(186) = ZE(3,2);
   pars(187) = ZE(3,3);
   pars(188) = ZE(3,4);
   pars(189) = ZE(3,5);
   pars(190) = ZE(4,0);
   pars(191) = ZE(4,1);
   pars(192) = ZE(4,2);
   pars(193) = ZE(4,3);
   pars(194) = ZE(4,4);
   pars(195) = ZE(4,5);
   pars(196) = ZE(5,0);
   pars(197) = ZE(5,1);
   pars(198) = ZE(5,2);
   pars(199) = ZE(5,3);
   pars(200) = ZE(5,4);
   pars(201) = ZE(5,5);
   pars(202) = ZH(0,0);
   pars(203) = ZH(0,1);
   pars(204) = ZH(0,2);
   pars(205) = ZH(1,0);
   pars(206) = ZH(1,1);
   pars(207) = ZH(1,2);
   pars(208) = ZH(2,0);
   pars(209) = ZH(2,1);
   pars(210) = ZH(2,2);
   pars(211) = ZA(0,0);
   pars(212) = ZA(0,1);
   pars(213) = ZA(0,2);
   pars(214) = ZA(1,0);
   pars(215) = ZA(1,1);
   pars(216) = ZA(1,2);
   pars(217) = ZA(2,0);
   pars(218) = ZA(2,1);
   pars(219) = ZA(2,2);
   pars(220) = ZP(0,0);
   pars(221) = ZP(0,1);
   pars(222) = ZP(1,0);
   pars(223) = ZP(1,1);
   pars(224) = Re(ZN(0,0));
   pars(225) = Im(ZN(0,0));
   pars(226) = Re(ZN(0,1));
   pars(227) = Im(ZN(0,1));
   pars(228) = Re(ZN(0,2));
   pars(229) = Im(ZN(0,2));
   pars(230) = Re(ZN(0,3));
   pars(231) = Im(ZN(0,3));
   pars(232) = Re(ZN(0,4));
   pars(233) = Im(ZN(0,4));
   pars(234) = Re(ZN(0,5));
   pars(235) = Im(ZN(0,5));
   pars(236) = Re(ZN(1,0));
   pars(237) = Im(ZN(1,0));
   pars(238) = Re(ZN(1,1));
   pars(239) = Im(ZN(1,1));
   pars(240) = Re(ZN(1,2));
   pars(241) = Im(ZN(1,2));
   pars(242) = Re(ZN(1,3));
   pars(243) = Im(ZN(1,3));
   pars(244) = Re(ZN(1,4));
   pars(245) = Im(ZN(1,4));
   pars(246) = Re(ZN(1,5));
   pars(247) = Im(ZN(1,5));
   pars(248) = Re(ZN(2,0));
   pars(249) = Im(ZN(2,0));
   pars(250) = Re(ZN(2,1));
   pars(251) = Im(ZN(2,1));
   pars(252) = Re(ZN(2,2));
   pars(253) = Im(ZN(2,2));
   pars(254) = Re(ZN(2,3));
   pars(255) = Im(ZN(2,3));
   pars(256) = Re(ZN(2,4));
   pars(257) = Im(ZN(2,4));
   pars(258) = Re(ZN(2,5));
   pars(259) = Im(ZN(2,5));
   pars(260) = Re(ZN(3,0));
   pars(261) = Im(ZN(3,0));
   pars(262) = Re(ZN(3,1));
   pars(263) = Im(ZN(3,1));
   pars(264) = Re(ZN(3,2));
   pars(265) = Im(ZN(3,2));
   pars(266) = Re(ZN(3,3));
   pars(267) = Im(ZN(3,3));
   pars(268) = Re(ZN(3,4));
   pars(269) = Im(ZN(3,4));
   pars(270) = Re(ZN(3,5));
   pars(271) = Im(ZN(3,5));
   pars(272) = Re(ZN(4,0));
   pars(273) = Im(ZN(4,0));
   pars(274) = Re(ZN(4,1));
   pars(275) = Im(ZN(4,1));
   pars(276) = Re(ZN(4,2));
   pars(277) = Im(ZN(4,2));
   pars(278) = Re(ZN(4,3));
   pars(279) = Im(ZN(4,3));
   pars(280) = Re(ZN(4,4));
   pars(281) = Im(ZN(4,4));
   pars(282) = Re(ZN(4,5));
   pars(283) = Im(ZN(4,5));
   pars(284) = Re(ZN(5,0));
   pars(285) = Im(ZN(5,0));
   pars(286) = Re(ZN(5,1));
   pars(287) = Im(ZN(5,1));
   pars(288) = Re(ZN(5,2));
   pars(289) = Im(ZN(5,2));
   pars(290) = Re(ZN(5,3));
   pars(291) = Im(ZN(5,3));
   pars(292) = Re(ZN(5,4));
   pars(293) = Im(ZN(5,4));
   pars(294) = Re(ZN(5,5));
   pars(295) = Im(ZN(5,5));
   pars(296) = Re(ZVL(0,0));
   pars(297) = Im(ZVL(0,0));
   pars(298) = Re(ZVL(0,1));
   pars(299) = Im(ZVL(0,1));
   pars(300) = Re(ZVL(0,2));
   pars(301) = Im(ZVL(0,2));
   pars(302) = Re(ZVL(1,0));
   pars(303) = Im(ZVL(1,0));
   pars(304) = Re(ZVL(1,1));
   pars(305) = Im(ZVL(1,1));
   pars(306) = Re(ZVL(1,2));
   pars(307) = Im(ZVL(1,2));
   pars(308) = Re(ZVL(2,0));
   pars(309) = Im(ZVL(2,0));
   pars(310) = Re(ZVL(2,1));
   pars(311) = Im(ZVL(2,1));
   pars(312) = Re(ZVL(2,2));
   pars(313) = Im(ZVL(2,2));
   pars(314) = Re(ZVR(0,0));
   pars(315) = Im(ZVR(0,0));
   pars(316) = Re(ZVR(0,1));
   pars(317) = Im(ZVR(0,1));
   pars(318) = Re(ZVR(0,2));
   pars(319) = Im(ZVR(0,2));
   pars(320) = Re(ZVR(1,0));
   pars(321) = Im(ZVR(1,0));
   pars(322) = Re(ZVR(1,1));
   pars(323) = Im(ZVR(1,1));
   pars(324) = Re(ZVR(1,2));
   pars(325) = Im(ZVR(1,2));
   pars(326) = Re(ZVR(2,0));
   pars(327) = Im(ZVR(2,0));
   pars(328) = Re(ZVR(2,1));
   pars(329) = Im(ZVR(2,1));
   pars(330) = Re(ZVR(2,2));
   pars(331) = Im(ZVR(2,2));
   pars(332) = Re(UM(0,0));
   pars(333) = Im(UM(0,0));
   pars(334) = Re(UM(0,1));
   pars(335) = Im(UM(0,1));
   pars(336) = Re(UM(1,0));
   pars(337) = Im(UM(1,0));
   pars(338) = Re(UM(1,1));
   pars(339) = Im(UM(1,1));
   pars(340) = Re(UP(0,0));
   pars(341) = Im(UP(0,0));
   pars(342) = Re(UP(0,1));
   pars(343) = Im(UP(0,1));
   pars(344) = Re(UP(1,0));
   pars(345) = Im(UP(1,0));
   pars(346) = Re(UP(1,1));
   pars(347) = Im(UP(1,1));
   pars(348) = Re(ZEL(0,0));
   pars(349) = Im(ZEL(0,0));
   pars(350) = Re(ZEL(0,1));
   pars(351) = Im(ZEL(0,1));
   pars(352) = Re(ZEL(0,2));
   pars(353) = Im(ZEL(0,2));
   pars(354) = Re(ZEL(1,0));
   pars(355) = Im(ZEL(1,0));
   pars(356) = Re(ZEL(1,1));
   pars(357) = Im(ZEL(1,1));
   pars(358) = Re(ZEL(1,2));
   pars(359) = Im(ZEL(1,2));
   pars(360) = Re(ZEL(2,0));
   pars(361) = Im(ZEL(2,0));
   pars(362) = Re(ZEL(2,1));
   pars(363) = Im(ZEL(2,1));
   pars(364) = Re(ZEL(2,2));
   pars(365) = Im(ZEL(2,2));
   pars(366) = Re(ZER(0,0));
   pars(367) = Im(ZER(0,0));
   pars(368) = Re(ZER(0,1));
   pars(369) = Im(ZER(0,1));
   pars(370) = Re(ZER(0,2));
   pars(371) = Im(ZER(0,2));
   pars(372) = Re(ZER(1,0));
   pars(373) = Im(ZER(1,0));
   pars(374) = Re(ZER(1,1));
   pars(375) = Im(ZER(1,1));
   pars(376) = Re(ZER(1,2));
   pars(377) = Im(ZER(1,2));
   pars(378) = Re(ZER(2,0));
   pars(379) = Im(ZER(2,0));
   pars(380) = Re(ZER(2,1));
   pars(381) = Im(ZER(2,1));
   pars(382) = Re(ZER(2,2));
   pars(383) = Im(ZER(2,2));
   pars(384) = Re(ZDL(0,0));
   pars(385) = Im(ZDL(0,0));
   pars(386) = Re(ZDL(0,1));
   pars(387) = Im(ZDL(0,1));
   pars(388) = Re(ZDL(0,2));
   pars(389) = Im(ZDL(0,2));
   pars(390) = Re(ZDL(1,0));
   pars(391) = Im(ZDL(1,0));
   pars(392) = Re(ZDL(1,1));
   pars(393) = Im(ZDL(1,1));
   pars(394) = Re(ZDL(1,2));
   pars(395) = Im(ZDL(1,2));
   pars(396) = Re(ZDL(2,0));
   pars(397) = Im(ZDL(2,0));
   pars(398) = Re(ZDL(2,1));
   pars(399) = Im(ZDL(2,1));
   pars(400) = Re(ZDL(2,2));
   pars(401) = Im(ZDL(2,2));
   pars(402) = Re(ZDR(0,0));
   pars(403) = Im(ZDR(0,0));
   pars(404) = Re(ZDR(0,1));
   pars(405) = Im(ZDR(0,1));
   pars(406) = Re(ZDR(0,2));
   pars(407) = Im(ZDR(0,2));
   pars(408) = Re(ZDR(1,0));
   pars(409) = Im(ZDR(1,0));
   pars(410) = Re(ZDR(1,1));
   pars(411) = Im(ZDR(1,1));
   pars(412) = Re(ZDR(1,2));
   pars(413) = Im(ZDR(1,2));
   pars(414) = Re(ZDR(2,0));
   pars(415) = Im(ZDR(2,0));
   pars(416) = Re(ZDR(2,1));
   pars(417) = Im(ZDR(2,1));
   pars(418) = Re(ZDR(2,2));
   pars(419) = Im(ZDR(2,2));
   pars(420) = Re(ZUL(0,0));
   pars(421) = Im(ZUL(0,0));
   pars(422) = Re(ZUL(0,1));
   pars(423) = Im(ZUL(0,1));
   pars(424) = Re(ZUL(0,2));
   pars(425) = Im(ZUL(0,2));
   pars(426) = Re(ZUL(1,0));
   pars(427) = Im(ZUL(1,0));
   pars(428) = Re(ZUL(1,1));
   pars(429) = Im(ZUL(1,1));
   pars(430) = Re(ZUL(1,2));
   pars(431) = Im(ZUL(1,2));
   pars(432) = Re(ZUL(2,0));
   pars(433) = Im(ZUL(2,0));
   pars(434) = Re(ZUL(2,1));
   pars(435) = Im(ZUL(2,1));
   pars(436) = Re(ZUL(2,2));
   pars(437) = Im(ZUL(2,2));
   pars(438) = Re(ZUR(0,0));
   pars(439) = Im(ZUR(0,0));
   pars(440) = Re(ZUR(0,1));
   pars(441) = Im(ZUR(0,1));
   pars(442) = Re(ZUR(0,2));
   pars(443) = Im(ZUR(0,2));
   pars(444) = Re(ZUR(1,0));
   pars(445) = Im(ZUR(1,0));
   pars(446) = Re(ZUR(1,1));
   pars(447) = Im(ZUR(1,1));
   pars(448) = Re(ZUR(1,2));
   pars(449) = Im(ZUR(1,2));
   pars(450) = Re(ZUR(2,0));
   pars(451) = Im(ZUR(2,0));
   pars(452) = Re(ZUR(2,1));
   pars(453) = Im(ZUR(2,1));
   pars(454) = Re(ZUR(2,2));
   pars(455) = Im(ZUR(2,2));
   pars(456) = ZZ(0,0);
   pars(457) = ZZ(0,1);
   pars(458) = ZZ(0,2);
   pars(459) = ZZ(1,0);
   pars(460) = ZZ(1,1);
   pars(461) = ZZ(1,2);
   pars(462) = ZZ(2,0);
   pars(463) = ZZ(2,1);
   pars(464) = ZZ(2,2);


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
   return "UMSSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   UMSSM_soft_parameters::run_to(scale, eps);
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

const UMSSM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

UMSSM_input_parameters& CLASSNAME::get_input_parameters()
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
double CLASSNAME::get_lsp(UMSSM_info::Particles& particle_type) const
{
   double lsp_mass = std::numeric_limits<double>::max();
   double tmp_mass;
   particle_type = UMSSM_info::NUMBER_OF_PARTICLES;

   tmp_mass = Abs(PHYSICAL(MChi(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Chi;
   }

   tmp_mass = Abs(PHYSICAL(MSv(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Sv;
   }

   tmp_mass = Abs(PHYSICAL(MSu(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Su;
   }

   tmp_mass = Abs(PHYSICAL(MSd(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Sd;
   }

   tmp_mass = Abs(PHYSICAL(MSe(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Se;
   }

   tmp_mass = Abs(PHYSICAL(MCha(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Cha;
   }

   tmp_mass = Abs(PHYSICAL(MGlu));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = UMSSM_info::Glu;
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

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qd = LOCALINPUT(Qd);

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)) +
      AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) + Conj(
      Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) + Conj(
      Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) - 0.5*vS*vu*Conj(
      Yd(0,0))*Lambdax;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) - 0.5*vS*vu*Conj(
      Yd(1,0))*Lambdax;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) - 0.5*vS*vu*Conj(
      Yd(2,0))*Lambdax;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) + Conj(
      Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) - 0.5*vS*vu*Conj(
      Yd(0,1))*Lambdax;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) - 0.5*vS*vu*Conj(
      Yd(1,1))*Lambdax;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) - 0.5*vS*vu*Conj(
      Yd(2,1))*Lambdax;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) - 0.5*vS*vu*Conj(
      Yd(0,2))*Lambdax;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) - 0.5*vS*vu*Conj(
      Yd(1,2))*Lambdax;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vS*vu*Conj(
      Yd(2,2))*Lambdax;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)) +
      AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) + Conj(
      Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) + Conj(
      Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) + Conj(
      Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*Sqr(
      gp)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Sd, eigenvalue_error > precision * Abs(
      MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sv() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qv = LOCALINPUT(Qv);

   Eigen::Matrix<double,6,6> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(0,0)) + AbsSqr(Yv(0,1)) + AbsSqr(Yv(0,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(1,0) + Conj(
      Yv(0,1))*Yv(1,1) + Conj(Yv(0,2))*Yv(1,2));
   mass_matrix_Sv(0,2) = ml2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(2,0) + Conj(
      Yv(0,1))*Yv(2,1) + Conj(Yv(0,2))*Yv(2,2));
   mass_matrix_Sv(0,3) = 0.7071067811865475*vu*Conj(TYv(0,0)) - 0.5*vd*vS*Conj(
      Yv(0,0))*Lambdax;
   mass_matrix_Sv(0,4) = 0.7071067811865475*vu*Conj(TYv(0,1)) - 0.5*vd*vS*Conj(
      Yv(0,1))*Lambdax;
   mass_matrix_Sv(0,5) = 0.7071067811865475*vu*Conj(TYv(0,2)) - 0.5*vd*vS*Conj(
      Yv(0,2))*Lambdax;
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(1,0)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(1,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(1,0))*Yv(2,0) + Conj(
      Yv(1,1))*Yv(2,1) + Conj(Yv(1,2))*Yv(2,2));
   mass_matrix_Sv(1,3) = 0.7071067811865475*vu*Conj(TYv(1,0)) - 0.5*vd*vS*Conj(
      Yv(1,0))*Lambdax;
   mass_matrix_Sv(1,4) = 0.7071067811865475*vu*Conj(TYv(1,1)) - 0.5*vd*vS*Conj(
      Yv(1,1))*Lambdax;
   mass_matrix_Sv(1,5) = 0.7071067811865475*vu*Conj(TYv(1,2)) - 0.5*vd*vS*Conj(
      Yv(1,2))*Lambdax;
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(2,0)) + AbsSqr(Yv(2,1)) + AbsSqr(Yv(2,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(2,3) = 0.7071067811865475*vu*Conj(TYv(2,0)) - 0.5*vd*vS*Conj(
      Yv(2,0))*Lambdax;
   mass_matrix_Sv(2,4) = 0.7071067811865475*vu*Conj(TYv(2,1)) - 0.5*vd*vS*Conj(
      Yv(2,1))*Lambdax;
   mass_matrix_Sv(2,5) = 0.7071067811865475*vu*Conj(TYv(2,2)) - 0.5*vd*vS*Conj(
      Yv(2,2))*Lambdax;
   mass_matrix_Sv(3,3) = mvR2(0,0) + 0.5*QHd*Qv*Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr
      (gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,0)) + AbsSqr(Yv(1,0)) + AbsSqr(Yv(2,0)))*
      Sqr(vu) + 0.5*QHu*Qv*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(3,4) = mvR2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,1))*Yv(0,0) + Conj(
      Yv(1,1))*Yv(1,0) + Conj(Yv(2,1))*Yv(2,0));
   mass_matrix_Sv(3,5) = mvR2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,2))*Yv(0,0) + Conj(
      Yv(1,2))*Yv(1,0) + Conj(Yv(2,2))*Yv(2,0));
   mass_matrix_Sv(4,4) = mvR2(1,1) + 0.5*QHd*Qv*Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr
      (gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,1)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(2,1)))*
      Sqr(vu) + 0.5*QHu*Qv*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(4,5) = mvR2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(0,2))*Yv(0,1) + Conj(
      Yv(1,2))*Yv(1,1) + Conj(Yv(2,2))*Yv(2,1));
   mass_matrix_Sv(5,5) = mvR2(2,2) + 0.5*QHd*Qv*Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr
      (gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,2)) + AbsSqr(Yv(1,2)) + AbsSqr(Yv(2,2)))*
      Sqr(vu) + 0.5*QHu*Qv*Sqr(gp)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Sv, eigenvalue_error > precision * Abs(
      MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qu = LOCALINPUT(Qu);

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) + Conj(
      Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) + Conj(
      Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) - 0.5*vd*vS*Conj(
      Yu(0,0))*Lambdax;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) - 0.5*vd*vS*Conj(
      Yu(1,0))*Lambdax;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) - 0.5*vd*vS*Conj(
      Yu(2,0))*Lambdax;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) + Conj(
      Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) - 0.5*vd*vS*Conj(
      Yu(0,1))*Lambdax;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vS*Conj(
      Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) - 0.5*vd*vS*Conj(
      Yu(2,1))*Lambdax;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) - 0.5*vd*vS*Conj(
      Yu(0,2))*Lambdax;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) - 0.5*vd*vS*Conj(
      Yu(1,2))*Lambdax;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vS*Conj(
      Yu(2,2))*Lambdax;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(gp)*
      Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(0,
      1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*Sqr(gp)
      *Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) + Conj(
      Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) + Conj(
      Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(gp)*
      Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(1,0)) + AbsSqr(Yu(1,
      1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*Sqr(gp)
      *Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) + Conj(
      Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(gp)*
      Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(2,0)) + AbsSqr(Yu(2,
      1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*Sqr(gp)
      *Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Su, eigenvalue_error > precision * Abs(
      MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qe = LOCALINPUT(Qe);

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)) +
      AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) + Conj(
      Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) + Conj(
      Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) - 0.5*vS*vu*Conj(
      Ye(0,0))*Lambdax;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) - 0.5*vS*vu*Conj(
      Ye(1,0))*Lambdax;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) - 0.5*vS*vu*Conj(
      Ye(2,0))*Lambdax;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) + Conj(
      Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) - 0.5*vS*vu*Conj(
      Ye(0,1))*Lambdax;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) - 0.5*vS*vu*Conj(
      Ye(1,1))*Lambdax;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) - 0.5*vS*vu*Conj(
      Ye(2,1))*Lambdax;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*
      Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) - 0.5*vS*vu*Conj(
      Ye(0,2))*Lambdax;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) - 0.5*vS*vu*Conj(
      Ye(1,2))*Lambdax;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vS*vu*Conj(
      Ye(2,2))*Lambdax;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)) +
      AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) + Conj(
      Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) + Conj(
      Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*Sqr(
      gp)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) + Conj(
      Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*Sqr(
      vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*Sqr(
      gp)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Se, eigenvalue_error > precision * Abs(
      MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      1.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*
      Sqr(gp)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vS*Conj(
      TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) + QHd*QHu*vd*vu*Sqr(
      gp) - 0.35355339059327373*vS*TLambdax;
   mass_matrix_hh(0,2) = vd*vS*AbsSqr(Lambdax) - 0.35355339059327373*vu*Conj(
      TLambdax) + QHd*Qs*vd*vS*Sqr(gp) - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(
      Lambdax)*Sqr(vS) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vS) + 0.225*Sqr(g1)*Sqr(vu) +
      0.375*Sqr(g2)*Sqr(vu) + 1.5*Sqr(gp)*Sqr(QHu)*Sqr(vu);
   mass_matrix_hh(1,2) = vS*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*Conj(
      TLambdax) + QHu*Qs*vS*vu*Sqr(gp) - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*Sqr(gp)
      *Sqr(vd) + 1.5*Sqr(gp)*Sqr(Qs)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) +
      0.5*QHu*Qs*Sqr(gp)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::hh, eigenvalue_error > precision * Abs(
      Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) +
      0.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*
      Sqr(gp)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu) + 0.3872983346207417*
      g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vd)*Sqr(Cos(ThetaWp())) + Sqr(gp)*
      Sqr(QHd)*Sqr(vd)*Sqr(Cos(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Cos(ThetaWp()))
      *Sqr(Sin(ThetaW())) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW()
      )*Sqr(vd)*Sqr(Sin(ThetaWp())) + Sqr(gp)*Sqr(QHd)*Sqr(vd)*Sqr(Sin(ThetaWp(
      ))) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*
      Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vS*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW())*Sqr(Cos(
      ThetaWp())) + QHd*QHu*vd*vu*Sqr(gp)*Sqr(Cos(ThetaWp())) - 0.25*vd*vu*Sqr(
      g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 0.15*vd*vu*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) - 0.3872983346207417*g1*g2*vd*vu*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + QHd*QHu*vd*vu*Sqr(gp)*Sqr(
      Sin(ThetaWp())) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()
      )) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vS*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) + QHd*Qs*vd*vS*
      Sqr(gp)*Sqr(Cos(ThetaWp())) + QHd*Qs*vd*vS*Sqr(gp)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(
      Lambdax)*Sqr(vS) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu) + 0.5*Sqr(gp)*Sqr(QHu)*Sqr(vu) + 0.3872983346207417
      *g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu)*Sqr(Cos(ThetaWp())) + Sqr(gp)*
      Sqr(QHu)*Sqr(vu)*Sqr(Cos(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Cos(ThetaWp()))
      *Sqr(Sin(ThetaW())) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW()
      )*Sqr(vu)*Sqr(Sin(ThetaWp())) + Sqr(gp)*Sqr(QHu)*Sqr(vu)*Sqr(Sin(ThetaWp(
      ))) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + 0.15*
      Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) + QHu*Qs*vS*vu*
      Sqr(gp)*Sqr(Cos(ThetaWp())) + QHu*Qs*vS*vu*Sqr(gp)*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vd*TLambdax;
   mass_matrix_Ah(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*Sqr(gp)
      *Sqr(vd) + 0.5*Sqr(gp)*Sqr(Qs)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) +
      0.5*QHu*Qs*Sqr(gp)*Sqr(vu) + Sqr(gp)*Sqr(Qs)*Sqr(vS)*Sqr(Cos(ThetaWp()))
      + Sqr(gp)*Sqr(Qs)*Sqr(vS)*Sqr(Sin(ThetaWp()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Ah, eigenvalue_error > precision * Abs(
      MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd)
      + 0.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs
      *Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*
      QHd*QHu*Sqr(gp)*Sqr(vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475*vS*
      Conj(TLambdax);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHu*Qs*
      Sqr(gp)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu) + 0.5*Sqr
      (gp)*Sqr(QHu)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Hpm, eigenvalue_error > precision * Abs(
      MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(UMSSM_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Chi() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   Eigen::Matrix<double,6,6> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassU;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = 0;
   mass_matrix_Chi(0,3) = gp*QHd*vd;
   mass_matrix_Chi(0,4) = gp*QHu*vu;
   mass_matrix_Chi(0,5) = gp*Qs*vS;
   mass_matrix_Chi(1,1) = MassB;
   mass_matrix_Chi(1,2) = 0;
   mass_matrix_Chi(1,3) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(1,4) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,5) = 0;
   mass_matrix_Chi(2,2) = MassWB;
   mass_matrix_Chi(2,3) = 0.5*g2*vd;
   mass_matrix_Chi(2,4) = -0.5*g2*vu;
   mass_matrix_Chi(2,5) = 0;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.7071067811865475*vS*Lambdax;
   mass_matrix_Chi(3,5) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Chi(4,4) = 0;
   mass_matrix_Chi(4,5) = -0.7071067811865475*vd*Lambdax;
   mass_matrix_Chi(5,5) = 0;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Chi, eigenvalue_error > precision * Abs(
      MChi(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
   normalize_to_interval(ZN);

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0.7071067811865475*vu*Yv(0,0);
   mass_matrix_Fv(0,1) = 0.7071067811865475*vu*Yv(0,1);
   mass_matrix_Fv(0,2) = 0.7071067811865475*vu*Yv(0,2);
   mass_matrix_Fv(1,0) = 0.7071067811865475*vu*Yv(1,0);
   mass_matrix_Fv(1,1) = 0.7071067811865475*vu*Yv(1,1);
   mass_matrix_Fv(1,2) = 0.7071067811865475*vu*Yv(1,2);
   mass_matrix_Fv(2,0) = 0.7071067811865475*vu*Yv(2,0);
   mass_matrix_Fv(2,1) = 0.7071067811865475*vu*Yv(2,1);
   mass_matrix_Fv(2,2) = 0.7071067811865475*vu*Yv(2,2);

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   const auto mass_matrix_Fv(get_mass_matrix_Fv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fv, MFv, ZVL, ZVR, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Fv, eigenvalue_error > precision * Abs(
      MFv(0)));
#else
   fs_svd(mass_matrix_Fv, MFv, ZVL, ZVR);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{

   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = 0.7071067811865475*vS*Lambdax;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Cha, eigenvalue_error > precision * Abs(
      MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*vd*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*vd*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*vd*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*vd*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*vd*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*vd*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*vd*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*vd*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*vd*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Fe, eigenvalue_error > precision * Abs(
      MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, ZEL, ZER);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*vd*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*vd*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*vd*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*vd*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*vd*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*vd*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*vd*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*vd*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*vd*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Fd, eigenvalue_error > precision * Abs(
      MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, ZDL, ZDR);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*vu*Yu(0,0);
   mass_matrix_Fu(0,1) = 0.7071067811865475*vu*Yu(1,0);
   mass_matrix_Fu(0,2) = 0.7071067811865475*vu*Yu(2,0);
   mass_matrix_Fu(1,0) = 0.7071067811865475*vu*Yu(0,1);
   mass_matrix_Fu(1,1) = 0.7071067811865475*vu*Yu(1,1);
   mass_matrix_Fu(1,2) = 0.7071067811865475*vu*Yu(2,1);
   mass_matrix_Fu(2,0) = 0.7071067811865475*vu*Yu(0,2);
   mass_matrix_Fu(2,1) = 0.7071067811865475*vu*Yu(1,2);
   mass_matrix_Fu(2,2) = 0.7071067811865475*vu*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Fu, eigenvalue_error > precision * Abs(
      MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

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
      problems.flag_running_tachyon(UMSSM_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_VPVZVZp() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   Eigen::Matrix<double,3,3> mass_matrix_VPVZVZp;

   mass_matrix_VPVZVZp(0,0) = 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_VPVZVZp(0,1) = -0.19364916731037085*g1*g2*Sqr(vd) -
      0.19364916731037085*g1*g2*Sqr(vu);
   mass_matrix_VPVZVZp(0,2) = -0.3872983346207417*g1*gp*QHd*Sqr(vd) +
      0.3872983346207417*g1*gp*QHu*Sqr(vu);
   mass_matrix_VPVZVZp(1,1) = 0.25*Sqr(g2)*Sqr(vd) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_VPVZVZp(1,2) = 0.5*g2*gp*QHd*Sqr(vd) - 0.5*g2*gp*QHu*Sqr(vu);
   mass_matrix_VPVZVZp(2,2) = Sqr(gp)*Sqr(QHd)*Sqr(vd) + Sqr(gp)*Sqr(Qs)*Sqr(vS
      ) + Sqr(gp)*Sqr(QHu)*Sqr(vu);

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
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   double result = Re(mHd2*vd - 0.35355339059327373*vS*vu*Conj(TLambdax) + 0.075*
      Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.5*Cube(vd)*Sqr(gp)*Sqr(QHd) +
      0.5*vd*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*vd*Sqr(gp)*Sqr(vS) + 0.5*vd*
      AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu
      ) + 0.5*QHd*QHu*vd*Sqr(gp)*Sqr(vu) - 0.35355339059327373*vS*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   double result = Re(mHu2*vu - 0.35355339059327373*vd*vS*Conj(TLambdax) + 0.075*
      Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) + 0.5*Cube(vu)*Sqr(gp)*Sqr(QHu) +
      0.5*vu*AbsSqr(Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*QHu*vu*Sqr(gp)*Sqr(vd) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vS) +
      0.5*QHu*Qs*vu*Sqr(gp)*Sqr(vS) - 0.35355339059327373*vd*vS*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   double result = Re(ms2*vS - 0.35355339059327373*vd*vu*Conj(TLambdax) + 0.5*Cube
      (vS)*Sqr(gp)*Sqr(Qs) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*vS*Sqr(gp
      )*Sqr(vd) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vu) + 0.5*QHu*Qs*vS*Sqr(gp)*Sqr(vu) -
      0.35355339059327373*vd*vu*TLambdax);

   return result;
}



std::complex<double> CLASSNAME::CpUSdconjUSdVZVZ(int gO1, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-2*g2*gp*Qq*Cos(ThetaW())*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1 < 3,-
      0.5163977794943222*g1*gp*Qq*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,0.2581988897471611*g1*g2*Cos(ThetaW
      ())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(gO1 <
      3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()))
      ,0) + IF(gO1 < 3,0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos
      (ThetaWp()))*Sqr(Sin(ThetaW())),0) + IF(gO1 < 3,2*KroneckerDelta(gO1,gO2)*
      Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())),0) - 1.0327955589886444*g1*gp*Qd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) + 0.13333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp(
      )))*Sqr(Sin(ThetaW()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(
      gO2,3 + j1)) + 2*Sqr(gp)*Sqr(Qd)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdVZpVZp(int gO1, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,g2*gp*Qq*Cos(ThetaW())*
      KroneckerDelta(gO1,gO2)*Sin(2*ThetaWp()),0) + IF(gO1 < 3,0.2581988897471611*
      g1*gp*Qq*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(2*ThetaWp()),0) + IF(gO1
      < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())),0) + IF(
      gO1 < 3,0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,
      0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
      Sin(ThetaWp())),0) + 1.0327955589886444*g1*gp*Qd*Cos(ThetaWp())*Sin(ThetaW()
      )*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
      + j1)) + 2*Sqr(gp)*Sqr(Qd)*Sqr(Cos(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.13333333333333333*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSdconjUSdconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSdconjHpmconjUSd(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(
      j1,gO2))*Yu(j1,gO1))*ZP(gI1,1)*ZP(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-(QHd*Qq
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(QHu*Qq
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)),0) + 0.1*Sqr(g1)*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(
      gI2,0) - Qd*QHd*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*
      SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1)))
      )*ZP(gI1,0)*ZP(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) - Qd*QHu*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZA(gI1,0))*
      Conj(ZA(gI2,0))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,gO1))),0),0) + IF(gO1 < 3,
      0.05*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF
      (gO1 < 3,-0.05*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,
      gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Qq*Conj(ZA(gI1,0))*
      Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Qq*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Qq*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*SUM(j1,
      0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1)),0) + IF(gO1 < 3,-0.5*Conj(Lambdax
      )*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(
      j1,gO1)),0) + IF(gO2 < 3,-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*Lambdax*SUM(j1
      ,0,2,Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1)),0) + IF(gO2 < 3,-0.5*Conj(
      ZA(gI1,1))*Conj(ZA(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*
      KroneckerDelta(gO1,3 + j1)),0) + 0.1*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.1*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) - Qd*QHd*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(
      gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Qd*
      QHu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1)) - Qd*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*
      Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) -
      Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2
      ,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZH(gI1,0))*
      Conj(ZH(gI2,0))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,gO1))),0),0) + IF(gO1 < 3,
      0.05*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF
      (gO1 < 3,-0.05*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,
      gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Qq*Conj(ZH(gI1,0))*
      Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Qq*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Qq*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))*SUM(j1,0
      ,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1)),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,
      gO1)),0) + IF(gO2 < 3,0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))*Lambdax*SUM(j1,0,2
      ,Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1)),0) + IF(gO2 < 3,0.5*Conj(ZH(
      gI1,1))*Conj(ZH(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*KroneckerDelta(
      gO1,3 + j1)),0) + 0.1*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.1*Conj(ZH(gI1,1))
      *Conj(ZH(gI2,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) - Qd*QHd*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Qd*QHu*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) - Qd*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Conj(
      ZH(gI1,0))*Conj(ZH(gI2,0))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpChaFuconjUSdPR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Yu(j1,gO2))*ZUR(
      gI1,j1))*UP(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpChaFuconjUSdPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UM(gI2,0))*Conj(ZUL(
      gI1,gO1))),0) + Conj(UM(gI2,1))*SUM(j2,0,2,Conj(ZUL(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPR(int gI2, int gI1, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO2))*
      ZDR(gI1,j1))*ZN(gI2,3)),0) - 0.09428090415820634*SUM(j1,0,2,KroneckerDelta(
      gO2,3 + j1)*ZDR(gI1,j1))*(15*gp*Qd*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPL(int gI2, int gI1, int gO1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Qq*Conj(
      ZDL(gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,-0.18257418583505536*g1*Conj(
      ZDL(gI1,gO1))*Conj(ZN(gI2,1)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZDL
      (gI1,gO1))*Conj(ZN(gI2,2)),0) - Conj(ZN(gI2,3))*SUM(j2,0,2,Conj(ZDL(gI1,j2))
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yd(j1,
      gO1)*ZD(gI1,3 + j1))*SUM(j3,0,2,Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*ZD
      (gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(
      gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-1.3333333333333333*Conj(ZD(gI2,gO2))
      *Sqr(g3)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZD(gI2,gO2))*Sqr(
      gp)*Sqr(Qq)*ZD(gI1,gO1)),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) +
      IF(gO1 < 3,-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-1.5*Qd*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 +
      j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.375*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-1.5*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,
      j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + IF(gO1 < 3,-1.5*Qd*Qq*KroneckerDelta(
      gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + IF(gO1
       < 3,-3*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1))*SUM(j4,0,2,SUM(j3,
      0,2,Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3)))*ZD(gI1,j4)),0) + IF(gO1 < 3,-
      0.016666666666666666*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j1,0
      ,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1
      < 3,-0.5*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3
      + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,-0.016666666666666666*Sqr(g1)*SUM(j2,0,2,
      Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO1 < 3
      ,0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(
      gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO1 < 3,-0.5*Qd*Qq*Sqr(gp)*SUM(j2,0,2,Conj(
      ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO2 < 3,-
      0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZD(gI2,gO2)
      )*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2
      < 3,-0.5*Qd*Qq*Conj(ZD(gI2,gO2))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*ZD(gI1,3 + j1)),0) + IF(gO2 < 3,-0.016666666666666666*Conj(ZD(gI2,gO2))*
      Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)),0) + IF(gO2 <
      3,0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1
      ,3 + j2)*ZD(gI1,3 + j2)),0) + IF(gO2 < 3,-0.5*Qd*Qq*Conj(ZD(gI2,gO2))*Sqr(gp
      )*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)),0) + IF(gO2 < 3,-3*
      SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1)))*SUM(j3,0,2
      ,Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3)),0) - 0.03333333333333333*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD
      (gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.6666666666666666*Sqr(g3)*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 +
      j2))*KroneckerDelta(gO2,3 + j2)) - 0.5*Sqr(gp)*Sqr(Qd)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*
      KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(
      gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) -
      1.5*Qd*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 1.5*Sqr(gp)*Sqr(Qd)*SUM(j1,0,2,Conj(ZD(gI2
      ,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) -
      1.5*Qd*Qq*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
       + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*ZD(gI1,3 + j2)) - 1.5*Sqr(gp)*Sqr(Qd)*SUM(j1,0,2,KroneckerDelta
      (gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(
      gI1,3 + j2)) - 0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 +
      j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 +
      j2)) - 0.5*Sqr(gp)*Sqr(Qd)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)) - SUM(j2,
      0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*SUM(
      j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZD(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yu(j1,
      gO1)*ZU(gI1,3 + j1))*SUM(j3,0,2,Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1),0),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1)
      )*ZU(gI1,j1)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0
      ,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,-1.5*KroneckerDelta(gO1,gO2)
      *Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.1
      *KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 +
      j1)),0) + IF(gO1 < 3,-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,
      Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(
      gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,
      0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)
      ),0) + IF(gO1 < 3,-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*SUM(j2,0,2,
      Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,-1.5*Qq
      *Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3
       + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2
      ,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 1.5*Qd*Qq*Sqr(gp)*
      SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2
      )*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      ZU(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 1.5*Qd*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM
      (j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0
      ,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) - 1.5*Qd*Qq*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,j2))*ZU(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2))
      - 1.5*Qd*Qu*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2
      ,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)) - SUM(j2,0,2,Conj(
      ZU(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*SUM(j4,0,2,SUM
      (j3,0,2,Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZU(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSdSeconjUSdconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1
      ,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-0.5*Qe*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(
      gp)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1
      ,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3
      ,j4))*Conj(ZE(gI1,3 + j3)))*ZE(gI2,j4))),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(
      ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1)))*SUM(j3,0,2,Conj(Yd(j3,gO2)
      )*KroneckerDelta(gO1,3 + j3))),0) + 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))
      *ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2
      )) - 0.5*Qd*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qd*Qe*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,3
       + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta
      (gO2,3 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 0.5*Qd
      *Ql*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.5*Qd*Qe*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(
      gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSdSvconjUSdconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1
      ,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-0.5*Qq*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)
      *SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) +
      IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZV(gI1,j2))
      *ZV(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(
      j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) + IF(gO1 < 3,-0.5*Qq*Qv*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*ZV(gI2,3 +
      j2)),0) + 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qd*Ql*Sqr(gp)*
      SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2
      )*KroneckerDelta(gO2,3 + j2)) - 0.5*Qd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 +
      j1))*ZV(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)) - 0.5*Qd
      *Ql*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)) - 0.5*Qd*Qv*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(
      gI1,3 + j2))*ZV(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpHpmSuconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZU(
      gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*vS*Lambdax*
      SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZP(gI2,0),0) + IF(gO2 < 3,
      0.7071067811865475*vd*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2)
      )*Yd(j1,j2)))*ZP(gI2,0),0) + IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZU(gI1,
      gO2))*Sqr(g2)*ZP(gI2,1),0) + IF(gO2 < 3,SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*Conj
      (TYu(j1,gO2)))*ZP(gI2,1),0) + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZP(gI2,1),0) + SUM(
      j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYd(j1,j2)))*
      ZP(gI2,0) + 0.7071067811865475*vu*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2
      ,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))))*ZP(gI2,0
      ) + 0.7071067811865475*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0
      ,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZP(gI2,1) + 0.7071067811865475*vd*
      SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1
      ,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))))*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhSdconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,0.5)*vS*
      Conj(ZA(gI2,1))*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1))),0)
      + IF(gO2 < 3,std::complex<double>(0,0.5)*vu*Conj(ZA(gI2,2))*Lambdax*SUM(j1,0
      ,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1))),0) + IF(gO2 < 3,std::complex<
      double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1
      ))*Conj(TYd(j1,gO2))),0) - std::complex<double>(0,0.5)*vS*Conj(Lambdax)*Conj
      (ZA(gI2,1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1
      )*Yd(j1,j2))) - std::complex<double>(0,0.5)*vu*Conj(Lambdax)*Conj(ZA(gI2,2))
      *SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)
      )) - std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CphhSdconjUSd(int gI2, int gI1, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO2 < 3,0.05*vd*Conj(ZD(gI1,gO2))*Conj(
      ZH(gI2,0))*Sqr(g1),0) + IF(gO2 < 3,-0.05*vu*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,1)
      )*Sqr(g1),0) + IF(gO2 < 3,0.25*vd*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(g2),
      0) + IF(gO2 < 3,-0.25*vu*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(g2),0) + IF(
      gO2 < 3,-(QHd*Qq*vd*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(gp)),0) + IF(gO2 <
      3,-(QHu*Qq*vu*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(gp)),0) + IF(gO2 < 3,-(
      Qq*Qs*vS*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,2))*Sqr(gp)),0) + IF(gO2 < 3,0.5*vS*
      Conj(ZH(gI2,1))*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1))),0)
      + IF(gO2 < 3,0.5*vu*Conj(ZH(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj
      (ZD(gI1,3 + j1))),0) + IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI2,0))*SUM(j1
      ,0,2,Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2))),0) + IF(gO2 < 3,-(vd*Conj(ZH(
      gI2,0))*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,j2))))
      ,0) + 0.1*vd*Conj(ZH(gI2,0))*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) - 0.1*vu*Conj(ZH(gI2,1))*Sqr(g1)*SUM(j1,0,2,Conj
      (ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - Qd*QHd*vd*Conj(ZH(gI2,0))*Sqr
      (gp)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - Qd*QHu*vu
      *Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,
      3 + j1)) - Qd*Qs*vS*Conj(ZH(gI2,2))*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) + 0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,1))*SUM(j2,0,
      2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2))) + 0.5*
      vu*Conj(Lambdax)*Conj(ZH(gI2,2))*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j2))) - 0.7071067811865475*Conj(ZH(gI2,0))*
      SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYd(j1,j2)
      )) - vd*Conj(ZH(gI2,0))*SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))));

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
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Conj(ZD(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZD(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,gp*Qq*Conj(ZD(gI2,gO2))*
      Sin(ThetaWp()),0) + 0.2581988897471611*g1*Cos(ThetaWp())*Sin(ThetaW())*SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) - gp*Qd*Sin(ThetaWp(
      ))*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVZp(int gI2, int gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO2 < 3,gp*Qq*Conj(ZD(gI2,gO2))*Cos(
      ThetaWp()),0) + IF(gO2 < 3,0.5*g2*Conj(ZD(gI2,gO2))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gO2 < 3,0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Sin(
      ThetaW())*Sin(ThetaWp()),0) - gp*Qd*Cos(ThetaWp())*SUM(j1,0,2,Conj(ZD(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1)) - 0.2581988897471611*g1*Sin(ThetaW())*Sin
      (ThetaWp())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

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
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,g2*gp*Ql*Cos(ThetaW())*
      KroneckerDelta(gO1,gO2)*Sin(2*ThetaWp()),0) + IF(gO1 < 3,0.7745966692414834*
      g1*gp*Ql*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(2*ThetaWp()),0) + IF(gO1
      < 3,0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,0.3*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())),0) +
      IF(gO1 < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())),0)
      + 2*Sqr(gp)*Sqr(Qv)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZpVZp(int gO1, int gO2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-2*g2*gp*Ql*Cos(ThetaW())*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1 < 3,-
      1.5491933384829668*g1*gp*Ql*Cos(ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(
      ThetaW())*Sin(ThetaWp()),0) + IF(gO1 < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*
      Sqr(Ql)*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,0.7745966692414834*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(
      gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp())),0) + IF(gO1 < 3,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())),0) + 2*Sqr(gp)*Sqr(Qv)*Sqr(Cos(ThetaWp()))*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSvconjHpmconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(
      j1,gO2))*Ye(j1,gO1))*ZP(gI1,0)*ZP(gI2,0)),0),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-(QHd*Ql
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(QHu*Ql
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)),0) - QHd*Qv*Sqr(gp)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*
      ZP(gI2,0) - QHu*Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) - SUM(j3,0,2,KroneckerDelta(
      gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yv(j1,j3))
      *Yv(j1,j2))))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFeconjUSvPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Ye(j1,gO2))*ZER(
      gI2,j1))*UM(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFeconjUSvPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UP(gI1,0))*Conj(ZEL(
      gI2,gO1))),0) + Conj(UP(gI1,1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*SUM(j1
      ,0,2,Conj(ZEL(gI2,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjHpmconjUSv(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZE(
      gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*
      Conj(TYe(j1,gO2)))*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,0,
      2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,j2)))*ZP(gI1,0),0) + IF
      (gO2 < 3,-0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) + IF
      (gO2 < 3,0.7071067811865475*vS*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(
      gI2,3 + j1)))*ZP(gI1,1),0) + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,
      Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Yv(gO2,j1))*Yv(j2,j1)))*ZP(gI1,1),0) +
      0.7071067811865475*vS*Conj(Lambdax)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*
      SUM(j1,0,2,Conj(ZE(gI2,j1))*Yv(j1,j2)))*ZP(gI1,0) + 0.7071067811865475*vu*
      SUM(j3,0,2,Conj(ZE(gI2,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1
      ,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))))*ZP(gI1,0) + SUM(j2,0,2,KroneckerDelta(gO2,
      3 + j2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*TYv(j1,j2)))*ZP(gI1,1) +
      0.7071067811865475*vd*SUM(j3,0,2,Conj(ZE(gI2,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))))*ZP(gI1,1)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZA(gI1,1))*
      Conj(ZA(gI2,1))*SUM(j1,0,2,Conj(Yv(gO2,j1))*Yv(gO1,j1))),0),0) + IF(gO1 < 3,
      -0.15*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) +
      IF(gO1 < 3,0.15*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,-0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1
      ,gO2)*Sqr(g2),0) + IF(gO1 < 3,0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Ql*Conj(ZA(gI1,0))*
      Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Ql*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Ql*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*SUM(j2,
      0,2,KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2)),0) + IF(gO1 < 3,-0.5*Conj(Lambdax
      )*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*Yv(
      gO1,j2)),0) + IF(gO2 < 3,-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*Lambdax*SUM(j2
      ,0,2,Conj(Yv(gO2,j2))*KroneckerDelta(gO1,3 + j2)),0) + IF(gO2 < 3,-0.5*Conj(
      ZA(gI1,0))*Conj(ZA(gI2,2))*Lambdax*SUM(j2,0,2,Conj(Yv(gO2,j2))*
      KroneckerDelta(gO1,3 + j2)),0) - QHd*Qv*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(
      gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - QHu*
      Qv*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1)) - Qs*Qv*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(
      gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Conj
      (ZA(gI1,1))*Conj(ZA(gI2,1))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2
      ,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZH(gI1,1))*
      Conj(ZH(gI2,1))*SUM(j1,0,2,Conj(Yv(gO2,j1))*Yv(gO1,j1))),0),0) + IF(gO1 < 3,
      -0.15*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) +
      IF(gO1 < 3,0.15*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,-0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1
      ,gO2)*Sqr(g2),0) + IF(gO1 < 3,0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Ql*Conj(ZH(gI1,0))*
      Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Ql*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Ql*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))*SUM(j2,0
      ,2,KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2)),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*
      Conj(ZH(gI1,0))*Conj(ZH(gI2,2))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*Yv(gO1
      ,j2)),0) + IF(gO2 < 3,0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))*Lambdax*SUM(j2,0,2
      ,Conj(Yv(gO2,j2))*KroneckerDelta(gO1,3 + j2)),0) + IF(gO2 < 3,0.5*Conj(ZH(
      gI1,0))*Conj(ZH(gI2,2))*Lambdax*SUM(j2,0,2,Conj(Yv(gO2,j2))*KroneckerDelta(
      gO1,3 + j2)),0) - QHd*Qv*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - QHu*Qv*Conj(ZH(gI1,
      1))*Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) - Qs*Qv*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Conj(ZH(
      gI1,1))*Conj(ZH(gI2,1))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiFvconjUSvPR(int gI2, int gI1, int gO2) const
{
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = IF(gO2 < 3,-(SUM(j2,0,2,Conj(Yv(gO2,j2))*
      ZVR(gI1,j2))*ZN(gI2,4)),0) - 1.4142135623730951*gp*Qv*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*ZVR(gI1,j1))*ZN(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFvconjUSvPL(int gI2, int gI1, int gO1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Ql*Conj(
      ZN(gI2,0))*Conj(ZVL(gI1,gO1)),0) + IF(gO1 < 3,0.5477225575051661*g1*Conj(ZN(
      gI2,1))*Conj(ZVL(gI1,gO1)),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZN(
      gI2,2))*Conj(ZVL(gI1,gO1)),0) - Conj(ZN(gI2,4))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*SUM(j1,0,2,Conj(ZVL(gI1,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSvconjSdconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1
      ,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)),0) + IF(gO1 < 3,-0.5*Qd*Ql*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) +
      IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI1,3 +
      j2))*ZD(gI2,3 + j2)),0) + IF(gO1 < 3,-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(
      gp)*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2)),0) - 0.5*Qq*Qv*Sqr(gp)*
      SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2
      )*KroneckerDelta(gO2,3 + j2)) - 0.5*Qd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,3 +
      j1))*ZD(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) - 0.5*Qq*Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) - 0.5*Qd
      *Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSvconjSeconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Ye(j1,
      gO1)*ZE(gI2,3 + j1))*SUM(j3,0,2,Conj(Ye(j3,gO2))*Conj(ZE(gI1,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZE(gI1,gO2))*Sqr(g2)*ZE(gI2,gO1),0),0) +
      IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1)
      )*ZE(gI2,j1)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0
      ,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)
      *Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,
      0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,
      3 + j1)),0) + IF(gO1 < 3,-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0
      ,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))
      *ZE(gI2,j2)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.15*KroneckerDelta(
      gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1
       < 3,-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) - 0.5*Ql*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(
      gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) -
      0.5*Qe*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2
      ,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Ql*Qv*Sqr(gp)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2
      ,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 0.5*Qe*Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(
      gI2,3 + j2)) - SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*Yv(j1,j2)))*SUM(j4,0,2,KroneckerDelta(gO1,3 + j4)*SUM(j3,0,2,Conj(Yv(j3
      ,j4))*ZE(gI2,j3)));

   return result;
}

std::complex<double> CLASSNAME::CpSuUSvconjSuconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qu = LOCALINPUT(Qu);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1
      ,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 3,-0.5*Ql*Qu*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) +
      IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 +
      j2))*ZU(gI2,3 + j2)),0) + IF(gO1 < 3,-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(
      gp)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j2
      ,0,2,KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3
      ,j4))*Conj(ZU(gI1,3 + j3)))*ZU(gI2,j4))),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(
      ZU(gI1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI2,3 + j1)))*SUM(j4,0,2,Conj(Yv(gO2,j4)
      )*KroneckerDelta(gO1,3 + j4))),0) - 0.5*Qq*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1
      ,j1))*ZU(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
       + j2)) - 0.5*Qu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qq*
      Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))
      *SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)) - 0.5*Qu*Qv*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI1,3 + j2))*ZU(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSvUSvconjSvconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j2,0,2,Yv(gO1,
      j2)*ZV(gI2,3 + j2))*SUM(j4,0,2,Conj(Yv(gO2,j4))*Conj(ZV(gI1,3 + j4)))),0),0)
      + IF(gO1 < 3,IF(gO2 < 3,-0.15*Conj(ZV(gI1,gO2))*Sqr(g1)*ZV(gI2,gO1),0),0) +
      IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZV(gI1,gO2))*Sqr(g2)*ZV(gI2,gO1),0),0) + IF
      (gO1 < 3,IF(gO2 < 3,-(Conj(ZV(gI1,gO2))*Sqr(gp)*Sqr(Ql)*ZV(gI2,gO1)),0),0) +
      IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gI1,j1)
      )*ZV(gI2,j1)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,
      0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2
      )*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-
      0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV
      (gI2,3 + j1)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,
      0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,
      gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) + IF(gO1 < 3,-0.5*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,
      j2)),0) + IF(gO1 < 3,-0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,
      Conj(ZV(gI1,3 + j2))*ZV(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2))*SUM(j4,0,2,Conj(ZV(gI1,3 + j4))*SUM(
      j3,0,2,Conj(Yv(j3,j4))*ZV(gI2,j3)))),0) + IF(gO1 < 3,-0.5*Ql*Qv*Sqr(gp)*SUM(
      j1,0,2,Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZV(gI2,gO1),0) + IF(
      gO1 < 3,-0.5*Ql*Qv*Sqr(gp)*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*KroneckerDelta(
      gO2,3 + j2))*ZV(gI2,gO1),0) + IF(gO2 < 3,-0.5*Ql*Qv*Conj(ZV(gI1,gO2))*Sqr(gp
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZV(gI2,3 + j1)),0) + IF(gO2 < 3,-0.5
      *Ql*Qv*Conj(ZV(gI1,gO2))*Sqr(gp)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZV(
      gI2,3 + j2)),0) + IF(gO2 < 3,-(SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,
      j2))*ZV(gI2,3 + j2))*SUM(j4,0,2,Conj(Yv(gO2,j4))*KroneckerDelta(gO1,3 + j4))
      ),0) - 0.5*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZV(gI2,3 +
      j1))*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.5*Ql*Qv
      *Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,
      Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.5*Ql*Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,
      j2)) - 0.5*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*ZV(gI2,3 + j2))
      - 0.5*Sqr(gp)*Sqr(Qv)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZV(gI2,3 + j2)) - SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,j2)))*SUM(j4,0,
      2,KroneckerDelta(gO1,3 + j4)*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gI2,j3)));

   return result;
}

std::complex<double> CLASSNAME::CpAhSvconjUSv(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,0.5)*vS*
      Conj(ZA(gI2,0))*Lambdax*SUM(j2,0,2,Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2))),0)
      + IF(gO2 < 3,std::complex<double>(0,0.5)*vd*Conj(ZA(gI2,2))*Lambdax*SUM(j2,0
      ,2,Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2))),0) + IF(gO2 < 3,std::complex<
      double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j2,0,2,Conj(ZV(gI1,3 + j2
      ))*Conj(TYv(gO2,j2))),0) - std::complex<double>(0,0.5)*vS*Conj(Lambdax)*Conj
      (ZA(gI2,0))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(ZV(gI1,j1)
      )*Yv(j1,j2))) - std::complex<double>(0,0.5)*vd*Conj(Lambdax)*Conj(ZA(gI2,2))
      *SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,j2)
      )) - std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*TYv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CphhSvconjUSv(int gI2, int gI1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO2 < 3,-0.15*vd*Conj(ZH(gI2,0))*Conj(ZV
      (gI1,gO2))*Sqr(g1),0) + IF(gO2 < 3,0.15*vu*Conj(ZH(gI2,1))*Conj(ZV(gI1,gO2))
      *Sqr(g1),0) + IF(gO2 < 3,-0.25*vd*Conj(ZH(gI2,0))*Conj(ZV(gI1,gO2))*Sqr(g2),
      0) + IF(gO2 < 3,0.25*vu*Conj(ZH(gI2,1))*Conj(ZV(gI1,gO2))*Sqr(g2),0) + IF(
      gO2 < 3,-(QHd*Ql*vd*Conj(ZH(gI2,0))*Conj(ZV(gI1,gO2))*Sqr(gp)),0) + IF(gO2 <
      3,-(QHu*Ql*vu*Conj(ZH(gI2,1))*Conj(ZV(gI1,gO2))*Sqr(gp)),0) + IF(gO2 < 3,-(
      Ql*Qs*vS*Conj(ZH(gI2,2))*Conj(ZV(gI1,gO2))*Sqr(gp)),0) + IF(gO2 < 3,0.5*vS*
      Conj(ZH(gI2,0))*Lambdax*SUM(j2,0,2,Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2))),0)
      + IF(gO2 < 3,0.5*vd*Conj(ZH(gI2,2))*Lambdax*SUM(j2,0,2,Conj(Yv(gO2,j2))*Conj
      (ZV(gI1,3 + j2))),0) + IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI2,1))*SUM(j2
      ,0,2,Conj(ZV(gI1,3 + j2))*Conj(TYv(gO2,j2))),0) + IF(gO2 < 3,-(vu*Conj(ZH(
      gI2,1))*SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Yv(gO2,j1))*Yv(j2,j1))))
      ,0) - QHd*Qv*vd*Conj(ZH(gI2,0))*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) - QHu*Qv*vu*Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,
      Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - Qs*Qv*vS*Conj(ZH(gI2,2))*
      Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) + 0.5*vS
      *Conj(Lambdax)*Conj(ZH(gI2,0))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,
      0,2,Conj(ZV(gI1,j1))*Yv(j1,j2))) + 0.5*vd*Conj(Lambdax)*Conj(ZH(gI2,2))*SUM(
      j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,j2))) -
      0.7071067811865475*Conj(ZH(gI2,1))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM
      (j1,0,2,Conj(ZV(gI1,j1))*TYv(j1,j2))) - vu*Conj(ZH(gI2,1))*SUM(j3,0,2,Conj(
      ZV(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yv(j1,
      j3))*Yv(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSvconjVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZE(
      gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSvVZ(int gI2, int gO2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO2 < 3,0.5*g2*Conj(ZV(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,0.3872983346207417*g1*Conj(ZV(gI2,
      gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,gp*Ql*Conj(ZV(gI2,gO2))*
      Sin(ThetaWp()),0) - gp*Qv*Sin(ThetaWp())*SUM(j1,0,2,Conj(ZV(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSvVZp(int gI2, int gO2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO2 < 3,gp*Ql*Conj(ZV(gI2,gO2))*Cos(
      ThetaWp()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gO2 < 3,-0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(
      ThetaW())*Sin(ThetaWp()),0) - gp*Qv*Cos(ThetaWp())*SUM(j1,0,2,Conj(ZV(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(int gO1, int gO2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-0.5163977794943222*g1*gp*Qq*Cos
      (ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp()),0) + IF(gO1
       < 3,g2*gp*Qq*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*Sin(2*ThetaWp()),0) + IF
      (gO1 < 3,-0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*Sin
      (ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,
      0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(
      Sin(ThetaW())),0) + IF(gO1 < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr
      (Sin(ThetaWp())),0) + 2.065591117977289*g1*gp*Qu*Cos(ThetaWp())*Sin(ThetaW()
      )*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
      + j1)) + 0.5333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 2*Sqr(gp
      )*Sqr(Qu)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZpVZp(int gO1, int gO2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-2*g2*gp*Qq*Cos(ThetaW())*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1 < 3,
      0.2581988897471611*g1*gp*Qq*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(2*
      ThetaWp()),0) + IF(gO1 < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos
      (ThetaWp())),0) + IF(gO1 < 3,-0.2581988897471611*g1*g2*Cos(ThetaW())*
      KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0
      ) + IF(gO1 < 3,0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp())),0) - 2.065591117977289*g1*gp*Qu*Cos(ThetaWp()
      )*Sin(ThetaW())*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 2*Sqr(gp)*Sqr(Qu)*Sqr(Cos(ThetaWp()))*SUM(j1,0
      ,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) +
      0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))*SUM(j1,0,2
      ,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSuconjUSuconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSuconjHpmconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(
      j1,gO2))*Yd(j1,gO1))*ZP(gI1,0)*ZP(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-(QHd*Qq
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(QHu*Qq
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)),0) - 0.2*Sqr(g1)*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(
      gI2,0) - QHd*Qu*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) - QHu*Qu*Sqr(gp)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)
      *ZP(gI2,1) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta
      (gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFdconjUSuPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Yd(j1,gO2))*ZDR(
      gI2,j1))*UM(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaFdconjUSuPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UP(gI1,0))*Conj(ZDL(
      gI2,gO1))),0) + Conj(UP(gI1,1))*SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjHpmconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZD(
      gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      Conj(TYd(j1,gO2)))*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,0,
      2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,j2)))*ZP(gI1,0),0) + IF
      (gO2 < 3,-0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) + IF
      (gO2 < 3,0.7071067811865475*vS*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(
      gI2,3 + j1)))*ZP(gI1,1),0) + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,
      Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZP(gI1,1),0) +
      0.7071067811865475*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZP(gI1,0) + 0.7071067811865475*vu*SUM
      (j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,
      2,Conj(Yd(j3,j1))*Yu(j2,j1))))*ZP(gI1,0) + SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(
      j1,0,2,KroneckerDelta(gO2,3 + j1)*TYu(j1,j2)))*ZP(gI1,1) +
      0.7071067811865475*vd*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))))*ZP(gI1,1)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZA(gI1,1))*
      Conj(ZA(gI2,1))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,gO1))),0),0) + IF(gO1 < 3,
      0.05*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF
      (gO1 < 3,-0.05*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,-0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1
      ,gO2)*Sqr(g2),0) + IF(gO1 < 3,0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Qq*Conj(ZA(gI1,0))*
      Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Qq*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Qq*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*SUM(j1,
      0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1)),0) + IF(gO1 < 3,-0.5*Conj(Lambdax
      )*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(
      j1,gO1)),0) + IF(gO2 < 3,-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*Lambdax*SUM(j1
      ,0,2,Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1)),0) + IF(gO2 < 3,-0.5*Conj(
      ZA(gI1,0))*Conj(ZA(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*
      KroneckerDelta(gO1,3 + j1)),0) - 0.2*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.2*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) - QHd*Qu*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(
      gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - QHu*
      Qu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1)) - Qs*Qu*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(
      gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Conj
      (ZA(gI1,1))*Conj(ZA(gI2,1))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2
      ,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZH(gI1,1))*
      Conj(ZH(gI2,1))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,gO1))),0),0) + IF(gO1 < 3,
      0.05*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF
      (gO1 < 3,-0.05*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,-0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1
      ,gO2)*Sqr(g2),0) + IF(gO1 < 3,0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Qq*Conj(ZH(gI1,0))*
      Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Qq*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Qq*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))*SUM(j1,0
      ,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1)),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*
      Conj(ZH(gI1,0))*Conj(ZH(gI2,2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,
      gO1)),0) + IF(gO2 < 3,0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))*Lambdax*SUM(j1,0,2
      ,Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1)),0) + IF(gO2 < 3,0.5*Conj(ZH(
      gI1,0))*Conj(ZH(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*KroneckerDelta(
      gO1,3 + j1)),0) - 0.2*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 0.2*Conj(ZH(gI1,1))
      *Conj(ZH(gI2,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) - QHd*Qu*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - QHu*Qu*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) - Qs*Qu*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Conj(
      ZH(gI1,1))*Conj(ZH(gI2,1))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPR(int gI2, int gI1, int gO2) const
{
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO2))*
      ZUR(gI1,j1))*ZN(gI2,4)),0) + 0.09428090415820634*SUM(j1,0,2,KroneckerDelta(
      gO2,3 + j1)*ZUR(gI1,j1))*(-15*gp*Qu*ZN(gI2,0) + 7.745966692414834*g1*ZN(gI2,
      1));

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPL(int gI2, int gI1, int gO1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Qq*Conj(
      ZN(gI2,0))*Conj(ZUL(gI1,gO1)),0) + IF(gO1 < 3,-0.18257418583505536*g1*Conj(
      ZN(gI2,1))*Conj(ZUL(gI1,gO1)),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZN
      (gI2,2))*Conj(ZUL(gI1,gO1)),0) - Conj(ZN(gI2,4))*SUM(j2,0,2,Conj(ZUL(gI1,j2)
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSuconjSeconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1
      ,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-0.5*Qe*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(
      gp)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.5*Ql*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)
      ) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qe*Qu*Sqr(gp)*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*ZE(gI2,j2)) - 0.5*Ql*Qu*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) +
      0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.5*Qe*Qu*Sqr(gp)*SUM(j1
      ,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(
      ZE(gI1,3 + j2))*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yd(j1,
      gO1)*ZD(gI1,3 + j1))*SUM(j3,0,2,Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1),0),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1)
      )*ZD(gI1,j1)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0
      ,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-1.5*KroneckerDelta(gO1,gO2)
      *Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-
      0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,
      3 + j1)),0) + IF(gO1 < 3,-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0
      ,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) +
      IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI2,j2))
      *ZD(gI1,j2)),0) + IF(gO1 < 3,-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*
      SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta
      (gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + IF(
      gO1 < 3,-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZD(gI2,3
      + j2))*ZD(gI1,3 + j2)),0) + 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,
      j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 1.5
      *Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 1.5*Qd*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI2,3
       + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta
      (gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) - 1.5*Qq
      *Qu*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*ZD(gI1,3 + j2)) - 1.5*Qd*Qu*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(
      gI1,3 + j2)) - SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 +
      j1)*Yu(j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 +
      j3))*ZD(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yu(j1,
      gO1)*ZU(gI1,3 + j1))*SUM(j3,0,2,Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666*Conj(ZU(gI2,gO2))*Sqr(g1)*ZU
      (gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(
      gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-1.3333333333333333*Conj(ZU(gI2,gO2))
      *Sqr(g3)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZU(gI2,gO2))*Sqr(
      gp)*Sqr(Qq)*ZU(gI1,gO1)),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) +
      IF(gO1 < 3,-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*SUM(j1,0,2,Conj(ZU(
      gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM
      (j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-1.5*Qq*Qu*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 +
      j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,-0.375*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,-1.5*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,
      j2)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,-1.5*Qq*Qu*KroneckerDelta(gO1,
      gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3
      ,-3*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,
      Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3)))*ZU(gI1,j4)),0) + IF(gO1 < 3,
      0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j1,0
      ,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1
      < 3,-0.5*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3
      + j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,0.03333333333333333*Sqr(g1)*SUM(j2,0,2,
      Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO1 < 3
      ,0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*KroneckerDelta(
      gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO1 < 3,-0.5*Qq*Qu*Sqr(gp)*SUM(j2,0,2,Conj(
      ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO2 < 3,
      0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*ZU(gI1,3 + j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZU(gI2,gO2))
      *Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1)),0) + IF(gO2 <
      3,-0.5*Qq*Qu*Conj(ZU(gI2,gO2))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)
      *ZU(gI1,3 + j1)),0) + IF(gO2 < 3,0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(
      g1)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 < 3,
      0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1,3
       + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 < 3,-0.5*Qq*Qu*Conj(ZU(gI2,gO2))*Sqr(gp)*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 < 3,-3*SUM
      (j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1)))*SUM(j3,0,2,
      Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3)),0) - 0.13333333333333333*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1
      ,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 +
      j2))*KroneckerDelta(gO2,3 + j2)) - 0.5*Sqr(gp)*Sqr(Qu)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*
      KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1
      ,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) -
      1.5*Qq*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.4*Sqr(g1)*SUM(j1,
      0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 1.5*Sqr(gp)*Sqr(Qu)*SUM(j1,0,2,Conj(ZU(gI2
      ,3 + j1))*ZU(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) -
      1.5*Qq*Qu*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
       + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) - 0.4*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*ZU(gI1,3 + j2)) - 1.5*Sqr(gp)*Sqr(Qu)*SUM(j1,0,2,KroneckerDelta
      (gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(
      gI1,3 + j2)) - 0.13333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 +
      j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 +
      j2)) - 0.5*Sqr(gp)*Sqr(Qu)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2)) - SUM(j2,
      0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*SUM(
      j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZU(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSuSvconjUSuconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1
      ,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-0.5*Qq*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)
      *SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) +
      IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZV(gI1,j2)
      )*ZV(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM
      (j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) + IF(gO1 < 3,-0.5*Qq*Qv*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*ZV(gI2,3 +
      j2)),0) + IF(gO1 < 3,-(SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1))*SUM
      (j4,0,2,Conj(ZV(gI1,3 + j4))*SUM(j3,0,2,Conj(Yv(j3,j4))*ZV(gI2,j3)))),0) +
      IF(gO2 < 3,-(SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,j2))*ZV(gI2,3 + j2
      ))*SUM(j3,0,2,Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3))),0) - 0.1*Sqr(g1)
      *SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Ql*Qu*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,
      j1))*ZV(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
      + j2)) - 0.5*Qu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2
      ,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)) - 0.5*Ql*Qu*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(
      gI1,j2))*ZV(gI2,j2)) - 0.5*Qu*Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*ZV(gI2,3 +
      j2));

   return result;
}

std::complex<double> CLASSNAME::CpAhSuconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,0.5)*vS*
      Conj(ZA(gI2,0))*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1))),0)
      + IF(gO2 < 3,std::complex<double>(0,0.5)*vd*Conj(ZA(gI2,2))*Lambdax*SUM(j1,0
      ,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1))),0) + IF(gO2 < 3,std::complex<
      double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1
      ))*Conj(TYu(j1,gO2))),0) - std::complex<double>(0,0.5)*vS*Conj(Lambdax)*Conj
      (ZA(gI2,0))*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1
      )*Yu(j1,j2))) - std::complex<double>(0,0.5)*vd*Conj(Lambdax)*Conj(ZA(gI2,2))
      *SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)
      )) - std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CphhSuconjUSu(int gI2, int gI1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO2 < 3,0.05*vd*Conj(ZH(gI2,0))*Conj(ZU(
      gI1,gO2))*Sqr(g1),0) + IF(gO2 < 3,-0.05*vu*Conj(ZH(gI2,1))*Conj(ZU(gI1,gO2))
      *Sqr(g1),0) + IF(gO2 < 3,-0.25*vd*Conj(ZH(gI2,0))*Conj(ZU(gI1,gO2))*Sqr(g2),
      0) + IF(gO2 < 3,0.25*vu*Conj(ZH(gI2,1))*Conj(ZU(gI1,gO2))*Sqr(g2),0) + IF(
      gO2 < 3,-(QHd*Qq*vd*Conj(ZH(gI2,0))*Conj(ZU(gI1,gO2))*Sqr(gp)),0) + IF(gO2 <
      3,-(QHu*Qq*vu*Conj(ZH(gI2,1))*Conj(ZU(gI1,gO2))*Sqr(gp)),0) + IF(gO2 < 3,-(
      Qq*Qs*vS*Conj(ZH(gI2,2))*Conj(ZU(gI1,gO2))*Sqr(gp)),0) + IF(gO2 < 3,0.5*vS*
      Conj(ZH(gI2,0))*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1))),0)
      + IF(gO2 < 3,0.5*vd*Conj(ZH(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj
      (ZU(gI1,3 + j1))),0) + IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI2,1))*SUM(j1
      ,0,2,Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2))),0) + IF(gO2 < 3,-(vu*Conj(ZH(
      gI2,1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,j2))))
      ,0) - 0.2*vd*Conj(ZH(gI2,0))*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) + 0.2*vu*Conj(ZH(gI2,1))*Sqr(g1)*SUM(j1,0,2,Conj
      (ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - QHd*Qu*vd*Conj(ZH(gI2,0))*Sqr
      (gp)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - QHu*Qu*vu
      *Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,
      3 + j1)) - Qs*Qu*vS*Conj(ZH(gI2,2))*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) + 0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,0))*SUM(j2,0,
      2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2))) + 0.5*
      vd*Conj(Lambdax)*Conj(ZH(gI2,2))*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j2))) - 0.7071067811865475*Conj(ZH(gI2,1))*
      SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYu(j1,j2)
      )) - vu*Conj(ZH(gI2,1))*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))));

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
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO2 < 3,0.5*g2*Conj(ZU(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZU(gI2
      ,gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,gp*Qq*Conj(ZU(gI2,gO2))*
      Sin(ThetaWp()),0) - 0.5163977794943222*g1*Cos(ThetaWp())*Sin(ThetaW())*SUM(
      j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) - gp*Qu*Sin(ThetaWp(
      ))*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVZp(int gI2, int gO2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO2 < 3,gp*Qq*Conj(ZU(gI2,gO2))*Cos(
      ThetaWp()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZU(gI2,gO2))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gO2 < 3,0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Sin(
      ThetaW())*Sin(ThetaWp()),0) - gp*Qu*Cos(ThetaWp())*SUM(j1,0,2,Conj(ZU(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1)) + 0.5163977794943222*g1*Sin(ThetaW())*Sin
      (ThetaWp())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(int gO1, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-2*g2*gp*Ql*Cos(ThetaW())*Cos(
      ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaWp()),0) + IF(gO1 < 3,
      0.7745966692414834*g1*gp*Ql*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(2*
      ThetaWp()),0) + IF(gO1 < 3,-0.7745966692414834*g1*g2*Cos(ThetaW())*
      KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp())),0) + IF(gO1 < 3,
      0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())),0
      ) + IF(gO1 < 3,0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(
      Sin(ThetaW())),0) + IF(gO1 < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*Sqr
      (Sin(ThetaWp())),0) - 3.0983866769659336*g1*gp*Qe*Cos(ThetaWp())*Sin(ThetaW(
      ))*Sin(ThetaWp())*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
       + j1)) + 1.2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 2*Sqr(gp)*Sqr(Qe)*
      Sqr(Sin(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2
      ,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZpVZp(int gO1, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-1.5491933384829668*g1*gp*Ql*Cos
      (ThetaWp())*KroneckerDelta(gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp()),0) + IF(gO1
       < 3,g2*gp*Ql*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*Sin(2*ThetaWp()),0) + IF
      (gO1 < 3,2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*Sqr(Cos(ThetaWp())),0) +
      IF(gO1 < 3,-0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,gO2)*
      Sin(ThetaW())*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2
      )*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp())),0) + IF(gO1 < 3,0.3*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())),0) +
      3.0983866769659336*g1*gp*Qe*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) + 2*Sqr(gp)*
      Sqr(Qe)*Sqr(Cos(ThetaWp()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) + 1.2*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp
      ()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSeconjUSeconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSeconjHpmconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yv(
      gO2,j1))*Yv(gO1,j1))*ZP(gI1,1)*ZP(gI2,1)),0),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-(QHd*Ql
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(QHu*Ql
      *KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)),0) + 0.3*Sqr(g1)*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(
      gI2,0) - Qe*QHd*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta
      (gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*
      SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1)))
      )*ZP(gI1,0)*ZP(gI2,0) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) - Qe*QHu*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZA(gI1,0))*
      Conj(ZA(gI2,0))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,gO1))),0),0) + IF(gO1 < 3,
      -0.15*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) +
      IF(gO1 < 3,0.15*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,
      gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Ql*Conj(ZA(gI1,0))*
      Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Ql*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Ql*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*SUM(j1,
      0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1)),0) + IF(gO1 < 3,-0.5*Conj(Lambdax
      )*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(
      j1,gO1)),0) + IF(gO2 < 3,-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*Lambdax*SUM(j1
      ,0,2,Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1)),0) + IF(gO2 < 3,-0.5*Conj(
      ZA(gI1,1))*Conj(ZA(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*
      KroneckerDelta(gO1,3 + j1)),0) + 0.3*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.3*
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) - Qe*QHd*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(
      gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Qe*
      QHu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1)) - Qe*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*
      Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) -
      Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2
      ,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZH(gI1,0))*
      Conj(ZH(gI2,0))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,gO1))),0),0) + IF(gO1 < 3,
      -0.15*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1),0) +
      IF(gO1 < 3,0.15*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(
      g1),0) + IF(gO1 < 3,0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,
      gO2)*Sqr(g2),0) + IF(gO1 < 3,-0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*
      KroneckerDelta(gO1,gO2)*Sqr(g2),0) + IF(gO1 < 3,-(QHd*Ql*Conj(ZH(gI1,0))*
      Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1 < 3,-(QHu*Ql*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(gp)),0) + IF(gO1
       < 3,-(Ql*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))*SUM(j1,0
      ,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1)),0) + IF(gO1 < 3,0.5*Conj(Lambdax)*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,
      gO1)),0) + IF(gO2 < 3,0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))*Lambdax*SUM(j1,0,2
      ,Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1)),0) + IF(gO2 < 3,0.5*Conj(ZH(
      gI1,1))*Conj(ZH(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*KroneckerDelta(
      gO1,3 + j1)),0) + 0.3*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - 0.3*Conj(ZH(gI1,1))
      *Conj(ZH(gI2,1))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1)) - Qe*QHd*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp)
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Qe*QHu*
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1
      )*KroneckerDelta(gO2,3 + j1)) - Qe*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - Conj(
      ZH(gI1,0))*Conj(ZH(gI2,0))*SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))));

   return result;
}

std::complex<double> CLASSNAME::CpChaFvconjUSePR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(Yv(gO2,j2))*ZVR(
      gI1,j2))*UP(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpChaFvconjUSePL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UM(gI2,0))*Conj(ZVL(
      gI1,gO1))),0) + Conj(UM(gI2,1))*SUM(j2,0,2,Conj(ZVL(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePR(int gI2, int gI1, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(j1,gO2))*
      ZER(gI1,j1))*ZN(gI2,3)),0) - 0.28284271247461906*SUM(j1,0,2,KroneckerDelta(
      gO2,3 + j1)*ZER(gI1,j1))*(5*gp*Qe*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePL(int gI2, int gI1, int gO1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Ql*Conj(
      ZEL(gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,0.5477225575051661*g1*Conj(ZEL
      (gI1,gO1))*Conj(ZN(gI2,1)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZEL(
      gI1,gO1))*Conj(ZN(gI2,2)),0) - Conj(ZN(gI2,3))*SUM(j2,0,2,Conj(ZEL(gI1,j2))*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSeconjSdconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1
      ,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)),0) + IF(gO1 < 3,-0.5*Qd*Ql*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(
      g2)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) +
      IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI1,3 +
      j2))*ZD(gI2,3 + j2)),0) + IF(gO1 < 3,-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(
      gp)*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1
      ,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3
      ,j4))*Conj(ZD(gI1,3 + j3)))*ZD(gI2,j4))),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(
      ZD(gI1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI2,3 + j1)))*SUM(j3,0,2,Conj(Ye(j3,gO2)
      )*KroneckerDelta(gO1,3 + j3))),0) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))
      *ZD(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2
      )) - 0.5*Qe*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qd*Qe*Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI1,3
       + j1))*ZD(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta
      (gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) - 0.5*Qe
      *Qq*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI1,3 + j2))*ZD(gI2,3 + j2)) - 0.5*Qd*Qe*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(
      gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSeconjSeconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Ye(j1,
      gO1)*ZE(gI2,3 + j1))*SUM(j3,0,2,Conj(Ye(j3,gO2))*Conj(ZE(gI1,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.15*Conj(ZE(gI1,gO2))*Sqr(g1)*ZE(gI2,gO1),0),0)
      + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZE(gI1,gO2))*Sqr(g2)*ZE(gI2,gO1),0),0) +
      IF(gO1 < 3,IF(gO2 < 3,-(Conj(ZE(gI1,gO2))*Sqr(gp)*Sqr(Ql)*ZE(gI2,gO1)),0),0)
      + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(
      j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,
      gO2)*Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3
      ,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2
      ,3 + j1)),0) + IF(gO1 < 3,-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,
      0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-0.075*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2)
      )*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.15*KroneckerDelta(
      gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1
       < 3,-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)
      *Ye(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gI1,3 + j3)))*ZE(
      gI2,j4))),0) + IF(gO1 < 3,0.15*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZE(gI2,gO1),0) + IF(gO1 < 3,-0.5*Qe*Ql*Sqr(gp)*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZE(gI2,gO1),0) +
      IF(gO1 < 3,0.15*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*KroneckerDelta(gO2,3
       + j2))*ZE(gI2,gO1),0) + IF(gO1 < 3,-0.5*Qe*Ql*Sqr(gp)*SUM(j2,0,2,Conj(ZE(
      gI1,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZE(gI2,gO1),0) + IF(gO2 < 3,0.15*
      Conj(ZE(gI1,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 +
      j1)),0) + IF(gO2 < 3,-0.5*Qe*Ql*Conj(ZE(gI1,gO2))*Sqr(gp)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1)),0) + IF(gO2 < 3,0.15*Conj(ZE(gI1,
      gO2))*Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)),0) + IF(
      gO2 < 3,-0.5*Qe*Ql*Conj(ZE(gI1,gO2))*Sqr(gp)*SUM(j2,0,2,KroneckerDelta(gO1,3
       + j2)*ZE(gI2,3 + j2)),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,
      0,2,Ye(j1,j2)*ZE(gI2,3 + j1)))*SUM(j3,0,2,Conj(Ye(j3,gO2))*KroneckerDelta(
      gO1,3 + j3))),0) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,
      3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.5*
      Sqr(gp)*Sqr(Qe)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1))*SUM(j2
      ,0,2,Conj(ZE(gI1,3 + j2))*KroneckerDelta(gO2,3 + j2)) + 0.15*Sqr(g1)*SUM(j1,
      0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.5*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)
      ) - 0.3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Sqr(gp)*Sqr(Qe)
      *SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.15*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*ZE(gI2,j2)) - 0.5*Qe*Ql*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)) -
      0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.5*Sqr(gp)*Sqr(Qe)*SUM(
      j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,
      Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3
      + j1))*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(
      gI2,3 + j2)) - 0.5*Sqr(gp)*Sqr(Qe)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 +
      j2)) - SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(
      j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3))*
      ZE(gI2,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSuconjUSeconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 3,-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1
      ,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(
      j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 3,-0.5*Ql*Qu*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 +
      j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(
      ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,-0.5*Ql*Qq*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) +
      IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 +
      j2))*ZU(gI2,3 + j2)),0) + IF(gO1 < 3,-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(
      gp)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 + j2)),0) - 0.05*Sqr(g1)*SUM(j1
      ,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) - 0.5*Qe*Qq*Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI1,j1))*
      ZU(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)
      ) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qe*Qu*Sqr(gp)*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(
      gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI1,j2))*ZU(gI2,j2)) - 0.5*Qe*Qq*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)) +
      0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)
      )*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 + j2)) - 0.5*Qe*Qu*Sqr(gp)*SUM(j1
      ,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(
      ZU(gI1,3 + j2))*ZU(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSvconjUSeconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j2,0,2,Yv(gO1,
      j2)*ZV(gI2,3 + j2))*SUM(j4,0,2,Conj(Yv(gO2,j4))*Conj(ZV(gI1,3 + j4)))),0),0)
      + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZV(gI1,gO2))*Sqr(g2)*ZV(gI2,gO1),0),0) +
      IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZV(gI1,j1)
      )*ZV(gI2,j1)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0
      ,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)
      *Sqr(gp)*Sqr(Ql)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)),0) + IF(gO1 < 3,-
      0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV
      (gI2,3 + j1)),0) + IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,
      0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta(gO1,
      gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)),0) + IF(gO1 < 3,-0.5*
      KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,
      j2)),0) + IF(gO1 < 3,-0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)*SUM(j2,0,2,
      Conj(ZV(gI1,3 + j2))*ZV(gI2,3 + j2)),0) + 0.15*Sqr(g1)*SUM(j1,0,2,Conj(ZV(
      gI1,j1))*ZV(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) - 0.5*Qe*Ql*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1))*SUM
      (j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.5*Qe*Qv*
      Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.15*Sqr(g1)*SUM(j1
      ,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(
      ZV(gI1,j2))*ZV(gI2,j2)) - 0.5*Qe*Ql*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZV(gI1,j2))*ZV(gI2,j2)) -
      0.5*Qe*Qv*Sqr(gp)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
       + j1))*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*ZV(gI2,3 + j2)) - SUM(j2,0,2,Conj(ZV
      (gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*SUM(j4,0,2,SUM(
      j3,0,2,Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZV(gI2,j4));

   return result;
}

std::complex<double> CLASSNAME::CpHpmSvconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZV(
      gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*vS*Lambdax*
      SUM(j2,0,2,Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2)))*ZP(gI2,0),0) + IF(gO2 < 3,
      0.7071067811865475*vd*SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2)
      )*Ye(j1,j2)))*ZP(gI2,0),0) + IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZV(gI1,
      gO2))*Sqr(g2)*ZP(gI2,1),0) + IF(gO2 < 3,SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*Conj
      (TYv(gO2,j2)))*ZP(gI2,1),0) + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,
      Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Yv(gO2,j1))*Yv(j2,j1)))*ZP(gI2,1),0) + SUM(
      j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYe(j1,j2)))*
      ZP(gI2,0) + 0.7071067811865475*vu*SUM(j3,0,2,Conj(ZV(gI1,3 + j3))*SUM(j2,0,2
      ,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))))*ZP(gI2,0
      ) + 0.7071067811865475*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0
      ,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZP(gI2,1) + 0.7071067811865475*vd*
      SUM(j3,0,2,Conj(ZV(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1
      ,0,2,Conj(Yv(j1,j3))*Ye(j2,j1))))*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhSeconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,0.5)*vS*
      Conj(ZA(gI2,1))*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1))),0)
      + IF(gO2 < 3,std::complex<double>(0,0.5)*vu*Conj(ZA(gI2,2))*Lambdax*SUM(j1,0
      ,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1))),0) + IF(gO2 < 3,std::complex<
      double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1
      ))*Conj(TYe(j1,gO2))),0) - std::complex<double>(0,0.5)*vS*Conj(Lambdax)*Conj
      (ZA(gI2,1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1
      )*Ye(j1,j2))) - std::complex<double>(0,0.5)*vu*Conj(Lambdax)*Conj(ZA(gI2,2))
      *SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)
      )) - std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYe(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CphhSeconjUSe(int gI2, int gI1, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO2 < 3,-0.15*vd*Conj(ZE(gI1,gO2))*Conj(
      ZH(gI2,0))*Sqr(g1),0) + IF(gO2 < 3,0.15*vu*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,1))
      *Sqr(g1),0) + IF(gO2 < 3,0.25*vd*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(g2),0
      ) + IF(gO2 < 3,-0.25*vu*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(g2),0) + IF(
      gO2 < 3,-(QHd*Ql*vd*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(gp)),0) + IF(gO2 <
      3,-(QHu*Ql*vu*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(gp)),0) + IF(gO2 < 3,-(
      Ql*Qs*vS*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,2))*Sqr(gp)),0) + IF(gO2 < 3,0.5*vS*
      Conj(ZH(gI2,1))*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1))),0)
      + IF(gO2 < 3,0.5*vu*Conj(ZH(gI2,2))*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj
      (ZE(gI1,3 + j1))),0) + IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI2,0))*SUM(j1
      ,0,2,Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2))),0) + IF(gO2 < 3,-(vd*Conj(ZH(
      gI2,0))*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,j2))))
      ,0) + 0.3*vd*Conj(ZH(gI2,0))*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) - 0.3*vu*Conj(ZH(gI2,1))*Sqr(g1)*SUM(j1,0,2,Conj
      (ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - Qe*QHd*vd*Conj(ZH(gI2,0))*Sqr
      (gp)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1)) - Qe*QHu*vu
      *Conj(ZH(gI2,1))*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,
      3 + j1)) - Qe*Qs*vS*Conj(ZH(gI2,2))*Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1)) + 0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,1))*SUM(j2,0,
      2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2))) + 0.5*
      vu*Conj(Lambdax)*Conj(ZH(gI2,2))*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j2))) - 0.7071067811865475*Conj(ZH(gI2,0))*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYe(j1,j2)
      )) - vd*Conj(ZH(gI2,0))*SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))));

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
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Conj(ZE(gI2,gO2))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gO2 < 3,0.3872983346207417*g1*Conj(ZE(gI2,
      gO2))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gO2 < 3,gp*Ql*Conj(ZE(gI2,gO2))*
      Sin(ThetaWp()),0) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())*SUM(
      j1,0,2,Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1)) - gp*Qe*Sin(ThetaWp(
      ))*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSeVZp(int gI2, int gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO2 < 3,gp*Ql*Conj(ZE(gI2,gO2))*Cos(
      ThetaWp()),0) + IF(gO2 < 3,0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gO2 < 3,-0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(
      ThetaW())*Sin(ThetaWp()),0) - gp*Qe*Cos(ThetaWp())*SUM(j1,0,2,Conj(ZE(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1)) - 0.7745966692414834*g1*Sin(ThetaW())*Sin
      (ThetaWp())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSeVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZV(
      gI2,gO2)),0);

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*vS*KroneckerDelta(2,gO1)*Sqr(gp)*
      Sqr(Qs)*Sqr(Sin(ThetaWp())) - vd*KroneckerDelta(0,gO1)*(20*g2*gp*QHd*Cos(
      ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr
      (Sin(ThetaWp()))) - vu*KroneckerDelta(1,gO1)*(-20*g2*gp*QHu*Cos(ThetaW())*
      Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgZUhh(int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*vS*Cos(ThetaWp())*KroneckerDelta(
      2,gO1)*Sin(ThetaWp())*Sqr(gp)*Sqr(Qs) + vd*KroneckerDelta(0,gO1)*(5*Cos(
      ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) - 7.745966692414834*g1*
      gp*QHd*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-
      20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) + 7.745966692414834*g1*
      gp*QHd*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 5*gp*QHd*
      Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(ThetaWp())))) + vu*KroneckerDelta(1,
      gO1)*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) +
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp(
      ))*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) -
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(
      ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())
      + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgZpUhh(int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*vS*KroneckerDelta(2,gO1)*Sqr(gp)*
      Sqr(Qs)*Sqr(Cos(ThetaWp())) - vd*KroneckerDelta(0,gO1)*(-2*gp*QHd*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp
      )*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos
      (ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(
      ThetaWp()))) - vu*KroneckerDelta(1,gO1)*(2*gp*QHu*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*vS*KroneckerDelta(2,gO2)*Sqr(gp)*
      Sqr(Qs)*Sqr(Sin(ThetaWp())) + vd*KroneckerDelta(0,gO2)*(20*g2*gp*QHd*Cos(
      ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr
      (Sin(ThetaWp()))) + vu*KroneckerDelta(1,gO2)*(-20*g2*gp*QHu*Cos(ThetaW())*
      Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZp(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*vS*Cos(ThetaWp())*KroneckerDelta(2,
      gO2)*Sin(ThetaWp())*Sqr(gp)*Sqr(Qs) - vd*KroneckerDelta(0,gO2)*(5*Cos(
      ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) - 7.745966692414834*g1*
      gp*QHd*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-
      20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) + 7.745966692414834*g1*
      gp*QHd*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 5*gp*QHd*
      Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(ThetaWp())))) - vu*KroneckerDelta(1,
      gO2)*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) +
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp(
      ))*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) -
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(
      ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())
      + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZpVZp(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*vS*KroneckerDelta(2,gO2)*Sqr(gp)*
      Sqr(Qs)*Sqr(Cos(ThetaWp())) + vd*KroneckerDelta(0,gO2)*(-2*gp*QHd*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp
      )*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos
      (ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(
      ThetaWp()))) + vu*KroneckerDelta(1,gO2)*(2*gp*QHu*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp())) + KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(20*g2*gp*QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(
      ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(-20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*Sin(
      ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZpVZp(int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp())) + KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(-2*gp*QHd*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(2*gp*QHu*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*((AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gI1,0)*ZP(gI2,0
      ) + (AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gI1,1)*ZP(gI2,1)) - KroneckerDelta
      (0,gO1)*(5*KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*
      ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr
      (g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) +
      20*QHd*QHu*Sqr(gp))*ZP(gI1,1)*ZP(gI2,1))) + KroneckerDelta(1,gO1)*(-5*
      KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) +
      ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*(Sqr(g2) + 4*
      QHd*QHu*Sqr(gp)))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*
      Sqr(QHu)))*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO2)*(ZP(gI1,0)*(
      vd*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI2,0) + 5*vu*(-2*AbsSqr
      (Lambdax) + Sqr(g2))*ZP(gI2,1)) + ZP(gI1,1)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(
      g2))*ZP(gI2,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI2,1))
      )) - 10*KroneckerDelta(2,gO2)*(ZP(gI1,0)*(2*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr
      (gp))*ZP(gI2,0) + 1.4142135623730951*Conj(TLambdax)*ZP(gI2,1)) + ZP(gI1,1)*(
      1.4142135623730951*TLambdax*ZP(gI2,0) + 2*vS*(AbsSqr(Lambdax) + QHu*Qs*Sqr(
      gp))*ZP(gI2,1))) + KroneckerDelta(1,gO2)*(ZP(gI1,0)*(vu*(3*Sqr(g1) - 5*(Sqr(
      g2) + 4*QHd*QHu*Sqr(gp)))*ZP(gI2,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP
      (gI2,1)) - ZP(gI1,1)*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vu*(3*
      Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZP(gI2,1))));

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

std::complex<double> CLASSNAME::CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1)
      - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*KroneckerDelta(2,gO1)*KroneckerDelta(
      2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))) +
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)
      *(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp
      )) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20
      *Sqr(gp)*Sqr(QHu))) - 20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*(KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp
      )) + KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-(Conj(ZH(gI1,0))*(20*Conj(ZH(gI2,2))
      *(KroneckerDelta(0,gO2)*KroneckerDelta(2,gO1) + KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) - Conj(ZH(gI2,1))*
      (KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*
      KroneckerDelta(1,gO2))*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd
      *QHu*Sqr(gp))) + Conj(ZH(gI2,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2
      )*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp
      )) + 3*KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) +
      4*Sqr(gp)*Sqr(QHd)))))) + Conj(ZH(gI1,1))*(-20*Conj(ZH(gI2,2))*(
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1) + KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gI2,0))*
      (KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*
      KroneckerDelta(1,gO2))*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd
      *QHu*Sqr(gp))) + Conj(ZH(gI2,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2
      )*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp
      )) - 3*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) +
      20*Sqr(gp)*Sqr(QHu)))) - 20*Conj(ZH(gI1,2))*(Conj(ZH(gI2,0))*(KroneckerDelta
      (0,gO2)*KroneckerDelta(2,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2))
      *(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZH(gI2,1))*(KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1) + KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2))*(
      AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gI2,2))*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhh(int gI1, int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-5*Conj(ZA(gI1,2))*(
      1.4142135623730951*Conj(TLambdax)*(Conj(ZA(gI2,1))*KroneckerDelta(0,gO2) +
      Conj(ZA(gI2,0))*KroneckerDelta(1,gO2)) + 4*Conj(ZA(gI2,2))*(vd*
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + vu*KroneckerDelta
      (1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + vS*KroneckerDelta(2,gO2)*Sqr(gp
      )*Sqr(Qs)) + 1.4142135623730951*(Conj(ZA(gI2,1))*KroneckerDelta(0,gO2) +
      Conj(ZA(gI2,0))*KroneckerDelta(1,gO2))*TLambdax) + Conj(ZA(gI1,1))*(Conj(ZA(
      gI2,1))*(vd*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(
      g2) - 20*QHd*QHu*Sqr(gp)) - 20*vS*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) +
      QHu*Qs*Sqr(gp)) - vu*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(
      gp)*Sqr(QHu))) - 7.0710678118654755*(Conj(ZA(gI2,2))*KroneckerDelta(0,gO2) +
      Conj(ZA(gI2,0))*KroneckerDelta(2,gO2))*(Conj(TLambdax) + TLambdax)) - Conj(
      ZA(gI1,0))*(Conj(ZA(gI2,0))*(vu*KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) -
      3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*vS*KroneckerDelta(2,gO2)*(
      AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + vd*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*
      (Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))) + 7.0710678118654755*(Conj(ZA(gI2,2))*
      KroneckerDelta(1,gO2) + Conj(ZA(gI2,1))*KroneckerDelta(2,gO2))*(Conj(
      TLambdax) + TLambdax)));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhUhh(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(Conj(ZA(gI2,2))*(Conj(ZH(gI1,1))*KroneckerDelta(0,gO2)
      + Conj(ZH(gI1,0))*KroneckerDelta(1,gO2)) + Conj(ZA(gI2,1))*(Conj(ZH(gI1,2))*
      KroneckerDelta(0,gO2) + Conj(ZH(gI1,0))*KroneckerDelta(2,gO2)) + Conj(ZA(gI2
      ,0))*(Conj(ZH(gI1,2))*KroneckerDelta(1,gO2) + Conj(ZH(gI1,1))*KroneckerDelta
      (2,gO2)))*(Conj(TLambdax) - TLambdax);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhh(int gI1, int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-5*Conj(ZH(gI1,2))*(-
      1.4142135623730951*Conj(TLambdax)*Conj(ZH(gI2,1))*KroneckerDelta(0,gO2) + 4*
      vd*AbsSqr(Lambdax)*Conj(ZH(gI2,2))*KroneckerDelta(0,gO2) + 4*vS*AbsSqr(
      Lambdax)*Conj(ZH(gI2,1))*KroneckerDelta(1,gO2) + 4*vu*AbsSqr(Lambdax)*Conj(
      ZH(gI2,2))*KroneckerDelta(1,gO2) + 4*vu*AbsSqr(Lambdax)*Conj(ZH(gI2,1))*
      KroneckerDelta(2,gO2) + 4*QHd*Qs*vd*Conj(ZH(gI2,2))*KroneckerDelta(0,gO2)*
      Sqr(gp) + 4*QHu*Qs*vS*Conj(ZH(gI2,1))*KroneckerDelta(1,gO2)*Sqr(gp) + 4*QHu*
      Qs*vu*Conj(ZH(gI2,2))*KroneckerDelta(1,gO2)*Sqr(gp) + 4*QHu*Qs*vu*Conj(ZH(
      gI2,1))*KroneckerDelta(2,gO2)*Sqr(gp) + 12*vS*Conj(ZH(gI2,2))*KroneckerDelta
      (2,gO2)*Sqr(gp)*Sqr(Qs) - 1.4142135623730951*Conj(ZH(gI2,1))*KroneckerDelta(
      0,gO2)*TLambdax + Conj(ZH(gI2,0))*(-1.4142135623730951*Conj(TLambdax)*
      KroneckerDelta(1,gO2) + 4*vd*AbsSqr(Lambdax)*KroneckerDelta(2,gO2) + 4*QHd*
      Qs*vd*KroneckerDelta(2,gO2)*Sqr(gp) + 4*vS*KroneckerDelta(0,gO2)*(AbsSqr(
      Lambdax) + QHd*Qs*Sqr(gp)) - 1.4142135623730951*KroneckerDelta(1,gO2)*
      TLambdax)) + Conj(ZH(gI1,1))*(Conj(ZH(gI2,1))*(vd*KroneckerDelta(0,gO2)*(-20
      *AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*vS*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - 3*vu*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))) + 5*
      Conj(ZH(gI2,2))*(1.4142135623730951*Conj(TLambdax)*KroneckerDelta(0,gO2) - 4
      *vu*AbsSqr(Lambdax)*KroneckerDelta(2,gO2) - 4*QHu*Qs*vu*KroneckerDelta(2,gO2
      )*Sqr(gp) - 4*vS*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) +
      1.4142135623730951*KroneckerDelta(0,gO2)*TLambdax) + Conj(ZH(gI2,0))*(vu*
      KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*
      QHu*Sqr(gp)) + vd*KroneckerDelta(1,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5
      *Sqr(g2) - 20*QHd*QHu*Sqr(gp)) + 7.0710678118654755*KroneckerDelta(2,gO2)*(
      Conj(TLambdax) + TLambdax))) - Conj(ZH(gI1,0))*(Conj(ZH(gI2,0))*(vu*
      KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*
      QHu*Sqr(gp)) + 20*vS*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)
      ) + 3*vd*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))
      )) + 5*Conj(ZH(gI2,2))*(-1.4142135623730951*Conj(TLambdax)*KroneckerDelta(1,
      gO2) + 4*vd*AbsSqr(Lambdax)*KroneckerDelta(2,gO2) + 4*QHd*Qs*vd*
      KroneckerDelta(2,gO2)*Sqr(gp) + 4*vS*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax)
      + QHd*Qs*Sqr(gp)) - 1.4142135623730951*KroneckerDelta(1,gO2)*TLambdax) -
      Conj(ZH(gI2,1))*(vu*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) +
      5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) + vd*KroneckerDelta(1,gO2)*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) + 7.0710678118654755*
      KroneckerDelta(2,gO2)*(Conj(TLambdax) + TLambdax))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZDL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO1)*
      SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZEL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO1)*
      SUM(j2,0,2,Conj(ZEL(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO1)*
      SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFvUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO2)*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZVL(gI1,j1))*ZVR(gI2,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFvUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO1)*
      SUM(j2,0,2,Conj(ZVR(gI1,j2))*SUM(j1,0,2,Conj(ZVL(gI2,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO1)*(-2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) +
      Qd*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) + KroneckerDelta(1,gO2)*
      (Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI2,3 +
      j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*
      ZD(gI2,j2)))) - KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((Sqr(g1) + 5*
      Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(
      Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))
      ) - 10*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(
      j1,0,2,Yd(j1,j2)*ZD(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1
      ,j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2)))) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((Sqr(g1) + 5*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2
      ,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(j1,0,2,
      Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI2,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI2
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO1)*(-2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) +
      Qe*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) + KroneckerDelta(1,gO2)*
      (Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 +
      j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*
      ZE(gI2,j2)))) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5
      *(Sqr(g2) + 4*QHu*Ql*Sqr(gp)))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 2*(
      3*Sqr(g1) + 10*Qe*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1
      ))) + 10*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*
      SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(
      Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*ZE(gI2,j2)))) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,
      2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + 20*Qe*QHd*Sqr(gp))*SUM(j1,0,2
      ,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI2
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO1)*(-2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) +
      Qu*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) + KroneckerDelta(0,gO2)*
      (Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI2,3 +
      j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*
      ZU(gI2,j2)))) + KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((Sqr(g1) - 5*(
      Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*(
      Sqr(g1) + 5*QHd*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)))
      + 10*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1
      ,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,
      j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,j2)))) - KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,
      Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0,2,
      Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI2,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI2
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO1)*(-2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)) +
      Qv*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))) + KroneckerDelta(0,gO2)*
      (Lambdax*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gI2,
      j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,j2))*ZV(
      gI2,3 + j2)))) - KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((3*Sqr(g1) +
      5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)) + 20*
      QHd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))) - 10*
      KroneckerDelta(2,gO2)*(Lambdax*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*SUM(j1,0,2,
      Conj(Yv(j1,j2))*ZV(gI2,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(
      gI1,j1))*Yv(j1,j2))*ZV(gI2,3 + j2)))) + KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)*((3*Sqr(g1) + 5*Sqr(g2) - 20*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gI1,
      j1))*ZV(gI2,j1)) - 20*(QHu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2
      ,3 + j1)) + SUM(j3,0,2,Conj(ZV(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1
      ,j3))*Yv(j1,j2))*ZV(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gI1,j2))*
      SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSdconjSd(int gO2, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO2)*(-2*Qq*Qs*vS
      *Sqr(gp)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) - 2*Qd*Qs*vS*Sqr(gp)*SUM(j1
      ,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)) + vu*Conj(Lambdax)*SUM(j2,0,2,Conj
      (ZD(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + vu*Lambdax*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2))) -
      KroneckerDelta(1,gO2)*(vu*(Sqr(g1) + 5*(Sqr(g2) + 4*QHu*Qq*Sqr(gp)))*SUM(j1,
      0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) + 2*(vu*(Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)) - 5*vS*(Conj(Lambdax)*SUM(j2,0,2
      ,Conj(ZD(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + Lambdax*SUM(j2,0,2
      ,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2))))) +
      KroneckerDelta(0,gO2)*(vd*(Sqr(g1) + 5*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))*SUM(j1,
      0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) + 2*(vd*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(
      j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)) - 5*(1.4142135623730951*SUM(j2,0
      ,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,ZD(gI1,3 + j1)*TYd(j1,j2))) +
      1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2
      )))*ZD(gI1,j2)) + 2*vd*(SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI1,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(
      ZD(gI2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI1,j3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSeconjSe(int gO2, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO2)*(-2*Ql*Qs*vS
      *Sqr(gp)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) - 2*Qe*Qs*vS*Sqr(gp)*SUM(j1
      ,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) + vu*Conj(Lambdax)*SUM(j2,0,2,Conj
      (ZE(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + vu*Lambdax*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2))) +
      KroneckerDelta(1,gO2)*(vu*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHu*Ql*Sqr(gp)))*SUM(
      j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) - 2*vu*(3*Sqr(g1) + 10*Qe*QHu*Sqr(gp))*
      SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) + 10*vS*(Conj(Lambdax)*SUM(
      j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + Lambdax*SUM(
      j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2)))) -
      KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1
      ,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) + (-6*vd*Sqr(g1) + 20*Qe*QHd*vd*Sqr(gp))*
      SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) + 10*(1.4142135623730951*SUM
      (j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,ZE(gI1,3 + j1)*TYe(j1,j2))) +
      1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2
      )))*ZE(gI1,j2)) + 2*vd*(SUM(j3,0,2,Conj(ZE(gI2,3 + j3))*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI1,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(
      ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI1,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSuconjSu(int gO2, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO2)*(-2*Qq*Qs*vS
      *Sqr(gp)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 2*Qs*Qu*vS*Sqr(gp)*SUM(j1
      ,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + vd*Conj(Lambdax)*SUM(j2,0,2,Conj
      (ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + vd*Lambdax*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2))) +
      KroneckerDelta(0,gO2)*(vd*(Sqr(g1) - 5*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,
      0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*vd*(Sqr(g1) + 5*QHd*Qu*Sqr(gp))*SUM(j1,
      0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + 10*vS*(Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + Lambdax*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2)))) -
      KroneckerDelta(1,gO2)*(vu*(Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0
      ,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*vu*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0
      ,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + 10*(1.4142135623730951*SUM(j2,0,2,
      Conj(ZU(gI2,j2))*SUM(j1,0,2,ZU(gI1,3 + j1)*TYu(j1,j2))) + 1.4142135623730951
      *SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2)))*ZU(gI1,j2)) +
      2*vu*(SUM(j3,0,2,Conj(ZU(gI2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*
      Yu(j2,j1))*ZU(gI1,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,
      0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI1,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSvconjSv(int gO2, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.05*(10*KroneckerDelta(2,gO2)*(-2*Ql*Qs*vS
      *Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI2,j1))*ZV(gI1,j1)) - 2*Qs*Qv*vS*Sqr(gp)*SUM(j1
      ,0,2,Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1)) + vd*Lambdax*SUM(j2,0,2,Conj(ZV(
      gI2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gI1,j1))) + vd*Conj(Lambdax)*SUM(
      j2,0,2,SUM(j1,0,2,Conj(ZV(gI2,j1))*Yv(j1,j2))*ZV(gI1,3 + j2))) -
      KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) + 5*(Sqr(g2) + 4*QHd*Ql*Sqr(gp)))*SUM(
      j1,0,2,Conj(ZV(gI2,j1))*ZV(gI1,j1)) + 20*QHd*Qv*vd*Sqr(gp)*SUM(j1,0,2,Conj(
      ZV(gI2,3 + j1))*ZV(gI1,3 + j1)) - 10*vS*(Lambdax*SUM(j2,0,2,Conj(ZV(gI2,3 +
      j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gI1,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(
      j1,0,2,Conj(ZV(gI2,j1))*Yv(j1,j2))*ZV(gI1,3 + j2)))) + KroneckerDelta(1,gO2)
      *(vu*(3*Sqr(g1) + 5*Sqr(g2) - 20*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gI2,j1))
      *ZV(gI1,j1)) - 10*(2*QHu*Qv*vu*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI2,3 + j1))*ZV(
      gI1,3 + j1)) + 1.4142135623730951*SUM(j2,0,2,Conj(ZV(gI2,3 + j2))*SUM(j1,0,2
      ,Conj(TYv(j1,j2))*ZV(gI1,j1))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,
      Conj(ZV(gI2,j1))*TYv(j1,j2))*ZV(gI1,3 + j2)) + 2*vu*SUM(j3,0,2,Conj(ZV(gI2,3
       + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j3))*Yv(j1,j2))*ZV(gI1,3 + j2))) +
      2*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gI2,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2
      ,j1)))*ZV(gI1,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUhhPR(int gI1, int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(5*KroneckerDelta(2,gO2)*(-2*gp*Qs*ZN(
      gI1,5)*ZN(gI2,0) + 1.4142135623730951*Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2,3) +
      ZN(gI1,3)*ZN(gI2,4)) - 2*gp*Qs*ZN(gI1,0)*ZN(gI2,5)) + KroneckerDelta(0,gO2)*
      (ZN(gI1,3)*(-10*gp*QHd*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(
      gI2,2)) - 10*gp*QHd*ZN(gI1,0)*ZN(gI2,3) + 3.872983346207417*g1*ZN(gI1,1)*ZN(
      gI2,3) - 5*g2*ZN(gI1,2)*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,
      5)*ZN(gI2,4) + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,5)) -
      KroneckerDelta(1,gO2)*(ZN(gI1,4)*(10*gp*QHu*ZN(gI2,0) + 3.872983346207417*g1
      *ZN(gI2,1) - 5*g2*ZN(gI2,2)) + (10*gp*QHu*ZN(gI1,0) + 3.872983346207417*g1*
      ZN(gI1,1) - 5*g2*ZN(gI1,2))*ZN(gI2,4) - 7.0710678118654755*Conj(Lambdax)*(ZN
      (gI1,5)*ZN(gI2,3) + ZN(gI1,3)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUhhPL(int gI1, int gI2, int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(3.872983346207417*g1*Conj(ZN(gI1,1))*
      Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) - 5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,3))
      *KroneckerDelta(0,gO1) - 10*gp*QHu*Conj(ZN(gI1,4))*Conj(ZN(gI2,0))*
      KroneckerDelta(1,gO1) - 3.872983346207417*g1*Conj(ZN(gI1,4))*Conj(ZN(gI2,1))
      *KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*KroneckerDelta
      (1,gO1) - 3.872983346207417*g1*Conj(ZN(gI1,1))*Conj(ZN(gI2,4))*
      KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,4))*KroneckerDelta(
      1,gO1) - 10*gp*Qs*Conj(ZN(gI1,5))*Conj(ZN(gI2,0))*KroneckerDelta(2,gO1) - 10
      *gp*Conj(ZN(gI1,0))*(QHd*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) + QHu*Conj(ZN
      (gI2,4))*KroneckerDelta(1,gO1) + Qs*Conj(ZN(gI2,5))*KroneckerDelta(2,gO1)) +
      7.0710678118654755*Conj(ZN(gI1,5))*Conj(ZN(gI2,4))*KroneckerDelta(0,gO1)*
      Lambdax + 7.0710678118654755*Conj(ZN(gI1,4))*Conj(ZN(gI2,5))*KroneckerDelta(
      0,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,5))*Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,4))*Conj(ZN(
      gI2,3))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,3))*(-10*gp*QHd*Conj(ZN(
      gI2,0))*KroneckerDelta(0,gO1) + 3.872983346207417*g1*Conj(ZN(gI2,1))*
      KroneckerDelta(0,gO1) - 5*g2*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) +
      7.0710678118654755*Conj(ZN(gI2,5))*KroneckerDelta(1,gO1)*Lambdax +
      7.0710678118654755*Conj(ZN(gI2,4))*KroneckerDelta(2,gO1)*Lambdax));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,0.1)*(10*gp*Qs*Conj(
      ZA(gI2,2))*KroneckerDelta(2,gO2)*Sin(ThetaWp()) + Conj(ZA(gI2,0))*
      KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417
      *g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) - Conj(ZA(gI2,1
      ))*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      ));

   return result;
}

std::complex<double> CLASSNAME::CpAhUhhVZp(int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,0.1)*(10*gp*Qs*Conj(
      ZA(gI2,2))*Cos(ThetaWp())*KroneckerDelta(2,gO2) + Conj(ZA(gI2,1))*
      KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(
      ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())) + Conj(ZA(
      gI2,0))*KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()
      ) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp())) + KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(20*g2*gp*QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(
      ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(-20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*Sin(
      ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZpVZp(int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp())) + KroneckerDelta(0
      ,gO1)*KroneckerDelta(0,gO2)*(-2*gp*QHd*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(2*gp*QHu*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*((AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gI1,0)*ZP(gI2,0
      ) + (AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gI1,1)*ZP(gI2,1)) - KroneckerDelta
      (0,gO1)*(-5*KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*
      ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr
      (g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) +
      20*QHd*QHu*Sqr(gp))*ZP(gI1,1)*ZP(gI2,1))) + KroneckerDelta(1,gO1)*(5*
      KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) +
      ZP(gI1,0)*ZP(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*(Sqr(g2) + 4*
      QHd*QHu*Sqr(gp)))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*
      Sqr(QHu)))*ZP(gI1,1)*ZP(gI2,1))));

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

std::complex<double> CLASSNAME::CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-(Conj(ZA(gI1,0))*(20*Conj(ZA(gI2,2))
      *(KroneckerDelta(0,gO2)*KroneckerDelta(2,gO1) + KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) - Conj(ZA(gI2,1))*
      (KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*
      KroneckerDelta(1,gO2))*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd
      *QHu*Sqr(gp))) + Conj(ZA(gI2,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2
      )*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp
      )) + 3*KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) +
      4*Sqr(gp)*Sqr(QHd)))))) + Conj(ZA(gI1,1))*(-20*Conj(ZA(gI2,2))*(
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1) + KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZA(gI2,0))*
      (KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*
      KroneckerDelta(1,gO2))*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd
      *QHu*Sqr(gp))) + Conj(ZA(gI2,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2
      )*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp
      )) - 3*KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) +
      20*Sqr(gp)*Sqr(QHu)))) - 20*Conj(ZA(gI1,2))*(Conj(ZA(gI2,0))*(KroneckerDelta
      (0,gO2)*KroneckerDelta(2,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2))
      *(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + Conj(ZA(gI2,1))*(KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1) + KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2))*(
      AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZA(gI2,2))*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-(Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*(
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1)
      - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*KroneckerDelta(2,gO1)*KroneckerDelta(
      2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))) +
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)
      *(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp
      )) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20
      *Sqr(gp)*Sqr(QHu))) - 20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*(KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp
      )) + KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.35355339059327373
      )*(Conj(ZA(gI1,2))*(Conj(ZA(gI2,1))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*
      KroneckerDelta(1,gO2)) + Conj(ZA(gI1,1))*(Conj(ZA(gI2,2))*KroneckerDelta(0,
      gO2) + Conj(ZA(gI2,0))*KroneckerDelta(2,gO2)) + Conj(ZA(gI1,0))*(Conj(ZA(gI2
      ,2))*KroneckerDelta(1,gO2) + Conj(ZA(gI2,1))*KroneckerDelta(2,gO2)))*(Conj(
      TLambdax) - TLambdax);

   return result;
}

std::complex<double> CLASSNAME::CpAhUAhhh(int gI2, int gO2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*vd*AbsSqr(Lambdax)*Conj(ZA(gI2,1)
      )*Conj(ZH(gI1,0))*KroneckerDelta(1,gO2) - 20*vS*AbsSqr(Lambdax)*Conj(ZA(gI2,
      1))*Conj(ZH(gI1,2))*KroneckerDelta(1,gO2) - 20*vd*AbsSqr(Lambdax)*Conj(ZA(
      gI2,2))*Conj(ZH(gI1,0))*KroneckerDelta(2,gO2) - 20*vu*AbsSqr(Lambdax)*Conj(
      ZA(gI2,2))*Conj(ZH(gI1,1))*KroneckerDelta(2,gO2) - 7.0710678118654755*Conj(
      TLambdax)*(Conj(ZA(gI2,2))*(Conj(ZH(gI1,1))*KroneckerDelta(0,gO2) + Conj(ZH(
      gI1,0))*KroneckerDelta(1,gO2)) + Conj(ZA(gI2,1))*(Conj(ZH(gI1,2))*
      KroneckerDelta(0,gO2) + Conj(ZH(gI1,0))*KroneckerDelta(2,gO2))) + 3*vd*Conj(
      ZA(gI2,1))*Conj(ZH(gI1,0))*KroneckerDelta(1,gO2)*Sqr(g1) - 3*vu*Conj(ZA(gI2,
      1))*Conj(ZH(gI1,1))*KroneckerDelta(1,gO2)*Sqr(g1) + 5*vd*Conj(ZA(gI2,1))*
      Conj(ZH(gI1,0))*KroneckerDelta(1,gO2)*Sqr(g2) - 5*vu*Conj(ZA(gI2,1))*Conj(ZH
      (gI1,1))*KroneckerDelta(1,gO2)*Sqr(g2) - 20*QHd*QHu*vd*Conj(ZA(gI2,1))*Conj(
      ZH(gI1,0))*KroneckerDelta(1,gO2)*Sqr(gp) - 20*QHu*Qs*vS*Conj(ZA(gI2,1))*Conj
      (ZH(gI1,2))*KroneckerDelta(1,gO2)*Sqr(gp) - 20*QHd*Qs*vd*Conj(ZA(gI2,2))*
      Conj(ZH(gI1,0))*KroneckerDelta(2,gO2)*Sqr(gp) - 20*QHu*Qs*vu*Conj(ZA(gI2,2))
      *Conj(ZH(gI1,1))*KroneckerDelta(2,gO2)*Sqr(gp) - 20*vu*Conj(ZA(gI2,1))*Conj(
      ZH(gI1,1))*KroneckerDelta(1,gO2)*Sqr(gp)*Sqr(QHu) - 20*vS*Conj(ZA(gI2,2))*
      Conj(ZH(gI1,2))*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs) - 7.0710678118654755*
      Conj(ZA(gI2,2))*Conj(ZH(gI1,1))*KroneckerDelta(0,gO2)*TLambdax -
      7.0710678118654755*Conj(ZA(gI2,1))*Conj(ZH(gI1,2))*KroneckerDelta(0,gO2)*
      TLambdax - 7.0710678118654755*Conj(ZA(gI2,2))*Conj(ZH(gI1,0))*KroneckerDelta
      (1,gO2)*TLambdax - 7.0710678118654755*Conj(ZA(gI2,1))*Conj(ZH(gI1,0))*
      KroneckerDelta(2,gO2)*TLambdax - Conj(ZA(gI2,0))*(vd*Conj(ZH(gI1,0))*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))) + 5*
      Conj(ZH(gI1,2))*(4*vS*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp
      )) + 1.4142135623730951*KroneckerDelta(1,gO2)*(Conj(TLambdax) + TLambdax)) +
      Conj(ZH(gI1,1))*(vu*KroneckerDelta(0,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) -
      5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 7.0710678118654755*KroneckerDelta(2,gO2)*(
      Conj(TLambdax) + TLambdax))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(Conj(ZH(gI1,2))*(Conj(ZH(gI2,1))*KroneckerDelta(0,gO2)
      + Conj(ZH(gI2,0))*KroneckerDelta(1,gO2)) + Conj(ZH(gI1,1))*(Conj(ZH(gI2,2))*
      KroneckerDelta(0,gO2) + Conj(ZH(gI2,0))*KroneckerDelta(2,gO2)) + Conj(ZH(gI1
      ,0))*(Conj(ZH(gI2,2))*KroneckerDelta(1,gO2) + Conj(ZH(gI2,1))*KroneckerDelta
      (2,gO2)))*(Conj(TLambdax) - TLambdax);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*
      ZDL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(0,gO1)*SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1
      ,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI2,j1))*
      ZEL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(0,gO1)*SUM(j2,0,2,Conj(ZEL(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1
      ,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(1,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*
      ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1
      ,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFvUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *KroneckerDelta(1,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZVL(gI1,j1))*
      ZVR(gI2,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFvUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(ZVR(gI1,j2))*SUM(j1,0,2,Conj(ZVL(gI2
      ,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.05*(-10*KroneckerDelta(2,gO1)*(2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) +
      Qd*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) + KroneckerDelta(1,gO2)*
      (Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI2,3 +
      j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*
      ZD(gI2,j2)))) - KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((Sqr(g1) + 5*
      Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(
      Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))
      ) + 10*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(
      j1,0,2,Yd(j1,j2)*ZD(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1
      ,j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2)))) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((Sqr(g1) + 5*(Sqr(g2) - 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2
      ,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(j1,0,2,
      Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI2,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI2
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.05*(-10*KroneckerDelta(2,gO1)*(2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) +
      Qe*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) + KroneckerDelta(1,gO2)*
      (Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 +
      j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*
      ZE(gI2,j2)))) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5
      *(Sqr(g2) + 4*QHu*Ql*Sqr(gp)))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 2*(
      3*Sqr(g1) + 10*Qe*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1
      ))) - 10*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*
      SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(
      Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*ZE(gI2,j2)))) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,
      2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + 20*Qe*QHd*Sqr(gp))*SUM(j1,0,2
      ,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI2
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.05*(-10*KroneckerDelta(2,gO1)*(2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Qq*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) +
      Qu*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) + KroneckerDelta(0,gO2)*
      (Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI2,3 +
      j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*
      ZU(gI2,j2)))) + KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((Sqr(g1) - 5*(
      Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*(
      Sqr(g1) + 5*QHd*Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)))
      - 10*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1
      ,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,
      j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,j2)))) - KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,
      Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0,2,
      Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI2,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI2
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.05*(-10*KroneckerDelta(2,gO1)*(2*Qs*
      KroneckerDelta(2,gO2)*Sqr(gp)*(Ql*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)) +
      Qv*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))) + KroneckerDelta(0,gO2)*
      (Lambdax*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gI2,
      j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gI1,j1))*Yv(j1,j2))*ZV(
      gI2,3 + j2)))) - KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((3*Sqr(g1) +
      5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)) + 20*
      QHd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1))) + 10*
      KroneckerDelta(2,gO2)*(Lambdax*SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*SUM(j1,0,2,
      Conj(Yv(j1,j2))*ZV(gI2,j1))) + Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(
      gI1,j1))*Yv(j1,j2))*ZV(gI2,3 + j2)))) + KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)*((3*Sqr(g1) + 5*Sqr(g2) - 20*QHu*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gI1,
      j1))*ZV(gI2,j1)) - 20*(QHu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2
      ,3 + j1)) + SUM(j3,0,2,Conj(ZV(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1
      ,j3))*Yv(j1,j2))*ZV(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gI1,j2))*
      SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSdconjSd(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(Lambdax)
      *(vS*KroneckerDelta(1,gO2) + vu*KroneckerDelta(2,gO2))*SUM(j2,0,2,Conj(ZD(
      gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) - (vS*KroneckerDelta(1,gO2) +
      vu*KroneckerDelta(2,gO2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj
      (ZD(gI2,3 + j1)))*ZD(gI1,j2)) + 1.4142135623730951*KroneckerDelta(0,gO2)*(
      SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,ZD(gI1,3 + j1)*TYd(j1,j2))) - SUM(j2,
      0,2,SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2)))*ZD(gI1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSeconjSe(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(Lambdax)
      *(vS*KroneckerDelta(1,gO2) + vu*KroneckerDelta(2,gO2))*SUM(j2,0,2,Conj(ZE(
      gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) - (vS*KroneckerDelta(1,gO2) +
      vu*KroneckerDelta(2,gO2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj
      (ZE(gI2,3 + j1)))*ZE(gI1,j2)) + 1.4142135623730951*KroneckerDelta(0,gO2)*(
      SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,ZE(gI1,3 + j1)*TYe(j1,j2))) - SUM(j2,
      0,2,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2)))*ZE(gI1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSuconjSu(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(Conj(Lambdax)
      *(vS*KroneckerDelta(0,gO2) + vd*KroneckerDelta(2,gO2))*SUM(j2,0,2,Conj(ZU(
      gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) - (vS*KroneckerDelta(0,gO2) +
      vd*KroneckerDelta(2,gO2))*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj
      (ZU(gI2,3 + j1)))*ZU(gI1,j2)) + 1.4142135623730951*KroneckerDelta(1,gO2)*(
      SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,ZU(gI1,3 + j1)*TYu(j1,j2))) - SUM(j2,
      0,2,SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2)))*ZU(gI1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSvconjSv(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(vS*
      KroneckerDelta(0,gO2)*(Lambdax*SUM(j2,0,2,Conj(ZV(gI2,3 + j2))*SUM(j1,0,2,
      Conj(Yv(j1,j2))*ZV(gI1,j1))) - Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(
      gI2,j1))*Yv(j1,j2))*ZV(gI1,3 + j2))) + vd*KroneckerDelta(2,gO2)*(Lambdax*SUM
      (j2,0,2,Conj(ZV(gI2,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gI1,j1))) - Conj(
      Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gI2,j1))*Yv(j1,j2))*ZV(gI1,3 + j2)))
      + 1.4142135623730951*KroneckerDelta(1,gO2)*(SUM(j2,0,2,Conj(ZV(gI2,3 + j2))*
      SUM(j1,0,2,Conj(TYv(j1,j2))*ZV(gI1,j1))) - SUM(j2,0,2,SUM(j1,0,2,Conj(ZV(gI2
      ,j1))*TYv(j1,j2))*ZV(gI1,3 + j2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUAhPR(int gI1, int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,-0.1)*(5*
      KroneckerDelta(2,gO2)*(2*gp*Qs*ZN(gI1,5)*ZN(gI2,0) + 1.4142135623730951*Conj
      (Lambdax)*(ZN(gI1,4)*ZN(gI2,3) + ZN(gI1,3)*ZN(gI2,4)) + 2*gp*Qs*ZN(gI1,0)*ZN
      (gI2,5)) + KroneckerDelta(0,gO2)*(ZN(gI1,3)*(10*gp*QHd*ZN(gI2,0) -
      3.872983346207417*g1*ZN(gI2,1) + 5*g2*ZN(gI2,2)) + 10*gp*QHd*ZN(gI1,0)*ZN(
      gI2,3) - 3.872983346207417*g1*ZN(gI1,1)*ZN(gI2,3) + 5*g2*ZN(gI1,2)*ZN(gI2,3)
      + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,5)*ZN(gI2,4) + 7.0710678118654755*
      Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,5)) + KroneckerDelta(1,gO2)*(ZN(gI1,4)*(10*gp
      *QHu*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(gI2,2)) + (10*gp*
      QHu*ZN(gI1,0) + 3.872983346207417*g1*ZN(gI1,1) - 5*g2*ZN(gI1,2))*ZN(gI2,4) +
      7.0710678118654755*Conj(Lambdax)*(ZN(gI1,5)*ZN(gI2,3) + ZN(gI1,3)*ZN(gI2,5))
      ));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUAhPL(int gI1, int gI2, int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,0.1)*(-
      3.872983346207417*g1*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) +
      5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) + 10*gp*QHu*Conj(
      ZN(gI1,4))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj
      (ZN(gI1,4))*Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,4))*
      Conj(ZN(gI2,2))*KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,1))
      *Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,4)
      )*KroneckerDelta(1,gO1) + 10*gp*Qs*Conj(ZN(gI1,5))*Conj(ZN(gI2,0))*
      KroneckerDelta(2,gO1) + 10*gp*Conj(ZN(gI1,0))*(QHd*Conj(ZN(gI2,3))*
      KroneckerDelta(0,gO1) + QHu*Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) + Qs*Conj(
      ZN(gI2,5))*KroneckerDelta(2,gO1)) + 7.0710678118654755*Conj(ZN(gI1,5))*Conj(
      ZN(gI2,4))*KroneckerDelta(0,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,4)
      )*Conj(ZN(gI2,5))*KroneckerDelta(0,gO1)*Lambdax + 7.0710678118654755*Conj(ZN
      (gI1,5))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)*Lambdax + 7.0710678118654755*
      Conj(ZN(gI1,4))*Conj(ZN(gI2,3))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,
      3))*(10*gp*QHd*Conj(ZN(gI2,0))*KroneckerDelta(0,gO1) - 3.872983346207417*g1*
      Conj(ZN(gI2,1))*KroneckerDelta(0,gO1) + 5*g2*Conj(ZN(gI2,2))*KroneckerDelta(
      0,gO1) + 7.0710678118654755*Conj(ZN(gI2,5))*KroneckerDelta(1,gO1)*Lambdax +
      7.0710678118654755*Conj(ZN(gI2,4))*KroneckerDelta(2,gO1)*Lambdax));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,0.1)*(10*gp*Qs*Conj(
      ZH(gI2,2))*KroneckerDelta(2,gO2)*Sin(ThetaWp()) + Conj(ZH(gI2,0))*
      KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417
      *g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) - Conj(ZH(gI2,1
      ))*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhVZp(int gO2, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,0.1)*(10*gp*Qs*Conj(
      ZH(gI2,2))*Cos(ThetaWp())*KroneckerDelta(2,gO2) + Conj(ZH(gI2,1))*
      KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(
      ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())) + Conj(ZH(
      gI2,0))*KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()
      ) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgZUHpm(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*g2*(vd*KroneckerDelta(0,gO2)*(5*g2*Cos
      (ThetaW())*Cos(ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()
      ) - 10*gp*QHd*Sin(ThetaWp())) + vu*KroneckerDelta(1,gO2)*(-5*g2*Cos(ThetaW()
      )*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp
      *QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmconjUHpm(int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*(-(g2*vd*KroneckerDelta(0,gO1)*(5*g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 10*gp*QHd*Sin(ThetaWp()))) + g2*vu*KroneckerDelta(1,gO1)*(5*g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) - 10*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgZconjUHpm(int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*g2*(vd*KroneckerDelta(0,gO1)*(5*g2*Cos
      (ThetaW())*Cos(ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()
      ) - 10*gp*QHd*Sin(ThetaWp())) + vu*KroneckerDelta(1,gO1)*(-5*g2*Cos(ThetaW()
      )*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp
      *QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmCUHpm(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*(-(g2*vd*KroneckerDelta(0,gO2)*(5*g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 10*gp*QHd*Sin(ThetaWp()))) + g2*vu*KroneckerDelta(1,gO2)*(5*g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) - 10*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgZpUHpm(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -0.05*g2*(vd*KroneckerDelta(0,gO2)*(10*gp*
      QHd*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()
      ))*Sin(ThetaWp())) + vu*KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (-
      5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgWmconjUHpm(int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*g2*(vd*KroneckerDelta(0,gO1)*(-10*gp*
      QHd*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()
      ))*Sin(ThetaWp())) - vu*KroneckerDelta(1,gO1)*(10*gp*QHu*Cos(ThetaWp()) + (5
      *g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgZpconjUHpm(int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -0.05*g2*(vd*KroneckerDelta(0,gO1)*(10*gp*
      QHd*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()
      ))*Sin(ThetaWp())) + vu*KroneckerDelta(1,gO1)*(10*gp*QHu*Cos(ThetaWp()) + (-
      5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZpgWmCUHpm(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*g2*(vd*KroneckerDelta(0,gO2)*(-10*gp*
      QHd*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()
      ))*Sin(ThetaWp())) - vu*KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (5
      *g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*g2*(vd*KroneckerDelta(0,gO2)*(
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())
      ) + vu*KroneckerDelta(1,gO2)*(-3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 10*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVZp(int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*g2*(vd*KroneckerDelta(0,gO2)*(10*gp*QHd
      *Cos(ThetaWp()) - 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())) + vu*
      KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + 3.872983346207417*g1*Sin(
      ThetaW())*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2)*(15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp
      ()) - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1*Cos(ThetaWp())
      *Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*
      Sqr(Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 20*
      Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)*(-15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1*Cos(
      ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(
      ThetaW())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZpVZp(int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2)*(2*gp*QHd*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*
      Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (-
      7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(ThetaWp()))) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-2*gp*QHu*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())
      + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(
      ThetaWp()))));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*(
      KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*
      QHu*Sqr(gp))*ZP(gI1,1)*ZP(gI2,0) + KroneckerDelta(0,gO2)*(2*(3*Sqr(g1) + 5*
      Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI1,0)*ZP(gI2,0) + (20*AbsSqr(Lambdax) - 3
      *Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI1,1)*ZP(gI2,1)))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZP(gI1,0)*ZP(gI2,1) + KroneckerDelta(1,
      gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZP(
      gI1,0)*ZP(gI2,0) - 2*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZP(gI1,1)
      *ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(vu*Conj(ZA(
      gI2,0))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(KroneckerDelta(1,gO2)*ZP(gI1,0) -
      KroneckerDelta(0,gO2)*ZP(gI1,1)) + vd*Conj(ZA(gI2,1))*(-2*AbsSqr(Lambdax) +
      Sqr(g2))*(KroneckerDelta(1,gO2)*ZP(gI1,0) - KroneckerDelta(0,gO2)*ZP(gI1,1))
      + 2.8284271247461903*Conj(ZA(gI2,2))*(-(KroneckerDelta(1,gO2)*TLambdax*ZP(
      gI1,0)) + Conj(TLambdax)*KroneckerDelta(0,gO2)*ZP(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-(Conj(ZH(gI2,0))*(KroneckerDelta(0,
      gO2)*(vd*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI1,0) + 5*vu*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) + KroneckerDelta(1,gO2)*(5*vu*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*
      QHu*Sqr(gp))*ZP(gI1,1)))) - 10*Conj(ZH(gI2,2))*(KroneckerDelta(0,gO2)*(2*vS*
      (AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gI1,0) + 1.4142135623730951*Conj(
      TLambdax)*ZP(gI1,1)) + KroneckerDelta(1,gO2)*(1.4142135623730951*TLambdax*ZP
      (gI1,0) + 2*vS*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gI1,1))) + Conj(ZH(gI2,
      1))*(KroneckerDelta(0,gO2)*(vu*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp)))
      *ZP(gI1,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) - KroneckerDelta
      (1,gO2)*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) + vu*(3*Sqr(g1) + 5*
      Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZP(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*(
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp
      )) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*
      Sqr(gp))) - Conj(ZA(gI1,0))*(-5*Conj(ZA(gI2,1))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*
      AbsSqr(Lambdax) + Sqr(g2)) + Conj(ZA(gI2,0))*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(
      gp)*Sqr(QHd)))) + Conj(ZA(gI1,1))*(5*Conj(ZA(gI2,0))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*
      AbsSqr(Lambdax) + Sqr(g2)) + Conj(ZA(gI2,1))*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp))) -
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(
      gp)*Sqr(QHu))))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.05*(-20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*(
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp
      )) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*
      Sqr(gp))) - Conj(ZH(gI1,0))*(5*Conj(ZH(gI2,1))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*
      AbsSqr(Lambdax) + Sqr(g2)) + Conj(ZH(gI2,0))*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(
      gp)*Sqr(QHd)))) + Conj(ZH(gI1,1))*(-5*Conj(ZH(gI2,0))*(KroneckerDelta(0,gO2)
      *KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*
      AbsSqr(Lambdax) + Sqr(g2)) + Conj(ZH(gI2,1))*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp))) -
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(
      gp)*Sqr(QHu))))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHpmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHpmPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(ZDL(
      gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjUHpmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZVL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjUHpmPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(ZVR(
      gI1,j2))*SUM(j1,0,2,Conj(ZEL(gI2,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSdconjUHpmconjSd(int gO1, int gI1, int gO2, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)*((Sqr(g1) - 5*(Sqr(g2) + 4*QHd*Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZD(gI1,j1
      ))*ZD(gI2,j1)) + 2*(Sqr(g1) - 10*Qd*QHd*Sqr(gp))*SUM(j1,0,2,Conj(ZD(gI1,3 +
      j1))*ZD(gI2,3 + j1)) - 20*SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,
      0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI2,3 + j2)))) - KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2) + 20*QHu*Qq*Sqr(gp))*SUM(j1,0,2,
      Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*(Sqr(g1) + 10*Qd*QHu*Sqr(gp))*SUM(j1,0,2,
      Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) + 20*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI1,
      j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gI2,j3))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSeconjUHpmconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*(Sqr(g2) + 4*QHd*Ql*Sqr(gp)))*SUM(j1,0
      ,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + (-6*Sqr(g1) + 20*Qe*QHd*Sqr(gp))*SUM(j1,0,
      2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) + 20*SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2))))) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*(Sqr(g2) - 4*QHu
      *Ql*Sqr(gp)))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 2*((3*Sqr(g1) + 10*
      Qe*QHu*Sqr(gp))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) + 10*SUM(j3,
      0,2,SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZE(
      gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSuconjUHpmconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.05*(-(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) + 5*(Sqr(g2) + 4*QHu*Qq*Sqr(gp)))*SUM(j1,0,2
      ,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*(Sqr(g1) - 5*QHu*Qu*Sqr(gp))*SUM(j1,0,2,
      Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) + 20*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI2,3 + j2))))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((Sqr(g1) + 5*(Sqr(g2) - 4*QHd*
      Qq*Sqr(gp)))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*((Sqr(g1) + 5*QHd*
      Qu*Sqr(gp))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) + 5*SUM(j3,0,2,
      SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gI2,j3
      )))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSvconjUHpmconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(
      1,gO2)*((3*Sqr(g1) - 5*(Sqr(g2) + 4*QHu*Ql*Sqr(gp)))*SUM(j1,0,2,Conj(ZV(gI1,
      j1))*ZV(gI2,j1)) - 20*(QHu*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2
      ,3 + j1)) + SUM(j3,0,2,Conj(ZV(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1
      ,j3))*Yv(j1,j2))*ZV(gI2,3 + j2))))) - KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2)*((3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*Ql*Sqr(gp))*SUM(j1,0,2,Conj(ZV(gI1,j1
      ))*ZV(gI2,j1)) + 20*(QHd*Qv*Sqr(gp)*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3
       + j1)) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*
      Ye(j1,j2)))*ZV(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjUHpmPR(int gI1, int gI2, int gO2) const
{
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -0.1*KroneckerDelta(1,gO2)*(
      1.4142135623730951*UP(gI2,1)*(10*gp*QHu*ZN(gI1,0) + 3.872983346207417*g1*ZN(
      gI1,1) + 5*g2*ZN(gI1,2)) + 10*g2*UP(gI2,0)*ZN(gI1,4)) - Conj(Lambdax)*
      KroneckerDelta(0,gO2)*UP(gI2,1)*ZN(gI1,5);

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjUHpmPL(int gI1, int gI2, int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Conj(ZN(gI1,3))*
      KroneckerDelta(0,gO1)) + Conj(UM(gI2,1))*(-1.4142135623730951*gp*QHd*Conj(ZN
      (gI1,0))*KroneckerDelta(0,gO1) + 0.5477225575051661*g1*Conj(ZN(gI1,1))*
      KroneckerDelta(0,gO1) + 0.7071067811865475*g2*Conj(ZN(gI1,2))*KroneckerDelta
      (0,gO1) - Conj(ZN(gI1,5))*KroneckerDelta(1,gO1)*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUHpmconjSu(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(KroneckerDelta(0,gO2)*(-
      1.4142135623730951*vd*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 2*(
      1.4142135623730951*vS*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      Yu(j1,j2)*ZU(gI1,3 + j1))) + 2*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      Conj(TYd(j1,j2)))*ZU(gI1,j2)) + 1.4142135623730951*vu*SUM(j3,0,2,Conj(ZD(gI2
      ,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2)))
      + 1.4142135623730951*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gI1,j3)))) + KroneckerDelta(1,gO2)*(-
      1.4142135623730951*vu*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 4*
      SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,ZU(gI1,3 + j1)*TYu(j1,j2))) +
      2.8284271247461903*(vS*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD
      (gI2,3 + j1)))*ZU(gI1,j2)) + vd*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2))) + vu*SUM(j3,0,2,SUM(
      j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI1,j3))))
      );

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUHpmconjSv(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(KroneckerDelta(0,gO2)*(-
      1.4142135623730951*vd*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZV(gI1,j1)) + 4*
      SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2)))*ZV(gI1,j2)) +
      2.8284271247461903*(vS*Conj(Lambdax)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gI2,j1))*
      Yv(j1,j2))*ZV(gI1,3 + j2)) + vu*SUM(j3,0,2,Conj(ZE(gI2,3 + j3))*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))*ZV(gI1,3 + j2))) + vd*SUM(j3,0,2,SUM(
      j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gI1,j3))))
      + KroneckerDelta(1,gO2)*(-1.4142135623730951*vu*Sqr(g2)*SUM(j1,0,2,Conj(ZE(
      gI2,j1))*ZV(gI1,j1)) + 2*(1.4142135623730951*vS*Lambdax*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZV(gI1,j2)) + 2*SUM(j2,0,2,SUM(j1,0,
      2,Conj(ZE(gI2,j1))*TYv(j1,j2))*ZV(gI1,3 + j2)) + 1.4142135623730951*vd*SUM(
      j3,0,2,Conj(ZE(gI2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Yv(j1,j2))
      *ZV(gI1,3 + j2))) + 1.4142135623730951*vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gI2,
      j2))*SUM(j1,0,2,Conj(Yv(j3,j1))*Yv(j2,j1)))*ZV(gI1,j3)))));

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
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*(KroneckerDelta(0,gO2)*(-5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 10*gp*QHd*Sin(ThetaWp()))*ZP(gI2,0) + KroneckerDelta(1,gO2)*(-5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      - 10*gp*QHu*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmVZp(int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*(KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*ZP(gI2,0) - KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (-
      5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZP(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjUHpmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(Conj(ZA(gI2
      ,0))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,1))*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CphhconjUHpmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.5*g2*(Conj(ZH(gI2,0))*KroneckerDelta(0,
      gO2) - Conj(ZH(gI2,1))*KroneckerDelta(1,gO2));

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

std::complex<double> CLASSNAME::CpHpmconjHpmVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + Cos(2*ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) + 5*Sqr(g2))*(ZP(
      gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpHpmconjHpmVP(int gI2, int gI1) const
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

std::complex<double> CLASSNAME::CpSdconjSdVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.03333333333333333*((-7.745966692414834*g1
      *g2*Cos(ThetaW())*Sin(ThetaW()) + Sqr(g1)*Sqr(Cos(ThetaW())) + 15*Sqr(g2)*
      Sqr(Sin(ThetaW())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*Sqr(g1)*Sqr(
      Cos(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

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

std::complex<double> CLASSNAME::CpHpmconjHpmVZVZ(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*((15.491933384829668*g1*gp*QHd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*
      (3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp()
      )) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp())))*ZP
      (gI1,0)*ZP(gI2,0) + (-15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW
      ())*Sin(ThetaWp()) - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1
      *Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(
      Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(
      ThetaW())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp())))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZ(int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*((-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())
      )*ZP(gI1,0)*ZP(gI2,0) + (-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      )*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZPL(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*Cos(
      ThetaWp())*UM(gI1,0) + 0.5*Conj(UM(gI2,1))*(g2*Cos(ThetaW())*Cos(ThetaWp())
      - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHd*Sin(ThetaWp(
      )))*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZPR(int gI1, int gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = g2*Conj(UP(gI1,0))*Cos(ThetaW())*Cos(
      ThetaWp())*UP(gI2,0) + 0.5*Conj(UP(gI1,1))*(g2*Cos(ThetaW())*Cos(ThetaWp())
      - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHu*Sin(ThetaWp(
      )))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr
      (gp)*Sqr(Qs)*Sqr(Sin(ThetaWp())) + Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(20*g2*gp
      *QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*
      QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr
      (Sin(ThetaWp()))) + Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(-20*g2*gp*QHu*Cos(
      ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr
      (Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr
      (gp)*Sqr(Qs)*Sqr(Sin(ThetaWp())) + Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*(20*g2*gp
      *QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*
      QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr
      (Sin(ThetaWp()))) + Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*(-20*g2*gp*QHu*Cos(
      ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr
      (Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ(int gI2, int gI1) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = std::complex<double>(0,0.5)*(2*gp*Qs*Conj(
      ZA(gI2,2))*Conj(ZH(gI1,2))*Sin(ThetaWp()) + Conj(ZA(gI2,0))*Conj(ZH(gI1,0))*
      (g2*Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 2*gp*QHd*Sin(ThetaWp())) - Conj(ZA(gI2,1))*Conj(ZH(gI1,1))*(g2*
      Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(
      ThetaW()) - 2*gp*QHu*Sin(ThetaWp())));

   return result;
}

double CLASSNAME::CpbarFdFdVZPL(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const double result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()
      ) - 6*gp*Qq*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdFdVZPR(int gI1, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   const double result = KroneckerDelta(gI1,gI2)*(-0.2581988897471611*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + gp*Qd*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPL(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW())*Cos(ThetaWp
      ()) - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*Ql*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPR(int gI1, int gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   const double result = -(KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) - gp*Qe*Sin(ThetaWp())));

   return result;
}

double CLASSNAME::CpbarFuFuVZPL(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(
      ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()
      ) + 6*gp*Qq*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPR(int gI1, int gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   const double result = KroneckerDelta(gI1,gI2)*(0.5163977794943222*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + gp*Qu*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZPL(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const double result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZPR(int gI1, int gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   const double result = gp*Qv*KroneckerDelta(gI1,gI2)*Sin(ThetaWp());

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZVZ(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.03333333333333333*((-60*g2*gp*Qq*Cos(
      ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) +
      15*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 60*Sqr(gp)*Sqr(Qq)*Sqr(
      Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*(-
      7.745966692414834*g1*gp*Qd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + Sqr
      (g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 15*Sqr(gp)*Sqr(Qd)*Sqr(Sin(
      ThetaWp())))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZVZ(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.1*((15.491933384829668*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*
      (3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*Ql*Sin(ThetaWp())
      ) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) + 20*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())))*SUM
      (j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 4*(-7.745966692414834*g1*gp*Qe*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(
      Sin(ThetaW())) + 5*Sqr(gp)*Sqr(Qe)*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZE(
      gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZVZ(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.03333333333333333*((-15.491933384829668*
      g1*gp*Qq*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 2*g2*Cos(ThetaW())*
      Cos(ThetaWp())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 30*gp*Qq
      *Sin(ThetaWp())) + 15*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + Sqr(
      g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 60*Sqr(gp)*Sqr(Qq)*Sqr(Sin(
      ThetaWp())))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + 4*(4*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) + gp*Qu*(7.745966692414834*g1*Sin(ThetaW())*
      Sin(2*ThetaWp()) + 15*gp*Qu*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZU(gI1,3 +
      j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZVZ(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.1*(15.491933384829668*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 10*g2*gp*Ql*Cos(ThetaW())*Sin(2*
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZV(gI1
      ,j1))*ZV(gI2,j1)) + 2*Sqr(gp)*Sqr(Qv)*Sqr(Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZV
      (gI1,3 + j1))*ZV(gI2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZ(int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.16666666666666666*(-((3*g2*Cos(ThetaW())*
      Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 6*gp*
      Qq*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))) + 2*(
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 3*gp*Qd*Sin(ThetaWp()))
      *SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZ(int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.1*((-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*Ql*Sin(ThetaWp()))
      *SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) + 2*(3.872983346207417*g1*Cos(
      ThetaWp())*Sin(ThetaW()) - 5*gp*Qe*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZE(gI2,3
      + j1))*ZE(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZ(int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.16666666666666666*((3*g2*Cos(ThetaW())*
      Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 6*gp*
      Qq*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 2*(
      1.5491933384829668*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3*gp*Qu*Sin(ThetaWp()))
      *SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZ(int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.5*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*Sin(ThetaWp()))
      *SUM(j1,0,2,Conj(ZV(gI2,j1))*ZV(gI1,j1)) - gp*Qv*Sin(ThetaWp())*SUM(j1,0,2,
      Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZPL(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.5*(-(Conj(ZN(gI2,3))*(g2*Cos(ThetaW())*
      Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*
      QHd*Sin(ThetaWp()))*ZN(gI1,3)) + Conj(ZN(gI2,4))*(g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*
      Sin(ThetaWp()))*ZN(gI1,4) - 2*gp*Qs*Conj(ZN(gI2,5))*Sin(ThetaWp())*ZN(gI1,5)
      );

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZPR(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.5*(Conj(ZN(gI1,3))*(g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*
      Sin(ThetaWp()))*ZN(gI2,3) - Conj(ZN(gI1,4))*(g2*Cos(ThetaW())*Cos(ThetaWp())
      + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*Sin(ThetaWp(
      )))*ZN(gI2,4) + 2*gp*Qs*Conj(ZN(gI1,5))*Sin(ThetaWp())*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjVWmVZ(int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*Sin(ThetaWp()))*ZP(gI2,0) + vu*(-
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHu*Sin(ThetaWp())
      )*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhVZVZ(int gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.5*(4*vS*Conj(ZH(gI2,2))*Sqr(gp)*Sqr(Qs)*
      Sqr(Sin(ThetaWp())) + vd*Conj(ZH(gI2,0))*Sqr(g2*Cos(ThetaW())*Cos(ThetaWp())
      + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*Sin(ThetaWp(
      ))) + vu*Conj(ZH(gI2,1))*Sqr(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*Sin(ThetaWp())
      ));

   return result;
}

std::complex<double> CLASSNAME::CphhVZVZp(int gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*(20*vS*Conj(ZH(gI2,2))*Cos(ThetaWp())*
      Sin(ThetaWp())*Sqr(gp)*Sqr(Qs) - vd*Conj(ZH(gI2,0))*(5*Cos(ThetaWp())*Sin(
      ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) - 7.745966692414834*g1*gp*QHd*Sin(
      ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*
      Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) + 7.745966692414834*g1*gp*QHd*Sin(
      ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos
      (ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*
      gp*QHd*Sqr(Sin(ThetaWp())))) - vu*Conj(ZH(gI2,1))*(5*Cos(ThetaWp())*Sin(
      ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) + 7.745966692414834*g1*gp*QHu*Sin(
      ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*
      Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) - 7.745966692414834*g1*gp*QHu*Sin(
      ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos
      (ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*
      gp*QHu*Sqr(Sin(ThetaWp())))));

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

std::complex<double> CLASSNAME::CpHpmconjHpmVZpVZp(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.1*((2*gp*QHd*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Cos(ThetaWp())) + (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())
      + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(
      ThetaWp())))*ZP(gI1,0)*ZP(gI2,0) + (-2*gp*QHu*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())
      + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(
      ThetaWp())))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZp(int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.5*((2*gp*QHd*Cos(ThetaWp()) + g2*Cos(
      ThetaW())*Sin(ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()
      ))*ZP(gI1,0)*ZP(gI2,0) + (-2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(
      ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI1,1)*
      ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZpPL(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Cos(ThetaW())*Sin(
      ThetaWp())*UM(gI1,0)) - 0.5*Conj(UM(gI2,1))*(2*gp*QHd*Cos(ThetaWp()) + g2*
      Cos(ThetaW())*Sin(ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(
      ThetaWp()))*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZpPR(int gI1, int gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*Cos(ThetaW())*Sin(
      ThetaWp())*UP(gI2,0)) + 0.1*Conj(UP(gI1,1))*(10*gp*QHu*Cos(ThetaWp()) + (-5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*UP(
      gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZpVZp(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr
      (gp)*Sqr(Qs)*Sqr(Cos(ThetaWp())) + Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(-2*gp*
      QHd*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp(
      )) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp()))) + Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(2*gp*
      QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp(
      )) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZpVZp(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr
      (gp)*Sqr(Qs)*Sqr(Cos(ThetaWp())) + Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*(-2*gp*
      QHd*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp(
      )) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp()))) + Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*(2*gp*
      QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp(
      )) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZp(int gI2, int gI1) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = std::complex<double>(0,0.5)*(2*gp*Qs*Conj(
      ZA(gI2,2))*Conj(ZH(gI1,2))*Cos(ThetaWp()) + Conj(ZA(gI2,1))*Conj(ZH(gI1,1))*
      (2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())) + Conj(ZA(gI2,0))*Conj(
      ZH(gI1,0))*(2*gp*QHd*Cos(ThetaWp()) - (g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

double CLASSNAME::CpbarFdFdVZpPL(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(6*gp*Qq*Cos
      (ThetaWp()) + 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(
      ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdFdVZpPR(int gI1, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   const double result = KroneckerDelta(gI1,gI2)*(gp*Qd*Cos(ThetaWp()) +
      0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZpPL(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const double result = -0.5*KroneckerDelta(gI1,gI2)*(2*gp*Ql*Cos(ThetaWp()) + g2
      *Cos(ThetaW())*Sin(ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(
      ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeFeVZpPR(int gI1, int gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   const double result = KroneckerDelta(gI1,gI2)*(gp*Qe*Cos(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZpPL(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(6*gp*Qq*Cos
      (ThetaWp()) - 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(
      ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuFuVZpPR(int gI1, int gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   const double result = KroneckerDelta(gI1,gI2)*(gp*Qu*Cos(ThetaWp()) -
      0.5163977794943222*g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZpPL(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const double result = -0.5*KroneckerDelta(gI1,gI2)*(2*gp*Ql*Cos(ThetaWp()) - (
      g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFvFvVZpPR(int gI1, int gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   const double result = gp*Qv*Cos(ThetaWp())*KroneckerDelta(gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZpVZp(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.03333333333333333*((2*gp*Qq*(15*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 60*Sqr(gp
      )*Sqr(Qq)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + g1*Sin(ThetaW())) + 15*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(
      ThetaWp())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*(15*Sqr(gp)*Sqr(Qd)
      *Sqr(Cos(ThetaWp())) + g1*Sin(ThetaW())*(3.872983346207417*gp*Qd*Sin(2*
      ThetaWp()) + g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZD(gI1,3
       + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZpVZp(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.1*((2*gp*Ql*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(Ql)*
      Sqr(Cos(ThetaWp())) + (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())
      + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(
      ThetaWp())))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 4*(5*Sqr(gp)*Sqr(Qe)*
      Sqr(Cos(ThetaWp())) + g1*Sin(ThetaW())*(3.872983346207417*gp*Qe*Sin(2*
      ThetaWp()) + 3*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZE(gI1
      ,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZpVZp(int gI1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.03333333333333333*((-2*gp*Qq*(15*g2*Cos(
      ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 60*Sqr(gp
      )*Sqr(Qq)*Sqr(Cos(ThetaWp())) + (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(
      ThetaW()) + 15*Sqr(g2)*Sqr(Cos(ThetaW())) + Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(
      Sin(ThetaWp())))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + 4*(15*Sqr(gp)*Sqr
      (Qu)*Sqr(Cos(ThetaWp())) + 2*g1*Sin(ThetaW())*(-3.872983346207417*gp*Qu*Sin(
      2*ThetaWp()) + 2*g1*Sin(ThetaW())*Sqr(Sin(ThetaWp()))))*SUM(j1,0,2,Conj(ZU(
      gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZpVZp(int gI1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = 0.1*(-2*gp*Ql*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(Ql)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp())))*
      SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(gI2,j1)) + 2*Sqr(gp)*Sqr(Qv)*Sqr(Cos(ThetaWp(
      )))*SUM(j1,0,2,Conj(ZV(gI1,3 + j1))*ZV(gI2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZp(int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = 0.16666666666666666*((6*gp*Qq*Cos(ThetaWp()
      ) + 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) - 2*(3*gp*Qd*Cos(
      ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM(j1,0,2,
      Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZp(int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = 0.5*(2*gp*Ql*Cos(ThetaWp()) + g2*Cos(ThetaW
      ())*Sin(ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM
      (j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) - (gp*Qe*Cos(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZE(gI2,3
       + j1))*ZE(gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZp(int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = 0.16666666666666666*((6*gp*Qq*Cos(ThetaWp()
      ) - 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) + 2*(-3*gp*Qu*Cos(
      ThetaWp()) + 1.5491933384829668*g1*Sin(ThetaW())*Sin(ThetaWp()))*SUM(j1,0,2,
      Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZp(int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = (gp*Ql*Cos(ThetaWp()) - 0.1*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*SUM(j1,0,2,
      Conj(ZV(gI2,j1))*ZV(gI1,j1)) - gp*Qv*Cos(ThetaWp())*SUM(j1,0,2,Conj(ZV(gI2,3
       + j1))*ZV(gI1,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZpPL(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.5*(-(Conj(ZN(gI2,3))*(2*gp*QHd*Cos(
      ThetaWp()) - (g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*ZN(gI1,3)) - Conj(ZN(gI2,4))*(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(
      ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()
      ))*ZN(gI1,4) - 2*gp*Qs*Conj(ZN(gI2,5))*Cos(ThetaWp())*ZN(gI1,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZpPR(int gI1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.5*(Conj(ZN(gI1,3))*(2*gp*QHd*Cos(ThetaWp(
      )) - (g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp())
      )*ZN(gI2,3) + Conj(ZN(gI1,4))*(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*
      Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZN(gI2,
      4) + 2*gp*Qs*Conj(ZN(gI1,5))*Cos(ThetaWp())*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjVWmVZp(int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp()) -
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI2,0) + vu*(2*gp*QHu
      *Cos(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI2
      ,1));

   return result;
}

std::complex<double> CLASSNAME::CphhVZpVZp(int gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = 0.5*(4*vS*Conj(ZH(gI2,2))*Sqr(gp)*Sqr(Qs)*
      Sqr(Cos(ThetaWp())) + vd*Conj(ZH(gI2,0))*Sqr(-2*gp*QHd*Cos(ThetaWp()) + g2*
      Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(
      ThetaWp())) + vu*Conj(ZH(gI2,1))*Sqr(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW
      ())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())));

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

std::complex<double> CLASSNAME::CpAhHpmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(Conj(ZA(gI2
      ,0))*ZP(gI1,0) + Conj(ZA(gI2,1))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5*g2*(Conj(ZH(gI2,0))*ZP(gI1,0) - Conj(
      ZH(gI2,1))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)) + Conj
      (ZA(gI1,1))*Conj(ZA(gI2,1)))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CphhhhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(Conj(ZH(gI1,0))*Conj(ZH(gI2,0)) + Conj
      (ZH(gI1,1))*Conj(ZH(gI2,1)))*Sqr(g2);

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
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZEL(
      gI2,j1))*ZVL(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFvFeconjVWmPR(int , int ) const
{
   
   const double result = 0;

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

std::complex<double> CLASSNAME::CpSvconjSvconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZV(
      gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,2) +
      1.4142135623730951*Conj(UM(gI2,1))*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjVWmPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(ZN(gI1,2))*UP(gI2,0)) +
      0.7071067811865475*g2*Conj(ZN(gI1,4))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSuconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*ZU(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSvconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZE(
      gI2,j1))*ZV(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CphhconjVWmVWm(int gI2) const
{
   
   const std::complex<double> result = 0.5*(vd*Conj(ZH(gI2,0)) + vu*Conj(ZH(gI2,1)
      ))*Sqr(g2);

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
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(4,gO2)*
      ZP(gI2,1)) - 0.1*Conj(UP(gI1,1))*(10*KroneckerDelta(5,gO2)*Lambdax*ZP(gI2,0)
      + 1.4142135623730951*(10*gp*QHu*KroneckerDelta(0,gO2) + 3.872983346207417*g1
      *KroneckerDelta(1,gO2) + 5*g2*KroneckerDelta(2,gO2))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiHpmPR(int gI1, int gO1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = -(g2*KroneckerDelta(3,gO1)*UM(gI1,0)*ZP(gI2
      ,0)) + 0.1*UM(gI1,1)*(-14.142135623730951*gp*QHd*KroneckerDelta(0,gO1)*ZP(
      gI2,0) + 5.477225575051661*g1*KroneckerDelta(1,gO1)*ZP(gI2,0) +
      7.0710678118654755*g2*KroneckerDelta(2,gO1)*ZP(gI2,0) - 10*Conj(Lambdax)*
      KroneckerDelta(5,gO1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjHpmPL(int gO2, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*KroneckerDelta(3,gO2)*
      ZP(gI1,0)) + Conj(UM(gI2,1))*(-1.4142135623730951*gp*QHd*KroneckerDelta(0,
      gO2)*ZP(gI1,0) + 0.5477225575051661*g1*KroneckerDelta(1,gO2)*ZP(gI1,0) +
      0.7071067811865475*g2*KroneckerDelta(2,gO2)*ZP(gI1,0) - KroneckerDelta(5,gO2
      )*Lambdax*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjHpmPR(int gO1, int gI2, int gI1) const
{
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -(Conj(Lambdax)*KroneckerDelta(5,gO1)*UP(
      gI2,1)*ZP(gI1,0)) - 0.1*(10*g2*KroneckerDelta(4,gO1)*UP(gI2,0) +
      1.4142135623730951*(10*gp*QHu*KroneckerDelta(0,gO1) + 3.872983346207417*g1*
      KroneckerDelta(1,gO1) + 5*g2*KroneckerDelta(2,gO1))*UP(gI2,1))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiVWmPL(int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*KroneckerDelta(2,gO2)*UM(gI1,0)
      + 1.4142135623730951*KroneckerDelta(3,gO2)*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiVWmPR(int gI1, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(2,gO1))
      + 0.7071067811865475*g2*Conj(UP(gI1,1))*KroneckerDelta(4,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdUChiSdPL(int gI1, int gO2, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = -1.4142135623730951*gp*Qd*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1))) - 0.3651483716701107
      *g1*KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1)))
      - KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1,
      j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdUChiSdPR(int gI1, int gO1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = -1.4142135623730951*gp*Qq*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1)) - 0.18257418583505536*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1)) +
      0.7071067811865475*g2*KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(
      gI1,j1)) - KroneckerDelta(3,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(
      ZD(gI2,3 + j1)))*ZDL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUChiSePL(int gI1, int gO2, int gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = -1.4142135623730951*gp*Qe*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(ZER(gI1,j1))) - 1.0954451150103321
      *g1*KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(ZER(gI1,j1)))
      - KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1,
      j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUChiSePR(int gI1, int gO1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = -1.4142135623730951*gp*Ql*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(gI1,j1)) + 0.5477225575051661*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(gI1,j1)) +
      0.7071067811865475*g2*KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(
      gI1,j1)) - KroneckerDelta(3,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(
      ZE(gI2,3 + j1)))*ZEL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuUChiSuPL(int gI1, int gO2, int gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = -1.4142135623730951*gp*Qu*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1))) + 0.7302967433402214
      *g1*KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1)))
      - KroneckerDelta(4,gO2)*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,
      j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuUChiSuPR(int gI1, int gO1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = -1.4142135623730951*gp*Qq*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1)) - 0.18257418583505536*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1)) -
      0.7071067811865475*g2*KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(
      gI1,j1)) - KroneckerDelta(4,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(
      ZU(gI2,3 + j1)))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvUChiSvPL(int gI1, int gO2, int gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = -1.4142135623730951*gp*Qv*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZV(gI2,3 + j1))*Conj(ZVR(gI1,j1))) - KroneckerDelta(4,
      gO2)*SUM(j2,0,2,Conj(ZVR(gI1,j2))*SUM(j1,0,2,Conj(ZV(gI2,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvUChiSvPR(int gI1, int gO1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = -1.4142135623730951*gp*Ql*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZV(gI2,j1))*ZVL(gI1,j1)) + 0.5477225575051661*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZV(gI2,j1))*ZVL(gI1,j1)) -
      0.7071067811865475*g2*KroneckerDelta(2,gO1)*SUM(j1,0,2,Conj(ZV(gI2,j1))*ZVL(
      gI1,j1)) - KroneckerDelta(4,gO1)*SUM(j2,0,2,Conj(ZV(gI2,3 + j2))*SUM(j1,0,2,
      Conj(Yv(j1,j2))*ZVL(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChihhPL(int gI2, int gO2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(5*Conj(ZH(gI1,2))*(-2*gp*Qs*Conj(ZN(
      gI2,5))*KroneckerDelta(0,gO2) - 2*gp*Qs*Conj(ZN(gI2,0))*KroneckerDelta(5,gO2
      ) + 1.4142135623730951*Conj(ZN(gI2,4))*KroneckerDelta(3,gO2)*Lambdax +
      1.4142135623730951*Conj(ZN(gI2,3))*KroneckerDelta(4,gO2)*Lambdax) - Conj(ZH(
      gI1,1))*(Conj(ZN(gI2,4))*(10*gp*QHu*KroneckerDelta(0,gO2) +
      3.872983346207417*g1*KroneckerDelta(1,gO2) - 5*g2*KroneckerDelta(2,gO2)) +
      10*gp*QHu*Conj(ZN(gI2,0))*KroneckerDelta(4,gO2) + 3.872983346207417*g1*Conj(
      ZN(gI2,1))*KroneckerDelta(4,gO2) - 5*g2*Conj(ZN(gI2,2))*KroneckerDelta(4,gO2
      ) - 7.0710678118654755*Conj(ZN(gI2,5))*KroneckerDelta(3,gO2)*Lambdax -
      7.0710678118654755*Conj(ZN(gI2,3))*KroneckerDelta(5,gO2)*Lambdax) + Conj(ZH(
      gI1,0))*(Conj(ZN(gI2,3))*(-10*gp*QHd*KroneckerDelta(0,gO2) +
      3.872983346207417*g1*KroneckerDelta(1,gO2) - 5*g2*KroneckerDelta(2,gO2)) -
      10*gp*QHd*Conj(ZN(gI2,0))*KroneckerDelta(3,gO2) + 3.872983346207417*g1*Conj(
      ZN(gI2,1))*KroneckerDelta(3,gO2) - 5*g2*Conj(ZN(gI2,2))*KroneckerDelta(3,gO2
      ) + 7.0710678118654755*Conj(ZN(gI2,5))*KroneckerDelta(4,gO2)*Lambdax +
      7.0710678118654755*Conj(ZN(gI2,4))*KroneckerDelta(5,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChihhPR(int gI2, int gO1, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(5*Conj(ZH(gI1,2))*(-2*gp*Qs*
      KroneckerDelta(5,gO1)*ZN(gI2,0) + 1.4142135623730951*Conj(Lambdax)*(
      KroneckerDelta(4,gO1)*ZN(gI2,3) + KroneckerDelta(3,gO1)*ZN(gI2,4)) - 2*gp*Qs
      *KroneckerDelta(0,gO1)*ZN(gI2,5)) + Conj(ZH(gI1,0))*(KroneckerDelta(3,gO1)*(
      -10*gp*QHd*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(gI2,2)) - 10
      *gp*QHd*KroneckerDelta(0,gO1)*ZN(gI2,3) + 3.872983346207417*g1*
      KroneckerDelta(1,gO1)*ZN(gI2,3) - 5*g2*KroneckerDelta(2,gO1)*ZN(gI2,3) +
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(5,gO1)*ZN(gI2,4) +
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZN(gI2,5)) - Conj(ZH(
      gI1,1))*(KroneckerDelta(4,gO1)*(10*gp*QHu*ZN(gI2,0) + 3.872983346207417*g1*
      ZN(gI2,1) - 5*g2*ZN(gI2,2)) + (10*gp*QHu*KroneckerDelta(0,gO1) +
      3.872983346207417*g1*KroneckerDelta(1,gO1) - 5*g2*KroneckerDelta(2,gO1))*ZN(
      gI2,4) - 7.0710678118654755*Conj(Lambdax)*(KroneckerDelta(5,gO1)*ZN(gI2,3) +
      KroneckerDelta(3,gO1)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiAhPL(int gI1, int gO2, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,0.1)*(5*Conj(ZA(gI2,
      2))*(2*gp*Qs*Conj(ZN(gI1,5))*KroneckerDelta(0,gO2) + 2*gp*Qs*Conj(ZN(gI1,0))
      *KroneckerDelta(5,gO2) + 1.4142135623730951*Conj(ZN(gI1,4))*KroneckerDelta(3
      ,gO2)*Lambdax + 1.4142135623730951*Conj(ZN(gI1,3))*KroneckerDelta(4,gO2)*
      Lambdax) + Conj(ZA(gI2,1))*(Conj(ZN(gI1,4))*(10*gp*QHu*KroneckerDelta(0,gO2)
      + 3.872983346207417*g1*KroneckerDelta(1,gO2) - 5*g2*KroneckerDelta(2,gO2)) +
      10*gp*QHu*Conj(ZN(gI1,0))*KroneckerDelta(4,gO2) + 3.872983346207417*g1*Conj(
      ZN(gI1,1))*KroneckerDelta(4,gO2) - 5*g2*Conj(ZN(gI1,2))*KroneckerDelta(4,gO2
      ) + 7.0710678118654755*Conj(ZN(gI1,5))*KroneckerDelta(3,gO2)*Lambdax +
      7.0710678118654755*Conj(ZN(gI1,3))*KroneckerDelta(5,gO2)*Lambdax) + Conj(ZA(
      gI2,0))*(Conj(ZN(gI1,3))*(10*gp*QHd*KroneckerDelta(0,gO2) -
      3.872983346207417*g1*KroneckerDelta(1,gO2) + 5*g2*KroneckerDelta(2,gO2)) +
      10*gp*QHd*Conj(ZN(gI1,0))*KroneckerDelta(3,gO2) - 3.872983346207417*g1*Conj(
      ZN(gI1,1))*KroneckerDelta(3,gO2) + 5*g2*Conj(ZN(gI1,2))*KroneckerDelta(3,gO2
      ) + 7.0710678118654755*Conj(ZN(gI1,5))*KroneckerDelta(4,gO2)*Lambdax +
      7.0710678118654755*Conj(ZN(gI1,4))*KroneckerDelta(5,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiAhPR(int gI1, int gO1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = std::complex<double>(0,-0.1)*(5*Conj(ZA(gI2
      ,2))*(2*gp*Qs*KroneckerDelta(5,gO1)*ZN(gI1,0) + 1.4142135623730951*Conj(
      Lambdax)*(KroneckerDelta(4,gO1)*ZN(gI1,3) + KroneckerDelta(3,gO1)*ZN(gI1,4))
      + 2*gp*Qs*KroneckerDelta(0,gO1)*ZN(gI1,5)) + Conj(ZA(gI2,0))*(KroneckerDelta
      (3,gO1)*(10*gp*QHd*ZN(gI1,0) - 3.872983346207417*g1*ZN(gI1,1) + 5*g2*ZN(gI1,
      2)) + 10*gp*QHd*KroneckerDelta(0,gO1)*ZN(gI1,3) - 3.872983346207417*g1*
      KroneckerDelta(1,gO1)*ZN(gI1,3) + 5*g2*KroneckerDelta(2,gO1)*ZN(gI1,3) +
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(5,gO1)*ZN(gI1,4) +
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZN(gI1,5)) + Conj(ZA(
      gI2,1))*(KroneckerDelta(4,gO1)*(10*gp*QHu*ZN(gI1,0) + 3.872983346207417*g1*
      ZN(gI1,1) - 5*g2*ZN(gI1,2)) + (10*gp*QHu*KroneckerDelta(0,gO1) +
      3.872983346207417*g1*KroneckerDelta(1,gO1) - 5*g2*KroneckerDelta(2,gO1))*ZN(
      gI1,4) + 7.0710678118654755*Conj(Lambdax)*(KroneckerDelta(5,gO1)*ZN(gI1,3) +
      KroneckerDelta(3,gO1)*ZN(gI1,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFdconjSdPL(int gO2, int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = -1.4142135623730951*gp*Qq*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(gI1,j1)) - 0.18257418583505536*g1*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(gI1,j1)) +
      0.7071067811865475*g2*KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(
      gI1,j1)) - KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,Yd(
      j1,j2)*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFdconjSdPR(int gO1, int gI2, int gI1) const
{
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = -1.4142135623730951*gp*Qd*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZD(gI1,3 + j1)*ZDR(gI2,j1)) - 0.3651483716701107*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,ZD(gI1,3 + j1)*ZDR(gI2,j1)) -
      KroneckerDelta(3,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZD(
      gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFeconjSePL(int gO2, int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = -1.4142135623730951*gp*Ql*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(gI1,j1)) + 0.5477225575051661*g1*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(gI1,j1)) +
      0.7071067811865475*g2*KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(
      gI1,j1)) - KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(ZEL(gI2,j2))*SUM(j1,0,2,Ye(
      j1,j2)*ZE(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFeconjSePR(int gO1, int gI2, int gI1) const
{
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = -1.4142135623730951*gp*Qe*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZE(gI1,3 + j1)*ZER(gI2,j1)) - 1.0954451150103321*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,ZE(gI1,3 + j1)*ZER(gI2,j1)) -
      KroneckerDelta(3,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZE(
      gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFuconjSuPL(int gO2, int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = -1.4142135623730951*gp*Qq*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(gI1,j1)) - 0.18257418583505536*g1*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(gI1,j1)) -
      0.7071067811865475*g2*KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(
      gI1,j1)) - KroneckerDelta(4,gO2)*SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Yu(
      j1,j2)*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFuconjSuPR(int gO1, int gI2, int gI1) const
{
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = -1.4142135623730951*gp*Qu*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZU(gI1,3 + j1)*ZUR(gI2,j1)) + 0.7302967433402214*g1*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,ZU(gI1,3 + j1)*ZUR(gI2,j1)) -
      KroneckerDelta(4,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZU(
      gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFvconjSvPL(int gO2, int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = -1.4142135623730951*gp*Ql*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZVL(gI2,j1))*ZV(gI1,j1)) + 0.5477225575051661*g1*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZVL(gI2,j1))*ZV(gI1,j1)) -
      0.7071067811865475*g2*KroneckerDelta(2,gO2)*SUM(j1,0,2,Conj(ZVL(gI2,j1))*ZV(
      gI1,j1)) - KroneckerDelta(4,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZVL(gI2,j1))*Yv(
      j1,j2))*ZV(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFvconjSvPR(int gO1, int gI2, int gI1) const
{
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = -1.4142135623730951*gp*Qv*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZV(gI1,3 + j1)*ZVR(gI2,j1)) - KroneckerDelta(4,gO1)*SUM(j2,0
      ,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZV(gI1,j1))*ZVR(gI2,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(2,gO2)*UP(gI2,0)) +
      0.7071067811865475*g2*KroneckerDelta(4,gO2)*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(2
      ,gO1) + 1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(3,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZPL(int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(-(KroneckerDelta(3,gO2)*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 10*gp*QHd*Sin(ThetaWp()))*ZN(gI2,3)) + KroneckerDelta(4,gO2)*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      - 10*gp*QHu*Sin(ThetaWp()))*ZN(gI2,4) - 10*gp*Qs*KroneckerDelta(5,gO2)*Sin(
      ThetaWp())*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZPR(int gI2, int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(10*gp*Qs*Conj(ZN(gI2,5))*
      KroneckerDelta(5,gO1)*Sin(ThetaWp()) + Conj(ZN(gI2,3))*KroneckerDelta(3,gO1)
      *(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) - Conj(ZN(gI2,4))*KroneckerDelta(4
      ,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp(
      ))*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZpPL(int gI2, int gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(KroneckerDelta(3,gO2)*(-10*gp*QHd*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*ZN(gI2,3) - KroneckerDelta(4,gO2)*(10*gp*QHu*Cos(ThetaWp()) + 5*
      g2*Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(
      ThetaWp()))*ZN(gI2,4) - 10*gp*Qs*Cos(ThetaWp())*KroneckerDelta(5,gO2)*ZN(gI2
      ,5));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZpPR(int gI2, int gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   const std::complex<double> result = 0.1*(10*gp*Qs*Conj(ZN(gI2,5))*Cos(ThetaWp()
      )*KroneckerDelta(5,gO1) + Conj(ZN(gI2,4))*KroneckerDelta(4,gO1)*(10*gp*QHu*
      Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*
      Sin(ThetaW())*Sin(ThetaWp())) + Conj(ZN(gI2,3))*KroneckerDelta(3,gO1)*(10*gp
      *QHd*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW(
      )))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFeconjHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(ZEL(gI2,j1))*Yv(
      j1,gO2))*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFeconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Ye(j1,gO1))*ZER(
      gI2,j1))*ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarUFvSePL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UP(gI1,1))*SUM(j1,0,2,Conj(
      ZE(gI2,j1))*Yv(j1,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarUFvSePR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZE(gI2,gO1))*UM(gI1,0)
      ),0) + IF(gO1 < 3,SUM(j1,0,2,Conj(Ye(j1,gO1))*Conj(ZE(gI2,3 + j1)))*UM(gI1,1
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j1,0,2,Conj(ZVL(gI1,j1))*Yv(j1,gO2))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j2,0,2,Conj(Yv(gO1,j2))*ZVR(gI1,j2))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI1,
      1))*SUM(j1,0,2,Conj(ZVL(gI2,j1))*Yv(j1,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(ZH(gI1,
      1))*SUM(j2,0,2,Conj(Yv(gO1,j2))*ZVR(gI2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvChiSvPL(int gO2, int gI2, int gI1) const
{
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = IF(gO2 < 3,-1.4142135623730951*gp*Qv*Conj(
      ZN(gI2,0))*Conj(ZV(gI1,3 + gO2)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,4))*SUM(j1,0,
      2,Conj(ZV(gI1,j1))*Yv(j1,gO2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvChiSvPR(int gO1, int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Ql*Conj(
      ZV(gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,0.5477225575051661*g1*Conj(ZV(gI1,gO1
      ))*ZN(gI2,1),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZV(gI1,gO1))*ZN(gI2
      ,2),0) + IF(gO1 < 3,-(SUM(j2,0,2,Conj(Yv(gO1,j2))*Conj(ZV(gI1,3 + j2)))*ZN(
      gI2,4)),0);

   return result;
}

double CLASSNAME::CpbarUFvFeconjVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFeconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZEL(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvVZPR(int gO2, int gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = IF(gI2 < 3,gp*Qv*Sin(ThetaWp())*ZVR(gI2,gO2
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvVZPL(int gO1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(ZVL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(ZVL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-(gp*Ql*Conj(ZVL(gI2,gO1
      ))*Sin(ThetaWp())),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvVZpPR(int gO2, int gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   const std::complex<double> result = IF(gI2 < 3,gp*Qv*Cos(ThetaWp())*ZVR(gI2,gO2
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvVZpPL(int gO1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gI2 < 3,-(gp*Ql*Conj(ZVL(gI2,gO1))*Cos(
      ThetaWp())),0) + IF(gI2 < 3,0.5*g2*Conj(ZVL(gI2,gO1))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gI2 < 3,0.3872983346207417*g1*Conj(ZVL(gI2,gO1))*Sin(
      ThetaW())*Sin(ThetaWp()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *(g2*Conj(UM(gI1,0))*Conj(ZA(gI2,1))*KroneckerDelta(1,gO2) + Conj(UM(gI1,1))
      *(g2*Conj(ZA(gI2,0))*KroneckerDelta(0,gO2) - Conj(ZA(gI2,2))*KroneckerDelta(
      1,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*(g2*Conj(ZA(gI2,0))*KroneckerDelta(1,gO1)*UP(gI1,0) + (g2*Conj(ZA(gI2,1))*
      KroneckerDelta(0,gO1) - Conj(Lambdax)*Conj(ZA(gI2,2))*KroneckerDelta(1,gO1))
      *UP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiHpmPL(int gO2, int gI2, int gI1) const
{
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -(Conj(ZN(gI2,5))*KroneckerDelta(1,gO2)*
      Lambdax*ZP(gI1,0)) - 0.1*(10*g2*Conj(ZN(gI2,4))*KroneckerDelta(0,gO2) +
      1.4142135623730951*(10*gp*QHu*Conj(ZN(gI2,0)) + 3.872983346207417*g1*Conj(ZN
      (gI2,1)) + 5*g2*Conj(ZN(gI2,2)))*KroneckerDelta(1,gO2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiHpmPR(int gO1, int gI2, int gI1) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*ZN(gI2,3)*ZP(gI1
      ,0)) + KroneckerDelta(1,gO1)*(-1.4142135623730951*gp*QHd*ZN(gI2,0)*ZP(gI1,0)
      + 0.5477225575051661*g1*ZN(gI2,1)*ZP(gI1,0) + 0.7071067811865475*g2*ZN(gI2,2
      )*ZP(gI1,0) - Conj(Lambdax)*ZN(gI2,5)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChahhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*
      Conj(ZH(gI1,1))*KroneckerDelta(1,gO2) + Conj(UM(gI2,1))*(g2*Conj(ZH(gI1,0))*
      KroneckerDelta(0,gO2) + Conj(ZH(gI1,2))*KroneckerDelta(1,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChahhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*(g2*Conj(ZH(gI1,0))*
      KroneckerDelta(1,gO1)*UP(gI2,0) + (g2*Conj(ZH(gI1,1))*KroneckerDelta(0,gO1)
      + Conj(Lambdax)*Conj(ZH(gI1,2))*KroneckerDelta(1,gO1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZD(
      gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*SUM(j1,0,2,Conj(
      ZD(gI2,j1))*ZUL(gI1,j1))) + KroneckerDelta(1,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj
      (Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZVR(
      gI1,j2))*SUM(j1,0,2,Conj(ZE(gI2,j1))*Yv(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*SUM(j1,0,2,Conj(
      ZE(gI2,j1))*ZVL(gI1,j1))) + KroneckerDelta(1,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj
      (Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZVL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFdconjSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(
      ZDL(gI2,j1))*ZU(gI1,j1))) + KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZDL(gI2,j2
      ))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFdconjSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZU(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFeconjSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(
      ZEL(gI2,j1))*ZV(gI1,j1))) + KroneckerDelta(1,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj
      (ZEL(gI2,j1))*Yv(j1,j2))*ZV(gI1,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaFeconjSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZV(gI1,j2));

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
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = g2*Cos(ThetaW())*Cos(ThetaWp())*
      KroneckerDelta(0,gO2)*UP(gI2,0) + 0.1*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW
      ())*Cos(ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*
      gp*QHu*Sin(ThetaWp()))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZPL(int gO1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*Cos(
      ThetaWp())*KroneckerDelta(0,gO1) + 0.1*Conj(UM(gI2,1))*KroneckerDelta(1,gO1)
      *(5*g2*Cos(ThetaW())*Cos(ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW()) - 10*gp*QHd*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZpPR(int gO2, int gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   const std::complex<double> result = -(g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*
      Sin(ThetaWp())*UP(gI2,0)) + 0.1*KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp
      ()) + (-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp
      ()))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZpPL(int gO1, int gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*Cos(ThetaW())*
      KroneckerDelta(0,gO1)*Sin(ThetaWp())) - 0.1*Conj(UM(gI2,1))*KroneckerDelta(1
      ,gO1)*(10*gp*QHd*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1
      *Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,2)) +
      0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZN(gI2,4);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(ZN(gI2,2))*KroneckerDelta(0
      ,gO1) + 1.4142135623730951*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(ZVL(gI2,j2))*Ye(
      gO2,j2))*ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,SUM(j2,0,2,Conj(Yv(gO1,j2))*ZVR(
      gI2,j2))*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j2,0,2,Conj(ZEL(gI1,j2))*Ye(gO2,j2))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j1,0,2,Conj(Ye(j1,gO1))*ZER(gI1,j1))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI1,
      0))*SUM(j2,0,2,Conj(ZEL(gI2,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(ZH(gI1,
      0))*SUM(j1,0,2,Conj(Ye(j1,gO1))*ZER(gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChaSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UM(gI2,1))*SUM(j2,0,2,Conj(
      ZV(gI1,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChaSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZV(gI1,gO1))*UP(gI2,0)
      ),0) + IF(gO1 < 3,SUM(j2,0,2,Conj(Yv(gO1,j2))*Conj(ZV(gI1,3 + j2)))*UP(gI2,1
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePL(int gO2, int gI2, int gI1) const
{
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = IF(gO2 < 3,-1.4142135623730951*gp*Qe*Conj(
      ZE(gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,-1.0954451150103321*g1*Conj
      (ZE(gI1,3 + gO2))*Conj(ZN(gI2,1)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,3))*SUM(j2,0
      ,2,Conj(ZE(gI1,j2))*Ye(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePR(int gO1, int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Ql*Conj(
      ZE(gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,0.5477225575051661*g1*Conj(ZE(gI1,gO1
      ))*ZN(gI2,1),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZE(gI1,gO1))*ZN(gI2,
      2),0) + IF(gO1 < 3,-(SUM(j1,0,2,Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1)))*ZN(
      gI2,3)),0);

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
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZER(gI2,gO2),0) + IF(gI2 < 3,gp*Qe*Sin(ThetaWp())*
      ZER(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPL(int gO1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZEL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(ZEL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-(gp*Ql*Conj(ZEL(gI2,gO1
      ))*Sin(ThetaWp())),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZpPR(int gO2, int gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = IF(gI2 < 3,gp*Qe*Cos(ThetaWp())*ZER(gI2,gO2
      ),0) + IF(gI2 < 3,0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())*ZER(gI2
      ,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZpPL(int gO1, int gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = IF(gI2 < 3,-(gp*Ql*Conj(ZEL(gI2,gO1))*Cos(
      ThetaWp())),0) + IF(gI2 < 3,-0.5*g2*Conj(ZEL(gI2,gO1))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gI2 < 3,0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Sin(
      ThetaW())*Sin(ThetaWp()),0);

   return result;
}

double CLASSNAME::CpbarUFeFvVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZVL(
      gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(ZUL(gI2,j2))*Yd(
      gO2,j2))*ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yu(j1,gO1))*ZUR(
      gI2,j1))*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j2,0,2,Conj(ZDL(gI1,j2))*Yd(gO2,j2))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZA(gI2,0))*SUM(j1,0,2,Conj(Yd(j1,gO1))*ZDR(gI1,j1))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI1,
      0))*SUM(j2,0,2,Conj(ZDL(gI2,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(ZH(gI1,
      0))*SUM(j1,0,2,Conj(Yd(j1,gO1))*ZDR(gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChaSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UM(gI2,1))*SUM(j2,0,2,Conj(
      ZU(gI1,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChaSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZU(gI1,gO1))*UP(gI2,0)
      ),0) + IF(gO1 < 3,SUM(j1,0,2,Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1)))*UP(gI2,1
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPL(int gO2, int gI2, int gI1) const
{
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = IF(gO2 < 3,-1.4142135623730951*gp*Qd*Conj(
      ZD(gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,-0.3651483716701107*g1*Conj
      (ZD(gI1,3 + gO2))*Conj(ZN(gI2,1)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,3))*SUM(j2,0
      ,2,Conj(ZD(gI1,j2))*Yd(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPR(int gO1, int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Qq*Conj(
      ZD(gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZD(gI1,
      gO1))*ZN(gI2,1),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZD(gI1,gO1))*ZN(
      gI2,2),0) + IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1)))*
      ZN(gI2,3)),0);

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
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZDR(gI2,gO2),0) + IF(gI2 < 3,gp*Qd*Sin(ThetaWp())*
      ZDR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPL(int gO1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZDL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZDL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-(gp*Qq*Conj(ZDL(gI2,gO1
      ))*Sin(ThetaWp())),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZpPR(int gO2, int gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = IF(gI2 < 3,gp*Qd*Cos(ThetaWp())*ZDR(gI2,gO2
      ),0) + IF(gI2 < 3,0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp())*ZDR(gI2
      ,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZpPL(int gO1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gI2 < 3,-(gp*Qq*Conj(ZDL(gI2,gO1))*Cos(
      ThetaWp())),0) + IF(gI2 < 3,-0.5*g2*Conj(ZDL(gI2,gO1))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gI2 < 3,-0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Sin(
      ThetaW())*Sin(ThetaWp()),0);

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
   
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(ZDL(gI2,j2))*Yu(
      gO2,j2))*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yd(j1,gO1))*ZDR(
      gI2,j1))*ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarUFuSdPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Conj(UP(gI1,1))*SUM(j2,0,2,Conj(
      ZD(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarUFuSdPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZD(gI2,gO1))*UM(gI1,0)
      ),0) + IF(gO1 < 3,SUM(j1,0,2,Conj(Yd(j1,gO1))*Conj(ZD(gI2,3 + j1)))*UM(gI1,1
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j2,0,2,Conj(ZUL(gI1,j2))*Yu(gO2,j2))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*Conj(ZA(gI2,1))*SUM(j1,0,2,Conj(Yu(j1,gO1))*ZUR(gI1,j1))
      ,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*Conj(ZH(gI1,
      1))*SUM(j2,0,2,Conj(ZUL(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*Conj(ZH(gI1,
      1))*SUM(j1,0,2,Conj(Yu(j1,gO1))*ZUR(gI2,j1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPL(int gO2, int gI2, int gI1) const
{
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = IF(gO2 < 3,-1.4142135623730951*gp*Qu*Conj(
      ZN(gI2,0))*Conj(ZU(gI1,3 + gO2)),0) + IF(gO2 < 3,0.7302967433402214*g1*Conj(
      ZN(gI2,1))*Conj(ZU(gI1,3 + gO2)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,4))*SUM(j2,0,
      2,Conj(ZU(gI1,j2))*Yu(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPR(int gO1, int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*gp*Qq*Conj(
      ZU(gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZU(gI1,
      gO1))*ZN(gI2,1),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZU(gI1,gO1))*ZN(
      gI2,2),0) + IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1)))*
      ZN(gI2,4)),0);

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
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Cos(
      ThetaWp())*Sin(ThetaW())*ZUR(gI2,gO2),0) + IF(gI2 < 3,gp*Qu*Sin(ThetaWp())*
      ZUR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(ZUL(gI2,gO1))*Cos(
      ThetaW())*Cos(ThetaWp()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZUL(gI2
      ,gO1))*Cos(ThetaWp())*Sin(ThetaW()),0) + IF(gI2 < 3,-(gp*Qq*Conj(ZUL(gI2,gO1
      ))*Sin(ThetaWp())),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZpPR(int gO2, int gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = IF(gI2 < 3,gp*Qu*Cos(ThetaWp())*ZUR(gI2,gO2
      ),0) + IF(gI2 < 3,-0.5163977794943222*g1*Sin(ThetaW())*Sin(ThetaWp())*ZUR(
      gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZpPL(int gO1, int gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = IF(gI2 < 3,-(gp*Qq*Conj(ZUL(gI2,gO1))*Cos(
      ThetaWp())),0) + IF(gI2 < 3,0.5*g2*Conj(ZUL(gI2,gO1))*Cos(ThetaW())*Sin(
      ThetaWp()),0) + IF(gI2 < 3,-0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Sin(
      ThetaW())*Sin(ThetaWp()),0);

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

std::complex<double> CLASSNAME::CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j2,0,2,Conj(ZVL(gI2,j2))*SUM(j1,0,2,
      Conj(ZER(gO2,j1))*Ye(j1,j2)))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yv(j1,j2))*ZEL(
      gO1,j1))*ZVR(gI2,j2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*Conj(ZA(gI2,0))*SUM(j2,0,2,Conj(ZEL(gI1,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*
      Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *Conj(ZA(gI2,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI1,j1))*ZEL(gO1,
      j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*Conj(ZH(gI1,0))*SUM(j2,
      0,2,Conj(ZEL(gI2,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*Conj(ZH(gI1,0))*SUM(j2,
      0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZEL(gO1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChaSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = Conj(UM(gI2,1))*SUM(j2,0,2,Conj(ZV(gI1,j2))
      *SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChaSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZV(gI1,j1))*ZEL(gO1,j1
      ))*UP(gI2,0)) + SUM(j2,0,2,Conj(ZV(gI1,3 + j2))*SUM(j1,0,2,Conj(Yv(j1,j2))*
      ZEL(gO1,j1)))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePL(int gO2, int gI2, int gI1) const
{
   const auto Qe = LOCALINPUT(Qe);

   const std::complex<double> result = -1.4142135623730951*gp*Qe*Conj(ZN(gI2,0))*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1))) - 1.0954451150103321*g1*
      Conj(ZN(gI2,1))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1))) - Conj(ZN
      (gI2,3))*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2))
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePR(int gO1, int gI2, int gI1) const
{
   const auto Ql = LOCALINPUT(Ql);

   const std::complex<double> result = 0.1414213562373095*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*ZEL(gO1,j1))*(-10*gp*Ql*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) + 5*
      g2*ZN(gI2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*
      ZEL(gO1,j2))*ZN(gI2,3);

   return result;
}

double CLASSNAME::CpbarFeFvVWmPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(ZVL(
      gI2,j1))*ZEL(gO1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,
      Conj(ZDR(gO2,j1))*Yd(j1,j2)))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(
      gI2,j1))*ZDL(gO1,j2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*Conj(ZA(gI2,0))*SUM(j2,0,2,Conj(ZDL(gI1,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*
      Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *Conj(ZA(gI2,0))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI1,j1))*ZDL(gO1,
      j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*Conj(ZH(gI1,0))*SUM(j2,
      0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*Conj(ZH(gI1,0))*SUM(j2,
      0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZDL(gO1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChaSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = Conj(UM(gI2,1))*SUM(j2,0,2,Conj(ZU(gI1,j2))
      *SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChaSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZDL(gO1,j1
      ))*UP(gI2,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*
      ZDL(gO1,j2))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPL(int gO2, int gI2, int gI1) const
{
   const auto Qd = LOCALINPUT(Qd);

   const std::complex<double> result = -1.4142135623730951*gp*Qd*Conj(ZN(gI2,0))*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1))) - 0.3651483716701107*g1*
      Conj(ZN(gI2,1))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1))) - Conj(ZN
      (gI2,3))*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2))
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPR(int gO1, int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = -0.2357022603955158*SUM(j1,0,2,Conj(ZD(gI1,
      j1))*ZDL(gO1,j1))*(6*gp*Qq*ZN(gI2,0) + 0.7745966692414834*g1*ZN(gI2,1) - 3*
      g2*ZN(gI2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*
      ZDL(gO1,j2))*ZN(gI2,3);

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
   
   const std::complex<double> result = SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,
      Conj(ZUR(gO2,j1))*Yu(j1,j2)))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(
      gI2,j1))*ZUL(gO1,j2))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarFuSdPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = Conj(UP(gI1,1))*SUM(j2,0,2,Conj(ZD(gI2,j2))
      *SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChabarFuSdPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZUL(gO1,j1
      ))*UM(gI1,0)) + SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*
      ZUL(gO1,j2))*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*Conj(ZA(gI2,1))*SUM(j2,0,2,Conj(ZUL(gI1,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*
      Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *Conj(ZA(gI2,1))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI1,j1))*ZUL(gO1,
      j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*Conj(ZH(gI1,1))*SUM(j2,
      0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*Conj(ZH(gI1,1))*SUM(j2,
      0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZUL(gO1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPL(int gO2, int gI2, int gI1) const
{
   const auto Qu = LOCALINPUT(Qu);

   const std::complex<double> result = -1.4142135623730951*gp*Qu*Conj(ZN(gI2,0))*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1))) + 0.7302967433402214*g1*
      Conj(ZN(gI2,1))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1))) - Conj(ZN
      (gI2,4))*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2))
      );

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPR(int gO1, int gI2, int gI1) const
{
   const auto Qq = LOCALINPUT(Qq);

   const std::complex<double> result = -0.2357022603955158*SUM(j1,0,2,Conj(ZU(gI1,
      j1))*ZUL(gO1,j1))*(6*gp*Qq*ZN(gI2,0) + 0.7745966692414834*g1*ZN(gI2,1) + 3*
      g2*ZN(gI2,2)) - SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*
      ZUL(gO1,j2))*ZN(gI2,4);

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
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSdconjUSd(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSdconjUSd(gI1,gI1,gO1,gO2))
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
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUSdSeconjUSdconjSe(gO1,gI1,gO2,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUSdSvconjUSdconjSv(gO1,gI1,gO2,gI1))
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
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*
      CpChiFvconjUSvPL(gI2,gI1,gO1) + Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*
      CpChiFvconjUSvPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*CpChiFvconjUSvPL(gI2,gI1,gO1) +
      Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*CpChiFvconjUSvPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSvconjSdconjUSv(gI1,gO1,gI1,gO2
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSvconjSeconjUSv(gI1,gO1,gI1,gO2))
      ;
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuUSvconjSuconjUSv(gI1,gO1,gI1,gO2
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpSvUSvconjSvconjUSv(gI1,gO1,gI1,gO2))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSvconjUSv(gI2,gI1,gO2))*CpAhSvconjUSv(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSvconjUSv(gI2,gI1,gO2))*CphhSvconjUSv(gI2,gI1,gO1)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSvconjVWm(gI2,gO2))*CpSeconjUSvconjVWm(gI2,
      gO1)*F0(Sqr(p),Sqr(MSe(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,5,Conj(CpSvconjUSvVZ(gI2,gO2))*CpSvconjUSvVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSvconjUSvVZp(gI2,gO2))*CpSvconjUSvVZp(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVZp)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,6,6> CLASSNAME::self_energy_Sv_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,6,6> self_energy;

   for (int i = 0; i < 6; i++)
      for (int k = i; k < 6; k++)
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
   result += SUM(gI1,0,2,SUM(gI2,0,5,(Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*
      CpChiFuconjUSuPL(gI2,gI1,gO1) + Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*
      CpChiFuconjUSuPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*CpChiFuconjUSuPL(gI2,gI1,gO1) +
      Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*CpChiFuconjUSuPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSuconjSeconjUSu(gI1,gO1,gI1,gO2))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUSuSvconjUSuconjSv(gO1,gI1,gO2,gI1))
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
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSeconjUSe(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSeconjUSe(gI1,gI1,gO1,gO2))
      ;
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
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSeconjSdconjUSe(gI1,gO1,gI1,gO2
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSeconjSeconjUSe(gI1,gO1,gI1,gO2))
      ;
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSeSuconjUSeconjSu(gO1,gI1,gO2,gI1
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUSeSvconjUSeconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MHpm(gI2)))*Conj(
      CpHpmSvconjUSe(gI2,gI1,gO2))*CpHpmSvconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSeconjUSe(gI2,gI1,gO2))*CpAhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSeconjUSe(gI2,gI1,gO2))*CphhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVP(gI2,gO2))*CpSeconjUSeVP(gI2,gO1)*F0(
      Sqr(p),Sqr(MSe(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVZ(gI2,gO2))*CpSeconjUSeVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSe(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVZp(gI2,gO2))*CpSeconjUSeVZp(gI2,gO1)*F0(
      Sqr(p),Sqr(MSe(gI2)),Sqr(MVZp)));
   result += SUM(gI2,0,5,Conj(CpSvconjUSeVWm(gI2,gO2))*CpSvconjUSeVWm(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVWm)));

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
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*Conj
      (CpUhhHpmconjHpm(gO2,gI2,gI1))*CpUhhHpmconjHpm(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaChaUhhPL(gI1,gI2,gO2))*
      CpbarChaChaUhhPL(gI1,gI2,gO1) + Conj(CpbarChaChaUhhPR(gI1,gI2,gO2))*
      CpbarChaChaUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpbarChaChaUhhPR(gI1,gI2,gO2))*CpbarChaChaUhhPL(gI1,gI2,
      gO1) + Conj(CpbarChaChaUhhPL(gI1,gI2,gO2))*CpbarChaChaUhhPR(gI1,gI2,gO1))*
      MCha(gI2)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhhUhh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhhUhh(gI1,gI1,gO1,gO2));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUhh(gI1,gI2,gO2))*CpAhAhUhh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhhhUhh(gI2,gI1,gO2))*CpAhhhUhh(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhhhUhh(gI1,gI2,gO2))*CphhhhUhh(gI1,gI2,gO1)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUhhPL(gI1,gI2,gO2))*
      CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*
      CpbarFdFdUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUhhPL(gI1,gI2,gO2))*
      CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*
      CpbarFeFeUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUhhPL(gI1,gI2,gO2))*
      CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*
      CpbarFuFuUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFvFvUhhPL(gI1,gI2,gO2))*
      CpbarFvFvUhhPL(gI1,gI2,gO1) + Conj(CpbarFvFvUhhPR(gI1,gI2,gO2))*
      CpbarFvFvUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUhhPL(gI1,gI2,gO2))*CpbarFdFdUhhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUhhPL(gI1,gI2,gO2))*CpbarFeFeUhhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUhhPL(gI1,gI2,gO2))*CpbarFuFuUhhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(
      gI2)))*(Conj(CpbarFvFvUhhPR(gI1,gI2,gO2))*CpbarFvFvUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFvFvUhhPL(gI1,gI2,gO2))*CpbarFvFvUhhPR(gI1,gI2,gO1))*MFv(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhUhhSdconjSd(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhUhhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhUhhSuconjSu(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUhhUhhSvconjSv(gO1,gO2,gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*Conj
      (CpUhhSdconjSd(gO2,gI2,gI1))*CpUhhSdconjSd(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*Conj(
      CpUhhSeconjSe(gO2,gI2,gI1))*CpUhhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*Conj
      (CpUhhSuconjSu(gO2,gI2,gI1))*CpUhhSuconjSu(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSv(gI2)))*Conj(
      CpUhhSvconjSv(gO2,gI2,gI1))*CpUhhSvconjSv(gO1,gI2,gI1)));
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
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*Conj
      (CpUAhHpmconjHpm(gO2,gI2,gI1))*CpUAhHpmconjHpm(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarChaChaUAhPL(gI1,gI2,gO2))*
      CpbarChaChaUAhPL(gI1,gI2,gO1) + Conj(CpbarChaChaUAhPR(gI1,gI2,gO2))*
      CpbarChaChaUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpbarChaChaUAhPR(gI1,gI2,gO2))*CpbarChaChaUAhPL(gI1,gI2,
      gO1) + Conj(CpbarChaChaUAhPL(gI1,gI2,gO2))*CpbarChaChaUAhPR(gI1,gI2,gO1))*
      MCha(gI2)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUAhUAh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUAh(gI1,gI2,gO2))*CpAhAhUAh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhUAhhh(gI2,gO2,gI1))*CpAhUAhhh(gI2,gO1,gI1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CpUAhhhhh(gO2,gI1,gI2))*CpUAhhhhh(gO1,gI1,gI2)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUAhPL(gI1,gI2,gO2))*
      CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*
      CpbarFdFdUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUAhPL(gI1,gI2,gO2))*
      CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*
      CpbarFeFeUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUAhPL(gI1,gI2,gO2))*
      CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*
      CpbarFuFuUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFvFvUAhPL(gI1,gI2,gO2))*
      CpbarFvFvUAhPL(gI1,gI2,gO1) + Conj(CpbarFvFvUAhPR(gI1,gI2,gO2))*
      CpbarFvFvUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUAhPL(gI1,gI2,gO2))*CpbarFdFdUAhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUAhPL(gI1,gI2,gO2))*CpbarFeFeUAhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUAhPL(gI1,gI2,gO2))*CpbarFuFuUAhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(
      gI2)))*(Conj(CpbarFvFvUAhPR(gI1,gI2,gO2))*CpbarFvFvUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFvFvUAhPL(gI1,gI2,gO2))*CpbarFvFvUAhPR(gI1,gI2,gO1))*MFv(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUAhUAhSdconjSd(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUAhUAhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUAhUAhSuconjSu(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUAhUAhSvconjSv(gO1,gO2,gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*Conj
      (CpUAhSdconjSd(gO2,gI2,gI1))*CpUAhSdconjSd(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*Conj(
      CpUAhSeconjSe(gO2,gI2,gI1))*CpUAhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*Conj
      (CpUAhSuconjSu(gO2,gI2,gI1))*CpUAhSuconjSu(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSv(gI2)))*Conj(
      CpUAhSvconjSv(gO2,gI2,gI1))*CpUAhSvconjSv(gO1,gI2,gI1)));
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
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhHpmconjUHpm(gI2,gI1,gO2))*CpAhHpmconjUHpm(gI2,gI1,gO1)));
   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhHpmconjUHpm(gI2,gI1,gO2))*CphhHpmconjUHpm(gI2,gI1,gO1)));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUHpmconjUHpm(gI1,gI1,gO1,gO2
      ));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUHpmconjUHpm(gI1,gI1,gO1,gO2
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
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUHpmSdconjUHpmconjSd(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUHpmSeconjUHpmconjSe(gO1,gI1,gO2,gI1
      ));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUHpmSuconjUHpmconjSu(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUHpmSvconjUHpmconjSv(gO1,gI1,gO2,gI1
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,1,(Conj(CpChiChaconjUHpmPL(gI1,gI2,gO2))*
      CpChiChaconjUHpmPL(gI1,gI2,gO1) + Conj(CpChiChaconjUHpmPR(gI1,gI2,gO2))*
      CpChiChaconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpChiChaconjUHpmPR(gI1,gI2,gO2))*CpChiChaconjUHpmPL(gI1,
      gI2,gO1) + Conj(CpChiChaconjUHpmPL(gI1,gI2,gO2))*CpChiChaconjUHpmPR(gI1,gI2,
      gO1))*MCha(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSd(gI2)))*Conj
      (CpSdconjUHpmconjSu(gI2,gO2,gI1))*CpSdconjUHpmconjSu(gI2,gO1,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSe(gI2)))*Conj(
      CpSeconjUHpmconjSv(gI2,gO2,gI1))*CpSeconjUHpmconjSv(gI2,gO1,gI1)));
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
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVGPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVGPL(gI1,
      gI2))*CpbarFuFuVGPR(gI1,gI2))));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVGVG(gI1,gI1));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVGVG(gI1,gI1));
   result += -2*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVG(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
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
   result += -2*AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm),
      Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVPVP(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHpmconjHpmVP(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaChaVPPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVPPL(gI1,gI2))*CpbarChaChaVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVPPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVPPL(gI1,
      gI2))*CpbarFdFdVPPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVPPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVPPL(gI1,
      gI2))*CpbarFeFeVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVPPL(gI1,gI2)) + AbsSqr(
      CpbarFuFuVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVPPL(gI1,
      gI2))*CpbarFuFuVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVPVP(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVPVP(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVPVP(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
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
   result += -2*AbsSqr(CpconjVWmVWmVZ())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm),
      Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHpmconjHpmVZ(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaChaVZPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVZPL(gI1,gI2))*CpbarChaChaVZPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhVZVZ(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpAhhhVZ(gI2,gI1))*B00(Sqr(p),Sqr(
      MAh(gI2)),Sqr(Mhh(gI1)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZPL(gI1,
      gI2))*CpbarFdFdVZPR(gI1,gI2))));
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
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVZVZ(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVZVZ(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVZVZ(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpSvconjSvVZVZ(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSvconjSvVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSv(gI1)),Sqr(MSv(gI2)))));
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
   result += -2*AbsSqr(CpconjVWmVWmVZp())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm),
      Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVZpVZp(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpHpmconjHpmVZp(gI2,gI1))*B00(Sqr(p
      ),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarChaChaVZpPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVZpPL(gI1,gI2))*CpbarChaChaVZpPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhVZpVZp(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhVZpVZp(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpAhhhVZp(gI2,gI1))*B00(Sqr(p),Sqr(
      MAh(gI2)),Sqr(Mhh(gI1)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZpPL(gI1,gI2)) + AbsSqr(
      CpbarFdFdVZpPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(
      p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZpPL(gI1
      ,gI2))*CpbarFdFdVZpPR(gI1,gI2))));
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
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVZpVZp(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVZpVZp(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVZpVZp(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpSvconjSvVZpVZp(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVZp(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVZp(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVZp(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSvconjSvVZp(gI2,gI1))*B00(Sqr(p),
      Sqr(MSv(gI1)),Sqr(MSv(gI2)))));
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
   result += -4*SUM(gI1,0,1,SUM(gI2,0,2,AbsSqr(CpAhHpmconjVWm(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MAh(gI2)),Sqr(MHpm(gI1)))));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,2,AbsSqr(CphhHpmconjVWm(gI2,gI1))*B00(Sqr(p)
      ,Sqr(Mhh(gI2)),Sqr(MHpm(gI1)))));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhconjVWmVWm(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFdconjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarFuFdconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))
      + 4*B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))*MFd(gI2)*MFu(gI1)*Re(Conj(
      CpbarFuFdconjVWmPL(gI1,gI2))*CpbarFuFdconjVWmPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFeconjVWmPL(gI1,gI2)) + AbsSqr
      (CpbarFvFeconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2))) + 4*B0
      (Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))*MFe(gI2)*MFv(gI1)*Re(Conj(
      CpbarFvFeconjVWmPL(gI1,gI2))*CpbarFvFeconjVWmPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpSvconjSvconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,(AbsSqr(CpChiChaconjVWmPL(gI1,gI2)) + AbsSqr(
      CpChiChaconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2))) + 4*B0
      (Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2)))*MCha(gI2)*MChi(gI1)*Re(Conj(
      CpChiChaconjVWmPL(gI1,gI2))*CpChiChaconjVWmPR(gI1,gI2))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSuconjVWm(gI2,gI1))*B00(
      Sqr(p),Sqr(MSd(gI2)),Sqr(MSu(gI1)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSvconjVWm(gI2,gI1))*B00(Sqr
      (p),Sqr(MSe(gI2)),Sqr(MSv(gI1)))));
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

   result += -4*SUM(gI1,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPL(gI1,gO2))*CpbarChaUChiVWmPR(gI1,gO1)*MCha(gI1));
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(
      gI2)))*Conj(CpbarChaUChiHpmPL(gI1,gO2,gI2))*CpbarChaUChiHpmPR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))*Conj
      (CpUChiChaconjHpmPL(gO2,gI2,gI1))*CpUChiChaconjHpmPR(gO1,gI2,gI1)*MCha(gI2))
      );
   result += 3*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarFdUChiSdPL(gI1,gO2,gI2))*CpbarFdUChiSdPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarFeUChiSePL(gI1,gO2,gI2))*CpbarFeUChiSePR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(
      gI2)))*Conj(CpbarFuUChiSuPL(gI1,gO2,gI2))*CpbarFuUChiSuPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)
      ))*Conj(CpbarFvUChiSvPL(gI1,gO2,gI2))*CpbarFvUChiSvPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpChiUChihhPL(gI2,gO2,gI1))*CpChiUChihhPR(gI2,gO1,gI1)*MChi(gI2)));
   result += SUM(gI1,0,5,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(
      gI2)))*Conj(CpChiUChiAhPL(gI1,gO2,gI2))*CpChiUChiAhPR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*Conj
      (CpUChiFdconjSdPL(gO2,gI2,gI1))*CpUChiFdconjSdPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*Conj(
      CpUChiFeconjSePL(gO2,gI2,gI1))*CpUChiFeconjSePR(gO1,gI2,gI1)*MFe(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*Conj
      (CpUChiFuconjSuPL(gO2,gI2,gI1))*CpUChiFuconjSuPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*Conj(
      CpUChiFvconjSvPL(gO2,gI2,gI1))*CpUChiFvconjSvPR(gO1,gI2,gI1)*MFv(gI2)));
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
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(gI2)))
      *Conj(CpbarChaUChiHpmPR(gI1,gO2,gI2))*CpbarChaUChiHpmPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpUChiChaconjHpmPR(gO2,gI2,gI1))*CpUChiChaconjHpmPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarFdUChiSdPR(gI1,gO2,gI2))*CpbarFdUChiSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarFeUChiSePR(gI1,gO2,gI2))*CpbarFeUChiSePR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpbarFuUChiSuPR(gI1,gO2,gI2))*CpbarFuUChiSuPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpbarFvUChiSvPR(gI1,gO2,gI2))*CpbarFvUChiSvPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpChiUChihhPR(gI2,gO2,gI1))*CpChiUChihhPR(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)))*
      Conj(CpChiUChiAhPR(gI1,gO2,gI2))*CpChiUChiAhPR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpUChiFdconjSdPR(gO2,gI2,gI1))*CpUChiFdconjSdPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpUChiFeconjSePR(gO2,gI2,gI1))*CpUChiFeconjSePR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpUChiFuconjSuPR(gO2,gI2,gI1))*CpUChiFuconjSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*
      Conj(CpUChiFvconjSvPR(gO2,gI2,gI1))*CpUChiFvconjSvPR(gO1,gI2,gI1)));
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
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(gI2)))
      *Conj(CpbarChaUChiHpmPL(gI1,gO2,gI2))*CpbarChaUChiHpmPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpUChiChaconjHpmPL(gO2,gI2,gI1))*CpUChiChaconjHpmPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarFdUChiSdPL(gI1,gO2,gI2))*CpbarFdUChiSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarFeUChiSePL(gI1,gO2,gI2))*CpbarFeUChiSePL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpbarFuUChiSuPL(gI1,gO2,gI2))*CpbarFuUChiSuPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpbarFvUChiSvPL(gI1,gO2,gI2))*CpbarFvUChiSvPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpChiUChihhPL(gI2,gO2,gI1))*CpChiUChihhPL(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)))*
      Conj(CpChiUChiAhPL(gI1,gO2,gI2))*CpChiUChiAhPL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpUChiFdconjSdPL(gO2,gI2,gI1))*CpUChiFdconjSdPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpUChiFeconjSePL(gO2,gI2,gI1))*CpUChiFeconjSePL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpUChiFuconjSuPL(gO2,gI2,gI1))*CpUChiFuconjSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*
      Conj(CpUChiFvconjSvPL(gO2,gI2,gI1))*CpUChiFvconjSvPL(gO1,gI2,gI1)));
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

std::complex<double> CLASSNAME::self_energy_Fv_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFvFeconjHpmPL(gO2,gI2,gI1))*CpbarUFvFeconjHpmPR(gO1,gI2,gI1)*MFe(gI2))
      );
   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(
      gI2)))*Conj(CpbarChabarUFvSePL(gI1,gO2,gI2))*CpbarChabarUFvSePR(gI1,gO1,gI2)
      ));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFvFvAhPL(gO2,gI1,gI2))*CpbarUFvFvAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFvFvhhPL(gO2,gI2,gI1))*CpbarUFvFvhhPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUFvChiSvPL(gO2,gI2,gI1))*CpbarUFvChiSvPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarUFvFeconjVWmPR(gO2,gI2))*CpbarUFvFeconjVWmPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarUFvFvVZPR(
      gO2,gI2))*CpbarUFvFvVZPL(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZp))*Conj(
      CpbarUFvFvVZpPR(gO2,gI2))*CpbarUFvFvVZpPL(gO1,gI2)*MFv(gI2));

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
      Conj(CpbarUFvFeconjHpmPR(gO2,gI2,gI1))*CpbarUFvFeconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarChabarUFvSePR(gI1,gO2,gI2))*CpbarChabarUFvSePR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFvFvAhPR(gO2,gI1,gI2))*CpbarUFvFvAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFvFvhhPR(gO2,gI2,gI1))*CpbarUFvFvhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFvChiSvPR(gO2,gI2,gI1))*CpbarUFvChiSvPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarUFvFeconjVWmPL(gO2,gI2))*CpbarUFvFeconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarUFvFvVZPL(
      gO2,gI2))*CpbarUFvFvVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZp))*Conj(CpbarUFvFvVZpPL(
      gO2,gI2))*CpbarUFvFvVZpPL(gO1,gI2));

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
      Conj(CpbarUFvFeconjHpmPL(gO2,gI2,gI1))*CpbarUFvFeconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarChabarUFvSePL(gI1,gO2,gI2))*CpbarChabarUFvSePL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFvFvAhPL(gO2,gI1,gI2))*CpbarUFvFvAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFvFvhhPL(gO2,gI2,gI1))*CpbarUFvFvhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFvChiSvPL(gO2,gI2,gI1))*CpbarUFvChiSvPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarUFvFeconjVWmPR(gO2,gI2))*CpbarUFvFeconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZp))*Conj(CpbarUFvFvVZpPR(
      gO2,gI2))*CpbarUFvFvVZpPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarUFvFvVZPR(
      gO2,gI2))*CpbarUFvFvVZPR(gO1,gI2));

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

std::complex<double> CLASSNAME::self_energy_Cha_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUChaChaAhPL(gO2,gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))*Conj
      (CpbarUChaChiHpmPL(gO2,gI2,gI1))*CpbarUChaChiHpmPR(gO1,gI2,gI1)*MChi(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUChaChahhPL(gO2,gI2,gI1))*CpbarUChaChahhPR(gO1,gI2,gI1)*MCha(gI2)));
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarUChabarFuSdPL(gO2,gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2)
      ));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarUChabarFvSePL(gO2,gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*Conj
      (CpbarUChaFdconjSuPL(gO2,gI2,gI1))*CpbarUChaFdconjSuPR(gO1,gI2,gI1)*MFd(gI2)
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUChaFeconjSvPL(gO2,gI2,gI1))*CpbarUChaFeconjSvPR(gO1,gI2,gI1)*MFe(gI2))
      );
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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUChaChaAhPR(gO2,gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpbarUChaChiHpmPR(gO2,gI2,gI1))*CpbarUChaChiHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChaChahhPR(gO2,gI2,gI1))*CpbarUChaChahhPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarUChabarFuSdPR(gO2,gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarUChabarFvSePR(gO2,gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChaFdconjSuPR(gO2,gI2,gI1))*CpbarUChaFdconjSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUChaFeconjSvPR(gO2,gI2,gI1))*CpbarUChaFeconjSvPR(gO1,gI2,gI1)));
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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUChaChaAhPL(gO2,gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpbarUChaChiHpmPL(gO2,gI2,gI1))*CpbarUChaChiHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChaChahhPL(gO2,gI2,gI1))*CpbarUChaChahhPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarUChabarFuSdPL(gO2,gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarUChabarFvSePL(gO2,gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChaFdconjSuPL(gO2,gI2,gI1))*CpbarUChaFdconjSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUChaFeconjSvPL(gO2,gI2,gI1))*CpbarUChaFeconjSvPL(gO1,gI2,gI1)));
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
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUFeChaSvPL(gO2,gI2,gI1))*CpbarUFeChaSvPR(gO1,gI2,gI1)*MCha(gI2)));
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
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPR(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPR(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFeChaSvPR(gO2,gI2,gI1))*CpbarUFeChaSvPR(gO1,gI2,gI1)));
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
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFeChaSvPL(gO2,gI2,gI1))*CpbarUFeChaSvPL(gO1,gI2,gI1)));
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

std::complex<double> CLASSNAME::self_energy_Glu_1loop_1(double p ) const
{
   std::complex<double> result;

   result += -12*MGlu*B0(Sqr(p),Sqr(MGlu),0)*Conj(CpGluGluVGPR())*CpGluGluVGPL();
   result += 0.5*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarFdGluSdPL(gI1,gI2))*CpbarFdGluSdPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(
      gI2)))*Conj(CpbarFuGluSuPL(gI1,gI2))*CpbarFuGluSuPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpGluFdconjSdPL(gI2,gI1))*CpGluFdconjSdPR(gI2,gI1)*MFd(gI2)));
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
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFuGluSuPR(gI1,gI2))*B1(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFdconjSdPR(gI2,gI1))*B1(Sqr
      (p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))));
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
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarFuGluSuPL(gI1,gI2))*B1(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFdconjSdPL(gI2,gI1))*B1(Sqr
      (p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpGluFuconjSuPL(gI2,gI1))*B1(Sqr
      (p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarFeChaSvPL(gO2,gI2,gI1))*CpbarFeChaSvPR(gO1,gI2,gI1)*MCha(gI2)));
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
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPR(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPR(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFeChaSvPR(gO2,gI2,gI1))*CpbarFeChaSvPR(gO1,gI2,gI1)));
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
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFeChaSvPL(gO2,gI2,gI1))*CpbarFeChaSvPL(gO1,gI2,gI1)));
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
   result += 2*SUM(gI1,0,1,A0(Sqr(MCha(gI1)))*(CpbarChaChaUhhPL(gI1,gI1,gO1) +
      CpbarChaChaUhhPR(gI1,gI1,gO1))*MCha(gI1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhh(gI1,gI1,gO1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhh(gI1,gI1,gO1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdUhhPL(gI1,gI1,gO1) +
      CpbarFdFdUhhPR(gI1,gI1,gO1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFeUhhPL(gI1,gI1,gO1) +
      CpbarFeFeUhhPR(gI1,gI1,gO1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuUhhPL(gI1,gI1,gO1) +
      CpbarFuFuUhhPR(gI1,gI1,gO1))*MFu(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFv(gI1)))*(CpbarFvFvUhhPL(gI1,gI1,gO1) +
      CpbarFvFvUhhPR(gI1,gI1,gO1))*MFv(gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhSdconjSd(gO1,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhSeconjSe(gO1,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhSuconjSu(gO1,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSv(gI1)))*CpUhhSvconjSv(gO1,gI1,gI1));
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   sf_data.mu  = Re(0.7071067811865475*vS*Lambdax);
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
   const double amu = Re(-0.7071067811865475*vS*Lambdax);
   const double mg = MGlu;
   const double mAsq = (0.7071067811865475*vS*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
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
   const double amu = Re(-0.7071067811865475*vS*Lambdax);
   const double mg = MGlu;
   const double mAsq = (0.7071067811865475*vS*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
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
   const double amu = Re(-0.7071067811865475*vS*Lambdax);
   const double mg = MGlu;
   const double mAsq = (0.7071067811865475*vS*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
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
   tadpole_2l(2) *= vS;

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

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(UMSSM_info::VZ);
   }

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MVZp_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::VZp))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZp));
   const double p = MVZp;
   const double self_energy = Re(self_energy_VZp_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(UMSSM_info::VZp);
   }

   PHYSICAL(MVZp) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::Sd))
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
         problems.flag_bad_mass(UMSSM_info::Sd, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(UMSSM_info::Sv))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sv());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSv(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Sv_1loop(p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZV,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Sv, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(UMSSM_info::Su))
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
         problems.flag_bad_mass(UMSSM_info::Su, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(UMSSM_info::Se))
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
         problems.flag_bad_mass(UMSSM_info::Se, eigenvalue_error > precision *
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

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::hh))
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
               problems.flag_bad_mass(UMSSM_info::hh);
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
            problems.flag_bad_mass(UMSSM_info::hh, eigenvalue_error > precision
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
      problems.flag_no_pole_mass_convergence(UMSSM_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(UMSSM_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::Ah))
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
               problems.flag_bad_mass(UMSSM_info::Ah);
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
            problems.flag_bad_mass(UMSSM_info::Ah, eigenvalue_error > precision
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
      problems.flag_no_pole_mass_convergence(UMSSM_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(UMSSM_info::Ah);
}

void CLASSNAME::calculate_MHpm_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::Hpm))
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
            problems.flag_bad_mass(UMSSM_info::Hpm, eigenvalue_error >
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
      problems.flag_no_pole_mass_convergence(UMSSM_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(UMSSM_info::Hpm);
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
         problems.flag_bad_mass(UMSSM_info::Chi, eigenvalue_error > precision *
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

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fv());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MFv(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fv_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fv_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fv_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZVL) mix_ZVL;
      decltype(ZVR) mix_ZVR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZVL, mix_ZVR, eigenvalue_error);
      problems.flag_bad_mass(UMSSM_info::Fv, eigenvalue_error > precision * Abs
         (eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZVL, mix_ZVR);
   #endif
      if (es == 0) {
         PHYSICAL(ZVL) = mix_ZVL;
         PHYSICAL(ZVR) = mix_ZVR;
      }
      PHYSICAL(MFv(es)) = Abs(eigen_values(es));
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
      problems.flag_bad_mass(UMSSM_info::Cha, eigenvalue_error > precision *
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
      problems.flag_bad_mass(UMSSM_info::Fe, eigenvalue_error > precision * Abs
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
      problems.flag_bad_mass(UMSSM_info::Fd, eigenvalue_error > precision * Abs
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
      problems.flag_bad_mass(UMSSM_info::Fu, eigenvalue_error > precision * Abs
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

void CLASSNAME::calculate_MVWm_pole()
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(UMSSM_info::VWm);
   }

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(UMSSM_info::VWm);
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(UMSSM_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(UMSSM_info::VZ);
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::calculate_MFv_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fv_1loop_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fv_1loop_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fv_1loop_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (self_energy_PL +
      self_energy_PR);

   return m_drbar;
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
      problems.flag_pole_tachyon(UMSSM_info::VZ);return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole) const
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(UMSSM_info::VWm);return m_pole;
   }

   return AbsSqrt(mass_sqr);
}



double CLASSNAME::v() const
{

   return Sqrt(Sqr(vd) + Sqr(vu));
}

double CLASSNAME::Betax() const
{

   return ArcSin(Abs(ZP(0,1)));
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
