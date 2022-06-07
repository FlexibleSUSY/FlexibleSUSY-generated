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
 * @file TMSSM_mass_eigenstates.cpp
 * @brief implementation of the TMSSM model class
 *
 * Contains the definition of the TMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "TMSSM_mass_eigenstates.hpp"
#include "TMSSM_ewsb_solver_interface.hpp"
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
#include "TMSSM_two_scale_ewsb_solver.hpp"
#endif






#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>

namespace flexiblesusy {

#define STRINGIFY(s) XSTRINGIFY(s)
#define XSTRINGIFY(s) #s
#define CLASSNAME TMSSM_mass_eigenstates

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

CLASSNAME::CLASSNAME(const TMSSM_input_parameters& input_)
   : TMSSM_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new TMSSM_ewsb_solver<Two_scale>())
#endif
{
}

std::unique_ptr<TMSSM_mass_eigenstates_interface> CLASSNAME::clone() const
{
   return std::make_unique<TMSSM_mass_eigenstates>(*this);
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

const TMSSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

TMSSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const TMSSM_physical& physical_)
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

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<TMSSM_ewsb_solver_interface>& solver)
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
   tadpole[2] /= vT;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;

   
   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_mT2 = mT2;

   mHd2 = Re((0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(BMu) + 20*MT*vT*vu
      *Conj(Lambdax) + 10*vT*vu*Conj(TLambdax) + 20*vT*vu*Conj(MT)*Lambdax - 20*vd
      *vT*Conj(Mu)*Lambdax - 20*vd*vT*Conj(Lambdax)*Mu - 3*Cube(vd)*Sqr(g1) - 5*
      Cube(vd)*Sqr(g2) - 10*vd*AbsSqr(Lambdax)*Sqr(vT) - 10*vd*AbsSqr(Lambdax)*Sqr
      (vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) + 10*vT*vu*TLambdax))/vd)
      ;
   mHu2 = Re((0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(BMu) + 20*MT*vd*vT
      *Conj(Lambdax) + 10*vd*vT*Conj(TLambdax) + 20*vd*vT*Conj(MT)*Lambdax - 20*vT
      *vu*Conj(Mu)*Lambdax - 20*vT*vu*Conj(Lambdax)*Mu - 3*Cube(vu)*Sqr(g1) - 5*
      Cube(vu)*Sqr(g2) - 10*vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) + 5*
      vu*Sqr(g2)*Sqr(vd) - 10*vu*AbsSqr(Lambdax)*Sqr(vT) + 10*vd*vT*TLambdax))/vu)
      ;
   mT2 = Re((0.25*(-16*vT*AbsSqr(MT) - 4*vT*BMT - 4*vT*Conj(BMT) + 2*MT*vd*vu*Conj
      (Lambdax) + vd*vu*Conj(TLambdax) + 2*vd*vu*Conj(MT)*Lambdax - vT*AbsSqr(
      Lambdax)*Sqr(vd) - Conj(Mu)*Lambdax*Sqr(vd) - Conj(Lambdax)*Mu*Sqr(vd) - vT*
      AbsSqr(Lambdax)*Sqr(vu) - Conj(Mu)*Lambdax*Sqr(vu) - Conj(Lambdax)*Mu*Sqr(vu
      ) + vd*vu*TLambdax))/vT);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(mT2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      mT2 = old_mT2;
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

   VERBOSE_MSG("\t\tSolving TMSSM EWSB at " << ewsb_loop_order << "-loop order");

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
           "TMSSM\n"
           "========================================\n";
   TMSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

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
   const auto save_mT2_raii = make_raii_save(mT2);

   solve_ewsb_tree_level_custom();

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha();
   calculate_MChi();
   calculate_MHpm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
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
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 18u));

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
      calculate_MChi_pole();
      calculate_MGlu_pole();
      calculate_Mhh_pole();
      calculate_MHpm_pole();
      calculate_MSd_pole();
      calculate_MSe_pole();
      calculate_MSu_pole();
      calculate_MSv_pole();
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
   move_goldstone_to(0, MVWm, PHYSICAL(MHpm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::Sd); }
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::Sv); }
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::Su); }
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::Se); }
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::hh); }
   if (PHYSICAL(MAh).tail<2>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::Ah); }
   if (PHYSICAL(MHpm).tail<3>().minCoeff() < 0.) { problems.flag_pole_tachyon(TMSSM_info::Hpm); }
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
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,4,1>::Zero();
   ZP = Eigen::Matrix<double,4,4>::Zero();
   MChi = Eigen::Matrix<double,5,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,5,5>::Zero();
   MCha = Eigen::Matrix<double,3,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   UP = Eigen::Matrix<std::complex<double>,3,3>::Zero();
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

   PhaseGlu = std::complex<double>(1.,0.);


}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   TMSSM_soft_parameters::clear();
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
   MSd(0) = pars(5);
   MSd(1) = pars(6);
   MSd(2) = pars(7);
   MSd(3) = pars(8);
   MSd(4) = pars(9);
   MSd(5) = pars(10);
   MSv(0) = pars(11);
   MSv(1) = pars(12);
   MSv(2) = pars(13);
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
   MHpm(2) = pars(34);
   MHpm(3) = pars(35);
   MChi(0) = pars(36);
   MChi(1) = pars(37);
   MChi(2) = pars(38);
   MChi(3) = pars(39);
   MChi(4) = pars(40);
   MCha(0) = pars(41);
   MCha(1) = pars(42);
   MCha(2) = pars(43);
   MFe(0) = pars(44);
   MFe(1) = pars(45);
   MFe(2) = pars(46);
   MFd(0) = pars(47);
   MFd(1) = pars(48);
   MFd(2) = pars(49);
   MFu(0) = pars(50);
   MFu(1) = pars(51);
   MFu(2) = pars(52);
   MVWm = pars(53);
   MVP = pars(54);
   MVZ = pars(55);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(56);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSd(0);
   pars(6) = MSd(1);
   pars(7) = MSd(2);
   pars(8) = MSd(3);
   pars(9) = MSd(4);
   pars(10) = MSd(5);
   pars(11) = MSv(0);
   pars(12) = MSv(1);
   pars(13) = MSv(2);
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
   pars(34) = MHpm(2);
   pars(35) = MHpm(3);
   pars(36) = MChi(0);
   pars(37) = MChi(1);
   pars(38) = MChi(2);
   pars(39) = MChi(3);
   pars(40) = MChi(4);
   pars(41) = MCha(0);
   pars(42) = MCha(1);
   pars(43) = MCha(2);
   pars(44) = MFe(0);
   pars(45) = MFe(1);
   pars(46) = MFe(2);
   pars(47) = MFd(0);
   pars(48) = MFd(1);
   pars(49) = MFd(2);
   pars(50) = MFu(0);
   pars(51) = MFu(1);
   pars(52) = MFu(2);
   pars(53) = MVWm;
   pars(54) = MVP;
   pars(55) = MVZ;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   ZD(0,0) = pars(56);
   ZD(0,1) = pars(57);
   ZD(0,2) = pars(58);
   ZD(0,3) = pars(59);
   ZD(0,4) = pars(60);
   ZD(0,5) = pars(61);
   ZD(1,0) = pars(62);
   ZD(1,1) = pars(63);
   ZD(1,2) = pars(64);
   ZD(1,3) = pars(65);
   ZD(1,4) = pars(66);
   ZD(1,5) = pars(67);
   ZD(2,0) = pars(68);
   ZD(2,1) = pars(69);
   ZD(2,2) = pars(70);
   ZD(2,3) = pars(71);
   ZD(2,4) = pars(72);
   ZD(2,5) = pars(73);
   ZD(3,0) = pars(74);
   ZD(3,1) = pars(75);
   ZD(3,2) = pars(76);
   ZD(3,3) = pars(77);
   ZD(3,4) = pars(78);
   ZD(3,5) = pars(79);
   ZD(4,0) = pars(80);
   ZD(4,1) = pars(81);
   ZD(4,2) = pars(82);
   ZD(4,3) = pars(83);
   ZD(4,4) = pars(84);
   ZD(4,5) = pars(85);
   ZD(5,0) = pars(86);
   ZD(5,1) = pars(87);
   ZD(5,2) = pars(88);
   ZD(5,3) = pars(89);
   ZD(5,4) = pars(90);
   ZD(5,5) = pars(91);
   ZV(0,0) = pars(92);
   ZV(0,1) = pars(93);
   ZV(0,2) = pars(94);
   ZV(1,0) = pars(95);
   ZV(1,1) = pars(96);
   ZV(1,2) = pars(97);
   ZV(2,0) = pars(98);
   ZV(2,1) = pars(99);
   ZV(2,2) = pars(100);
   ZU(0,0) = pars(101);
   ZU(0,1) = pars(102);
   ZU(0,2) = pars(103);
   ZU(0,3) = pars(104);
   ZU(0,4) = pars(105);
   ZU(0,5) = pars(106);
   ZU(1,0) = pars(107);
   ZU(1,1) = pars(108);
   ZU(1,2) = pars(109);
   ZU(1,3) = pars(110);
   ZU(1,4) = pars(111);
   ZU(1,5) = pars(112);
   ZU(2,0) = pars(113);
   ZU(2,1) = pars(114);
   ZU(2,2) = pars(115);
   ZU(2,3) = pars(116);
   ZU(2,4) = pars(117);
   ZU(2,5) = pars(118);
   ZU(3,0) = pars(119);
   ZU(3,1) = pars(120);
   ZU(3,2) = pars(121);
   ZU(3,3) = pars(122);
   ZU(3,4) = pars(123);
   ZU(3,5) = pars(124);
   ZU(4,0) = pars(125);
   ZU(4,1) = pars(126);
   ZU(4,2) = pars(127);
   ZU(4,3) = pars(128);
   ZU(4,4) = pars(129);
   ZU(4,5) = pars(130);
   ZU(5,0) = pars(131);
   ZU(5,1) = pars(132);
   ZU(5,2) = pars(133);
   ZU(5,3) = pars(134);
   ZU(5,4) = pars(135);
   ZU(5,5) = pars(136);
   ZE(0,0) = pars(137);
   ZE(0,1) = pars(138);
   ZE(0,2) = pars(139);
   ZE(0,3) = pars(140);
   ZE(0,4) = pars(141);
   ZE(0,5) = pars(142);
   ZE(1,0) = pars(143);
   ZE(1,1) = pars(144);
   ZE(1,2) = pars(145);
   ZE(1,3) = pars(146);
   ZE(1,4) = pars(147);
   ZE(1,5) = pars(148);
   ZE(2,0) = pars(149);
   ZE(2,1) = pars(150);
   ZE(2,2) = pars(151);
   ZE(2,3) = pars(152);
   ZE(2,4) = pars(153);
   ZE(2,5) = pars(154);
   ZE(3,0) = pars(155);
   ZE(3,1) = pars(156);
   ZE(3,2) = pars(157);
   ZE(3,3) = pars(158);
   ZE(3,4) = pars(159);
   ZE(3,5) = pars(160);
   ZE(4,0) = pars(161);
   ZE(4,1) = pars(162);
   ZE(4,2) = pars(163);
   ZE(4,3) = pars(164);
   ZE(4,4) = pars(165);
   ZE(4,5) = pars(166);
   ZE(5,0) = pars(167);
   ZE(5,1) = pars(168);
   ZE(5,2) = pars(169);
   ZE(5,3) = pars(170);
   ZE(5,4) = pars(171);
   ZE(5,5) = pars(172);
   ZH(0,0) = pars(173);
   ZH(0,1) = pars(174);
   ZH(0,2) = pars(175);
   ZH(1,0) = pars(176);
   ZH(1,1) = pars(177);
   ZH(1,2) = pars(178);
   ZH(2,0) = pars(179);
   ZH(2,1) = pars(180);
   ZH(2,2) = pars(181);
   ZA(0,0) = pars(182);
   ZA(0,1) = pars(183);
   ZA(0,2) = pars(184);
   ZA(1,0) = pars(185);
   ZA(1,1) = pars(186);
   ZA(1,2) = pars(187);
   ZA(2,0) = pars(188);
   ZA(2,1) = pars(189);
   ZA(2,2) = pars(190);
   ZP(0,0) = pars(191);
   ZP(0,1) = pars(192);
   ZP(0,2) = pars(193);
   ZP(0,3) = pars(194);
   ZP(1,0) = pars(195);
   ZP(1,1) = pars(196);
   ZP(1,2) = pars(197);
   ZP(1,3) = pars(198);
   ZP(2,0) = pars(199);
   ZP(2,1) = pars(200);
   ZP(2,2) = pars(201);
   ZP(2,3) = pars(202);
   ZP(3,0) = pars(203);
   ZP(3,1) = pars(204);
   ZP(3,2) = pars(205);
   ZP(3,3) = pars(206);
   ZN(0,0) = std::complex<double>(pars(207), pars(208));
   ZN(0,1) = std::complex<double>(pars(209), pars(210));
   ZN(0,2) = std::complex<double>(pars(211), pars(212));
   ZN(0,3) = std::complex<double>(pars(213), pars(214));
   ZN(0,4) = std::complex<double>(pars(215), pars(216));
   ZN(1,0) = std::complex<double>(pars(217), pars(218));
   ZN(1,1) = std::complex<double>(pars(219), pars(220));
   ZN(1,2) = std::complex<double>(pars(221), pars(222));
   ZN(1,3) = std::complex<double>(pars(223), pars(224));
   ZN(1,4) = std::complex<double>(pars(225), pars(226));
   ZN(2,0) = std::complex<double>(pars(227), pars(228));
   ZN(2,1) = std::complex<double>(pars(229), pars(230));
   ZN(2,2) = std::complex<double>(pars(231), pars(232));
   ZN(2,3) = std::complex<double>(pars(233), pars(234));
   ZN(2,4) = std::complex<double>(pars(235), pars(236));
   ZN(3,0) = std::complex<double>(pars(237), pars(238));
   ZN(3,1) = std::complex<double>(pars(239), pars(240));
   ZN(3,2) = std::complex<double>(pars(241), pars(242));
   ZN(3,3) = std::complex<double>(pars(243), pars(244));
   ZN(3,4) = std::complex<double>(pars(245), pars(246));
   ZN(4,0) = std::complex<double>(pars(247), pars(248));
   ZN(4,1) = std::complex<double>(pars(249), pars(250));
   ZN(4,2) = std::complex<double>(pars(251), pars(252));
   ZN(4,3) = std::complex<double>(pars(253), pars(254));
   ZN(4,4) = std::complex<double>(pars(255), pars(256));
   UM(0,0) = std::complex<double>(pars(257), pars(258));
   UM(0,1) = std::complex<double>(pars(259), pars(260));
   UM(0,2) = std::complex<double>(pars(261), pars(262));
   UM(1,0) = std::complex<double>(pars(263), pars(264));
   UM(1,1) = std::complex<double>(pars(265), pars(266));
   UM(1,2) = std::complex<double>(pars(267), pars(268));
   UM(2,0) = std::complex<double>(pars(269), pars(270));
   UM(2,1) = std::complex<double>(pars(271), pars(272));
   UM(2,2) = std::complex<double>(pars(273), pars(274));
   UP(0,0) = std::complex<double>(pars(275), pars(276));
   UP(0,1) = std::complex<double>(pars(277), pars(278));
   UP(0,2) = std::complex<double>(pars(279), pars(280));
   UP(1,0) = std::complex<double>(pars(281), pars(282));
   UP(1,1) = std::complex<double>(pars(283), pars(284));
   UP(1,2) = std::complex<double>(pars(285), pars(286));
   UP(2,0) = std::complex<double>(pars(287), pars(288));
   UP(2,1) = std::complex<double>(pars(289), pars(290));
   UP(2,2) = std::complex<double>(pars(291), pars(292));
   ZEL(0,0) = std::complex<double>(pars(293), pars(294));
   ZEL(0,1) = std::complex<double>(pars(295), pars(296));
   ZEL(0,2) = std::complex<double>(pars(297), pars(298));
   ZEL(1,0) = std::complex<double>(pars(299), pars(300));
   ZEL(1,1) = std::complex<double>(pars(301), pars(302));
   ZEL(1,2) = std::complex<double>(pars(303), pars(304));
   ZEL(2,0) = std::complex<double>(pars(305), pars(306));
   ZEL(2,1) = std::complex<double>(pars(307), pars(308));
   ZEL(2,2) = std::complex<double>(pars(309), pars(310));
   ZER(0,0) = std::complex<double>(pars(311), pars(312));
   ZER(0,1) = std::complex<double>(pars(313), pars(314));
   ZER(0,2) = std::complex<double>(pars(315), pars(316));
   ZER(1,0) = std::complex<double>(pars(317), pars(318));
   ZER(1,1) = std::complex<double>(pars(319), pars(320));
   ZER(1,2) = std::complex<double>(pars(321), pars(322));
   ZER(2,0) = std::complex<double>(pars(323), pars(324));
   ZER(2,1) = std::complex<double>(pars(325), pars(326));
   ZER(2,2) = std::complex<double>(pars(327), pars(328));
   ZDL(0,0) = std::complex<double>(pars(329), pars(330));
   ZDL(0,1) = std::complex<double>(pars(331), pars(332));
   ZDL(0,2) = std::complex<double>(pars(333), pars(334));
   ZDL(1,0) = std::complex<double>(pars(335), pars(336));
   ZDL(1,1) = std::complex<double>(pars(337), pars(338));
   ZDL(1,2) = std::complex<double>(pars(339), pars(340));
   ZDL(2,0) = std::complex<double>(pars(341), pars(342));
   ZDL(2,1) = std::complex<double>(pars(343), pars(344));
   ZDL(2,2) = std::complex<double>(pars(345), pars(346));
   ZDR(0,0) = std::complex<double>(pars(347), pars(348));
   ZDR(0,1) = std::complex<double>(pars(349), pars(350));
   ZDR(0,2) = std::complex<double>(pars(351), pars(352));
   ZDR(1,0) = std::complex<double>(pars(353), pars(354));
   ZDR(1,1) = std::complex<double>(pars(355), pars(356));
   ZDR(1,2) = std::complex<double>(pars(357), pars(358));
   ZDR(2,0) = std::complex<double>(pars(359), pars(360));
   ZDR(2,1) = std::complex<double>(pars(361), pars(362));
   ZDR(2,2) = std::complex<double>(pars(363), pars(364));
   ZUL(0,0) = std::complex<double>(pars(365), pars(366));
   ZUL(0,1) = std::complex<double>(pars(367), pars(368));
   ZUL(0,2) = std::complex<double>(pars(369), pars(370));
   ZUL(1,0) = std::complex<double>(pars(371), pars(372));
   ZUL(1,1) = std::complex<double>(pars(373), pars(374));
   ZUL(1,2) = std::complex<double>(pars(375), pars(376));
   ZUL(2,0) = std::complex<double>(pars(377), pars(378));
   ZUL(2,1) = std::complex<double>(pars(379), pars(380));
   ZUL(2,2) = std::complex<double>(pars(381), pars(382));
   ZUR(0,0) = std::complex<double>(pars(383), pars(384));
   ZUR(0,1) = std::complex<double>(pars(385), pars(386));
   ZUR(0,2) = std::complex<double>(pars(387), pars(388));
   ZUR(1,0) = std::complex<double>(pars(389), pars(390));
   ZUR(1,1) = std::complex<double>(pars(391), pars(392));
   ZUR(1,2) = std::complex<double>(pars(393), pars(394));
   ZUR(2,0) = std::complex<double>(pars(395), pars(396));
   ZUR(2,1) = std::complex<double>(pars(397), pars(398));
   ZUR(2,2) = std::complex<double>(pars(399), pars(400));
   ZZ(0,0) = pars(401);
   ZZ(0,1) = pars(402);
   ZZ(1,0) = pars(403);
   ZZ(1,1) = pars(404);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(405);

   pars(56) = ZD(0,0);
   pars(57) = ZD(0,1);
   pars(58) = ZD(0,2);
   pars(59) = ZD(0,3);
   pars(60) = ZD(0,4);
   pars(61) = ZD(0,5);
   pars(62) = ZD(1,0);
   pars(63) = ZD(1,1);
   pars(64) = ZD(1,2);
   pars(65) = ZD(1,3);
   pars(66) = ZD(1,4);
   pars(67) = ZD(1,5);
   pars(68) = ZD(2,0);
   pars(69) = ZD(2,1);
   pars(70) = ZD(2,2);
   pars(71) = ZD(2,3);
   pars(72) = ZD(2,4);
   pars(73) = ZD(2,5);
   pars(74) = ZD(3,0);
   pars(75) = ZD(3,1);
   pars(76) = ZD(3,2);
   pars(77) = ZD(3,3);
   pars(78) = ZD(3,4);
   pars(79) = ZD(3,5);
   pars(80) = ZD(4,0);
   pars(81) = ZD(4,1);
   pars(82) = ZD(4,2);
   pars(83) = ZD(4,3);
   pars(84) = ZD(4,4);
   pars(85) = ZD(4,5);
   pars(86) = ZD(5,0);
   pars(87) = ZD(5,1);
   pars(88) = ZD(5,2);
   pars(89) = ZD(5,3);
   pars(90) = ZD(5,4);
   pars(91) = ZD(5,5);
   pars(92) = ZV(0,0);
   pars(93) = ZV(0,1);
   pars(94) = ZV(0,2);
   pars(95) = ZV(1,0);
   pars(96) = ZV(1,1);
   pars(97) = ZV(1,2);
   pars(98) = ZV(2,0);
   pars(99) = ZV(2,1);
   pars(100) = ZV(2,2);
   pars(101) = ZU(0,0);
   pars(102) = ZU(0,1);
   pars(103) = ZU(0,2);
   pars(104) = ZU(0,3);
   pars(105) = ZU(0,4);
   pars(106) = ZU(0,5);
   pars(107) = ZU(1,0);
   pars(108) = ZU(1,1);
   pars(109) = ZU(1,2);
   pars(110) = ZU(1,3);
   pars(111) = ZU(1,4);
   pars(112) = ZU(1,5);
   pars(113) = ZU(2,0);
   pars(114) = ZU(2,1);
   pars(115) = ZU(2,2);
   pars(116) = ZU(2,3);
   pars(117) = ZU(2,4);
   pars(118) = ZU(2,5);
   pars(119) = ZU(3,0);
   pars(120) = ZU(3,1);
   pars(121) = ZU(3,2);
   pars(122) = ZU(3,3);
   pars(123) = ZU(3,4);
   pars(124) = ZU(3,5);
   pars(125) = ZU(4,0);
   pars(126) = ZU(4,1);
   pars(127) = ZU(4,2);
   pars(128) = ZU(4,3);
   pars(129) = ZU(4,4);
   pars(130) = ZU(4,5);
   pars(131) = ZU(5,0);
   pars(132) = ZU(5,1);
   pars(133) = ZU(5,2);
   pars(134) = ZU(5,3);
   pars(135) = ZU(5,4);
   pars(136) = ZU(5,5);
   pars(137) = ZE(0,0);
   pars(138) = ZE(0,1);
   pars(139) = ZE(0,2);
   pars(140) = ZE(0,3);
   pars(141) = ZE(0,4);
   pars(142) = ZE(0,5);
   pars(143) = ZE(1,0);
   pars(144) = ZE(1,1);
   pars(145) = ZE(1,2);
   pars(146) = ZE(1,3);
   pars(147) = ZE(1,4);
   pars(148) = ZE(1,5);
   pars(149) = ZE(2,0);
   pars(150) = ZE(2,1);
   pars(151) = ZE(2,2);
   pars(152) = ZE(2,3);
   pars(153) = ZE(2,4);
   pars(154) = ZE(2,5);
   pars(155) = ZE(3,0);
   pars(156) = ZE(3,1);
   pars(157) = ZE(3,2);
   pars(158) = ZE(3,3);
   pars(159) = ZE(3,4);
   pars(160) = ZE(3,5);
   pars(161) = ZE(4,0);
   pars(162) = ZE(4,1);
   pars(163) = ZE(4,2);
   pars(164) = ZE(4,3);
   pars(165) = ZE(4,4);
   pars(166) = ZE(4,5);
   pars(167) = ZE(5,0);
   pars(168) = ZE(5,1);
   pars(169) = ZE(5,2);
   pars(170) = ZE(5,3);
   pars(171) = ZE(5,4);
   pars(172) = ZE(5,5);
   pars(173) = ZH(0,0);
   pars(174) = ZH(0,1);
   pars(175) = ZH(0,2);
   pars(176) = ZH(1,0);
   pars(177) = ZH(1,1);
   pars(178) = ZH(1,2);
   pars(179) = ZH(2,0);
   pars(180) = ZH(2,1);
   pars(181) = ZH(2,2);
   pars(182) = ZA(0,0);
   pars(183) = ZA(0,1);
   pars(184) = ZA(0,2);
   pars(185) = ZA(1,0);
   pars(186) = ZA(1,1);
   pars(187) = ZA(1,2);
   pars(188) = ZA(2,0);
   pars(189) = ZA(2,1);
   pars(190) = ZA(2,2);
   pars(191) = ZP(0,0);
   pars(192) = ZP(0,1);
   pars(193) = ZP(0,2);
   pars(194) = ZP(0,3);
   pars(195) = ZP(1,0);
   pars(196) = ZP(1,1);
   pars(197) = ZP(1,2);
   pars(198) = ZP(1,3);
   pars(199) = ZP(2,0);
   pars(200) = ZP(2,1);
   pars(201) = ZP(2,2);
   pars(202) = ZP(2,3);
   pars(203) = ZP(3,0);
   pars(204) = ZP(3,1);
   pars(205) = ZP(3,2);
   pars(206) = ZP(3,3);
   pars(207) = Re(ZN(0,0));
   pars(208) = Im(ZN(0,0));
   pars(209) = Re(ZN(0,1));
   pars(210) = Im(ZN(0,1));
   pars(211) = Re(ZN(0,2));
   pars(212) = Im(ZN(0,2));
   pars(213) = Re(ZN(0,3));
   pars(214) = Im(ZN(0,3));
   pars(215) = Re(ZN(0,4));
   pars(216) = Im(ZN(0,4));
   pars(217) = Re(ZN(1,0));
   pars(218) = Im(ZN(1,0));
   pars(219) = Re(ZN(1,1));
   pars(220) = Im(ZN(1,1));
   pars(221) = Re(ZN(1,2));
   pars(222) = Im(ZN(1,2));
   pars(223) = Re(ZN(1,3));
   pars(224) = Im(ZN(1,3));
   pars(225) = Re(ZN(1,4));
   pars(226) = Im(ZN(1,4));
   pars(227) = Re(ZN(2,0));
   pars(228) = Im(ZN(2,0));
   pars(229) = Re(ZN(2,1));
   pars(230) = Im(ZN(2,1));
   pars(231) = Re(ZN(2,2));
   pars(232) = Im(ZN(2,2));
   pars(233) = Re(ZN(2,3));
   pars(234) = Im(ZN(2,3));
   pars(235) = Re(ZN(2,4));
   pars(236) = Im(ZN(2,4));
   pars(237) = Re(ZN(3,0));
   pars(238) = Im(ZN(3,0));
   pars(239) = Re(ZN(3,1));
   pars(240) = Im(ZN(3,1));
   pars(241) = Re(ZN(3,2));
   pars(242) = Im(ZN(3,2));
   pars(243) = Re(ZN(3,3));
   pars(244) = Im(ZN(3,3));
   pars(245) = Re(ZN(3,4));
   pars(246) = Im(ZN(3,4));
   pars(247) = Re(ZN(4,0));
   pars(248) = Im(ZN(4,0));
   pars(249) = Re(ZN(4,1));
   pars(250) = Im(ZN(4,1));
   pars(251) = Re(ZN(4,2));
   pars(252) = Im(ZN(4,2));
   pars(253) = Re(ZN(4,3));
   pars(254) = Im(ZN(4,3));
   pars(255) = Re(ZN(4,4));
   pars(256) = Im(ZN(4,4));
   pars(257) = Re(UM(0,0));
   pars(258) = Im(UM(0,0));
   pars(259) = Re(UM(0,1));
   pars(260) = Im(UM(0,1));
   pars(261) = Re(UM(0,2));
   pars(262) = Im(UM(0,2));
   pars(263) = Re(UM(1,0));
   pars(264) = Im(UM(1,0));
   pars(265) = Re(UM(1,1));
   pars(266) = Im(UM(1,1));
   pars(267) = Re(UM(1,2));
   pars(268) = Im(UM(1,2));
   pars(269) = Re(UM(2,0));
   pars(270) = Im(UM(2,0));
   pars(271) = Re(UM(2,1));
   pars(272) = Im(UM(2,1));
   pars(273) = Re(UM(2,2));
   pars(274) = Im(UM(2,2));
   pars(275) = Re(UP(0,0));
   pars(276) = Im(UP(0,0));
   pars(277) = Re(UP(0,1));
   pars(278) = Im(UP(0,1));
   pars(279) = Re(UP(0,2));
   pars(280) = Im(UP(0,2));
   pars(281) = Re(UP(1,0));
   pars(282) = Im(UP(1,0));
   pars(283) = Re(UP(1,1));
   pars(284) = Im(UP(1,1));
   pars(285) = Re(UP(1,2));
   pars(286) = Im(UP(1,2));
   pars(287) = Re(UP(2,0));
   pars(288) = Im(UP(2,0));
   pars(289) = Re(UP(2,1));
   pars(290) = Im(UP(2,1));
   pars(291) = Re(UP(2,2));
   pars(292) = Im(UP(2,2));
   pars(293) = Re(ZEL(0,0));
   pars(294) = Im(ZEL(0,0));
   pars(295) = Re(ZEL(0,1));
   pars(296) = Im(ZEL(0,1));
   pars(297) = Re(ZEL(0,2));
   pars(298) = Im(ZEL(0,2));
   pars(299) = Re(ZEL(1,0));
   pars(300) = Im(ZEL(1,0));
   pars(301) = Re(ZEL(1,1));
   pars(302) = Im(ZEL(1,1));
   pars(303) = Re(ZEL(1,2));
   pars(304) = Im(ZEL(1,2));
   pars(305) = Re(ZEL(2,0));
   pars(306) = Im(ZEL(2,0));
   pars(307) = Re(ZEL(2,1));
   pars(308) = Im(ZEL(2,1));
   pars(309) = Re(ZEL(2,2));
   pars(310) = Im(ZEL(2,2));
   pars(311) = Re(ZER(0,0));
   pars(312) = Im(ZER(0,0));
   pars(313) = Re(ZER(0,1));
   pars(314) = Im(ZER(0,1));
   pars(315) = Re(ZER(0,2));
   pars(316) = Im(ZER(0,2));
   pars(317) = Re(ZER(1,0));
   pars(318) = Im(ZER(1,0));
   pars(319) = Re(ZER(1,1));
   pars(320) = Im(ZER(1,1));
   pars(321) = Re(ZER(1,2));
   pars(322) = Im(ZER(1,2));
   pars(323) = Re(ZER(2,0));
   pars(324) = Im(ZER(2,0));
   pars(325) = Re(ZER(2,1));
   pars(326) = Im(ZER(2,1));
   pars(327) = Re(ZER(2,2));
   pars(328) = Im(ZER(2,2));
   pars(329) = Re(ZDL(0,0));
   pars(330) = Im(ZDL(0,0));
   pars(331) = Re(ZDL(0,1));
   pars(332) = Im(ZDL(0,1));
   pars(333) = Re(ZDL(0,2));
   pars(334) = Im(ZDL(0,2));
   pars(335) = Re(ZDL(1,0));
   pars(336) = Im(ZDL(1,0));
   pars(337) = Re(ZDL(1,1));
   pars(338) = Im(ZDL(1,1));
   pars(339) = Re(ZDL(1,2));
   pars(340) = Im(ZDL(1,2));
   pars(341) = Re(ZDL(2,0));
   pars(342) = Im(ZDL(2,0));
   pars(343) = Re(ZDL(2,1));
   pars(344) = Im(ZDL(2,1));
   pars(345) = Re(ZDL(2,2));
   pars(346) = Im(ZDL(2,2));
   pars(347) = Re(ZDR(0,0));
   pars(348) = Im(ZDR(0,0));
   pars(349) = Re(ZDR(0,1));
   pars(350) = Im(ZDR(0,1));
   pars(351) = Re(ZDR(0,2));
   pars(352) = Im(ZDR(0,2));
   pars(353) = Re(ZDR(1,0));
   pars(354) = Im(ZDR(1,0));
   pars(355) = Re(ZDR(1,1));
   pars(356) = Im(ZDR(1,1));
   pars(357) = Re(ZDR(1,2));
   pars(358) = Im(ZDR(1,2));
   pars(359) = Re(ZDR(2,0));
   pars(360) = Im(ZDR(2,0));
   pars(361) = Re(ZDR(2,1));
   pars(362) = Im(ZDR(2,1));
   pars(363) = Re(ZDR(2,2));
   pars(364) = Im(ZDR(2,2));
   pars(365) = Re(ZUL(0,0));
   pars(366) = Im(ZUL(0,0));
   pars(367) = Re(ZUL(0,1));
   pars(368) = Im(ZUL(0,1));
   pars(369) = Re(ZUL(0,2));
   pars(370) = Im(ZUL(0,2));
   pars(371) = Re(ZUL(1,0));
   pars(372) = Im(ZUL(1,0));
   pars(373) = Re(ZUL(1,1));
   pars(374) = Im(ZUL(1,1));
   pars(375) = Re(ZUL(1,2));
   pars(376) = Im(ZUL(1,2));
   pars(377) = Re(ZUL(2,0));
   pars(378) = Im(ZUL(2,0));
   pars(379) = Re(ZUL(2,1));
   pars(380) = Im(ZUL(2,1));
   pars(381) = Re(ZUL(2,2));
   pars(382) = Im(ZUL(2,2));
   pars(383) = Re(ZUR(0,0));
   pars(384) = Im(ZUR(0,0));
   pars(385) = Re(ZUR(0,1));
   pars(386) = Im(ZUR(0,1));
   pars(387) = Re(ZUR(0,2));
   pars(388) = Im(ZUR(0,2));
   pars(389) = Re(ZUR(1,0));
   pars(390) = Im(ZUR(1,0));
   pars(391) = Re(ZUR(1,1));
   pars(392) = Im(ZUR(1,1));
   pars(393) = Re(ZUR(1,2));
   pars(394) = Im(ZUR(1,2));
   pars(395) = Re(ZUR(2,0));
   pars(396) = Im(ZUR(2,0));
   pars(397) = Re(ZUR(2,1));
   pars(398) = Im(ZUR(2,1));
   pars(399) = Re(ZUR(2,2));
   pars(400) = Im(ZUR(2,2));
   pars(401) = ZZ(0,0);
   pars(402) = ZZ(0,1);
   pars(403) = ZZ(1,0);
   pars(404) = ZZ(1,1);


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
   return "TMSSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   TMSSM_soft_parameters::run_to(scale, eps);
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

const TMSSM_input_parameters& CLASSNAME::get_input_parameters() const
{
   return get_input();
}

TMSSM_input_parameters& CLASSNAME::get_input_parameters()
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


Eigen::Array<double,3,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_goldstone;
   MHpm_goldstone(0) = MVWm;

   return remove_if_equal(MHpm, MHpm_goldstone);
}

Eigen::Array<double,2,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;

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
double CLASSNAME::get_lsp(TMSSM_info::Particles& particle_type) const
{
   double lsp_mass = std::numeric_limits<double>::max();
   double tmp_mass;
   particle_type = TMSSM_info::NUMBER_OF_PARTICLES;

   tmp_mass = Abs(PHYSICAL(MChi(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Chi;
   }

   tmp_mass = Abs(PHYSICAL(MSv(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Sv;
   }

   tmp_mass = Abs(PHYSICAL(MSu(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Su;
   }

   tmp_mass = Abs(PHYSICAL(MSd(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Sd;
   }

   tmp_mass = Abs(PHYSICAL(MSe(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Se;
   }

   tmp_mass = Abs(PHYSICAL(MCha(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Cha;
   }

   tmp_mass = Abs(PHYSICAL(MGlu));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = TMSSM_info::Glu;
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

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)) +
      AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) + Conj(
      Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) + Conj(
      Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) -
      0.35355339059327373*vT*vu*Conj(Yd(0,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) -
      0.35355339059327373*vT*vu*Conj(Yd(1,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(1,0))*Mu;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) -
      0.35355339059327373*vT*vu*Conj(Yd(2,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(2,0))*Mu;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) + Conj(
      Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) -
      0.35355339059327373*vT*vu*Conj(Yd(0,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(0,1))*Mu;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) -
      0.35355339059327373*vT*vu*Conj(Yd(1,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(1,1))*Mu;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) -
      0.35355339059327373*vT*vu*Conj(Yd(2,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(2,1))*Mu;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) -
      0.35355339059327373*vT*vu*Conj(Yd(0,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) -
      0.35355339059327373*vT*vu*Conj(Yd(1,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) -
      0.35355339059327373*vT*vu*Conj(Yd(2,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)) +
      AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) + Conj(
      Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) + Conj(
      Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)) +
      AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) + Conj(
      Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)) +
      AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Sd, eigenvalue_error > precision * Abs(
      MSd(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Sv, eigenvalue_error > precision * Abs(
      MSv(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) + Conj(
      Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) + Conj(
      Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) -
      0.35355339059327373*vd*vT*Conj(Yu(0,0))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(0,0))*Mu;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) -
      0.35355339059327373*vd*vT*Conj(Yu(1,0))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(1,0))*Mu;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) -
      0.35355339059327373*vd*vT*Conj(Yu(2,0))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(2,0))*Mu;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) + Conj(
      Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) -
      0.35355339059327373*vd*vT*Conj(Yu(0,1))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(0,1))*Mu;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) -
      0.35355339059327373*vd*vT*Conj(Yu(1,1))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(1,1))*Mu;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) -
      0.35355339059327373*vd*vT*Conj(Yu(2,1))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(2,1))*Mu;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*(AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) +
      0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) -
      0.35355339059327373*vd*vT*Conj(Yu(0,2))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) -
      0.35355339059327373*vd*vT*Conj(Yu(1,2))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) -
      0.35355339059327373*vd*vT*Conj(Yu(2,2))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,0))
      + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) + Conj(
      Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) + Conj(
      Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(1,0))
      + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) + Conj(
      Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(2,0))
      + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Su, eigenvalue_error > precision * Abs(
      MSu(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{

   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)) +
      AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) + Conj(
      Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) + Conj(
      Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) -
      0.35355339059327373*vT*vu*Conj(Ye(0,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(0,0))*Mu;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) -
      0.35355339059327373*vT*vu*Conj(Ye(1,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(1,0))*Mu;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) -
      0.35355339059327373*vT*vu*Conj(Ye(2,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(2,0))*Mu;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) + Conj(
      Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) -
      0.35355339059327373*vT*vu*Conj(Ye(0,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(0,1))*Mu;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) -
      0.35355339059327373*vT*vu*Conj(Ye(1,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(1,1))*Mu;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) -
      0.35355339059327373*vT*vu*Conj(Ye(2,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(2,1))*Mu;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd)
      - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) -
      0.35355339059327373*vT*vu*Conj(Ye(0,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) -
      0.35355339059327373*vT*vu*Conj(Ye(1,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) -
      0.35355339059327373*vT*vu*Conj(Ye(2,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)) +
      AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) + Conj(
      Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) + Conj(
      Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)) +
      AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) + Conj(
      Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)) +
      AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Se, eigenvalue_error > precision * Abs(
      MSe(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax + 0.5*vT*
      Conj(Lambdax)*Mu + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) + 0.25*
      AbsSqr(Lambdax)*Sqr(vT) + 0.25*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*
      Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = 0.5*vd*vu*AbsSqr(Lambdax) - 0.5*BMu - 0.5*Conj(BMu) -
      0.5*MT*vT*Conj(Lambdax) - 0.25*vT*Conj(TLambdax) - 0.5*vT*Conj(MT)*
      Lambdax - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) - 0.25*vT*TLambdax;
   mass_matrix_hh(0,2) = 0.5*vd*vT*AbsSqr(Lambdax) - 0.5*MT*vu*Conj(Lambdax) -
      0.25*vu*Conj(TLambdax) - 0.5*vu*Conj(MT)*Lambdax + 0.5*vd*Conj(Mu)*
      Lambdax + 0.5*vd*Conj(Lambdax)*Mu - 0.25*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax + 0.5*vT*
      Conj(Lambdax)*Mu + 0.25*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vd) -
      0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.225*Sqr(g1)*Sqr(
      vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(1,2) = 0.5*vT*vu*AbsSqr(Lambdax) - 0.5*MT*vd*Conj(Lambdax) -
      0.25*vd*Conj(TLambdax) - 0.5*vd*Conj(MT)*Lambdax + 0.5*vu*Conj(Mu)*
      Lambdax + 0.5*vu*Conj(Lambdax)*Mu - 0.25*vd*TLambdax;
   mass_matrix_hh(2,2) = mT2 + 4*AbsSqr(MT) + BMT + Conj(BMT) + 0.25*AbsSqr(
      Lambdax)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::hh, eigenvalue_error > precision * Abs(
      Mhh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax + 0.5*vT*
      Conj(Lambdax)*Mu + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*
      Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(
      Lambdax)*Sqr(vT) + 0.25*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW())) + 0.15*
      Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*Conj(BMu) + 0.5*MT*vT*Conj(Lambdax) +
      0.25*vT*Conj(TLambdax) + 0.5*vT*Conj(MT)*Lambdax - 0.3872983346207417*g1*
      g2*vd*vu*Cos(ThetaW())*Sin(ThetaW()) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW(
      ))) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW())) + 0.25*vT*TLambdax;
   mass_matrix_Ah(0,2) = -0.5*MT*vu*Conj(Lambdax) + 0.25*vu*Conj(TLambdax) -
      0.5*vu*Conj(MT)*Lambdax + 0.25*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax + 0.5*vT*
      Conj(Lambdax)*Mu + 0.25*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vd) -
      0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.3872983346207417
      *g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu) + 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW())) + 0.15*
      Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = -0.5*MT*vd*Conj(Lambdax) + 0.25*vd*Conj(TLambdax) -
      0.5*vd*Conj(MT)*Lambdax + 0.25*vd*TLambdax;
   mass_matrix_Ah(2,2) = mT2 + 4*AbsSqr(MT) - BMT - Conj(BMT) + 0.25*AbsSqr(
      Lambdax)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vu);

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Ah, eigenvalue_error > precision * Abs(
      MAh(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Hpm() const
{

   Eigen::Matrix<double,4,4> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + AbsSqr(Mu) - 0.5*vT*Conj(Mu)*Lambdax - 0.5*vT*
      Conj(Lambdax)*Mu + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) + 0.25*
      AbsSqr(Lambdax)*Sqr(vT) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr
      (vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = 0.25*vd*vu*AbsSqr(Lambdax) + Conj(BMu) - MT*vT*Conj(
      Lambdax) - 0.5*vT*Conj(TLambdax);
   mass_matrix_Hpm(0,2) = -0.35355339059327373*vd*vT*AbsSqr(Lambdax) -
      1.4142135623730951*MT*vu*Conj(Lambdax) + 0.7071067811865475*vd*Conj(Mu)*
      Lambdax + 0.7071067811865475*vd*vT*Sqr(g2);
   mass_matrix_Hpm(0,3) = 0.35355339059327373*vd*vT*AbsSqr(Lambdax) -
      0.7071067811865475*vu*Conj(TLambdax) + 0.7071067811865475*vd*Conj(Lambdax
      )*Mu;
   mass_matrix_Hpm(1,1) = mHu2 + AbsSqr(Mu) - 0.5*vT*Conj(Mu)*Lambdax - 0.5*vT*
      Conj(Lambdax)*Mu + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vd) +
      0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.075*Sqr(g1)*Sqr(
      vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(1,2) = -0.35355339059327373*vT*vu*AbsSqr(Lambdax) -
      0.7071067811865475*vu*Conj(Mu)*Lambdax + 0.7071067811865475*vd*TLambdax;
   mass_matrix_Hpm(1,3) = 0.35355339059327373*vT*vu*AbsSqr(Lambdax) +
      1.4142135623730951*vd*Conj(MT)*Lambdax - 0.7071067811865475*vu*Conj(
      Lambdax)*Mu - 0.7071067811865475*vT*vu*Sqr(g2);
   mass_matrix_Hpm(2,2) = mT2 + 4*AbsSqr(MT) + 0.5*AbsSqr(Lambdax)*Sqr(vd) -
      0.25*Sqr(g2)*Sqr(vd) + Sqr(g2)*Sqr(vT) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(2,3) = 2*Conj(BMT);
   mass_matrix_Hpm(3,3) = mT2 + 4*AbsSqr(MT) + 0.25*Sqr(g2)*Sqr(vd) + Sqr(g2)*
      Sqr(vT) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.25*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Hpm, eigenvalue_error > precision * Abs(
      MHpm(0)));
#else

   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(TMSSM_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,5,5> CLASSNAME::get_mass_matrix_Chi() const
{

   Eigen::Matrix<double,5,5> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(0,4) = 0;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(1,4) = 0;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -0.5*vT*Lambdax - Mu;
   mass_matrix_Chi(2,4) = -0.5*vu*Lambdax;
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.5*vd*Lambdax;
   mass_matrix_Chi(4,4) = 2*MT;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Chi, eigenvalue_error > precision * Abs(
      MChi(0)));
#else

   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
   normalize_to_interval(ZN);

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Cha() const
{

   Eigen::Matrix<double,3,3> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(0,2) = -(g2*vT);
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = -0.5*vT*Lambdax + Mu;
   mass_matrix_Cha(1,2) = -0.7071067811865475*vu*Lambdax;
   mass_matrix_Cha(2,0) = g2*vT;
   mass_matrix_Cha(2,1) = 0.7071067811865475*vd*Lambdax;
   mass_matrix_Cha(2,2) = 2*MT;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Cha, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(TMSSM_info::Fe, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(TMSSM_info::Fd, eigenvalue_error > precision * Abs(
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
   problems.flag_bad_mass(TMSSM_info::Fu, eigenvalue_error > precision * Abs(
      MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{

   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(vd) + 4*Sqr(vT) + Sqr(vu
      )));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{

   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(TMSSM_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{

   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(vd) -
      0.19364916731037085*g1*g2*Sqr(vu);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(vd) + 0.25*Sqr(g2)*Sqr(vu);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void CLASSNAME::calculate_MVPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());
   Eigen::Array<double,2,1> MVPVZ;


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error);
#else

   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}



double CLASSNAME::get_ewsb_eq_hh_1() const
{
   
   double result = Re(mHd2*vd + vd*AbsSqr(Mu) - 0.5*vu*BMu - 0.5*vu*Conj(BMu) -
      0.5*MT*vT*vu*Conj(Lambdax) - 0.25*vT*vu*Conj(TLambdax) - 0.5*vT*vu*Conj(MT)*
      Lambdax + 0.5*vd*vT*Conj(Mu)*Lambdax + 0.5*vd*vT*Conj(Lambdax)*Mu + 0.075*
      Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.25*vd*AbsSqr(Lambdax)*Sqr(vT)
      + 0.25*vd*AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr(
      g2)*Sqr(vu) - 0.25*vT*vu*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   
   double result = Re(mHu2*vu + vu*AbsSqr(Mu) - 0.5*vd*BMu - 0.5*vd*Conj(BMu) -
      0.5*MT*vd*vT*Conj(Lambdax) - 0.25*vd*vT*Conj(TLambdax) - 0.5*vd*vT*Conj(MT)*
      Lambdax + 0.5*vT*vu*Conj(Mu)*Lambdax + 0.5*vT*vu*Conj(Lambdax)*Mu + 0.075*
      Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) + 0.25*vu*AbsSqr(Lambdax)*Sqr(vd)
      - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.25*vu*AbsSqr(
      Lambdax)*Sqr(vT) - 0.25*vd*vT*TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   
   double result = Re(mT2*vT + 4*vT*AbsSqr(MT) + vT*BMT + vT*Conj(BMT) - 0.5*MT*vd
      *vu*Conj(Lambdax) - 0.25*vd*vu*Conj(TLambdax) - 0.5*vd*vu*Conj(MT)*Lambdax +
      0.25*vT*AbsSqr(Lambdax)*Sqr(vd) + 0.25*Conj(Mu)*Lambdax*Sqr(vd) + 0.25*Conj(
      Lambdax)*Mu*Sqr(vd) + 0.25*vT*AbsSqr(Lambdax)*Sqr(vu) + 0.25*Conj(Mu)*
      Lambdax*Sqr(vu) + 0.25*Conj(Lambdax)*Mu*Sqr(vu) - 0.25*vd*vu*TLambdax);

   return result;
}



std::complex<double> CLASSNAME::CpUSdconjUSdVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.2581988897471611*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW()),0) + IF(gO1 < 3,0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW())),0) + IF(gO1 < 3,
      0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW())),0) +
      0.13333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSdconjUSdconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(
      j1,gO2))*Yd(j1,gO1))*ZA(gI1,0)*ZA(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-
      0.35355339059327373*Conj(Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(
      j1,gO1))*ZA(gI1,2)*ZA(gI2,1),0) + IF(gO1 < 3,-0.35355339059327373*Conj(
      Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1))*ZA(gI1,1)*ZA(gI2,
      2),0) + IF(gO2 < 3,-0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*
      KroneckerDelta(gO1,3 + j1))*ZA(gI1,2)*ZA(gI2,1),0) + IF(gO2 < 3,-
      0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3
       + j1))*ZA(gI1,1)*ZA(gI2,2),0) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) - SUM(j3,0,2,
      KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,
      Conj(Yd(j3,j1))*Yd(j2,j1))))*ZA(gI1,0)*ZA(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(
      j1,gO2))*Yd(j1,gO1))*ZH(gI1,0)*ZH(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,
      0.35355339059327373*Conj(Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(
      j1,gO1))*ZH(gI1,2)*ZH(gI2,1),0) + IF(gO1 < 3,0.35355339059327373*Conj(
      Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1))*ZH(gI1,1)*ZH(gI2,
      2),0) + IF(gO2 < 3,0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*
      KroneckerDelta(gO1,3 + j1))*ZH(gI1,2)*ZH(gI2,1),0) + IF(gO2 < 3,
      0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3
       + j1))*ZH(gI1,1)*ZH(gI2,2),0) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) - SUM(j3,0,2,
      KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,
      Conj(Yd(j3,j1))*Yd(j2,j1))))*ZH(gI1,0)*ZH(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSdSvconjUSdconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gI1,gI2)
      *KroneckerDelta(gO1,gO2)*Sqr(g2),0) + 0.1*KroneckerDelta(gI1,gI2)*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

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
   
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO2))*
      ZDR(gI1,j1))*ZN(gI2,2)),0) - 0.3651483716701107*g1*SUM(j1,0,2,KroneckerDelta
      (gO2,3 + j1)*ZDR(gI1,j1))*ZN(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZDL
      (gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZDL(
      gI1,gO1))*Conj(ZN(gI2,1)),0) - Conj(ZN(gI2,2))*SUM(j2,0,2,Conj(ZDL(gI1,j2))*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSdconjHpmconjUSd(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(
      j1,gO2))*Yu(j1,gO1))*ZP(gI1,1)*ZP(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(Conj(
      Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1))*ZP(gI1,1)*ZP(gI2,
      2)),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2)
      ,0) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3),0)
      + IF(gO2 < 3,-(Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1
      ))*ZP(gI1,2)*ZP(gI2,1)),0) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 +
      j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) - SUM(j3,0,2,
      KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,
      Conj(Yd(j3,j1))*Yd(j2,j1))))*ZP(gI1,0)*ZP(gI2,0) - 0.1*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yd(j1,
      gO1)*ZD(gI1,3 + j1))*SUM(j3,0,2,Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*ZD
      (gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(
      gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-1.3333333333333333*Conj(ZD(gI2,gO2))
      *Sqr(g3)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 +
      j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,3 +
      j2))*ZD(gI1,3 + j2)),0) + IF(gO1 < 3,-3*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1
      )*Yd(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3)))*ZD
      (gI1,j4)),0) + IF(gO1 < 3,-0.016666666666666666*Sqr(g1)*SUM(j1,0,2,Conj(ZD(
      gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,
      0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,-0.016666666666666666*Sqr(g1)*SUM(
      j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(
      gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*
      KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO2 < 3,-
      0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZD(gI2,gO2)
      )*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2
      < 3,-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)),0) + IF(gO2 < 3,
      0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1,3
       + j2)*ZD(gI1,3 + j2)),0) + IF(gO2 < 3,-3*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1
      ,0,2,Yd(j1,j2)*ZD(gI1,3 + j1)))*SUM(j3,0,2,Conj(Yd(j3,gO2))*KroneckerDelta(
      gO1,3 + j3)),0) - 0.03333333333333333*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3
       + j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      ZD(gI1,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))
      - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2))
      - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)) - 0.03333333333333333*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*SUM(j2,0
      ,2,KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)) - 0.6666666666666666*Sqr(g3)*
      SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2)) - SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM
      (j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd
      (j3,j4))*KroneckerDelta(gO1,3 + j3))*ZD(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yu(j1,
      gO1)*ZU(gI1,3 + j1))*SUM(j3,0,2,Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1),0),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,j1)
      )*ZU(gI1,j1)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0
      ,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-
      0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)
      ),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(
      gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM
      (j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2,
      Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU
      (gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2
      ,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*ZU(gI1,3 + j2)) - SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))
      *KroneckerDelta(gO1,3 + j3))*ZU(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSdSeconjUSdconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 +
      j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)
      *Yd(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gI1,3 + j3)))*ZE(
      gI2,j4))),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2)
      *ZE(gI2,3 + j1)))*SUM(j3,0,2,Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3))),0
      ) + 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2))
      - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpAhSdconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2)))*ZA(
      gI2,0),0) + IF(gO2 < 3,std::complex<double>(0.,0.35355339059327373)*vT*
      Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1)))*ZA(gI2,1),0) + IF(
      gO2 < 3,std::complex<double>(0.,0.7071067811865475)*Mu*SUM(j1,0,2,Conj(Yd(j1
      ,gO2))*Conj(ZD(gI1,3 + j1)))*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(
      0.,0.35355339059327373)*vu*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3
       + j1)))*ZA(gI2,2),0) - std::complex<double>(0.,0.7071067811865475)*SUM(j2,0
      ,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYd(j1,j2)))*ZA(
      gI2,0) - std::complex<double>(0.,0.35355339059327373)*vT*Conj(Lambdax)*SUM(
      j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZA
      (gI2,1) - std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZA(gI2,1)
      - std::complex<double>(0.,0.35355339059327373)*vu*Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZA(gI2,2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CphhSdconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1
      )*ZH(gI2,0),0) + IF(gO2 < 3,0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0),0) +
      IF(gO2 < 3,-0.7071067811865475*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,
      gO2)))*ZH(gI2,0),0) + IF(gO2 < 3,-(vd*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2
      ,Conj(Yd(j1,gO2))*Yd(j1,j2)))*ZH(gI2,0)),0) + IF(gO2 < 3,-0.05*vu*Conj(ZD(
      gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gO2 < 3,-0.25*vu*Conj(ZD(gI1,gO2))*Sqr(
      g2)*ZH(gI2,1),0) + IF(gO2 < 3,0.35355339059327373*vT*Lambdax*SUM(j1,0,2,Conj
      (Yd(j1,gO2))*Conj(ZD(gI1,3 + j1)))*ZH(gI2,1),0) + IF(gO2 < 3,
      0.7071067811865475*Mu*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1)))*ZH(
      gI2,1),0) + IF(gO2 < 3,0.35355339059327373*vu*Lambdax*SUM(j1,0,2,Conj(Yd(j1,
      gO2))*Conj(ZD(gI1,3 + j1)))*ZH(gI2,2),0) + 0.1*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZD
      (gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) - 0.7071067811865475*SUM
      (j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYd(j1,j2)))*
      ZH(gI2,0) - vd*SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2
      ,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))))*ZH(gI2,0) - 0.1*vu*Sqr(g1)*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) +
      0.35355339059327373*vT*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZH(gI2,1) + 0.7071067811865475*Conj(
      Mu)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,
      j2)))*ZH(gI2,1) + 0.35355339059327373*vu*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(
      gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmSuconjUSd(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZU(
      gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,-0.5*vT*Lambdax*SUM(j1,0,2,Conj(
      Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZP(gI2,0),0) + IF(gO2 < 3,Mu*SUM(j1,0,2,
      Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZP(gI2,0),0) + IF(gO2 < 3,
      0.7071067811865475*vd*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2)
      )*Yd(j1,j2)))*ZP(gI2,0),0) + IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZU(gI1,
      gO2))*Sqr(g2)*ZP(gI2,1),0) + IF(gO2 < 3,SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*Conj
      (TYu(j1,gO2)))*ZP(gI2,1),0) + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZP(gI2,1),0) + IF(
      gO2 < 3,-0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,2),0) + IF(gO2 < 3,
      0.7071067811865475*vd*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1
      )))*ZP(gI2,2),0) + IF(gO2 < 3,0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,3),0)
      + SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYd(j1,
      j2)))*ZP(gI2,0) + 0.7071067811865475*vu*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(
      j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))))*ZP
      (gI2,0) - 0.5*vT*Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZP(gI2,1) + Conj(Mu)*SUM(j2,0,2,Conj(
      ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZP(gI2,1) +
      0.7071067811865475*vd*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))))*ZP(gI2,1)
      - 0.7071067811865475*vu*Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2
      ,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZP(gI2,3);

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
      ThetaW()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Sin(
      ThetaW()),0) + 0.2581988897471611*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

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
   
   const std::complex<double> result = 0.1*KroneckerDelta(gO1,gO2)*(g1*Sin(ThetaW(
      ))*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr
      (Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) +
      5*Sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) +
      5*Sqr(g2))*(ZH(gI1,0)*ZH(gI2,0) - ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSvUSvconjSvconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,IF(gI2 < 3,-0.15*Conj(ZV(gI1,gO2
      ))*Sqr(g1)*ZV(gI2,gO1),0),0) + IF(gI1 < 3,IF(gI2 < 3,-0.25*Conj(ZV(gI1,gO2))
      *Sqr(g2)*ZV(gI2,gO1),0),0) - 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1
      ,gO2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CphhSvconjUSv(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-0.15*vd*Conj(ZV(gI1,gO2))*Sqr(
      g1)*ZH(gI2,0),0) + IF(gI1 < 3,-0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,0),0
      ) + IF(gI1 < 3,0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gI1 < 3,
      0.25*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0);

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
      gI2,gO1))),0);

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
      ZN(gI2,1))*KroneckerDelta(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSvconjHpmconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(
      j1,gO2))*Ye(j1,gO1))*ZP(gI1,0)*ZP(gI2,0)),0),0) + 0.05*KroneckerDelta(gO1,
      gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (3*Sqr(g1) - 5*Sqr(g2))
      *ZP(gI1,1)*ZP(gI2,1) + 10*Sqr(g2)*(ZP(gI1,2)*ZP(gI2,2) - ZP(gI1,3)*ZP(gI2,3)
      ));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjHpmconjUSv(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZE(
      gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*
      Conj(TYe(j1,gO2)))*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,0,
      2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,j2)))*ZP(gI1,0),0) + IF
      (gO2 < 3,-0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) + IF
      (gO2 < 3,-0.5*vT*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1)))*
      ZP(gI1,1),0) + IF(gO2 < 3,Mu*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1)
      ))*ZP(gI1,1),0) + IF(gO2 < 3,-0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,2),0)
      + IF(gO2 < 3,0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,3),0) + IF(gO2 < 3,-
      0.7071067811865475*vu*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1
      )))*ZP(gI1,3),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdUSvconjSdconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((Sqr(g1) + 5*
      Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj
      (ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSvconjSeconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Ye(j1,
      gO1)*ZE(gI2,3 + j1))*SUM(j3,0,2,Conj(Ye(j3,gO2))*Conj(ZE(gI1,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZE(gI1,gO2))*Sqr(g2)*ZE(gI2,gO1),0),0) +
      0.05*KroneckerDelta(gO1,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(
      gI1,j1))*ZE(gI2,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 +
      j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuUSvconjSuconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((Sqr(g1) - 5*
      Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj
      (ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSvVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZV(gI2,gO2))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW
      ()),0);

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
   
   const std::complex<double> result = IF(gO1 < 3,-0.2581988897471611*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW()),0) + IF(gO1 < 3,0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW())),0) + IF(gO1 < 3,
      0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW())),0) +
      0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

double CLASSNAME::CpUSuconjUSuconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(
      j1,gO2))*Yu(j1,gO1))*ZA(gI1,1)*ZA(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-
      0.35355339059327373*Conj(Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(
      j1,gO1))*ZA(gI1,2)*ZA(gI2,0),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr
      (g2)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.35355339059327373*Conj(Lambdax)*
      SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1))*ZA(gI1,0)*ZA(gI2,2),0) +
      IF(gO2 < 3,-0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*
      KroneckerDelta(gO1,3 + j1))*ZA(gI1,2)*ZA(gI2,0),0) + IF(gO2 < 3,-
      0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3
       + j1))*ZA(gI1,0)*ZA(gI2,2),0) - 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) + 0.2*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,
      1) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZA(gI1,1)*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(
      j1,gO2))*Yu(j1,gO1))*ZH(gI1,1)*ZH(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,
      0.35355339059327373*Conj(Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(
      j1,gO1))*ZH(gI1,2)*ZH(gI2,0),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr
      (g2)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.35355339059327373*Conj(Lambdax)*
      SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1))*ZH(gI1,0)*ZH(gI2,2),0) +
      IF(gO2 < 3,0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*
      KroneckerDelta(gO1,3 + j1))*ZH(gI1,2)*ZH(gI2,0),0) + IF(gO2 < 3,
      0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3
       + j1))*ZH(gI1,0)*ZH(gI2,2),0) - 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) + 0.2*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,
      1) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZH(gI1,1)*ZH(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSuSvconjUSuconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gI1,gI2
      )*KroneckerDelta(gO1,gO2)*Sqr(g2),0) - 0.2*KroneckerDelta(gI1,gI2)*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

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

std::complex<double> CLASSNAME::CpChiFuconjUSuPR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO2))*
      ZUR(gI1,j1))*ZN(gI2,3)),0) + 0.7302967433402214*g1*SUM(j1,0,2,KroneckerDelta
      (gO2,3 + j1)*ZUR(gI1,j1))*ZN(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZN(
      gI2,0))*Conj(ZUL(gI1,gO1)),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZN(
      gI2,1))*Conj(ZUL(gI1,gO1)),0) - Conj(ZN(gI2,3))*SUM(j2,0,2,Conj(ZUL(gI1,j2))
      *SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSuconjHpmconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(
      j1,gO2))*Yd(j1,gO1))*ZP(gI1,0)*ZP(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,Conj(
      Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1))*ZP(gI1,3)*ZP(gI2,
      0),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)
      ,0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0
      ) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2),0) +
      IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3),0) + IF(
      gO2 < 3,Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1))*ZP(
      gI1,0)*ZP(gI2,3),0) - 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1) -
      SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*
      SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjHpmconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZD(
      gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      Conj(TYd(j1,gO2)))*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,0,
      2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,j2)))*ZP(gI1,0),0) + IF
      (gO2 < 3,-0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) + IF
      (gO2 < 3,-0.5*vT*Lambdax*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1)))*
      ZP(gI1,1),0) + IF(gO2 < 3,Mu*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1)
      ))*ZP(gI1,1),0) + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,Conj(ZD(gI2,j2
      ))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZP(gI1,1),0) + IF(gO2 < 3,-0.5*vT
      *Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,2),0) + IF(gO2 < 3,0.5*vT*Conj(ZD(gI2,gO2)
      )*Sqr(g2)*ZP(gI1,3),0) + IF(gO2 < 3,-0.7071067811865475*vu*Lambdax*SUM(j1,0,
      2,Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1)))*ZP(gI1,3),0) - 0.5*vT*Conj(Lambdax)
      *SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)
      ))*ZP(gI1,0) + Conj(Mu)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZP(gI1,0) + 0.7071067811865475*vu*SUM
      (j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,
      2,Conj(Yd(j3,j1))*Yu(j2,j1))))*ZP(gI1,0) + SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(
      j1,0,2,KroneckerDelta(gO2,3 + j1)*TYu(j1,j2)))*ZP(gI1,1) +
      0.7071067811865475*vd*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))))*ZP(gI1,1)
      + 0.7071067811865475*vd*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2
      ,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZP(gI1,2);

   return result;
}

std::complex<double> CLASSNAME::CpSeUSuconjSeconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 +
      j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta
      (gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-
      0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,
      3 + j2)),0) - 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2
      ,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1
      ,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2))
      + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yd(j1,
      gO1)*ZD(gI1,3 + j1))*SUM(j3,0,2,Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1),0),0) +
      IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,j1)
      )*ZD(gI1,j1)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0
      ,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-
      0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)
      ),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(
      gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + 0.1*Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
      KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD
      (gI1,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 +
      j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,
      3 + j1))*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(
      gI2,3 + j2))*ZD(gI1,3 + j2)) - SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))
      *KroneckerDelta(gO1,3 + j3))*ZD(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yu(j1,
      gO1)*ZU(gI1,3 + j1))*SUM(j3,0,2,Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666*Conj(ZU(gI2,gO2))*Sqr(g1)*ZU
      (gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(
      gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-1.3333333333333333*Conj(ZU(gI2,gO2))
      *Sqr(g3)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1
      ))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2
      ))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,-3*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*
      Yu(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3)))*ZU(
      gI1,j4)),0) + IF(gO1 < 3,0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,
      3 + j1))*KroneckerDelta(gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,
      0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,0.03333333333333333*Sqr(g1)*SUM(j2,
      0,2,Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO1
       < 3,0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*
      KroneckerDelta(gO2,3 + j2))*ZU(gI1,gO1),0) + IF(gO2 < 3,0.03333333333333333*
      Conj(ZU(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 +
      j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)*SUM(j1,0,2
      ,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1)),0) + IF(gO2 < 3,
      0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1,
      3 + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZU(gI2,gO2))
      *Sqr(g3)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2)),0) + IF(gO2 <
      3,-3*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1)))*SUM(
      j3,0,2,Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3)),0) - 0.13333333333333333
      *Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1))*SUM(j2,0,2,
      Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)) - 0.6666666666666666*Sqr(g3
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1))*SUM(j2,0,2,Conj(ZU(
      gI2,3 + j2))*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZU(
      gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(
      gO2,3 + j2)) - 0.4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))*
      SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2
      ,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)) - 0.4*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(
      gI1,3 + j2)) - 0.13333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 +
      j2)) - 0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 +
      j2)) - SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(
      j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3))*
      ZU(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpAhSuconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.35355339059327373)*vT*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 +
      j1)))*ZA(gI2,0),0) + IF(gO2 < 3,std::complex<double>(0.,0.7071067811865475)*
      Mu*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZA(gI2,0),0) + IF(gO2 <
      3,std::complex<double>(0.,0.7071067811865475)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1)
      )*Conj(TYu(j1,gO2)))*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(0.,
      0.35355339059327373)*vd*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 +
      j1)))*ZA(gI2,2),0) - std::complex<double>(0.,0.35355339059327373)*vT*Conj(
      Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*
      Yu(j1,j2)))*ZA(gI2,0) - std::complex<double>(0.,0.7071067811865475)*Conj(Mu)
      *SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)
      ))*ZA(gI2,0) - std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(
      ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYu(j1,j2)))*ZA(gI2,1) -
      std::complex<double>(0.,0.35355339059327373)*vd*Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZA(gI2,2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CphhSuconjUSu(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1
      )*ZH(gI2,0),0) + IF(gO2 < 3,-0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0),0)
      + IF(gO2 < 3,0.35355339059327373*vT*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj
      (ZU(gI1,3 + j1)))*ZH(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*Mu*SUM(j1,0,2
      ,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZH(gI2,0),0) + IF(gO2 < 3,-0.05*vu*
      Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gO2 < 3,0.25*vu*Conj(ZU(gI1,gO2)
      )*Sqr(g2)*ZH(gI2,1),0) + IF(gO2 < 3,-0.7071067811865475*SUM(j1,0,2,Conj(ZU(
      gI1,3 + j1))*Conj(TYu(j1,gO2)))*ZH(gI2,1),0) + IF(gO2 < 3,-(vu*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZH(gI2,1)),0) + IF(
      gO2 < 3,0.35355339059327373*vd*Lambdax*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(
      gI1,3 + j1)))*ZH(gI2,2),0) - 0.2*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) + 0.35355339059327373*vT*Conj(Lambdax)
      *SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)
      ))*ZH(gI2,0) + 0.7071067811865475*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(
      j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZH(gI2,0) + 0.2*vu*Sqr(g1)*SUM
      (j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) -
      0.7071067811865475*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2
      ,3 + j1)*TYu(j1,j2)))*ZH(gI2,1) - vu*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,
      0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZH(
      gI2,1) + 0.35355339059327373*vd*Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*
      SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZH(gI2,2);

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
      ThetaW()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Sin(
      ThetaW()),0) - 0.5163977794943222*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7745966692414834*g1*g2*Cos(
      ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW()),0) + IF(gO1 < 3,0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW())),0) + IF(gO1 < 3,0.3*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW())),0) + 1.2*Sqr(g1)*Sqr(Sin(
      ThetaW()))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))
      ;

   return result;
}

double CLASSNAME::CpUSeconjUSeconjVWmVWm(int gO1, int gO2) const
{
   
   const double result = IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(
      j1,gO2))*Ye(j1,gO1))*ZA(gI1,0)*ZA(gI2,0)),0),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-
      0.35355339059327373*Conj(Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(
      j1,gO1))*ZA(gI1,2)*ZA(gI2,1),0) + IF(gO1 < 3,-0.35355339059327373*Conj(
      Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1))*ZA(gI1,1)*ZA(gI2,
      2),0) + IF(gO2 < 3,-0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*
      KroneckerDelta(gO1,3 + j1))*ZA(gI1,2)*ZA(gI2,1),0) + IF(gO2 < 3,-
      0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3
       + j1))*ZA(gI1,1)*ZA(gI2,2),0) + 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,0) - SUM(j3,0,2,
      KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,
      Conj(Ye(j3,j1))*Ye(j2,j1))))*ZA(gI1,0)*ZA(gI2,0) - 0.3*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,1)*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(
      j1,gO2))*Ye(j1,gO1))*ZH(gI1,0)*ZH(gI2,0)),0),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,
      0.35355339059327373*Conj(Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(
      j1,gO1))*ZH(gI1,2)*ZH(gI2,1),0) + IF(gO1 < 3,0.35355339059327373*Conj(
      Lambdax)*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1))*ZH(gI1,1)*ZH(gI2,
      2),0) + IF(gO2 < 3,0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*
      KroneckerDelta(gO1,3 + j1))*ZH(gI1,2)*ZH(gI2,1),0) + IF(gO2 < 3,
      0.35355339059327373*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3
       + j1))*ZH(gI1,1)*ZH(gI2,2),0) + 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
       + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,0) - SUM(j3,0,2,
      KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,
      Conj(Ye(j3,j1))*Ye(j2,j1))))*ZH(gI1,0)*ZH(gI2,0) - 0.3*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,1)*ZH(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUSeSvconjUSeconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZV(gI1,gO2)
      )*Sqr(g2)*ZV(gI2,gO1),0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gI1,gI2)
      *KroneckerDelta(gO1,gO2)*Sqr(g2),0) + 0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - SUM(j2,0
      ,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*SUM(j4
      ,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZV(gI2,j4));

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
      KroneckerDelta(gI1,gO1)),0) + Conj(UM(gI2,1))*SUM(j1,0,2,KroneckerDelta(gO1,
      3 + j1)*Ye(j1,gI1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmSvconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(ZV(
      gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,0,2
      ,Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,j2)))*ZP(gI2,0),0) + IF(
      gO2 < 3,-0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1),0) + IF(
      gO2 < 3,-0.5*vT*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,2),0) + IF(gO2 < 3,0.5*vT*
      Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,3),0) + SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,
      0,2,KroneckerDelta(gO2,3 + j1)*TYe(j1,j2)))*ZP(gI2,0) - 0.5*vT*Conj(Lambdax)
      *SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)
      ))*ZP(gI2,1) + Conj(Mu)*SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZP(gI2,1) - 0.7071067811865475*vu*
      Conj(Lambdax)*SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 +
      j1)*Ye(j1,j2)))*ZP(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePR(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(j1,gO2))*
      ZER(gI1,j1))*ZN(gI2,2)),0) - 1.0954451150103321*g1*SUM(j1,0,2,KroneckerDelta
      (gO2,3 + j1)*ZER(gI1,j1))*ZN(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePL(int gI2, int gI1, int gO1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZEL(
      gI1,gO1))*Conj(ZN(gI2,0)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZEL(gI1
      ,gO1))*Conj(ZN(gI2,1)),0) - Conj(ZN(gI2,2))*SUM(j2,0,2,Conj(ZEL(gI1,j2))*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSeconjHpmconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-(Conj(Lambdax)*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1))*ZP(gI1,1)*ZP(gI2,2)),0) + IF(gO1 < 3,
      -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2),0) + IF(gO1 < 3,0.5
      *KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3),0) + IF(gO2 < 3,-(
      Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1))*ZP(gI1,2)*ZP
      (gI2,1)),0) + 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) - SUM(j3,0,2,KroneckerDelta(
      gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))
      *Ye(j2,j1))))*ZP(gI1,0)*ZP(gI2,0) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSdUSeconjSdconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) +
      IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 +
      j1))*ZD(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2)),0) +
      IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI1,3 +
      j2))*ZD(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)
      *Ye(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gI1,3 + j3)))*ZD(
      gI2,j4))),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Yd(j1,j2)
      *ZD(gI2,3 + j1)))*SUM(j3,0,2,Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3))),0
      ) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZD(gI1,j2))*ZD(gI2,j2))
      - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZD(gI1,3 + j2))*ZD(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSeconjSeconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Ye(j1,
      gO1)*ZE(gI2,3 + j1))*SUM(j3,0,2,Conj(Ye(j3,gO2))*Conj(ZE(gI1,3 + j3)))),0),0
      ) + IF(gO1 < 3,IF(gO2 < 3,-0.15*Conj(ZE(gI1,gO2))*Sqr(g1)*ZE(gI2,gO1),0),0)
      + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZE(gI1,gO2))*Sqr(g2)*ZE(gI2,gO1),0),0) +
      IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1)
      )*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,
      0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2
      )*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,-
      0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)
      ),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1,0,2
      ,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4)
      )*Conj(ZE(gI1,3 + j3)))*ZE(gI2,j4))),0) + IF(gO1 < 3,0.15*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZE(gI2,gO1),0) + IF(gO1 <
      3,0.15*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*KroneckerDelta(gO2,3 + j2))*
      ZE(gI2,gO1),0) + IF(gO2 < 3,0.15*Conj(ZE(gI1,gO2))*Sqr(g1)*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1)),0) + IF(gO2 < 3,0.15*Conj(ZE(gI1,
      gO2))*Sqr(g1)*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)),0) + IF(
      gO2 < 3,-(SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1)))*
      SUM(j3,0,2,Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3))),0) - 0.3*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZE(gI2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1
      ,3 + j2))*KroneckerDelta(gO2,3 + j2)) + 0.15*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*ZE(gI2,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3
      + j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,
      0,2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.15*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2
      ,Conj(ZE(gI1,j2))*ZE(gI2,j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3
      + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 +
      j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1
      ))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)) - SUM(j2,0,2,Conj(
      ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*SUM(j4,0,2,SUM
      (j3,0,2,Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZE(gI2,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSuconjUSeconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 +
      j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta
      (gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,-
      0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3
       + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,2
      ,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(j1
      ,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2))
      + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpAhSeconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2)))*ZA(
      gI2,0),0) + IF(gO2 < 3,std::complex<double>(0.,0.35355339059327373)*vT*
      Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1)))*ZA(gI2,1),0) + IF(
      gO2 < 3,std::complex<double>(0.,0.7071067811865475)*Mu*SUM(j1,0,2,Conj(Ye(j1
      ,gO2))*Conj(ZE(gI1,3 + j1)))*ZA(gI2,1),0) + IF(gO2 < 3,std::complex<double>(
      0.,0.35355339059327373)*vu*Lambdax*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3
       + j1)))*ZA(gI2,2),0) - std::complex<double>(0.,0.7071067811865475)*SUM(j2,0
      ,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYe(j1,j2)))*ZA(
      gI2,0) - std::complex<double>(0.,0.35355339059327373)*vT*Conj(Lambdax)*SUM(
      j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZA
      (gI2,1) - std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZA(gI2,1)
      - std::complex<double>(0.,0.35355339059327373)*vu*Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZA(gI2,2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CphhSeconjUSe(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.15*vd*Conj(ZE(gI1,gO2))*Sqr(
      g1)*ZH(gI2,0),0) + IF(gO2 < 3,0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0),0)
      + IF(gO2 < 3,-0.7071067811865475*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(TYe(j1
      ,gO2)))*ZH(gI2,0),0) + IF(gO2 < 3,-(vd*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,
      2,Conj(Ye(j1,gO2))*Ye(j1,j2)))*ZH(gI2,0)),0) + IF(gO2 < 3,0.15*vu*Conj(ZE(
      gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gO2 < 3,-0.25*vu*Conj(ZE(gI1,gO2))*Sqr(
      g2)*ZH(gI2,1),0) + IF(gO2 < 3,0.35355339059327373*vT*Lambdax*SUM(j1,0,2,Conj
      (Ye(j1,gO2))*Conj(ZE(gI1,3 + j1)))*ZH(gI2,1),0) + IF(gO2 < 3,
      0.7071067811865475*Mu*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1)))*ZH(
      gI2,1),0) + IF(gO2 < 3,0.35355339059327373*vu*Lambdax*SUM(j1,0,2,Conj(Ye(j1,
      gO2))*Conj(ZE(gI1,3 + j1)))*ZH(gI2,2),0) + 0.3*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZE
      (gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) - 0.7071067811865475*SUM
      (j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*TYe(j1,j2)))*
      ZH(gI2,0) - vd*SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2
      ,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))))*ZH(gI2,0) - 0.3*vu*Sqr(g1)*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) +
      0.35355339059327373*vT*Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZH(gI2,1) + 0.7071067811865475*Conj(
      Mu)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,
      j2)))*ZH(gI2,1) + 0.35355339059327373*vu*Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZH(gI2,2);

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
      ThetaW()),0) + IF(gO2 < 3,0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(ThetaW
      ()),0) + 0.7745966692414834*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))
      *KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgWmUhh(int gO1) const
{
   
   const std::complex<double> result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1) + 4*vT*KroneckerDelta(2,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgWmCUhh(int gO1) const
{
   
   const std::complex<double> result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1) + 4*vT*KroneckerDelta(2,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargZgZUhh(int gO1) const
{
   
   const std::complex<double> result = -0.025*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1)
      + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(int gO2) const
{
   
   const std::complex<double> result = 0.05*(vd*KroneckerDelta(0,gO2) + vu*
      KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1)
      + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(int gO2) const
{
   
   const std::complex<double> result = 0.5*(vd*KroneckerDelta(0,gO2) + vu*
      KroneckerDelta(1,gO2) + 4*vT*KroneckerDelta(2,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*
      g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5
      *Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) + 4*KroneckerDelta(2,gO1
      )*KroneckerDelta(2,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,
      gO1)*KroneckerDelta(2,gO2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(
      gI2,1) - 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2)) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (10*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 10*AbsSqr(
      Lambdax)*ZA(gI1,2)*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,
      gO1)*(KroneckerDelta(2,gO2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) +
      KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) +
      KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 5*Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) - 10*AbsSqr(
      Lambdax)*KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) +
      KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,
      0)*ZH(gI2,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) - 10*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2))) - KroneckerDelta(0,gO1)*(-(KroneckerDelta(1,
      gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH
      (gI1,0)*ZH(gI2,1))) + 10*AbsSqr(Lambdax)*KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH
      (gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + KroneckerDelta(0,gO2)*(3*(3*Sqr(g1) + 5*Sqr
      (g2))*ZH(gI1,0)*ZH(gI2,0) + (10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(
      gI1,1)*ZH(gI2,1) + 10*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,
      gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-5*KroneckerDelta(2,gO2)*((Conj(
      TLambdax) + 2*Conj(MT)*Lambdax + TLambdax)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*
      ZA(gI2,1)) + 2*Conj(Mu)*Lambdax*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1))
      + 2*Conj(Lambdax)*(ZA(gI1,0)*((vT*Lambdax + Mu)*ZA(gI2,0) + MT*ZA(gI2,1)) +
      ZA(gI1,1)*(MT*ZA(gI2,0) + (vT*Lambdax + Mu)*ZA(gI2,1)))) + KroneckerDelta(1,
      gO2)*(10*MT*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,0) - 5*Conj(TLambdax)*ZA(gI1,2)*
      ZA(gI2,0) + 10*Conj(MT)*Lambdax*ZA(gI1,2)*ZA(gI2,0) - 5*TLambdax*ZA(gI1,2)*
      ZA(gI2,0) - 3*vu*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1) - 5*vu*Sqr(g2)*ZA(gI1,1)*ZA(gI2
      ,1) - 10*vu*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2) + ZA(gI1,0)*(vu*(-10*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI2,0) + 5*(2*MT*Conj(Lambdax) - Conj(
      TLambdax) + 2*Conj(MT)*Lambdax - TLambdax)*ZA(gI2,2))) - KroneckerDelta(0,
      gO2)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*(vd*(10*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI2,1) + 5*(-2*MT*Conj(Lambdax)
      + Conj(TLambdax) - 2*Conj(MT)*Lambdax + TLambdax)*ZA(gI2,2)) + 5*ZA(gI1,2)*(
      -2*Conj(MT)*Lambdax*ZA(gI2,1) + (Conj(TLambdax) + TLambdax)*ZA(gI2,1) + Conj
      (Lambdax)*(-2*MT*ZA(gI2,1) + 2*vd*Lambdax*ZA(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhUhh(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(-(Conj(
      TLambdax)*KroneckerDelta(2,gO2)*ZA(gI2,1)*ZH(gI1,0)) + KroneckerDelta(2,gO2)
      *TLambdax*ZA(gI2,1)*ZH(gI1,0) - Conj(TLambdax)*KroneckerDelta(1,gO2)*ZA(gI2,
      2)*ZH(gI1,0) - 2*Conj(Mu)*KroneckerDelta(0,gO2)*Lambdax*ZA(gI2,2)*ZH(gI1,0)
      + KroneckerDelta(1,gO2)*TLambdax*ZA(gI2,2)*ZH(gI1,0) - Conj(TLambdax)*
      KroneckerDelta(2,gO2)*ZA(gI2,0)*ZH(gI1,1) + KroneckerDelta(2,gO2)*TLambdax*
      ZA(gI2,0)*ZH(gI1,1) - Conj(TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZH(gI1,
      1) - 2*Conj(Mu)*KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,2)*ZH(gI1,1) +
      KroneckerDelta(0,gO2)*TLambdax*ZA(gI2,2)*ZH(gI1,1) - Conj(TLambdax)*
      KroneckerDelta(1,gO2)*ZA(gI2,0)*ZH(gI1,2) + KroneckerDelta(1,gO2)*TLambdax*
      ZA(gI2,0)*ZH(gI1,2) - Conj(TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,1)*ZH(gI1,
      2) + KroneckerDelta(0,gO2)*TLambdax*ZA(gI2,1)*ZH(gI1,2) - 2*Conj(Lambdax)*(-
      (KroneckerDelta(0,gO2)*ZA(gI2,2)*(Mu*ZH(gI1,0) + MT*ZH(gI1,1))) -
      KroneckerDelta(1,gO2)*ZA(gI2,2)*(MT*ZH(gI1,0) + Mu*ZH(gI1,1)) + MT*
      KroneckerDelta(2,gO2)*(ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + MT*
      KroneckerDelta(1,gO2)*ZA(gI2,0)*ZH(gI1,2) + MT*KroneckerDelta(0,gO2)*ZA(gI2,
      1)*ZH(gI1,2)) + 2*Conj(MT)*Lambdax*(KroneckerDelta(2,gO2)*(ZA(gI2,1)*ZH(gI1,
      0) + ZA(gI2,0)*ZH(gI1,1)) + KroneckerDelta(1,gO2)*(-(ZA(gI2,2)*ZH(gI1,0)) +
      ZA(gI2,0)*ZH(gI1,2)) + KroneckerDelta(0,gO2)*(-(ZA(gI2,2)*ZH(gI1,1)) + ZA(
      gI2,1)*ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO2)*(-(ZH(gI1,0)*(3
      *vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,0) + vu*(10*AbsSqr(Lambdax) - 3*Sqr(g1) -
      5*Sqr(g2))*ZH(gI2,1) + 10*(Conj(Mu)*Lambdax + Conj(Lambdax)*(vT*Lambdax + Mu
      ))*ZH(gI2,2))) + ZH(gI1,1)*(vu*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))
      *ZH(gI2,0) + vd*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) + 5*
      (2*MT*Conj(Lambdax) + Conj(TLambdax) + 2*Conj(MT)*Lambdax + TLambdax)*ZH(gI2
      ,2)) - 5*ZH(gI1,2)*(2*Conj(Mu)*Lambdax*ZH(gI2,0) - (Conj(TLambdax) + 2*Conj(
      MT)*Lambdax + TLambdax)*ZH(gI2,1) + 2*Conj(Lambdax)*((vT*Lambdax + Mu)*ZH(
      gI2,0) - MT*ZH(gI2,1) + vd*Lambdax*ZH(gI2,2)))) + KroneckerDelta(1,gO2)*(ZH(
      gI1,1)*(vd*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,0) - 3*vu*(3
      *Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) - 10*(Conj(Mu)*Lambdax + Conj(Lambdax)*(vT*
      Lambdax + Mu))*ZH(gI2,2)) + ZH(gI1,0)*(vu*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) +
      5*Sqr(g2))*ZH(gI2,0) + vd*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(
      gI2,1) + 5*(2*MT*Conj(Lambdax) + Conj(TLambdax) + 2*Conj(MT)*Lambdax +
      TLambdax)*ZH(gI2,2)) + 5*ZH(gI1,2)*(Conj(TLambdax)*ZH(gI2,0) + 2*Conj(MT)*
      Lambdax*ZH(gI2,0) + TLambdax*ZH(gI2,0) - 2*Conj(Mu)*Lambdax*ZH(gI2,1) + 2*
      Conj(Lambdax)*(MT*ZH(gI2,0) - (vT*Lambdax + Mu)*ZH(gI2,1) - vu*Lambdax*ZH(
      gI2,2)))) - 5*KroneckerDelta(2,gO2)*(-((Conj(TLambdax) + 2*Conj(MT)*Lambdax
      + TLambdax)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1))) + 2*Conj(Mu)*
      Lambdax*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) + 2*Conj(Lambdax)*(
      Lambdax*ZH(gI1,2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1)) + ZH(gI1,0)*((vT*Lambdax +
      Mu)*ZH(gI2,0) - MT*ZH(gI2,1) + vd*Lambdax*ZH(gI2,2)) + ZH(gI1,1)*(-(MT*ZH(
      gI2,0)) + (vT*Lambdax + Mu)*ZH(gI2,1) + vu*Lambdax*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSvconjSv(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.05*(vd*KroneckerDelta(0,gO2) - vu*
      KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.5*(-1.4142135623730951*KroneckerDelta(0,
      gO2)*(g2*UM(gI1,1)*UP(gI2,0) + Conj(Lambdax)*UM(gI1,2)*UP(gI2,1)) +
      KroneckerDelta(2,gO2)*(-2*g2*UM(gI1,2)*UP(gI2,0) + Conj(Lambdax)*UM(gI1,1)*
      UP(gI2,1) + 2*g2*UM(gI1,0)*UP(gI2,2)) + 1.4142135623730951*KroneckerDelta(1,
      gO2)*(-(g2*UM(gI1,0)*UP(gI2,1)) + Conj(Lambdax)*UM(gI1,1)*UP(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.5*(g2*Conj(UM(gI2,0))*(-
      1.4142135623730951*Conj(UP(gI1,1))*KroneckerDelta(1,gO1) + 2*Conj(UP(gI1,2))
      *KroneckerDelta(2,gO1)) - Conj(UM(gI2,2))*(2*g2*Conj(UP(gI1,0))*
      KroneckerDelta(2,gO1) + 1.4142135623730951*Conj(UP(gI1,1))*KroneckerDelta(0,
      gO1)*Lambdax) + Conj(UM(gI2,1))*(-1.4142135623730951*g2*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) + 1.4142135623730951*Conj(UP(gI1,2))*KroneckerDelta(1,
      gO1)*Lambdax + Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

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

std::complex<double> CLASSNAME::CpUhhUhhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*(5*(
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(
      gI1,0)*ZP(gI2,1)) + 1.4142135623730951*KroneckerDelta(2,gO2)*(-AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,2)*ZP(gI2,0) - ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(
      ZP(gI2,2) - ZP(gI2,3)))) + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP
      (gI1,0)*ZP(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*
      ZP(gI2,1) - 10*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2) + 10*Sqr(
      g2)*ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(1,gO1)*(-5*(KroneckerDelta(0,gO2
      )*(AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      1.4142135623730951*KroneckerDelta(2,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZP(
      gI1,2)*ZP(gI2,1) - ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) - ZP(gI2,3))))
      + KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZP(
      gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*Sqr(g2)*
      ZP(gI1,2)*ZP(gI2,2) + 10*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3))
      ) - 5*KroneckerDelta(2,gO1)*(Sqr(g2)*(1.4142135623730951*KroneckerDelta(0,
      gO2)*(ZP(gI1,2)*ZP(gI2,0) - ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) - ZP(
      gI2,3))) + 1.4142135623730951*KroneckerDelta(1,gO2)*(ZP(gI1,2)*ZP(gI2,1) -
      ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) - ZP(gI2,3))) + 4*KroneckerDelta(
      2,gO2)*(ZP(gI1,2) - ZP(gI1,3))*(ZP(gI2,2) - ZP(gI2,3))) + AbsSqr(Lambdax)*(2
      *KroneckerDelta(2,gO2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1)) +
      1.4142135623730951*(KroneckerDelta(0,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,3
      )*ZP(gI2,0) + ZP(gI1,0)*(-ZP(gI2,2) + ZP(gI2,3))) + KroneckerDelta(1,gO2)*(-
      (ZP(gI1,2)*ZP(gI2,1)) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(-ZP(gI2,2) + ZP(gI2
      ,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO2)*(5*(
      5.656854249492381*Conj(MT)*Lambdax*ZP(gI1,2)*ZP(gI2,0) + 2.8284271247461903*
      TLambdax*ZP(gI1,3)*ZP(gI2,0) + 1.4142135623730951*vT*AbsSqr(Lambdax)*ZP(gI1,
      2)*ZP(gI2,1) + 2.8284271247461903*Conj(Lambdax)*Mu*ZP(gI1,2)*ZP(gI2,1) -
      1.4142135623730951*vT*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951*vT*
      AbsSqr(Lambdax)*ZP(gI1,3)*ZP(gI2,1) + 2.8284271247461903*Conj(Mu)*Lambdax*ZP
      (gI1,3)*ZP(gI2,1) + 1.4142135623730951*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) - 2*vu
      *Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 4*vu*AbsSqr(Lambdax)*ZP(gI1,3)*ZP(gI2,3) + 2*
      vu*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)) + ZP(gI1,0)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr
      (g1) - 5*Sqr(g2))*ZP(gI2,0) - 5*vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,1) +
      14.142135623730951*(2*MT*Conj(Lambdax)*ZP(gI2,2) + Conj(TLambdax)*ZP(gI2,3))
      ) - ZP(gI1,1)*(5*vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vu*(3*Sqr(g1) +
      5*Sqr(g2))*ZP(gI2,1) + 7.0710678118654755*((-(vT*AbsSqr(Lambdax)) - 2*Conj(
      Mu)*Lambdax + vT*Sqr(g2))*ZP(gI2,2) + (vT*AbsSqr(Lambdax) - 2*Conj(Lambdax)*
      Mu - vT*Sqr(g2))*ZP(gI2,3)))) - KroneckerDelta(0,gO2)*(ZP(gI1,1)*(5*vu*(
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vd*(20*AbsSqr(Lambdax) - 3*Sqr(g1) +
      5*Sqr(g2))*ZP(gI2,1) + 14.142135623730951*(TLambdax*ZP(gI2,2) + 2*Conj(MT)*
      Lambdax*ZP(gI2,3))) + 5*(ZP(gI1,2)*(1.4142135623730951*(Conj(Lambdax)*(-(vT*
      Lambdax) + 2*Mu) + vT*Sqr(g2))*ZP(gI2,0) + 2.8284271247461903*Conj(TLambdax)
      *ZP(gI2,1) - 2*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,2)) + ZP(gI1,3)*(
      1.4142135623730951*(vT*AbsSqr(Lambdax) + 2*Conj(Mu)*Lambdax - vT*Sqr(g2))*ZP
      (gI2,0) + 5.656854249492381*MT*Conj(Lambdax)*ZP(gI2,1) + 2*vd*Sqr(g2)*ZP(gI2
      ,3))) + ZP(gI1,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,0) + 5*vu*(AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gI2,1) + 7.0710678118654755*((-(vT*AbsSqr(Lambdax)) +
      2*Conj(Mu)*Lambdax + vT*Sqr(g2))*ZP(gI2,2) + (Conj(Lambdax)*(vT*Lambdax + 2*
      Mu) - vT*Sqr(g2))*ZP(gI2,3)))) - 5*KroneckerDelta(2,gO2)*(-2*Conj(Mu)*
      Lambdax*ZP(gI1,0)*ZP(gI2,0) - 4*Conj(MT)*Lambdax*ZP(gI1,1)*ZP(gI2,0) - 2*
      TLambdax*ZP(gI1,1)*ZP(gI2,0) + 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,2)*ZP(
      gI2,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) - 2*Conj(TLambdax
      )*ZP(gI1,0)*ZP(gI2,1) - 2*Conj(Mu)*Lambdax*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*vu*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951*vu*
      Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) + 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,0)*ZP(gI2
      ,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) + 4*vT*Sqr(g2)*ZP(
      gI1,2)*ZP(gI2,2) - 4*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) - 1.4142135623730951*vd*
      Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) - 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2
      ,3) - 4*vT*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) + 4*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3) +
      Conj(Lambdax)*(1.4142135623730951*vd*Lambdax*ZP(gI1,3)*ZP(gI2,0) + 2*vT*
      Lambdax*ZP(gI1,1)*ZP(gI2,1) - 2*Mu*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*
      vu*Lambdax*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*Lambdax*ZP(gI1,2)*(vd*ZP
      (gI2,0) + vu*ZP(gI2,1)) - 1.4142135623730951*vu*Lambdax*ZP(gI1,1)*ZP(gI2,2)
      + 1.4142135623730951*vu*Lambdax*ZP(gI1,1)*ZP(gI2,3) + ZP(gI1,0)*(2*(vT*
      Lambdax - Mu)*ZP(gI2,0) - 4*MT*ZP(gI2,1) + 1.4142135623730951*vd*Lambdax*(-
      ZP(gI2,2) + ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUhhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.1*(5*Conj(Lambdax)*KroneckerDelta(2,gO2)*
      (ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,3)) + KroneckerDelta(0,gO2)*(ZN(gI1,
      2)*(3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + 3.872983346207417*g1*
      ZN(gI1,0)*ZN(gI2,2) - 5*g2*ZN(gI1,1)*ZN(gI2,2) + 5*Conj(Lambdax)*ZN(gI1,4)*
      ZN(gI2,3) + 5*Conj(Lambdax)*ZN(gI1,3)*ZN(gI2,4)) + KroneckerDelta(1,gO2)*(ZN
      (gI1,3)*(-3.872983346207417*g1*ZN(gI2,0) + 5*g2*ZN(gI2,1)) + (-
      3.872983346207417*g1*ZN(gI1,0) + 5*g2*ZN(gI1,1))*ZN(gI2,3) + 5*Conj(Lambdax)
      *(ZN(gI1,4)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUhhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.1*(-5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2))*
      KroneckerDelta(0,gO1) - 3.872983346207417*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))
      *KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta
      (1,gO1) + 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) +
      3.872983346207417*g1*Conj(ZN(gI1,0))*(Conj(ZN(gI2,2))*KroneckerDelta(0,gO1)
      - Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)) + 5*Conj(ZN(gI1,4))*Conj(ZN(gI2,3))
      *KroneckerDelta(0,gO1)*Lambdax + 5*Conj(ZN(gI1,3))*Conj(ZN(gI2,4))*
      KroneckerDelta(0,gO1)*Lambdax + 5*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*
      KroneckerDelta(1,gO1)*Lambdax + 5*Conj(ZN(gI1,3))*Conj(ZN(gI2,2))*
      KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(3.872983346207417*g1*Conj(
      ZN(gI2,0))*KroneckerDelta(0,gO1) - 5*g2*Conj(ZN(gI2,1))*KroneckerDelta(0,gO1
      ) + 5*Conj(ZN(gI2,4))*KroneckerDelta(1,gO1)*Lambdax + 5*Conj(ZN(gI2,3))*
      KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.0710678118654755*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1
      ,0,2,Yd(j1,j2)*ZD(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,
      j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2))) - KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) -
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(
      gI1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2)))) + KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1)
      )*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) -
      20*(SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd
      (j2,j1))*ZD(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,
      2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.0710678118654755*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1
      ,0,2,Ye(j1,j2)*ZE(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,
      j2))*Conj(ZE(gI1,3 + j1)))*ZE(gI2,j2))) + KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) +
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*ZE(gI2,j2)))) + KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*ZE(gI2,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))
      - 20*(SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*
      Ye(j2,j1))*ZE(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,
      0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.0710678118654755*KroneckerDelta(0,
      gO2)*KroneckerDelta(2,gO1)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1
      ,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,
      j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,j2))) + KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(
      gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) +
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(
      gI1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,j2)))) - KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) +
      20*(SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu
      (j2,j1))*ZU(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,
      2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSdconjSd(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(7.0710678118654755*vu*KroneckerDelta(
      2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(
      gI1,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3
      + j1)))*ZD(gI1,j2))) + KroneckerDelta(1,gO2)*(-(vu*(Sqr(g1) + 5*Sqr(g2))*SUM
      (j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))) - 2*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,
      3 + j1))*ZD(gI1,3 + j1)) + 7.0710678118654755*(vT*Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZD(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + 2*Conj(Mu)*SUM(j2,0
      ,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + (vT*Lambdax + 2*
      Mu)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2)))
      ) + KroneckerDelta(0,gO2)*(vd*(Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI2,
      j1))*ZD(gI1,j1)) + 2*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 +
      j1)) - 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,ZD(gI1,
      3 + j1)*TYd(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gI2,
      3 + j1))*Conj(TYd(j1,j2)))*ZD(gI1,j2)) + 2*vd*(SUM(j3,0,2,Conj(ZD(gI2,3 + j3
      ))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI1,3 + j2))) + SUM(
      j3,0,2,SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD
      (gI1,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSeconjSe(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(7.0710678118654755*vu*KroneckerDelta(
      2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(
      gI1,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3
      + j1)))*ZE(gI1,j2))) + KroneckerDelta(1,gO2)*(vu*(3*Sqr(g1) - 5*Sqr(g2))*SUM
      (j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) - 6*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI2,3
       + j1))*ZE(gI1,3 + j1)) + 7.0710678118654755*(vT*Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZE(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + 2*Conj(Mu)*SUM(j2,0
      ,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + (vT*Lambdax + 2*
      Mu)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2)))
      ) - KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI2
      ,j1))*ZE(gI1,j1)) - 6*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 +
      j1)) + 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,ZE(gI1,
      3 + j1)*TYe(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gI2,
      3 + j1))*Conj(TYe(j1,j2)))*ZE(gI1,j2)) + 2*vd*(SUM(j3,0,2,Conj(ZE(gI2,3 + j3
      ))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI1,3 + j2))) + SUM(
      j3,0,2,SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE
      (gI1,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSuconjSu(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(7.0710678118654755*vd*KroneckerDelta(
      2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(
      gI1,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3
      + j1)))*ZU(gI1,j2))) + KroneckerDelta(0,gO2)*(vd*(Sqr(g1) - 5*Sqr(g2))*SUM(
      j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3
      + j1))*ZU(gI1,3 + j1)) + 7.0710678118654755*(vT*Conj(Lambdax)*SUM(j2,0,2,
      Conj(ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + 2*Conj(Mu)*SUM(j2,0
      ,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + (vT*Lambdax + 2*
      Mu)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2)))
      ) - KroneckerDelta(1,gO2)*(vu*(Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI2,
      j1))*ZU(gI1,j1)) - 4*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 +
      j1)) + 10*(1.4142135623730951*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,ZU(gI1,
      3 + j1)*TYu(j1,j2))) + 1.4142135623730951*SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gI2,
      3 + j1))*Conj(TYu(j1,j2)))*ZU(gI1,j2)) + 2*vu*(SUM(j3,0,2,Conj(ZU(gI2,3 + j3
      ))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2))) + SUM(
      j3,0,2,SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU
      (gI1,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpAhUhhVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZA(
      gI2,0) - KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjVWm(int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) -
      KroneckerDelta(1,gO2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(2,gO2)*(
      ZP(gI2,2) + ZP(gI2,3)));

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
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*
      g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5
      *Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) + 4*KroneckerDelta(2,gO1
      )*KroneckerDelta(2,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,
      gO1)*(KroneckerDelta(2,gO2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1)) +
      KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) +
      KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 5*Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) - 10*AbsSqr(
      Lambdax)*KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)) +
      KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,
      0)*ZA(gI2,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) - 10*AbsSqr(
      Lambdax)*ZA(gI1,2)*ZA(gI2,2))) - KroneckerDelta(0,gO1)*(-(KroneckerDelta(1,
      gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,1))) + 10*AbsSqr(Lambdax)*KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA
      (gI2,0) + ZA(gI1,0)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(3*(3*Sqr(g1) + 5*Sqr
      (g2))*ZA(gI1,0)*ZA(gI2,0) + (10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZA(
      gI1,1)*ZA(gI2,1) + 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,
      gO1)*KroneckerDelta(2,gO2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(
      gI2,1) - 10*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2)) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (10*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 10*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,
      gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAh(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*((Conj(
      TLambdax) + 2*Conj(MT)*Lambdax - TLambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,1)
      *ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,
      0) + ZA(gI1,0)*ZA(gI2,2)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(
      gI1,1)*ZA(gI2,2))) - 2*Conj(Mu)*Lambdax*(KroneckerDelta(2,gO2)*(ZA(gI1,0)*ZA
      (gI2,0) + ZA(gI1,1)*ZA(gI2,1)) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,0)
      + ZA(gI1,0)*ZA(gI2,2)) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1
      ,1)*ZA(gI2,2))) + 2*Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,0)*(Mu*ZA(
      gI2,0) - MT*ZA(gI2,1)) + ZA(gI1,1)*(-(MT*ZA(gI2,0)) + Mu*ZA(gI2,1))) +
      KroneckerDelta(0,gO2)*(ZA(gI1,2)*(Mu*ZA(gI2,0) - MT*ZA(gI2,1)) + (Mu*ZA(gI1,
      0) - MT*ZA(gI1,1))*ZA(gI2,2)) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*(-(MT*ZA(
      gI2,0)) + Mu*ZA(gI2,1)) + (-(MT*ZA(gI1,0)) + Mu*ZA(gI1,1))*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpAhUAhhh(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.05*(5*KroneckerDelta(2,gO2)*(2*Conj(MT)*
      Lambdax*(ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) - (Conj(TLambdax) +
      TLambdax)*(ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + 2*Conj(Lambdax)*(MT*
      ZA(gI2,1)*ZH(gI1,0) + MT*ZA(gI2,0)*ZH(gI1,1) - Lambdax*ZA(gI2,2)*(vd*ZH(gI1,
      0) + vu*ZH(gI1,1)))) + KroneckerDelta(1,gO2)*(ZA(gI2,1)*(vd*(-10*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(
      gI1,1) - 10*(Conj(Mu)*Lambdax + Conj(Lambdax)*(vT*Lambdax + Mu))*ZH(gI1,2))
      + 5*(2*MT*Conj(Lambdax)*(ZA(gI2,2)*ZH(gI1,0) - ZA(gI2,0)*ZH(gI1,2)) + 2*Conj
      (MT)*Lambdax*(ZA(gI2,2)*ZH(gI1,0) - ZA(gI2,0)*ZH(gI1,2)) - (Conj(TLambdax) +
      TLambdax)*(ZA(gI2,2)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2)))) - KroneckerDelta(0,
      gO2)*(ZA(gI2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) + vu*(10*AbsSqr(
      Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,1) + 10*(Conj(Mu)*Lambdax + Conj(
      Lambdax)*(vT*Lambdax + Mu))*ZH(gI1,2)) + 5*((Conj(TLambdax) + TLambdax)*(ZA(
      gI2,2)*ZH(gI1,1) + ZA(gI2,1)*ZH(gI1,2)) + Conj(Lambdax)*(-2*MT*ZA(gI2,2)*ZH(
      gI1,1) + 2*MT*ZA(gI2,1)*ZH(gI1,2)) + Conj(MT)*(-2*Lambdax*ZA(gI2,2)*ZH(gI1,1
      ) + 2*Lambdax*ZA(gI2,1)*ZH(gI1,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(-(Conj(
      TLambdax)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,0)) - 2*Conj(MT)*
      KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1)*ZH(gI2,0) + KroneckerDelta(2,gO2)*
      TLambdax*ZH(gI1,1)*ZH(gI2,0) - Conj(TLambdax)*KroneckerDelta(1,gO2)*ZH(gI1,2
      )*ZH(gI2,0) + 2*Conj(MT)*KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,2)*ZH(gI2,0) +
      KroneckerDelta(1,gO2)*TLambdax*ZH(gI1,2)*ZH(gI2,0) - Conj(TLambdax)*
      KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,1) - 2*Conj(MT)*KroneckerDelta(2,gO2)
      *Lambdax*ZH(gI1,0)*ZH(gI2,1) + KroneckerDelta(2,gO2)*TLambdax*ZH(gI1,0)*ZH(
      gI2,1) - Conj(TLambdax)*KroneckerDelta(0,gO2)*ZH(gI1,2)*ZH(gI2,1) + 2*Conj(
      MT)*KroneckerDelta(0,gO2)*Lambdax*ZH(gI1,2)*ZH(gI2,1) + KroneckerDelta(0,gO2
      )*TLambdax*ZH(gI1,2)*ZH(gI2,1) - 2*Conj(Mu)*KroneckerDelta(2,gO2)*Lambdax*(
      ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) - Conj(TLambdax)*KroneckerDelta(1
      ,gO2)*ZH(gI1,0)*ZH(gI2,2) + 2*Conj(MT)*KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,
      0)*ZH(gI2,2) + KroneckerDelta(1,gO2)*TLambdax*ZH(gI1,0)*ZH(gI2,2) - Conj(
      TLambdax)*KroneckerDelta(0,gO2)*ZH(gI1,1)*ZH(gI2,2) + 2*Conj(MT)*
      KroneckerDelta(0,gO2)*Lambdax*ZH(gI1,1)*ZH(gI2,2) + KroneckerDelta(0,gO2)*
      TLambdax*ZH(gI1,1)*ZH(gI2,2) + 2*Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZH(
      gI1,0)*(Mu*ZH(gI2,0) + MT*ZH(gI2,1)) + ZH(gI1,1)*(MT*ZH(gI2,0) + Mu*ZH(gI2,1
      ))) - MT*(KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2))
      + KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(
      1.4142135623730951*KroneckerDelta(0,gO2)*(g2*UM(gI1,1)*UP(gI2,0) - Conj(
      Lambdax)*UM(gI1,2)*UP(gI2,1)) + KroneckerDelta(2,gO2)*(2*g2*UM(gI1,2)*UP(gI2
      ,0) + Conj(Lambdax)*UM(gI1,1)*UP(gI2,1) - 2*g2*UM(gI1,0)*UP(gI2,2)) +
      1.4142135623730951*KroneckerDelta(1,gO2)*(g2*UM(gI1,0)*UP(gI2,1) + Conj(
      Lambdax)*UM(gI1,1)*UP(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Conj(UM(gI2
      ,0))*(1.4142135623730951*Conj(UP(gI1,1))*KroneckerDelta(1,gO1) - 2*Conj(UP(
      gI1,2))*KroneckerDelta(2,gO1)) + Conj(UM(gI2,2))*(2*g2*Conj(UP(gI1,0))*
      KroneckerDelta(2,gO1) - 1.4142135623730951*Conj(UP(gI1,1))*KroneckerDelta(0,
      gO1)*Lambdax) + Conj(UM(gI2,1))*(1.4142135623730951*g2*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) + 1.4142135623730951*Conj(UP(gI1,2))*KroneckerDelta(1,
      gO1)*Lambdax + Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

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

std::complex<double> CLASSNAME::CpUAhUAhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*(5*
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(
      gI1,0)*ZP(gI2,1)) - KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)
      *ZP(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1
      ) - 10*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2) + 10*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,3)) - 7.0710678118654755*KroneckerDelta(2,gO2)*(-AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(
      ZP(gI2,2) + ZP(gI2,3)))) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-
      20*AbsSqr(Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1)
      + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) + 10*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3)) + 5*(KroneckerDelta(0,gO2)*(
      AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      1.4142135623730951*KroneckerDelta(2,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZP(
      gI1,2)*ZP(gI2,1) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3))))
      ) - 5*KroneckerDelta(2,gO1)*(Sqr(g2)*(4*KroneckerDelta(2,gO2)*(ZP(gI1,2) +
      ZP(gI1,3))*(ZP(gI2,2) + ZP(gI2,3)) + 1.4142135623730951*KroneckerDelta(0,gO2
      )*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2
      ,3))) - 1.4142135623730951*KroneckerDelta(1,gO2)*(ZP(gI1,2)*ZP(gI2,1) + ZP(
      gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3)))) + AbsSqr(Lambdax)*(2*
      KroneckerDelta(2,gO2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1)) +
      1.4142135623730951*(-(KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,3)
      *ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2,3)))) + KroneckerDelta(1,gO2)*(ZP
      (gI1,2)*ZP(gI2,1) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3)))
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.25)*(4*Conj(MT)*
      KroneckerDelta(2,gO2)*Lambdax*ZP(gI1,1)*ZP(gI2,0) - vu*KroneckerDelta(0,gO2)
      *Sqr(g2)*ZP(gI1,1)*ZP(gI2,0) - vd*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP
      (gI2,0) - 2*KroneckerDelta(2,gO2)*TLambdax*ZP(gI1,1)*ZP(gI2,0) -
      5.656854249492381*Conj(MT)*KroneckerDelta(1,gO2)*Lambdax*ZP(gI1,2)*ZP(gI2,0)
      - 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) +
      1.4142135623730951*vd*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) +
      1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) +
      1.4142135623730951*vd*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) -
      2.8284271247461903*KroneckerDelta(1,gO2)*TLambdax*ZP(gI1,3)*ZP(gI2,0) + 2*
      Conj(TLambdax)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,1) + vu*KroneckerDelta
      (0,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,1) + vd*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,1) - 2.8284271247461903*Conj(TLambdax)*KroneckerDelta(0,gO2)*
      ZP(gI1,2)*ZP(gI2,1) + 1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP
      (gI1,2)*ZP(gI2,1) + 1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(
      gI1,2)*ZP(gI2,1) - 1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,1) + 1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,1) + 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,2) - 1.4142135623730951*vd*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,2) - 1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(
      gI1,1)*ZP(gI2,2) - 1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(
      gI1,1)*ZP(gI2,2) + 2.8284271247461903*KroneckerDelta(0,gO2)*TLambdax*ZP(gI1,
      1)*ZP(gI2,2) + 4*vT*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) - 2*
      Conj(Mu)*Lambdax*(KroneckerDelta(2,gO2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(
      gI2,1)) + 1.4142135623730951*(KroneckerDelta(0,gO2)*(ZP(gI1,3)*ZP(gI2,0) -
      ZP(gI1,0)*ZP(gI2,2)) + KroneckerDelta(1,gO2)*(ZP(gI1,3)*ZP(gI2,1) - ZP(gI1,1
      )*ZP(gI2,2)))) + 2.8284271247461903*Conj(TLambdax)*KroneckerDelta(1,gO2)*ZP(
      gI1,0)*ZP(gI2,3) - 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,3) - 1.4142135623730951*vd*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,3) + 5.656854249492381*Conj(MT)*KroneckerDelta(0,gO2)*Lambdax*
      ZP(gI1,1)*ZP(gI2,3) + 1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP
      (gI1,1)*ZP(gI2,3) - 1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(
      gI1,1)*ZP(gI2,3) - 4*vT*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) +
      Conj(Lambdax)*(KroneckerDelta(0,gO2)*(-(vu*Lambdax*ZP(gI1,1)*ZP(gI2,0)) +
      1.4142135623730951*(vT*Lambdax - 2*Mu)*ZP(gI1,2)*ZP(gI2,0) -
      1.4142135623730951*vT*Lambdax*ZP(gI1,3)*ZP(gI2,0) + vu*Lambdax*ZP(gI1,0)*ZP(
      gI2,1) - 5.656854249492381*MT*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*vT*
      Lambdax*ZP(gI1,0)*ZP(gI2,2) + 1.4142135623730951*vT*Lambdax*ZP(gI1,0)*ZP(gI2
      ,3) + 2.8284271247461903*Mu*ZP(gI1,0)*ZP(gI2,3)) + KroneckerDelta(2,gO2)*(-
      1.4142135623730951*vd*Lambdax*ZP(gI1,3)*ZP(gI2,0) + 2*Mu*ZP(gI1,1)*ZP(gI2,1)
      - 1.4142135623730951*vu*Lambdax*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*
      Lambdax*ZP(gI1,2)*(vd*ZP(gI2,0) + vu*ZP(gI2,1)) + 1.4142135623730951*vu*
      Lambdax*ZP(gI1,1)*ZP(gI2,2) + 1.4142135623730951*vu*Lambdax*ZP(gI1,1)*ZP(gI2
      ,3) + ZP(gI1,0)*(2*Mu*ZP(gI2,0) - 4*MT*ZP(gI2,1) + 1.4142135623730951*vd*
      Lambdax*(ZP(gI2,2) + ZP(gI2,3)))) + KroneckerDelta(1,gO2)*(
      1.4142135623730951*(-((vT*Lambdax + 2*Mu)*ZP(gI1,2)) + vT*Lambdax*ZP(gI1,3))
      *ZP(gI2,1) + ZP(gI1,0)*(vd*Lambdax*ZP(gI2,1) + 5.656854249492381*MT*ZP(gI2,2
      )) + ZP(gI1,1)*(-(vd*Lambdax*ZP(gI2,0)) + 1.4142135623730951*(vT*Lambdax*ZP(
      gI2,2) + (-(vT*Lambdax) + 2*Mu)*ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUAhPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(-5*Conj(
      Lambdax)*KroneckerDelta(2,gO2)*(ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,3)) -
      KroneckerDelta(1,gO2)*(ZN(gI1,3)*(3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(
      gI2,1)) + (3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3) + 5*
      Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,4))) + KroneckerDelta(
      0,gO2)*(ZN(gI1,2)*(3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) +
      3.872983346207417*g1*ZN(gI1,0)*ZN(gI2,2) - 5*(g2*ZN(gI1,1)*ZN(gI2,2) + Conj(
      Lambdax)*ZN(gI1,4)*ZN(gI2,3) + Conj(Lambdax)*ZN(gI1,3)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiUAhPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Conj(ZN(
      gI1,1))*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) + 3.872983346207417*g1*Conj(ZN
      (gI1,3))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,3))*Conj(
      ZN(gI2,1))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,0))*(-(Conj(ZN(gI2,
      2))*KroneckerDelta(0,gO1)) + Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)) + 5*Conj
      (ZN(gI1,4))*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1)*Lambdax + 5*Conj(ZN(gI1,3)
      )*Conj(ZN(gI2,4))*KroneckerDelta(0,gO1)*Lambdax + 5*Conj(ZN(gI1,4))*Conj(ZN(
      gI2,2))*KroneckerDelta(1,gO1)*Lambdax + 5*Conj(ZN(gI1,3))*Conj(ZN(gI2,2))*
      KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(-3.872983346207417*g1*Conj(
      ZN(gI2,0))*KroneckerDelta(0,gO1) + 5*(g2*Conj(ZN(gI2,1))*KroneckerDelta(0,
      gO1) + Conj(ZN(gI2,4))*KroneckerDelta(1,gO1)*Lambdax + Conj(ZN(gI2,3))*
      KroneckerDelta(2,gO1)*Lambdax)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-7.0710678118654755*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1
      ,0,2,Yd(j1,j2)*ZD(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,
      j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2))) - KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1))) +
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(
      gI1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2)))) + KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1)
      )*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) -
      20*(SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd
      (j2,j1))*ZD(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,
      2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-7.0710678118654755*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1
      ,0,2,Ye(j1,j2)*ZE(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,
      j2))*Conj(ZE(gI1,3 + j1)))*ZE(gI2,j2))) + KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) -
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZE(
      gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*ZE(gI2,j2)))) + KroneckerDelta(0,
      gO1)*KroneckerDelta(0,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*ZE(gI2,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))
      - 20*(SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*
      Ye(j2,j1))*ZE(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,
      0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(-7.0710678118654755*KroneckerDelta(0,
      gO2)*KroneckerDelta(2,gO1)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1
      ,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,
      j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,j2))) + KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(
      gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) -
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(Lambdax)*SUM(j2,0,2,Conj(ZU(
      gI1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))) + Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,j2)))) - KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) +
      20*(SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu
      (j2,j1))*ZU(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,
      2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSdconjSd(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(2*Conj(Mu)*KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZD(
      gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + Conj(Lambdax)*(vT*
      KroneckerDelta(1,gO2) + vu*KroneckerDelta(2,gO2))*SUM(j2,0,2,Conj(ZD(gI2,j2)
      )*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + 2*KroneckerDelta(0,gO2)*SUM(j2,0,2
      ,Conj(ZD(gI2,j2))*SUM(j1,0,2,ZD(gI1,3 + j1)*TYd(j1,j2))) - vT*KroneckerDelta
      (1,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*
      ZD(gI1,j2)) - vu*KroneckerDelta(2,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd
      (j1,j2))*Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2)) - 2*KroneckerDelta(1,gO2)*Mu*SUM(
      j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2)) - 2*
      KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(TYd(j1
      ,j2)))*ZD(gI1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSeconjSe(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(2*Conj(Mu)*KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZE(
      gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + Conj(Lambdax)*(vT*
      KroneckerDelta(1,gO2) + vu*KroneckerDelta(2,gO2))*SUM(j2,0,2,Conj(ZE(gI2,j2)
      )*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + 2*KroneckerDelta(0,gO2)*SUM(j2,0,2
      ,Conj(ZE(gI2,j2))*SUM(j1,0,2,ZE(gI1,3 + j1)*TYe(j1,j2))) - vT*KroneckerDelta
      (1,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*
      ZE(gI1,j2)) - vu*KroneckerDelta(2,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye
      (j1,j2))*Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2)) - 2*KroneckerDelta(1,gO2)*Mu*SUM(
      j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2)) - 2*
      KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1
      ,j2)))*ZE(gI1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSuconjSu(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-
      0.35355339059327373)*(2*Conj(Mu)*KroneckerDelta(0,gO2)*SUM(j2,0,2,Conj(ZU(
      gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + Conj(Lambdax)*(vT*
      KroneckerDelta(0,gO2) + vd*KroneckerDelta(2,gO2))*SUM(j2,0,2,Conj(ZU(gI2,j2)
      )*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + 2*KroneckerDelta(1,gO2)*SUM(j2,0,2
      ,Conj(ZU(gI2,j2))*SUM(j1,0,2,ZU(gI1,3 + j1)*TYu(j1,j2))) - vT*KroneckerDelta
      (0,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1)))*
      ZU(gI1,j2)) - vd*KroneckerDelta(2,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu
      (j1,j2))*Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2)) - 2*KroneckerDelta(0,gO2)*Mu*SUM(
      j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2)) - 2*
      KroneckerDelta(1,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(TYu(j1
      ,j2)))*ZU(gI1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhVZ(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZH(
      gI2,0) - KroneckerDelta(1,gO2)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHpmconjVWm(int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZP(gI2,0) + KroneckerDelta(1,gO2)*ZP(gI2,1) +
      1.4142135623730951*KroneckerDelta(2,gO2)*(ZP(gI2,2) - ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgZUHpm(int gO2) const
{
   
   const std::complex<double> result = 0.05*g2*(14.142135623730951*g2*vT*Cos(
      ThetaW())*(KroneckerDelta(2,gO2) + KroneckerDelta(3,gO2)) + vd*
      KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      ))) + KroneckerDelta(1,gO2)*(-5*g2*vu*Cos(ThetaW()) + 3.872983346207417*g1*
      vu*Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmconjUHpm(int gO1) const
{
   
   const std::complex<double> result = -0.05*g2*(vd*KroneckerDelta(0,gO1) - vu*
      KroneckerDelta(1,gO1))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgZconjUHpm(int gO1) const
{
   
   const std::complex<double> result = 0.05*g2*(14.142135623730951*g2*vT*Cos(
      ThetaW())*(KroneckerDelta(2,gO1) + KroneckerDelta(3,gO1)) + vd*
      KroneckerDelta(0,gO1)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      ))) + KroneckerDelta(1,gO1)*(-5*g2*vu*Cos(ThetaW()) + 3.872983346207417*g1*
      vu*Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpbargZgWmCUHpm(int gO2) const
{
   
   const std::complex<double> result = -0.05*g2*(vd*KroneckerDelta(0,gO2) - vu*
      KroneckerDelta(1,gO2))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVPVWm(int gO2) const
{
   
   const std::complex<double> result = -0.1*g2*(3.872983346207417*g1*vd*Cos(ThetaW
      ())*KroneckerDelta(0,gO2) - 3.872983346207417*g1*vu*Cos(ThetaW())*
      KroneckerDelta(1,gO2) + 7.0710678118654755*g2*vT*(KroneckerDelta(2,gO2) +
      KroneckerDelta(3,gO2))*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVZ(int gO2) const
{
   
   const std::complex<double> result = -0.1*g2*(7.0710678118654755*g2*vT*Cos(
      ThetaW())*KroneckerDelta(2,gO2) + 7.0710678118654755*g2*vT*Cos(ThetaW())*
      KroneckerDelta(3,gO2) + 3.872983346207417*g1*(-(vd*KroneckerDelta(0,gO2)) +
      vu*KroneckerDelta(1,gO2))*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.1*(20*(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2) + KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2)
      *Sqr(Cos(ThetaW())) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-
      7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())
      + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0
      ,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) + 2*KroneckerDelta(2,gO1
      )*KroneckerDelta(2,gO2) + 2*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr
      (g2);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*(5*
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA(
      gI1,0)*ZA(gI2,1)) - 7.0710678118654755*(KroneckerDelta(2,gO2) +
      KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,2)) - KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1
      ,0)*ZA(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(
      gI2,1) + 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2))) - 5*KroneckerDelta(3,gO1)*
      (1.4142135623730951*AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1
      ) - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,1) +
      1.4142135623730951*AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2)
      - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,2) + 4*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,2) + 1.4142135623730951*
      KroneckerDelta(0,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,2)) + 2*KroneckerDelta(3,gO2)*(Sqr(g2)*ZA(gI1,0)*ZA(gI2,0) -
      (-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 2*Sqr(g2)*ZA(gI1,2)*ZA(
      gI2,2))) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax
      ) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(
      gI1,1)*ZA(gI2,1) - 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2)) + 5*(
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA(
      gI1,0)*ZA(gI2,1)) + 1.4142135623730951*(KroneckerDelta(2,gO2) +
      KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,1) +
      ZA(gI1,1)*ZA(gI2,2)))) + 5*KroneckerDelta(2,gO1)*(-1.4142135623730951*AbsSqr
      (Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1) + 1.4142135623730951*
      KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,1) - 1.4142135623730951*
      AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2) +
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,2) - 4*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,2) - 1.4142135623730951*
      KroneckerDelta(0,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,2)) + 2*KroneckerDelta(2,gO2)*((-2*AbsSqr(Lambdax) + Sqr(g2))
      *ZA(gI1,0)*ZA(gI2,0) - Sqr(g2)*(ZA(gI1,1)*ZA(gI2,1) + 2*ZA(gI1,2)*ZA(gI2,2))
      )));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (20*
      AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 10*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2)) + 5*(KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) +
      Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + 1.4142135623730951*(
      KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*
      (ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2))))) + KroneckerDelta(1,gO1)*(
      KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,
      0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) - 10*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2)) - 5*(KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) +
      Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + 1.4142135623730951*(
      KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*
      (ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))) + 5*(KroneckerDelta(3,gO1)*(-
      1.4142135623730951*AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,1)
      + 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,2)*ZH(gI2,1) -
      1.4142135623730951*AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,2)
      + 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,2) + 4*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZH(gI1,2)*ZH(gI2,2) + 1.4142135623730951*
      KroneckerDelta(0,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,2)*ZH(gI2,0) + ZH
      (gI1,0)*ZH(gI2,2)) - 2*KroneckerDelta(3,gO2)*(Sqr(g2)*ZH(gI1,0)*ZH(gI2,0) -
      (-2*AbsSqr(Lambdax) + Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 2*Sqr(g2)*ZH(gI1,2)*ZH(
      gI2,2))) + KroneckerDelta(2,gO1)*(1.4142135623730951*AbsSqr(Lambdax)*
      KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,1) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,2)*ZH(gI2,1) + 1.4142135623730951*
      AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,2) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,2) + 4*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZH(gI1,2)*ZH(gI2,2) - 1.4142135623730951*
      KroneckerDelta(0,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,2)*ZH(gI2,0) + ZH
      (gI1,0)*ZH(gI2,2)) + 2*KroneckerDelta(2,gO2)*((-2*AbsSqr(Lambdax) + Sqr(g2))
      *ZH(gI1,0)*ZH(gI2,0) - Sqr(g2)*(ZH(gI1,1)*ZH(gI2,1) + 2*ZH(gI1,2)*ZH(gI2,2))
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSvconjUHpmconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(KroneckerDelta(gI1,gI2)*(
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) - 5*Sqr(g2)) + 10*(
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2) - KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2))*Sqr(g2)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2
      )*(KroneckerDelta(gI1,gI2)*(-3*Sqr(g1) + 5*Sqr(g2)) - 20*SUM(j3,0,2,SUM(j2,0
      ,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gI2,j3))));

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
   
   const std::complex<double> result = KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(Ye(j1
      ,gI1))*ZER(gI2,j1));

   return result;
}

double CLASSNAME::CpbarFvFeconjUHpmPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUHpmconjSv(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(-(KroneckerDelta(1,gO2)*(
      1.4142135623730951*vu*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZV(gI1,j1)) + 2*(
      vT*Lambdax - 2*Mu)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)
      ))*ZV(gI1,j2)))) - 2*(vT*KroneckerDelta(2,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(
      gI2,j1))*ZV(gI1,j1)) + KroneckerDelta(3,gO2)*(-(vT*Sqr(g2)*SUM(j1,0,2,Conj(
      ZE(gI2,j1))*ZV(gI1,j1))) + 1.4142135623730951*vu*Lambdax*SUM(j2,0,2,SUM(j1,0
      ,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZV(gI1,j2)))) + KroneckerDelta(0,
      gO2)*(-1.4142135623730951*vd*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZV(gI1,j1))
      + 4*SUM(j2,0,2,SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2)))*ZV(gI1,j2))
      + 2.8284271247461903*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,
      Conj(Ye(j1,j3))*Ye(j1,j2)))*ZV(gI1,j3))));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUHpmconjHpmconjUHpm(int gI1, int gO1, int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*((
      KroneckerDelta(1,gO2)*(10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,1)
      + 10*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2) - 10*KroneckerDelta(3,gO2)*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,3))*ZP(gI2,0) + KroneckerDelta(0,gO2)*(2*(
      3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (10*AbsSqr(Lambdax) - 3*Sqr(g1)
      - 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 10*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(1,gO1)*((
      KroneckerDelta(0,gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0
      ) + 10*KroneckerDelta(2,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2) - 10*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3))*ZP(gI2,1) + KroneckerDelta(1,gO2)*(
      (-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - 2*((3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 5*(-2*AbsSqr(Lambdax) + Sqr(g2))*
      ZP(gI1,2)*ZP(gI2,2) + 5*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)))) - 10*(KroneckerDelta(
      2,gO1)*((KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,0) - KroneckerDelta(1,gO2)*(-2
      *AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1) - 2*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(
      gI1,3))*ZP(gI2,2) + KroneckerDelta(2,gO2)*(Sqr(g2)*ZP(gI1,0)*ZP(gI2,0) - (-2
      *AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 2*Sqr(g2)*(2*ZP(gI1,2)*ZP(
      gI2,2) - ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(3,gO1)*((-(KroneckerDelta(0
      ,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0)) + Sqr(g2)*(KroneckerDelta(1,
      gO2)*ZP(gI1,1) - 2*KroneckerDelta(2,gO2)*ZP(gI1,2)))*ZP(gI2,3) +
      KroneckerDelta(3,gO2)*(-((-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0)*ZP(gI2,0))
      + Sqr(g2)*(ZP(gI1,1)*ZP(gI2,1) - 2*ZP(gI1,2)*ZP(gI2,2) + 4*ZP(gI1,3)*ZP(gI2,
      3))))));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.25)*(
      1.4142135623730951*vT*AbsSqr(Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1
      ,0) + 2.8284271247461903*Conj(Mu)*KroneckerDelta(3,gO2)*Lambdax*ZA(gI2,0)*ZP
      (gI1,0) - 1.4142135623730951*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(
      gI1,0) + 2.8284271247461903*KroneckerDelta(3,gO2)*TLambdax*ZA(gI2,1)*ZP(gI1,
      0) + 1.4142135623730951*vd*AbsSqr(Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,2)*
      ZP(gI1,0) + 2*Conj(Mu)*KroneckerDelta(0,gO2)*Lambdax*ZA(gI2,2)*ZP(gI1,0) - 2
      *Conj(Lambdax)*KroneckerDelta(0,gO2)*Mu*ZA(gI2,2)*ZP(gI1,0) -
      1.4142135623730951*vd*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,0) - vu
      *AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,1) +
      5.656854249492381*MT*Conj(Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,1)
      - vu*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,1) - vd*AbsSqr(Lambdax)*
      KroneckerDelta(0,gO2)*ZA(gI2,1)*ZP(gI1,1) - 1.4142135623730951*vT*AbsSqr(
      Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,1)*ZP(gI1,1) + 2.8284271247461903*Conj
      (Mu)*KroneckerDelta(3,gO2)*Lambdax*ZA(gI2,1)*ZP(gI1,1) - vd*KroneckerDelta(0
      ,gO2)*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) + 1.4142135623730951*vT*KroneckerDelta(3,
      gO2)*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) + 4*MT*Conj(Lambdax)*KroneckerDelta(0,gO2)*
      ZA(gI2,2)*ZP(gI1,1) - 2*Conj(TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(
      gI1,1) + 1.4142135623730951*vu*AbsSqr(Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,
      2)*ZP(gI1,1) - 1.4142135623730951*vu*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,2)
      *ZP(gI1,1) + 1.4142135623730951*vT*AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA(
      gI2,0)*ZP(gI1,2) - 2.8284271247461903*Conj(Mu)*KroneckerDelta(0,gO2)*Lambdax
      *ZA(gI2,0)*ZP(gI1,2) - 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*
      ZA(gI2,0)*ZP(gI1,2) - 5.656854249492381*MT*Conj(Lambdax)*KroneckerDelta(0,
      gO2)*ZA(gI2,1)*ZP(gI1,2) - 1.4142135623730951*vd*AbsSqr(Lambdax)*
      KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,2) + 1.4142135623730951*vd*
      KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,2) - 4*vT*KroneckerDelta(3,
      gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,2) - 1.4142135623730951*vT*AbsSqr(Lambdax)*
      KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,3) - 2.8284271247461903*Conj(Lambdax)
      *KroneckerDelta(0,gO2)*Mu*ZA(gI2,0)*ZP(gI1,3) + 1.4142135623730951*vT*
      KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,3) - 2.8284271247461903*Conj(
      TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,1)*ZP(gI1,3) - 1.4142135623730951*vd*
      AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,3) +
      1.4142135623730951*vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,3) +
      KroneckerDelta(2,gO2)*(5.656854249492381*Conj(MT)*Lambdax*ZA(gI2,1)*ZP(gI1,0
      ) + 1.4142135623730951*vd*AbsSqr(Lambdax)*ZA(gI2,2)*ZP(gI1,0) -
      1.4142135623730951*vd*Sqr(g2)*ZA(gI2,2)*ZP(gI1,0) + 1.4142135623730951*vT*
      AbsSqr(Lambdax)*ZA(gI2,1)*ZP(gI1,1) + 2.8284271247461903*Conj(Lambdax)*Mu*ZA
      (gI2,1)*ZP(gI1,1) - 1.4142135623730951*vT*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) +
      1.4142135623730951*vu*AbsSqr(Lambdax)*ZA(gI2,2)*ZP(gI1,1) -
      1.4142135623730951*vu*Sqr(g2)*ZA(gI2,2)*ZP(gI1,1) + 1.4142135623730951*ZA(
      gI2,0)*((-(vT*AbsSqr(Lambdax)) + 2*Conj(Lambdax)*Mu + vT*Sqr(g2))*ZP(gI1,0)
      + 2*Conj(TLambdax)*ZP(gI1,1)) + 4*vT*Sqr(g2)*ZA(gI2,2)*ZP(gI1,3)) +
      KroneckerDelta(1,gO2)*(ZA(gI2,2)*(-4*Conj(MT)*Lambdax*ZP(gI1,0) + 2*TLambdax
      *ZP(gI1,0) + 2*Conj(Mu)*Lambdax*ZP(gI1,1) - 2*Conj(Lambdax)*Mu*ZP(gI1,1) -
      1.4142135623730951*vu*AbsSqr(Lambdax)*ZP(gI1,2) + 1.4142135623730951*vu*Sqr(
      g2)*ZP(gI1,2) - 1.4142135623730951*vu*AbsSqr(Lambdax)*ZP(gI1,3) +
      1.4142135623730951*vu*Sqr(g2)*ZP(gI1,3)) + ZA(gI2,0)*(vu*(AbsSqr(Lambdax) +
      Sqr(g2))*ZP(gI1,0) - 2.8284271247461903*(TLambdax*ZP(gI1,2) + 2*Conj(MT)*
      Lambdax*ZP(gI1,3))) + ZA(gI2,1)*(vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) +
      1.4142135623730951*((-(vT*AbsSqr(Lambdax)) - 2*Conj(Mu)*Lambdax + vT*Sqr(g2)
      )*ZP(gI1,2) + (vT*AbsSqr(Lambdax) - 2*Conj(Lambdax)*Mu - vT*Sqr(g2))*ZP(gI1,
      3)))));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   
   const std::complex<double> result = 0.05*(-5*(KroneckerDelta(2,gO2)*(-
      5.656854249492381*Conj(MT)*Lambdax*ZH(gI2,1)*ZP(gI1,0) - 1.4142135623730951*
      vd*AbsSqr(Lambdax)*ZH(gI2,2)*ZP(gI1,0) + 1.4142135623730951*vd*Sqr(g2)*ZH(
      gI2,2)*ZP(gI1,0) - 1.4142135623730951*vT*AbsSqr(Lambdax)*ZH(gI2,1)*ZP(gI1,1)
      - 2.8284271247461903*Conj(Lambdax)*Mu*ZH(gI2,1)*ZP(gI1,1) +
      1.4142135623730951*vT*Sqr(g2)*ZH(gI2,1)*ZP(gI1,1) - 1.4142135623730951*vu*
      AbsSqr(Lambdax)*ZH(gI2,2)*ZP(gI1,1) + 1.4142135623730951*vu*Sqr(g2)*ZH(gI2,2
      )*ZP(gI1,1) + 2*vu*Sqr(g2)*ZH(gI2,1)*ZP(gI1,2) + 4*vT*Sqr(g2)*ZH(gI2,2)*ZP(
      gI1,2) + ZH(gI2,0)*(1.4142135623730951*(Conj(Lambdax)*(-(vT*Lambdax) + 2*Mu)
      + vT*Sqr(g2))*ZP(gI1,0) + 2.8284271247461903*Conj(TLambdax)*ZP(gI1,1) - 2*vd
      *(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2)) - 4*vT*Sqr(g2)*ZH(gI2,2)*ZP(gI1,3
      )) + KroneckerDelta(3,gO2)*(-2.8284271247461903*TLambdax*ZH(gI2,1)*ZP(gI1,0)
      + 1.4142135623730951*vd*AbsSqr(Lambdax)*ZH(gI2,2)*ZP(gI1,0) -
      1.4142135623730951*vd*Sqr(g2)*ZH(gI2,2)*ZP(gI1,0) + 1.4142135623730951*vT*
      AbsSqr(Lambdax)*ZH(gI2,1)*ZP(gI1,1) - 2.8284271247461903*Conj(Mu)*Lambdax*ZH
      (gI2,1)*ZP(gI1,1) - 1.4142135623730951*vT*Sqr(g2)*ZH(gI2,1)*ZP(gI1,1) +
      1.4142135623730951*vu*AbsSqr(Lambdax)*ZH(gI2,2)*ZP(gI1,1) -
      1.4142135623730951*vu*Sqr(g2)*ZH(gI2,2)*ZP(gI1,1) - 4*vT*Sqr(g2)*ZH(gI2,2)*
      ZP(gI1,2) + 4*vu*AbsSqr(Lambdax)*ZH(gI2,1)*ZP(gI1,3) - 2*vu*Sqr(g2)*ZH(gI2,1
      )*ZP(gI1,3) + 4*vT*Sqr(g2)*ZH(gI2,2)*ZP(gI1,3) + ZH(gI2,0)*(
      1.4142135623730951*(vT*AbsSqr(Lambdax) + 2*Conj(Mu)*Lambdax - vT*Sqr(g2))*ZP
      (gI1,0) + 5.656854249492381*MT*Conj(Lambdax)*ZP(gI1,1) + 2*vd*Sqr(g2)*ZP(gI1
      ,3)))) - KroneckerDelta(1,gO2)*(-5*ZH(gI2,2)*(4*Conj(MT)*Lambdax*ZP(gI1,0) +
      2*TLambdax*ZP(gI1,0) - 2*vT*AbsSqr(Lambdax)*ZP(gI1,1) + 2*Conj(Mu)*Lambdax*
      ZP(gI1,1) + 2*Conj(Lambdax)*Mu*ZP(gI1,1) + 1.4142135623730951*vu*AbsSqr(
      Lambdax)*ZP(gI1,2) - 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,2) -
      1.4142135623730951*vu*AbsSqr(Lambdax)*ZP(gI1,3) + 1.4142135623730951*vu*Sqr(
      g2)*ZP(gI1,3)) + ZH(gI2,0)*(5*vu*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) + vd*
      (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1) + 14.142135623730951*
      (TLambdax*ZP(gI1,2) + 2*Conj(MT)*Lambdax*ZP(gI1,3))) + ZH(gI2,1)*(5*vd*(
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)
      + 7.0710678118654755*((-(vT*AbsSqr(Lambdax)) - 2*Conj(Mu)*Lambdax + vT*Sqr(
      g2))*ZP(gI1,2) + (vT*AbsSqr(Lambdax) - 2*Conj(Lambdax)*Mu - vT*Sqr(g2))*ZP(
      gI1,3)))) - KroneckerDelta(0,gO2)*(ZH(gI2,1)*(vu*(20*AbsSqr(Lambdax) - 3*Sqr
      (g1) + 5*Sqr(g2))*ZP(gI1,0) + 5*vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1) -
      14.142135623730951*(2*MT*Conj(Lambdax)*ZP(gI1,2) + Conj(TLambdax)*ZP(gI1,3))
      ) + ZH(gI2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0) + 5*vu*(AbsSqr(Lambdax)
      + Sqr(g2))*ZP(gI1,1) + 7.0710678118654755*((-(vT*AbsSqr(Lambdax)) + 2*Conj(
      Mu)*Lambdax + vT*Sqr(g2))*ZP(gI1,2) + (Conj(Lambdax)*(vT*Lambdax + 2*Mu) -
      vT*Sqr(g2))*ZP(gI1,3))) + 5*ZH(gI2,2)*(-2*Conj(Mu)*Lambdax*ZP(gI1,0) - 2*
      Conj(TLambdax)*ZP(gI1,1) + 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,2) -
      1.4142135623730951*vd*Sqr(g2)*ZP(gI1,3) + Conj(Lambdax)*(2*(vT*Lambdax - Mu)
      *ZP(gI1,0) - 4*MT*ZP(gI1,1) + 1.4142135623730951*vd*Lambdax*(-ZP(gI1,2) + ZP
      (gI1,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjUHpmPR(int gI1, int gI2, int gO2) const
{
   
   const std::complex<double> result = -(Conj(Lambdax)*KroneckerDelta(2,gO2)*UP(
      gI2,1)*ZN(gI1,2)) + Conj(Lambdax)*KroneckerDelta(0,gO2)*UP(gI2,2)*ZN(gI1,3)
      - 0.1*KroneckerDelta(1,gO2)*(UP(gI2,1)*(5.477225575051661*g1*ZN(gI1,0) +
      7.0710678118654755*g2*ZN(gI1,1)) + 10*g2*UP(gI2,0)*ZN(gI1,3)) +
      0.7071067811865475*Conj(Lambdax)*KroneckerDelta(0,gO2)*UP(gI2,1)*ZN(gI1,4) +
      1.4142135623730951*g2*KroneckerDelta(3,gO2)*(-(UP(gI2,2)*ZN(gI1,1)) + UP(gI2
      ,0)*ZN(gI1,4));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjUHpmPL(int gI1, int gI2, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UM(gI2,0))*(Conj(ZN(gI1,2))*
      KroneckerDelta(0,gO1) + 1.4142135623730951*Conj(ZN(gI1,4))*KroneckerDelta(2,
      gO1))) + Conj(UM(gI2,2))*(1.4142135623730951*g2*Conj(ZN(gI1,1))*
      KroneckerDelta(2,gO1) - Conj(ZN(gI1,2))*KroneckerDelta(1,gO1)*Lambdax) +
      Conj(UM(gI2,1))*(0.5477225575051661*g1*Conj(ZN(gI1,0))*KroneckerDelta(0,gO1)
      + 0.7071067811865475*g2*Conj(ZN(gI1,1))*KroneckerDelta(0,gO1) +
      0.7071067811865475*Conj(ZN(gI1,4))*KroneckerDelta(1,gO1)*Lambdax + Conj(ZN(
      gI1,3))*KroneckerDelta(3,gO1)*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSdconjUHpmconjSd(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(10*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) - 10*
      KroneckerDelta(2,gO1)*(KroneckerDelta(2,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,
      j1))*ZD(gI2,j1)) + 2*KroneckerDelta(1,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*ZD(gI2,j2))) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) - 20*
      SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,
      j1))*ZD(gI2,3 + j2)))) - KroneckerDelta(1,gO1)*(20*Conj(Lambdax)*
      KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(
      gI2,3 + j1))) + KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj
      (ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,
      3 + j1)) + 20*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,
      j3))*Yu(j1,j2)))*ZD(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSeconjUHpmconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(10*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1
      ,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1)
      )*ZE(gI2,3 + j1))) - 20*Conj(Lambdax)*KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(
      ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI2,3 + j1)))) - 10*KroneckerDelta(2,gO1
      )*(KroneckerDelta(2,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 2
      *KroneckerDelta(1,gO2)*Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE
      (gI1,3 + j1)))*ZE(gI2,j2))) - KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((
      3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 6*Sqr(g1)*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)) + 20*SUM(j3,0,2,Conj(ZE(gI1,
      3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2)))))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSuconjUHpmconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.05*(10*(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) +
      KroneckerDelta(3,gO1)*(-(KroneckerDelta(3,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(
      gI1,j1))*ZU(gI2,j1))) + 2*Conj(Lambdax)*KroneckerDelta(0,gO2)*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI2,3 + j1))))) - KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1)
      )*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) +
      20*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(
      j2,j1))*ZU(gI2,3 + j2)))) + KroneckerDelta(0,gO1)*(20*KroneckerDelta(3,gO2)*
      Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*ZU(gI2,
      j2)) + KroneckerDelta(0,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,
      j1))*ZU(gI2,j1)) - 4*(Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)
      ) + 5*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(
      j1,j2)))*ZU(gI2,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUHpmconjSu(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.25*(-2*vT*KroneckerDelta(2,gO2)*Sqr(g2)*
      SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 2*vT*KroneckerDelta(3,gO2)*Sqr(g2)
      *SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 2.8284271247461903*vd*Conj(
      Lambdax)*KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Yu(j1,
      j2)*ZU(gI1,3 + j1))) - 2.8284271247461903*vu*KroneckerDelta(3,gO2)*Lambdax*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZU(gI1,j2)) +
      KroneckerDelta(0,gO2)*(-1.4142135623730951*vd*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2
      ,j1))*ZU(gI1,j1)) - 2*vT*Conj(Lambdax)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,
      2,Yu(j1,j2)*ZU(gI1,3 + j1))) + 4*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1
      ,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + 4*SUM(j2,0,2,SUM(j1,0,2,Conj(ZD(gI2,3 + j1
      ))*Conj(TYd(j1,j2)))*ZU(gI1,j2)) + 2.8284271247461903*vu*SUM(j3,0,2,Conj(ZD(
      gI2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2)
      )) + 2.8284271247461903*vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2
      ,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gI1,j3))) + KroneckerDelta(1,gO2)*(-
      1.4142135623730951*vu*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 4*
      SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,ZU(gI1,3 + j1)*TYu(j1,j2))) - 2*vT*
      Lambdax*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZU(gI1,
      j2)) + 4*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZU(
      gI1,j2)) + 2.8284271247461903*vd*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2))) + 2.8284271247461903*
      vu*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,
      j2)))*ZU(gI1,j3))));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjUHpmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZA(gI2,0) + KroneckerDelta(1,gO2)*ZA(gI2,1) +
      1.4142135623730951*(KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2))*ZA(gI2,2)
      );

   return result;
}

std::complex<double> CLASSNAME::CphhconjUHpmVWm(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) -
      KroneckerDelta(1,gO2)*ZH(gI2,1) + 1.4142135623730951*(KroneckerDelta(2,gO2)
      + KroneckerDelta(3,gO2))*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmVP(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.1*(-(KroneckerDelta(0,gO2)*(
      3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW()))*ZP(gI2,0)) -
      KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW(
      )))*ZP(gI2,1) - 10*g2*Sin(ThetaW())*(KroneckerDelta(2,gO2)*ZP(gI2,2) +
      KroneckerDelta(3,gO2)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmVZ(int gI2, int gO2) const
{
   
   const std::complex<double> result = 0.1*(-(KroneckerDelta(0,gO2)*(5*g2*Cos(
      ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*ZP(gI2,0)) - KroneckerDelta(
      1,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*ZP(gI2,1) -
      10*g2*Cos(ThetaW())*(KroneckerDelta(2,gO2)*ZP(gI2,2) + KroneckerDelta(3,gO2)
      *ZP(gI2,3)));

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

std::complex<double> CLASSNAME::CpbarChaChaVPPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*Sin(ThetaW())*UM(gI1,0)
      + 0.5*Conj(UM(gI2,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()
      ))*UM(gI1,1) + g2*Conj(UM(gI2,2))*Sin(ThetaW())*UM(gI1,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVPPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UP(gI1,0))*Sin(ThetaW())*UP(gI2,0)
      + 0.5*Conj(UP(gI1,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()
      ))*UP(gI2,1) + g2*Conj(UP(gI1,2))*Sin(ThetaW())*UP(gI2,2);

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

std::complex<double> CLASSNAME::CpHpmconjHpmVPVP(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*((g2*Sin(ThetaW())*(7.745966692414834*
      g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())))*ZP(
      gI1,0)*ZP(gI2,0) + (g2*Sin(ThetaW())*(7.745966692414834*g1*Cos(ThetaW()) + 5
      *g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())))*ZP(gI1,1)*ZP(gI2,1) + 20*
      Sqr(g2)*Sqr(Sin(ThetaW()))*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVP(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5*(-((0.7745966692414834*g1*Cos(ThetaW())
      + g2*Sin(ThetaW()))*ZP(gI1,0)*ZP(gI2,0)) - (0.7745966692414834*g1*Cos(ThetaW
      ()) + g2*Sin(ThetaW()))*ZP(gI1,1)*ZP(gI2,1) - 2*g2*Sin(ThetaW())*(ZP(gI1,2)*
      ZP(gI2,2) + ZP(gI1,3)*ZP(gI2,3)));

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
   
   const std::complex<double> result = -0.5*g2*(0.7745966692414834*g1*vd*Cos(
      ThetaW())*ZP(gI2,0) - 0.7745966692414834*g1*vu*Cos(ThetaW())*ZP(gI2,1) +
      1.4142135623730951*g2*vT*Sin(ThetaW())*(ZP(gI2,2) + ZP(gI2,3)));

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
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargWmCgWmCVZ() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmVWmVZ() const
{
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZA(
      gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW(
      )) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZH(
      gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*KroneckerDelta(gI1,gI2)*(g1*Sin(ThetaW(
      ))*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr
      (Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(ThetaW(
      )) + 0.7745966692414834*g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) - ZA(gI2,1)*
      ZH(gI1,1));

   return result;
}

double CLASSNAME::CpSvconjSvVZ(int gI2, int gI1) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*UM(gI1,0)
      + 0.5*Conj(UM(gI2,1))*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()
      ))*UM(gI1,1) + g2*Conj(UM(gI2,2))*Cos(ThetaW())*UM(gI1,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarChaChaVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UP(gI1,0))*Cos(ThetaW())*UP(gI2,0)
      + 0.5*Conj(UP(gI1,1))*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()
      ))*UP(gI2,1) + g2*Conj(UP(gI1,2))*Cos(ThetaW())*UP(gI2,2);

   return result;
}

double CLASSNAME::CpbarFdFdVZPL(int gI1, int gI2) const
{
   
   const double result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFdFdVZPR(int gI1, int gI2) const
{
   
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW
      ());

   return result;
}

double CLASSNAME::CpbarFeFeVZPL(int gI1, int gI2) const
{
   
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) -
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVZPR(int gI1, int gI2) const
{
   
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW
      ());

   return result;
}

double CLASSNAME::CpbarFuFuVZPL(int gI1, int gI2) const
{
   
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPR(int gI1, int gI2) const
{
   
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW(
      ));

   return result;
}

double CLASSNAME::CpbarFvFvVZPL(int gI1, int gI2) const
{
   
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFvFvVZPR(int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*((-7.745966692414834*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())
      ))*ZP(gI1,0)*ZP(gI2,0) + (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW(
      )) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*ZP(gI1,1)*
      ZP(gI2,1) + 20*Sqr(g2)*Sqr(Cos(ThetaW()))*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)*
      ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*((-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*ZP(gI1,0)*ZP(gI2,0) + (-5*g2*Cos(ThetaW(
      )) + 3.872983346207417*g1*Sin(ThetaW()))*ZP(gI1,1)*ZP(gI2,1) - 10*g2*Cos(
      ThetaW())*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.5*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*(Conj(ZN(gI2,2))*ZN(gI1,2) - Conj(ZN(gI2,3))*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiVZPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*
      g1*Sin(ThetaW()))*(Conj(ZN(gI1,2))*ZN(gI2,2) - Conj(ZN(gI1,3))*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.03333333333333333*((g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + g1*Sin(ThetaW())) + 15*Sqr(g2)*Sqr(Cos(
      ThetaW())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*Sqr(g1)*Sqr(Sin(
      ThetaW()))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.1*((-7.745966692414834*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())
      ))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 12*Sqr(g1)*Sqr(Sin(ThetaW()))*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZVZ(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.03333333333333333*((-7.745966692414834*g1
      *g2*Cos(ThetaW())*Sin(ThetaW()) + 15*Sqr(g2)*Sqr(Cos(ThetaW())) + Sqr(g1)*
      Sqr(Sin(ThetaW())))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + 16*Sqr(g1)*Sqr
      (Sin(ThetaW()))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.16666666666666666*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))
      + 0.2581988897471611*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1
      ,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.1*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1))
      + 0.7745966692414834*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1
      ,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZ(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.16666666666666666*(3*g2*Cos(ThetaW()) -
      0.7745966692414834*g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))
      - 0.5163977794943222*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1
      ,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CphhVZVZ(int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*(vd*ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjVWmVZ(int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(-0.7745966692414834*g1*vd*Sin(
      ThetaW())*ZP(gI2,0) + 0.7745966692414834*g1*vu*Sin(ThetaW())*ZP(gI2,1) +
      1.4142135623730951*g2*vT*Cos(ThetaW())*(ZP(gI2,2) + ZP(gI2,3)));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ1() const
{
   
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ2() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVWmVZVZ3() const
{
   
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

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
   
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpbargZgWmconjVWm() const
{
   
   const double result = g2*Cos(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpAhAhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1
      )*ZA(gI2,1) + 4*ZA(gI1,2)*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CphhhhconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1
      )*ZH(gI2,1) + 4*ZH(gI1,2)*ZH(gI2,2));

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

std::complex<double> CLASSNAME::CpHpmconjHpmconjVWmVWm(int gI1, int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1
      )*ZP(gI2,1) + 2*ZP(gI1,2)*ZP(gI2,2) + 2*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(ZA(gI2,0)*
      ZP(gI1,0) + ZA(gI2,1)*ZP(gI1,1) + 1.4142135623730951*ZA(gI2,2)*(ZP(gI1,2) -
      ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP
      (gI1,1) + 1.4142135623730951*ZH(gI2,2)*(ZP(gI1,2) + ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjVWmPL(int gI1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,1) +
      1.4142135623730951*Conj(UM(gI2,1))*ZN(gI1,2) + 2*Conj(UM(gI2,2))*ZN(gI1,4));

   return result;
}

std::complex<double> CLASSNAME::CpChiChaconjVWmPR(int gI1, int gI2) const
{
   
   const std::complex<double> result = -(g2*Conj(ZN(gI1,1))*UP(gI2,0)) +
      0.7071067811865475*g2*Conj(ZN(gI1,3))*UP(gI2,1) - g2*Conj(ZN(gI1,4))*UP(gI2,
      2);

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

std::complex<double> CLASSNAME::CpSdconjSuconjVWm(int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*ZU(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CphhconjVWmVWm(int gI2) const
{
   
   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1) +
      4*vT*ZH(gI2,2));

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

double CLASSNAME::CpbarFvUChiSvPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvUChiSvPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,0.5477225575051661*g1*Conj(ZV(
      gI2,gI1))*KroneckerDelta(0,gO1),0) + IF(gI1 < 3,-0.7071067811865475*g2*Conj(
      ZV(gI2,gI1))*KroneckerDelta(1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiFvconjSvPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5477225575051661*g1*
      KroneckerDelta(0,gO2)*ZV(gI1,gI2),0) + IF(gI2 < 3,-0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZV(gI1,gI2),0);

   return result;
}

double CLASSNAME::CpUChiFvconjSvPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiHpmPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = Conj(UP(gI1,1))*(0.7071067811865475*
      KroneckerDelta(4,gO2)*Lambdax*ZP(gI2,0) - 0.5477225575051661*g1*
      KroneckerDelta(0,gO2)*ZP(gI2,1) - 0.7071067811865475*g2*KroneckerDelta(1,gO2
      )*ZP(gI2,1) - KroneckerDelta(2,gO2)*Lambdax*ZP(gI2,2)) + Conj(UP(gI1,2))*(
      KroneckerDelta(3,gO2)*Lambdax*ZP(gI2,0) - 1.4142135623730951*g2*
      KroneckerDelta(1,gO2)*ZP(gI2,3)) + g2*Conj(UP(gI1,0))*(-(KroneckerDelta(3,
      gO2)*ZP(gI2,1)) + 1.4142135623730951*KroneckerDelta(4,gO2)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiHpmPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO1)
      *UM(gI1,1)*ZP(gI2,0) + 0.7071067811865475*g2*KroneckerDelta(1,gO1)*UM(gI1,1)
      *ZP(gI2,0) + 0.7071067811865475*Conj(Lambdax)*KroneckerDelta(4,gO1)*UM(gI1,1
      )*ZP(gI2,1) - KroneckerDelta(2,gO1)*(g2*UM(gI1,0)*ZP(gI2,0) + Conj(Lambdax)*
      UM(gI1,2)*ZP(gI2,1)) - 1.4142135623730951*g2*KroneckerDelta(4,gO1)*UM(gI1,0)
      *ZP(gI2,2) + 1.4142135623730951*g2*KroneckerDelta(1,gO1)*UM(gI1,2)*ZP(gI2,2)
      + Conj(Lambdax)*KroneckerDelta(3,gO1)*UM(gI1,1)*ZP(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpChiUChihhPL(int gI2, int gO2, int gI1) const
{
   
   const std::complex<double> result = 0.1*(-5*g2*Conj(ZN(gI2,1))*KroneckerDelta(2
      ,gO2)*ZH(gI1,0) + 5*Conj(ZN(gI2,4))*KroneckerDelta(3,gO2)*Lambdax*ZH(gI1,0)
      + 5*Conj(ZN(gI2,3))*KroneckerDelta(4,gO2)*Lambdax*ZH(gI1,0) -
      3.872983346207417*g1*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2)*ZH(gI1,1) + 5*g2*
      Conj(ZN(gI2,3))*KroneckerDelta(1,gO2)*ZH(gI1,1) + 5*g2*Conj(ZN(gI2,1))*
      KroneckerDelta(3,gO2)*ZH(gI1,1) + 5*Conj(ZN(gI2,4))*KroneckerDelta(2,gO2)*
      Lambdax*ZH(gI1,1) + 3.872983346207417*g1*Conj(ZN(gI2,0))*(KroneckerDelta(2,
      gO2)*ZH(gI1,0) - KroneckerDelta(3,gO2)*ZH(gI1,1)) + 5*Conj(ZN(gI2,3))*
      KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,2) + Conj(ZN(gI2,2))*(3.872983346207417
      *g1*KroneckerDelta(0,gO2)*ZH(gI1,0) - 5*g2*KroneckerDelta(1,gO2)*ZH(gI1,0) +
      5*KroneckerDelta(4,gO2)*Lambdax*ZH(gI1,1) + 5*KroneckerDelta(3,gO2)*Lambdax*
      ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChihhPR(int gI2, int gO1, int gI1) const
{
   
   const std::complex<double> result = 0.1*(3.872983346207417*g1*KroneckerDelta(0,
      gO1)*ZH(gI1,0)*ZN(gI2,2) - 5*g2*KroneckerDelta(1,gO1)*ZH(gI1,0)*ZN(gI2,2) +
      5*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZH(gI1,1)*ZN(gI2,2) + 5*Conj(Lambdax)*
      KroneckerDelta(4,gO1)*ZH(gI1,0)*ZN(gI2,3) - 3.872983346207417*g1*
      KroneckerDelta(0,gO1)*ZH(gI1,1)*ZN(gI2,3) + 5*g2*KroneckerDelta(1,gO1)*ZH(
      gI1,1)*ZN(gI2,3) + KroneckerDelta(3,gO1)*(ZH(gI1,1)*(-3.872983346207417*g1*
      ZN(gI2,0) + 5*g2*ZN(gI2,1)) + 5*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,2) + ZH(gI1,
      0)*ZN(gI2,4))) + KroneckerDelta(2,gO1)*(ZH(gI1,0)*(3.872983346207417*g1*ZN(
      gI2,0) - 5*g2*ZN(gI2,1)) + 5*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,3) + ZH(gI1,1)*
      ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdUChiSdPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -0.3651483716701107*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*Conj(ZDR(gI1,j1))) - KroneckerDelta(2,gO2)
      *SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdUChiSdPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZDL(gI1,j1)) -
      KroneckerDelta(2,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 +
      j1)))*ZDL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUChiSePL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = -1.0954451150103321*g1*KroneckerDelta(0,gO2
      )*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*Conj(ZER(gI1,j1))) - KroneckerDelta(2,gO2)
      *SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeUChiSePR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO1)
      *SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZEL(gI1,j1)) -
      KroneckerDelta(2,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 +
      j1)))*ZEL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuUChiSuPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = 0.7302967433402214*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*Conj(ZUR(gI1,j1))) - KroneckerDelta(3,gO2)*
      SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuUChiSuPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1)) - 0.7071067811865475*g2*
      KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZUL(gI1,j1)) -
      KroneckerDelta(3,gO1)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI2,3 +
      j1)))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiVWmPL(int gI1, int gO2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*KroneckerDelta(1,gO2)*UM(gI1,0)
      + 1.4142135623730951*KroneckerDelta(2,gO2)*UM(gI1,1) + 2*KroneckerDelta(4,
      gO2)*UM(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpbarChaUChiVWmPR(int gI1, int gO1) const
{
   
   const std::complex<double> result = -(g2*Conj(UP(gI1,0))*KroneckerDelta(1,gO1))
      + 0.7071067811865475*g2*Conj(UP(gI1,1))*KroneckerDelta(3,gO1) - g2*Conj(UP(
      gI1,2))*KroneckerDelta(4,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = Conj(UM(gI2,2))*(-(KroneckerDelta(2,gO2)*
      Lambdax*ZP(gI1,1)) + 1.4142135623730951*g2*KroneckerDelta(1,gO2)*ZP(gI1,2))
      - g2*Conj(UM(gI2,0))*(KroneckerDelta(2,gO2)*ZP(gI1,0) + 1.4142135623730951*
      KroneckerDelta(4,gO2)*ZP(gI1,2)) + Conj(UM(gI2,1))*(0.5477225575051661*g1*
      KroneckerDelta(0,gO2)*ZP(gI1,0) + 0.7071067811865475*g2*KroneckerDelta(1,gO2
      )*ZP(gI1,0) + 0.7071067811865475*KroneckerDelta(4,gO2)*Lambdax*ZP(gI1,1) +
      KroneckerDelta(3,gO2)*Lambdax*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(3,gO1)*UP(gI2,0)*ZP(gI1
      ,1)) + Conj(Lambdax)*(0.7071067811865475*KroneckerDelta(4,gO1)*UP(gI2,1)*ZP(
      gI1,0) + KroneckerDelta(3,gO1)*UP(gI2,2)*ZP(gI1,0) - KroneckerDelta(2,gO1)*
      UP(gI2,1)*ZP(gI1,2)) + 0.1414213562373095*(-3.872983346207417*g1*
      KroneckerDelta(0,gO1)*UP(gI2,1)*ZP(gI1,1) + 10*g2*KroneckerDelta(4,gO1)*UP(
      gI2,0)*ZP(gI1,3) - 5*g2*KroneckerDelta(1,gO1)*(UP(gI2,1)*ZP(gI1,1) + 2*UP(
      gI2,2)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiAhPL(int gI1, int gO2, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Conj(ZN(
      gI1,1))*KroneckerDelta(2,gO2)*ZA(gI2,0) + 5*Conj(ZN(gI1,4))*KroneckerDelta(3
      ,gO2)*Lambdax*ZA(gI2,0) + 5*Conj(ZN(gI1,3))*KroneckerDelta(4,gO2)*Lambdax*ZA
      (gI2,0) + 3.872983346207417*g1*Conj(ZN(gI1,3))*KroneckerDelta(0,gO2)*ZA(gI2,
      1) - 5*g2*Conj(ZN(gI1,3))*KroneckerDelta(1,gO2)*ZA(gI2,1) - 5*g2*Conj(ZN(gI1
      ,1))*KroneckerDelta(3,gO2)*ZA(gI2,1) + 5*Conj(ZN(gI1,4))*KroneckerDelta(2,
      gO2)*Lambdax*ZA(gI2,1) + 3.872983346207417*g1*Conj(ZN(gI1,0))*(-(
      KroneckerDelta(2,gO2)*ZA(gI2,0)) + KroneckerDelta(3,gO2)*ZA(gI2,1)) + 5*Conj
      (ZN(gI1,3))*KroneckerDelta(2,gO2)*Lambdax*ZA(gI2,2) + Conj(ZN(gI1,2))*(-
      3.872983346207417*g1*KroneckerDelta(0,gO2)*ZA(gI2,0) + 5*(g2*KroneckerDelta(
      1,gO2)*ZA(gI2,0) + KroneckerDelta(4,gO2)*Lambdax*ZA(gI2,1) + KroneckerDelta(
      3,gO2)*Lambdax*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiAhPR(int gI1, int gO1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.1)*(
      3.872983346207417*g1*KroneckerDelta(0,gO1)*ZA(gI2,0)*ZN(gI1,2) - 5*g2*
      KroneckerDelta(1,gO1)*ZA(gI2,0)*ZN(gI1,2) - 5*Conj(Lambdax)*KroneckerDelta(4
      ,gO1)*ZA(gI2,1)*ZN(gI1,2) - 5*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZA(gI2,0)*
      ZN(gI1,3) - 3.872983346207417*g1*KroneckerDelta(0,gO1)*ZA(gI2,1)*ZN(gI1,3) +
      5*g2*KroneckerDelta(1,gO1)*ZA(gI2,1)*ZN(gI1,3) - KroneckerDelta(3,gO1)*(ZA(
      gI2,1)*(3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1)) + 5*Conj(Lambdax)*(
      ZA(gI2,2)*ZN(gI1,2) + ZA(gI2,0)*ZN(gI1,4))) + KroneckerDelta(2,gO1)*(ZA(gI2,
      0)*(3.872983346207417*g1*ZN(gI1,0) - 5*g2*ZN(gI1,1)) - 5*Conj(Lambdax)*(ZA(
      gI2,2)*ZN(gI1,3) + ZA(gI2,1)*ZN(gI1,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFdconjSdPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZDL(gI2,j1))*ZD(gI1,j1)) -
      KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(
      gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFdconjSdPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.3651483716701107*g1*KroneckerDelta(0,gO1
      )*SUM(j1,0,2,ZD(gI1,3 + j1)*ZDR(gI2,j1)) - KroneckerDelta(2,gO1)*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZD(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFeconjSePL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,gO2)
      *SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(gI1,j1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZEL(gI2,j1))*ZE(gI1,j1)) -
      KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(ZEL(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(
      gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFeconjSePR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.0954451150103321*g1*KroneckerDelta(0,gO1
      )*SUM(j1,0,2,ZE(gI1,3 + j1)*ZER(gI2,j1)) - KroneckerDelta(2,gO1)*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZE(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFuconjSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.18257418583505536*g1*KroneckerDelta(0,
      gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(gI1,j1)) - 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*SUM(j1,0,2,Conj(ZUL(gI2,j1))*ZU(gI1,j1)) -
      KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(
      gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiFuconjSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7302967433402214*g1*KroneckerDelta(0,gO1)
      *SUM(j1,0,2,ZU(gI1,3 + j1)*ZUR(gI2,j1)) - KroneckerDelta(3,gO1)*SUM(j2,0,2,
      SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZU(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)) +
      0.7071067811865475*g2*KroneckerDelta(3,gO2)*UP(gI2,1) - g2*KroneckerDelta(4,
      gO2)*UP(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChaconjVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(1
      ,gO1) + 1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(2,gO1) + 2*Conj(UM
      (gI2,2))*KroneckerDelta(4,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZPL(int gI2, int gO2) const
{
   
   const std::complex<double> result = -0.1*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(2,gO2)*ZN(gI2,2) -
      KroneckerDelta(3,gO2)*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpChiUChiVZPR(int gI2, int gO1) const
{
   
   const std::complex<double> result = 0.1*(Conj(ZN(gI2,2))*KroneckerDelta(2,gO1)
      - Conj(ZN(gI2,3))*KroneckerDelta(3,gO1))*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,0.5)*(Conj(UM(gI1,2)
      )*(-1.4142135623730951*KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,0) + 2*g2*
      KroneckerDelta(0,gO2)*ZA(gI2,2)) + g2*Conj(UM(gI1,0))*(1.4142135623730951*
      KroneckerDelta(1,gO2)*ZA(gI2,1) - 2*KroneckerDelta(2,gO2)*ZA(gI2,2)) + Conj(
      UM(gI1,1))*(1.4142135623730951*g2*KroneckerDelta(0,gO2)*ZA(gI2,0) +
      1.4142135623730951*KroneckerDelta(2,gO2)*Lambdax*ZA(gI2,1) + KroneckerDelta(
      1,gO2)*Lambdax*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0,-0.5)*(-
      1.4142135623730951*Conj(Lambdax)*KroneckerDelta(2,gO1)*UP(gI1,1)*ZA(gI2,0) +
      KroneckerDelta(1,gO1)*(1.4142135623730951*g2*UP(gI1,0)*ZA(gI2,0) + Conj(
      Lambdax)*(1.4142135623730951*UP(gI1,2)*ZA(gI2,1) + UP(gI1,1)*ZA(gI2,2))) +
      g2*(2*KroneckerDelta(2,gO1)*UP(gI1,0)*ZA(gI2,2) + KroneckerDelta(0,gO1)*(
      1.4142135623730951*UP(gI1,1)*ZA(gI2,1) - 2*UP(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChahhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5*(-(Conj(UM(gI2,2))*(1.4142135623730951*
      KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,0) + 2*g2*KroneckerDelta(0,gO2)*ZH(gI1,
      2))) + g2*Conj(UM(gI2,0))*(-1.4142135623730951*KroneckerDelta(1,gO2)*ZH(gI1,
      1) + 2*KroneckerDelta(2,gO2)*ZH(gI1,2)) + Conj(UM(gI2,1))*(-
      1.4142135623730951*g2*KroneckerDelta(0,gO2)*ZH(gI1,0) + 1.4142135623730951*
      KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1) + KroneckerDelta(1,gO2)*Lambdax*ZH(
      gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChahhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.5*(-1.4142135623730951*Conj(Lambdax)*
      KroneckerDelta(2,gO1)*UP(gI2,1)*ZH(gI1,0) + KroneckerDelta(1,gO1)*(-
      1.4142135623730951*g2*UP(gI2,0)*ZH(gI1,0) + Conj(Lambdax)*(
      1.4142135623730951*UP(gI2,2)*ZH(gI1,1) + UP(gI2,1)*ZH(gI1,2))) - g2*(2*
      KroneckerDelta(2,gO1)*UP(gI2,0)*ZH(gI1,2) + KroneckerDelta(0,gO1)*(
      1.4142135623730951*UP(gI2,1)*ZH(gI1,1) - 2*UP(gI2,2)*ZH(gI1,2))));

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
   
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j2,0,2,SUM(j1,0,2
      ,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZV(gI1,j2));

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

double CLASSNAME::CpbarUChabarFvSePL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gI1 < 3,-(g2*Conj(ZE(gI2,gI1))*
      KroneckerDelta(0,gO1)),0) + KroneckerDelta(1,gO1)*SUM(j1,0,2,Conj(Ye(j1,gI1)
      )*Conj(ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.5477225575051661*g1*Conj(ZN(gI2,0))*
      KroneckerDelta(1,gO2)*ZP(gI1,1) - 0.7071067811865475*g2*Conj(ZN(gI2,1))*
      KroneckerDelta(1,gO2)*ZP(gI1,1) + Conj(ZN(gI2,3))*(KroneckerDelta(2,gO2)*
      Lambdax*ZP(gI1,0) - g2*KroneckerDelta(0,gO2)*ZP(gI1,1)) - Conj(ZN(gI2,2))*
      KroneckerDelta(1,gO2)*Lambdax*ZP(gI1,2) - 1.4142135623730951*g2*Conj(ZN(gI2,
      1))*KroneckerDelta(2,gO2)*ZP(gI1,3) + 0.7071067811865475*Conj(ZN(gI2,4))*(
      KroneckerDelta(1,gO2)*Lambdax*ZP(gI1,0) + 2*g2*KroneckerDelta(0,gO2)*ZP(gI1,
      3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = KroneckerDelta(2,gO1)*(-(Conj(Lambdax)*ZN(
      gI2,2)*ZP(gI1,1)) + 1.4142135623730951*g2*ZN(gI2,1)*ZP(gI1,2)) - g2*
      KroneckerDelta(0,gO1)*(ZN(gI2,2)*ZP(gI1,0) + 1.4142135623730951*ZN(gI2,4)*ZP
      (gI1,2)) + KroneckerDelta(1,gO1)*(0.5477225575051661*g1*ZN(gI2,0)*ZP(gI1,0)
      + 0.7071067811865475*g2*ZN(gI2,1)*ZP(gI1,0) + 0.7071067811865475*Conj(
      Lambdax)*ZN(gI2,4)*ZP(gI1,1) + Conj(Lambdax)*ZN(gI2,3)*ZP(gI1,3));

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

std::complex<double> CLASSNAME::CpbarUChaChaVPPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*UP(
      gI2,0) + 0.1*KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*
      g2*Sin(ThetaW()))*UP(gI2,1) + g2*KroneckerDelta(2,gO2)*Sin(ThetaW())*UP(gI2,
      2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVPPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = g2*(Conj(UM(gI2,0))*KroneckerDelta(0,gO1) +
      Conj(UM(gI2,2))*KroneckerDelta(2,gO1))*Sin(ThetaW()) + 0.1*Conj(UM(gI2,1))*
      KroneckerDelta(1,gO1)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW(
      )));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*UP(
      gI2,0) + 0.1*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*
      g1*Sin(ThetaW()))*UP(gI2,1) + g2*Cos(ThetaW())*KroneckerDelta(2,gO2)*UP(gI2,
      2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*
      KroneckerDelta(0,gO1) + g2*Conj(UM(gI2,2))*Cos(ThetaW())*KroneckerDelta(2,
      gO1) + 0.1*Conj(UM(gI2,1))*KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiVWmPR(int gO2, int gI2) const
{
   
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,1)) +
      0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZN(gI2,3) - g2*KroneckerDelta(2,
      gO2)*ZN(gI2,4);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChiVWmPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = -0.5*g2*(2*Conj(ZN(gI2,1))*KroneckerDelta(0
      ,gO1) + 1.4142135623730951*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1) + 2*Conj(ZN
      (gI2,4))*KroneckerDelta(2,gO1));

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
      ),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gI1,j2))*Ye(gO2,j2))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,gO1))*ZER(gI1,j1))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(ZEL(gI2,j2))*Ye(gO2,j2))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Ye(j1,gO1))*ZER(gI2,j1))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,Ye(gO2,gI2)*ZP(gI1,0),0);

   return result;
}

double CLASSNAME::CpbarUFeFvHpmPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-1.0954451150103321*g1*Conj(ZE(
      gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,2))*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*Ye(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZE(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZE(gI1,gO1))*
      ZN(gI2,1),0) + IF(gO1 < 3,-(SUM(j1,0,2,Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1))
      )*ZN(gI2,2)),0);

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
   
   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Sin(
      ThetaW())*ZER(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZEL(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Sin(
      ThetaW()),0);

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

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j2,0,2,Conj(ZDL(gI1,j2))*Yd(gO2,j2))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,gO1))*ZDR(gI1,j1))*ZA(gI2,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(ZDL(gI2,j2))*Yd(gO2,j2))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Yd(j1,gO1))*ZDR(gI2,j1))*ZH(gI1,0),0);

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
   
   const std::complex<double> result = IF(gO2 < 3,-0.3651483716701107*g1*Conj(ZD(
      gI1,3 + gO2))*Conj(ZN(gI2,0)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,2))*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*Yd(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZD(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZD(gI1,gO1))*
      ZN(gI2,1),0) + IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1))
      )*ZN(gI2,2)),0);

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
   
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Sin(
      ThetaW())*ZDR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZDL(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Sin(
      ThetaW()),0);

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

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,-
      0.7071067811865475)*SUM(j2,0,2,Conj(ZUL(gI1,j2))*Yu(gO2,j2))*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*ZUR(gI1,j1))*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,2,
      Conj(ZUL(gI2,j2))*Yu(gO2,j2))*ZH(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,2,
      Conj(Yu(j1,gO1))*ZUR(gI2,j1))*ZH(gI1,1),0);

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

std::complex<double> CLASSNAME::CpbarUFuChiSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO2 < 3,0.7302967433402214*g1*Conj(ZN(
      gI2,0))*Conj(ZU(gI1,3 + gO2)),0) + IF(gO2 < 3,-(Conj(ZN(gI2,3))*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*Yu(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(ZU(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZU(gI1,gO1))
      *ZN(gI2,1),0) + IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1)
      ))*ZN(gI2,3)),0);

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
   
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Sin(ThetaW
      ())*ZUR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(ZUL(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Sin(
      ThetaW()),0);

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

double CLASSNAME::CpbarFvChiSvPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvChiSvPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZV(
      gI1,gO1))*ZN(gI2,0),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZV(gI1,gO1))
      *ZN(gI2,1),0);

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
      ),0) + SUM(j1,0,2,Conj(Ye(j1,gO1))*Conj(ZE(gI2,3 + j1)))*UM(gI1,1);

   return result;
}

double CLASSNAME::CpbarFvFeconjHpmPL(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjHpmPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(Ye(j1,gO1))*ZER(gI2,j1))*ZP
      (gI1,0);

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
      ))*UP(gI2,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(ZEL(gI1,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)))*ZA(
      gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI1,j1))*ZEL(gO1,j2))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(gI2
      ,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Ye(j1,j2))*ZER(gI2,j1))*ZEL(gO1,j2))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,gI2))*ZP
      (gI1,0);

   return result;
}

double CLASSNAME::CpbarFeFvHpmPR(int , int , int ) const
{
   
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -1.0954451150103321*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1))) - Conj(ZN(gI2,2))*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(gI1,
      j1))*ZEL(gO1,j1))*(0.7745966692414834*g1*ZN(gI2,0) + g2*ZN(gI2,1)) - SUM(j2,
      0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1)))*ZEL(gO1,j2))*ZN(gI2,2);

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

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(ZDL(gI1,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)))*ZA(
      gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI1,j1))*ZDL(gO1,j2))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZDL(gI2
      ,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZDL(gO1,j2))*ZH(gI1,0);

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
   
   const std::complex<double> result = -0.3651483716701107*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1))) - Conj(ZN(gI2,2))*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(ZD(gI1,j1))*ZDL(gO1,j1))*(-
      0.18257418583505536*g1*ZN(gI2,0) + 0.7071067811865475*g2*ZN(gI2,1)) - SUM(j2
      ,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1)))*ZDL(gO1,j2))*ZN(gI2,2)
      ;

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

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,-0.7071067811865475
      )*SUM(j2,0,2,Conj(ZUL(gI1,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)))*ZA(
      gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gO1, int gI1, int gI2) const
{
   
   const std::complex<double> result = std::complex<double>(0.,0.7071067811865475)
      *SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI1,j1))*ZUL(gO1,j2))*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZUL(gI2
      ,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)))*ZH(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZUL(gO1,j2))*ZH(gI1,1);

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

std::complex<double> CLASSNAME::CpbarFuChiSuPL(int gO2, int gI2, int gI1) const
{
   
   const std::complex<double> result = 0.7302967433402214*g1*Conj(ZN(gI2,0))*SUM(
      j1,0,2,Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1))) - Conj(ZN(gI2,3))*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPR(int gO1, int gI2, int gI1) const
{
   
   const std::complex<double> result = SUM(j1,0,2,Conj(ZU(gI1,j1))*ZUL(gO1,j1))*(-
      0.18257418583505536*g1*ZN(gI2,0) - 0.7071067811865475*g2*ZN(gI2,1)) - SUM(j2
      ,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1)))*ZUL(gO1,j2))*ZN(gI2,3)
      ;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSdconjUSdVZVZ(gO1,gO2);
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSdconjUSd(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSdconjUSd(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSdSvconjUSdconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpChaFuconjUSdPL(gI2,gI1,gO2))*
      CpChaFuconjUSdPL(gI2,gI1,gO1) + Conj(CpChaFuconjUSdPR(gI2,gI1,gO2))*
      CpChaFuconjUSdPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MCha(
      gI2)))*(Conj(CpChaFuconjUSdPR(gI2,gI1,gO2))*CpChaFuconjUSdPL(gI2,gI1,gO1) +
      Conj(CpChaFuconjUSdPL(gI2,gI1,gO2))*CpChaFuconjUSdPR(gI2,gI1,gO1))*MCha(gI2)
      ));
   result += SUM(gI1,0,2,SUM(gI2,0,4,(Conj(CpChiFdconjUSdPL(gI2,gI1,gO2))*
      CpChiFdconjUSdPL(gI2,gI1,gO1) + Conj(CpChiFdconjUSdPR(gI2,gI1,gO2))*
      CpChiFdconjUSdPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,4,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFdconjUSdPR(gI2,gI1,gO2))*CpChiFdconjUSdPL(gI2,gI1,gO1) +
      Conj(CpChiFdconjUSdPL(gI2,gI1,gO2))*CpChiFdconjUSdPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSdconjHpmconjUSd(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUSdSeconjUSdconjSe(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhSdconjUSd(gI2,gI1,gO2))*CpAhSdconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSdconjUSd(gI2,gI1,gO2))*CphhSdconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MHpm(gI2)))*Conj(
      CpHpmSuconjUSd(gI2,gI1,gO2))*CpHpmSuconjUSd(gI2,gI1,gO1)));
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
   result += 2*A0(Sqr(MVZ))*CpUSvconjUSvVZVZ(gO1,gO2);
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSvconjUSv(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSvconjUSv(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvUSvconjSvconjUSv(gI1,gO1,gI1,gO2))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhSvconjUSv(gI2,gI1,gO2))*CphhSvconjUSv(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarChaFeconjUSvPL(gI1,gI2,gO2))*
      CpbarChaFeconjUSvPL(gI1,gI2,gO1) + Conj(CpbarChaFeconjUSvPR(gI1,gI2,gO2))*
      CpbarChaFeconjUSvPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFe(gI2)))));
   result += -2*SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFe
      (gI2)))*(Conj(CpbarChaFeconjUSvPR(gI1,gI2,gO2))*CpbarChaFeconjUSvPL(gI1,gI2,
      gO1) + Conj(CpbarChaFeconjUSvPL(gI1,gI2,gO2))*CpbarChaFeconjUSvPR(gI1,gI2,
      gO1))*MFe(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,4,(Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*
      CpChiFvconjUSvPL(gI2,gI1,gO1) + Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*
      CpChiFvconjUSvPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,4,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*CpChiFvconjUSvPL(gI2,gI1,gO1) +
      Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*CpChiFvconjUSvPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSvconjHpmconjUSv(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSe(gI2)))*Conj(
      CpSeconjHpmconjUSv(gI2,gI1,gO2))*CpSeconjHpmconjUSv(gI2,gI1,gO1)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSvconjSdconjUSv(gI1,gO1,gI1,gO2
      ));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSvconjSeconjUSv(gI1,gO1,gI1,gO2))
      ;
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuUSvconjSuconjUSv(gI1,gO1,gI1,gO2
      ));
   result += SUM(gI2,0,2,Conj(CpSvconjUSvVZ(gI2,gO2))*CpSvconjUSvVZ(gI2,gO1)*F0(
      Sqr(p),Sqr(MSv(gI2)),Sqr(MVZ)));
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
   result += 2*A0(Sqr(MVZ))*CpUSuconjUSuVZVZ(gO1,gO2);
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSuconjUSu(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSuconjUSu(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSuSvconjUSuconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarChaFdconjUSuPL(gI1,gI2,gO2))*
      CpbarChaFdconjUSuPL(gI1,gI2,gO1) + Conj(CpbarChaFdconjUSuPR(gI1,gI2,gO2))*
      CpbarChaFdconjUSuPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFd(gI2)))));
   result += -2*SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MFd
      (gI2)))*(Conj(CpbarChaFdconjUSuPR(gI1,gI2,gO2))*CpbarChaFdconjUSuPL(gI1,gI2,
      gO1) + Conj(CpbarChaFdconjUSuPL(gI1,gI2,gO2))*CpbarChaFdconjUSuPR(gI1,gI2,
      gO1))*MFd(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,4,(Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*
      CpChiFuconjUSuPL(gI2,gI1,gO1) + Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*
      CpChiFuconjUSuPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,4,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*CpChiFuconjUSuPL(gI2,gI1,gO1) +
      Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*CpChiFuconjUSuPR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSuconjHpmconjUSu(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSd(gI2)))*Conj(
      CpSdconjHpmconjUSu(gI2,gI1,gO2))*CpSdconjHpmconjUSu(gI2,gI1,gO1)));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSuconjSeconjUSu(gI1,gO1,gI1,gO2))
      ;
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1))
      ;
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
   result += 2*A0(Sqr(MVZ))*CpUSeconjUSeVZVZ(gO1,gO2);
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUSeconjUSe(gI1,gI1,gO1,gO2))
      ;
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUSeconjUSe(gI1,gI1,gO1,gO2))
      ;
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSeSvconjUSeconjSv(gO1,gI1,gO2,gI1))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpChaFvconjUSePL(gI2,gI1,gO2))*
      CpChaFvconjUSePL(gI2,gI1,gO1) + Conj(CpChaFvconjUSePR(gI2,gI1,gO2))*
      CpChaFvconjUSePR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MCha(
      gI2)))*(Conj(CpChaFvconjUSePR(gI2,gI1,gO2))*CpChaFvconjUSePL(gI2,gI1,gO1) +
      Conj(CpChaFvconjUSePL(gI2,gI1,gO2))*CpChaFvconjUSePR(gI2,gI1,gO1))*MCha(gI2)
      ));
   result += SUM(gI1,0,2,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MHpm(gI2)))*Conj(
      CpHpmSvconjUSe(gI2,gI1,gO2))*CpHpmSvconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,4,(Conj(CpChiFeconjUSePL(gI2,gI1,gO2))*
      CpChiFeconjUSePL(gI2,gI1,gO1) + Conj(CpChiFeconjUSePR(gI2,gI1,gO2))*
      CpChiFeconjUSePR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,4,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiFeconjUSePR(gI2,gI1,gO2))*CpChiFeconjUSePL(gI2,gI1,gO1) +
      Conj(CpChiFeconjUSePL(gI2,gI1,gO2))*CpChiFeconjUSePR(gI2,gI1,gO1))*MChi(gI2)
      ));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSeconjHpmconjUSe(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSeconjSdconjUSe(gI1,gO1,gI1,gO2
      ));
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
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*Conj(CpUhhconjVWmVWm(gO2))*
      CpUhhconjVWmVWm(gO1);
   result += 4*A0(Sqr(MVWm))*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUhhUhhVZVZ(gO1,gO2);
   result += 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
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
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarChaChaUhhPL(gI1,gI2,gO2))*
      CpbarChaChaUhhPL(gI1,gI2,gO1) + Conj(CpbarChaChaUhhPR(gI1,gI2,gO2))*
      CpbarChaChaUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUhhPL(gI1,gI2,gO2))*
      CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*
      CpbarFdFdUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUhhPL(gI1,gI2,gO2))*
      CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*
      CpbarFeFeUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUhhPL(gI1,gI2,gO2))*
      CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*
      CpbarFuFuUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -2*SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpbarChaChaUhhPR(gI1,gI2,gO2))*CpbarChaChaUhhPL(gI1,gI2,
      gO1) + Conj(CpbarChaChaUhhPL(gI1,gI2,gO2))*CpbarChaChaUhhPR(gI1,gI2,gO1))*
      MCha(gI2)));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUhhPL(gI1,gI2,gO2))*CpbarFdFdUhhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUhhPL(gI1,gI2,gO2))*CpbarFeFeUhhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUhhPL(gI1,gI2,gO2))*CpbarFuFuUhhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpUhhUhhHpmconjHpm(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*Conj
      (CpUhhHpmconjHpm(gO2,gI2,gI1))*CpUhhHpmconjHpm(gO1,gI2,gI1)));
   result += 0.5*SUM(gI1,0,4,SUM(gI2,0,4,(Conj(CpChiChiUhhPL(gI1,gI2,gO2))*
      CpChiChiUhhPL(gI1,gI2,gO1) + Conj(CpChiChiUhhPR(gI1,gI2,gO2))*CpChiChiUhhPR(
      gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -SUM(gI1,0,4,MChi(gI1)*SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiChiUhhPR(gI1,gI2,gO2))*CpChiChiUhhPL(gI1,gI2,gO1) + Conj(
      CpChiChiUhhPL(gI1,gI2,gO2))*CpChiChiUhhPR(gI1,gI2,gO1))*MChi(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhUhhSdconjSd(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhUhhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhUhhSuconjSu(gO1,gO2,gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*Conj
      (CpUhhSdconjSd(gO2,gI2,gI1))*CpUhhSdconjSd(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*Conj(
      CpUhhSeconjSe(gO2,gI2,gI1))*CpUhhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*Conj
      (CpUhhSuconjSu(gO2,gI2,gI1))*CpUhhSuconjSu(gO1,gI2,gI1)));
   result += SUM(gI2,0,2,Conj(CpAhUhhVZ(gI2,gO2))*CpAhUhhVZ(gI2,gO1)*F0(Sqr(p),Sqr
      (MAh(gI2)),Sqr(MVZ)));
   result += 2*SUM(gI2,0,3,Conj(CpUhhHpmconjVWm(gO2,gI2))*CpUhhHpmconjVWm(gO1,gI2)
      *F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVWm)));

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
   result += 2*A0(Sqr(MVZ))*CpUAhUAhVZVZ(gO1,gO2);
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUAhUAh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUAhUAhSvconjSv(gO1,gO2,gI1,gI1));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhAhUAh(gI1,gI2,gO2))*CpAhAhUAh(gI1,gI2,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhUAhhh(gI2,gO2,gI1))*CpAhUAhhh(gI2,gO1,gI1)));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CpUAhhhhh(gO2,gI1,gI2))*CpUAhhhhh(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarChaChaUAhPL(gI1,gI2,gO2))*
      CpbarChaChaUAhPL(gI1,gI2,gO1) + Conj(CpbarChaChaUAhPR(gI1,gI2,gO2))*
      CpbarChaChaUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUAhPL(gI1,gI2,gO2))*
      CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*
      CpbarFdFdUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUAhPL(gI1,gI2,gO2))*
      CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*
      CpbarFeFeUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUAhPL(gI1,gI2,gO2))*
      CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*
      CpbarFuFuUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -2*SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpbarChaChaUAhPR(gI1,gI2,gO2))*CpbarChaChaUAhPL(gI1,gI2,
      gO1) + Conj(CpbarChaChaUAhPL(gI1,gI2,gO2))*CpbarChaChaUAhPR(gI1,gI2,gO1))*
      MCha(gI2)));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFdFdUAhPL(gI1,gI2,gO2))*CpbarFdFdUAhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFeFeUAhPL(gI1,gI2,gO2))*CpbarFeFeUAhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj
      (CpbarFuFuUAhPL(gI1,gI2,gO2))*CpbarFuFuUAhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpUAhUAhHpmconjHpm(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*Conj
      (CpUAhHpmconjHpm(gO2,gI2,gI1))*CpUAhHpmconjHpm(gO1,gI2,gI1)));
   result += 0.5*SUM(gI1,0,4,SUM(gI2,0,4,(Conj(CpChiChiUAhPL(gI1,gI2,gO2))*
      CpChiChiUAhPL(gI1,gI2,gO1) + Conj(CpChiChiUAhPR(gI1,gI2,gO2))*CpChiChiUAhPR(
      gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -SUM(gI1,0,4,MChi(gI1)*SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(
      gI2)))*(Conj(CpChiChiUAhPR(gI1,gI2,gO2))*CpChiChiUAhPL(gI1,gI2,gO1) + Conj(
      CpChiChiUAhPL(gI1,gI2,gO2))*CpChiChiUAhPR(gI1,gI2,gO1))*MChi(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUAhUAhSdconjSd(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUAhUAhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUAhUAhSuconjSu(gO1,gO2,gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*Conj
      (CpUAhSdconjSd(gO2,gI2,gI1))*CpUAhSdconjSd(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*Conj(
      CpUAhSeconjSe(gO2,gI2,gI1))*CpUAhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*Conj
      (CpUAhSuconjSu(gO2,gI2,gI1))*CpUAhSuconjSu(gO1,gI2,gI1)));
   result += SUM(gI2,0,2,Conj(CpUAhhhVZ(gO2,gI2))*CpUAhhhVZ(gO1,gI2)*F0(Sqr(p),Sqr
      (Mhh(gI2)),Sqr(MVZ)));
   result += 2*SUM(gI2,0,3,Conj(CpUAhHpmconjVWm(gO2,gI2))*CpUAhHpmconjVWm(gO1,gI2)
      *F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVWm)));

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
   result += 4*B0(Sqr(p),0,Sqr(MVWm))*Conj(CpconjUHpmVPVWm(gO2))*CpconjUHpmVPVWm(
      gO1);
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVZ))*Conj(CpconjUHpmVWmVZ(gO2))*
      CpconjUHpmVWmVZ(gO1);
   result += 4*A0(Sqr(MVWm))*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUHpmconjUHpmVZVZ(gO1,gO2);
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
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUHpmconjHpmconjUHpm(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MAh(gI2)))*Conj(
      CpAhHpmconjUHpm(gI2,gI1,gO2))*CpAhHpmconjUHpm(gI2,gI1,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(Mhh(gI2)))*Conj(
      CphhHpmconjUHpm(gI2,gI1,gO2))*CphhHpmconjUHpm(gI2,gI1,gO1)));
   result += SUM(gI1,0,4,SUM(gI2,0,2,(Conj(CpChiChaconjUHpmPL(gI1,gI2,gO2))*
      CpChiChaconjUHpmPL(gI1,gI2,gO1) + Conj(CpChiChaconjUHpmPR(gI1,gI2,gO2))*
      CpChiChaconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2)))));
   result += -2*SUM(gI1,0,4,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha(gI2)))*(Conj(CpChiChaconjUHpmPR(gI1,gI2,gO2))*CpChiChaconjUHpmPL(gI1,
      gI2,gO1) + Conj(CpChiChaconjUHpmPL(gI1,gI2,gO2))*CpChiChaconjUHpmPR(gI1,gI2,
      gO1))*MCha(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUHpmSdconjUHpmconjSd(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUHpmSeconjUHpmconjSe(gO1,gI1,gO2,gI1
      ));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUHpmSuconjUHpmconjSu(gO1,gI1,gO2,
      gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSd(gI2)))*Conj
      (CpSdconjUHpmconjSu(gI2,gO2,gI1))*CpSdconjUHpmconjSu(gI2,gO1,gI1)));
   result += SUM(gI2,0,2,Conj(CpAhconjUHpmVWm(gI2,gO2))*CpAhconjUHpmVWm(gI2,gO1)*
      F0(Sqr(p),Sqr(MAh(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,2,Conj(CphhconjUHpmVWm(gI2,gO2))*CphhconjUHpmVWm(gI2,gO1)*
      F0(Sqr(p),Sqr(Mhh(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,3,Conj(CpHpmconjUHpmVP(gI2,gO2))*CpHpmconjUHpmVP(gI2,gO1)*
      F0(Sqr(p),Sqr(MHpm(gI2)),0));
   result += SUM(gI2,0,3,Conj(CpHpmconjUHpmVZ(gI2,gO2))*CpHpmconjUHpmVZ(gI2,gO1)*
      F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVZ)));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_Hpm_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = i; k < 4; k++)
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
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarChaChaVPPL(gI1,gI2)) + AbsSqr(
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
   result += SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVPVP(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpHpmconjHpmVP(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVPVP(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVPVP(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVPVP(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVP(gI2,gI1))*B00(Sqr(p),
      Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += 2*SUM(gI2,0,3,AbsSqr(CpHpmconjVWmVP(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(
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
   result += 0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhVZVZ(gI1,gI1));
   result += 0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhVZVZ(gI1,gI1));
   result += SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvconjSvVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpAhhhVZ(gI2,gI1))*B00(Sqr(p),Sqr(
      MAh(gI2)),Sqr(Mhh(gI1)))));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpSvconjSvVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSv(gI1)),Sqr(MSv(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarChaChaVZPL(gI1,gI2)) + AbsSqr(
      CpbarChaChaVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2))) + 4*B0(
      Sqr(p),Sqr(MCha(gI1)),Sqr(MCha(gI2)))*MCha(gI1)*MCha(gI2)*Re(Conj(
      CpbarChaChaVZPL(gI1,gI2))*CpbarChaChaVZPR(gI1,gI2))));
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
   result += SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpHpmconjHpmVZ(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += 0.5*SUM(gI1,0,4,SUM(gI2,0,4,(AbsSqr(CpChiChiVZPL(gI1,gI2)) + AbsSqr(
      CpChiChiVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2))) + 4*B0(Sqr(
      p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))*MChi(gI1)*MChi(gI2)*Re(Conj(CpChiChiVZPL(
      gI1,gI2))*CpChiChiVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVZVZ(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVZVZ(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVZVZ(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += SUM(gI2,0,2,AbsSqr(CphhVZVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh(gI2))));
   result += 2*SUM(gI2,0,3,AbsSqr(CpHpmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(
      MHpm(gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_VWm_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargPgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVP));
   result += AbsSqr(CpbargWmCgPconjVWm())*B00(Sqr(p),Sqr(MVP),Sqr(MVWm));
   result += AbsSqr(CpbargWmCgZconjVWm())*B00(Sqr(p),Sqr(MVZ),Sqr(MVWm));
   result += AbsSqr(CpbargZgWmconjVWm())*B00(Sqr(p),Sqr(MVWm),Sqr(MVZ));
   result += -(A0(Sqr(MVWm))*(CpconjVWmconjVWmVWmVWm1() + 4*
      CpconjVWmconjVWmVWmVWm2() + CpconjVWmconjVWmVWmVWm3()));
   result += 0;
   result += -0.5*A0(Sqr(MVZ))*(4*CpconjVWmVWmVZVZ1() + CpconjVWmVWmVZVZ2() +
      CpconjVWmVWmVZVZ3());
   result += -(AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 10*B00(Sqr(p),Sqr(MVWm),0
      ) + B0(Sqr(p),Sqr(MVWm),0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += AbsSqr(CpconjVWmVWmVZ())*(-A0(Sqr(MVWm)) - A0(Sqr(MVZ)) - 10*B00(Sqr(
      p),Sqr(MVZ),Sqr(MVWm)) - B0(Sqr(p),Sqr(MVZ),Sqr(MVWm))*(Sqr(MVWm) + Sqr(MVZ)
      + 4*Sqr(p)));
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
   result += SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmconjVWmVWm(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,2,AbsSqr(CpAhHpmconjVWm(gI2,gI1))*B00(Sqr(p)
      ,Sqr(MAh(gI2)),Sqr(MHpm(gI1)))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,2,AbsSqr(CphhHpmconjVWm(gI2,gI1))*B00(Sqr(p)
      ,Sqr(Mhh(gI2)),Sqr(MHpm(gI1)))));
   result += SUM(gI1,0,4,SUM(gI2,0,2,(AbsSqr(CpChiChaconjVWmPL(gI1,gI2)) + AbsSqr(
      CpChiChaconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2))) + 4*B0
      (Sqr(p),Sqr(MChi(gI1)),Sqr(MCha(gI2)))*MCha(gI2)*MChi(gI1)*Re(Conj(
      CpChiChaconjVWmPL(gI1,gI2))*CpChiChaconjVWmPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuconjVWmVWm(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSuconjVWm(gI2,gI1))*B00(
      Sqr(p),Sqr(MSd(gI2)),Sqr(MSu(gI1)))));
   result += SUM(gI2,0,2,AbsSqr(CphhconjVWmVWm(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(Mhh(
      gI2))));
   result += SUM(gI2,0,3,AbsSqr(CpHpmconjVWmVP(gI2))*B0(Sqr(p),0,Sqr(MHpm(gI2))));
   result += SUM(gI2,0,3,AbsSqr(CpHpmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(MHpm(
      gI2))));

   return result * oneOver16PiSqr;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -4*SUM(gI1,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPL(gI1,gO2))*CpbarChaUChiVWmPR(gI1,gO1)*MCha(gI1));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)
      ))*Conj(CpbarFvUChiSvPL(gI1,gO2,gI2))*CpbarFvUChiSvPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*Conj(
      CpUChiFvconjSvPL(gO2,gI2,gI1))*CpUChiFvconjSvPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(
      gI2)))*Conj(CpbarChaUChiHpmPL(gI1,gO2,gI2))*CpbarChaUChiHpmPR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,2,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpChiUChihhPL(gI2,gO2,gI1))*CpChiUChihhPR(gI2,gO1,gI1)*MChi(gI2)));
   result += 3*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarFdUChiSdPL(gI1,gO2,gI2))*CpbarFdUChiSdPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarFeUChiSePL(gI1,gO2,gI2))*CpbarFeUChiSePR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(
      gI2)))*Conj(CpbarFuUChiSuPL(gI1,gO2,gI2))*CpbarFuUChiSuPR(gI1,gO1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))*Conj
      (CpUChiChaconjHpmPL(gO2,gI2,gI1))*CpUChiChaconjHpmPR(gO1,gI2,gI1)*MCha(gI2))
      );
   result += SUM(gI1,0,4,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(
      gI2)))*Conj(CpChiUChiAhPL(gI1,gO2,gI2))*CpChiUChiAhPR(gI1,gO1,gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*Conj
      (CpUChiFdconjSdPL(gO2,gI2,gI1))*CpUChiFdconjSdPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*Conj(
      CpUChiFeconjSePL(gO2,gI2,gI1))*CpUChiFeconjSePR(gO1,gI2,gI1)*MFe(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*Conj
      (CpUChiFuconjSuPL(gO2,gI2,gI1))*CpUChiFuconjSuPR(gO1,gI2,gI1)*MFu(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MVWm))*Conj(
      CpUChiChaconjVWmPR(gO2,gI2))*CpUChiChaconjVWmPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(CpChiUChiVZPL(
      gI2,gO2))*CpChiUChiVZPR(gI2,gO1)*MChi(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,5,5> CLASSNAME::self_energy_Chi_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,5,5> self_energy;

   for (int i = 0; i < 5; i++)
      for (int k = 0; k < 5; k++)
         self_energy(i, k) = self_energy_Chi_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -SUM(gI1,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPR(gI1,gO2))*CpbarChaUChiVWmPR(gI1,gO1));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpbarFvUChiSvPR(gI1,gO2,gI2))*CpbarFvUChiSvPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*
      Conj(CpUChiFvconjSvPR(gO2,gI2,gI1))*CpUChiFvconjSvPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(gI2)))
      *Conj(CpbarChaUChiHpmPR(gI1,gO2,gI2))*CpbarChaUChiHpmPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpChiUChihhPR(gI2,gO2,gI1))*CpChiUChihhPR(gI2,gO1,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarFdUChiSdPR(gI1,gO2,gI2))*CpbarFdUChiSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarFeUChiSePR(gI1,gO2,gI2))*CpbarFeUChiSePR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpbarFuUChiSuPR(gI1,gO2,gI2))*CpbarFuUChiSuPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpUChiChaconjHpmPR(gO2,gI2,gI1))*CpUChiChaconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,4,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)))*
      Conj(CpChiUChiAhPR(gI1,gO2,gI2))*CpChiUChiAhPR(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpUChiFdconjSdPR(gO2,gI2,gI1))*CpUChiFdconjSdPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpUChiFeconjSePR(gO2,gI2,gI1))*CpUChiFeconjSePR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpUChiFuconjSuPR(gO2,gI2,gI1))*CpUChiFuconjSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVWm))*Conj(
      CpUChiChaconjVWmPL(gO2,gI2))*CpUChiChaconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(CpChiUChiVZPR(
      gI2,gO2))*CpChiUChiVZPR(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,5,5> CLASSNAME::self_energy_Chi_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,5,5> self_energy;

   for (int i = 0; i < 5; i++)
      for (int k = 0; k < 5; k++)
         self_energy(i, k) = self_energy_Chi_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -SUM(gI1,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MVWm))*Conj(
      CpbarChaUChiVWmPL(gI1,gO2))*CpbarChaUChiVWmPL(gI1,gO1));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpbarFvUChiSvPL(gI1,gO2,gI2))*CpbarFvUChiSvPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSv(gI1)))*
      Conj(CpUChiFvconjSvPL(gO2,gI2,gI1))*CpUChiFvconjSvPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MHpm(gI2)))
      *Conj(CpbarChaUChiHpmPL(gI1,gO2,gI2))*CpbarChaUChiHpmPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpChiUChihhPL(gI2,gO2,gI1))*CpChiUChihhPL(gI2,gO1,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarFdUChiSdPL(gI1,gO2,gI2))*CpbarFdUChiSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarFeUChiSePL(gI1,gO2,gI2))*CpbarFeUChiSePL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpbarFuUChiSuPL(gI1,gO2,gI2))*CpbarFuUChiSuPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpUChiChaconjHpmPL(gO2,gI2,gI1))*CpUChiChaconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,4,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)))*
      Conj(CpChiUChiAhPL(gI1,gO2,gI2))*CpChiUChiAhPL(gI1,gO1,gI2)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpUChiFdconjSdPL(gO2,gI2,gI1))*CpUChiFdconjSdPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpUChiFeconjSePL(gO2,gI2,gI1))*CpUChiFeconjSePL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpUChiFuconjSuPL(gO2,gI2,gI1))*CpUChiFuconjSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVWm))*Conj(
      CpUChiChaconjVWmPR(gO2,gI2))*CpUChiChaconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(CpChiUChiVZPL(
      gI2,gO2))*CpChiUChiVZPL(gI2,gO1));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,5,5> CLASSNAME::self_energy_Chi_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,5,5> self_energy;

   for (int i = 0; i < 5; i++)
      for (int k = 0; k < 5; k++)
         self_energy(i, k) = self_energy_Chi_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUChaChaAhPL(gO2,gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUChaChahhPL(gO2,gI2,gI1))*CpbarUChaChahhPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUChaFeconjSvPL(gO2,gI2,gI1))*CpbarUChaFeconjSvPR(gO1,gI2,gI1)*MFe(gI2))
      );
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarUChabarFuSdPL(gO2,gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2)
      ));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarUChabarFvSePL(gO2,gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))*Conj
      (CpbarUChaChiHpmPL(gO2,gI2,gI1))*CpbarUChaChiHpmPR(gO1,gI2,gI1)*MChi(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*Conj
      (CpbarUChaFdconjSuPL(gO2,gI2,gI1))*CpbarUChaFdconjSuPR(gO1,gI2,gI1)*MFd(gI2)
      ));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),0)*Conj(CpbarUChaChaVPPR(gO2,
      gI2))*CpbarUChaChaVPPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZ))*Conj(
      CpbarUChaChaVZPR(gO2,gI2))*CpbarUChaChaVZPL(gO1,gI2)*MCha(gI2));
   result += -4*SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaChiVWmPR(gO2,gI2))*CpbarUChaChiVWmPL(gO1,gI2)*MChi(gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Cha_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Cha_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUChaChaAhPR(gO2,gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChaChahhPR(gO2,gI2,gI1))*CpbarUChaChahhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUChaFeconjSvPR(gO2,gI2,gI1))*CpbarUChaFeconjSvPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarUChabarFuSdPR(gO2,gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarUChabarFvSePR(gO2,gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpbarUChaChiHpmPR(gO2,gI2,gI1))*CpbarUChaChiHpmPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChaFdconjSuPR(gO2,gI2,gI1))*CpbarUChaFdconjSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),0)*Conj(CpbarUChaChaVPPL(gO2,
      gI2))*CpbarUChaChaVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZ))*Conj(CpbarUChaChaVZPL
      (gO2,gI2))*CpbarUChaChaVZPL(gO1,gI2));
   result += -SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaChiVWmPL(gO2,gI2))*CpbarUChaChiVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Cha_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Cha_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUChaChaAhPL(gO2,gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChaChahhPL(gO2,gI2,gI1))*CpbarUChaChahhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUChaFeconjSvPL(gO2,gI2,gI1))*CpbarUChaFeconjSvPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarUChabarFuSdPL(gO2,gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarUChabarFvSePL(gO2,gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))
      *Conj(CpbarUChaChiHpmPL(gO2,gI2,gI1))*CpbarUChaChiHpmPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChaFdconjSuPL(gO2,gI2,gI1))*CpbarUChaFdconjSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),0)*Conj(CpbarUChaChaVPPR(gO2,
      gI2))*CpbarUChaChaVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MVZ))*Conj(CpbarUChaChaVZPR
      (gO2,gI2))*CpbarUChaChaVZPR(gO1,gI2));
   result += -SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUChaChiVWmPR(gO2,gI2))*CpbarUChaChiVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;
}

Eigen::Matrix<std::complex<double>,3,3> CLASSNAME::self_energy_Cha_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Cha_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarUFeChaSvPL(gO2,gI2,gI1))*CpbarUFeChaSvPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFeFvHpmPL(gO2,gI2,gI1))*CpbarUFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*Conj(
      CpbarUFeChiSePL(gO2,gI2,gI1))*CpbarUFeChiSePR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPR(gO2,gI2
      ))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPR(
      gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFeChaSvPR(gO2,gI2,gI1))*CpbarUFeChaSvPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPR(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPR(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFeFvHpmPR(gO2,gI2,gI1))*CpbarUFeFvHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUFeChiSePR(gO2,gI2,gI1))*CpbarUFeChiSePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPL(gO2,gI2))
      *CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPL(
      gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarUFeChaSvPL(gO2,gI2,gI1))*CpbarUFeChaSvPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFeFvHpmPL(gO2,gI2,gI1))*CpbarUFeFvHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUFeChiSePL(gO2,gI2,gI1))*CpbarUFeChiSePL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPR(gO2,gI2))
      *CpbarUFeFeVPPR(gO1,gI2));
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

   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFdFuHpmPL(gO2,gI2,gI1))*CpbarUFdFuHpmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1))
      )*Conj(CpbarUFdGluSdPL(gO2,gI1))*CpbarUFdGluSdPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarUFdChaSuPL(gO2,gI2,gI1))*CpbarUFdChaSuPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*Conj(
      CpbarUFdChiSdPL(gO2,gI2,gI1))*CpbarUFdChiSdPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPR(gO2,gI2
      ))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPR(
      gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFdFdAhPR(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPR(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFdFuHpmPR(gO2,gI2,gI1))*CpbarUFdFuHpmPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarUFdGluSdPR(gO2,gI1))*CpbarUFdGluSdPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFdChaSuPR(gO2,gI2,gI1))*CpbarUFdChaSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFdChiSdPR(gO2,gI2,gI1))*CpbarUFdChiSdPR(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPL(gO2,gI2))
      *CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPL(
      gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFdFuHpmPL(gO2,gI2,gI1))*CpbarUFdFuHpmPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarUFdGluSdPL(gO2,gI1))*CpbarUFdGluSdPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFdChaSuPL(gO2,gI2,gI1))*CpbarUFdChaSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFdChiSdPL(gO2,gI2,gI1))*CpbarUFdChiSdPL(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPR(gO2,gI2))
      *CpbarUFdFdVPPR(gO1,gI2));
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

   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)
      ));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)*MFd(gI2))
      );
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1))
      )*Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2
      ))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPR(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPR(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPR(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPR(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPL(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPL(gO2,gI2))
      *CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPL(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2))
      *CpbarUFuFuVPPR(gO1,gI2));
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

std::complex<double> CLASSNAME::self_energy_Fv_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarFvChiSvPL(gO2,gI2,gI1))*CpbarFvChiSvPR(gO1,gI2,gI1)*MChi(gI2)));
   result += SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(
      gI2)))*Conj(CpbarChabarFvSePL(gI1,gO2,gI2))*CpbarChabarFvSePR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFvFeconjHpmPL(gO2,gI2,gI1))*CpbarFvFeconjHpmPR(gO1,gI2,gI1)*MFe(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPR(
      gO2,gI2))*CpbarFvFvVZPL(gO1,gI2)*MFv(gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFvChiSvPR(gO2,gI2,gI1))*CpbarFvChiSvPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarChabarFvSePR(gI1,gO2,gI2))*CpbarChabarFvSePR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFvFeconjHpmPR(gO2,gI2,gI1))*CpbarFvFeconjHpmPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPL(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPL(gO2
      ,gI2))*CpbarFvFvVZPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFvChiSvPL(gO2,gI2,gI1))*CpbarFvChiSvPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSe(gI2)))*
      Conj(CpbarChabarFvSePL(gI1,gO2,gI2))*CpbarChabarFvSePL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFvFeconjHpmPL(gO2,gI2,gI1))*CpbarFvFeconjHpmPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPR(gO1,gI2));
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

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1_heavy_rotated(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*Conj(
      CpbarFeChaSvPL(gO2,gI2,gI1))*CpbarFeChaSvPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*Conj(
      CpbarFeChiSePL(gO2,gI2,gI1))*CpbarFeChiSePR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPR(
      gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFeChaSvPR(gO2,gI2,gI1))*CpbarFeChaSvPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPR(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPR(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFeFvHpmPR(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarFeChiSePR(gO2,gI2,gI1))*CpbarFeChiSePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPL(gO2
      ,gI2))*CpbarFeFeVZPL(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSv(gI1)))*
      Conj(CpbarFeChaSvPL(gO2,gI2,gI1))*CpbarFeChaSvPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarFeChiSePL(gO2,gI2,gI1))*CpbarFeChiSePL(gO1,gI2,gI1)));
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

   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFdFuHpmPL(gO2,gI2,gI1))*CpbarFdFuHpmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1))
      )*Conj(CpbarFdGluSdPL(gO2,gI1))*CpbarFdGluSdPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarFdChaSuPL(gO2,gI2,gI1))*CpbarFdChaSuPR(gO1,gI2,gI1)*MCha(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*Conj(
      CpbarFdChiSdPL(gO2,gI2,gI1))*CpbarFdChiSdPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPR(
      gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFdFdAhPR(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPR(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFdFuHpmPR(gO2,gI2,gI1))*CpbarFdFuHpmPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarFdGluSdPR(gO2,gI1))*CpbarFdGluSdPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFdChaSuPR(gO2,gI2,gI1))*CpbarFdChaSuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarFdChiSdPR(gO2,gI2,gI1))*CpbarFdChiSdPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPL(gO2
      ,gI2))*CpbarFdFdVZPL(gO1,gI2));
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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFdFuHpmPL(gO2,gI2,gI1))*CpbarFdFuHpmPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))*
      Conj(CpbarFdGluSdPL(gO2,gI1))*CpbarFdGluSdPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFdChaSuPL(gO2,gI2,gI1))*CpbarFdChaSuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarFdChiSdPL(gO2,gI2,gI1))*CpbarFdChiSdPL(gO1,gI2,gI1)));
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

   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarChabarFuSdPL(gI1,gO2,gI2))*CpbarChabarFuSdPR(gI1,gO1,gI2)))
      ;
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarFuFdconjHpmPL(gO2,gI2,gI1))*CpbarFuFdconjHpmPR(gO1,gI2,gI1)*MFd(gI2)));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1))
      )*Conj(CpbarFuGluSuPL(gO2,gI1))*CpbarFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarFuChiSuPL(gO2,gI2,gI1))*CpbarFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPR(gO2,gI2)
      )*CpbarFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPR(
      gO2,gI2))*CpbarFuFuVZPL(gO1,gI2)*MFu(gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFuFuAhPR(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPR(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarFuSdPR(gI1,gO2,gI2))*CpbarChabarFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFuFdconjHpmPR(gO2,gI2,gI1))*CpbarFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarFuGluSuPR(gO2,gI1))*CpbarFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFuChiSuPR(gO2,gI2,gI1))*CpbarFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPL(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPL(gO2,gI2))*
      CpbarFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPL(gO2
      ,gI2))*CpbarFuFuVZPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarFuSdPL(gI1,gO2,gI2))*CpbarChabarFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFuFdconjHpmPL(gO2,gI2,gI1))*CpbarFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarFuGluSuPL(gO2,gI1))*CpbarFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFuChiSuPL(gO2,gI2,gI1))*CpbarFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPR(gO2,gI2))*
      CpbarFuFuVPPR(gO1,gI2));
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

   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*Conj(
      CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,2,MCha(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(
      gI2)))*Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)
      ));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*Conj(
      CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)*MFd(gI2))
      );
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1))
      )*Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,4,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*Conj(
      CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2
      ))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPR(gI1,gO2,gI2))*CpbarChabarUFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPR(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPR(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPR(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPL(gO2,gI2))
      *CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPL(
      gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha(gI1)),Sqr(MSd(gI2)))*
      Conj(CpbarChabarUFuSdPL(gI1,gO2,gI2))*CpbarChabarUFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))*
      Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,4,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,gI2))
      *CpbarUFuFuVPPR(gO1,gI2));
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
   result += 4*A0(Sqr(MVWm))*CpUhhconjVWmVWm(gO1);
   result += 2*A0(Sqr(MVZ))*CpUhhVZVZ(gO1);
   result += -0.5*SUM(gI1,0,2,A0(Sqr(MAh(gI1)))*CpAhAhUhh(gI1,gI1,gO1));
   result += -0.5*SUM(gI1,0,2,A0(Sqr(Mhh(gI1)))*CphhhhUhh(gI1,gI1,gO1));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUhhSvconjSv(gO1,gI1,gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MCha(gI1)))*(CpbarChaChaUhhPL(gI1,gI1,gO1) +
      CpbarChaChaUhhPR(gI1,gI1,gO1))*MCha(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdUhhPL(gI1,gI1,gO1) +
      CpbarFdFdUhhPR(gI1,gI1,gO1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFeUhhPL(gI1,gI1,gO1) +
      CpbarFeFeUhhPR(gI1,gI1,gO1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuUhhPL(gI1,gI1,gO1) +
      CpbarFuFuUhhPR(gI1,gI1,gO1))*MFu(gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpUhhHpmconjHpm(gO1,gI1,gI1));
   result += SUM(gI1,0,4,A0(Sqr(MChi(gI1)))*(CpChiChiUhhPL(gI1,gI1,gO1) +
      CpChiChiUhhPR(gI1,gI1,gO1))*MChi(gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhSdconjSd(gO1,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhSeconjSe(gO1,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhSuconjSu(gO1,gI1,gI1));

   return result * oneOver16PiSqr;
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

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(TMSSM_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(TMSSM_info::VZ);
   }

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_running_tachyon(TMSSM_info::Sd))
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
         problems.flag_bad_mass(TMSSM_info::Sd, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(TMSSM_info::Sv))
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
         problems.flag_bad_mass(TMSSM_info::Sv, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(TMSSM_info::Su))
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
         problems.flag_bad_mass(TMSSM_info::Su, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(TMSSM_info::Se))
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
         problems.flag_bad_mass(TMSSM_info::Se, eigenvalue_error > precision *
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
   if (!force_output && problems.is_running_tachyon(TMSSM_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_hh());

      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_Mhh(es));
         Eigen::Matrix<double,3,3> self_energy = Re(self_energy_hh_1loop(p));
         const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZH,
               eigenvalue_error);
            problems.flag_bad_mass(TMSSM_info::hh, eigenvalue_error > precision
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
      problems.flag_no_pole_mass_convergence(TMSSM_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(TMSSM_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_running_tachyon(TMSSM_info::Ah))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   do {
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Ah());

      for (int es = 0; es < 3; ++es) {
         const double p = Abs(old_MAh(es));
         Eigen::Matrix<double,3,3> self_energy = Re(self_energy_Ah_1loop(p));
         const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA,
               eigenvalue_error);
            problems.flag_bad_mass(TMSSM_info::Ah, eigenvalue_error > precision
                * Abs(eigen_values(0)));
         #else

            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA);
         #endif
            normalize_to_interval(mix_ZA);

         PHYSICAL(MAh(es)) = SignedAbsSqrt(eigen_values(es));
         if (es == 1)
            PHYSICAL(ZA) = mix_ZA;
      }

      new_MAh = PHYSICAL(MAh);
      diff = MaxRelDiff(new_MAh, old_MAh);
      old_MAh = new_MAh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(TMSSM_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(TMSSM_info::Ah);
}

void CLASSNAME::calculate_MHpm_pole()
{
   if (!force_output && problems.is_running_tachyon(TMSSM_info::Hpm))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations();
   int iteration = 0;
   double diff = 0.0;
   decltype(MHpm) old_MHpm(MHpm), new_MHpm(MHpm);

   do {
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Hpm());

      for (int es = 0; es < 4; ++es) {
         const double p = Abs(old_MHpm(es));
         Eigen::Matrix<double,4,4> self_energy = Re(self_energy_Hpm_1loop(p));
         const Eigen::Matrix<double,4,4> M_loop(M_tree - self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP,
               eigenvalue_error);
            problems.flag_bad_mass(TMSSM_info::Hpm, eigenvalue_error >
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
      problems.flag_no_pole_mass_convergence(TMSSM_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(TMSSM_info::Hpm);
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,5,5> M_tree(get_mass_matrix_Chi());
   for (int es = 0; es < 5; ++es) {
      const double p = Abs(MChi(es));
      const Eigen::Matrix<double,5,5> self_energy_1  = Re(
         self_energy_Chi_1loop_1(p));
      const Eigen::Matrix<double,5,5> self_energy_PL = Re(
         self_energy_Chi_1loop_PL(p));
      const Eigen::Matrix<double,5,5> self_energy_PR = Re(
         self_energy_Chi_1loop_PR(p));
      const Eigen::Matrix<double,5,5> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,5,5> M_loop(M_tree + 0.5 * (delta_M + delta_M.
         transpose()));
      Eigen::Array<double,5,1> eigen_values;
      decltype(ZN) mix_ZN;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZN,
            eigenvalue_error);
         problems.flag_bad_mass(TMSSM_info::Chi, eigenvalue_error > precision *
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
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Cha());
   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MCha(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Cha_1loop_1(p));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Cha_1loop_PL(p));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Cha_1loop_PR(p));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(UM) mix_UM;
      decltype(UP) mix_UP;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_UM, mix_UP, eigenvalue_error);
      problems.flag_bad_mass(TMSSM_info::Cha, eigenvalue_error > precision *
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
      problems.flag_bad_mass(TMSSM_info::Fe, eigenvalue_error > precision * Abs
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
      problems.flag_bad_mass(TMSSM_info::Fd, eigenvalue_error > precision * Abs
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
      problems.flag_bad_mass(TMSSM_info::Fu, eigenvalue_error > precision * Abs
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
   if (!force_output && problems.is_running_tachyon(TMSSM_info::VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(TMSSM_info::VWm);
   }

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(TMSSM_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(TMSSM_info::VWm);
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(TMSSM_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(TMSSM_info::VZ);
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
      problems.flag_pole_tachyon(TMSSM_info::VZ);return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole) const
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(TMSSM_info::VWm);return m_pole;
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
