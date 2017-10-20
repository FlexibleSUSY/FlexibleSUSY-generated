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

// File generated at Fri 20 Oct 2017 08:49:46

/**
 * @file MRSSM_mass_eigenstates.cpp
 * @brief implementation of the MRSSM model class
 *
 * Contains the definition of the MRSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Fri 20 Oct 2017 08:49:46 with FlexibleSUSY
 * 2.0.1 (git commit: 5296739235bd0ef7020eda218da9c069270c3f45) and SARAH 4.12.0 .
 */

#include "MRSSM_mass_eigenstates.hpp"
#include "MRSSM_ewsb_solver_interface.hpp"
#include "eigen_utils.hpp"
#include "ewsb_solver.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "pv.hpp"
#include "raii.hpp"
#include "thread_pool.hpp"
#include "functors.hpp"

#include "config.h"

#ifdef ENABLE_TWO_SCALE_SOLVER
#include "MRSSM_two_scale_ewsb_solver.hpp"
#endif





#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>
#include <stdexcept>

namespace flexiblesusy {

#define CLASSNAME MRSSM_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model.get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model.get_##parameter()
#define EXTRAPARAMETER(parameter) model.get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS     loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION          loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS  loop_corrections.higgs_at_as_as
#define HIGGS_3LOOP_CORRECTION_AB_AS_AS  loop_corrections.higgs_ab_as_as
#define HIGGS_3LOOP_MDR_SCHEME           loop_corrections.higgs_3L_mdr_scheme
#define HIGGS_3LOOP_CORRECTION_AT_AT_AS  loop_corrections.higgs_at_at_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AT  loop_corrections.higgs_at_at_at

CLASSNAME::MRSSM_mass_eigenstates(const MRSSM_input_parameters& input_)
   : MRSSM_soft_parameters(input_)
#if defined(ENABLE_TWO_SCALE_SOLVER)
   , ewsb_solver(new MRSSM_ewsb_solver<Two_scale>())
#endif
{
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

const MRSSM_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MRSSM_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MRSSM_physical& physical_)
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

void CLASSNAME::set_ewsb_solver(const std::shared_ptr<MRSSM_ewsb_solver_interface>& solver)
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
   tadpole[3] = get_ewsb_eq_hh_4();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop(0));
      tadpole[1] -= Re(tadpole_hh_1loop(1));
      tadpole[2] -= Re(tadpole_hh_1loop(3));
      tadpole[3] -= Re(tadpole_hh_1loop(2));

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
   tadpole[2] /= vS;
   tadpole[3] /= vT;


   return tadpole;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = EWSB_solver::SUCCESS;

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_mT2 = mT2;
   const double old_mS2 = mS2;

   mHd2 = Re((0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*
      vd*AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS
      *Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) - 3*Cube(vd)*Sqr(g1) - 5*Cube(vd)*Sqr(g2) - 20*vd*AbsSqr(LamSD)*
      Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*
      Sqr(vu)))/vd);
   mHu2 = Re((0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40
      *vu*AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*
      vu*Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*
      vu*Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) - 3*Cube(vu)*Sqr(g1) - 5*Cube(vu)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr(vd) +
      5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu*AbsSqr(LamTU)*
      Sqr(vT)))/vu);
   mT2 = Re((0.125*(-32*vT*Sqr(MDWBT) - 2*g2*MDWBT*Sqr(vd) - 2*vT*AbsSqr(LamTD)
      *Sqr(vd) - 1.4142135623730951*LamTD*vS*Conj(LamSD)*Sqr(vd) - 2*MuD*Conj(
      LamTD)*Sqr(vd) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*Sqr(vd) - 2*g2*Conj
      (MDWBT)*Sqr(vd) - 2*LamTD*Conj(MuD)*Sqr(vd) + 2*g2*MDWBT*Sqr(vu) - 2*vT*
      AbsSqr(LamTU)*Sqr(vu) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*Sqr(vu) + 2*
      MuU*Conj(LamTU)*Sqr(vu) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*Sqr(vu) +
      2*g2*Conj(MDWBT)*Sqr(vu) + 2*LamTU*Conj(MuU)*Sqr(vu)))/vT);
   mS2 = Re((0.025*(-160*vS*Sqr(MDBS) + 7.745966692414834*g1*MDBS*Sqr(vd) - 20*
      vS*AbsSqr(LamSD)*Sqr(vd) - 14.142135623730951*MuD*Conj(LamSD)*Sqr(vd) -
      7.0710678118654755*LamTD*vT*Conj(LamSD)*Sqr(vd) - 7.0710678118654755*LamSD*
      vT*Conj(LamTD)*Sqr(vd) + 7.745966692414834*g1*Conj(MDBS)*Sqr(vd) -
      14.142135623730951*LamSD*Conj(MuD)*Sqr(vd) - 7.745966692414834*g1*MDBS*Sqr(
      vu) - 20*vS*AbsSqr(LamSU)*Sqr(vu) - 14.142135623730951*MuU*Conj(LamSU)*Sqr(
      vu) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*Sqr(vu) + 7.0710678118654755*
      LamSU*vT*Conj(LamTU)*Sqr(vu) - 7.745966692414834*g1*Conj(MDBS)*Sqr(vu) -
      14.142135623730951*LamSU*Conj(MuU)*Sqr(vu)))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(mT2) &&
      IsFinite(mS2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      mT2 = old_mT2;
      mS2 = old_mS2;
      error = EWSB_solver::FAIL;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   if (!ewsb_solver) {
      throw SetupError("MRSSM_mass_eigenstates::solve_ewsb_tree_level: "
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
      throw SetupError("MRSSM_mass_eigenstates::solve_ewsb_one_loop: "
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
      throw SetupError("MRSSM_mass_eigenstates::solve_ewsb: "
                       "no EWSB solver set");
   }

   VERBOSE_MSG("\t\tSolving MRSSM EWSB at " << ewsb_loop_order << "-loop order");

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
           "MRSSM\n"
           "========================================\n";
   MRSSM_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSRdp = " << MSRdp << '\n';
   ostr << "MSRum = " << MSRum << '\n';
   ostr << "MsigmaO = " << MsigmaO << '\n';
   ostr << "MphiO = " << MphiO << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MRh = " << MRh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha1 = " << MCha1.transpose() << '\n';
   ostr << "MCha2 = " << MCha2.transpose() << '\n';
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
   ostr << "ZHR = " << ZHR << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN1 = " << ZN1 << '\n';
   ostr << "ZN2 = " << ZN2 << '\n';
   ostr << "UM1 = " << UM1 << '\n';
   ostr << "UP1 = " << UP1 << '\n';
   ostr << "UM2 = " << UM2 << '\n';
   ostr << "UP2 = " << UP2 << '\n';
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
   return passarino_veltman::ReA0(m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB1(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB00(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB22(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReH0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReF0(p, m1, m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReG0(p, m1, m2, Sqr(get_scale()));
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto save_mHd2_raii = make_raii_save(mHd2);
   const auto save_mHu2_raii = make_raii_save(mHu2);
   const auto save_mS2_raii = make_raii_save(mS2);
   const auto save_mT2_raii = make_raii_save(mT2);

   solve_ewsb_tree_level_custom();

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MFu();
   calculate_MFd();
   calculate_MFe();
   calculate_MCha2();
   calculate_MCha1();
   calculate_MChi();
   calculate_MHpm();
   calculate_MRh();
   calculate_MAh();
   calculate_Mhh();
   calculate_MSe();
   calculate_MSu();
   calculate_MSv();
   calculate_MSd();
   calculate_MphiO();
   calculate_MsigmaO();
   calculate_MSRum();
   calculate_MSRdp();
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
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 24u));

   if (calculate_bsm_pole_masses) {
      tp.run_task([this] () { calculate_MAh_pole(); });
      tp.run_task([this] () { calculate_MCha1_pole(); });
      tp.run_task([this] () { calculate_MCha2_pole(); });
      tp.run_task([this] () { calculate_MChi_pole(); });
      tp.run_task([this] () { calculate_MGlu_pole(); });
      tp.run_task([this] () { calculate_Mhh_pole(); });
      tp.run_task([this] () { calculate_MHpm_pole(); });
      tp.run_task([this] () { calculate_MphiO_pole(); });
      tp.run_task([this] () { calculate_MRh_pole(); });
      tp.run_task([this] () { calculate_MSd_pole(); });
      tp.run_task([this] () { calculate_MSe_pole(); });
      tp.run_task([this] () { calculate_MsigmaO_pole(); });
      tp.run_task([this] () { calculate_MSRdp_pole(); });
      tp.run_task([this] () { calculate_MSRum_pole(); });
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
      calculate_MCha1_pole();
      calculate_MCha2_pole();
      calculate_MChi_pole();
      calculate_MGlu_pole();
      calculate_Mhh_pole();
      calculate_MHpm_pole();
      calculate_MphiO_pole();
      calculate_MRh_pole();
      calculate_MSd_pole();
      calculate_MSe_pole();
      calculate_MsigmaO_pole();
      calculate_MSRdp_pole();
      calculate_MSRum_pole();
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
   PHYSICAL(MSRdp) = MSRdp;
   PHYSICAL(MSRum) = MSRum;
   PHYSICAL(MsigmaO) = MsigmaO;
   PHYSICAL(MphiO) = MphiO;
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
   PHYSICAL(MRh) = MRh;
   PHYSICAL(ZHR) = ZHR;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN1) = ZN1;
   PHYSICAL(ZN2) = ZN2;
   PHYSICAL(MCha1) = MCha1;
   PHYSICAL(UM1) = UM1;
   PHYSICAL(UP1) = UP1;
   PHYSICAL(MCha2) = MCha2;
   PHYSICAL(UM2) = UM2;
   PHYSICAL(UP2) = UP2;
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
   if (PHYSICAL(MSRdp) < 0.) problems.flag_pole_tachyon(MRSSM_info::SRdp);
   if (PHYSICAL(MSRum) < 0.) problems.flag_pole_tachyon(MRSSM_info::SRum);
   if (PHYSICAL(MsigmaO) < 0.) problems.flag_pole_tachyon(MRSSM_info::sigmaO);
   if (PHYSICAL(MphiO) < 0.) problems.flag_pole_tachyon(MRSSM_info::phiO);
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Sd);
   if (PHYSICAL(MSv).tail<3>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Sv);
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Su);
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Se);
   if (PHYSICAL(Mhh).tail<4>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::hh);
   if (PHYSICAL(MAh).tail<3>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Ah);
   if (PHYSICAL(MRh).tail<2>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Rh);
   if (PHYSICAL(MHpm).tail<3>().minCoeff() < 0.) problems.flag_pole_tachyon(MRSSM_info::Hpm);

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
   MSRdp = 0.;
   MSRum = 0.;
   MsigmaO = 0.;
   MphiO = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,4,1>::Zero();
   ZH = Eigen::Matrix<double,4,4>::Zero();
   MAh = Eigen::Matrix<double,4,1>::Zero();
   ZA = Eigen::Matrix<double,4,4>::Zero();
   MRh = Eigen::Matrix<double,2,1>::Zero();
   ZHR = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,4,1>::Zero();
   ZP = Eigen::Matrix<double,4,4>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN1 = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   ZN2 = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha1 = Eigen::Matrix<double,2,1>::Zero();
   UM1 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP1 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MCha2 = Eigen::Matrix<double,2,1>::Zero();
   UM2 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP2 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
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



}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   MRSSM_soft_parameters::clear();
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
   MSRdp = pars(5);
   MSRum = pars(6);
   MsigmaO = pars(7);
   MphiO = pars(8);
   MSd(0) = pars(9);
   MSd(1) = pars(10);
   MSd(2) = pars(11);
   MSd(3) = pars(12);
   MSd(4) = pars(13);
   MSd(5) = pars(14);
   MSv(0) = pars(15);
   MSv(1) = pars(16);
   MSv(2) = pars(17);
   MSu(0) = pars(18);
   MSu(1) = pars(19);
   MSu(2) = pars(20);
   MSu(3) = pars(21);
   MSu(4) = pars(22);
   MSu(5) = pars(23);
   MSe(0) = pars(24);
   MSe(1) = pars(25);
   MSe(2) = pars(26);
   MSe(3) = pars(27);
   MSe(4) = pars(28);
   MSe(5) = pars(29);
   Mhh(0) = pars(30);
   Mhh(1) = pars(31);
   Mhh(2) = pars(32);
   Mhh(3) = pars(33);
   MAh(0) = pars(34);
   MAh(1) = pars(35);
   MAh(2) = pars(36);
   MAh(3) = pars(37);
   MRh(0) = pars(38);
   MRh(1) = pars(39);
   MHpm(0) = pars(40);
   MHpm(1) = pars(41);
   MHpm(2) = pars(42);
   MHpm(3) = pars(43);
   MChi(0) = pars(44);
   MChi(1) = pars(45);
   MChi(2) = pars(46);
   MChi(3) = pars(47);
   MCha1(0) = pars(48);
   MCha1(1) = pars(49);
   MCha2(0) = pars(50);
   MCha2(1) = pars(51);
   MFe(0) = pars(52);
   MFe(1) = pars(53);
   MFe(2) = pars(54);
   MFd(0) = pars(55);
   MFd(1) = pars(56);
   MFd(2) = pars(57);
   MFu(0) = pars(58);
   MFu(1) = pars(59);
   MFu(2) = pars(60);
   MVWm = pars(61);
   MVP = pars(62);
   MVZ = pars(63);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(64);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSRdp;
   pars(6) = MSRum;
   pars(7) = MsigmaO;
   pars(8) = MphiO;
   pars(9) = MSd(0);
   pars(10) = MSd(1);
   pars(11) = MSd(2);
   pars(12) = MSd(3);
   pars(13) = MSd(4);
   pars(14) = MSd(5);
   pars(15) = MSv(0);
   pars(16) = MSv(1);
   pars(17) = MSv(2);
   pars(18) = MSu(0);
   pars(19) = MSu(1);
   pars(20) = MSu(2);
   pars(21) = MSu(3);
   pars(22) = MSu(4);
   pars(23) = MSu(5);
   pars(24) = MSe(0);
   pars(25) = MSe(1);
   pars(26) = MSe(2);
   pars(27) = MSe(3);
   pars(28) = MSe(4);
   pars(29) = MSe(5);
   pars(30) = Mhh(0);
   pars(31) = Mhh(1);
   pars(32) = Mhh(2);
   pars(33) = Mhh(3);
   pars(34) = MAh(0);
   pars(35) = MAh(1);
   pars(36) = MAh(2);
   pars(37) = MAh(3);
   pars(38) = MRh(0);
   pars(39) = MRh(1);
   pars(40) = MHpm(0);
   pars(41) = MHpm(1);
   pars(42) = MHpm(2);
   pars(43) = MHpm(3);
   pars(44) = MChi(0);
   pars(45) = MChi(1);
   pars(46) = MChi(2);
   pars(47) = MChi(3);
   pars(48) = MCha1(0);
   pars(49) = MCha1(1);
   pars(50) = MCha2(0);
   pars(51) = MCha2(1);
   pars(52) = MFe(0);
   pars(53) = MFe(1);
   pars(54) = MFe(2);
   pars(55) = MFd(0);
   pars(56) = MFd(1);
   pars(57) = MFd(2);
   pars(58) = MFu(0);
   pars(59) = MFu(1);
   pars(60) = MFu(2);
   pars(61) = MVWm;
   pars(62) = MVP;
   pars(63) = MVZ;

   return pars;
}

void CLASSNAME::set_DRbar_masses_and_mixings(const Eigen::ArrayXd& pars)
{
   set_DRbar_masses(pars);

   ZD(0,0) = pars(64);
   ZD(0,1) = pars(65);
   ZD(0,2) = pars(66);
   ZD(0,3) = pars(67);
   ZD(0,4) = pars(68);
   ZD(0,5) = pars(69);
   ZD(1,0) = pars(70);
   ZD(1,1) = pars(71);
   ZD(1,2) = pars(72);
   ZD(1,3) = pars(73);
   ZD(1,4) = pars(74);
   ZD(1,5) = pars(75);
   ZD(2,0) = pars(76);
   ZD(2,1) = pars(77);
   ZD(2,2) = pars(78);
   ZD(2,3) = pars(79);
   ZD(2,4) = pars(80);
   ZD(2,5) = pars(81);
   ZD(3,0) = pars(82);
   ZD(3,1) = pars(83);
   ZD(3,2) = pars(84);
   ZD(3,3) = pars(85);
   ZD(3,4) = pars(86);
   ZD(3,5) = pars(87);
   ZD(4,0) = pars(88);
   ZD(4,1) = pars(89);
   ZD(4,2) = pars(90);
   ZD(4,3) = pars(91);
   ZD(4,4) = pars(92);
   ZD(4,5) = pars(93);
   ZD(5,0) = pars(94);
   ZD(5,1) = pars(95);
   ZD(5,2) = pars(96);
   ZD(5,3) = pars(97);
   ZD(5,4) = pars(98);
   ZD(5,5) = pars(99);
   ZV(0,0) = pars(100);
   ZV(0,1) = pars(101);
   ZV(0,2) = pars(102);
   ZV(1,0) = pars(103);
   ZV(1,1) = pars(104);
   ZV(1,2) = pars(105);
   ZV(2,0) = pars(106);
   ZV(2,1) = pars(107);
   ZV(2,2) = pars(108);
   ZU(0,0) = pars(109);
   ZU(0,1) = pars(110);
   ZU(0,2) = pars(111);
   ZU(0,3) = pars(112);
   ZU(0,4) = pars(113);
   ZU(0,5) = pars(114);
   ZU(1,0) = pars(115);
   ZU(1,1) = pars(116);
   ZU(1,2) = pars(117);
   ZU(1,3) = pars(118);
   ZU(1,4) = pars(119);
   ZU(1,5) = pars(120);
   ZU(2,0) = pars(121);
   ZU(2,1) = pars(122);
   ZU(2,2) = pars(123);
   ZU(2,3) = pars(124);
   ZU(2,4) = pars(125);
   ZU(2,5) = pars(126);
   ZU(3,0) = pars(127);
   ZU(3,1) = pars(128);
   ZU(3,2) = pars(129);
   ZU(3,3) = pars(130);
   ZU(3,4) = pars(131);
   ZU(3,5) = pars(132);
   ZU(4,0) = pars(133);
   ZU(4,1) = pars(134);
   ZU(4,2) = pars(135);
   ZU(4,3) = pars(136);
   ZU(4,4) = pars(137);
   ZU(4,5) = pars(138);
   ZU(5,0) = pars(139);
   ZU(5,1) = pars(140);
   ZU(5,2) = pars(141);
   ZU(5,3) = pars(142);
   ZU(5,4) = pars(143);
   ZU(5,5) = pars(144);
   ZE(0,0) = pars(145);
   ZE(0,1) = pars(146);
   ZE(0,2) = pars(147);
   ZE(0,3) = pars(148);
   ZE(0,4) = pars(149);
   ZE(0,5) = pars(150);
   ZE(1,0) = pars(151);
   ZE(1,1) = pars(152);
   ZE(1,2) = pars(153);
   ZE(1,3) = pars(154);
   ZE(1,4) = pars(155);
   ZE(1,5) = pars(156);
   ZE(2,0) = pars(157);
   ZE(2,1) = pars(158);
   ZE(2,2) = pars(159);
   ZE(2,3) = pars(160);
   ZE(2,4) = pars(161);
   ZE(2,5) = pars(162);
   ZE(3,0) = pars(163);
   ZE(3,1) = pars(164);
   ZE(3,2) = pars(165);
   ZE(3,3) = pars(166);
   ZE(3,4) = pars(167);
   ZE(3,5) = pars(168);
   ZE(4,0) = pars(169);
   ZE(4,1) = pars(170);
   ZE(4,2) = pars(171);
   ZE(4,3) = pars(172);
   ZE(4,4) = pars(173);
   ZE(4,5) = pars(174);
   ZE(5,0) = pars(175);
   ZE(5,1) = pars(176);
   ZE(5,2) = pars(177);
   ZE(5,3) = pars(178);
   ZE(5,4) = pars(179);
   ZE(5,5) = pars(180);
   ZH(0,0) = pars(181);
   ZH(0,1) = pars(182);
   ZH(0,2) = pars(183);
   ZH(0,3) = pars(184);
   ZH(1,0) = pars(185);
   ZH(1,1) = pars(186);
   ZH(1,2) = pars(187);
   ZH(1,3) = pars(188);
   ZH(2,0) = pars(189);
   ZH(2,1) = pars(190);
   ZH(2,2) = pars(191);
   ZH(2,3) = pars(192);
   ZH(3,0) = pars(193);
   ZH(3,1) = pars(194);
   ZH(3,2) = pars(195);
   ZH(3,3) = pars(196);
   ZA(0,0) = pars(197);
   ZA(0,1) = pars(198);
   ZA(0,2) = pars(199);
   ZA(0,3) = pars(200);
   ZA(1,0) = pars(201);
   ZA(1,1) = pars(202);
   ZA(1,2) = pars(203);
   ZA(1,3) = pars(204);
   ZA(2,0) = pars(205);
   ZA(2,1) = pars(206);
   ZA(2,2) = pars(207);
   ZA(2,3) = pars(208);
   ZA(3,0) = pars(209);
   ZA(3,1) = pars(210);
   ZA(3,2) = pars(211);
   ZA(3,3) = pars(212);
   ZHR(0,0) = pars(213);
   ZHR(0,1) = pars(214);
   ZHR(1,0) = pars(215);
   ZHR(1,1) = pars(216);
   ZP(0,0) = pars(217);
   ZP(0,1) = pars(218);
   ZP(0,2) = pars(219);
   ZP(0,3) = pars(220);
   ZP(1,0) = pars(221);
   ZP(1,1) = pars(222);
   ZP(1,2) = pars(223);
   ZP(1,3) = pars(224);
   ZP(2,0) = pars(225);
   ZP(2,1) = pars(226);
   ZP(2,2) = pars(227);
   ZP(2,3) = pars(228);
   ZP(3,0) = pars(229);
   ZP(3,1) = pars(230);
   ZP(3,2) = pars(231);
   ZP(3,3) = pars(232);
   ZN1(0,0) = std::complex<double>(pars(233), pars(234));
   ZN1(0,1) = std::complex<double>(pars(235), pars(236));
   ZN1(0,2) = std::complex<double>(pars(237), pars(238));
   ZN1(0,3) = std::complex<double>(pars(239), pars(240));
   ZN1(1,0) = std::complex<double>(pars(241), pars(242));
   ZN1(1,1) = std::complex<double>(pars(243), pars(244));
   ZN1(1,2) = std::complex<double>(pars(245), pars(246));
   ZN1(1,3) = std::complex<double>(pars(247), pars(248));
   ZN1(2,0) = std::complex<double>(pars(249), pars(250));
   ZN1(2,1) = std::complex<double>(pars(251), pars(252));
   ZN1(2,2) = std::complex<double>(pars(253), pars(254));
   ZN1(2,3) = std::complex<double>(pars(255), pars(256));
   ZN1(3,0) = std::complex<double>(pars(257), pars(258));
   ZN1(3,1) = std::complex<double>(pars(259), pars(260));
   ZN1(3,2) = std::complex<double>(pars(261), pars(262));
   ZN1(3,3) = std::complex<double>(pars(263), pars(264));
   ZN2(0,0) = std::complex<double>(pars(265), pars(266));
   ZN2(0,1) = std::complex<double>(pars(267), pars(268));
   ZN2(0,2) = std::complex<double>(pars(269), pars(270));
   ZN2(0,3) = std::complex<double>(pars(271), pars(272));
   ZN2(1,0) = std::complex<double>(pars(273), pars(274));
   ZN2(1,1) = std::complex<double>(pars(275), pars(276));
   ZN2(1,2) = std::complex<double>(pars(277), pars(278));
   ZN2(1,3) = std::complex<double>(pars(279), pars(280));
   ZN2(2,0) = std::complex<double>(pars(281), pars(282));
   ZN2(2,1) = std::complex<double>(pars(283), pars(284));
   ZN2(2,2) = std::complex<double>(pars(285), pars(286));
   ZN2(2,3) = std::complex<double>(pars(287), pars(288));
   ZN2(3,0) = std::complex<double>(pars(289), pars(290));
   ZN2(3,1) = std::complex<double>(pars(291), pars(292));
   ZN2(3,2) = std::complex<double>(pars(293), pars(294));
   ZN2(3,3) = std::complex<double>(pars(295), pars(296));
   UM1(0,0) = std::complex<double>(pars(297), pars(298));
   UM1(0,1) = std::complex<double>(pars(299), pars(300));
   UM1(1,0) = std::complex<double>(pars(301), pars(302));
   UM1(1,1) = std::complex<double>(pars(303), pars(304));
   UP1(0,0) = std::complex<double>(pars(305), pars(306));
   UP1(0,1) = std::complex<double>(pars(307), pars(308));
   UP1(1,0) = std::complex<double>(pars(309), pars(310));
   UP1(1,1) = std::complex<double>(pars(311), pars(312));
   UM2(0,0) = std::complex<double>(pars(313), pars(314));
   UM2(0,1) = std::complex<double>(pars(315), pars(316));
   UM2(1,0) = std::complex<double>(pars(317), pars(318));
   UM2(1,1) = std::complex<double>(pars(319), pars(320));
   UP2(0,0) = std::complex<double>(pars(321), pars(322));
   UP2(0,1) = std::complex<double>(pars(323), pars(324));
   UP2(1,0) = std::complex<double>(pars(325), pars(326));
   UP2(1,1) = std::complex<double>(pars(327), pars(328));
   ZEL(0,0) = std::complex<double>(pars(329), pars(330));
   ZEL(0,1) = std::complex<double>(pars(331), pars(332));
   ZEL(0,2) = std::complex<double>(pars(333), pars(334));
   ZEL(1,0) = std::complex<double>(pars(335), pars(336));
   ZEL(1,1) = std::complex<double>(pars(337), pars(338));
   ZEL(1,2) = std::complex<double>(pars(339), pars(340));
   ZEL(2,0) = std::complex<double>(pars(341), pars(342));
   ZEL(2,1) = std::complex<double>(pars(343), pars(344));
   ZEL(2,2) = std::complex<double>(pars(345), pars(346));
   ZER(0,0) = std::complex<double>(pars(347), pars(348));
   ZER(0,1) = std::complex<double>(pars(349), pars(350));
   ZER(0,2) = std::complex<double>(pars(351), pars(352));
   ZER(1,0) = std::complex<double>(pars(353), pars(354));
   ZER(1,1) = std::complex<double>(pars(355), pars(356));
   ZER(1,2) = std::complex<double>(pars(357), pars(358));
   ZER(2,0) = std::complex<double>(pars(359), pars(360));
   ZER(2,1) = std::complex<double>(pars(361), pars(362));
   ZER(2,2) = std::complex<double>(pars(363), pars(364));
   ZDL(0,0) = std::complex<double>(pars(365), pars(366));
   ZDL(0,1) = std::complex<double>(pars(367), pars(368));
   ZDL(0,2) = std::complex<double>(pars(369), pars(370));
   ZDL(1,0) = std::complex<double>(pars(371), pars(372));
   ZDL(1,1) = std::complex<double>(pars(373), pars(374));
   ZDL(1,2) = std::complex<double>(pars(375), pars(376));
   ZDL(2,0) = std::complex<double>(pars(377), pars(378));
   ZDL(2,1) = std::complex<double>(pars(379), pars(380));
   ZDL(2,2) = std::complex<double>(pars(381), pars(382));
   ZDR(0,0) = std::complex<double>(pars(383), pars(384));
   ZDR(0,1) = std::complex<double>(pars(385), pars(386));
   ZDR(0,2) = std::complex<double>(pars(387), pars(388));
   ZDR(1,0) = std::complex<double>(pars(389), pars(390));
   ZDR(1,1) = std::complex<double>(pars(391), pars(392));
   ZDR(1,2) = std::complex<double>(pars(393), pars(394));
   ZDR(2,0) = std::complex<double>(pars(395), pars(396));
   ZDR(2,1) = std::complex<double>(pars(397), pars(398));
   ZDR(2,2) = std::complex<double>(pars(399), pars(400));
   ZUL(0,0) = std::complex<double>(pars(401), pars(402));
   ZUL(0,1) = std::complex<double>(pars(403), pars(404));
   ZUL(0,2) = std::complex<double>(pars(405), pars(406));
   ZUL(1,0) = std::complex<double>(pars(407), pars(408));
   ZUL(1,1) = std::complex<double>(pars(409), pars(410));
   ZUL(1,2) = std::complex<double>(pars(411), pars(412));
   ZUL(2,0) = std::complex<double>(pars(413), pars(414));
   ZUL(2,1) = std::complex<double>(pars(415), pars(416));
   ZUL(2,2) = std::complex<double>(pars(417), pars(418));
   ZUR(0,0) = std::complex<double>(pars(419), pars(420));
   ZUR(0,1) = std::complex<double>(pars(421), pars(422));
   ZUR(0,2) = std::complex<double>(pars(423), pars(424));
   ZUR(1,0) = std::complex<double>(pars(425), pars(426));
   ZUR(1,1) = std::complex<double>(pars(427), pars(428));
   ZUR(1,2) = std::complex<double>(pars(429), pars(430));
   ZUR(2,0) = std::complex<double>(pars(431), pars(432));
   ZUR(2,1) = std::complex<double>(pars(433), pars(434));
   ZUR(2,2) = std::complex<double>(pars(435), pars(436));
   ZZ(0,0) = pars(437);
   ZZ(0,1) = pars(438);
   ZZ(1,0) = pars(439);
   ZZ(1,1) = pars(440);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses_and_mixings() const
{
   Eigen::ArrayXd pars(get_DRbar_masses());

   pars.conservativeResize(441);

   pars(64) = ZD(0,0);
   pars(65) = ZD(0,1);
   pars(66) = ZD(0,2);
   pars(67) = ZD(0,3);
   pars(68) = ZD(0,4);
   pars(69) = ZD(0,5);
   pars(70) = ZD(1,0);
   pars(71) = ZD(1,1);
   pars(72) = ZD(1,2);
   pars(73) = ZD(1,3);
   pars(74) = ZD(1,4);
   pars(75) = ZD(1,5);
   pars(76) = ZD(2,0);
   pars(77) = ZD(2,1);
   pars(78) = ZD(2,2);
   pars(79) = ZD(2,3);
   pars(80) = ZD(2,4);
   pars(81) = ZD(2,5);
   pars(82) = ZD(3,0);
   pars(83) = ZD(3,1);
   pars(84) = ZD(3,2);
   pars(85) = ZD(3,3);
   pars(86) = ZD(3,4);
   pars(87) = ZD(3,5);
   pars(88) = ZD(4,0);
   pars(89) = ZD(4,1);
   pars(90) = ZD(4,2);
   pars(91) = ZD(4,3);
   pars(92) = ZD(4,4);
   pars(93) = ZD(4,5);
   pars(94) = ZD(5,0);
   pars(95) = ZD(5,1);
   pars(96) = ZD(5,2);
   pars(97) = ZD(5,3);
   pars(98) = ZD(5,4);
   pars(99) = ZD(5,5);
   pars(100) = ZV(0,0);
   pars(101) = ZV(0,1);
   pars(102) = ZV(0,2);
   pars(103) = ZV(1,0);
   pars(104) = ZV(1,1);
   pars(105) = ZV(1,2);
   pars(106) = ZV(2,0);
   pars(107) = ZV(2,1);
   pars(108) = ZV(2,2);
   pars(109) = ZU(0,0);
   pars(110) = ZU(0,1);
   pars(111) = ZU(0,2);
   pars(112) = ZU(0,3);
   pars(113) = ZU(0,4);
   pars(114) = ZU(0,5);
   pars(115) = ZU(1,0);
   pars(116) = ZU(1,1);
   pars(117) = ZU(1,2);
   pars(118) = ZU(1,3);
   pars(119) = ZU(1,4);
   pars(120) = ZU(1,5);
   pars(121) = ZU(2,0);
   pars(122) = ZU(2,1);
   pars(123) = ZU(2,2);
   pars(124) = ZU(2,3);
   pars(125) = ZU(2,4);
   pars(126) = ZU(2,5);
   pars(127) = ZU(3,0);
   pars(128) = ZU(3,1);
   pars(129) = ZU(3,2);
   pars(130) = ZU(3,3);
   pars(131) = ZU(3,4);
   pars(132) = ZU(3,5);
   pars(133) = ZU(4,0);
   pars(134) = ZU(4,1);
   pars(135) = ZU(4,2);
   pars(136) = ZU(4,3);
   pars(137) = ZU(4,4);
   pars(138) = ZU(4,5);
   pars(139) = ZU(5,0);
   pars(140) = ZU(5,1);
   pars(141) = ZU(5,2);
   pars(142) = ZU(5,3);
   pars(143) = ZU(5,4);
   pars(144) = ZU(5,5);
   pars(145) = ZE(0,0);
   pars(146) = ZE(0,1);
   pars(147) = ZE(0,2);
   pars(148) = ZE(0,3);
   pars(149) = ZE(0,4);
   pars(150) = ZE(0,5);
   pars(151) = ZE(1,0);
   pars(152) = ZE(1,1);
   pars(153) = ZE(1,2);
   pars(154) = ZE(1,3);
   pars(155) = ZE(1,4);
   pars(156) = ZE(1,5);
   pars(157) = ZE(2,0);
   pars(158) = ZE(2,1);
   pars(159) = ZE(2,2);
   pars(160) = ZE(2,3);
   pars(161) = ZE(2,4);
   pars(162) = ZE(2,5);
   pars(163) = ZE(3,0);
   pars(164) = ZE(3,1);
   pars(165) = ZE(3,2);
   pars(166) = ZE(3,3);
   pars(167) = ZE(3,4);
   pars(168) = ZE(3,5);
   pars(169) = ZE(4,0);
   pars(170) = ZE(4,1);
   pars(171) = ZE(4,2);
   pars(172) = ZE(4,3);
   pars(173) = ZE(4,4);
   pars(174) = ZE(4,5);
   pars(175) = ZE(5,0);
   pars(176) = ZE(5,1);
   pars(177) = ZE(5,2);
   pars(178) = ZE(5,3);
   pars(179) = ZE(5,4);
   pars(180) = ZE(5,5);
   pars(181) = ZH(0,0);
   pars(182) = ZH(0,1);
   pars(183) = ZH(0,2);
   pars(184) = ZH(0,3);
   pars(185) = ZH(1,0);
   pars(186) = ZH(1,1);
   pars(187) = ZH(1,2);
   pars(188) = ZH(1,3);
   pars(189) = ZH(2,0);
   pars(190) = ZH(2,1);
   pars(191) = ZH(2,2);
   pars(192) = ZH(2,3);
   pars(193) = ZH(3,0);
   pars(194) = ZH(3,1);
   pars(195) = ZH(3,2);
   pars(196) = ZH(3,3);
   pars(197) = ZA(0,0);
   pars(198) = ZA(0,1);
   pars(199) = ZA(0,2);
   pars(200) = ZA(0,3);
   pars(201) = ZA(1,0);
   pars(202) = ZA(1,1);
   pars(203) = ZA(1,2);
   pars(204) = ZA(1,3);
   pars(205) = ZA(2,0);
   pars(206) = ZA(2,1);
   pars(207) = ZA(2,2);
   pars(208) = ZA(2,3);
   pars(209) = ZA(3,0);
   pars(210) = ZA(3,1);
   pars(211) = ZA(3,2);
   pars(212) = ZA(3,3);
   pars(213) = ZHR(0,0);
   pars(214) = ZHR(0,1);
   pars(215) = ZHR(1,0);
   pars(216) = ZHR(1,1);
   pars(217) = ZP(0,0);
   pars(218) = ZP(0,1);
   pars(219) = ZP(0,2);
   pars(220) = ZP(0,3);
   pars(221) = ZP(1,0);
   pars(222) = ZP(1,1);
   pars(223) = ZP(1,2);
   pars(224) = ZP(1,3);
   pars(225) = ZP(2,0);
   pars(226) = ZP(2,1);
   pars(227) = ZP(2,2);
   pars(228) = ZP(2,3);
   pars(229) = ZP(3,0);
   pars(230) = ZP(3,1);
   pars(231) = ZP(3,2);
   pars(232) = ZP(3,3);
   pars(233) = Re(ZN1(0,0));
   pars(234) = Im(ZN1(0,0));
   pars(235) = Re(ZN1(0,1));
   pars(236) = Im(ZN1(0,1));
   pars(237) = Re(ZN1(0,2));
   pars(238) = Im(ZN1(0,2));
   pars(239) = Re(ZN1(0,3));
   pars(240) = Im(ZN1(0,3));
   pars(241) = Re(ZN1(1,0));
   pars(242) = Im(ZN1(1,0));
   pars(243) = Re(ZN1(1,1));
   pars(244) = Im(ZN1(1,1));
   pars(245) = Re(ZN1(1,2));
   pars(246) = Im(ZN1(1,2));
   pars(247) = Re(ZN1(1,3));
   pars(248) = Im(ZN1(1,3));
   pars(249) = Re(ZN1(2,0));
   pars(250) = Im(ZN1(2,0));
   pars(251) = Re(ZN1(2,1));
   pars(252) = Im(ZN1(2,1));
   pars(253) = Re(ZN1(2,2));
   pars(254) = Im(ZN1(2,2));
   pars(255) = Re(ZN1(2,3));
   pars(256) = Im(ZN1(2,3));
   pars(257) = Re(ZN1(3,0));
   pars(258) = Im(ZN1(3,0));
   pars(259) = Re(ZN1(3,1));
   pars(260) = Im(ZN1(3,1));
   pars(261) = Re(ZN1(3,2));
   pars(262) = Im(ZN1(3,2));
   pars(263) = Re(ZN1(3,3));
   pars(264) = Im(ZN1(3,3));
   pars(265) = Re(ZN2(0,0));
   pars(266) = Im(ZN2(0,0));
   pars(267) = Re(ZN2(0,1));
   pars(268) = Im(ZN2(0,1));
   pars(269) = Re(ZN2(0,2));
   pars(270) = Im(ZN2(0,2));
   pars(271) = Re(ZN2(0,3));
   pars(272) = Im(ZN2(0,3));
   pars(273) = Re(ZN2(1,0));
   pars(274) = Im(ZN2(1,0));
   pars(275) = Re(ZN2(1,1));
   pars(276) = Im(ZN2(1,1));
   pars(277) = Re(ZN2(1,2));
   pars(278) = Im(ZN2(1,2));
   pars(279) = Re(ZN2(1,3));
   pars(280) = Im(ZN2(1,3));
   pars(281) = Re(ZN2(2,0));
   pars(282) = Im(ZN2(2,0));
   pars(283) = Re(ZN2(2,1));
   pars(284) = Im(ZN2(2,1));
   pars(285) = Re(ZN2(2,2));
   pars(286) = Im(ZN2(2,2));
   pars(287) = Re(ZN2(2,3));
   pars(288) = Im(ZN2(2,3));
   pars(289) = Re(ZN2(3,0));
   pars(290) = Im(ZN2(3,0));
   pars(291) = Re(ZN2(3,1));
   pars(292) = Im(ZN2(3,1));
   pars(293) = Re(ZN2(3,2));
   pars(294) = Im(ZN2(3,2));
   pars(295) = Re(ZN2(3,3));
   pars(296) = Im(ZN2(3,3));
   pars(297) = Re(UM1(0,0));
   pars(298) = Im(UM1(0,0));
   pars(299) = Re(UM1(0,1));
   pars(300) = Im(UM1(0,1));
   pars(301) = Re(UM1(1,0));
   pars(302) = Im(UM1(1,0));
   pars(303) = Re(UM1(1,1));
   pars(304) = Im(UM1(1,1));
   pars(305) = Re(UP1(0,0));
   pars(306) = Im(UP1(0,0));
   pars(307) = Re(UP1(0,1));
   pars(308) = Im(UP1(0,1));
   pars(309) = Re(UP1(1,0));
   pars(310) = Im(UP1(1,0));
   pars(311) = Re(UP1(1,1));
   pars(312) = Im(UP1(1,1));
   pars(313) = Re(UM2(0,0));
   pars(314) = Im(UM2(0,0));
   pars(315) = Re(UM2(0,1));
   pars(316) = Im(UM2(0,1));
   pars(317) = Re(UM2(1,0));
   pars(318) = Im(UM2(1,0));
   pars(319) = Re(UM2(1,1));
   pars(320) = Im(UM2(1,1));
   pars(321) = Re(UP2(0,0));
   pars(322) = Im(UP2(0,0));
   pars(323) = Re(UP2(0,1));
   pars(324) = Im(UP2(0,1));
   pars(325) = Re(UP2(1,0));
   pars(326) = Im(UP2(1,0));
   pars(327) = Re(UP2(1,1));
   pars(328) = Im(UP2(1,1));
   pars(329) = Re(ZEL(0,0));
   pars(330) = Im(ZEL(0,0));
   pars(331) = Re(ZEL(0,1));
   pars(332) = Im(ZEL(0,1));
   pars(333) = Re(ZEL(0,2));
   pars(334) = Im(ZEL(0,2));
   pars(335) = Re(ZEL(1,0));
   pars(336) = Im(ZEL(1,0));
   pars(337) = Re(ZEL(1,1));
   pars(338) = Im(ZEL(1,1));
   pars(339) = Re(ZEL(1,2));
   pars(340) = Im(ZEL(1,2));
   pars(341) = Re(ZEL(2,0));
   pars(342) = Im(ZEL(2,0));
   pars(343) = Re(ZEL(2,1));
   pars(344) = Im(ZEL(2,1));
   pars(345) = Re(ZEL(2,2));
   pars(346) = Im(ZEL(2,2));
   pars(347) = Re(ZER(0,0));
   pars(348) = Im(ZER(0,0));
   pars(349) = Re(ZER(0,1));
   pars(350) = Im(ZER(0,1));
   pars(351) = Re(ZER(0,2));
   pars(352) = Im(ZER(0,2));
   pars(353) = Re(ZER(1,0));
   pars(354) = Im(ZER(1,0));
   pars(355) = Re(ZER(1,1));
   pars(356) = Im(ZER(1,1));
   pars(357) = Re(ZER(1,2));
   pars(358) = Im(ZER(1,2));
   pars(359) = Re(ZER(2,0));
   pars(360) = Im(ZER(2,0));
   pars(361) = Re(ZER(2,1));
   pars(362) = Im(ZER(2,1));
   pars(363) = Re(ZER(2,2));
   pars(364) = Im(ZER(2,2));
   pars(365) = Re(ZDL(0,0));
   pars(366) = Im(ZDL(0,0));
   pars(367) = Re(ZDL(0,1));
   pars(368) = Im(ZDL(0,1));
   pars(369) = Re(ZDL(0,2));
   pars(370) = Im(ZDL(0,2));
   pars(371) = Re(ZDL(1,0));
   pars(372) = Im(ZDL(1,0));
   pars(373) = Re(ZDL(1,1));
   pars(374) = Im(ZDL(1,1));
   pars(375) = Re(ZDL(1,2));
   pars(376) = Im(ZDL(1,2));
   pars(377) = Re(ZDL(2,0));
   pars(378) = Im(ZDL(2,0));
   pars(379) = Re(ZDL(2,1));
   pars(380) = Im(ZDL(2,1));
   pars(381) = Re(ZDL(2,2));
   pars(382) = Im(ZDL(2,2));
   pars(383) = Re(ZDR(0,0));
   pars(384) = Im(ZDR(0,0));
   pars(385) = Re(ZDR(0,1));
   pars(386) = Im(ZDR(0,1));
   pars(387) = Re(ZDR(0,2));
   pars(388) = Im(ZDR(0,2));
   pars(389) = Re(ZDR(1,0));
   pars(390) = Im(ZDR(1,0));
   pars(391) = Re(ZDR(1,1));
   pars(392) = Im(ZDR(1,1));
   pars(393) = Re(ZDR(1,2));
   pars(394) = Im(ZDR(1,2));
   pars(395) = Re(ZDR(2,0));
   pars(396) = Im(ZDR(2,0));
   pars(397) = Re(ZDR(2,1));
   pars(398) = Im(ZDR(2,1));
   pars(399) = Re(ZDR(2,2));
   pars(400) = Im(ZDR(2,2));
   pars(401) = Re(ZUL(0,0));
   pars(402) = Im(ZUL(0,0));
   pars(403) = Re(ZUL(0,1));
   pars(404) = Im(ZUL(0,1));
   pars(405) = Re(ZUL(0,2));
   pars(406) = Im(ZUL(0,2));
   pars(407) = Re(ZUL(1,0));
   pars(408) = Im(ZUL(1,0));
   pars(409) = Re(ZUL(1,1));
   pars(410) = Im(ZUL(1,1));
   pars(411) = Re(ZUL(1,2));
   pars(412) = Im(ZUL(1,2));
   pars(413) = Re(ZUL(2,0));
   pars(414) = Im(ZUL(2,0));
   pars(415) = Re(ZUL(2,1));
   pars(416) = Im(ZUL(2,1));
   pars(417) = Re(ZUL(2,2));
   pars(418) = Im(ZUL(2,2));
   pars(419) = Re(ZUR(0,0));
   pars(420) = Im(ZUR(0,0));
   pars(421) = Re(ZUR(0,1));
   pars(422) = Im(ZUR(0,1));
   pars(423) = Re(ZUR(0,2));
   pars(424) = Im(ZUR(0,2));
   pars(425) = Re(ZUR(1,0));
   pars(426) = Im(ZUR(1,0));
   pars(427) = Re(ZUR(1,1));
   pars(428) = Im(ZUR(1,1));
   pars(429) = Re(ZUR(1,2));
   pars(430) = Im(ZUR(1,2));
   pars(431) = Re(ZUR(2,0));
   pars(432) = Im(ZUR(2,0));
   pars(433) = Re(ZUR(2,1));
   pars(434) = Im(ZUR(2,1));
   pars(435) = Re(ZUR(2,2));
   pars(436) = Im(ZUR(2,2));
   pars(437) = ZZ(0,0);
   pars(438) = ZZ(0,1);
   pars(439) = ZZ(1,0);
   pars(440) = ZZ(1,1);


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
   return "MRSSM";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   MRSSM_soft_parameters::run_to(scale, eps);
}


Eigen::Array<double,3,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_goldstone;
   MHpm_goldstone(0) = MVWm;

   return remove_if_equal(MHpm, MHpm_goldstone);
}

Eigen::Array<double,3,1> CLASSNAME::get_MPseudoscalarHiggs() const
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
 * {Chi, Su, Sd, Se, Cha1, Cha2}
 *
 * @param particle_type particle type
 * @return mass of LSP
 */
double CLASSNAME::get_lsp(MRSSM_info::Particles& particle_type) const
{
   double lsp_mass = std::numeric_limits<double>::max();
   double tmp_mass;
   particle_type = MRSSM_info::NUMBER_OF_PARTICLES;

   tmp_mass = Abs(PHYSICAL(MChi(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = MRSSM_info::Chi;
   }

   tmp_mass = Abs(PHYSICAL(MSu(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = MRSSM_info::Su;
   }

   tmp_mass = Abs(PHYSICAL(MSd(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = MRSSM_info::Sd;
   }

   tmp_mass = Abs(PHYSICAL(MSe(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = MRSSM_info::Se;
   }

   tmp_mass = Abs(PHYSICAL(MCha1(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = MRSSM_info::Cha1;
   }

   tmp_mass = Abs(PHYSICAL(MCha2(0)));
   if (tmp_mass < lsp_mass) {
      lsp_mass = tmp_mass;
      particle_type = MRSSM_info::Cha2;
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
   const double mass_matrix_Glu = Re(MDGoc);

   return mass_matrix_Glu;
}

void CLASSNAME::calculate_MGlu()
{
   const auto mass_matrix_Glu = get_mass_matrix_Glu();
   MGlu = calculate_singlet_mass(mass_matrix_Glu);
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

double CLASSNAME::get_mass_matrix_SRdp() const
{
   const double mass_matrix_SRdp = Re(0.125*(8*mRd2 + 3.0983866769659336*
      g1*MDBS*vS + 4*g2*MDWBT*vT + 8*AbsSqr(MuD) - 2*vS*(-2.8284271247461903*
      MuD - 2*LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamSD) +
      3.0983866769659336*g1*vS*Conj(MDBS) + 4*g2*vT*Conj(MDWBT) +
      5.656854249492381*LamSD*vS*Conj(MuD) - 4*LamTD*vT*Conj(MuD) - 0.6*Sqr(g1)
      *Sqr(vd) + Sqr(g2)*Sqr(vd) - Conj(LamTD)*(4*MuD*vT + 2.8284271247461903*
      LamSD*vS*vT - 4*LamTD*Sqr(vd) - 2*LamTD*Sqr(vT)) + 0.6*Sqr(g1)*Sqr(vu) -
      Sqr(g2)*Sqr(vu)));

   return mass_matrix_SRdp;
}

void CLASSNAME::calculate_MSRdp()
{
   const auto mass_matrix_SRdp = get_mass_matrix_SRdp();
   MSRdp = mass_matrix_SRdp;

   if (MSRdp < 0.) {
      problems.flag_running_tachyon(MRSSM_info::SRdp);
   }

   MSRdp = AbsSqrt(MSRdp);
}

double CLASSNAME::get_mass_matrix_SRum() const
{
   const double mass_matrix_SRum = Re(0.125*(8*mRu2 - 3.0983866769659336*
      g1*MDBS*vS - 4*g2*MDWBT*vT + 8*AbsSqr(MuU) + 2*vS*(2.8284271247461903*MuU
      + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamSU) -
      3.0983866769659336*g1*vS*Conj(MDBS) - 4*g2*vT*Conj(MDWBT) +
      5.656854249492381*LamSU*vS*Conj(MuU) + 4*LamTU*vT*Conj(MuU) + 0.6*Sqr(g1)
      *Sqr(vd) - Sqr(g2)*Sqr(vd) - 0.6*Sqr(g1)*Sqr(vu) + Sqr(g2)*Sqr(vu) + 2*
      Conj(LamTU)*(2*MuU*vT + 1.4142135623730951*LamSU*vS*vT + LamTU*Sqr(vT) +
      2*LamTU*Sqr(vu))));

   return mass_matrix_SRum;
}

void CLASSNAME::calculate_MSRum()
{
   const auto mass_matrix_SRum = get_mass_matrix_SRum();
   MSRum = mass_matrix_SRum;

   if (MSRum < 0.) {
      problems.flag_running_tachyon(MRSSM_info::SRum);
   }

   MSRum = AbsSqrt(MSRum);
}

double CLASSNAME::get_mass_matrix_sigmaO() const
{
   const double mass_matrix_sigmaO = Re(moc2);

   return mass_matrix_sigmaO;
}

void CLASSNAME::calculate_MsigmaO()
{
   const auto mass_matrix_sigmaO = get_mass_matrix_sigmaO();
   MsigmaO = mass_matrix_sigmaO;

   if (MsigmaO < 0.) {
      problems.flag_running_tachyon(MRSSM_info::sigmaO);
   }

   MsigmaO = AbsSqrt(MsigmaO);
}

double CLASSNAME::get_mass_matrix_phiO() const
{
   const double mass_matrix_phiO = Re(moc2 + 4*Sqr(MDGoc));

   return mass_matrix_phiO;
}

void CLASSNAME::calculate_MphiO()
{
   const auto mass_matrix_phiO = get_mass_matrix_phiO();
   MphiO = mass_matrix_phiO;

   if (MphiO < 0.) {
      problems.flag_running_tachyon(MRSSM_info::phiO);
   }

   MphiO = AbsSqrt(MphiO);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(0,0
      ) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)) + AbsSqr(Yd(2,0)))*Sqr(vd) -
      0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) +
      Conj(Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) +
      Conj(Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = -0.7071067811865475*vu*Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(0,4) = -0.7071067811865475*vu*Conj(Yd(1,0))*Mu;
   mass_matrix_Sd(0,5) = -0.7071067811865475*vu*Conj(Yd(2,0))*Mu;
   mass_matrix_Sd(1,1) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(1,1
      ) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)) + AbsSqr(Yd(2,1)))*Sqr(vd) -
      0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) +
      Conj(Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = -0.7071067811865475*vu*Conj(Yd(0,1))*Mu;
   mass_matrix_Sd(1,4) = -0.7071067811865475*vu*Conj(Yd(1,1))*Mu;
   mass_matrix_Sd(1,5) = -0.7071067811865475*vu*Conj(Yd(2,1))*Mu;
   mass_matrix_Sd(2,2) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(2,2
      ) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) + AbsSqr(Yd(2,2)))*Sqr(vd) -
      0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = -0.7071067811865475*vu*Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = -0.7071067811865475*vu*Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = -0.7071067811865475*vu*Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,3) = 0.2581988897471611*g1*MDBS*vS +
      0.2581988897471611*g1*vS*Conj(MDBS) + md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) +
      AbsSqr(Yd(0,1)) + AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) +
      Conj(Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) +
      Conj(Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = 0.2581988897471611*g1*MDBS*vS +
      0.2581988897471611*g1*vS*Conj(MDBS) + md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) +
      AbsSqr(Yd(1,1)) + AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) +
      Conj(Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = 0.2581988897471611*g1*MDBS*vS +
      0.2581988897471611*g1*vS*Conj(MDBS) + md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) +
      AbsSqr(Yd(2,1)) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*
      Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Sd, eigenvalue_error > precision *
      Abs(MSd(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif
   normalize_to_interval(ZD);


   if (MSd.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Sd);
   }

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(0,0)
      + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu)
      - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(1,1)
      + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu)
      - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(2,2)
      + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu)
      - 0.125*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Sv, eigenvalue_error > precision *
      Abs(MSv(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif
   normalize_to_interval(ZV);


   if (MSv.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Sv);
   }

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(0,0
      ) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,0))
      + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) +
      Conj(Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) +
      Conj(Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = -0.7071067811865475*vd*Conj(Yu(0,0))*Mu;
   mass_matrix_Su(0,4) = -0.7071067811865475*vd*Conj(Yu(1,0))*Mu;
   mass_matrix_Su(0,5) = -0.7071067811865475*vd*Conj(Yu(2,0))*Mu;
   mass_matrix_Su(1,1) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(1,1
      ) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,1))
      + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) +
      Conj(Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = -0.7071067811865475*vd*Conj(Yu(0,1))*Mu;
   mass_matrix_Su(1,4) = -0.7071067811865475*vd*Conj(Yu(1,1))*Mu;
   mass_matrix_Su(1,5) = -0.7071067811865475*vd*Conj(Yu(2,1))*Mu;
   mass_matrix_Su(2,2) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(2,2
      ) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,2))
      + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = -0.7071067811865475*vd*Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = -0.7071067811865475*vd*Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = -0.7071067811865475*vd*Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,3) = -0.5163977794943222*g1*MDBS*vS -
      0.5163977794943222*g1*vS*Conj(MDBS) + mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) +
      0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) +
      Conj(Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) +
      Conj(Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = -0.5163977794943222*g1*MDBS*vS -
      0.5163977794943222*g1*vS*Conj(MDBS) + mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) +
      0.5*(AbsSqr(Yu(1,0)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) +
      Conj(Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = -0.5163977794943222*g1*MDBS*vS -
      0.5163977794943222*g1*vS*Conj(MDBS) + mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) +
      0.5*(AbsSqr(Yu(2,0)) + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*
      Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Su, eigenvalue_error > precision *
      Abs(MSu(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif
   normalize_to_interval(ZU);


   if (MSu.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Su);
   }

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(0,0)
      + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)) + AbsSqr(Ye(2,0)))*Sqr(vd) +
      0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) +
      Conj(Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) +
      Conj(Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = -0.7071067811865475*vu*Conj(Ye(0,0))*Mu;
   mass_matrix_Se(0,4) = -0.7071067811865475*vu*Conj(Ye(1,0))*Mu;
   mass_matrix_Se(0,5) = -0.7071067811865475*vu*Conj(Ye(2,0))*Mu;
   mass_matrix_Se(1,1) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(1,1)
      + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)) + AbsSqr(Ye(2,1)))*Sqr(vd) +
      0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) +
      Conj(Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = -0.7071067811865475*vu*Conj(Ye(0,1))*Mu;
   mass_matrix_Se(1,4) = -0.7071067811865475*vu*Conj(Ye(1,1))*Mu;
   mass_matrix_Se(1,5) = -0.7071067811865475*vu*Conj(Ye(2,1))*Mu;
   mass_matrix_Se(2,2) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(2,2)
      + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) + AbsSqr(Ye(2,2)))*Sqr(vd) +
      0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = -0.7071067811865475*vu*Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = -0.7071067811865475*vu*Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = -0.7071067811865475*vu*Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,3) = 0.7745966692414834*g1*MDBS*vS +
      0.7745966692414834*g1*vS*Conj(MDBS) + me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) +
      AbsSqr(Ye(0,1)) + AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) +
      Conj(Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) +
      Conj(Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = 0.7745966692414834*g1*MDBS*vS +
      0.7745966692414834*g1*vS*Conj(MDBS) + me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) +
      AbsSqr(Ye(1,1)) + AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) +
      Conj(Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = 0.7745966692414834*g1*MDBS*vS +
      0.7745966692414834*g1*vS*Conj(MDBS) + me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) +
      AbsSqr(Ye(2,1)) + AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*
      Sqr(g1)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Se, eigenvalue_error > precision *
      Abs(MSe(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif
   normalize_to_interval(ZE);


   if (MSe.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Se);
   }

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,4,4> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*
      MDWBT*vT + AbsSqr(MuD) + AbsSqr(Mu) + 0.7071067811865475*MuD*vS*Conj(
      LamSD) + 0.35355339059327373*LamTD*vS*vT*Conj(LamSD) + 0.5*MuD*vT*Conj(
      LamTD) + 0.35355339059327373*LamSD*vS*vT*Conj(LamTD) - 0.3872983346207417
      *g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*
      Conj(MuD) + 0.5*LamTD*vT*Conj(MuD) + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2
      )*Sqr(vd) + 0.5*AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) -
      0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = -0.5*BMu - 0.5*Conj(BMu) - 0.15*vd*vu*Sqr(g1) -
      0.25*vd*vu*Sqr(g2);
   mass_matrix_hh(0,2) = -0.3872983346207417*g1*MDBS*vd + vd*vS*AbsSqr(
      LamSD) + 0.7071067811865475*MuD*vd*Conj(LamSD) + 0.35355339059327373*
      LamTD*vd*vT*Conj(LamSD) + 0.35355339059327373*LamSD*vd*vT*Conj(LamTD) -
      0.3872983346207417*g1*vd*Conj(MDBS) + 0.7071067811865475*LamSD*vd*Conj(
      MuD);
   mass_matrix_hh(0,3) = 0.5*g2*MDWBT*vd + 0.5*vd*vT*AbsSqr(LamTD) +
      0.35355339059327373*LamTD*vd*vS*Conj(LamSD) + 0.5*MuD*vd*Conj(LamTD) +
      0.35355339059327373*LamSD*vd*vS*Conj(LamTD) + 0.5*g2*vd*Conj(MDWBT) + 0.5
      *LamTD*vd*Conj(MuD);
   mass_matrix_hh(1,1) = mHu2 + 0.3872983346207417*g1*MDBS*vS - 0.5*g2*
      MDWBT*vT + AbsSqr(MuU) + AbsSqr(Mu) + 0.7071067811865475*MuU*vS*Conj(
      LamSU) - 0.35355339059327373*LamTU*vS*vT*Conj(LamSU) - 0.5*MuU*vT*Conj(
      LamTU) - 0.35355339059327373*LamSU*vS*vT*Conj(LamTU) + 0.3872983346207417
      *g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*
      Conj(MuU) - 0.5*LamTU*vT*Conj(MuU) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2
      )*Sqr(vd) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) +
      0.225*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(1,2) = 0.3872983346207417*g1*MDBS*vu + vS*vu*AbsSqr(
      LamSU) + 0.7071067811865475*MuU*vu*Conj(LamSU) - 0.35355339059327373*
      LamTU*vT*vu*Conj(LamSU) - 0.35355339059327373*LamSU*vT*vu*Conj(LamTU) +
      0.3872983346207417*g1*vu*Conj(MDBS) + 0.7071067811865475*LamSU*vu*Conj(
      MuU);
   mass_matrix_hh(1,3) = -0.5*g2*MDWBT*vu + 0.5*vT*vu*AbsSqr(LamTU) -
      0.35355339059327373*LamTU*vS*vu*Conj(LamSU) - 0.5*MuU*vu*Conj(LamTU) -
      0.35355339059327373*LamSU*vS*vu*Conj(LamTU) - 0.5*g2*vu*Conj(MDWBT) - 0.5
      *LamTU*vu*Conj(MuU);
   mass_matrix_hh(2,2) = mS2 + 4*Sqr(MDBS) + 0.5*AbsSqr(LamSD)*Sqr(vd) +
      0.5*AbsSqr(LamSU)*Sqr(vu);
   mass_matrix_hh(2,3) = 0.17677669529663687*LamTD*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*Conj(LamTD)*Sqr(vd) - 0.17677669529663687*LamTU
      *Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*Conj(LamTU)*Sqr(vu);
   mass_matrix_hh(3,3) = mT2 + 4*Sqr(MDWBT) + 0.25*AbsSqr(LamTD)*Sqr(vd)
      + 0.25*AbsSqr(LamTU)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::hh, eigenvalue_error > precision *
      Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif
   normalize_to_interval(ZH);


   if (Mhh.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*
      MDWBT*vT + AbsSqr(MuD) + AbsSqr(Mu) + 0.7071067811865475*MuD*vS*Conj(
      LamSD) + 0.35355339059327373*LamTD*vS*vT*Conj(LamSD) + 0.5*MuD*vT*Conj(
      LamTD) + 0.35355339059327373*LamSD*vS*vT*Conj(LamTD) - 0.3872983346207417
      *g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*
      Conj(MuD) + 0.5*LamTD*vT*Conj(MuD) + 0.3872983346207417*g1*g2*Cos(ThetaW(
      ))*Sin(ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd)
      + 0.5*AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) - 0.075*Sqr(g1)*
      Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))
      + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*Conj(BMu) - 0.3872983346207417*g1*
      g2*vd*vu*Cos(ThetaW())*Sin(ThetaW()) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW(
      ))) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,2) = 0;
   mass_matrix_Ah(0,3) = 0;
   mass_matrix_Ah(1,1) = mHu2 + 0.3872983346207417*g1*MDBS*vS - 0.5*g2*
      MDWBT*vT + AbsSqr(MuU) + AbsSqr(Mu) + 0.7071067811865475*MuU*vS*Conj(
      LamSU) - 0.35355339059327373*LamTU*vS*vT*Conj(LamSU) - 0.5*MuU*vT*Conj(
      LamTU) - 0.35355339059327373*LamSU*vS*vT*Conj(LamTU) + 0.3872983346207417
      *g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*
      Conj(MuU) - 0.5*LamTU*vT*Conj(MuU) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2
      )*Sqr(vd) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu) + 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW
      ())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = 0;
   mass_matrix_Ah(1,3) = 0;
   mass_matrix_Ah(2,2) = mS2 + 0.5*AbsSqr(LamSD)*Sqr(vd) + 0.5*AbsSqr(
      LamSU)*Sqr(vu);
   mass_matrix_Ah(2,3) = 0.17677669529663687*LamTD*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*Conj(LamTD)*Sqr(vd) - 0.17677669529663687*LamTU
      *Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*Conj(LamTU)*Sqr(vu);
   mass_matrix_Ah(3,3) = mT2 + 0.25*AbsSqr(LamTD)*Sqr(vd) + 0.25*AbsSqr(
      LamTU)*Sqr(vu);

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Ah, eigenvalue_error > precision *
      Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif
   normalize_to_interval(ZA);


   if (MAh.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Rh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Rh;

   mass_matrix_Rh(0,0) = mRd2 + 0.3872983346207417*g1*MDBS*vS - 0.5*g2*
      MDWBT*vT + AbsSqr(MuD) + 0.7071067811865475*MuD*vS*Conj(LamSD) +
      0.35355339059327373*LamTD*vS*vT*Conj(LamSD) + 0.5*MuD*vT*Conj(LamTD) +
      0.35355339059327373*LamSD*vS*vT*Conj(LamTD) + 0.3872983346207417*g1*vS*
      Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*Conj(MuD
      ) + 0.5*LamTD*vT*Conj(MuD) + 0.5*AbsSqr(LamSD)*Sqr(vd) + 0.25*AbsSqr(
      LamTD)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*
      AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) + 0.075*Sqr(g1)*Sqr(vu
      ) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Rh(0,1) = -0.5*LamSU*vd*vu*Conj(LamSD) + 0.25*LamTU*vd*vu*
      Conj(LamTD);
   mass_matrix_Rh(1,1) = mRu2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*
      MDWBT*vT + AbsSqr(MuU) + 0.7071067811865475*MuU*vS*Conj(LamSU) -
      0.35355339059327373*LamTU*vS*vT*Conj(LamSU) - 0.5*MuU*vT*Conj(LamTU) -
      0.35355339059327373*LamSU*vS*vT*Conj(LamTU) - 0.3872983346207417*g1*vS*
      Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*Conj(MuU
      ) - 0.5*LamTU*vT*Conj(MuU) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) + 0.5*AbsSqr(
      LamSU)*Sqr(vu) + 0.25*AbsSqr(LamTU)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);

   Hermitianize(mass_matrix_Rh);

   return mass_matrix_Rh;
}

void CLASSNAME::calculate_MRh()
{
   const auto mass_matrix_Rh(get_mass_matrix_Rh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Rh, MRh, ZHR, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Rh, eigenvalue_error > precision *
      Abs(MRh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Rh, MRh, ZHR);
#endif
   normalize_to_interval(ZHR);


   if (MRh.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Rh);
   }

   MRh = AbsSqrt(MRh);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Hpm() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 - 0.3872983346207417*g1*MDBS*vS - 0.5*g2*
      MDWBT*vT + AbsSqr(MuD) + AbsSqr(Mu) + 0.7071067811865475*MuD*vS*Conj(
      LamSD) - 0.35355339059327373*LamTD*vS*vT*Conj(LamSD) - 0.5*MuD*vT*Conj(
      LamTD) - 0.35355339059327373*LamSD*vS*vT*Conj(LamTD) - 0.3872983346207417
      *g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*
      Conj(MuD) - 0.5*LamTD*vT*Conj(MuD) + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2
      )*Sqr(vd) + 0.5*AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr(LamTD)*Sqr(vT) -
      0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = Conj(BMu);
   mass_matrix_Hpm(0,2) = 0.7071067811865475*g2*MDWBT*vd -
      0.35355339059327373*vd*vT*AbsSqr(LamTD) + 0.5*LamTD*vd*vS*Conj(LamSD) +
      0.7071067811865475*LamTD*vd*Conj(MuD) + 0.7071067811865475*vd*vT*Sqr(g2);
   mass_matrix_Hpm(0,3) = 0.35355339059327373*vd*vT*AbsSqr(LamTD) +
      0.7071067811865475*MuD*vd*Conj(LamTD) + 0.5*LamSD*vd*vS*Conj(LamTD) +
      0.7071067811865475*g2*vd*Conj(MDWBT);
   mass_matrix_Hpm(1,0) = BMu;
   mass_matrix_Hpm(1,1) = mHu2 + 0.3872983346207417*g1*MDBS*vS + 0.5*g2*
      MDWBT*vT + AbsSqr(MuU) + AbsSqr(Mu) + 0.7071067811865475*MuU*vS*Conj(
      LamSU) + 0.35355339059327373*LamTU*vS*vT*Conj(LamSU) + 0.5*MuU*vT*Conj(
      LamTU) + 0.35355339059327373*LamSU*vS*vT*Conj(LamTU) + 0.3872983346207417
      *g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*
      Conj(MuU) + 0.5*LamTU*vT*Conj(MuU) - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2
      )*Sqr(vd) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) +
      0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(1,2) = 0.7071067811865475*g2*MDWBT*vu -
      0.35355339059327373*vT*vu*AbsSqr(LamTU) + 0.5*LamTU*vS*vu*Conj(LamSU) +
      0.7071067811865475*LamTU*vu*Conj(MuU);
   mass_matrix_Hpm(1,3) = 0.35355339059327373*vT*vu*AbsSqr(LamTU) +
      0.7071067811865475*MuU*vu*Conj(LamTU) + 0.5*LamSU*vS*vu*Conj(LamTU) +
      0.7071067811865475*g2*vu*Conj(MDWBT) - 0.7071067811865475*vT*vu*Sqr(g2);
   mass_matrix_Hpm(2,0) = -0.35355339059327373*vd*vT*AbsSqr(LamTD) +
      0.7071067811865475*MuD*vd*Conj(LamTD) + 0.5*LamSD*vd*vS*Conj(LamTD) +
      0.7071067811865475*g2*vd*Conj(MDWBT) + 0.7071067811865475*vd*vT*Sqr(g2);
   mass_matrix_Hpm(2,1) = -0.35355339059327373*vT*vu*AbsSqr(LamTU) +
      0.7071067811865475*MuU*vu*Conj(LamTU) + 0.5*LamSU*vS*vu*Conj(LamTU) +
      0.7071067811865475*g2*vu*Conj(MDWBT);
   mass_matrix_Hpm(2,2) = mT2 + 2*Sqr(MDWBT) + 0.5*AbsSqr(LamTD)*Sqr(vd)
      - 0.25*Sqr(g2)*Sqr(vd) + Sqr(g2)*Sqr(vT) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(2,3) = 2*Sqr(MDWBT);
   mass_matrix_Hpm(3,0) = 0.7071067811865475*g2*MDWBT*vd +
      0.35355339059327373*vd*vT*AbsSqr(LamTD) + 0.5*LamTD*vd*vS*Conj(LamSD) +
      0.7071067811865475*LamTD*vd*Conj(MuD);
   mass_matrix_Hpm(3,1) = 0.7071067811865475*g2*MDWBT*vu +
      0.35355339059327373*vT*vu*AbsSqr(LamTU) + 0.5*LamTU*vS*vu*Conj(LamSU) +
      0.7071067811865475*LamTU*vu*Conj(MuU) - 0.7071067811865475*vT*vu*Sqr(g2);
   mass_matrix_Hpm(3,2) = 2*Sqr(MDWBT);
   mass_matrix_Hpm(3,3) = mT2 + 2*Sqr(MDWBT) + 0.25*Sqr(g2)*Sqr(vd) + Sqr
      (g2)*Sqr(vT) + 0.5*AbsSqr(LamTU)*Sqr(vu) - 0.25*Sqr(g2)*Sqr(vu);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Hpm, eigenvalue_error > precision *
      Abs(MHpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif
   normalize_to_interval(ZP);


   if (MHpm.minCoeff() < 0.) {
      problems.flag_running_tachyon(MRSSM_info::Hpm);
   }

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Chi() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MDBS;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,0) = 0;
   mass_matrix_Chi(1,1) = MDWBT;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(2,0) = -0.7071067811865475*LamSD*vd;
   mass_matrix_Chi(2,1) = -0.5*LamTD*vd;
   mass_matrix_Chi(2,2) = -MuD - 0.7071067811865475*LamSD*vS - 0.5*LamTD*
      vT;
   mass_matrix_Chi(2,3) = 0;
   mass_matrix_Chi(3,0) = 0.7071067811865475*LamSU*vu;
   mass_matrix_Chi(3,1) = -0.5*LamTU*vu;
   mass_matrix_Chi(3,2) = 0;
   mass_matrix_Chi(3,3) = MuU + 0.7071067811865475*LamSU*vS - 0.5*LamTU*
      vT;

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Chi, MChi, ZN1, ZN2, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Chi, eigenvalue_error > precision *
      Abs(MChi(0)));
#else
   fs_svd(mass_matrix_Chi, MChi, ZN1, ZN2);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha1() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha1;

   mass_matrix_Cha1(0,0) = MDWBT + g2*vT;
   mass_matrix_Cha1(0,1) = 0.7071067811865475*LamTD*vd;
   mass_matrix_Cha1(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha1(1,1) = MuD + 0.7071067811865475*LamSD*vS - 0.5*LamTD*
      vT;

   return mass_matrix_Cha1;
}

void CLASSNAME::calculate_MCha1()
{
   const auto mass_matrix_Cha1(get_mass_matrix_Cha1());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha1, MCha1, UM1, UP1, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Cha1, eigenvalue_error > precision
      * Abs(MCha1(0)));
#else
   fs_svd(mass_matrix_Cha1, MCha1, UM1, UP1);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha2() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha2;

   mass_matrix_Cha2(0,0) = MDWBT - g2*vT;
   mass_matrix_Cha2(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha2(1,0) = -0.7071067811865475*LamTU*vu;
   mass_matrix_Cha2(1,1) = -MuU - 0.7071067811865475*LamSU*vS - 0.5*LamTU
      *vT;

   return mass_matrix_Cha2;
}

void CLASSNAME::calculate_MCha2()
{
   const auto mass_matrix_Cha2(get_mass_matrix_Cha2());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha2, MCha2, UM2, UP2, eigenvalue_error);
   problems.flag_bad_mass(MRSSM_info::Cha2, eigenvalue_error > precision
      * Abs(MCha2(0)));
#else
   fs_svd(mass_matrix_Cha2, MCha2, UM2, UP2);
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
   problems.flag_bad_mass(MRSSM_info::Fe, eigenvalue_error > precision *
      Abs(MFe(0)));
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
   problems.flag_bad_mass(MRSSM_info::Fd, eigenvalue_error > precision *
      Abs(MFd(0)));
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
   problems.flag_bad_mass(MRSSM_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(vd) + 4*Sqr(vT) +
      Sqr(vu)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_running_tachyon(MRSSM_info::VWm);
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
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ, eigenvalue_error
      );
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
   double result = Re(mHd2*vd - 0.3872983346207417*g1*MDBS*vd*vS + 0.5*g2*MDWBT
      *vd*vT + vd*AbsSqr(MuD) + vd*AbsSqr(Mu) - 0.5*vu*BMu + 0.7071067811865475*
      MuD*vd*vS*Conj(LamSD) + 0.35355339059327373*LamTD*vd*vS*vT*Conj(LamSD) + 0.5
      *MuD*vd*vT*Conj(LamTD) + 0.35355339059327373*LamSD*vd*vS*vT*Conj(LamTD) -
      0.3872983346207417*g1*vd*vS*Conj(MDBS) + 0.5*g2*vd*vT*Conj(MDWBT) +
      0.7071067811865475*LamSD*vd*vS*Conj(MuD) + 0.5*LamTD*vd*vT*Conj(MuD) - 0.5*
      vu*Conj(BMu) + 0.075*Cube(vd)*Sqr(g1) + 0.125*Cube(vd)*Sqr(g2) + 0.5*vd*
      AbsSqr(LamSD)*Sqr(vS) + 0.25*vd*AbsSqr(LamTD)*Sqr(vT) - 0.075*vd*Sqr(g1)*Sqr
      (vu) - 0.125*vd*Sqr(g2)*Sqr(vu));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = Re(mHu2*vu + 0.3872983346207417*g1*MDBS*vS*vu - 0.5*g2*MDWBT
      *vT*vu + vu*AbsSqr(MuU) + vu*AbsSqr(Mu) - 0.5*vd*BMu + 0.7071067811865475*
      MuU*vS*vu*Conj(LamSU) - 0.35355339059327373*LamTU*vS*vT*vu*Conj(LamSU) - 0.5
      *MuU*vT*vu*Conj(LamTU) - 0.35355339059327373*LamSU*vS*vT*vu*Conj(LamTU) +
      0.3872983346207417*g1*vS*vu*Conj(MDBS) - 0.5*g2*vT*vu*Conj(MDWBT) +
      0.7071067811865475*LamSU*vS*vu*Conj(MuU) - 0.5*LamTU*vT*vu*Conj(MuU) - 0.5*
      vd*Conj(BMu) + 0.075*Cube(vu)*Sqr(g1) + 0.125*Cube(vu)*Sqr(g2) - 0.075*vu*
      Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*vu*AbsSqr(LamSU)*Sqr(vS) +
      0.25*vu*AbsSqr(LamTU)*Sqr(vT));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   double result = Re(mT2*vT + 4*vT*Sqr(MDWBT) + 0.25*g2*MDWBT*Sqr(vd) + 0.25*
      vT*AbsSqr(LamTD)*Sqr(vd) + 0.17677669529663687*LamTD*vS*Conj(LamSD)*Sqr(vd)
      + 0.25*MuD*Conj(LamTD)*Sqr(vd) + 0.17677669529663687*LamSD*vS*Conj(LamTD)*
      Sqr(vd) + 0.25*g2*Conj(MDWBT)*Sqr(vd) + 0.25*LamTD*Conj(MuD)*Sqr(vd) - 0.25*
      g2*MDWBT*Sqr(vu) + 0.25*vT*AbsSqr(LamTU)*Sqr(vu) - 0.17677669529663687*LamTU
      *vS*Conj(LamSU)*Sqr(vu) - 0.25*MuU*Conj(LamTU)*Sqr(vu) - 0.17677669529663687
      *LamSU*vS*Conj(LamTU)*Sqr(vu) - 0.25*g2*Conj(MDWBT)*Sqr(vu) - 0.25*LamTU*
      Conj(MuU)*Sqr(vu));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_4() const
{
   double result = Re(mS2*vS + 4*vS*Sqr(MDBS) - 0.19364916731037085*g1*MDBS*Sqr
      (vd) + 0.5*vS*AbsSqr(LamSD)*Sqr(vd) + 0.35355339059327373*MuD*Conj(LamSD)*
      Sqr(vd) + 0.17677669529663687*LamTD*vT*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*vT*Conj(LamTD)*Sqr(vd) - 0.19364916731037085*g1*
      Conj(MDBS)*Sqr(vd) + 0.35355339059327373*LamSD*Conj(MuD)*Sqr(vd) +
      0.19364916731037085*g1*MDBS*Sqr(vu) + 0.5*vS*AbsSqr(LamSU)*Sqr(vu) +
      0.35355339059327373*MuU*Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamTU*vT*
      Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*vT*Conj(LamTU)*Sqr(vu) +
      0.19364916731037085*g1*Conj(MDBS)*Sqr(vu) + 0.35355339059327373*LamSU*Conj(
      MuU)*Sqr(vu));

   return result;
}



std::complex<double> CLASSNAME::CpSRdpUSdconjSRdpconjUSd(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)
      *Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2),0) - 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUSdconjSRumconjUSd(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2),0) + 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

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

std::complex<double> CLASSNAME::CpRhUSdconjRhconjUSd(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)
      *Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1),0) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZHR(gI1,0)*ZHR(gI2,0) + 0.1*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZHR(gI1,1)*
      ZHR(gI2,1);

   return result;
}

double CLASSNAME::CpbarCha1FuconjUSdPR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1FuconjUSdPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = Conj(UM1(gI1,1))*SUM(j2,0,2,Conj(ZUL(gI2
      ,j2))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yd(j1,j2)));

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

std::complex<double> CLASSNAME::CpCha2FuconjUSdPR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Yu(j1,gO2))*
      ZUR(gI1,j1))*UP2(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpCha2FuconjUSdPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UM2(gI2,0))*Conj(
      ZUL(gI1,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO2))
      *ZDR(gI1,j1))*ZN2(gI2,2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFdconjUSdPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(
      ZDL(gI1,gO1))*Conj(ZN1(gI2,0)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(
      ZDL(gI1,gO1))*Conj(ZN1(gI2,1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Yd(j1,gO2))*Yd(j1,gO1))*ZA(gI1,0)*ZA(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + 0.1*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,
      0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))))*ZA(gI1,0)*ZA(gI2,0) - 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1
      ,1)*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSdconjUSd(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Yd(j1,gO2))*Yd(j1,gO1))*ZH(gI1,0)*ZH(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + 0.1*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,
      0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))))*ZH(gI1,0)*ZH(gI2,0) - 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1
      ,1)*ZH(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSdconjHpmconjUSd(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Yu(j1,gO2))*Yu(j1,gO1))*ZP(gI1,1)*ZP(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2),0) + IF(gO1 < 3,0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3),0) + 0.1*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,
      0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))))*ZP(gI1,0)*ZP(gI2,0) - 0.1*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1
      ,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiFdconjUSdPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -0.3651483716701107*g1*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1))*ZN1(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiFdconjUSdPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -(Conj(ZN2(gI1,2))*SUM(j2,0,2,Conj(ZDL(
      gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yd(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpphiOSdconjUSd(int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(g3*MDGoc*Conj(ZD(gI1,gO2)))
      ,0) + IF(gO2 < 3,-(g3*Conj(MDGoc)*Conj(ZD(gI1,gO2))),0) + g3*(MDGoc + Conj(
      MDGoc))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSRumSuconjUSd(int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(MuU*SUM(j1,0,2,Conj(Yu(j1,
      gO2))*Conj(ZU(gI1,3 + j1)))),0) + IF(gO2 < 3,-0.7071067811865475*LamSU*vS*
      SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1))),0) + IF(gO2 < 3,-0.5*LamTU
      *vT*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yd(j1
      ,gO1)*ZD(gI1,3 + j1))*SUM(j3,0,2,Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3)))),0),
      0) + IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*
      ZD(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD
      (gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-1.3333333333333333*Conj(ZD(gI2,gO2)
      )*Sqr(g3)*ZD(gI1,gO1),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 +
      j1))*ZD(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,3 +
      j2))*ZD(gI1,3 + j2)),0) + IF(gO1 < 3,-3*SUM(j1,0,2,KroneckerDelta(gO2,3 +
      j1)*Yd(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3)))*
      ZD(gI1,j4)),0) + IF(gO1 < 3,-0.016666666666666666*Sqr(g1)*SUM(j1,0,2,Conj(ZD
      (gI2,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,
      0.6666666666666666*Sqr(g3)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(
      gO2,3 + j1))*ZD(gI1,gO1),0) + IF(gO1 < 3,-0.016666666666666666*Sqr(g1)*SUM(
      j2,0,2,Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(
      gO1 < 3,0.6666666666666666*Sqr(g3)*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*
      KroneckerDelta(gO2,3 + j2))*ZD(gI1,gO1),0) + IF(gO2 < 3,
      -0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*ZD(gI1,3 + j1)),0) + IF(gO2 < 3,0.6666666666666666*Conj(ZD(gI2,
      gO2))*Sqr(g3)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1)),0) + IF(
      gO2 < 3,-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*SUM(j2,0,2,
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
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yu(j1
      ,gO1)*ZU(gI1,3 + j1))*SUM(j3,0,2,Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3)))),0),
      0) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1),0),0)
      + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,
      j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(
      j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,
      gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3
      ,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,
      j2)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(
      ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZU(gI2,j1))*ZU(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
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
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 +
      j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,-0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 +
      j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1,0,2,KroneckerDelta(gO2,3 + j1
      )*Yd(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*Conj(ZE(gI1,3 + j3)))*ZE
      (gI2,j4))),0) + IF(gO2 < 3,-(SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Ye(j1,j2
      )*ZE(gI2,3 + j1)))*SUM(j3,0,2,Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3))),
      0) + 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,2,
      KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,
      0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 +
      j2)*KroneckerDelta(gO2,3 + j2)) + 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2))
      - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpRhSdconjUSd(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,MuD*SUM(j1,0,2,Conj(Yd(j1,gO2
      ))*Conj(ZD(gI1,3 + j1)))*ZHR(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*LamSD
      *vS*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1)))*ZHR(gI2,0),0) + IF(gO2
      < 3,0.5*LamTD*vT*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1)))*ZHR(gI2,
      0),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhSdconjUSd(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*Mu*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1)))*ZA(
      gI2,1),0) + IF(gO2 < 3,std::complex<double>(0.,-0.12909944487358055)*g1*MDBS
      *Conj(ZD(gI1,gO2))*ZA(gI2,2),0) + IF(gO2 < 3,std::complex<double>(0.,
      0.12909944487358055)*g1*Conj(MDBS)*Conj(ZD(gI1,gO2))*ZA(gI2,2),0) + IF(gO2 <
      3,std::complex<double>(0,0.5)*g2*MDWBT*Conj(ZD(gI1,gO2))*ZA(gI2,3),0) + IF(
      gO2 < 3,std::complex<double>(0,-0.5)*g2*Conj(MDWBT)*Conj(ZD(gI1,gO2))*ZA(gI2
      ,3),0) - std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*SUM(j2,0,2,
      Conj(ZD(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZA(gI2,1)
      - std::complex<double>(0.,0.2581988897471611)*g1*MDBS*SUM(j1,0,2,Conj(ZD(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZA(gI2,2) + std::complex<double>(0.
      ,0.2581988897471611)*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSdconjUSd(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.05*vd*Conj(ZD(gI1,gO2))*Sqr
      (g1)*ZH(gI2,0),0) + IF(gO2 < 3,0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0),0
      ) + IF(gO2 < 3,-(vd*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2))*
      Yd(j1,j2)))*ZH(gI2,0)),0) + IF(gO2 < 3,-0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH
      (gI2,1),0) + IF(gO2 < 3,-0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF
      (gO2 < 3,0.7071067811865475*Mu*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 +
      j1)))*ZH(gI2,1),0) + IF(gO2 < 3,-0.12909944487358055*g1*MDBS*Conj(ZD(gI1,gO2
      ))*ZH(gI2,2),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(MDBS)*Conj(ZD(gI1,
      gO2))*ZH(gI2,2),0) + IF(gO2 < 3,0.5*g2*MDWBT*Conj(ZD(gI1,gO2))*ZH(gI2,3),0)
      + IF(gO2 < 3,0.5*g2*Conj(MDWBT)*Conj(ZD(gI1,gO2))*ZH(gI2,3),0) + 0.1*vd*Sqr(
      g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) -
      vd*SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM
      (j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))))*ZH(gI2,0) - 0.1*vu*Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) +
      0.7071067811865475*Conj(Mu)*SUM(j2,0,2,Conj(ZD(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZH(gI2,1) - 0.2581988897471611*g1*
      MDBS*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,2) -
      0.2581988897471611*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmSuconjUSd(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(
      ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,Mu*SUM(j1,0,2,Conj(Yu(j1,gO2)
      )*Conj(ZU(gI1,3 + j1)))*ZP(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(
      j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,j2)))*ZP(gI2,0),0)
      + IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1),0)
      + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,
      Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZP(gI2,1),0) + IF(gO2 < 3,-(g2*MDWBT*Conj(ZU(
      gI1,gO2))*ZP(gI2,2)),0) + IF(gO2 < 3,-0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(
      gI2,2),0) + IF(gO2 < 3,-(g2*Conj(MDWBT)*Conj(ZU(gI1,gO2))*ZP(gI2,3)),0) + IF
      (gO2 < 3,0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,3),0) + 0.7071067811865475*
      vu*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM
      (j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))))*ZP(gI2,0) + Conj(Mu)*SUM(j2,0,2,Conj(ZU
      (gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)))*ZP(gI2,1) +
      0.7071067811865475*vd*SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yd(j2,j1))))*ZP(gI2,1)
      ;

   return result;
}

double CLASSNAME::CpGluFdconjUSdPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpGluFdconjUSdPL(int gI2, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      ZDL(gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarGluFdconjUSdPR(int gI2, int gO2) const
{
   const std::complex<double> result = 1.4142135623730951*g3*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1));

   return result;
}

double CLASSNAME::CpbarGluFdconjUSdPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOSdconjUSd(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,-1)*g3
      *MDGoc*Conj(ZD(gI2,gO2)),0) + IF(gO2 < 3,std::complex<double>(0,1)*g3*Conj(
      MDGoc)*Conj(ZD(gI2,gO2)),0) + std::complex<double>(0,1)*g3*(MDGoc - Conj(
      MDGoc))*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSRdpconjUSd(int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) -
      vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVG(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gI2 < 6,g3*Conj(ZD(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVP(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.12909944487358055*g1*Conj(
      ZD(gI2,gO2))*Cos(ThetaW()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZD(gI2,gO2))*Sin(
      ThetaW()),0) - 0.2581988897471611*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSdVZ(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Conj(ZD(gI2,gO2))*Cos
      (ThetaW()),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Sin(
      ThetaW()),0) + 0.2581988897471611*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSdVWm(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZU
      (gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUSvconjSRdpconjUSv(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1)
      - 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUSvconjSRumconjUSv(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*(-3*Sqr(g1)
      + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = 0.1*KroneckerDelta(gO1,gO2)*(g1*Sin(
      ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(int gO1, int gO2) const
{
   const double result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpRhUSvconjRhconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1)
      + 5*Sqr(g2))*(ZHR(gI1,0)*ZHR(gI2,0) - ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSvUSvconjSvconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gI1 < 3,IF(gI2 < 3,-0.15*Conj(ZV(gI1,
      gO2))*Sqr(g1)*ZV(gI2,gO1),0),0) + IF(gI1 < 3,IF(gI2 < 3,-0.25*Conj(ZV(gI1,
      gO2))*Sqr(g2)*ZV(gI2,gO1),0),0) - 0.05*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpCha1FeconjUSvPR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Ye(j1,gO2))*
      ZER(gI1,j1))*UM1(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpCha1FeconjUSvPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UP1(gI2,0))*Conj(
      ZEL(gI1,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhSvconjUSv(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gI1 < 3,std::complex<double>(0.,
      0.3872983346207417)*g1*MDBS*Conj(ZV(gI1,gO2))*ZA(gI2,2),0) + IF(gI1 < 3,
      std::complex<double>(0.,-0.3872983346207417)*g1*Conj(MDBS)*Conj(ZV(gI1,gO2))
      *ZA(gI2,2),0) + IF(gI1 < 3,std::complex<double>(0,-0.5)*g2*MDWBT*Conj(ZV(gI1
      ,gO2))*ZA(gI2,3),0) + IF(gI1 < 3,std::complex<double>(0,0.5)*g2*Conj(MDWBT)*
      Conj(ZV(gI1,gO2))*ZA(gI2,3),0);

   return result;
}

std::complex<double> CLASSNAME::CphhSvconjUSv(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gI1 < 3,-0.15*vd*Conj(ZV(gI1,gO2))*
      Sqr(g1)*ZH(gI2,0),0) + IF(gI1 < 3,-0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,
      0),0) + IF(gI1 < 3,0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1),0) + IF(gI1 <
      3,0.25*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) + IF(gI1 < 3,
      0.3872983346207417*g1*MDBS*Conj(ZV(gI1,gO2))*ZH(gI2,2),0) + IF(gI1 < 3,
      0.3872983346207417*g1*Conj(MDBS)*Conj(ZV(gI1,gO2))*ZH(gI2,2),0) + IF(gI1 < 3
      ,-0.5*g2*MDWBT*Conj(ZV(gI1,gO2))*ZH(gI2,3),0) + IF(gI1 < 3,-0.5*g2*Conj(
      MDWBT)*Conj(ZV(gI1,gO2))*ZH(gI2,3),0);

   return result;
}

double CLASSNAME::CpChiFvconjUSvPR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpChiFvconjUSvPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gI1 < 3,0.5477225575051661*g1*Conj(
      ZN1(gI2,0))*KroneckerDelta(gI1,gO1),0) + IF(gI1 < 3,-0.7071067811865475*g2*
      Conj(ZN1(gI2,1))*KroneckerDelta(gI1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1)
      + 5*Sqr(g2))*(ZA(gI1,0)*ZA(gI2,0) - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSvconjUSv(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1)
      + 5*Sqr(g2))*(ZH(gI1,0)*ZH(gI2,0) - ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSvconjHpmconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*(20*IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,
      0,2,Conj(Ye(j1,gO2))*Ye(j1,gO1))*ZP(gI1,0)*ZP(gI2,0)),0),0) + KroneckerDelta
      (gO1,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (3*Sqr(g1) - 5*Sqr
      (g2))*ZP(gI1,1)*ZP(gI2,1) + 10*Sqr(g2)*(ZP(gI1,2)*ZP(gI2,2) - ZP(gI1,3)*ZP(
      gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjHpmconjUSv(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(
      ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,
      0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,j2)))*ZP(gI1,0),0) +
      IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) +
      IF(gO2 < 3,Mu*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1)))*ZP(gI1,1),0)
      + IF(gO2 < 3,-(g2*Conj(MDWBT)*Conj(ZE(gI2,gO2))*ZP(gI1,2)),0) + IF(gO2 < 3,
      -0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,2),0) + IF(gO2 < 3,-(g2*MDWBT*Conj(
      ZE(gI2,gO2))*ZP(gI1,3)),0) + IF(gO2 < 3,0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(
      gI1,3),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdUSvconjSdconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((Sqr(g1) +
      5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,
      Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSvconjSeconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Ye(j1
      ,gO1)*ZE(gI2,3 + j1))*SUM(j3,0,2,Conj(Ye(j3,gO2))*Conj(ZE(gI1,3 + j3)))),0),
      0) + 0.05*(20*IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZE(gI1,gO2))*Sqr(g2)*ZE(gI2,
      gO1),0),0) + KroneckerDelta(gO1,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,
      Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(
      gI2,3 + j1))));

   return result;
}

std::complex<double> CLASSNAME::CpSuUSvconjSuconjUSv(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*KroneckerDelta(gO1,gO2)*((Sqr(g1) -
      5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,
      Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSvVZ(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZV(gI2,gO2))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> CLASSNAME::CpSRdpSeconjUSv(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,MuD*SUM(j1,0,2,Conj(Ye(j1,gO2
      ))*Conj(ZE(gI2,3 + j1))),0) + IF(gO2 < 3,0.7071067811865475*LamSD*vS*SUM(j1,
      0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1))),0) + IF(gO2 < 3,-0.5*LamTD*vT*SUM
      (j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSvconjVWm(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZE
      (gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUSuconjSRdpconjUSu(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)
      *Sqr(g1),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2),0) + 0.2*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUSuconjSRumconjUSu(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2),0) - 0.2*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.2581988897471611*g1*g2*Cos
      (ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW()),0) + IF(gO1 < 3,0.5*
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

std::complex<double> CLASSNAME::CpRhUSuconjRhconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0),0) + IF(gO1 < 3,0.05*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1),0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZHR(gI1,0)*ZHR(gI2,0) - 0.2*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZHR(gI1,1)*
      ZHR(gI2,1);

   return result;
}

double CLASSNAME::CpbarCha2FdconjUSuPR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2FdconjUSuPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = Conj(UP2(gI1,1))*SUM(j2,0,2,Conj(ZDL(gI2
      ,j2))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yu(j1,j2)));

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

std::complex<double> CLASSNAME::CpCha1FdconjUSuPR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,SUM(j1,0,2,Conj(Yd(j1,gO2))*
      ZDR(gI1,j1))*UM1(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpCha1FdconjUSuPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(UP1(gI2,0))*Conj(
      ZDL(gI1,gO1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO2))
      *ZUR(gI1,j1))*ZN2(gI2,3)),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFuconjUSuPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(
      ZN1(gI2,0))*Conj(ZUL(gI1,gO1)),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      ZN1(gI2,1))*Conj(ZUL(gI1,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Yu(j1,gO2))*Yu(j1,gO1))*ZA(gI1,1)*ZA(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) - 0.2*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,
      0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
      + j1))*ZA(gI1,1)*ZA(gI2,1) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,
      2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZA(gI1,
      1)*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSuconjUSu(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Yu(j1,gO2))*Yu(j1,gO1))*ZH(gI1,1)*ZH(gI2,1)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) - 0.2*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,
      0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
      + j1))*ZH(gI1,1)*ZH(gI2,1) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,
      2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZH(gI1,
      1)*ZH(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSuconjHpmconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Yd(j1,gO2))*Yd(j1,gO1))*ZP(gI1,0)*ZP(gI2,0)),0),0) + IF(gO1 < 3,0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.05*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2),0) + IF(gO1 < 3,-0.5*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3),0) - 0.2*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,
      0) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3
      + j1))*ZP(gI1,1)*ZP(gI2,1) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,
      2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))))*ZP(gI1,
      1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiFuconjUSuPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.7302967433402214*g1*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1))*ZN1(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiFuconjUSuPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -(Conj(ZN2(gI1,3))*SUM(j2,0,2,Conj(ZUL(
      gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Yu(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjHpmconjUSu(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(
      ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,
      0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yd(j1,gO2))*Yd(j1,j2)))*ZP(gI1,0),0) +
      IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1),0) +
      IF(gO2 < 3,Mu*SUM(j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1)))*ZP(gI1,1),0)
      + IF(gO2 < 3,0.7071067811865475*vu*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      Conj(Yu(j1,gO2))*Yu(j1,j2)))*ZP(gI1,1),0) + IF(gO2 < 3,-(g2*Conj(MDWBT)*Conj
      (ZD(gI2,gO2))*ZP(gI1,2)),0) + IF(gO2 < 3,-0.5*vT*Conj(ZD(gI2,gO2))*Sqr(g2)*
      ZP(gI1,2),0) + IF(gO2 < 3,-(g2*MDWBT*Conj(ZD(gI2,gO2))*ZP(gI1,3)),0) + IF(
      gO2 < 3,0.5*vT*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,3),0) + Conj(Mu)*SUM(j2,0,2,
      Conj(ZD(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZP(gI1,0)
      + 0.7071067811865475*vu*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))))*ZP(gI1,0)
      + 0.7071067811865475*vd*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,
      KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Yd(j3,j1))*Yu(j2,j1))))*ZP(gI1,1)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpphiOSuconjUSu(int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(g3*MDGoc*Conj(ZU(gI1,gO2)))
      ,0) + IF(gO2 < 3,-(g3*Conj(MDGoc)*Conj(ZU(gI1,gO2))),0) + g3*(MDGoc + Conj(
      MDGoc))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOSuconjUSu(int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0,-1)*g3
      *MDGoc*Conj(ZU(gI1,gO2)),0) + IF(gO2 < 3,std::complex<double>(0,1)*g3*Conj(
      MDGoc)*Conj(ZU(gI1,gO2)),0) + std::complex<double>(0,1)*g3*(MDGoc - Conj(
      MDGoc))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeUSuconjSeconjUSu(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) +
      IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 +
      j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta
      (gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,
      -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2
      ,3 + j2)),0) - 0.1*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))*SUM(j2,0,
      2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(
      j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) - 0.1*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,
      j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,
      3 + j1))*SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yd(j1
      ,gO1)*ZD(gI1,3 + j1))*SUM(j3,0,2,Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3)))),0),
      0) + IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1),0),0)
      + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,
      j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(
      j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1
      ,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)),0) + IF(gO1 <
      3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZD(gI2,j2))*ZD(gI1,
      j2)),0) + IF(gO1 < 3,0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(
      ZD(gI2,j2))*ZD(gI1,j2)),0) + IF(gO1 < 3,-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1
      )*SUM(j2,0,2,Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2)),0) + 0.1*Sqr(g1)*SUM(j1,0,
      2,Conj(ZD(gI2,j1))*ZD(gI1,j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*
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
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Yu(j1
      ,gO1)*ZU(gI1,3 + j1))*SUM(j3,0,2,Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3)))),0),
      0) + IF(gO1 < 3,IF(gO2 < 3,-0.016666666666666666*Conj(ZU(gI2,gO2))*Sqr(g1)*
      ZU(gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU
      (gI1,gO1),0),0) + IF(gO1 < 3,IF(gO2 < 3,-1.3333333333333333*Conj(ZU(gI2,gO2)
      )*Sqr(g3)*ZU(gI1,gO1),0),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(
      g1)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 +
      j1))*ZU(gI1,3 + j1)),0) + IF(gO1 < 3,-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) + IF(gO1 < 3,-0.375*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI2,j2))*ZU(gI1,j2)),0) +
      IF(gO1 < 3,0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI2,3 +
      j2))*ZU(gI1,3 + j2)),0) + IF(gO1 < 3,-3*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1
      )*Yu(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3)))*ZU
      (gI1,j4)),0) + IF(gO1 < 3,0.03333333333333333*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2
      ,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZU(gI1,gO1),0) + IF(gO1 < 3,
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
      j2)) - SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu
      (j1,j2)))*SUM(j4,0,2,SUM(j3,0,2,Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3))*
      ZU(gI1,j4));

   return result;
}

std::complex<double> CLASSNAME::CpRhSuconjUSu(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(MuU*SUM(j1,0,2,Conj(Yu(j1,
      gO2))*Conj(ZU(gI1,3 + j1)))*ZHR(gI2,1)),0) + IF(gO2 < 3,-0.7071067811865475*
      LamSU*vS*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZHR(gI2,1),0) +
      IF(gO2 < 3,0.5*LamTU*vT*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*
      ZHR(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhSuconjUSu(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*Mu*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1)))*ZA(
      gI2,0),0) + IF(gO2 < 3,std::complex<double>(0.,-0.12909944487358055)*g1*MDBS
      *Conj(ZU(gI1,gO2))*ZA(gI2,2),0) + IF(gO2 < 3,std::complex<double>(0.,
      0.12909944487358055)*g1*Conj(MDBS)*Conj(ZU(gI1,gO2))*ZA(gI2,2),0) + IF(gO2 <
      3,std::complex<double>(0,-0.5)*g2*MDWBT*Conj(ZU(gI1,gO2))*ZA(gI2,3),0) + IF
      (gO2 < 3,std::complex<double>(0,0.5)*g2*Conj(MDWBT)*Conj(ZU(gI1,gO2))*ZA(gI2
      ,3),0) - std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*SUM(j2,0,2,
      Conj(ZU(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZA(gI2,0)
      + std::complex<double>(0.,0.5163977794943222)*g1*MDBS*SUM(j1,0,2,Conj(ZU(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZA(gI2,2) - std::complex<double>(0.
      ,0.5163977794943222)*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSuconjUSu(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.05*vd*Conj(ZU(gI1,gO2))*Sqr
      (g1)*ZH(gI2,0),0) + IF(gO2 < 3,-0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0),
      0) + IF(gO2 < 3,0.7071067811865475*Mu*SUM(j1,0,2,Conj(Yu(j1,gO2))*Conj(ZU(
      gI1,3 + j1)))*ZH(gI2,0),0) + IF(gO2 < 3,-0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*
      ZH(gI2,1),0) + IF(gO2 < 3,0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) +
      IF(gO2 < 3,-(vu*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,gO2))*Yu(
      j1,j2)))*ZH(gI2,1)),0) + IF(gO2 < 3,-0.12909944487358055*g1*MDBS*Conj(ZU(gI1
      ,gO2))*ZH(gI2,2),0) + IF(gO2 < 3,-0.12909944487358055*g1*Conj(MDBS)*Conj(ZU(
      gI1,gO2))*ZH(gI2,2),0) + IF(gO2 < 3,-0.5*g2*MDWBT*Conj(ZU(gI1,gO2))*ZH(gI2,3
      ),0) + IF(gO2 < 3,-0.5*g2*Conj(MDWBT)*Conj(ZU(gI1,gO2))*ZH(gI2,3),0) - 0.2*
      vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(
      gI2,0) + 0.7071067811865475*Conj(Mu)*SUM(j2,0,2,Conj(ZU(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)))*ZH(gI2,0) + 0.2*vu*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) - vu*SUM(j3,0,2,
      Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(
      Yu(j3,j1))*Yu(j2,j1))))*ZH(gI2,1) + 0.5163977794943222*g1*MDBS*SUM(j1,0,2,
      Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,2) +
      0.5163977794943222*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,2);

   return result;
}

double CLASSNAME::CpGluFuconjUSuPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpGluFuconjUSuPL(int gI2, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      ZUL(gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarGluFuconjUSuPR(int gI2, int gO2) const
{
   const std::complex<double> result = 1.4142135623730951*g3*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1));

   return result;
}

double CLASSNAME::CpbarGluFuconjUSuPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSRdpSdconjUSu(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,MuD*SUM(j1,0,2,Conj(Yd(j1,gO2
      ))*Conj(ZD(gI2,3 + j1))),0) + IF(gO2 < 3,0.7071067811865475*LamSD*vS*SUM(j1,
      0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1))),0) + IF(gO2 < 3,-0.5*LamTD*vT*SUM
      (j1,0,2,Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1))),0);

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSRumconjUSu(int gI2, int gO2) const
{
   const std::complex<double> result = -0.5*(1.4142135623730951*vS*Conj(LamSU)
      + vT*Conj(LamTU) + 2*Conj(MuU))*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUSuconjVWm(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZD
      (gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVG(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gI2 < 6,g3*Conj(ZU(gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSuconjUSuVP(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.12909944487358055*g1*Conj(
      ZU(gI2,gO2))*Cos(ThetaW()),0) + IF(gO2 < 3,0.5*g2*Conj(ZU(gI2,gO2))*Sin(
      ThetaW()),0) + 0.5163977794943222*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

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

std::complex<double> CLASSNAME::CpSRdpUSeconjSRdpconjUSe(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2),0) - 0.3*Sqr(g1
      )*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUSeconjSRumconjUSe(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)
      *Sqr(g1),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*Sqr(g2),0) + 0.3*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7745966692414834*g1*g2*Cos
      (ThetaW())*KroneckerDelta(gO1,gO2)*Sin(ThetaW()),0) + IF(gO1 < 3,0.5*
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

std::complex<double> CLASSNAME::CpRhUSeconjRhconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)*
      Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1),0) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1
      ,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZHR(gI1,0)*ZHR(gI2,0) + 0.3*Sqr(g1)*SUM
      (j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZHR(gI1,1)*
      ZHR(gI2,1);

   return result;
}

double CLASSNAME::CpbarCha1FvconjUSePR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1FvconjUSePL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = Conj(UM1(gI1,1))*SUM(j1,0,2,
      KroneckerDelta(gO1,3 + j1)*Ye(j1,gI2));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSvconjUSeconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-0.5*Conj(ZV(gI1,
      gO2))*Sqr(g2)*ZV(gI2,gO1),0),0) + IF(gO1 < 3,-0.15*KroneckerDelta(gI1,gI2)*
      KroneckerDelta(gO1,gO2)*Sqr(g1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gI1,gI2)
      *KroneckerDelta(gO1,gO2)*Sqr(g2),0) + 0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)*
      SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1)) - SUM(j2,0
      ,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*SUM(j4
      ,0,2,SUM(j3,0,2,Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZV(gI2,j4));

   return result;
}

double CLASSNAME::CpCha2FvconjUSePR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpCha2FvconjUSePL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gI1 < 3,-(g2*Conj(UM2(gI2,0))*
      KroneckerDelta(gI1,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpHpmSvconjUSe(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.35355339059327373*vd*Conj(
      ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*vd*SUM(j2,
      0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2))*Ye(j1,j2)))*ZP(gI2,0),0) +
      IF(gO2 < 3,-0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1),0) +
      IF(gO2 < 3,-(g2*MDWBT*Conj(ZV(gI1,gO2))*ZP(gI2,2)),0) + IF(gO2 < 3,-0.5*vT*
      Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,2),0) + IF(gO2 < 3,-(g2*Conj(MDWBT)*Conj(ZV
      (gI1,gO2))*ZP(gI2,3)),0) + IF(gO2 < 3,0.5*vT*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(
      gI2,3),0) + Conj(Mu)*SUM(j2,0,2,Conj(ZV(gI1,j2))*SUM(j1,0,2,KroneckerDelta(
      gO2,3 + j1)*Ye(j1,j2)))*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(SUM(j1,0,2,Conj(Ye(j1,gO2))
      *ZER(gI1,j1))*ZN2(gI2,2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpChiFeconjUSePL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(
      ZEL(gI1,gO1))*Conj(ZN1(gI2,0)),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(
      ZEL(gI1,gO1))*Conj(ZN1(gI2,1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Ye(j1,gO2))*Ye(j1,gO1))*ZA(gI1,0)*ZA(gI2,0)),0),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1),0) + 0.3*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1,0)*ZA(gI2,
      0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))))*ZA(gI1,0)*ZA(gI2,0) - 0.3*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZA(gI1
      ,1)*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CphhhhUSeconjUSe(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Conj(
      Ye(j1,gO2))*Ye(j1,gO1))*ZH(gI1,0)*ZH(gI2,0)),0),0) + IF(gO1 < 3,-0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0),0) + IF(gO1 < 3,0.15*
      KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1),0) + IF(gO1 < 3,-0.25*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1),0) + 0.3*Sqr(g1)*SUM(j1,
      0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1,0)*ZH(gI2,
      0) - SUM(j3,0,2,KroneckerDelta(gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 +
      j2)*SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))))*ZH(gI1,0)*ZH(gI2,0) - 0.3*Sqr(
      g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZH(gI1
      ,1)*ZH(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpHpmUSeconjHpmconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.15*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,-0.25*KroneckerDelta(gO1,gO2)*
      Sqr(g2)*ZP(gI1,0)*ZP(gI2,0),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr
      (g1)*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,0.25*KroneckerDelta(gO1,gO2)*Sqr(g2
      )*ZP(gI1,1)*ZP(gI2,1),0) + IF(gO1 < 3,-0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*
      ZP(gI1,2)*ZP(gI2,2),0) + IF(gO1 < 3,0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,3),0) + 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*
      KroneckerDelta(gO2,3 + j1))*ZP(gI1,0)*ZP(gI2,0) - SUM(j3,0,2,KroneckerDelta(
      gO1,3 + j3)*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(j1,0,2,Conj(Ye(j3,j1))
      *Ye(j2,j1))))*ZP(gI1,0)*ZP(gI2,0) - 0.3*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*ZP(gI1,1)*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiFeconjUSePR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -1.0954451150103321*g1*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*ZER(gI2,j1))*ZN1(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiFeconjUSePL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -(Conj(ZN2(gI1,2))*SUM(j2,0,2,Conj(ZEL(
      gI2,j2))*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*Ye(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpSdUSeconjSdconjUSe(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)),0) + IF(gO1 < 3,-0.125*
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
   const std::complex<double> result = IF(gO1 < 3,IF(gO2 < 3,-(SUM(j1,0,2,Ye(j1
      ,gO1)*ZE(gI2,3 + j1))*SUM(j3,0,2,Conj(Ye(j3,gO2))*Conj(ZE(gI1,3 + j3)))),0),
      0) + IF(gO1 < 3,IF(gO2 < 3,-0.15*Conj(ZE(gI1,gO2))*Sqr(g1)*ZE(gI2,gO1),0),0)
      + IF(gO1 < 3,IF(gO2 < 3,-0.25*Conj(ZE(gI1,gO2))*Sqr(g2)*ZE(gI2,gO1),0),0) +
      IF(gO1 < 3,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,j1
      ))*ZE(gI2,j1)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1
      ,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,
      gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)),0) + IF(gO1 < 3
      ,-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZE(gI1,j2))*ZE(gI2,
      j2)),0) + IF(gO1 < 3,-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(
      ZE(gI1,j2))*ZE(gI2,j2)),0) + IF(gO1 < 3,0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)
      *SUM(j2,0,2,Conj(ZE(gI1,3 + j2))*ZE(gI2,3 + j2)),0) + IF(gO1 < 3,-(SUM(j1,0,
      2,KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1))*SUM(j4,0,2,SUM(j3,0,2,Conj(Ye(j3,j4
      ))*Conj(ZE(gI1,3 + j3)))*ZE(gI2,j4))),0) + IF(gO1 < 3,0.15*Sqr(g1)*SUM(j1,0,
      2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZE(gI2,gO1),0) + IF(gO1 <
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
      j2)) - 0.3*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 +
      j1))*SUM(j2,0,2,KroneckerDelta(gO1,3 + j2)*ZE(gI2,3 + j2)) - SUM(j2,0,2,Conj
      (ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*SUM(j4,0,2,
      SUM(j3,0,2,Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3))*ZE(gI2,j4));

   return result;
}

std::complex<double> CLASSNAME::CpUSeSuconjUSeconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)
      *Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) + IF(gO1 < 3,0.125*
      KroneckerDelta(gO1,gO2)*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)),0) +
      IF(gO1 < 3,-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 +
      j1))*ZU(gI2,3 + j1)),0) + IF(gO1 < 3,0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)*
      SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,0.125*KroneckerDelta
      (gO1,gO2)*Sqr(g2)*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,j2)),0) + IF(gO1 < 3,
      -0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,
      3 + j2)),0) - 0.05*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1))*SUM(j2,0,
      2,KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2)) + 0.2*Sqr(g1)*SUM(
      j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))*SUM(j2,0,2,KroneckerDelta(gO1,3
      + j2)*KroneckerDelta(gO2,3 + j2)) - 0.05*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(
      gO1,3 + j1)*KroneckerDelta(gO2,3 + j1))*SUM(j2,0,2,Conj(ZU(gI1,j2))*ZU(gI2,
      j2)) + 0.2*Sqr(g1)*SUM(j1,0,2,KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,
      3 + j1))*SUM(j2,0,2,Conj(ZU(gI1,3 + j2))*ZU(gI2,3 + j2));

   return result;
}

std::complex<double> CLASSNAME::CpRhSeconjUSe(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,MuD*SUM(j1,0,2,Conj(Ye(j1,gO2
      ))*Conj(ZE(gI1,3 + j1)))*ZHR(gI2,0),0) + IF(gO2 < 3,0.7071067811865475*LamSD
      *vS*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1)))*ZHR(gI2,0),0) + IF(gO2
      < 3,0.5*LamTD*vT*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1)))*ZHR(gI2,
      0),0);

   return result;
}

std::complex<double> CLASSNAME::CpAhSeconjUSe(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*Mu*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1)))*ZA(
      gI2,1),0) + IF(gO2 < 3,std::complex<double>(0.,0.3872983346207417)*g1*MDBS*
      Conj(ZE(gI1,gO2))*ZA(gI2,2),0) + IF(gO2 < 3,std::complex<double>(0.,
      -0.3872983346207417)*g1*Conj(MDBS)*Conj(ZE(gI1,gO2))*ZA(gI2,2),0) + IF(gO2 <
      3,std::complex<double>(0,0.5)*g2*MDWBT*Conj(ZE(gI1,gO2))*ZA(gI2,3),0) + IF(
      gO2 < 3,std::complex<double>(0,-0.5)*g2*Conj(MDWBT)*Conj(ZE(gI1,gO2))*ZA(gI2
      ,3),0) - std::complex<double>(0.,0.7071067811865475)*Conj(Mu)*SUM(j2,0,2,
      Conj(ZE(gI1,j2))*SUM(j1,0,2,KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZA(gI2,1)
      - std::complex<double>(0.,0.7745966692414834)*g1*MDBS*SUM(j1,0,2,Conj(ZE(
      gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZA(gI2,2) + std::complex<double>(0.
      ,0.7745966692414834)*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZA(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CphhSeconjUSe(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.15*vd*Conj(ZE(gI1,gO2))*
      Sqr(g1)*ZH(gI2,0),0) + IF(gO2 < 3,0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0
      ),0) + IF(gO2 < 3,-(vd*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,gO2
      ))*Ye(j1,j2)))*ZH(gI2,0)),0) + IF(gO2 < 3,0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*
      ZH(gI2,1),0) + IF(gO2 < 3,-0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1),0) +
      IF(gO2 < 3,0.7071067811865475*Mu*SUM(j1,0,2,Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 +
      j1)))*ZH(gI2,1),0) + IF(gO2 < 3,0.3872983346207417*g1*MDBS*Conj(ZE(gI1,gO2)
      )*ZH(gI2,2),0) + IF(gO2 < 3,0.3872983346207417*g1*Conj(MDBS)*Conj(ZE(gI1,gO2
      ))*ZH(gI2,2),0) + IF(gO2 < 3,0.5*g2*MDWBT*Conj(ZE(gI1,gO2))*ZH(gI2,3),0) +
      IF(gO2 < 3,0.5*g2*Conj(MDWBT)*Conj(ZE(gI1,gO2))*ZH(gI2,3),0) + 0.3*vd*Sqr(g1
      )*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,0) - vd
      *SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,KroneckerDelta(gO2,3 + j2)*SUM(
      j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))))*ZH(gI2,0) - 0.3*vu*Sqr(g1)*SUM(j1,0,2,
      Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,1) +
      0.7071067811865475*Conj(Mu)*SUM(j2,0,2,Conj(ZE(gI1,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)))*ZH(gI2,1) - 0.7745966692414834*g1*
      MDBS*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1))*ZH(gI2,2) -
      0.7745966692414834*g1*Conj(MDBS)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*
      KroneckerDelta(gO2,3 + j1))*ZH(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSRdpconjUSe(int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) -
      vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZV(gI2,j2))*SUM(j1,0,2,
      KroneckerDelta(gO2,3 + j1)*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjUSeVWm(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,0.7071067811865475*g2*Conj(ZV
      (gI2,gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSeVP(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.3872983346207417*g1*Conj(
      ZE(gI2,gO2))*Cos(ThetaW()),0) + IF(gO2 < 3,-0.5*g2*Conj(ZE(gI2,gO2))*Sin(
      ThetaW()),0) - 0.7745966692414834*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUSeVZ(int gI2, int gO2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.5*g2*Conj(ZE(gI2,gO2))*Cos
      (ThetaW()),0) + IF(gO2 < 3,0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(
      ThetaW()),0) + 0.7745966692414834*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3
      + j1))*KroneckerDelta(gO2,3 + j1));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUhhconjSRdp(int gO2) const
{
   const std::complex<double> result = 0.05*(-7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2) - 20*vS*AbsSqr(LamSD)*KroneckerDelta(2,gO2) -
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(2,gO2) +
      7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(2,gO2) +
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2) -
      7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2) - 14.142135623730951*
      LamSD*Conj(MuD)*KroneckerDelta(2,gO2) - 10*g2*MDWBT*KroneckerDelta(3,gO2) -
      10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2) + 7.0710678118654755*LamTD*vS*Conj
      (LamSD)*KroneckerDelta(3,gO2) + 10*MuD*Conj(LamTD)*KroneckerDelta(3,gO2) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2) - 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2) + 10*LamTD*Conj(MuD)*KroneckerDelta(3,gO2) + vd
      *KroneckerDelta(0,gO2)*(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(1,gO2)*(-3*vu*Sqr(g1) + 5*vu*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUhhconjSRum(int gO2) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2) - 20*vS*AbsSqr(LamSU)*KroneckerDelta(2,gO2) -
      14.142135623730951*MuU*Conj(LamSU)*KroneckerDelta(2,gO2) -
      7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(2,gO2) -
      7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(2,gO2) +
      7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2) - 14.142135623730951*
      LamSU*Conj(MuU)*KroneckerDelta(2,gO2) + 10*g2*MDWBT*KroneckerDelta(3,gO2) -
      10*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2) - 7.0710678118654755*LamTU*vS*Conj
      (LamSU)*KroneckerDelta(3,gO2) - 10*MuU*Conj(LamTU)*KroneckerDelta(3,gO2) -
      7.0710678118654755*LamSU*vS*Conj(LamTU)*KroneckerDelta(3,gO2) + 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2) - 10*LamTU*Conj(MuU)*KroneckerDelta(3,gO2) + vu
      *KroneckerDelta(1,gO2)*(-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(0,gO2)*(-3*vd*Sqr(g1) + 5*vd*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpbargWmgWmUhh(int gO1) const
{
   const std::complex<double> result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1) + 4*vT*KroneckerDelta(3,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbargWmCgWmCUhh(int gO1) const
{
   const std::complex<double> result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*
      KroneckerDelta(1,gO1) + 4*vT*KroneckerDelta(3,gO1))*Sqr(g2);

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
      KroneckerDelta(1,gO2) + 4*vT*KroneckerDelta(3,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUhhUhhconjSRdp(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(5*(Conj(LamTD)*(1.4142135623730951
      *LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*
      LamSD*KroneckerDelta(2,gO1) - 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(
      3,gO2)) + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*LamSD*KroneckerDelta(2,gO2
      ) + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1
      )*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUhhUhhconjSRum(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-5*(Conj(LamTU)*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + 2*LamTU*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2)) + Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSU
      *KroneckerDelta(2,gO2) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)))) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-20*AbsSqr(LamTU) + 3*Sqr(g1)
      - 5*Sqr(g2)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-3*Sqr(g1) + 5*
      Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(int gO1, int gO2) const
{
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhRhconjRh(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + 2*LamTD*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))*ZHR(gI1,0)*ZHR(gI2,0) + Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTD*KroneckerDelta(3,gO2)))*ZHR(gI1,0)*ZHR(gI2,0) + (-(Conj(LamTU)*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*LamTU*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*
      LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)
      )))*ZHR(gI1,1)*ZHR(gI2,1)) - KroneckerDelta(1,gO1)*(5*KroneckerDelta(0,gO2)*
      (-2*LamSU*Conj(LamSD)*ZHR(gI1,1)*ZHR(gI2,0) + LamTU*Conj(LamTD)*ZHR(gI1,1)*
      ZHR(gI2,0) + (-2*LamSD*Conj(LamSU) + LamTD*Conj(LamTU))*ZHR(gI1,0)*ZHR(gI2,1
      )) + KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) +
      (20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,1)*ZHR
      (gI2,1))) + KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*(2*LamSU*Conj(
      LamSD)*ZHR(gI1,1)*ZHR(gI2,0) - LamTU*Conj(LamTD)*ZHR(gI1,1)*ZHR(gI2,0) + (2*
      LamSD*Conj(LamSU) - LamTD*Conj(LamTU))*ZHR(gI1,0)*ZHR(gI2,1)) +
      KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5
      *Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2
      ,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhRhconjRh(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.05*(-7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 20*vS*AbsSqr(LamSD)*
      KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 14.142135623730951*MuD*Conj(
      LamSD)*KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 7.0710678118654755*
      LamTD*vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(gI2,0) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR
      (gI2,0) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZHR(gI1,0)*
      ZHR(gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)*ZHR(
      gI1,0)*ZHR(gI2,0) + 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR(gI2,0)
      - 10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR(gI2,0) -
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR
      (gI2,0) - 10*MuD*Conj(LamTD)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR(gI2,0) -
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR
      (gI2,0) + 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 10
      *LamTD*Conj(MuD)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR(gI2,0) +
      7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZHR(gI1,1)*ZHR(gI2,1) - 20*
      vS*AbsSqr(LamSU)*KroneckerDelta(2,gO2)*ZHR(gI1,1)*ZHR(gI2,1) -
      14.142135623730951*MuU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZHR(gI1,1)*ZHR(gI2,
      1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(2,gO2)*ZHR(gI1,1
      )*ZHR(gI2,1) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(2,gO2)
      *ZHR(gI1,1)*ZHR(gI2,1) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,
      gO2)*ZHR(gI1,1)*ZHR(gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*
      KroneckerDelta(2,gO2)*ZHR(gI1,1)*ZHR(gI2,1) - 10*g2*MDWBT*KroneckerDelta(3,
      gO2)*ZHR(gI1,1)*ZHR(gI2,1) - 10*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZHR(
      gI1,1)*ZHR(gI2,1) + 7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(3
      ,gO2)*ZHR(gI1,1)*ZHR(gI2,1) + 10*MuU*Conj(LamTU)*KroneckerDelta(3,gO2)*ZHR(
      gI1,1)*ZHR(gI2,1) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*KroneckerDelta(3
      ,gO2)*ZHR(gI1,1)*ZHR(gI2,1) - 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZHR(
      gI1,1)*ZHR(gI2,1) + 10*LamTU*Conj(MuU)*KroneckerDelta(3,gO2)*ZHR(gI1,1)*ZHR(
      gI2,1) + KroneckerDelta(0,gO2)*(ZHR(gI1,0)*(vd*(-20*AbsSqr(LamSD) - 10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI2,0) + 5*vu*(2*LamSU*Conj(LamSD
      ) - LamTU*Conj(LamTD))*ZHR(gI2,1)) + ZHR(gI1,1)*(10*LamSD*vu*Conj(LamSU)*ZHR
      (gI2,0) - 5*LamTD*vu*Conj(LamTU)*ZHR(gI2,0) - vd*(3*Sqr(g1) + 5*Sqr(g2))*ZHR
      (gI2,1))) - KroneckerDelta(1,gO2)*(ZHR(gI1,0)*(vu*(3*Sqr(g1) + 5*Sqr(g2))*
      ZHR(gI2,0) + 5*vd*(-2*LamSU*Conj(LamSD) + LamTU*Conj(LamTD))*ZHR(gI2,1)) +
      ZHR(gI1,1)*(-(vu*(3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI2,1)) - 10*Conj(LamSU)*(LamSD
      *vd*ZHR(gI2,0) - 2*LamSU*vu*ZHR(gI2,1)) + 5*Conj(LamTU)*(LamTD*vd*ZHR(gI2,0)
      + 2*LamTU*vu*ZHR(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1UhhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*(KroneckerDelta(3,gO2)*(-2*g2*UM1(
      gI2,0)*UP1(gI1,0) + Conj(LamTD)*UM1(gI2,1)*UP1(gI1,1)) - 1.4142135623730951*
      (Conj(LamSD)*KroneckerDelta(2,gO2)*UM1(gI2,1)*UP1(gI1,1) + KroneckerDelta(0,
      gO2)*(g2*UM1(gI2,1)*UP1(gI1,0) + Conj(LamTD)*UM1(gI2,0)*UP1(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1UhhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = 0.5*(-(Conj(UM1(gI1,0))*(
      1.4142135623730951*LamTD*Conj(UP1(gI2,1))*KroneckerDelta(0,gO1) + 2*g2*Conj(
      UP1(gI2,0))*KroneckerDelta(3,gO1))) - Conj(UM1(gI1,1))*(1.4142135623730951*
      g2*Conj(UP1(gI2,0))*KroneckerDelta(0,gO1) + Conj(UP1(gI2,1))*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1) - LamTD*KroneckerDelta(3,gO1)
      )));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2UhhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*(KroneckerDelta(3,gO2)*(2*g2*UM2(gI1
      ,0)*UP2(gI2,0) + Conj(LamTU)*UM2(gI1,1)*UP2(gI2,1)) + 1.4142135623730951*(
      Conj(LamTU)*KroneckerDelta(1,gO2)*UM2(gI1,1)*UP2(gI2,0) + (-(g2*
      KroneckerDelta(1,gO2)*UM2(gI1,0)) + Conj(LamSU)*KroneckerDelta(2,gO2)*UM2(
      gI1,1))*UP2(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2UhhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = 0.5*(g2*Conj(UM2(gI2,0))*(
      -1.4142135623730951*Conj(UP2(gI1,1))*KroneckerDelta(1,gO1) + 2*Conj(UP2(gI1,
      0))*KroneckerDelta(3,gO1)) + Conj(UM2(gI2,1))*(1.4142135623730951*LamTU*Conj
      (UP2(gI1,0))*KroneckerDelta(1,gO1) + Conj(UP2(gI1,1))*(1.4142135623730951*
      LamSU*KroneckerDelta(2,gO1) + LamTU*KroneckerDelta(3,gO1))));

   return result;
}

std::complex<double> CLASSNAME::CpAhUhhconjRh(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,-0.25)*Mu*(2*Conj
      (LamSD)*(KroneckerDelta(2,gO2)*ZA(gI2,1) - KroneckerDelta(1,gO2)*ZA(gI2,2))*
      ZHR(gI1,0) + 1.4142135623730951*Conj(LamTD)*(KroneckerDelta(3,gO2)*ZA(gI2,1)
      - KroneckerDelta(1,gO2)*ZA(gI2,3))*ZHR(gI1,0) + (Conj(LamSU)*(-2*
      KroneckerDelta(2,gO2)*ZA(gI2,0) + 2*KroneckerDelta(0,gO2)*ZA(gI2,2)) +
      1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZA(gI2,0) -
      KroneckerDelta(0,gO2)*ZA(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CphhUhhconjRh(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = -0.25*Mu*(2*Conj(LamSD)*(KroneckerDelta(
      2,gO2)*ZH(gI2,1) + KroneckerDelta(1,gO2)*ZH(gI2,2))*ZHR(gI1,0) +
      1.4142135623730951*Conj(LamTD)*(KroneckerDelta(3,gO2)*ZH(gI2,1) +
      KroneckerDelta(1,gO2)*ZH(gI2,3))*ZHR(gI1,0) + (-2*Conj(LamSU)*(
      KroneckerDelta(2,gO2)*ZH(gI2,0) + KroneckerDelta(0,gO2)*ZH(gI2,2)) +
      1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZH(gI2,0) +
      KroneckerDelta(0,gO2)*ZH(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSvconjSv(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.05*KroneckerDelta(gI1,gI2)*(2*(
      3.872983346207417*g1*(MDBS + Conj(MDBS))*KroneckerDelta(2,gO2) - 5*g2*(MDWBT
      + Conj(MDWBT))*KroneckerDelta(3,gO2)) - vd*KroneckerDelta(0,gO2)*(3*Sqr(g1)
      + 5*Sqr(g2)) + vu*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO2
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZDL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUhhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO1
      )*SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO2
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZEL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUhhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(0,gO1
      )*SUM(j2,0,2,Conj(ZEL(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO2
      )*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUhhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -0.7071067811865475*KroneckerDelta(1,gO1
      )*SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + 2*LamTD*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))*ZA(gI1,0)*ZA(gI2,0) + Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTD*KroneckerDelta(3,gO2)))*ZA(gI1,0)*ZA(gI2,0) + (-(Conj(LamTU)*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*LamTU*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*
      LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)
      )))*ZA(gI1,1)*ZA(gI2,1)) - KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*
      Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)
      *ZA(gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) +
      ZA(gI1,2)*(4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))) + Conj(
      LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(
      1.4142135623730951*LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3))))) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) -
      (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTU)*(
      1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSU*ZA(gI2,2) - 2*LamTU*ZA(gI2,3))) + Conj(LamSU)*(1.4142135623730951*
      LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSU*ZA(gI2,2) +
      1.4142135623730951*LamTU*ZA(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhhUhh(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-5*(4*AbsSqr(LamSU)*KroneckerDelta
      (2,gO1)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU
      *Conj(LamSU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZH(gI1,1)*ZH(gI2,1)
      - 1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*KroneckerDelta
      (3,gO1)*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO1)*KroneckerDelta(3,
      gO2)*ZH(gI1,1)*ZH(gI2,1) + 2*AbsSqr(LamTU)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1) + 4*AbsSqr(LamSU)*KroneckerDelta(1
      ,gO2)*KroneckerDelta(2,gO1)*ZH(gI1,2)*ZH(gI2,1) - 1.4142135623730951*LamTU*
      Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*ZH(gI1,2)*ZH(gI2,1)
      - 1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      3,gO1)*ZH(gI1,2)*ZH(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZH(gI1,3)*ZH(gI2,1) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,
      gO1)*ZH(gI1,3)*ZH(gI2,1) + 2*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZH(gI1,3)*ZH(gI2,1) + 4*AbsSqr(LamSU)*KroneckerDelta(1
      ,gO2)*KroneckerDelta(2,gO1)*ZH(gI1,1)*ZH(gI2,2) - 1.4142135623730951*LamTU*
      Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*ZH(gI1,1)*ZH(gI2,2)
      - 1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      3,gO1)*ZH(gI1,1)*ZH(gI2,2) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZH(gI1,1)*ZH(gI2,3) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,
      gO1)*ZH(gI1,1)*ZH(gI2,3) + 2*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZH(gI1,1)*ZH(gI2,3) + Conj(LamSD)*(1.4142135623730951*
      LamTD*KroneckerDelta(3,gO1)*(KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) +
      KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2))) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) +
      1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) +
      KroneckerDelta(0,gO2)*(4*LamSD*ZH(gI1,2)*ZH(gI2,0) + 1.4142135623730951*
      LamTD*ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*(4*LamSD*ZH(gI2,2) +
      1.4142135623730951*LamTD*ZH(gI2,3))))) + Conj(LamTD)*(1.4142135623730951*
      LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZH(gI1,0)*ZH(gI2,0) +
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1)*(KroneckerDelta(3,gO2)*ZH(gI1
      ,0)*ZH(gI2,0) + KroneckerDelta(0,gO2)*(ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*ZH(
      gI2,3))) + KroneckerDelta(3,gO1)*(2*LamTD*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH
      (gI2,0) + KroneckerDelta(0,gO2)*(1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,0
      ) + 2*LamTD*ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*(1.4142135623730951*LamSD*ZH(gI2
      ,2) + 2*LamTD*ZH(gI2,3)))))) - KroneckerDelta(0,gO1)*(-(KroneckerDelta(1,gO2
      )*(3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1))) + 5*(
      Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*(ZH(gI1,2)*ZH(
      gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + KroneckerDelta(2,gO2)*(4*LamSD*ZH(gI1,2)*ZH(
      gI2,0) + 1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*(4*LamSD*
      ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3)))) + Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*(ZH(gI1,3)*ZH(gI2,0) + ZH(gI1
      ,0)*ZH(gI2,3)) + KroneckerDelta(3,gO2)*(1.4142135623730951*LamSD*ZH(gI1,2)*
      ZH(gI2,0) + 2*LamTD*ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*(1.4142135623730951*
      LamSD*ZH(gI2,2) + 2*LamTD*ZH(gI2,3))))) + KroneckerDelta(0,gO2)*(3*(3*Sqr(g1
      ) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(
      gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(
      gI1,2)*(4*LamSD*ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3))) + Conj(
      LamTD)*(1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(
      1.4142135623730951*LamSD*ZH(gI2,2) + 2*LamTD*ZH(gI2,3)))))) + KroneckerDelta
      (1,gO1)*(KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0)
      + ZH(gI1,0)*ZH(gI2,1)) + 5*(Conj(LamTU)*(1.4142135623730951*LamSU*
      KroneckerDelta(2,gO2)*(ZH(gI1,3)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,3)) +
      KroneckerDelta(3,gO2)*(1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,1) - 2*
      LamTU*ZH(gI1,3)*ZH(gI2,1) + ZH(gI1,1)*(1.4142135623730951*LamSU*ZH(gI2,2) -
      2*LamTU*ZH(gI2,3)))) + Conj(LamSU)*(1.4142135623730951*LamTU*KroneckerDelta(
      3,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) + KroneckerDelta(2,gO2)*(
      -4*LamSU*ZH(gI1,2)*ZH(gI2,1) + 1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,1)
      + ZH(gI1,1)*(-4*LamSU*ZH(gI2,2) + 1.4142135623730951*LamTU*ZH(gI2,3))))) +
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - 3*(3*
      Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 5*(Conj(LamTU)*(
      1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951
      *LamSU*ZH(gI2,2) - 2*LamTU*ZH(gI2,3))) + Conj(LamSU)*(1.4142135623730951*
      LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSU*ZH(gI2,2) +
      1.4142135623730951*LamTU*ZH(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*(5*(
      -1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,2)*ZP(gI2,0)
      + 1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 2*
      LamSD*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,0) +
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,3)*ZP(gI2,0) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      2*LamSD*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,3) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) + 2*
      LamTD*Conj(LamSD)*KroneckerDelta(2,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,0)*ZP(
      gI2,3))) + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0
      ) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*AbsSqr(LamTD) +
      Sqr(g2))*ZP(gI1,2)*ZP(gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3))))) +
      KroneckerDelta(1,gO1)*(-5*(-1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(
      3,gO2)*ZP(gI1,2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2
      )*ZP(gI1,2)*ZP(gI2,1) + 2*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*ZP(gI1,3)*
      ZP(gI2,1) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,3)
      *ZP(gI2,1) - 1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(
      gI2,1) + KroneckerDelta(0,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(
      gI2,1)) + 2*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,2) -
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,3) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) + 2*
      LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*(ZP(gI1,2)*ZP(gI2,1) + ZP(gI1,1)*ZP(
      gI2,3))) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0
      ) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(
      gI2,2) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3)))) - 5*(
      1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,2)*ZP(gI2,0) - 1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3
      ,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) + 4*AbsSqr(LamSU)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamTU*Conj(
      LamSU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,
      gO1)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO1)*KroneckerDelta(3,
      gO2)*ZP(gI1,1)*ZP(gI2,1) + 2*AbsSqr(LamTU)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) + 2*LamTU*Conj(LamSU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,2)*ZP(gI2,1) -
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)
      *ZP(gI1,2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) + 2*LamSU*Conj(LamTU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)
      *ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) + 1.4142135623730951*
      KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) + 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,1)*ZP(
      gI2,2) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,2) + 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) -
      1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,3) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      2,gO1)*ZP(gI1,1)*ZP(gI2,3) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta
      (1,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,3) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) - 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3) +
      Conj(LamTD)*(-1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(
      3,gO1)*ZP(gI1,0)*ZP(gI2,0) + LamSD*KroneckerDelta(2,gO1)*(
      -1.4142135623730951*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 2*
      KroneckerDelta(0,gO2)*(ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,2))) + LamTD*
      KroneckerDelta(3,gO1)*(2*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) +
      1.4142135623730951*KroneckerDelta(0,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,3)
      *ZP(gI2,0) + ZP(gI1,0)*(-ZP(gI2,2) + ZP(gI2,3))))) + Conj(LamSD)*(
      -1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZP(gI1
      ,0)*ZP(gI2,0) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)*ZP(gI1,
      0)*ZP(gI2,0) + LamTD*(-1.4142135623730951*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP
      (gI2,0) + 2*KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,3)
      ))))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUhh(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) - 20*vS*AbsSqr(LamSD)*
      KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) - 14.142135623730951*MuD*Conj(
      LamSD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) - 7.0710678118654755*LamTD*
      vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(
      gI2,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(
      gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*
      ZA(gI2,0) - 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 10*vT*
      AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 7.0710678118654755
      *LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 10*MuD*
      Conj(LamTD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 7.0710678118654755*
      LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 10*LamTD*Conj(MuD)*
      KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) - 20*vS*AbsSqr(LamSU)*
      KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) - 14.142135623730951*MuU*Conj(
      LamSU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) + 7.0710678118654755*LamTU*
      vT*Conj(LamSU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) +
      7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(
      gI2,1) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(
      gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*
      ZA(gI2,1) + 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) - 10*vT*
      AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 7.0710678118654755
      *LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 10*MuU*
      Conj(LamTU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 7.0710678118654755*
      LamSU*vS*Conj(LamTU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 10*LamTU*Conj(MuU)*
      KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) - vd*KroneckerDelta(0,gO2)*((3*Sqr
      (g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA
      (gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(
      gI1,2)*(4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))) + Conj(
      LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(
      1.4142135623730951*LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3))))) + vu*
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(
      g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTU)*(1.4142135623730951*
      LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSU*ZA(gI2,2) -
      2*LamTU*ZA(gI2,3))) + Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2
      ,2) + ZA(gI1,2)*(-4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3)))))
      );

   return result;
}

std::complex<double> CLASSNAME::CpAhhhUhh(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      -7.0710678118654755*LamSD*vd*Conj(LamTD)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZH(
      gI1,0) + 7.0710678118654755*LamSD*vd*Conj(LamTD)*KroneckerDelta(2,gO2)*ZA(
      gI2,3)*ZH(gI1,0) + 7.0710678118654755*LamTD*vd*Conj(LamSD)*(KroneckerDelta(3
      ,gO2)*ZA(gI2,2) - KroneckerDelta(2,gO2)*ZA(gI2,3))*ZH(gI1,0) -
      7.745966692414834*g1*MDBS*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,1) +
      14.142135623730951*MuU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,1)
      - 7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*
      ZH(gI1,1) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA
      (gI2,2)*ZH(gI1,1) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(1,gO2)*ZA
      (gI2,2)*ZH(gI1,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(1,gO2)
      *ZA(gI2,2)*ZH(gI1,1) - 7.0710678118654755*LamTU*vu*Conj(LamSU)*
      KroneckerDelta(3,gO2)*ZA(gI2,2)*ZH(gI1,1) + 7.0710678118654755*LamSU*vu*Conj
      (LamTU)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZH(gI1,1) + 10*g2*MDWBT*
      KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) + 7.0710678118654755*LamTU*vS*Conj
      (LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) - 10*MuU*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) - 7.0710678118654755*LamSU*vS*Conj
      (LamTU)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) - 10*g2*Conj(MDWBT)*
      KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) + 10*LamTU*Conj(MuU)*
      KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) + 7.0710678118654755*LamTU*vu*Conj
      (LamSU)*KroneckerDelta(2,gO2)*ZA(gI2,3)*ZH(gI1,1) - 7.0710678118654755*LamSU
      *vu*Conj(LamTU)*KroneckerDelta(2,gO2)*ZA(gI2,3)*ZH(gI1,1) +
      7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(
      gI1,2) - 7.0710678118654755*LamSU*vu*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(
      gI2,3)*ZH(gI1,2) - 7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(1,
      gO2)*ZA(gI2,2)*ZH(gI1,3) + 7.0710678118654755*LamSU*vu*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,3) + KroneckerDelta(0,gO2)*(-5*ZA(gI2
      ,3)*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) - 2*MuD*Conj(
      LamTD) - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) + 2*
      LamTD*Conj(MuD))*ZH(gI1,0) + 1.4142135623730951*vd*(LamTD*Conj(LamSD) -
      LamSD*Conj(LamTD))*ZH(gI1,2)) + ZA(gI2,2)*((7.745966692414834*g1*MDBS +
      7.0710678118654755*(2*MuD + LamTD*vT)*Conj(LamSD) - 7.0710678118654755*LamSD
      *vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSD
      *Conj(MuD))*ZH(gI1,0) + 7.0710678118654755*vd*(LamTD*Conj(LamSD) - LamSD*
      Conj(LamTD))*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUhh(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) - 20*vS*AbsSqr(LamSD)*
      KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) - 14.142135623730951*MuD*Conj(
      LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) - 7.0710678118654755*LamTD*
      vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(
      gI2,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(
      gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*
      ZH(gI2,0) - 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 10*vT*
      AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 7.0710678118654755
      *LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 10*MuD*
      Conj(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 7.0710678118654755*
      LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 10*LamTD*Conj(MuD)*
      KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 20*vd*AbsSqr(LamSD)*
      KroneckerDelta(2,gO2)*ZH(gI1,2)*ZH(gI2,0) - 7.0710678118654755*LamTD*vd*Conj
      (LamSD)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(gI2,0) - 7.0710678118654755*LamSD
      *vd*Conj(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(gI2,0) -
      7.0710678118654755*LamTD*vd*Conj(LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,3)*ZH(
      gI2,0) - 7.0710678118654755*LamSD*vd*Conj(LamTD)*KroneckerDelta(2,gO2)*ZH(
      gI1,3)*ZH(gI2,0) - 10*vd*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,3)*ZH(
      gI2,0) - 7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,1)
      - 20*vS*AbsSqr(LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,1) -
      14.142135623730951*MuU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,1)
      + 7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*
      ZH(gI2,1) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(2,gO2)*ZH
      (gI1,1)*ZH(gI2,1) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZH
      (gI1,1)*ZH(gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(2,gO2)
      *ZH(gI1,1)*ZH(gI2,1) + 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1)
      - 10*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1) +
      7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(
      gI2,1) + 10*MuU*Conj(LamTU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1) +
      7.0710678118654755*LamSU*vS*Conj(LamTU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(
      gI2,1) + 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1) + 10*
      LamTU*Conj(MuU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,1) - 20*vu*AbsSqr(
      LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,2)*ZH(gI2,1) + 7.0710678118654755*LamTU*
      vu*Conj(LamSU)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(gI2,1) +
      7.0710678118654755*LamSU*vu*Conj(LamTU)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(
      gI2,1) + 7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(2,gO2)*ZH(
      gI1,3)*ZH(gI2,1) + 7.0710678118654755*LamSU*vu*Conj(LamTU)*KroneckerDelta(2,
      gO2)*ZH(gI1,3)*ZH(gI2,1) - 10*vu*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZH(gI1,
      3)*ZH(gI2,1) - 20*vd*AbsSqr(LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,2)
      - 7.0710678118654755*LamTD*vd*Conj(LamSD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*
      ZH(gI2,2) - 7.0710678118654755*LamSD*vd*Conj(LamTD)*KroneckerDelta(3,gO2)*ZH
      (gI1,0)*ZH(gI2,2) - 20*vu*AbsSqr(LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(
      gI2,2) + 7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(3,gO2)*ZH(
      gI1,1)*ZH(gI2,2) + 7.0710678118654755*LamSU*vu*Conj(LamTU)*KroneckerDelta(3,
      gO2)*ZH(gI1,1)*ZH(gI2,2) - 7.0710678118654755*LamTD*vd*Conj(LamSD)*
      KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,3) - 7.0710678118654755*LamSD*vd*Conj
      (LamTD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,3) - 10*vd*AbsSqr(LamTD)*
      KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,3) + 7.0710678118654755*LamTU*vu*Conj
      (LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,3) + 7.0710678118654755*LamSU
      *vu*Conj(LamTU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,3) - 10*vu*AbsSqr(
      LamTU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,3) - KroneckerDelta(0,gO2)*(
      -7.745966692414834*g1*MDBS*ZH(gI1,2)*ZH(gI2,0) + 20*vS*AbsSqr(LamSD)*ZH(gI1,
      2)*ZH(gI2,0) + 14.142135623730951*MuD*Conj(LamSD)*ZH(gI1,2)*ZH(gI2,0) +
      7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gI1,2)*ZH(gI2,0) +
      7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gI1,2)*ZH(gI2,0) -
      7.745966692414834*g1*Conj(MDBS)*ZH(gI1,2)*ZH(gI2,0) + 14.142135623730951*
      LamSD*Conj(MuD)*ZH(gI1,2)*ZH(gI2,0) + 10*g2*MDWBT*ZH(gI1,3)*ZH(gI2,0) + 10*
      vT*AbsSqr(LamTD)*ZH(gI1,3)*ZH(gI2,0) + 7.0710678118654755*LamTD*vS*Conj(
      LamSD)*ZH(gI1,3)*ZH(gI2,0) + 10*MuD*Conj(LamTD)*ZH(gI1,3)*ZH(gI2,0) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gI1,3)*ZH(gI2,0) + 10*g2*Conj(
      MDWBT)*ZH(gI1,3)*ZH(gI2,0) + 10*LamTD*Conj(MuD)*ZH(gI1,3)*ZH(gI2,0) - (3*Sqr
      (g1) + 5*Sqr(g2))*ZH(gI1,1)*(vu*ZH(gI2,0) + vd*ZH(gI2,1)) + 20*vd*AbsSqr(
      LamSD)*ZH(gI1,2)*ZH(gI2,2) + 7.0710678118654755*LamTD*vd*Conj(LamSD)*ZH(gI1,
      3)*ZH(gI2,2) + 7.0710678118654755*LamSD*vd*Conj(LamTD)*ZH(gI1,3)*ZH(gI2,2) +
      7.0710678118654755*LamTD*vd*Conj(LamSD)*ZH(gI1,2)*ZH(gI2,3) +
      7.0710678118654755*LamSD*vd*Conj(LamTD)*ZH(gI1,2)*ZH(gI2,3) + 10*vd*AbsSqr(
      LamTD)*ZH(gI1,3)*ZH(gI2,3) + ZH(gI1,0)*(3*vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,
      0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) - 7.745966692414834*g1*MDBS*ZH(gI2
      ,2) + 20*vS*AbsSqr(LamSD)*ZH(gI2,2) + 14.142135623730951*MuD*Conj(LamSD)*ZH(
      gI2,2) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gI2,2) +
      7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gI2,2) - 7.745966692414834*g1*
      Conj(MDBS)*ZH(gI2,2) + 14.142135623730951*LamSD*Conj(MuD)*ZH(gI2,2) + 10*g2*
      MDWBT*ZH(gI2,3) + 10*vT*AbsSqr(LamTD)*ZH(gI2,3) + 7.0710678118654755*LamTD*
      vS*Conj(LamSD)*ZH(gI2,3) + 10*MuD*Conj(LamTD)*ZH(gI2,3) + 7.0710678118654755
      *LamSD*vS*Conj(LamTD)*ZH(gI2,3) + 10*g2*Conj(MDWBT)*ZH(gI2,3) + 10*LamTD*
      Conj(MuD)*ZH(gI2,3))) + KroneckerDelta(1,gO2)*(-7.745966692414834*g1*MDBS*ZH
      (gI1,2)*ZH(gI2,1) - 20*vS*AbsSqr(LamSU)*ZH(gI1,2)*ZH(gI2,1) -
      14.142135623730951*MuU*Conj(LamSU)*ZH(gI1,2)*ZH(gI2,1) + 7.0710678118654755*
      LamTU*vT*Conj(LamSU)*ZH(gI1,2)*ZH(gI2,1) + 7.0710678118654755*LamSU*vT*Conj(
      LamTU)*ZH(gI1,2)*ZH(gI2,1) - 7.745966692414834*g1*Conj(MDBS)*ZH(gI1,2)*ZH(
      gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gI1,2)*ZH(gI2,1) + 10*g2*
      MDWBT*ZH(gI1,3)*ZH(gI2,1) - 10*vT*AbsSqr(LamTU)*ZH(gI1,3)*ZH(gI2,1) +
      7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gI1,3)*ZH(gI2,1) + 10*MuU*Conj(
      LamTU)*ZH(gI1,3)*ZH(gI2,1) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gI1,
      3)*ZH(gI2,1) + 10*g2*Conj(MDWBT)*ZH(gI1,3)*ZH(gI2,1) + 10*LamTU*Conj(MuU)*ZH
      (gI1,3)*ZH(gI2,1) + (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*(vu*ZH(gI2,0) + vd*ZH(
      gI2,1)) - 20*vu*AbsSqr(LamSU)*ZH(gI1,2)*ZH(gI2,2) + 7.0710678118654755*LamTU
      *vu*Conj(LamSU)*ZH(gI1,3)*ZH(gI2,2) + 7.0710678118654755*LamSU*vu*Conj(LamTU
      )*ZH(gI1,3)*ZH(gI2,2) + 7.0710678118654755*LamTU*vu*Conj(LamSU)*ZH(gI1,2)*ZH
      (gI2,3) + 7.0710678118654755*LamSU*vu*Conj(LamTU)*ZH(gI1,2)*ZH(gI2,3) - 10*
      vu*AbsSqr(LamTU)*ZH(gI1,3)*ZH(gI2,3) + ZH(gI1,1)*(vd*(3*Sqr(g1) + 5*Sqr(g2))
      *ZH(gI2,0) - 3*vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) - 7.745966692414834*g1*
      MDBS*ZH(gI2,2) - 20*vS*AbsSqr(LamSU)*ZH(gI2,2) - 14.142135623730951*MuU*Conj
      (LamSU)*ZH(gI2,2) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gI2,2) +
      7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gI2,2) - 7.745966692414834*g1*
      Conj(MDBS)*ZH(gI2,2) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gI2,2) + 10*g2*
      MDWBT*ZH(gI2,3) - 10*vT*AbsSqr(LamTU)*ZH(gI2,3) + 7.0710678118654755*LamTU*
      vS*Conj(LamSU)*ZH(gI2,3) + 10*MuU*Conj(LamTU)*ZH(gI2,3) + 7.0710678118654755
      *LamSU*vS*Conj(LamTU)*ZH(gI2,3) + 10*g2*Conj(MDWBT)*ZH(gI2,3) + 10*LamTU*
      Conj(MuU)*ZH(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*MDBS*
      KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0) - 20*vS*AbsSqr(LamSD)*
      KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0) - 14.142135623730951*MuD*Conj(
      LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0) + 7.0710678118654755*LamTD*
      vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0) +
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(
      gI2,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(
      gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*
      ZP(gI2,0) + 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) - 10*vT*
      AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 7.0710678118654755
      *LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 10*MuD*
      Conj(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 7.0710678118654755*
      LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 10*LamTD*Conj(MuD)*
      KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) - 10*LamSD*vd*Conj(LamTD)*
      KroneckerDelta(2,gO2)*ZP(gI1,2)*ZP(gI2,0) + 7.0710678118654755*vd*AbsSqr(
      LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,2)*ZP(gI2,0) - 7.0710678118654755*vd*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) - 10*LamTD*vd*Conj(LamSD)*
      KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,0) - 7.0710678118654755*vd*AbsSqr(
      LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,3)*ZP(gI2,0) + 7.0710678118654755*vd*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) - 7.745966692414834*g1*
      MDBS*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) - 20*vS*AbsSqr(LamSU)*
      KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) - 14.142135623730951*MuU*Conj(
      LamSU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) - 7.0710678118654755*LamTU*
      vT*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) -
      7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(
      gI2,1) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(
      gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*
      ZP(gI2,1) - 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 10*vT*
      AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 7.0710678118654755
      *LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 10*MuU*
      Conj(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 7.0710678118654755*
      LamSU*vS*Conj(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 10*g2*Conj(
      MDWBT)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 10*LamTU*Conj(MuU)*
      KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 10*LamSU*vu*Conj(LamTU)*
      KroneckerDelta(2,gO2)*ZP(gI1,2)*ZP(gI2,1) + 7.0710678118654755*vu*AbsSqr(
      LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,2)*ZP(gI2,1) - 7.0710678118654755*vu*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) - 10*LamTU*vu*Conj(LamSU)*
      KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,1) - 7.0710678118654755*vu*AbsSqr(
      LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,3)*ZP(gI2,1) + 7.0710678118654755*vu*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) - 10*LamTD*vd*Conj(LamSD)*
      KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,2) + 7.0710678118654755*vd*AbsSqr(
      LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,2) - 7.0710678118654755*vd*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) - 10*LamTU*vu*Conj(LamSU)*
      KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,2) + 7.0710678118654755*vu*AbsSqr(
      LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,2) - 7.0710678118654755*vu*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) - 20*vT*KroneckerDelta(3,
      gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) + 20*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,2) - 10*LamSD*vd*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*
      ZP(gI2,3) - 7.0710678118654755*vd*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1
      ,0)*ZP(gI2,3) + 7.0710678118654755*vd*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0
      )*ZP(gI2,3) - 10*LamSU*vu*Conj(LamTU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2
      ,3) - 7.0710678118654755*vu*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP
      (gI2,3) + 7.0710678118654755*vu*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(
      gI2,3) + 20*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) - 20*vT*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3) - KroneckerDelta(0,gO2)*(
      ZP(gI1,1)*(5*vu*Sqr(g2)*ZP(gI2,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,1)) +
      ZP(gI1,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,0) + 5*(vu*Sqr(g2)*ZP(gI2,1) +
      (2*LamTD*vS*Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD)
      + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*ZP(gI2,2) + ((2.8284271247461903*MuD + 2
      *LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2
      *(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gI2,3))) + 5*(ZP(gI1,2)*(((
      2.8284271247461903*MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(
      LamTD) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gI2,0) - 2*vd*(-2
      *AbsSqr(LamTD) + Sqr(g2))*ZP(gI2,2)) + ZP(gI1,3)*((2*LamTD*vS*Conj(LamSD) +
      1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) - vT*
      Sqr(g2)))*ZP(gI2,0) + 2*vd*Sqr(g2)*ZP(gI2,3)))) - KroneckerDelta(1,gO2)*(ZP(
      gI1,0)*((-3*vu*Sqr(g1) + 5*vu*Sqr(g2))*ZP(gI2,0) + 5*vd*Sqr(g2)*ZP(gI2,1)) +
      ZP(gI1,1)*(5*vd*Sqr(g2)*ZP(gI2,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,1) +
      5*((2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(
      LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gI2,2) + ((2.8284271247461903*
      MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(LamTU) +
      1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gI2,3))) + 5*(ZP(gI1,2)
      *(((2.8284271247461903*MuU + 2*LamSU*vS - 1.4142135623730951*LamTU*vT)*Conj(
      LamTU) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gI2,1) + 2*vu*Sqr
      (g2)*ZP(gI2,2)) + ZP(gI1,3)*((2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2
      *g2*MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gI2,1) -
      2*vu*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiChiUhhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.1*(5*Conj(LamTD)*KroneckerDelta(0,gO2)
      *ZN1(gI1,2)*ZN2(gI2,1) + 5*Conj(LamTU)*KroneckerDelta(1,gO2)*ZN1(gI1,3)*ZN2(
      gI2,1) + 3.872983346207417*g1*KroneckerDelta(0,gO2)*ZN1(gI1,0)*ZN2(gI2,2) -
      5*g2*KroneckerDelta(0,gO2)*ZN1(gI1,1)*ZN2(gI2,2) + 5*Conj(LamTD)*
      KroneckerDelta(3,gO2)*ZN1(gI1,2)*ZN2(gI2,2) + 7.0710678118654755*Conj(LamSD)
      *ZN1(gI1,2)*(KroneckerDelta(0,gO2)*ZN2(gI2,0) + KroneckerDelta(2,gO2)*ZN2(
      gI2,2)) - 3.872983346207417*g1*KroneckerDelta(1,gO2)*ZN1(gI1,0)*ZN2(gI2,3) +
      5*g2*KroneckerDelta(1,gO2)*ZN1(gI1,1)*ZN2(gI2,3) + 5*Conj(LamTU)*
      KroneckerDelta(3,gO2)*ZN1(gI1,3)*ZN2(gI2,3) - 7.0710678118654755*Conj(LamSU)
      *ZN1(gI1,3)*(KroneckerDelta(1,gO2)*ZN2(gI2,0) + KroneckerDelta(2,gO2)*ZN2(
      gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiChiUhhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = 0.1*(3.872983346207417*g1*Conj(ZN1(gI2,0
      ))*(Conj(ZN2(gI1,2))*KroneckerDelta(0,gO1) - Conj(ZN2(gI1,3))*KroneckerDelta
      (1,gO1)) + 5*Conj(ZN1(gI2,2))*(1.4142135623730951*LamSD*Conj(ZN2(gI1,0))*
      KroneckerDelta(0,gO1) + LamTD*Conj(ZN2(gI1,1))*KroneckerDelta(0,gO1) + Conj(
      ZN2(gI1,2))*(1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + LamTD*
      KroneckerDelta(3,gO1))) + 5*(Conj(ZN1(gI2,1))*(-(g2*Conj(ZN2(gI1,2))*
      KroneckerDelta(0,gO1)) + g2*Conj(ZN2(gI1,3))*KroneckerDelta(1,gO1)) + Conj(
      ZN1(gI2,3))*(-1.4142135623730951*LamSU*Conj(ZN2(gI1,0))*KroneckerDelta(1,gO1
      ) + LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(1,gO1) + Conj(ZN2(gI1,3))*(
      -1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + LamTU*KroneckerDelta(3,gO1
      )))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,
      2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      ZD(gI2,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(
      gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1
      ,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1)
      )*ZE(gI2,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(
      ZE(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(
      gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) -
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,
      2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      ZU(gI2,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(
      gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSdconjSd(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.016666666666666666*(-2*(-15*g2*(MDWBT
      + Conj(MDWBT))*KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1))
      + 3.872983346207417*g1*(MDBS + Conj(MDBS))*KroneckerDelta(2,gO2)*(SUM(j1,0,
      2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) + 2*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3
      + j1)))) - 3*KroneckerDelta(1,gO2)*(vu*(Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj
      (ZD(gI2,j1))*ZD(gI1,j1)) + 2*(vu*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(
      gI1,3 + j1)) - 7.0710678118654755*(Conj(Mu)*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(
      j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))
      *Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2))))) + 3*vd*KroneckerDelta(0,gO2)*((Sqr(g1)
      + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) + 2*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI1,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI1
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSeconjSe(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.05*(2*(5*g2*(MDWBT + Conj(MDWBT))*
      KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) +
      3.872983346207417*g1*(MDBS + Conj(MDBS))*KroneckerDelta(2,gO2)*(SUM(j1,0,2,
      Conj(ZE(gI2,j1))*ZE(gI1,j1)) - 2*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 +
      j1)))) + KroneckerDelta(1,gO2)*(vu*(3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(
      ZE(gI2,j1))*ZE(gI1,j1)) - 6*vu*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(
      gI1,3 + j1)) + 14.142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(
      j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))
      *Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2)))) + vd*KroneckerDelta(0,gO2)*((-3*Sqr(g1)
      + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1)) + 6*Sqr(g1)*SUM(j1,0,2
      ,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZE(gI2,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI1,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI1
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhSuconjSu(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.016666666666666666*(-2*(15*g2*(MDWBT +
      Conj(MDWBT))*KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1))
      + 3.872983346207417*g1*(MDBS + Conj(MDBS))*KroneckerDelta(2,gO2)*(SUM(j1,0,2
      ,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 +
      j1)))) + 3*KroneckerDelta(0,gO2)*(vd*(Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(
      ZU(gI2,j1))*ZU(gI1,j1)) - 4*vd*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(
      gI1,3 + j1)) + 14.142135623730951*(Conj(Mu)*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(
      j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))
      *Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2)))) - 3*vu*KroneckerDelta(1,gO2)*((Sqr(g1)
      - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*Sqr(g1)*SUM(j1,0,2,
      Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gI2,3 + j3))*
      SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2))) + SUM(j3,0
      ,2,SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI1
      ,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUhhHpm(int gO2, int gI2) const
{
   const std::complex<double> result = 0.5*Conj(Mu)*(-1.4142135623730951*LamSD*
      KroneckerDelta(2,gO2)*ZP(gI2,1) + LamTD*KroneckerDelta(3,gO2)*ZP(gI2,1) +
      1.4142135623730951*LamTD*KroneckerDelta(1,gO2)*ZP(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUhhHpmconjSRum(int gO2, int gI2) const
{
   const std::complex<double> result = 0.5*Mu*(1.4142135623730951*Conj(LamSU)*
      KroneckerDelta(2,gO2)*ZP(gI2,0) + Conj(LamTU)*(KroneckerDelta(3,gO2)*ZP(gI2,
      0) - 1.4142135623730951*KroneckerDelta(0,gO2)*ZP(gI2,3)));

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
   const std::complex<double> result = -0.5*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0)
      - KroneckerDelta(1,gO2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(3,gO2
      )*(ZP(gI2,2) + ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUAhconjSRdp(int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.05)*((
      7.745966692414834*g1*MDBS + 7.0710678118654755*(-2*MuD + LamTD*vT)*Conj(
      LamSD) - 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj
      (MDBS) + 14.142135623730951*LamSD*Conj(MuD))*KroneckerDelta(2,gO2) + 5*(2*g2
      *MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) +
      1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(
      MuD))*KroneckerDelta(3,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUAhconjSRum(int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*((
      7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU
      ) - 7.0710678118654755*LamSU*vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS
      ) - 14.142135623730951*LamSU*Conj(MuU))*KroneckerDelta(2,gO2) + 5*(2*g2*
      MDWBT - 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) +
      1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(
      MuU))*KroneckerDelta(3,gO2));

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

std::complex<double> CLASSNAME::CpSRdpUAhUAhconjSRdp(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(5*(Conj(LamTD)*(1.4142135623730951
      *LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*
      LamSD*KroneckerDelta(2,gO1) - 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(
      3,gO2)) + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*LamSD*KroneckerDelta(2,gO2
      ) + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1
      )*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUAhUAhconjSRum(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-5*(Conj(LamTU)*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + 2*LamTU*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2)) + Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSU
      *KroneckerDelta(2,gO2) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)))) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-20*AbsSqr(LamTU) + 3*Sqr(g1)
      - 5*Sqr(g2)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-3*Sqr(g1) + 5*
      Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(int gO1, int gO2) const
{
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhRhconjRh(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + 2*LamTD*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))*ZHR(gI1,0)*ZHR(gI2,0) + Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTD*KroneckerDelta(3,gO2)))*ZHR(gI1,0)*ZHR(gI2,0) + (-(Conj(LamTU)*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*LamTU*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*
      LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)
      )))*ZHR(gI1,1)*ZHR(gI2,1)) - KroneckerDelta(1,gO1)*(5*KroneckerDelta(0,gO2)*
      (-2*LamSU*Conj(LamSD)*ZHR(gI1,1)*ZHR(gI2,0) + LamTU*Conj(LamTD)*ZHR(gI1,1)*
      ZHR(gI2,0) + (-2*LamSD*Conj(LamSU) + LamTD*Conj(LamTU))*ZHR(gI1,0)*ZHR(gI2,1
      )) + KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) +
      (20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,1)*ZHR
      (gI2,1))) + KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*(2*LamSU*Conj(
      LamSD)*ZHR(gI1,1)*ZHR(gI2,0) - LamTU*Conj(LamTD)*ZHR(gI1,1)*ZHR(gI2,0) + (2*
      LamSD*Conj(LamSU) - LamTD*Conj(LamTU))*ZHR(gI1,0)*ZHR(gI2,1)) +
      KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5
      *Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2
      ,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhRhconjRh(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,-0.05)*(
      KroneckerDelta(2,gO2)*((7.745966692414834*g1*MDBS - 7.0710678118654755*(2*
      MuD + LamTD*vT)*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) -
      7.745966692414834*g1*Conj(MDBS) + 14.142135623730951*LamSD*Conj(MuD))*ZHR(
      gI1,0)*ZHR(gI2,0) + (-7.745966692414834*g1*MDBS + 7.0710678118654755*(-2*MuU
      + LamTU*vT)*Conj(LamSU) - 7.0710678118654755*LamSU*vT*Conj(LamTU) +
      7.745966692414834*g1*Conj(MDBS) + 14.142135623730951*LamSU*Conj(MuU))*ZHR(
      gI1,1)*ZHR(gI2,1)) - 5*((vu*KroneckerDelta(0,gO2) - vd*KroneckerDelta(1,gO2)
      )*(2*LamSD*Conj(LamSU)*ZHR(gI1,1)*ZHR(gI2,0) - LamTD*Conj(LamTU)*ZHR(gI1,1)*
      ZHR(gI2,0) + (-2*LamSU*Conj(LamSD) + LamTU*Conj(LamTD))*ZHR(gI1,0)*ZHR(gI2,1
      )) + KroneckerDelta(3,gO2)*((2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(
      LamSD) + (2*MuD + 1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(MDWBT
      ) - 2*LamTD*Conj(MuD))*ZHR(gI1,0)*ZHR(gI2,0) + (-2*g2*MDWBT +
      1.4142135623730951*LamTU*vS*Conj(LamSU) - (2*MuU + 1.4142135623730951*LamSU*
      vS)*Conj(LamTU) + 2*g2*Conj(MDWBT) + 2*LamTU*Conj(MuU))*ZHR(gI1,1)*ZHR(gI2,1
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1UAhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*(
      KroneckerDelta(3,gO2)*(2*g2*UM1(gI2,0)*UP1(gI1,0) + Conj(LamTD)*UM1(gI2,1)*
      UP1(gI1,1)) + 1.4142135623730951*(-(Conj(LamSD)*KroneckerDelta(2,gO2)*UM1(
      gI2,1)*UP1(gI1,1)) + KroneckerDelta(0,gO2)*(g2*UM1(gI2,1)*UP1(gI1,0) - Conj(
      LamTD)*UM1(gI2,0)*UP1(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1UAhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(Conj(UM1(
      gI1,0))*(-1.4142135623730951*LamTD*Conj(UP1(gI2,1))*KroneckerDelta(0,gO1) +
      2*g2*Conj(UP1(gI2,0))*KroneckerDelta(3,gO1)) + Conj(UM1(gI1,1))*(
      1.4142135623730951*g2*Conj(UP1(gI2,0))*KroneckerDelta(0,gO1) + Conj(UP1(gI2,
      1))*(-1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + LamTD*KroneckerDelta(
      3,gO1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2UAhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(
      KroneckerDelta(3,gO2)*(2*g2*UM2(gI1,0)*UP2(gI2,0) - Conj(LamTU)*UM2(gI1,1)*
      UP2(gI2,1)) - 1.4142135623730951*(Conj(LamTU)*KroneckerDelta(1,gO2)*UM2(gI1,
      1)*UP2(gI2,0) + (g2*KroneckerDelta(1,gO2)*UM2(gI1,0) + Conj(LamSU)*
      KroneckerDelta(2,gO2)*UM2(gI1,1))*UP2(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2UAhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Conj(UM2
      (gI2,0))*(1.4142135623730951*Conj(UP2(gI1,1))*KroneckerDelta(1,gO1) - 2*Conj
      (UP2(gI1,0))*KroneckerDelta(3,gO1)) + Conj(UM2(gI2,1))*(1.4142135623730951*
      LamTU*Conj(UP2(gI1,0))*KroneckerDelta(1,gO1) + Conj(UP2(gI1,1))*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + LamTU*KroneckerDelta(3,gO1)
      )));

   return result;
}

std::complex<double> CLASSNAME::CpAhUAhconjRh(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = -0.25*Mu*(2*Conj(LamSD)*(KroneckerDelta(
      2,gO2)*ZA(gI2,1) + KroneckerDelta(1,gO2)*ZA(gI2,2))*ZHR(gI1,0) +
      1.4142135623730951*Conj(LamTD)*(KroneckerDelta(3,gO2)*ZA(gI2,1) +
      KroneckerDelta(1,gO2)*ZA(gI2,3))*ZHR(gI1,0) + (-2*Conj(LamSU)*(
      KroneckerDelta(2,gO2)*ZA(gI2,0) + KroneckerDelta(0,gO2)*ZA(gI2,2)) +
      1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZA(gI2,0) +
      KroneckerDelta(0,gO2)*ZA(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhconjRh(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.25)*Mu*(2*Conj(
      LamSD)*(KroneckerDelta(2,gO2)*ZH(gI2,1) - KroneckerDelta(1,gO2)*ZH(gI2,2))*
      ZHR(gI1,0) + 1.4142135623730951*Conj(LamTD)*(KroneckerDelta(3,gO2)*ZH(gI2,1)
      - KroneckerDelta(1,gO2)*ZH(gI2,3))*ZHR(gI1,0) + (Conj(LamSU)*(-2*
      KroneckerDelta(2,gO2)*ZH(gI2,0) + 2*KroneckerDelta(0,gO2)*ZH(gI2,2)) +
      1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZH(gI2,0) -
      KroneckerDelta(0,gO2)*ZH(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSvconjSv(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSvconjSv(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(
      3.872983346207417*g1*(MDBS - Conj(MDBS))*KroneckerDelta(2,gO2) + 5*g2*(
      -MDWBT + Conj(MDWBT))*KroneckerDelta(3,gO2))*KroneckerDelta(gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,
      j2))*ZDR(gI2,j1))*ZDL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdUAhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*KroneckerDelta(0,gO1)*SUM(j2,0,2,Conj(ZDL(gI2,j2))*SUM(
      j1,0,2,Conj(ZDR(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,
      j2))*ZER(gI2,j1))*ZEL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeUAhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*KroneckerDelta(0,gO1)*SUM(j2,0,2,Conj(ZEL(gI2,j2))*SUM(
      j1,0,2,Conj(ZER(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*KroneckerDelta(1,gO2)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,
      j2))*ZUR(gI2,j1))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuUAhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(
      j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAhUAh(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-5*(4*AbsSqr(LamSU)*KroneckerDelta
      (2,gO1)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) - 1.4142135623730951*LamTU
      *Conj(LamSU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZA(gI1,1)*ZA(gI2,1)
      - 1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*KroneckerDelta
      (3,gO1)*ZA(gI1,1)*ZA(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO1)*KroneckerDelta(3,
      gO2)*ZA(gI1,1)*ZA(gI2,1) + 2*AbsSqr(LamTU)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 4*AbsSqr(LamSU)*KroneckerDelta(1
      ,gO2)*KroneckerDelta(2,gO1)*ZA(gI1,2)*ZA(gI2,1) - 1.4142135623730951*LamTU*
      Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*ZA(gI1,2)*ZA(gI2,1)
      - 1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      3,gO1)*ZA(gI1,2)*ZA(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZA(gI1,3)*ZA(gI2,1) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,
      gO1)*ZA(gI1,3)*ZA(gI2,1) + 2*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZA(gI1,3)*ZA(gI2,1) + 4*AbsSqr(LamSU)*KroneckerDelta(1
      ,gO2)*KroneckerDelta(2,gO1)*ZA(gI1,1)*ZA(gI2,2) - 1.4142135623730951*LamTU*
      Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*ZA(gI1,1)*ZA(gI2,2)
      - 1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      3,gO1)*ZA(gI1,1)*ZA(gI2,2) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZA(gI1,1)*ZA(gI2,3) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,
      gO1)*ZA(gI1,1)*ZA(gI2,3) + 2*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZA(gI1,1)*ZA(gI2,3) + Conj(LamSD)*(1.4142135623730951*
      LamTD*KroneckerDelta(3,gO1)*(KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) +
      KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2))) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) +
      1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) +
      KroneckerDelta(0,gO2)*(4*LamSD*ZA(gI1,2)*ZA(gI2,0) + 1.4142135623730951*
      LamTD*ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*(4*LamSD*ZA(gI2,2) +
      1.4142135623730951*LamTD*ZA(gI2,3))))) + Conj(LamTD)*(1.4142135623730951*
      LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZA(gI1,0)*ZA(gI2,0) +
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1)*(KroneckerDelta(3,gO2)*ZA(gI1
      ,0)*ZA(gI2,0) + KroneckerDelta(0,gO2)*(ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*ZA(
      gI2,3))) + KroneckerDelta(3,gO1)*(2*LamTD*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA
      (gI2,0) + KroneckerDelta(0,gO2)*(1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,0
      ) + 2*LamTD*ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*(1.4142135623730951*LamSD*ZA(gI2
      ,2) + 2*LamTD*ZA(gI2,3)))))) - KroneckerDelta(0,gO1)*(-(KroneckerDelta(1,gO2
      )*(3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1))) + 5*(
      Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*(ZA(gI1,2)*ZA(
      gI2,0) + ZA(gI1,0)*ZA(gI2,2)) + KroneckerDelta(2,gO2)*(4*LamSD*ZA(gI1,2)*ZA(
      gI2,0) + 1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*(4*LamSD*
      ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3)))) + Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*(ZA(gI1,3)*ZA(gI2,0) + ZA(gI1
      ,0)*ZA(gI2,3)) + KroneckerDelta(3,gO2)*(1.4142135623730951*LamSD*ZA(gI1,2)*
      ZA(gI2,0) + 2*LamTD*ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*(1.4142135623730951*
      LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3))))) + KroneckerDelta(0,gO2)*(3*(3*Sqr(g1
      ) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(
      gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(
      gI1,2)*(4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))) + Conj(
      LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(
      1.4142135623730951*LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3)))))) + KroneckerDelta
      (1,gO1)*(KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0)
      + ZA(gI1,0)*ZA(gI2,1)) + 5*(Conj(LamTU)*(1.4142135623730951*LamSU*
      KroneckerDelta(2,gO2)*(ZA(gI1,3)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,3)) +
      KroneckerDelta(3,gO2)*(1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,1) - 2*
      LamTU*ZA(gI1,3)*ZA(gI2,1) + ZA(gI1,1)*(1.4142135623730951*LamSU*ZA(gI2,2) -
      2*LamTU*ZA(gI2,3)))) + Conj(LamSU)*(1.4142135623730951*LamTU*KroneckerDelta(
      3,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)) + KroneckerDelta(2,gO2)*(
      -4*LamSU*ZA(gI1,2)*ZA(gI2,1) + 1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,1)
      + ZA(gI1,1)*(-4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3))))) +
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - 3*(3*
      Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTU)*(
      1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSU*ZA(gI2,2) - 2*LamTU*ZA(gI2,3))) + Conj(LamSU)*(1.4142135623730951*
      LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSU*ZA(gI2,2) +
      1.4142135623730951*LamTU*ZA(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + 2*LamTD*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))*ZH(gI1,0)*ZH(gI2,0) + Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTD*KroneckerDelta(3,gO2)))*ZH(gI1,0)*ZH(gI2,0) + (-(Conj(LamTU)*(
      1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (
      1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*LamTU*KroneckerDelta(3,
      gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*
      LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)
      )))*ZH(gI1,1)*ZH(gI2,1)) - KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*
      Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)
      *ZH(gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,2) +
      ZH(gI1,2)*(4*LamSD*ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3))) + Conj(
      LamTD)*(1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(
      1.4142135623730951*LamSD*ZH(gI2,2) + 2*LamTD*ZH(gI2,3))))) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) -
      (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 5*(Conj(LamTU)*(
      1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951
      *LamSU*ZH(gI2,2) - 2*LamTU*ZH(gI2,3))) + Conj(LamSU)*(1.4142135623730951*
      LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSU*ZH(gI2,2) +
      1.4142135623730951*LamTU*ZH(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhHpmconjHpm(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*(5*(
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,2)*ZP(gI2,0) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 2*
      LamSD*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,0) +
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,3)*ZP(gI2,0) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) -
      2*LamSD*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,3) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) + 2*
      LamTD*Conj(LamSD)*KroneckerDelta(2,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,0)*
      ZP(gI2,3))) - KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(
      gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*AbsSqr(
      LamTD) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)))) +
      KroneckerDelta(1,gO1)*(5*(-1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3
      ,gO2)*ZP(gI1,2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)
      *ZP(gI1,2)*ZP(gI2,1) - 2*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*ZP(gI1,3)*
      ZP(gI2,1) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,3)
      *ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(
      gI2,1) + KroneckerDelta(0,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(
      gI2,1)) + 2*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,2) -
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) -
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,3) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) + 2*
      LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*(ZP(gI1,2)*ZP(gI2,1) - ZP(gI1,1)*ZP(
      gI2,3))) + KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0
      ) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(
      gI2,2) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3)))) - 5*(
      1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,2)*ZP(gI2,0) + 1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3
      ,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) + 4*AbsSqr(LamSU)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamTU*Conj(
      LamSU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,
      gO1)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO1)*KroneckerDelta(3,
      gO2)*ZP(gI1,1)*ZP(gI2,1) + 2*AbsSqr(LamTU)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 2*LamTU*Conj(LamSU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,2)*ZP(gI2,1) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)
      *ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) + 2*LamSU*Conj(LamTU)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)
      *ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) + 1.4142135623730951*
      KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) - 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,1)*ZP(
      gI2,2) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,2) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,3) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      2,gO1)*ZP(gI1,1)*ZP(gI2,3) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta
      (1,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,3) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3) -
      Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*KroneckerDelta(3
      ,gO1)*ZP(gI1,0)*ZP(gI2,0) + LamSD*KroneckerDelta(2,gO1)*(1.4142135623730951*
      KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 2*KroneckerDelta(0,gO2)*(ZP(gI1,
      3)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,2))) + LamTD*KroneckerDelta(3,gO1)*(-2*
      KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 1.4142135623730951*
      KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)
      *(ZP(gI2,2) + ZP(gI2,3))))) + Conj(LamSD)*(-1.4142135623730951*LamTD*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,0)*ZP(gI2,0) +
      KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0) -
      LamTD*(1.4142135623730951*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 2*
      KroneckerDelta(0,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,0)*ZP(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUAh(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      7.745966692414834*g1*MDBS*KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(gI2,0) +
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(gI2,0)
      + 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,2)*
      ZA(gI2,0) - 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA
      (gI1,2)*ZA(gI2,0) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(0,gO2)*ZA
      (gI1,2)*ZA(gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(0,gO2)
      *ZA(gI1,2)*ZA(gI2,0) - 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(gI2,0)
      - 7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,3)*
      ZA(gI2,0) + 10*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(gI2,0) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(
      gI2,0) + 10*g2*Conj(MDWBT)*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(gI2,0) - 10*
      LamTD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(gI2,0) -
      7.745966692414834*g1*MDBS*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1) +
      14.142135623730951*MuU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1)
      - 7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,2)*
      ZA(gI2,1) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA
      (gI1,2)*ZA(gI2,1) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(1,gO2)*ZA
      (gI1,2)*ZA(gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(1,gO2)
      *ZA(gI1,2)*ZA(gI2,1) + 10*g2*MDWBT*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA(gI2,1)
      + 7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,3)*
      ZA(gI2,1) - 10*MuU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA(gI2,1) -
      7.0710678118654755*LamSU*vS*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA(
      gI2,1) - 10*g2*Conj(MDWBT)*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA(gI2,1) + 10*
      LamTU*Conj(MuU)*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA(gI2,1) + KroneckerDelta(2
      ,gO2)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD + LamTD*vT)*
      Conj(LamSD) - 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1
      *Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZA(gI1,0)*ZA(gI2,0) + (
      -7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuU - LamTU*vT)*Conj(
      LamSU) + 7.0710678118654755*LamSU*vT*Conj(LamTU) + 7.745966692414834*g1*Conj
      (MDBS) - 14.142135623730951*LamSU*Conj(MuU))*ZA(gI1,1)*ZA(gI2,1)) - 5*
      KroneckerDelta(3,gO2)*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD)
      - 2*MuD*Conj(LamTD) - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(
      MDWBT) + 2*LamTD*Conj(MuD))*ZA(gI1,0)*ZA(gI2,0) + (-2*g2*MDWBT -
      1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) +
      1.4142135623730951*LamSU*vS*Conj(LamTU) + 2*g2*Conj(MDWBT) - 2*LamTU*Conj(
      MuU))*ZA(gI1,1)*ZA(gI2,1)) + 7.745966692414834*g1*MDBS*KroneckerDelta(0,gO2)
      *ZA(gI1,0)*ZA(gI2,2) + 14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(0,
      gO2)*ZA(gI1,0)*ZA(gI2,2) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*
      KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,2) - 7.0710678118654755*LamSD*vT*Conj
      (LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,2) - 7.745966692414834*g1*
      Conj(MDBS)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,2) - 14.142135623730951*
      LamSD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,2) -
      7.745966692414834*g1*MDBS*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2) +
      14.142135623730951*MuU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2)
      - 7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*
      ZA(gI2,2) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA
      (gI1,1)*ZA(gI2,2) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(1,gO2)*ZA
      (gI1,1)*ZA(gI2,2) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(1,gO2)
      *ZA(gI1,1)*ZA(gI2,2) - 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,3)
      - 7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*
      ZA(gI2,3) + 10*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,3) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(
      gI2,3) + 10*g2*Conj(MDWBT)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,3) - 10*
      LamTD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2,3) + 10*g2*MDWBT*
      KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3) + 7.0710678118654755*LamTU*vS*Conj
      (LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3) - 10*MuU*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3) - 7.0710678118654755*LamSU*vS*Conj
      (LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3) - 10*g2*Conj(MDWBT)*
      KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3) + 10*LamTU*Conj(MuU)*
      KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpAhUAhhh(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = 0.05*(-5*(vd*Conj(LamSD)*(
      1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*ZA(gI2,2) + KroneckerDelta(2,
      gO2)*(4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3)))*ZH(gI1,0) +
      vd*Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*ZA(gI2,3) +
      KroneckerDelta(3,gO2)*(1.4142135623730951*LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3
      )))*ZH(gI1,0) + vu*(-(Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2
      ,gO2)*ZA(gI2,3) + KroneckerDelta(3,gO2)*(1.4142135623730951*LamSU*ZA(gI2,2)
      - 2*LamTU*ZA(gI2,3)))) - Conj(LamSU)*(1.4142135623730951*LamTU*
      KroneckerDelta(3,gO2)*ZA(gI2,2) + KroneckerDelta(2,gO2)*(-4*LamSU*ZA(gI2,2)
      + 1.4142135623730951*LamTU*ZA(gI2,3))))*ZH(gI1,1)) - KroneckerDelta(0,gO2)*
      ZA(gI2,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))
      *ZH(gI1,1) - 7.745966692414834*g1*MDBS*ZH(gI1,2) + 20*vS*AbsSqr(LamSD)*ZH(
      gI1,2) + 14.142135623730951*MuD*Conj(LamSD)*ZH(gI1,2) + 7.0710678118654755*
      LamTD*vT*Conj(LamSD)*ZH(gI1,2) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(
      gI1,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gI1,2) + 14.142135623730951*
      LamSD*Conj(MuD)*ZH(gI1,2) + 10*g2*MDWBT*ZH(gI1,3) + 10*vT*AbsSqr(LamTD)*ZH(
      gI1,3) + 7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gI1,3) + 10*MuD*Conj(
      LamTD)*ZH(gI1,3) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gI1,3) + 10*g2
      *Conj(MDWBT)*ZH(gI1,3) + 10*LamTD*Conj(MuD)*ZH(gI1,3)) + KroneckerDelta(1,
      gO2)*ZA(gI2,1)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr
      (g2))*ZH(gI1,1) - 7.745966692414834*g1*MDBS*ZH(gI1,2) - 20*vS*AbsSqr(LamSU)*
      ZH(gI1,2) - 14.142135623730951*MuU*Conj(LamSU)*ZH(gI1,2) +
      7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gI1,2) + 7.0710678118654755*LamSU
      *vT*Conj(LamTU)*ZH(gI1,2) - 7.745966692414834*g1*Conj(MDBS)*ZH(gI1,2) -
      14.142135623730951*LamSU*Conj(MuU)*ZH(gI1,2) + 10*g2*MDWBT*ZH(gI1,3) - 10*vT
      *AbsSqr(LamTU)*ZH(gI1,3) + 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gI1,3)
      + 10*MuU*Conj(LamTU)*ZH(gI1,3) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH
      (gI1,3) + 10*g2*Conj(MDWBT)*ZH(gI1,3) + 10*LamTU*Conj(MuU)*ZH(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*(-5*
      KroneckerDelta(3,gO2)*(1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gI1,2)*ZH(
      gI2,0) - 1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gI1,2)*ZH(gI2,0) - 2*g2*
      MDWBT*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gI1,1
      )*ZH(gI2,1) + 2*MuU*Conj(LamTU)*ZH(gI1,1)*ZH(gI2,1) + 1.4142135623730951*
      LamSU*vS*Conj(LamTU)*ZH(gI1,1)*ZH(gI2,1) + 2*g2*Conj(MDWBT)*ZH(gI1,1)*ZH(gI2
      ,1) - 2*LamTU*Conj(MuU)*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*vu*
      Conj(LamSU)*ZH(gI1,2)*ZH(gI2,1) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH
      (gI1,2)*ZH(gI2,1) - 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gI1,1)*ZH(gI2
      ,2) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gI1,1)*ZH(gI2,2) + ZH(gI1,0
      )*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) - 2*MuD*Conj(LamTD)
      - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) + 2*LamTD*Conj
      (MuD))*ZH(gI2,0) + 1.4142135623730951*vd*(LamTD*Conj(LamSD) - LamSD*Conj(
      LamTD))*ZH(gI2,2))) + KroneckerDelta(2,gO2)*(7.0710678118654755*LamTD*vd*
      Conj(LamSD)*ZH(gI1,3)*ZH(gI2,0) - 7.0710678118654755*LamSD*vd*Conj(LamTD)*ZH
      (gI1,3)*ZH(gI2,0) - 7.745966692414834*g1*MDBS*ZH(gI1,1)*ZH(gI2,1) +
      14.142135623730951*MuU*Conj(LamSU)*ZH(gI1,1)*ZH(gI2,1) - 7.0710678118654755*
      LamTU*vT*Conj(LamSU)*ZH(gI1,1)*ZH(gI2,1) + 7.0710678118654755*LamSU*vT*Conj(
      LamTU)*ZH(gI1,1)*ZH(gI2,1) + 7.745966692414834*g1*Conj(MDBS)*ZH(gI1,1)*ZH(
      gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*ZH(gI1,1)*ZH(gI2,1) -
      7.0710678118654755*LamTU*vu*Conj(LamSU)*ZH(gI1,3)*ZH(gI2,1) +
      7.0710678118654755*LamSU*vu*Conj(LamTU)*ZH(gI1,3)*ZH(gI2,1) -
      7.0710678118654755*LamTU*vu*Conj(LamSU)*ZH(gI1,1)*ZH(gI2,3) +
      7.0710678118654755*LamSU*vu*Conj(LamTU)*ZH(gI1,1)*ZH(gI2,3) + ZH(gI1,0)*((
      7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD + LamTD*vT)*Conj(LamSD
      ) - 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS
      ) - 14.142135623730951*LamSD*Conj(MuD))*ZH(gI2,0) + 7.0710678118654755*vd*(
      LamTD*Conj(LamSD) - LamSD*Conj(LamTD))*ZH(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHpmconjHpm(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      KroneckerDelta(2,gO2)*(-10*LamSD*vd*Conj(LamTD)*ZP(gI1,2)*ZP(gI2,0) + 10*
      LamTD*vd*Conj(LamSD)*ZP(gI1,3)*ZP(gI2,0) - 7.745966692414834*g1*MDBS*ZP(gI1,
      1)*ZP(gI2,1) + 14.142135623730951*MuU*Conj(LamSU)*ZP(gI1,1)*ZP(gI2,1) +
      7.0710678118654755*LamTU*vT*Conj(LamSU)*ZP(gI1,1)*ZP(gI2,1) -
      7.0710678118654755*LamSU*vT*Conj(LamTU)*ZP(gI1,1)*ZP(gI2,1) +
      7.745966692414834*g1*Conj(MDBS)*ZP(gI1,1)*ZP(gI2,1) - 14.142135623730951*
      LamSU*Conj(MuU)*ZP(gI1,1)*ZP(gI2,1) - 10*LamSU*vu*Conj(LamTU)*ZP(gI1,2)*ZP(
      gI2,1) + 10*LamTU*vu*Conj(LamSU)*ZP(gI1,3)*ZP(gI2,1) + 10*LamTU*vu*Conj(
      LamSU)*ZP(gI1,1)*ZP(gI2,2) - 10*LamSU*vu*Conj(LamTU)*ZP(gI1,1)*ZP(gI2,3) +
      ZP(gI1,0)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT
      )*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*
      g1*Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZP(gI2,0) + 10*vd*(LamTD
      *Conj(LamSD)*ZP(gI2,2) - LamSD*Conj(LamTD)*ZP(gI2,3)))) + 5*(KroneckerDelta(
      0,gO2)*(vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,0) + ((2.8284271247461903*MuD + 2*LamSD*
      vS - 1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(g2*vT
      + 2*Conj(MDWBT)))*ZP(gI1,2)*ZP(gI2,0) + 2.8284271247461903*g2*MDWBT*ZP(gI1,
      3)*ZP(gI2,0) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gI1,3)*ZP(gI2,0) + 2*
      LamTD*vS*Conj(LamSD)*ZP(gI1,3)*ZP(gI2,0) + 2.8284271247461903*LamTD*Conj(MuD
      )*ZP(gI1,3)*ZP(gI2,0) - 1.4142135623730951*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) -
      vu*Sqr(g2)*ZP(gI1,0)*ZP(gI2,1) - 2.8284271247461903*g2*MDWBT*ZP(gI1,0)*ZP(
      gI2,2) + 1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gI1,0)*ZP(gI2,2) - 2*LamTD*
      vS*Conj(LamSD)*ZP(gI1,0)*ZP(gI2,2) - 2.8284271247461903*LamTD*Conj(MuD)*ZP(
      gI1,0)*ZP(gI2,2) - 1.4142135623730951*vT*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gI1,0)*ZP(gI2,3) - 2.8284271247461903
      *MuD*Conj(LamTD)*ZP(gI1,0)*ZP(gI2,3) - 2*LamSD*vS*Conj(LamTD)*ZP(gI1,0)*ZP(
      gI2,3) - 2.8284271247461903*g2*Conj(MDWBT)*ZP(gI1,0)*ZP(gI2,3) +
      1.4142135623730951*vT*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3)) + KroneckerDelta(1,gO2)*(
      -((vd*Sqr(g2)*ZP(gI1,0) + ((2.8284271247461903*MuU + 2*LamSU*vS -
      1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(g2*vT + 2*
      Conj(MDWBT)))*ZP(gI1,2) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2
      *MDWBT + vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) - vT*Sqr(g2)))*ZP(gI1,3))*ZP(
      gI2,1)) + ZP(gI1,1)*(vd*Sqr(g2)*ZP(gI2,0) + (2*LamTU*vS*Conj(LamSU) +
      1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*
      Sqr(g2)))*ZP(gI2,2) + ((2.8284271247461903*MuU + 2*LamSU*vS +
      1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) +
      2*Conj(MDWBT)))*ZP(gI2,3))) + KroneckerDelta(3,gO2)*(1.4142135623730951*vd*
      AbsSqr(LamTD)*ZP(gI1,3)*ZP(gI2,0) - 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,3)*
      ZP(gI2,0) - 2*g2*MDWBT*ZP(gI1,1)*ZP(gI2,1) - 1.4142135623730951*LamTU*vS*
      Conj(LamSU)*ZP(gI1,1)*ZP(gI2,1) + 2*MuU*Conj(LamTU)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamSU*vS*Conj(LamTU)*ZP(gI1,1)*ZP(gI2,1) + 2*g2*Conj(
      MDWBT)*ZP(gI1,1)*ZP(gI2,1) - 2*LamTU*Conj(MuU)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951
      *vu*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gI1
      ,1)*ZP(gI2,2) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) - 4*vT*Sqr
      (g2)*ZP(gI1,3)*ZP(gI2,2) - 1.4142135623730951*vu*AbsSqr(LamTU)*ZP(gI1,1)*ZP(
      gI2,3) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) + ZP(gI1,2)*(
      1.4142135623730951*vd*(AbsSqr(LamTD) - Sqr(g2))*ZP(gI2,0) +
      1.4142135623730951*vu*(AbsSqr(LamTU) - Sqr(g2))*ZP(gI2,1) + 4*vT*Sqr(g2)*ZP(
      gI2,3)) + ZP(gI1,0)*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) -
      2*MuD*Conj(LamTD) - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(
      MDWBT) + 2*LamTD*Conj(MuD))*ZP(gI2,0) + 1.4142135623730951*vd*(-AbsSqr(LamTD
      ) + Sqr(g2))*(ZP(gI2,2) + ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiChiUAhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.1)*(5*Conj(
      LamTD)*KroneckerDelta(0,gO2)*ZN1(gI1,2)*ZN2(gI2,1) + 5*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZN1(gI1,3)*ZN2(gI2,1) - 3.872983346207417*g1*
      KroneckerDelta(0,gO2)*ZN1(gI1,0)*ZN2(gI2,2) + 5*g2*KroneckerDelta(0,gO2)*ZN1
      (gI1,1)*ZN2(gI2,2) + 5*Conj(LamTD)*KroneckerDelta(3,gO2)*ZN1(gI1,2)*ZN2(gI2,
      2) + 7.0710678118654755*Conj(LamSD)*ZN1(gI1,2)*(KroneckerDelta(0,gO2)*ZN2(
      gI2,0) + KroneckerDelta(2,gO2)*ZN2(gI2,2)) + 3.872983346207417*g1*
      KroneckerDelta(1,gO2)*ZN1(gI1,0)*ZN2(gI2,3) - 5*g2*KroneckerDelta(1,gO2)*ZN1
      (gI1,1)*ZN2(gI2,3) + 5*Conj(LamTU)*KroneckerDelta(3,gO2)*ZN1(gI1,3)*ZN2(gI2,
      3) - 7.0710678118654755*Conj(LamSU)*ZN1(gI1,3)*(KroneckerDelta(1,gO2)*ZN2(
      gI2,0) + KroneckerDelta(2,gO2)*ZN2(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiChiUAhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(
      3.872983346207417*g1*Conj(ZN1(gI2,0))*(-(Conj(ZN2(gI1,2))*KroneckerDelta(0,
      gO1)) + Conj(ZN2(gI1,3))*KroneckerDelta(1,gO1)) + 5*Conj(ZN1(gI2,2))*(
      1.4142135623730951*LamSD*Conj(ZN2(gI1,0))*KroneckerDelta(0,gO1) + LamTD*Conj
      (ZN2(gI1,1))*KroneckerDelta(0,gO1) + Conj(ZN2(gI1,2))*(1.4142135623730951*
      LamSD*KroneckerDelta(2,gO1) + LamTD*KroneckerDelta(3,gO1))) + 5*(g2*Conj(ZN1
      (gI2,1))*(Conj(ZN2(gI1,2))*KroneckerDelta(0,gO1) - Conj(ZN2(gI1,3))*
      KroneckerDelta(1,gO1)) + Conj(ZN1(gI2,3))*(-1.4142135623730951*LamSU*Conj(
      ZN2(gI1,0))*KroneckerDelta(1,gO1) + LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(1,
      gO1) + Conj(ZN2(gI1,3))*(-1.4142135623730951*LamSU*KroneckerDelta(2,gO1) +
      LamTU*KroneckerDelta(3,gO1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSdconjSd(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,
      2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      ZD(gI2,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yd(j3,j1))*Yd(j2,j1))*ZD(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(
      gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZD(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSeconjSe(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*SUM(j1
      ,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1)
      )*ZE(gI2,3 + j1)) - 20*(SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(
      ZE(gI1,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,j2)))*ZE(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhSuconjSu(int gO1, int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(
      gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1))) -
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,
      2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      ZU(gI2,3 + j1)) + 20*(SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Yu(j3,j1))*Yu(j2,j1))*ZU(gI2,3 + j2))) + SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(
      gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSdconjSd(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,
      -0.03333333333333333)*(3.872983346207417*g1*(MDBS - Conj(MDBS))*
      KroneckerDelta(2,gO2)*(SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) + 2*SUM(j1,0,
      2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1))) - 15*(g2*(MDWBT - Conj(MDWBT))*
      KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) +
      1.4142135623730951*KroneckerDelta(1,gO2)*(-(Conj(Mu)*SUM(j2,0,2,Conj(ZD(gI2,
      j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1)))) + Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(
      Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZD(gI1,j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSeconjSe(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(
      3.872983346207417*g1*(MDBS - Conj(MDBS))*KroneckerDelta(2,gO2)*(SUM(j1,0,2,
      Conj(ZE(gI2,j1))*ZE(gI1,j1)) - 2*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3 +
      j1))) + 5*(g2*(MDWBT - Conj(MDWBT))*KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZE
      (gI2,j1))*ZE(gI1,j1)) + 1.4142135623730951*KroneckerDelta(1,gO2)*(-(Conj(Mu)
      *SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1)))) + Mu*SUM
      (j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZE(gI1,j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhSuconjSu(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,
      -0.03333333333333333)*(3.872983346207417*g1*(MDBS - Conj(MDBS))*
      KroneckerDelta(2,gO2)*(SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - 4*SUM(j1,0,
      2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1))) + 15*(g2*(MDWBT - Conj(MDWBT))*
      KroneckerDelta(3,gO2)*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) +
      1.4142135623730951*KroneckerDelta(0,gO2)*(Conj(Mu)*SUM(j2,0,2,Conj(ZU(gI2,j2
      ))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) - Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(
      j1,j2))*Conj(ZU(gI2,3 + j1)))*ZU(gI1,j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpUAhHpm(int gO2, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*Conj(Mu)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*ZP(gI2,1) - LamTD*
      KroneckerDelta(3,gO2)*ZP(gI2,1) + 1.4142135623730951*LamTD*KroneckerDelta(1,
      gO2)*ZP(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUAhHpmconjSRum(int gO2, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*Mu*(
      1.4142135623730951*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI2,0) + Conj(LamTU)
      *(KroneckerDelta(3,gO2)*ZP(gI2,0) + 1.4142135623730951*KroneckerDelta(0,gO2)
      *ZP(gI2,3)));

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
      1.4142135623730951*KroneckerDelta(3,gO2)*(ZP(gI2,2) - ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpURhconjSRdpconjURh(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) - 5*Sqr(g2)) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumURhconjSRumconjURh(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*Sqr(g2)) - KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjVWmVWm(int gO1, int gO2) const
{
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*Sqr(g2)
      ;

   return result;
}

double CLASSNAME::CpSRdpconjURhVWm(int gO2) const
{
   const double result = 0.7071067811865475*g2*KroneckerDelta(0,gO2);

   return result;
}

double CLASSNAME::CpSRumconjURhconjVWm(int gO2) const
{
   const double result = 0.7071067811865475*g2*KroneckerDelta(1,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpRhURhconjRhconjURh(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*(3*Sqr(g1) + 5*Sqr(g2))*(
      KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*ZHR(gI1,0)*ZHR(gI2,1) +
      KroneckerDelta(1,gO2)*(ZHR(gI1,0)*ZHR(gI2,0) - 2*ZHR(gI1,1)*ZHR(gI2,1))) +
      KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*ZHR(gI1,1)*ZHR(gI2,0) +
      KroneckerDelta(0,gO2)*(-2*ZHR(gI1,0)*ZHR(gI2,0) + ZHR(gI1,1)*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpCha1Cha2conjURhPR(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = Conj(LamTD)*KroneckerDelta(0,gO2)*UM1(
      gI2,1)*UP2(gI1,0) - Conj(LamTU)*KroneckerDelta(1,gO2)*UM1(gI2,0)*UP2(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpCha1Cha2conjURhPL(int gI2, int gI1, int gO1) const
{
   const std::complex<double> result = -(g2*(Conj(UM2(gI1,0))*Conj(UP1(gI2,1))*
      KroneckerDelta(0,gO1) + Conj(UM2(gI1,1))*Conj(UP1(gI2,0))*KroneckerDelta(1,
      gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpAhRhconjURh(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      -7.745966692414834*g1*MDBS*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR(gI1,0) +
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR(gI1,0
      ) + 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*
      ZHR(gI1,0) - 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(0,gO2)*
      ZA(gI2,2)*ZHR(gI1,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(0,gO2)
      *ZA(gI2,2)*ZHR(gI1,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(0,
      gO2)*ZA(gI2,2)*ZHR(gI1,0) + 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(
      gI1,0) - 7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(
      gI2,3)*ZHR(gI1,0) + 10*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(
      gI1,0) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(
      gI2,3)*ZHR(gI1,0) - 10*g2*Conj(MDWBT)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(
      gI1,0) - 10*LamTD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(gI1,0) - 10*
      LamSU*vu*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZHR(gI1,1) + 5*LamTU*vu
      *Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZHR(gI1,1) + 10*LamSU*vd*Conj(
      LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,1)*ZHR(gI1,1) - 5*LamTU*vd*Conj(LamTD)*
      KroneckerDelta(0,gO2)*ZA(gI2,1)*ZHR(gI1,1) + 7.745966692414834*g1*MDBS*
      KroneckerDelta(1,gO2)*ZA(gI2,2)*ZHR(gI1,1) - 7.745966692414834*g1*Conj(MDBS)
      *KroneckerDelta(1,gO2)*ZA(gI2,2)*ZHR(gI1,1) - 14.142135623730951*LamSU*Conj(
      MuU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZHR(gI1,1) - 10*g2*MDWBT*KroneckerDelta
      (1,gO2)*ZA(gI2,3)*ZHR(gI1,1) + 10*g2*Conj(MDWBT)*KroneckerDelta(1,gO2)*ZA(
      gI2,3)*ZHR(gI1,1) + 10*LamTU*Conj(MuU)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZHR(
      gI1,1) + 5*Conj(LamSU)*KroneckerDelta(1,gO2)*(2*LamSD*vu*ZA(gI2,0)*ZHR(gI1,0
      ) - 2*LamSD*vd*ZA(gI2,1)*ZHR(gI1,0) + 1.4142135623730951*((2*MuU - LamTU*vT)
      *ZA(gI2,2) + LamTU*vS*ZA(gI2,3))*ZHR(gI1,1)) - 5*Conj(LamTU)*KroneckerDelta(
      1,gO2)*(LamTD*vu*ZA(gI2,0)*ZHR(gI1,0) - LamTD*vd*ZA(gI2,1)*ZHR(gI1,0) + (
      -1.4142135623730951*LamSU*vT*ZA(gI2,2) + (2*MuU + 1.4142135623730951*LamSU*
      vS)*ZA(gI2,3))*ZHR(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CphhRhconjURh(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO2)*(-(((5*(
      2.8284271247461903*MuD + 4*LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(
      LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) + 2*(3.872983346207417*g1*
      MDBS + 3.872983346207417*g1*Conj(MDBS) + 7.0710678118654755*LamSD*Conj(MuD))
      )*ZH(gI2,2) + 5*(1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD +
      1.4142135623730951*LamSD*vS + 2*LamTD*vT)*Conj(LamTD) - 2*(g2*MDWBT + g2*
      Conj(MDWBT) - LamTD*Conj(MuD)))*ZH(gI2,3))*ZHR(gI1,0)) + ZH(gI2,0)*(vd*(-20*
      AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0) + 5*vu*
      (2*LamSU*Conj(LamSD) - LamTU*Conj(LamTD))*ZHR(gI1,1)) - ZH(gI2,1)*(vu*(3*Sqr
      (g1) + 5*Sqr(g2))*ZHR(gI1,0) + 5*vd*(-2*LamSU*Conj(LamSD) + LamTU*Conj(LamTD
      ))*ZHR(gI1,1))) + KroneckerDelta(1,gO2)*((-(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(
      gI2,0)) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) + 2*((3.872983346207417*g1*
      MDBS + 3.872983346207417*g1*Conj(MDBS) - 7.0710678118654755*LamSU*Conj(MuU))
      *ZH(gI2,2) - 5*(g2*MDWBT + g2*Conj(MDWBT) - LamTU*Conj(MuU))*ZH(gI2,3)))*ZHR
      (gI1,1) + 5*Conj(LamSU)*(2*LamSD*vu*ZH(gI2,0)*ZHR(gI1,0) + ((
      -2.8284271247461903*MuU - 4*LamSU*vS + 1.4142135623730951*LamTU*vT)*ZH(gI2,2
      ) + 1.4142135623730951*LamTU*vS*ZH(gI2,3))*ZHR(gI1,1) + 2*ZH(gI2,1)*(LamSD*
      vd*ZHR(gI1,0) - 2*LamSU*vu*ZHR(gI1,1))) - 5*Conj(LamTU)*(LamTD*vu*ZH(gI2,0)*
      ZHR(gI1,0) - (1.4142135623730951*LamSU*vT*ZH(gI2,2) + (2*MuU +
      1.4142135623730951*LamSU*vS - 2*LamTU*vT)*ZH(gI2,3))*ZHR(gI1,1) + ZH(gI2,1)*
      (LamTD*vd*ZHR(gI1,0) + 2*LamTU*vu*ZHR(gI1,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpURhSvconjURhconjSv(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjHpmconjURh(int gI1, int gO2) const
{
   const std::complex<double> result = 0.25*(-2.8284271247461903*LamSU*vd*Conj(
      LamSD)*KroneckerDelta(0,gO2)*ZP(gI1,1) - 1.4142135623730951*LamTU*Conj(LamTD
      )*KroneckerDelta(0,gO2)*(2*vu*ZP(gI1,0) + vd*ZP(gI1,1)) - KroneckerDelta(1,
      gO2)*(1.4142135623730951*vd*Sqr(g2)*ZP(gI1,0) + 1.4142135623730951*vu*(-2*
      AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,1) + 2*((-((2*MuU +
      1.4142135623730951*LamSU*vS + LamTU*vT)*Conj(LamTU)) + g2*(g2*vT + 2*Conj(
      MDWBT)))*ZP(gI1,2) - (-2*g2*MDWBT - vT*AbsSqr(LamTU) + 1.4142135623730951*
      LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2))*ZP(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhURhconjURh(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*(5*(2*LamSD*
      Conj(LamSU) - LamTD*Conj(LamTU))*KroneckerDelta(1,gO2)*(ZA(gI1,1)*ZA(gI2,0)
      + ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*
      Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) - 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(
      gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*
      ZA(gI2,3))) + Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA
      (gI1,3)*(1.4142135623730951*LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3)))))) +
      KroneckerDelta(1,gO1)*(5*(2*LamSU*Conj(LamSD) - LamTU*Conj(LamTD))*
      KroneckerDelta(0,gO2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) +
      KroneckerDelta(1,gO2)*(-((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0)) + (-20
      *AbsSqr(LamSU) - 10*AbsSqr(LamTU) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,
      1) + 5*(Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3
      )*(1.4142135623730951*LamSU*ZA(gI2,2) - 2*LamTU*ZA(gI2,3))) + Conj(LamSU)*(
      1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSU*ZA(gI2,2)
      + 1.4142135623730951*LamTU*ZA(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhURhconjURh(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*(5*(2*LamSD*
      Conj(LamSU) - LamTD*Conj(LamTU))*KroneckerDelta(1,gO2)*(ZH(gI1,1)*ZH(gI2,0)
      + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*
      Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) - 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZH(
      gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(4*LamSD*ZH(gI2,2) + 1.4142135623730951*LamTD*
      ZH(gI2,3))) + Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH
      (gI1,3)*(1.4142135623730951*LamSD*ZH(gI2,2) + 2*LamTD*ZH(gI2,3)))))) +
      KroneckerDelta(1,gO1)*(5*(2*LamSU*Conj(LamSD) - LamTU*Conj(LamTD))*
      KroneckerDelta(0,gO2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) +
      KroneckerDelta(1,gO2)*(-((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0)) + (-20
      *AbsSqr(LamSU) - 10*AbsSqr(LamTU) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,
      1) + 5*(Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3
      )*(1.4142135623730951*LamSU*ZH(gI2,2) - 2*LamTU*ZH(gI2,3))) + Conj(LamSU)*(
      1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSU*ZH(gI2,2)
      + 1.4142135623730951*LamTU*ZH(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpHpmURhconjHpmconjURh(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (20*
      AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*
      AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3))
      )) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-20*AbsSqr(LamTD) + 3*Sqr
      (g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*
      ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - (-2*AbsSqr(LamTD) + Sqr(g2))*
      ZP(gI1,3)*ZP(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhconjURh(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.25*Mu*(2*Conj(LamSU)*KroneckerDelta(1,
      gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) - 1.4142135623730951*Conj(
      LamTU)*KroneckerDelta(1,gO2)*(ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,3)) -
      KroneckerDelta(0,gO2)*(2*Conj(LamSD)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2
      ,2)) + 1.4142135623730951*Conj(LamTD)*(ZA(gI1,3)*ZA(gI2,1) + ZA(gI1,1)*ZA(
      gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhconjURh(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.25)*Mu*(2*Conj
      (LamSU)*KroneckerDelta(1,gO2)*(ZA(gI2,2)*ZH(gI1,0) - ZA(gI2,0)*ZH(gI1,2)) +
      1.4142135623730951*Conj(LamTU)*KroneckerDelta(1,gO2)*(-(ZA(gI2,3)*ZH(gI1,0))
      + ZA(gI2,0)*ZH(gI1,3)) + KroneckerDelta(0,gO2)*(Conj(LamSD)*(-2*ZA(gI2,2)*
      ZH(gI1,1) + 2*ZA(gI2,1)*ZH(gI1,2)) + 1.4142135623730951*Conj(LamTD)*(-(ZA(
      gI2,3)*ZH(gI1,1)) + ZA(gI2,1)*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhconjURh(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.25*Mu*(2*Conj(LamSU)*KroneckerDelta(1,
      gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) - 1.4142135623730951*Conj(
      LamTU)*KroneckerDelta(1,gO2)*(ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,3)) -
      KroneckerDelta(0,gO2)*(2*Conj(LamSD)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2
      ,2)) + 1.4142135623730951*Conj(LamTD)*(ZH(gI1,3)*ZH(gI2,1) + ZH(gI1,1)*ZH(
      gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmconjURh(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = -(Conj(LamTU)*KroneckerDelta(1,gO2)*Mu*
      ZP(gI1,2)*ZP(gI2,0)) + Conj(LamTD)*KroneckerDelta(0,gO2)*Mu*ZP(gI1,1)*ZP(gI2
      ,3);

   return result;
}

std::complex<double> CLASSNAME::CpChiChiconjURhPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = Conj(LamSD)*KroneckerDelta(0,gO2)*(ZN2(
      gI1,2)*ZN2(gI2,0) + ZN2(gI1,0)*ZN2(gI2,2)) + 0.5*(-2*Conj(LamSU)*
      KroneckerDelta(1,gO2)*(ZN2(gI1,3)*ZN2(gI2,0) + ZN2(gI1,0)*ZN2(gI2,3)) +
      1.4142135623730951*(Conj(LamTD)*KroneckerDelta(0,gO2)*(ZN2(gI1,2)*ZN2(gI2,1)
      + ZN2(gI1,1)*ZN2(gI2,2)) + Conj(LamTU)*KroneckerDelta(1,gO2)*(ZN2(gI1,3)*
      ZN2(gI2,1) + ZN2(gI1,1)*ZN2(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpChiChiconjURhPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = 0.1414213562373095*(Conj(ZN1(gI1,2))*(
      -3.872983346207417*g1*Conj(ZN1(gI2,0)) + 5*g2*Conj(ZN1(gI2,1)))*
      KroneckerDelta(0,gO1) + 5*g2*Conj(ZN1(gI1,1))*Conj(ZN1(gI2,2))*
      KroneckerDelta(0,gO1) + 3.872983346207417*g1*Conj(ZN1(gI1,3))*Conj(ZN1(gI2,0
      ))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN1(gI1,3))*Conj(ZN1(gI2,1))*
      KroneckerDelta(1,gO1) - 5*g2*Conj(ZN1(gI1,1))*Conj(ZN1(gI2,3))*
      KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN1(gI1,0))*(-(Conj(ZN1(
      gI2,2))*KroneckerDelta(0,gO1)) + Conj(ZN1(gI2,3))*KroneckerDelta(1,gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpURhSdconjURhconjSd(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((Sqr(
      g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpURhSeconjURhconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((3*Sqr
      (g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1
      ,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpURhSuconjURhconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = -0.05*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2))*((Sqr(
      g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,
      0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjURhconjSd(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) +
      vT*Conj(LamTD) + 2*Conj(MuD))*KroneckerDelta(0,gO2)*SUM(j2,0,2,Conj(ZD(gI2,
      j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjURhconjSe(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) +
      vT*Conj(LamTD) + 2*Conj(MuD))*KroneckerDelta(0,gO2)*SUM(j2,0,2,Conj(ZE(gI2,
      j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjURhconjSu(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = -0.5*(1.4142135623730951*vS*Conj(LamSU)
      - vT*Conj(LamTU) + 2*Conj(MuU))*KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZU(gI2
      ,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpRhconjURhVZ(int gI2, int gO2) const
{
   const std::complex<double> result = -0.1*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZHR(gI2,0) -
      KroneckerDelta(1,gO2)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpHpmconjURh(int gI2, int gO2) const
{
   const std::complex<double> result = 0.25*(-1.4142135623730951*KroneckerDelta
      (1,gO2)*(2*LamSD*vu*Conj(LamSU)*ZP(gI2,0) + LamTD*Conj(LamTU)*(vu*ZP(gI2,0)
      + 2*vd*ZP(gI2,1))) - KroneckerDelta(0,gO2)*(1.4142135623730951*vd*(-2*AbsSqr
      (LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZP(gI2,0) + 1.4142135623730951*vu*Sqr(g2)
      *ZP(gI2,1) + 4*g2*MDWBT*ZP(gI2,2) - 2*vT*AbsSqr(LamTD)*ZP(gI2,2) -
      2.8284271247461903*LamTD*vS*Conj(LamSD)*ZP(gI2,2) - 4*LamTD*Conj(MuD)*ZP(gI2
      ,2) + 2*vT*Sqr(g2)*ZP(gI2,2) + 2*vT*AbsSqr(LamTD)*ZP(gI2,3) - 4*MuD*Conj(
      LamTD)*ZP(gI2,3) - 2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gI2,3) + 4*g2*
      Conj(MDWBT)*ZP(gI2,3) - 2*vT*Sqr(g2)*ZP(gI2,3)));

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
   const std::complex<double> result = -0.1*g2*(3.872983346207417*g1*vd*Cos(
      ThetaW())*KroneckerDelta(0,gO2) - 3.872983346207417*g1*vu*Cos(ThetaW())*
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

std::complex<double> CLASSNAME::CpSRdpUHpmconjSRdpconjUHpm(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5*
      Sqr(g2)) + 10*(-(KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)) +
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(-2*AbsSqr(LamTD) + Sqr(g2))));

   return result;
}

std::complex<double> CLASSNAME::CpSRumUHpmconjSRumconjUHpm(int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2))) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-20*AbsSqr(LamSU) - 10*AbsSqr(LamTU) + 3*Sqr(g1) + 5*
      Sqr(g2)) - 10*(KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(g2) -
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*(-2*AbsSqr(LamTU) + Sqr(g2))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(int gO1, int gO2) const
{
   const std::complex<double> result = 0.1*(20*(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2) + KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2)
      *Sqr(Cos(ThetaW())) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(
      -7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW())
      + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(int gO1, int gO2) const
{
   const std::complex<double> result = 0.5*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2) + 2*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2) + 2*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjUHpmconjRh(int gO2, int gI1) const
{
   const std::complex<double> result = 0.25*(-2.8284271247461903*LamSU*vd*Conj(
      LamSD)*KroneckerDelta(1,gO2)*ZHR(gI1,0) - 1.4142135623730951*LamTU*Conj(
      LamTD)*(2*vu*KroneckerDelta(0,gO2) + vd*KroneckerDelta(1,gO2))*ZHR(gI1,0) -
      (1.4142135623730951*vd*KroneckerDelta(0,gO2)*Sqr(g2) + 1.4142135623730951*vu
      *KroneckerDelta(1,gO2)*(-2*AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2)) + 2*((-(
      (2*MuU + 1.4142135623730951*LamSU*vS + LamTU*vT)*Conj(LamTU)) + g2*(g2*vT +
      2*Conj(MDWBT)))*KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2)*(-2*g2*MDWBT -
      vT*AbsSqr(LamTU) + 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(
      MuU) + vT*Sqr(g2))))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmRhconjUHpmconjRh(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((-3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) + (-20
      *AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1)) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-20*AbsSqr(LamTD) + 3*Sqr(g1)
      - 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,1)*ZHR
      (gI2,1)) - 10*(KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*(-((-2*AbsSqr(
      LamTD) + Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0)) + Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1)) +
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0) -
      (-2*AbsSqr(LamTU) + Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjUHpmconjRh(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = -(Conj(LamTU)*KroneckerDelta(2,gO2)*Mu*
      ZHR(gI1,1)*ZP(gI2,0)) + Conj(LamTD)*KroneckerDelta(1,gO2)*Mu*ZHR(gI1,0)*ZP(
      gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1ChiconjUHpmPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -(Conj(LamSD)*KroneckerDelta(0,gO2)*UP1(
      gI1,1)*ZN2(gI2,0)) + 1.4142135623730951*g2*KroneckerDelta(3,gO2)*UP1(gI1,0)*
      ZN2(gI2,1) + 0.7071067811865475*Conj(LamTD)*KroneckerDelta(0,gO2)*UP1(gI1,1)
      *ZN2(gI2,1) - Conj(LamTD)*KroneckerDelta(2,gO2)*UP1(gI1,1)*ZN2(gI2,2) - g2*
      KroneckerDelta(1,gO2)*UP1(gI1,0)*ZN2(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1ChiconjUHpmPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = Conj(UM1(gI1,0))*(-(LamTU*Conj(ZN1(gI2,3
      ))*KroneckerDelta(1,gO1)) + 1.4142135623730951*g2*Conj(ZN1(gI2,1))*
      KroneckerDelta(2,gO1)) + Conj(UM1(gI1,1))*(0.5477225575051661*g1*Conj(ZN1(
      gI2,0))*KroneckerDelta(0,gO1) + 0.7071067811865475*g2*Conj(ZN1(gI2,1))*
      KroneckerDelta(0,gO1) + LamTD*Conj(ZN1(gI2,2))*KroneckerDelta(3,gO1));

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
   const std::complex<double> result = KroneckerDelta(0,gO2)*SUM(j2,0,2,SUM(j1,
      0,2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZUL(gI1,j2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjUHpmPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = KroneckerDelta(1,gO1)*SUM(j2,0,2,Conj(
      ZDL(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjUHpmPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = KroneckerDelta(0,gO2)*SUM(j1,0,2,Conj(Ye
      (j1,gI1))*ZER(gI2,j1));

   return result;
}

double CLASSNAME::CpbarFvFeconjUHpmPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSeconjUHpmconjSv(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = 0.25*(-2*g2*((g2*vT + 2*Conj(MDWBT))*
      KroneckerDelta(2,gO2) + (2*MDWBT - g2*vT)*KroneckerDelta(3,gO2))*SUM(j1,0,2,
      Conj(ZE(gI2,j1))*ZV(gI1,j1)) + KroneckerDelta(1,gO2)*(-1.4142135623730951*vu
      *Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZV(gI1,j1)) + 4*Mu*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1)))*ZV(gI1,j2))) + 1.4142135623730951*vd
      *KroneckerDelta(0,gO2)*(-(Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZV(gI1,j1))) +
      2*SUM(j3,0,2,SUM(j2,0,2,Conj(ZE(gI2,j2))*SUM(j1,0,2,Conj(Ye(j1,j3))*Ye(j1,
      j2)))*ZV(gI1,j3))));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(5*(-(KroneckerDelta(3,gO1)*(
      1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,0) + 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA(gI2,1) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,1) + 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2) +
      1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,3) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,3) + 4*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,3) + 2*KroneckerDelta(3,gO2)*
      (Sqr(g2)*ZA(gI1,0)*ZA(gI2,0) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZA(gI1,1)*ZA(gI2
      ,1) + 2*Sqr(g2)*ZA(gI1,3)*ZA(gI2,3)) - Conj(LamTD)*KroneckerDelta(0,gO2)*(2*
      LamSD*ZA(gI1,2)*ZA(gI2,0) + 1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,0) +
      ZA(gI1,0)*(2*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))))) +
      KroneckerDelta(2,gO1)*(1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2
      )*ZA(gI1,3)*ZA(gI2,0) - 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(
      gI1,3)*ZA(gI2,0) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(
      gI2,1) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,3)*ZA
      (gI2,1) + 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,
      1) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2) - 2*LamTD
      *Conj(LamSD)*KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2
      )) + 1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(gI2
      ,3) - 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,3) -
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,3)
      + 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,3) - 4*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,3) + 2*KroneckerDelta(2,gO2)*
      ((-2*AbsSqr(LamTD) + Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - Sqr(g2)*(ZA(gI1,1)*ZA(
      gI2,1) + 2*ZA(gI1,3)*ZA(gI2,3))))) + KroneckerDelta(0,gO1)*(KroneckerDelta(0
      ,gO2)*(-((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0)) + (3*Sqr(g1) - 5*Sqr(
      g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2
      )*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSD*ZA(gI2,2) - 2*LamTD*ZA(gI2
      ,3))) + Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2
      )*(-4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))))) + 5*(2*LamTD*
      Conj(LamSD)*KroneckerDelta(3,gO2)*ZA(gI1,2)*ZA(gI2,0) - 1.4142135623730951*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,0) - 1.4142135623730951*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,0) + KroneckerDelta(1,gO2)*
      Sqr(g2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + 2*LamTD*Conj(LamSD)*
      KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,2) - 1.4142135623730951*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,3) - 1.4142135623730951*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,3) + Conj(LamTD)*(
      1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*(ZA(gI1,3)*ZA(gI2,0) + ZA(gI1
      ,0)*ZA(gI2,3)) + KroneckerDelta(2,gO2)*(-2*LamSD*ZA(gI1,2)*ZA(gI2,0) +
      1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*(-2*LamSD*ZA(gI2,2)
      + 1.4142135623730951*LamTD*ZA(gI2,3)))))) + KroneckerDelta(1,gO1)*(5*(-2*
      LamTU*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA(gI1,2)*ZA(gI2,1) +
      1.4142135623730951*KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,1) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,1) +
      KroneckerDelta(0,gO2)*Sqr(g2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) -
      2*LamTU*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,2) +
      1.4142135623730951*KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,3) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,3) + Conj(
      LamTU)*(-1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*(ZA(gI1,3)*ZA(gI2,1)
      + ZA(gI1,1)*ZA(gI2,3)) + KroneckerDelta(2,gO2)*(2*LamSU*ZA(gI1,2)*ZA(gI2,1)
      - 1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,1) + ZA(gI1,1)*(2*LamSU*ZA(gI2,
      2) - 1.4142135623730951*LamTU*ZA(gI2,3))))) + KroneckerDelta(1,gO2)*((3*Sqr(
      g1) - 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(
      gI2,1) - 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(
      gI1,2)*(4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3))) + Conj(
      LamTU)*(1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(
      1.4142135623730951*LamSU*ZA(gI2,2) + 2*LamTU*ZA(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CphhhhUHpmconjUHpm(int gI1, int gI2, int gO1, int gO2) const
{
   const std::complex<double> result = 0.05*(5*(-(KroneckerDelta(3,gO1)*(
      -1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,0) + 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,1) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,3)*ZH(gI2,1) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,1) + 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,2) -
      1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,3) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,3) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,3) - 4*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,3) + 2*KroneckerDelta(3,gO2)*
      (Sqr(g2)*ZH(gI1,0)*ZH(gI2,0) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZH(gI1,1)*ZH(gI2
      ,1) + 2*Sqr(g2)*ZH(gI1,3)*ZH(gI2,3)) + Conj(LamTD)*KroneckerDelta(0,gO2)*(2*
      LamSD*ZH(gI1,2)*ZH(gI2,0) + 1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,0) +
      ZH(gI1,0)*(2*LamSD*ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3))))) +
      KroneckerDelta(2,gO1)*(1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2
      )*ZH(gI1,3)*ZH(gI2,0) - 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZH(
      gI1,3)*ZH(gI2,0) - 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(
      gI2,1) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,3)*ZH
      (gI2,1) - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,
      1) - 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,2) - 2*LamTD
      *Conj(LamSD)*KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2
      )) + 1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZH(gI1,0)*ZH(gI2
      ,3) - 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,3) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,3)
      - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,3) + 4*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,3) + 2*KroneckerDelta(2,gO2)*
      ((-2*AbsSqr(LamTD) + Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - Sqr(g2)*(ZH(gI1,1)*ZH(
      gI2,1) + 2*ZH(gI1,3)*ZH(gI2,3))))) - KroneckerDelta(0,gO1)*(5*(2*LamTD*Conj(
      LamSD)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(gI2,0) + 1.4142135623730951*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,0) - 1.4142135623730951*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,0) + KroneckerDelta(1,gO2)*
      Sqr(g2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + 2*LamTD*Conj(LamSD)*
      KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,2) + 1.4142135623730951*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,3) - 1.4142135623730951*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,3) + Conj(LamTD)*(
      1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*(ZH(gI1,3)*ZH(gI2,0) + ZH(gI1
      ,0)*ZH(gI2,3)) + KroneckerDelta(2,gO2)*(2*LamSD*ZH(gI1,2)*ZH(gI2,0) -
      1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*(2*LamSD*ZH(gI2,2)
      - 1.4142135623730951*LamTD*ZH(gI2,3))))) + KroneckerDelta(0,gO2)*((3*Sqr(g1)
      + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(
      gI2,1) + 5*(-(Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH
      (gI1,3)*(1.4142135623730951*LamSD*ZH(gI2,2) - 2*LamTD*ZH(gI2,3)))) - Conj(
      LamSD)*(1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSD*
      ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3)))))) + KroneckerDelta(1,gO1)*
      (-5*(2*LamTU*Conj(LamSU)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(gI2,1) +
      1.4142135623730951*KroneckerDelta(2,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,1) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,1) +
      KroneckerDelta(0,gO2)*Sqr(g2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) +
      2*LamTU*Conj(LamSU)*KroneckerDelta(3,gO2)*ZH(gI1,1)*ZH(gI2,2) +
      1.4142135623730951*KroneckerDelta(2,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,3) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,3) + Conj(
      LamTU)*(1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*(ZH(gI1,3)*ZH(gI2,1)
      + ZH(gI1,1)*ZH(gI2,3)) + KroneckerDelta(2,gO2)*(2*LamSU*ZH(gI1,2)*ZH(gI2,1)
      - 1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,1) + ZH(gI1,1)*(2*LamSU*ZH(gI2,2
      ) - 1.4142135623730951*LamTU*ZH(gI2,3))))) + KroneckerDelta(1,gO2)*((3*Sqr(
      g1) - 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(
      gI2,1) - 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(
      gI1,2)*(4*LamSU*ZH(gI2,2) + 1.4142135623730951*LamTU*ZH(gI2,3))) + Conj(
      LamTU)*(1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(
      1.4142135623730951*LamSU*ZH(gI2,2) + 2*LamTU*ZH(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpHpmUHpmconjHpmconjUHpm(int gI1, int gO1, int gI2, int gO2) const
{
   const std::complex<double> result = 0.05*(KroneckerDelta(1,gO1)*((
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0) + 10*(KroneckerDelta
      (2,gO2)*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,2) - KroneckerDelta(3,gO2)*Sqr(
      g2)*ZP(gI1,3)))*ZP(gI2,1) + KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*
      ZP(gI1,0)*ZP(gI2,0) - 2*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 5*(-2
      *AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2) + 5*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3
      )))) + KroneckerDelta(0,gO1)*((KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))
      *ZP(gI1,1) - 10*(KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2) - KroneckerDelta(3,
      gO2)*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,3)))*ZP(gI2,0) + KroneckerDelta(0,
      gO2)*(-2*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (3*Sqr(g1) + 5*Sqr(g2
      ))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - (-2*AbsSqr(LamTD)
      + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3)))) - 10*(KroneckerDelta(2,gO1)*((
      KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,0) - KroneckerDelta(1,gO2)*(-2*AbsSqr(
      LamTU) + Sqr(g2))*ZP(gI1,1) - 2*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3))*ZP(
      gI2,2) + KroneckerDelta(2,gO2)*(Sqr(g2)*ZP(gI1,0)*ZP(gI2,0) - (-2*AbsSqr(
      LamTU) + Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 2*Sqr(g2)*(2*ZP(gI1,2)*ZP(gI2,2) -
      ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(3,gO1)*((-(KroneckerDelta(0,gO2)*(-2
      *AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,0)) + Sqr(g2)*(KroneckerDelta(1,gO2)*ZP(gI1
      ,1) - 2*KroneckerDelta(2,gO2)*ZP(gI1,2)))*ZP(gI2,3) + KroneckerDelta(3,gO2)*
      (-((-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,0)*ZP(gI2,0)) + Sqr(g2)*(ZP(gI1,1)*ZP
      (gI2,1) - 2*ZP(gI1,2)*ZP(gI2,2) + 4*ZP(gI1,3)*ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiCha2conjUHpmPR(int gI1, int gI2, int gO2) const
{
   const std::complex<double> result = -1.4142135623730951*g2*KroneckerDelta(3,
      gO2)*UP2(gI2,0)*ZN1(gI1,1) - 0.1414213562373095*KroneckerDelta(1,gO2)*UP2(
      gI2,1)*(3.872983346207417*g1*ZN1(gI1,0) + 5*g2*ZN1(gI1,1)) + Conj(LamTD)*
      KroneckerDelta(0,gO2)*UP2(gI2,0)*ZN1(gI1,2) - Conj(LamTU)*KroneckerDelta(2,
      gO2)*UP2(gI2,1)*ZN1(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarChiCha2conjUHpmPL(int gI1, int gI2, int gO1) const
{
   const std::complex<double> result = -(g2*Conj(UM2(gI2,0))*(Conj(ZN2(gI1,2))*
      KroneckerDelta(0,gO1) + 1.4142135623730951*Conj(ZN2(gI1,1))*KroneckerDelta(2
      ,gO1))) + 0.5*Conj(UM2(gI2,1))*(2*LamSU*Conj(ZN2(gI1,0))*KroneckerDelta(1,
      gO1) + 1.4142135623730951*LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(1,gO1) + 2*
      LamTU*Conj(ZN2(gI1,3))*KroneckerDelta(3,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*(
      14.142135623730951*g2*MDWBT*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) +
      7.0710678118654755*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0
      ) + 10*LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) +
      14.142135623730951*LamTD*Conj(MuD)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0)
      - 7.0710678118654755*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,0) +
      7.745966692414834*g1*MDBS*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,0) +
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,0)
      - 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*
      ZP(gI1,0) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA
      (gI2,2)*ZP(gI1,0) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(0,gO2)*ZA
      (gI2,2)*ZP(gI1,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(0,gO2)
      *ZA(gI2,2)*ZP(gI1,0) + 10*LamTD*vd*Conj(LamSD)*KroneckerDelta(3,gO2)*ZA(gI2,
      2)*ZP(gI1,0) + 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,0) +
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(
      gI1,0) - 10*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,0) -
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(
      gI1,0) - 10*g2*Conj(MDWBT)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,0) + 10*
      LamTD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,0) +
      7.0710678118654755*vd*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZA(gI2,3)*ZP(gI1,0
      ) - 7.0710678118654755*vd*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,3)*ZP(gI1,0)
      - 5*vu*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,1) -
      14.142135623730951*g2*MDWBT*KroneckerDelta(3,gO2)*ZA(gI2,1)*ZP(gI1,1) -
      7.0710678118654755*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZA(gI2,1)*ZP(gI1,1
      ) - 10*LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA(gI2,1)*ZP(gI1,1) -
      14.142135623730951*LamTU*Conj(MuU)*KroneckerDelta(3,gO2)*ZA(gI2,1)*ZP(gI1,1)
      - 5*vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) +
      7.0710678118654755*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) + 10
      *LamTU*vu*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZP(gI1,1) +
      7.0710678118654755*vu*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZA(gI2,3)*ZP(gI1,1
      ) - 7.0710678118654755*vu*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,3)*ZP(gI1,1)
      - 14.142135623730951*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,2) +
      7.0710678118654755*vT*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,2
      ) - 10*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,2) -
      14.142135623730951*LamTD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,2)
      - 7.0710678118654755*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,2) +
      10*LamTD*vd*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,2) -
      7.0710678118654755*vd*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,2
      ) + 7.0710678118654755*vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,3)*ZP(gI1,2)
      - 20*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,3)*ZP(gI1,2) -
      7.0710678118654755*vT*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,3
      ) - 14.142135623730951*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(
      gI1,3) - 10*LamSD*vS*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,3) -
      14.142135623730951*g2*Conj(MDWBT)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,3)
      + 7.0710678118654755*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,3) -
      10*LamSD*vd*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,3) -
      7.0710678118654755*vd*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,3
      ) + 7.0710678118654755*vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,3)*ZP(gI1,3)
      + 5*KroneckerDelta(2,gO2)*(((2.8284271247461903*MuD + 2*LamSD*vS -
      1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(g2*vT + 2*
      Conj(MDWBT)))*ZA(gI2,0)*ZP(gI1,0) - 1.4142135623730951*vd*Sqr(g2)*ZA(gI2,3)*
      ZP(gI1,0) + vd*Conj(LamTD)*(-2*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA
      (gI2,3))*ZP(gI1,0) + 1.4142135623730951*vT*AbsSqr(LamTU)*ZA(gI2,1)*ZP(gI1,1)
      - 2.8284271247461903*MuU*Conj(LamTU)*ZA(gI2,1)*ZP(gI1,1) - 2*LamSU*vS*Conj(
      LamTU)*ZA(gI2,1)*ZP(gI1,1) - 2.8284271247461903*g2*Conj(MDWBT)*ZA(gI2,1)*ZP(
      gI1,1) - 1.4142135623730951*vT*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) - 2*LamSU*vu*Conj
      (LamTU)*ZA(gI2,2)*ZP(gI1,1) + 1.4142135623730951*vu*AbsSqr(LamTU)*ZA(gI2,3)*
      ZP(gI1,1) - 1.4142135623730951*vu*Sqr(g2)*ZA(gI2,3)*ZP(gI1,1) + 4*vT*Sqr(g2)
      *ZA(gI2,3)*ZP(gI1,3)) + KroneckerDelta(1,gO2)*(5*vu*Sqr(g2)*ZA(gI2,0)*ZP(gI1
      ,0) - 7.745966692414834*g1*MDBS*ZA(gI2,2)*ZP(gI1,1) + 14.142135623730951*MuU
      *Conj(LamSU)*ZA(gI2,2)*ZP(gI1,1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*
      ZA(gI2,2)*ZP(gI1,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZA(gI2,2)*ZP(
      gI1,1) + 7.745966692414834*g1*Conj(MDBS)*ZA(gI2,2)*ZP(gI1,1) -
      14.142135623730951*LamSU*Conj(MuU)*ZA(gI2,2)*ZP(gI1,1) - 10*g2*MDWBT*ZA(gI2,
      3)*ZP(gI1,1) - 7.0710678118654755*LamTU*vS*Conj(LamSU)*ZA(gI2,3)*ZP(gI1,1) +
      10*MuU*Conj(LamTU)*ZA(gI2,3)*ZP(gI1,1) + 7.0710678118654755*LamSU*vS*Conj(
      LamTU)*ZA(gI2,3)*ZP(gI1,1) + 10*g2*Conj(MDWBT)*ZA(gI2,3)*ZP(gI1,1) - 10*
      LamTU*Conj(MuU)*ZA(gI2,3)*ZP(gI1,1) + 10*LamTU*vu*Conj(LamSU)*ZA(gI2,2)*ZP(
      gI1,2) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZA(gI2,3)*ZP(gI1,2) +
      7.0710678118654755*vu*Sqr(g2)*ZA(gI2,3)*ZP(gI1,2) - 10*LamSU*vu*Conj(LamTU)*
      ZA(gI2,2)*ZP(gI1,3) - 7.0710678118654755*vu*AbsSqr(LamTU)*ZA(gI2,3)*ZP(gI1,3
      ) + 7.0710678118654755*vu*Sqr(g2)*ZA(gI2,3)*ZP(gI1,3) + 5*ZA(gI2,1)*(vd*Sqr(
      g2)*ZP(gI1,0) + (2*LamTU*vS*Conj(LamSU) + 1.4142135623730951*(2*g2*MDWBT -
      vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2)))*ZP(gI1,2) + ((
      2.8284271247461903*MuU + 2*LamSU*vS + 1.4142135623730951*LamTU*vT)*Conj(
      LamTU) + 1.4142135623730951*g2*(-(g2*vT) + 2*Conj(MDWBT)))*ZP(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjUHpm(int gI2, int gI1, int gO2) const
{
   const std::complex<double> result = 0.05*(-(KroneckerDelta(0,gO2)*(-3*vu*Sqr
      (g1)*ZH(gI2,1)*ZP(gI1,0) + 5*vu*Sqr(g2)*ZH(gI2,1)*ZP(gI1,0) -
      7.745966692414834*g1*MDBS*ZH(gI2,2)*ZP(gI1,0) + 20*vS*AbsSqr(LamSD)*ZH(gI2,2
      )*ZP(gI1,0) + 14.142135623730951*MuD*Conj(LamSD)*ZH(gI2,2)*ZP(gI1,0) -
      7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH(gI2,2)*ZP(gI1,0) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gI2,2)*ZP(gI1,0) -
      7.745966692414834*g1*Conj(MDBS)*ZH(gI2,2)*ZP(gI1,0) + 14.142135623730951*
      LamSD*Conj(MuD)*ZH(gI2,2)*ZP(gI1,0) - 10*g2*MDWBT*ZH(gI2,3)*ZP(gI1,0) + 10*
      vT*AbsSqr(LamTD)*ZH(gI2,3)*ZP(gI1,0) - 7.0710678118654755*LamTD*vS*Conj(
      LamSD)*ZH(gI2,3)*ZP(gI1,0) - 10*MuD*Conj(LamTD)*ZH(gI2,3)*ZP(gI1,0) -
      7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gI2,3)*ZP(gI1,0) - 10*g2*Conj(
      MDWBT)*ZH(gI2,3)*ZP(gI1,0) - 10*LamTD*Conj(MuD)*ZH(gI2,3)*ZP(gI1,0) + 5*vd*
      Sqr(g2)*ZH(gI2,1)*ZP(gI1,1) + 10*LamTD*vd*Conj(LamSD)*ZH(gI2,2)*ZP(gI1,2) -
      7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gI2,3)*ZP(gI1,2) + 7.0710678118654755
      *vd*Sqr(g2)*ZH(gI2,3)*ZP(gI1,2) + 10*LamSD*vd*Conj(LamTD)*ZH(gI2,2)*ZP(gI1,3
      ) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gI2,3)*ZP(gI1,3) -
      7.0710678118654755*vd*Sqr(g2)*ZH(gI2,3)*ZP(gI1,3) + ZH(gI2,0)*(vd*(3*Sqr(g1)
      + 5*Sqr(g2))*ZP(gI1,0) + 5*(vu*Sqr(g2)*ZP(gI1,1) + (2*LamTD*vS*Conj(LamSD)
      + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(MuD) + vT
      *Sqr(g2)))*ZP(gI1,2) + ((2.8284271247461903*MuD + 2*LamSD*vS +
      1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(-(g2*vT) +
      2*Conj(MDWBT)))*ZP(gI1,3))))) - KroneckerDelta(1,gO2)*(7.745966692414834*g1
      *MDBS*ZH(gI2,2)*ZP(gI1,1) + 20*vS*AbsSqr(LamSU)*ZH(gI2,2)*ZP(gI1,1) +
      14.142135623730951*MuU*Conj(LamSU)*ZH(gI2,2)*ZP(gI1,1) + 7.0710678118654755*
      LamTU*vT*Conj(LamSU)*ZH(gI2,2)*ZP(gI1,1) + 7.0710678118654755*LamSU*vT*Conj(
      LamTU)*ZH(gI2,2)*ZP(gI1,1) + 7.745966692414834*g1*Conj(MDBS)*ZH(gI2,2)*ZP(
      gI1,1) + 14.142135623730951*LamSU*Conj(MuU)*ZH(gI2,2)*ZP(gI1,1) + 10*g2*
      MDWBT*ZH(gI2,3)*ZP(gI1,1) + 10*vT*AbsSqr(LamTU)*ZH(gI2,3)*ZP(gI1,1) +
      7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gI2,3)*ZP(gI1,1) + 10*MuU*Conj(
      LamTU)*ZH(gI2,3)*ZP(gI1,1) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gI2,
      3)*ZP(gI1,1) + 10*g2*Conj(MDWBT)*ZH(gI2,3)*ZP(gI1,1) + 10*LamTU*Conj(MuU)*ZH
      (gI2,3)*ZP(gI1,1) + ZH(gI2,0)*(5*vu*Sqr(g2)*ZP(gI1,0) + vd*(-3*Sqr(g1) + 5*
      Sqr(g2))*ZP(gI1,1)) + 10*LamTU*vu*Conj(LamSU)*ZH(gI2,2)*ZP(gI1,2) -
      7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gI2,3)*ZP(gI1,2) + 7.0710678118654755
      *vu*Sqr(g2)*ZH(gI2,3)*ZP(gI1,2) + 10*LamSU*vu*Conj(LamTU)*ZH(gI2,2)*ZP(gI1,3
      ) + 7.0710678118654755*vu*AbsSqr(LamTU)*ZH(gI2,3)*ZP(gI1,3) -
      7.0710678118654755*vu*Sqr(g2)*ZH(gI2,3)*ZP(gI1,3) + ZH(gI2,1)*(5*vd*Sqr(g2)*
      ZP(gI1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1) + 5*((2*LamTU*vS*Conj(LamSU
      ) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTU) + 2*LamTU*Conj(MuU) +
      vT*Sqr(g2)))*ZP(gI1,2) + ((2.8284271247461903*MuU + 2*LamSU*vS +
      1.4142135623730951*LamTU*vT)*Conj(LamTU) + 1.4142135623730951*g2*(-(g2*vT) +
      2*Conj(MDWBT)))*ZP(gI1,3)))) - 5*(KroneckerDelta(2,gO2)*(1.4142135623730951
      *vd*Sqr(g2)*ZH(gI2,3)*ZP(gI1,0) + vd*Conj(LamTD)*(2*LamSD*ZH(gI2,2) -
      1.4142135623730951*LamTD*ZH(gI2,3))*ZP(gI1,0) - 1.4142135623730951*vT*AbsSqr
      (LamTU)*ZH(gI2,1)*ZP(gI1,1) + 2.8284271247461903*MuU*Conj(LamTU)*ZH(gI2,1)*
      ZP(gI1,1) + 2*LamSU*vS*Conj(LamTU)*ZH(gI2,1)*ZP(gI1,1) + 2.8284271247461903*
      g2*Conj(MDWBT)*ZH(gI2,1)*ZP(gI1,1) + 1.4142135623730951*vT*Sqr(g2)*ZH(gI2,1)
      *ZP(gI1,1) + 2*LamSU*vu*Conj(LamTU)*ZH(gI2,2)*ZP(gI1,1) - 1.4142135623730951
      *vu*AbsSqr(LamTU)*ZH(gI2,3)*ZP(gI1,1) + 1.4142135623730951*vu*Sqr(g2)*ZH(gI2
      ,3)*ZP(gI1,1) + 2*vu*Sqr(g2)*ZH(gI2,1)*ZP(gI1,2) + 4*vT*Sqr(g2)*ZH(gI2,3)*ZP
      (gI1,2) + ZH(gI2,0)*(((2.8284271247461903*MuD + 2*LamSD*vS -
      1.4142135623730951*LamTD*vT)*Conj(LamTD) + 1.4142135623730951*g2*(g2*vT + 2*
      Conj(MDWBT)))*ZP(gI1,0) - 2*vd*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,2)) - 4*
      vT*Sqr(g2)*ZH(gI2,3)*ZP(gI1,3)) + KroneckerDelta(3,gO2)*(2*LamTD*vd*Conj(
      LamSD)*ZH(gI2,2)*ZP(gI1,0) + 1.4142135623730951*vd*AbsSqr(LamTD)*ZH(gI2,3)*
      ZP(gI1,0) - 1.4142135623730951*vd*Sqr(g2)*ZH(gI2,3)*ZP(gI1,0) +
      2.8284271247461903*g2*MDWBT*ZH(gI2,1)*ZP(gI1,1) + 1.4142135623730951*vT*
      AbsSqr(LamTU)*ZH(gI2,1)*ZP(gI1,1) + 2*LamTU*vS*Conj(LamSU)*ZH(gI2,1)*ZP(gI1,
      1) + 2.8284271247461903*LamTU*Conj(MuU)*ZH(gI2,1)*ZP(gI1,1) -
      1.4142135623730951*vT*Sqr(g2)*ZH(gI2,1)*ZP(gI1,1) + 2*LamTU*vu*Conj(LamSU)*
      ZH(gI2,2)*ZP(gI1,1) + 1.4142135623730951*vu*AbsSqr(LamTU)*ZH(gI2,3)*ZP(gI1,1
      ) - 1.4142135623730951*vu*Sqr(g2)*ZH(gI2,3)*ZP(gI1,1) - 4*vT*Sqr(g2)*ZH(gI2,
      3)*ZP(gI1,2) + 4*vu*AbsSqr(LamTU)*ZH(gI2,1)*ZP(gI1,3) - 2*vu*Sqr(g2)*ZH(gI2,
      1)*ZP(gI1,3) + 4*vT*Sqr(g2)*ZH(gI2,3)*ZP(gI1,3) + ZH(gI2,0)*((2*LamTD*vS*
      Conj(LamSD) + 1.4142135623730951*(2*g2*MDWBT + vT*AbsSqr(LamTD) + 2*LamTD*
      Conj(MuD) - vT*Sqr(g2)))*ZP(gI1,0) + 2*vd*Sqr(g2)*ZP(gI1,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSdconjUHpmconjSd(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.05*(10*(-(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)) + KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2
      )*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(
      gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) - 20*
      SUM(j3,0,2,Conj(ZD(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yd(j2,
      j1))*ZD(gI2,3 + j2)))) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((Sqr(
      g1) - 5*Sqr(g2))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,
      0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)) + 20*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(
      gI1,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZD(gI2,j3))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSeconjUHpmconjSe(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.05*(10*(-(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)) + KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2
      )*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1))) -
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,
      0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))
      *ZE(gI2,3 + j1)) + 20*SUM(j3,0,2,Conj(ZE(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,
      Conj(Ye(j3,j1))*Ye(j2,j1))*ZE(gI2,3 + j2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmSuconjUHpmconjSu(int gO1, int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.05*(10*(KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2) - KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2)
      *SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(
      gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) + 20*
      SUM(j3,0,2,Conj(ZU(gI1,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j3,j1))*Yu(j2,
      j1))*ZU(gI2,3 + j2)))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((Sqr(
      g1) + 5*Sqr(g2))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*(Sqr(g1)*SUM(j1
      ,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)) + 5*SUM(j3,0,2,SUM(j2,0,2,Conj(ZU(
      gI1,j2))*SUM(j1,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gI2,j3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjUHpmconjSu(int gI2, int gO2, int gI1) const
{
   const std::complex<double> result = 0.25*(-2*g2*((g2*vT + 2*Conj(MDWBT))*
      KroneckerDelta(2,gO2) + (2*MDWBT - g2*vT)*KroneckerDelta(3,gO2))*SUM(j1,0,2,
      Conj(ZD(gI2,j1))*ZU(gI1,j1)) + KroneckerDelta(0,gO2)*(-1.4142135623730951*vd
      *Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 2*(2*Conj(Mu)*SUM(j2,0,2,
      Conj(ZD(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))) + 1.4142135623730951*
      (vu*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j3,j1))*Yu
      (j2,j1))*ZU(gI1,3 + j2))) + vd*SUM(j3,0,2,SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1
      ,0,2,Conj(Yd(j1,j3))*Yd(j1,j2)))*ZU(gI1,j3))))) + KroneckerDelta(1,gO2)*(
      -1.4142135623730951*vu*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZU(gI1,j1)) + 2*(
      2*Mu*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1)))*ZU(gI1,j2))
      + 1.4142135623730951*(vd*SUM(j3,0,2,Conj(ZD(gI2,3 + j3))*SUM(j2,0,2,SUM(j1,
      0,2,Conj(Yd(j3,j1))*Yu(j2,j1))*ZU(gI1,3 + j2))) + vu*SUM(j3,0,2,SUM(j2,0,2,
      Conj(ZD(gI2,j2))*SUM(j1,0,2,Conj(Yu(j1,j3))*Yu(j1,j2)))*ZU(gI1,j3))))));

   return result;
}

std::complex<double> CLASSNAME::CpRhconjSRdpconjUHpm(int gI2, int gO2) const
{
   const std::complex<double> result = 0.25*(2*(((2*MuD + 1.4142135623730951*
      LamSD*vS + LamTD*vT)*Conj(LamTD) - g2*(g2*vT + 2*Conj(MDWBT)))*
      KroneckerDelta(2,gO2) + KroneckerDelta(3,gO2)*(-2*g2*MDWBT - vT*AbsSqr(LamTD
      ) + 1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)
      ))*ZHR(gI2,0) - 1.4142135623730951*KroneckerDelta(1,gO2)*(vu*Sqr(g2)*ZHR(gI2
      ,0) + 2*LamTU*vd*Conj(LamTD)*ZHR(gI2,1)) - 1.4142135623730951*KroneckerDelta
      (0,gO2)*(vd*(-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZHR(gI2,0) + vu*(2*
      LamSU*Conj(LamSD) + LamTU*Conj(LamTD))*ZHR(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumAhconjUHpm(int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*Conj(Mu)*(
      1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*ZA(gI2,0) + KroneckerDelta(0,
      gO2)*(1.4142135623730951*LamSU*ZA(gI2,2) + LamTU*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumhhconjUHpm(int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*Conj(Mu)*(-1.4142135623730951*LamTU*
      KroneckerDelta(3,gO2)*ZH(gI2,0) + KroneckerDelta(0,gO2)*(1.4142135623730951*
      LamSU*ZH(gI2,2) + LamTU*ZH(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjSRdpconjUHpm(int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*Mu*(
      1.4142135623730951*Conj(LamSD)*KroneckerDelta(1,gO2)*ZA(gI2,2) + Conj(LamTD)
      *(1.4142135623730951*KroneckerDelta(2,gO2)*ZA(gI2,1) - KroneckerDelta(1,gO2)
      *ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CphhconjSRdpconjUHpm(int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*Mu*(-1.4142135623730951*Conj(LamSD)*
      KroneckerDelta(1,gO2)*ZH(gI2,2) + Conj(LamTD)*(1.4142135623730951*
      KroneckerDelta(2,gO2)*ZH(gI2,1) + KroneckerDelta(1,gO2)*ZH(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjUHpmVWm(int gI2, int gO2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(
      KroneckerDelta(0,gO2)*ZA(gI2,0) + KroneckerDelta(1,gO2)*ZA(gI2,1) +
      1.4142135623730951*(KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2))*ZA(gI2,3)
      );

   return result;
}

std::complex<double> CLASSNAME::CphhconjUHpmVWm(int gI2, int gO2) const
{
   const std::complex<double> result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0)
      - KroneckerDelta(1,gO2)*ZH(gI2,1) + 1.4142135623730951*(KroneckerDelta(2,gO2
      ) + KroneckerDelta(3,gO2))*ZH(gI2,3));

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
      10*g2*Cos(ThetaW())*(KroneckerDelta(2,gO2)*ZP(gI2,2) + KroneckerDelta(3,gO2
      )*ZP(gI2,3)));

   return result;
}

double CLASSNAME::CpSRdpSRdpconjSRdpconjSRdp() const
{
   const double result = 0.1*(-3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

double CLASSNAME::CpSRdpSRumconjSRdpconjSRum() const
{
   const double result = 0.25*(0.6*Sqr(g1) + Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpconjSRdpVZVZ() const
{
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW
      ())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW(
      ))));

   return result;
}

double CLASSNAME::CpSRdpconjSRdpconjVWmVWm() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpSRdpconjSRdpVP() const
{
   const double result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpSRdpconjSRdpVZ() const
{
   const double result = 0.5*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpRhconjSRdpconjRh(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-((3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,
      0)*ZHR(gI2,0)) + (3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

double CLASSNAME::CpSRdpSvconjSRdpconjSv(int gI1, int gI2) const
{
   const double result = 0.05*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpAhAhconjSRdp(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5
      *Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1)
      + 5*(Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*
      (1.4142135623730951*LamSD*ZA(gI2,2) - 2*LamTD*ZA(gI2,3))) + Conj(LamSD)*(
      1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSD*ZA(gI2,2)
      + 1.4142135623730951*LamTD*ZA(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRdphhhhconjSRdp(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5
      *Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1)
      + 5*(Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*
      (1.4142135623730951*LamSD*ZH(gI2,2) - 2*LamTD*ZH(gI2,3))) + Conj(LamSD)*(
      1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSD*ZH(gI2,2)
      + 1.4142135623730951*LamTD*ZH(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpHpmconjSRdpconjHpm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.25*((-4*AbsSqr(LamSD) - 2*AbsSqr(LamTD
      ) + 0.6*Sqr(g1) + Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - (0.6*Sqr(g1) + Sqr(g2))*ZP(
      gI1,1)*ZP(gI2,1) - 4*AbsSqr(LamTD)*ZP(gI1,2)*ZP(gI2,2) + 2*Sqr(g2)*ZP(gI1,2)
      *ZP(gI2,2) - 2*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpRhconjSRdpconjHpm(int gI2, int gI1) const
{
   const std::complex<double> result = 0.25*(-1.4142135623730951*ZHR(gI2,1)*(2*
      LamSU*vu*Conj(LamSD)*ZP(gI1,0) + LamTU*Conj(LamTD)*(vu*ZP(gI1,0) + 2*vd*ZP(
      gI1,1))) - ZHR(gI2,0)*(1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(
      LamTD) + Sqr(g2))*ZP(gI1,0) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1) - 2*vT
      *AbsSqr(LamTD)*ZP(gI1,2) - 4*MuD*Conj(LamTD)*ZP(gI1,2) - 2.8284271247461903*
      LamSD*vS*Conj(LamTD)*ZP(gI1,2) + 4*g2*Conj(MDWBT)*ZP(gI1,2) + 2*vT*Sqr(g2)*
      ZP(gI1,2) + 4*g2*MDWBT*ZP(gI1,3) + 2*vT*AbsSqr(LamTD)*ZP(gI1,3) -
      2.8284271247461903*LamTD*vS*Conj(LamSD)*ZP(gI1,3) - 4*LamTD*Conj(MuD)*ZP(gI1
      ,3) - 2*vT*Sqr(g2)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpCha1ChiconjSRdpPR(int gI2, int gI1) const
{
   const std::complex<double> result = -(Conj(LamSD)*UM1(gI2,1)*ZN2(gI1,0)) +
      Conj(LamTD)*(0.7071067811865475*UM1(gI2,1)*ZN2(gI1,1) - UM1(gI2,0)*ZN2(gI1,2
      ));

   return result;
}

std::complex<double> CLASSNAME::CpCha1ChiconjSRdpPL(int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*Conj(UP1(gI2,1))*(
      0.7745966692414834*g1*Conj(ZN1(gI1,0)) + g2*Conj(ZN1(gI1,1))) - g2*Conj(UP1(
      gI2,0))*Conj(ZN1(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpAhconjSRdpconjHpm(int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*Mu*(
      1.4142135623730951*Conj(LamSD)*ZA(gI2,2)*ZP(gI1,1) + Conj(LamTD)*(-(ZA(gI2,3
      )*ZP(gI1,1)) + 1.4142135623730951*ZA(gI2,1)*ZP(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CphhconjSRdpconjHpm(int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*Mu*(-1.4142135623730951*Conj(LamSD)*
      ZH(gI2,2)*ZP(gI1,1) + Conj(LamTD)*(ZH(gI2,3)*ZP(gI1,1) + 1.4142135623730951*
      ZH(gI2,1)*ZP(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpSdconjSRdpconjSd(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2
      ,Conj(ZD(gI1,j1))*ZD(gI2,j1))) - 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*
      ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpSeconjSRdpconjSe(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2
      ,Conj(ZE(gI1,j1))*ZE(gI2,j1)) - 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE
      (gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSRdpSuconjSRdpconjSu(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2
      ,Conj(ZU(gI1,j1))*ZU(gI2,j1))) + 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*
      ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSRdpconjSe(int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) -
      vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZV(gI2,j2))*SUM(j1,0,2,Ye(j1,
      j2)*ZE(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSRdpconjSd(int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*vS*Conj(LamSD) -
      vT*Conj(LamTD) + 2*Conj(MuD))*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Yd(j1,
      j2)*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpRhconjSRdpconjVWm(int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*g2*ZHR(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpSRdpAhconjSRdp(int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.25)*((
      1.5491933384829668*g1*MDBS + 1.4142135623730951*(-2*MuD + LamTD*vT)*Conj(
      LamSD) - 1.4142135623730951*LamSD*vT*Conj(LamTD) - 1.5491933384829668*g1*
      Conj(MDBS) + 2.8284271247461903*LamSD*Conj(MuD))*ZA(gI2,2) + (2*g2*MDWBT -
      1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD + 1.4142135623730951*LamSD*
      vS)*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(MuD))*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpSRdphhconjSRdp(int gI2) const
{
   const std::complex<double> result = 0.25*(vd*(-4*AbsSqr(LamTD) + 0.6*Sqr(g1)
      - Sqr(g2))*ZH(gI2,0) + vu*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gI2,1) -
      1.5491933384829668*g1*MDBS*ZH(gI2,2) - 4*vS*AbsSqr(LamSD)*ZH(gI2,2) -
      2.8284271247461903*MuD*Conj(LamSD)*ZH(gI2,2) + 1.4142135623730951*LamTD*vT*
      Conj(LamSD)*ZH(gI2,2) + 1.4142135623730951*LamSD*vT*Conj(LamTD)*ZH(gI2,2) -
      1.5491933384829668*g1*Conj(MDBS)*ZH(gI2,2) - 2.8284271247461903*LamSD*Conj(
      MuD)*ZH(gI2,2) - 2*g2*MDWBT*ZH(gI2,3) - 2*vT*AbsSqr(LamTD)*ZH(gI2,3) +
      1.4142135623730951*LamTD*vS*Conj(LamSD)*ZH(gI2,3) + 2*MuD*Conj(LamTD)*ZH(gI2
      ,3) + 1.4142135623730951*LamSD*vS*Conj(LamTD)*ZH(gI2,3) - 2*g2*Conj(MDWBT)*
      ZH(gI2,3) + 2*LamTD*Conj(MuD)*ZH(gI2,3));

   return result;
}

double CLASSNAME::CpSRumSRumconjSRumconjSRum() const
{
   const double result = 0.1*(-3*Sqr(g1) - 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumVZVZ() const
{
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW
      ())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW(
      ))));

   return result;
}

double CLASSNAME::CpSRumconjSRumconjVWmVWm() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

double CLASSNAME::CpSRumconjSRumVP() const
{
   const double result = 0.1*(-3.872983346207417*g1*Cos(ThetaW()) - 5*g2*Sin(
      ThetaW()));

   return result;
}

double CLASSNAME::CpSRumconjSRumVZ() const
{
   const double result = 0.1*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(
      ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpSRumRhconjSRumconjRh(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,0)
      *ZHR(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmRhconjSRum(int gI2, int gI1) const
{
   const std::complex<double> result = 0.25*(-2.8284271247461903*LamSD*vd*Conj(
      LamSU)*ZHR(gI1,0)*ZP(gI2,1) - ZHR(gI1,1)*(1.4142135623730951*vd*Sqr(g2)*ZP(
      gI2,0) + 1.4142135623730951*vu*(-2*AbsSqr(LamSU) + Sqr(g2))*ZP(gI2,1) + 4*g2
      *MDWBT*ZP(gI2,2) - 2.8284271247461903*LamTU*vS*Conj(LamSU)*ZP(gI2,2) - 4*
      LamTU*Conj(MuU)*ZP(gI2,2) + 2*vT*Sqr(g2)*ZP(gI2,2) + 4*g2*Conj(MDWBT)*ZP(gI2
      ,3) - 2*vT*Sqr(g2)*ZP(gI2,3)) - Conj(LamTU)*(1.4142135623730951*LamTD*ZHR(
      gI1,0)*(2*vu*ZP(gI2,0) + vd*ZP(gI2,1)) + ZHR(gI1,1)*(1.4142135623730951*
      LamTU*vu*ZP(gI2,1) - 2*(LamTU*vT*ZP(gI2,2) + (2*MuU + 1.4142135623730951*
      LamSU*vS - LamTU*vT)*ZP(gI2,3)))));

   return result;
}

double CLASSNAME::CpSRumSvconjSRumconjSv(int gI1, int gI2) const
{
   const double result = 0.25*KroneckerDelta(gI1,gI2)*(-0.6*Sqr(g1) + Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpSRumAhAhconjSRum(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((-3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)
      *ZA(gI2,0) + (-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1)
      - 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*
      (4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3))) + Conj(LamTU)*(
      1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSU*ZA(gI2,2) + 2*LamTU*ZA(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRumhhhhconjSRum(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((-3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)
      *ZH(gI2,0) + (-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1)
      - 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*
      (4*LamSU*ZH(gI2,2) + 1.4142135623730951*LamTU*ZH(gI2,3))) + Conj(LamTU)*(
      1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951
      *LamSU*ZH(gI2,2) + 2*LamTU*ZH(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpSRumHpmconjSRumconjHpm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.25*(-((0.6*Sqr(g1) + Sqr(g2))*ZP(gI1,0
      )*ZP(gI2,0)) + (-4*AbsSqr(LamSU) - 2*AbsSqr(LamTU) + 0.6*Sqr(g1) + Sqr(g2))*
      ZP(gI1,1)*ZP(gI2,1) - 2*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 4*AbsSqr(LamTU)*ZP(gI1
      ,3)*ZP(gI2,3) + 2*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpCha2ChiconjSRumPR(int gI2, int gI1) const
{
   const std::complex<double> result = Conj(LamSU)*UP2(gI2,1)*ZN2(gI1,0) + Conj
      (LamTU)*(0.7071067811865475*UP2(gI2,1)*ZN2(gI1,1) + UP2(gI2,0)*ZN2(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpCha2ChiconjSRumPL(int gI2, int gI1) const
{
   const std::complex<double> result = 0.7071067811865475*Conj(UM2(gI2,1))*(
      0.7745966692414834*g1*Conj(ZN1(gI1,0)) + g2*Conj(ZN1(gI1,1))) - g2*Conj(UM2(
      gI2,0))*Conj(ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjSRum(int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*Mu*(
      1.4142135623730951*Conj(LamSU)*ZA(gI2,2)*ZP(gI1,0) + Conj(LamTU)*(ZA(gI2,3)*
      ZP(gI1,0) + 1.4142135623730951*ZA(gI2,0)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjSRum(int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*Mu*(1.4142135623730951*Conj(LamSU)*
      ZH(gI2,2)*ZP(gI1,0) + Conj(LamTU)*(ZH(gI2,3)*ZP(gI1,0) - 1.4142135623730951*
      ZH(gI2,0)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumSdconjSRumconjSd(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((Sqr(g1) - 5*Sqr(g2))*SUM(j1,0,2,
      Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 2*Sqr(g1)*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(
      gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumSeconjSRumconjSe(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(-((3*Sqr(g1) + 5*Sqr(g2))*SUM(j1,0
      ,2,Conj(ZE(gI1,j1))*ZE(gI2,j1))) + 6*Sqr(g1)*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))
      *ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSRumSuconjSRumconjSu(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*((Sqr(g1) + 5*Sqr(g2))*SUM(j1,0,2,
      Conj(ZU(gI1,j1))*ZU(gI2,j1)) - 4*Sqr(g1)*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(
      gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSRumconjSu(int gI2, int gI1) const
{
   const std::complex<double> result = -0.5*(1.4142135623730951*vS*Conj(LamSU)
      + vT*Conj(LamTU) + 2*Conj(MuU))*SUM(j2,0,2,Conj(ZD(gI2,j2))*SUM(j1,0,2,Yu(j1
      ,j2)*ZU(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpRhconjSRumVWm(int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*g2*ZHR(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpSRumAhconjSRum(int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,0.05)*((
      7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU
      ) - 7.0710678118654755*LamSU*vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS
      ) - 14.142135623730951*LamSU*Conj(MuU))*ZA(gI2,2) + 5*(2*g2*MDWBT -
      1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) +
      1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(
      MuU))*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpSRumhhconjSRum(int gI2) const
{
   const std::complex<double> result = 0.25*(vd*(-0.6*Sqr(g1) + Sqr(g2))*ZH(gI2
      ,0) + vu*(-4*AbsSqr(LamTU) + 0.6*Sqr(g1) - Sqr(g2))*ZH(gI2,1) +
      1.5491933384829668*g1*MDBS*ZH(gI2,2) - 4*vS*AbsSqr(LamSU)*ZH(gI2,2) -
      2.8284271247461903*MuU*Conj(LamSU)*ZH(gI2,2) - 1.4142135623730951*LamTU*vT*
      Conj(LamSU)*ZH(gI2,2) - 1.4142135623730951*LamSU*vT*Conj(LamTU)*ZH(gI2,2) +
      1.5491933384829668*g1*Conj(MDBS)*ZH(gI2,2) - 2.8284271247461903*LamSU*Conj(
      MuU)*ZH(gI2,2) + 2*g2*MDWBT*ZH(gI2,3) - 2*vT*AbsSqr(LamTU)*ZH(gI2,3) -
      1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gI2,3) - 2*MuU*Conj(LamTU)*ZH(gI2
      ,3) - 1.4142135623730951*LamSU*vS*Conj(LamTU)*ZH(gI2,3) + 2*g2*Conj(MDWBT)*
      ZH(gI2,3) - 2*LamTU*Conj(MuU)*ZH(gI2,3));

   return result;
}

double CLASSNAME::CpsigmaOsigmaOphiOphiO() const
{
   const double result = -6*Sqr(g3);

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOsigmaOVG() const
{
   const std::complex<double> result = std::complex<double>(0,1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOSdconjSd(int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,-1)*g3*(MDGoc -
      Conj(MDGoc))*(SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,j1)) - SUM(j1,0,2,Conj(ZD(
      gI2,3 + j1))*ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpsigmaOSuconjSu(int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,-1)*g3*(MDGoc -
      Conj(MDGoc))*(SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,j1)) - SUM(j1,0,2,Conj(ZU(
      gI2,3 + j1))*ZU(gI1,3 + j1)));

   return result;
}

double CLASSNAME::CpbarGluGlusigmaOPR() const
{
   const double result = -g3;

   return result;
}

double CLASSNAME::CpbarGluGlusigmaOPL() const
{
   const double result = g3;

   return result;
}

double CLASSNAME::CpphiOphiOsigmaOsigmaO() const
{
   const double result = -6*Sqr(g3);

   return result;
}

std::complex<double> CLASSNAME::CpphiOphiOVG() const
{
   const std::complex<double> result = std::complex<double>(0,1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpphiOSdconjSd(int gI2, int gI1) const
{
   const std::complex<double> result = -(g3*(MDGoc + Conj(MDGoc))*(SUM(j1,0,2,
      Conj(ZD(gI2,j1))*ZD(gI1,j1)) - SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpphiOSuconjSu(int gI2, int gI1) const
{
   const std::complex<double> result = -(g3*(MDGoc + Conj(MDGoc))*(SUM(j1,0,2,
      Conj(ZU(gI2,j1))*ZU(gI1,j1)) - SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpbarGluGluphiOPR() const
{
   const std::complex<double> result = std::complex<double>(0,1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluGluphiOPL() const
{
   const std::complex<double> result = std::complex<double>(0,1)*g3;

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

double CLASSNAME::CpphiOphiOVGVG() const
{
   const double result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpsigmaOsigmaOVGVG() const
{
   const double result = 16*Sqr(g3);

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

std::complex<double> CLASSNAME::CpbarGluGluVGPL() const
{
   const std::complex<double> result = std::complex<double>(0,1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluGluVGPR() const
{
   const std::complex<double> result = std::complex<double>(0,1)*g3;

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

std::complex<double> CLASSNAME::CpSRdpconjSRdpVPVP() const
{
   const std::complex<double> result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834
      *g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpSRumconjSRumVPVP() const
{
   const std::complex<double> result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834
      *g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpconjVWmVPVWm() const
{
   const double result = g2*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1VPPL(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(-2*g2*Conj(UP1(gI2,0))*Sin(ThetaW()
      )*UP1(gI1,0) - Conj(UP1(gI2,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*
      Sin(ThetaW()))*UP1(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1VPPR(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(-2*g2*Conj(UM1(gI1,0))*Sin(ThetaW()
      )*UM1(gI2,0) - Conj(UM1(gI1,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*
      Sin(ThetaW()))*UM1(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2VPPL(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(2*g2*Conj(UM2(gI2,0))*Sin(ThetaW())
      *UM2(gI1,0) + Conj(UM2(gI2,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin
      (ThetaW()))*UM2(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2VPPR(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(2*g2*Conj(UP2(gI1,0))*Sin(ThetaW())
      *UP2(gI2,0) + Conj(UP2(gI1,1))*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin
      (ThetaW()))*UP2(gI2,1));

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
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1
      ,gI2);

   return result;
}

double CLASSNAME::CpbarFeFeVPPL(int gI1, int gI2) const
{
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos
      (ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFeFeVPPR(int gI1, int gI2) const
{
   const double result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1
      ,gI2);

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
   const double result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(
      gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVPVP(int gI1, int gI2) const
{
   const std::complex<double> result = 0.1*((g2*Sin(ThetaW())*(
      7.745966692414834*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos
      (ThetaW())))*ZP(gI1,0)*ZP(gI2,0) + (g2*Sin(ThetaW())*(7.745966692414834*g1*
      Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())))*ZP(gI1,1
      )*ZP(gI2,1) + 20*Sqr(g2)*Sqr(Sin(ThetaW()))*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)
      *ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVP(int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(-((0.7745966692414834*g1*Cos(ThetaW
      ()) + g2*Sin(ThetaW()))*ZP(gI1,0)*ZP(gI2,0)) - (0.7745966692414834*g1*Cos(
      ThetaW()) + g2*Sin(ThetaW()))*ZP(gI1,1)*ZP(gI2,1) - 2*g2*Sin(ThetaW())*(ZP(
      gI1,2)*ZP(gI2,2) + ZP(gI1,3)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVPVP(int gI1, int gI2) const
{
   const std::complex<double> result = 0.03333333333333333*((-7.745966692414834
      *g1*g2*Cos(ThetaW())*Sin(ThetaW()) + Sqr(g1)*Sqr(Cos(ThetaW())) + 15*Sqr(g2)
      *Sqr(Sin(ThetaW())))*SUM(j1,0,2,Conj(ZD(gI1,j1))*ZD(gI2,j1)) + 4*Sqr(g1)*Sqr
      (Cos(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*ZD(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVPVP(int gI1, int gI2) const
{
   const std::complex<double> result = 0.1*((g2*Sin(ThetaW())*(
      7.745966692414834*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos
      (ThetaW())))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 12*Sqr(g1)*Sqr(Cos(
      ThetaW()))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

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
   const std::complex<double> result = 0.16666666666666666*((0.7745966692414834
      *g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1,
      j1)) - 1.5491933384829668*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))*
      ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVP(int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(-((0.7745966692414834*g1*Cos(ThetaW
      ()) + g2*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1))) -
      1.5491933384829668*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,3
      + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVP(int gI2, int gI1) const
{
   const std::complex<double> result = 0.16666666666666666*((0.7745966692414834
      *g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,
      j1)) + 3.0983866769659336*g1*Cos(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      ZU(gI1,3 + j1)));

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

std::complex<double> CLASSNAME::CpRhconjRhVZVZ(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*
      ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2))
      )*(ZHR(gI1,0)*ZHR(gI2,0) + ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpRhconjRhVZ(int gI2, int gI1) const
{
   const std::complex<double> result = -0.5*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*(ZHR(gI1,0)*ZHR(gI2,0) - ZHR(gI1,1)*ZHR
      (gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1VZPL(int gI1, int gI2) const
{
   const std::complex<double> result = 0.1*(-10*g2*Conj(UP1(gI2,0))*Cos(ThetaW(
      ))*UP1(gI1,0) + Conj(UP1(gI2,1))*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1
      *Sin(ThetaW()))*UP1(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1Cha1VZPR(int gI1, int gI2) const
{
   const std::complex<double> result = 0.1*(-10*g2*Conj(UM1(gI1,0))*Cos(ThetaW(
      ))*UM1(gI2,0) + Conj(UM1(gI1,1))*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1
      *Sin(ThetaW()))*UM1(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2VZPL(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(2*g2*Conj(UM2(gI2,0))*Cos(ThetaW())
      *UM2(gI1,0) + Conj(UM2(gI2,1))*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin
      (ThetaW()))*UM2(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2Cha2VZPR(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(2*g2*Conj(UP2(gI1,0))*Cos(ThetaW())
      *UP2(gI2,0) + Conj(UP2(gI1,1))*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin
      (ThetaW()))*UP2(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpSvconjSvVZVZ(int gI1, int gI2) const
{
   const std::complex<double> result = 0.1*KroneckerDelta(gI1,gI2)*(g1*Sin(
      ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpSvconjSvVZ(int gI2, int gI1) const
{
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

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
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

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
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double CLASSNAME::CpbarFuFuVZPL(int gI1, int gI2) const
{
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*
      Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuFuVZPR(int gI1, int gI2) const
{
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

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

std::complex<double> CLASSNAME::CpAhAhVZVZ(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*
      ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2))
      )*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CphhhhVZVZ(int gI1, int gI2) const
{
   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Sin(2*
      ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2))
      )*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmVZVZ(int gI1, int gI2) const
{
   const std::complex<double> result = 0.1*((-7.745966692414834*g1*g2*Cos(
      ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(
      ThetaW())))*ZP(gI1,0)*ZP(gI2,0) + (-7.745966692414834*g1*g2*Cos(ThetaW())*
      Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))
      *ZP(gI1,1)*ZP(gI2,1) + 20*Sqr(g2)*Sqr(Cos(ThetaW()))*(ZP(gI1,2)*ZP(gI2,2) +
      ZP(gI1,3)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpAhhhVZ(int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) - ZA(
      gI2,1)*ZH(gI1,1));

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

std::complex<double> CLASSNAME::CpbarChiChiVZPL(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*(Conj(ZN1(gI2,2))*ZN1(gI1,2) - Conj(ZN1
      (gI2,3))*ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiChiVZPR(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*(Conj(ZN2(gI1,2))*ZN2(gI2,2) - Conj(ZN2
      (gI1,3))*ZN2(gI2,3));

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
   const std::complex<double> result = 0.1*((-7.745966692414834*g1*g2*Cos(
      ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(
      ThetaW())))*SUM(j1,0,2,Conj(ZE(gI1,j1))*ZE(gI2,j1)) + 12*Sqr(g1)*Sqr(Sin(
      ThetaW()))*SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*ZE(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZVZ(int gI1, int gI2) const
{
   const std::complex<double> result = 0.03333333333333333*((-7.745966692414834
      *g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 15*Sqr(g2)*Sqr(Cos(ThetaW())) + Sqr(g1)
      *Sqr(Sin(ThetaW())))*SUM(j1,0,2,Conj(ZU(gI1,j1))*ZU(gI2,j1)) + 16*Sqr(g1)*
      Sqr(Sin(ThetaW()))*SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*ZU(gI2,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdVZ(int gI2, int gI1) const
{
   const std::complex<double> result = 0.16666666666666666*(-((3*g2*Cos(ThetaW(
      )) + 0.7745966692414834*g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZD(gI1
      ,j1))) + 1.5491933384829668*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZD(gI2,3 + j1))
      *ZD(gI1,3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeVZ(int gI2, int gI1) const
{
   const std::complex<double> result = 0.1*((-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZE(gI2,j1))*ZE(gI1,j1))
      + 7.745966692414834*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZE(gI2,3 + j1))*ZE(gI1,
      3 + j1)));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuVZ(int gI2, int gI1) const
{
   const std::complex<double> result = 0.16666666666666666*((3*g2*Cos(ThetaW())
      - 0.7745966692414834*g1*Sin(ThetaW()))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZU(gI1,
      j1)) - 3.0983866769659336*g1*Sin(ThetaW())*SUM(j1,0,2,Conj(ZU(gI2,3 + j1))*
      ZU(gI1,3 + j1)));

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

std::complex<double> CLASSNAME::CpSRumconjRhconjVWm(int gI1) const
{
   const std::complex<double> result = 0.7071067811865475*g2*ZHR(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpRhconjRhconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*(ZHR(gI1,0)*ZHR(gI2,0) + ZHR
      (gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1ChiconjVWmPL(int gI1, int gI2) const
{
   const std::complex<double> result = g2*Conj(ZN1(gI2,1))*UP1(gI1,0) -
      0.7071067811865475*g2*Conj(ZN1(gI2,2))*UP1(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1ChiconjVWmPR(int gI1, int gI2) const
{
   const std::complex<double> result = g2*Conj(UM1(gI1,0))*ZN2(gI2,1) +
      0.7071067811865475*g2*Conj(UM1(gI1,1))*ZN2(gI2,2);

   return result;
}

double CLASSNAME::CpSvconjSvconjVWmVWm(int gI1, int gI2) const
{
   const double result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFdconjVWmPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      ZDL(gI2,j1))*ZUL(gI1,j1));

   return result;
}

double CLASSNAME::CpbarFuFdconjVWmPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjVWmPL(int gI1, int gI2) const
{
   const std::complex<double> result = IF(gI1 < 3,-0.7071067811865475*g2*Conj(
      ZEL(gI2,gI1)),0);

   return result;
}

double CLASSNAME::CpbarFvFeconjVWmPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSvconjVWm(int gI2, int gI1) const
{
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZE
      (gI2,j1))*ZV(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpAhAhconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(
      gI1,1)*ZA(gI2,1) + 4*ZA(gI1,3)*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CphhhhconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(
      gI1,1)*ZH(gI2,1) + 4*ZH(gI1,3)*ZH(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpHpmconjHpmconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(
      gI1,1)*ZP(gI2,1) + 2*ZP(gI1,2)*ZP(gI2,2) + 2*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiCha2conjVWmPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.5*g2*(2*Conj(UM2(gI2,0))*ZN1(gI1,1) +
      1.4142135623730951*Conj(UM2(gI2,1))*ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarChiCha2conjVWmPR(int gI1, int gI2) const
{
   const std::complex<double> result = -(g2*Conj(ZN2(gI1,1))*UP2(gI2,0)) +
      0.7071067811865475*g2*Conj(ZN2(gI1,3))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpAhHpmconjVWm(int gI2, int gI1) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*g2*(ZA(gI2,0
      )*ZP(gI1,0) + ZA(gI2,1)*ZP(gI1,1) + 1.4142135623730951*ZA(gI2,3)*(ZP(gI1,2)
      - ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CphhHpmconjVWm(int gI2, int gI1) const
{
   const std::complex<double> result = -0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)
      *ZP(gI1,1) + 1.4142135623730951*ZH(gI2,3)*(ZP(gI1,2) + ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSdconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZD(gI1,j1))*
      ZD(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSeconjSeconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZE(gI1,j1))*
      ZE(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSuconjSuconjVWmVWm(int gI1, int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*SUM(j1,0,2,Conj(ZU(gI1,j1))*
      ZU(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpSdconjSuconjVWm(int gI2, int gI1) const
{
   const std::complex<double> result = 0.7071067811865475*g2*SUM(j1,0,2,Conj(ZD
      (gI2,j1))*ZU(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CphhconjVWmVWm(int gI2) const
{
   const std::complex<double> result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1)
      + 4*vT*ZH(gI2,3));

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

std::complex<double> CLASSNAME::CpbarCha1barUChiSRdpPL(int gI1, int gO2) const
{
   const std::complex<double> result = Conj(UM1(gI1,1))*(-(LamSD*KroneckerDelta
      (0,gO2)) + 0.7071067811865475*LamTD*KroneckerDelta(1,gO2)) - LamTD*Conj(UM1(
      gI1,0))*KroneckerDelta(2,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUChiSRdpPR(int gI1, int gO1) const
{
   const std::complex<double> result = -(g2*KroneckerDelta(2,gO1)*UP1(gI1,0)) -
      0.1414213562373095*(3.872983346207417*g1*KroneckerDelta(0,gO1) + 5*g2*
      KroneckerDelta(1,gO1))*UP1(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barUChiSRumPL(int gI1, int gO2) const
{
   const std::complex<double> result = Conj(UP2(gI1,1))*(LamSU*KroneckerDelta(0
      ,gO2) + 0.7071067811865475*LamTU*KroneckerDelta(1,gO2)) + LamTU*Conj(UP2(gI1
      ,0))*KroneckerDelta(3,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barUChiSRumPR(int gI1, int gO1) const
{
   const std::complex<double> result = -(g2*KroneckerDelta(3,gO1)*UM2(gI1,0)) +
      0.5477225575051661*g1*KroneckerDelta(0,gO1)*UM2(gI1,1) + 0.7071067811865475
      *g2*KroneckerDelta(1,gO1)*UM2(gI1,1);

   return result;
}

double CLASSNAME::CpbarUChibarFvSvPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFvSvPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gI1 < 3,0.5477225575051661*g1*Conj(ZV
      (gI2,gI1))*KroneckerDelta(0,gO1),0) + IF(gI1 < 3,-0.7071067811865475*g2*Conj
      (ZV(gI2,gI1))*KroneckerDelta(1,gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFdSdPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = -(KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(
      ZD(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1,j1))*Yd(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFdSdPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = -0.04714045207910316*(3.872983346207417*
      g1*KroneckerDelta(0,gO1) - 15*g2*KroneckerDelta(1,gO1))*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*ZDL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFeSePL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = -(KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(
      ZE(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1,j1))*Ye(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFeSePR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = 0.1414213562373095*(3.872983346207417*g1
      *KroneckerDelta(0,gO1) + 5*g2*KroneckerDelta(1,gO1))*SUM(j1,0,2,Conj(ZE(gI2,
      j1))*ZEL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFuSuPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = -(KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(
      ZU(gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFuSuPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = -0.04714045207910316*(3.872983346207417*
      g1*KroneckerDelta(0,gO1) + 15*g2*KroneckerDelta(1,gO1))*SUM(j1,0,2,Conj(ZU(
      gI2,j1))*ZUL(gI1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUChiRhPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*LamTD*Conj(ZN2(gI1,1)
      )*KroneckerDelta(2,gO2)*ZHR(gI2,0) + Conj(ZN2(gI1,2))*(LamSD*KroneckerDelta(
      0,gO2)*ZHR(gI2,0) + 0.7071067811865475*LamTD*KroneckerDelta(1,gO2)*ZHR(gI2,0
      )) - LamSU*Conj(ZN2(gI1,3))*KroneckerDelta(0,gO2)*ZHR(gI2,1) +
      0.7071067811865475*LamTU*Conj(ZN2(gI1,3))*KroneckerDelta(1,gO2)*ZHR(gI2,1) +
      0.7071067811865475*LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(3,gO2)*ZHR(gI2,1)
      + Conj(ZN2(gI1,0))*(LamSD*KroneckerDelta(2,gO2)*ZHR(gI2,0) - LamSU*
      KroneckerDelta(3,gO2)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUChiRhPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = 0.1414213562373095*(KroneckerDelta(3,gO1
      )*ZHR(gI2,1)*(3.872983346207417*g1*ZN1(gI1,0) - 5*g2*ZN1(gI1,1)) +
      KroneckerDelta(2,gO1)*ZHR(gI2,0)*(-3.872983346207417*g1*ZN1(gI1,0) + 5*g2*
      ZN1(gI1,1)) - (3.872983346207417*g1*KroneckerDelta(0,gO1) - 5*g2*
      KroneckerDelta(1,gO1))*(ZHR(gI2,0)*ZN1(gI1,2) - ZHR(gI2,1)*ZN1(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha1HpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = Conj(UP1(gI2,1))*(-(LamSD*KroneckerDelta
      (0,gO2)*ZP(gI1,0)) + 0.7071067811865475*LamTD*KroneckerDelta(1,gO2)*ZP(gI1,0
      ) - LamTD*KroneckerDelta(2,gO2)*ZP(gI1,2)) + g2*Conj(UP1(gI2,0))*(-(
      KroneckerDelta(3,gO2)*ZP(gI1,1)) + 1.4142135623730951*KroneckerDelta(1,gO2)*
      ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha1HpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = 0.5477225575051661*g1*KroneckerDelta(0,
      gO1)*UM1(gI2,1)*ZP(gI1,0) - Conj(LamTU)*KroneckerDelta(3,gO1)*UM1(gI2,0)*ZP(
      gI1,1) + 0.7071067811865475*g2*KroneckerDelta(1,gO1)*(UM1(gI2,1)*ZP(gI1,0) +
      2*UM1(gI2,0)*ZP(gI1,2)) + Conj(LamTD)*KroneckerDelta(2,gO1)*UM1(gI2,1)*ZP(
      gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha2conjHpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -(g2*Conj(UM2(gI2,0))*(KroneckerDelta(2,
      gO2)*ZP(gI1,0) + 1.4142135623730951*KroneckerDelta(1,gO2)*ZP(gI1,2))) + 0.5*
      Conj(UM2(gI2,1))*(2*LamSU*KroneckerDelta(0,gO2)*ZP(gI1,1) +
      1.4142135623730951*LamTU*KroneckerDelta(1,gO2)*ZP(gI1,1) + 2*LamTU*
      KroneckerDelta(3,gO2)*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha2conjHpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = Conj(LamTD)*KroneckerDelta(2,gO1)*UP2(
      gI2,0)*ZP(gI1,0) - 0.5477225575051661*g1*KroneckerDelta(0,gO1)*UP2(gI2,1)*ZP
      (gI1,1) - 0.7071067811865475*g2*KroneckerDelta(1,gO1)*UP2(gI2,1)*ZP(gI1,1) -
      Conj(LamTU)*KroneckerDelta(3,gO1)*UP2(gI2,1)*ZP(gI1,2) - 1.4142135623730951
      *g2*KroneckerDelta(1,gO1)*UP2(gI2,0)*ZP(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChiAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,0.1)*(
      3.872983346207417*g1*Conj(ZN1(gI1,0))*(-(KroneckerDelta(2,gO2)*ZA(gI2,0)) +
      KroneckerDelta(3,gO2)*ZA(gI2,1)) + 5*Conj(ZN1(gI1,2))*(1.4142135623730951*
      LamSD*KroneckerDelta(0,gO2)*ZA(gI2,0) + LamTD*KroneckerDelta(1,gO2)*ZA(gI2,0
      ) + KroneckerDelta(2,gO2)*(1.4142135623730951*LamSD*ZA(gI2,2) + LamTD*ZA(gI2
      ,3))) + 5*(g2*Conj(ZN1(gI1,1))*(KroneckerDelta(2,gO2)*ZA(gI2,0) -
      KroneckerDelta(3,gO2)*ZA(gI2,1)) + Conj(ZN1(gI1,3))*(-1.4142135623730951*
      LamSU*KroneckerDelta(0,gO2)*ZA(gI2,1) + LamTU*KroneckerDelta(1,gO2)*ZA(gI2,1
      ) + KroneckerDelta(3,gO2)*(-1.4142135623730951*LamSU*ZA(gI2,2) + LamTU*ZA(
      gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChiAhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.1)*(5*Conj(
      LamTD)*KroneckerDelta(2,gO1)*ZA(gI2,0)*ZN2(gI1,1) + 5*Conj(LamTU)*
      KroneckerDelta(3,gO1)*ZA(gI2,1)*ZN2(gI1,1) - 3.872983346207417*g1*
      KroneckerDelta(0,gO1)*ZA(gI2,0)*ZN2(gI1,2) + 5*g2*KroneckerDelta(1,gO1)*ZA(
      gI2,0)*ZN2(gI1,2) + 5*Conj(LamTD)*KroneckerDelta(2,gO1)*ZA(gI2,3)*ZN2(gI1,2)
      + 7.0710678118654755*Conj(LamSD)*KroneckerDelta(2,gO1)*(ZA(gI2,0)*ZN2(gI1,0
      ) + ZA(gI2,2)*ZN2(gI1,2)) + 3.872983346207417*g1*KroneckerDelta(0,gO1)*ZA(
      gI2,1)*ZN2(gI1,3) - 5*g2*KroneckerDelta(1,gO1)*ZA(gI2,1)*ZN2(gI1,3) + 5*Conj
      (LamTU)*KroneckerDelta(3,gO1)*ZA(gI2,3)*ZN2(gI1,3) - 7.0710678118654755*Conj
      (LamSU)*KroneckerDelta(3,gO1)*(ZA(gI2,1)*ZN2(gI1,0) + ZA(gI2,2)*ZN2(gI1,3)))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChihhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.1*(3.872983346207417*g1*Conj(ZN1(gI2,0
      ))*(KroneckerDelta(2,gO2)*ZH(gI1,0) - KroneckerDelta(3,gO2)*ZH(gI1,1)) + 5*
      Conj(ZN1(gI2,2))*(1.4142135623730951*LamSD*KroneckerDelta(0,gO2)*ZH(gI1,0) +
      LamTD*KroneckerDelta(1,gO2)*ZH(gI1,0) + KroneckerDelta(2,gO2)*(
      1.4142135623730951*LamSD*ZH(gI1,2) + LamTD*ZH(gI1,3))) + 5*(Conj(ZN1(gI2,1))
      *(-(g2*KroneckerDelta(2,gO2)*ZH(gI1,0)) + g2*KroneckerDelta(3,gO2)*ZH(gI1,1)
      ) + Conj(ZN1(gI2,3))*(-1.4142135623730951*LamSU*KroneckerDelta(0,gO2)*ZH(gI1
      ,1) + LamTU*KroneckerDelta(1,gO2)*ZH(gI1,1) + KroneckerDelta(3,gO2)*(
      -1.4142135623730951*LamSU*ZH(gI1,2) + LamTU*ZH(gI1,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChihhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = 0.1*(5*Conj(LamTD)*KroneckerDelta(2,gO1)
      *ZH(gI1,0)*ZN2(gI2,1) + 5*Conj(LamTU)*KroneckerDelta(3,gO1)*ZH(gI1,1)*ZN2(
      gI2,1) + 3.872983346207417*g1*KroneckerDelta(0,gO1)*ZH(gI1,0)*ZN2(gI2,2) - 5
      *g2*KroneckerDelta(1,gO1)*ZH(gI1,0)*ZN2(gI2,2) + 5*Conj(LamTD)*
      KroneckerDelta(2,gO1)*ZH(gI1,3)*ZN2(gI2,2) + 7.0710678118654755*Conj(LamSD)*
      KroneckerDelta(2,gO1)*(ZH(gI1,0)*ZN2(gI2,0) + ZH(gI1,2)*ZN2(gI2,2)) -
      3.872983346207417*g1*KroneckerDelta(0,gO1)*ZH(gI1,1)*ZN2(gI2,3) + 5*g2*
      KroneckerDelta(1,gO1)*ZH(gI1,1)*ZN2(gI2,3) + 5*Conj(LamTU)*KroneckerDelta(3,
      gO1)*ZH(gI1,3)*ZN2(gI2,3) - 7.0710678118654755*Conj(LamSU)*KroneckerDelta(3,
      gO1)*(ZH(gI1,1)*ZN2(gI2,0) + ZH(gI1,2)*ZN2(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiFdconjSdPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -(KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(
      ZDL(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiFdconjSdPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -0.3651483716701107*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZD(gI1,3 + j1)*ZDR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiFeconjSePL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -(KroneckerDelta(2,gO2)*SUM(j2,0,2,Conj(
      ZEL(gI2,j2))*SUM(j1,0,2,Ye(j1,j2)*ZE(gI1,3 + j1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiFeconjSePR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -1.0954451150103321*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZE(gI1,3 + j1)*ZER(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiFuconjSuPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -(KroneckerDelta(3,gO2)*SUM(j2,0,2,Conj(
      ZUL(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiFuconjSuPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = 0.7302967433402214*g1*KroneckerDelta(0,
      gO1)*SUM(j1,0,2,ZU(gI1,3 + j1)*ZUR(gI2,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha1VWmPR(int gO2, int gI2) const
{
   const std::complex<double> result = g2*KroneckerDelta(1,gO2)*UM1(gI2,0) +
      0.7071067811865475*g2*KroneckerDelta(2,gO2)*UM1(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha1VWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = g2*Conj(UP1(gI2,0))*KroneckerDelta(1,gO1
      ) - 0.7071067811865475*g2*Conj(UP1(gI2,1))*KroneckerDelta(2,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha2conjVWmPR(int gO2, int gI2) const
{
   const std::complex<double> result = -(g2*KroneckerDelta(1,gO2)*UP2(gI2,0)) +
      0.7071067811865475*g2*KroneckerDelta(3,gO2)*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiCha2conjVWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = -0.5*g2*(2*Conj(UM2(gI2,0))*
      KroneckerDelta(1,gO1) + 1.4142135623730951*Conj(UM2(gI2,1))*KroneckerDelta(3
      ,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChiVZPR(int gO2, int gI2) const
{
   const std::complex<double> result = 0.1*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(2,gO2)*ZN2(gI2,2) -
      KroneckerDelta(3,gO2)*ZN2(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChiVZPL(int gO1, int gI2) const
{
   const std::complex<double> result = 0.1*(Conj(ZN1(gI2,2))*KroneckerDelta(2,
      gO1) - Conj(ZN1(gI2,3))*KroneckerDelta(3,gO1))*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barCha2RhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = LamTD*Conj(UP2(gI1,0))*SUM(gl1385,0,1,
      Conj(UM1(gl1385,1))*UP1(gl1385,gO2))*ZHR(gI2,0) - LamTU*Conj(UP2(gI1,1))*SUM
      (gl1385,0,1,Conj(UM1(gl1385,0))*UP1(gl1385,gO2))*ZHR(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barCha2RhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = -(g2*(SUM(gl1388,0,1,Conj(UM1(gl1388,gO1
      ))*UP1(gl1388,1))*UM2(gI1,0)*ZHR(gI2,0) + SUM(gl1388,0,1,Conj(UM1(gl1388,gO1
      ))*UP1(gl1388,0))*UM2(gI1,1)*ZHR(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1AhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Conj(UP1
      (gI1,0))*(1.4142135623730951*SUM(gl1391,0,1,Conj(UM1(gl1391,1))*UP1(gl1391,
      gO2))*ZA(gI2,0) + 2*SUM(gl1391,0,1,Conj(UM1(gl1391,0))*UP1(gl1391,gO2))*ZA(
      gI2,3)) - Conj(UP1(gI1,1))*(1.4142135623730951*LamTD*SUM(gl1391,0,1,Conj(UM1
      (gl1391,0))*UP1(gl1391,gO2))*ZA(gI2,0) + SUM(gl1391,0,1,Conj(UM1(gl1391,1))*
      UP1(gl1391,gO2))*(1.4142135623730951*LamSD*ZA(gI2,2) - LamTD*ZA(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1AhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(
      1.4142135623730951*Conj(LamSD)*SUM(gl1394,0,1,Conj(UM1(gl1394,gO1))*UP1(
      gl1394,1))*UM1(gI1,1)*ZA(gI2,2) - g2*SUM(gl1394,0,1,Conj(UM1(gl1394,gO1))*
      UP1(gl1394,0))*(1.4142135623730951*UM1(gI1,1)*ZA(gI2,0) + 2*UM1(gI1,0)*ZA(
      gI2,3)) + Conj(LamTD)*SUM(gl1394,0,1,Conj(UM1(gl1394,gO1))*UP1(gl1394,1))*(
      1.4142135623730951*UM1(gI1,0)*ZA(gI2,0) - UM1(gI1,1)*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFeSvPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = SUM(gl1397,0,1,Conj(UM1(gl1397,1))*UP1(
      gl1397,gO2))*SUM(j2,0,2,Conj(ZV(gI2,j2))*SUM(j1,0,2,Conj(ZER(gI1,j1))*Ye(j1,
      j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFeSvPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = -(g2*SUM(gl1400,0,1,Conj(UM1(gl1400,gO1)
      )*UP1(gl1400,0))*SUM(j1,0,2,Conj(ZV(gI2,j1))*ZEL(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFdSuPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = SUM(gl1403,0,1,Conj(UM1(gl1403,1))*UP1(
      gl1403,gO2))*SUM(j2,0,2,Conj(ZU(gI2,j2))*SUM(j1,0,2,Conj(ZDR(gI1,j1))*Yd(j1,
      j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFdSuPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = -(g2*SUM(gl1406,0,1,Conj(UM1(gl1406,gO1)
      )*UP1(gl1406,0))*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZDL(gI1,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1hhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(-(g2*Conj(UP1(gI2,0))*(
      1.4142135623730951*SUM(gl1409,0,1,Conj(UM1(gl1409,1))*UP1(gl1409,gO2))*ZH(
      gI1,0) + 2*SUM(gl1409,0,1,Conj(UM1(gl1409,0))*UP1(gl1409,gO2))*ZH(gI1,3))) -
      Conj(UP1(gI2,1))*(1.4142135623730951*LamTD*SUM(gl1409,0,1,Conj(UM1(gl1409,0
      ))*UP1(gl1409,gO2))*ZH(gI1,0) + SUM(gl1409,0,1,Conj(UM1(gl1409,1))*UP1(
      gl1409,gO2))*(1.4142135623730951*LamSD*ZH(gI1,2) - LamTD*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1hhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(-1.4142135623730951*Conj(LamSD)*SUM
      (gl1412,0,1,Conj(UM1(gl1412,gO1))*UP1(gl1412,1))*UM1(gI2,1)*ZH(gI1,2) - g2*
      SUM(gl1412,0,1,Conj(UM1(gl1412,gO1))*UP1(gl1412,0))*(1.4142135623730951*UM1(
      gI2,1)*ZH(gI1,0) + 2*UM1(gI2,0)*ZH(gI1,3)) + Conj(LamTD)*SUM(gl1412,0,1,Conj
      (UM1(gl1412,gO1))*UP1(gl1412,1))*(-1.4142135623730951*UM1(gI2,0)*ZH(gI1,0) +
      UM1(gI2,1)*ZH(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1ChiconjHpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.5477225575051661*g1*Conj(ZN1(gI2,0))*
      SUM(gl1415,0,1,Conj(UM1(gl1415,1))*UP1(gl1415,gO2))*ZP(gI1,0) - LamTU*Conj(
      ZN1(gI2,3))*SUM(gl1415,0,1,Conj(UM1(gl1415,0))*UP1(gl1415,gO2))*ZP(gI1,1) +
      0.7071067811865475*g2*Conj(ZN1(gI2,1))*(SUM(gl1415,0,1,Conj(UM1(gl1415,1))*
      UP1(gl1415,gO2))*ZP(gI1,0) + 2*SUM(gl1415,0,1,Conj(UM1(gl1415,0))*UP1(gl1415
      ,gO2))*ZP(gI1,2)) + LamTD*Conj(ZN1(gI2,2))*SUM(gl1415,0,1,Conj(UM1(gl1415,1)
      )*UP1(gl1415,gO2))*ZP(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1ChiconjHpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -(Conj(LamSD)*SUM(gl1418,0,1,Conj(UM1(
      gl1418,gO1))*UP1(gl1418,1))*ZN2(gI2,0)*ZP(gI1,0)) + 0.5*Conj(LamTD)*SUM(
      gl1418,0,1,Conj(UM1(gl1418,gO1))*UP1(gl1418,1))*(1.4142135623730951*ZN2(gI2,
      1)*ZP(gI1,0) - 2*ZN2(gI2,2)*ZP(gI1,2)) + g2*SUM(gl1418,0,1,Conj(UM1(gl1418,
      gO1))*UP1(gl1418,0))*(-(ZN2(gI2,3)*ZP(gI1,1)) + 1.4142135623730951*ZN2(gI2,1
      )*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barChiSRdpPL(int gO2, int gI1) const
{
   const std::complex<double> result = -(LamTD*Conj(ZN2(gI1,2))*SUM(gl1421,0,1,
      Conj(UM1(gl1421,0))*UP1(gl1421,gO2))) + 0.5*(-2*LamSD*Conj(ZN2(gI1,0)) +
      1.4142135623730951*LamTD*Conj(ZN2(gI1,1)))*SUM(gl1421,0,1,Conj(UM1(gl1421,1)
      )*UP1(gl1421,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barChiSRdpPR(int gO1, int gI1) const
{
   const std::complex<double> result = -0.1*SUM(gl1424,0,1,Conj(UM1(gl1424,gO1)
      )*UP1(gl1424,1))*(5.477225575051661*g1*ZN1(gI1,0) + 7.0710678118654755*g2*
      ZN1(gI1,1)) - g2*SUM(gl1424,0,1,Conj(UM1(gl1424,gO1))*UP1(gl1424,0))*ZN1(gI1
      ,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1FuconjSdPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(gl1427,0,1,Conj(UM1(gl1427,1))*UP1(
      gl1427,gO2))*SUM(j2,0,2,Conj(ZUL(gI2,j2))*SUM(j1,0,2,Yd(j1,j2)*ZD(gI1,3 + j1
      )));

   return result;
}

double CLASSNAME::CpbarUCha1FuconjSdPR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1FvconjSePL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(gl1433,0,1,Conj(UM1(gl1433,1))*UP1(
      gl1433,gO2))*SUM(j1,0,2,Ye(j1,gI2)*ZE(gI1,3 + j1));

   return result;
}

double CLASSNAME::CpbarUCha1FvconjSePR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1VPPR(int gO2, int gI2) const
{
   const std::complex<double> result = -(g2*Sin(ThetaW())*SUM(gl1439,0,1,Conj(
      UM1(gl1439,0))*UP1(gl1439,gO2))*UM1(gI2,0)) - 0.1*(3.872983346207417*g1*Cos(
      ThetaW()) + 5*g2*Sin(ThetaW()))*SUM(gl1439,0,1,Conj(UM1(gl1439,1))*UP1(
      gl1439,gO2))*UM1(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1VPPL(int gO1, int gI2) const
{
   const std::complex<double> result = -(g2*Conj(UP1(gI2,0))*Sin(ThetaW())*SUM(
      gl1442,0,1,Conj(UM1(gl1442,gO1))*UP1(gl1442,0))) - 0.1*Conj(UP1(gI2,1))*(
      3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW()))*SUM(gl1442,0,1,Conj
      (UM1(gl1442,gO1))*UP1(gl1442,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1VZPR(int gO2, int gI2) const
{
   const std::complex<double> result = 0.1*(-10*g2*Cos(ThetaW())*SUM(gl1445,0,1
      ,Conj(UM1(gl1445,0))*UP1(gl1445,gO2))*UM1(gI2,0) + (-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*SUM(gl1445,0,1,Conj(UM1(gl1445,1))*UP1(
      gl1445,gO2))*UM1(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1VZPL(int gO1, int gI2) const
{
   const std::complex<double> result = 0.1*(-10*g2*Conj(UP1(gI2,0))*Cos(ThetaW(
      ))*SUM(gl1448,0,1,Conj(UM1(gl1448,gO1))*UP1(gl1448,0)) + Conj(UP1(gI2,1))*(
      -5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*SUM(gl1448,0,1,
      Conj(UM1(gl1448,gO1))*UP1(gl1448,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1ChiconjVWmPR(int gO2, int gI2) const
{
   const std::complex<double> result = g2*SUM(gl1451,0,1,Conj(UM1(gl1451,0))*
      UP1(gl1451,gO2))*ZN2(gI2,1) + 0.7071067811865475*g2*SUM(gl1451,0,1,Conj(UM1(
      gl1451,1))*UP1(gl1451,gO2))*ZN2(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1ChiconjVWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = g2*Conj(ZN1(gI2,1))*SUM(gl1454,0,1,Conj(
      UM1(gl1454,gO1))*UP1(gl1454,0)) - 0.7071067811865475*g2*Conj(ZN1(gI2,2))*SUM
      (gl1454,0,1,Conj(UM1(gl1454,gO1))*UP1(gl1454,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUCha2RhPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = LamTD*Conj(UM1(gI1,1))*KroneckerDelta(0,
      gO2)*ZHR(gI2,0) - LamTU*Conj(UM1(gI1,0))*KroneckerDelta(1,gO2)*ZHR(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUCha2RhPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = -(g2*(KroneckerDelta(0,gO1)*UP1(gI1,1)*
      ZHR(gI2,0) + KroneckerDelta(1,gO1)*UP1(gI1,0)*ZHR(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2AhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*(g2*Conj(UM2
      (gI1,0))*(1.4142135623730951*KroneckerDelta(1,gO2)*ZA(gI2,1) - 2*
      KroneckerDelta(0,gO2)*ZA(gI2,3)) + Conj(UM2(gI1,1))*(1.4142135623730951*
      LamTU*KroneckerDelta(0,gO2)*ZA(gI2,1) + KroneckerDelta(1,gO2)*(
      1.4142135623730951*LamSU*ZA(gI2,2) + LamTU*ZA(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2AhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0,-0.5)*(
      1.4142135623730951*Conj(LamSU)*KroneckerDelta(1,gO1)*UP2(gI1,1)*ZA(gI2,2) +
      g2*KroneckerDelta(0,gO1)*(1.4142135623730951*UP2(gI1,1)*ZA(gI2,1) - 2*UP2(
      gI1,0)*ZA(gI2,3)) + Conj(LamTU)*KroneckerDelta(1,gO1)*(1.4142135623730951*
      UP2(gI1,0)*ZA(gI2,1) + UP2(gI1,1)*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFuSdPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(ZD
      (gI2,j2))*SUM(j1,0,2,Conj(ZUR(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFuSdPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*SUM(j1,0,2,
      Conj(ZD(gI2,j1))*ZUL(gI1,j1)));

   return result;
}

double CLASSNAME::CpbarUCha2barFvSePL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFvSePR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gI1 < 3,-(g2*Conj(ZE(gI2,gI1))*
      KroneckerDelta(0,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2hhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(g2*Conj(UM2(gI2,0))*(
      -1.4142135623730951*KroneckerDelta(1,gO2)*ZH(gI1,1) + 2*KroneckerDelta(0,gO2
      )*ZH(gI1,3)) + Conj(UM2(gI2,1))*(1.4142135623730951*LamTU*KroneckerDelta(0,
      gO2)*ZH(gI1,1) + KroneckerDelta(1,gO2)*(1.4142135623730951*LamSU*ZH(gI1,2) +
      LamTU*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2hhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = 0.5*(1.4142135623730951*Conj(LamSU)*
      KroneckerDelta(1,gO1)*UP2(gI2,1)*ZH(gI1,2) + KroneckerDelta(0,gO1)*(
      -1.4142135623730951*g2*UP2(gI2,1)*ZH(gI1,1) + 2*g2*UP2(gI2,0)*ZH(gI1,3)) +
      Conj(LamTU)*KroneckerDelta(1,gO1)*(1.4142135623730951*UP2(gI2,0)*ZH(gI1,1) +
      UP2(gI2,1)*ZH(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2ChiHpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = LamTD*Conj(ZN1(gI2,2))*KroneckerDelta(0,
      gO2)*ZP(gI1,0) - 0.5477225575051661*g1*Conj(ZN1(gI2,0))*KroneckerDelta(1,gO2
      )*ZP(gI1,1) - 0.7071067811865475*g2*Conj(ZN1(gI2,1))*KroneckerDelta(1,gO2)*
      ZP(gI1,1) - LamTU*Conj(ZN1(gI2,3))*KroneckerDelta(1,gO2)*ZP(gI1,2) -
      1.4142135623730951*g2*Conj(ZN1(gI2,1))*KroneckerDelta(0,gO2)*ZP(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2ChiHpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO1)*(ZN2(gI2,2)*
      ZP(gI1,0) + 1.4142135623730951*ZN2(gI2,1)*ZP(gI1,2))) + 0.5*KroneckerDelta(1
      ,gO1)*(2*Conj(LamSU)*ZN2(gI2,0)*ZP(gI1,1) + Conj(LamTU)*(1.4142135623730951*
      ZN2(gI2,1)*ZP(gI1,1) + 2*ZN2(gI2,3)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barChiSRumPL(int gO2, int gI1) const
{
   const std::complex<double> result = LamTU*Conj(ZN2(gI1,3))*KroneckerDelta(0,
      gO2) + LamSU*Conj(ZN2(gI1,0))*KroneckerDelta(1,gO2) + 0.7071067811865475*
      LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(1,gO2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barChiSRumPR(int gO1, int gI1) const
{
   const std::complex<double> result = KroneckerDelta(1,gO1)*(
      0.5477225575051661*g1*ZN1(gI1,0) + 0.7071067811865475*g2*ZN1(gI1,1)) - g2*
      KroneckerDelta(0,gO1)*ZN1(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2FdconjSuPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = KroneckerDelta(1,gO2)*SUM(j2,0,2,Conj(
      ZDL(gI2,j2))*SUM(j1,0,2,Yu(j1,j2)*ZU(gI1,3 + j1)));

   return result;
}

double CLASSNAME::CpbarUCha2FdconjSuPR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2VPPR(int gO2, int gI2) const
{
   const std::complex<double> result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*
      UP2(gI2,0) + 0.1*KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) +
      5*g2*Sin(ThetaW()))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2VPPL(int gO1, int gI2) const
{
   const std::complex<double> result = g2*Conj(UM2(gI2,0))*KroneckerDelta(0,gO1
      )*Sin(ThetaW()) + 0.1*Conj(UM2(gI2,1))*KroneckerDelta(1,gO1)*(
      3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2VZPR(int gO2, int gI2) const
{
   const std::complex<double> result = g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*
      UP2(gI2,0) + 0.1*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2VZPL(int gO1, int gI2) const
{
   const std::complex<double> result = g2*Conj(UM2(gI2,0))*Cos(ThetaW())*
      KroneckerDelta(0,gO1) + 0.1*Conj(UM2(gI2,1))*KroneckerDelta(1,gO1)*(5*g2*Cos
      (ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2ChiVWmPR(int gO2, int gI2) const
{
   const std::complex<double> result = -(g2*KroneckerDelta(0,gO2)*ZN2(gI2,1)) +
      0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZN2(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2ChiVWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = -0.5*g2*(2*Conj(ZN1(gI2,1))*
      KroneckerDelta(0,gO1) + 1.4142135623730951*Conj(ZN1(gI2,3))*KroneckerDelta(1
      ,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUFeSvPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,Conj(UM1(gI1,1))*SUM(j2,0,2,
      Conj(ZV(gI2,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUFeSvPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZV(gI2,gO1))*UP1(
      gI1,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gI1,j2))*Ye(gO2,j2))*ZA(gI2,0),0);

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
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(ZEL(gI2,j2))*Ye(gO2,j2))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFehhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Ye(j1,gO1))*ZER(gI2,j1))*ZH(gI1,0),0);

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

std::complex<double> CLASSNAME::CpbarChibarUFeSePL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(Conj(ZN2(gI1,2))*SUM(j2,0,2
      ,Conj(ZE(gI2,j2))*Ye(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUFeSePR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZE
      (gI2,gO1))*ZN1(gI1,0),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZE(gI2,gO1)
      )*ZN1(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,-1.0954451150103321*g1*Conj(
      ZE(gI1,3 + gO2))*Conj(ZN1(gI2,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeChiSePR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,-(SUM(j1,0,2,Conj(Ye(j1,gO1))
      *Conj(ZE(gI1,3 + j1)))*ZN2(gI2,2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.7745966692414834*g1*Cos(
      ThetaW())*ZER(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeVPPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.3872983346207417*g1*Conj(
      ZEL(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(ZEL(gI2,gO1))*Sin(
      ThetaW()),0);

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
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZEL(gI2,gO1))*Cos
      (ThetaW()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Sin(
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
   const double result = IF(gI2 < 3,-0.7071067811865475*g2*KroneckerDelta(gI2,
      gO1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUFdSuPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,Conj(UM1(gI1,1))*SUM(j2,0,2,
      Conj(ZU(gI2,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barUFdSuPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZU(gI2,gO1))*UP1(
      gI1,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(ZDL(gI1,j2))*Yd(gO2,j2))*ZA(gI2,0),0);

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
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(ZDL(gI2,j2))*Yd(gO2,j2))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdhhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Yd(j1,gO1))*ZDR(gI2,j1))*ZH(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(ZUL(gI2,j2))*
      Yd(gO2,j2))*ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFuHpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yu(j1,gO1))*
      ZUR(gI2,j1))*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUFdSdPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(Conj(ZN2(gI1,2))*SUM(j2,0,2
      ,Conj(ZD(gI2,j2))*Yd(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUFdSdPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(
      ZD(gI2,gO1))*ZN1(gI1,0),0) + IF(gO1 < 3,0.7071067811865475*g2*Conj(ZD(gI2,
      gO1))*ZN1(gI1,1),0);

   return result;
}

double CLASSNAME::CpbarUFdCha2SuPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdCha2SuPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yu(j1,gO1))*
      Conj(ZU(gI1,3 + j1)))*UP2(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.3651483716701107*g1*Conj(
      ZD(gI1,3 + gO2))*Conj(ZN1(gI2,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdChiSdPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yd(j1,gO1))
      *Conj(ZD(gI1,3 + j1)))*ZN2(gI2,2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdGluSdPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,1.4142135623730951*g3*Conj(ZD
      (gI1,3 + gO2)),0);

   return result;
}

double CLASSNAME::CpbarUFdGluSdPR(int , int ) const
{
   const double result = 0;

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
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(
      ThetaW())*ZDR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdVPPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(
      ZDL(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(ZDL(gI2,gO1))*Sin(
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
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(ZDL(gI2,gO1))*Cos
      (ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Sin(
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
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      ZUL(gI2,gO1)),0);

   return result;
}

double CLASSNAME::CpbarGlubarUFdSdPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGlubarUFdSdPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      ZD(gI2,gO1)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barUFuSdPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,Conj(UP2(gI1,1))*SUM(j2,0,2,
      Conj(ZD(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barUFuSdPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZD(gI2,gO1))*UM2(
      gI1,0)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(ZUL(gI1,j2))*Yu(gO2,j2))*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*ZUR(gI1,j1))*ZA(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(ZDL(gI2,j2))*
      Yu(gO2,j2))*ZP(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjHpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yd(j1,gO1))*
      ZDR(gI2,j1))*ZP(gI1,0),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(ZUL(gI2,j2))*Yu(gO2,j2))*ZH(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuhhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Yu(j1,gO1))*ZUR(gI2,j1))*ZH(gI1,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUFuSuPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-(Conj(ZN2(gI1,3))*SUM(j2,0,2
      ,Conj(ZU(gI2,j2))*Yu(gO2,j2))),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarUFuSuPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.18257418583505536*g1*Conj(
      ZU(gI2,gO1))*ZN1(gI1,0),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZU(gI2,
      gO1))*ZN1(gI1,1),0);

   return result;
}

double CLASSNAME::CpbarUFuCha1SdPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuCha1SdPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yd(j1,gO1))*
      Conj(ZD(gI1,3 + j1)))*UM1(gI2,1),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,0.7302967433402214*g1*Conj(
      ZN1(gI2,0))*Conj(ZU(gI1,3 + gO2)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuChiSuPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,-(SUM(j1,0,2,Conj(Yu(j1,gO1))
      *Conj(ZU(gI1,3 + j1)))*ZN2(gI2,3)),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuGluSuPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,1.4142135623730951*g3*Conj(ZU
      (gI1,3 + gO2)),0);

   return result;
}

double CLASSNAME::CpbarUFuGluSuPR(int , int ) const
{
   const double result = 0;

   return result;
}

double CLASSNAME::CpbarUFuFdconjVWmPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFdconjVWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      ZDL(gI2,gO1)),0);

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
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(
      ZUL(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,-0.5*g2*Conj(ZUL(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Sin(
      ThetaW())*ZUR(gI2,gO2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(ZUL(gI2,gO1))*
      Cos(ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

double CLASSNAME::CpbarGlubarUFuSuPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGlubarUFuSuPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-1.4142135623730951*g3*Conj(
      ZU(gI2,gO1)),0);

   return result;
}

double CLASSNAME::CpbarGlubarFdSdPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGlubarFdSdPR(int gI1, int gI2) const
{
   const std::complex<double> result = -1.4142135623730951*g3*SUM(j1,0,2,Conj(
      ZD(gI2,j1))*ZDL(gI1,j1));

   return result;
}

double CLASSNAME::CpbarGlubarFuSuPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGlubarFuSuPR(int gI1, int gI2) const
{
   const std::complex<double> result = -1.4142135623730951*g3*SUM(j1,0,2,Conj(
      ZU(gI2,j1))*ZUL(gI1,j1));

   return result;
}

double CLASSNAME::CpbarGluFdconjSdPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluFdconjSdPR(int gI2, int gI1) const
{
   const std::complex<double> result = 1.4142135623730951*g3*SUM(j1,0,2,ZD(gI1,
      3 + j1)*ZDR(gI2,j1));

   return result;
}

double CLASSNAME::CpbarGluFuconjSuPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluFuconjSuPR(int gI2, int gI1) const
{
   const std::complex<double> result = 1.4142135623730951*g3*SUM(j1,0,2,ZU(gI1,
      3 + j1)*ZUR(gI2,j1));

   return result;
}

double CLASSNAME::CpbarCha2barFvSePL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barFvSePR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-(g2*Conj(ZE(gI2,gO1))*UM2(
      gI1,0)),0);

   return result;
}

double CLASSNAME::CpbarChibarFvSvPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFvSvPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,0.5477225575051661*g1*Conj(ZV
      (gI2,gO1))*ZN1(gI1,0),0) + IF(gO1 < 3,-0.7071067811865475*g2*Conj(ZV(gI2,gO1
      ))*ZN1(gI1,1),0);

   return result;
}

double CLASSNAME::CpbarFvFeconjHpmPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvFeconjHpmPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(j1,0,2,Conj(Ye(j1,gO1))*ZER(gI2,j1))
      *ZP(gI1,0);

   return result;
}

double CLASSNAME::CpbarFvCha1SePL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvCha1SePR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(j1,0,2,Conj(Ye(j1,gO1))*Conj(ZE(gI1,
      3 + j1)))*UM1(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barFeSvPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = Conj(UM1(gI1,1))*SUM(j2,0,2,Conj(ZV(gI2,
      j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barFeSvPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZV(gI2,j1))*ZEL(gO1
      ,j1))*UP1(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(ZEL(gI1,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1
      ))*Ye(j1,j2)))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*ZER(gI1,j1))*ZEL(
      gO1,j2))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZEL(
      gI2,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFehhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j1,j2))*ZER(gI2,j1))*ZEL(gO1,j2))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvHpmPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,gI2))
      *ZP(gI1,0);

   return result;
}

double CLASSNAME::CpbarFeFvHpmPR(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFeSePL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = -(Conj(ZN2(gI1,2))*SUM(j2,0,2,Conj(ZE(
      gI2,j2))*SUM(j1,0,2,Conj(ZER(gO2,j1))*Ye(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFeSePR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = 0.7071067811865475*SUM(j1,0,2,Conj(ZE(
      gI2,j1))*ZEL(gO1,j1))*(0.7745966692414834*g1*ZN1(gI1,0) + g2*ZN1(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -1.0954451150103321*g1*Conj(ZN1(gI2,0))*
      SUM(j1,0,2,Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFeChiSePR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*
      Conj(ZE(gI1,3 + j1)))*ZEL(gO1,j2))*ZN2(gI2,2));

   return result;
}

double CLASSNAME::CpbarFeFvVWmPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFvVWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*ZEL(
      gO1,gI2),0);

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barFdSuPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = Conj(UM1(gI1,1))*SUM(j2,0,2,Conj(ZU(gI2,
      j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha1barFdSuPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZU(gI2,j1))*ZDL(gO1
      ,j1))*UP1(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(ZDL(gI1,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1
      ))*Yd(j1,j2)))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*ZDR(gI1,j1))*ZDL(
      gO1,j2))*ZA(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZDL(
      gI2,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2)))*ZH(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdhhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yd(j1,j2))*ZDR(gI2,j1))*ZDL(gO1,j2))*ZH(gI1,0);

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
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*
      ZUR(gI2,j1))*ZDL(gO1,j2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFdSdPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = -(Conj(ZN2(gI1,2))*SUM(j2,0,2,Conj(ZD(
      gI2,j2))*SUM(j1,0,2,Conj(ZDR(gO2,j1))*Yd(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFdSdPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = -0.2357022603955158*SUM(j1,0,2,Conj(ZD(
      gI2,j1))*ZDL(gO1,j1))*(0.7745966692414834*g1*ZN1(gI1,0) - 3*g2*ZN1(gI1,1));

   return result;
}

double CLASSNAME::CpbarFdCha2SuPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdCha2SuPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*
      Conj(ZU(gI1,3 + j1)))*ZDL(gO1,j2))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -0.3651483716701107*g1*Conj(ZN1(gI2,0))*
      SUM(j1,0,2,Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdChiSdPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*
      Conj(ZD(gI1,3 + j1)))*ZDL(gO1,j2))*ZN2(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdGluSdPL(int gO2, int gI1) const
{
   const std::complex<double> result = 1.4142135623730951*g3*SUM(j1,0,2,Conj(ZD
      (gI1,3 + j1))*Conj(ZDR(gO2,j1)));

   return result;
}

double CLASSNAME::CpbarFdGluSdPR(int , int ) const
{
   const double result = 0;

   return result;
}

double CLASSNAME::CpbarFdFuVWmPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFuVWmPL(int gO1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      ZUL(gI2,j1))*ZDL(gO1,j1));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barFuSdPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = Conj(UP2(gI1,1))*SUM(j2,0,2,Conj(ZD(gI2,
      j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarCha2barFuSdPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = -(g2*SUM(j1,0,2,Conj(ZD(gI2,j1))*ZUL(gO1
      ,j1))*UM2(gI1,0));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(int gO2, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(ZUL(gI1,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1
      ))*Yu(j1,j2)))*ZA(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(int gO1, int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*ZUR(gI1,j1))*ZUL(
      gO1,j2))*ZA(gI2,1);

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
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*
      ZDR(gI2,j1))*ZUL(gO1,j2))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(ZUL(
      gI2,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2)))*ZH(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuhhPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yu(j1,j2))*ZUR(gI2,j1))*ZUL(gO1,j2))*ZH(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFuSuPL(int gI1, int gO2, int gI2) const
{
   const std::complex<double> result = -(Conj(ZN2(gI1,3))*SUM(j2,0,2,Conj(ZU(
      gI2,j2))*SUM(j1,0,2,Conj(ZUR(gO2,j1))*Yu(j1,j2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarChibarFuSuPR(int gI1, int gO1, int gI2) const
{
   const std::complex<double> result = -0.2357022603955158*SUM(j1,0,2,Conj(ZU(
      gI2,j1))*ZUL(gO1,j1))*(0.7745966692414834*g1*ZN1(gI1,0) + 3*g2*ZN1(gI1,1));

   return result;
}

double CLASSNAME::CpbarFuCha1SdPL(int , int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuCha1SdPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*
      Conj(ZD(gI1,3 + j1)))*ZUL(gO1,j2))*UM1(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPL(int gO2, int gI2, int gI1) const
{
   const std::complex<double> result = 0.7302967433402214*g1*Conj(ZN1(gI2,0))*
      SUM(j1,0,2,Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuChiSuPR(int gO1, int gI2, int gI1) const
{
   const std::complex<double> result = -(SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*
      Conj(ZU(gI1,3 + j1)))*ZUL(gO1,j2))*ZN2(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuGluSuPL(int gO2, int gI1) const
{
   const std::complex<double> result = 1.4142135623730951*g3*SUM(j1,0,2,Conj(ZU
      (gI1,3 + j1))*Conj(ZUR(gO2,j1)));

   return result;
}

double CLASSNAME::CpbarFuGluSuPR(int , int ) const
{
   const double result = 0;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(A0(Sqr(MSRdp))*CpSRdpUSdconjSRdpconjUSd(gO1,gO2));
   result += -(A0(Sqr(MSRum))*CpSRumUSdconjSRumconjUSd(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSdconjUSdVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhUSdconjRhconjUSd(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,2,(Conj(CpbarCha1FuconjUSdPL(gI1,gI2,gO2))*
      CpbarCha1FuconjUSdPL(gI1,gI2,gO1) + Conj(CpbarCha1FuconjUSdPR(gI1,gI2,gO2))*
      CpbarCha1FuconjUSdPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MFu(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha1(gI1)),
      Sqr(MFu(gI2)))*(Conj(CpbarCha1FuconjUSdPR(gI1,gI2,gO2))*CpbarCha1FuconjUSdPL
      (gI1,gI2,gO1) + Conj(CpbarCha1FuconjUSdPL(gI1,gI2,gO2))*CpbarCha1FuconjUSdPR
      (gI1,gI2,gO1))*MFu(gI2)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSdSvconjUSdconjSv(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,2,SUM(gI2,0,1,(Conj(CpCha2FuconjUSdPL(gI2,gI1,gO2))*
      CpCha2FuconjUSdPL(gI2,gI1,gO1) + Conj(CpCha2FuconjUSdPR(gI2,gI1,gO2))*
      CpCha2FuconjUSdPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MCha2(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(
      MCha2(gI2)))*(Conj(CpCha2FuconjUSdPR(gI2,gI1,gO2))*CpCha2FuconjUSdPL(gI2,gI1
      ,gO1) + Conj(CpCha2FuconjUSdPL(gI2,gI1,gO2))*CpCha2FuconjUSdPR(gI2,gI1,gO1))
      *MCha2(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,(Conj(CpChiFdconjUSdPL(gI2,gI1,gO2))*
      CpChiFdconjUSdPL(gI2,gI1,gO1) + Conj(CpChiFdconjUSdPR(gI2,gI1,gO2))*
      CpChiFdconjUSdPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiFdconjUSdPR(gI2,gI1,gO2))*CpChiFdconjUSdPL(gI2,gI1,
      gO1) + Conj(CpChiFdconjUSdPL(gI2,gI1,gO2))*CpChiFdconjUSdPR(gI2,gI1,gO1))*
      MChi(gI2)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUSdconjUSd(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUSdconjUSd(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSdconjHpmconjUSd(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,2,(Conj(CpbarChiFdconjUSdPL(gI1,gI2,gO2))*
      CpbarChiFdconjUSdPL(gI1,gI2,gO1) + Conj(CpbarChiFdconjUSdPR(gI1,gI2,gO2))*
      CpbarChiFdconjUSdPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MFd(gI2)))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MFd(gI2)))*(Conj(CpbarChiFdconjUSdPR(gI1,gI2,gO2))*CpbarChiFdconjUSdPL(gI1,
      gI2,gO1) + Conj(CpbarChiFdconjUSdPL(gI1,gI2,gO2))*CpbarChiFdconjUSdPR(gI1,
      gI2,gO1))*MFd(gI2)));
   result += 1.3333333333333333*SUM(gI1,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MphiO))
      *Conj(CpphiOSdconjUSd(gI1,gO2))*CpphiOSdconjUSd(gI1,gO1));
   result += SUM(gI1,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSRum))*Conj(
      CpSRumSuconjUSd(gI1,gO2))*CpSRumSuconjUSd(gI1,gO1));
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUSdSeconjUSdconjSe(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MRh(gI2)))*
      Conj(CpRhSdconjUSd(gI2,gI1,gO2))*CpRhSdconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhSdconjUSd(gI2,gI1,gO2))*CpAhSdconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhSdconjUSd(gI2,gI1,gO2))*CphhSdconjUSd(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpHpmSuconjUSd(gI2,gI1,gO2))*CpHpmSuconjUSd(gI2,gI1,gO1)));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpbarGluFdconjUSdPL(gI2,gO2))
      *CpbarGluFdconjUSdPL(gI2,gO1) + Conj(CpbarGluFdconjUSdPR(gI2,gO2))*
      CpbarGluFdconjUSdPR(gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFd(gI2))));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpGluFdconjUSdPL(gI2,gO2))*
      CpGluFdconjUSdPL(gI2,gO1) + Conj(CpGluFdconjUSdPR(gI2,gO2))*CpGluFdconjUSdPR
      (gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFd(gI2))));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFd(
      gI2)))*(Conj(CpbarGluFdconjUSdPR(gI2,gO2))*CpbarGluFdconjUSdPL(gI2,gO1) +
      Conj(CpbarGluFdconjUSdPL(gI2,gO2))*CpbarGluFdconjUSdPR(gI2,gO1))*MFd(gI2));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFd(
      gI2)))*(Conj(CpGluFdconjUSdPR(gI2,gO2))*CpGluFdconjUSdPL(gI2,gO1) + Conj(
      CpGluFdconjUSdPL(gI2,gO2))*CpGluFdconjUSdPR(gI2,gO1))*MFd(gI2));
   result += 1.3333333333333333*SUM(gI2,0,5,B0(Sqr(p),Sqr(MsigmaO),Sqr(MSd(gI2)
      ))*Conj(CpsigmaOSdconjUSd(gI2,gO2))*CpsigmaOSdconjUSd(gI2,gO1));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MSRdp),Sqr(MSu(gI2)))*Conj(
      CpSuconjSRdpconjUSd(gI2,gO2))*CpSuconjSRdpconjUSd(gI2,gO1));
   result += 1.3333333333333333*SUM(gI2,0,5,Conj(CpSdconjUSdVG(gI2,gO2))*
      CpSdconjUSdVG(gI2,gO1)*F0(Sqr(p),Sqr(MSd(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSdconjUSdVP(gI2,gO2))*CpSdconjUSdVP(gI2,gO1)*F0
      (Sqr(p),Sqr(MSd(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSdconjUSdVZ(gI2,gO2))*CpSdconjUSdVZ(gI2,gO1)*F0
      (Sqr(p),Sqr(MSd(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,Conj(CpSuconjUSdVWm(gI2,gO2))*CpSuconjUSdVWm(gI2,gO1)*
      F0(Sqr(p),Sqr(MSu(gI2)),Sqr(MVWm)));

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

   result += -(A0(Sqr(MSRdp))*CpSRdpUSvconjSRdpconjUSv(gO1,gO2));
   result += -(A0(Sqr(MSRum))*CpSRumUSvconjSRumconjUSv(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSvconjUSvVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhUSvconjRhconjUSv(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvUSvconjSvconjUSv(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,2,SUM(gI2,0,1,(Conj(CpCha1FeconjUSvPL(gI2,gI1,gO2))*
      CpCha1FeconjUSvPL(gI2,gI1,gO1) + Conj(CpCha1FeconjUSvPR(gI2,gI1,gO2))*
      CpCha1FeconjUSvPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MCha1(gI2)))));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(
      MCha1(gI2)))*(Conj(CpCha1FeconjUSvPR(gI2,gI1,gO2))*CpCha1FeconjUSvPL(gI2,gI1
      ,gO1) + Conj(CpCha1FeconjUSvPL(gI2,gI1,gO2))*CpCha1FeconjUSvPR(gI2,gI1,gO1))
      *MCha1(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhSvconjUSv(gI2,gI1,gO2))*CpAhSvconjUSv(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhSvconjUSv(gI2,gI1,gO2))*CphhSvconjUSv(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,(Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*
      CpChiFvconjUSvPL(gI2,gI1,gO1) + Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*
      CpChiFvconjUSvPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiFvconjUSvPR(gI2,gI1,gO2))*CpChiFvconjUSvPL(gI2,gI1,
      gO1) + Conj(CpChiFvconjUSvPL(gI2,gI1,gO2))*CpChiFvconjUSvPR(gI2,gI1,gO1))*
      MChi(gI2)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUSvconjUSv(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUSvconjUSv(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSvconjHpmconjUSv(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSe(gI2)))*
      Conj(CpSeconjHpmconjUSv(gI2,gI1,gO2))*CpSeconjHpmconjUSv(gI2,gI1,gO1)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSvconjSdconjUSv(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSvconjSeconjUSv(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuUSvconjSuconjUSv(gI1,gO1,gI1,
      gO2));
   result += SUM(gI2,0,2,Conj(CpSvconjUSvVZ(gI2,gO2))*CpSvconjUSvVZ(gI2,gO1)*F0
      (Sqr(p),Sqr(MSv(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MSRdp),Sqr(MSe(gI2)))*Conj(
      CpSRdpSeconjUSv(gI2,gO2))*CpSRdpSeconjUSv(gI2,gO1));
   result += SUM(gI2,0,5,Conj(CpSeconjUSvconjVWm(gI2,gO2))*CpSeconjUSvconjVWm(
      gI2,gO1)*F0(Sqr(p),Sqr(MSe(gI2)),Sqr(MVWm)));

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

   result += -(A0(Sqr(MSRdp))*CpSRdpUSuconjSRdpconjUSu(gO1,gO2));
   result += -(A0(Sqr(MSRum))*CpSRumUSuconjSRumconjUSu(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSuconjUSuVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhUSuconjRhconjUSu(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,2,(Conj(CpbarCha2FdconjUSuPL(gI1,gI2,gO2))*
      CpbarCha2FdconjUSuPL(gI1,gI2,gO1) + Conj(CpbarCha2FdconjUSuPR(gI1,gI2,gO2))*
      CpbarCha2FdconjUSuPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MFd(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha2(gI1)),
      Sqr(MFd(gI2)))*(Conj(CpbarCha2FdconjUSuPR(gI1,gI2,gO2))*CpbarCha2FdconjUSuPL
      (gI1,gI2,gO1) + Conj(CpbarCha2FdconjUSuPL(gI1,gI2,gO2))*CpbarCha2FdconjUSuPR
      (gI1,gI2,gO1))*MFd(gI2)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSuSvconjUSuconjSv(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,2,SUM(gI2,0,1,(Conj(CpCha1FdconjUSuPL(gI2,gI1,gO2))*
      CpCha1FdconjUSuPL(gI2,gI1,gO1) + Conj(CpCha1FdconjUSuPR(gI2,gI1,gO2))*
      CpCha1FdconjUSuPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MCha1(gI2)))));
   result += -2*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(
      MCha1(gI2)))*(Conj(CpCha1FdconjUSuPR(gI2,gI1,gO2))*CpCha1FdconjUSuPL(gI2,gI1
      ,gO1) + Conj(CpCha1FdconjUSuPL(gI2,gI1,gO2))*CpCha1FdconjUSuPR(gI2,gI1,gO1))
      *MCha1(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,(Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*
      CpChiFuconjUSuPL(gI2,gI1,gO1) + Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*
      CpChiFuconjUSuPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiFuconjUSuPR(gI2,gI1,gO2))*CpChiFuconjUSuPL(gI2,gI1,
      gO1) + Conj(CpChiFuconjUSuPL(gI2,gI1,gO2))*CpChiFuconjUSuPR(gI2,gI1,gO1))*
      MChi(gI2)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUSuconjUSu(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUSuconjUSu(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSuconjHpmconjUSu(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,2,(Conj(CpbarChiFuconjUSuPL(gI1,gI2,gO2))*
      CpbarChiFuconjUSuPL(gI1,gI2,gO1) + Conj(CpbarChiFuconjUSuPR(gI1,gI2,gO2))*
      CpbarChiFuconjUSuPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MFu(gI2)))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MFu(gI2)))*(Conj(CpbarChiFuconjUSuPR(gI1,gI2,gO2))*CpbarChiFuconjUSuPL(gI1,
      gI2,gO1) + Conj(CpbarChiFuconjUSuPL(gI1,gI2,gO2))*CpbarChiFuconjUSuPR(gI1,
      gI2,gO1))*MFu(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,5,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSd(gI2)))*
      Conj(CpSdconjHpmconjUSu(gI2,gI1,gO2))*CpSdconjHpmconjUSu(gI2,gI1,gO1)));
   result += 1.3333333333333333*SUM(gI1,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MphiO))
      *Conj(CpphiOSuconjUSu(gI1,gO2))*CpphiOSuconjUSu(gI1,gO1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSuconjSeconjUSu(gI1,gO1,gI1,
      gO2));
   result += 1.3333333333333333*SUM(gI1,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MsigmaO
      ))*Conj(CpsigmaOSuconjUSu(gI1,gO2))*CpsigmaOSuconjUSu(gI1,gO1));
   result += -SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,
      gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MRh(gI2)))*
      Conj(CpRhSuconjUSu(gI2,gI1,gO2))*CpRhSuconjUSu(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhSuconjUSu(gI2,gI1,gO2))*CpAhSuconjUSu(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhSuconjUSu(gI2,gI1,gO2))*CphhSuconjUSu(gI2,gI1,gO1)));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpbarGluFuconjUSuPL(gI2,gO2))
      *CpbarGluFuconjUSuPL(gI2,gO1) + Conj(CpbarGluFuconjUSuPR(gI2,gO2))*
      CpbarGluFuconjUSuPR(gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFu(gI2))));
   result += 1.3333333333333333*SUM(gI2,0,2,(Conj(CpGluFuconjUSuPL(gI2,gO2))*
      CpGluFuconjUSuPL(gI2,gO1) + Conj(CpGluFuconjUSuPR(gI2,gO2))*CpGluFuconjUSuPR
      (gI2,gO1))*G0(Sqr(p),Sqr(MGlu),Sqr(MFu(gI2))));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFu(
      gI2)))*(Conj(CpbarGluFuconjUSuPR(gI2,gO2))*CpbarGluFuconjUSuPL(gI2,gO1) +
      Conj(CpbarGluFuconjUSuPL(gI2,gO2))*CpbarGluFuconjUSuPR(gI2,gO1))*MFu(gI2));
   result += -2.6666666666666665*MGlu*SUM(gI2,0,2,B0(Sqr(p),Sqr(MGlu),Sqr(MFu(
      gI2)))*(Conj(CpGluFuconjUSuPR(gI2,gO2))*CpGluFuconjUSuPL(gI2,gO1) + Conj(
      CpGluFuconjUSuPL(gI2,gO2))*CpGluFuconjUSuPR(gI2,gO1))*MFu(gI2));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MSRum),Sqr(MSd(gI2)))*Conj(
      CpSdconjSRumconjUSu(gI2,gO2))*CpSdconjSRumconjUSu(gI2,gO1));
   result += SUM(gI2,0,5,B0(Sqr(p),Sqr(MSRdp),Sqr(MSd(gI2)))*Conj(
      CpSRdpSdconjUSu(gI2,gO2))*CpSRdpSdconjUSu(gI2,gO1));
   result += SUM(gI2,0,5,Conj(CpSdconjUSuconjVWm(gI2,gO2))*CpSdconjUSuconjVWm(
      gI2,gO1)*F0(Sqr(p),Sqr(MSd(gI2)),Sqr(MVWm)));
   result += 1.3333333333333333*SUM(gI2,0,5,Conj(CpSuconjUSuVG(gI2,gO2))*
      CpSuconjUSuVG(gI2,gO1)*F0(Sqr(p),Sqr(MSu(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSuconjUSuVP(gI2,gO2))*CpSuconjUSuVP(gI2,gO1)*F0
      (Sqr(p),Sqr(MSu(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSuconjUSuVZ(gI2,gO2))*CpSuconjUSuVZ(gI2,gO1)*F0
      (Sqr(p),Sqr(MSu(gI2)),Sqr(MVZ)));

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

   result += -(A0(Sqr(MSRdp))*CpSRdpUSeconjSRdpconjUSe(gO1,gO2));
   result += -(A0(Sqr(MSRum))*CpSRumUSeconjSRumconjUSe(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUSeconjUSeVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhUSeconjRhconjUSe(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,2,(Conj(CpbarCha1FvconjUSePL(gI1,gI2,gO2))*
      CpbarCha1FvconjUSePL(gI1,gI2,gO1) + Conj(CpbarCha1FvconjUSePR(gI1,gI2,gO2))*
      CpbarCha1FvconjUSePR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MFv(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha1(gI1)),
      Sqr(MFv(gI2)))*(Conj(CpbarCha1FvconjUSePR(gI1,gI2,gO2))*CpbarCha1FvconjUSePL
      (gI1,gI2,gO1) + Conj(CpbarCha1FvconjUSePL(gI1,gI2,gO2))*CpbarCha1FvconjUSePR
      (gI1,gI2,gO1))*MFv(gI2)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUSeSvconjUSeconjSv(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,2,SUM(gI2,0,1,(Conj(CpCha2FvconjUSePL(gI2,gI1,gO2))*
      CpCha2FvconjUSePL(gI2,gI1,gO1) + Conj(CpCha2FvconjUSePR(gI2,gI1,gO2))*
      CpCha2FvconjUSePR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MCha2(gI2)))));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(
      MCha2(gI2)))*(Conj(CpCha2FvconjUSePR(gI2,gI1,gO2))*CpCha2FvconjUSePL(gI2,gI1
      ,gO1) + Conj(CpCha2FvconjUSePL(gI2,gI1,gO2))*CpCha2FvconjUSePR(gI2,gI1,gO1))
      *MCha2(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpHpmSvconjUSe(gI2,gI1,gO2))*CpHpmSvconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,2,SUM(gI2,0,3,(Conj(CpChiFeconjUSePL(gI2,gI1,gO2))*
      CpChiFeconjUSePL(gI2,gI1,gO1) + Conj(CpChiFeconjUSePR(gI2,gI1,gO2))*
      CpChiFeconjUSePR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiFeconjUSePR(gI2,gI1,gO2))*CpChiFeconjUSePL(gI2,gI1,
      gO1) + Conj(CpChiFeconjUSePL(gI2,gI1,gO2))*CpChiFeconjUSePR(gI2,gI1,gO1))*
      MChi(gI2)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUSeconjUSe(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUSeconjUSe(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUSeconjHpmconjUSe(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,2,(Conj(CpbarChiFeconjUSePL(gI1,gI2,gO2))*
      CpbarChiFeconjUSePL(gI1,gI2,gO1) + Conj(CpbarChiFeconjUSePR(gI1,gI2,gO2))*
      CpbarChiFeconjUSePR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MFe(gI2)))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MFe(gI2)))*(Conj(CpbarChiFeconjUSePR(gI1,gI2,gO2))*CpbarChiFeconjUSePL(gI1,
      gI2,gO1) + Conj(CpbarChiFeconjUSePL(gI1,gI2,gO2))*CpbarChiFeconjUSePR(gI1,
      gI2,gO1))*MFe(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdUSeconjSdconjUSe(gI1,gO1,gI1,
      gO2));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeUSeconjSeconjUSe(gI1,gO1,gI1,
      gO2));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUSeSuconjUSeconjSu(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MRh(gI2)))*
      Conj(CpRhSeconjUSe(gI2,gI1,gO2))*CpRhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhSeconjUSe(gI2,gI1,gO2))*CpAhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhSeconjUSe(gI2,gI1,gO2))*CphhSeconjUSe(gI2,gI1,gO1)));
   result += SUM(gI2,0,2,B0(Sqr(p),Sqr(MSRdp),Sqr(MSv(gI2)))*Conj(
      CpSvconjSRdpconjUSe(gI2,gO2))*CpSvconjSRdpconjUSe(gI2,gO1));
   result += SUM(gI2,0,2,Conj(CpSvconjUSeVWm(gI2,gO2))*CpSvconjUSeVWm(gI2,gO1)*
      F0(Sqr(p),Sqr(MSv(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVP(gI2,gO2))*CpSeconjUSeVP(gI2,gO1)*F0
      (Sqr(p),Sqr(MSe(gI2)),0));
   result += SUM(gI2,0,5,Conj(CpSeconjUSeVZ(gI2,gO2))*CpSeconjUSeVZ(gI2,gO1)*F0
      (Sqr(p),Sqr(MSe(gI2)),Sqr(MVZ)));

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
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmgWmUhh(gO1)*
      CpbargWmgWmUhh(gO2));
   result += -(B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*CpbargZgZUhh(gO1)*CpbargZgZUhh(gO2)
      );
   result += B0(Sqr(p),Sqr(MSRdp),Sqr(MSRdp))*Conj(CpSRdpUhhconjSRdp(gO2))*
      CpSRdpUhhconjSRdp(gO1);
   result += -(A0(Sqr(MSRdp))*CpSRdpUhhUhhconjSRdp(gO1,gO2));
   result += B0(Sqr(p),Sqr(MSRum),Sqr(MSRum))*Conj(CpSRumUhhconjSRum(gO2))*
      CpSRumUhhconjSRum(gO1);
   result += -(A0(Sqr(MSRum))*CpSRumUhhUhhconjSRum(gO1,gO2));
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*Conj(CpUhhconjVWmVWm(gO2))*
      CpUhhconjVWmVWm(gO1);
   result += 4*A0(Sqr(MVWm))*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUhhUhhVZVZ(gO1,gO2);
   result += 2*B0(Sqr(p),Sqr(MVZ),Sqr(MVZ))*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1)
      ;
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpUhhUhhRhconjRh(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MRh(gI2)))*
      Conj(CpUhhRhconjRh(gO2,gI2,gI1))*CpUhhRhconjRh(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarCha1Cha1UhhPL(gI1,gI2,gO2))*
      CpbarCha1Cha1UhhPL(gI1,gI2,gO1) + Conj(CpbarCha1Cha1UhhPR(gI1,gI2,gO2))*
      CpbarCha1Cha1UhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MCha1(gI2))))
      );
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarCha2Cha2UhhPL(gI1,gI2,gO2))*
      CpbarCha2Cha2UhhPL(gI1,gI2,gO1) + Conj(CpbarCha2Cha2UhhPR(gI1,gI2,gO2))*
      CpbarCha2Cha2UhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha2(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI1)),
      Sqr(MCha1(gI2)))*(Conj(CpbarCha1Cha1UhhPR(gI1,gI2,gO2))*CpbarCha1Cha1UhhPL(
      gI1,gI2,gO1) + Conj(CpbarCha1Cha1UhhPL(gI1,gI2,gO2))*CpbarCha1Cha1UhhPR(gI1,
      gI2,gO1))*MCha1(gI2)));
   result += -2*SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI1)),
      Sqr(MCha2(gI2)))*(Conj(CpbarCha2Cha2UhhPR(gI1,gI2,gO2))*CpbarCha2Cha2UhhPL(
      gI1,gI2,gO1) + Conj(CpbarCha2Cha2UhhPL(gI1,gI2,gO2))*CpbarCha2Cha2UhhPR(gI1,
      gI2,gO1))*MCha2(gI2)));
   result += 2*SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhUhhconjRh(gI2,gO2,gI1))*CpAhUhhconjRh(gI2,gO1,gI1)));
   result += 2*SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhUhhconjRh(gI2,gO2,gI1))*CphhUhhconjRh(gI2,gO1,gI1)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUhhUhhSvconjSv(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpUhhSvconjSv(gO2,gI2,gI1))*CpUhhSvconjSv(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUhhPL(gI1,gI2,gO2))*
      CpbarFdFdUhhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*
      CpbarFdFdUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUhhPL(gI1,gI2,gO2))*
      CpbarFeFeUhhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*
      CpbarFeFeUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUhhPL(gI1,gI2,gO2))*
      CpbarFuFuUhhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*
      CpbarFuFuUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(
      MFd(gI2)))*(Conj(CpbarFdFdUhhPR(gI1,gI2,gO2))*CpbarFdFdUhhPL(gI1,gI2,gO1) +
      Conj(CpbarFdFdUhhPL(gI1,gI2,gO2))*CpbarFdFdUhhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(
      MFe(gI2)))*(Conj(CpbarFeFeUhhPR(gI1,gI2,gO2))*CpbarFeFeUhhPL(gI1,gI2,gO1) +
      Conj(CpbarFeFeUhhPL(gI1,gI2,gO2))*CpbarFeFeUhhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(
      MFu(gI2)))*(Conj(CpbarFuFuUhhPR(gI1,gI2,gO2))*CpbarFuFuUhhPL(gI1,gI2,gO1) +
      Conj(CpbarFuFuUhhPL(gI1,gI2,gO2))*CpbarFuFuUhhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUhhUhh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUhhUhh(gI1,gI1,gO1,gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpUhhUhhHpmconjHpm(gO1,gO2,gI1,gI1
      ));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))
      *Conj(CpAhAhUhh(gI1,gI2,gO2))*CpAhAhUhh(gI1,gI2,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhhhUhh(gI2,gI1,gO2))*CpAhhhUhh(gI2,gI1,gO1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))
      *Conj(CphhhhUhh(gI1,gI2,gO2))*CphhhhUhh(gI1,gI2,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpUhhHpmconjHpm(gO2,gI2,gI1))*CpUhhHpmconjHpm(gO1,gI2,gI1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,(Conj(CpbarChiChiUhhPL(gI1,gI2,gO2))*
      CpbarChiChiUhhPL(gI1,gI2,gO1) + Conj(CpbarChiChiUhhPR(gI1,gI2,gO2))*
      CpbarChiChiUhhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpbarChiChiUhhPR(gI1,gI2,gO2))*CpbarChiChiUhhPL(gI1,gI2,
      gO1) + Conj(CpbarChiChiUhhPL(gI1,gI2,gO2))*CpbarChiChiUhhPR(gI1,gI2,gO1))*
      MChi(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhUhhSdconjSd(gO1,gO2,gI1,gI1)
      );
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhUhhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhUhhSuconjSu(gO1,gO2,gI1,gI1)
      );
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpUhhSdconjSd(gO2,gI2,gI1))*CpUhhSdconjSd(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpUhhSeconjSe(gO2,gI2,gI1))*CpUhhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpUhhSuconjSu(gO2,gI2,gI1))*CpUhhSuconjSu(gO1,gI2,gI1)));
   result += 2*SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRdp),Sqr(MHpm(gI2)))*Conj(
      CpSRdpUhhHpm(gO2,gI2))*CpSRdpUhhHpm(gO1,gI2));
   result += 2*SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRum),Sqr(MHpm(gI2)))*Conj(
      CpUhhHpmconjSRum(gO2,gI2))*CpUhhHpmconjSRum(gO1,gI2));
   result += SUM(gI2,0,3,Conj(CpAhUhhVZ(gI2,gO2))*CpAhUhhVZ(gI2,gO1)*F0(Sqr(p),
      Sqr(MAh(gI2)),Sqr(MVZ)));
   result += 2*SUM(gI2,0,3,Conj(CpUhhHpmconjVWm(gO2,gI2))*CpUhhHpmconjVWm(gO1,
      gI2)*F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_hh_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = i; k < 4; k++)
         self_energy(i, k) = self_energy_hh_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Ah_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmCgWmCUAh(gO1)*
      CpbargWmCgWmCUAh(gO2));
   result += -(B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*CpbargWmgWmUAh(gO1)*
      CpbargWmgWmUAh(gO2));
   result += B0(Sqr(p),Sqr(MSRdp),Sqr(MSRdp))*Conj(CpSRdpUAhconjSRdp(gO2))*
      CpSRdpUAhconjSRdp(gO1);
   result += -(A0(Sqr(MSRdp))*CpSRdpUAhUAhconjSRdp(gO1,gO2));
   result += B0(Sqr(p),Sqr(MSRum),Sqr(MSRum))*Conj(CpSRumUAhconjSRum(gO2))*
      CpSRumUAhconjSRum(gO1);
   result += -(A0(Sqr(MSRum))*CpSRumUAhUAhconjSRum(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUAhUAhVZVZ(gO1,gO2);
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpUAhUAhRhconjRh(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,1,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MRh(gI2)))*
      Conj(CpUAhRhconjRh(gO2,gI2,gI1))*CpUAhRhconjRh(gO1,gI2,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarCha1Cha1UAhPL(gI1,gI2,gO2))*
      CpbarCha1Cha1UAhPL(gI1,gI2,gO1) + Conj(CpbarCha1Cha1UAhPR(gI1,gI2,gO2))*
      CpbarCha1Cha1UAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MCha1(gI2))))
      );
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpbarCha2Cha2UAhPL(gI1,gI2,gO2))*
      CpbarCha2Cha2UAhPL(gI1,gI2,gO1) + Conj(CpbarCha2Cha2UAhPR(gI1,gI2,gO2))*
      CpbarCha2Cha2UAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha2(gI2))))
      );
   result += -2*SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI1)),
      Sqr(MCha1(gI2)))*(Conj(CpbarCha1Cha1UAhPR(gI1,gI2,gO2))*CpbarCha1Cha1UAhPL(
      gI1,gI2,gO1) + Conj(CpbarCha1Cha1UAhPL(gI1,gI2,gO2))*CpbarCha1Cha1UAhPR(gI1,
      gI2,gO1))*MCha1(gI2)));
   result += -2*SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI1)),
      Sqr(MCha2(gI2)))*(Conj(CpbarCha2Cha2UAhPR(gI1,gI2,gO2))*CpbarCha2Cha2UAhPL(
      gI1,gI2,gO1) + Conj(CpbarCha2Cha2UAhPL(gI1,gI2,gO2))*CpbarCha2Cha2UAhPR(gI1,
      gI2,gO1))*MCha2(gI2)));
   result += 2*SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhUAhconjRh(gI2,gO2,gI1))*CpAhUAhconjRh(gI2,gO1,gI1)));
   result += 2*SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CpUAhhhconjRh(gO2,gI2,gI1))*CpUAhhhconjRh(gO1,gI2,gI1)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUAhUAhSvconjSv(gO1,gO2,gI1,gI1));
   result += SUM(gI1,0,2,SUM(gI2,0,2,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSv(gI2)))*
      Conj(CpUAhSvconjSv(gO2,gI2,gI1))*CpUAhSvconjSv(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFdFdUAhPL(gI1,gI2,gO2))*
      CpbarFdFdUAhPL(gI1,gI2,gO1) + Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*
      CpbarFdFdUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFeFeUAhPL(gI1,gI2,gO2))*
      CpbarFeFeUAhPL(gI1,gI2,gO1) + Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*
      CpbarFeFeUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFuUAhPL(gI1,gI2,gO2))*
      CpbarFuFuUAhPL(gI1,gI2,gO1) + Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*
      CpbarFuFuUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(
      MFd(gI2)))*(Conj(CpbarFdFdUAhPR(gI1,gI2,gO2))*CpbarFdFdUAhPL(gI1,gI2,gO1) +
      Conj(CpbarFdFdUAhPL(gI1,gI2,gO2))*CpbarFdFdUAhPR(gI1,gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(
      MFe(gI2)))*(Conj(CpbarFeFeUAhPR(gI1,gI2,gO2))*CpbarFeFeUAhPL(gI1,gI2,gO1) +
      Conj(CpbarFeFeUAhPL(gI1,gI2,gO2))*CpbarFeFeUAhPR(gI1,gI2,gO1))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(
      MFu(gI2)))*(Conj(CpbarFuFuUAhPR(gI1,gI2,gO2))*CpbarFuFuUAhPL(gI1,gI2,gO1) +
      Conj(CpbarFuFuUAhPL(gI1,gI2,gO2))*CpbarFuFuUAhPR(gI1,gI2,gO1))*MFu(gI2)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUAhUAh(gI1,gI1,gO1,gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpUAhUAhHpmconjHpm(gO1,gO2,gI1,gI1
      ));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2)))
      *Conj(CpAhAhUAh(gI1,gI2,gO2))*CpAhAhUAh(gI1,gI2,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhUAhhh(gI2,gO2,gI1))*CpAhUAhhh(gI2,gO1,gI1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2)))
      *Conj(CpUAhhhhh(gO2,gI1,gI2))*CpUAhhhhh(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpUAhHpmconjHpm(gO2,gI2,gI1))*CpUAhHpmconjHpm(gO1,gI2,gI1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,(Conj(CpbarChiChiUAhPL(gI1,gI2,gO2))*
      CpbarChiChiUAhPL(gI1,gI2,gO1) + Conj(CpbarChiChiUAhPR(gI1,gI2,gO2))*
      CpbarChiChiUAhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpbarChiChiUAhPR(gI1,gI2,gO2))*CpbarChiChiUAhPL(gI1,gI2,
      gO1) + Conj(CpbarChiChiUAhPL(gI1,gI2,gO2))*CpbarChiChiUAhPR(gI1,gI2,gO1))*
      MChi(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUAhUAhSdconjSd(gO1,gO2,gI1,gI1)
      );
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUAhUAhSeconjSe(gO1,gO2,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUAhUAhSuconjSu(gO1,gO2,gI1,gI1)
      );
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))*
      Conj(CpUAhSdconjSd(gO2,gI2,gI1))*CpUAhSdconjSd(gO1,gI2,gI1)));
   result += SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))*
      Conj(CpUAhSeconjSe(gO2,gI2,gI1))*CpUAhSeconjSe(gO1,gI2,gI1)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))*
      Conj(CpUAhSuconjSu(gO2,gI2,gI1))*CpUAhSuconjSu(gO1,gI2,gI1)));
   result += 2*SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRdp),Sqr(MHpm(gI2)))*Conj(
      CpSRdpUAhHpm(gO2,gI2))*CpSRdpUAhHpm(gO1,gI2));
   result += 2*SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRum),Sqr(MHpm(gI2)))*Conj(
      CpUAhHpmconjSRum(gO2,gI2))*CpUAhHpmconjSRum(gO1,gI2));
   result += SUM(gI2,0,3,Conj(CpUAhhhVZ(gO2,gI2))*CpUAhhhVZ(gO1,gI2)*F0(Sqr(p),
      Sqr(Mhh(gI2)),Sqr(MVZ)));
   result += 2*SUM(gI2,0,3,Conj(CpUAhHpmconjVWm(gO2,gI2))*CpUAhHpmconjVWm(gO1,
      gI2)*F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVWm)));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_Ah_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = i; k < 4; k++)
         self_energy(i, k) = self_energy_Ah_1loop(p, i, k);

   Hermitianize(self_energy);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Rh_1loop(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -(A0(Sqr(MSRdp))*CpSRdpURhconjSRdpconjURh(gO1,gO2));
   result += -(A0(Sqr(MSRum))*CpSRumURhconjSRumconjURh(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpURhconjURhconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpURhconjURhVZVZ(gO1,gO2);
   result += Conj(CpSRdpconjURhVWm(gO2))*CpSRdpconjURhVWm(gO1)*F0(Sqr(p),Sqr(
      MSRdp),Sqr(MVWm));
   result += Conj(CpSRumconjURhconjVWm(gO2))*CpSRumconjURhconjVWm(gO1)*F0(Sqr(p
      ),Sqr(MSRum),Sqr(MVWm));
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhURhconjRhconjURh(gI1,gO1,gI1,
      gO2));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(Conj(CpCha1Cha2conjURhPL(gI2,gI1,gO2))*
      CpCha1Cha2conjURhPL(gI2,gI1,gO1) + Conj(CpCha1Cha2conjURhPR(gI2,gI1,gO2))*
      CpCha1Cha2conjURhPR(gI2,gI1,gO1))*G0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha1(gI2)))
      ));
   result += -2*SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI1)),
      Sqr(MCha1(gI2)))*(Conj(CpCha1Cha2conjURhPR(gI2,gI1,gO2))*CpCha1Cha2conjURhPL
      (gI2,gI1,gO1) + Conj(CpCha1Cha2conjURhPL(gI2,gI1,gO2))*CpCha1Cha2conjURhPR(
      gI2,gI1,gO1))*MCha1(gI2)));
   result += SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhRhconjURh(gI2,gI1,gO2))*CpAhRhconjURh(gI2,gI1,gO1)));
   result += SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhRhconjURh(gI2,gI1,gO2))*CphhRhconjURh(gI2,gI1,gO1)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpURhSvconjURhconjSv(gO1,gI1,gO2,
      gI1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhURhconjURh(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhURhconjURh(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmURhconjHpmconjURh(gI1,gO1,gI1
      ,gO2));
   result += SUM(gI1,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MSRum))*Conj(
      CpSRumconjHpmconjURh(gI1,gO2))*CpSRumconjHpmconjURh(gI1,gO1));
   result += 0.25*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MAh(gI1)),Sqr(MAh(gI2))
      )*Conj(CpAhAhconjURh(gI1,gI2,gO2))*CpAhAhconjURh(gI1,gI2,gO1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(MAh(gI2)))
      *Conj(CpAhhhconjURh(gI2,gI1,gO2))*CpAhhhconjURh(gI2,gI1,gO1)));
   result += 0.25*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(Mhh(gI1)),Sqr(Mhh(gI2))
      )*Conj(CphhhhconjURh(gI1,gI2,gO2))*CphhhhconjURh(gI1,gI2,gO1)));
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)
      ))*Conj(CpHpmconjHpmconjURh(gI2,gI1,gO2))*CpHpmconjHpmconjURh(gI2,gI1,gO1)))
      ;
   result += 0.5*SUM(gI1,0,3,SUM(gI2,0,3,(Conj(CpChiChiconjURhPL(gI1,gI2,gO2))*
      CpChiChiconjURhPL(gI1,gI2,gO1) + Conj(CpChiChiconjURhPR(gI1,gI2,gO2))*
      CpChiChiconjURhPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))));
   result += -SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MChi(gI2)))*(Conj(CpChiChiconjURhPR(gI1,gI2,gO2))*CpChiChiconjURhPL(gI1,gI2,
      gO1) + Conj(CpChiChiconjURhPL(gI1,gI2,gO2))*CpChiChiconjURhPR(gI1,gI2,gO1))*
      MChi(gI2)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpURhSdconjURhconjSd(gO1,gI1,gO2,
      gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpURhSeconjURhconjSe(gO1,gI1,gO2,
      gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpURhSuconjURhconjSu(gO1,gI1,gO2,
      gI1));
   result += 1.5*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))
      *Conj(CpSdconjURhconjSd(gI2,gO2,gI1))*CpSdconjURhconjSd(gI2,gO1,gI1)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSe(gI1)),Sqr(MSe(gI2)))
      *Conj(CpSeconjURhconjSe(gI2,gO2,gI1))*CpSeconjURhconjSe(gI2,gO1,gI1)));
   result += 1.5*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))
      *Conj(CpSuconjURhconjSu(gI2,gO2,gI1))*CpSuconjURhconjSu(gI2,gO1,gI1)));
   result += SUM(gI2,0,1,Conj(CpRhconjURhVZ(gI2,gO2))*CpRhconjURhVZ(gI2,gO1)*F0
      (Sqr(p),Sqr(MRh(gI2)),Sqr(MVZ)));
   result += SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRdp),Sqr(MHpm(gI2)))*Conj(
      CpSRdpHpmconjURh(gI2,gO2))*CpSRdpHpmconjURh(gI2,gO1));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Rh_1loop(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = i; k < 2; k++)
         self_energy(i, k) = self_energy_Rh_1loop(p, i, k);

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
   result += 4*B0(Sqr(p),0,Sqr(MVWm))*Conj(CpconjUHpmVPVWm(gO2))*
      CpconjUHpmVPVWm(gO1);
   result += 4*B0(Sqr(p),Sqr(MVWm),Sqr(MVZ))*Conj(CpconjUHpmVWmVZ(gO2))*
      CpconjUHpmVWmVZ(gO1);
   result += -(A0(Sqr(MSRdp))*CpSRdpUHpmconjSRdpconjUHpm(gO1,gO2));
   result += -(A0(Sqr(MSRum))*CpSRumUHpmconjSRumconjUHpm(gO1,gO2));
   result += 4*A0(Sqr(MVWm))*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(Sqr(MVZ))*CpUHpmconjUHpmVZVZ(gO1,gO2);
   result += SUM(gI1,0,1,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MSRum))*Conj(
      CpSRumconjUHpmconjRh(gO2,gI1))*CpSRumconjUHpmconjRh(gO1,gI1));
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpUHpmRhconjUHpmconjRh(gO1,gI1,gO2,
      gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,3,B0(Sqr(p),Sqr(MRh(gI1)),Sqr(MHpm(gI2)))*
      Conj(CpHpmconjUHpmconjRh(gI2,gO2,gI1))*CpHpmconjUHpmconjRh(gI2,gO1,gI1)));
   result += SUM(gI1,0,1,SUM(gI2,0,3,(Conj(CpbarCha1ChiconjUHpmPL(gI1,gI2,gO2))
      *CpbarCha1ChiconjUHpmPL(gI1,gI2,gO1) + Conj(CpbarCha1ChiconjUHpmPR(gI1,gI2,
      gO2))*CpbarCha1ChiconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MChi(gI2)))));
   result += -2*SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MCha1(gI1)),
      Sqr(MChi(gI2)))*(Conj(CpbarCha1ChiconjUHpmPR(gI1,gI2,gO2))*
      CpbarCha1ChiconjUHpmPL(gI1,gI2,gO1) + Conj(CpbarCha1ChiconjUHpmPL(gI1,gI2,
      gO2))*CpbarCha1ChiconjUHpmPR(gI1,gI2,gO1))*MChi(gI2)));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUHpmSvconjUHpmconjSv(gO1,gI1,gO2,
      gI1));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFuFdconjUHpmPL(gI1,gI2,gO2))*
      CpbarFuFdconjUHpmPL(gI1,gI2,gO1) + Conj(CpbarFuFdconjUHpmPR(gI1,gI2,gO2))*
      CpbarFuFdconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(Conj(CpbarFvFeconjUHpmPL(gI1,gI2,gO2))*
      CpbarFvFeconjUHpmPL(gI1,gI2,gO1) + Conj(CpbarFvFeconjUHpmPR(gI1,gI2,gO2))*
      CpbarFvFeconjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(
      MFd(gI2)))*(Conj(CpbarFuFdconjUHpmPR(gI1,gI2,gO2))*CpbarFuFdconjUHpmPL(gI1,
      gI2,gO1) + Conj(CpbarFuFdconjUHpmPL(gI1,gI2,gO2))*CpbarFuFdconjUHpmPR(gI1,
      gI2,gO1))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(
      MFe(gI2)))*(Conj(CpbarFvFeconjUHpmPR(gI1,gI2,gO2))*CpbarFvFeconjUHpmPL(gI1,
      gI2,gO1) + Conj(CpbarFvFeconjUHpmPL(gI1,gI2,gO2))*CpbarFvFeconjUHpmPR(gI1,
      gI2,gO1))*MFe(gI2)));
   result += SUM(gI1,0,2,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSv(gI1)),Sqr(MSe(gI2)))*
      Conj(CpSeconjUHpmconjSv(gI2,gO2,gI1))*CpSeconjUHpmconjSv(gI2,gO1,gI1)));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUHpmconjUHpm(gI1,gI1,gO1,
      gO2));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUHpmconjUHpm(gI1,gI1,gO1,
      gO2));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmUHpmconjHpmconjUHpm(gI1,gO1,
      gI1,gO2));
   result += SUM(gI1,0,3,SUM(gI2,0,1,(Conj(CpbarChiCha2conjUHpmPL(gI1,gI2,gO2))
      *CpbarChiCha2conjUHpmPL(gI1,gI2,gO1) + Conj(CpbarChiCha2conjUHpmPR(gI1,gI2,
      gO2))*CpbarChiCha2conjUHpmPR(gI1,gI2,gO1))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha2(gI2)))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha2(gI2)))*(Conj(CpbarChiCha2conjUHpmPR(gI1,gI2,gO2))*
      CpbarChiCha2conjUHpmPL(gI1,gI2,gO1) + Conj(CpbarChiCha2conjUHpmPL(gI1,gI2,
      gO2))*CpbarChiCha2conjUHpmPR(gI1,gI2,gO1))*MCha2(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(MAh(gI2)))*
      Conj(CpAhHpmconjUHpm(gI2,gI1,gO2))*CpAhHpmconjUHpm(gI2,gI1,gO1)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MHpm(gI1)),Sqr(Mhh(gI2)))*
      Conj(CphhHpmconjUHpm(gI2,gI1,gO2))*CphhHpmconjUHpm(gI2,gI1,gO1)));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUHpmSdconjUHpmconjSd(gO1,gI1,
      gO2,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUHpmSeconjUHpmconjSe(gO1,gI1,gO2,
      gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUHpmSuconjUHpmconjSu(gO1,gI1,
      gO2,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,B0(Sqr(p),Sqr(MSu(gI1)),Sqr(MSd(gI2)))*
      Conj(CpSdconjUHpmconjSu(gI2,gO2,gI1))*CpSdconjUHpmconjSu(gI2,gO1,gI1)));
   result += SUM(gI2,0,1,B0(Sqr(p),Sqr(MSRdp),Sqr(MRh(gI2)))*Conj(
      CpRhconjSRdpconjUHpm(gI2,gO2))*CpRhconjSRdpconjUHpm(gI2,gO1));
   result += SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRdp),Sqr(MAh(gI2)))*Conj(
      CpAhconjSRdpconjUHpm(gI2,gO2))*CpAhconjSRdpconjUHpm(gI2,gO1));
   result += SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRdp),Sqr(Mhh(gI2)))*Conj(
      CphhconjSRdpconjUHpm(gI2,gO2))*CphhconjSRdpconjUHpm(gI2,gO1));
   result += SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRum),Sqr(MAh(gI2)))*Conj(
      CpSRumAhconjUHpm(gI2,gO2))*CpSRumAhconjUHpm(gI2,gO1));
   result += SUM(gI2,0,3,B0(Sqr(p),Sqr(MSRum),Sqr(Mhh(gI2)))*Conj(
      CpSRumhhconjUHpm(gI2,gO2))*CpSRumhhconjUHpm(gI2,gO1));
   result += SUM(gI2,0,3,Conj(CpAhconjUHpmVWm(gI2,gO2))*CpAhconjUHpmVWm(gI2,gO1
      )*F0(Sqr(p),Sqr(MAh(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,3,Conj(CphhconjUHpmVWm(gI2,gO2))*CphhconjUHpmVWm(gI2,gO1
      )*F0(Sqr(p),Sqr(Mhh(gI2)),Sqr(MVWm)));
   result += SUM(gI2,0,3,Conj(CpHpmconjUHpmVP(gI2,gO2))*CpHpmconjUHpmVP(gI2,gO1
      )*F0(Sqr(p),Sqr(MHpm(gI2)),0));
   result += SUM(gI2,0,3,Conj(CpHpmconjUHpmVZ(gI2,gO2))*CpHpmconjUHpmVZ(gI2,gO1
      )*F0(Sqr(p),Sqr(MHpm(gI2)),Sqr(MVZ)));

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

std::complex<double> CLASSNAME::self_energy_SRdp_1loop(double p ) const
{
   std::complex<double> result;

   result += 4*A0(Sqr(MVWm))*CpSRdpconjSRdpconjVWmVWm();
   result += 2*A0(Sqr(MVZ))*CpSRdpconjSRdpVZVZ();
   result += -(A0(Sqr(MSRdp))*CpSRdpSRdpconjSRdpconjSRdp());
   result += -(A0(Sqr(MSRum))*CpSRdpSRumconjSRdpconjSRum());
   result += AbsSqr(CpSRdpconjSRdpVP())*F0(Sqr(p),Sqr(MSRdp),0);
   result += AbsSqr(CpSRdpconjSRdpVZ())*F0(Sqr(p),Sqr(MSRdp),Sqr(MVZ));
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpSRdpRhconjSRdpconjRh(gI1,gI1));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSRdpSvconjSRdpconjSv(gI1,gI1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpSRdpAhAhconjSRdp(gI1,gI1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CpSRdphhhhconjSRdp(gI1,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpSRdpHpmconjSRdpconjHpm(gI1,gI1))
      ;
   result += SUM(gI1,0,3,SUM(gI2,0,1,AbsSqr(CpRhconjSRdpconjHpm(gI2,gI1))*B0(
      Sqr(p),Sqr(MHpm(gI1)),Sqr(MRh(gI2)))));
   result += SUM(gI1,0,3,SUM(gI2,0,1,(AbsSqr(CpCha1ChiconjSRdpPL(gI2,gI1)) +
      AbsSqr(CpCha1ChiconjSRdpPR(gI2,gI1)))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha1(gI2
      )))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha1(gI2)))*(Conj(CpCha1ChiconjSRdpPR(gI2,gI1))*CpCha1ChiconjSRdpPL(gI2,gI1
      ) + Conj(CpCha1ChiconjSRdpPL(gI2,gI1))*CpCha1ChiconjSRdpPR(gI2,gI1))*MCha1(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpAhconjSRdpconjHpm(gI2,gI1))*B0(
      Sqr(p),Sqr(MHpm(gI1)),Sqr(MAh(gI2)))));
   result += SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CphhconjSRdpconjHpm(gI2,gI1))*B0(
      Sqr(p),Sqr(MHpm(gI1)),Sqr(Mhh(gI2)))));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSRdpSdconjSRdpconjSd(gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSRdpSeconjSRdpconjSe(gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSRdpSuconjSRdpconjSu(gI1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpSvconjSRdpconjSe(gI2,gI1))*B0(Sqr
      (p),Sqr(MSe(gI1)),Sqr(MSv(gI2)))));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSRdpconjSd(gI2,gI1))*B0(
      Sqr(p),Sqr(MSd(gI1)),Sqr(MSu(gI2)))));
   result += SUM(gI2,0,1,AbsSqr(CpRhconjSRdpconjVWm(gI2))*F0(Sqr(p),Sqr(MRh(gI2
      )),Sqr(MVWm)));
   result += SUM(gI2,0,3,AbsSqr(CpSRdpAhconjSRdp(gI2))*B0(Sqr(p),Sqr(MSRdp),Sqr
      (MAh(gI2))));
   result += SUM(gI2,0,3,AbsSqr(CpSRdphhconjSRdp(gI2))*B0(Sqr(p),Sqr(MSRdp),Sqr
      (Mhh(gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_SRum_1loop(double p ) const
{
   std::complex<double> result;

   result += -(A0(Sqr(MSRdp))*CpSRdpSRumconjSRdpconjSRum());
   result += 4*A0(Sqr(MVWm))*CpSRumconjSRumconjVWmVWm();
   result += 2*A0(Sqr(MVZ))*CpSRumconjSRumVZVZ();
   result += -(A0(Sqr(MSRum))*CpSRumSRumconjSRumconjSRum());
   result += AbsSqr(CpSRumconjSRumVP())*F0(Sqr(p),Sqr(MSRum),0);
   result += AbsSqr(CpSRumconjSRumVZ())*F0(Sqr(p),Sqr(MSRum),Sqr(MVZ));
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpSRumRhconjSRumconjRh(gI1,gI1));
   result += SUM(gI1,0,1,SUM(gI2,0,3,AbsSqr(CpHpmRhconjSRum(gI2,gI1))*B0(Sqr(p)
      ,Sqr(MRh(gI1)),Sqr(MHpm(gI2)))));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSRumSvconjSRumconjSv(gI1,gI1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpSRumAhAhconjSRum(gI1,gI1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CpSRumhhhhconjSRum(gI1,gI1));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpSRumHpmconjSRumconjHpm(gI1,gI1))
      ;
   result += SUM(gI1,0,3,SUM(gI2,0,1,(AbsSqr(CpCha2ChiconjSRumPL(gI2,gI1)) +
      AbsSqr(CpCha2ChiconjSRumPR(gI2,gI1)))*G0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha2(gI2
      )))));
   result += -2*SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(
      MCha2(gI2)))*(Conj(CpCha2ChiconjSRumPR(gI2,gI1))*CpCha2ChiconjSRumPL(gI2,gI1
      ) + Conj(CpCha2ChiconjSRumPL(gI2,gI1))*CpCha2ChiconjSRumPR(gI2,gI1))*MCha2(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpAhHpmconjSRum(gI2,gI1))*B0(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(MAh(gI2)))));
   result += SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CphhHpmconjSRum(gI2,gI1))*B0(Sqr(p)
      ,Sqr(MHpm(gI1)),Sqr(Mhh(gI2)))));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSRumSdconjSRumconjSd(gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSRumSeconjSRumconjSe(gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSRumSuconjSRumconjSu(gI1,gI1));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSRumconjSu(gI2,gI1))*B0(
      Sqr(p),Sqr(MSu(gI1)),Sqr(MSd(gI2)))));
   result += SUM(gI2,0,1,AbsSqr(CpRhconjSRumVWm(gI2))*F0(Sqr(p),Sqr(MRh(gI2)),
      Sqr(MVWm)));
   result += SUM(gI2,0,3,AbsSqr(CpSRumAhconjSRum(gI2))*B0(Sqr(p),Sqr(MSRum),Sqr
      (MAh(gI2))));
   result += SUM(gI2,0,3,AbsSqr(CpSRumhhconjSRum(gI2))*B0(Sqr(p),Sqr(MSRum),Sqr
      (Mhh(gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_sigmaO_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MphiO))*CpsigmaOsigmaOphiOphiO();
   result += 3*AbsSqr(CpsigmaOsigmaOVG())*F0(Sqr(p),Sqr(MsigmaO),0);
   result += 3*(AbsSqr(CpbarGluGlusigmaOPL()) + AbsSqr(CpbarGluGlusigmaOPR()))*
      G0(Sqr(p),Sqr(MGlu),Sqr(MGlu));
   result += -6*B0(Sqr(p),Sqr(MGlu),Sqr(MGlu))*(Conj(CpbarGluGlusigmaOPR())*
      CpbarGluGlusigmaOPL() + Conj(CpbarGluGlusigmaOPL())*CpbarGluGlusigmaOPR())*
      Sqr(MGlu);
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpsigmaOSdconjSd(gI2,gI1))*B0(
      Sqr(p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpsigmaOSuconjSu(gI2,gI1))*B0(
      Sqr(p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_phiO_1loop(double p ) const
{
   std::complex<double> result;

   result += -0.5*A0(Sqr(MsigmaO))*CpphiOphiOsigmaOsigmaO();
   result += 3*AbsSqr(CpphiOphiOVG())*F0(Sqr(p),Sqr(MphiO),0);
   result += 3*(AbsSqr(CpbarGluGluphiOPL()) + AbsSqr(CpbarGluGluphiOPR()))*G0(
      Sqr(p),Sqr(MGlu),Sqr(MGlu));
   result += -6*B0(Sqr(p),Sqr(MGlu),Sqr(MGlu))*(Conj(CpbarGluGluphiOPR())*
      CpbarGluGluphiOPL() + Conj(CpbarGluGluphiOPL())*CpbarGluGluphiOPR())*Sqr(
      MGlu);
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpphiOSdconjSd(gI2,gI1))*B0(Sqr
      (p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpphiOSuconjSu(gI2,gI1))*B0(Sqr
      (p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VG_1loop(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpbargGgGVG())*B00(Sqr(p),Sqr(MVG),Sqr(MVG));
   result += -6*AbsSqr(CpphiOphiOVG())*B00(Sqr(p),Sqr(MphiO),Sqr(MphiO));
   result += 499.5*A0(Sqr(MphiO))*CpphiOphiOVGVG();
   result += -6*AbsSqr(CpsigmaOsigmaOVG())*B00(Sqr(p),Sqr(MsigmaO),Sqr(MsigmaO)
      );
   result += 499.5*A0(Sqr(MsigmaO))*CpsigmaOsigmaOVGVG();
   result += -3*AbsSqr(CpVGVGVG())*(5*B00(Sqr(p),0,0) + 2*B0(Sqr(p),0,0)*Sqr(p)
      );
   result += 0;
   result += 3*((AbsSqr(CpbarGluGluVGPL()) + AbsSqr(CpbarGluGluVGPR()))*H0(Sqr(
      p),Sqr(MGlu),Sqr(MGlu)) + 4*B0(Sqr(p),Sqr(MGlu),Sqr(MGlu))*Re(Conj(
      CpbarGluGluVGPL())*CpbarGluGluVGPR())*Sqr(MGlu));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVGPL(gI1,gI2)) +
      AbsSqr(CpbarFdFdVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*
      B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(
      CpbarFdFdVGPL(gI1,gI2))*CpbarFdFdVGPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVGPL(gI1,gI2)) +
      AbsSqr(CpbarFuFuVGPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*
      B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(
      CpbarFuFuVGPL(gI1,gI2))*CpbarFuFuVGPR(gI1,gI2))));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVGVG(gI1,gI1));
   result += 999*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVGVG(gI1,gI1));
   result += -2*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVG(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -2*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVG(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSu(gI1)),Sqr(MSu(gI2)))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VP_1loop(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWmCgWmCVP())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += AbsSqr(CpbargWmgWmVP())*B00(Sqr(p),Sqr(MVWm),Sqr(MVWm));
   result += -(A0(Sqr(MVWm))*(CpconjVWmVPVPVWm1() + CpconjVWmVPVPVWm2() + 4*
      CpconjVWmVPVPVWm3()));
   result += -4*AbsSqr(CpSRdpconjSRdpVP())*B00(Sqr(p),Sqr(MSRdp),Sqr(MSRdp));
   result += A0(Sqr(MSRdp))*CpSRdpconjSRdpVPVP();
   result += -4*AbsSqr(CpSRumconjSRumVP())*B00(Sqr(p),Sqr(MSRum),Sqr(MSRum));
   result += A0(Sqr(MSRum))*CpSRumconjSRumVPVP();
   result += -2*AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm
      ),Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarCha1Cha1VPPL(gI1,gI2)) +
      AbsSqr(CpbarCha1Cha1VPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MCha1(gI2)
      )) + 4*B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MCha1(gI2)))*MCha1(gI1)*MCha1(gI2)*Re(
      Conj(CpbarCha1Cha1VPPL(gI1,gI2))*CpbarCha1Cha1VPPR(gI1,gI2))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarCha2Cha2VPPL(gI1,gI2)) +
      AbsSqr(CpbarCha2Cha2VPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha2(gI2)
      )) + 4*B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha2(gI2)))*MCha2(gI1)*MCha2(gI2)*Re(
      Conj(CpbarCha2Cha2VPPL(gI1,gI2))*CpbarCha2Cha2VPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVPPL(gI1,gI2)) + AbsSqr
      (CpbarFdFdVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(
      p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVPPL(gI1,
      gI2))*CpbarFdFdVPPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVPPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVPPL(gI1,
      gI2))*CpbarFeFeVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVPPL(gI1,gI2)) + AbsSqr
      (CpbarFuFuVPPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVPPL(gI1,
      gI2))*CpbarFuFuVPPR(gI1,gI2))));
   result += SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVPVP(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpHpmconjHpmVP(gI2,gI1))*B00(Sqr
      (p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVPVP(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVPVP(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVPVP(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVP(gI2,gI1))*B00(Sqr(
      p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVP(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVP(gI2,gI1))*B00(Sqr(
      p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
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
   result += -4*AbsSqr(CpSRdpconjSRdpVZ())*B00(Sqr(p),Sqr(MSRdp),Sqr(MSRdp));
   result += A0(Sqr(MSRdp))*CpSRdpconjSRdpVZVZ();
   result += -4*AbsSqr(CpSRumconjSRumVZ())*B00(Sqr(p),Sqr(MSRum),Sqr(MSRum));
   result += A0(Sqr(MSRum))*CpSRumconjSRumVZVZ();
   result += -2*AbsSqr(CpconjVWmVWmVZ())*(A0(Sqr(MVWm)) + 5*B00(Sqr(p),Sqr(MVWm
      ),Sqr(MVWm)) + B0(Sqr(p),Sqr(MVWm),Sqr(MVWm))*(Sqr(MVWm) + 2*Sqr(p)));
   result += SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhconjRhVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,1,SUM(gI2,0,1,AbsSqr(CpRhconjRhVZ(gI2,gI1))*B00(Sqr(p
      ),Sqr(MRh(gI1)),Sqr(MRh(gI2)))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarCha1Cha1VZPL(gI1,gI2)) +
      AbsSqr(CpbarCha1Cha1VZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MCha1(gI2)
      )) + 4*B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MCha1(gI2)))*MCha1(gI1)*MCha1(gI2)*Re(
      Conj(CpbarCha1Cha1VZPL(gI1,gI2))*CpbarCha1Cha1VZPR(gI1,gI2))));
   result += SUM(gI1,0,1,SUM(gI2,0,1,(AbsSqr(CpbarCha2Cha2VZPL(gI1,gI2)) +
      AbsSqr(CpbarCha2Cha2VZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha2(gI2)
      )) + 4*B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MCha2(gI2)))*MCha2(gI1)*MCha2(gI2)*Re(
      Conj(CpbarCha2Cha2VZPL(gI1,gI2))*CpbarCha2Cha2VZPR(gI1,gI2))));
   result += SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvconjSvVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,2,AbsSqr(CpSvconjSvVZ(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSv(gI1)),Sqr(MSv(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZPL(gI1,gI2)) + AbsSqr
      (CpbarFdFdVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(Sqr(
      p),Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZPL(gI1,
      gI2))*CpbarFdFdVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVZPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVZPL(gI1,
      gI2))*CpbarFeFeVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVZPL(gI1,gI2)) + AbsSqr
      (CpbarFuFuVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(Sqr(
      p),Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVZPL(gI1,
      gI2))*CpbarFuFuVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFvVZPL(gI1,gI2)) + AbsSqr(
      CpbarFvFvVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFv(gI2))) + 4*B0(Sqr(p
      ),Sqr(MFv(gI1)),Sqr(MFv(gI2)))*MFv(gI1)*MFv(gI2)*Re(Conj(CpbarFvFvVZPL(gI1,
      gI2))*CpbarFvFvVZPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhVZVZ(gI1,gI1));
   result += 0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhVZVZ(gI1,gI1));
   result += SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmVZVZ(gI1,gI1));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpAhhhVZ(gI2,gI1))*B00(Sqr(p),
      Sqr(MAh(gI2)),Sqr(Mhh(gI1)))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpHpmconjHpmVZ(gI2,gI1))*B00(Sqr
      (p),Sqr(MHpm(gI1)),Sqr(MHpm(gI2)))));
   result += SUM(gI1,0,3,SUM(gI2,0,3,(AbsSqr(CpbarChiChiVZPL(gI1,gI2)) + AbsSqr
      (CpbarChiChiVZPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2))) + 4*B0(
      Sqr(p),Sqr(MChi(gI1)),Sqr(MChi(gI2)))*MChi(gI1)*MChi(gI2)*Re(Conj(
      CpbarChiChiVZPL(gI1,gI2))*CpbarChiChiVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdVZVZ(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeVZVZ(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuVZVZ(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSdVZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSd(gI1)),Sqr(MSd(gI2)))));
   result += -4*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSeconjSeVZ(gI2,gI1))*B00(Sqr(p
      ),Sqr(MSe(gI1)),Sqr(MSe(gI2)))));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSuconjSuVZ(gI2,gI1))*B00(Sqr(
      p),Sqr(MSu(gI1)),Sqr(MSu(gI2)))));
   result += SUM(gI2,0,3,AbsSqr(CphhVZVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(Mhh(gI2))
      ));
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
   result += A0(Sqr(MSRdp))*CpSRdpconjSRdpconjVWmVWm();
   result += A0(Sqr(MSRum))*CpSRumconjSRumconjVWmVWm();
   result += -(AbsSqr(CpconjVWmVPVWm())*(A0(Sqr(MVWm)) + 10*B00(Sqr(p),Sqr(MVWm
      ),0) + B0(Sqr(p),Sqr(MVWm),0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVWmVZ())*(A0(Sqr(MVWm)) + A0(Sqr(MVZ)) + 10*B00(
      Sqr(p),Sqr(MVZ),Sqr(MVWm)) + B0(Sqr(p),Sqr(MVZ),Sqr(MVWm))*(Sqr(MVWm) + Sqr(
      MVZ) + 4*Sqr(p))));
   result += SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpRhconjRhconjVWmVWm(gI1,gI1));
   result += -4*SUM(gI1,0,1,AbsSqr(CpSRumconjRhconjVWm(gI1))*B00(Sqr(p),Sqr(
      MSRum),Sqr(MRh(gI1))));
   result += SUM(gI1,0,1,SUM(gI2,0,3,(AbsSqr(CpbarCha1ChiconjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarCha1ChiconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MChi(
      gI2))) + 4*B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MChi(gI2)))*MCha1(gI1)*MChi(gI2)*Re
      (Conj(CpbarCha1ChiconjVWmPL(gI1,gI2))*CpbarCha1ChiconjVWmPR(gI1,gI2))));
   result += SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpSvconjSvconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFdconjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarFuFdconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))
      + 4*B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MFd(gI2)))*MFd(gI2)*MFu(gI1)*Re(Conj(
      CpbarFuFdconjVWmPL(gI1,gI2))*CpbarFuFdconjVWmPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFeconjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarFvFeconjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))
      + 4*B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MFe(gI2)))*MFe(gI2)*MFv(gI1)*Re(Conj(
      CpbarFvFeconjVWmPL(gI1,gI2))*CpbarFvFeconjVWmPR(gI1,gI2))));
   result += -4*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpSeconjSvconjVWm(gI2,gI1))*B00(
      Sqr(p),Sqr(MSe(gI2)),Sqr(MSv(gI1)))));
   result += 0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhconjVWmVWm(gI1,gI1));
   result += 0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpHpmconjHpmconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,1,(AbsSqr(CpbarChiCha2conjVWmPL(gI1,gI2)) +
      AbsSqr(CpbarChiCha2conjVWmPR(gI1,gI2)))*H0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha2(
      gI2))) + 4*B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MCha2(gI2)))*MCha2(gI2)*MChi(gI1)*Re
      (Conj(CpbarChiCha2conjVWmPL(gI1,gI2))*CpbarChiCha2conjVWmPR(gI1,gI2))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CpAhHpmconjVWm(gI2,gI1))*B00(Sqr
      (p),Sqr(MAh(gI2)),Sqr(MHpm(gI1)))));
   result += -4*SUM(gI1,0,3,SUM(gI2,0,3,AbsSqr(CphhHpmconjVWm(gI2,gI1))*B00(Sqr
      (p),Sqr(Mhh(gI2)),Sqr(MHpm(gI1)))));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpSdconjSdconjVWmVWm(gI1,gI1));
   result += SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpSeconjSeconjVWmVWm(gI1,gI1));
   result += 3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpSuconjSuconjVWmVWm(gI1,gI1));
   result += -12*SUM(gI1,0,5,SUM(gI2,0,5,AbsSqr(CpSdconjSuconjVWm(gI2,gI1))*B00
      (Sqr(p),Sqr(MSd(gI2)),Sqr(MSu(gI1)))));
   result += -4*SUM(gI2,0,1,AbsSqr(CpRhconjSRdpconjVWm(gI2))*B00(Sqr(p),Sqr(MRh
      (gI2)),Sqr(MSRdp)));
   result += SUM(gI2,0,3,AbsSqr(CphhconjVWmVWm(gI2))*B0(Sqr(p),Sqr(MVWm),Sqr(
      Mhh(gI2))));
   result += SUM(gI2,0,3,AbsSqr(CpHpmconjVWmVP(gI2))*B0(Sqr(p),0,Sqr(MHpm(gI2))
      ));
   result += SUM(gI2,0,3,AbsSqr(CpHpmconjVWmVZ(gI2))*B0(Sqr(p),Sqr(MVZ),Sqr(
      MHpm(gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSRdp))*Conj(
      CpbarCha1barUChiSRdpPL(gI1,gO2))*CpbarCha1barUChiSRdpPR(gI1,gO1)*MCha1(gI1))
      ;
   result += SUM(gI1,0,1,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSRum))*Conj(
      CpbarCha2barUChiSRumPL(gI1,gO2))*CpbarCha2barUChiSRumPR(gI1,gO1)*MCha2(gI1))
      ;
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(
      gI2)))*Conj(CpbarUChibarFvSvPL(gO2,gI1,gI2))*CpbarUChibarFvSvPR(gO1,gI1,gI2)
      ));
   result += 3*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd
      (gI2)))*Conj(CpbarUChibarFdSdPL(gO2,gI1,gI2))*CpbarUChibarFdSdPR(gO1,gI1,gI2
      )));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(
      gI2)))*Conj(CpbarUChibarFeSePL(gO2,gI1,gI2))*CpbarUChibarFeSePR(gO1,gI1,gI2)
      ));
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu
      (gI2)))*Conj(CpbarUChibarFuSuPL(gO2,gI1,gI2))*CpbarUChibarFuSuPR(gO1,gI1,gI2
      )));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MRh
      (gI2)))*Conj(CpbarChibarUChiRhPL(gI1,gO2,gI2))*CpbarChibarUChiRhPR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUChiCha1HpmPL(gO2,gI2,gI1))*CpbarUChiCha1HpmPR(gO1,gI2,gI1)*MCha1(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUChiCha2conjHpmPL(gO2,gI2,gI1))*CpbarUChiCha2conjHpmPR(gO1,gI2,gI1
      )*MCha2(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh
      (gI2)))*Conj(CpbarUChiChiAhPL(gO2,gI1,gI2))*CpbarUChiChiAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUChiChihhPL(gO2,gI2,gI1))*CpbarUChiChihhPR(gO1,gI2,gI1)*MChi(gI2))
      );
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUChiFdconjSdPL(gO2,gI2,gI1))*CpbarUChiFdconjSdPR(gO1,gI2,gI1)*MFd(
      gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUChiFeconjSePL(gO2,gI2,gI1))*CpbarUChiFeconjSePR(gO1,gI2,gI1)*MFe(
      gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUChiFuconjSuPL(gO2,gI2,gI1))*CpbarUChiFuconjSuPR(gO1,gI2,gI1)*MFu(
      gI2)));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MVWm))*Conj(
      CpbarUChiCha1VWmPR(gO2,gI2))*CpbarUChiCha1VWmPL(gO1,gI2)*MCha1(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),Sqr(MVWm))*Conj(
      CpbarUChiCha2conjVWmPR(gO2,gI2))*CpbarUChiCha2conjVWmPL(gO1,gI2)*MCha2(gI2))
      ;
   result += -4*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(
      CpbarUChiChiVZPR(gO2,gI2))*CpbarUChiChiVZPL(gO1,gI2)*MChi(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_Chi_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++)
         self_energy(i, k) = self_energy_Chi_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSRdp))*Conj(
      CpbarCha1barUChiSRdpPR(gI1,gO2))*CpbarCha1barUChiSRdpPR(gI1,gO1));
   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSRum))*Conj(
      CpbarCha2barUChiSRumPR(gI1,gO2))*CpbarCha2barUChiSRumPR(gI1,gO1));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2))
      )*Conj(CpbarUChibarFvSvPR(gO2,gI1,gI2))*CpbarUChibarFvSvPR(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2))
      )*Conj(CpbarUChibarFdSdPR(gO2,gI1,gI2))*CpbarUChibarFdSdPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2))
      )*Conj(CpbarUChibarFeSePR(gO2,gI1,gI2))*CpbarUChibarFeSePR(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2))
      )*Conj(CpbarUChibarFuSuPR(gO2,gI1,gI2))*CpbarUChibarFuSuPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MRh(gI2)
      ))*Conj(CpbarChibarUChiRhPR(gI1,gO2,gI2))*CpbarChibarUChiRhPR(gI1,gO1,gI2)))
      ;
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MHpm(
      gI1)))*Conj(CpbarUChiCha1HpmPR(gO2,gI2,gI1))*CpbarUChiCha1HpmPR(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MHpm(
      gI1)))*Conj(CpbarUChiCha2conjHpmPR(gO2,gI2,gI1))*CpbarUChiCha2conjHpmPR(gO1,
      gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUChiChiAhPR(gO2,gI1,gI2))*CpbarUChiChiAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)
      ))*Conj(CpbarUChiChihhPR(gO2,gI2,gI1))*CpbarUChiChihhPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1))
      )*Conj(CpbarUChiFdconjSdPR(gO2,gI2,gI1))*CpbarUChiFdconjSdPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1))
      )*Conj(CpbarUChiFeconjSePR(gO2,gI2,gI1))*CpbarUChiFeconjSePR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1))
      )*Conj(CpbarUChiFuconjSuPR(gO2,gI2,gI1))*CpbarUChiFuconjSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MVWm))*Conj(
      CpbarUChiCha1VWmPL(gO2,gI2))*CpbarUChiCha1VWmPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MVWm))*Conj(
      CpbarUChiCha2conjVWmPL(gO2,gI2))*CpbarUChiCha2conjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(
      CpbarUChiChiVZPL(gO2,gI2))*CpbarUChiChiVZPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_Chi_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++)
         self_energy(i, k) = self_energy_Chi_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Chi_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSRdp))*Conj(
      CpbarCha1barUChiSRdpPL(gI1,gO2))*CpbarCha1barUChiSRdpPL(gI1,gO1));
   result += -0.5*SUM(gI1,0,1,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSRum))*Conj(
      CpbarCha2barUChiSRumPL(gI1,gO2))*CpbarCha2barUChiSRumPL(gI1,gO1));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSv(gI2))
      )*Conj(CpbarUChibarFvSvPL(gO2,gI1,gI2))*CpbarUChibarFvSvPL(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2))
      )*Conj(CpbarUChibarFdSdPL(gO2,gI1,gI2))*CpbarUChibarFdSdPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSe(gI2))
      )*Conj(CpbarUChibarFeSePL(gO2,gI1,gI2))*CpbarUChibarFeSePL(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2))
      )*Conj(CpbarUChibarFuSuPL(gO2,gI1,gI2))*CpbarUChibarFuSuPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MRh(gI2)
      ))*Conj(CpbarChibarUChiRhPL(gI1,gO2,gI2))*CpbarChibarUChiRhPL(gI1,gO1,gI2)))
      ;
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MHpm(
      gI1)))*Conj(CpbarUChiCha1HpmPL(gO2,gI2,gI1))*CpbarUChiCha1HpmPL(gO1,gI2,gI1)
      ));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MHpm(
      gI1)))*Conj(CpbarUChiCha2conjHpmPL(gO2,gI2,gI1))*CpbarUChiCha2conjHpmPL(gO1,
      gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MAh(gI2)
      ))*Conj(CpbarUChiChiAhPL(gO2,gI1,gI2))*CpbarUChiChiAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(Mhh(gI1)
      ))*Conj(CpbarUChiChihhPL(gO2,gI2,gI1))*CpbarUChiChihhPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1))
      )*Conj(CpbarUChiFdconjSdPL(gO2,gI2,gI1))*CpbarUChiFdconjSdPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MSe(gI1))
      )*Conj(CpbarUChiFeconjSePL(gO2,gI2,gI1))*CpbarUChiFeconjSePL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1))
      )*Conj(CpbarUChiFuconjSuPL(gO2,gI2,gI1))*CpbarUChiFuconjSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MVWm))*Conj(
      CpbarUChiCha1VWmPR(gO2,gI2))*CpbarUChiCha1VWmPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MVWm))*Conj(
      CpbarUChiCha2conjVWmPR(gO2,gI2))*CpbarUChiCha2conjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVZ))*Conj(
      CpbarUChiChiVZPR(gO2,gI2))*CpbarUChiChiVZPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,4,4> CLASSNAME::self_energy_Chi_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,4,4> self_energy;

   for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++)
         self_energy(i, k) = self_energy_Chi_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha1_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(
      MRh(gI2)))*Conj(CpbarUCha1barCha2RhPL(gO2,gI1,gI2))*CpbarUCha1barCha2RhPR(
      gO1,gI1,gI2)));
   result += SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MAh(gI2)))*Conj(CpbarUCha1Cha1AhPL(gO2,gI1,gI2))*CpbarUCha1Cha1AhPR(gO1,gI1,
      gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MSv(
      gI2)))*Conj(CpbarUCha1barFeSvPL(gO2,gI1,gI2))*CpbarUCha1barFeSvPR(gO1,gI1,
      gI2)));
   result += 3*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MSu
      (gI2)))*Conj(CpbarUCha1barFdSuPL(gO2,gI1,gI2))*CpbarUCha1barFdSuPR(gO1,gI1,
      gI2)));
   result += SUM(gI1,0,3,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSRdp))*Conj(
      CpbarUCha1barChiSRdpPL(gO2,gI1))*CpbarUCha1barChiSRdpPR(gO1,gI1)*MChi(gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUCha1Cha1hhPL(gO2,gI2,gI1))*CpbarUCha1Cha1hhPR(gO1,gI2,gI1)*MCha1(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUCha1ChiconjHpmPL(gO2,gI2,gI1))*CpbarUCha1ChiconjHpmPR(gO1,gI2,gI1
      )*MChi(gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUCha1FuconjSdPL(gO2,gI2,gI1))*CpbarUCha1FuconjSdPR(gO1,gI2,gI1)*
      MFu(gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUCha1FvconjSePL(gO2,gI2,gI1))*CpbarUCha1FvconjSePR(gO1,gI2,gI1)*
      MFv(gI2)));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),0)*Conj(
      CpbarUCha1Cha1VPPR(gO2,gI2))*CpbarUCha1Cha1VPPL(gO1,gI2)*MCha1(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MVZ))*Conj(
      CpbarUCha1Cha1VZPR(gO2,gI2))*CpbarUCha1Cha1VZPL(gO1,gI2)*MCha1(gI2));
   result += -4*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUCha1ChiconjVWmPR(gO2,gI2))*CpbarUCha1ChiconjVWmPL(gO1,gI2)*MChi(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha1_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha1_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha1_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MRh(gI2
      )))*Conj(CpbarUCha1barCha2RhPR(gO2,gI1,gI2))*CpbarUCha1barCha2RhPR(gO1,gI1,
      gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MAh(gI2
      )))*Conj(CpbarUCha1Cha1AhPR(gO2,gI1,gI2))*CpbarUCha1Cha1AhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSv(gI2))
      )*Conj(CpbarUCha1barFeSvPR(gO2,gI1,gI2))*CpbarUCha1barFeSvPR(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSu(gI2))
      )*Conj(CpbarUCha1barFdSuPR(gO2,gI1,gI2))*CpbarUCha1barFdSuPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSRdp))*Conj(
      CpbarUCha1barChiSRdpPR(gO2,gI1))*CpbarUCha1barChiSRdpPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(Mhh(gI1
      )))*Conj(CpbarUCha1Cha1hhPR(gO2,gI2,gI1))*CpbarUCha1Cha1hhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1
      )))*Conj(CpbarUCha1ChiconjHpmPR(gO2,gI2,gI1))*CpbarUCha1ChiconjHpmPR(gO1,gI2
      ,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSd(gI1))
      )*Conj(CpbarUCha1FuconjSdPR(gO2,gI2,gI1))*CpbarUCha1FuconjSdPR(gO1,gI2,gI1))
      );
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSe(gI1))
      )*Conj(CpbarUCha1FvconjSePR(gO2,gI2,gI1))*CpbarUCha1FvconjSePR(gO1,gI2,gI1))
      );
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),0)*Conj(CpbarUCha1Cha1VPPL(
      gO2,gI2))*CpbarUCha1Cha1VPPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MVZ))*Conj(
      CpbarUCha1Cha1VZPL(gO2,gI2))*CpbarUCha1Cha1VZPL(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUCha1ChiconjVWmPL(gO2,gI2))*CpbarUCha1ChiconjVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha1_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha1_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha1_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MRh(gI2
      )))*Conj(CpbarUCha1barCha2RhPL(gO2,gI1,gI2))*CpbarUCha1barCha2RhPL(gO1,gI1,
      gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MAh(gI2
      )))*Conj(CpbarUCha1Cha1AhPL(gO2,gI1,gI2))*CpbarUCha1Cha1AhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MSv(gI2))
      )*Conj(CpbarUCha1barFeSvPL(gO2,gI1,gI2))*CpbarUCha1barFeSvPL(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSu(gI2))
      )*Conj(CpbarUCha1barFdSuPL(gO2,gI1,gI2))*CpbarUCha1barFdSuPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSRdp))*Conj(
      CpbarUCha1barChiSRdpPL(gO2,gI1))*CpbarUCha1barChiSRdpPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(Mhh(gI1
      )))*Conj(CpbarUCha1Cha1hhPL(gO2,gI2,gI1))*CpbarUCha1Cha1hhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1
      )))*Conj(CpbarUCha1ChiconjHpmPL(gO2,gI2,gI1))*CpbarUCha1ChiconjHpmPL(gO1,gI2
      ,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSd(gI1))
      )*Conj(CpbarUCha1FuconjSdPL(gO2,gI2,gI1))*CpbarUCha1FuconjSdPL(gO1,gI2,gI1))
      );
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MSe(gI1))
      )*Conj(CpbarUCha1FvconjSePL(gO2,gI2,gI1))*CpbarUCha1FvconjSePL(gO1,gI2,gI1))
      );
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),0)*Conj(CpbarUCha1Cha1VPPR(
      gO2,gI2))*CpbarUCha1Cha1VPPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MVZ))*Conj(
      CpbarUCha1Cha1VZPR(gO2,gI2))*CpbarUCha1Cha1VZPR(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUCha1ChiconjVWmPR(gO2,gI2))*CpbarUCha1ChiconjVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha1_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha1_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha2_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MRh(gI2)))*Conj(CpbarCha1barUCha2RhPL(gI1,gO2,gI2))*CpbarCha1barUCha2RhPR(
      gI1,gO1,gI2)));
   result += SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(
      MAh(gI2)))*Conj(CpbarUCha2Cha2AhPL(gO2,gI1,gI2))*CpbarUCha2Cha2AhPR(gO1,gI1,
      gI2)));
   result += 3*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd
      (gI2)))*Conj(CpbarUCha2barFuSdPL(gO2,gI1,gI2))*CpbarUCha2barFuSdPR(gO1,gI1,
      gI2)));
   result += SUM(gI1,0,2,MFv(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(
      gI2)))*Conj(CpbarUCha2barFvSePL(gO2,gI1,gI2))*CpbarUCha2barFvSePR(gO1,gI1,
      gI2)));
   result += SUM(gI1,0,3,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSRum))*Conj(
      CpbarUCha2barChiSRumPL(gO2,gI1))*CpbarUCha2barChiSRumPR(gO1,gI1)*MChi(gI1));
   result += SUM(gI1,0,3,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUCha2Cha2hhPL(gO2,gI2,gI1))*CpbarUCha2Cha2hhPR(gO1,gI2,gI1)*MCha2(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUCha2ChiHpmPL(gO2,gI2,gI1))*CpbarUCha2ChiHpmPR(gO1,gI2,gI1)*MChi(
      gI2)));
   result += 3*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUCha2FdconjSuPL(gO2,gI2,gI1))*CpbarUCha2FdconjSuPR(gO1,gI2,gI1)*
      MFd(gI2)));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),0)*Conj(
      CpbarUCha2Cha2VPPR(gO2,gI2))*CpbarUCha2Cha2VPPL(gO1,gI2)*MCha2(gI2));
   result += -4*SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),Sqr(MVZ))*Conj(
      CpbarUCha2Cha2VZPR(gO2,gI2))*CpbarUCha2Cha2VZPL(gO1,gI2)*MCha2(gI2));
   result += -4*SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUCha2ChiVWmPR(gO2,gI2))*CpbarUCha2ChiVWmPL(gO1,gI2)*MChi(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha2_1loop_1(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha2_1loop_1(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha2_1loop_PR(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MRh(gI2
      )))*Conj(CpbarCha1barUCha2RhPR(gI1,gO2,gI2))*CpbarCha1barUCha2RhPR(gI1,gO1,
      gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MAh(gI2
      )))*Conj(CpbarUCha2Cha2AhPR(gO2,gI1,gI2))*CpbarUCha2Cha2AhPR(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2))
      )*Conj(CpbarUCha2barFuSdPR(gO2,gI1,gI2))*CpbarUCha2barFuSdPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2))
      )*Conj(CpbarUCha2barFvSePR(gO2,gI1,gI2))*CpbarUCha2barFvSePR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSRum))*Conj(
      CpbarUCha2barChiSRumPR(gO2,gI1))*CpbarUCha2barChiSRumPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(Mhh(gI1
      )))*Conj(CpbarUCha2Cha2hhPR(gO2,gI2,gI1))*CpbarUCha2Cha2hhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1
      )))*Conj(CpbarUCha2ChiHpmPR(gO2,gI2,gI1))*CpbarUCha2ChiHpmPR(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1))
      )*Conj(CpbarUCha2FdconjSuPR(gO2,gI2,gI1))*CpbarUCha2FdconjSuPR(gO1,gI2,gI1))
      );
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),0)*Conj(CpbarUCha2Cha2VPPL(
      gO2,gI2))*CpbarUCha2Cha2VPPL(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MVZ))*Conj(
      CpbarUCha2Cha2VZPL(gO2,gI2))*CpbarUCha2Cha2VZPL(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUCha2ChiVWmPL(gO2,gI2))*CpbarUCha2ChiVWmPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha2_1loop_PR(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha2_1loop_PR(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Cha2_1loop_PL(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MRh(gI2
      )))*Conj(CpbarCha1barUCha2RhPL(gI1,gO2,gI2))*CpbarCha1barUCha2RhPL(gI1,gO1,
      gI2)));
   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,3,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MAh(gI2
      )))*Conj(CpbarUCha2Cha2AhPL(gO2,gI1,gI2))*CpbarUCha2Cha2AhPL(gO1,gI1,gI2)));
   result += -1.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSd(gI2))
      )*Conj(CpbarUCha2barFuSdPL(gO2,gI1,gI2))*CpbarUCha2barFuSdPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,5,B1(Sqr(p),Sqr(MFv(gI1)),Sqr(MSe(gI2))
      )*Conj(CpbarUCha2barFvSePL(gO2,gI1,gI2))*CpbarUCha2barFvSePL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSRum))*Conj(
      CpbarUCha2barChiSRumPL(gO2,gI1))*CpbarUCha2barChiSRumPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(Mhh(gI1
      )))*Conj(CpbarUCha2Cha2hhPL(gO2,gI2,gI1))*CpbarUCha2Cha2hhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MHpm(gI1
      )))*Conj(CpbarUCha2ChiHpmPL(gO2,gI2,gI1))*CpbarUCha2ChiHpmPL(gO1,gI2,gI1)));
   result += -1.5*SUM(gI1,0,5,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSu(gI1))
      )*Conj(CpbarUCha2FdconjSuPL(gO2,gI2,gI1))*CpbarUCha2FdconjSuPL(gO1,gI2,gI1))
      );
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),0)*Conj(CpbarUCha2Cha2VPPR(
      gO2,gI2))*CpbarUCha2Cha2VPPR(gO1,gI2));
   result += -SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MVZ))*Conj(
      CpbarUCha2Cha2VZPR(gO2,gI2))*CpbarUCha2Cha2VZPR(gO1,gI2));
   result += -SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MVWm))*Conj(
      CpbarUCha2ChiVWmPR(gO2,gI2))*CpbarUCha2ChiVWmPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,2,2> CLASSNAME::self_energy_Cha2_1loop_PL(double p) const
{
   Eigen::Matrix<std::complex<double>,2,2> self_energy;

   for (int i = 0; i < 2; i++)
      for (int k = 0; k < 2; k++)
         self_energy(i, k) = self_energy_Cha2_1loop_PL(p, i, k);

   return self_energy;
}

std::complex<double> CLASSNAME::self_energy_Fe_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MSv(gI2)))*Conj(CpbarCha1barUFeSvPL(gI1,gO2,gI2))*CpbarCha1barUFeSvPR(gI1,
      gO1,gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFeFvHpmPL(gO2,gI2,gI1))*CpbarUFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSe
      (gI2)))*Conj(CpbarChibarUFeSePL(gI1,gO2,gI2))*CpbarChibarUFeSePR(gI1,gO1,gI2
      )));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarUFeChiSePL(gO2,gI2,gI1))*CpbarUFeChiSePR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPR(gO2,
      gI2))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSv(gI2
      )))*Conj(CpbarCha1barUFeSvPR(gI1,gO2,gI2))*CpbarCha1barUFeSvPR(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFeFeAhPR(gO2,gI1,gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFeFehhPR(gO2,gI2,gI1))*CpbarUFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFeFvHpmPR(gO2,gI2,gI1))*CpbarUFeFvHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarChibarUFeSePR(gI1,gO2,gI2))*CpbarChibarUFeSePR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)
      ))*Conj(CpbarUFeChiSePR(gO2,gI2,gI1))*CpbarUFeChiSePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPL(gO2,
      gI2))*CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPL
      (gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(
      CpbarUFeFvVWmPL(gO2,gI2))*CpbarUFeFvVWmPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSv(gI2
      )))*Conj(CpbarCha1barUFeSvPL(gI1,gO2,gI2))*CpbarCha1barUFeSvPL(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFeFeAhPL(gO2,gI1,gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFeFehhPL(gO2,gI2,gI1))*CpbarUFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFeFvHpmPL(gO2,gI2,gI1))*CpbarUFeFvHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarChibarUFeSePL(gI1,gO2,gI2))*CpbarChibarUFeSePL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)
      ))*Conj(CpbarUFeChiSePL(gO2,gI2,gI1))*CpbarUFeChiSePL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),0)*Conj(CpbarUFeFeVPPR(gO2,
      gI2))*CpbarUFeFeVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarUFeFeVZPR
      (gO2,gI2))*CpbarUFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(
      CpbarUFeFvVWmPR(gO2,gI2))*CpbarUFeFvVWmPR(gO1,gI2));

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

   result += SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MSu(gI2)))*Conj(CpbarCha1barUFdSuPL(gI1,gO2,gI2))*CpbarCha1barUFdSuPR(gI1,
      gO1,gI2)));
   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFdFuHpmPL(gO2,gI2,gI1))*CpbarUFdFuHpmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSd
      (gI2)))*Conj(CpbarChibarUFdSdPL(gI1,gO2,gI2))*CpbarChibarUFdSdPR(gI1,gO1,gI2
      )));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(
      gI1)))*Conj(CpbarUFdGluSdPL(gO2,gI1))*CpbarUFdGluSdPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFdCha2SuPL(gO2,gI2,gI1))*CpbarUFdCha2SuPR(gO1,gI2,gI1)*MCha2(gI2)
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFdChiSdPL(gO2,gI2,gI1))*CpbarUFdChiSdPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPR(gO2,
      gI2))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarUFdFuVWmPR(gO2,gI2))*CpbarUFdFuVWmPL(gO1,gI2)*MFu(gI2));
   result += 1.3333333333333333*MGlu*SUM(gI2,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(
      gI2)))*Conj(CpbarGlubarUFdSdPL(gO2,gI2))*CpbarGlubarUFdSdPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSu(gI2
      )))*Conj(CpbarCha1barUFdSuPR(gI1,gO2,gI2))*CpbarCha1barUFdSuPR(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFdFdAhPR(gO2,gI1,gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFdFdhhPR(gO2,gI2,gI1))*CpbarUFdFdhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFdFuHpmPR(gO2,gI2,gI1))*CpbarUFdFuHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSd(gI2)
      ))*Conj(CpbarChibarUFdSdPR(gI1,gO2,gI2))*CpbarChibarUFdSdPR(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))
      *Conj(CpbarUFdGluSdPR(gO2,gI1))*CpbarUFdGluSdPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MSu(gI1
      )))*Conj(CpbarUFdCha2SuPR(gO2,gI2,gI1))*CpbarUFdCha2SuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)
      ))*Conj(CpbarUFdChiSdPR(gO2,gI2,gI1))*CpbarUFdChiSdPR(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPL(gO2,
      gI2))*CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPL
      (gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarUFdFuVWmPL(gO2,gI2))*CpbarUFdFuVWmPL(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI2)))
      *Conj(CpbarGlubarUFdSdPR(gO2,gI2))*CpbarGlubarUFdSdPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSu(gI2
      )))*Conj(CpbarCha1barUFdSuPL(gI1,gO2,gI2))*CpbarCha1barUFdSuPL(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFdFdAhPL(gO2,gI1,gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFdFdhhPL(gO2,gI2,gI1))*CpbarUFdFdhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFdFuHpmPL(gO2,gI2,gI1))*CpbarUFdFuHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSd(gI2)
      ))*Conj(CpbarChibarUFdSdPL(gI1,gO2,gI2))*CpbarChibarUFdSdPL(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))
      *Conj(CpbarUFdGluSdPL(gO2,gI1))*CpbarUFdGluSdPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MSu(gI1
      )))*Conj(CpbarUFdCha2SuPL(gO2,gI2,gI1))*CpbarUFdCha2SuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)
      ))*Conj(CpbarUFdChiSdPL(gO2,gI2,gI1))*CpbarUFdChiSdPL(gO1,gI2,gI1)));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(
      CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),0)*Conj(CpbarUFdFdVPPR(gO2,
      gI2))*CpbarUFdFdVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarUFdFdVZPR
      (gO2,gI2))*CpbarUFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarUFdFuVWmPR(gO2,gI2))*CpbarUFdFuVWmPR(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI2)))
      *Conj(CpbarGlubarUFdSdPL(gO2,gI2))*CpbarGlubarUFdSdPL(gO1,gI2));

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

   result += SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(
      MSd(gI2)))*Conj(CpbarCha2barUFuSdPL(gI1,gO2,gI2))*CpbarCha2barUFuSdPR(gI1,
      gO1,gI2)));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)*MFd(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu
      (gI2)))*Conj(CpbarChibarUFuSuPL(gI1,gO2,gI2))*CpbarChibarUFuSuPR(gI1,gO1,gI2
      )));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(
      gI1)))*Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFuCha1SdPL(gO2,gI2,gI1))*CpbarUFuCha1SdPR(gO1,gI2,gI1)*MCha1(gI2)
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,
      gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));
   result += 1.3333333333333333*MGlu*SUM(gI2,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(
      gI2)))*Conj(CpbarGlubarUFuSuPL(gO2,gI2))*CpbarGlubarUFuSuPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSd(gI2
      )))*Conj(CpbarCha2barUFuSdPR(gI1,gO2,gI2))*CpbarCha2barUFuSdPR(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFuFdconjHpmPR(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)))
      ;
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu(gI2)
      ))*Conj(CpbarChibarUFuSuPR(gI1,gO2,gI2))*CpbarChibarUFuSuPR(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))
      *Conj(CpbarUFuGluSuPR(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1
      )))*Conj(CpbarUFuCha1SdPR(gO2,gI2,gI1))*CpbarUFuCha1SdPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)
      ))*Conj(CpbarUFuChiSuPR(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPL(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPL(gO2,
      gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPL
      (gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI2)))
      *Conj(CpbarGlubarUFuSuPR(gO2,gI2))*CpbarGlubarUFuSuPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSd(gI2
      )))*Conj(CpbarCha2barUFuSdPL(gI1,gO2,gI2))*CpbarCha2barUFuSdPL(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPL(gO1,gI2,gI1)))
      ;
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu(gI2)
      ))*Conj(CpbarChibarUFuSuPL(gI1,gO2,gI2))*CpbarChibarUFuSuPL(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))
      *Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1
      )))*Conj(CpbarUFuCha1SdPL(gO2,gI2,gI1))*CpbarUFuCha1SdPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)
      ))*Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(
      CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,
      gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR
      (gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI2)))
      *Conj(CpbarGlubarUFuSuPL(gO2,gI2))*CpbarGlubarUFuSuPL(gO1,gI2));

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

   result += 3*MGlu*B0(Sqr(p),Sqr(MGlu),Sqr(MphiO))*Conj(CpbarGluGluphiOPL())*
      CpbarGluGluphiOPR();
   result += 3*MGlu*B0(Sqr(p),Sqr(MGlu),Sqr(MsigmaO))*Conj(CpbarGluGlusigmaOPL(
      ))*CpbarGluGlusigmaOPR();
   result += -12*MGlu*B0(Sqr(p),Sqr(MGlu),0)*Conj(CpbarGluGluVGPR())*
      CpbarGluGluVGPL();
   result += 0.5*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(
      MSd(gI2)))*Conj(CpbarGlubarFdSdPL(gI1,gI2))*CpbarGlubarFdSdPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(
      MSu(gI2)))*Conj(CpbarGlubarFuSuPL(gI1,gI2))*CpbarGlubarFuSuPR(gI1,gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))
      *Conj(CpbarGluFdconjSdPL(gI2,gI1))*CpbarGluFdconjSdPR(gI2,gI1)*MFd(gI2)));
   result += 0.5*SUM(gI1,0,5,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))
      *Conj(CpbarGluFuconjSuPL(gI2,gI1))*CpbarGluFuconjSuPR(gI2,gI1)*MFu(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1loop_PR(double p ) const
{
   std::complex<double> result;

   result += -1.5*AbsSqr(CpbarGluGluphiOPR())*B1(Sqr(p),Sqr(MGlu),Sqr(MphiO));
   result += -1.5*AbsSqr(CpbarGluGlusigmaOPR())*B1(Sqr(p),Sqr(MGlu),Sqr(MsigmaO
      ));
   result += -3*AbsSqr(CpbarGluGluVGPL())*B1(Sqr(p),Sqr(MGlu),0);
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarGlubarFdSdPR(gI1,gI2))*
      B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))));
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarGlubarFuSuPR(gI1,gI2))*
      B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpbarGluFdconjSdPR(gI2,gI1))*
      B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpbarGluFuconjSuPR(gI2,gI1))*
      B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1loop_PL(double p ) const
{
   std::complex<double> result;

   result += -1.5*AbsSqr(CpbarGluGluphiOPL())*B1(Sqr(p),Sqr(MGlu),Sqr(MphiO));
   result += -1.5*AbsSqr(CpbarGluGlusigmaOPL())*B1(Sqr(p),Sqr(MGlu),Sqr(MsigmaO
      ));
   result += -3*AbsSqr(CpbarGluGluVGPR())*B1(Sqr(p),Sqr(MGlu),0);
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarGlubarFdSdPL(gI1,gI2))*
      B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MSd(gI2)))));
   result += -0.25*SUM(gI1,0,2,SUM(gI2,0,5,AbsSqr(CpbarGlubarFuSuPL(gI1,gI2))*
      B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MSu(gI2)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpbarGluFdconjSdPL(gI2,gI1))*
      B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MSd(gI1)))));
   result += -0.25*SUM(gI1,0,5,SUM(gI2,0,2,AbsSqr(CpbarGluFuconjSuPL(gI2,gI1))*
      B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MSu(gI1)))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1loop_1(double p , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(
      MSe(gI2)))*Conj(CpbarCha2barFvSePL(gI1,gO2,gI2))*CpbarCha2barFvSePR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSv
      (gI2)))*Conj(CpbarChibarFvSvPL(gI1,gO2,gI2))*CpbarChibarFvSvPR(gI1,gO1,gI2))
      );
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFvFeconjHpmPL(gO2,gI2,gI1))*CpbarFvFeconjHpmPR(gO1,gI2,gI1)*MFe(
      gI2)));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarFvCha1SePL(gO2,gI2,gI1))*CpbarFvCha1SePR(gO1,gI2,gI1)*MCha1(gI2)))
      ;
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(
      CpbarFvFvVZPR(gO2,gI2))*CpbarFvFvVZPL(gO1,gI2)*MFv(gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSe(gI2
      )))*Conj(CpbarCha2barFvSePR(gI1,gO2,gI2))*CpbarCha2barFvSePR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSv(gI2)
      ))*Conj(CpbarChibarFvSvPR(gI1,gO2,gI2))*CpbarChibarFvSvPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFvFeconjHpmPR(gO2,gI2,gI1))*CpbarFvFeconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSe(gI1
      )))*Conj(CpbarFvCha1SePR(gO2,gI2,gI1))*CpbarFvCha1SePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPL(gO2,gI2))*CpbarFvFeconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPL(
      gO2,gI2))*CpbarFvFvVZPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSe(gI2
      )))*Conj(CpbarCha2barFvSePL(gI1,gO2,gI2))*CpbarCha2barFvSePL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSv(gI2)
      ))*Conj(CpbarChibarFvSvPL(gI1,gO2,gI2))*CpbarChibarFvSvPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFvFeconjHpmPL(gO2,gI2,gI1))*CpbarFvFeconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSe(gI1
      )))*Conj(CpbarFvCha1SePL(gO2,gI2,gI1))*CpbarFvCha1SePL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVWm))*Conj(
      CpbarFvFeconjVWmPR(gO2,gI2))*CpbarFvFeconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVZ))*Conj(CpbarFvFvVZPR(
      gO2,gI2))*CpbarFvFvVZPR(gO1,gI2));

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

   result += SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,2,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MSv(gI2)))*Conj(CpbarCha1barFeSvPL(gI1,gO2,gI2))*CpbarCha1barFeSvPR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)*MFe(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)*MFv(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSe
      (gI2)))*Conj(CpbarChibarFeSePL(gI1,gO2,gI2))*CpbarChibarFeSePR(gI1,gO1,gI2))
      );
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)))*
      Conj(CpbarFeChiSePL(gO2,gI2,gI1))*CpbarFeChiSePR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(
      CpbarFeFvVWmPR(gO2,gI2))*CpbarFeFvVWmPL(gO1,gI2)*MFv(gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSv(gI2
      )))*Conj(CpbarCha1barFeSvPR(gI1,gO2,gI2))*CpbarCha1barFeSvPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarFeFeAhPR(gO2,gI1,gI2))*CpbarFeFeAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarFeFehhPR(gO2,gI2,gI1))*CpbarFeFehhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFeFvHpmPR(gO2,gI2,gI1))*CpbarFeFvHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarChibarFeSePR(gI1,gO2,gI2))*CpbarChibarFeSePR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)
      ))*Conj(CpbarFeChiSePR(gO2,gI2,gI1))*CpbarFeChiSePR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPL(
      gO2,gI2))*CpbarFeFeVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(
      CpbarFeFvVWmPL(gO2,gI2))*CpbarFeFvVWmPL(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,2,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSv(gI2
      )))*Conj(CpbarCha1barFeSvPL(gI1,gO2,gI2))*CpbarCha1barFeSvPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFe(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarFeFeAhPL(gO2,gI1,gI2))*CpbarFeFeAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarFeFehhPL(gO2,gI2,gI1))*CpbarFeFehhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFeFvHpmPL(gO2,gI2,gI1))*CpbarFeFvHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSe(gI2)
      ))*Conj(CpbarChibarFeSePL(gI1,gO2,gI2))*CpbarChibarFeSePL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSe(gI1)
      ))*Conj(CpbarFeChiSePL(gO2,gI2,gI1))*CpbarFeChiSePL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFe(gI2)),Sqr(MVZ))*Conj(CpbarFeFeVZPR(
      gO2,gI2))*CpbarFeFeVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFv(gI2)),Sqr(MVWm))*Conj(
      CpbarFeFvVWmPR(gO2,gI2))*CpbarFeFvVWmPR(gO1,gI2));

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

   result += SUM(gI1,0,1,MCha1(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha1(gI1)),Sqr(
      MSu(gI2)))*Conj(CpbarCha1barFdSuPL(gI1,gO2,gI2))*CpbarCha1barFdSuPR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)*MFd(gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFdFuHpmPL(gO2,gI2,gI1))*CpbarFdFuHpmPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSd
      (gI2)))*Conj(CpbarChibarFdSdPL(gI1,gO2,gI2))*CpbarChibarFdSdPR(gI1,gO1,gI2))
      );
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(
      gI1)))*Conj(CpbarFdGluSdPL(gO2,gI1))*CpbarFdGluSdPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha2(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFdCha2SuPL(gO2,gI2,gI1))*CpbarFdCha2SuPR(gO1,gI2,gI1)*MCha2(gI2)))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarFdChiSdPL(gO2,gI2,gI1))*CpbarFdChiSdPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarFdFuVWmPR(gO2,gI2))*CpbarFdFuVWmPL(gO1,gI2)*MFu(gI2));
   result += 1.3333333333333333*MGlu*SUM(gI2,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSd(
      gI2)))*Conj(CpbarGlubarFdSdPL(gO2,gI2))*CpbarGlubarFdSdPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSu(gI2
      )))*Conj(CpbarCha1barFdSuPR(gI1,gO2,gI2))*CpbarCha1barFdSuPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarFdFdAhPR(gO2,gI1,gI2))*CpbarFdFdAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarFdFdhhPR(gO2,gI2,gI1))*CpbarFdFdhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFdFuHpmPR(gO2,gI2,gI1))*CpbarFdFuHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSd(gI2)
      ))*Conj(CpbarChibarFdSdPR(gI1,gO2,gI2))*CpbarChibarFdSdPR(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))
      *Conj(CpbarFdGluSdPR(gO2,gI1))*CpbarFdGluSdPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MSu(gI1
      )))*Conj(CpbarFdCha2SuPR(gO2,gI2,gI1))*CpbarFdCha2SuPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)
      ))*Conj(CpbarFdChiSdPR(gO2,gI2,gI1))*CpbarFdChiSdPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPL(
      gO2,gI2))*CpbarFdFdVZPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarFdFuVWmPL(gO2,gI2))*CpbarFdFuVWmPL(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI2)))
      *Conj(CpbarGlubarFdSdPR(gO2,gI2))*CpbarGlubarFdSdPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha1(gI1)),Sqr(MSu(gI2
      )))*Conj(CpbarCha1barFdSuPL(gI1,gO2,gI2))*CpbarCha1barFdSuPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFd(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarFdFdAhPL(gO2,gI1,gI2))*CpbarFdFdAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarFdFdhhPL(gO2,gI2,gI1))*CpbarFdFdhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFdFuHpmPL(gO2,gI2,gI1))*CpbarFdFuHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSd(gI2)
      ))*Conj(CpbarChibarFdSdPL(gI1,gO2,gI2))*CpbarChibarFdSdPL(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI1)))
      *Conj(CpbarFdGluSdPL(gO2,gI1))*CpbarFdGluSdPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha2(gI2)),Sqr(MSu(gI1
      )))*Conj(CpbarFdCha2SuPL(gO2,gI2,gI1))*CpbarFdCha2SuPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSd(gI1)
      ))*Conj(CpbarFdChiSdPL(gO2,gI2,gI1))*CpbarFdChiSdPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVZ))*Conj(CpbarFdFdVZPR(
      gO2,gI2))*CpbarFdFdVZPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVWm))*Conj(
      CpbarFdFuVWmPR(gO2,gI2))*CpbarFdFuVWmPR(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSd(gI2)))
      *Conj(CpbarGlubarFdSdPL(gO2,gI2))*CpbarGlubarFdSdPL(gO1,gI2));

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

   result += SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(
      MSd(gI2)))*Conj(CpbarCha2barFuSdPL(gI1,gO2,gI2))*CpbarCha2barFuSdPR(gI1,gO1,
      gI2)));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarFuFdconjHpmPL(gO2,gI2,gI1))*CpbarFuFdconjHpmPR(gO1,gI2,gI1)*MFd(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu
      (gI2)))*Conj(CpbarChibarFuSuPL(gI1,gO2,gI2))*CpbarChibarFuSuPR(gI1,gO1,gI2))
      );
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(
      gI1)))*Conj(CpbarFuGluSuPL(gO2,gI1))*CpbarFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarFuCha1SdPL(gO2,gI2,gI1))*CpbarFuCha1SdPR(gO1,gI2,gI1)*MCha1(gI2)))
      ;
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarFuChiSuPL(gO2,gI2,gI1))*CpbarFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPR(gO2,
      gI2))*CpbarFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(
      CpbarFuFuVZPR(gO2,gI2))*CpbarFuFuVZPL(gO1,gI2)*MFu(gI2));
   result += 1.3333333333333333*MGlu*SUM(gI2,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(
      gI2)))*Conj(CpbarGlubarFuSuPL(gO2,gI2))*CpbarGlubarFuSuPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSd(gI2
      )))*Conj(CpbarCha2barFuSdPR(gI1,gO2,gI2))*CpbarCha2barFuSdPR(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarFuFuAhPR(gO2,gI1,gI2))*CpbarFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFuFdconjHpmPR(gO2,gI2,gI1))*CpbarFuFdconjHpmPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarFuFuhhPR(gO2,gI2,gI1))*CpbarFuFuhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu(gI2)
      ))*Conj(CpbarChibarFuSuPR(gI1,gO2,gI2))*CpbarChibarFuSuPR(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))
      *Conj(CpbarFuGluSuPR(gO2,gI1))*CpbarFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1
      )))*Conj(CpbarFuCha1SdPR(gO2,gI2,gI1))*CpbarFuCha1SdPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)
      ))*Conj(CpbarFuChiSuPR(gO2,gI2,gI1))*CpbarFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPL(gO2,gI2))*CpbarFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPL(gO2,gI2
      ))*CpbarFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPL(
      gO2,gI2))*CpbarFuFuVZPL(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI2)))
      *Conj(CpbarGlubarFuSuPR(gO2,gI2))*CpbarGlubarFuSuPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSd(gI2
      )))*Conj(CpbarCha2barFuSdPL(gI1,gO2,gI2))*CpbarCha2barFuSdPL(gI1,gO1,gI2)));
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarFuFuAhPL(gO2,gI1,gI2))*CpbarFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarFuFdconjHpmPL(gO2,gI2,gI1))*CpbarFuFdconjHpmPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarFuFuhhPL(gO2,gI2,gI1))*CpbarFuFuhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu(gI2)
      ))*Conj(CpbarChibarFuSuPL(gI1,gO2,gI2))*CpbarChibarFuSuPL(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))
      *Conj(CpbarFuGluSuPL(gO2,gI1))*CpbarFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1
      )))*Conj(CpbarFuCha1SdPL(gO2,gI2,gI1))*CpbarFuCha1SdPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)
      ))*Conj(CpbarFuChiSuPL(gO2,gI2,gI1))*CpbarFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarFuFdconjVWmPR(gO2,gI2))*CpbarFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarFuFuVPPR(gO2,gI2
      ))*CpbarFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarFuFuVZPR(
      gO2,gI2))*CpbarFuFuVZPR(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI2)))
      *Conj(CpbarGlubarFuSuPL(gO2,gI2))*CpbarGlubarFuSuPL(gO1,gI2));

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

   result += SUM(gI1,0,1,MCha2(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MCha2(gI1)),Sqr(
      MSd(gI2)))*Conj(CpbarCha2barUFuSdPL(gI1,gO2,gI2))*CpbarCha2barUFuSdPR(gI1,
      gO1,gI2)));
   result += SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,3,B0(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(
      gI2)))*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)))*
      Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)*MFd(
      gI2)));
   result += SUM(gI1,0,3,SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1)))*
      Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)*MFu(gI2)));
   result += SUM(gI1,0,3,MChi(gI1)*SUM(gI2,0,5,B0(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu
      (gI2)))*Conj(CpbarChibarUFuSuPL(gI1,gO2,gI2))*CpbarChibarUFuSuPR(gI1,gO1,gI2
      )));
   result += 1.3333333333333333*MGlu*SUM(gI1,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(
      gI1)))*Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += SUM(gI1,0,5,SUM(gI2,0,1,B0(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1)))*
      Conj(CpbarUFuCha1SdPL(gO2,gI2,gI1))*CpbarUFuCha1SdPR(gO1,gI2,gI1)*MCha1(gI2)
      ));
   result += SUM(gI1,0,5,SUM(gI2,0,3,B0(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)))*
      Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)*MChi(gI2)));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,
      gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,B0(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));
   result += 1.3333333333333333*MGlu*SUM(gI2,0,5,B0(Sqr(p),Sqr(MGlu),Sqr(MSu(
      gI2)))*Conj(CpbarGlubarUFuSuPL(gO2,gI2))*CpbarGlubarUFuSuPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSd(gI2
      )))*Conj(CpbarCha2barUFuSdPR(gI1,gO2,gI2))*CpbarCha2barUFuSdPR(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFuFuAhPR(gO2,gI1,gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFuFdconjHpmPR(gO2,gI2,gI1))*CpbarUFuFdconjHpmPR(gO1,gI2,gI1)))
      ;
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFuFuhhPR(gO2,gI2,gI1))*CpbarUFuFuhhPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu(gI2)
      ))*Conj(CpbarChibarUFuSuPR(gI1,gO2,gI2))*CpbarChibarUFuSuPR(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))
      *Conj(CpbarUFuGluSuPR(gO2,gI1))*CpbarUFuGluSuPR(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1
      )))*Conj(CpbarUFuCha1SdPR(gO2,gI2,gI1))*CpbarUFuCha1SdPR(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)
      ))*Conj(CpbarUFuChiSuPR(gO2,gI2,gI1))*CpbarUFuChiSuPR(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPL(gO2,gI2))*CpbarUFuFdconjVWmPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPL(gO2,
      gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPL
      (gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI2)))
      *Conj(CpbarGlubarUFuSuPR(gO2,gI2))*CpbarGlubarUFuSuPR(gO1,gI2));

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

   result += -0.5*SUM(gI1,0,1,SUM(gI2,0,5,B1(Sqr(p),Sqr(MCha2(gI1)),Sqr(MSd(gI2
      )))*Conj(CpbarCha2barUFuSdPL(gI1,gO2,gI2))*CpbarCha2barUFuSdPL(gI1,gO1,gI2))
      );
   result += -0.5*SUM(gI1,0,2,SUM(gI2,0,3,B1(Sqr(p),Sqr(MFu(gI1)),Sqr(MAh(gI2))
      )*Conj(CpbarUFuFuAhPL(gO2,gI1,gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MHpm(gI1)
      ))*Conj(CpbarUFuFdconjHpmPL(gO2,gI2,gI1))*CpbarUFuFdconjHpmPL(gO1,gI2,gI1)))
      ;
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(Mhh(gI1))
      )*Conj(CpbarUFuFuhhPL(gO2,gI2,gI1))*CpbarUFuFuhhPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,3,SUM(gI2,0,5,B1(Sqr(p),Sqr(MChi(gI1)),Sqr(MSu(gI2)
      ))*Conj(CpbarChibarUFuSuPL(gI1,gO2,gI2))*CpbarChibarUFuSuPL(gI1,gO1,gI2)));
   result += -0.6666666666666666*SUM(gI1,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI1)))
      *Conj(CpbarUFuGluSuPL(gO2,gI1))*CpbarUFuGluSuPL(gO1,gI1));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,1,B1(Sqr(p),Sqr(MCha1(gI2)),Sqr(MSd(gI1
      )))*Conj(CpbarUFuCha1SdPL(gO2,gI2,gI1))*CpbarUFuCha1SdPL(gO1,gI2,gI1)));
   result += -0.5*SUM(gI1,0,5,SUM(gI2,0,3,B1(Sqr(p),Sqr(MChi(gI2)),Sqr(MSu(gI1)
      ))*Conj(CpbarUFuChiSuPL(gO2,gI2,gI1))*CpbarUFuChiSuPL(gO1,gI2,gI1)));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFd(gI2)),Sqr(MVWm))*Conj(
      CpbarUFuFdconjVWmPR(gO2,gI2))*CpbarUFuFdconjVWmPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),0)*Conj(CpbarUFuFuVPPR(gO2,
      gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,B1(Sqr(p),Sqr(MFu(gI2)),Sqr(MVZ))*Conj(CpbarUFuFuVZPR
      (gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));
   result += -0.6666666666666666*SUM(gI2,0,5,B1(Sqr(p),Sqr(MGlu),Sqr(MSu(gI2)))
      *Conj(CpbarGlubarUFuSuPL(gO2,gI2))*CpbarGlubarUFuSuPL(gO1,gI2));

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
   result += -(A0(Sqr(MSRdp))*CpSRdpUhhconjSRdp(gO1));
   result += -(A0(Sqr(MSRum))*CpSRumUhhconjSRum(gO1));
   result += 4*A0(Sqr(MVWm))*CpUhhconjVWmVWm(gO1);
   result += 2*A0(Sqr(MVZ))*CpUhhVZVZ(gO1);
   result += -SUM(gI1,0,1,A0(Sqr(MRh(gI1)))*CpUhhRhconjRh(gO1,gI1,gI1));
   result += 2*SUM(gI1,0,1,A0(Sqr(MCha1(gI1)))*(CpbarCha1Cha1UhhPL(gI1,gI1,gO1)
      + CpbarCha1Cha1UhhPR(gI1,gI1,gO1))*MCha1(gI1));
   result += 2*SUM(gI1,0,1,A0(Sqr(MCha2(gI1)))*(CpbarCha2Cha2UhhPL(gI1,gI1,gO1)
      + CpbarCha2Cha2UhhPR(gI1,gI1,gO1))*MCha2(gI1));
   result += -SUM(gI1,0,2,A0(Sqr(MSv(gI1)))*CpUhhSvconjSv(gO1,gI1,gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdUhhPL(gI1,gI1,gO1) +
      CpbarFdFdUhhPR(gI1,gI1,gO1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFeUhhPL(gI1,gI1,gO1) +
      CpbarFeFeUhhPR(gI1,gI1,gO1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuUhhPL(gI1,gI1,gO1) +
      CpbarFuFuUhhPR(gI1,gI1,gO1))*MFu(gI1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(MAh(gI1)))*CpAhAhUhh(gI1,gI1,gO1));
   result += -0.5*SUM(gI1,0,3,A0(Sqr(Mhh(gI1)))*CphhhhUhh(gI1,gI1,gO1));
   result += -SUM(gI1,0,3,A0(Sqr(MHpm(gI1)))*CpUhhHpmconjHpm(gO1,gI1,gI1));
   result += 2*SUM(gI1,0,3,A0(Sqr(MChi(gI1)))*(CpbarChiChiUhhPL(gI1,gI1,gO1) +
      CpbarChiChiUhhPR(gI1,gI1,gO1))*MChi(gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpUhhSdconjSd(gO1,gI1,gI1));
   result += -SUM(gI1,0,5,A0(Sqr(MSe(gI1)))*CpUhhSeconjSe(gO1,gI1,gI1));
   result += -3*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpUhhSuconjSu(gO1,gI1,gI1));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_phiO_1loop() const
{
   std::complex<double> result;

   result += 6*MGlu*A0(Sqr(MGlu))*(CpbarGluGluphiOPL() + CpbarGluGluphiOPR());
   result += -0.5*SUM(gI1,0,5,A0(Sqr(MSd(gI1)))*CpphiOSdconjSd(gI1,gI1));
   result += -0.5*SUM(gI1,0,5,A0(Sqr(MSu(gI1)))*CpphiOSuconjSu(gI1,gI1));

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
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MGlu) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void CLASSNAME::calculate_MSRdp_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::SRdp))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_SRdp());
   const double p = MSRdp;
   double self_energy = Re(self_energy_SRdp_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   PHYSICAL(MSRdp) = SignedAbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MSRum_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::SRum))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_SRum());
   const double p = MSRum;
   double self_energy = Re(self_energy_SRum_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   PHYSICAL(MSRum) = SignedAbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MsigmaO_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::sigmaO))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_sigmaO());
   const double p = MsigmaO;
   double self_energy = Re(self_energy_sigmaO_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   PHYSICAL(MsigmaO) = SignedAbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MphiO_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::phiO))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_phiO());
   const double p = MphiO;
   double self_energy = Re(self_energy_phiO_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   PHYSICAL(MphiO) = SignedAbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(MRSSM_info::VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Sd))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sd());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSd(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Sd_1loop(
         p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZD;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZD,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Sd, eigenvalue_error >
            precision * Abs(eigen_values(0)));
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
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Sv))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Sv());

   for (int es = 0; es < 3; ++es) {
      const double p = Abs(MSv(es));
      Eigen::Matrix<double,3,3> self_energy = Re(self_energy_Sv_1loop(
         p));
      const Eigen::Matrix<double,3,3> M_loop(M_tree - self_energy);
      Eigen::Array<double,3,1> eigen_values;
      Eigen::Matrix<double,3,3> mix_ZV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZV,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Sv, eigenvalue_error >
            precision * Abs(eigen_values(0)));
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
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Su))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Su());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSu(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Su_1loop(
         p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZU;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZU,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Su, eigenvalue_error >
            precision * Abs(eigen_values(0)));
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
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Se))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Se());

   for (int es = 0; es < 6; ++es) {
      const double p = Abs(MSe(es));
      Eigen::Matrix<double,6,6> self_energy = Re(self_energy_Se_1loop(
         p));
      const Eigen::Matrix<double,6,6> M_loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZE;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZE,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Se, eigenvalue_error >
            precision * Abs(eigen_values(0)));
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
   if (!force_output && problems.is_running_tachyon(MRSSM_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations()
      ;
   int iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_hh());

      for (int es = 0; es < 4; ++es) {
         const double p = Abs(old_Mhh(es));
         Eigen::Matrix<double,4,4> self_energy = Re(
            self_energy_hh_1loop(p));
         const Eigen::Matrix<double,4,4> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(MRSSM_info::hh,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH);
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
      problems.flag_no_pole_mass_convergence(MRSSM_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(MRSSM_info::hh);
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Ah))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations()
      ;
   int iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   do {
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Ah());

      for (int es = 0; es < 4; ++es) {
         const double p = Abs(old_MAh(es));
         Eigen::Matrix<double,4,4> self_energy = Re(
            self_energy_Ah_1loop(p));
         const Eigen::Matrix<double,4,4> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(MRSSM_info::Ah,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZA);
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
      problems.flag_no_pole_mass_convergence(MRSSM_info::Ah);
   else
      problems.unflag_no_pole_mass_convergence(MRSSM_info::Ah);
}

void CLASSNAME::calculate_MRh_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Rh))
      return;

   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Rh());

   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MRh(es));
      Eigen::Matrix<double,2,2> self_energy = Re(self_energy_Rh_1loop(
         p));
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZHR;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZHR,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Rh, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZHR);
      #endif
         normalize_to_interval(mix_ZHR);

      PHYSICAL(MRh(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZHR) = mix_ZHR;
   }
}

void CLASSNAME::calculate_MHpm_pole()
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::Hpm))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations()
      ;
   int iteration = 0;
   double diff = 0.0;
   decltype(MHpm) old_MHpm(MHpm), new_MHpm(MHpm);

   do {
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Hpm());

      for (int es = 0; es < 4; ++es) {
         const double p = Abs(old_MHpm(es));
         Eigen::Matrix<double,4,4> self_energy = Re(
            self_energy_Hpm_1loop(p));
         const Eigen::Matrix<double,4,4> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZP, eigenvalue_error);
            problems.flag_bad_mass(MRSSM_info::Hpm,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZP);
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
      problems.flag_no_pole_mass_convergence(MRSSM_info::Hpm);
   else
      problems.unflag_no_pole_mass_convergence(MRSSM_info::Hpm);
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Chi());
   for (int es = 0; es < 4; ++es) {
      const double p = Abs(MChi(es));
      const Eigen::Matrix<double,4,4> self_energy_1  = Re(
         self_energy_Chi_1loop_1(p));
      const Eigen::Matrix<double,4,4> self_energy_PL = Re(
         self_energy_Chi_1loop_PL(p));
      const Eigen::Matrix<double,4,4> self_energy_PR = Re(
         self_energy_Chi_1loop_PR(p));
      const Eigen::Matrix<double,4,4> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,4,4> M_loop(M_tree + delta_M);
      Eigen::Array<double,4,1> eigen_values;
      decltype(ZN1) mix_ZN1;
      decltype(ZN2) mix_ZN2;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZN1, mix_ZN2, eigenvalue_error)
         ;
      problems.flag_bad_mass(MRSSM_info::Chi, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_ZN1, mix_ZN2);
   #endif
      if (es == 0) {
         PHYSICAL(ZN1) = mix_ZN1;
         PHYSICAL(ZN2) = mix_ZN2;
      }
      PHYSICAL(MChi(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MCha1_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha1());
   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MCha1(es));
      const Eigen::Matrix<double,2,2> self_energy_1  = Re(
         self_energy_Cha1_1loop_1(p));
      const Eigen::Matrix<double,2,2> self_energy_PL = Re(
         self_energy_Cha1_1loop_PL(p));
      const Eigen::Matrix<double,2,2> self_energy_PR = Re(
         self_energy_Cha1_1loop_PR(p));
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM1) mix_UM1;
      decltype(UP1) mix_UP1;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_UM1, mix_UP1, eigenvalue_error)
         ;
      problems.flag_bad_mass(MRSSM_info::Cha1, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_UM1, mix_UP1);
   #endif
      if (es == 0) {
         PHYSICAL(UM1) = mix_UM1;
         PHYSICAL(UP1) = mix_UP1;
      }
      PHYSICAL(MCha1(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MCha2_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha2());
   for (int es = 0; es < 2; ++es) {
      const double p = Abs(MCha2(es));
      const Eigen::Matrix<double,2,2> self_energy_1  = Re(
         self_energy_Cha2_1loop_1(p));
      const Eigen::Matrix<double,2,2> self_energy_PL = Re(
         self_energy_Cha2_1loop_PL(p));
      const Eigen::Matrix<double,2,2> self_energy_PR = Re(
         self_energy_Cha2_1loop_PR(p));
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM2) mix_UM2;
      decltype(UP2) mix_UP2;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_UM2, mix_UP2, eigenvalue_error)
         ;
      problems.flag_bad_mass(MRSSM_info::Cha2, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_UM2, mix_UP2);
   #endif
      if (es == 0) {
         PHYSICAL(UM2) = mix_UM2;
         PHYSICAL(UP2) = mix_UP2;
      }
      PHYSICAL(MCha2(es)) = Abs(eigen_values(es));
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
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZEL) mix_ZEL;
      decltype(ZER) mix_ZER;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZEL, mix_ZER, eigenvalue_error)
         ;
      problems.flag_bad_mass(MRSSM_info::Fe, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZDL) mix_ZDL;
      decltype(ZDR) mix_ZDR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZDL, mix_ZDR, eigenvalue_error)
         ;
      problems.flag_bad_mass(MRSSM_info::Fd, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
      qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(2))/Sqr(
         currentScale)))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = -0.005191204615668296*Quad(g3) - 0.0032883224409535764*
         Log(Sqr(currentScale)/Sqr(MFu(2)))*Quad(g3) - 0.0008822328500119351*
         Quad(g3)*Sqr(Log(Sqr(currentScale)/Sqr(MFu(2))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = 0;
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
               self_energy_1(i1,i2)  = Re(
                  self_energy_Fu_1loop_1_heavy(p,i1,i2));
               self_energy_PL(i1,i2) = Re(
                  self_energy_Fu_1loop_PL_heavy(p,i1,i2));
               self_energy_PR(i1,i2) = Re(
                  self_energy_Fu_1loop_PR_heavy(p,i1,i2));
            } else {
               self_energy_1(i1,i2)  = Re(
                  self_energy_Fu_1loop_1(p,i1,i2));
               self_energy_PL(i1,i2) = Re(
                  self_energy_Fu_1loop_PL(p,i1,i2));
               self_energy_PR(i1,i2) = Re(
                  self_energy_Fu_1loop_PR(p,i1,i2));
            }
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZUL) mix_ZUL;
      decltype(ZUR) mix_ZUR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_ZUL, mix_ZUR, eigenvalue_error)
         ;
      problems.flag_bad_mass(MRSSM_info::Fu, eigenvalue_error >
         precision * Abs(eigen_values(0)));
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
   if (!force_output && problems.is_running_tachyon(MRSSM_info::VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(MRSSM_info::VWm);

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(MRSSM_info::VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_running_tachyon(MRSSM_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(MRSSM_info::VZ);

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1loop_1_heavy_rotated(
      p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL_heavy_rotated
      (p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR_heavy_rotated
      (p, idx, idx));
   const double drbar_conversion = 1 - 0.0023747152416172916*(0.6*Sqr(g1)
      - Sqr(g2));
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mf_1loop = - self_energy_1/m_sm_drbar -
      self_energy_PL - self_energy_PR;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mf_1loop);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1loop_1_heavy_rotated(
      p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fu_1loop_PL_heavy_rotated
      (p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_1loop_PR_heavy_rotated
      (p, idx, idx));

   const double currentScale = get_scale();
   double qcd_1l = 0., qcd_2l = 0., qcd_3l = 0.;

   qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(idx))/Sqr(
      currentScale)))*Sqr(g3);

   if (get_thresholds() > 1 && threshold_corrections.mt > 1) {
      const double q_2l = 0.005191204615668296*Quad(g3) +
         0.0032883224409535764*Log(Sqr(currentScale)/Sqr(MFu(idx)))*Quad(g3) +
         0.0008822328500119351*Quad(g3)*Sqr(Log(Sqr(currentScale)/Sqr(MFu(idx))
         ));

      qcd_2l = -q_2l + qcd_1l * qcd_1l;
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1loop_1_heavy_rotated(
      p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL_heavy_rotated
      (p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR_heavy_rotated
      (p, idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1 + 0.0006860288475783287*Sqr(g1) +
      0.0023747152416172916*Sqr(g2) - 0.008443431970194815*Sqr(g3);
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mb_1loop = - self_energy_1/m_tree - self_energy_PL
      - self_energy_PR;
   double qcd_2l = 0.;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mb_1loop +
      qcd_2l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(MRSSM_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm_1loop(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(MRSSM_info::VWm);
      return m_pole;
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


std::ostream& operator<<(std::ostream& ostr, const MRSSM_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
