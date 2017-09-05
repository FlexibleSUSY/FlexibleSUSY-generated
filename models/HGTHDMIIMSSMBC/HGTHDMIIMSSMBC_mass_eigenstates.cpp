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

// File generated at Tue 5 Sep 2017 10:24:36

/**
 * @file HGTHDMIIMSSMBC_mass_eigenstates.cpp
 * @brief implementation of the HGTHDMIIMSSMBC model class
 *
 * Contains the definition of the HGTHDMIIMSSMBC model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 5 Sep 2017 10:24:36 with FlexibleSUSY
 * 1.7.5 (git commit: c98e024e1e74ea3309b68f7006d5f91f8df6c678) and SARAH 4.12.0 .
 */

#include "HGTHDMIIMSSMBC_mass_eigenstates.hpp"
#include "eigen_utils.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "parallel.hpp"
#include "pv.hpp"
#include "functors.hpp"




#include <cmath>
#include <iostream>
#include <memory>
#include <algorithm>

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

#define CLASSNAME HGTHDMIIMSSMBC_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model->get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS     two_loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     two_loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     two_loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU two_loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION          two_loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS  1

CLASSNAME::HGTHDMIIMSSMBC_mass_eigenstates(const HGTHDMIIMSSMBC_input_parameters& input_)
   : HGTHDMIIMSSMBC_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , calculate_bsm_pole_masses(true)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(HGTHDMIIMSSMBC_info::particle_names)
   , two_loop_corrections()
   , MVG(0), MFv(Eigen::Array<double,3,1>::Zero()), MGlu(0), Mhh(Eigen::Array<
      double,2,1>::Zero()), MAh(Eigen::Array<double,2,1>::Zero()), MHm(
      Eigen::Array<double,2,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()),
      MFu(Eigen::Array<double,3,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero())
      , MChi(Eigen::Array<double,4,1>::Zero()), MCha(Eigen::Array<double,2,1>
      ::Zero()), MVWm(0), MVP(0), MVZ(0)

   , ZH(Eigen::Matrix<double,2,2>::Zero()), ZA(Eigen::Matrix<double,2,2>::Zero(
      )), ZP(Eigen::Matrix<double,2,2>::Zero()), Vd(Eigen::Matrix<std::complex<
      double>,3,3>::Zero()), Ud(Eigen::Matrix<std::complex<double>,3,3>::Zero()),
      Vu(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Uu(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), Ve(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), Ue(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZN(
      Eigen::Matrix<std::complex<double>,4,4>::Zero()), UM(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), UP(Eigen::Matrix<std::complex<double>,2,
      2>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~HGTHDMIIMSSMBC_mass_eigenstates()
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

void CLASSNAME::set_ewsb_loop_order(unsigned loop_order)
{
   ewsb_loop_order = loop_order;
}

void CLASSNAME::set_two_loop_corrections(const Two_loop_corrections& two_loop_corrections_)
{
   two_loop_corrections = two_loop_corrections_;
}

const Two_loop_corrections& CLASSNAME::get_two_loop_corrections() const
{
   return two_loop_corrections;
}

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

std::size_t CLASSNAME::get_number_of_ewsb_iterations() const
{
   return number_of_ewsb_iterations;
}

void CLASSNAME::set_number_of_mass_iterations(std::size_t iterations)
{
   number_of_mass_iterations = iterations;
}

std::size_t CLASSNAME::get_number_of_mass_iterations() const
{
   return number_of_mass_iterations;
}

void CLASSNAME::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
}

void CLASSNAME::set_pole_mass_loop_order(unsigned loop_order)
{
   pole_mass_loop_order = loop_order;
}

unsigned CLASSNAME::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
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

const HGTHDMIIMSSMBC_physical& CLASSNAME::get_physical() const
{
   return physical;
}

HGTHDMIIMSSMBC_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const HGTHDMIIMSSMBC_physical& physical_)
{
   physical = physical_;
}

const Problems<HGTHDMIIMSSMBC_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<HGTHDMIIMSSMBC_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
{
   return problems;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::tadpole_equations() const
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0));
      tadpole[1] -= Re(tadpole_hh(1));

      if (ewsb_loop_order > 1) {

      }
   }

   return tadpole;
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
 * Method which calculates the tadpoles at loop order specified in the
 * pointer to the CLASSNAME::EWSB_args struct.
 *
 * @param x GSL vector of EWSB output parameters
 * @param params pointer to CLASSNAME::EWSB_args struct
 * @param f GSL vector with tadpoles
 *
 * @return GSL_EDOM if x contains Nans, GSL_SUCCESS otherwise.
 */
int CLASSNAME::tadpole_equations(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   HGTHDMIIMSSMBC_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_M112(gsl_vector_get(x, 0));
   model->set_M222(gsl_vector_get(x, 1));


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   const auto tadpole(model->tadpole_equations());

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return IsFinite(tadpole) ? GSL_SUCCESS : GSL_EDOM;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   EWSB_args params = {this, ewsb_loop_order};

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, fixed_point_iterator::Convergence_tester_relative(ewsb_iteration_precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrids)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_broyden))
   };

   const std::size_t number_of_solvers = sizeof(solvers)/sizeof(*solvers);
   const auto x_init(ewsb_initial_guess());

   VERBOSE_MSG("Solving EWSB equations ...");
   VERBOSE_MSG("\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (std::size_t i = 0; i < number_of_solvers; ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i].get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\tSolver " << i << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\tSolver " << i << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int CLASSNAME::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const Eigen::Matrix<double, number_of_ewsb_equations, 1>& x_init
)
{
   const int status = solver->solve(&x_init[0]);

   M112 = solver->get_solution(0);
   M222 = solver->get_solution(1);


   return status;
}

int CLASSNAME::solve_ewsb_iteratively(unsigned loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const unsigned old_loop_order = ewsb_loop_order;
   ewsb_loop_order = loop_order;
   const int status = solve_ewsb_iteratively();
   ewsb_loop_order = old_loop_order;
   return status;
}


int CLASSNAME::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_M112 = M112;
   const double old_M222 = M222;

   M112 = Re((0.25*(-4*Lambda1*Power(v1,3) + 2*M122*v2 - Lambda7*Power(v2,3) -
      Power(v2,3)*Conj(Lambda7) + 2*v2*Conj(M122) - 3*Lambda6*v2*Sqr(v1) - 3*v2*
      Conj(Lambda6)*Sqr(v1) - 2*Lambda3*v1*Sqr(v2) - 2*Lambda4*v1*Sqr(v2) -
      Lambda5*v1*Sqr(v2) - v1*Conj(Lambda5)*Sqr(v2)))/v1);
   M222 = Re((0.25*(2*M122*v1 - Lambda6*Power(v1,3) - 4*Lambda2*Power(v2,3) -
      Power(v1,3)*Conj(Lambda6) + 2*v1*Conj(M122) - 2*Lambda3*v2*Sqr(v1) - 2*
      Lambda4*v2*Sqr(v1) - Lambda5*v2*Sqr(v1) - v2*Conj(Lambda5)*Sqr(v1) - 3*
      Lambda7*v1*Sqr(v2) - 3*v1*Conj(Lambda7)*Sqr(v2)))/v2);

   const bool is_finite = IsFinite(M112) && IsFinite(M222);

   if (!is_finite) {
      M112 = old_M112;
      M222 = old_M222;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_custom()
{
   int error = 0;



   return error;
}

int CLASSNAME::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int CLASSNAME::solve_ewsb()
{
   VERBOSE_MSG("\tSolving EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::ewsb_initial_guess()
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> x_init(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   x_init[0] = M112;
   x_init[1] = M222;


   return x_init;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @return new set of EWSB output parameters
 */
Eigen::Matrix<double, CLASSNAME::number_of_ewsb_equations, 1> CLASSNAME::ewsb_step() const
{
   double tadpole[number_of_ewsb_equations] = { 0. };
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh(0));
      tadpole[1] += Re(tadpole_hh(1));

      if (ewsb_loop_order > 1) {

      }
   }

   double M112;
   double M222;

   M112 = Re((0.25*(-4*Lambda1*Power(v1,3) + 2*M122*v2 - Lambda7*Power(v2,3) -
      Power(v2,3)*Conj(Lambda7) + 2*v2*Conj(M122) + 4*tadpole[0] - 3*Lambda6*v2*
      Sqr(v1) - 3*v2*Conj(Lambda6)*Sqr(v1) - 2*Lambda3*v1*Sqr(v2) - 2*Lambda4*v1*
      Sqr(v2) - Lambda5*v1*Sqr(v2) - v1*Conj(Lambda5)*Sqr(v2)))/v1);
   M222 = Re((0.25*(2*M122*v1 - Lambda6*Power(v1,3) - 4*Lambda2*Power(v2,3) -
      Power(v1,3)*Conj(Lambda6) + 2*v1*Conj(M122) + 4*tadpole[1] - 2*Lambda3*v2*
      Sqr(v1) - 2*Lambda4*v2*Sqr(v1) - Lambda5*v2*Sqr(v1) - v2*Conj(Lambda5)*Sqr(
      v1) - 3*Lambda7*v1*Sqr(v2) - 3*v1*Conj(Lambda7)*Sqr(v2)))/v2);

   const bool is_finite = IsFinite(M112) && IsFinite(M222);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = M112;
   ewsb_parameters[1] = M222;


   return ewsb_parameters;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param x old EWSB output parameters
 * @param params further function parameters
 * @param f new EWSB output parameters
 *
 * @return Returns status of CLASSNAME::ewsb_step
 */
int CLASSNAME::ewsb_step(const gsl_vector* x, void* params, gsl_vector* f)
{
   if (!is_finite(x)) {
      gsl_vector_set_all(f, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   HGTHDMIIMSSMBC_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double M112 = gsl_vector_get(x, 0);
   const double M222 = gsl_vector_get(x, 1);

   model->set_M112(M112);
   model->set_M222(M222);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters;
   ewsb_parameters[0] = M112;
   ewsb_parameters[1] = M222;


   int status = GSL_SUCCESS;

   try {
      ewsb_parameters = model->ewsb_step();
      status = GSL_SUCCESS;
   } catch (...) {
      status = GSL_EDOM;
   }

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "HGTHDMIIMSSMBC\n"
           "========================================\n";
   HGTHDMIIMSSMBC_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHm = " << MHm.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 */

double CLASSNAME::A0(double m) const
{
   return passarino_veltman::ReA0(m*m, Sqr(get_scale()));
}

double CLASSNAME::B0(double p, double m1, double m2) const
{
   return passarino_veltman::ReB0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B1(double p, double m1, double m2) const
{
   return passarino_veltman::ReB1(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B00(double p, double m1, double m2) const
{
   return passarino_veltman::ReB00(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::B22(double p, double m1, double m2) const
{
   return passarino_veltman::ReB22(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::H0(double p, double m1, double m2) const
{
   return passarino_veltman::ReH0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::F0(double p, double m1, double m2) const
{
   return passarino_veltman::ReF0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

double CLASSNAME::G0(double p, double m1, double m2) const
{
   return passarino_veltman::ReG0(p*p, m1*m1, m2*m2, Sqr(get_scale()));
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto old_M112 = M112;
   const auto old_M222 = M222;

   solve_ewsb_tree_level();

   calculate_MVPVZ();
   calculate_MVWm();
   calculate_MCha();
   calculate_MChi();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_MHm();
   calculate_MAh();
   calculate_Mhh();
   calculate_MGlu();
   calculate_MFv();
   calculate_MVG();

   M112 = old_M112;
   M222 = old_M222;

}

/**
 * Backward compatibility routine which finds the DRbar mass
 * eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_parameters()
{
   calculate_DRbar_masses();
}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void CLASSNAME::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   auto obj_ptr = this;

   std::future<void> fut_MVG;
   std::future<void> fut_MFv;
   std::future<void> fut_MGlu;
   std::future<void> fut_MVP;
   std::future<void> fut_MVZ;
   std::future<void> fut_Mhh;
   std::future<void> fut_MAh;
   std::future<void> fut_MHm;
   std::future<void> fut_MFd;
   std::future<void> fut_MFu;
   std::future<void> fut_MFe;
   std::future<void> fut_MChi;
   std::future<void> fut_MCha;
   std::future<void> fut_MVWm;

   if (calculate_bsm_pole_masses) {
      fut_MAh = run_async([obj_ptr] () { obj_ptr->calculate_MAh_pole(); });
      fut_MCha = run_async([obj_ptr] () { obj_ptr->calculate_MCha_pole(); });
      fut_MChi = run_async([obj_ptr] () { obj_ptr->calculate_MChi_pole(); });
      fut_MGlu = run_async([obj_ptr] () { obj_ptr->calculate_MGlu_pole(); });
      fut_Mhh = run_async([obj_ptr] () { obj_ptr->calculate_Mhh_pole(); });
      fut_MHm = run_async([obj_ptr] () { obj_ptr->calculate_MHm_pole(); });
   }

   if (calculate_sm_pole_masses) {
      fut_MVG = run_async([obj_ptr] () { obj_ptr->calculate_MVG_pole(); });
      fut_MFv = run_async([obj_ptr] () { obj_ptr->calculate_MFv_pole(); });
      fut_MVP = run_async([obj_ptr] () { obj_ptr->calculate_MVP_pole(); });
      fut_MVZ = run_async([obj_ptr] () { obj_ptr->calculate_MVZ_pole(); });
      fut_MFd = run_async([obj_ptr] () { obj_ptr->calculate_MFd_pole(); });
      fut_MFu = run_async([obj_ptr] () { obj_ptr->calculate_MFu_pole(); });
      fut_MFe = run_async([obj_ptr] () { obj_ptr->calculate_MFe_pole(); });
      fut_MVWm = run_async([obj_ptr] () { obj_ptr->calculate_MVWm_pole(); });
   }

   if (fut_MAh.valid()) fut_MAh.get();
   if (fut_MCha.valid()) fut_MCha.get();
   if (fut_MChi.valid()) fut_MChi.get();
   if (fut_MGlu.valid()) fut_MGlu.get();
   if (fut_Mhh.valid()) fut_Mhh.get();
   if (fut_MHm.valid()) fut_MHm.get();
   if (fut_MVG.valid()) fut_MVG.get();
   if (fut_MFv.valid()) fut_MFv.get();
   if (fut_MVP.valid()) fut_MVP.get();
   if (fut_MVZ.valid()) fut_MVZ.get();
   if (fut_MFd.valid()) fut_MFd.get();
   if (fut_MFu.valid()) fut_MFu.get();
   if (fut_MFe.valid()) fut_MFe.get();
   if (fut_MVWm.valid()) fut_MVWm.get();

#else
   if (calculate_bsm_pole_masses) {
      calculate_MAh_pole();
      calculate_MCha_pole();
      calculate_MChi_pole();
      calculate_MGlu_pole();
      calculate_Mhh_pole();
      calculate_MHm_pole();
   }

   if (calculate_sm_pole_masses) {
      calculate_MVG_pole();
      calculate_MFv_pole();
      calculate_MVP_pole();
      calculate_MVZ_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MFe_pole();
      calculate_MVWm_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHm) = MHm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MCha) = MCha;
   PHYSICAL(UM) = UM;
   PHYSICAL(UP) = UP;
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
   move_goldstone_to(0, MVWm, MHm, ZP);

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
   move_goldstone_to(0, MVWm, PHYSICAL(MHm), PHYSICAL(ZP));

}

/**
 * Checks the pole masses for tachyons
 */
void CLASSNAME::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(Mhh).tail<2>().minCoeff() < 0.) problems.flag_tachyon(HGTHDMIIMSSMBC_info::hh);
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) problems.flag_tachyon(HGTHDMIIMSSMBC_info::Ah);
   if (PHYSICAL(MHm).tail<1>().minCoeff() < 0.) problems.flag_tachyon(HGTHDMIIMSSMBC_info::Hm);

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
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MGlu = 0.;
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;

   PhaseGlu = std::complex<double>(1.,0.);

}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   HGTHDMIIMSSMBC_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

void CLASSNAME::set_DRbar_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MFv(0) = pars(1);
   MFv(1) = pars(2);
   MFv(2) = pars(3);
   MGlu = pars(4);
   Mhh(0) = pars(5);
   Mhh(1) = pars(6);
   MAh(0) = pars(7);
   MAh(1) = pars(8);
   MHm(0) = pars(9);
   MHm(1) = pars(10);
   MFd(0) = pars(11);
   MFd(1) = pars(12);
   MFd(2) = pars(13);
   MFu(0) = pars(14);
   MFu(1) = pars(15);
   MFu(2) = pars(16);
   MFe(0) = pars(17);
   MFe(1) = pars(18);
   MFe(2) = pars(19);
   MChi(0) = pars(20);
   MChi(1) = pars(21);
   MChi(2) = pars(22);
   MChi(3) = pars(23);
   MCha(0) = pars(24);
   MCha(1) = pars(25);
   MVWm = pars(26);
   MVP = pars(27);
   MVZ = pars(28);

}

Eigen::ArrayXd CLASSNAME::get_DRbar_masses() const
{
   Eigen::ArrayXd pars(29);

   pars(0) = MVG;
   pars(1) = MFv(0);
   pars(2) = MFv(1);
   pars(3) = MFv(2);
   pars(4) = MGlu;
   pars(5) = Mhh(0);
   pars(6) = Mhh(1);
   pars(7) = MAh(0);
   pars(8) = MAh(1);
   pars(9) = MHm(0);
   pars(10) = MHm(1);
   pars(11) = MFd(0);
   pars(12) = MFd(1);
   pars(13) = MFd(2);
   pars(14) = MFu(0);
   pars(15) = MFu(1);
   pars(16) = MFu(2);
   pars(17) = MFe(0);
   pars(18) = MFe(1);
   pars(19) = MFe(2);
   pars(20) = MChi(0);
   pars(21) = MChi(1);
   pars(22) = MChi(2);
   pars(23) = MChi(3);
   pars(24) = MCha(0);
   pars(25) = MCha(1);
   pars(26) = MVWm;
   pars(27) = MVP;
   pars(28) = MVZ;

   return pars;
}

std::string CLASSNAME::name() const
{
   return "HGTHDMIIMSSMBC";
}

void CLASSNAME::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   HGTHDMIIMSSMBC_soft_parameters::run_to(scale, eps);
}


Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHm_goldstone;
   MHm_goldstone(0) = MVWm;

   return remove_if_equal(MHm, MHm_goldstone);
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_goldstone;
   MAh_goldstone(0) = MVZ;

   return remove_if_equal(MAh, MAh_goldstone);
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

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = M112 + 1.5*Lambda6*v1*v2 + 1.5*v1*v2*Conj(
      Lambda6) + 3*Lambda1*Sqr(v1) + 0.5*Lambda3*Sqr(v2) + 0.5*Lambda4*Sqr(v2)
      + 0.25*Lambda5*Sqr(v2) + 0.25*Conj(Lambda5)*Sqr(v2);
   mass_matrix_hh(0,1) = -0.5*M122 + Lambda3*v1*v2 + Lambda4*v1*v2 + 0.5*
      Lambda5*v1*v2 + 0.5*v1*v2*Conj(Lambda5) - 0.5*Conj(M122) + 0.75*Lambda6*
      Sqr(v1) + 0.75*Conj(Lambda6)*Sqr(v1) + 0.75*Lambda7*Sqr(v2) + 0.75*Conj(
      Lambda7)*Sqr(v2);
   mass_matrix_hh(1,1) = M222 + 1.5*Lambda7*v1*v2 + 1.5*v1*v2*Conj(
      Lambda7) + 0.5*Lambda3*Sqr(v1) + 0.5*Lambda4*Sqr(v1) + 0.25*Lambda5*Sqr(
      v1) + 0.25*Conj(Lambda5)*Sqr(v1) + 3*Lambda2*Sqr(v2);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::hh, eigenvalue_error >
      precision * Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif


   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::hh);
   }

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = M112 + 0.5*Lambda6*v1*v2 + 0.5*v1*v2*Conj(
      Lambda6) + Lambda1*Sqr(v1) + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(
      ThetaW())*Sqr(v1) + 0.5*Lambda3*Sqr(v2) + 0.5*Lambda4*Sqr(v2) - 0.25*
      Lambda5*Sqr(v2) - 0.25*Conj(Lambda5)*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v1)*Sqr(
      Cos(ThetaW())) + 0.15*Sqr(g1)*Sqr(v1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = -0.5*M122 + 0.5*Lambda5*v1*v2 + 0.5*v1*v2*Conj(
      Lambda5) - 0.5*Conj(M122) + 0.3872983346207417*g1*g2*v1*v2*Cos(ThetaW())*
      Sin(ThetaW()) + 0.25*Lambda6*Sqr(v1) + 0.25*Conj(Lambda6)*Sqr(v1) + 0.25*
      Lambda7*Sqr(v2) + 0.25*Conj(Lambda7)*Sqr(v2) + 0.25*v1*v2*Sqr(g2)*Sqr(Cos
      (ThetaW())) + 0.15*v1*v2*Sqr(g1)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,1) = M222 + 0.5*Lambda7*v1*v2 + 0.5*v1*v2*Conj(
      Lambda7) + 0.5*Lambda3*Sqr(v1) + 0.5*Lambda4*Sqr(v1) - 0.25*Lambda5*Sqr(
      v1) - 0.25*Conj(Lambda5)*Sqr(v1) + Lambda2*Sqr(v2) + 0.3872983346207417*
      g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v2)*Sqr(Cos(
      ThetaW())) + 0.15*Sqr(g1)*Sqr(v2)*Sqr(Sin(ThetaW()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Ah, eigenvalue_error >
      precision * Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif


   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::Ah);
   }

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hm;

   mass_matrix_Hm(0,0) = M112 + 0.5*Lambda6*v1*v2 + 0.5*v1*v2*Conj(
      Lambda6) + Lambda1*Sqr(v1) + 0.25*Sqr(g2)*Sqr(v1) + 0.5*Lambda3*Sqr(v2);
   mass_matrix_Hm(0,1) = 0.5*Lambda4*v1*v2 + 0.5*Lambda5*v1*v2 - Conj(
      M122) + 0.25*v1*v2*Sqr(g2) + 0.5*Conj(Lambda6)*Sqr(v1) + 0.5*Conj(Lambda7
      )*Sqr(v2);
   mass_matrix_Hm(1,0) = -M122 + 0.5*Lambda4*v1*v2 + 0.5*v1*v2*Conj(
      Lambda5) + 0.25*v1*v2*Sqr(g2) + 0.5*Lambda6*Sqr(v1) + 0.5*Lambda7*Sqr(v2)
      ;
   mass_matrix_Hm(1,1) = M222 + 0.5*Lambda7*v1*v2 + 0.5*v1*v2*Conj(
      Lambda7) + 0.5*Lambda3*Sqr(v1) + Lambda2*Sqr(v2) + 0.25*Sqr(g2)*Sqr(v2);

   return mass_matrix_Hm;
}

void CLASSNAME::calculate_MHm()
{
   const auto mass_matrix_Hm(get_mass_matrix_Hm());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hm, MHm, ZP, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Hm, eigenvalue_error >
      precision * Abs(MHm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hm, MHm, ZP);
#endif


   if (MHm.minCoeff() < 0.) {
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::Hm);
   }

   MHm = AbsSqrt(MHm);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v1*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v1*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v1*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v1*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v1*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v1*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v1*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v1*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v1*Yd(2,2);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Fd, eigenvalue_error >
      precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = -0.7071067811865475*v2*Yu(0,0);
   mass_matrix_Fu(0,1) = -0.7071067811865475*v2*Yu(1,0);
   mass_matrix_Fu(0,2) = -0.7071067811865475*v2*Yu(2,0);
   mass_matrix_Fu(1,0) = -0.7071067811865475*v2*Yu(0,1);
   mass_matrix_Fu(1,1) = -0.7071067811865475*v2*Yu(1,1);
   mass_matrix_Fu(1,2) = -0.7071067811865475*v2*Yu(2,1);
   mass_matrix_Fu(2,0) = -0.7071067811865475*v2*Yu(0,2);
   mass_matrix_Fu(2,1) = -0.7071067811865475*v2*Yu(1,2);
   mass_matrix_Fu(2,2) = -0.7071067811865475*v2*Yu(2,2);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Fu, eigenvalue_error >
      precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v1*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v1*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v1*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v1*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v1*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v1*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v1*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v1*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v1*Ye(2,2);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Fe, eigenvalue_error >
      precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Chi() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.5*g1dp*v1;
   mass_matrix_Chi(0,3) = 0.5*g2up*v2;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g1d*v1;
   mass_matrix_Chi(1,3) = -0.5*g2u*v2;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -Mu;
   mass_matrix_Chi(3,3) = 0;

   Symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Chi, eigenvalue_error >
      precision * Abs(MChi(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif

}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2u*v2;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g1d*v1;
   mass_matrix_Cha(1,1) = Mu;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Cha, MCha, UM, UP, eigenvalue_error);
   problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Cha, eigenvalue_error >
      precision * Abs(MCha(0)));
#else
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
#endif

}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = Re(0.25*Sqr(g2)*(Sqr(v1) + Sqr(v2)));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = mass_matrix_VWm;

   if (MVWm < 0.) {
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VWm);
   }

   MVWm = AbsSqrt(MVWm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_VPVZ() const
{
   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v1) + 0.15*Sqr(g1)*Sqr(v2);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v1) -
      0.19364916731037085*g1*g2*Sqr(v2);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v1) + 0.25*Sqr(g2)*Sqr(v2);

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
   ZZ.transposeInPlace();
#else
   fs_diagonalize_hermitian(mass_matrix_VPVZ, MVPVZ, ZZ);
   ZZ.transposeInPlace();
#endif


   MVPVZ = AbsSqrt(MVPVZ);

   MVP = 0.;
   MVZ = MVPVZ(1);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = Re(M112*v1 + Lambda1*Power(v1,3) - 0.5*M122*v2 + 0.25*
      Lambda7*Power(v2,3) + 0.25*Power(v2,3)*Conj(Lambda7) - 0.5*v2*Conj(M122) +
      0.75*Lambda6*v2*Sqr(v1) + 0.75*v2*Conj(Lambda6)*Sqr(v1) + 0.5*Lambda3*v1*Sqr
      (v2) + 0.5*Lambda4*v1*Sqr(v2) + 0.25*Lambda5*v1*Sqr(v2) + 0.25*v1*Conj(
      Lambda5)*Sqr(v2));

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = Re(-0.5*M122*v1 + 0.25*Lambda6*Power(v1,3) + M222*v2 +
      Lambda2*Power(v2,3) + 0.25*Power(v1,3)*Conj(Lambda6) - 0.5*v1*Conj(M122) +
      0.5*Lambda3*v2*Sqr(v1) + 0.5*Lambda4*v2*Sqr(v1) + 0.25*Lambda5*v2*Sqr(v1) +
      0.25*v2*Conj(Lambda5)*Sqr(v1) + 0.75*Lambda7*v1*Sqr(v2) + 0.75*v1*Conj(
      Lambda7)*Sqr(v2));

   return result;
}



std::complex<double> CLASSNAME::CpUhhVZVZ(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(v1*KroneckerDelta(0,gO2) + v2*KroneckerDelta(1,gO2))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(v1*KroneckerDelta(0,gO2) + v2*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(v1*KroneckerDelta(0,gO1) + v2*KroneckerDelta(1,gO1))*Sqr(g2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(v1*KroneckerDelta(0,gO1) + v2*KroneckerDelta(1,gO1))*Sqr(g2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargZgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.025*(v1*KroneckerDelta(0,gO1) + v2*KroneckerDelta(1,gO1))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*(ZA(gI1,1)*((
      Lambda7 + Conj(Lambda7))*ZA(gI2,0) + 4*Lambda2*ZA(gI2,1)) + ZA(gI1,0)*((2*
      Lambda3 + 2*Lambda4 - Lambda5 - Conj(Lambda5))*ZA(gI2,0) + (Lambda7 + Conj(
      Lambda7))*ZA(gI2,1))) + KroneckerDelta(0,gO2)*(ZA(gI1,0)*((Lambda6 + Conj(
      Lambda6))*ZA(gI2,0) + (Lambda5 + Conj(Lambda5))*ZA(gI2,1)) + ZA(gI1,1)*((
      Lambda5 + Conj(Lambda5))*ZA(gI2,0) + (Lambda7 + Conj(Lambda7))*ZA(gI2,1)))))
      - KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*(ZA(gI1,1)*((Lambda6 + Conj(
      Lambda6))*ZA(gI2,0) + (2*Lambda3 + 2*Lambda4 - Lambda5 - Conj(Lambda5))*ZA(
      gI2,1)) + ZA(gI1,0)*(4*Lambda1*ZA(gI2,0) + (Lambda6 + Conj(Lambda6))*ZA(gI2,
      1))) + KroneckerDelta(1,gO2)*(ZA(gI1,0)*((Lambda6 + Conj(Lambda6))*ZA(gI2,0)
      + (Lambda5 + Conj(Lambda5))*ZA(gI2,1)) + ZA(gI1,1)*((Lambda5 + Conj(Lambda5
      ))*ZA(gI2,0) + (Lambda7 + Conj(Lambda7))*ZA(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjHmHm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO1)*(2*KroneckerDelta(1,gO2)*(ZP(gI1,1)*(
      Lambda7*ZP(gI2,0) + 2*Lambda2*ZP(gI2,1)) + ZP(gI1,0)*(Lambda3*ZP(gI2,0) +
      Conj(Lambda7)*ZP(gI2,1))) + KroneckerDelta(0,gO2)*(ZP(gI1,0)*((Lambda6 +
      Conj(Lambda6))*ZP(gI2,0) + (Lambda4 + Lambda5)*ZP(gI2,1)) + ZP(gI1,1)*((
      Lambda4 + Conj(Lambda5))*ZP(gI2,0) + (Lambda7 + Conj(Lambda7))*ZP(gI2,1)))))
      - KroneckerDelta(0,gO1)*(2*KroneckerDelta(0,gO2)*(ZP(gI1,1)*(Lambda6*ZP(gI2
      ,0) + Lambda3*ZP(gI2,1)) + ZP(gI1,0)*(2*Lambda1*ZP(gI2,0) + Conj(Lambda6)*ZP
      (gI2,1))) + KroneckerDelta(1,gO2)*(ZP(gI1,0)*((Lambda6 + Conj(Lambda6))*ZP(
      gI2,0) + (Lambda4 + Lambda5)*ZP(gI2,1)) + ZP(gI1,1)*((Lambda4 + Conj(Lambda5
      ))*ZP(gI2,0) + (Lambda7 + Conj(Lambda7))*ZP(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*(3*ZH(gI1,1)*((
      Lambda7 + Conj(Lambda7))*ZH(gI2,0) + 4*Lambda2*ZH(gI2,1)) + ZH(gI1,0)*((2*
      Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gI2,0) + 3*(Lambda7 + Conj
      (Lambda7))*ZH(gI2,1))) + KroneckerDelta(0,gO2)*(ZH(gI1,0)*(3*(Lambda6 + Conj
      (Lambda6))*ZH(gI2,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(
      gI2,1)) + ZH(gI1,1)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(
      gI2,0) + 3*(Lambda7 + Conj(Lambda7))*ZH(gI2,1))))) - KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*(ZH(gI1,1)*(3*(Lambda6 + Conj(Lambda6))*ZH(gI2,0) + (2
      *Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gI2,1)) + 3*ZH(gI1,0)*(4*
      Lambda1*ZH(gI2,0) + (Lambda6 + Conj(Lambda6))*ZH(gI2,1))) + KroneckerDelta(1
      ,gO2)*(ZH(gI1,0)*(3*(Lambda6 + Conj(Lambda6))*ZH(gI2,0) + (2*Lambda3 + 2*
      Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gI2,1)) + ZH(gI1,1)*((2*Lambda3 + 2*
      Lambda4 + Lambda5 + Conj(Lambda5))*ZH(gI2,0) + 3*(Lambda7 + Conj(Lambda7))*
      ZH(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(0,gO2)*(ZA(gI1,0)*((4*Lambda1*v1 + Lambda6*v2
      + v2*Conj(Lambda6))*ZA(gI2,0) + (Lambda6*v1 + Lambda5*v2 + v2*Conj(Lambda5)
      + v1*Conj(Lambda6))*ZA(gI2,1)) + ZA(gI1,1)*((Lambda6*v1 + Lambda5*v2 + v2*
      Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gI2,0) + (2*Lambda3*v1 + 2*Lambda4*v1 -
      Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZA(gI2,1))))
      - KroneckerDelta(1,gO2)*(ZA(gI1,1)*((Lambda5*v1 + Lambda7*v2 + v1*Conj(
      Lambda5) + v2*Conj(Lambda7))*ZA(gI2,0) + (Lambda7*v1 + 4*Lambda2*v2 + v1*
      Conj(Lambda7))*ZA(gI2,1)) + ZA(gI1,0)*((Lambda6*v1 + 2*Lambda3*v2 + 2*
      Lambda4*v2 - Lambda5*v2 - v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gI2,0) + (
      Lambda5*v1 + Lambda7*v2 + v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZA(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjHmHm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(0,gO2)*(ZP(gI1,0)*((4*Lambda1*v1 + Lambda6*v2
      + v2*Conj(Lambda6))*ZP(gI2,0) + ((Lambda4 + Lambda5)*v2 + 2*v1*Conj(Lambda6
      ))*ZP(gI2,1)) + ZP(gI1,1)*((2*Lambda6*v1 + Lambda4*v2 + v2*Conj(Lambda5))*ZP
      (gI2,0) + (2*Lambda3*v1 + Lambda7*v2 + v2*Conj(Lambda7))*ZP(gI2,1)))) -
      KroneckerDelta(1,gO2)*(ZP(gI1,1)*((Lambda4*v1 + 2*Lambda7*v2 + v1*Conj(
      Lambda5))*ZP(gI2,0) + (Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZP(gI2,
      1)) + ZP(gI1,0)*((Lambda6*v1 + 2*Lambda3*v2 + v1*Conj(Lambda6))*ZP(gI2,0) +
      ((Lambda4 + Lambda5)*v1 + 2*v2*Conj(Lambda7))*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(1,gO2)*(ZA(gI2,1)*((
      Lambda5*v1 - Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gI1,0) +
      v1*(-Lambda7 + Conj(Lambda7))*ZH(gI1,1)) + ZA(gI2,0)*((Lambda6*v1 - Lambda5*
      v2 + v2*Conj(Lambda5) - v1*Conj(Lambda6))*ZH(gI1,0) + (-(Lambda5*v1) + 3*
      Lambda7*v2 + v1*Conj(Lambda5) - 3*v2*Conj(Lambda7))*ZH(gI1,1))) +
      KroneckerDelta(0,gO2)*(ZA(gI2,0)*(v2*(Lambda6 - Conj(Lambda6))*ZH(gI1,0) + (
      Lambda6*v1 - Lambda5*v2 + v2*Conj(Lambda5) - v1*Conj(Lambda6))*ZH(gI1,1)) +
      ZA(gI2,1)*((-3*Lambda6*v1 + Lambda5*v2 - v2*Conj(Lambda5) + 3*v1*Conj(
      Lambda6))*ZH(gI1,0) + (Lambda5*v1 - Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(
      Lambda7))*ZH(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO2)*(ZH(gI1,1)*((2*Lambda3*v1 + 2*Lambda4*
      v1 + Lambda5*v1 + 3*Lambda7*v2 + v1*Conj(Lambda5) + 3*v2*Conj(Lambda7))*ZH(
      gI2,0) + 3*(Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZH(gI2,1)) + ZH(
      gI1,0)*((3*Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 + Lambda5*v2 + v2*Conj(
      Lambda5) + 3*v1*Conj(Lambda6))*ZH(gI2,0) + (2*Lambda3*v1 + 2*Lambda4*v1 +
      Lambda5*v1 + 3*Lambda7*v2 + v1*Conj(Lambda5) + 3*v2*Conj(Lambda7))*ZH(gI2,1)
      ))) - KroneckerDelta(0,gO2)*(ZH(gI1,0)*(3*(4*Lambda1*v1 + Lambda6*v2 + v2*
      Conj(Lambda6))*ZH(gI2,0) + (3*Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 +
      Lambda5*v2 + v2*Conj(Lambda5) + 3*v1*Conj(Lambda6))*ZH(gI2,1)) + ZH(gI1,1)*(
      (3*Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 + Lambda5*v2 + v2*Conj(Lambda5)
      + 3*v1*Conj(Lambda6))*ZH(gI2,0) + (2*Lambda3*v1 + 2*Lambda4*v1 + Lambda5*v1
      + 3*Lambda7*v2 + v1*Conj(Lambda5) + 3*v2*Conj(Lambda7))*ZH(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g1d*KroneckerDelta(0,gO2)*UM(gI1,1)*UP(gI2,0)
      + g2u*KroneckerDelta(1,gO2)*UM(gI1,0)*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g1d*Conj(UM(gI2,1))*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) + g2u*Conj(UM(gI2,0))*Conj(UP(gI1,1))*KroneckerDelta(1
      ,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_0;
   std::complex<double> tmp_1;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2;
      std::complex<double> tmp_3;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_2 += tmp_3;
      tmp_1 += (Vd(gI1,j2)) * tmp_2;
   }
   tmp_0 += tmp_1;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_0;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_6;
      std::complex<double> tmp_7;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_7 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_6 += tmp_7;
      tmp_5 += (Conj(Vd(gI2,j2))) * tmp_6;
   }
   tmp_4 += tmp_5;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_4;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_8;
   std::complex<double> tmp_9;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_10;
      std::complex<double> tmp_11;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_11 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_10 += tmp_11;
      tmp_9 += (Ve(gI1,j2)) * tmp_10;
   }
   tmp_8 += tmp_9;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_8;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_12;
   std::complex<double> tmp_13;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_14;
      std::complex<double> tmp_15;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_15 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_14 += tmp_15;
      tmp_13 += (Conj(Ve(gI2,j2))) * tmp_14;
   }
   tmp_12 += tmp_13;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_12;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_16;
   std::complex<double> tmp_17;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_18;
      std::complex<double> tmp_19;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_19 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_18 += tmp_19;
      tmp_17 += (Vu(gI1,j2)) * tmp_18;
   }
   tmp_16 += tmp_17;
   result += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_16;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_20;
   std::complex<double> tmp_21;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_22;
      std::complex<double> tmp_23;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_23 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_22 += tmp_23;
      tmp_21 += (Conj(Vu(gI2,j2))) * tmp_22;
   }
   tmp_20 += tmp_21;
   result += (0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_20;

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO2)*(ZN(gI1,2)*(g1dp*ZN(gI2,0) - g1d*ZN(gI2,
      1)) + (g1dp*ZN(gI1,0) - g1d*ZN(gI1,1))*ZN(gI2,2)) + KroneckerDelta(1,gO2)*(
      ZN(gI1,3)*(-(g2up*ZN(gI2,0)) + g2u*ZN(gI2,1)) + (-(g2up*ZN(gI1,0)) + g2u*ZN(
      gI1,1))*ZN(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZN(gI1,2))*(g1dp*Conj(ZN(gI2,0)) - g1d*Conj(ZN(gI2,1)))*
      KroneckerDelta(0,gO1) - g1d*Conj(ZN(gI1,1))*Conj(ZN(gI2,2))*KroneckerDelta(0
      ,gO1) - g2up*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) + g2u*
      Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) + g2u*Conj(ZN(gI1,1))*
      Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) + Conj(ZN(gI1,0))*(g1dp*Conj(ZN(gI2,2)
      )*KroneckerDelta(0,gO1) - g2up*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZA(gI2,0) + KroneckerDelta(1,gO2)*
      ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmHm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) + KroneckerDelta(1,gO2)*ZP
      (gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(v1*KroneckerDelta(0,gO1) + v2*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(v1*KroneckerDelta(0,gO1) + v2*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*(3*ZA(gI1,1)*((
      Lambda7 + Conj(Lambda7))*ZA(gI2,0) + 4*Lambda2*ZA(gI2,1)) + ZA(gI1,0)*((2*
      Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gI2,0) + 3*(Lambda7 + Conj
      (Lambda7))*ZA(gI2,1))) + KroneckerDelta(0,gO2)*(ZA(gI1,0)*(3*(Lambda6 + Conj
      (Lambda6))*ZA(gI2,0) + (2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(
      gI2,1)) + ZA(gI1,1)*((2*Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(
      gI2,0) + 3*(Lambda7 + Conj(Lambda7))*ZA(gI2,1))))) - KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*(ZA(gI1,1)*(3*(Lambda6 + Conj(Lambda6))*ZA(gI2,0) + (2
      *Lambda3 + 2*Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gI2,1)) + 3*ZA(gI1,0)*(4*
      Lambda1*ZA(gI2,0) + (Lambda6 + Conj(Lambda6))*ZA(gI2,1))) + KroneckerDelta(1
      ,gO2)*(ZA(gI1,0)*(3*(Lambda6 + Conj(Lambda6))*ZA(gI2,0) + (2*Lambda3 + 2*
      Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gI2,1)) + ZA(gI1,1)*((2*Lambda3 + 2*
      Lambda4 + Lambda5 + Conj(Lambda5))*ZA(gI2,0) + 3*(Lambda7 + Conj(Lambda7))*
      ZA(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjHmHm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO1)*(2*KroneckerDelta(1,gO2)*(ZP(gI1,1)*(
      Lambda7*ZP(gI2,0) + 2*Lambda2*ZP(gI2,1)) + ZP(gI1,0)*(Lambda3*ZP(gI2,0) +
      Conj(Lambda7)*ZP(gI2,1))) + KroneckerDelta(0,gO2)*(ZP(gI1,0)*((Lambda6 +
      Conj(Lambda6))*ZP(gI2,0) + (Lambda4 + Lambda5)*ZP(gI2,1)) + ZP(gI1,1)*((
      Lambda4 + Conj(Lambda5))*ZP(gI2,0) + (Lambda7 + Conj(Lambda7))*ZP(gI2,1)))))
      - KroneckerDelta(0,gO1)*(2*KroneckerDelta(0,gO2)*(ZP(gI1,1)*(Lambda6*ZP(gI2
      ,0) + Lambda3*ZP(gI2,1)) + ZP(gI1,0)*(2*Lambda1*ZP(gI2,0) + Conj(Lambda6)*ZP
      (gI2,1))) + KroneckerDelta(1,gO2)*(ZP(gI1,0)*((Lambda6 + Conj(Lambda6))*ZP(
      gI2,0) + (Lambda4 + Lambda5)*ZP(gI2,1)) + ZP(gI1,1)*((Lambda4 + Conj(Lambda5
      ))*ZP(gI2,0) + (Lambda7 + Conj(Lambda7))*ZP(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*(ZH(gI1,1)*((
      Lambda7 + Conj(Lambda7))*ZH(gI2,0) + 4*Lambda2*ZH(gI2,1)) + ZH(gI1,0)*((2*
      Lambda3 + 2*Lambda4 - Lambda5 - Conj(Lambda5))*ZH(gI2,0) + (Lambda7 + Conj(
      Lambda7))*ZH(gI2,1))) + KroneckerDelta(0,gO2)*(ZH(gI1,0)*((Lambda6 + Conj(
      Lambda6))*ZH(gI2,0) + (Lambda5 + Conj(Lambda5))*ZH(gI2,1)) + ZH(gI1,1)*((
      Lambda5 + Conj(Lambda5))*ZH(gI2,0) + (Lambda7 + Conj(Lambda7))*ZH(gI2,1)))))
      - KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*(ZH(gI1,1)*((Lambda6 + Conj(
      Lambda6))*ZH(gI2,0) + (2*Lambda3 + 2*Lambda4 - Lambda5 - Conj(Lambda5))*ZH(
      gI2,1)) + ZH(gI1,0)*(4*Lambda1*ZH(gI2,0) + (Lambda6 + Conj(Lambda6))*ZH(gI2,
      1))) + KroneckerDelta(1,gO2)*(ZH(gI1,0)*((Lambda6 + Conj(Lambda6))*ZH(gI2,0)
      + (Lambda5 + Conj(Lambda5))*ZH(gI2,1)) + ZH(gI1,1)*((Lambda5 + Conj(Lambda5
      ))*ZH(gI2,0) + (Lambda7 + Conj(Lambda7))*ZH(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(1,gO2)*(ZA(gI1,1)*((
      Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) - v2*Conj(Lambda7))*ZA(gI2,0) + 3
      *v1*(-Lambda7 + Conj(Lambda7))*ZA(gI2,1)) + ZA(gI1,0)*((-(Lambda6*v1) -
      Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gI2,0) + (Lambda5*v1 +
      Lambda7*v2 - v1*Conj(Lambda5) - v2*Conj(Lambda7))*ZA(gI2,1))) +
      KroneckerDelta(0,gO2)*(ZA(gI1,0)*(3*v2*(Lambda6 - Conj(Lambda6))*ZA(gI2,0) +
      (-(Lambda6*v1) - Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZA(gI2,1
      )) + ZA(gI1,1)*((-(Lambda6*v1) - Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(
      Lambda6))*ZA(gI2,0) + (Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) - v2*Conj(
      Lambda7))*ZA(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjHmHm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(v2*KroneckerDelta(0,gO2) - v1*
      KroneckerDelta(1,gO2))*(ZP(gI1,0)*((Lambda6 - Conj(Lambda6))*ZP(gI2,0) + (
      Lambda4 - Lambda5)*ZP(gI2,1)) + ZP(gI1,1)*((-Lambda4 + Conj(Lambda5))*ZP(gI2
      ,0) + (Lambda7 - Conj(Lambda7))*ZP(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO2)*(ZA(gI2,1)*((2*Lambda3*v1 + 2*Lambda4*
      v1 - Lambda5*v1 + Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gI1,0
      ) + (Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(Lambda7))*ZH(gI1,1)) + ZA(gI2,0)*((
      Lambda6*v1 + Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZH(gI1,0) + (
      Lambda5*v1 + Lambda7*v2 + v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gI1,1))))
      - KroneckerDelta(0,gO2)*(ZA(gI2,0)*((4*Lambda1*v1 + Lambda6*v2 + v2*Conj(
      Lambda6))*ZH(gI1,0) + (Lambda6*v1 + 2*Lambda3*v2 + 2*Lambda4*v2 - Lambda5*v2
      - v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZH(gI1,1)) + ZA(gI2,1)*((Lambda6*v1
      + Lambda5*v2 + v2*Conj(Lambda5) + v1*Conj(Lambda6))*ZH(gI1,0) + (Lambda5*v1
      + Lambda7*v2 + v1*Conj(Lambda5) + v2*Conj(Lambda7))*ZH(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(0,gO2)*(ZH(gI1,0)*(v2*(
      Lambda6 - Conj(Lambda6))*ZH(gI2,0) + (Lambda6*v1 - Lambda5*v2 + v2*Conj(
      Lambda5) - v1*Conj(Lambda6))*ZH(gI2,1)) + ZH(gI1,1)*((Lambda6*v1 - Lambda5*
      v2 + v2*Conj(Lambda5) - v1*Conj(Lambda6))*ZH(gI2,0) + (-(Lambda5*v1) + 3*
      Lambda7*v2 + v1*Conj(Lambda5) - 3*v2*Conj(Lambda7))*ZH(gI2,1))) +
      KroneckerDelta(1,gO2)*(ZH(gI1,1)*((Lambda5*v1 - Lambda7*v2 - v1*Conj(Lambda5
      ) + v2*Conj(Lambda7))*ZH(gI2,0) + v1*(-Lambda7 + Conj(Lambda7))*ZH(gI2,1)) +
      ZH(gI1,0)*((-3*Lambda6*v1 + Lambda5*v2 - v2*Conj(Lambda5) + 3*v1*Conj(
      Lambda6))*ZH(gI2,0) + (Lambda5*v1 - Lambda7*v2 - v1*Conj(Lambda5) + v2*Conj(
      Lambda7))*ZH(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g1d*KroneckerDelta(0,
      gO2)*UM(gI1,1)*UP(gI2,0) - g2u*KroneckerDelta(1,gO2)*UM(gI1,0)*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g1d*Conj(UM(gI2,1))*
      Conj(UP(gI1,0))*KroneckerDelta(0,gO1) - g2u*Conj(UM(gI2,0))*Conj(UP(gI1,1))*
      KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_24;
   std::complex<double> tmp_25;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_26;
      std::complex<double> tmp_27;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_27 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_26 += tmp_27;
      tmp_25 += (Vd(gI1,j2)) * tmp_26;
   }
   tmp_24 += tmp_25;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO2
      )) * tmp_24;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_28;
   std::complex<double> tmp_29;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_30;
      std::complex<double> tmp_31;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_31 += Conj(Ud(gI1,j1))*Yd(j1,j2);
      }
      tmp_30 += tmp_31;
      tmp_29 += (Conj(Vd(gI2,j2))) * tmp_30;
   }
   tmp_28 += tmp_29;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_28;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_32;
   std::complex<double> tmp_33;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_34;
      std::complex<double> tmp_35;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_35 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_34 += tmp_35;
      tmp_33 += (Ve(gI1,j2)) * tmp_34;
   }
   tmp_32 += tmp_33;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO2
      )) * tmp_32;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_36;
   std::complex<double> tmp_37;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_38;
      std::complex<double> tmp_39;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_39 += Conj(Ue(gI1,j1))*Ye(j1,j2);
      }
      tmp_38 += tmp_39;
      tmp_37 += (Conj(Ve(gI2,j2))) * tmp_38;
   }
   tmp_36 += tmp_37;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_36;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_40;
   std::complex<double> tmp_41;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_42;
      std::complex<double> tmp_43;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_43 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_42 += tmp_43;
      tmp_41 += (Vu(gI1,j2)) * tmp_42;
   }
   tmp_40 += tmp_41;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO2
      )) * tmp_40;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_44;
   std::complex<double> tmp_45;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_46;
      std::complex<double> tmp_47;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_47 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_46 += tmp_47;
      tmp_45 += (Conj(Vu(gI2,j2))) * tmp_46;
   }
   tmp_44 += tmp_45;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO1)
      ) * tmp_44;

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(0,gO2)*(ZN(gI1,2)*(-(
      g1dp*ZN(gI2,0)) + g1d*ZN(gI2,1)) + (-(g1dp*ZN(gI1,0)) + g1d*ZN(gI1,1))*ZN(
      gI2,2)) + KroneckerDelta(1,gO2)*(ZN(gI1,3)*(-(g2up*ZN(gI2,0)) + g2u*ZN(gI2,1
      )) + (-(g2up*ZN(gI1,0)) + g2u*ZN(gI1,1))*ZN(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(ZN(gI1,2))*(g1dp*Conj(ZN(gI2,0))
      - g1d*Conj(ZN(gI2,1)))*KroneckerDelta(0,gO1) - g1d*Conj(ZN(gI1,1))*Conj(ZN(
      gI2,2))*KroneckerDelta(0,gO1) + g2up*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))*
      KroneckerDelta(1,gO1) - g2u*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta(1
      ,gO1) - g2u*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) + Conj(ZN(
      gI1,0))*(g1dp*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) + g2up*Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhVZhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZH(gI2,0) + KroneckerDelta(1,gO2)*
      ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjVWmHm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVWmVP(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*Cos(ThetaW())*(v1*KroneckerDelta(0,gO2) +
      v2*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVZVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*(v1*KroneckerDelta(0,gO2) + v2*
      KroneckerDelta(1,gO2))*Sin(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmbargWmCgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.05*g2*(v1*KroneckerDelta(0,gO1) + v2*KroneckerDelta(1,gO1))*(5*
      g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHmgWmCbargZ(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*g2*(v1*KroneckerDelta(0,gO2) + v2*KroneckerDelta(1,gO2))*(5*g2
      *Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmbargZgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = 0.05*g2*(v1*KroneckerDelta(0,gO1) + v2*KroneckerDelta(1,gO1))*(5*g2
      *Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHmgZbargWm(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.05*g2*(v1*KroneckerDelta(0,gO2) + v2*KroneckerDelta(1,gO2))*(5*
      g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*(ZA(gI1,1)*((
      Lambda4 + Conj(Lambda5))*ZA(gI2,0) + 2*Lambda7*ZA(gI2,1)) + ZA(gI1,0)*(2*
      Lambda6*ZA(gI2,0) + (Lambda4 + Conj(Lambda5))*ZA(gI2,1))) + KroneckerDelta(0
      ,gO2)*(ZA(gI1,1)*((Lambda6 + Conj(Lambda6))*ZA(gI2,0) + 2*Lambda3*ZA(gI2,1))
      + ZA(gI1,0)*(4*Lambda1*ZA(gI2,0) + (Lambda6 + Conj(Lambda6))*ZA(gI2,1)))))
      - KroneckerDelta(1,gO1)*(2*Conj(Lambda6)*KroneckerDelta(0,gO2)*ZA(gI1,0)*ZA(
      gI2,0) + KroneckerDelta(0,gO2)*((Lambda4 + Lambda5)*ZA(gI1,0)*ZA(gI2,1) + ZA
      (gI1,1)*((Lambda4 + Lambda5)*ZA(gI2,0) + 2*Conj(Lambda7)*ZA(gI2,1))) +
      KroneckerDelta(1,gO2)*(ZA(gI1,1)*((Lambda7 + Conj(Lambda7))*ZA(gI2,0) + 4*
      Lambda2*ZA(gI2,1)) + ZA(gI1,0)*(2*Lambda3*ZA(gI2,0) + (Lambda7 + Conj(
      Lambda7))*ZA(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmconjHmHm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*(ZP(gI1,0)*(2*
      Lambda6*ZP(gI2,0) + (Lambda3 + Lambda4)*ZP(gI2,1)) + 2*ZP(gI1,1)*(Conj(
      Lambda5)*ZP(gI2,0) + Lambda7*ZP(gI2,1))) + KroneckerDelta(0,gO2)*(ZP(gI1,1)*
      (2*Lambda6*ZP(gI2,0) + (Lambda3 + Lambda4)*ZP(gI2,1)) + 2*ZP(gI1,0)*(2*
      Lambda1*ZP(gI2,0) + Conj(Lambda6)*ZP(gI2,1))))) - KroneckerDelta(1,gO1)*(2*
      Conj(Lambda6)*KroneckerDelta(0,gO2)*ZP(gI1,0)*ZP(gI2,0) + KroneckerDelta(1,
      gO2)*(2*ZP(gI1,1)*(Lambda7*ZP(gI2,0) + 2*Lambda2*ZP(gI2,1)) + ZP(gI1,0)*((
      Lambda3 + Lambda4)*ZP(gI2,0) + 2*Conj(Lambda7)*ZP(gI2,1))) + KroneckerDelta(
      0,gO2)*(2*Lambda5*ZP(gI1,0)*ZP(gI2,1) + ZP(gI1,1)*((Lambda3 + Lambda4)*ZP(
      gI2,0) + 2*Conj(Lambda7)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHmconjUHmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*(ZH(gI1,1)*((
      Lambda4 + Conj(Lambda5))*ZH(gI2,0) + 2*Lambda7*ZH(gI2,1)) + ZH(gI1,0)*(2*
      Lambda6*ZH(gI2,0) + (Lambda4 + Conj(Lambda5))*ZH(gI2,1))) + KroneckerDelta(0
      ,gO2)*(ZH(gI1,1)*((Lambda6 + Conj(Lambda6))*ZH(gI2,0) + 2*Lambda3*ZH(gI2,1))
      + ZH(gI1,0)*(4*Lambda1*ZH(gI2,0) + (Lambda6 + Conj(Lambda6))*ZH(gI2,1)))))
      - KroneckerDelta(1,gO1)*(2*Conj(Lambda6)*KroneckerDelta(0,gO2)*ZH(gI1,0)*ZH(
      gI2,0) + KroneckerDelta(0,gO2)*((Lambda4 + Lambda5)*ZH(gI1,0)*ZH(gI2,1) + ZH
      (gI1,1)*((Lambda4 + Lambda5)*ZH(gI2,0) + 2*Conj(Lambda7)*ZH(gI2,1))) +
      KroneckerDelta(1,gO2)*(ZH(gI1,1)*((Lambda7 + Conj(Lambda7))*ZH(gI2,0) + 4*
      Lambda2*ZH(gI2,1)) + ZH(gI1,0)*(2*Lambda3*ZH(gI2,0) + (Lambda7 + Conj(
      Lambda7))*ZH(gI2,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmHmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(v2*ZA(gI2,0) - v1*ZA(gI2,1))*(
      KroneckerDelta(0,gO2)*((Lambda6 - Conj(Lambda6))*ZP(gI1,0) + (Lambda4 -
      Lambda5)*ZP(gI1,1)) + KroneckerDelta(1,gO2)*((-Lambda4 + Conj(Lambda5))*ZP(
      gI1,0) + (Lambda7 - Conj(Lambda7))*ZP(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmHmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(KroneckerDelta(1,gO2)*(ZH(gI2,1)*((Lambda4*v1 + 2*Lambda7*v2
      + v1*Conj(Lambda5))*ZP(gI1,0) + (Lambda7*v1 + 4*Lambda2*v2 + v1*Conj(
      Lambda7))*ZP(gI1,1)) + ZH(gI2,0)*((2*Lambda6*v1 + Lambda4*v2 + v2*Conj(
      Lambda5))*ZP(gI1,0) + (2*Lambda3*v1 + Lambda7*v2 + v2*Conj(Lambda7))*ZP(gI1,
      1)))) - KroneckerDelta(0,gO2)*(ZH(gI2,0)*((4*Lambda1*v1 + Lambda6*v2 + v2*
      Conj(Lambda6))*ZP(gI1,0) + ((Lambda4 + Lambda5)*v2 + 2*v1*Conj(Lambda6))*ZP(
      gI1,1)) + ZH(gI2,1)*((Lambda6*v1 + 2*Lambda3*v2 + v1*Conj(Lambda6))*ZP(gI1,0
      ) + ((Lambda4 + Lambda5)*v1 + 2*v2*Conj(Lambda7))*ZP(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_48;
   std::complex<double> tmp_49;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_50;
      std::complex<double> tmp_51;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_51 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_50 += tmp_51;
      tmp_49 += (Vu(gI1,j2)) * tmp_50;
   }
   tmp_48 += tmp_49;
   result += (-KroneckerDelta(0,gO2)) * tmp_48;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_52;
   std::complex<double> tmp_53;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_54;
      std::complex<double> tmp_55;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_55 += Conj(Uu(gI1,j1))*Yu(j1,j2);
      }
      tmp_54 += tmp_55;
      tmp_53 += (Conj(Vd(gI2,j2))) * tmp_54;
   }
   tmp_52 += tmp_53;
   result += (-KroneckerDelta(1,gO1)) * tmp_52;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_56;
   std::complex<double> tmp_57;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_57 += Conj(Ye(j1,gI1))*Ue(gI2,j1);
   }
   tmp_56 += tmp_57;
   result += (-KroneckerDelta(0,gO2)) * tmp_56;

   return result;
}

double CLASSNAME::CpconjUHmbarFvFePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmChiChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*KroneckerDelta(1,gO2)*(1.4142135623730951*UP(gI2,1)*(g2up*ZN(
      gI1,0) + g2u*ZN(gI1,1)) + 2*g2u*UP(gI2,0)*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmChiChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(1.4142135623730951*Conj(UM(gI2,1))*(g1dp*Conj(ZN(gI1,0)) +
      g1d*Conj(ZN(gI1,1))) - 2*g1d*Conj(UM(gI2,0))*Conj(ZN(gI1,2)))*KroneckerDelta
      (0,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVWmAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZA(gI2,0) +
      KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVWmhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) + KroneckerDelta(1,gO2)*ZH
      (gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVPHm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 2) {
      result += -0.3872983346207417*g1*Cos(ThetaW())*ZP(gI2,gO2);
   }
   if (gI2 < 2) {
      result += -0.5*g2*Sin(ThetaW())*ZP(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUHmVZHm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 2) {
      result += -0.5*g2*Cos(ThetaW())*ZP(gI2,gO2);
   }
   if (gI2 < 2) {
      result += 0.3872983346207417*g1*Sin(ThetaW())*ZP(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpVGVGVG() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpVGbargGgG() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpVGbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double CLASSNAME::CpVGGluGluPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGGluGluPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG1() const
{
   double result = 0.0;

   result = -16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVGVGVGVG2() const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpVGVGVGVG3() const
{
   double result = 0.0;

   result = 16*Sqr(g3);

   return result;
}

double CLASSNAME::CpVPbargWmgWm() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbargWmCgWmC() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVPconjVWmVWm() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVPVPconjHmHm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + Cos(2*
      ThetaW())*(3*Sqr(g1) - 5*Sqr(g2)) + 5*Sqr(g2))*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1
      ,1)*ZP(gI2,1));

   return result;
}

double CLASSNAME::CpVPconjHmHm(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaW()) +
      g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM(gI2,0))*Sin(ThetaW())*UM(gI1,0) + Conj(UM(gI2,1))
      *(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVPbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP(gI1,0))*Sin(ThetaW())*UP(gI2,0) + Conj(UP(gI1,1))
      *(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(ThetaW()))*UP(gI2,1));

   return result;
}

double CLASSNAME::CpVPbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1
      *Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaW()) +
      g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

double CLASSNAME::CpVPbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpVPconjVWmHm(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*Cos(ThetaW())*(v1*ZP(gI2,0) + v2*ZP(gI2,1)
      );

   return result;
}

double CLASSNAME::CpVPVPconjVWmVWm1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPVPconjVWmVWm2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVPVPconjVWmVWm3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWmgWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbargWmCgWmC() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpVZconjVWmVWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZA(gI1,0)*ZA(gI2,0) + ZA(
      gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjHmHm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZP(gI1,0)*ZP(gI2,0) + ZP(
      gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZH(gI1,0)*ZH(gI2,0) + ZH(
      gI1,1)*ZH(gI2,1));

   return result;
}

double CLASSNAME::CpVZconjHmHm(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.1*KroneckerDelta(gI1,gI2)*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW()) + 3.872983346207417
      *g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) + ZA(gI2,1)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM(gI2,0))*Cos(ThetaW())*UM(gI1,0) + Conj(UM(gI2,1))
      *(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP(gI1,0))*Cos(ThetaW())*UP(gI2,0) + Conj(UP(gI1,1))
      *(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UP(gI2,1));

   return result;
}

double CLASSNAME::CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFeFePL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFeFePR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVZbarFvFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpVZChiChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(Conj
      (ZN(gI2,2))*ZN(gI1,2) - Conj(ZN(gI2,3))*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpVZChiChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(Conj(
      ZN(gI1,2))*ZN(gI2,2) - Conj(ZN(gI1,3))*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(v1
      *ZH(gI2,0) + v2*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZconjVWmHm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Sin(ThetaW())*(v1*ZP(gI2,0) + v2*ZP(gI2,1
      ));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmbargPgWm() const
{
   double result = 0.0;

   result = g2*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmbargWmCgP() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmbargWmCgZ() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmbargZgWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

double CLASSNAME::CpconjVWmVWmVP() const
{
   double result = 0.0;

   result = -(g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpconjVWmVZVWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW());

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjHmHm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(ZA(gI2,0)*ZP(gI1,0) + ZA(gI2,1)*ZP
      (gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(ZH(gI2,0)*ZP(gI1,0) + ZH(gI2,1)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_59 += Conj(Vd(gI2,j1))*Vu(gI1,j1);
   }
   tmp_58 += tmp_59;
   result += (-0.7071067811865475*g2) * tmp_58;

   return result;
}

double CLASSNAME::CpconjVWmbarFuFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFvFePL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(Ve(gI2,gI1));
   }

   return result;
}

double CLASSNAME::CpconjVWmbarFvFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,1) + 1.4142135623730951*Conj(UM(
      gI2,1))*ZN(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(ZN(gI1,1))*UP(gI2,0)) + 0.7071067811865475*g2*Conj(ZN(gI1
      ,3))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVPHm(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.3872983346207417*g1*g2*Cos(ThetaW())*(v1*ZP(gI2,0) + v2*ZP(gI2,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVWmhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(v1*ZH(gI2,0) + v2*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVZHm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Sin(ThetaW())*(v1*ZP(gI2,0) + v2*ZP(gI2,1
      ));

   return result;
}

double CLASSNAME::CpVWmconjVWmVPVP1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVPVP2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVPVP3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double CLASSNAME::CpVWmconjVWmconjVWmVWm1() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjVWmVWm2() const
{
   double result = 0.0;

   result = -Sqr(g2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjVWmVWm3() const
{
   double result = 0.0;

   result = 2*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_60;
      std::complex<double> tmp_61;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_61 += Conj(Vd(gI2,j2))*Yd(gO2,j2);
      }
      tmp_60 += tmp_61;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_60;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_62;
      std::complex<double> tmp_63;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_63 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_62 += tmp_63;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_62;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_64;
      std::complex<double> tmp_65;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_65 += Conj(Vu(gI2,j2))*Yd(gO2,j2);
      }
      tmp_64 += tmp_65;
      result += (-ZP(gI1,0)) * tmp_64;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_66;
      std::complex<double> tmp_67;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_67 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      tmp_66 += tmp_67;
      result += (-ZP(gI1,1)) * tmp_66;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_68;
      std::complex<double> tmp_69;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_69 += Conj(Vd(gI1,j2))*Yd(gO2,j2);
      }
      tmp_68 += tmp_69;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_68;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_70;
      std::complex<double> tmp_71;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_71 += Conj(Yd(j1,gO1))*Ud(gI1,j1);
      }
      tmp_70 += tmp_71;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_70;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Ud(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(Vd(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.2581988897471611*g1*Cos(ThetaW())*Ud(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(Vd(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFdVWmFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVWmFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vu(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.2581988897471611*g1*Sin(ThetaW())*Ud(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(Vd(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_72;
      std::complex<double> tmp_73;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_73 += Conj(Vd(gI2,j2))*Yu(gO2,j2);
      }
      tmp_72 += tmp_73;
      result += (-ZP(gI1,1)) * tmp_72;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_74;
      std::complex<double> tmp_75;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_75 += Conj(Yd(j1,gO1))*Ud(gI2,j1);
      }
      tmp_74 += tmp_75;
      result += (-ZP(gI1,0)) * tmp_74;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_76;
      std::complex<double> tmp_77;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_77 += Conj(Vu(gI2,j2))*Yu(gO2,j2);
      }
      tmp_76 += tmp_77;
      result += (0.7071067811865475*ZH(gI1,1)) * tmp_76;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_78;
      std::complex<double> tmp_79;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_79 += Conj(Yu(j1,gO1))*Uu(gI2,j1);
      }
      tmp_78 += tmp_79;
      result += (0.7071067811865475*ZH(gI1,1)) * tmp_78;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_80;
      std::complex<double> tmp_81;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_81 += Conj(Vu(gI1,j2))*Yu(gO2,j2);
      }
      tmp_80 += tmp_81;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) *
         tmp_80;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_82;
      std::complex<double> tmp_83;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_83 += Conj(Yu(j1,gO1))*Uu(gI1,j1);
      }
      tmp_82 += tmp_83;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_82;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Uu(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(Vu(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5163977794943222*g1*Cos(ThetaW())*Uu(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(Vu(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5163977794943222*g1*Sin(ThetaW())*Uu(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(Vu(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFuconjVWmFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjVWmFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Vd(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_84;
      std::complex<double> tmp_85;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_85 += Conj(Ve(gI2,j2))*Ye(gO2,j2);
      }
      tmp_84 += tmp_85;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_84;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_86;
      std::complex<double> tmp_87;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_87 += Conj(Ye(j1,gO1))*Ue(gI2,j1);
      }
      tmp_86 += tmp_87;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_86;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeHmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -(Ye(gO2,gI2)*ZP(gI1,0));
   }

   return result;
}

double CLASSNAME::CpbarUFeHmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_88;
      std::complex<double> tmp_89;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_89 += Conj(Ve(gI1,j2))*Ye(gO2,j2);
      }
      tmp_88 += tmp_89;
      result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) *
         tmp_88;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_90;
      std::complex<double> tmp_91;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_91 += Conj(Ye(j1,gO1))*Ue(gI1,j1);
      }
      tmp_90 += tmp_91;
      result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_90;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Cos(ThetaW())*Ue(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpbarUFeVWmFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarUFeVWmFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(gI2,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Sin(ThetaW())*Ue(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(Ve(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaHmPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(1.4142135623730951*Conj(UP(gI1,1))*(g2up*KroneckerDelta(0,gO2
      ) + g2u*KroneckerDelta(1,gO2)) + 2*g2u*Conj(UP(gI1,0))*KroneckerDelta(3,gO2)
      )*ZP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaHmPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g1d*KroneckerDelta(2,gO1)*UM(gI1,0) - 1.4142135623730951*(
      g1dp*KroneckerDelta(0,gO1) + g1d*KroneckerDelta(1,gO1))*UM(gI1,1))*ZP(gI2,0)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHmChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(1.4142135623730951*Conj(UM(gI2,1))*(g1dp*KroneckerDelta(0,gO2
      ) + g1d*KroneckerDelta(1,gO2)) - 2*g1d*Conj(UM(gI2,0))*KroneckerDelta(2,gO2)
      )*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHmChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(2*g2u*KroneckerDelta(3,gO1)*UP(gI2,0) + 1.4142135623730951*(
      g2up*KroneckerDelta(0,gO1) + g2u*KroneckerDelta(1,gO1))*UP(gI2,1))*ZP(gI1,1)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZN(gI2,2))*(g1dp*KroneckerDelta(0,gO2) - g1d*
      KroneckerDelta(1,gO2))*ZH(gI1,0) - g1d*Conj(ZN(gI2,1))*KroneckerDelta(2,gO2)
      *ZH(gI1,0) - g2up*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2)*ZH(gI1,1) + g2u*Conj
      (ZN(gI2,3))*KroneckerDelta(1,gO2)*ZH(gI1,1) + g2u*Conj(ZN(gI2,1))*
      KroneckerDelta(3,gO2)*ZH(gI1,1) + Conj(ZN(gI2,0))*(g1dp*KroneckerDelta(2,gO2
      )*ZH(gI1,0) - g2up*KroneckerDelta(3,gO2)*ZH(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(2,gO1)*ZH(gI1,0)*(g1dp*ZN(gI2,0) - g1d*ZN(gI2,1
      )) + KroneckerDelta(3,gO1)*ZH(gI1,1)*(-(g2up*ZN(gI2,0)) + g2u*ZN(gI2,1)) +
      g1dp*KroneckerDelta(0,gO1)*ZH(gI1,0)*ZN(gI2,2) - g1d*KroneckerDelta(1,gO1)*
      ZH(gI1,0)*ZN(gI2,2) - g2up*KroneckerDelta(0,gO1)*ZH(gI1,1)*ZN(gI2,3) + g2u*
      KroneckerDelta(1,gO1)*ZH(gI1,1)*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaVWmPR(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(1,gO2)*UM(gI1,0) + 0.7071067811865475*g2*
      KroneckerDelta(2,gO2)*UM(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChibarChaVWmPL(unsigned gO1, unsigned gI1) const
{
   std::complex<double> result;

   result = g2*Conj(UP(gI1,0))*KroneckerDelta(1,gO1) - 0.7071067811865475*g2*
      Conj(UP(gI1,1))*KroneckerDelta(3,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(ZN(gI1,2))*(g1dp*KroneckerDelta(0
      ,gO2) - g1d*KroneckerDelta(1,gO2))*ZA(gI2,0) - g1d*Conj(ZN(gI1,1))*
      KroneckerDelta(2,gO2)*ZA(gI2,0) + g2up*Conj(ZN(gI1,3))*KroneckerDelta(0,gO2)
      *ZA(gI2,1) - g2u*Conj(ZN(gI1,3))*KroneckerDelta(1,gO2)*ZA(gI2,1) - g2u*Conj(
      ZN(gI1,1))*KroneckerDelta(3,gO2)*ZA(gI2,1) + Conj(ZN(gI1,0))*(g1dp*
      KroneckerDelta(2,gO2)*ZA(gI2,0) + g2up*KroneckerDelta(3,gO2)*ZA(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(2,gO1)*ZA(gI2,0)*(-(
      g1dp*ZN(gI1,0)) + g1d*ZN(gI1,1)) + KroneckerDelta(3,gO1)*ZA(gI2,1)*(-(g2up*
      ZN(gI1,0)) + g2u*ZN(gI1,1)) - g1dp*KroneckerDelta(0,gO1)*ZA(gI2,0)*ZN(gI1,2)
      + g1d*KroneckerDelta(1,gO1)*ZA(gI2,0)*ZN(gI1,2) - g2up*KroneckerDelta(0,gO1
      )*ZA(gI2,1)*ZN(gI1,3) + g2u*KroneckerDelta(1,gO1)*ZA(gI2,1)*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)) + 0.7071067811865475*g2*
      KroneckerDelta(3,gO2)*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(1,gO1) +
      1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(2,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(
      KroneckerDelta(2,gO2)*ZN(gI2,2) - KroneckerDelta(3,gO2)*ZN(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(Conj(ZN(gI2,2))*KroneckerDelta(2,gO1) - Conj(ZN(gI2,3))*
      KroneckerDelta(3,gO1))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g1d*Conj(UM(gI1,1))*
      KroneckerDelta(0,gO2)*ZA(gI2,0) - g2u*Conj(UM(gI1,0))*KroneckerDelta(1,gO2)*
      ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g1d*KroneckerDelta(1,
      gO1)*UP(gI1,0)*ZA(gI2,0) - g2u*KroneckerDelta(0,gO1)*UP(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g1d*Conj(UM(gI2,1))*KroneckerDelta(0,gO2)*ZH(
      gI1,0) + g2u*Conj(UM(gI2,0))*KroneckerDelta(1,gO2)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g1d*KroneckerDelta(1,gO1)*UP(gI2,0)*ZH(gI1,0)
      + g2u*KroneckerDelta(0,gO1)*UP(gI2,1)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(2*g2u*Conj(ZN(gI2,3))*KroneckerDelta(0,gO2) +
      1.4142135623730951*(g2up*Conj(ZN(gI2,0)) + g2u*Conj(ZN(gI2,1)))*
      KroneckerDelta(1,gO2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(1.4142135623730951*KroneckerDelta(1,gO1)*(g1dp*ZN(gI2,0) +
      g1d*ZN(gI2,1)) - 2*g1d*KroneckerDelta(0,gO1)*ZN(gI2,2))*ZP(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVPChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*UP(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW(
      )))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVPChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM(gI2,0))*KroneckerDelta(0,gO1)*Sin(ThetaW()) + 0.1*Conj(
      UM(gI2,1))*KroneckerDelta(1,gO1)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*
      Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*UP(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      )))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*KroneckerDelta(0,gO1) + 0.1*Conj(
      UM(gI2,1))*KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZN(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(ZN(gI2,1))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1));

   return result;
}

double CLASSNAME::CpbarFvconjHmFePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvconjHmFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_92;
   std::complex<double> tmp_93;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_93 += Conj(Ye(j1,gO1))*Ue(gI2,j1);
   }
   tmp_92 += tmp_93;
   result += (-ZP(gI1,0)) * tmp_92;

   return result;
}

double CLASSNAME::CpbarFvVZFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpbarFvVZFvPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFvconjVWmFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFvconjVWmFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(Ve(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpGluVGGluPR() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3*AbsSqr(PhaseGlu);

   return result;
}

std::complex<double> CLASSNAME::CpGluVGGluPL() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3*AbsSqr(PhaseGlu);

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_94;
   std::complex<double> tmp_95;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_96;
      std::complex<double> tmp_97;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_97 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_96 += tmp_97;
      tmp_95 += (Conj(Vd(gI2,j2))) * tmp_96;
   }
   tmp_94 += tmp_95;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_94;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_98;
   std::complex<double> tmp_99;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_100;
      std::complex<double> tmp_101;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_101 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_100 += tmp_101;
      tmp_99 += (Vd(gO1,j2)) * tmp_100;
   }
   tmp_98 += tmp_99;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_98;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_102;
   std::complex<double> tmp_103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_104;
      std::complex<double> tmp_105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_105 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_104 += tmp_105;
      tmp_103 += (Conj(Vu(gI2,j2))) * tmp_104;
   }
   tmp_102 += tmp_103;
   result += (-ZP(gI1,0)) * tmp_102;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_106;
   std::complex<double> tmp_107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_108;
      std::complex<double> tmp_109;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_109 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_108 += tmp_109;
      tmp_107 += (Vd(gO1,j2)) * tmp_108;
   }
   tmp_106 += tmp_107;
   result += (-ZP(gI1,1)) * tmp_106;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_110;
   std::complex<double> tmp_111;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_113;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_113 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_112 += tmp_113;
      tmp_111 += (Conj(Vd(gI1,j2))) * tmp_112;
   }
   tmp_110 += tmp_111;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_110;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_114;
   std::complex<double> tmp_115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_116;
      std::complex<double> tmp_117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_117 += Conj(Yd(j1,j2))*Ud(gI1,j1);
      }
      tmp_116 += tmp_117;
      tmp_115 += (Vd(gO1,j2)) * tmp_116;
   }
   tmp_114 += tmp_115;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_114
      ;

   return result;
}

double CLASSNAME::CpbarFdVWmFuPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdVWmFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_118;
   std::complex<double> tmp_119;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_119 += Conj(Vu(gI2,j1))*Vd(gO1,j1);
   }
   tmp_118 += tmp_119;
   result += (-0.7071067811865475*g2) * tmp_118;

   return result;
}

double CLASSNAME::CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.2581988897471611*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI2,gO1)*(3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_120;
   std::complex<double> tmp_121;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_122;
      std::complex<double> tmp_123;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_123 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_122 += tmp_123;
      tmp_121 += (Conj(Ve(gI2,j2))) * tmp_122;
   }
   tmp_120 += tmp_121;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_120;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_124;
   std::complex<double> tmp_125;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_126;
      std::complex<double> tmp_127;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_127 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_126 += tmp_127;
      tmp_125 += (Ve(gO1,j2)) * tmp_126;
   }
   tmp_124 += tmp_125;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_124;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_128;
   std::complex<double> tmp_129;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_129 += Conj(Ue(gO2,j1))*Ye(j1,gI2);
   }
   tmp_128 += tmp_129;
   result += (-ZP(gI1,0)) * tmp_128;

   return result;
}

double CLASSNAME::CpbarFeHmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_130;
   std::complex<double> tmp_131;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_132;
      std::complex<double> tmp_133;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_133 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_132 += tmp_133;
      tmp_131 += (Conj(Ve(gI1,j2))) * tmp_132;
   }
   tmp_130 += tmp_131;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_130;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_134;
   std::complex<double> tmp_135;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_136;
      std::complex<double> tmp_137;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_137 += Conj(Ye(j1,j2))*Ue(gI1,j1);
      }
      tmp_136 += tmp_137;
      tmp_135 += (Ve(gO1,j2)) * tmp_136;
   }
   tmp_134 += tmp_135;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_134
      ;

   return result;
}

double CLASSNAME::CpbarFeVWmFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeVWmFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7071067811865475*g2*Ve(gO1,gI2);
   }

   return result;
}

double CLASSNAME::CpbarFeVZFePR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.7745966692414834*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFeVZFePL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW()) - 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_138;
   std::complex<double> tmp_139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_140;
      std::complex<double> tmp_141;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_141 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_140 += tmp_141;
      tmp_139 += (Conj(Vd(gI2,j2))) * tmp_140;
   }
   tmp_138 += tmp_139;
   result += (-ZP(gI1,1)) * tmp_138;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_142;
   std::complex<double> tmp_143;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_144;
      std::complex<double> tmp_145;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_145 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_144 += tmp_145;
      tmp_143 += (Vu(gO1,j2)) * tmp_144;
   }
   tmp_142 += tmp_143;
   result += (-ZP(gI1,0)) * tmp_142;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_146;
   std::complex<double> tmp_147;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_148;
      std::complex<double> tmp_149;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_149 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_148 += tmp_149;
      tmp_147 += (Conj(Vu(gI2,j2))) * tmp_148;
   }
   tmp_146 += tmp_147;
   result += (0.7071067811865475*ZH(gI1,1)) * tmp_146;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_150;
   std::complex<double> tmp_151;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_152;
      std::complex<double> tmp_153;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_153 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_152 += tmp_153;
      tmp_151 += (Vu(gO1,j2)) * tmp_152;
   }
   tmp_150 += tmp_151;
   result += (0.7071067811865475*ZH(gI1,1)) * tmp_150;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_154;
   std::complex<double> tmp_155;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_156;
      std::complex<double> tmp_157;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_157 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_156 += tmp_157;
      tmp_155 += (Conj(Vu(gI1,j2))) * tmp_156;
   }
   tmp_154 += tmp_155;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) * tmp_154;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_158;
   std::complex<double> tmp_159;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_160;
      std::complex<double> tmp_161;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_161 += Conj(Yu(j1,j2))*Uu(gI1,j1);
      }
      tmp_160 += tmp_161;
      tmp_159 += (Vu(gO1,j2)) * tmp_160;
   }
   tmp_158 += tmp_159;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) * tmp_158
      ;

   return result;
}

double CLASSNAME::CpbarFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(gI2,gO2);

   return result;
}

double CLASSNAME::CpbarFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(0.7745966692414834*g1
      *Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5163977794943222*g1*KroneckerDelta(gI2,gO2)*Sin(ThetaW());

   return result;
}

double CLASSNAME::CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.03333333333333333*KroneckerDelta(gI2,gO1)*(-15*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double CLASSNAME::CpbarFuconjVWmFdPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjVWmFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_162;
   std::complex<double> tmp_163;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_163 += Conj(Vd(gI2,j1))*Vu(gO1,j1);
   }
   tmp_162 += tmp_163;
   result += (-0.7071067811865475*g2) * tmp_162;

   return result;
}


std::complex<double> CLASSNAME::self_energy_hh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmCgWmC(gO1)*CpUhhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmgWm(gO1)*CpUhhbargWmgWm(gO2));
   result += -(B0(p,MVZ,MVZ)*CpUhhbargZgZ(gO1)*CpUhhbargZgZ(gO2));
   result += 4*(-0.5 + B0(p,MVWm,MVWm))*Conj(CpUhhconjVWmVWm(gO2))*
      CpUhhconjVWmVWm(gO1);
   result += 2*(-0.5 + B0(p,MVZ,MVZ))*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
   std::complex<double> tmp_164;
   std::complex<double> tmp_165;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_165 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_164 += tmp_165;
   result += (-0.5) * tmp_164;
   std::complex<double> tmp_166;
   std::complex<double> tmp_167;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_167 += A0(MHm(gI1))*CpUhhUhhconjHmHm(gO1,gO2,gI1,gI1);
   }
   tmp_166 += tmp_167;
   result += (-1) * tmp_166;
   std::complex<double> tmp_168;
   std::complex<double> tmp_169;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_169 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_168 += tmp_169;
   result += (-0.5) * tmp_168;
   std::complex<double> tmp_170;
   std::complex<double> tmp_171;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_172;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_172 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))*
            CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_171 += tmp_172;
   }
   tmp_170 += tmp_171;
   result += (0.5) * tmp_170;
   std::complex<double> tmp_173;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_174;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_174 += B0(p,MHm(gI1),MHm(gI2))*Conj(CpUhhconjHmHm(gO2,gI1,
            gI2))*CpUhhconjHmHm(gO1,gI1,gI2);
      }
      tmp_173 += tmp_174;
   }
   result += tmp_173;
   std::complex<double> tmp_175;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_176;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_176 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))*
            CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_175 += tmp_176;
   }
   result += tmp_175;
   std::complex<double> tmp_177;
   std::complex<double> tmp_178;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_179;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_179 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))*
            CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_178 += tmp_179;
   }
   tmp_177 += tmp_178;
   result += (0.5) * tmp_177;
   std::complex<double> tmp_180;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_181;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_181 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*CpUhhbarChaChaPL
            (gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*CpUhhbarChaChaPR(
            gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_180 += tmp_181;
   }
   result += tmp_180;
   std::complex<double> tmp_182;
   std::complex<double> tmp_183;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_184;
      std::complex<double> tmp_185;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_185 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2,
            gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_184 += tmp_185;
      tmp_183 += (MCha(gI1)) * tmp_184;
   }
   tmp_182 += tmp_183;
   result += (-2) * tmp_182;
   std::complex<double> tmp_186;
   std::complex<double> tmp_187;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_188;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_188 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(gO1
            ,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,gI1,
            gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_187 += tmp_188;
   }
   tmp_186 += tmp_187;
   result += (3) * tmp_186;
   std::complex<double> tmp_189;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_190;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_190 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(gO1
            ,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,gI1,
            gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_189 += tmp_190;
   }
   result += tmp_189;
   std::complex<double> tmp_191;
   std::complex<double> tmp_192;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_193;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_193 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(gO1
            ,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,gI1,
            gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_192 += tmp_193;
   }
   tmp_191 += tmp_192;
   result += (3) * tmp_191;
   std::complex<double> tmp_194;
   std::complex<double> tmp_195;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_196;
      std::complex<double> tmp_197;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_197 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1,
            gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*
            CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_196 += tmp_197;
      tmp_195 += (MFd(gI1)) * tmp_196;
   }
   tmp_194 += tmp_195;
   result += (-6) * tmp_194;
   std::complex<double> tmp_198;
   std::complex<double> tmp_199;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_200;
      std::complex<double> tmp_201;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_201 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1,
            gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*
            CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_200 += tmp_201;
      tmp_199 += (MFe(gI1)) * tmp_200;
   }
   tmp_198 += tmp_199;
   result += (-2) * tmp_198;
   std::complex<double> tmp_202;
   std::complex<double> tmp_203;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_204;
      std::complex<double> tmp_205;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_205 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1,
            gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*
            CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_204 += tmp_205;
      tmp_203 += (MFu(gI1)) * tmp_204;
   }
   tmp_202 += tmp_203;
   result += (-6) * tmp_202;
   std::complex<double> tmp_206;
   std::complex<double> tmp_207;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_208;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_208 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_207 += tmp_208;
   }
   tmp_206 += tmp_207;
   result += (0.5) * tmp_206;
   std::complex<double> tmp_209;
   std::complex<double> tmp_210;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_211;
      std::complex<double> tmp_212;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_212 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,gI1
            ,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2))*
            CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_211 += tmp_212;
      tmp_210 += (MChi(gI1)) * tmp_211;
   }
   tmp_209 += tmp_210;
   result += (-1) * tmp_209;
   std::complex<double> tmp_213;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_213 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_213;
   std::complex<double> tmp_214;
   std::complex<double> tmp_215;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_215 += Conj(CpUhhconjVWmHm(gO2,gI2))*CpUhhconjVWmHm(gO1,gI2)*F0(p,
         MHm(gI2),MVWm);
   }
   tmp_214 += tmp_215;
   result += (2) * tmp_214;
   result += 4*CpUhhUhhconjVWmVWm(gO1,gO2)*(A0(MVWm) - 0.5*Sqr(MVWm));
   result += 2*CpUhhUhhVZVZ(gO1,gO2)*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   std::complex<double> tmp_216;
   std::complex<double> tmp_217;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_217 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_216 += tmp_217;
   result += (-0.5) * tmp_216;
   std::complex<double> tmp_218;
   std::complex<double> tmp_219;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_219 += A0(MHm(gI1))*CpUAhUAhconjHmHm(gO1,gO2,gI1,gI1);
   }
   tmp_218 += tmp_219;
   result += (-1) * tmp_218;
   std::complex<double> tmp_220;
   std::complex<double> tmp_221;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_221 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_220 += tmp_221;
   result += (-0.5) * tmp_220;
   std::complex<double> tmp_222;
   std::complex<double> tmp_223;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_224;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_224 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))*
            CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_223 += tmp_224;
   }
   tmp_222 += tmp_223;
   result += (0.5) * tmp_222;
   std::complex<double> tmp_225;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_226;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_226 += B0(p,MHm(gI1),MHm(gI2))*Conj(CpUAhconjHmHm(gO2,gI1,
            gI2))*CpUAhconjHmHm(gO1,gI1,gI2);
      }
      tmp_225 += tmp_226;
   }
   result += tmp_225;
   std::complex<double> tmp_227;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_228;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_228 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))*
            CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_227 += tmp_228;
   }
   result += tmp_227;
   std::complex<double> tmp_229;
   std::complex<double> tmp_230;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_231;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_231 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))*
            CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_230 += tmp_231;
   }
   tmp_229 += tmp_230;
   result += (0.5) * tmp_229;
   std::complex<double> tmp_232;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_233;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_233 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*CpUAhbarChaChaPL
            (gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*CpUAhbarChaChaPR(
            gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_232 += tmp_233;
   }
   result += tmp_232;
   std::complex<double> tmp_234;
   std::complex<double> tmp_235;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_236;
      std::complex<double> tmp_237;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_237 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2,
            gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_236 += tmp_237;
      tmp_235 += (MCha(gI1)) * tmp_236;
   }
   tmp_234 += tmp_235;
   result += (-2) * tmp_234;
   std::complex<double> tmp_238;
   std::complex<double> tmp_239;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_240;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_240 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(gO1
            ,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,gI1,
            gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_239 += tmp_240;
   }
   tmp_238 += tmp_239;
   result += (3) * tmp_238;
   std::complex<double> tmp_241;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_242;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_242 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(gO1
            ,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,gI1,
            gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_241 += tmp_242;
   }
   result += tmp_241;
   std::complex<double> tmp_243;
   std::complex<double> tmp_244;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_245;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_245 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(gO1
            ,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,gI1,
            gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_244 += tmp_245;
   }
   tmp_243 += tmp_244;
   result += (3) * tmp_243;
   std::complex<double> tmp_246;
   std::complex<double> tmp_247;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_248;
      std::complex<double> tmp_249;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_249 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1,
            gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*
            CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_248 += tmp_249;
      tmp_247 += (MFd(gI1)) * tmp_248;
   }
   tmp_246 += tmp_247;
   result += (-6) * tmp_246;
   std::complex<double> tmp_250;
   std::complex<double> tmp_251;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_252;
      std::complex<double> tmp_253;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_253 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1,
            gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*
            CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_252 += tmp_253;
      tmp_251 += (MFe(gI1)) * tmp_252;
   }
   tmp_250 += tmp_251;
   result += (-2) * tmp_250;
   std::complex<double> tmp_254;
   std::complex<double> tmp_255;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_256;
      std::complex<double> tmp_257;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_257 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1,
            gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*
            CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_256 += tmp_257;
      tmp_255 += (MFu(gI1)) * tmp_256;
   }
   tmp_254 += tmp_255;
   result += (-6) * tmp_254;
   std::complex<double> tmp_258;
   std::complex<double> tmp_259;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_260;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_260 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_259 += tmp_260;
   }
   tmp_258 += tmp_259;
   result += (0.5) * tmp_258;
   std::complex<double> tmp_261;
   std::complex<double> tmp_262;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_263;
      std::complex<double> tmp_264;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_264 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,gI1
            ,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2))*
            CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_263 += tmp_264;
      tmp_262 += (MChi(gI1)) * tmp_263;
   }
   tmp_261 += tmp_262;
   result += (-1) * tmp_261;
   std::complex<double> tmp_265;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_265 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_265;
   std::complex<double> tmp_266;
   std::complex<double> tmp_267;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_267 += Conj(CpUAhconjVWmHm(gO2,gI2))*CpUAhconjVWmHm(gO1,gI2)*F0(p,
         MHm(gI2),MVWm);
   }
   tmp_266 += tmp_267;
   result += (2) * tmp_266;
   result += 4*CpUAhUAhconjVWmVWm(gO1,gO2)*(A0(MVWm) - 0.5*Sqr(MVWm));
   result += 2*CpUAhUAhVZVZ(gO1,gO2)*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Hm(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*(-0.5 + B0(p,0,MVWm))*Conj(CpconjUHmVWmVP(gO2))*CpconjUHmVWmVP(
      gO1);
   result += 4*(-0.5 + B0(p,MVWm,MVZ))*Conj(CpconjUHmVZVWm(gO2))*CpconjUHmVZVWm
      (gO1);
   result += -(B0(p,MVZ,MVWm)*CpconjUHmbargWmCgZ(gO1)*CpUHmgWmCbargZ(gO2));
   result += -(B0(p,MVWm,MVZ)*CpconjUHmbargZgWm(gO1)*CpUHmgZbargWm(gO2));
   std::complex<double> tmp_268;
   std::complex<double> tmp_269;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_269 += A0(MAh(gI1))*CpUHmconjUHmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_268 += tmp_269;
   result += (-0.5) * tmp_268;
   std::complex<double> tmp_270;
   std::complex<double> tmp_271;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_271 += A0(MHm(gI1))*CpUHmconjUHmconjHmHm(gO1,gO2,gI1,gI1);
   }
   tmp_270 += tmp_271;
   result += (-1) * tmp_270;
   std::complex<double> tmp_272;
   std::complex<double> tmp_273;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_273 += A0(Mhh(gI1))*CpUHmconjUHmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_272 += tmp_273;
   result += (-0.5) * tmp_272;
   std::complex<double> tmp_274;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_275;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_275 += B0(p,MHm(gI1),MAh(gI2))*Conj(CpconjUHmHmAh(gO2,gI1,
            gI2))*CpconjUHmHmAh(gO1,gI1,gI2);
      }
      tmp_274 += tmp_275;
   }
   result += tmp_274;
   std::complex<double> tmp_276;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_277;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_277 += B0(p,MHm(gI1),Mhh(gI2))*Conj(CpconjUHmHmhh(gO2,gI1,
            gI2))*CpconjUHmHmhh(gO1,gI1,gI2);
      }
      tmp_276 += tmp_277;
   }
   result += tmp_276;
   std::complex<double> tmp_278;
   std::complex<double> tmp_279;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_280;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_280 += (Conj(CpconjUHmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHmbarFuFdPR(gO2,gI1,gI2)
            )*CpconjUHmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_279 += tmp_280;
   }
   tmp_278 += tmp_279;
   result += (3) * tmp_278;
   std::complex<double> tmp_281;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_282;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_282 += (Conj(CpconjUHmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHmbarFvFePR(gO2,gI1,gI2)
            )*CpconjUHmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_281 += tmp_282;
   }
   result += tmp_281;
   std::complex<double> tmp_283;
   std::complex<double> tmp_284;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_285;
      std::complex<double> tmp_286;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_286 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHmbarFuFdPR(gO2,
            gI1,gI2))*CpconjUHmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHmbarFuFdPL(
            gO2,gI1,gI2))*CpconjUHmbarFuFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_285 += tmp_286;
      tmp_284 += (MFu(gI1)) * tmp_285;
   }
   tmp_283 += tmp_284;
   result += (-6) * tmp_283;
   std::complex<double> tmp_287;
   std::complex<double> tmp_288;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_289;
      std::complex<double> tmp_290;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_290 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHmbarFvFePR(gO2,
            gI1,gI2))*CpconjUHmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHmbarFvFePL(
            gO2,gI1,gI2))*CpconjUHmbarFvFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_289 += tmp_290;
      tmp_288 += (MFv(gI1)) * tmp_289;
   }
   tmp_287 += tmp_288;
   result += (-2) * tmp_287;
   std::complex<double> tmp_291;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_292;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_292 += (Conj(CpconjUHmChiChaPL(gO2,gI1,gI2))*
            CpconjUHmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHmChiChaPR(gO2,gI1,gI2))*
            CpconjUHmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_291 += tmp_292;
   }
   result += tmp_291;
   std::complex<double> tmp_293;
   std::complex<double> tmp_294;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_295;
      std::complex<double> tmp_296;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_296 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHmChiChaPR(gO2
            ,gI1,gI2))*CpconjUHmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHmChiChaPL(gO2
            ,gI1,gI2))*CpconjUHmChiChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_295 += tmp_296;
      tmp_294 += (MChi(gI1)) * tmp_295;
   }
   tmp_293 += tmp_294;
   result += (-2) * tmp_293;
   std::complex<double> tmp_297;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_297 += Conj(CpconjUHmVWmAh(gO2,gI2))*CpconjUHmVWmAh(gO1,gI2)*F0(p,
         MAh(gI2),MVWm);
   }
   result += tmp_297;
   std::complex<double> tmp_298;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_298 += Conj(CpconjUHmVWmhh(gO2,gI2))*CpconjUHmVWmhh(gO1,gI2)*F0(p,
         Mhh(gI2),MVWm);
   }
   result += tmp_298;
   std::complex<double> tmp_299;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_299 += Conj(CpconjUHmVPHm(gO2,gI2))*CpconjUHmVPHm(gO1,gI2)*F0(p,
         MHm(gI2),0);
   }
   result += tmp_299;
   std::complex<double> tmp_300;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_300 += Conj(CpconjUHmVZHm(gO2,gI2))*CpconjUHmVZHm(gO1,gI2)*F0(p,
         MHm(gI2),MVZ);
   }
   result += tmp_300;
   result += 4*CpUHmconjUHmconjVWmVWm(gO1,gO2)*(A0(MVWm) - 0.5*Sqr(MVWm));
   result += 2*CpUHmconjUHmVZVZ(gO1,gO2)*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VG(double p ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpVGbargGgG())*B00(p,MVG,MVG);
   result += -1.5*AbsSqr(CpVGVGVG())*(10*B00(p,0,0) + 0.6666666666666666*Sqr(p)
      + 4*B0(p,0,0)*Sqr(p));
   result += 0;
   std::complex<double> tmp_301;
   std::complex<double> tmp_302;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_303;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_303 += (AbsSqr(CpVGbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVGbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_303 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVGbarFdFdPL(gI1,gI2))*CpVGbarFdFdPR(gI1,gI2));
      }
      tmp_302 += tmp_303;
   }
   tmp_301 += tmp_302;
   result += (0.5) * tmp_301;
   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_306;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_306 += (AbsSqr(CpVGbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVGbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_306 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVGbarFuFuPL(gI1,gI2))*CpVGbarFuFuPR(gI1,gI2));
      }
      tmp_305 += tmp_306;
   }
   tmp_304 += tmp_305;
   result += (0.5) * tmp_304;
   result += 1.5*((AbsSqr(CpVGGluGluPL(0,0)) + AbsSqr(CpVGGluGluPR(0,0)))*H0(p,
      MGlu,MGlu) + 4*B0(p,MGlu,MGlu)*Re(Conj(CpVGGluGluPL(0,0))*CpVGGluGluPR(0,0))
      *Sqr(MGlu));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VP(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVPbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVPbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVPVPconjVWmVWm1() + CpVPVPconjVWmVWm2() +
      CpVPVPconjVWmVWm3()));
   std::complex<double> tmp_307;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_307 += A0(MHm(gI1))*CpVPVPconjHmHm(gI1,gI1);
   }
   result += tmp_307;
   std::complex<double> tmp_308;
   std::complex<double> tmp_309;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_310;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_310 += AbsSqr(CpVPconjHmHm(gI1,gI2))*B00(p,MHm(gI1),MHm(gI2)
            );
      }
      tmp_309 += tmp_310;
   }
   tmp_308 += tmp_309;
   result += (-4) * tmp_308;
   std::complex<double> tmp_311;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_312;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_312 += (AbsSqr(CpVPbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVPbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_312 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVPbarChaChaPL(gI1,gI2))*CpVPbarChaChaPR(gI1,gI2));
      }
      tmp_311 += tmp_312;
   }
   result += tmp_311;
   std::complex<double> tmp_313;
   std::complex<double> tmp_314;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_315;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_315 += (AbsSqr(CpVPbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVPbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_315 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVPbarFdFdPL(gI1,gI2))*CpVPbarFdFdPR(gI1,gI2));
      }
      tmp_314 += tmp_315;
   }
   tmp_313 += tmp_314;
   result += (3) * tmp_313;
   std::complex<double> tmp_316;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_317;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_317 += (AbsSqr(CpVPbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVPbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_317 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVPbarFeFePL(gI1,gI2))*CpVPbarFeFePR(gI1,gI2));
      }
      tmp_316 += tmp_317;
   }
   result += tmp_316;
   std::complex<double> tmp_318;
   std::complex<double> tmp_319;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_320;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_320 += (AbsSqr(CpVPbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVPbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_320 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVPbarFuFuPL(gI1,gI2))*CpVPbarFuFuPR(gI1,gI2));
      }
      tmp_319 += tmp_320;
   }
   tmp_318 += tmp_319;
   result += (3) * tmp_318;
   std::complex<double> tmp_321;
   std::complex<double> tmp_322;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_322 += AbsSqr(CpVPconjVWmHm(gI2))*B0(p,MVWm,MHm(gI2));
   }
   tmp_321 += tmp_322;
   result += (2) * tmp_321;
   result += 2*CpVPVPconjVWmVWm1()*Sqr(MVWm);
   result += -(AbsSqr(CpVPconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) - 2*(
      2*Sqr(MVWm) - 0.3333333333333333*Sqr(p)) + B0(p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*
      Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_323;
   std::complex<double> tmp_324;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_324 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_323 += tmp_324;
   result += (0.5) * tmp_323;
   std::complex<double> tmp_325;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_325 += A0(MHm(gI1))*CpVZVZconjHmHm(gI1,gI1);
   }
   result += tmp_325;
   std::complex<double> tmp_326;
   std::complex<double> tmp_327;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_327 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_326 += tmp_327;
   result += (0.5) * tmp_326;
   std::complex<double> tmp_328;
   std::complex<double> tmp_329;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_330;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_330 += AbsSqr(CpVZconjHmHm(gI1,gI2))*B00(p,MHm(gI1),MHm(gI2)
            );
      }
      tmp_329 += tmp_330;
   }
   tmp_328 += tmp_329;
   result += (-4) * tmp_328;
   std::complex<double> tmp_331;
   std::complex<double> tmp_332;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_333;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_333 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_332 += tmp_333;
   }
   tmp_331 += tmp_332;
   result += (-4) * tmp_331;
   std::complex<double> tmp_334;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_335;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_335 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_335 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_334 += tmp_335;
   }
   result += tmp_334;
   std::complex<double> tmp_336;
   std::complex<double> tmp_337;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_338;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_338 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_338 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_337 += tmp_338;
   }
   tmp_336 += tmp_337;
   result += (3) * tmp_336;
   std::complex<double> tmp_339;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_340;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_340 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_340 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_339 += tmp_340;
   }
   result += tmp_339;
   std::complex<double> tmp_341;
   std::complex<double> tmp_342;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_343;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_343 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_343 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_342 += tmp_343;
   }
   tmp_341 += tmp_342;
   result += (3) * tmp_341;
   std::complex<double> tmp_344;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_345;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_345 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_345 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_344 += tmp_345;
   }
   result += tmp_344;
   std::complex<double> tmp_346;
   std::complex<double> tmp_347;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_348;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_348 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR(
            gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_348 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_347 += tmp_348;
   }
   tmp_346 += tmp_347;
   result += (0.5) * tmp_346;
   std::complex<double> tmp_349;
   std::complex<double> tmp_350;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_350 += AbsSqr(CpVZconjVWmHm(gI2))*B0(p,MVWm,MHm(gI2));
   }
   tmp_349 += tmp_350;
   result += (2) * tmp_349;
   std::complex<double> tmp_351;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_351 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_351;
   result += 2*CpVZVZconjVWmVWm1()*Sqr(MVWm);
   result += -(AbsSqr(CpVZconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) - 2*(
      2*Sqr(MVWm) - 0.3333333333333333*Sqr(p)) + B0(p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*
      Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWmbargPgWm())*B00(p,MVWm,MVP);
   result += AbsSqr(CpconjVWmbargWmCgP())*B00(p,MVP,MVWm);
   result += AbsSqr(CpconjVWmbargWmCgZ())*B00(p,MVZ,MVWm);
   result += AbsSqr(CpconjVWmbargZgWm())*B00(p,MVWm,MVZ);
   result += -(A0(MVWm)*(4*CpVWmconjVWmconjVWmVWm1() + CpVWmconjVWmconjVWmVWm2(
      ) + CpVWmconjVWmconjVWmVWm3()));
   result += 0;
   std::complex<double> tmp_352;
   std::complex<double> tmp_353;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_353 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_352 += tmp_353;
   result += (0.5) * tmp_352;
   std::complex<double> tmp_354;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_354 += A0(MHm(gI1))*CpVWmconjVWmconjHmHm(gI1,gI1);
   }
   result += tmp_354;
   std::complex<double> tmp_355;
   std::complex<double> tmp_356;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_356 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_355 += tmp_356;
   result += (0.5) * tmp_355;
   std::complex<double> tmp_357;
   std::complex<double> tmp_358;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_359;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_359 += AbsSqr(CpconjVWmHmAh(gI1,gI2))*B00(p,MAh(gI2),MHm(gI1
            ));
      }
      tmp_358 += tmp_359;
   }
   tmp_357 += tmp_358;
   result += (-4) * tmp_357;
   std::complex<double> tmp_360;
   std::complex<double> tmp_361;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_362;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_362 += AbsSqr(CpconjVWmHmhh(gI1,gI2))*B00(p,Mhh(gI2),MHm(gI1
            ));
      }
      tmp_361 += tmp_362;
   }
   tmp_360 += tmp_361;
   result += (-4) * tmp_360;
   std::complex<double> tmp_363;
   std::complex<double> tmp_364;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_365;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_365 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_365 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_364 += tmp_365;
   }
   tmp_363 += tmp_364;
   result += (3) * tmp_363;
   std::complex<double> tmp_366;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_367;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_367 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_367 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_366 += tmp_367;
   }
   result += tmp_366;
   std::complex<double> tmp_368;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_369;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_369 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_369 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_368 += tmp_369;
   }
   result += tmp_368;
   std::complex<double> tmp_370;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_370 += AbsSqr(CpconjVWmVPHm(gI2))*B0(p,0,MHm(gI2));
   }
   result += tmp_370;
   std::complex<double> tmp_371;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_371 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_371;
   std::complex<double> tmp_372;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_372 += AbsSqr(CpconjVWmVZHm(gI2))*B0(p,MVZ,MHm(gI2));
   }
   result += tmp_372;
   result += 2*CpVWmconjVWmconjVWmVWm1()*Sqr(MVWm);
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) - 2*(Sqr(
      MVWm) - 0.3333333333333333*Sqr(p)) + B0(p,MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += 0.5*(-(A0(MVZ)*(4*CpVWmconjVWmVZVZ1() + CpVWmconjVWmVZVZ2() +
      CpVWmconjVWmVZVZ3())) + 2*CpVWmconjVWmVZVZ1()*Sqr(MVZ));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) - 2*(Sqr(MVWm) + Sqr(MVZ) - 0.3333333333333333*Sqr(p)) + B0(p,MVZ,MVWm)*(
      Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_373;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_374;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_374 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_373 += tmp_374;
   }
   result += tmp_373;
   std::complex<double> tmp_375;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_376;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_376 += B0(p,MFu(gI2),MHm(gI1))*Conj(CpbarUFdHmFuPL(gO2,gI1,
            gI2))*CpbarUFdHmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_375 += tmp_376;
   }
   result += tmp_375;
   std::complex<double> tmp_377;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_378;
      std::complex<double> tmp_379;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_379 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_378 += tmp_379;
      tmp_377 += (MFd(gI1)) * tmp_378;
   }
   result += tmp_377;
   std::complex<double> tmp_380;
   std::complex<double> tmp_381;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_381 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_380 += tmp_381;
   result += (-5.333333333333333) * tmp_380;
   std::complex<double> tmp_382;
   std::complex<double> tmp_383;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_383 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_382 += tmp_383;
   result += (-4) * tmp_382;
   std::complex<double> tmp_384;
   std::complex<double> tmp_385;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_385 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_384 += tmp_385;
   result += (-4) * tmp_384;
   std::complex<double> tmp_386;
   std::complex<double> tmp_387;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_387 += (-0.5 + B0(p,MFu(gI2),MVWm))*Conj(CpbarUFdVWmFuPR(gO2,gI2))
         *CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_386 += tmp_387;
   result += (-4) * tmp_386;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_388;
   std::complex<double> tmp_389;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_390;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_390 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_389 += tmp_390;
   }
   tmp_388 += tmp_389;
   result += (-0.5) * tmp_388;
   std::complex<double> tmp_391;
   std::complex<double> tmp_392;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_393;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_393 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarUFdHmFuPR(gO2,gI1,
            gI2))*CpbarUFdHmFuPR(gO1,gI1,gI2);
      }
      tmp_392 += tmp_393;
   }
   tmp_391 += tmp_392;
   result += (-0.5) * tmp_391;
   std::complex<double> tmp_394;
   std::complex<double> tmp_395;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_396;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_396 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_395 += tmp_396;
   }
   tmp_394 += tmp_395;
   result += (-0.5) * tmp_394;
   std::complex<double> tmp_397;
   std::complex<double> tmp_398;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_398 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_397 += tmp_398;
   result += (-1.3333333333333333) * tmp_397;
   std::complex<double> tmp_399;
   std::complex<double> tmp_400;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_400 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_399 += tmp_400;
   result += (-1) * tmp_399;
   std::complex<double> tmp_401;
   std::complex<double> tmp_402;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_402 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_401 += tmp_402;
   result += (-1) * tmp_401;
   std::complex<double> tmp_403;
   std::complex<double> tmp_404;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_404 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_403 += tmp_404;
   result += (-1) * tmp_403;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_405;
   std::complex<double> tmp_406;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_407;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_407 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_406 += tmp_407;
   }
   tmp_405 += tmp_406;
   result += (-0.5) * tmp_405;
   std::complex<double> tmp_408;
   std::complex<double> tmp_409;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_410;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_410 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarUFdHmFuPL(gO2,gI1,
            gI2))*CpbarUFdHmFuPL(gO1,gI1,gI2);
      }
      tmp_409 += tmp_410;
   }
   tmp_408 += tmp_409;
   result += (-0.5) * tmp_408;
   std::complex<double> tmp_411;
   std::complex<double> tmp_412;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_413;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_413 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_412 += tmp_413;
   }
   tmp_411 += tmp_412;
   result += (-0.5) * tmp_411;
   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_415 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_414 += tmp_415;
   result += (-1.3333333333333333) * tmp_414;
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_417 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_416 += tmp_417;
   result += (-1) * tmp_416;
   std::complex<double> tmp_418;
   std::complex<double> tmp_419;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_419 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_418 += tmp_419;
   result += (-1) * tmp_418;
   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_421 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_420 += tmp_421;
   result += (-1) * tmp_420;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_422;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_423;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_423 += B0(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_422 += tmp_423;
   }
   result += tmp_422;
   std::complex<double> tmp_424;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_425;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_425 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_424 += tmp_425;
   }
   result += tmp_424;
   std::complex<double> tmp_426;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_427;
      std::complex<double> tmp_428;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_428 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_427 += tmp_428;
      tmp_426 += (MFu(gI1)) * tmp_427;
   }
   result += tmp_426;
   std::complex<double> tmp_429;
   std::complex<double> tmp_430;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_430 += (-0.5 + B0(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_429 += tmp_430;
   result += (-4) * tmp_429;
   std::complex<double> tmp_431;
   std::complex<double> tmp_432;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_432 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_431 += tmp_432;
   result += (-5.333333333333333) * tmp_431;
   std::complex<double> tmp_433;
   std::complex<double> tmp_434;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_434 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_433 += tmp_434;
   result += (-4) * tmp_433;
   std::complex<double> tmp_435;
   std::complex<double> tmp_436;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_436 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_435 += tmp_436;
   result += (-4) * tmp_435;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_437;
   std::complex<double> tmp_438;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_439;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_439 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPR(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2);
      }
      tmp_438 += tmp_439;
   }
   tmp_437 += tmp_438;
   result += (-0.5) * tmp_437;
   std::complex<double> tmp_440;
   std::complex<double> tmp_441;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_442;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_442 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_441 += tmp_442;
   }
   tmp_440 += tmp_441;
   result += (-0.5) * tmp_440;
   std::complex<double> tmp_443;
   std::complex<double> tmp_444;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_445;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_445 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_444 += tmp_445;
   }
   tmp_443 += tmp_444;
   result += (-0.5) * tmp_443;
   std::complex<double> tmp_446;
   std::complex<double> tmp_447;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_447 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPL(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_446 += tmp_447;
   result += (-1) * tmp_446;
   std::complex<double> tmp_448;
   std::complex<double> tmp_449;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_449 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_448 += tmp_449;
   result += (-1.3333333333333333) * tmp_448;
   std::complex<double> tmp_450;
   std::complex<double> tmp_451;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_451 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_450 += tmp_451;
   result += (-1) * tmp_450;
   std::complex<double> tmp_452;
   std::complex<double> tmp_453;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_453 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_452 += tmp_453;
   result += (-1) * tmp_452;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_454;
   std::complex<double> tmp_455;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_456;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_456 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPL(gO1,gI1,gI2);
      }
      tmp_455 += tmp_456;
   }
   tmp_454 += tmp_455;
   result += (-0.5) * tmp_454;
   std::complex<double> tmp_457;
   std::complex<double> tmp_458;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_459;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_459 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_458 += tmp_459;
   }
   tmp_457 += tmp_458;
   result += (-0.5) * tmp_457;
   std::complex<double> tmp_460;
   std::complex<double> tmp_461;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_462;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_462 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_461 += tmp_462;
   }
   tmp_460 += tmp_461;
   result += (-0.5) * tmp_460;
   std::complex<double> tmp_463;
   std::complex<double> tmp_464;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_464 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_463 += tmp_464;
   result += (-1) * tmp_463;
   std::complex<double> tmp_465;
   std::complex<double> tmp_466;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_466 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_465 += tmp_466;
   result += (-1.3333333333333333) * tmp_465;
   std::complex<double> tmp_467;
   std::complex<double> tmp_468;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_468 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_467 += tmp_468;
   result += (-1) * tmp_467;
   std::complex<double> tmp_469;
   std::complex<double> tmp_470;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_470 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_469 += tmp_470;
   result += (-1) * tmp_469;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_471;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_472;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_472 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_471 += tmp_472;
   }
   result += tmp_471;
   std::complex<double> tmp_473;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_474;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_474 += B0(p,MFv(gI2),MHm(gI1))*Conj(CpbarUFeHmFvPL(gO2,gI1,
            gI2))*CpbarUFeHmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_473 += tmp_474;
   }
   result += tmp_473;
   std::complex<double> tmp_475;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_476;
      std::complex<double> tmp_477;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_477 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_476 += tmp_477;
      tmp_475 += (MFe(gI1)) * tmp_476;
   }
   result += tmp_475;
   std::complex<double> tmp_478;
   std::complex<double> tmp_479;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_479 += (-0.5 + B0(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_478 += tmp_479;
   result += (-4) * tmp_478;
   std::complex<double> tmp_480;
   std::complex<double> tmp_481;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_481 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_480 += tmp_481;
   result += (-4) * tmp_480;
   std::complex<double> tmp_482;
   std::complex<double> tmp_483;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_483 += (-0.5 + B0(p,MFv(gI2),MVWm))*Conj(CpbarUFeVWmFvPR(gO2,gI2))
         *CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_482 += tmp_483;
   result += (-4) * tmp_482;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_484;
   std::complex<double> tmp_485;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_486;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_486 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_485 += tmp_486;
   }
   tmp_484 += tmp_485;
   result += (-0.5) * tmp_484;
   std::complex<double> tmp_487;
   std::complex<double> tmp_488;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_489;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_489 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarUFeHmFvPR(gO2,gI1,
            gI2))*CpbarUFeHmFvPR(gO1,gI1,gI2);
      }
      tmp_488 += tmp_489;
   }
   tmp_487 += tmp_488;
   result += (-0.5) * tmp_487;
   std::complex<double> tmp_490;
   std::complex<double> tmp_491;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_492;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_492 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_491 += tmp_492;
   }
   tmp_490 += tmp_491;
   result += (-0.5) * tmp_490;
   std::complex<double> tmp_493;
   std::complex<double> tmp_494;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_494 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_493 += tmp_494;
   result += (-1) * tmp_493;
   std::complex<double> tmp_495;
   std::complex<double> tmp_496;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_496 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_495 += tmp_496;
   result += (-1) * tmp_495;
   std::complex<double> tmp_497;
   std::complex<double> tmp_498;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_498 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_497 += tmp_498;
   result += (-1) * tmp_497;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_499;
   std::complex<double> tmp_500;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_501;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_501 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_500 += tmp_501;
   }
   tmp_499 += tmp_500;
   result += (-0.5) * tmp_499;
   std::complex<double> tmp_502;
   std::complex<double> tmp_503;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_504;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_504 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarUFeHmFvPL(gO2,gI1,
            gI2))*CpbarUFeHmFvPL(gO1,gI1,gI2);
      }
      tmp_503 += tmp_504;
   }
   tmp_502 += tmp_503;
   result += (-0.5) * tmp_502;
   std::complex<double> tmp_505;
   std::complex<double> tmp_506;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_507;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_507 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_506 += tmp_507;
   }
   tmp_505 += tmp_506;
   result += (-0.5) * tmp_505;
   std::complex<double> tmp_508;
   std::complex<double> tmp_509;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_509 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_508 += tmp_509;
   result += (-1) * tmp_508;
   std::complex<double> tmp_510;
   std::complex<double> tmp_511;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_511 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_510 += tmp_511;
   result += (-1) * tmp_510;
   std::complex<double> tmp_512;
   std::complex<double> tmp_513;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_513 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_512 += tmp_513;
   result += (-1) * tmp_512;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_514;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_515;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_515 += B0(p,MCha(gI2),MHm(gI1))*Conj(CpUChiconjHmChaPL(gO2,
            gI1,gI2))*CpUChiconjHmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_514 += tmp_515;
   }
   result += tmp_514;
   std::complex<double> tmp_516;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_517;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_517 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_516 += tmp_517;
   }
   result += tmp_516;
   std::complex<double> tmp_518;
   std::complex<double> tmp_519;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_519 += (-0.5 + B0(p,MCha(gI1),MVWm))*Conj(CpUChibarChaVWmPR(gO2,
         gI1))*CpUChibarChaVWmPL(gO1,gI1)*MCha(gI1);
   }
   tmp_518 += tmp_519;
   result += (-4) * tmp_518;
   std::complex<double> tmp_520;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_521;
      std::complex<double> tmp_522;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_522 += B0(p,MCha(gI1),MHm(gI2))*Conj(CpUChibarChaHmPL(gO2,
            gI1,gI2))*CpUChibarChaHmPR(gO1,gI1,gI2);
      }
      tmp_521 += tmp_522;
      tmp_520 += (MCha(gI1)) * tmp_521;
   }
   result += tmp_520;
   std::complex<double> tmp_523;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_524;
      std::complex<double> tmp_525;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_525 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_524 += tmp_525;
      tmp_523 += (MChi(gI1)) * tmp_524;
   }
   result += tmp_523;
   std::complex<double> tmp_526;
   std::complex<double> tmp_527;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_527 += (-0.5 + B0(p,MCha(gI2),MVWm))*Conj(CpUChiconjVWmChaPR(gO2,
         gI2))*CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_526 += tmp_527;
   result += (-4) * tmp_526;
   std::complex<double> tmp_528;
   std::complex<double> tmp_529;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_529 += (-0.5 + B0(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_528 += tmp_529;
   result += (-4) * tmp_528;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_530;
   std::complex<double> tmp_531;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_531 += (0.5 + B1(p,MCha(gI1),MVWm))*Conj(CpUChibarChaVWmPL(gO2,gI1
         ))*CpUChibarChaVWmPL(gO1,gI1);
   }
   tmp_530 += tmp_531;
   result += (-1) * tmp_530;
   std::complex<double> tmp_532;
   std::complex<double> tmp_533;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_534;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_534 += B1(p,MCha(gI1),MHm(gI2))*Conj(CpUChibarChaHmPR(gO2,
            gI1,gI2))*CpUChibarChaHmPR(gO1,gI1,gI2);
      }
      tmp_533 += tmp_534;
   }
   tmp_532 += tmp_533;
   result += (-0.5) * tmp_532;
   std::complex<double> tmp_535;
   std::complex<double> tmp_536;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_537;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_537 += B1(p,MCha(gI2),MHm(gI1))*Conj(CpUChiconjHmChaPR(gO2,
            gI1,gI2))*CpUChiconjHmChaPR(gO1,gI1,gI2);
      }
      tmp_536 += tmp_537;
   }
   tmp_535 += tmp_536;
   result += (-0.5) * tmp_535;
   std::complex<double> tmp_538;
   std::complex<double> tmp_539;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_540;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_540 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_539 += tmp_540;
   }
   tmp_538 += tmp_539;
   result += (-0.5) * tmp_538;
   std::complex<double> tmp_541;
   std::complex<double> tmp_542;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_543;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_543 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_542 += tmp_543;
   }
   tmp_541 += tmp_542;
   result += (-0.5) * tmp_541;
   std::complex<double> tmp_544;
   std::complex<double> tmp_545;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_545 += (0.5 + B1(p,MCha(gI2),MVWm))*Conj(CpUChiconjVWmChaPL(gO2,
         gI2))*CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_544 += tmp_545;
   result += (-1) * tmp_544;
   std::complex<double> tmp_546;
   std::complex<double> tmp_547;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_547 += (0.5 + B1(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_546 += tmp_547;
   result += (-1) * tmp_546;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_548;
   std::complex<double> tmp_549;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_549 += (0.5 + B1(p,MCha(gI1),MVWm))*Conj(CpUChibarChaVWmPR(gO2,gI1
         ))*CpUChibarChaVWmPR(gO1,gI1);
   }
   tmp_548 += tmp_549;
   result += (-1) * tmp_548;
   std::complex<double> tmp_550;
   std::complex<double> tmp_551;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_552;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_552 += B1(p,MCha(gI1),MHm(gI2))*Conj(CpUChibarChaHmPL(gO2,
            gI1,gI2))*CpUChibarChaHmPL(gO1,gI1,gI2);
      }
      tmp_551 += tmp_552;
   }
   tmp_550 += tmp_551;
   result += (-0.5) * tmp_550;
   std::complex<double> tmp_553;
   std::complex<double> tmp_554;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_555;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_555 += B1(p,MCha(gI2),MHm(gI1))*Conj(CpUChiconjHmChaPL(gO2,
            gI1,gI2))*CpUChiconjHmChaPL(gO1,gI1,gI2);
      }
      tmp_554 += tmp_555;
   }
   tmp_553 += tmp_554;
   result += (-0.5) * tmp_553;
   std::complex<double> tmp_556;
   std::complex<double> tmp_557;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_558;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_558 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_557 += tmp_558;
   }
   tmp_556 += tmp_557;
   result += (-0.5) * tmp_556;
   std::complex<double> tmp_559;
   std::complex<double> tmp_560;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_561;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_561 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_560 += tmp_561;
   }
   tmp_559 += tmp_560;
   result += (-0.5) * tmp_559;
   std::complex<double> tmp_562;
   std::complex<double> tmp_563;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_563 += (0.5 + B1(p,MCha(gI2),MVWm))*Conj(CpUChiconjVWmChaPR(gO2,
         gI2))*CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_562 += tmp_563;
   result += (-1) * tmp_562;
   std::complex<double> tmp_564;
   std::complex<double> tmp_565;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_565 += (0.5 + B1(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_564 += tmp_565;
   result += (-1) * tmp_564;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_566;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_567;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_567 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_566 += tmp_567;
   }
   result += tmp_566;
   std::complex<double> tmp_568;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_569;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_569 += B0(p,MChi(gI2),MHm(gI1))*Conj(CpbarUChaHmChiPL(gO2,
            gI1,gI2))*CpbarUChaHmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_568 += tmp_569;
   }
   result += tmp_568;
   std::complex<double> tmp_570;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_571;
      std::complex<double> tmp_572;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_572 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_571 += tmp_572;
      tmp_570 += (MCha(gI1)) * tmp_571;
   }
   result += tmp_570;
   std::complex<double> tmp_573;
   std::complex<double> tmp_574;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_574 += (-0.5 + B0(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_573 += tmp_574;
   result += (-4) * tmp_573;
   std::complex<double> tmp_575;
   std::complex<double> tmp_576;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_576 += (-0.5 + B0(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPR(gO2,gI2)
         )*CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_575 += tmp_576;
   result += (-4) * tmp_575;
   std::complex<double> tmp_577;
   std::complex<double> tmp_578;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_578 += (-0.5 + B0(p,MChi(gI2),MVWm))*Conj(CpbarUChaVWmChiPR(gO2,
         gI2))*CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_577 += tmp_578;
   result += (-4) * tmp_577;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_579;
   std::complex<double> tmp_580;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_581;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_581 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_580 += tmp_581;
   }
   tmp_579 += tmp_580;
   result += (-0.5) * tmp_579;
   std::complex<double> tmp_582;
   std::complex<double> tmp_583;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_584;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_584 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_583 += tmp_584;
   }
   tmp_582 += tmp_583;
   result += (-0.5) * tmp_582;
   std::complex<double> tmp_585;
   std::complex<double> tmp_586;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_587;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_587 += B1(p,MChi(gI2),MHm(gI1))*Conj(CpbarUChaHmChiPR(gO2,
            gI1,gI2))*CpbarUChaHmChiPR(gO1,gI1,gI2);
      }
      tmp_586 += tmp_587;
   }
   tmp_585 += tmp_586;
   result += (-0.5) * tmp_585;
   std::complex<double> tmp_588;
   std::complex<double> tmp_589;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_589 += (0.5 + B1(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_588 += tmp_589;
   result += (-1) * tmp_588;
   std::complex<double> tmp_590;
   std::complex<double> tmp_591;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_591 += (0.5 + B1(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPL(gO2,gI2))
         *CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_590 += tmp_591;
   result += (-1) * tmp_590;
   std::complex<double> tmp_592;
   std::complex<double> tmp_593;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_593 += (0.5 + B1(p,MChi(gI2),MVWm))*Conj(CpbarUChaVWmChiPL(gO2,gI2
         ))*CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_592 += tmp_593;
   result += (-1) * tmp_592;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_594;
   std::complex<double> tmp_595;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_596;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_596 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_595 += tmp_596;
   }
   tmp_594 += tmp_595;
   result += (-0.5) * tmp_594;
   std::complex<double> tmp_597;
   std::complex<double> tmp_598;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_599;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_599 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_598 += tmp_599;
   }
   tmp_597 += tmp_598;
   result += (-0.5) * tmp_597;
   std::complex<double> tmp_600;
   std::complex<double> tmp_601;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_602;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_602 += B1(p,MChi(gI2),MHm(gI1))*Conj(CpbarUChaHmChiPL(gO2,
            gI1,gI2))*CpbarUChaHmChiPL(gO1,gI1,gI2);
      }
      tmp_601 += tmp_602;
   }
   tmp_600 += tmp_601;
   result += (-0.5) * tmp_600;
   std::complex<double> tmp_603;
   std::complex<double> tmp_604;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_604 += (0.5 + B1(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_603 += tmp_604;
   result += (-1) * tmp_603;
   std::complex<double> tmp_605;
   std::complex<double> tmp_606;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_606 += (0.5 + B1(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPR(gO2,gI2))
         *CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_605 += tmp_606;
   result += (-1) * tmp_605;
   std::complex<double> tmp_607;
   std::complex<double> tmp_608;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_608 += (0.5 + B1(p,MChi(gI2),MVWm))*Conj(CpbarUChaVWmChiPR(gO2,gI2
         ))*CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_607 += tmp_608;
   result += (-1) * tmp_607;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_609;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_610;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_610 += B0(p,MFe(gI2),MHm(gI1))*Conj(CpbarFvconjHmFePL(gO2,
            gI1,gI2))*CpbarFvconjHmFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_609 += tmp_610;
   }
   result += tmp_609;
   std::complex<double> tmp_611;
   std::complex<double> tmp_612;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_612 += (-0.5 + B0(p,MFe(gI2),MVWm))*Conj(CpbarFvconjVWmFePR(gO2,
         gI2))*CpbarFvconjVWmFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_611 += tmp_612;
   result += (-4) * tmp_611;
   std::complex<double> tmp_613;
   std::complex<double> tmp_614;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_614 += (-0.5 + B0(p,MFv(gI2),MVZ))*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_613 += tmp_614;
   result += (-4) * tmp_613;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_615;
   std::complex<double> tmp_616;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_617;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_617 += B1(p,MFe(gI2),MHm(gI1))*Conj(CpbarFvconjHmFePR(gO2,
            gI1,gI2))*CpbarFvconjHmFePR(gO1,gI1,gI2);
      }
      tmp_616 += tmp_617;
   }
   tmp_615 += tmp_616;
   result += (-0.5) * tmp_615;
   std::complex<double> tmp_618;
   std::complex<double> tmp_619;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_619 += (0.5 + B1(p,MFe(gI2),MVWm))*Conj(CpbarFvconjVWmFePL(gO2,gI2
         ))*CpbarFvconjVWmFePL(gO1,gI2);
   }
   tmp_618 += tmp_619;
   result += (-1) * tmp_618;
   std::complex<double> tmp_620;
   std::complex<double> tmp_621;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_621 += (0.5 + B1(p,MFv(gI2),MVZ))*Conj(CpbarFvVZFvPL(gO2,gI2))*
         CpbarFvVZFvPL(gO1,gI2);
   }
   tmp_620 += tmp_621;
   result += (-1) * tmp_620;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_622;
   std::complex<double> tmp_623;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_624;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_624 += B1(p,MFe(gI2),MHm(gI1))*Conj(CpbarFvconjHmFePL(gO2,
            gI1,gI2))*CpbarFvconjHmFePL(gO1,gI1,gI2);
      }
      tmp_623 += tmp_624;
   }
   tmp_622 += tmp_623;
   result += (-0.5) * tmp_622;
   std::complex<double> tmp_625;
   std::complex<double> tmp_626;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_626 += (0.5 + B1(p,MFe(gI2),MVWm))*Conj(CpbarFvconjVWmFePR(gO2,gI2
         ))*CpbarFvconjVWmFePR(gO1,gI2);
   }
   tmp_625 += tmp_626;
   result += (-1) * tmp_625;
   std::complex<double> tmp_627;
   std::complex<double> tmp_628;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_628 += (0.5 + B1(p,MFv(gI2),MVZ))*Conj(CpbarFvVZFvPR(gO2,gI2))*
         CpbarFvVZFvPR(gO1,gI2);
   }
   tmp_627 += tmp_628;
   result += (-1) * tmp_627;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   result += -12*MGlu*(-0.5 + B0(p,MGlu,0))*Conj(CpGluVGGluPR())*CpGluVGGluPL()
      ;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPL())*(0.5 + B1(p,MGlu,0));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPR())*(0.5 + B1(p,MGlu,0));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_629;
   std::complex<double> tmp_630;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_630 += AbsSqr(CpVZhhAh(gI1,1))*B00(p,MAh(1),Mhh(gI1));
   }
   tmp_629 += tmp_630;
   result += (-4) * tmp_629;
   std::complex<double> tmp_631;
   std::complex<double> tmp_632;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_632 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_631 += tmp_632;
   result += (0.5) * tmp_631;
   std::complex<double> tmp_633;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_634;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_634 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_634 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_633 += tmp_634;
   }
   result += tmp_633;
   std::complex<double> tmp_635;
   std::complex<double> tmp_636;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_637;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_637 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR(
            gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_637 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_636 += tmp_637;
   }
   tmp_635 += tmp_636;
   result += (0.5) * tmp_635;
   std::complex<double> tmp_638;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_638 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_638;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_639;
   std::complex<double> tmp_640;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_640 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_639 += tmp_640;
   result += (0.5) * tmp_639;
   std::complex<double> tmp_641;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_642;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_642 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_642 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_641 += tmp_642;
   }
   result += tmp_641;
   std::complex<double> tmp_643;
   std::complex<double> tmp_644;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_644 += AbsSqr(CpconjVWmHmhh(1,gI2))*B00(p,Mhh(gI2),MHm(1));
   }
   tmp_643 += tmp_644;
   result += (-4) * tmp_643;
   std::complex<double> tmp_645;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_645 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_645;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_646;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_647;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_647 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_646 += tmp_647;
   }
   result += tmp_646;
   std::complex<double> tmp_648;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_649;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_649 += B0(p,MFu(gI2),MHm(gI1))*Conj(CpbarFdHmFuPL(gO2,gI1,
            gI2))*CpbarFdHmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_648 += tmp_649;
   }
   result += tmp_648;
   std::complex<double> tmp_650;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_651;
      std::complex<double> tmp_652;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_652 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_651 += tmp_652;
      tmp_650 += (MFd(gI1)) * tmp_651;
   }
   result += tmp_650;
   std::complex<double> tmp_653;
   std::complex<double> tmp_654;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_654 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_653 += tmp_654;
   result += (-4) * tmp_653;
   std::complex<double> tmp_655;
   std::complex<double> tmp_656;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_656 += (-0.5 + B0(p,MFu(gI2),MVWm))*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_655 += tmp_656;
   result += (-4) * tmp_655;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_657;
   std::complex<double> tmp_658;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_659;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_659 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_658 += tmp_659;
   }
   tmp_657 += tmp_658;
   result += (-0.5) * tmp_657;
   std::complex<double> tmp_660;
   std::complex<double> tmp_661;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_662;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_662 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarFdHmFuPR(gO2,gI1,
            gI2))*CpbarFdHmFuPR(gO1,gI1,gI2);
      }
      tmp_661 += tmp_662;
   }
   tmp_660 += tmp_661;
   result += (-0.5) * tmp_660;
   std::complex<double> tmp_663;
   std::complex<double> tmp_664;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_665;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_665 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_664 += tmp_665;
   }
   tmp_663 += tmp_664;
   result += (-0.5) * tmp_663;
   std::complex<double> tmp_666;
   std::complex<double> tmp_667;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_667 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_666 += tmp_667;
   result += (-1) * tmp_666;
   std::complex<double> tmp_668;
   std::complex<double> tmp_669;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_669 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_668 += tmp_669;
   result += (-1) * tmp_668;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_670;
   std::complex<double> tmp_671;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_672;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_672 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_671 += tmp_672;
   }
   tmp_670 += tmp_671;
   result += (-0.5) * tmp_670;
   std::complex<double> tmp_673;
   std::complex<double> tmp_674;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_675;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_675 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarFdHmFuPL(gO2,gI1,
            gI2))*CpbarFdHmFuPL(gO1,gI1,gI2);
      }
      tmp_674 += tmp_675;
   }
   tmp_673 += tmp_674;
   result += (-0.5) * tmp_673;
   std::complex<double> tmp_676;
   std::complex<double> tmp_677;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_678;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_678 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_677 += tmp_678;
   }
   tmp_676 += tmp_677;
   result += (-0.5) * tmp_676;
   std::complex<double> tmp_679;
   std::complex<double> tmp_680;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_680 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_679 += tmp_680;
   result += (-1) * tmp_679;
   std::complex<double> tmp_681;
   std::complex<double> tmp_682;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_682 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_681 += tmp_682;
   result += (-1) * tmp_681;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_683;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_684;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_684 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_683 += tmp_684;
   }
   result += tmp_683;
   std::complex<double> tmp_685;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_686;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_686 += B0(p,MFv(gI2),MHm(gI1))*Conj(CpbarFeHmFvPL(gO2,gI1,
            gI2))*CpbarFeHmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_685 += tmp_686;
   }
   result += tmp_685;
   std::complex<double> tmp_687;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_688;
      std::complex<double> tmp_689;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_689 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_688 += tmp_689;
      tmp_687 += (MFe(gI1)) * tmp_688;
   }
   result += tmp_687;
   std::complex<double> tmp_690;
   std::complex<double> tmp_691;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_691 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_690 += tmp_691;
   result += (-4) * tmp_690;
   std::complex<double> tmp_692;
   std::complex<double> tmp_693;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_693 += (-0.5 + B0(p,MFv(gI2),MVWm))*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_692 += tmp_693;
   result += (-4) * tmp_692;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_694;
   std::complex<double> tmp_695;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_696;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_696 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_695 += tmp_696;
   }
   tmp_694 += tmp_695;
   result += (-0.5) * tmp_694;
   std::complex<double> tmp_697;
   std::complex<double> tmp_698;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_699;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_699 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarFeHmFvPR(gO2,gI1,
            gI2))*CpbarFeHmFvPR(gO1,gI1,gI2);
      }
      tmp_698 += tmp_699;
   }
   tmp_697 += tmp_698;
   result += (-0.5) * tmp_697;
   std::complex<double> tmp_700;
   std::complex<double> tmp_701;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_702;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_702 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_701 += tmp_702;
   }
   tmp_700 += tmp_701;
   result += (-0.5) * tmp_700;
   std::complex<double> tmp_703;
   std::complex<double> tmp_704;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_704 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_703 += tmp_704;
   result += (-1) * tmp_703;
   std::complex<double> tmp_705;
   std::complex<double> tmp_706;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_706 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_705 += tmp_706;
   result += (-1) * tmp_705;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_707;
   std::complex<double> tmp_708;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_709;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_709 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_708 += tmp_709;
   }
   tmp_707 += tmp_708;
   result += (-0.5) * tmp_707;
   std::complex<double> tmp_710;
   std::complex<double> tmp_711;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_712;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_712 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarFeHmFvPL(gO2,gI1,
            gI2))*CpbarFeHmFvPL(gO1,gI1,gI2);
      }
      tmp_711 += tmp_712;
   }
   tmp_710 += tmp_711;
   result += (-0.5) * tmp_710;
   std::complex<double> tmp_713;
   std::complex<double> tmp_714;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_715;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_715 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_714 += tmp_715;
   }
   tmp_713 += tmp_714;
   result += (-0.5) * tmp_713;
   std::complex<double> tmp_716;
   std::complex<double> tmp_717;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_717 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_716 += tmp_717;
   result += (-1) * tmp_716;
   std::complex<double> tmp_718;
   std::complex<double> tmp_719;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_719 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_718 += tmp_719;
   result += (-1) * tmp_718;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_720;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_721;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_721 += B0(p,MFd(gI2),MHm(gI1))*Conj(CpbarFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarFuconjHmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_720 += tmp_721;
   }
   result += tmp_720;
   std::complex<double> tmp_722;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_723;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_723 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_722 += tmp_723;
   }
   result += tmp_722;
   std::complex<double> tmp_724;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_725;
      std::complex<double> tmp_726;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_726 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_725 += tmp_726;
      tmp_724 += (MFu(gI1)) * tmp_725;
   }
   result += tmp_724;
   std::complex<double> tmp_727;
   std::complex<double> tmp_728;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_728 += (-0.5 + B0(p,MFd(gI2),MVWm))*Conj(CpbarFuconjVWmFdPR(gO2,
         gI2))*CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_727 += tmp_728;
   result += (-4) * tmp_727;
   std::complex<double> tmp_729;
   std::complex<double> tmp_730;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_730 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_729 += tmp_730;
   result += (-4) * tmp_729;
   std::complex<double> tmp_731;
   std::complex<double> tmp_732;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_732 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_731 += tmp_732;
   result += (-4) * tmp_731;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_733;
   std::complex<double> tmp_734;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_735;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_735 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarFuconjHmFdPR(gO2,
            gI1,gI2))*CpbarFuconjHmFdPR(gO1,gI1,gI2);
      }
      tmp_734 += tmp_735;
   }
   tmp_733 += tmp_734;
   result += (-0.5) * tmp_733;
   std::complex<double> tmp_736;
   std::complex<double> tmp_737;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_738;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_738 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_737 += tmp_738;
   }
   tmp_736 += tmp_737;
   result += (-0.5) * tmp_736;
   std::complex<double> tmp_739;
   std::complex<double> tmp_740;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_741;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_741 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_740 += tmp_741;
   }
   tmp_739 += tmp_740;
   result += (-0.5) * tmp_739;
   std::complex<double> tmp_742;
   std::complex<double> tmp_743;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_743 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarFuconjVWmFdPL(gO2,gI2
         ))*CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_742 += tmp_743;
   result += (-1) * tmp_742;
   std::complex<double> tmp_744;
   std::complex<double> tmp_745;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_745 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_744 += tmp_745;
   result += (-1) * tmp_744;
   std::complex<double> tmp_746;
   std::complex<double> tmp_747;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_747 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_746 += tmp_747;
   result += (-1) * tmp_746;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_748;
   std::complex<double> tmp_749;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_750;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_750 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarFuconjHmFdPL(gO1,gI1,gI2);
      }
      tmp_749 += tmp_750;
   }
   tmp_748 += tmp_749;
   result += (-0.5) * tmp_748;
   std::complex<double> tmp_751;
   std::complex<double> tmp_752;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_753;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_753 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_752 += tmp_753;
   }
   tmp_751 += tmp_752;
   result += (-0.5) * tmp_751;
   std::complex<double> tmp_754;
   std::complex<double> tmp_755;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_756;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_756 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_755 += tmp_756;
   }
   tmp_754 += tmp_755;
   result += (-0.5) * tmp_754;
   std::complex<double> tmp_757;
   std::complex<double> tmp_758;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_758 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarFuconjVWmFdPR(gO2,gI2
         ))*CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_757 += tmp_758;
   result += (-1) * tmp_757;
   std::complex<double> tmp_759;
   std::complex<double> tmp_760;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_760 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_759 += tmp_760;
   result += (-1) * tmp_759;
   std::complex<double> tmp_761;
   std::complex<double> tmp_762;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_762 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_761 += tmp_762;
   result += (-1) * tmp_761;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_763;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_764;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_764 += B0(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_763 += tmp_764;
   }
   result += tmp_763;
   std::complex<double> tmp_765;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_766;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_766 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_765 += tmp_766;
   }
   result += tmp_765;
   std::complex<double> tmp_767;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_768;
      std::complex<double> tmp_769;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_769 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_768 += tmp_769;
      tmp_767 += (MFu(gI1)) * tmp_768;
   }
   result += tmp_767;
   std::complex<double> tmp_770;
   std::complex<double> tmp_771;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_771 += (-0.5 + B0(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_770 += tmp_771;
   result += (-4) * tmp_770;
   std::complex<double> tmp_772;
   std::complex<double> tmp_773;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_773 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_772 += tmp_773;
   result += (-4) * tmp_772;
   std::complex<double> tmp_774;
   std::complex<double> tmp_775;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_775 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_774 += tmp_775;
   result += (-4) * tmp_774;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_776;
   std::complex<double> tmp_777;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_778;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_778 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPR(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2);
      }
      tmp_777 += tmp_778;
   }
   tmp_776 += tmp_777;
   result += (-0.5) * tmp_776;
   std::complex<double> tmp_779;
   std::complex<double> tmp_780;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_781;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_781 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_780 += tmp_781;
   }
   tmp_779 += tmp_780;
   result += (-0.5) * tmp_779;
   std::complex<double> tmp_782;
   std::complex<double> tmp_783;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_784;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_784 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_783 += tmp_784;
   }
   tmp_782 += tmp_783;
   result += (-0.5) * tmp_782;
   std::complex<double> tmp_785;
   std::complex<double> tmp_786;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_786 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPL(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_785 += tmp_786;
   result += (-1) * tmp_785;
   std::complex<double> tmp_787;
   std::complex<double> tmp_788;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_788 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_787 += tmp_788;
   result += (-1) * tmp_787;
   std::complex<double> tmp_789;
   std::complex<double> tmp_790;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_790 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_789 += tmp_790;
   result += (-1) * tmp_789;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_791;
   std::complex<double> tmp_792;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_793;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_793 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPL(gO1,gI1,gI2);
      }
      tmp_792 += tmp_793;
   }
   tmp_791 += tmp_792;
   result += (-0.5) * tmp_791;
   std::complex<double> tmp_794;
   std::complex<double> tmp_795;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_796;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_796 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_795 += tmp_796;
   }
   tmp_794 += tmp_795;
   result += (-0.5) * tmp_794;
   std::complex<double> tmp_797;
   std::complex<double> tmp_798;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_799;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_799 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_798 += tmp_799;
   }
   tmp_797 += tmp_798;
   result += (-0.5) * tmp_797;
   std::complex<double> tmp_800;
   std::complex<double> tmp_801;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_801 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_800 += tmp_801;
   result += (-1) * tmp_800;
   std::complex<double> tmp_802;
   std::complex<double> tmp_803;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_803 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_802 += tmp_803;
   result += (-1) * tmp_802;
   std::complex<double> tmp_804;
   std::complex<double> tmp_805;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_805 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_804 += tmp_805;
   result += (-1) * tmp_804;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh(unsigned gO1) const
{
   std::complex<double> result;

   result += A0(MVWm)*CpUhhbargWmCgWmC(gO1);
   result += A0(MVWm)*CpUhhbargWmgWm(gO1);
   result += A0(MVZ)*CpUhhbargZgZ(gO1);
   std::complex<double> tmp_806;
   std::complex<double> tmp_807;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_807 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_806 += tmp_807;
   result += (-0.5) * tmp_806;
   std::complex<double> tmp_808;
   std::complex<double> tmp_809;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_809 += A0(MHm(gI1))*CpUhhconjHmHm(gO1,gI1,gI1);
   }
   tmp_808 += tmp_809;
   result += (-1) * tmp_808;
   std::complex<double> tmp_810;
   std::complex<double> tmp_811;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_811 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_810 += tmp_811;
   result += (-0.5) * tmp_810;
   std::complex<double> tmp_812;
   std::complex<double> tmp_813;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_813 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_812 += tmp_813;
   result += (2) * tmp_812;
   std::complex<double> tmp_814;
   std::complex<double> tmp_815;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_815 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR(
         gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_814 += tmp_815;
   result += (6) * tmp_814;
   std::complex<double> tmp_816;
   std::complex<double> tmp_817;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_817 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR(
         gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_816 += tmp_817;
   result += (2) * tmp_816;
   std::complex<double> tmp_818;
   std::complex<double> tmp_819;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_819 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR(
         gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_818 += tmp_819;
   result += (6) * tmp_818;
   std::complex<double> tmp_820;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_820 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_820;
   result += 4*CpUhhconjVWmVWm(gO1)*(A0(MVWm) - 0.5*Sqr(MVWm));
   result += 2*CpUhhVZVZ(gO1)*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}










void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void CLASSNAME::calculate_MGlu_pole()
{
   // diagonalization with medium precision
   const double M_tree(MGlu);
   const double p = MGlu;
   const double self_energy_1  = Re(self_energy_Glu_1(p));
   const double self_energy_PL = Re(self_energy_Glu_PL(p));
   const double self_energy_PR = Re(self_energy_Glu_PR(p));
   const auto M_loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MGlu) = calculate_singlet_mass(M_loop);
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
}

void CLASSNAME::calculate_MVZ_pole()
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::hh))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(Mhh) old_Mhh(Mhh), new_Mhh(Mhh);

   do {
      Eigen::Matrix<double,2,2> self_energy;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_hh());

      for (unsigned es = 0; es < 2; ++es) {
         const double p = Abs(old_Mhh(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = i1; i2 < 2; ++i2) {
               self_energy(i1,i2) = Re(self_energy_hh(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,2,2> M_loop(M_tree -
            self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(HGTHDMIIMSSMBC_info::hh,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_loop, eigen_values,
               mix_ZH);
         #endif

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
      problems.flag_no_pole_mass_convergence(HGTHDMIIMSSMBC_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(HGTHDMIIMSSMBC_info::hh
         );
}

void CLASSNAME::calculate_MAh_pole()
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::Ah))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,2,2> self_energy;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Ah());

   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MAh(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = i1; i2 < 2; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Ah(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZA;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA,
            eigenvalue_error);
         problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Ah,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZA);
      #endif

      PHYSICAL(MAh(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 1)
         PHYSICAL(ZA) = mix_ZA;
   }
}

void CLASSNAME::calculate_MHm_pole()
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::Hm))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,2,2> self_energy;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Hm());

   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MHm(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = i1; i2 < 2; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Hm(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,2,2> M_loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZP;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP,
            eigenvalue_error);
         problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Hm,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_loop, eigen_values, mix_ZP);
      #endif

      PHYSICAL(MHm(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 1)
         PHYSICAL(ZP) = mix_ZP;
   }
}

void CLASSNAME::calculate_MFd_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFd(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fd_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fd_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fd_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vd) mix_Vd;
      decltype(Ud) mix_Ud;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud, eigenvalue_error);
      problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Fd, eigenvalue_error
         > precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud);
   #endif
      if (es == 0) {
         PHYSICAL(Vd) = mix_Vd;
         PHYSICAL(Ud) = mix_Ud;
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
      qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(2))/Sqr(
         currentScale)))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = -0.005284774766427138*Power(g3,4) -
         0.0032348537833770956*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(2))) -
         0.0008822328500119351*Power(g3,4)*Sqr(Log(Power(currentScale,2)/Sqr(
         MFu(2))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = -0.00003352082872926087*Power(g3,6)*(35.702577217116016
         + 15.387410814884797*Log(Sqr(currentScale)/Sqr(MFu(2))) + 1.*Power(
         Log(Sqr(currentScale)/Sqr(MFu(2))),3) + 5.378787878787879*Sqr(Log(
         Power(currentScale,2)/Sqr(MFu(2)))));
   }

   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            if (i1 == 2 && i2 == 2) {
               self_energy_1(i1,i2)  = Re(
                  self_energy_Fu_1_heavy(p,i1,i2));
               self_energy_PL(i1,i2) = Re(
                  self_energy_Fu_PL_heavy(p,i1,i2));
               self_energy_PR(i1,i2) = Re(
                  self_energy_Fu_PR_heavy(p,i1,i2));
            } else {
               self_energy_1(i1,i2)  = Re(self_energy_Fu_1(p,
                  i1,i2));
               self_energy_PL(i1,i2) = Re(self_energy_Fu_PL(p
                  ,i1,i2));
               self_energy_PR(i1,i2) = Re(self_energy_Fu_PR(p
                  ,i1,i2));
            }
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vu) mix_Vu;
      decltype(Uu) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Fu, eigenvalue_error
         > precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Vu) = mix_Vu;
         PHYSICAL(Uu) = mix_Uu;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MFe_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFe(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fe_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fe_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fe_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Ve;
      decltype(Ue) mix_Ue;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue, eigenvalue_error);
      problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Fe, eigenvalue_error
         > precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Ve;
         PHYSICAL(Ue) = mix_Ue;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,4,4> self_energy_1;
   Eigen::Matrix<double,4,4> self_energy_PL;
   Eigen::Matrix<double,4,4> self_energy_PR;
   const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Chi());
   for (unsigned es = 0; es < 4; ++es) {
      const double p = Abs(MChi(es));
      for (unsigned i1 = 0; i1 < 4; ++i1) {
         for (unsigned i2 = 0; i2 < 4; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Chi_1(p,i1,i2
               ));
            self_energy_PL(i1,i2) = Re(self_energy_Chi_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Chi_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,4,4> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,4,4> M_loop(M_tree + 0.5 * (delta_M +
         delta_M.transpose()));
      Eigen::Array<double,4,1> eigen_values;
      decltype(ZN) mix_ZN;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZN,
            eigenvalue_error);
         problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Chi,
            eigenvalue_error > precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_symmetric(M_loop, eigen_values, mix_ZN);
      #endif
      if (es == 0)
         PHYSICAL(ZN) = mix_ZN;
      PHYSICAL(MChi(es)) = Abs(eigen_values(es));
   }
}

void CLASSNAME::calculate_MCha_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,2,2> self_energy_1;
   Eigen::Matrix<double,2,2> self_energy_PL;
   Eigen::Matrix<double,2,2> self_energy_PR;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha());
   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MCha(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = 0; i2 < 2; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Cha_1(p,i1,i2
               ));
            self_energy_PL(i1,i2) = Re(self_energy_Cha_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Cha_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM) mix_UM;
      decltype(UP) mix_UP;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_UM, mix_UP, eigenvalue_error);
      problems.flag_bad_mass(HGTHDMIIMSSMBC_info::Cha,
         eigenvalue_error > precision * Abs(eigen_values(0)));
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

void CLASSNAME::calculate_MVWm_pole()
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VWm);

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_tachyon(HGTHDMIIMSSMBC_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VZ);

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double CLASSNAME::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fe_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fe_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_PR_heavy_rotated(p,
      idx, idx));
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar + self_energy_1 + m_sm_drbar *
      (self_energy_PL + self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fu_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fu_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_PR_heavy_rotated(p,
      idx, idx));

   const double currentScale = get_scale();
   const double qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(idx))
      /Sqr(currentScale)))*Sqr(g3);
   double qcd_2l = 0., qcd_3l = 0.;

   if (get_thresholds() > 1) {
      qcd_2l = -0.0041441100714622115*Power(g3,4) -
         0.0015238567409297061*Power(g3,4)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(currentScale,2)
         /Sqr(MFu(idx))));
   }

   if (get_thresholds() > 2) {
      qcd_3l = -0.0008783313853540776*Power(g3,6) -
         0.0004114970933517977*Power(g3,6)*Log(Sqr(currentScale)/Sqr(MFu(idx)))
         - 5.078913443827405e-6*Power(g3,6)*Power(Log(Sqr(currentScale)/Sqr(
         MFu(idx))),3) - 0.0002952541682011665*Power(g3,6)*Log(Sqr(MFu(idx))
         /Sqr(currentScale)) + 0.00005282069981580501*Power(g3,6)*Sqr(Log(Power
         (MFu(idx),2)/Sqr(currentScale))) - 0.00007466002762426286*Power(g3,6)*
         Sqr(Log(Power(currentScale,2)/Sqr(MFu(idx))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p = m_sm_msbar;
   const double self_energy_1  = Re(self_energy_Fd_1_heavy_rotated(p, idx
      , idx));
   const double self_energy_PL = Re(self_energy_Fd_PL_heavy_rotated(p,
      idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_PR_heavy_rotated(p,
      idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree -
      self_energy_PL - self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MVP_DRbar(double)
{
   return 0.0;
}

double CLASSNAME::calculate_MVZ_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VZ);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_DRbar(double m_pole)
{
   const double p = m_pole;
   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = Sqr(m_pole) + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_tachyon(HGTHDMIIMSSMBC_info::VWm);
      return m_pole;
   }

   return AbsSqrt(mass_sqr);
}


double CLASSNAME::Betax() const
{
   return ArcSin(Abs(ZP(0,1)));
}

double CLASSNAME::Alpha() const
{
   return ArcCos(ZH(0,1));
}

double CLASSNAME::ThetaW() const
{
   return ArcCos(Abs(ZZ(0,0)));
}


std::ostream& operator<<(std::ostream& ostr, const HGTHDMIIMSSMBC_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
