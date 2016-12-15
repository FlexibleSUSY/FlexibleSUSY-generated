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

// File generated at Thu 15 Dec 2016 12:40:55

/**
 * @file HGTHDMIIMSSMBC_mass_eigenstates.cpp
 * @brief implementation of the HGTHDMIIMSSMBC model class
 *
 * Contains the definition of the HGTHDMIIMSSMBC model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Thu 15 Dec 2016 12:40:55 with FlexibleSUSY
 * 1.7.2 (git commit: 0d19299fef514160cb7541a03abb9b2c3365f927) and SARAH 4.9.1 .
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

using namespace HGTHDMIIMSSMBC_info;

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
   if (PHYSICAL(Mhh).tail<2>().minCoeff() < 0.) problems.flag_tachyon(hh);
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) problems.flag_tachyon(Ah);
   if (PHYSICAL(MHm).tail<1>().minCoeff() < 0.) problems.flag_tachyon(Hm);

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
   Eigen::Array<double,1,1> MHm_ChargedHiggs;
   Eigen::Array<double,1,1> MHm_goldstone;

   MHm_goldstone(0) = MVWm;

   remove_if_equal(MHm, MHm_goldstone, MHm_ChargedHiggs);

   return MHm_ChargedHiggs;
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_PseudoscalarHiggs;
   Eigen::Array<double,1,1> MAh_goldstone;

   MAh_goldstone(0) = MVZ;

   remove_if_equal(MAh, MAh_goldstone, MAh_PseudoscalarHiggs);

   return MAh_PseudoscalarHiggs;
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

   std::complex<double> tmp_92;
   std::complex<double> tmp_93;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_94;
      std::complex<double> tmp_95;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_95 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_94 += tmp_95;
      tmp_93 += (Conj(Vd(gI2,j2))) * tmp_94;
   }
   tmp_92 += tmp_93;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_92;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_96;
   std::complex<double> tmp_97;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_98;
      std::complex<double> tmp_99;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_99 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_98 += tmp_99;
      tmp_97 += (Vd(gO1,j2)) * tmp_98;
   }
   tmp_96 += tmp_97;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_96;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_100;
   std::complex<double> tmp_101;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_102;
      std::complex<double> tmp_103;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_103 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_102 += tmp_103;
      tmp_101 += (Conj(Vu(gI2,j2))) * tmp_102;
   }
   tmp_100 += tmp_101;
   result += (-ZP(gI1,0)) * tmp_100;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_104;
   std::complex<double> tmp_105;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_106;
      std::complex<double> tmp_107;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_107 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_106 += tmp_107;
      tmp_105 += (Vd(gO1,j2)) * tmp_106;
   }
   tmp_104 += tmp_105;
   result += (-ZP(gI1,1)) * tmp_104;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_108;
   std::complex<double> tmp_109;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_110;
      std::complex<double> tmp_111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_111 += Conj(Ud(gO2,j1))*Yd(j1,j2);
      }
      tmp_110 += tmp_111;
      tmp_109 += (Conj(Vd(gI1,j2))) * tmp_110;
   }
   tmp_108 += tmp_109;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_108;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_112;
   std::complex<double> tmp_113;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_114;
      std::complex<double> tmp_115;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_115 += Conj(Yd(j1,j2))*Ud(gI1,j1);
      }
      tmp_114 += tmp_115;
      tmp_113 += (Vd(gO1,j2)) * tmp_114;
   }
   tmp_112 += tmp_113;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_112
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

   std::complex<double> tmp_116;
   std::complex<double> tmp_117;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_117 += Conj(Vu(gI2,j1))*Vd(gO1,j1);
   }
   tmp_116 += tmp_117;
   result += (-0.7071067811865475*g2) * tmp_116;

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

   std::complex<double> tmp_118;
   std::complex<double> tmp_119;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_120;
      std::complex<double> tmp_121;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_121 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_120 += tmp_121;
      tmp_119 += (Conj(Ve(gI2,j2))) * tmp_120;
   }
   tmp_118 += tmp_119;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_118;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_122;
   std::complex<double> tmp_123;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_124;
      std::complex<double> tmp_125;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_125 += Conj(Ye(j1,j2))*Ue(gI2,j1);
      }
      tmp_124 += tmp_125;
      tmp_123 += (Ve(gO1,j2)) * tmp_124;
   }
   tmp_122 += tmp_123;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_122;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_126;
   std::complex<double> tmp_127;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_127 += Conj(Ue(gO2,j1))*Ye(j1,gI2);
   }
   tmp_126 += tmp_127;
   result += (-ZP(gI1,0)) * tmp_126;

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

   std::complex<double> tmp_128;
   std::complex<double> tmp_129;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_130;
      std::complex<double> tmp_131;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_131 += Conj(Ue(gO2,j1))*Ye(j1,j2);
      }
      tmp_130 += tmp_131;
      tmp_129 += (Conj(Ve(gI1,j2))) * tmp_130;
   }
   tmp_128 += tmp_129;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,0)) * tmp_128;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_132;
   std::complex<double> tmp_133;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_135 += Conj(Ye(j1,j2))*Ue(gI1,j1);
      }
      tmp_134 += tmp_135;
      tmp_133 += (Ve(gO1,j2)) * tmp_134;
   }
   tmp_132 += tmp_133;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,0)) * tmp_132
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

   std::complex<double> tmp_136;
   std::complex<double> tmp_137;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_138;
      std::complex<double> tmp_139;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_139 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_138 += tmp_139;
      tmp_137 += (Conj(Vd(gI2,j2))) * tmp_138;
   }
   tmp_136 += tmp_137;
   result += (-ZP(gI1,1)) * tmp_136;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_140;
   std::complex<double> tmp_141;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_142;
      std::complex<double> tmp_143;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_143 += Conj(Yd(j1,j2))*Ud(gI2,j1);
      }
      tmp_142 += tmp_143;
      tmp_141 += (Vu(gO1,j2)) * tmp_142;
   }
   tmp_140 += tmp_141;
   result += (-ZP(gI1,0)) * tmp_140;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_144;
   std::complex<double> tmp_145;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_146;
      std::complex<double> tmp_147;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_147 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_146 += tmp_147;
      tmp_145 += (Conj(Vu(gI2,j2))) * tmp_146;
   }
   tmp_144 += tmp_145;
   result += (0.7071067811865475*ZH(gI1,1)) * tmp_144;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_148;
   std::complex<double> tmp_149;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_150;
      std::complex<double> tmp_151;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_151 += Conj(Yu(j1,j2))*Uu(gI2,j1);
      }
      tmp_150 += tmp_151;
      tmp_149 += (Vu(gO1,j2)) * tmp_150;
   }
   tmp_148 += tmp_149;
   result += (0.7071067811865475*ZH(gI1,1)) * tmp_148;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_152;
   std::complex<double> tmp_153;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_154;
      std::complex<double> tmp_155;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_155 += Conj(Uu(gO2,j1))*Yu(j1,j2);
      }
      tmp_154 += tmp_155;
      tmp_153 += (Conj(Vu(gI1,j2))) * tmp_154;
   }
   tmp_152 += tmp_153;
   result += (std::complex<double>(0.,0.7071067811865475)*ZA(gI2,1)) * tmp_152;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_156;
   std::complex<double> tmp_157;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_158;
      std::complex<double> tmp_159;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_159 += Conj(Yu(j1,j2))*Uu(gI1,j1);
      }
      tmp_158 += tmp_159;
      tmp_157 += (Vu(gO1,j2)) * tmp_158;
   }
   tmp_156 += tmp_157;
   result += (std::complex<double>(0.,-0.7071067811865475)*ZA(gI2,1)) * tmp_156
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

   std::complex<double> tmp_160;
   std::complex<double> tmp_161;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_161 += Conj(Vd(gI2,j1))*Vu(gO1,j1);
   }
   tmp_160 += tmp_161;
   result += (-0.7071067811865475*g2) * tmp_160;

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
   std::complex<double> tmp_162;
   std::complex<double> tmp_163;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_163 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_162 += tmp_163;
   result += (-0.5) * tmp_162;
   std::complex<double> tmp_164;
   std::complex<double> tmp_165;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_165 += A0(MHm(gI1))*CpUhhUhhconjHmHm(gO1,gO2,gI1,gI1);
   }
   tmp_164 += tmp_165;
   result += (-1) * tmp_164;
   std::complex<double> tmp_166;
   std::complex<double> tmp_167;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_167 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_166 += tmp_167;
   result += (-0.5) * tmp_166;
   std::complex<double> tmp_168;
   std::complex<double> tmp_169;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_170;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_170 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))*
            CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_169 += tmp_170;
   }
   tmp_168 += tmp_169;
   result += (0.5) * tmp_168;
   std::complex<double> tmp_171;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_172;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_172 += B0(p,MHm(gI1),MHm(gI2))*Conj(CpUhhconjHmHm(gO2,gI1,
            gI2))*CpUhhconjHmHm(gO1,gI1,gI2);
      }
      tmp_171 += tmp_172;
   }
   result += tmp_171;
   std::complex<double> tmp_173;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_174;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_174 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))*
            CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_173 += tmp_174;
   }
   result += tmp_173;
   std::complex<double> tmp_175;
   std::complex<double> tmp_176;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_177;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_177 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))*
            CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_176 += tmp_177;
   }
   tmp_175 += tmp_176;
   result += (0.5) * tmp_175;
   std::complex<double> tmp_178;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_179;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_179 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*CpUhhbarChaChaPL
            (gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*CpUhhbarChaChaPR(
            gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_178 += tmp_179;
   }
   result += tmp_178;
   std::complex<double> tmp_180;
   std::complex<double> tmp_181;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_182;
      std::complex<double> tmp_183;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_183 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2,
            gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_182 += tmp_183;
      tmp_181 += (MCha(gI1)) * tmp_182;
   }
   tmp_180 += tmp_181;
   result += (-2) * tmp_180;
   std::complex<double> tmp_184;
   std::complex<double> tmp_185;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_186;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_186 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(gO1
            ,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,gI1,
            gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_185 += tmp_186;
   }
   tmp_184 += tmp_185;
   result += (3) * tmp_184;
   std::complex<double> tmp_187;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_188;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_188 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(gO1
            ,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,gI1,
            gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_187 += tmp_188;
   }
   result += tmp_187;
   std::complex<double> tmp_189;
   std::complex<double> tmp_190;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_191;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_191 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(gO1
            ,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,gI1,
            gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_190 += tmp_191;
   }
   tmp_189 += tmp_190;
   result += (3) * tmp_189;
   std::complex<double> tmp_192;
   std::complex<double> tmp_193;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_194;
      std::complex<double> tmp_195;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_195 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1,
            gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*
            CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_194 += tmp_195;
      tmp_193 += (MFd(gI1)) * tmp_194;
   }
   tmp_192 += tmp_193;
   result += (-6) * tmp_192;
   std::complex<double> tmp_196;
   std::complex<double> tmp_197;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_198;
      std::complex<double> tmp_199;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_199 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1,
            gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*
            CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_198 += tmp_199;
      tmp_197 += (MFe(gI1)) * tmp_198;
   }
   tmp_196 += tmp_197;
   result += (-2) * tmp_196;
   std::complex<double> tmp_200;
   std::complex<double> tmp_201;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_202;
      std::complex<double> tmp_203;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_203 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1,
            gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*
            CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_202 += tmp_203;
      tmp_201 += (MFu(gI1)) * tmp_202;
   }
   tmp_200 += tmp_201;
   result += (-6) * tmp_200;
   std::complex<double> tmp_204;
   std::complex<double> tmp_205;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_206;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_206 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_205 += tmp_206;
   }
   tmp_204 += tmp_205;
   result += (0.5) * tmp_204;
   std::complex<double> tmp_207;
   std::complex<double> tmp_208;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_209;
      std::complex<double> tmp_210;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_210 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,gI1
            ,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2))*
            CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_209 += tmp_210;
      tmp_208 += (MChi(gI1)) * tmp_209;
   }
   tmp_207 += tmp_208;
   result += (-1) * tmp_207;
   std::complex<double> tmp_211;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_211 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_211;
   std::complex<double> tmp_212;
   std::complex<double> tmp_213;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_213 += Conj(CpUhhconjVWmHm(gO2,gI2))*CpUhhconjVWmHm(gO1,gI2)*F0(p,
         MHm(gI2),MVWm);
   }
   tmp_212 += tmp_213;
   result += (2) * tmp_212;
   result += 4*CpUhhUhhconjVWmVWm(gO1,gO2)*(A0(MVWm) - 0.5*Sqr(MVWm));
   result += 2*CpUhhUhhVZVZ(gO1,gO2)*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   std::complex<double> tmp_214;
   std::complex<double> tmp_215;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_215 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_214 += tmp_215;
   result += (-0.5) * tmp_214;
   std::complex<double> tmp_216;
   std::complex<double> tmp_217;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_217 += A0(MHm(gI1))*CpUAhUAhconjHmHm(gO1,gO2,gI1,gI1);
   }
   tmp_216 += tmp_217;
   result += (-1) * tmp_216;
   std::complex<double> tmp_218;
   std::complex<double> tmp_219;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_219 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_218 += tmp_219;
   result += (-0.5) * tmp_218;
   std::complex<double> tmp_220;
   std::complex<double> tmp_221;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_222;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_222 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))*
            CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_221 += tmp_222;
   }
   tmp_220 += tmp_221;
   result += (0.5) * tmp_220;
   std::complex<double> tmp_223;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_224;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_224 += B0(p,MHm(gI1),MHm(gI2))*Conj(CpUAhconjHmHm(gO2,gI1,
            gI2))*CpUAhconjHmHm(gO1,gI1,gI2);
      }
      tmp_223 += tmp_224;
   }
   result += tmp_223;
   std::complex<double> tmp_225;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_226;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_226 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))*
            CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_225 += tmp_226;
   }
   result += tmp_225;
   std::complex<double> tmp_227;
   std::complex<double> tmp_228;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_229;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_229 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))*
            CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_228 += tmp_229;
   }
   tmp_227 += tmp_228;
   result += (0.5) * tmp_227;
   std::complex<double> tmp_230;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_231;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_231 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*CpUAhbarChaChaPL
            (gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*CpUAhbarChaChaPR(
            gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_230 += tmp_231;
   }
   result += tmp_230;
   std::complex<double> tmp_232;
   std::complex<double> tmp_233;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_234;
      std::complex<double> tmp_235;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_235 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2,
            gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_234 += tmp_235;
      tmp_233 += (MCha(gI1)) * tmp_234;
   }
   tmp_232 += tmp_233;
   result += (-2) * tmp_232;
   std::complex<double> tmp_236;
   std::complex<double> tmp_237;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_238;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_238 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(gO1
            ,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,gI1,
            gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_237 += tmp_238;
   }
   tmp_236 += tmp_237;
   result += (3) * tmp_236;
   std::complex<double> tmp_239;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_240;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_240 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(gO1
            ,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,gI1,
            gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_239 += tmp_240;
   }
   result += tmp_239;
   std::complex<double> tmp_241;
   std::complex<double> tmp_242;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_243;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_243 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(gO1
            ,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,gI1,
            gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_242 += tmp_243;
   }
   tmp_241 += tmp_242;
   result += (3) * tmp_241;
   std::complex<double> tmp_244;
   std::complex<double> tmp_245;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_246;
      std::complex<double> tmp_247;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_247 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1,
            gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*
            CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_246 += tmp_247;
      tmp_245 += (MFd(gI1)) * tmp_246;
   }
   tmp_244 += tmp_245;
   result += (-6) * tmp_244;
   std::complex<double> tmp_248;
   std::complex<double> tmp_249;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_250;
      std::complex<double> tmp_251;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_251 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1,
            gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*
            CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_250 += tmp_251;
      tmp_249 += (MFe(gI1)) * tmp_250;
   }
   tmp_248 += tmp_249;
   result += (-2) * tmp_248;
   std::complex<double> tmp_252;
   std::complex<double> tmp_253;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_254;
      std::complex<double> tmp_255;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_255 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1,
            gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*
            CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_254 += tmp_255;
      tmp_253 += (MFu(gI1)) * tmp_254;
   }
   tmp_252 += tmp_253;
   result += (-6) * tmp_252;
   std::complex<double> tmp_256;
   std::complex<double> tmp_257;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_258;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_258 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_257 += tmp_258;
   }
   tmp_256 += tmp_257;
   result += (0.5) * tmp_256;
   std::complex<double> tmp_259;
   std::complex<double> tmp_260;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_261;
      std::complex<double> tmp_262;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_262 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,gI1
            ,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2))*
            CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_261 += tmp_262;
      tmp_260 += (MChi(gI1)) * tmp_261;
   }
   tmp_259 += tmp_260;
   result += (-1) * tmp_259;
   std::complex<double> tmp_263;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_263 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_263;
   std::complex<double> tmp_264;
   std::complex<double> tmp_265;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_265 += Conj(CpUAhconjVWmHm(gO2,gI2))*CpUAhconjVWmHm(gO1,gI2)*F0(p,
         MHm(gI2),MVWm);
   }
   tmp_264 += tmp_265;
   result += (2) * tmp_264;
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
   std::complex<double> tmp_266;
   std::complex<double> tmp_267;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_267 += A0(MAh(gI1))*CpUHmconjUHmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_266 += tmp_267;
   result += (-0.5) * tmp_266;
   std::complex<double> tmp_268;
   std::complex<double> tmp_269;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_269 += A0(MHm(gI1))*CpUHmconjUHmconjHmHm(gO1,gO2,gI1,gI1);
   }
   tmp_268 += tmp_269;
   result += (-1) * tmp_268;
   std::complex<double> tmp_270;
   std::complex<double> tmp_271;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_271 += A0(Mhh(gI1))*CpUHmconjUHmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_270 += tmp_271;
   result += (-0.5) * tmp_270;
   std::complex<double> tmp_272;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_273;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_273 += B0(p,MHm(gI1),MAh(gI2))*Conj(CpconjUHmHmAh(gO2,gI1,
            gI2))*CpconjUHmHmAh(gO1,gI1,gI2);
      }
      tmp_272 += tmp_273;
   }
   result += tmp_272;
   std::complex<double> tmp_274;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_275;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_275 += B0(p,MHm(gI1),Mhh(gI2))*Conj(CpconjUHmHmhh(gO2,gI1,
            gI2))*CpconjUHmHmhh(gO1,gI1,gI2);
      }
      tmp_274 += tmp_275;
   }
   result += tmp_274;
   std::complex<double> tmp_276;
   std::complex<double> tmp_277;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_278;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_278 += (Conj(CpconjUHmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHmbarFuFdPR(gO2,gI1,gI2)
            )*CpconjUHmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_277 += tmp_278;
   }
   tmp_276 += tmp_277;
   result += (3) * tmp_276;
   std::complex<double> tmp_279;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_280;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_280 += (Conj(CpconjUHmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHmbarFvFePR(gO2,gI1,gI2)
            )*CpconjUHmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_279 += tmp_280;
   }
   result += tmp_279;
   std::complex<double> tmp_281;
   std::complex<double> tmp_282;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_283;
      std::complex<double> tmp_284;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_284 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHmbarFuFdPR(gO2,
            gI1,gI2))*CpconjUHmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHmbarFuFdPL(
            gO2,gI1,gI2))*CpconjUHmbarFuFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_283 += tmp_284;
      tmp_282 += (MFu(gI1)) * tmp_283;
   }
   tmp_281 += tmp_282;
   result += (-6) * tmp_281;
   std::complex<double> tmp_285;
   std::complex<double> tmp_286;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_287;
      std::complex<double> tmp_288;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_288 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHmbarFvFePR(gO2,
            gI1,gI2))*CpconjUHmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHmbarFvFePL(
            gO2,gI1,gI2))*CpconjUHmbarFvFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_287 += tmp_288;
      tmp_286 += (MFv(gI1)) * tmp_287;
   }
   tmp_285 += tmp_286;
   result += (-2) * tmp_285;
   std::complex<double> tmp_289;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_290;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_290 += (Conj(CpconjUHmChiChaPL(gO2,gI1,gI2))*
            CpconjUHmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHmChiChaPR(gO2,gI1,gI2))*
            CpconjUHmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_289 += tmp_290;
   }
   result += tmp_289;
   std::complex<double> tmp_291;
   std::complex<double> tmp_292;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_293;
      std::complex<double> tmp_294;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_294 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHmChiChaPR(gO2
            ,gI1,gI2))*CpconjUHmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHmChiChaPL(gO2
            ,gI1,gI2))*CpconjUHmChiChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_293 += tmp_294;
      tmp_292 += (MChi(gI1)) * tmp_293;
   }
   tmp_291 += tmp_292;
   result += (-2) * tmp_291;
   std::complex<double> tmp_295;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_295 += Conj(CpconjUHmVWmAh(gO2,gI2))*CpconjUHmVWmAh(gO1,gI2)*F0(p,
         MAh(gI2),MVWm);
   }
   result += tmp_295;
   std::complex<double> tmp_296;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_296 += Conj(CpconjUHmVWmhh(gO2,gI2))*CpconjUHmVWmhh(gO1,gI2)*F0(p,
         Mhh(gI2),MVWm);
   }
   result += tmp_296;
   std::complex<double> tmp_297;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_297 += Conj(CpconjUHmVPHm(gO2,gI2))*CpconjUHmVPHm(gO1,gI2)*F0(p,
         MHm(gI2),0);
   }
   result += tmp_297;
   std::complex<double> tmp_298;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_298 += Conj(CpconjUHmVZHm(gO2,gI2))*CpconjUHmVZHm(gO1,gI2)*F0(p,
         MHm(gI2),MVZ);
   }
   result += tmp_298;
   result += 4*CpUHmconjUHmconjVWmVWm(gO1,gO2)*(A0(MVWm) - 0.5*Sqr(MVWm));
   result += 2*CpUHmconjUHmVZVZ(gO1,gO2)*(A0(MVZ) - 0.5*Sqr(MVZ));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_299;
   std::complex<double> tmp_300;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_300 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_299 += tmp_300;
   result += (0.5) * tmp_299;
   std::complex<double> tmp_301;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_301 += A0(MHm(gI1))*CpVZVZconjHmHm(gI1,gI1);
   }
   result += tmp_301;
   std::complex<double> tmp_302;
   std::complex<double> tmp_303;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_303 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_302 += tmp_303;
   result += (0.5) * tmp_302;
   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_306;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_306 += AbsSqr(CpVZconjHmHm(gI1,gI2))*B00(p,MHm(gI1),MHm(gI2)
            );
      }
      tmp_305 += tmp_306;
   }
   tmp_304 += tmp_305;
   result += (-4) * tmp_304;
   std::complex<double> tmp_307;
   std::complex<double> tmp_308;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_309;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_309 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_308 += tmp_309;
   }
   tmp_307 += tmp_308;
   result += (-4) * tmp_307;
   std::complex<double> tmp_310;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_311;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_311 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_311 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_310 += tmp_311;
   }
   result += tmp_310;
   std::complex<double> tmp_312;
   std::complex<double> tmp_313;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_314;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_314 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_314 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_313 += tmp_314;
   }
   tmp_312 += tmp_313;
   result += (3) * tmp_312;
   std::complex<double> tmp_315;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_316;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_316 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_316 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_315 += tmp_316;
   }
   result += tmp_315;
   std::complex<double> tmp_317;
   std::complex<double> tmp_318;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_319;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_319 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_319 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_318 += tmp_319;
   }
   tmp_317 += tmp_318;
   result += (3) * tmp_317;
   std::complex<double> tmp_320;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_321;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_321 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_321 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_320 += tmp_321;
   }
   result += tmp_320;
   std::complex<double> tmp_322;
   std::complex<double> tmp_323;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_324;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_324 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR(
            gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_324 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_323 += tmp_324;
   }
   tmp_322 += tmp_323;
   result += (0.5) * tmp_322;
   std::complex<double> tmp_325;
   std::complex<double> tmp_326;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_326 += AbsSqr(CpVZconjVWmHm(gI2))*B0(p,MVWm,MHm(gI2));
   }
   tmp_325 += tmp_326;
   result += (2) * tmp_325;
   std::complex<double> tmp_327;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_327 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_327;
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
   std::complex<double> tmp_328;
   std::complex<double> tmp_329;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_329 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_328 += tmp_329;
   result += (0.5) * tmp_328;
   std::complex<double> tmp_330;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_330 += A0(MHm(gI1))*CpVWmconjVWmconjHmHm(gI1,gI1);
   }
   result += tmp_330;
   std::complex<double> tmp_331;
   std::complex<double> tmp_332;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_332 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_331 += tmp_332;
   result += (0.5) * tmp_331;
   std::complex<double> tmp_333;
   std::complex<double> tmp_334;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_335;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_335 += AbsSqr(CpconjVWmHmAh(gI1,gI2))*B00(p,MAh(gI2),MHm(gI1
            ));
      }
      tmp_334 += tmp_335;
   }
   tmp_333 += tmp_334;
   result += (-4) * tmp_333;
   std::complex<double> tmp_336;
   std::complex<double> tmp_337;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_338;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_338 += AbsSqr(CpconjVWmHmhh(gI1,gI2))*B00(p,Mhh(gI2),MHm(gI1
            ));
      }
      tmp_337 += tmp_338;
   }
   tmp_336 += tmp_337;
   result += (-4) * tmp_336;
   std::complex<double> tmp_339;
   std::complex<double> tmp_340;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_341;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_341 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_341 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_340 += tmp_341;
   }
   tmp_339 += tmp_340;
   result += (3) * tmp_339;
   std::complex<double> tmp_342;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_343;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_343 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_343 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_342 += tmp_343;
   }
   result += tmp_342;
   std::complex<double> tmp_344;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_345;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_345 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_345 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_344 += tmp_345;
   }
   result += tmp_344;
   std::complex<double> tmp_346;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_346 += AbsSqr(CpconjVWmVPHm(gI2))*B0(p,0,MHm(gI2));
   }
   result += tmp_346;
   std::complex<double> tmp_347;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_347 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_347;
   std::complex<double> tmp_348;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_348 += AbsSqr(CpconjVWmVZHm(gI2))*B0(p,MVZ,MHm(gI2));
   }
   result += tmp_348;
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

   std::complex<double> tmp_349;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_350;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_350 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_349 += tmp_350;
   }
   result += tmp_349;
   std::complex<double> tmp_351;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_352;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_352 += B0(p,MFu(gI2),MHm(gI1))*Conj(CpbarUFdHmFuPL(gO2,gI1,
            gI2))*CpbarUFdHmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_351 += tmp_352;
   }
   result += tmp_351;
   std::complex<double> tmp_353;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_354;
      std::complex<double> tmp_355;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_355 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_354 += tmp_355;
      tmp_353 += (MFd(gI1)) * tmp_354;
   }
   result += tmp_353;
   std::complex<double> tmp_356;
   std::complex<double> tmp_357;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_357 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_356 += tmp_357;
   result += (-5.333333333333333) * tmp_356;
   std::complex<double> tmp_358;
   std::complex<double> tmp_359;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_359 += (-0.5 + B0(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_358 += tmp_359;
   result += (-4) * tmp_358;
   std::complex<double> tmp_360;
   std::complex<double> tmp_361;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_361 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_360 += tmp_361;
   result += (-4) * tmp_360;
   std::complex<double> tmp_362;
   std::complex<double> tmp_363;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_363 += (-0.5 + B0(p,MFu(gI2),MVWm))*Conj(CpbarUFdVWmFuPR(gO2,gI2))
         *CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_362 += tmp_363;
   result += (-4) * tmp_362;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_364;
   std::complex<double> tmp_365;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_366;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_366 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_365 += tmp_366;
   }
   tmp_364 += tmp_365;
   result += (-0.5) * tmp_364;
   std::complex<double> tmp_367;
   std::complex<double> tmp_368;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_369;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_369 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarUFdHmFuPR(gO2,gI1,
            gI2))*CpbarUFdHmFuPR(gO1,gI1,gI2);
      }
      tmp_368 += tmp_369;
   }
   tmp_367 += tmp_368;
   result += (-0.5) * tmp_367;
   std::complex<double> tmp_370;
   std::complex<double> tmp_371;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_372;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_372 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_371 += tmp_372;
   }
   tmp_370 += tmp_371;
   result += (-0.5) * tmp_370;
   std::complex<double> tmp_373;
   std::complex<double> tmp_374;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_374 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_373 += tmp_374;
   result += (-1.3333333333333333) * tmp_373;
   std::complex<double> tmp_375;
   std::complex<double> tmp_376;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_376 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_375 += tmp_376;
   result += (-1) * tmp_375;
   std::complex<double> tmp_377;
   std::complex<double> tmp_378;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_378 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_377 += tmp_378;
   result += (-1) * tmp_377;
   std::complex<double> tmp_379;
   std::complex<double> tmp_380;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_380 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_379 += tmp_380;
   result += (-1) * tmp_379;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_381;
   std::complex<double> tmp_382;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_383;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_383 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_382 += tmp_383;
   }
   tmp_381 += tmp_382;
   result += (-0.5) * tmp_381;
   std::complex<double> tmp_384;
   std::complex<double> tmp_385;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_386;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_386 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarUFdHmFuPL(gO2,gI1,
            gI2))*CpbarUFdHmFuPL(gO1,gI1,gI2);
      }
      tmp_385 += tmp_386;
   }
   tmp_384 += tmp_385;
   result += (-0.5) * tmp_384;
   std::complex<double> tmp_387;
   std::complex<double> tmp_388;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_389;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_389 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_388 += tmp_389;
   }
   tmp_387 += tmp_388;
   result += (-0.5) * tmp_387;
   std::complex<double> tmp_390;
   std::complex<double> tmp_391;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_391 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_390 += tmp_391;
   result += (-1.3333333333333333) * tmp_390;
   std::complex<double> tmp_392;
   std::complex<double> tmp_393;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_393 += (0.5 + B1(p,MFd(gI2),0))*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_392 += tmp_393;
   result += (-1) * tmp_392;
   std::complex<double> tmp_394;
   std::complex<double> tmp_395;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_395 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_394 += tmp_395;
   result += (-1) * tmp_394;
   std::complex<double> tmp_396;
   std::complex<double> tmp_397;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_397 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_396 += tmp_397;
   result += (-1) * tmp_396;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_398;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_399;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_399 += B0(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_398 += tmp_399;
   }
   result += tmp_398;
   std::complex<double> tmp_400;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_401;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_401 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_400 += tmp_401;
   }
   result += tmp_400;
   std::complex<double> tmp_402;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_403;
      std::complex<double> tmp_404;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_404 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_403 += tmp_404;
      tmp_402 += (MFu(gI1)) * tmp_403;
   }
   result += tmp_402;
   std::complex<double> tmp_405;
   std::complex<double> tmp_406;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_406 += (-0.5 + B0(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_405 += tmp_406;
   result += (-4) * tmp_405;
   std::complex<double> tmp_407;
   std::complex<double> tmp_408;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_408 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_407 += tmp_408;
   result += (-5.333333333333333) * tmp_407;
   std::complex<double> tmp_409;
   std::complex<double> tmp_410;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_410 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_409 += tmp_410;
   result += (-4) * tmp_409;
   std::complex<double> tmp_411;
   std::complex<double> tmp_412;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_412 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_411 += tmp_412;
   result += (-4) * tmp_411;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_413;
   std::complex<double> tmp_414;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_415;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_415 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPR(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2);
      }
      tmp_414 += tmp_415;
   }
   tmp_413 += tmp_414;
   result += (-0.5) * tmp_413;
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_418;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_418 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_417 += tmp_418;
   }
   tmp_416 += tmp_417;
   result += (-0.5) * tmp_416;
   std::complex<double> tmp_419;
   std::complex<double> tmp_420;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_421;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_421 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_420 += tmp_421;
   }
   tmp_419 += tmp_420;
   result += (-0.5) * tmp_419;
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_423 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPL(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_422 += tmp_423;
   result += (-1) * tmp_422;
   std::complex<double> tmp_424;
   std::complex<double> tmp_425;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_425 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_424 += tmp_425;
   result += (-1.3333333333333333) * tmp_424;
   std::complex<double> tmp_426;
   std::complex<double> tmp_427;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_427 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_426 += tmp_427;
   result += (-1) * tmp_426;
   std::complex<double> tmp_428;
   std::complex<double> tmp_429;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_429 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_428 += tmp_429;
   result += (-1) * tmp_428;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_430;
   std::complex<double> tmp_431;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_432;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_432 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPL(gO1,gI1,gI2);
      }
      tmp_431 += tmp_432;
   }
   tmp_430 += tmp_431;
   result += (-0.5) * tmp_430;
   std::complex<double> tmp_433;
   std::complex<double> tmp_434;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_435;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_435 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_434 += tmp_435;
   }
   tmp_433 += tmp_434;
   result += (-0.5) * tmp_433;
   std::complex<double> tmp_436;
   std::complex<double> tmp_437;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_438;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_438 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_437 += tmp_438;
   }
   tmp_436 += tmp_437;
   result += (-0.5) * tmp_436;
   std::complex<double> tmp_439;
   std::complex<double> tmp_440;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_440 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_439 += tmp_440;
   result += (-1) * tmp_439;
   std::complex<double> tmp_441;
   std::complex<double> tmp_442;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_442 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_441 += tmp_442;
   result += (-1.3333333333333333) * tmp_441;
   std::complex<double> tmp_443;
   std::complex<double> tmp_444;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_444 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_443 += tmp_444;
   result += (-1) * tmp_443;
   std::complex<double> tmp_445;
   std::complex<double> tmp_446;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_446 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_445 += tmp_446;
   result += (-1) * tmp_445;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_447;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_448;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_448 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_447 += tmp_448;
   }
   result += tmp_447;
   std::complex<double> tmp_449;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_450;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_450 += B0(p,MFv(gI2),MHm(gI1))*Conj(CpbarUFeHmFvPL(gO2,gI1,
            gI2))*CpbarUFeHmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_449 += tmp_450;
   }
   result += tmp_449;
   std::complex<double> tmp_451;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_452;
      std::complex<double> tmp_453;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_453 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_452 += tmp_453;
      tmp_451 += (MFe(gI1)) * tmp_452;
   }
   result += tmp_451;
   std::complex<double> tmp_454;
   std::complex<double> tmp_455;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_455 += (-0.5 + B0(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_454 += tmp_455;
   result += (-4) * tmp_454;
   std::complex<double> tmp_456;
   std::complex<double> tmp_457;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_457 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_456 += tmp_457;
   result += (-4) * tmp_456;
   std::complex<double> tmp_458;
   std::complex<double> tmp_459;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_459 += (-0.5 + B0(p,MFv(gI2),MVWm))*Conj(CpbarUFeVWmFvPR(gO2,gI2))
         *CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_458 += tmp_459;
   result += (-4) * tmp_458;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_460;
   std::complex<double> tmp_461;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_462;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_462 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_461 += tmp_462;
   }
   tmp_460 += tmp_461;
   result += (-0.5) * tmp_460;
   std::complex<double> tmp_463;
   std::complex<double> tmp_464;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_465;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_465 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarUFeHmFvPR(gO2,gI1,
            gI2))*CpbarUFeHmFvPR(gO1,gI1,gI2);
      }
      tmp_464 += tmp_465;
   }
   tmp_463 += tmp_464;
   result += (-0.5) * tmp_463;
   std::complex<double> tmp_466;
   std::complex<double> tmp_467;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_468;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_468 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_467 += tmp_468;
   }
   tmp_466 += tmp_467;
   result += (-0.5) * tmp_466;
   std::complex<double> tmp_469;
   std::complex<double> tmp_470;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_470 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_469 += tmp_470;
   result += (-1) * tmp_469;
   std::complex<double> tmp_471;
   std::complex<double> tmp_472;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_472 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_471 += tmp_472;
   result += (-1) * tmp_471;
   std::complex<double> tmp_473;
   std::complex<double> tmp_474;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_474 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_473 += tmp_474;
   result += (-1) * tmp_473;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_475;
   std::complex<double> tmp_476;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_477;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_477 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_476 += tmp_477;
   }
   tmp_475 += tmp_476;
   result += (-0.5) * tmp_475;
   std::complex<double> tmp_478;
   std::complex<double> tmp_479;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_480;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_480 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarUFeHmFvPL(gO2,gI1,
            gI2))*CpbarUFeHmFvPL(gO1,gI1,gI2);
      }
      tmp_479 += tmp_480;
   }
   tmp_478 += tmp_479;
   result += (-0.5) * tmp_478;
   std::complex<double> tmp_481;
   std::complex<double> tmp_482;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_483;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_483 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_482 += tmp_483;
   }
   tmp_481 += tmp_482;
   result += (-0.5) * tmp_481;
   std::complex<double> tmp_484;
   std::complex<double> tmp_485;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_485 += (0.5 + B1(p,MFe(gI2),0))*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_484 += tmp_485;
   result += (-1) * tmp_484;
   std::complex<double> tmp_486;
   std::complex<double> tmp_487;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_487 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_486 += tmp_487;
   result += (-1) * tmp_486;
   std::complex<double> tmp_488;
   std::complex<double> tmp_489;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_489 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_488 += tmp_489;
   result += (-1) * tmp_488;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_490;
   std::complex<double> tmp_491;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_492;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_492 += B0(p,MCha(gI2),MHm(gI1))*Conj(CpUChiconjHmChaPL(gO2,
            gI1,gI2))*CpUChiconjHmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_491 += tmp_492;
   }
   tmp_490 += tmp_491;
   result += (2) * tmp_490;
   std::complex<double> tmp_493;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_494;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_494 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_493 += tmp_494;
   }
   result += tmp_493;
   std::complex<double> tmp_495;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_496;
      std::complex<double> tmp_497;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_497 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_496 += tmp_497;
      tmp_495 += (MChi(gI1)) * tmp_496;
   }
   result += tmp_495;
   std::complex<double> tmp_498;
   std::complex<double> tmp_499;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_499 += (-0.5 + B0(p,MCha(gI2),MVWm))*Conj(CpUChiconjVWmChaPR(gO2,
         gI2))*CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_498 += tmp_499;
   result += (-8) * tmp_498;
   std::complex<double> tmp_500;
   std::complex<double> tmp_501;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_501 += (-0.5 + B0(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_500 += tmp_501;
   result += (-4) * tmp_500;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_502;
   std::complex<double> tmp_503;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_504;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_504 += B1(p,MCha(gI2),MHm(gI1))*Conj(CpUChiconjHmChaPR(gO2,
            gI1,gI2))*CpUChiconjHmChaPR(gO1,gI1,gI2);
      }
      tmp_503 += tmp_504;
   }
   tmp_502 += tmp_503;
   result += (-1) * tmp_502;
   std::complex<double> tmp_505;
   std::complex<double> tmp_506;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_507;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_507 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_506 += tmp_507;
   }
   tmp_505 += tmp_506;
   result += (-0.5) * tmp_505;
   std::complex<double> tmp_508;
   std::complex<double> tmp_509;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_510;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_510 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_509 += tmp_510;
   }
   tmp_508 += tmp_509;
   result += (-0.5) * tmp_508;
   std::complex<double> tmp_511;
   std::complex<double> tmp_512;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_512 += (0.5 + B1(p,MCha(gI2),MVWm))*Conj(CpUChiconjVWmChaPL(gO2,
         gI2))*CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_511 += tmp_512;
   result += (-2) * tmp_511;
   std::complex<double> tmp_513;
   std::complex<double> tmp_514;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_514 += (0.5 + B1(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_513 += tmp_514;
   result += (-1) * tmp_513;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_515;
   std::complex<double> tmp_516;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_517;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_517 += B1(p,MCha(gI2),MHm(gI1))*Conj(CpUChiconjHmChaPL(gO2,
            gI1,gI2))*CpUChiconjHmChaPL(gO1,gI1,gI2);
      }
      tmp_516 += tmp_517;
   }
   tmp_515 += tmp_516;
   result += (-1) * tmp_515;
   std::complex<double> tmp_518;
   std::complex<double> tmp_519;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_520;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_520 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_519 += tmp_520;
   }
   tmp_518 += tmp_519;
   result += (-0.5) * tmp_518;
   std::complex<double> tmp_521;
   std::complex<double> tmp_522;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_523;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_523 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_522 += tmp_523;
   }
   tmp_521 += tmp_522;
   result += (-0.5) * tmp_521;
   std::complex<double> tmp_524;
   std::complex<double> tmp_525;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_525 += (0.5 + B1(p,MCha(gI2),MVWm))*Conj(CpUChiconjVWmChaPR(gO2,
         gI2))*CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_524 += tmp_525;
   result += (-2) * tmp_524;
   std::complex<double> tmp_526;
   std::complex<double> tmp_527;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_527 += (0.5 + B1(p,MChi(gI2),MVZ))*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_526 += tmp_527;
   result += (-1) * tmp_526;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_528;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_529;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_529 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_528 += tmp_529;
   }
   result += tmp_528;
   std::complex<double> tmp_530;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_531;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_531 += B0(p,MChi(gI2),MHm(gI1))*Conj(CpbarUChaHmChiPL(gO2,
            gI1,gI2))*CpbarUChaHmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_530 += tmp_531;
   }
   result += tmp_530;
   std::complex<double> tmp_532;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_533;
      std::complex<double> tmp_534;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_534 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_533 += tmp_534;
      tmp_532 += (MCha(gI1)) * tmp_533;
   }
   result += tmp_532;
   std::complex<double> tmp_535;
   std::complex<double> tmp_536;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_536 += (-0.5 + B0(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_535 += tmp_536;
   result += (-4) * tmp_535;
   std::complex<double> tmp_537;
   std::complex<double> tmp_538;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_538 += (-0.5 + B0(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPR(gO2,gI2)
         )*CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_537 += tmp_538;
   result += (-4) * tmp_537;
   std::complex<double> tmp_539;
   std::complex<double> tmp_540;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_540 += (-0.5 + B0(p,MChi(gI2),MVWm))*Conj(CpbarUChaVWmChiPR(gO2,
         gI2))*CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_539 += tmp_540;
   result += (-4) * tmp_539;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_541;
   std::complex<double> tmp_542;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_543;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_543 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_542 += tmp_543;
   }
   tmp_541 += tmp_542;
   result += (-0.5) * tmp_541;
   std::complex<double> tmp_544;
   std::complex<double> tmp_545;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_546;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_546 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_545 += tmp_546;
   }
   tmp_544 += tmp_545;
   result += (-0.5) * tmp_544;
   std::complex<double> tmp_547;
   std::complex<double> tmp_548;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_549;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_549 += B1(p,MChi(gI2),MHm(gI1))*Conj(CpbarUChaHmChiPR(gO2,
            gI1,gI2))*CpbarUChaHmChiPR(gO1,gI1,gI2);
      }
      tmp_548 += tmp_549;
   }
   tmp_547 += tmp_548;
   result += (-0.5) * tmp_547;
   std::complex<double> tmp_550;
   std::complex<double> tmp_551;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_551 += (0.5 + B1(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_550 += tmp_551;
   result += (-1) * tmp_550;
   std::complex<double> tmp_552;
   std::complex<double> tmp_553;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_553 += (0.5 + B1(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPL(gO2,gI2))
         *CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_552 += tmp_553;
   result += (-1) * tmp_552;
   std::complex<double> tmp_554;
   std::complex<double> tmp_555;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_555 += (0.5 + B1(p,MChi(gI2),MVWm))*Conj(CpbarUChaVWmChiPL(gO2,gI2
         ))*CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_554 += tmp_555;
   result += (-1) * tmp_554;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_556;
   std::complex<double> tmp_557;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_558;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_558 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_557 += tmp_558;
   }
   tmp_556 += tmp_557;
   result += (-0.5) * tmp_556;
   std::complex<double> tmp_559;
   std::complex<double> tmp_560;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_561;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_561 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_560 += tmp_561;
   }
   tmp_559 += tmp_560;
   result += (-0.5) * tmp_559;
   std::complex<double> tmp_562;
   std::complex<double> tmp_563;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_564;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_564 += B1(p,MChi(gI2),MHm(gI1))*Conj(CpbarUChaHmChiPL(gO2,
            gI1,gI2))*CpbarUChaHmChiPL(gO1,gI1,gI2);
      }
      tmp_563 += tmp_564;
   }
   tmp_562 += tmp_563;
   result += (-0.5) * tmp_562;
   std::complex<double> tmp_565;
   std::complex<double> tmp_566;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_566 += (0.5 + B1(p,MCha(gI2),0))*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_565 += tmp_566;
   result += (-1) * tmp_565;
   std::complex<double> tmp_567;
   std::complex<double> tmp_568;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_568 += (0.5 + B1(p,MCha(gI2),MVZ))*Conj(CpbarUChaVZChaPR(gO2,gI2))
         *CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_567 += tmp_568;
   result += (-1) * tmp_567;
   std::complex<double> tmp_569;
   std::complex<double> tmp_570;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_570 += (0.5 + B1(p,MChi(gI2),MVWm))*Conj(CpbarUChaVWmChiPR(gO2,gI2
         ))*CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_569 += tmp_570;
   result += (-1) * tmp_569;

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

   std::complex<double> tmp_571;
   std::complex<double> tmp_572;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_572 += AbsSqr(CpVZhhAh(gI1,1))*B00(p,MAh(1),Mhh(gI1));
   }
   tmp_571 += tmp_572;
   result += (-4) * tmp_571;
   std::complex<double> tmp_573;
   std::complex<double> tmp_574;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_574 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_573 += tmp_574;
   result += (0.5) * tmp_573;
   std::complex<double> tmp_575;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_576;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_576 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_576 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_575 += tmp_576;
   }
   result += tmp_575;
   std::complex<double> tmp_577;
   std::complex<double> tmp_578;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_579;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_579 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR(
            gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_579 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_578 += tmp_579;
   }
   tmp_577 += tmp_578;
   result += (0.5) * tmp_577;
   std::complex<double> tmp_580;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_580 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_580;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_581;
   std::complex<double> tmp_582;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_582 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_581 += tmp_582;
   result += (0.5) * tmp_581;
   std::complex<double> tmp_583;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_584;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_584 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_584 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_583 += tmp_584;
   }
   result += tmp_583;
   std::complex<double> tmp_585;
   std::complex<double> tmp_586;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_586 += AbsSqr(CpconjVWmHmhh(1,gI2))*B00(p,Mhh(gI2),MHm(1));
   }
   tmp_585 += tmp_586;
   result += (-4) * tmp_585;
   std::complex<double> tmp_587;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_587 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_587;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_588;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_589;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_589 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_588 += tmp_589;
   }
   result += tmp_588;
   std::complex<double> tmp_590;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_591;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_591 += B0(p,MFu(gI2),MHm(gI1))*Conj(CpbarFdHmFuPL(gO2,gI1,
            gI2))*CpbarFdHmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_590 += tmp_591;
   }
   result += tmp_590;
   std::complex<double> tmp_592;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_593;
      std::complex<double> tmp_594;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_594 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_593 += tmp_594;
      tmp_592 += (MFd(gI1)) * tmp_593;
   }
   result += tmp_592;
   std::complex<double> tmp_595;
   std::complex<double> tmp_596;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_596 += (-0.5 + B0(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_595 += tmp_596;
   result += (-4) * tmp_595;
   std::complex<double> tmp_597;
   std::complex<double> tmp_598;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_598 += (-0.5 + B0(p,MFu(gI2),MVWm))*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_597 += tmp_598;
   result += (-4) * tmp_597;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_599;
   std::complex<double> tmp_600;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_601;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_601 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_600 += tmp_601;
   }
   tmp_599 += tmp_600;
   result += (-0.5) * tmp_599;
   std::complex<double> tmp_602;
   std::complex<double> tmp_603;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_604;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_604 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarFdHmFuPR(gO2,gI1,
            gI2))*CpbarFdHmFuPR(gO1,gI1,gI2);
      }
      tmp_603 += tmp_604;
   }
   tmp_602 += tmp_603;
   result += (-0.5) * tmp_602;
   std::complex<double> tmp_605;
   std::complex<double> tmp_606;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_607;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_607 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_606 += tmp_607;
   }
   tmp_605 += tmp_606;
   result += (-0.5) * tmp_605;
   std::complex<double> tmp_608;
   std::complex<double> tmp_609;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_609 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_608 += tmp_609;
   result += (-1) * tmp_608;
   std::complex<double> tmp_610;
   std::complex<double> tmp_611;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_611 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_610 += tmp_611;
   result += (-1) * tmp_610;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_612;
   std::complex<double> tmp_613;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_614;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_614 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_613 += tmp_614;
   }
   tmp_612 += tmp_613;
   result += (-0.5) * tmp_612;
   std::complex<double> tmp_615;
   std::complex<double> tmp_616;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_617;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_617 += B1(p,MFu(gI2),MHm(gI1))*Conj(CpbarFdHmFuPL(gO2,gI1,
            gI2))*CpbarFdHmFuPL(gO1,gI1,gI2);
      }
      tmp_616 += tmp_617;
   }
   tmp_615 += tmp_616;
   result += (-0.5) * tmp_615;
   std::complex<double> tmp_618;
   std::complex<double> tmp_619;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_620;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_620 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_619 += tmp_620;
   }
   tmp_618 += tmp_619;
   result += (-0.5) * tmp_618;
   std::complex<double> tmp_621;
   std::complex<double> tmp_622;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_622 += (0.5 + B1(p,MFu(gI2),MVWm))*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_621 += tmp_622;
   result += (-1) * tmp_621;
   std::complex<double> tmp_623;
   std::complex<double> tmp_624;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_624 += (0.5 + B1(p,MFd(gI2),MVZ))*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_623 += tmp_624;
   result += (-1) * tmp_623;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_625;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_626;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_626 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_625 += tmp_626;
   }
   result += tmp_625;
   std::complex<double> tmp_627;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_628;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_628 += B0(p,MFv(gI2),MHm(gI1))*Conj(CpbarFeHmFvPL(gO2,gI1,
            gI2))*CpbarFeHmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_627 += tmp_628;
   }
   result += tmp_627;
   std::complex<double> tmp_629;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_630;
      std::complex<double> tmp_631;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_631 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_630 += tmp_631;
      tmp_629 += (MFe(gI1)) * tmp_630;
   }
   result += tmp_629;
   std::complex<double> tmp_632;
   std::complex<double> tmp_633;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_633 += (-0.5 + B0(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_632 += tmp_633;
   result += (-4) * tmp_632;
   std::complex<double> tmp_634;
   std::complex<double> tmp_635;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_635 += (-0.5 + B0(p,MFv(gI2),MVWm))*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_634 += tmp_635;
   result += (-4) * tmp_634;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_636;
   std::complex<double> tmp_637;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_638;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_638 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_637 += tmp_638;
   }
   tmp_636 += tmp_637;
   result += (-0.5) * tmp_636;
   std::complex<double> tmp_639;
   std::complex<double> tmp_640;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_641;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_641 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarFeHmFvPR(gO2,gI1,
            gI2))*CpbarFeHmFvPR(gO1,gI1,gI2);
      }
      tmp_640 += tmp_641;
   }
   tmp_639 += tmp_640;
   result += (-0.5) * tmp_639;
   std::complex<double> tmp_642;
   std::complex<double> tmp_643;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_644;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_644 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_643 += tmp_644;
   }
   tmp_642 += tmp_643;
   result += (-0.5) * tmp_642;
   std::complex<double> tmp_645;
   std::complex<double> tmp_646;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_646 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_645 += tmp_646;
   result += (-1) * tmp_645;
   std::complex<double> tmp_647;
   std::complex<double> tmp_648;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_648 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_647 += tmp_648;
   result += (-1) * tmp_647;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_649;
   std::complex<double> tmp_650;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_651;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_651 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_650 += tmp_651;
   }
   tmp_649 += tmp_650;
   result += (-0.5) * tmp_649;
   std::complex<double> tmp_652;
   std::complex<double> tmp_653;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_654;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_654 += B1(p,MFv(gI2),MHm(gI1))*Conj(CpbarFeHmFvPL(gO2,gI1,
            gI2))*CpbarFeHmFvPL(gO1,gI1,gI2);
      }
      tmp_653 += tmp_654;
   }
   tmp_652 += tmp_653;
   result += (-0.5) * tmp_652;
   std::complex<double> tmp_655;
   std::complex<double> tmp_656;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_657;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_657 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_656 += tmp_657;
   }
   tmp_655 += tmp_656;
   result += (-0.5) * tmp_655;
   std::complex<double> tmp_658;
   std::complex<double> tmp_659;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_659 += (0.5 + B1(p,MFv(gI2),MVWm))*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_658 += tmp_659;
   result += (-1) * tmp_658;
   std::complex<double> tmp_660;
   std::complex<double> tmp_661;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_661 += (0.5 + B1(p,MFe(gI2),MVZ))*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_660 += tmp_661;
   result += (-1) * tmp_660;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_662;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_663;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_663 += B0(p,MFd(gI2),MHm(gI1))*Conj(CpbarFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarFuconjHmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_662 += tmp_663;
   }
   result += tmp_662;
   std::complex<double> tmp_664;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_665;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_665 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_664 += tmp_665;
   }
   result += tmp_664;
   std::complex<double> tmp_666;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_667;
      std::complex<double> tmp_668;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_668 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_667 += tmp_668;
      tmp_666 += (MFu(gI1)) * tmp_667;
   }
   result += tmp_666;
   std::complex<double> tmp_669;
   std::complex<double> tmp_670;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_670 += (-0.5 + B0(p,MFd(gI2),MVWm))*Conj(CpbarFuconjVWmFdPR(gO2,
         gI2))*CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_669 += tmp_670;
   result += (-4) * tmp_669;
   std::complex<double> tmp_671;
   std::complex<double> tmp_672;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_672 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_671 += tmp_672;
   result += (-4) * tmp_671;
   std::complex<double> tmp_673;
   std::complex<double> tmp_674;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_674 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_673 += tmp_674;
   result += (-4) * tmp_673;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_675;
   std::complex<double> tmp_676;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_677;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_677 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarFuconjHmFdPR(gO2,
            gI1,gI2))*CpbarFuconjHmFdPR(gO1,gI1,gI2);
      }
      tmp_676 += tmp_677;
   }
   tmp_675 += tmp_676;
   result += (-0.5) * tmp_675;
   std::complex<double> tmp_678;
   std::complex<double> tmp_679;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_680;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_680 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_679 += tmp_680;
   }
   tmp_678 += tmp_679;
   result += (-0.5) * tmp_678;
   std::complex<double> tmp_681;
   std::complex<double> tmp_682;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_683;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_683 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_682 += tmp_683;
   }
   tmp_681 += tmp_682;
   result += (-0.5) * tmp_681;
   std::complex<double> tmp_684;
   std::complex<double> tmp_685;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_685 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarFuconjVWmFdPL(gO2,gI2
         ))*CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_684 += tmp_685;
   result += (-1) * tmp_684;
   std::complex<double> tmp_686;
   std::complex<double> tmp_687;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_687 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_686 += tmp_687;
   result += (-1) * tmp_686;
   std::complex<double> tmp_688;
   std::complex<double> tmp_689;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_689 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_688 += tmp_689;
   result += (-1) * tmp_688;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_690;
   std::complex<double> tmp_691;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_692;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_692 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarFuconjHmFdPL(gO1,gI1,gI2);
      }
      tmp_691 += tmp_692;
   }
   tmp_690 += tmp_691;
   result += (-0.5) * tmp_690;
   std::complex<double> tmp_693;
   std::complex<double> tmp_694;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_695;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_695 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_694 += tmp_695;
   }
   tmp_693 += tmp_694;
   result += (-0.5) * tmp_693;
   std::complex<double> tmp_696;
   std::complex<double> tmp_697;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_698;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_698 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_697 += tmp_698;
   }
   tmp_696 += tmp_697;
   result += (-0.5) * tmp_696;
   std::complex<double> tmp_699;
   std::complex<double> tmp_700;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_700 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarFuconjVWmFdPR(gO2,gI2
         ))*CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_699 += tmp_700;
   result += (-1) * tmp_699;
   std::complex<double> tmp_701;
   std::complex<double> tmp_702;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_702 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_701 += tmp_702;
   result += (-1) * tmp_701;
   std::complex<double> tmp_703;
   std::complex<double> tmp_704;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_704 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_703 += tmp_704;
   result += (-1) * tmp_703;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_705;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_706;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_706 += B0(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_705 += tmp_706;
   }
   result += tmp_705;
   std::complex<double> tmp_707;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_708;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_708 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_707 += tmp_708;
   }
   result += tmp_707;
   std::complex<double> tmp_709;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_710;
      std::complex<double> tmp_711;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_711 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_710 += tmp_711;
      tmp_709 += (MFu(gI1)) * tmp_710;
   }
   result += tmp_709;
   std::complex<double> tmp_712;
   std::complex<double> tmp_713;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_713 += (-0.5 + B0(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_712 += tmp_713;
   result += (-4) * tmp_712;
   std::complex<double> tmp_714;
   std::complex<double> tmp_715;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_715 += (-0.5 + B0(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_714 += tmp_715;
   result += (-4) * tmp_714;
   std::complex<double> tmp_716;
   std::complex<double> tmp_717;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_717 += (-0.5 + B0(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_716 += tmp_717;
   result += (-4) * tmp_716;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_718;
   std::complex<double> tmp_719;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_720;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_720 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPR(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPR(gO1,gI1,gI2);
      }
      tmp_719 += tmp_720;
   }
   tmp_718 += tmp_719;
   result += (-0.5) * tmp_718;
   std::complex<double> tmp_721;
   std::complex<double> tmp_722;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_723;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_723 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_722 += tmp_723;
   }
   tmp_721 += tmp_722;
   result += (-0.5) * tmp_721;
   std::complex<double> tmp_724;
   std::complex<double> tmp_725;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_726;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_726 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_725 += tmp_726;
   }
   tmp_724 += tmp_725;
   result += (-0.5) * tmp_724;
   std::complex<double> tmp_727;
   std::complex<double> tmp_728;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_728 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPL(gO2,
         gI2))*CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_727 += tmp_728;
   result += (-1) * tmp_727;
   std::complex<double> tmp_729;
   std::complex<double> tmp_730;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_730 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_729 += tmp_730;
   result += (-1) * tmp_729;
   std::complex<double> tmp_731;
   std::complex<double> tmp_732;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_732 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_731 += tmp_732;
   result += (-1) * tmp_731;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_733;
   std::complex<double> tmp_734;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_735;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_735 += B1(p,MFd(gI2),MHm(gI1))*Conj(CpbarUFuconjHmFdPL(gO2,
            gI1,gI2))*CpbarUFuconjHmFdPL(gO1,gI1,gI2);
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
         tmp_738 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
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
         tmp_741 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_740 += tmp_741;
   }
   tmp_739 += tmp_740;
   result += (-0.5) * tmp_739;
   std::complex<double> tmp_742;
   std::complex<double> tmp_743;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_743 += (0.5 + B1(p,MFd(gI2),MVWm))*Conj(CpbarUFuconjVWmFdPR(gO2,
         gI2))*CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_742 += tmp_743;
   result += (-1) * tmp_742;
   std::complex<double> tmp_744;
   std::complex<double> tmp_745;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_745 += (0.5 + B1(p,MFu(gI2),0))*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_744 += tmp_745;
   result += (-1) * tmp_744;
   std::complex<double> tmp_746;
   std::complex<double> tmp_747;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_747 += (0.5 + B1(p,MFu(gI2),MVZ))*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_746 += tmp_747;
   result += (-1) * tmp_746;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh(unsigned gO1) const
{
   std::complex<double> result;

   result += A0(MVWm)*CpUhhbargWmCgWmC(gO1);
   result += A0(MVWm)*CpUhhbargWmgWm(gO1);
   result += A0(MVZ)*CpUhhbargZgZ(gO1);
   std::complex<double> tmp_748;
   std::complex<double> tmp_749;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_749 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_748 += tmp_749;
   result += (-0.5) * tmp_748;
   std::complex<double> tmp_750;
   std::complex<double> tmp_751;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_751 += A0(MHm(gI1))*CpUhhconjHmHm(gO1,gI1,gI1);
   }
   tmp_750 += tmp_751;
   result += (-1) * tmp_750;
   std::complex<double> tmp_752;
   std::complex<double> tmp_753;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_753 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_752 += tmp_753;
   result += (-0.5) * tmp_752;
   std::complex<double> tmp_754;
   std::complex<double> tmp_755;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_755 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_754 += tmp_755;
   result += (2) * tmp_754;
   std::complex<double> tmp_756;
   std::complex<double> tmp_757;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_757 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR(
         gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_756 += tmp_757;
   result += (6) * tmp_756;
   std::complex<double> tmp_758;
   std::complex<double> tmp_759;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_759 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR(
         gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_758 += tmp_759;
   result += (2) * tmp_758;
   std::complex<double> tmp_760;
   std::complex<double> tmp_761;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_761 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR(
         gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_760 += tmp_761;
   result += (6) * tmp_760;
   std::complex<double> tmp_762;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_762 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_762;
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
   if (!force_output && problems.is_tachyon(VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVZ));
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_Mhh_pole()
{
   if (!force_output && problems.is_tachyon(hh))
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
   if (!force_output && problems.is_tachyon(Ah))
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
   if (!force_output && problems.is_tachyon(Hm))
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
   if (!force_output && problems.is_tachyon(VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(Sqr(MVWm));
   const double p = MVWm;
   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VWm);

   PHYSICAL(MVWm) = AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVWm_pole(double p)
{
   if (!force_output && problems.is_tachyon(VWm))
      return 0.;

   const double self_energy = Re(self_energy_VWm(p));
   const double mass_sqr = Sqr(MVWm) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_tachyon(VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = Sqr(MVZ) - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VZ);

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
      problems.flag_tachyon(VZ);
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
      problems.flag_tachyon(VWm);
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
