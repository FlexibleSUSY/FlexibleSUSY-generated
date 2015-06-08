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

// File generated at Mon 8 Jun 2015 17:54:38

/**
 * @file UMSSM_mass_eigenstates.cpp
 * @brief implementation of the UMSSM model class
 *
 * Contains the definition of the UMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 8 Jun 2015 17:54:38 with FlexibleSUSY
 * 1.1.1 (git commit: v1.1.1) and SARAH 4.5.6 .
 */

#include "UMSSM_mass_eigenstates.hpp"
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
#include "pv.hpp"
#include "functors.hpp"

#include "sfermions.hpp"
#include "nmssm_twoloophiggs.h"


#include <cmath>
#include <iostream>
#include <algorithm>

#ifdef ENABLE_THREADS
#include <thread>
#endif

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace UMSSM_info;

#define CLASSNAME UMSSM_mass_eigenstates

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter
#define MODELPARAMETER(parameter) model->get_##parameter()

#define HIGGS_2LOOP_CORRECTION_AT_AS     two_loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     two_loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     two_loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU two_loop_corrections.higgs_atau_atau
#define TOP_2LOOP_CORRECTION_QCD         two_loop_corrections.top_qcd

#ifdef ENABLE_THREADS
   std::mutex CLASSNAME::mtx_fortran;
   #define LOCK_MUTEX() mtx_fortran.lock()
   #define UNLOCK_MUTEX() mtx_fortran.unlock()
#else
   #define LOCK_MUTEX()
   #define UNLOCK_MUTEX()
#endif

CLASSNAME::UMSSM_mass_eigenstates(const UMSSM_input_parameters& input_)
   : UMSSM_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(UMSSM_info::particle_names)
   , two_loop_corrections()
#ifdef ENABLE_THREADS
   , thread_exception()
#endif
   , MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MVP(0), MVZ(0),
      MVZp(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>
      ::Zero()), MSu(Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,
      1>::Zero()), Mhh(Eigen::Array<double,3,1>::Zero()), MAh(Eigen::Array<double,
      3,1>::Zero()), MHpm(Eigen::Array<double,2,1>::Zero()), MChi(Eigen::Array<
      double,6,1>::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(
      Eigen::Array<double,3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()),
      MFu(Eigen::Array<double,3,1>::Zero()), MVWm(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,6,6>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,
      2>::Zero()), UP(Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~UMSSM_mass_eigenstates()
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

void CLASSNAME::set_number_of_ewsb_iterations(std::size_t iterations)
{
   number_of_ewsb_iterations = iterations;
}

void CLASSNAME::set_number_of_mass_iterations(std::size_t iterations)
{
   number_of_mass_iterations = iterations;
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

void CLASSNAME::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double CLASSNAME::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
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

const Problems<UMSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<UMSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
{
   return problems;
}

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @param tadpole array of tadpole
 */
void CLASSNAME::tadpole_equations(double tadpole[number_of_ewsb_equations]) const
{
   tadpole[0] = get_ewsb_eq_hh_1();
   tadpole[1] = get_ewsb_eq_hh_2();
   tadpole[2] = get_ewsb_eq_hh_3();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0));
      tadpole[1] -= Re(tadpole_hh(1));
      tadpole[2] -= Re(tadpole_hh(2));

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] -= two_loop_tadpole[0];
         tadpole[1] -= two_loop_tadpole[1];
         tadpole[2] -= two_loop_tadpole[2];

      }
   }
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
   UMSSM_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_mHd2(gsl_vector_get(x, 0));
   model->set_mHu2(gsl_vector_get(x, 1));
   model->set_ms2(gsl_vector_get(x, 2));


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double tadpole[number_of_ewsb_equations] = { 0. };

   model->tadpole_equations(tadpole);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   return is_finite<number_of_ewsb_equations>(tadpole) ? GSL_SUCCESS : GSL_EDOM;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int CLASSNAME::solve_ewsb_iteratively()
{
   EWSB_args params = {this, ewsb_loop_order};

   EWSB_solver* solvers[] = {
      new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(CLASSNAME::ewsb_step, &params, number_of_ewsb_iterations, ewsb_iteration_precision),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_hybrids),
      new Root_finder<number_of_ewsb_equations>(CLASSNAME::tadpole_equations, &params, number_of_ewsb_iterations, ewsb_iteration_precision, gsl_multiroot_fsolver_broyden)
   };

   const std::size_t number_of_solvers = sizeof(solvers)/sizeof(*solvers);
   double x_init[number_of_ewsb_equations] = { 0. };
   ewsb_initial_guess(x_init);

#ifdef ENABLE_VERBOSE
   std::cout << "Solving EWSB equations ...\n"
      "\tInitial guess: x_init =";
   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      std::cout << ' ' << x_init[i];
   std::cout << '\n';
#endif

   int status;
   for (std::size_t i = 0; i < number_of_solvers; ++i) {
      VERBOSE_MSG("\tStarting EWSB iteration using solver " << i);
      status = solve_ewsb_iteratively_with(solvers[i], x_init);
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

   for_each(solvers, solvers + number_of_solvers, Delete_object());

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
   const double x_init[number_of_ewsb_equations]
)
{
   const int status = solver->solve(x_init);

   mHd2 = solver->get_solution(0);
   mHu2 = solver->get_solution(1);
   ms2 = solver->get_solution(2);


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

   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_ms2 = ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) - 3*Power(vd,3)*
      Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*Power(vd,3)*Sqr(gp)*Sqr(QHd) - 20*vd*
      AbsSqr(Lambdax)*Sqr(vS) - 20*QHd*Qs*vd*Sqr(gp)*Sqr(vS) - 20*vd*AbsSqr(
      Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) - 20*QHd*QHu*
      vd*Sqr(gp)*Sqr(vu) + 14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) - 3*Power(vu,3)*
      Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 20*Power(vu,3)*Sqr(gp)*Sqr(QHu) - 20*vu*
      AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*
      QHd*QHu*vu*Sqr(gp)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) - 20*QHu*Qs*vu*
      Sqr(gp)*Sqr(vS) + 14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(1.4142135623730951*vd*vu*Conj(TLambdax) - 2*Power(vS,3)*Sqr(
      gp)*Sqr(Qs) - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 2*QHd*Qs*vS*Sqr(gp)*Sqr(vd) - 2
      *vS*AbsSqr(Lambdax)*Sqr(vu) - 2*QHu*Qs*vS*Sqr(gp)*Sqr(vu) +
      1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      ms2 = old_ms2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_via_soft_higgs_masses()
{
   int error = 0;

   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const double new_mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) -
      3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*Power(vd,3)*Sqr(gp)*Sqr(
      QHd) - 20*vd*AbsSqr(Lambdax)*Sqr(vS) - 20*QHd*Qs*vd*Sqr(gp)*Sqr(vS) - 20*vd*
      AbsSqr(Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) - 20*
      QHd*QHu*vd*Sqr(gp)*Sqr(vu) + 14.142135623730951*vS*vu*TLambdax))/vd);
   const double new_mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) -
      3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 20*Power(vu,3)*Sqr(gp)*Sqr(
      QHu) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*
      Sqr(vd) - 20*QHd*QHu*vu*Sqr(gp)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) - 20
      *QHu*Qs*vu*Sqr(gp)*Sqr(vS) + 14.142135623730951*vd*vS*TLambdax))/vu);
   const double new_ms2 = Re((0.25*(1.4142135623730951*vd*vu*Conj(TLambdax) - 2
      *Power(vS,3)*Sqr(gp)*Sqr(Qs) - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 2*QHd*Qs*vS*
      Sqr(gp)*Sqr(vd) - 2*vS*AbsSqr(Lambdax)*Sqr(vu) - 2*QHu*Qs*vS*Sqr(gp)*Sqr(vu)
      + 1.4142135623730951*vd*vu*TLambdax))/vS);

   if (IsFinite(new_mHd2))
      mHd2 = new_mHd2;
   else
      error = 1;

   if (IsFinite(new_mHu2))
      mHu2 = new_mHu2;
   else
      error = 1;

   if (IsFinite(new_ms2))
      ms2 = new_ms2;
   else
      error = 1;


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

void CLASSNAME::ewsb_initial_guess(double x_init[number_of_ewsb_equations])
{
   x_init[0] = mHd2;
   x_init[1] = mHu2;
   x_init[2] = ms2;

}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * @param ewsb_parameters new EWSB output parameters.  \a
 * ewsb_parameters is only modified if all new parameters are finite.
 *
 * @return GSL_SUCCESS if new EWSB output parameters are finite,
 * GSL_EDOM otherwise.
 */
int CLASSNAME::ewsb_step(double ewsb_parameters[number_of_ewsb_equations]) const
{
   int error;
   double tadpole[number_of_ewsb_equations] = { 0. };

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh(0));
      tadpole[1] += Re(tadpole_hh(1));
      tadpole[2] += Re(tadpole_hh(2));

      if (ewsb_loop_order > 1) {
         double two_loop_tadpole[3];
         tadpole_hh_2loop(two_loop_tadpole);
         tadpole[0] += two_loop_tadpole[0];
         tadpole[1] += two_loop_tadpole[1];
         tadpole[2] += two_loop_tadpole[2];

      }
   }

   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   double mHd2;
   double mHu2;
   double ms2;

   mHd2 = Re((0.025*(14.142135623730951*vS*vu*Conj(TLambdax) + 40*tadpole[0] -
      3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*Power(vd,3)*Sqr(gp)*Sqr(
      QHd) - 20*vd*AbsSqr(Lambdax)*Sqr(vS) - 20*QHd*Qs*vd*Sqr(gp)*Sqr(vS) - 20*vd*
      AbsSqr(Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu) - 20*
      QHd*QHu*vd*Sqr(gp)*Sqr(vu) + 14.142135623730951*vS*vu*TLambdax))/vd);
   mHu2 = Re((0.025*(14.142135623730951*vd*vS*Conj(TLambdax) + 40*tadpole[1] -
      3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 20*Power(vu,3)*Sqr(gp)*Sqr(
      QHu) - 20*vu*AbsSqr(Lambdax)*Sqr(vd) + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*
      Sqr(vd) - 20*QHd*QHu*vu*Sqr(gp)*Sqr(vd) - 20*vu*AbsSqr(Lambdax)*Sqr(vS) - 20
      *QHu*Qs*vu*Sqr(gp)*Sqr(vS) + 14.142135623730951*vd*vS*TLambdax))/vu);
   ms2 = Re((0.25*(1.4142135623730951*vd*vu*Conj(TLambdax) + 4*tadpole[2] - 2*
      Power(vS,3)*Sqr(gp)*Sqr(Qs) - 2*vS*AbsSqr(Lambdax)*Sqr(vd) - 2*QHd*Qs*vS*Sqr
      (gp)*Sqr(vd) - 2*vS*AbsSqr(Lambdax)*Sqr(vu) - 2*QHu*Qs*vS*Sqr(gp)*Sqr(vu) +
      1.4142135623730951*vd*vu*TLambdax))/vS);

   const bool is_finite = IsFinite(mHd2) && IsFinite(mHu2) && IsFinite(ms2);


   if (is_finite) {
      error = GSL_SUCCESS;
      ewsb_parameters[0] = mHd2;
      ewsb_parameters[1] = mHu2;
      ewsb_parameters[2] = ms2;

   } else {
      error = GSL_EDOM;
   }

   return error;
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
   UMSSM_mass_eigenstates* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double mHd2 = gsl_vector_get(x, 0);
   const double mHu2 = gsl_vector_get(x, 1);
   const double ms2 = gsl_vector_get(x, 2);

   model->set_mHd2(mHd2);
   model->set_mHu2(mHu2);
   model->set_ms2(ms2);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { mHd2, mHu2, ms2 };

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
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
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MVZp = " << MVZp << '\n';
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
   const auto old_mHd2 = mHd2;
   const auto old_mHu2 = mHu2;
   const auto old_ms2 = ms2;

   solve_ewsb_tree_level_via_soft_higgs_masses();

   calculate_MVG();
   calculate_MVP();
   calculate_MVZ();
   calculate_MVZp();
   calculate_MVWm();
   calculate_MGlu();
   calculate_MFv();
   calculate_MSd();
   calculate_MSv();
   calculate_MSu();
   calculate_MSe();
   calculate_Mhh();
   calculate_MAh();
   calculate_MHpm();
   calculate_MChi();
   calculate_MCha();
   calculate_MFe();
   calculate_MFd();
   calculate_MFu();

   mHd2 = old_mHd2;
   mHu2 = old_mHu2;
   ms2 = old_ms2;

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
   thread_exception = 0;

   std::thread thread_MAh(Thread(this, &CLASSNAME::calculate_MAh_pole));
   std::thread thread_MCha(Thread(this, &CLASSNAME::calculate_MCha_pole));
   std::thread thread_MChi(Thread(this, &CLASSNAME::calculate_MChi_pole));
   std::thread thread_MGlu(Thread(this, &CLASSNAME::calculate_MGlu_pole));
   std::thread thread_Mhh(Thread(this, &CLASSNAME::calculate_Mhh_pole));
   std::thread thread_MHpm(Thread(this, &CLASSNAME::calculate_MHpm_pole));
   std::thread thread_MSd(Thread(this, &CLASSNAME::calculate_MSd_pole));
   std::thread thread_MSe(Thread(this, &CLASSNAME::calculate_MSe_pole));
   std::thread thread_MSu(Thread(this, &CLASSNAME::calculate_MSu_pole));
   std::thread thread_MSv(Thread(this, &CLASSNAME::calculate_MSv_pole));
   std::thread thread_MVZp(Thread(this, &CLASSNAME::calculate_MVZp_pole));

   if (calculate_sm_pole_masses) {
      std::thread thread_MVG(Thread(this, &CLASSNAME::calculate_MVG_pole));
      std::thread thread_MFv(Thread(this, &CLASSNAME::calculate_MFv_pole));
      std::thread thread_MVP(Thread(this, &CLASSNAME::calculate_MVP_pole));
      std::thread thread_MVZ(Thread(this, &CLASSNAME::calculate_MVZ_pole));
      std::thread thread_MFe(Thread(this, &CLASSNAME::calculate_MFe_pole));
      std::thread thread_MFd(Thread(this, &CLASSNAME::calculate_MFd_pole));
      std::thread thread_MFu(Thread(this, &CLASSNAME::calculate_MFu_pole));
      std::thread thread_MVWm(Thread(this, &CLASSNAME::calculate_MVWm_pole));
      thread_MVG.join();
      thread_MFv.join();
      thread_MVP.join();
      thread_MVZ.join();
      thread_MFe.join();
      thread_MFd.join();
      thread_MFu.join();
      thread_MVWm.join();
   }

   thread_MAh.join();
   thread_MCha.join();
   thread_MChi.join();
   thread_MGlu.join();
   thread_Mhh.join();
   thread_MHpm.join();
   thread_MSd.join();
   thread_MSe.join();
   thread_MSu.join();
   thread_MSv.join();
   thread_MVZp.join();


   if (thread_exception != 0)
      std::rethrow_exception(thread_exception);
#else
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
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(MVZp) = MVZp;
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

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associuated
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
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;
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
   MVWm = 0.;

   PhaseGlu = std::complex<double>(1.,0.);

}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   UMSSM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
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


Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_ChargedHiggs;
   Eigen::Array<double,1,1> MHpm_goldstone;

   MHpm_goldstone(0) = MVWm;

   remove_if_equal(MHpm, MHpm_goldstone, MHpm_ChargedHiggs);

   return MHpm_ChargedHiggs;
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_PseudoscalarHiggs;
   Eigen::Array<double,2,1> MAh_goldstone;

   MAh_goldstone(0) = MVZ;
   MAh_goldstone(1) = MVZp;

   remove_if_equal(MAh, MAh_goldstone, MAh_PseudoscalarHiggs);

   return MAh_PseudoscalarHiggs;
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
   MVG = calculate_singlet_mass(mass_matrix_VG);
}

double CLASSNAME::get_mass_matrix_Glu() const
{
   const double mass_matrix_Glu = Re(MassG);

   return mass_matrix_Glu;
}

void CLASSNAME::calculate_MGlu()
{
   const auto mass_matrix_Glu = get_mass_matrix_Glu();
   MGlu = calculate_singlet_mass(mass_matrix_Glu, PhaseGlu);
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

double CLASSNAME::get_mass_matrix_VP() const
{
   const double mass_matrix_VP = Re(0);

   return mass_matrix_VP;
}

void CLASSNAME::calculate_MVP()
{
   const auto mass_matrix_VP = get_mass_matrix_VP();
   MVP = calculate_singlet_mass(mass_matrix_VP);
}

double CLASSNAME::get_mass_matrix_VZ() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const double mass_matrix_VZ = Re(0.25*(2*g2*gp*Cos(ThetaW())*Sin(2*
      ThetaWp())*(QHd*Sqr(vd) - QHu*Sqr(vu)) + 1.5491933384829668*g1*gp*Sin(
      ThetaW())*Sin(2*ThetaWp())*(QHd*Sqr(vd) - QHu*Sqr(vu)) +
      0.7745966692414834*g1*Sin(ThetaW())*(2*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*(Sqr(vd) + Sqr(vu))*Sqr(Cos(ThetaWp(
      ))) + Sqr(g2)*(Sqr(vd) + Sqr(vu))*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 4*Sqr(gp)*(Sqr(QHd)*Sqr(vd) + Sqr(Qs)*Sqr(vS) + Sqr(QHu)*Sqr(vu))*Sqr(
      Sin(ThetaWp()))));

   return mass_matrix_VZ;
}

void CLASSNAME::calculate_MVZ()
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const auto mass_matrix_VZ = get_mass_matrix_VZ();
   MVZ = calculate_singlet_mass(mass_matrix_VZ);

   if (MVZ < 0.)
      problems.flag_tachyon(UMSSM_info::VZ);

   MVZ = AbsSqrt(MVZ);
}

double CLASSNAME::get_mass_matrix_VZp() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const double mass_matrix_VZp = Re(0.25*(-2*gp*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sin(2*ThetaWp())*(QHd*Sqr(vd) - QHu*
      Sqr(vu)) + 4*Sqr(gp)*(Sqr(QHd)*Sqr(vd) + Sqr(Qs)*Sqr(vS) + Sqr(QHu)*Sqr(
      vu))*Sqr(Cos(ThetaWp())) + (Sqr(vd) + Sqr(vu))*Sqr(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sqr(Sin(ThetaWp()))));

   return mass_matrix_VZp;
}

void CLASSNAME::calculate_MVZp()
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   const auto mass_matrix_VZp = get_mass_matrix_VZp();
   MVZp = calculate_singlet_mass(mass_matrix_VZp);

   if (MVZp < 0.)
      problems.flag_tachyon(UMSSM_info::VZp);

   MVZp = AbsSqrt(MVZp);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qd = LOCALINPUT(Qd);

   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)
      ) + AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) +
      Conj(Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) +
      Conj(Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) - 0.5*vS*vu
      *Conj(Yd(0,0))*Lambdax;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) - 0.5*vS*vu
      *Conj(Yd(1,0))*Lambdax;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) - 0.5*vS*vu
      *Conj(Yd(2,0))*Lambdax;
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) +
      Conj(Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) - 0.5*vS*vu
      *Conj(Yd(0,1))*Lambdax;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) - 0.5*vS*vu
      *Conj(Yd(1,1))*Lambdax;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) - 0.5*vS*vu
      *Conj(Yd(2,1))*Lambdax;
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.025*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) - 0.5*vS*vu
      *Conj(Yd(0,2))*Lambdax;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) - 0.5*vS*vu
      *Conj(Yd(1,2))*Lambdax;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) - 0.5*vS*vu
      *Conj(Yd(2,2))*Lambdax;
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)
      ) + AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*
      Sqr(vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*
      Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) +
      Conj(Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) +
      Conj(Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*
      Sqr(vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*
      Sqr(gp)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) +
      Conj(Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.5*Qd*QHd*Sqr(gp)*
      Sqr(vd) + 0.5*Qd*Qs*Sqr(gp)*Sqr(vS) + 0.05*Sqr(g1)*Sqr(vu) + 0.5*Qd*QHu*
      Sqr(gp)*Sqr(vu);

   Hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Sd, eigenvalue_error > precision *
      Abs(MSd(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif

   if (MSd.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::Sd);

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075
      *Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075
      *Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075
      *Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);

   Hermitianize(mass_matrix_Sv);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Sv, eigenvalue_error > precision *
      Abs(MSv(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif

   if (MSv.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::Sv);

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

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) +
      Conj(Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) +
      Conj(Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) - 0.5*vd*vS
      *Conj(Yu(0,0))*Lambdax;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) - 0.5*vd*vS
      *Conj(Yu(1,0))*Lambdax;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) - 0.5*vd*vS
      *Conj(Yu(2,0))*Lambdax;
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) +
      Conj(Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) - 0.5*vd*vS
      *Conj(Yu(0,1))*Lambdax;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) - 0.5*vd*vS
      *Conj(Yu(1,1))*Lambdax;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) - 0.5*vd*vS
      *Conj(Yu(2,1))*Lambdax;
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Qq*Sqr(gp)*Sqr(vd) + 0.5*Qq*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) + 0.025*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Qq*Sqr(gp)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) - 0.5*vd*vS
      *Conj(Yu(0,2))*Lambdax;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) - 0.5*vd*vS
      *Conj(Yu(1,2))*Lambdax;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) - 0.5*vd*vS
      *Conj(Yu(2,2))*Lambdax;
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(
      gp)*Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(
      Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*
      Sqr(gp)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) +
      Conj(Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) +
      Conj(Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(
      gp)*Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(1,0)) + AbsSqr(
      Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*
      Sqr(gp)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) +
      Conj(Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qu*Sqr(
      gp)*Sqr(vd) + 0.5*Qs*Qu*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yu(2,0)) + AbsSqr(
      Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu) + 0.5*QHu*Qu*
      Sqr(gp)*Sqr(vu);

   Hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Su, eigenvalue_error > precision *
      Abs(MSu(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif

   if (MSu.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::Su);

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

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)
      ) + AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) +
      Conj(Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) +
      Conj(Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) - 0.5*vS*vu
      *Conj(Ye(0,0))*Lambdax;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) - 0.5*vS*vu
      *Conj(Ye(1,0))*Lambdax;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) - 0.5*vS*vu
      *Conj(Ye(2,0))*Lambdax;
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) +
      Conj(Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) - 0.5*vS*vu
      *Conj(Ye(0,1))*Lambdax;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) - 0.5*vS*vu
      *Conj(Ye(1,1))*Lambdax;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) - 0.5*vS*vu
      *Conj(Ye(2,1))*Lambdax;
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) - 0.5*vS*vu
      *Conj(Ye(0,2))*Lambdax;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) - 0.5*vS*vu
      *Conj(Ye(1,2))*Lambdax;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) - 0.5*vS*vu
      *Conj(Ye(2,2))*Lambdax;
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)
      ) + AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*
      Sqr(vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*
      Sqr(gp)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) +
      Conj(Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) +
      Conj(Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*
      Sqr(vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*
      Sqr(gp)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) +
      Conj(Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*Qe*QHd*Sqr(gp)*
      Sqr(vd) + 0.5*Qe*Qs*Sqr(gp)*Sqr(vS) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*Qe*QHu*
      Sqr(gp)*Sqr(vu);

   Hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Se, eigenvalue_error > precision *
      Abs(MSe(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif

   if (MSe.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::Se);

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr
      (vd) + 1.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*
      QHd*Qs*Sqr(gp)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu);
   mass_matrix_hh(0,1) = vd*vu*AbsSqr(Lambdax) - 0.35355339059327373*vS*
      Conj(TLambdax) - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) + QHd*QHu*vd*vu*
      Sqr(gp) - 0.35355339059327373*vS*TLambdax;
   mass_matrix_hh(0,2) = vd*vS*AbsSqr(Lambdax) - 0.35355339059327373*vu*
      Conj(TLambdax) + QHd*Qs*vd*vS*Sqr(gp) - 0.35355339059327373*vu*TLambdax;
   mass_matrix_hh(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(
      g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*
      AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vS) + 0.225*Sqr(g1)*Sqr(
      vu) + 0.375*Sqr(g2)*Sqr(vu) + 1.5*Sqr(gp)*Sqr(QHu)*Sqr(vu);
   mass_matrix_hh(1,2) = vS*vu*AbsSqr(Lambdax) - 0.35355339059327373*vd*
      Conj(TLambdax) + QHu*Qs*vS*vu*Sqr(gp) - 0.35355339059327373*vd*TLambdax;
   mass_matrix_hh(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*
      Sqr(gp)*Sqr(vd) + 1.5*Sqr(gp)*Sqr(Qs)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(
      vu) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vu);

   Symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::hh, eigenvalue_error > precision *
      Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif

   if (Mhh.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::hh);

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr
      (vd) + 0.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5*
      QHd*Qs*Sqr(gp)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(
      vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vd)*Sqr(Cos(
      ThetaWp())) + Sqr(gp)*Sqr(QHd)*Sqr(vd)*Sqr(Cos(ThetaWp())) + 0.25*Sqr(g2)
      *Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vd)*
      Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 0.3872983346207417*g1*g2*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(vd)*Sqr(Sin(ThetaWp())) + Sqr(gp)*Sqr(QHd)*
      Sqr(vd)*Sqr(Sin(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))*Sqr
      (Sin(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()));
   mass_matrix_Ah(0,1) = 0.35355339059327373*vS*Conj(TLambdax) -
      0.3872983346207417*g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW())*Sqr(Cos(
      ThetaWp())) + QHd*QHu*vd*vu*Sqr(gp)*Sqr(Cos(ThetaWp())) - 0.25*vd*vu*Sqr(
      g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) - 0.15*vd*vu*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW())) - 0.3872983346207417*g1*g2*vd*vu*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + QHd*QHu*vd*vu*Sqr(gp)*Sqr(
      Sin(ThetaWp())) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()
      )) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())) +
      0.35355339059327373*vS*TLambdax;
   mass_matrix_Ah(0,2) = 0.35355339059327373*vu*Conj(TLambdax) + QHd*Qs*
      vd*vS*Sqr(gp)*Sqr(Cos(ThetaWp())) + QHd*Qs*vd*vS*Sqr(gp)*Sqr(Sin(ThetaWp(
      ))) + 0.35355339059327373*vu*TLambdax;
   mass_matrix_Ah(1,1) = mHu2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(
      g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*
      AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(
      vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.5*Sqr(gp)*Sqr(QHu)*Sqr(vu) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu)*Sqr(Cos(
      ThetaWp())) + Sqr(gp)*Sqr(QHu)*Sqr(vu)*Sqr(Cos(ThetaWp())) + 0.25*Sqr(g2)
      *Sqr(vu)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vu)*
      Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 0.3872983346207417*g1*g2*Cos(
      ThetaW())*Sin(ThetaW())*Sqr(vu)*Sqr(Sin(ThetaWp())) + Sqr(gp)*Sqr(QHu)*
      Sqr(vu)*Sqr(Sin(ThetaWp())) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW()))*Sqr
      (Sin(ThetaWp())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()));
   mass_matrix_Ah(1,2) = 0.35355339059327373*vd*Conj(TLambdax) + QHu*Qs*
      vS*vu*Sqr(gp)*Sqr(Cos(ThetaWp())) + QHu*Qs*vS*vu*Sqr(gp)*Sqr(Sin(ThetaWp(
      ))) + 0.35355339059327373*vd*TLambdax;
   mass_matrix_Ah(2,2) = ms2 + 0.5*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*
      Sqr(gp)*Sqr(vd) + 0.5*Sqr(gp)*Sqr(Qs)*Sqr(vS) + 0.5*AbsSqr(Lambdax)*Sqr(
      vu) + 0.5*QHu*Qs*Sqr(gp)*Sqr(vu) + Sqr(gp)*Sqr(Qs)*Sqr(vS)*Sqr(Cos(
      ThetaWp())) + Sqr(gp)*Sqr(Qs)*Sqr(vS)*Sqr(Sin(ThetaWp()));

   Symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Ah, eigenvalue_error > precision *
      Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif

   if (MAh.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::Ah);

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*
      Sqr(vd) + 0.5*Sqr(gp)*Sqr(QHd)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) +
      0.5*QHd*Qs*Sqr(gp)*Sqr(vS) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu
      ) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vu);
   mass_matrix_Hpm(0,1) = -0.5*vd*vu*AbsSqr(Lambdax) + 0.7071067811865475
      *vS*Conj(TLambdax);
   mass_matrix_Hpm(1,1) = mHu2 - 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*
      Sqr(vd) + 0.5*QHd*QHu*Sqr(gp)*Sqr(vd) + 0.5*AbsSqr(Lambdax)*Sqr(vS) + 0.5
      *QHu*Qs*Sqr(gp)*Sqr(vS) + 0.075*Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu) +
      0.5*Sqr(gp)*Sqr(QHu)*Sqr(vu);

   Hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(UMSSM_info::Hpm, eigenvalue_error > precision *
      Abs(MHpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif

   if (MHpm.minCoeff() < 0.)
      problems.flag_tachyon(UMSSM_info::Hpm);

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
   problems.flag_bad_mass(UMSSM_info::Chi, eigenvalue_error > precision *
      Abs(MChi(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
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
   problems.flag_bad_mass(UMSSM_info::Cha, eigenvalue_error > precision *
      Abs(MCha(0)));
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
   problems.flag_bad_mass(UMSSM_info::Fe, eigenvalue_error > precision *
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
   problems.flag_bad_mass(UMSSM_info::Fd, eigenvalue_error > precision *
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
   problems.flag_bad_mass(UMSSM_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
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
   MVWm = calculate_singlet_mass(mass_matrix_VWm);

   if (MVWm < 0.)
      problems.flag_tachyon(UMSSM_info::VWm);

   MVWm = AbsSqrt(MVWm);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   double result = Re(mHd2*vd - 0.35355339059327373*vS*vu*Conj(TLambdax) +
      0.075*Power(vd,3)*Sqr(g1) + 0.125*Power(vd,3)*Sqr(g2) + 0.5*Power(vd,3)*Sqr(
      gp)*Sqr(QHd) + 0.5*vd*AbsSqr(Lambdax)*Sqr(vS) + 0.5*QHd*Qs*vd*Sqr(gp)*Sqr(vS
      ) + 0.5*vd*AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125*vd*Sqr
      (g2)*Sqr(vu) + 0.5*QHd*QHu*vd*Sqr(gp)*Sqr(vu) - 0.35355339059327373*vS*vu*
      TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   double result = Re(mHu2*vu - 0.35355339059327373*vd*vS*Conj(TLambdax) +
      0.075*Power(vu,3)*Sqr(g1) + 0.125*Power(vu,3)*Sqr(g2) + 0.5*Power(vu,3)*Sqr(
      gp)*Sqr(QHu) + 0.5*vu*AbsSqr(Lambdax)*Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) -
      0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*QHd*QHu*vu*Sqr(gp)*Sqr(vd) + 0.5*vu*AbsSqr(
      Lambdax)*Sqr(vS) + 0.5*QHu*Qs*vu*Sqr(gp)*Sqr(vS) - 0.35355339059327373*vd*vS
      *TLambdax);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);

   double result = Re(ms2*vS - 0.35355339059327373*vd*vu*Conj(TLambdax) + 0.5*
      Power(vS,3)*Sqr(gp)*Sqr(Qs) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vd) + 0.5*QHd*Qs*vS
      *Sqr(gp)*Sqr(vd) + 0.5*vS*AbsSqr(Lambdax)*Sqr(vu) + 0.5*QHu*Qs*vS*Sqr(gp)*
      Sqr(vu) - 0.35355339059327373*vd*vu*TLambdax);

   return result;
}



std::complex<double> CLASSNAME::CpUSdconjUSdVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_0;
   std::complex<double> tmp_1;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_0 += tmp_1;
   result += (0.13333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())
      )) * tmp_0;
   std::complex<double> tmp_2;
   std::complex<double> tmp_3;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_2 += tmp_3;
   result += (-1.0327955589886444*g1*gp*Qd*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_2;
   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_4 += tmp_5;
   result += (2*Sqr(gp)*Sqr(Qd)*Sqr(Sin(ThetaWp()))) * tmp_4;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(
         ThetaWp()))*Sqr(Sin(ThetaW()));
   }
   if (gO1 < 3) {
      result += -2*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += -0.5163977794943222*g1*gp*Qq*Cos(ThetaWp())*KroneckerDelta(
         gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_6;
   std::complex<double> tmp_7;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_6 += tmp_7;
   result += (2*Sqr(gp)*Sqr(Qd)*Sqr(Cos(ThetaWp()))) * tmp_6;
   std::complex<double> tmp_8;
   std::complex<double> tmp_9;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_9 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_8 += tmp_9;
   result += (1.0327955589886444*g1*gp*Qd*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_8;
   std::complex<double> tmp_10;
   std::complex<double> tmp_11;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_11 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_10 += tmp_11;
   result += (0.13333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())
      )) * tmp_10;
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())
         );
   }
   if (gO1 < 3) {
      result += 2*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 0.5163977794943222*g1*gp*Qq*Cos(ThetaWp())*KroneckerDelta(
         gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
         ThetaW()))*Sqr(Sin(ThetaWp()));
   }

   return result;
}

double CLASSNAME::CpUSdconjUSdconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   std::complex<double> tmp_12;
   std::complex<double> tmp_13;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_13 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_12 += tmp_13;
   result += (0.1*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_12;
   std::complex<double> tmp_14;
   std::complex<double> tmp_15;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_15 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_14 += tmp_15;
   result += (-(Qd*QHd*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0))) * tmp_14;
   std::complex<double> tmp_16;
   std::complex<double> tmp_17;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_18;
      std::complex<double> tmp_19;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_20;
         std::complex<double> tmp_21;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_21 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_20 += tmp_21;
         tmp_19 += (KroneckerDelta(gO2,3 + j2)) * tmp_20;
      }
      tmp_18 += tmp_19;
      tmp_17 += (KroneckerDelta(gO1,3 + j3)) * tmp_18;
   }
   tmp_16 += tmp_17;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_16;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)
         );
   }
   std::complex<double> tmp_22;
   std::complex<double> tmp_23;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_23 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_22 += tmp_23;
   result += (-0.1*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_22;
   std::complex<double> tmp_24;
   std::complex<double> tmp_25;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_25 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_24 += tmp_25;
   result += (-(Qd*QHu*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1))) * tmp_24;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)
         );
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_26;
      std::complex<double> tmp_27;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_27 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_26 += tmp_27;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_26;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_28;
   std::complex<double> tmp_29;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_29 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_28 += tmp_29;
   result += (0.1*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)) * tmp_28;
   std::complex<double> tmp_30;
   std::complex<double> tmp_31;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_31 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_30 += tmp_31;
   result += (-(Qd*QHd*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(gp))) * tmp_30;
   std::complex<double> tmp_32;
   std::complex<double> tmp_33;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_33 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_32 += tmp_33;
   result += (-0.1*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)) * tmp_32;
   std::complex<double> tmp_34;
   std::complex<double> tmp_35;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_35 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_34 += tmp_35;
   result += (-(Qd*QHu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp))) * tmp_34;
   std::complex<double> tmp_36;
   std::complex<double> tmp_37;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_37 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_36 += tmp_37;
   result += (-(Qd*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp))) * tmp_36;
   std::complex<double> tmp_38;
   std::complex<double> tmp_39;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_40;
      std::complex<double> tmp_41;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_42;
         std::complex<double> tmp_43;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_43 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_42 += tmp_43;
         tmp_41 += (KroneckerDelta(gO2,3 + j2)) * tmp_42;
      }
      tmp_40 += tmp_41;
      tmp_39 += (KroneckerDelta(gO1,3 + j3)) * tmp_40;
   }
   tmp_38 += tmp_39;
   result += (-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)))) * tmp_38;
   if (gO1 < 3) {
      result += 0.05*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -0.05*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Qq*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_44;
      std::complex<double> tmp_45;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_45 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_44 += tmp_45;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))) *
         tmp_44;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_46;
      std::complex<double> tmp_47;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_47 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_46 += tmp_47;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))) *
         tmp_46;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_48;
      std::complex<double> tmp_49;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_49 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_48 += tmp_49;
      result += (-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*Lambdax) * tmp_48;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_50;
      std::complex<double> tmp_51;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_51 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_50 += tmp_51;
      result += (-0.5*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))*Lambdax) * tmp_50;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_52;
      std::complex<double> tmp_53;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_53 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_52 += tmp_53;
      result += (-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)))) * tmp_52;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_54;
   std::complex<double> tmp_55;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_55 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_54 += tmp_55;
   result += (0.1*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_54;
   std::complex<double> tmp_56;
   std::complex<double> tmp_57;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_57 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_56 += tmp_57;
   result += (-(Qd*Ql*KroneckerDelta(gI1,gI2)*Sqr(gp))) * tmp_56;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1)
         ;
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2)
         ;
   }
   if (gO1 < 3) {
      result += -(Ql*Qq*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(
         gp));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_59 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_58 += tmp_59;
   result += (0.1*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_58;
   std::complex<double> tmp_60;
   std::complex<double> tmp_61;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_61 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_60 += tmp_61;
   result += (-(Qd*QHd*Sqr(gp)*ZH(gI1,0)*ZH(gI2,0))) * tmp_60;
   std::complex<double> tmp_62;
   std::complex<double> tmp_63;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_64;
      std::complex<double> tmp_65;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_66;
         std::complex<double> tmp_67;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_67 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_66 += tmp_67;
         tmp_65 += (KroneckerDelta(gO2,3 + j2)) * tmp_66;
      }
      tmp_64 += tmp_65;
      tmp_63 += (KroneckerDelta(gO1,3 + j3)) * tmp_64;
   }
   tmp_62 += tmp_63;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_62;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,0)*ZH(gI2,0)
         );
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_68;
      std::complex<double> tmp_69;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_69 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_68 += tmp_69;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_68;
   }
   std::complex<double> tmp_70;
   std::complex<double> tmp_71;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_71 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_70 += tmp_71;
   result += (-0.1*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_70;
   std::complex<double> tmp_72;
   std::complex<double> tmp_73;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_73 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_72 += tmp_73;
   result += (-(Qd*QHu*Sqr(gp)*ZH(gI1,1)*ZH(gI2,1))) * tmp_72;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,1)*ZH(gI2,1)
         );
   }
   if (gO1 < 3) {
      std::complex<double> tmp_74;
      std::complex<double> tmp_75;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_75 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_74 += tmp_75;
      result += (0.5*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,1)) * tmp_74;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_76;
      std::complex<double> tmp_77;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_77 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_76 += tmp_77;
      result += (0.5*Lambdax*ZH(gI1,2)*ZH(gI2,1)) * tmp_76;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_78;
      std::complex<double> tmp_79;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_79 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_78 += tmp_79;
      result += (0.5*Conj(Lambdax)*ZH(gI1,1)*ZH(gI2,2)) * tmp_78;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_80;
      std::complex<double> tmp_81;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_81 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_80 += tmp_81;
      result += (0.5*Lambdax*ZH(gI1,1)*ZH(gI2,2)) * tmp_80;
   }
   std::complex<double> tmp_82;
   std::complex<double> tmp_83;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_83 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_82 += tmp_83;
   result += (-(Qd*Qs*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2))) * tmp_82;
   if (gO1 < 3) {
      result += -(Qq*Qs*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2))
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_84;
      std::complex<double> tmp_85;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_85 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_84 += tmp_85;
      result += (UP(gI2,1)) * tmp_84;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_86;
   std::complex<double> tmp_87;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_88;
      std::complex<double> tmp_89;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_89 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_88 += tmp_89;
      tmp_87 += (Conj(ZUL(gI1,j2))) * tmp_88;
   }
   tmp_86 += tmp_87;
   result += (Conj(UM(gI2,1))) * tmp_86;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(ZUL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_90;
   std::complex<double> tmp_91;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_91 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_90 += tmp_91;
   result += (-1.4142135623730951*gp*Qd*ZN(gI2,0)) * tmp_90;
   std::complex<double> tmp_92;
   std::complex<double> tmp_93;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_93 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_92 += tmp_93;
   result += (-0.3651483716701107*g1*ZN(gI2,1)) * tmp_92;
   if (gO2 < 3) {
      std::complex<double> tmp_94;
      std::complex<double> tmp_95;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_95 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_94 += tmp_95;
      result += (-ZN(gI2,3)) * tmp_94;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_96;
   std::complex<double> tmp_97;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_98;
      std::complex<double> tmp_99;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_99 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_98 += tmp_99;
      tmp_97 += (Conj(ZDL(gI1,j2))) * tmp_98;
   }
   tmp_96 += tmp_97;
   result += (-Conj(ZN(gI2,3))) * tmp_96;
   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Qq*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,0))
         ;
   }
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,1));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_100;
   std::complex<double> tmp_102;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_102 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_100 += tmp_102;
   std::complex<double> tmp_101;
   std::complex<double> tmp_103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_103 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_101 += tmp_103;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_100 * tmp_101;
   std::complex<double> tmp_104;
   std::complex<double> tmp_106;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_106 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_104 += tmp_106;
   std::complex<double> tmp_105;
   std::complex<double> tmp_107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_107 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_105 += tmp_107;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_104 * tmp_105;
   std::complex<double> tmp_108;
   std::complex<double> tmp_110;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_110 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_108 += tmp_110;
   std::complex<double> tmp_109;
   std::complex<double> tmp_111;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_111 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_109 += tmp_111;
   result += (-0.5*Sqr(gp)*Sqr(Qd)) * tmp_108 * tmp_109;
   std::complex<double> tmp_112;
   std::complex<double> tmp_114;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_114 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_112 += tmp_114;
   std::complex<double> tmp_113;
   std::complex<double> tmp_115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_115 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_113 += tmp_115;
   result += (-0.05*Sqr(g1)) * tmp_112 * tmp_113;
   std::complex<double> tmp_116;
   std::complex<double> tmp_118;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_118 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_116 += tmp_118;
   std::complex<double> tmp_117;
   std::complex<double> tmp_119;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_119 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_117 += tmp_119;
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_116 * tmp_117;
   std::complex<double> tmp_120;
   std::complex<double> tmp_122;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_122 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_120 += tmp_122;
   std::complex<double> tmp_121;
   std::complex<double> tmp_123;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_123 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_121 += tmp_123;
   result += (-0.1*Sqr(g1)) * tmp_120 * tmp_121;
   std::complex<double> tmp_124;
   std::complex<double> tmp_126;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_126 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_124 += tmp_126;
   std::complex<double> tmp_125;
   std::complex<double> tmp_127;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_127 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_125 += tmp_127;
   result += (-1.5*Sqr(gp)*Sqr(Qd)) * tmp_124 * tmp_125;
   std::complex<double> tmp_128;
   std::complex<double> tmp_130;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_130 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_128 += tmp_130;
   std::complex<double> tmp_129;
   std::complex<double> tmp_131;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_131 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_129 += tmp_131;
   result += (-0.05*Sqr(g1)) * tmp_128 * tmp_129;
   std::complex<double> tmp_132;
   std::complex<double> tmp_134;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_134 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_132 += tmp_134;
   std::complex<double> tmp_133;
   std::complex<double> tmp_135;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_135 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_133 += tmp_135;
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_132 * tmp_133;
   std::complex<double> tmp_136;
   std::complex<double> tmp_138;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_138 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_136 += tmp_138;
   std::complex<double> tmp_137;
   std::complex<double> tmp_139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_139 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_137 += tmp_139;
   result += (-0.1*Sqr(g1)) * tmp_136 * tmp_137;
   std::complex<double> tmp_140;
   std::complex<double> tmp_142;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_142 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_140 += tmp_142;
   std::complex<double> tmp_141;
   std::complex<double> tmp_143;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_143 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_141 += tmp_143;
   result += (-1.5*Sqr(gp)*Sqr(Qd)) * tmp_140 * tmp_141;
   std::complex<double> tmp_144;
   std::complex<double> tmp_146;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_146 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_144 += tmp_146;
   std::complex<double> tmp_145;
   std::complex<double> tmp_147;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_147 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_145 += tmp_147;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_144 * tmp_145;
   std::complex<double> tmp_148;
   std::complex<double> tmp_150;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_150 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_148 += tmp_150;
   std::complex<double> tmp_149;
   std::complex<double> tmp_151;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_151 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_149 += tmp_151;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_148 * tmp_149;
   std::complex<double> tmp_152;
   std::complex<double> tmp_154;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_154 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_152 += tmp_154;
   std::complex<double> tmp_153;
   std::complex<double> tmp_155;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_155 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_153 += tmp_155;
   result += (-0.5*Sqr(gp)*Sqr(Qd)) * tmp_152 * tmp_153;
   std::complex<double> tmp_156;
   std::complex<double> tmp_158;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_159;
      std::complex<double> tmp_160;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_160 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_159 += tmp_160;
      tmp_158 += (Conj(ZD(gI2,j2))) * tmp_159;
   }
   tmp_156 += tmp_158;
   std::complex<double> tmp_157;
   std::complex<double> tmp_161;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_162;
      std::complex<double> tmp_163;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_163 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_162 += tmp_163;
      tmp_161 += (ZD(gI1,j4)) * tmp_162;
   }
   tmp_157 += tmp_161;
   result += (-1) * tmp_156 * tmp_157;
   if (gO1 < 3) {
      std::complex<double> tmp_164;
      std::complex<double> tmp_165;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_165 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_164 += tmp_165;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_164;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_166;
      std::complex<double> tmp_167;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_167 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_166 += tmp_167;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_166;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_168;
      std::complex<double> tmp_169;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_169 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_168 += tmp_169;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_168;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_170;
      std::complex<double> tmp_171;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_171 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_170 += tmp_171;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_170;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_172;
      std::complex<double> tmp_173;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_173 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_172 += tmp_173;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_172;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_174;
      std::complex<double> tmp_175;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_175 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_174 += tmp_175;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_174;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_176;
      std::complex<double> tmp_177;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_177 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_176 += tmp_177;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_176;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_178;
      std::complex<double> tmp_179;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_179 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_178 += tmp_179;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_178;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_181;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_181 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_180 += tmp_181;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_180;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_182;
      std::complex<double> tmp_183;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_183 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_182 += tmp_183;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_182;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_184;
      std::complex<double> tmp_186;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_186 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_184 += tmp_186;
      std::complex<double> tmp_185;
      std::complex<double> tmp_187;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_188;
         std::complex<double> tmp_189;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_189 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_188 += tmp_189;
         tmp_187 += (ZD(gI1,j4)) * tmp_188;
      }
      tmp_185 += tmp_187;
      result += (-3) * tmp_184 * tmp_185;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_190;
      std::complex<double> tmp_191;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_191 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_190 += tmp_191;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_190;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_192;
      std::complex<double> tmp_193;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_193 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_192 += tmp_193;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_192;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_194;
      std::complex<double> tmp_195;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_195 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_194 += tmp_195;
      result += (-0.5*Qd*Qq*Conj(ZD(gI2,gO2))*Sqr(gp)) * tmp_194;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_196;
      std::complex<double> tmp_197;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_197 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_196 += tmp_197;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_196;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_198;
      std::complex<double> tmp_199;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_199 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_198 += tmp_199;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_198;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_200;
      std::complex<double> tmp_201;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_201 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_200 += tmp_201;
      result += (-0.5*Qd*Qq*Conj(ZD(gI2,gO2))*Sqr(gp)) * tmp_200;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_202;
      std::complex<double> tmp_204;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_205;
         std::complex<double> tmp_206;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_206 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_205 += tmp_206;
         tmp_204 += (Conj(ZD(gI2,j2))) * tmp_205;
      }
      tmp_202 += tmp_204;
      std::complex<double> tmp_203;
      std::complex<double> tmp_207;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_207 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_203 += tmp_207;
      result += (-3) * tmp_202 * tmp_203;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_208;
      std::complex<double> tmp_210;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_210 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_208 += tmp_210;
      std::complex<double> tmp_209;
      std::complex<double> tmp_211;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_211 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_209 += tmp_211;
      result += (-1) * tmp_208 * tmp_209;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_212;
      std::complex<double> tmp_213;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_213 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_212 += tmp_213;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_212;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_214;
      std::complex<double> tmp_215;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_215 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_214 += tmp_215;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_214;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_216;
      std::complex<double> tmp_217;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_217 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_216 += tmp_217;
      result += (-0.5*Qd*Qq*Sqr(gp)*ZD(gI1,gO1)) * tmp_216;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_218;
      std::complex<double> tmp_219;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_219 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_218 += tmp_219;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_218;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_220;
      std::complex<double> tmp_221;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_221 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_220 += tmp_221;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_220;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_222;
      std::complex<double> tmp_223;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_223 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_222 += tmp_223;
      result += (-0.5*Qd*Qq*Sqr(gp)*ZD(gI1,gO1)) * tmp_222;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)*ZD(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -1.3333333333333333*Conj(ZD(gI2,gO2))*Sqr(g3)*ZD(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -(Conj(ZD(gI2,gO2))*Sqr(gp)*Sqr(Qq)*ZD(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_224;
   std::complex<double> tmp_226;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_226 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_224 += tmp_226;
   std::complex<double> tmp_225;
   std::complex<double> tmp_227;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_227 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_225 += tmp_227;
   result += (0.05*Sqr(g1)) * tmp_224 * tmp_225;
   std::complex<double> tmp_228;
   std::complex<double> tmp_230;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_230 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_228 += tmp_230;
   std::complex<double> tmp_229;
   std::complex<double> tmp_231;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_231 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_229 += tmp_231;
   result += (-0.5*Qd*Ql*Sqr(gp)) * tmp_228 * tmp_229;
   std::complex<double> tmp_232;
   std::complex<double> tmp_234;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_234 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_232 += tmp_234;
   std::complex<double> tmp_233;
   std::complex<double> tmp_235;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_235 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_233 += tmp_235;
   result += (-0.1*Sqr(g1)) * tmp_232 * tmp_233;
   std::complex<double> tmp_236;
   std::complex<double> tmp_238;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_238 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_236 += tmp_238;
   std::complex<double> tmp_237;
   std::complex<double> tmp_239;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_239 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_237 += tmp_239;
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_236 * tmp_237;
   std::complex<double> tmp_240;
   std::complex<double> tmp_242;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_242 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_240 += tmp_242;
   std::complex<double> tmp_241;
   std::complex<double> tmp_243;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_243 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_241 += tmp_243;
   result += (0.05*Sqr(g1)) * tmp_240 * tmp_241;
   std::complex<double> tmp_244;
   std::complex<double> tmp_246;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_246 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_244 += tmp_246;
   std::complex<double> tmp_245;
   std::complex<double> tmp_247;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_247 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_245 += tmp_247;
   result += (-0.5*Qd*Ql*Sqr(gp)) * tmp_244 * tmp_245;
   std::complex<double> tmp_248;
   std::complex<double> tmp_250;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_250 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_248 += tmp_250;
   std::complex<double> tmp_249;
   std::complex<double> tmp_251;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_251 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_249 += tmp_251;
   result += (-0.1*Sqr(g1)) * tmp_248 * tmp_249;
   std::complex<double> tmp_252;
   std::complex<double> tmp_254;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_254 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_252 += tmp_254;
   std::complex<double> tmp_253;
   std::complex<double> tmp_255;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_255 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_253 += tmp_255;
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_252 * tmp_253;
   if (gO1 < 3) {
      std::complex<double> tmp_256;
      std::complex<double> tmp_257;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_257 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_256 += tmp_257;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_256;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_258;
      std::complex<double> tmp_259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_259 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_258 += tmp_259;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_258;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_260;
      std::complex<double> tmp_261;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_261 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_260 += tmp_261;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_260;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_262;
      std::complex<double> tmp_263;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_263 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_262 += tmp_263;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_262;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_264;
      std::complex<double> tmp_265;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_265 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_264 += tmp_265;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_264;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_266;
      std::complex<double> tmp_267;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_267 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_266 += tmp_267;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_266;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_268;
      std::complex<double> tmp_269;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_269 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_268 += tmp_269;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_268;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_270;
      std::complex<double> tmp_271;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_271 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_270 += tmp_271;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_270;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_272;
      std::complex<double> tmp_273;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_273 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_272 += tmp_273;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_272;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_274;
      std::complex<double> tmp_275;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_275 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_274 += tmp_275;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_274;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_276;
      std::complex<double> tmp_278;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_278 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_276 += tmp_278;
      std::complex<double> tmp_277;
      std::complex<double> tmp_279;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_280;
         std::complex<double> tmp_281;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_281 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_280 += tmp_281;
         tmp_279 += (ZE(gI1,j4)) * tmp_280;
      }
      tmp_277 += tmp_279;
      result += (-1) * tmp_276 * tmp_277;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_282;
      std::complex<double> tmp_284;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_285;
         std::complex<double> tmp_286;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_286 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_285 += tmp_286;
         tmp_284 += (Conj(ZE(gI2,j2))) * tmp_285;
      }
      tmp_282 += tmp_284;
      std::complex<double> tmp_283;
      std::complex<double> tmp_287;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_287 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_283 += tmp_287;
      result += (-1) * tmp_282 * tmp_283;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_288;
   std::complex<double> tmp_290;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_290 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_288 += tmp_290;
   std::complex<double> tmp_289;
   std::complex<double> tmp_291;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_291 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_289 += tmp_291;
   result += (-0.05*Sqr(g1)) * tmp_288 * tmp_289;
   std::complex<double> tmp_292;
   std::complex<double> tmp_294;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_294 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_292 += tmp_294;
   std::complex<double> tmp_293;
   std::complex<double> tmp_295;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_295 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_293 += tmp_295;
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_292 * tmp_293;
   std::complex<double> tmp_296;
   std::complex<double> tmp_298;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_298 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_296 += tmp_298;
   std::complex<double> tmp_297;
   std::complex<double> tmp_299;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_299 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_297 += tmp_299;
   result += (0.2*Sqr(g1)) * tmp_296 * tmp_297;
   std::complex<double> tmp_300;
   std::complex<double> tmp_302;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_302 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_300 += tmp_302;
   std::complex<double> tmp_301;
   std::complex<double> tmp_303;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_303 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_301 += tmp_303;
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_300 * tmp_301;
   std::complex<double> tmp_304;
   std::complex<double> tmp_306;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_306 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_304 += tmp_306;
   std::complex<double> tmp_305;
   std::complex<double> tmp_307;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_307 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_305 += tmp_307;
   result += (-0.05*Sqr(g1)) * tmp_304 * tmp_305;
   std::complex<double> tmp_308;
   std::complex<double> tmp_310;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_310 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_308 += tmp_310;
   std::complex<double> tmp_309;
   std::complex<double> tmp_311;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_311 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_309 += tmp_311;
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_308 * tmp_309;
   std::complex<double> tmp_312;
   std::complex<double> tmp_314;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_314 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_312 += tmp_314;
   std::complex<double> tmp_313;
   std::complex<double> tmp_315;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_315 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_313 += tmp_315;
   result += (0.2*Sqr(g1)) * tmp_312 * tmp_313;
   std::complex<double> tmp_316;
   std::complex<double> tmp_318;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_318 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_316 += tmp_318;
   std::complex<double> tmp_317;
   std::complex<double> tmp_319;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_319 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_317 += tmp_319;
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_316 * tmp_317;
   std::complex<double> tmp_320;
   std::complex<double> tmp_322;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_323;
      std::complex<double> tmp_324;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_324 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_323 += tmp_324;
      tmp_322 += (Conj(ZU(gI2,j2))) * tmp_323;
   }
   tmp_320 += tmp_322;
   std::complex<double> tmp_321;
   std::complex<double> tmp_325;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_326;
      std::complex<double> tmp_327;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_327 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_326 += tmp_327;
      tmp_325 += (ZU(gI1,j4)) * tmp_326;
   }
   tmp_321 += tmp_325;
   result += (-1) * tmp_320 * tmp_321;
   if (gO1 < 3) {
      std::complex<double> tmp_328;
      std::complex<double> tmp_329;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_329 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_328 += tmp_329;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_328;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_330;
      std::complex<double> tmp_331;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_331 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_330 += tmp_331;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_330;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_332;
      std::complex<double> tmp_333;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_333 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_332 += tmp_333;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_332;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_334;
      std::complex<double> tmp_335;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_335 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_334 += tmp_335;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_334;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_336;
      std::complex<double> tmp_337;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_337 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_336 += tmp_337;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_336;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_338;
      std::complex<double> tmp_339;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_339 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_338 += tmp_339;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_338;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_340;
      std::complex<double> tmp_341;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_341 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_340 += tmp_341;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_340;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_342;
      std::complex<double> tmp_343;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_343 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_342 += tmp_343;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_342;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_344;
      std::complex<double> tmp_345;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_345 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_344 += tmp_345;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_344;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_346;
      std::complex<double> tmp_347;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_347 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_346 += tmp_347;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_346;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_348;
      std::complex<double> tmp_350;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_350 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_348 += tmp_350;
      std::complex<double> tmp_349;
      std::complex<double> tmp_351;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_351 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_349 += tmp_351;
      result += (-1) * tmp_348 * tmp_349;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_352;
   std::complex<double> tmp_353;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_354;
      std::complex<double> tmp_355;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_355 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_354 += tmp_355;
      tmp_353 += (Conj(ZU(gI1,j2))) * tmp_354;
   }
   tmp_352 += tmp_353;
   result += (ZP(gI2,0)) * tmp_352;
   std::complex<double> tmp_356;
   std::complex<double> tmp_357;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_358;
      std::complex<double> tmp_359;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_360;
         std::complex<double> tmp_361;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_361 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_360 += tmp_361;
         tmp_359 += (KroneckerDelta(gO2,3 + j2)) * tmp_360;
      }
      tmp_358 += tmp_359;
      tmp_357 += (Conj(ZU(gI1,3 + j3))) * tmp_358;
   }
   tmp_356 += tmp_357;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_356;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_362;
      std::complex<double> tmp_363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_363 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_362 += tmp_363;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI2,0)) * tmp_362;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_364;
      std::complex<double> tmp_365;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_366;
         std::complex<double> tmp_367;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_367 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_366 += tmp_367;
         tmp_365 += (Conj(ZU(gI1,j2))) * tmp_366;
      }
      tmp_364 += tmp_365;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_364;
   }
   std::complex<double> tmp_368;
   std::complex<double> tmp_369;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_370;
      std::complex<double> tmp_371;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_371 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_370 += tmp_371;
      tmp_369 += (Conj(ZU(gI1,j2))) * tmp_370;
   }
   tmp_368 += tmp_369;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI2,1)) * tmp_368;
   std::complex<double> tmp_372;
   std::complex<double> tmp_373;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_374;
      std::complex<double> tmp_375;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_376;
         std::complex<double> tmp_377;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_377 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_376 += tmp_377;
         tmp_375 += (KroneckerDelta(gO2,3 + j2)) * tmp_376;
      }
      tmp_374 += tmp_375;
      tmp_373 += (Conj(ZU(gI1,3 + j3))) * tmp_374;
   }
   tmp_372 += tmp_373;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_372;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_378;
      std::complex<double> tmp_379;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_379 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_378 += tmp_379;
      result += (ZP(gI2,1)) * tmp_378;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_380;
      std::complex<double> tmp_381;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_382;
         std::complex<double> tmp_383;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_383 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_382 += tmp_383;
         tmp_381 += (Conj(ZU(gI1,j2))) * tmp_382;
      }
      tmp_380 += tmp_381;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_380;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_384;
   std::complex<double> tmp_385;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_386;
      std::complex<double> tmp_387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_387 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_386 += tmp_387;
      tmp_385 += (Conj(ZD(gI1,j2))) * tmp_386;
   }
   tmp_384 += tmp_385;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,1))) *
      tmp_384;
   std::complex<double> tmp_388;
   std::complex<double> tmp_389;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_390;
      std::complex<double> tmp_391;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_391 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_390 += tmp_391;
      tmp_389 += (Conj(ZD(gI1,j2))) * tmp_390;
   }
   tmp_388 += tmp_389;
   result += (std::complex<double>(0,-0.5)*vu*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_388;
   std::complex<double> tmp_392;
   std::complex<double> tmp_393;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_394;
      std::complex<double> tmp_395;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_395 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_394 += tmp_395;
      tmp_393 += (Conj(ZD(gI1,j2))) * tmp_394;
   }
   tmp_392 += tmp_393;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_392;
   if (gO2 < 3) {
      std::complex<double> tmp_396;
      std::complex<double> tmp_397;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_397 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_396 += tmp_397;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,1))*Lambdax) *
         tmp_396;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_398;
      std::complex<double> tmp_399;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_399 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_398 += tmp_399;
      result += (std::complex<double>(0,0.5)*vu*Conj(ZA(gI2,2))*Lambdax) *
         tmp_398;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_400;
      std::complex<double> tmp_401;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_401 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_400 += tmp_401;
      result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,0)))
         * tmp_400;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   std::complex<double> tmp_402;
   std::complex<double> tmp_403;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_403 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_402 += tmp_403;
   result += (0.1*vd*Sqr(g1)*ZH(gI2,0)) * tmp_402;
   std::complex<double> tmp_404;
   std::complex<double> tmp_405;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_405 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_404 += tmp_405;
   result += (-(Qd*QHd*vd*Sqr(gp)*ZH(gI2,0))) * tmp_404;
   std::complex<double> tmp_406;
   std::complex<double> tmp_407;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_408;
      std::complex<double> tmp_409;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_409 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_408 += tmp_409;
      tmp_407 += (Conj(ZD(gI1,j2))) * tmp_408;
   }
   tmp_406 += tmp_407;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_406;
   std::complex<double> tmp_410;
   std::complex<double> tmp_411;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_412;
      std::complex<double> tmp_413;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_414;
         std::complex<double> tmp_415;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_415 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_414 += tmp_415;
         tmp_413 += (KroneckerDelta(gO2,3 + j2)) * tmp_414;
      }
      tmp_412 += tmp_413;
      tmp_411 += (Conj(ZD(gI1,3 + j3))) * tmp_412;
   }
   tmp_410 += tmp_411;
   result += (-(vd*ZH(gI2,0))) * tmp_410;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -(QHd*Qq*vd*Conj(ZD(gI1,gO2))*Sqr(gp)*ZH(gI2,0));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_416;
      std::complex<double> tmp_417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_417 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_416 += tmp_417;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_416;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_418;
      std::complex<double> tmp_419;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_420;
         std::complex<double> tmp_421;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_421 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_420 += tmp_421;
         tmp_419 += (Conj(ZD(gI1,j2))) * tmp_420;
      }
      tmp_418 += tmp_419;
      result += (-(vd*ZH(gI2,0))) * tmp_418;
   }
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_423 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_422 += tmp_423;
   result += (-0.1*vu*Sqr(g1)*ZH(gI2,1)) * tmp_422;
   std::complex<double> tmp_424;
   std::complex<double> tmp_425;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_425 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_424 += tmp_425;
   result += (-(Qd*QHu*vu*Sqr(gp)*ZH(gI2,1))) * tmp_424;
   std::complex<double> tmp_426;
   std::complex<double> tmp_427;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_428;
      std::complex<double> tmp_429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_429 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_428 += tmp_429;
      tmp_427 += (Conj(ZD(gI1,j2))) * tmp_428;
   }
   tmp_426 += tmp_427;
   result += (0.5*vS*Conj(Lambdax)*ZH(gI2,1)) * tmp_426;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -(QHu*Qq*vu*Conj(ZD(gI1,gO2))*Sqr(gp)*ZH(gI2,1));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_430;
      std::complex<double> tmp_431;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_431 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_430 += tmp_431;
      result += (0.5*vS*Lambdax*ZH(gI2,1)) * tmp_430;
   }
   std::complex<double> tmp_432;
   std::complex<double> tmp_433;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_433 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_432 += tmp_433;
   result += (-(Qd*Qs*vS*Sqr(gp)*ZH(gI2,2))) * tmp_432;
   std::complex<double> tmp_434;
   std::complex<double> tmp_435;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_436;
      std::complex<double> tmp_437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_437 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_436 += tmp_437;
      tmp_435 += (Conj(ZD(gI1,j2))) * tmp_436;
   }
   tmp_434 += tmp_435;
   result += (0.5*vu*Conj(Lambdax)*ZH(gI2,2)) * tmp_434;
   if (gO2 < 3) {
      result += -(Qq*Qs*vS*Conj(ZD(gI1,gO2))*Sqr(gp)*ZH(gI2,2));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_438;
      std::complex<double> tmp_439;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_439 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_438 += tmp_439;
      result += (0.5*vu*Lambdax*ZH(gI2,2)) * tmp_438;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_440;
   std::complex<double> tmp_441;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_441 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_440 += tmp_441;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_440;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPL(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*PhaseGlu*Conj(ZDL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVGSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 6) {
      result += g3*Conj(ZD(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVPSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_442;
   std::complex<double> tmp_443;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_443 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_442 += tmp_443;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_442;
   if (gO2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZD(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVZSd(unsigned gO2, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_444;
   std::complex<double> tmp_445;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_445 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_444 += tmp_445;
   result += (0.2581988897471611*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_444;
   std::complex<double> tmp_446;
   std::complex<double> tmp_447;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_447 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_446 += tmp_447;
   result += (-(gp*Qd*Sin(ThetaWp()))) * tmp_446;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZD(gI2,gO2))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Cos(ThetaWp())*Sin
         (ThetaW());
   }
   if (gO2 < 3) {
      result += gp*Qq*Conj(ZD(gI2,gO2))*Sin(ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVZpSd(unsigned gO2, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_448;
   std::complex<double> tmp_449;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_449 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_448 += tmp_449;
   result += (-(gp*Qd*Cos(ThetaWp()))) * tmp_448;
   std::complex<double> tmp_450;
   std::complex<double> tmp_451;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_451 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_450 += tmp_451;
   result += (-0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_450;
   if (gO2 < 3) {
      result += gp*Qq*Conj(ZD(gI2,gO2))*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZD(gI2,gO2))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gO2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdVWmSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZU(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   result = 0.1*KroneckerDelta(gO1,gO2)*(20*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp(
      ))*Sin(ThetaWp()) + 15.491933384829668*g1*gp*Ql*Cos(ThetaWp())*Sin(ThetaW())
      *Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*
      g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos
      (ThetaWp())) + 20*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   result = 0.1*KroneckerDelta(gO1,gO2)*(-2*gp*Ql*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(Ql)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp())));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   result += -(QHd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0));
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_452;
      std::complex<double> tmp_453;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_453 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_452 += tmp_453;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_452;
   }
   result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   result += -(QHu*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_454;
      std::complex<double> tmp_455;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_455 += Conj(Ye(j1,gO2))*ZER(gI2,j1);
      }
      tmp_454 += tmp_455;
      result += (UM(gI1,1)) * tmp_454;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZEL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_456;
      std::complex<double> tmp_457;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_457 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_456 += tmp_457;
      result += (ZP(gI1,0)) * tmp_456;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_458;
      std::complex<double> tmp_459;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_460;
         std::complex<double> tmp_461;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_461 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_460 += tmp_461;
         tmp_459 += (Conj(ZE(gI2,j2))) * tmp_460;
      }
      tmp_458 += tmp_459;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_458;
   }
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_462;
      std::complex<double> tmp_463;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_463 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_462 += tmp_463;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI1,1)) * tmp_462;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(20*Ql*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2
      ,2))*Sqr(gp) + Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(-3*Sqr(g1) - 5*Sqr(g2) + 20*
      QHu*Ql*Sqr(gp)) + Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(3*Sqr(g1) + 5*(Sqr(g2) +
      4*QHd*Ql*Sqr(gp))));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   result += -0.15*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1);
   result += -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2);
   result += -(KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql))
      ;
   if (gI1 < 3 && gI2 < 3) {
      result += -0.15*Conj(ZV(gI2,gO2))*Sqr(g1)*ZV(gI1,gO1);
   }
   if (gI1 < 3 && gI2 < 3) {
      result += -0.25*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }
   if (gI1 < 3 && gI2 < 3) {
      result += -(Conj(ZV(gI2,gO2))*Sqr(gp)*Sqr(Ql)*ZV(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*((3*Sqr(g1) + 5*(Sqr(g2) + 4*QHd*Ql*
      Sqr(gp)))*ZH(gI1,0)*ZH(gI2,0) + (-3*Sqr(g1) - 5*Sqr(g2) + 20*QHu*Ql*Sqr(gp))
      *ZH(gI1,1)*ZH(gI2,1) + 20*Ql*Qs*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   if (gI1 < 3) {
      result += -0.15*vd*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gI1 < 3) {
      result += -0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gI1 < 3) {
      result += -(QHd*Ql*vd*Conj(ZV(gI1,gO2))*Sqr(gp)*ZH(gI2,0));
   }
   if (gI1 < 3) {
      result += 0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gI1 < 3) {
      result += 0.25*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gI1 < 3) {
      result += -(QHu*Ql*vu*Conj(ZV(gI1,gO2))*Sqr(gp)*ZH(gI2,1));
   }
   if (gI1 < 3) {
      result += -(Ql*Qs*vS*Conj(ZV(gI1,gO2))*Sqr(gp)*ZH(gI2,2));
   }

   return result;
}

double CLASSNAME::CpconjUSvFvChiPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI1 < 3) {
      result += -1.4142135623730951*gp*Ql*Conj(ZN(gI2,0))*KroneckerDelta(gI1
         ,gO1);
   }
   if (gI1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZN(gI2,1))*KroneckerDelta(gI1,gO1
         );
   }
   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,2))*KroneckerDelta(gI1,
         gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_464;
   std::complex<double> tmp_465;
   std::complex<double> tmp_466;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_466 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_465 += tmp_466;
   tmp_464 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_465;
   std::complex<double> tmp_467;
   std::complex<double> tmp_468;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_468 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_467 += tmp_468;
   tmp_464 += (std::complex<double>(0,0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_467;
   std::complex<double> tmp_469;
   std::complex<double> tmp_470;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_470 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_469 += tmp_470;
   tmp_464 += (std::complex<double>(0,-1)*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ) * tmp_469;
   std::complex<double> tmp_471;
   std::complex<double> tmp_472;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_472 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_471 += tmp_472;
   tmp_464 += (std::complex<double>(0,0.1)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_471;
   std::complex<double> tmp_473;
   std::complex<double> tmp_474;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_474 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_473 += tmp_474;
   tmp_464 += (std::complex<double>(0,-1)*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ) * tmp_473;
   result += (std::complex<double>(0,-1)) * tmp_464;

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_475;
   std::complex<double> tmp_476;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_476 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_475 += tmp_476;
   result += (-0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_475;
   std::complex<double> tmp_477;
   std::complex<double> tmp_478;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_478 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_477 += tmp_478;
   result += (0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_477;
   std::complex<double> tmp_479;
   std::complex<double> tmp_480;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_480 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_479 += tmp_480;
   result += (-(KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql))) * tmp_479;
   std::complex<double> tmp_481;
   std::complex<double> tmp_482;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_482 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_481 += tmp_482;
   result += (0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_481;
   std::complex<double> tmp_483;
   std::complex<double> tmp_484;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_484 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_483 += tmp_484;
   result += (-(Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp))) * tmp_483;
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_485;
      std::complex<double> tmp_487;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_487 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_485 += tmp_487;
      std::complex<double> tmp_486;
      std::complex<double> tmp_488;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_488 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_486 += tmp_488;
      result += (-1) * tmp_485 * tmp_486;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_489;
   std::complex<double> tmp_490;
   std::complex<double> tmp_491;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_491 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_490 += tmp_491;
   tmp_489 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_490;
   std::complex<double> tmp_492;
   std::complex<double> tmp_493;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_493 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_492 += tmp_493;
   tmp_489 += (std::complex<double>(0,-0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_492;
   std::complex<double> tmp_494;
   std::complex<double> tmp_495;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_495 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_494 += tmp_495;
   tmp_489 += (std::complex<double>(0,-1)*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ) * tmp_494;
   std::complex<double> tmp_496;
   std::complex<double> tmp_497;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_497 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_496 += tmp_497;
   tmp_489 += (std::complex<double>(0,-0.2)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_496;
   std::complex<double> tmp_498;
   std::complex<double> tmp_499;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_499 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_498 += tmp_499;
   tmp_489 += (std::complex<double>(0,-1)*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)
      ) * tmp_498;
   result += (std::complex<double>(0,-1)) * tmp_489;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvVZSv(unsigned gO2, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Cos(ThetaWp())*Sin(
         ThetaW());
   }
   if (gI2 < 3) {
      result += gp*Ql*Conj(ZV(gI2,gO2))*Sin(ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvVZpSv(unsigned gO2, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += gp*Ql*Conj(ZV(gI2,gO2))*Cos(ThetaWp());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjVWmSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZE(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_500;
   std::complex<double> tmp_501;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_501 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_500 += tmp_501;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))
      ) * tmp_500;
   std::complex<double> tmp_502;
   std::complex<double> tmp_503;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_503 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_502 += tmp_503;
   result += (2.065591117977289*g1*gp*Qu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_502;
   std::complex<double> tmp_504;
   std::complex<double> tmp_505;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_505 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_504 += tmp_505;
   result += (2*Sqr(gp)*Sqr(Qu)*Sqr(Sin(ThetaWp()))) * tmp_504;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += -0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(
         ThetaWp()))*Sqr(Sin(ThetaW()));
   }
   if (gO1 < 3) {
      result += 2*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += -0.5163977794943222*g1*gp*Qq*Cos(ThetaWp())*KroneckerDelta(
         gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_506;
   std::complex<double> tmp_507;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_507 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_506 += tmp_507;
   result += (2*Sqr(gp)*Sqr(Qu)*Sqr(Cos(ThetaWp()))) * tmp_506;
   std::complex<double> tmp_508;
   std::complex<double> tmp_509;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_509 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_508 += tmp_509;
   result += (-2.065591117977289*g1*gp*Qu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_508;
   std::complex<double> tmp_510;
   std::complex<double> tmp_511;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_511 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_510 += tmp_511;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))
      ) * tmp_510;
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())
         );
   }
   if (gO1 < 3) {
      result += -2*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 0.5163977794943222*g1*gp*Qq*Cos(ThetaWp())*KroneckerDelta(
         gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += -0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
         ThetaW()))*Sqr(Sin(ThetaWp()));
   }

   return result;
}

double CLASSNAME::CpUSuconjUSuconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   std::complex<double> tmp_512;
   std::complex<double> tmp_513;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_513 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_512 += tmp_513;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_512;
   std::complex<double> tmp_514;
   std::complex<double> tmp_515;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_515 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_514 += tmp_515;
   result += (-(QHd*Qu*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0))) * tmp_514;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)
         );
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_516;
      std::complex<double> tmp_517;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_517 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_516 += tmp_517;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_516;
   }
   std::complex<double> tmp_518;
   std::complex<double> tmp_519;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_519 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_518 += tmp_519;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_518;
   std::complex<double> tmp_520;
   std::complex<double> tmp_521;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_521 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_520 += tmp_521;
   result += (-(QHu*Qu*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1))) * tmp_520;
   std::complex<double> tmp_522;
   std::complex<double> tmp_523;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_524;
      std::complex<double> tmp_525;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_526;
         std::complex<double> tmp_527;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_527 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_526 += tmp_527;
         tmp_525 += (KroneckerDelta(gO2,3 + j2)) * tmp_526;
      }
      tmp_524 += tmp_525;
      tmp_523 += (KroneckerDelta(gO1,3 + j3)) * tmp_524;
   }
   tmp_522 += tmp_523;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_522;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_528;
      std::complex<double> tmp_529;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_529 += Conj(Yd(j1,gO2))*ZDR(gI2,j1);
      }
      tmp_528 += tmp_529;
      result += (UM(gI1,1)) * tmp_528;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_530;
   std::complex<double> tmp_531;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_532;
      std::complex<double> tmp_533;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_533 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_532 += tmp_533;
      tmp_531 += (Conj(ZDL(gI2,j2))) * tmp_532;
   }
   tmp_530 += tmp_531;
   result += (Conj(UP(gI1,1))) * tmp_530;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_534;
   std::complex<double> tmp_535;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_536;
      std::complex<double> tmp_537;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_537 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_536 += tmp_537;
      tmp_535 += (Conj(ZD(gI2,j2))) * tmp_536;
   }
   tmp_534 += tmp_535;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI1,0)) * tmp_534;
   std::complex<double> tmp_538;
   std::complex<double> tmp_539;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_540;
      std::complex<double> tmp_541;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_542;
         std::complex<double> tmp_543;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_543 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_542 += tmp_543;
         tmp_541 += (KroneckerDelta(gO2,3 + j2)) * tmp_542;
      }
      tmp_540 += tmp_541;
      tmp_539 += (Conj(ZD(gI2,3 + j3))) * tmp_540;
   }
   tmp_538 += tmp_539;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_538;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_544;
      std::complex<double> tmp_545;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_545 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_544 += tmp_545;
      result += (ZP(gI1,0)) * tmp_544;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_546;
      std::complex<double> tmp_547;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_548;
         std::complex<double> tmp_549;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_549 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_548 += tmp_549;
         tmp_547 += (Conj(ZD(gI2,j2))) * tmp_548;
      }
      tmp_546 += tmp_547;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_546;
   }
   std::complex<double> tmp_550;
   std::complex<double> tmp_551;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_552;
      std::complex<double> tmp_553;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_553 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_552 += tmp_553;
      tmp_551 += (Conj(ZD(gI2,j2))) * tmp_552;
   }
   tmp_550 += tmp_551;
   result += (ZP(gI1,1)) * tmp_550;
   std::complex<double> tmp_554;
   std::complex<double> tmp_555;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_556;
      std::complex<double> tmp_557;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_558;
         std::complex<double> tmp_559;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_559 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_558 += tmp_559;
         tmp_557 += (KroneckerDelta(gO2,3 + j2)) * tmp_558;
      }
      tmp_556 += tmp_557;
      tmp_555 += (Conj(ZD(gI2,3 + j3))) * tmp_556;
   }
   tmp_554 += tmp_555;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_554;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_560;
      std::complex<double> tmp_561;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_561 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_560 += tmp_561;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI1,1)) * tmp_560;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_562;
      std::complex<double> tmp_563;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_564;
         std::complex<double> tmp_565;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_565 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_564 += tmp_565;
         tmp_563 += (Conj(ZD(gI2,j2))) * tmp_564;
      }
      tmp_562 += tmp_563;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_562;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_566;
   std::complex<double> tmp_567;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_567 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_566 += tmp_567;
   result += (-0.2*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)) * tmp_566;
   std::complex<double> tmp_568;
   std::complex<double> tmp_569;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_569 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_568 += tmp_569;
   result += (-(QHd*Qu*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(gp))) * tmp_568;
   std::complex<double> tmp_570;
   std::complex<double> tmp_571;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_571 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_570 += tmp_571;
   result += (0.2*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)) * tmp_570;
   std::complex<double> tmp_572;
   std::complex<double> tmp_573;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_573 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_572 += tmp_573;
   result += (-(QHu*Qu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp))) * tmp_572;
   std::complex<double> tmp_574;
   std::complex<double> tmp_575;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_575 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_574 += tmp_575;
   result += (-(Qs*Qu*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp))) * tmp_574;
   std::complex<double> tmp_576;
   std::complex<double> tmp_577;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_578;
      std::complex<double> tmp_579;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_580;
         std::complex<double> tmp_581;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_581 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_580 += tmp_581;
         tmp_579 += (KroneckerDelta(gO2,3 + j2)) * tmp_580;
      }
      tmp_578 += tmp_579;
      tmp_577 += (KroneckerDelta(gO1,3 + j3)) * tmp_578;
   }
   tmp_576 += tmp_577;
   result += (-(Conj(ZA(gI1,1))*Conj(ZA(gI2,1)))) * tmp_576;
   if (gO1 < 3) {
      result += 0.05*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -0.05*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Qq*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_582;
      std::complex<double> tmp_583;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_583 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_582 += tmp_583;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))) *
         tmp_582;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_584;
      std::complex<double> tmp_585;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_585 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_584 += tmp_585;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))) *
         tmp_584;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_586;
      std::complex<double> tmp_587;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_587 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_586 += tmp_587;
      result += (-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*Lambdax) * tmp_586;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_588;
      std::complex<double> tmp_589;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_589 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_588 += tmp_589;
      result += (-0.5*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))*Lambdax) * tmp_588;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_590;
      std::complex<double> tmp_591;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_591 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_590 += tmp_591;
      result += (-(Conj(ZA(gI1,1))*Conj(ZA(gI2,1)))) * tmp_590;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_592;
   std::complex<double> tmp_593;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_593 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_592 += tmp_593;
   result += (-0.2*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_592;
   std::complex<double> tmp_594;
   std::complex<double> tmp_595;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_595 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_594 += tmp_595;
   result += (-(Ql*Qu*KroneckerDelta(gI1,gI2)*Sqr(gp))) * tmp_594;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1)
         ;
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2
         );
   }
   if (gO1 < 3) {
      result += -(Ql*Qq*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(
         gp));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   std::complex<double> tmp_596;
   std::complex<double> tmp_597;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_597 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_596 += tmp_597;
   result += (-0.2*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_596;
   std::complex<double> tmp_598;
   std::complex<double> tmp_599;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_599 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_598 += tmp_599;
   result += (-(QHd*Qu*Sqr(gp)*ZH(gI1,0)*ZH(gI2,0))) * tmp_598;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,0)*ZH(gI2,0)
         );
   }
   if (gO1 < 3) {
      std::complex<double> tmp_600;
      std::complex<double> tmp_601;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_601 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_600 += tmp_601;
      result += (0.5*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,0)) * tmp_600;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_602;
      std::complex<double> tmp_603;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_603 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_602 += tmp_603;
      result += (0.5*Lambdax*ZH(gI1,2)*ZH(gI2,0)) * tmp_602;
   }
   std::complex<double> tmp_604;
   std::complex<double> tmp_605;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_605 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_604 += tmp_605;
   result += (0.2*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_604;
   std::complex<double> tmp_606;
   std::complex<double> tmp_607;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_607 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_606 += tmp_607;
   result += (-(QHu*Qu*Sqr(gp)*ZH(gI1,1)*ZH(gI2,1))) * tmp_606;
   std::complex<double> tmp_608;
   std::complex<double> tmp_609;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_610;
      std::complex<double> tmp_611;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_612;
         std::complex<double> tmp_613;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_613 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_612 += tmp_613;
         tmp_611 += (KroneckerDelta(gO2,3 + j2)) * tmp_612;
      }
      tmp_610 += tmp_611;
      tmp_609 += (KroneckerDelta(gO1,3 + j3)) * tmp_610;
   }
   tmp_608 += tmp_609;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_608;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,1)*ZH(gI2,1)
         );
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_614;
      std::complex<double> tmp_615;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_615 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_614 += tmp_615;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_614;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_616;
      std::complex<double> tmp_617;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_617 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_616 += tmp_617;
      result += (0.5*Conj(Lambdax)*ZH(gI1,0)*ZH(gI2,2)) * tmp_616;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_618;
      std::complex<double> tmp_619;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_619 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_618 += tmp_619;
      result += (0.5*Lambdax*ZH(gI1,0)*ZH(gI2,2)) * tmp_618;
   }
   std::complex<double> tmp_620;
   std::complex<double> tmp_621;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_621 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_620 += tmp_621;
   result += (-(Qs*Qu*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2))) * tmp_620;
   if (gO1 < 3) {
      result += -(Qq*Qs*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2))
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_622;
   std::complex<double> tmp_623;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_623 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_622 += tmp_623;
   result += (-1.4142135623730951*gp*Qu*ZN(gI2,0)) * tmp_622;
   std::complex<double> tmp_624;
   std::complex<double> tmp_625;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_625 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_624 += tmp_625;
   result += (0.7302967433402214*g1*ZN(gI2,1)) * tmp_624;
   if (gO2 < 3) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_627 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_626 += tmp_627;
      result += (-ZN(gI2,4)) * tmp_626;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_628;
   std::complex<double> tmp_629;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_630;
      std::complex<double> tmp_631;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_631 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_630 += tmp_631;
      tmp_629 += (Conj(ZUL(gI1,j2))) * tmp_630;
   }
   tmp_628 += tmp_629;
   result += (-Conj(ZN(gI2,4))) * tmp_628;
   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Qq*Conj(ZN(gI2,0))*Conj(ZUL(gI1,gO1))
         ;
   }
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZN(gI2,1))*Conj(ZUL(gI1,gO1));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,2))*Conj(ZUL(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_632;
   std::complex<double> tmp_634;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_634 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_632 += tmp_634;
   std::complex<double> tmp_633;
   std::complex<double> tmp_635;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_635 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_633 += tmp_635;
   result += (0.1*Sqr(g1)) * tmp_632 * tmp_633;
   std::complex<double> tmp_636;
   std::complex<double> tmp_638;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_638 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_636 += tmp_638;
   std::complex<double> tmp_637;
   std::complex<double> tmp_639;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_639 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_637 += tmp_639;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_636 * tmp_637;
   std::complex<double> tmp_640;
   std::complex<double> tmp_642;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_642 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_640 += tmp_642;
   std::complex<double> tmp_641;
   std::complex<double> tmp_643;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_643 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_641 += tmp_643;
   result += (0.2*Sqr(g1)) * tmp_640 * tmp_641;
   std::complex<double> tmp_644;
   std::complex<double> tmp_646;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_646 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_644 += tmp_646;
   std::complex<double> tmp_645;
   std::complex<double> tmp_647;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_647 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_645 += tmp_647;
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_644 * tmp_645;
   std::complex<double> tmp_648;
   std::complex<double> tmp_650;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_650 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_648 += tmp_650;
   std::complex<double> tmp_649;
   std::complex<double> tmp_651;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_651 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_649 += tmp_651;
   result += (0.1*Sqr(g1)) * tmp_648 * tmp_649;
   std::complex<double> tmp_652;
   std::complex<double> tmp_654;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_654 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_652 += tmp_654;
   std::complex<double> tmp_653;
   std::complex<double> tmp_655;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_655 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_653 += tmp_655;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_652 * tmp_653;
   std::complex<double> tmp_656;
   std::complex<double> tmp_658;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_658 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_656 += tmp_658;
   std::complex<double> tmp_657;
   std::complex<double> tmp_659;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_659 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_657 += tmp_659;
   result += (0.2*Sqr(g1)) * tmp_656 * tmp_657;
   std::complex<double> tmp_660;
   std::complex<double> tmp_662;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_662 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_660 += tmp_662;
   std::complex<double> tmp_661;
   std::complex<double> tmp_663;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_663 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_661 += tmp_663;
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_660 * tmp_661;
   std::complex<double> tmp_664;
   std::complex<double> tmp_666;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_667;
      std::complex<double> tmp_668;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_668 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_667 += tmp_668;
      tmp_666 += (Conj(ZD(gI2,j2))) * tmp_667;
   }
   tmp_664 += tmp_666;
   std::complex<double> tmp_665;
   std::complex<double> tmp_669;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_670;
      std::complex<double> tmp_671;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_671 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_670 += tmp_671;
      tmp_669 += (ZD(gI1,j4)) * tmp_670;
   }
   tmp_665 += tmp_669;
   result += (-1) * tmp_664 * tmp_665;
   if (gO1 < 3) {
      std::complex<double> tmp_672;
      std::complex<double> tmp_673;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_673 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_672 += tmp_673;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_672;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_674;
      std::complex<double> tmp_675;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_675 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_674 += tmp_675;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_674;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_676;
      std::complex<double> tmp_677;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_677 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_676 += tmp_677;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_676;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_678;
      std::complex<double> tmp_679;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_679 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_678 += tmp_679;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_678;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_680;
      std::complex<double> tmp_681;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_681 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_680 += tmp_681;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_680;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_682;
      std::complex<double> tmp_683;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_683 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_682 += tmp_683;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_682;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_684;
      std::complex<double> tmp_685;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_685 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_684 += tmp_685;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_684;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_686;
      std::complex<double> tmp_687;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_687 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_686 += tmp_687;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_686;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_688;
      std::complex<double> tmp_689;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_689 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_688 += tmp_689;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_688;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_690;
      std::complex<double> tmp_691;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_691 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_690 += tmp_691;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_690;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_692;
      std::complex<double> tmp_694;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_694 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_692 += tmp_694;
      std::complex<double> tmp_693;
      std::complex<double> tmp_695;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_695 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_693 += tmp_695;
      result += (-1) * tmp_692 * tmp_693;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_696;
   std::complex<double> tmp_698;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_698 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_696 += tmp_698;
   std::complex<double> tmp_697;
   std::complex<double> tmp_699;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_699 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_697 += tmp_699;
   result += (-0.1*Sqr(g1)) * tmp_696 * tmp_697;
   std::complex<double> tmp_700;
   std::complex<double> tmp_702;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_702 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_700 += tmp_702;
   std::complex<double> tmp_701;
   std::complex<double> tmp_703;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_703 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_701 += tmp_703;
   result += (-0.5*Ql*Qu*Sqr(gp)) * tmp_700 * tmp_701;
   std::complex<double> tmp_704;
   std::complex<double> tmp_706;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_706 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_704 += tmp_706;
   std::complex<double> tmp_705;
   std::complex<double> tmp_707;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_707 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_705 += tmp_707;
   result += (0.2*Sqr(g1)) * tmp_704 * tmp_705;
   std::complex<double> tmp_708;
   std::complex<double> tmp_710;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_710 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_708 += tmp_710;
   std::complex<double> tmp_709;
   std::complex<double> tmp_711;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_711 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_709 += tmp_711;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_708 * tmp_709;
   std::complex<double> tmp_712;
   std::complex<double> tmp_714;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_714 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_712 += tmp_714;
   std::complex<double> tmp_713;
   std::complex<double> tmp_715;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_715 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_713 += tmp_715;
   result += (-0.1*Sqr(g1)) * tmp_712 * tmp_713;
   std::complex<double> tmp_716;
   std::complex<double> tmp_718;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_718 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_716 += tmp_718;
   std::complex<double> tmp_717;
   std::complex<double> tmp_719;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_719 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_717 += tmp_719;
   result += (-0.5*Ql*Qu*Sqr(gp)) * tmp_716 * tmp_717;
   std::complex<double> tmp_720;
   std::complex<double> tmp_722;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_722 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_720 += tmp_722;
   std::complex<double> tmp_721;
   std::complex<double> tmp_723;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_723 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_721 += tmp_723;
   result += (0.2*Sqr(g1)) * tmp_720 * tmp_721;
   std::complex<double> tmp_724;
   std::complex<double> tmp_726;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_726 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_724 += tmp_726;
   std::complex<double> tmp_725;
   std::complex<double> tmp_727;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_727 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_725 += tmp_727;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_724 * tmp_725;
   if (gO1 < 3) {
      std::complex<double> tmp_728;
      std::complex<double> tmp_729;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_729 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_728 += tmp_729;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_728;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_730;
      std::complex<double> tmp_731;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_731 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_730 += tmp_731;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_730;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_732;
      std::complex<double> tmp_733;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_733 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_732 += tmp_733;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_732;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_734;
      std::complex<double> tmp_735;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_735 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_734 += tmp_735;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_734;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_736;
      std::complex<double> tmp_737;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_737 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_736 += tmp_737;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_736;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_738;
      std::complex<double> tmp_739;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_739 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_738 += tmp_739;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_738;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_740;
      std::complex<double> tmp_741;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_741 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_740 += tmp_741;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_740;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_742;
      std::complex<double> tmp_743;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_743 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_742 += tmp_743;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_742;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_744;
      std::complex<double> tmp_745;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_745 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_744 += tmp_745;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_744;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_746;
      std::complex<double> tmp_747;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_747 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_746 += tmp_747;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_746;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_748;
   std::complex<double> tmp_750;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_750 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_748 += tmp_750;
   std::complex<double> tmp_749;
   std::complex<double> tmp_751;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_751 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_749 += tmp_751;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_748 * tmp_749;
   std::complex<double> tmp_752;
   std::complex<double> tmp_754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_754 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_752 += tmp_754;
   std::complex<double> tmp_753;
   std::complex<double> tmp_755;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_755 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_753 += tmp_755;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_752 * tmp_753;
   std::complex<double> tmp_756;
   std::complex<double> tmp_758;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_758 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_756 += tmp_758;
   std::complex<double> tmp_757;
   std::complex<double> tmp_759;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_759 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_757 += tmp_759;
   result += (-0.5*Sqr(gp)*Sqr(Qu)) * tmp_756 * tmp_757;
   std::complex<double> tmp_760;
   std::complex<double> tmp_762;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_762 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_760 += tmp_762;
   std::complex<double> tmp_761;
   std::complex<double> tmp_763;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_763 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_761 += tmp_763;
   result += (0.1*Sqr(g1)) * tmp_760 * tmp_761;
   std::complex<double> tmp_764;
   std::complex<double> tmp_766;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_766 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_764 += tmp_766;
   std::complex<double> tmp_765;
   std::complex<double> tmp_767;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_767 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_765 += tmp_767;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_764 * tmp_765;
   std::complex<double> tmp_768;
   std::complex<double> tmp_770;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_770 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_768 += tmp_770;
   std::complex<double> tmp_769;
   std::complex<double> tmp_771;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_771 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_769 += tmp_771;
   result += (-0.4*Sqr(g1)) * tmp_768 * tmp_769;
   std::complex<double> tmp_772;
   std::complex<double> tmp_774;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_774 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_772 += tmp_774;
   std::complex<double> tmp_773;
   std::complex<double> tmp_775;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_775 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_773 += tmp_775;
   result += (-1.5*Sqr(gp)*Sqr(Qu)) * tmp_772 * tmp_773;
   std::complex<double> tmp_776;
   std::complex<double> tmp_778;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_778 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_776 += tmp_778;
   std::complex<double> tmp_777;
   std::complex<double> tmp_779;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_779 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_777 += tmp_779;
   result += (0.1*Sqr(g1)) * tmp_776 * tmp_777;
   std::complex<double> tmp_780;
   std::complex<double> tmp_782;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_782 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_780 += tmp_782;
   std::complex<double> tmp_781;
   std::complex<double> tmp_783;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_783 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_781 += tmp_783;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_780 * tmp_781;
   std::complex<double> tmp_784;
   std::complex<double> tmp_786;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_786 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_784 += tmp_786;
   std::complex<double> tmp_785;
   std::complex<double> tmp_787;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_787 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_785 += tmp_787;
   result += (-0.4*Sqr(g1)) * tmp_784 * tmp_785;
   std::complex<double> tmp_788;
   std::complex<double> tmp_790;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_790 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_788 += tmp_790;
   std::complex<double> tmp_789;
   std::complex<double> tmp_791;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_791 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_789 += tmp_791;
   result += (-1.5*Sqr(gp)*Sqr(Qu)) * tmp_788 * tmp_789;
   std::complex<double> tmp_792;
   std::complex<double> tmp_794;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_794 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_792 += tmp_794;
   std::complex<double> tmp_793;
   std::complex<double> tmp_795;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_795 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_793 += tmp_795;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_792 * tmp_793;
   std::complex<double> tmp_796;
   std::complex<double> tmp_798;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_798 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_796 += tmp_798;
   std::complex<double> tmp_797;
   std::complex<double> tmp_799;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_799 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_797 += tmp_799;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_796 * tmp_797;
   std::complex<double> tmp_800;
   std::complex<double> tmp_802;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_802 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_800 += tmp_802;
   std::complex<double> tmp_801;
   std::complex<double> tmp_803;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_803 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_801 += tmp_803;
   result += (-0.5*Sqr(gp)*Sqr(Qu)) * tmp_800 * tmp_801;
   std::complex<double> tmp_804;
   std::complex<double> tmp_806;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_807;
      std::complex<double> tmp_808;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_808 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_807 += tmp_808;
      tmp_806 += (Conj(ZU(gI2,j2))) * tmp_807;
   }
   tmp_804 += tmp_806;
   std::complex<double> tmp_805;
   std::complex<double> tmp_809;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_810;
      std::complex<double> tmp_811;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_811 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_810 += tmp_811;
      tmp_809 += (ZU(gI1,j4)) * tmp_810;
   }
   tmp_805 += tmp_809;
   result += (-1) * tmp_804 * tmp_805;
   if (gO1 < 3) {
      std::complex<double> tmp_812;
      std::complex<double> tmp_813;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_813 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_812 += tmp_813;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_812;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_814;
      std::complex<double> tmp_815;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_815 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_814 += tmp_815;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_814;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_816;
      std::complex<double> tmp_817;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_817 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_816 += tmp_817;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_816;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_818;
      std::complex<double> tmp_819;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_819 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_818 += tmp_819;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_818;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_820;
      std::complex<double> tmp_821;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_821 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_820 += tmp_821;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_820;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_822;
      std::complex<double> tmp_823;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_823 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_822 += tmp_823;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_822;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_824;
      std::complex<double> tmp_825;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_825 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_824 += tmp_825;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_824;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_826;
      std::complex<double> tmp_827;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_827 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_826 += tmp_827;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_826;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_828;
      std::complex<double> tmp_829;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_829 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_828 += tmp_829;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_828;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_830;
      std::complex<double> tmp_831;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_831 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_830 += tmp_831;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_830;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_832;
      std::complex<double> tmp_834;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_834 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_832 += tmp_834;
      std::complex<double> tmp_833;
      std::complex<double> tmp_835;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_836;
         std::complex<double> tmp_837;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_837 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_836 += tmp_837;
         tmp_835 += (ZU(gI1,j4)) * tmp_836;
      }
      tmp_833 += tmp_835;
      result += (-3) * tmp_832 * tmp_833;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_838;
      std::complex<double> tmp_839;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_839 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_838 += tmp_839;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_838;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_840;
      std::complex<double> tmp_841;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_841 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_840 += tmp_841;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_840;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_842;
      std::complex<double> tmp_843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_843 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_842 += tmp_843;
      result += (-0.5*Qq*Qu*Conj(ZU(gI2,gO2))*Sqr(gp)) * tmp_842;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_844;
      std::complex<double> tmp_845;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_845 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_844 += tmp_845;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_844;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_846;
      std::complex<double> tmp_847;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_847 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_846 += tmp_847;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_846;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_848;
      std::complex<double> tmp_849;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_849 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_848 += tmp_849;
      result += (-0.5*Qq*Qu*Conj(ZU(gI2,gO2))*Sqr(gp)) * tmp_848;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_850;
      std::complex<double> tmp_852;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_853;
         std::complex<double> tmp_854;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_854 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_853 += tmp_854;
         tmp_852 += (Conj(ZU(gI2,j2))) * tmp_853;
      }
      tmp_850 += tmp_852;
      std::complex<double> tmp_851;
      std::complex<double> tmp_855;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_855 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_851 += tmp_855;
      result += (-3) * tmp_850 * tmp_851;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_856;
      std::complex<double> tmp_858;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_858 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_856 += tmp_858;
      std::complex<double> tmp_857;
      std::complex<double> tmp_859;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_859 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_857 += tmp_859;
      result += (-1) * tmp_856 * tmp_857;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_860;
      std::complex<double> tmp_861;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_861 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_860 += tmp_861;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_860;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_862;
      std::complex<double> tmp_863;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_863 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_862 += tmp_863;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_862;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_864;
      std::complex<double> tmp_865;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_865 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_864 += tmp_865;
      result += (-0.5*Qq*Qu*Sqr(gp)*ZU(gI1,gO1)) * tmp_864;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_866;
      std::complex<double> tmp_867;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_867 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_866 += tmp_867;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_866;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_868;
      std::complex<double> tmp_869;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_869 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_868 += tmp_869;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_868;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_870;
      std::complex<double> tmp_871;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_871 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_870 += tmp_871;
      result += (-0.5*Qq*Qu*Sqr(gp)*ZU(gI1,gO1)) * tmp_870;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.016666666666666666*Conj(ZU(gI2,gO2))*Sqr(g1)*ZU(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -1.3333333333333333*Conj(ZU(gI2,gO2))*Sqr(g3)*ZU(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -(Conj(ZU(gI2,gO2))*Sqr(gp)*Sqr(Qq)*ZU(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_872;
   std::complex<double> tmp_873;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_874;
      std::complex<double> tmp_875;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_875 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_874 += tmp_875;
      tmp_873 += (Conj(ZU(gI1,j2))) * tmp_874;
   }
   tmp_872 += tmp_873;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,0))) *
      tmp_872;
   std::complex<double> tmp_876;
   std::complex<double> tmp_877;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_878;
      std::complex<double> tmp_879;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_879 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_878 += tmp_879;
      tmp_877 += (Conj(ZU(gI1,j2))) * tmp_878;
   }
   tmp_876 += tmp_877;
   result += (std::complex<double>(0,-0.5)*vd*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_876;
   std::complex<double> tmp_880;
   std::complex<double> tmp_881;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_882;
      std::complex<double> tmp_883;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_883 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_882 += tmp_883;
      tmp_881 += (Conj(ZU(gI1,j2))) * tmp_882;
   }
   tmp_880 += tmp_881;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_880;
   if (gO2 < 3) {
      std::complex<double> tmp_884;
      std::complex<double> tmp_885;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_885 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_884 += tmp_885;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,0))*Lambdax) *
         tmp_884;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_886;
      std::complex<double> tmp_887;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_887 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_886 += tmp_887;
      result += (std::complex<double>(0,0.5)*vd*Conj(ZA(gI2,2))*Lambdax) *
         tmp_886;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_888;
      std::complex<double> tmp_889;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_889 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_888 += tmp_889;
      result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,1)))
         * tmp_888;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   std::complex<double> tmp_890;
   std::complex<double> tmp_891;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_891 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_890 += tmp_891;
   result += (-0.2*vd*Sqr(g1)*ZH(gI2,0)) * tmp_890;
   std::complex<double> tmp_892;
   std::complex<double> tmp_893;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_893 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_892 += tmp_893;
   result += (-(QHd*Qu*vd*Sqr(gp)*ZH(gI2,0))) * tmp_892;
   std::complex<double> tmp_894;
   std::complex<double> tmp_895;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_896;
      std::complex<double> tmp_897;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_897 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_896 += tmp_897;
      tmp_895 += (Conj(ZU(gI1,j2))) * tmp_896;
   }
   tmp_894 += tmp_895;
   result += (0.5*vS*Conj(Lambdax)*ZH(gI2,0)) * tmp_894;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -(QHd*Qq*vd*Conj(ZU(gI1,gO2))*Sqr(gp)*ZH(gI2,0));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_898;
      std::complex<double> tmp_899;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_899 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_898 += tmp_899;
      result += (0.5*vS*Lambdax*ZH(gI2,0)) * tmp_898;
   }
   std::complex<double> tmp_900;
   std::complex<double> tmp_901;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_901 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_900 += tmp_901;
   result += (0.2*vu*Sqr(g1)*ZH(gI2,1)) * tmp_900;
   std::complex<double> tmp_902;
   std::complex<double> tmp_903;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_903 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_902 += tmp_903;
   result += (-(QHu*Qu*vu*Sqr(gp)*ZH(gI2,1))) * tmp_902;
   std::complex<double> tmp_904;
   std::complex<double> tmp_905;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_906;
      std::complex<double> tmp_907;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_907 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_906 += tmp_907;
      tmp_905 += (Conj(ZU(gI1,j2))) * tmp_906;
   }
   tmp_904 += tmp_905;
   result += (-0.7071067811865475*ZH(gI2,1)) * tmp_904;
   std::complex<double> tmp_908;
   std::complex<double> tmp_909;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_910;
      std::complex<double> tmp_911;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_912;
         std::complex<double> tmp_913;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_913 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_912 += tmp_913;
         tmp_911 += (KroneckerDelta(gO2,3 + j2)) * tmp_912;
      }
      tmp_910 += tmp_911;
      tmp_909 += (Conj(ZU(gI1,3 + j3))) * tmp_910;
   }
   tmp_908 += tmp_909;
   result += (-(vu*ZH(gI2,1))) * tmp_908;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -(QHu*Qq*vu*Conj(ZU(gI1,gO2))*Sqr(gp)*ZH(gI2,1));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_914;
      std::complex<double> tmp_915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_915 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_914 += tmp_915;
      result += (-0.7071067811865475*ZH(gI2,1)) * tmp_914;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_916;
      std::complex<double> tmp_917;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_918;
         std::complex<double> tmp_919;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_919 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_918 += tmp_919;
         tmp_917 += (Conj(ZU(gI1,j2))) * tmp_918;
      }
      tmp_916 += tmp_917;
      result += (-(vu*ZH(gI2,1))) * tmp_916;
   }
   std::complex<double> tmp_920;
   std::complex<double> tmp_921;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_921 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_920 += tmp_921;
   result += (-(Qs*Qu*vS*Sqr(gp)*ZH(gI2,2))) * tmp_920;
   std::complex<double> tmp_922;
   std::complex<double> tmp_923;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_924;
      std::complex<double> tmp_925;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_925 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_924 += tmp_925;
      tmp_923 += (Conj(ZU(gI1,j2))) * tmp_924;
   }
   tmp_922 += tmp_923;
   result += (0.5*vd*Conj(Lambdax)*ZH(gI2,2)) * tmp_922;
   if (gO2 < 3) {
      result += -(Qq*Qs*vS*Conj(ZU(gI1,gO2))*Sqr(gp)*ZH(gI2,2));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_926;
      std::complex<double> tmp_927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_927 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_926 += tmp_927;
      result += (0.5*vd*Lambdax*ZH(gI2,2)) * tmp_926;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_928;
   std::complex<double> tmp_929;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_929 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_928 += tmp_929;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_928;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPL(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*PhaseGlu*Conj(ZUL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjVWmSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZD(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVGSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 6) {
      result += g3*Conj(ZU(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVPSu(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_930;
   std::complex<double> tmp_931;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_931 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_930 += tmp_931;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_930;
   if (gO2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZU(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVZSu(unsigned gO2, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_932;
   std::complex<double> tmp_933;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_933 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_932 += tmp_933;
   result += (-0.5163977794943222*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_932;
   std::complex<double> tmp_934;
   std::complex<double> tmp_935;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_935 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_934 += tmp_935;
   result += (-(gp*Qu*Sin(ThetaWp()))) * tmp_934;
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZU(gI2,gO2))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Cos(ThetaWp())*Sin
         (ThetaW());
   }
   if (gO2 < 3) {
      result += gp*Qq*Conj(ZU(gI2,gO2))*Sin(ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuVZpSu(unsigned gO2, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_936;
   std::complex<double> tmp_937;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_937 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_936 += tmp_937;
   result += (-(gp*Qu*Cos(ThetaWp()))) * tmp_936;
   std::complex<double> tmp_938;
   std::complex<double> tmp_939;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_939 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_938 += tmp_939;
   result += (0.5163977794943222*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_938;
   if (gO2 < 3) {
      result += gp*Qq*Conj(ZU(gI2,gO2))*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZU(gI2,gO2))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gO2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_940;
   std::complex<double> tmp_941;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_941 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_940 += tmp_941;
   result += (1.2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_940;
   std::complex<double> tmp_942;
   std::complex<double> tmp_943;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_943 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_942 += tmp_943;
   result += (-3.0983866769659336*g1*gp*Qe*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_942;
   std::complex<double> tmp_944;
   std::complex<double> tmp_945;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_945 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_944 += tmp_945;
   result += (2*Sqr(gp)*Sqr(Qe)*Sqr(Sin(ThetaWp()))) * tmp_944;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += -0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(
         Sin(ThetaW()));
   }
   if (gO1 < 3) {
      result += -2*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 1.5491933384829668*g1*gp*Ql*Cos(ThetaWp())*KroneckerDelta(
         gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_946;
   std::complex<double> tmp_947;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_947 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_946 += tmp_947;
   result += (2*Sqr(gp)*Sqr(Qe)*Sqr(Cos(ThetaWp()))) * tmp_946;
   std::complex<double> tmp_948;
   std::complex<double> tmp_949;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_949 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_948 += tmp_949;
   result += (3.0983866769659336*g1*gp*Qe*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_948;
   std::complex<double> tmp_950;
   std::complex<double> tmp_951;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_951 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_950 += tmp_951;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_950;
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*Sqr(Cos(ThetaWp())
         );
   }
   if (gO1 < 3) {
      result += 2*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += -1.5491933384829668*g1*gp*Ql*Cos(ThetaWp())*KroneckerDelta(
         gO1,gO2)*Sin(ThetaW())*Sin(ThetaWp());
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += -0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
         Sin(ThetaWp()));
   }

   return result;
}

double CLASSNAME::CpUSeconjUSeconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   std::complex<double> tmp_952;
   std::complex<double> tmp_953;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_953 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_952 += tmp_953;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_952;
   std::complex<double> tmp_954;
   std::complex<double> tmp_955;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_955 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_954 += tmp_955;
   result += (-(Qe*QHd*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0))) * tmp_954;
   std::complex<double> tmp_956;
   std::complex<double> tmp_957;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_958;
      std::complex<double> tmp_959;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_960;
         std::complex<double> tmp_961;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_961 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_960 += tmp_961;
         tmp_959 += (KroneckerDelta(gO2,3 + j2)) * tmp_960;
      }
      tmp_958 += tmp_959;
      tmp_957 += (KroneckerDelta(gO1,3 + j3)) * tmp_958;
   }
   tmp_956 += tmp_957;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_956;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)
         );
   }
   std::complex<double> tmp_962;
   std::complex<double> tmp_963;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_963 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_962 += tmp_963;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_962;
   std::complex<double> tmp_964;
   std::complex<double> tmp_965;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_965 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_964 += tmp_965;
   result += (-(Qe*QHu*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1))) * tmp_964;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_966;
   std::complex<double> tmp_967;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_967 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_966 += tmp_967;
   result += (0.3*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)) * tmp_966;
   std::complex<double> tmp_968;
   std::complex<double> tmp_969;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_969 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_968 += tmp_969;
   result += (-(Qe*QHd*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(gp))) * tmp_968;
   std::complex<double> tmp_970;
   std::complex<double> tmp_971;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_971 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_970 += tmp_971;
   result += (-0.3*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)) * tmp_970;
   std::complex<double> tmp_972;
   std::complex<double> tmp_973;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_973 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_972 += tmp_973;
   result += (-(Qe*QHu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp))) * tmp_972;
   std::complex<double> tmp_974;
   std::complex<double> tmp_975;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_975 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_974 += tmp_975;
   result += (-(Qe*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp))) * tmp_974;
   std::complex<double> tmp_976;
   std::complex<double> tmp_977;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_978;
      std::complex<double> tmp_979;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_980;
         std::complex<double> tmp_981;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_981 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_980 += tmp_981;
         tmp_979 += (KroneckerDelta(gO2,3 + j2)) * tmp_980;
      }
      tmp_978 += tmp_979;
      tmp_977 += (KroneckerDelta(gO1,3 + j3)) * tmp_978;
   }
   tmp_976 += tmp_977;
   result += (-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)))) * tmp_976;
   if (gO1 < 3) {
      result += -0.15*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Ql*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += 0.15*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Ql*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Ql*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_982;
      std::complex<double> tmp_983;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_983 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_982 += tmp_983;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))) *
         tmp_982;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_984;
      std::complex<double> tmp_985;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_985 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_984 += tmp_985;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))) *
         tmp_984;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_986;
      std::complex<double> tmp_987;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_987 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_986 += tmp_987;
      result += (-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*Lambdax) * tmp_986;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_988;
      std::complex<double> tmp_989;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_989 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_988 += tmp_989;
      result += (-0.5*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))*Lambdax) * tmp_988;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_990;
      std::complex<double> tmp_991;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_991 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_990 += tmp_991;
      result += (-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)))) * tmp_990;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_992;
   std::complex<double> tmp_993;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_993 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_992 += tmp_993;
   result += (0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_992;
   std::complex<double> tmp_994;
   std::complex<double> tmp_995;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_995 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_994 += tmp_995;
   result += (-(Qe*Ql*KroneckerDelta(gI1,gI2)*Sqr(gp))) * tmp_994;
   std::complex<double> tmp_996;
   std::complex<double> tmp_998;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_999;
      std::complex<double> tmp_1000;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1000 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_999 += tmp_1000;
      tmp_998 += (Conj(ZV(gI2,j2))) * tmp_999;
   }
   tmp_996 += tmp_998;
   std::complex<double> tmp_997;
   std::complex<double> tmp_1001;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1002;
      std::complex<double> tmp_1003;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1003 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1002 += tmp_1003;
      tmp_1001 += (ZV(gI1,j4)) * tmp_1002;
   }
   tmp_997 += tmp_1001;
   result += (-1) * tmp_996 * tmp_997;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1
         );
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2)
         ;
   }
   if (gO1 < 3) {
      result += -(KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(gp)*
         Sqr(Ql));
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   std::complex<double> tmp_1004;
   std::complex<double> tmp_1005;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1005 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1004 += tmp_1005;
   result += (0.3*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_1004;
   std::complex<double> tmp_1006;
   std::complex<double> tmp_1007;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1007 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1006 += tmp_1007;
   result += (-(Qe*QHd*Sqr(gp)*ZH(gI1,0)*ZH(gI2,0))) * tmp_1006;
   std::complex<double> tmp_1008;
   std::complex<double> tmp_1009;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1010;
      std::complex<double> tmp_1011;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1012;
         std::complex<double> tmp_1013;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1013 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1012 += tmp_1013;
         tmp_1011 += (KroneckerDelta(gO2,3 + j2)) * tmp_1012;
      }
      tmp_1010 += tmp_1011;
      tmp_1009 += (KroneckerDelta(gO1,3 + j3)) * tmp_1010;
   }
   tmp_1008 += tmp_1009;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_1008;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,0)*ZH(gI2,0)
         );
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1014;
      std::complex<double> tmp_1015;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1015 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_1014 += tmp_1015;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_1014;
   }
   std::complex<double> tmp_1016;
   std::complex<double> tmp_1017;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1017 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1016 += tmp_1017;
   result += (-0.3*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_1016;
   std::complex<double> tmp_1018;
   std::complex<double> tmp_1019;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1019 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1018 += tmp_1019;
   result += (-(Qe*QHu*Sqr(gp)*ZH(gI1,1)*ZH(gI2,1))) * tmp_1018;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,1)*ZH(gI2,1)
         );
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1020;
      std::complex<double> tmp_1021;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1021 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1020 += tmp_1021;
      result += (0.5*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,1)) * tmp_1020;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1022;
      std::complex<double> tmp_1023;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1023 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1022 += tmp_1023;
      result += (0.5*Lambdax*ZH(gI1,2)*ZH(gI2,1)) * tmp_1022;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1024;
      std::complex<double> tmp_1025;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1025 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1024 += tmp_1025;
      result += (0.5*Conj(Lambdax)*ZH(gI1,1)*ZH(gI2,2)) * tmp_1024;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1026;
      std::complex<double> tmp_1027;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1027 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1026 += tmp_1027;
      result += (0.5*Lambdax*ZH(gI1,1)*ZH(gI2,2)) * tmp_1026;
   }
   std::complex<double> tmp_1028;
   std::complex<double> tmp_1029;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1029 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1028 += tmp_1029;
   result += (-(Qe*Qs*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2))) * tmp_1028;
   if (gO1 < 3) {
      result += -(Ql*Qs*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZH(gI1,2)*ZH(gI2,2))
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1030;
   std::complex<double> tmp_1031;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1032;
      std::complex<double> tmp_1033;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1033 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_1032 += tmp_1033;
      tmp_1031 += (Conj(ZV(gI1,j2))) * tmp_1032;
   }
   tmp_1030 += tmp_1031;
   result += (ZP(gI2,0)) * tmp_1030;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1034;
      std::complex<double> tmp_1035;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1036;
         std::complex<double> tmp_1037;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1037 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_1036 += tmp_1037;
         tmp_1035 += (Conj(ZV(gI1,j2))) * tmp_1036;
      }
      tmp_1034 += tmp_1035;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_1034;
   }
   std::complex<double> tmp_1038;
   std::complex<double> tmp_1039;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1040;
      std::complex<double> tmp_1041;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1041 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1040 += tmp_1041;
      tmp_1039 += (Conj(ZV(gI1,j2))) * tmp_1040;
   }
   tmp_1038 += tmp_1039;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI2,1)) * tmp_1038;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }

   return result;
}

double CLASSNAME::CpconjUSeFvChaPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1042;
   std::complex<double> tmp_1043;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1043 += KroneckerDelta(gO1,3 + j1)*Ye(j1,gI1);
   }
   tmp_1042 += tmp_1043;
   result += (Conj(UM(gI2,1))) * tmp_1042;
   if (gI1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*KroneckerDelta(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_1044;
   std::complex<double> tmp_1045;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1045 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_1044 += tmp_1045;
   result += (-1.4142135623730951*gp*Qe*ZN(gI2,0)) * tmp_1044;
   std::complex<double> tmp_1046;
   std::complex<double> tmp_1047;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1047 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_1046 += tmp_1047;
   result += (-1.0954451150103321*g1*ZN(gI2,1)) * tmp_1046;
   if (gO2 < 3) {
      std::complex<double> tmp_1048;
      std::complex<double> tmp_1049;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1049 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_1048 += tmp_1049;
      result += (-ZN(gI2,3)) * tmp_1048;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1050;
   std::complex<double> tmp_1051;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1052;
      std::complex<double> tmp_1053;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1053 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_1052 += tmp_1053;
      tmp_1051 += (Conj(ZEL(gI1,j2))) * tmp_1052;
   }
   tmp_1050 += tmp_1051;
   result += (-Conj(ZN(gI2,3))) * tmp_1050;
   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Ql*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,0))
         ;
   }
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,1));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1054;
   std::complex<double> tmp_1056;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1056 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1054 += tmp_1056;
   std::complex<double> tmp_1055;
   std::complex<double> tmp_1057;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1057 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1055 += tmp_1057;
   result += (-0.05*Sqr(g1)) * tmp_1054 * tmp_1055;
   std::complex<double> tmp_1058;
   std::complex<double> tmp_1060;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1060 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1058 += tmp_1060;
   std::complex<double> tmp_1059;
   std::complex<double> tmp_1061;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1061 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1059 += tmp_1061;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1058 * tmp_1059;
   std::complex<double> tmp_1062;
   std::complex<double> tmp_1064;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1064 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1062 += tmp_1064;
   std::complex<double> tmp_1063;
   std::complex<double> tmp_1065;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1065 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1063 += tmp_1065;
   result += (-0.1*Sqr(g1)) * tmp_1062 * tmp_1063;
   std::complex<double> tmp_1066;
   std::complex<double> tmp_1068;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1068 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1066 += tmp_1068;
   std::complex<double> tmp_1067;
   std::complex<double> tmp_1069;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1069 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1067 += tmp_1069;
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_1066 * tmp_1067;
   std::complex<double> tmp_1070;
   std::complex<double> tmp_1072;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1072 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1070 += tmp_1072;
   std::complex<double> tmp_1071;
   std::complex<double> tmp_1073;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1073 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_1071 += tmp_1073;
   result += (-0.05*Sqr(g1)) * tmp_1070 * tmp_1071;
   std::complex<double> tmp_1074;
   std::complex<double> tmp_1076;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1076 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1074 += tmp_1076;
   std::complex<double> tmp_1075;
   std::complex<double> tmp_1077;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1077 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_1075 += tmp_1077;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1074 * tmp_1075;
   std::complex<double> tmp_1078;
   std::complex<double> tmp_1080;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1080 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1078 += tmp_1080;
   std::complex<double> tmp_1079;
   std::complex<double> tmp_1081;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1081 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_1079 += tmp_1081;
   result += (-0.1*Sqr(g1)) * tmp_1078 * tmp_1079;
   std::complex<double> tmp_1082;
   std::complex<double> tmp_1084;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1084 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1082 += tmp_1084;
   std::complex<double> tmp_1083;
   std::complex<double> tmp_1085;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1085 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_1083 += tmp_1085;
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_1082 * tmp_1083;
   if (gO1 < 3) {
      std::complex<double> tmp_1086;
      std::complex<double> tmp_1087;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1087 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1086 += tmp_1087;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1086;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1088;
      std::complex<double> tmp_1089;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1089 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1088 += tmp_1089;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1088;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1090;
      std::complex<double> tmp_1091;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1091 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1090 += tmp_1091;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1090;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1092;
      std::complex<double> tmp_1093;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1093 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_1092 += tmp_1093;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1092;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1094;
      std::complex<double> tmp_1095;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1095 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_1094 += tmp_1095;
      result += (-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1094;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1096;
      std::complex<double> tmp_1097;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1097 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1096 += tmp_1097;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1096;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1098;
      std::complex<double> tmp_1099;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1099 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1098 += tmp_1099;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1098;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1100;
      std::complex<double> tmp_1101;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1101 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1100 += tmp_1101;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1100;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1102;
      std::complex<double> tmp_1103;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1103 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_1102 += tmp_1103;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1102;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1104;
      std::complex<double> tmp_1105;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1105 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_1104 += tmp_1105;
      result += (-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1104;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1106;
      std::complex<double> tmp_1108;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1108 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1106 += tmp_1108;
      std::complex<double> tmp_1107;
      std::complex<double> tmp_1109;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1110;
         std::complex<double> tmp_1111;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1111 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_1110 += tmp_1111;
         tmp_1109 += (ZD(gI1,j4)) * tmp_1110;
      }
      tmp_1107 += tmp_1109;
      result += (-1) * tmp_1106 * tmp_1107;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1112;
      std::complex<double> tmp_1114;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1115;
         std::complex<double> tmp_1116;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1116 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_1115 += tmp_1116;
         tmp_1114 += (Conj(ZD(gI2,j2))) * tmp_1115;
      }
      tmp_1112 += tmp_1114;
      std::complex<double> tmp_1113;
      std::complex<double> tmp_1117;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1117 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1113 += tmp_1117;
      result += (-1) * tmp_1112 * tmp_1113;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1118;
   std::complex<double> tmp_1120;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1120 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_1118 += tmp_1120;
   std::complex<double> tmp_1119;
   std::complex<double> tmp_1121;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1121 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1119 += tmp_1121;
   result += (-0.3*Sqr(g1)) * tmp_1118 * tmp_1119;
   std::complex<double> tmp_1122;
   std::complex<double> tmp_1124;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1124 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_1122 += tmp_1124;
   std::complex<double> tmp_1123;
   std::complex<double> tmp_1125;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1125 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1123 += tmp_1125;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1122 * tmp_1123;
   std::complex<double> tmp_1126;
   std::complex<double> tmp_1128;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1128 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1126 += tmp_1128;
   std::complex<double> tmp_1127;
   std::complex<double> tmp_1129;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1129 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1127 += tmp_1129;
   result += (0.15*Sqr(g1)) * tmp_1126 * tmp_1127;
   std::complex<double> tmp_1130;
   std::complex<double> tmp_1132;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1132 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1130 += tmp_1132;
   std::complex<double> tmp_1131;
   std::complex<double> tmp_1133;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1133 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1131 += tmp_1133;
   result += (-0.5*Qe*Ql*Sqr(gp)) * tmp_1130 * tmp_1131;
   std::complex<double> tmp_1134;
   std::complex<double> tmp_1136;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1136 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1134 += tmp_1136;
   std::complex<double> tmp_1135;
   std::complex<double> tmp_1137;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1137 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1135 += tmp_1137;
   result += (-0.3*Sqr(g1)) * tmp_1134 * tmp_1135;
   std::complex<double> tmp_1138;
   std::complex<double> tmp_1140;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1140 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1138 += tmp_1140;
   std::complex<double> tmp_1139;
   std::complex<double> tmp_1141;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1141 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1139 += tmp_1141;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1138 * tmp_1139;
   std::complex<double> tmp_1142;
   std::complex<double> tmp_1144;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1144 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1142 += tmp_1144;
   std::complex<double> tmp_1143;
   std::complex<double> tmp_1145;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1145 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_1143 += tmp_1145;
   result += (0.15*Sqr(g1)) * tmp_1142 * tmp_1143;
   std::complex<double> tmp_1146;
   std::complex<double> tmp_1148;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1148 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1146 += tmp_1148;
   std::complex<double> tmp_1147;
   std::complex<double> tmp_1149;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1149 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_1147 += tmp_1149;
   result += (-0.5*Qe*Ql*Sqr(gp)) * tmp_1146 * tmp_1147;
   std::complex<double> tmp_1150;
   std::complex<double> tmp_1152;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1152 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1150 += tmp_1152;
   std::complex<double> tmp_1151;
   std::complex<double> tmp_1153;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1153 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_1151 += tmp_1153;
   result += (-0.3*Sqr(g1)) * tmp_1150 * tmp_1151;
   std::complex<double> tmp_1154;
   std::complex<double> tmp_1156;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1156 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1154 += tmp_1156;
   std::complex<double> tmp_1155;
   std::complex<double> tmp_1157;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1157 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_1155 += tmp_1157;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1154 * tmp_1155;
   std::complex<double> tmp_1158;
   std::complex<double> tmp_1160;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1160 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1158 += tmp_1160;
   std::complex<double> tmp_1159;
   std::complex<double> tmp_1161;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1161 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_1159 += tmp_1161;
   result += (-0.3*Sqr(g1)) * tmp_1158 * tmp_1159;
   std::complex<double> tmp_1162;
   std::complex<double> tmp_1164;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1164 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1162 += tmp_1164;
   std::complex<double> tmp_1163;
   std::complex<double> tmp_1165;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1165 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_1163 += tmp_1165;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1162 * tmp_1163;
   std::complex<double> tmp_1166;
   std::complex<double> tmp_1168;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1169;
      std::complex<double> tmp_1170;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1170 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1169 += tmp_1170;
      tmp_1168 += (Conj(ZE(gI2,j2))) * tmp_1169;
   }
   tmp_1166 += tmp_1168;
   std::complex<double> tmp_1167;
   std::complex<double> tmp_1171;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1172;
      std::complex<double> tmp_1173;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1173 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1172 += tmp_1173;
      tmp_1171 += (ZE(gI1,j4)) * tmp_1172;
   }
   tmp_1167 += tmp_1171;
   result += (-1) * tmp_1166 * tmp_1167;
   if (gO1 < 3) {
      std::complex<double> tmp_1174;
      std::complex<double> tmp_1175;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1175 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1174 += tmp_1175;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1174;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1176;
      std::complex<double> tmp_1177;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1177 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1176 += tmp_1177;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1176;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1178;
      std::complex<double> tmp_1179;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1179 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1178 += tmp_1179;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_1178;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1180;
      std::complex<double> tmp_1181;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1181 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1180 += tmp_1181;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1180;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1182;
      std::complex<double> tmp_1183;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1183 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1182 += tmp_1183;
      result += (-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1182;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1184;
      std::complex<double> tmp_1185;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1185 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1184 += tmp_1185;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1184;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1186;
      std::complex<double> tmp_1187;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1187 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1186 += tmp_1187;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1186;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1188;
      std::complex<double> tmp_1189;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1189 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1188 += tmp_1189;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_1188;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1190;
      std::complex<double> tmp_1191;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1191 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1190 += tmp_1191;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1190;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1192;
      std::complex<double> tmp_1193;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1193 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1192 += tmp_1193;
      result += (-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1192;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1194;
      std::complex<double> tmp_1196;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1196 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1194 += tmp_1196;
      std::complex<double> tmp_1195;
      std::complex<double> tmp_1197;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1198;
         std::complex<double> tmp_1199;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1199 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_1198 += tmp_1199;
         tmp_1197 += (ZE(gI1,j4)) * tmp_1198;
      }
      tmp_1195 += tmp_1197;
      result += (-1) * tmp_1194 * tmp_1195;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1200;
      std::complex<double> tmp_1201;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1201 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_1200 += tmp_1201;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_1200;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1202;
      std::complex<double> tmp_1203;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1203 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_1202 += tmp_1203;
      result += (-0.5*Qe*Ql*Conj(ZE(gI2,gO2))*Sqr(gp)) * tmp_1202;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1204;
      std::complex<double> tmp_1205;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1205 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_1204 += tmp_1205;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_1204;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1206;
      std::complex<double> tmp_1207;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1207 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_1206 += tmp_1207;
      result += (-0.5*Qe*Ql*Conj(ZE(gI2,gO2))*Sqr(gp)) * tmp_1206;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1208;
      std::complex<double> tmp_1210;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1211;
         std::complex<double> tmp_1212;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1212 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_1211 += tmp_1212;
         tmp_1210 += (Conj(ZE(gI2,j2))) * tmp_1211;
      }
      tmp_1208 += tmp_1210;
      std::complex<double> tmp_1209;
      std::complex<double> tmp_1213;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1213 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1209 += tmp_1213;
      result += (-1) * tmp_1208 * tmp_1209;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1214;
      std::complex<double> tmp_1216;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1216 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_1214 += tmp_1216;
      std::complex<double> tmp_1215;
      std::complex<double> tmp_1217;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1217 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_1215 += tmp_1217;
      result += (-1) * tmp_1214 * tmp_1215;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1218;
      std::complex<double> tmp_1219;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1219 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1218 += tmp_1219;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_1218;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1220;
      std::complex<double> tmp_1221;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1221 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1220 += tmp_1221;
      result += (-0.5*Qe*Ql*Sqr(gp)*ZE(gI1,gO1)) * tmp_1220;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1222;
      std::complex<double> tmp_1223;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1223 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1222 += tmp_1223;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_1222;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1224;
      std::complex<double> tmp_1225;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1225 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1224 += tmp_1225;
      result += (-0.5*Qe*Ql*Sqr(gp)*ZE(gI1,gO1)) * tmp_1224;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.15*Conj(ZE(gI2,gO2))*Sqr(g1)*ZE(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -(Conj(ZE(gI2,gO2))*Sqr(gp)*Sqr(Ql)*ZE(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1226;
   std::complex<double> tmp_1228;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1228 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1226 += tmp_1228;
   std::complex<double> tmp_1227;
   std::complex<double> tmp_1229;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1229 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1227 += tmp_1229;
   result += (-0.05*Sqr(g1)) * tmp_1226 * tmp_1227;
   std::complex<double> tmp_1230;
   std::complex<double> tmp_1232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1232 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1230 += tmp_1232;
   std::complex<double> tmp_1231;
   std::complex<double> tmp_1233;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1233 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1231 += tmp_1233;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1230 * tmp_1231;
   std::complex<double> tmp_1234;
   std::complex<double> tmp_1236;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1236 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1234 += tmp_1236;
   std::complex<double> tmp_1235;
   std::complex<double> tmp_1237;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1237 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1235 += tmp_1237;
   result += (0.2*Sqr(g1)) * tmp_1234 * tmp_1235;
   std::complex<double> tmp_1238;
   std::complex<double> tmp_1240;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1240 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1238 += tmp_1240;
   std::complex<double> tmp_1239;
   std::complex<double> tmp_1241;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1241 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1239 += tmp_1241;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_1238 * tmp_1239;
   std::complex<double> tmp_1242;
   std::complex<double> tmp_1244;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1244 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1242 += tmp_1244;
   std::complex<double> tmp_1243;
   std::complex<double> tmp_1245;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1245 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_1243 += tmp_1245;
   result += (-0.05*Sqr(g1)) * tmp_1242 * tmp_1243;
   std::complex<double> tmp_1246;
   std::complex<double> tmp_1248;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1248 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1246 += tmp_1248;
   std::complex<double> tmp_1247;
   std::complex<double> tmp_1249;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1249 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_1247 += tmp_1249;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1246 * tmp_1247;
   std::complex<double> tmp_1250;
   std::complex<double> tmp_1252;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1252 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1250 += tmp_1252;
   std::complex<double> tmp_1251;
   std::complex<double> tmp_1253;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1253 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_1251 += tmp_1253;
   result += (0.2*Sqr(g1)) * tmp_1250 * tmp_1251;
   std::complex<double> tmp_1254;
   std::complex<double> tmp_1256;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1256 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1254 += tmp_1256;
   std::complex<double> tmp_1255;
   std::complex<double> tmp_1257;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1257 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_1255 += tmp_1257;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_1254 * tmp_1255;
   if (gO1 < 3) {
      std::complex<double> tmp_1258;
      std::complex<double> tmp_1259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1259 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1258 += tmp_1259;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1258;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1260;
      std::complex<double> tmp_1261;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1261 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1260 += tmp_1261;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1260;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1262;
      std::complex<double> tmp_1263;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1263 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1262 += tmp_1263;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1262;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1264;
      std::complex<double> tmp_1265;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1265 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1264 += tmp_1265;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1264;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1266;
      std::complex<double> tmp_1267;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1267 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1266 += tmp_1267;
      result += (-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1266;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1268;
      std::complex<double> tmp_1269;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1269 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1268 += tmp_1269;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1268;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1270;
      std::complex<double> tmp_1271;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1271 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1270 += tmp_1271;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1270;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1272;
      std::complex<double> tmp_1273;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1273 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1272 += tmp_1273;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1272;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1274;
      std::complex<double> tmp_1275;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1275 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1274 += tmp_1275;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1274;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1276;
      std::complex<double> tmp_1277;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1277 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1276 += tmp_1277;
      result += (-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1276;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1278;
   std::complex<double> tmp_1279;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1280;
      std::complex<double> tmp_1281;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1281 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1280 += tmp_1281;
      tmp_1279 += (Conj(ZE(gI1,j2))) * tmp_1280;
   }
   tmp_1278 += tmp_1279;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,1))) *
      tmp_1278;
   std::complex<double> tmp_1282;
   std::complex<double> tmp_1283;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1284;
      std::complex<double> tmp_1285;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1285 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1284 += tmp_1285;
      tmp_1283 += (Conj(ZE(gI1,j2))) * tmp_1284;
   }
   tmp_1282 += tmp_1283;
   result += (std::complex<double>(0,-0.5)*vu*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_1282;
   std::complex<double> tmp_1286;
   std::complex<double> tmp_1287;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1288;
      std::complex<double> tmp_1289;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1289 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_1288 += tmp_1289;
      tmp_1287 += (Conj(ZE(gI1,j2))) * tmp_1288;
   }
   tmp_1286 += tmp_1287;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_1286;
   if (gO2 < 3) {
      std::complex<double> tmp_1290;
      std::complex<double> tmp_1291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1291 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1290 += tmp_1291;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,1))*Lambdax) *
         tmp_1290;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1292;
      std::complex<double> tmp_1293;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1293 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1292 += tmp_1293;
      result += (std::complex<double>(0,0.5)*vu*Conj(ZA(gI2,2))*Lambdax) *
         tmp_1292;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1294;
      std::complex<double> tmp_1295;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1295 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_1294 += tmp_1295;
      result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,0)))
         * tmp_1294;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   std::complex<double> tmp_1296;
   std::complex<double> tmp_1297;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1297 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1296 += tmp_1297;
   result += (0.3*vd*Sqr(g1)*ZH(gI2,0)) * tmp_1296;
   std::complex<double> tmp_1298;
   std::complex<double> tmp_1299;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1299 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1298 += tmp_1299;
   result += (-(Qe*QHd*vd*Sqr(gp)*ZH(gI2,0))) * tmp_1298;
   std::complex<double> tmp_1300;
   std::complex<double> tmp_1301;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1302;
      std::complex<double> tmp_1303;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1303 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_1302 += tmp_1303;
      tmp_1301 += (Conj(ZE(gI1,j2))) * tmp_1302;
   }
   tmp_1300 += tmp_1301;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_1300;
   std::complex<double> tmp_1304;
   std::complex<double> tmp_1305;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1306;
      std::complex<double> tmp_1307;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1308;
         std::complex<double> tmp_1309;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1309 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1308 += tmp_1309;
         tmp_1307 += (KroneckerDelta(gO2,3 + j2)) * tmp_1308;
      }
      tmp_1306 += tmp_1307;
      tmp_1305 += (Conj(ZE(gI1,3 + j3))) * tmp_1306;
   }
   tmp_1304 += tmp_1305;
   result += (-(vd*ZH(gI2,0))) * tmp_1304;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -(QHd*Ql*vd*Conj(ZE(gI1,gO2))*Sqr(gp)*ZH(gI2,0));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1310;
      std::complex<double> tmp_1311;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1311 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_1310 += tmp_1311;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_1310;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1312;
      std::complex<double> tmp_1313;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1314;
         std::complex<double> tmp_1315;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1315 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_1314 += tmp_1315;
         tmp_1313 += (Conj(ZE(gI1,j2))) * tmp_1314;
      }
      tmp_1312 += tmp_1313;
      result += (-(vd*ZH(gI2,0))) * tmp_1312;
   }
   std::complex<double> tmp_1316;
   std::complex<double> tmp_1317;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1317 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1316 += tmp_1317;
   result += (-0.3*vu*Sqr(g1)*ZH(gI2,1)) * tmp_1316;
   std::complex<double> tmp_1318;
   std::complex<double> tmp_1319;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1319 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1318 += tmp_1319;
   result += (-(Qe*QHu*vu*Sqr(gp)*ZH(gI2,1))) * tmp_1318;
   std::complex<double> tmp_1320;
   std::complex<double> tmp_1321;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1322;
      std::complex<double> tmp_1323;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1323 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1322 += tmp_1323;
      tmp_1321 += (Conj(ZE(gI1,j2))) * tmp_1322;
   }
   tmp_1320 += tmp_1321;
   result += (0.5*vS*Conj(Lambdax)*ZH(gI2,1)) * tmp_1320;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -(QHu*Ql*vu*Conj(ZE(gI1,gO2))*Sqr(gp)*ZH(gI2,1));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1324;
      std::complex<double> tmp_1325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1325 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1324 += tmp_1325;
      result += (0.5*vS*Lambdax*ZH(gI2,1)) * tmp_1324;
   }
   std::complex<double> tmp_1326;
   std::complex<double> tmp_1327;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1327 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1326 += tmp_1327;
   result += (-(Qe*Qs*vS*Sqr(gp)*ZH(gI2,2))) * tmp_1326;
   std::complex<double> tmp_1328;
   std::complex<double> tmp_1329;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1330;
      std::complex<double> tmp_1331;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1331 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1330 += tmp_1331;
      tmp_1329 += (Conj(ZE(gI1,j2))) * tmp_1330;
   }
   tmp_1328 += tmp_1329;
   result += (0.5*vu*Conj(Lambdax)*ZH(gI2,2)) * tmp_1328;
   if (gO2 < 3) {
      result += -(Ql*Qs*vS*Conj(ZE(gI1,gO2))*Sqr(gp)*ZH(gI2,2));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1332;
      std::complex<double> tmp_1333;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1333 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1332 += tmp_1333;
      result += (0.5*vu*Lambdax*ZH(gI2,2)) * tmp_1332;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVWmSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZV(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVPSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1334;
   std::complex<double> tmp_1335;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1335 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1334 += tmp_1335;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_1334;
   if (gO2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZE(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVZSe(unsigned gO2, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1336;
   std::complex<double> tmp_1337;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1337 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1336 += tmp_1337;
   result += (0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_1336;
   std::complex<double> tmp_1338;
   std::complex<double> tmp_1339;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1339 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1338 += tmp_1339;
   result += (-(gp*Qe*Sin(ThetaWp()))) * tmp_1338;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Cos(ThetaWp())*Sin(
         ThetaW());
   }
   if (gO2 < 3) {
      result += gp*Ql*Conj(ZE(gI2,gO2))*Sin(ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVZpSe(unsigned gO2, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1340;
   std::complex<double> tmp_1341;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1341 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1340 += tmp_1341;
   result += (-(gp*Qe*Cos(ThetaWp()))) * tmp_1340;
   std::complex<double> tmp_1342;
   std::complex<double> tmp_1343;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1343 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1342 += tmp_1343;
   result += (-0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_1342;
   if (gO2 < 3) {
      result += gp*Ql*Conj(ZE(gI2,gO2))*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gO2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*vS*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp())
      ) + vd*KroneckerDelta(0,gO2)*(20*g2*gp*QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(
      ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin
      (ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) + vu*KroneckerDelta(1
      ,gO2)*(-20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) -
      15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) +
      g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*
      Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*
      Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZpVZ(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*vS*Cos(ThetaWp())*KroneckerDelta(2,gO2)*Sin(ThetaWp())*Sqr(
      gp)*Sqr(Qs) - vd*KroneckerDelta(0,gO2)*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(
      g2)*Sqr(Cos(ThetaW())) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Cos(
      ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1
      )*Sqr(Sin(ThetaW()))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Sin(
      ThetaWp())) + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW())*Sin(ThetaWp()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(
      ThetaWp())))) - vu*KroneckerDelta(1,gO2)*(5*Cos(ThetaWp())*Sin(ThetaWp())*
      Sqr(g2)*Sqr(Cos(ThetaW())) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(
      Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*
      Sqr(g1)*Sqr(Sin(ThetaW()))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(
      Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(
      Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZpVZp(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*vS*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp())
      ) + vd*KroneckerDelta(0,gO2)*(-2*gp*QHd*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))) +
      vu*KroneckerDelta(1,gO2)*(2*gp*QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*
      g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp()))
      + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()
      )) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(vd*KroneckerDelta(0,gO2) + vu*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1))*Sqr(g2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1))*Sqr(g2)
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargZgZ(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*vS*KroneckerDelta(2,gO1)*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp(
      ))) - vd*KroneckerDelta(0,gO1)*(20*g2*gp*QHd*Cos(ThetaW())*Cos(ThetaWp())*
      Sin(ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1
      *Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) - vu*KroneckerDelta(1
      ,gO1)*(-20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) -
      15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) +
      g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*
      Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*
      Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargZpgZ(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*vS*Cos(ThetaWp())*KroneckerDelta(2,gO1)*Sin(ThetaWp())*
      Sqr(gp)*Sqr(Qs) + vd*KroneckerDelta(0,gO1)*(5*Cos(ThetaWp())*Sin(ThetaWp())*
      Sqr(g2)*Sqr(Cos(ThetaW())) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(
      Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHd) + 3*
      Sqr(g1)*Sqr(Sin(ThetaW()))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(
      Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(
      Sin(ThetaWp())))) + vu*KroneckerDelta(1,gO1)*(5*Cos(ThetaWp())*Sin(ThetaWp()
      )*Sqr(g2)*Sqr(Cos(ThetaW())) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr
      (Cos(ThetaWp())) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*
      Sqr(g1)*Sqr(Sin(ThetaW()))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(
      Sin(ThetaWp())) + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(
      Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargZpgZp(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*vS*KroneckerDelta(2,gO1)*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp(
      ))) - vd*KroneckerDelta(0,gO1)*(-2*gp*QHd*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))) -
      vu*KroneckerDelta(1,gO1)*(2*gp*QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*
      g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp()))
      + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()
      )) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)
      *Sqr(Sin(ThetaWp())) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(20*g2*gp
      *QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*
      QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-20*g2*
      gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*
      gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)
      *Sqr(Cos(ThetaWp())) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-2*gp*
      QHd*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp(
      )) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(2*gp*QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin
      (2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*
      (7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(
      Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*((AbsSqr(
      Lambdax) + QHd*Qs*Sqr(gp))*ZP(gI1,0)*ZP(gI2,0) + (AbsSqr(Lambdax) + QHu*Qs*
      Sqr(gp))*ZP(gI1,1)*ZP(gI2,1)) - KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,
      gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,
      1)) + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*
      ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI1,1
      )*ZP(gI2,1))) + KroneckerDelta(1,gO1)*(-5*KroneckerDelta(0,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp)))*ZP(gI1,
      0)*ZP(gI2,0) - (3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHu)))*ZP(gI1,1)*ZP(
      gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO2)*(ZP(gI1,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2)
      + 20*Sqr(gp)*Sqr(QHd))*ZP(gI2,0) + 5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(
      gI2,1)) + ZP(gI1,1)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vd*(-3*
      Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI2,1)))) - 10*KroneckerDelta(2
      ,gO2)*(ZP(gI1,0)*(2*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZP(gI2,0) +
      1.4142135623730951*Conj(TLambdax)*ZP(gI2,1)) + ZP(gI1,1)*(1.4142135623730951
      *TLambdax*ZP(gI2,0) + 2*vS*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZP(gI2,1))) +
      KroneckerDelta(1,gO2)*(ZP(gI1,0)*(vu*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr
      (gp)))*ZP(gI2,0) - 5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,1)) - ZP(gI1,1
      )*(5*vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vu*(3*Sqr(g1) + 5*Sqr(g2)
      + 20*Sqr(gp)*Sqr(QHu))*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*KroneckerDelta(0,gO2)*UM(gI1,1)*UP(gI2,0) +
      (g2*KroneckerDelta(1,gO2)*UM(gI1,0) + Conj(Lambdax)*KroneckerDelta(2,gO2)*
      UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*Conj(UP(gI1,1))*
      KroneckerDelta(1,gO1) + Conj(UM(gI2,1))*(g2*Conj(UP(gI1,0))*KroneckerDelta(0
      ,gO1) + Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*
      QHu*Sqr(gp)) + 20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(
      Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*
      Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))) + Conj(ZA(gI1,1))*Conj(ZA(gI2,
      1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*
      Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))) -
      20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,
      gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + KroneckerDelta(2,
      gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = -0.05*KroneckerDelta(gI1,gI2)*(20*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*
      (-3*Sqr(g1) - 5*Sqr(g2) + 20*QHu*Ql*Sqr(gp)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*QHd*Ql*Sqr(gp))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(-(KroneckerDelta(1,gO2)*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd*QHu*Sqr(gp)))*(ZH(gI1,1)*ZH(gI2,0)
      + ZH(gI1,0)*ZH(gI2,1))) + 20*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*
      Qs*Sqr(gp))*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + KroneckerDelta(0,
      gO2)*(3*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))*ZH(gI1,0)*ZH(gI2,0) +
      (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZH(gI1,1)
      *ZH(gI2,1) + 20*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,2)*ZH(gI2,2)))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1
      ) + 5*(Sqr(g2) - 4*QHd*QHu*Sqr(gp)))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2
      ,1)) - 20*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*(ZH(gI1,2
      )*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*((-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZH(gI1,0)*ZH(gI2,0) -
      3*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZH(gI1,1)*ZH(gI2,1) - 20*(
      AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZH(gI1,2)*ZH(gI2,2))) - 20*KroneckerDelta(
      2,gO1)*(KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*(ZH(gI1,2)*
      ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) +
      QHu*Qs*Sqr(gp))*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) + KroneckerDelta
      (2,gO2)*((AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,0)*ZH(gI2,0) + (AbsSqr(
      Lambdax) + QHu*Qs*Sqr(gp))*ZH(gI1,1)*ZH(gI2,1) + 3*Sqr(gp)*Sqr(Qs)*ZH(gI1,2)
      *ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-5*Conj(ZA(gI1,2))*(1.4142135623730951*Conj(TLambdax)*(Conj(
      ZA(gI2,1))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*KroneckerDelta(1,gO2)) +
      4*Conj(ZA(gI2,2))*(vd*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp
      )) + vu*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + vS*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)) + 1.4142135623730951*(Conj(ZA(gI2,1))
      *KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*KroneckerDelta(1,gO2))*TLambdax) +
      Conj(ZA(gI1,1))*(Conj(ZA(gI2,1))*(vd*KroneckerDelta(0,gO2)*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*vS*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - vu*KroneckerDelta
      (1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))) - 7.0710678118654755*
      (Conj(ZA(gI2,2))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*KroneckerDelta(2,
      gO2))*(Conj(TLambdax) + TLambdax)) - Conj(ZA(gI1,0))*(Conj(ZA(gI2,0))*(vu*
      KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*
      QHu*Sqr(gp)) + 20*vS*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)
      ) + vd*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))
      + 7.0710678118654755*(Conj(ZA(gI2,2))*KroneckerDelta(1,gO2) + Conj(ZA(gI2,1
      ))*KroneckerDelta(2,gO2))*(Conj(TLambdax) + TLambdax)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = -0.05*KroneckerDelta(gI1,gI2)*(20*Ql*Qs*vS*KroneckerDelta(2,gO2)*
      Sqr(gp) + vu*KroneckerDelta(1,gO2)*(-3*Sqr(g1) - 5*Sqr(g2) + 20*QHu*Ql*Sqr(
      gp)) + vd*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*QHd*Ql*Sqr(gp)))
      );

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.35355339059327373)*(Conj(TLambdax) -
      TLambdax)*(Conj(ZA(gI2,2))*(KroneckerDelta(1,gO2)*ZH(gI1,0) + KroneckerDelta
      (0,gO2)*ZH(gI1,1)) + Conj(ZA(gI2,1))*(KroneckerDelta(2,gO2)*ZH(gI1,0) +
      KroneckerDelta(0,gO2)*ZH(gI1,2)) + Conj(ZA(gI2,0))*(KroneckerDelta(2,gO2)*ZH
      (gI1,1) + KroneckerDelta(1,gO2)*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-5*KroneckerDelta(2,gO2)*(-1.4142135623730951*Conj(TLambdax)*
      ZH(gI1,1)*ZH(gI2,0) - 1.4142135623730951*TLambdax*ZH(gI1,1)*ZH(gI2,0) + 4*vd
      *AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,0) + 4*QHd*Qs*vd*Sqr(gp)*ZH(gI1,2)*ZH(gI2,
      0) + 4*vS*AbsSqr(Lambdax)*ZH(gI1,1)*ZH(gI2,1) + 4*QHu*Qs*vS*Sqr(gp)*ZH(gI1,1
      )*ZH(gI2,1) + 4*vu*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,1) + 4*QHu*Qs*vu*Sqr(gp)
      *ZH(gI1,2)*ZH(gI2,1) + 4*vu*AbsSqr(Lambdax)*ZH(gI1,1)*ZH(gI2,2) + 4*QHu*Qs*
      vu*Sqr(gp)*ZH(gI1,1)*ZH(gI2,2) + 12*vS*Sqr(gp)*Sqr(Qs)*ZH(gI1,2)*ZH(gI2,2) +
      ZH(gI1,0)*(4*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI2,0) -
      1.4142135623730951*Conj(TLambdax)*ZH(gI2,1) - 1.4142135623730951*TLambdax*ZH
      (gI2,1) + 4*vd*AbsSqr(Lambdax)*ZH(gI2,2) + 4*QHd*Qs*vd*Sqr(gp)*ZH(gI2,2))) -
      KroneckerDelta(0,gO2)*(5*ZH(gI1,2)*(4*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))
      *ZH(gI2,0) - 1.4142135623730951*Conj(TLambdax)*ZH(gI2,1) -
      1.4142135623730951*TLambdax*ZH(gI2,1) + 4*vd*AbsSqr(Lambdax)*ZH(gI2,2) + 4*
      QHd*Qs*vd*Sqr(gp)*ZH(gI2,2)) + ZH(gI1,0)*(3*vd*(3*Sqr(g1) + 5*(Sqr(g2) + 4*
      Sqr(gp)*Sqr(QHd)))*ZH(gI2,0) + vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2
      ) + 20*QHd*QHu*Sqr(gp))*ZH(gI2,1) + 20*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))
      *ZH(gI2,2)) - ZH(gI1,1)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) -
      20*QHd*QHu*Sqr(gp))*ZH(gI2,0) + vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(
      g2) - 20*QHd*QHu*Sqr(gp))*ZH(gI2,1) + 7.0710678118654755*(Conj(TLambdax) +
      TLambdax)*ZH(gI2,2))) + KroneckerDelta(1,gO2)*(ZH(gI1,1)*(vd*(-20*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZH(gI2,0) - 3*vu*(3*
      Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZH(gI2,1) - 20*vS*(AbsSqr(Lambdax
      ) + QHu*Qs*Sqr(gp))*ZH(gI2,2)) + ZH(gI1,0)*(vu*(-20*AbsSqr(Lambdax) + 3*Sqr(
      g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZH(gI2,0) + vd*(-20*AbsSqr(Lambdax) +
      3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZH(gI2,1) + 7.0710678118654755*(
      Conj(TLambdax) + TLambdax)*ZH(gI2,2)) + 5*ZH(gI1,2)*(1.4142135623730951*Conj
      (TLambdax)*ZH(gI2,0) + 1.4142135623730951*TLambdax*ZH(gI2,0) - 4*(AbsSqr(
      Lambdax) + QHu*Qs*Sqr(gp))*(vS*ZH(gI2,1) + vu*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1344;
   std::complex<double> tmp_1345;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1346;
      std::complex<double> tmp_1347;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1347 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1346 += tmp_1347;
      tmp_1345 += (ZDL(gI1,j2)) * tmp_1346;
   }
   tmp_1344 += tmp_1345;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1344;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1348;
   std::complex<double> tmp_1349;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1350;
      std::complex<double> tmp_1351;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1351 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1350 += tmp_1351;
      tmp_1349 += (Conj(ZDL(gI2,j2))) * tmp_1350;
   }
   tmp_1348 += tmp_1349;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1348;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1352;
   std::complex<double> tmp_1353;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1354;
      std::complex<double> tmp_1355;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1355 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1354 += tmp_1355;
      tmp_1353 += (ZEL(gI1,j2)) * tmp_1354;
   }
   tmp_1352 += tmp_1353;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1352;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1356;
   std::complex<double> tmp_1357;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1358;
      std::complex<double> tmp_1359;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1359 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1358 += tmp_1359;
      tmp_1357 += (Conj(ZEL(gI2,j2))) * tmp_1358;
   }
   tmp_1356 += tmp_1357;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1356;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1360;
   std::complex<double> tmp_1361;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1362;
      std::complex<double> tmp_1363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1363 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1362 += tmp_1363;
      tmp_1361 += (ZUL(gI1,j2)) * tmp_1362;
   }
   tmp_1360 += tmp_1361;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1360;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1364;
   std::complex<double> tmp_1365;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1366;
      std::complex<double> tmp_1367;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1367 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1366 += tmp_1367;
      tmp_1365 += (Conj(ZUL(gI2,j2))) * tmp_1366;
   }
   tmp_1364 += tmp_1365;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_1364;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_1368;
   std::complex<double> tmp_1369;
   std::complex<double> tmp_1370;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1370 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1369 += tmp_1370;
   tmp_1368 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1369;
   std::complex<double> tmp_1371;
   std::complex<double> tmp_1372;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1372 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1371 += tmp_1372;
   tmp_1368 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1371;
   std::complex<double> tmp_1373;
   std::complex<double> tmp_1374;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1374 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1373 += tmp_1374;
   tmp_1368 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1373;
   std::complex<double> tmp_1375;
   std::complex<double> tmp_1376;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1376 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1375 += tmp_1376;
   tmp_1368 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1375;
   std::complex<double> tmp_1377;
   std::complex<double> tmp_1378;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1378 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1377 += tmp_1378;
   tmp_1368 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1377;
   std::complex<double> tmp_1379;
   std::complex<double> tmp_1380;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1380 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1379 += tmp_1380;
   tmp_1368 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1379;
   std::complex<double> tmp_1381;
   std::complex<double> tmp_1382;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1382 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1381 += tmp_1382;
   tmp_1368 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1381;
   std::complex<double> tmp_1383;
   std::complex<double> tmp_1384;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1384 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1383 += tmp_1384;
   tmp_1368 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1383;
   std::complex<double> tmp_1385;
   std::complex<double> tmp_1386;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1386 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1385 += tmp_1386;
   tmp_1368 += (std::complex<double>(0,-1)*Qd*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1385;
   std::complex<double> tmp_1387;
   std::complex<double> tmp_1388;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1388 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1387 += tmp_1388;
   tmp_1368 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1387;
   std::complex<double> tmp_1389;
   std::complex<double> tmp_1390;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1390 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1389 += tmp_1390;
   tmp_1368 += (std::complex<double>(0,-1)*Qd*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1389;
   std::complex<double> tmp_1391;
   std::complex<double> tmp_1392;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1392 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1391 += tmp_1392;
   tmp_1368 += (std::complex<double>(0,-1)*Qd*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1391;
   std::complex<double> tmp_1393;
   std::complex<double> tmp_1394;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1395;
      std::complex<double> tmp_1396;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1396 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1395 += tmp_1396;
      tmp_1394 += (Conj(ZD(gI2,j2))) * tmp_1395;
   }
   tmp_1393 += tmp_1394;
   tmp_1368 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1)) * tmp_1393;
   std::complex<double> tmp_1397;
   std::complex<double> tmp_1398;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1399;
      std::complex<double> tmp_1400;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1400 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1399 += tmp_1400;
      tmp_1398 += (Conj(ZD(gI2,j2))) * tmp_1399;
   }
   tmp_1397 += tmp_1398;
   tmp_1368 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1)
      *KroneckerDelta(2,gO2)) * tmp_1397;
   std::complex<double> tmp_1401;
   std::complex<double> tmp_1402;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1403;
      std::complex<double> tmp_1404;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1404 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1403 += tmp_1404;
      tmp_1402 += (ZD(gI1,j2)) * tmp_1403;
   }
   tmp_1401 += tmp_1402;
   tmp_1368 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1401;
   std::complex<double> tmp_1405;
   std::complex<double> tmp_1406;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1407;
      std::complex<double> tmp_1408;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1408 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1407 += tmp_1408;
      tmp_1406 += (ZD(gI1,j2)) * tmp_1407;
   }
   tmp_1405 += tmp_1406;
   tmp_1368 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1405;
   std::complex<double> tmp_1409;
   std::complex<double> tmp_1410;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1411;
      std::complex<double> tmp_1412;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1413;
         std::complex<double> tmp_1414;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1414 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1413 += tmp_1414;
         tmp_1412 += (ZD(gI1,3 + j2)) * tmp_1413;
      }
      tmp_1411 += tmp_1412;
      tmp_1410 += (Conj(ZD(gI2,3 + j3))) * tmp_1411;
   }
   tmp_1409 += tmp_1410;
   tmp_1368 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1409;
   std::complex<double> tmp_1415;
   std::complex<double> tmp_1416;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1417;
      std::complex<double> tmp_1418;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1419;
         std::complex<double> tmp_1420;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1420 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1419 += tmp_1420;
         tmp_1418 += (Conj(ZD(gI2,j2))) * tmp_1419;
      }
      tmp_1417 += tmp_1418;
      tmp_1416 += (ZD(gI1,j3)) * tmp_1417;
   }
   tmp_1415 += tmp_1416;
   tmp_1368 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1415;
   result += (std::complex<double>(0,-1)) * tmp_1368;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_1421;
   std::complex<double> tmp_1422;
   std::complex<double> tmp_1423;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1423 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1422 += tmp_1423;
   tmp_1421 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1422;
   std::complex<double> tmp_1424;
   std::complex<double> tmp_1425;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1425 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1424 += tmp_1425;
   tmp_1421 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1424;
   std::complex<double> tmp_1426;
   std::complex<double> tmp_1427;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1427 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1426 += tmp_1427;
   tmp_1421 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1426;
   std::complex<double> tmp_1428;
   std::complex<double> tmp_1429;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1429 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1428 += tmp_1429;
   tmp_1421 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1428;
   std::complex<double> tmp_1430;
   std::complex<double> tmp_1431;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1431 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1430 += tmp_1431;
   tmp_1421 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1430;
   std::complex<double> tmp_1432;
   std::complex<double> tmp_1433;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1433 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1432 += tmp_1433;
   tmp_1421 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1432;
   std::complex<double> tmp_1434;
   std::complex<double> tmp_1435;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1435 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1434 += tmp_1435;
   tmp_1421 += (std::complex<double>(0,-1)*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1434;
   std::complex<double> tmp_1436;
   std::complex<double> tmp_1437;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1437 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1436 += tmp_1437;
   tmp_1421 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1436;
   std::complex<double> tmp_1438;
   std::complex<double> tmp_1439;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1439 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1438 += tmp_1439;
   tmp_1421 += (std::complex<double>(0,-1)*Qe*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1438;
   std::complex<double> tmp_1440;
   std::complex<double> tmp_1441;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1441 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1440 += tmp_1441;
   tmp_1421 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1440;
   std::complex<double> tmp_1442;
   std::complex<double> tmp_1443;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1443 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1442 += tmp_1443;
   tmp_1421 += (std::complex<double>(0,-1)*Qe*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1442;
   std::complex<double> tmp_1444;
   std::complex<double> tmp_1445;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1445 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1444 += tmp_1445;
   tmp_1421 += (std::complex<double>(0,-1)*Qe*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1444;
   std::complex<double> tmp_1446;
   std::complex<double> tmp_1447;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1448;
      std::complex<double> tmp_1449;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1449 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1448 += tmp_1449;
      tmp_1447 += (Conj(ZE(gI2,j2))) * tmp_1448;
   }
   tmp_1446 += tmp_1447;
   tmp_1421 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1)) * tmp_1446;
   std::complex<double> tmp_1450;
   std::complex<double> tmp_1451;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1452;
      std::complex<double> tmp_1453;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1453 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1452 += tmp_1453;
      tmp_1451 += (Conj(ZE(gI2,j2))) * tmp_1452;
   }
   tmp_1450 += tmp_1451;
   tmp_1421 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1)
      *KroneckerDelta(2,gO2)) * tmp_1450;
   std::complex<double> tmp_1454;
   std::complex<double> tmp_1455;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1456;
      std::complex<double> tmp_1457;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1457 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1456 += tmp_1457;
      tmp_1455 += (ZE(gI1,j2)) * tmp_1456;
   }
   tmp_1454 += tmp_1455;
   tmp_1421 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1454;
   std::complex<double> tmp_1458;
   std::complex<double> tmp_1459;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1460;
      std::complex<double> tmp_1461;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1461 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1460 += tmp_1461;
      tmp_1459 += (ZE(gI1,j2)) * tmp_1460;
   }
   tmp_1458 += tmp_1459;
   tmp_1421 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1458;
   std::complex<double> tmp_1462;
   std::complex<double> tmp_1463;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1464;
      std::complex<double> tmp_1465;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1466;
         std::complex<double> tmp_1467;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1467 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1466 += tmp_1467;
         tmp_1465 += (ZE(gI1,3 + j2)) * tmp_1466;
      }
      tmp_1464 += tmp_1465;
      tmp_1463 += (Conj(ZE(gI2,3 + j3))) * tmp_1464;
   }
   tmp_1462 += tmp_1463;
   tmp_1421 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1462;
   std::complex<double> tmp_1468;
   std::complex<double> tmp_1469;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1470;
      std::complex<double> tmp_1471;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1472;
         std::complex<double> tmp_1473;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1473 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1472 += tmp_1473;
         tmp_1471 += (Conj(ZE(gI2,j2))) * tmp_1472;
      }
      tmp_1470 += tmp_1471;
      tmp_1469 += (ZE(gI1,j3)) * tmp_1470;
   }
   tmp_1468 += tmp_1469;
   tmp_1421 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1468;
   result += (std::complex<double>(0,-1)) * tmp_1421;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_1474;
   std::complex<double> tmp_1475;
   std::complex<double> tmp_1476;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1476 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1475 += tmp_1476;
   tmp_1474 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1475;
   std::complex<double> tmp_1477;
   std::complex<double> tmp_1478;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1478 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1477 += tmp_1478;
   tmp_1474 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1477;
   std::complex<double> tmp_1479;
   std::complex<double> tmp_1480;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1480 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1479 += tmp_1480;
   tmp_1474 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1479;
   std::complex<double> tmp_1481;
   std::complex<double> tmp_1482;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1482 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1481 += tmp_1482;
   tmp_1474 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1481;
   std::complex<double> tmp_1483;
   std::complex<double> tmp_1484;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1484 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1483 += tmp_1484;
   tmp_1474 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1483;
   std::complex<double> tmp_1485;
   std::complex<double> tmp_1486;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1486 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1485 += tmp_1486;
   tmp_1474 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1485;
   std::complex<double> tmp_1487;
   std::complex<double> tmp_1488;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1488 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1487 += tmp_1488;
   tmp_1474 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1487;
   std::complex<double> tmp_1489;
   std::complex<double> tmp_1490;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1490 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1489 += tmp_1490;
   tmp_1474 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1489;
   std::complex<double> tmp_1491;
   std::complex<double> tmp_1492;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1492 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1491 += tmp_1492;
   tmp_1474 += (std::complex<double>(0,-1)*QHd*Qu*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1491;
   std::complex<double> tmp_1493;
   std::complex<double> tmp_1494;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1494 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1493 += tmp_1494;
   tmp_1474 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1493;
   std::complex<double> tmp_1495;
   std::complex<double> tmp_1496;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1496 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1495 += tmp_1496;
   tmp_1474 += (std::complex<double>(0,-1)*QHu*Qu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1495;
   std::complex<double> tmp_1497;
   std::complex<double> tmp_1498;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1498 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1497 += tmp_1498;
   tmp_1474 += (std::complex<double>(0,-1)*Qs*Qu*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1497;
   std::complex<double> tmp_1499;
   std::complex<double> tmp_1500;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1501;
      std::complex<double> tmp_1502;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1502 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1501 += tmp_1502;
      tmp_1500 += (Conj(ZU(gI2,j2))) * tmp_1501;
   }
   tmp_1499 += tmp_1500;
   tmp_1474 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2)
      *KroneckerDelta(2,gO1)) * tmp_1499;
   std::complex<double> tmp_1503;
   std::complex<double> tmp_1504;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1505;
      std::complex<double> tmp_1506;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1506 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1505 += tmp_1506;
      tmp_1504 += (Conj(ZU(gI2,j2))) * tmp_1505;
   }
   tmp_1503 += tmp_1504;
   tmp_1474 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1)
      *KroneckerDelta(2,gO2)) * tmp_1503;
   std::complex<double> tmp_1507;
   std::complex<double> tmp_1508;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1509;
      std::complex<double> tmp_1510;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1510 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1509 += tmp_1510;
      tmp_1508 += (ZU(gI1,j2)) * tmp_1509;
   }
   tmp_1507 += tmp_1508;
   tmp_1474 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1507;
   std::complex<double> tmp_1511;
   std::complex<double> tmp_1512;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1513;
      std::complex<double> tmp_1514;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1514 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1513 += tmp_1514;
      tmp_1512 += (ZU(gI1,j2)) * tmp_1513;
   }
   tmp_1511 += tmp_1512;
   tmp_1474 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1511;
   std::complex<double> tmp_1515;
   std::complex<double> tmp_1516;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1517;
      std::complex<double> tmp_1518;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1519;
         std::complex<double> tmp_1520;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1520 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1519 += tmp_1520;
         tmp_1518 += (ZU(gI1,3 + j2)) * tmp_1519;
      }
      tmp_1517 += tmp_1518;
      tmp_1516 += (Conj(ZU(gI2,3 + j3))) * tmp_1517;
   }
   tmp_1515 += tmp_1516;
   tmp_1474 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1515;
   std::complex<double> tmp_1521;
   std::complex<double> tmp_1522;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1523;
      std::complex<double> tmp_1524;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1525;
         std::complex<double> tmp_1526;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1526 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1525 += tmp_1526;
         tmp_1524 += (Conj(ZU(gI2,j2))) * tmp_1525;
      }
      tmp_1523 += tmp_1524;
      tmp_1522 += (ZU(gI1,j3)) * tmp_1523;
   }
   tmp_1521 += tmp_1522;
   tmp_1474 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1521;
   result += (std::complex<double>(0,-1)) * tmp_1474;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_1527;
   std::complex<double> tmp_1528;
   std::complex<double> tmp_1529;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1529 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1528 += tmp_1529;
   tmp_1527 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1528;
   std::complex<double> tmp_1530;
   std::complex<double> tmp_1531;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1531 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1530 += tmp_1531;
   tmp_1527 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1530;
   std::complex<double> tmp_1532;
   std::complex<double> tmp_1533;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1533 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1532 += tmp_1533;
   tmp_1527 += (std::complex<double>(0,-1)*QHd*Qq*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_1532;
   std::complex<double> tmp_1534;
   std::complex<double> tmp_1535;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1535 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1534 += tmp_1535;
   tmp_1527 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1534;
   std::complex<double> tmp_1536;
   std::complex<double> tmp_1537;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1537 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1536 += tmp_1537;
   tmp_1527 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1536;
   std::complex<double> tmp_1538;
   std::complex<double> tmp_1539;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1539 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1538 += tmp_1539;
   tmp_1527 += (std::complex<double>(0,-1)*QHu*Qq*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_1538;
   std::complex<double> tmp_1540;
   std::complex<double> tmp_1541;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1541 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1540 += tmp_1541;
   tmp_1527 += (std::complex<double>(0,-1)*Qq*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_1540;
   std::complex<double> tmp_1542;
   std::complex<double> tmp_1543;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1543 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1542 += tmp_1543;
   tmp_1527 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1542;
   std::complex<double> tmp_1544;
   std::complex<double> tmp_1545;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1545 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1544 += tmp_1545;
   tmp_1527 += (std::complex<double>(0,-1)*Qd*QHd*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_1544;
   std::complex<double> tmp_1546;
   std::complex<double> tmp_1547;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1547 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1546 += tmp_1547;
   tmp_1527 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1546;
   std::complex<double> tmp_1548;
   std::complex<double> tmp_1549;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1549 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1548 += tmp_1549;
   tmp_1527 += (std::complex<double>(0,-1)*Qd*QHu*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_1548;
   std::complex<double> tmp_1550;
   std::complex<double> tmp_1551;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1551 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1550 += tmp_1551;
   tmp_1527 += (std::complex<double>(0,-1)*Qd*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_1550;
   std::complex<double> tmp_1552;
   std::complex<double> tmp_1553;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1554;
      std::complex<double> tmp_1555;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1555 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1554 += tmp_1555;
      tmp_1553 += (Conj(ZD(gI2,j2))) * tmp_1554;
   }
   tmp_1552 += tmp_1553;
   tmp_1527 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(1,
      gO2)) * tmp_1552;
   std::complex<double> tmp_1556;
   std::complex<double> tmp_1557;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1558;
      std::complex<double> tmp_1559;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1559 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1558 += tmp_1559;
      tmp_1557 += (Conj(ZD(gI2,j2))) * tmp_1558;
   }
   tmp_1556 += tmp_1557;
   tmp_1527 += (std::complex<double>(0,0.5)*vu*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_1556;
   std::complex<double> tmp_1560;
   std::complex<double> tmp_1561;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1562;
      std::complex<double> tmp_1563;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1563 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1562 += tmp_1563;
      tmp_1561 += (Conj(ZD(gI2,j2))) * tmp_1562;
   }
   tmp_1560 += tmp_1561;
   tmp_1527 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1560;
   std::complex<double> tmp_1564;
   std::complex<double> tmp_1565;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1566;
      std::complex<double> tmp_1567;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1567 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1566 += tmp_1567;
      tmp_1565 += (ZD(gI1,j2)) * tmp_1566;
   }
   tmp_1564 += tmp_1565;
   tmp_1527 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_1564;
   std::complex<double> tmp_1568;
   std::complex<double> tmp_1569;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1570;
      std::complex<double> tmp_1571;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1571 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1570 += tmp_1571;
      tmp_1569 += (ZD(gI1,j2)) * tmp_1570;
   }
   tmp_1568 += tmp_1569;
   tmp_1527 += (std::complex<double>(0,0.5)*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1568;
   std::complex<double> tmp_1572;
   std::complex<double> tmp_1573;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1574;
      std::complex<double> tmp_1575;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1575 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1574 += tmp_1575;
      tmp_1573 += (ZD(gI1,j2)) * tmp_1574;
   }
   tmp_1572 += tmp_1573;
   tmp_1527 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1572;
   std::complex<double> tmp_1576;
   std::complex<double> tmp_1577;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1578;
      std::complex<double> tmp_1579;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1580;
         std::complex<double> tmp_1581;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1581 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1580 += tmp_1581;
         tmp_1579 += (ZD(gI1,3 + j2)) * tmp_1580;
      }
      tmp_1578 += tmp_1579;
      tmp_1577 += (Conj(ZD(gI2,3 + j3))) * tmp_1578;
   }
   tmp_1576 += tmp_1577;
   tmp_1527 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1576
      ;
   std::complex<double> tmp_1582;
   std::complex<double> tmp_1583;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1584;
      std::complex<double> tmp_1585;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1586;
         std::complex<double> tmp_1587;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1587 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1586 += tmp_1587;
         tmp_1585 += (Conj(ZD(gI2,j2))) * tmp_1586;
      }
      tmp_1584 += tmp_1585;
      tmp_1583 += (ZD(gI1,j3)) * tmp_1584;
   }
   tmp_1582 += tmp_1583;
   tmp_1527 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1582
      ;
   result += (std::complex<double>(0,-1)) * tmp_1527;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_1588;
   std::complex<double> tmp_1589;
   std::complex<double> tmp_1590;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1590 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1589 += tmp_1590;
   tmp_1588 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1589;
   std::complex<double> tmp_1591;
   std::complex<double> tmp_1592;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1592 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1591 += tmp_1592;
   tmp_1588 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1591;
   std::complex<double> tmp_1593;
   std::complex<double> tmp_1594;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1594 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1593 += tmp_1594;
   tmp_1588 += (std::complex<double>(0,-1)*QHd*Ql*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_1593;
   std::complex<double> tmp_1595;
   std::complex<double> tmp_1596;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1596 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1595 += tmp_1596;
   tmp_1588 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1595;
   std::complex<double> tmp_1597;
   std::complex<double> tmp_1598;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1598 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1597 += tmp_1598;
   tmp_1588 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1597;
   std::complex<double> tmp_1599;
   std::complex<double> tmp_1600;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1600 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1599 += tmp_1600;
   tmp_1588 += (std::complex<double>(0,-1)*QHu*Ql*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_1599;
   std::complex<double> tmp_1601;
   std::complex<double> tmp_1602;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1602 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1601 += tmp_1602;
   tmp_1588 += (std::complex<double>(0,-1)*Ql*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_1601;
   std::complex<double> tmp_1603;
   std::complex<double> tmp_1604;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1604 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1603 += tmp_1604;
   tmp_1588 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1603;
   std::complex<double> tmp_1605;
   std::complex<double> tmp_1606;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1606 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1605 += tmp_1606;
   tmp_1588 += (std::complex<double>(0,-1)*Qe*QHd*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_1605;
   std::complex<double> tmp_1607;
   std::complex<double> tmp_1608;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1608 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1607 += tmp_1608;
   tmp_1588 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1607;
   std::complex<double> tmp_1609;
   std::complex<double> tmp_1610;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1610 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1609 += tmp_1610;
   tmp_1588 += (std::complex<double>(0,-1)*Qe*QHu*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_1609;
   std::complex<double> tmp_1611;
   std::complex<double> tmp_1612;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1612 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1611 += tmp_1612;
   tmp_1588 += (std::complex<double>(0,-1)*Qe*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_1611;
   std::complex<double> tmp_1613;
   std::complex<double> tmp_1614;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1615;
      std::complex<double> tmp_1616;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1616 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1615 += tmp_1616;
      tmp_1614 += (Conj(ZE(gI2,j2))) * tmp_1615;
   }
   tmp_1613 += tmp_1614;
   tmp_1588 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(1,
      gO2)) * tmp_1613;
   std::complex<double> tmp_1617;
   std::complex<double> tmp_1618;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1619;
      std::complex<double> tmp_1620;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1620 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1619 += tmp_1620;
      tmp_1618 += (Conj(ZE(gI2,j2))) * tmp_1619;
   }
   tmp_1617 += tmp_1618;
   tmp_1588 += (std::complex<double>(0,0.5)*vu*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_1617;
   std::complex<double> tmp_1621;
   std::complex<double> tmp_1622;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1623;
      std::complex<double> tmp_1624;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1624 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1623 += tmp_1624;
      tmp_1622 += (Conj(ZE(gI2,j2))) * tmp_1623;
   }
   tmp_1621 += tmp_1622;
   tmp_1588 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1621;
   std::complex<double> tmp_1625;
   std::complex<double> tmp_1626;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1627;
      std::complex<double> tmp_1628;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1628 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1627 += tmp_1628;
      tmp_1626 += (ZE(gI1,j2)) * tmp_1627;
   }
   tmp_1625 += tmp_1626;
   tmp_1588 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_1625;
   std::complex<double> tmp_1629;
   std::complex<double> tmp_1630;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1631;
      std::complex<double> tmp_1632;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1632 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1631 += tmp_1632;
      tmp_1630 += (ZE(gI1,j2)) * tmp_1631;
   }
   tmp_1629 += tmp_1630;
   tmp_1588 += (std::complex<double>(0,0.5)*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1629;
   std::complex<double> tmp_1633;
   std::complex<double> tmp_1634;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1635;
      std::complex<double> tmp_1636;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1636 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1635 += tmp_1636;
      tmp_1634 += (ZE(gI1,j2)) * tmp_1635;
   }
   tmp_1633 += tmp_1634;
   tmp_1588 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1633;
   std::complex<double> tmp_1637;
   std::complex<double> tmp_1638;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1639;
      std::complex<double> tmp_1640;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1641;
         std::complex<double> tmp_1642;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1642 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1641 += tmp_1642;
         tmp_1640 += (ZE(gI1,3 + j2)) * tmp_1641;
      }
      tmp_1639 += tmp_1640;
      tmp_1638 += (Conj(ZE(gI2,3 + j3))) * tmp_1639;
   }
   tmp_1637 += tmp_1638;
   tmp_1588 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1637
      ;
   std::complex<double> tmp_1643;
   std::complex<double> tmp_1644;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1645;
      std::complex<double> tmp_1646;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1647;
         std::complex<double> tmp_1648;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1648 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1647 += tmp_1648;
         tmp_1646 += (Conj(ZE(gI2,j2))) * tmp_1647;
      }
      tmp_1645 += tmp_1646;
      tmp_1644 += (ZE(gI1,j3)) * tmp_1645;
   }
   tmp_1643 += tmp_1644;
   tmp_1588 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1643
      ;
   result += (std::complex<double>(0,-1)) * tmp_1588;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_1649;
   std::complex<double> tmp_1650;
   std::complex<double> tmp_1651;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1651 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1650 += tmp_1651;
   tmp_1649 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1650;
   std::complex<double> tmp_1652;
   std::complex<double> tmp_1653;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1653 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1652 += tmp_1653;
   tmp_1649 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1652;
   std::complex<double> tmp_1654;
   std::complex<double> tmp_1655;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1655 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1654 += tmp_1655;
   tmp_1649 += (std::complex<double>(0,-1)*QHd*Qq*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_1654;
   std::complex<double> tmp_1656;
   std::complex<double> tmp_1657;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1657 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1656 += tmp_1657;
   tmp_1649 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1656;
   std::complex<double> tmp_1658;
   std::complex<double> tmp_1659;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1659 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1658 += tmp_1659;
   tmp_1649 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1658;
   std::complex<double> tmp_1660;
   std::complex<double> tmp_1661;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1661 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1660 += tmp_1661;
   tmp_1649 += (std::complex<double>(0,-1)*QHu*Qq*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_1660;
   std::complex<double> tmp_1662;
   std::complex<double> tmp_1663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1663 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1662 += tmp_1663;
   tmp_1649 += (std::complex<double>(0,-1)*Qq*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_1662;
   std::complex<double> tmp_1664;
   std::complex<double> tmp_1665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1665 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1664 += tmp_1665;
   tmp_1649 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1664;
   std::complex<double> tmp_1666;
   std::complex<double> tmp_1667;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1667 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1666 += tmp_1667;
   tmp_1649 += (std::complex<double>(0,-1)*QHd*Qu*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_1666;
   std::complex<double> tmp_1668;
   std::complex<double> tmp_1669;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1669 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1668 += tmp_1669;
   tmp_1649 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_1668;
   std::complex<double> tmp_1670;
   std::complex<double> tmp_1671;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1671 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1670 += tmp_1671;
   tmp_1649 += (std::complex<double>(0,-1)*QHu*Qu*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_1670;
   std::complex<double> tmp_1672;
   std::complex<double> tmp_1673;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1673 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1672 += tmp_1673;
   tmp_1649 += (std::complex<double>(0,-1)*Qs*Qu*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_1672;
   std::complex<double> tmp_1674;
   std::complex<double> tmp_1675;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1676;
      std::complex<double> tmp_1677;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1677 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1676 += tmp_1677;
      tmp_1675 += (Conj(ZU(gI2,j2))) * tmp_1676;
   }
   tmp_1674 += tmp_1675;
   tmp_1649 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(0,
      gO2)) * tmp_1674;
   std::complex<double> tmp_1678;
   std::complex<double> tmp_1679;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1680;
      std::complex<double> tmp_1681;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1681 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1680 += tmp_1681;
      tmp_1679 += (Conj(ZU(gI2,j2))) * tmp_1680;
   }
   tmp_1678 += tmp_1679;
   tmp_1649 += (std::complex<double>(0,0.5)*vd*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_1678;
   std::complex<double> tmp_1682;
   std::complex<double> tmp_1683;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1684;
      std::complex<double> tmp_1685;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1685 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1684 += tmp_1685;
      tmp_1683 += (Conj(ZU(gI2,j2))) * tmp_1684;
   }
   tmp_1682 += tmp_1683;
   tmp_1649 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1682;
   std::complex<double> tmp_1686;
   std::complex<double> tmp_1687;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1688;
      std::complex<double> tmp_1689;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1689 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1688 += tmp_1689;
      tmp_1687 += (ZU(gI1,j2)) * tmp_1688;
   }
   tmp_1686 += tmp_1687;
   tmp_1649 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(0,gO2)*Lambdax) *
      tmp_1686;
   std::complex<double> tmp_1690;
   std::complex<double> tmp_1691;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1692;
      std::complex<double> tmp_1693;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1693 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1692 += tmp_1693;
      tmp_1691 += (ZU(gI1,j2)) * tmp_1692;
   }
   tmp_1690 += tmp_1691;
   tmp_1649 += (std::complex<double>(0,0.5)*vd*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1690;
   std::complex<double> tmp_1694;
   std::complex<double> tmp_1695;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1696;
      std::complex<double> tmp_1697;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1697 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1696 += tmp_1697;
      tmp_1695 += (ZU(gI1,j2)) * tmp_1696;
   }
   tmp_1694 += tmp_1695;
   tmp_1649 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1694;
   std::complex<double> tmp_1698;
   std::complex<double> tmp_1699;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1700;
      std::complex<double> tmp_1701;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1702;
         std::complex<double> tmp_1703;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1703 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1702 += tmp_1703;
         tmp_1701 += (ZU(gI1,3 + j2)) * tmp_1702;
      }
      tmp_1700 += tmp_1701;
      tmp_1699 += (Conj(ZU(gI2,3 + j3))) * tmp_1700;
   }
   tmp_1698 += tmp_1699;
   tmp_1649 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1698
      ;
   std::complex<double> tmp_1704;
   std::complex<double> tmp_1705;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1706;
      std::complex<double> tmp_1707;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1708;
         std::complex<double> tmp_1709;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1709 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1708 += tmp_1709;
         tmp_1707 += (Conj(ZU(gI2,j2))) * tmp_1708;
      }
      tmp_1706 += tmp_1707;
      tmp_1705 += (ZU(gI1,j3)) * tmp_1706;
   }
   tmp_1704 += tmp_1705;
   tmp_1649 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1704
      ;
   result += (std::complex<double>(0,-1)) * tmp_1649;

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(5*KroneckerDelta(2,gO2)*(-2*gp*Qs*ZN(gI1,5)*ZN(gI2,0) +
      1.4142135623730951*Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2,3) + ZN(gI1,3)*ZN(gI2,4))
      - 2*gp*Qs*ZN(gI1,0)*ZN(gI2,5)) + KroneckerDelta(0,gO2)*(ZN(gI1,3)*(-10*gp*
      QHd*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(gI2,2)) - 10*gp*QHd
      *ZN(gI1,0)*ZN(gI2,3) + 3.872983346207417*g1*ZN(gI1,1)*ZN(gI2,3) - 5*g2*ZN(
      gI1,2)*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,5)*ZN(gI2,4) +
      7.0710678118654755*Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,5)) - KroneckerDelta(1,gO2
      )*(ZN(gI1,4)*(10*gp*QHu*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN
      (gI2,2)) + (10*gp*QHu*ZN(gI1,0) + 3.872983346207417*g1*ZN(gI1,1) - 5*g2*ZN(
      gI1,2))*ZN(gI2,4) - 7.0710678118654755*Conj(Lambdax)*(ZN(gI1,5)*ZN(gI2,3) +
      ZN(gI1,3)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(3.872983346207417*g1*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*
      KroneckerDelta(0,gO1) - 5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,3))*KroneckerDelta(
      0,gO1) - 10*gp*QHu*Conj(ZN(gI1,4))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) -
      3.872983346207417*g1*Conj(ZN(gI1,4))*Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) +
      5*g2*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1) -
      3.872983346207417*g1*Conj(ZN(gI1,1))*Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) +
      5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,4))*KroneckerDelta(1,gO1) - 10*gp*Qs*Conj(
      ZN(gI1,5))*Conj(ZN(gI2,0))*KroneckerDelta(2,gO1) - 10*gp*Conj(ZN(gI1,0))*(
      QHd*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) + QHu*Conj(ZN(gI2,4))*
      KroneckerDelta(1,gO1) + Qs*Conj(ZN(gI2,5))*KroneckerDelta(2,gO1)) +
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

std::complex<double> CLASSNAME::CpUhhconjVWmHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) - KroneckerDelta(1,gO2)*ZP(
      gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZAh(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(10*gp*Qs*Conj(ZA(gI2,2))*
      KroneckerDelta(2,gO2)*Sin(ThetaWp()) + Conj(ZA(gI2,0))*KroneckerDelta(0,gO2)
      *(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) - Conj(ZA(gI2,1))*KroneckerDelta(1
      ,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp(
      ))*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZpAh(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(10*gp*Qs*Conj(ZA(gI2,2))*Cos(ThetaWp(
      ))*KroneckerDelta(2,gO2) + Conj(ZA(gI2,1))*KroneckerDelta(1,gO2)*(10*gp*QHu*
      Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*
      Sin(ThetaW())*Sin(ThetaWp())) + Conj(ZA(gI2,0))*KroneckerDelta(0,gO2)*(10*gp
      *QHd*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW(
      )))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(vd*KroneckerDelta(0,gO1) - vu*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(vd*KroneckerDelta(0,gO1) - vu*
      KroneckerDelta(1,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)
      *Sqr(Sin(ThetaWp())) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(20*g2*gp
      *QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*
      QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*
      Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-20*g2*
      gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*
      gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)
      *Sqr(Cos(ThetaWp())) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-2*gp*
      QHd*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp(
      )) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(2*gp*QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin
      (2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*
      (7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(
      Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*((AbsSqr(
      Lambdax) + QHd*Qs*Sqr(gp))*ZP(gI1,0)*ZP(gI2,0) + (AbsSqr(Lambdax) + QHu*Qs*
      Sqr(gp))*ZP(gI1,1)*ZP(gI2,1)) - KroneckerDelta(0,gO1)*(-5*KroneckerDelta(1,
      gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,
      1)) + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*
      ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI1,1
      )*ZP(gI2,1))) + KroneckerDelta(1,gO1)*(5*KroneckerDelta(0,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp)))*ZP(gI1,
      0)*ZP(gI2,0) - (3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHu)))*ZP(gI1,1)*ZP(
      gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(vu*KroneckerDelta(0,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) - ZP(gI1,0)*ZP(gI2,1)) + vd*
      KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) -
      ZP(gI1,0)*ZP(gI2,1)) + 2.8284271247461903*KroneckerDelta(2,gO2)*(-(TLambdax*
      ZP(gI1,1)*ZP(gI2,0)) + Conj(TLambdax)*ZP(gI1,0)*ZP(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.7071067811865475)*(g2*KroneckerDelta(0,
      gO2)*UM(gI1,1)*UP(gI2,0) + (g2*KroneckerDelta(1,gO2)*UM(gI1,0) - Conj(
      Lambdax)*KroneckerDelta(2,gO2)*UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.7071067811865475)*(g2*Conj(UM(gI2,0))*Conj
      (UP(gI1,1))*KroneckerDelta(1,gO1) + Conj(UM(gI2,1))*(g2*Conj(UP(gI1,0))*
      KroneckerDelta(0,gO1) - Conj(UP(gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(Conj(ZA(gI1,0))*(20*Conj(ZA(gI2,2))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2))*(AbsSqr
      (Lambdax) + QHd*Qs*Sqr(gp)) - Conj(ZA(gI2,1))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + Conj(ZA(gI2
      ,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*
      Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 3*KroneckerDelta(
      0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))
      )) + Conj(ZA(gI1,1))*(-20*Conj(ZA(gI2,2))*(KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1) + KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2))*(AbsSqr
      (Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZA(gI2,0))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + Conj(ZA(gI2
      ,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*
      Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - 3*KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))))
      - 20*Conj(ZA(gI1,2))*(Conj(ZA(gI2,0))*(KroneckerDelta(0,gO2)*KroneckerDelta
      (2,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) +
      QHd*Qs*Sqr(gp)) + Conj(ZA(gI2,1))*(KroneckerDelta(1,gO2)*KroneckerDelta(2,
      gO1) + KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) + QHu*
      Qs*Sqr(gp)) + Conj(ZA(gI2,2))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(
      AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = -0.05*KroneckerDelta(gI1,gI2)*(20*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*
      (-3*Sqr(g1) - 5*Sqr(g2) + 20*QHu*Ql*Sqr(gp)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*QHd*Ql*Sqr(gp))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5
      *(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))*ZH(gI1,0)*ZH(gI2,0) + (20*AbsSqr(Lambdax) -
      3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZH(gI1,1)*ZH(gI2,1) + 20*(
      AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,2)*ZH(gI2,2))) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) -
      20*QHd*QHu*Sqr(gp))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp
      )*Sqr(QHu))*ZH(gI1,1)*ZH(gI2,1) - 20*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZH(
      gI1,2)*ZH(gI2,2)) - 20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*((AbsSqr(
      Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,0)*ZH(gI2,0) + (AbsSqr(Lambdax) + QHu*Qs*
      Sqr(gp))*ZH(gI1,1)*ZH(gI2,1) + Sqr(gp)*Sqr(Qs)*ZH(gI1,2)*ZH(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.35355339059327373)*(Conj(ZA(gI1,2))*(Conj(
      ZA(gI2,1))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*KroneckerDelta(1,gO2)) +
      Conj(ZA(gI1,1))*(Conj(ZA(gI2,2))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*
      KroneckerDelta(2,gO2)) + Conj(ZA(gI1,0))*(Conj(ZA(gI2,2))*KroneckerDelta(1,
      gO2) + Conj(ZA(gI2,1))*KroneckerDelta(2,gO2)))*(Conj(TLambdax) - TLambdax);

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(Conj(ZA(gI2,0))*(7.0710678118654755*(Conj(TLambdax) +
      TLambdax)*(KroneckerDelta(2,gO2)*ZH(gI1,1) + KroneckerDelta(1,gO2)*ZH(gI1,2)
      ) + KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))
      *ZH(gI1,0) + vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr
      (gp))*ZH(gI1,1) + 20*vS*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,2)))) +
      Conj(ZA(gI2,1))*(-7.0710678118654755*(Conj(TLambdax) + TLambdax)*(
      KroneckerDelta(2,gO2)*ZH(gI1,0) + KroneckerDelta(0,gO2)*ZH(gI1,2)) +
      KroneckerDelta(1,gO2)*(vd*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*
      QHd*QHu*Sqr(gp))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu)
      )*ZH(gI1,1) - 20*vS*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZH(gI1,2))) - 5*Conj(
      ZA(gI2,2))*(1.4142135623730951*Conj(TLambdax)*(KroneckerDelta(1,gO2)*ZH(gI1,
      0) + KroneckerDelta(0,gO2)*ZH(gI1,1)) + 1.4142135623730951*TLambdax*(
      KroneckerDelta(1,gO2)*ZH(gI1,0) + KroneckerDelta(0,gO2)*ZH(gI1,1)) + 4*
      KroneckerDelta(2,gO2)*(vd*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,0) + vu*
      (AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))*ZH(gI1,1) + vS*Sqr(gp)*Sqr(Qs)*ZH(gI1,2))
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.35355339059327373)*(Conj(TLambdax) -
      TLambdax)*(KroneckerDelta(2,gO2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1))
      + KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) +
      KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1710;
   std::complex<double> tmp_1711;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1712;
      std::complex<double> tmp_1713;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1713 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1712 += tmp_1713;
      tmp_1711 += (ZDL(gI1,j2)) * tmp_1712;
   }
   tmp_1710 += tmp_1711;
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2))
      * tmp_1710;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1714;
   std::complex<double> tmp_1715;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1716;
      std::complex<double> tmp_1717;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1717 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1716 += tmp_1717;
      tmp_1715 += (Conj(ZDL(gI2,j2))) * tmp_1716;
   }
   tmp_1714 += tmp_1715;
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_1714;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1718;
   std::complex<double> tmp_1719;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1720;
      std::complex<double> tmp_1721;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1721 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1720 += tmp_1721;
      tmp_1719 += (ZEL(gI1,j2)) * tmp_1720;
   }
   tmp_1718 += tmp_1719;
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2))
      * tmp_1718;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1722;
   std::complex<double> tmp_1723;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1724;
      std::complex<double> tmp_1725;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1725 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1724 += tmp_1725;
      tmp_1723 += (Conj(ZEL(gI2,j2))) * tmp_1724;
   }
   tmp_1722 += tmp_1723;
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_1722;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1726;
   std::complex<double> tmp_1727;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1728;
      std::complex<double> tmp_1729;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1729 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1728 += tmp_1729;
      tmp_1727 += (ZUL(gI1,j2)) * tmp_1728;
   }
   tmp_1726 += tmp_1727;
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2))
      * tmp_1726;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1730;
   std::complex<double> tmp_1731;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1732;
      std::complex<double> tmp_1733;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1733 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1732 += tmp_1733;
      tmp_1731 += (Conj(ZUL(gI2,j2))) * tmp_1732;
   }
   tmp_1730 += tmp_1731;
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,gO1)
      ) * tmp_1730;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_1734;
   std::complex<double> tmp_1735;
   std::complex<double> tmp_1736;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1736 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1735 += tmp_1736;
   tmp_1734 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1735;
   std::complex<double> tmp_1737;
   std::complex<double> tmp_1738;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1738 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1737 += tmp_1738;
   tmp_1734 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1737;
   std::complex<double> tmp_1739;
   std::complex<double> tmp_1740;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1740 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1739 += tmp_1740;
   tmp_1734 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1739;
   std::complex<double> tmp_1741;
   std::complex<double> tmp_1742;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1742 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1741 += tmp_1742;
   tmp_1734 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1741;
   std::complex<double> tmp_1743;
   std::complex<double> tmp_1744;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1744 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1743 += tmp_1744;
   tmp_1734 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1743;
   std::complex<double> tmp_1745;
   std::complex<double> tmp_1746;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1746 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1745 += tmp_1746;
   tmp_1734 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1745;
   std::complex<double> tmp_1747;
   std::complex<double> tmp_1748;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1748 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1747 += tmp_1748;
   tmp_1734 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1747;
   std::complex<double> tmp_1749;
   std::complex<double> tmp_1750;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1750 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1749 += tmp_1750;
   tmp_1734 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1749;
   std::complex<double> tmp_1751;
   std::complex<double> tmp_1752;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1752 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1751 += tmp_1752;
   tmp_1734 += (std::complex<double>(0,-1)*Qd*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1751;
   std::complex<double> tmp_1753;
   std::complex<double> tmp_1754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1754 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1753 += tmp_1754;
   tmp_1734 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1753;
   std::complex<double> tmp_1755;
   std::complex<double> tmp_1756;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1756 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1755 += tmp_1756;
   tmp_1734 += (std::complex<double>(0,-1)*Qd*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1755;
   std::complex<double> tmp_1757;
   std::complex<double> tmp_1758;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1758 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1757 += tmp_1758;
   tmp_1734 += (std::complex<double>(0,-1)*Qd*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1757;
   std::complex<double> tmp_1759;
   std::complex<double> tmp_1760;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1761;
      std::complex<double> tmp_1762;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1762 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1761 += tmp_1762;
      tmp_1760 += (Conj(ZD(gI2,j2))) * tmp_1761;
   }
   tmp_1759 += tmp_1760;
   tmp_1734 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)) * tmp_1759;
   std::complex<double> tmp_1763;
   std::complex<double> tmp_1764;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1765;
      std::complex<double> tmp_1766;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1766 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1765 += tmp_1766;
      tmp_1764 += (Conj(ZD(gI2,j2))) * tmp_1765;
   }
   tmp_1763 += tmp_1764;
   tmp_1734 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1
      )*KroneckerDelta(2,gO2)) * tmp_1763;
   std::complex<double> tmp_1767;
   std::complex<double> tmp_1768;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1769;
      std::complex<double> tmp_1770;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1770 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1769 += tmp_1770;
      tmp_1768 += (ZD(gI1,j2)) * tmp_1769;
   }
   tmp_1767 += tmp_1768;
   tmp_1734 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1767;
   std::complex<double> tmp_1771;
   std::complex<double> tmp_1772;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1773;
      std::complex<double> tmp_1774;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1774 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1773 += tmp_1774;
      tmp_1772 += (ZD(gI1,j2)) * tmp_1773;
   }
   tmp_1771 += tmp_1772;
   tmp_1734 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1771;
   std::complex<double> tmp_1775;
   std::complex<double> tmp_1776;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1777;
      std::complex<double> tmp_1778;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1779;
         std::complex<double> tmp_1780;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1780 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1779 += tmp_1780;
         tmp_1778 += (ZD(gI1,3 + j2)) * tmp_1779;
      }
      tmp_1777 += tmp_1778;
      tmp_1776 += (Conj(ZD(gI2,3 + j3))) * tmp_1777;
   }
   tmp_1775 += tmp_1776;
   tmp_1734 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1775;
   std::complex<double> tmp_1781;
   std::complex<double> tmp_1782;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1783;
      std::complex<double> tmp_1784;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1785;
         std::complex<double> tmp_1786;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1786 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1785 += tmp_1786;
         tmp_1784 += (Conj(ZD(gI2,j2))) * tmp_1785;
      }
      tmp_1783 += tmp_1784;
      tmp_1782 += (ZD(gI1,j3)) * tmp_1783;
   }
   tmp_1781 += tmp_1782;
   tmp_1734 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1781;
   result += (std::complex<double>(0,-1)) * tmp_1734;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_1787;
   std::complex<double> tmp_1788;
   std::complex<double> tmp_1789;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1789 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1788 += tmp_1789;
   tmp_1787 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1788;
   std::complex<double> tmp_1790;
   std::complex<double> tmp_1791;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1791 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1790 += tmp_1791;
   tmp_1787 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1790;
   std::complex<double> tmp_1792;
   std::complex<double> tmp_1793;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1793 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1792 += tmp_1793;
   tmp_1787 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1792;
   std::complex<double> tmp_1794;
   std::complex<double> tmp_1795;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1795 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1794 += tmp_1795;
   tmp_1787 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1794;
   std::complex<double> tmp_1796;
   std::complex<double> tmp_1797;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1797 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1796 += tmp_1797;
   tmp_1787 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1796;
   std::complex<double> tmp_1798;
   std::complex<double> tmp_1799;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1799 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1798 += tmp_1799;
   tmp_1787 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1798;
   std::complex<double> tmp_1800;
   std::complex<double> tmp_1801;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1801 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1800 += tmp_1801;
   tmp_1787 += (std::complex<double>(0,-1)*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1800;
   std::complex<double> tmp_1802;
   std::complex<double> tmp_1803;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1803 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1802 += tmp_1803;
   tmp_1787 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1802;
   std::complex<double> tmp_1804;
   std::complex<double> tmp_1805;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1805 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1804 += tmp_1805;
   tmp_1787 += (std::complex<double>(0,-1)*Qe*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1804;
   std::complex<double> tmp_1806;
   std::complex<double> tmp_1807;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1807 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1806 += tmp_1807;
   tmp_1787 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1806;
   std::complex<double> tmp_1808;
   std::complex<double> tmp_1809;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1809 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1808 += tmp_1809;
   tmp_1787 += (std::complex<double>(0,-1)*Qe*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1808;
   std::complex<double> tmp_1810;
   std::complex<double> tmp_1811;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1811 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1810 += tmp_1811;
   tmp_1787 += (std::complex<double>(0,-1)*Qe*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1810;
   std::complex<double> tmp_1812;
   std::complex<double> tmp_1813;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1814;
      std::complex<double> tmp_1815;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1815 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1814 += tmp_1815;
      tmp_1813 += (Conj(ZE(gI2,j2))) * tmp_1814;
   }
   tmp_1812 += tmp_1813;
   tmp_1787 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)) * tmp_1812;
   std::complex<double> tmp_1816;
   std::complex<double> tmp_1817;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1818;
      std::complex<double> tmp_1819;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1819 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1818 += tmp_1819;
      tmp_1817 += (Conj(ZE(gI2,j2))) * tmp_1818;
   }
   tmp_1816 += tmp_1817;
   tmp_1787 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1
      )*KroneckerDelta(2,gO2)) * tmp_1816;
   std::complex<double> tmp_1820;
   std::complex<double> tmp_1821;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1822;
      std::complex<double> tmp_1823;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1823 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1822 += tmp_1823;
      tmp_1821 += (ZE(gI1,j2)) * tmp_1822;
   }
   tmp_1820 += tmp_1821;
   tmp_1787 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1820;
   std::complex<double> tmp_1824;
   std::complex<double> tmp_1825;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1826;
      std::complex<double> tmp_1827;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1827 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1826 += tmp_1827;
      tmp_1825 += (ZE(gI1,j2)) * tmp_1826;
   }
   tmp_1824 += tmp_1825;
   tmp_1787 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1824;
   std::complex<double> tmp_1828;
   std::complex<double> tmp_1829;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1830;
      std::complex<double> tmp_1831;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1832;
         std::complex<double> tmp_1833;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1833 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1832 += tmp_1833;
         tmp_1831 += (ZE(gI1,3 + j2)) * tmp_1832;
      }
      tmp_1830 += tmp_1831;
      tmp_1829 += (Conj(ZE(gI2,3 + j3))) * tmp_1830;
   }
   tmp_1828 += tmp_1829;
   tmp_1787 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1828;
   std::complex<double> tmp_1834;
   std::complex<double> tmp_1835;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1836;
      std::complex<double> tmp_1837;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1838;
         std::complex<double> tmp_1839;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1839 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1838 += tmp_1839;
         tmp_1837 += (Conj(ZE(gI2,j2))) * tmp_1838;
      }
      tmp_1836 += tmp_1837;
      tmp_1835 += (ZE(gI1,j3)) * tmp_1836;
   }
   tmp_1834 += tmp_1835;
   tmp_1787 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1834;
   result += (std::complex<double>(0,-1)) * tmp_1787;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_1840;
   std::complex<double> tmp_1841;
   std::complex<double> tmp_1842;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1842 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1841 += tmp_1842;
   tmp_1840 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1841;
   std::complex<double> tmp_1843;
   std::complex<double> tmp_1844;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1844 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1843 += tmp_1844;
   tmp_1840 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1843;
   std::complex<double> tmp_1845;
   std::complex<double> tmp_1846;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1846 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1845 += tmp_1846;
   tmp_1840 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1845;
   std::complex<double> tmp_1847;
   std::complex<double> tmp_1848;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1848 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1847 += tmp_1848;
   tmp_1840 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1847;
   std::complex<double> tmp_1849;
   std::complex<double> tmp_1850;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1850 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1849 += tmp_1850;
   tmp_1840 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1849;
   std::complex<double> tmp_1851;
   std::complex<double> tmp_1852;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1852 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1851 += tmp_1852;
   tmp_1840 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1851;
   std::complex<double> tmp_1853;
   std::complex<double> tmp_1854;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1854 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1853 += tmp_1854;
   tmp_1840 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1853;
   std::complex<double> tmp_1855;
   std::complex<double> tmp_1856;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1856 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1855 += tmp_1856;
   tmp_1840 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1855;
   std::complex<double> tmp_1857;
   std::complex<double> tmp_1858;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1858 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1857 += tmp_1858;
   tmp_1840 += (std::complex<double>(0,-1)*QHd*Qu*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_1857;
   std::complex<double> tmp_1859;
   std::complex<double> tmp_1860;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1860 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1859 += tmp_1860;
   tmp_1840 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1859;
   std::complex<double> tmp_1861;
   std::complex<double> tmp_1862;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1862 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1861 += tmp_1862;
   tmp_1840 += (std::complex<double>(0,-1)*QHu*Qu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_1861;
   std::complex<double> tmp_1863;
   std::complex<double> tmp_1864;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1864 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1863 += tmp_1864;
   tmp_1840 += (std::complex<double>(0,-1)*Qs*Qu*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_1863;
   std::complex<double> tmp_1865;
   std::complex<double> tmp_1866;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1867;
      std::complex<double> tmp_1868;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1868 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1867 += tmp_1868;
      tmp_1866 += (Conj(ZU(gI2,j2))) * tmp_1867;
   }
   tmp_1865 += tmp_1866;
   tmp_1840 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2
      )*KroneckerDelta(2,gO1)) * tmp_1865;
   std::complex<double> tmp_1869;
   std::complex<double> tmp_1870;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1871;
      std::complex<double> tmp_1872;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1872 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1871 += tmp_1872;
      tmp_1870 += (Conj(ZU(gI2,j2))) * tmp_1871;
   }
   tmp_1869 += tmp_1870;
   tmp_1840 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1
      )*KroneckerDelta(2,gO2)) * tmp_1869;
   std::complex<double> tmp_1873;
   std::complex<double> tmp_1874;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1875;
      std::complex<double> tmp_1876;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1876 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1875 += tmp_1876;
      tmp_1874 += (ZU(gI1,j2)) * tmp_1875;
   }
   tmp_1873 += tmp_1874;
   tmp_1840 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_1873;
   std::complex<double> tmp_1877;
   std::complex<double> tmp_1878;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1879;
      std::complex<double> tmp_1880;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1880 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1879 += tmp_1880;
      tmp_1878 += (ZU(gI1,j2)) * tmp_1879;
   }
   tmp_1877 += tmp_1878;
   tmp_1840 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_1877;
   std::complex<double> tmp_1881;
   std::complex<double> tmp_1882;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1883;
      std::complex<double> tmp_1884;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1885;
         std::complex<double> tmp_1886;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1886 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1885 += tmp_1886;
         tmp_1884 += (ZU(gI1,3 + j2)) * tmp_1885;
      }
      tmp_1883 += tmp_1884;
      tmp_1882 += (Conj(ZU(gI2,3 + j3))) * tmp_1883;
   }
   tmp_1881 += tmp_1882;
   tmp_1840 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1881;
   std::complex<double> tmp_1887;
   std::complex<double> tmp_1888;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1889;
      std::complex<double> tmp_1890;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1891;
         std::complex<double> tmp_1892;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1892 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1891 += tmp_1892;
         tmp_1890 += (Conj(ZU(gI2,j2))) * tmp_1891;
      }
      tmp_1889 += tmp_1890;
      tmp_1888 += (ZU(gI1,j3)) * tmp_1889;
   }
   tmp_1887 += tmp_1888;
   tmp_1840 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1887;
   result += (std::complex<double>(0,-1)) * tmp_1840;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1893;
   std::complex<double> tmp_1894;
   std::complex<double> tmp_1895;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1896;
      std::complex<double> tmp_1897;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1897 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1896 += tmp_1897;
      tmp_1895 += (Conj(ZD(gI2,j2))) * tmp_1896;
   }
   tmp_1894 += tmp_1895;
   tmp_1893 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(1,gO2)) * tmp_1894;
   std::complex<double> tmp_1898;
   std::complex<double> tmp_1899;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1900;
      std::complex<double> tmp_1901;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1901 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1900 += tmp_1901;
      tmp_1899 += (Conj(ZD(gI2,j2))) * tmp_1900;
   }
   tmp_1898 += tmp_1899;
   tmp_1893 += (0.5*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_1898;
   std::complex<double> tmp_1902;
   std::complex<double> tmp_1903;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1904;
      std::complex<double> tmp_1905;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1905 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1904 += tmp_1905;
      tmp_1903 += (Conj(ZD(gI2,j2))) * tmp_1904;
   }
   tmp_1902 += tmp_1903;
   tmp_1893 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1902;
   std::complex<double> tmp_1906;
   std::complex<double> tmp_1907;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1908;
      std::complex<double> tmp_1909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1909 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1908 += tmp_1909;
      tmp_1907 += (ZD(gI1,j2)) * tmp_1908;
   }
   tmp_1906 += tmp_1907;
   tmp_1893 += (-0.5*vS*KroneckerDelta(1,gO2)*Lambdax) * tmp_1906;
   std::complex<double> tmp_1910;
   std::complex<double> tmp_1911;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1912;
      std::complex<double> tmp_1913;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1913 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1912 += tmp_1913;
      tmp_1911 += (ZD(gI1,j2)) * tmp_1912;
   }
   tmp_1910 += tmp_1911;
   tmp_1893 += (-0.5*vu*KroneckerDelta(2,gO2)*Lambdax) * tmp_1910;
   std::complex<double> tmp_1914;
   std::complex<double> tmp_1915;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1916;
      std::complex<double> tmp_1917;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1917 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1916 += tmp_1917;
      tmp_1915 += (ZD(gI1,j2)) * tmp_1916;
   }
   tmp_1914 += tmp_1915;
   tmp_1893 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1914;
   result += (std::complex<double>(0,-1)) * tmp_1893;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1918;
   std::complex<double> tmp_1919;
   std::complex<double> tmp_1920;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1921;
      std::complex<double> tmp_1922;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1922 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1921 += tmp_1922;
      tmp_1920 += (Conj(ZE(gI2,j2))) * tmp_1921;
   }
   tmp_1919 += tmp_1920;
   tmp_1918 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(1,gO2)) * tmp_1919;
   std::complex<double> tmp_1923;
   std::complex<double> tmp_1924;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1925;
      std::complex<double> tmp_1926;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1926 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1925 += tmp_1926;
      tmp_1924 += (Conj(ZE(gI2,j2))) * tmp_1925;
   }
   tmp_1923 += tmp_1924;
   tmp_1918 += (0.5*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_1923;
   std::complex<double> tmp_1927;
   std::complex<double> tmp_1928;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1929;
      std::complex<double> tmp_1930;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1930 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1929 += tmp_1930;
      tmp_1928 += (Conj(ZE(gI2,j2))) * tmp_1929;
   }
   tmp_1927 += tmp_1928;
   tmp_1918 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1927;
   std::complex<double> tmp_1931;
   std::complex<double> tmp_1932;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1933;
      std::complex<double> tmp_1934;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1934 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1933 += tmp_1934;
      tmp_1932 += (ZE(gI1,j2)) * tmp_1933;
   }
   tmp_1931 += tmp_1932;
   tmp_1918 += (-0.5*vS*KroneckerDelta(1,gO2)*Lambdax) * tmp_1931;
   std::complex<double> tmp_1935;
   std::complex<double> tmp_1936;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1937;
      std::complex<double> tmp_1938;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1938 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1937 += tmp_1938;
      tmp_1936 += (ZE(gI1,j2)) * tmp_1937;
   }
   tmp_1935 += tmp_1936;
   tmp_1918 += (-0.5*vu*KroneckerDelta(2,gO2)*Lambdax) * tmp_1935;
   std::complex<double> tmp_1939;
   std::complex<double> tmp_1940;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1941;
      std::complex<double> tmp_1942;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1942 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1941 += tmp_1942;
      tmp_1940 += (ZE(gI1,j2)) * tmp_1941;
   }
   tmp_1939 += tmp_1940;
   tmp_1918 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1939;
   result += (std::complex<double>(0,-1)) * tmp_1918;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1943;
   std::complex<double> tmp_1944;
   std::complex<double> tmp_1945;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1946;
      std::complex<double> tmp_1947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1947 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1946 += tmp_1947;
      tmp_1945 += (Conj(ZU(gI2,j2))) * tmp_1946;
   }
   tmp_1944 += tmp_1945;
   tmp_1943 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(0,gO2)) * tmp_1944;
   std::complex<double> tmp_1948;
   std::complex<double> tmp_1949;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1950;
      std::complex<double> tmp_1951;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1951 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1950 += tmp_1951;
      tmp_1949 += (Conj(ZU(gI2,j2))) * tmp_1950;
   }
   tmp_1948 += tmp_1949;
   tmp_1943 += (0.5*vd*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_1948;
   std::complex<double> tmp_1952;
   std::complex<double> tmp_1953;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1954;
      std::complex<double> tmp_1955;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1955 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1954 += tmp_1955;
      tmp_1953 += (Conj(ZU(gI2,j2))) * tmp_1954;
   }
   tmp_1952 += tmp_1953;
   tmp_1943 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1952;
   std::complex<double> tmp_1956;
   std::complex<double> tmp_1957;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1958;
      std::complex<double> tmp_1959;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1959 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1958 += tmp_1959;
      tmp_1957 += (ZU(gI1,j2)) * tmp_1958;
   }
   tmp_1956 += tmp_1957;
   tmp_1943 += (-0.5*vS*KroneckerDelta(0,gO2)*Lambdax) * tmp_1956;
   std::complex<double> tmp_1960;
   std::complex<double> tmp_1961;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1962;
      std::complex<double> tmp_1963;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1963 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1962 += tmp_1963;
      tmp_1961 += (ZU(gI1,j2)) * tmp_1962;
   }
   tmp_1960 += tmp_1961;
   tmp_1943 += (-0.5*vd*KroneckerDelta(2,gO2)*Lambdax) * tmp_1960;
   std::complex<double> tmp_1964;
   std::complex<double> tmp_1965;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1966;
      std::complex<double> tmp_1967;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1967 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1966 += tmp_1967;
      tmp_1965 += (ZU(gI1,j2)) * tmp_1966;
   }
   tmp_1964 += tmp_1965;
   tmp_1943 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1964;
   result += (std::complex<double>(0,-1)) * tmp_1943;

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*KroneckerDelta(2,gO2)*(2*gp*Qs*ZN(
      gI1,5)*ZN(gI2,0) + 1.4142135623730951*Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2,3) +
      ZN(gI1,3)*ZN(gI2,4)) + 2*gp*Qs*ZN(gI1,0)*ZN(gI2,5)) + KroneckerDelta(0,gO2)*
      (ZN(gI1,3)*(10*gp*QHd*ZN(gI2,0) - 3.872983346207417*g1*ZN(gI2,1) + 5*g2*ZN(
      gI2,2)) + 10*gp*QHd*ZN(gI1,0)*ZN(gI2,3) - 3.872983346207417*g1*ZN(gI1,1)*ZN(
      gI2,3) + 5*g2*ZN(gI1,2)*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,
      5)*ZN(gI2,4) + 7.0710678118654755*Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,5)) +
      KroneckerDelta(1,gO2)*(ZN(gI1,4)*(10*gp*QHu*ZN(gI2,0) + 3.872983346207417*g1
      *ZN(gI2,1) - 5*g2*ZN(gI2,2)) + (10*gp*QHu*ZN(gI1,0) + 3.872983346207417*g1*
      ZN(gI1,1) - 5*g2*ZN(gI1,2))*ZN(gI2,4) + 7.0710678118654755*Conj(Lambdax)*(ZN
      (gI1,5)*ZN(gI2,3) + ZN(gI1,3)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(-3.872983346207417*g1*Conj(ZN(gI1,1))*
      Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) + 5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,3))
      *KroneckerDelta(0,gO1) + 10*gp*QHu*Conj(ZN(gI1,4))*Conj(ZN(gI2,0))*
      KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,4))*Conj(ZN(gI2,1))
      *KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*KroneckerDelta
      (1,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,1))*Conj(ZN(gI2,4))*
      KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,2))*Conj(ZN(gI2,4))*KroneckerDelta(
      1,gO1) + 10*gp*Qs*Conj(ZN(gI1,5))*Conj(ZN(gI2,0))*KroneckerDelta(2,gO1) + 10
      *gp*Conj(ZN(gI1,0))*(QHd*Conj(ZN(gI2,3))*KroneckerDelta(0,gO1) + QHu*Conj(ZN
      (gI2,4))*KroneckerDelta(1,gO1) + Qs*Conj(ZN(gI2,5))*KroneckerDelta(2,gO1)) +
      7.0710678118654755*Conj(ZN(gI1,5))*Conj(ZN(gI2,4))*KroneckerDelta(0,gO1)*
      Lambdax + 7.0710678118654755*Conj(ZN(gI1,4))*Conj(ZN(gI2,5))*KroneckerDelta(
      0,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,5))*Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1)*Lambdax + 7.0710678118654755*Conj(ZN(gI1,4))*Conj(ZN(
      gI2,3))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,3))*(10*gp*QHd*Conj(ZN(
      gI2,0))*KroneckerDelta(0,gO1) - 3.872983346207417*g1*Conj(ZN(gI2,1))*
      KroneckerDelta(0,gO1) + 5*g2*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) +
      7.0710678118654755*Conj(ZN(gI2,5))*KroneckerDelta(1,gO1)*Lambdax +
      7.0710678118654755*Conj(ZN(gI2,4))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjVWmHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhVZhh(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(KroneckerDelta(0,gO2)*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      + 10*gp*QHd*Sin(ThetaWp()))*ZH(gI2,0) - KroneckerDelta(1,gO2)*(5*g2*Cos(
      ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())
      - 10*gp*QHu*Sin(ThetaWp()))*ZH(gI2,1) + 10*gp*Qs*KroneckerDelta(2,gO2)*Sin(
      ThetaWp())*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpUAhVZphh(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(
      ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()))*ZH(gI2,0) + KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + 5*
      g2*Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(
      ThetaWp()))*ZH(gI2,1) + 10*gp*Qs*Cos(ThetaWp())*KroneckerDelta(2,gO2)*ZH(gI2
      ,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVP(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd*KroneckerDelta(0,gO2) -
      vu*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZVWm(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*g2*(vd*KroneckerDelta(0,gO2)*(3.872983346207417*g1*Cos(ThetaWp(
      ))*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) + vu*KroneckerDelta(1,gO2)*(
      -3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHu*Sin(ThetaWp()
      )));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZpVWm(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*g2*(vd*KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(ThetaWp()) -
      3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp())) + vu*KroneckerDelta(1,gO2
      )*(10*gp*QHu*Cos(ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp
      ())));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargWmCgZ(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()
      ) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHd*Sin(
      ThetaWp())) + vu*KroneckerDelta(1,gO1)*(-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp()
      )));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgWmCbargZ(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp(
      )) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(
      ThetaWp())) - vu*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargWmCgZp(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO1)*(10*gp*QHd*Cos(ThetaWp()) + (5*
      g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) + vu*
      KroneckerDelta(1,gO1)*(10*gp*QHu*Cos(ThetaWp()) + (-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgWmCbargZp(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO2)*(-10*gp*QHd*Cos(ThetaWp()) + (5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) - vu*
      KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargZgWm(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp(
      )) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(
      ThetaWp())) - vu*KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgZbargWm(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()
      ) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHd*Sin(
      ThetaWp())) + vu*KroneckerDelta(1,gO2)*(-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp()
      )));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargZpgWm(unsigned gO1) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.05*g2*(vd*KroneckerDelta(0,gO1)*(-10*gp*QHd*Cos(ThetaWp()) + (5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) - vu*
      KroneckerDelta(1,gO1)*(10*gp*QHu*Cos(ThetaWp()) + (5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgZpbargWm(unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(ThetaWp()) + (5*
      g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())) + vu*
      KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(
      15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) - 2
      *g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos
      (ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 20*Sqr(gp)
      *Sqr(QHd)*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)
      *(-15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())
      - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1*Cos(ThetaWp())*Sin
      (ThetaW()) - 10*gp*QHu*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
      Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 20*Sqr(
      gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(2*gp*QHd*(5*g2*
      Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*
      Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (-7.745966692414834*g1*g2*Cos(ThetaW(
      ))*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()
      )))*Sqr(Sin(ThetaWp()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-2*
      gp*QHu*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*
      ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (-7.745966692414834*
      g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)
      *Sqr(Sin(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*(20*AbsSqr(
      Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI1,0)*ZP(gI2,1) +
      KroneckerDelta(0,gO2)*(2*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(
      gI1,0)*ZP(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*
      Sqr(gp))*ZP(gI1,1)*ZP(gI2,1)))) + KroneckerDelta(1,gO1)*(KroneckerDelta(0,
      gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZP(
      gI1,1)*ZP(gI2,0) + KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) +
      5*Sqr(g2) - 20*QHd*QHu*Sqr(gp))*ZP(gI1,0)*ZP(gI2,0) - 2*(3*Sqr(g1) + 5*Sqr(
      g2) + 20*Sqr(gp)*Sqr(QHu))*ZP(gI1,1)*ZP(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(vu*Conj(ZA(gI2,0))*(-2*AbsSqr(Lambdax
      ) + Sqr(g2))*(KroneckerDelta(1,gO2)*ZP(gI1,0) - KroneckerDelta(0,gO2)*ZP(gI1
      ,1)) + vd*Conj(ZA(gI2,1))*(-2*AbsSqr(Lambdax) + Sqr(g2))*(KroneckerDelta(1,
      gO2)*ZP(gI1,0) - KroneckerDelta(0,gO2)*ZP(gI1,1)) + 2.8284271247461903*Conj(
      ZA(gI2,2))*(-(KroneckerDelta(1,gO2)*TLambdax*ZP(gI1,0)) + Conj(TLambdax)*
      KroneckerDelta(0,gO2)*ZP(gI1,1)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO2)*(10*ZH(gI2,2)*(2*vS*(AbsSqr(Lambdax)
      + QHd*Qs*Sqr(gp))*ZP(gI1,0) + 1.4142135623730951*Conj(TLambdax)*ZP(gI1,1)) +
      ZH(gI2,1)*(vu*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI1,0) + 5*
      vd*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) + ZH(gI2,0)*(vd*(3*Sqr(g1) + 5*
      Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI1,0) + 5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2
      ))*ZP(gI1,1)))) - KroneckerDelta(1,gO2)*(ZH(gI2,0)*(5*vu*(-2*AbsSqr(Lambdax)
      + Sqr(g2))*ZP(gI1,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(
      gI1,1)) + 10*ZH(gI2,2)*(1.4142135623730951*TLambdax*ZP(gI1,0) + 2*vS*(AbsSqr
      (Lambdax) + QHu*Qs*Sqr(gp))*ZP(gI1,1)) + ZH(gI2,1)*(5*vd*(-2*AbsSqr(Lambdax)
      + Sqr(g2))*ZP(gI1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZP(
      gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))) - Conj(ZA(gI1
      ,0))*(-5*Conj(ZA(gI2,1))*(KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) +
      KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*AbsSqr(Lambdax) + Sqr(g2))
      + Conj(ZA(gI2,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-3*Sqr(g1) +
      5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,
      gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd)))) + Conj(ZA(gI1,1))*(5*
      Conj(ZA(gI2,0))*(KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) +
      KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*AbsSqr(Lambdax) + Sqr(g2))
      + Conj(ZA(gI2,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) -
      5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp))) - KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   std::complex<double> tmp_1968;
   tmp_1968 += std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1968 += std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1968 += std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(gp);
   tmp_1968 += std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1968 += std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1968 += std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(gp);
   std::complex<double> tmp_1969;
   std::complex<double> tmp_1970;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1971;
      std::complex<double> tmp_1972;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1973;
         std::complex<double> tmp_1974;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1974 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1973 += tmp_1974;
         tmp_1972 += (Conj(ZV(gI2,j2))) * tmp_1973;
      }
      tmp_1971 += tmp_1972;
      tmp_1970 += (ZV(gI1,j3)) * tmp_1971;
   }
   tmp_1969 += tmp_1970;
   tmp_1968 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1969;
   result += (std::complex<double>(0,-1)) * tmp_1968;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(5*KroneckerDelta(1,gO2)*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) +
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))*ZH(gI1
      ,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZH(gI1,1)*ZH(
      gI2,1) + 20*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))*ZH(gI1,2)*ZH(gI2,2)))) +
      KroneckerDelta(1,gO1)*(-5*KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2
      ))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(1,gO2)*((3*
      Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp)))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1)
      + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))*ZH(gI1,1)*ZH(gI2,1) - 20*(AbsSqr(Lambdax)
      + QHu*Qs*Sqr(gp))*ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1975;
   std::complex<double> tmp_1976;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1977;
      std::complex<double> tmp_1978;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1978 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1977 += tmp_1978;
      tmp_1976 += (ZUL(gI1,j2)) * tmp_1977;
   }
   tmp_1975 += tmp_1976;
   result += (KroneckerDelta(0,gO2)) * tmp_1975;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1979;
   std::complex<double> tmp_1980;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1981;
      std::complex<double> tmp_1982;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1982 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1981 += tmp_1982;
      tmp_1980 += (Conj(ZDL(gI2,j2))) * tmp_1981;
   }
   tmp_1979 += tmp_1980;
   result += (KroneckerDelta(1,gO1)) * tmp_1979;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1983;
   std::complex<double> tmp_1984;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1984 += Conj(Ye(j1,gI1))*ZER(gI2,j1);
   }
   tmp_1983 += tmp_1984;
   result += (KroneckerDelta(0,gO2)) * tmp_1983;

   return result;
}

double CLASSNAME::CpconjUHpmbarFvFePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1985;
   std::complex<double> tmp_1986;
   std::complex<double> tmp_1987;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1987 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1986 += tmp_1987;
   tmp_1985 += (std::complex<double>(0,-0.35355339059327373)*vd*KroneckerDelta(
      0,gO2)*Sqr(g2)) * tmp_1986;
   std::complex<double> tmp_1988;
   std::complex<double> tmp_1989;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1989 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1988 += tmp_1989;
   tmp_1985 += (std::complex<double>(0,-0.35355339059327373)*vu*KroneckerDelta(
      1,gO2)*Sqr(g2)) * tmp_1988;
   std::complex<double> tmp_1990;
   std::complex<double> tmp_1991;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1992;
      std::complex<double> tmp_1993;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1993 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1992 += tmp_1993;
      tmp_1991 += (ZV(gI1,j2)) * tmp_1992;
   }
   tmp_1990 += tmp_1991;
   tmp_1985 += (std::complex<double>(0,0.7071067811865475)*vS*KroneckerDelta(1,
      gO2)*Lambdax) * tmp_1990;
   std::complex<double> tmp_1994;
   std::complex<double> tmp_1995;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1996;
      std::complex<double> tmp_1997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1997 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1996 += tmp_1997;
      tmp_1995 += (ZV(gI1,j2)) * tmp_1996;
   }
   tmp_1994 += tmp_1995;
   tmp_1985 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1994;
   std::complex<double> tmp_1998;
   std::complex<double> tmp_1999;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2000;
      std::complex<double> tmp_2001;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2002;
         std::complex<double> tmp_2003;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2003 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_2002 += tmp_2003;
         tmp_2001 += (Conj(ZE(gI2,j2))) * tmp_2002;
      }
      tmp_2000 += tmp_2001;
      tmp_1999 += (ZV(gI1,j3)) * tmp_2000;
   }
   tmp_1998 += tmp_1999;
   tmp_1985 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(0,
      gO2)) * tmp_1998;
   result += (std::complex<double>(0,-1)) * tmp_1985;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2004;
   std::complex<double> tmp_2005;
   std::complex<double> tmp_2006;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2006 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2005 += tmp_2006;
   tmp_2004 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2005;
   std::complex<double> tmp_2007;
   std::complex<double> tmp_2008;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2008 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2007 += tmp_2008;
   tmp_2004 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2007;
   std::complex<double> tmp_2009;
   std::complex<double> tmp_2010;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2010 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2009 += tmp_2010;
   tmp_2004 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2009;
   std::complex<double> tmp_2011;
   std::complex<double> tmp_2012;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2012 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2011 += tmp_2012;
   tmp_2004 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2011;
   std::complex<double> tmp_2013;
   std::complex<double> tmp_2014;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2014 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2013 += tmp_2014;
   tmp_2004 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2013;
   std::complex<double> tmp_2015;
   std::complex<double> tmp_2016;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2016 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2015 += tmp_2016;
   tmp_2004 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2015;
   std::complex<double> tmp_2017;
   std::complex<double> tmp_2018;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2018 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2017 += tmp_2018;
   tmp_2004 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2017;
   std::complex<double> tmp_2019;
   std::complex<double> tmp_2020;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2020 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2019 += tmp_2020;
   tmp_2004 += (std::complex<double>(0,-1)*Qd*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2019;
   std::complex<double> tmp_2021;
   std::complex<double> tmp_2022;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2022 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2021 += tmp_2022;
   tmp_2004 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2021;
   std::complex<double> tmp_2023;
   std::complex<double> tmp_2024;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2024 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2023 += tmp_2024;
   tmp_2004 += (std::complex<double>(0,-1)*Qd*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2023;
   std::complex<double> tmp_2025;
   std::complex<double> tmp_2026;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2027;
      std::complex<double> tmp_2028;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2029;
         std::complex<double> tmp_2030;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2030 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_2029 += tmp_2030;
         tmp_2028 += (ZD(gI1,3 + j2)) * tmp_2029;
      }
      tmp_2027 += tmp_2028;
      tmp_2026 += (Conj(ZD(gI2,3 + j3))) * tmp_2027;
   }
   tmp_2025 += tmp_2026;
   tmp_2004 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2025;
   std::complex<double> tmp_2031;
   std::complex<double> tmp_2032;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2033;
      std::complex<double> tmp_2034;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2035;
         std::complex<double> tmp_2036;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2036 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2035 += tmp_2036;
         tmp_2034 += (Conj(ZD(gI2,j2))) * tmp_2035;
      }
      tmp_2033 += tmp_2034;
      tmp_2032 += (ZD(gI1,j3)) * tmp_2033;
   }
   tmp_2031 += tmp_2032;
   tmp_2004 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2031;
   result += (std::complex<double>(0,-1)) * tmp_2004;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2037;
   std::complex<double> tmp_2038;
   std::complex<double> tmp_2039;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2039 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2038 += tmp_2039;
   tmp_2037 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2038;
   std::complex<double> tmp_2040;
   std::complex<double> tmp_2041;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2041 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2040 += tmp_2041;
   tmp_2037 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2040;
   std::complex<double> tmp_2042;
   std::complex<double> tmp_2043;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2043 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2042 += tmp_2043;
   tmp_2037 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2042;
   std::complex<double> tmp_2044;
   std::complex<double> tmp_2045;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2045 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2044 += tmp_2045;
   tmp_2037 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2044;
   std::complex<double> tmp_2046;
   std::complex<double> tmp_2047;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2047 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2046 += tmp_2047;
   tmp_2037 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2046;
   std::complex<double> tmp_2048;
   std::complex<double> tmp_2049;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2049 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2048 += tmp_2049;
   tmp_2037 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2048;
   std::complex<double> tmp_2050;
   std::complex<double> tmp_2051;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2051 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2050 += tmp_2051;
   tmp_2037 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2050;
   std::complex<double> tmp_2052;
   std::complex<double> tmp_2053;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2053 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2052 += tmp_2053;
   tmp_2037 += (std::complex<double>(0,-1)*Qe*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2052;
   std::complex<double> tmp_2054;
   std::complex<double> tmp_2055;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2055 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2054 += tmp_2055;
   tmp_2037 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2054;
   std::complex<double> tmp_2056;
   std::complex<double> tmp_2057;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2057 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2056 += tmp_2057;
   tmp_2037 += (std::complex<double>(0,-1)*Qe*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2056;
   std::complex<double> tmp_2058;
   std::complex<double> tmp_2059;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2060;
      std::complex<double> tmp_2061;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2062;
         std::complex<double> tmp_2063;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2063 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_2062 += tmp_2063;
         tmp_2061 += (ZE(gI1,3 + j2)) * tmp_2062;
      }
      tmp_2060 += tmp_2061;
      tmp_2059 += (Conj(ZE(gI2,3 + j3))) * tmp_2060;
   }
   tmp_2058 += tmp_2059;
   tmp_2037 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2058;
   result += (std::complex<double>(0,-1)) * tmp_2037;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2064;
   std::complex<double> tmp_2065;
   std::complex<double> tmp_2066;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2066 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2065 += tmp_2066;
   tmp_2064 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2065;
   std::complex<double> tmp_2067;
   std::complex<double> tmp_2068;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2068 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2067 += tmp_2068;
   tmp_2064 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2067;
   std::complex<double> tmp_2069;
   std::complex<double> tmp_2070;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2070 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2069 += tmp_2070;
   tmp_2064 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2069;
   std::complex<double> tmp_2071;
   std::complex<double> tmp_2072;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2072 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2071 += tmp_2072;
   tmp_2064 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2071;
   std::complex<double> tmp_2073;
   std::complex<double> tmp_2074;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2074 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2073 += tmp_2074;
   tmp_2064 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2073;
   std::complex<double> tmp_2075;
   std::complex<double> tmp_2076;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2076 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2075 += tmp_2076;
   tmp_2064 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2075;
   std::complex<double> tmp_2077;
   std::complex<double> tmp_2078;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2078 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2077 += tmp_2078;
   tmp_2064 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2077;
   std::complex<double> tmp_2079;
   std::complex<double> tmp_2080;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2080 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2079 += tmp_2080;
   tmp_2064 += (std::complex<double>(0,-1)*QHd*Qu*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2079;
   std::complex<double> tmp_2081;
   std::complex<double> tmp_2082;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2082 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2081 += tmp_2082;
   tmp_2064 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2081;
   std::complex<double> tmp_2083;
   std::complex<double> tmp_2084;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2084 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2083 += tmp_2084;
   tmp_2064 += (std::complex<double>(0,-1)*QHu*Qu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2083;
   std::complex<double> tmp_2085;
   std::complex<double> tmp_2086;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2087;
      std::complex<double> tmp_2088;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2089;
         std::complex<double> tmp_2090;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2090 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_2089 += tmp_2090;
         tmp_2088 += (ZU(gI1,3 + j2)) * tmp_2089;
      }
      tmp_2087 += tmp_2088;
      tmp_2086 += (Conj(ZU(gI2,3 + j3))) * tmp_2087;
   }
   tmp_2085 += tmp_2086;
   tmp_2064 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2085;
   std::complex<double> tmp_2091;
   std::complex<double> tmp_2092;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2093;
      std::complex<double> tmp_2094;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2095;
         std::complex<double> tmp_2096;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2096 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2095 += tmp_2096;
         tmp_2094 += (Conj(ZU(gI2,j2))) * tmp_2095;
      }
      tmp_2093 += tmp_2094;
      tmp_2092 += (ZU(gI1,j3)) * tmp_2093;
   }
   tmp_2091 += tmp_2092;
   tmp_2064 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2091;
   result += (std::complex<double>(0,-1)) * tmp_2064;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -0.1*KroneckerDelta(1,gO2)*(1.4142135623730951*UP(gI2,1)*(10*gp*QHu
      *ZN(gI1,0) + 3.872983346207417*g1*ZN(gI1,1) + 5*g2*ZN(gI1,2)) + 10*g2*UP(gI2
      ,0)*ZN(gI1,4)) - Conj(Lambdax)*KroneckerDelta(0,gO2)*UP(gI2,1)*ZN(gI1,5);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = -(g2*Conj(UM(gI2,0))*Conj(ZN(gI1,3))*KroneckerDelta(0,gO1)) + Conj(
      UM(gI2,1))*(-1.4142135623730951*gp*QHd*Conj(ZN(gI1,0))*KroneckerDelta(0,gO1)
      + 0.5477225575051661*g1*Conj(ZN(gI1,1))*KroneckerDelta(0,gO1) +
      0.7071067811865475*g2*Conj(ZN(gI1,2))*KroneckerDelta(0,gO1) - Conj(ZN(gI1,5)
      )*KroneckerDelta(1,gO1)*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2097;
   std::complex<double> tmp_2098;
   std::complex<double> tmp_2099;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2099 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2098 += tmp_2099;
   tmp_2097 += (std::complex<double>(0,-0.35355339059327373)*vd*KroneckerDelta(
      0,gO2)*Sqr(g2)) * tmp_2098;
   std::complex<double> tmp_2100;
   std::complex<double> tmp_2101;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2101 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2100 += tmp_2101;
   tmp_2097 += (std::complex<double>(0,-0.35355339059327373)*vu*KroneckerDelta(
      1,gO2)*Sqr(g2)) * tmp_2100;
   std::complex<double> tmp_2102;
   std::complex<double> tmp_2103;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2104;
      std::complex<double> tmp_2105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2105 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2104 += tmp_2105;
      tmp_2103 += (Conj(ZD(gI2,j2))) * tmp_2104;
   }
   tmp_2102 += tmp_2103;
   tmp_2097 += (std::complex<double>(0,0.7071067811865475)*vS*Conj(Lambdax)*
      KroneckerDelta(0,gO2)) * tmp_2102;
   std::complex<double> tmp_2106;
   std::complex<double> tmp_2107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2108;
      std::complex<double> tmp_2109;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2109 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_2108 += tmp_2109;
      tmp_2107 += (Conj(ZD(gI2,j2))) * tmp_2108;
   }
   tmp_2106 += tmp_2107;
   tmp_2097 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_2106;
   std::complex<double> tmp_2110;
   std::complex<double> tmp_2111;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2112;
      std::complex<double> tmp_2113;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2113 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2112 += tmp_2113;
      tmp_2111 += (ZU(gI1,j2)) * tmp_2112;
   }
   tmp_2110 += tmp_2111;
   tmp_2097 += (std::complex<double>(0,0.7071067811865475)*vS*KroneckerDelta(1,
      gO2)*Lambdax) * tmp_2110;
   std::complex<double> tmp_2114;
   std::complex<double> tmp_2115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2116;
      std::complex<double> tmp_2117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2117 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_2116 += tmp_2117;
      tmp_2115 += (ZU(gI1,j2)) * tmp_2116;
   }
   tmp_2114 += tmp_2115;
   tmp_2097 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_2114;
   std::complex<double> tmp_2118;
   std::complex<double> tmp_2119;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2120;
      std::complex<double> tmp_2121;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2122;
         std::complex<double> tmp_2123;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2123 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_2122 += tmp_2123;
         tmp_2121 += (ZU(gI1,3 + j2)) * tmp_2122;
      }
      tmp_2120 += tmp_2121;
      tmp_2119 += (Conj(ZD(gI2,3 + j3))) * tmp_2120;
   }
   tmp_2118 += tmp_2119;
   tmp_2097 += (std::complex<double>(0,0.7071067811865475)*vu*KroneckerDelta(0,
      gO2)) * tmp_2118;
   std::complex<double> tmp_2124;
   std::complex<double> tmp_2125;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2126;
      std::complex<double> tmp_2127;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2128;
         std::complex<double> tmp_2129;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2129 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_2128 += tmp_2129;
         tmp_2127 += (ZU(gI1,3 + j2)) * tmp_2128;
      }
      tmp_2126 += tmp_2127;
      tmp_2125 += (Conj(ZD(gI2,3 + j3))) * tmp_2126;
   }
   tmp_2124 += tmp_2125;
   tmp_2097 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(1,
      gO2)) * tmp_2124;
   std::complex<double> tmp_2130;
   std::complex<double> tmp_2131;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2132;
      std::complex<double> tmp_2133;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2134;
         std::complex<double> tmp_2135;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2135 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2134 += tmp_2135;
         tmp_2133 += (Conj(ZD(gI2,j2))) * tmp_2134;
      }
      tmp_2132 += tmp_2133;
      tmp_2131 += (ZU(gI1,j3)) * tmp_2132;
   }
   tmp_2130 += tmp_2131;
   tmp_2097 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(0,
      gO2)) * tmp_2130;
   std::complex<double> tmp_2136;
   std::complex<double> tmp_2137;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2138;
      std::complex<double> tmp_2139;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2140;
         std::complex<double> tmp_2141;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2141 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2140 += tmp_2141;
         tmp_2139 += (Conj(ZD(gI2,j2))) * tmp_2140;
      }
      tmp_2138 += tmp_2139;
      tmp_2137 += (ZU(gI1,j3)) * tmp_2138;
   }
   tmp_2136 += tmp_2137;
   tmp_2097 += (std::complex<double>(0,0.7071067811865475)*vu*KroneckerDelta(1,
      gO2)) * tmp_2136;
   result += (std::complex<double>(0,-1)) * tmp_2097;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVPHpm(unsigned gO2, unsigned gI2) const
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

std::complex<double> CLASSNAME::CpconjUHpmVZHpm(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(KroneckerDelta(0,gO2)*(-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())
      )*ZP(gI2,0) + KroneckerDelta(1,gO2)*(-5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      )*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZpHpm(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(KroneckerDelta(0,gO2)*(10*gp*QHd*Cos(ThetaWp()) + (5*g2*Cos(
      ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZP(gI2,0) -
      KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gI2,0))*KroneckerDelta(0,
      gO2) + Conj(ZA(gI2,1))*KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) - KroneckerDelta(1,gO2)*ZH(
      gI2,1));

   return result;
}

double CLASSNAME::CpVZbargWmgWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW())*Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZbargWmCgWmC() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW())*Cos(ThetaWp());

   return result;
}

double CLASSNAME::CpVZconjVWmVWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW())*Cos(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*((15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin
      (ThetaWp()) - 2*g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(
      ThetaW())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp())))*ZP(gI1,0)*ZP(gI2,0) +
      (-15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) -
      2*g2*Cos(ThetaW())*Cos(ThetaWp())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(
      ThetaW()) - 10*gp*QHu*Sin(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos
      (ThetaWp())) + 3*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW())) + 20*Sqr(gp)
      *Sqr(QHu)*Sqr(Sin(ThetaWp())))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*((-5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos
      (ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp()))*ZP(gI1,0)*ZP(gI2,0) +
      (-5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp()))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM(gI2,0))*Cos(ThetaW())*Cos(ThetaWp())*UM(gI1,0) +
      Conj(UM(gI2,1))*(g2*Cos(ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos
      (ThetaWp())*Sin(ThetaW()) - 2*gp*QHd*Sin(ThetaWp()))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP(gI1,0))*Cos(ThetaW())*Cos(ThetaWp())*UP(gI2,0) +
      Conj(UP(gI1,1))*(g2*Cos(ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos
      (ThetaWp())*Sin(ThetaW()) + 2*gp*QHu*Sin(ThetaWp()))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZAhAh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(
      ThetaWp())) + Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(20*g2*gp*QHd*Cos(ThetaW())*
      Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) +
      Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(-20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*
      Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1
      *Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   result = 0.1*KroneckerDelta(gI1,gI2)*(20*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp(
      ))*Sin(ThetaWp()) + 15.491933384829668*g1*gp*Ql*Cos(ThetaWp())*Sin(ThetaW())
      *Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*
      g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos
      (ThetaWp())) + 20*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhhhh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*((20*g2*gp*QHd*Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) +
      15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) +
      g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*
      Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*
      Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp())))*ZH(gI1,0)*ZH(gI2,0) + (-20*g2*gp*QHu*
      Cos(ThetaW())*Cos(ThetaWp())*Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*
      Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp()))
      + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Sin(ThetaWp())))*ZH(gI1,1)*ZH(gI2,1) + 20*Sqr(gp)*Sqr(Qs)*Sqr(Sin(
      ThetaWp()))*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*Sin(ThetaWp()))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(Conj(ZA(gI2,0))*(g2*Cos(ThetaW())*Cos
      (ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*
      Sin(ThetaWp()))*ZH(gI1,0) - Conj(ZA(gI2,1))*(g2*Cos(ThetaW())*Cos(ThetaWp())
      + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*Sin(ThetaWp
      ()))*ZH(gI1,1) + 2*gp*Qs*Conj(ZA(gI2,2))*Sin(ThetaWp())*ZH(gI1,2));

   return result;
}

double CLASSNAME::CpVZbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(ThetaW())*Cos
      (ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 6*gp*Qq*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   double result = 0.0;

   result = KroneckerDelta(gI1,gI2)*(-0.2581988897471611*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + gp*Qd*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZbarFeFePL(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW())*Cos(ThetaWp()) -
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*Ql*Sin(ThetaWp()))
      ;

   return result;
}

double CLASSNAME::CpVZbarFeFePR(unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   double result = 0.0;

   result = -(KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaWp())*Sin
      (ThetaW()) - gp*Qe*Sin(ThetaWp())));

   return result;
}

double CLASSNAME::CpVZbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(ThetaW())*
      Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 6*gp*
      Qq*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   double result = 0.0;

   result = KroneckerDelta(gI1,gI2)*(0.5163977794943222*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + gp*Qu*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZbarFvFvPL(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*Sin(ThetaWp()))
      ;

   return result;
}

double CLASSNAME::CpVZbarFvFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2142;
   std::complex<double> tmp_2143;
   std::complex<double> tmp_2144;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2144 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2143 += tmp_2144;
   tmp_2142 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_2143;
   std::complex<double> tmp_2145;
   std::complex<double> tmp_2146;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2146 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2145 += tmp_2146;
   tmp_2142 += (std::complex<double>(0,0.2581988897471611)*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_2145;
   std::complex<double> tmp_2147;
   std::complex<double> tmp_2148;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2148 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2147 += tmp_2148;
   tmp_2142 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_2147;
   std::complex<double> tmp_2149;
   std::complex<double> tmp_2150;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2150 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2149 += tmp_2150;
   tmp_2142 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_2149;
   std::complex<double> tmp_2151;
   std::complex<double> tmp_2152;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2152 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2151 += tmp_2152;
   tmp_2142 += (std::complex<double>(0,-2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_2151;
   std::complex<double> tmp_2153;
   std::complex<double> tmp_2154;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2154 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2153 += tmp_2154;
   tmp_2142 += (std::complex<double>(0,-0.5163977794943222)*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2153;
   std::complex<double> tmp_2155;
   std::complex<double> tmp_2156;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2156 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2155 += tmp_2156;
   tmp_2142 += (std::complex<double>(0,-1.0327955589886444)*g1*gp*Qd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2155;
   std::complex<double> tmp_2157;
   std::complex<double> tmp_2158;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2158 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2157 += tmp_2158;
   tmp_2142 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())))
      * tmp_2157;
   std::complex<double> tmp_2159;
   std::complex<double> tmp_2160;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2160 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2159 += tmp_2160;
   tmp_2142 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qd)*Sqr(Sin(ThetaWp())))
      * tmp_2159;
   result += (std::complex<double>(0,-1)) * tmp_2142;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2161;
   std::complex<double> tmp_2162;
   std::complex<double> tmp_2163;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2163 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2162 += tmp_2163;
   tmp_2161 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_2162;
   std::complex<double> tmp_2164;
   std::complex<double> tmp_2165;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2165 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2164 += tmp_2165;
   tmp_2161 += (std::complex<double>(0,-0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_2164;
   std::complex<double> tmp_2166;
   std::complex<double> tmp_2167;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2167 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2166 += tmp_2167;
   tmp_2161 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin
      (ThetaW()))) * tmp_2166;
   std::complex<double> tmp_2168;
   std::complex<double> tmp_2169;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2169 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2168 += tmp_2169;
   tmp_2161 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin
      (ThetaW()))) * tmp_2168;
   std::complex<double> tmp_2170;
   std::complex<double> tmp_2171;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2171 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2170 += tmp_2171;
   tmp_2161 += (std::complex<double>(0,-2)*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_2170;
   std::complex<double> tmp_2172;
   std::complex<double> tmp_2173;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2173 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2172 += tmp_2173;
   tmp_2161 += (std::complex<double>(0,1.5491933384829668)*g1*gp*Ql*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2172;
   std::complex<double> tmp_2174;
   std::complex<double> tmp_2175;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2175 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2174 += tmp_2175;
   tmp_2161 += (std::complex<double>(0,-3.0983866769659336)*g1*gp*Qe*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2174;
   std::complex<double> tmp_2176;
   std::complex<double> tmp_2177;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2177 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2176 += tmp_2177;
   tmp_2161 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())))
      * tmp_2176;
   std::complex<double> tmp_2178;
   std::complex<double> tmp_2179;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2179 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2178 += tmp_2179;
   tmp_2161 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qe)*Sqr(Sin(ThetaWp())))
      * tmp_2178;
   result += (std::complex<double>(0,-1)) * tmp_2161;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2180;
   std::complex<double> tmp_2181;
   std::complex<double> tmp_2182;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2182 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2181 += tmp_2182;
   tmp_2180 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_2181;
   std::complex<double> tmp_2183;
   std::complex<double> tmp_2184;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2184 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2183 += tmp_2184;
   tmp_2180 += (std::complex<double>(0,-0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_2183;
   std::complex<double> tmp_2185;
   std::complex<double> tmp_2186;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2186 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2185 += tmp_2186;
   tmp_2180 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_2185;
   std::complex<double> tmp_2187;
   std::complex<double> tmp_2188;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2188 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2187 += tmp_2188;
   tmp_2180 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_2187;
   std::complex<double> tmp_2189;
   std::complex<double> tmp_2190;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2190 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2189 += tmp_2190;
   tmp_2180 += (std::complex<double>(0,2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_2189;
   std::complex<double> tmp_2191;
   std::complex<double> tmp_2192;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2192 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2191 += tmp_2192;
   tmp_2180 += (std::complex<double>(0,-0.5163977794943222)*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2191;
   std::complex<double> tmp_2193;
   std::complex<double> tmp_2194;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2194 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2193 += tmp_2194;
   tmp_2180 += (std::complex<double>(0,2.065591117977289)*g1*gp*Qu*Cos(ThetaWp(
      ))*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2193;
   std::complex<double> tmp_2195;
   std::complex<double> tmp_2196;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2196 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2195 += tmp_2196;
   tmp_2180 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())))
      * tmp_2195;
   std::complex<double> tmp_2197;
   std::complex<double> tmp_2198;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2198 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2197 += tmp_2198;
   tmp_2180 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qu)*Sqr(Sin(ThetaWp())))
      * tmp_2197;
   result += (std::complex<double>(0,-1)) * tmp_2180;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2199;
   std::complex<double> tmp_2200;
   std::complex<double> tmp_2201;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2201 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2200 += tmp_2201;
   tmp_2199 += (2*(0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 3*gp*Qd
      *Sin(ThetaWp()))) * tmp_2200;
   std::complex<double> tmp_2202;
   std::complex<double> tmp_2203;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2203 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2202 += tmp_2203;
   tmp_2199 += (-3*g2*Cos(ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 6*gp*Qq*Sin(ThetaWp())) * tmp_2202;
   result += (0.16666666666666666) * tmp_2199;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2204;
   std::complex<double> tmp_2205;
   std::complex<double> tmp_2206;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2206 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2205 += tmp_2206;
   tmp_2204 += (-2*(-0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + gp*Qe
      *Sin(ThetaWp()))) * tmp_2205;
   std::complex<double> tmp_2207;
   std::complex<double> tmp_2208;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2208 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2207 += tmp_2208;
   tmp_2204 += (-(g2*Cos(ThetaW())*Cos(ThetaWp())) + 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*Sin(ThetaWp())) * tmp_2207;
   result += (0.5) * tmp_2204;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2209;
   std::complex<double> tmp_2210;
   std::complex<double> tmp_2211;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2211 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2210 += tmp_2211;
   tmp_2209 += (3*g2*Cos(ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 6*gp*Qq*Sin(ThetaWp())) * tmp_2210;
   std::complex<double> tmp_2212;
   std::complex<double> tmp_2213;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2213 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2212 += tmp_2213;
   tmp_2209 += (-2*(1.5491933384829668*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3*gp*
      Qu*Sin(ThetaWp()))) * tmp_2212;
   result += (0.16666666666666666) * tmp_2209;

   return result;
}

std::complex<double> CLASSNAME::CpVZChiChiPL(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.5*(-(Conj(ZN(gI2,3))*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*Sin(ThetaWp())
      )*ZN(gI1,3)) + Conj(ZN(gI2,4))*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*Sin(ThetaWp())
      )*ZN(gI1,4) - 2*gp*Qs*Conj(ZN(gI2,5))*Sin(ThetaWp())*ZN(gI1,5));

   return result;
}

std::complex<double> CLASSNAME::CpVZChiChiPR(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.5*(Conj(ZN(gI1,3))*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*Sin(ThetaWp())
      )*ZN(gI2,3) - Conj(ZN(gI1,4))*(g2*Cos(ThetaW())*Cos(ThetaWp()) +
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*Sin(ThetaWp())
      )*ZN(gI2,4) + 2*gp*Qs*Conj(ZN(gI1,5))*Sin(ThetaWp())*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpVZconjVWmHpm(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*
      gp*QHd*Sin(ThetaWp()))*ZP(gI2,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp())
      *Sin(ThetaW()) + 2*gp*QHu*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhh(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.5*(vd*Sqr(g2*Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1
      *Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*Sin(ThetaWp()))*ZH(gI2,0) + vu*Sqr(
      g2*Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(
      ThetaW()) - 2*gp*QHu*Sin(ThetaWp()))*ZH(gI2,1) + 4*vS*Sqr(gp)*Sqr(Qs)*Sqr(
      Sin(ThetaWp()))*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZphh(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(-(vd*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW()
      )) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(
      ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW()
      ))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*
      Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp
      ()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(ThetaWp()))))*ZH(gI2,0
      )) - vu*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) +
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp(
      ))*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) -
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(
      ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())
      + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(Sin(ThetaWp()))))*ZH(gI2,1) +
      20*vS*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(gp)*Sqr(Qs)*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZVZconjVWmVWm3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbargWmgWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW())*Sin(ThetaWp());

   return result;
}

double CLASSNAME::CpVZpbargWmCgWmC() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpconjVWmVWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW())*Sin(ThetaWp());

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*((2*gp*QHd*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(
      ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (
      -7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(ThetaWp())))*ZP(gI1,0)*ZP
      (gI2,0) + (-2*gp*QHu*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()
      ))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (
      -7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*Sqr(Sin(ThetaWp())))*ZP(gI1,1)*ZP
      (gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*((2*gp*QHd*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) -
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI1,0)*ZP(gI2,0) + (
      -2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) -
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZpbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = 0.5*(-2*g2*Conj(UM(gI2,0))*Cos(ThetaW())*Sin(ThetaWp())*UM(gI1,0) -
      Conj(UM(gI2,1))*(2*gp*QHd*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp())
      - 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*UM(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZpbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(-10*g2*Conj(UP(gI1,0))*Cos(ThetaW())*Sin(ThetaWp())*UP(gI2,0)
      + Conj(UP(gI1,1))*(10*gp*QHu*Cos(ThetaWp()) + (-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpAhAh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(
      ThetaWp())) + Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*(-2*gp*QHd*(5*g2*Cos(ThetaW())
      + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd
      )*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()
      ) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp())))
      + Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(2*gp*QHu*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   result = 0.1*KroneckerDelta(gI1,gI2)*(-2*gp*Ql*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(Ql)*
      Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW())
      + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZphhhh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*((-2*gp*QHd*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(
      ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Cos(ThetaWp())) + (g1*
      Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*
      Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp())))*ZH(gI1,0)*ZH(gI2,0) + (2*gp
      *QHu*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp
      ()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(
      7.745966692414834*g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos
      (ThetaW())))*Sqr(Sin(ThetaWp())))*ZH(gI1,1)*ZH(gI2,1) + 20*Sqr(gp)*Sqr(Qs)*
      Sqr(Cos(ThetaWp()))*ZH(gI1,2)*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpVZpconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(2*gp*Ql*Cos(ThetaWp()) - (g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpVZphhAh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(Conj(ZA(gI2,0))*(2*gp*QHd*Cos(ThetaWp
      ()) - (g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp()
      ))*ZH(gI1,0) + Conj(ZA(gI2,1))*(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*
      Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZH(gI1,
      1) + 2*gp*Qs*Conj(ZA(gI2,2))*Cos(ThetaWp())*ZH(gI1,2));

   return result;
}

double CLASSNAME::CpVZpbarFdFdPL(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(6*gp*Qq*Cos(ThetaWp()
      ) + 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbarFdFdPR(unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   double result = 0.0;

   result = KroneckerDelta(gI1,gI2)*(gp*Qd*Cos(ThetaWp()) + 0.2581988897471611*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbarFeFePL(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(2*gp*Ql*Cos(ThetaWp()) + g2*Cos(
      ThetaW())*Sin(ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()
      ));

   return result;
}

double CLASSNAME::CpVZpbarFeFePR(unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   double result = 0.0;

   result = KroneckerDelta(gI1,gI2)*(gp*Qe*Cos(ThetaWp()) + 0.7745966692414834*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbarFuFuPL(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(6*gp*Qq*Cos(ThetaWp()
      ) - 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbarFuFuPR(unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   double result = 0.0;

   result = KroneckerDelta(gI1,gI2)*(gp*Qu*Cos(ThetaWp()) - 0.5163977794943222*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbarFvFvPL(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = -0.5*KroneckerDelta(gI1,gI2)*(2*gp*Ql*Cos(ThetaWp()) - (g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpbarFvFvPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2214;
   std::complex<double> tmp_2215;
   std::complex<double> tmp_2216;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2216 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2215 += tmp_2216;
   tmp_2214 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())))
      * tmp_2215;
   std::complex<double> tmp_2217;
   std::complex<double> tmp_2218;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2218 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2217 += tmp_2218;
   tmp_2214 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qd)*Sqr(Cos(ThetaWp())))
      * tmp_2217;
   std::complex<double> tmp_2219;
   std::complex<double> tmp_2220;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2220 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2219 += tmp_2220;
   tmp_2214 += (std::complex<double>(0,2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_2219;
   std::complex<double> tmp_2221;
   std::complex<double> tmp_2222;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2222 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2221 += tmp_2222;
   tmp_2214 += (std::complex<double>(0,0.5163977794943222)*g1*gp*Qq*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2221;
   std::complex<double> tmp_2223;
   std::complex<double> tmp_2224;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2224 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2223 += tmp_2224;
   tmp_2214 += (std::complex<double>(0,1.0327955589886444)*g1*gp*Qd*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2223;
   std::complex<double> tmp_2225;
   std::complex<double> tmp_2226;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2226 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2225 += tmp_2226;
   tmp_2214 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_2225;
   std::complex<double> tmp_2227;
   std::complex<double> tmp_2228;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2228 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2227 += tmp_2228;
   tmp_2214 += (std::complex<double>(0,0.2581988897471611)*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_2227;
   std::complex<double> tmp_2229;
   std::complex<double> tmp_2230;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2230 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2229 += tmp_2230;
   tmp_2214 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_2229;
   std::complex<double> tmp_2231;
   std::complex<double> tmp_2232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2232 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2231 += tmp_2232;
   tmp_2214 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_2231;
   result += (std::complex<double>(0,-1)) * tmp_2214;

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2233;
   std::complex<double> tmp_2234;
   std::complex<double> tmp_2235;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2235 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2234 += tmp_2235;
   tmp_2233 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Ql)*Sqr(Cos(ThetaWp())))
      * tmp_2234;
   std::complex<double> tmp_2236;
   std::complex<double> tmp_2237;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2237 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2236 += tmp_2237;
   tmp_2233 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qe)*Sqr(Cos(ThetaWp())))
      * tmp_2236;
   std::complex<double> tmp_2238;
   std::complex<double> tmp_2239;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2239 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2238 += tmp_2239;
   tmp_2233 += (std::complex<double>(0,2)*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_2238;
   std::complex<double> tmp_2240;
   std::complex<double> tmp_2241;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2241 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2240 += tmp_2241;
   tmp_2233 += (std::complex<double>(0,-1.5491933384829668)*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2240;
   std::complex<double> tmp_2242;
   std::complex<double> tmp_2243;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2243 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2242 += tmp_2243;
   tmp_2233 += (std::complex<double>(0,3.0983866769659336)*g1*gp*Qe*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2242;
   std::complex<double> tmp_2244;
   std::complex<double> tmp_2245;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2245 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2244 += tmp_2245;
   tmp_2233 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_2244;
   std::complex<double> tmp_2246;
   std::complex<double> tmp_2247;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2247 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2246 += tmp_2247;
   tmp_2233 += (std::complex<double>(0,-0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_2246;
   std::complex<double> tmp_2248;
   std::complex<double> tmp_2249;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2249 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2248 += tmp_2249;
   tmp_2233 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_2248;
   std::complex<double> tmp_2250;
   std::complex<double> tmp_2251;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2251 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2250 += tmp_2251;
   tmp_2233 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_2250;
   result += (std::complex<double>(0,-1)) * tmp_2233;

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2252;
   std::complex<double> tmp_2253;
   std::complex<double> tmp_2254;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2254 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2253 += tmp_2254;
   tmp_2252 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())))
      * tmp_2253;
   std::complex<double> tmp_2255;
   std::complex<double> tmp_2256;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2256 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2255 += tmp_2256;
   tmp_2252 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qu)*Sqr(Cos(ThetaWp())))
      * tmp_2255;
   std::complex<double> tmp_2257;
   std::complex<double> tmp_2258;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2258 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2257 += tmp_2258;
   tmp_2252 += (std::complex<double>(0,-2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_2257;
   std::complex<double> tmp_2259;
   std::complex<double> tmp_2260;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2260 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2259 += tmp_2260;
   tmp_2252 += (std::complex<double>(0,0.5163977794943222)*g1*gp*Qq*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2259;
   std::complex<double> tmp_2261;
   std::complex<double> tmp_2262;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2262 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2261 += tmp_2262;
   tmp_2252 += (std::complex<double>(0,-2.065591117977289)*g1*gp*Qu*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_2261;
   std::complex<double> tmp_2263;
   std::complex<double> tmp_2264;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2264 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2263 += tmp_2264;
   tmp_2252 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_2263;
   std::complex<double> tmp_2265;
   std::complex<double> tmp_2266;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2266 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2265 += tmp_2266;
   tmp_2252 += (std::complex<double>(0,-0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_2265;
   std::complex<double> tmp_2267;
   std::complex<double> tmp_2268;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2268 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2267 += tmp_2268;
   tmp_2252 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_2267;
   std::complex<double> tmp_2269;
   std::complex<double> tmp_2270;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2270 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2269 += tmp_2270;
   tmp_2252 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_2269;
   result += (std::complex<double>(0,-1)) * tmp_2252;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2271;
   std::complex<double> tmp_2272;
   std::complex<double> tmp_2273;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2273 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2272 += tmp_2273;
   tmp_2271 += (-2*(3*gp*Qd*Cos(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW()
      )*Sin(ThetaWp()))) * tmp_2272;
   std::complex<double> tmp_2274;
   std::complex<double> tmp_2275;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2275 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2274 += tmp_2275;
   tmp_2271 += (6*gp*Qq*Cos(ThetaWp()) + (3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_2274;
   result += (0.16666666666666666) * tmp_2271;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2276;
   std::complex<double> tmp_2277;
   std::complex<double> tmp_2278;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2278 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2277 += tmp_2278;
   tmp_2276 += (-2*(gp*Qe*Cos(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()))) * tmp_2277;
   std::complex<double> tmp_2279;
   std::complex<double> tmp_2280;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2280 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2279 += tmp_2280;
   tmp_2276 += (2*gp*Ql*Cos(ThetaWp()) + (g2*Cos(ThetaW()) - 0.7745966692414834
      *g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_2279;
   result += (0.5) * tmp_2276;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2281;
   std::complex<double> tmp_2282;
   std::complex<double> tmp_2283;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2283 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2282 += tmp_2283;
   tmp_2281 += (2*(-3*gp*Qu*Cos(ThetaWp()) + 1.5491933384829668*g1*Sin(ThetaW()
      )*Sin(ThetaWp()))) * tmp_2282;
   std::complex<double> tmp_2284;
   std::complex<double> tmp_2285;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2285 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2284 += tmp_2285;
   tmp_2281 += (6*gp*Qq*Cos(ThetaWp()) + (-3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_2284;
   result += (0.16666666666666666) * tmp_2281;

   return result;
}

std::complex<double> CLASSNAME::CpVZpChiChiPL(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.5*(-(Conj(ZN(gI2,3))*(2*gp*QHd*Cos(ThetaWp()) - (g2*Cos(ThetaW())
      + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZN(gI1,3)) - Conj(ZN
      (gI2,4))*(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZN(gI1,4) - 2*gp*Qs*Conj
      (ZN(gI2,5))*Cos(ThetaWp())*ZN(gI1,5));

   return result;
}

std::complex<double> CLASSNAME::CpVZpChiChiPR(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.5*(Conj(ZN(gI1,3))*(2*gp*QHd*Cos(ThetaWp()) - (g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZN(gI2,3) + Conj(ZN(
      gI1,4))*(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZN(gI2,4) + 2*gp*Qs*Conj
      (ZN(gI1,5))*Cos(ThetaWp())*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjVWmHpm(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp()) - 0.7745966692414834*g1*Sin(
      ThetaW())*Sin(ThetaWp()))*ZP(gI2,0) + vu*(2*gp*QHu*Cos(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZhh(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(-(vd*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW()
      )) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(
      ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin(ThetaW()
      ))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*
      Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp
      ()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(ThetaWp()))))*ZH(gI2,0
      )) - vu*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW())) +
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(ThetaWp(
      ))*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))) -
      7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*Cos(
      ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())
      + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(Sin(ThetaWp()))))*ZH(gI2,1) +
      20*vS*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(gp)*Sqr(Qs)*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZphh(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.5*(vd*Sqr(-2*gp*QHd*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp
      ()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZH(gI2,0) + vu*Sqr
      (2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZH(gI2,1) + 4*vS*Sqr(gp)
      *Sqr(Qs)*Sqr(Cos(ThetaWp()))*ZH(gI2,2));

   return result;
}

double CLASSNAME::CpVZpVZpconjVWmVWm1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpVZpconjVWmVWm2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVZpVZpconjVWmVWm3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

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

   result = -(g2*Cos(ThetaW())*Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpconjVWmbargWmCgZp() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW())*Sin(ThetaWp());

   return result;
}

double CLASSNAME::CpconjVWmbargZgWm() const
{
   double result = 0.0;

   result = g2*Cos(ThetaW())*Cos(ThetaWp());

   return result;
}

double CLASSNAME::CpconjVWmbargZpgWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW())*Sin(ThetaWp()));

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

   result = g2*Cos(ThetaW())*Cos(ThetaWp());

   return result;
}

double CLASSNAME::CpconjVWmVZpVWm() const
{
   double result = 0.0;

   result = -(g2*Cos(ThetaW())*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(Conj(ZA(gI2,0))*ZP(gI1,0) + Conj(
      ZA(gI2,1))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)) + Conj(ZA(gI1,1))*Conj(ZA(gI2,
      1)))*Sqr(g2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2286;
   std::complex<double> tmp_2287;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2287 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_2286 += tmp_2287;
   result += (-0.7071067811865475*g2) * tmp_2286;

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
      result += -0.7071067811865475*g2*Conj(ZEL(gI2,gI1));
   }

   return result;
}

double CLASSNAME::CpconjVWmbarFvFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSvSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2288;
   std::complex<double> tmp_2289;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2289 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2288 += tmp_2289;
   result += (0.7071067811865475*g2) * tmp_2288;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2290;
   std::complex<double> tmp_2291;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2291 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2290 += tmp_2291;
   result += (0.5*Sqr(g2)) * tmp_2290;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2292;
   std::complex<double> tmp_2293;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2293 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2292 += tmp_2293;
   result += (0.5*Sqr(g2)) * tmp_2292;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2294;
   std::complex<double> tmp_2295;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2295 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2294 += tmp_2295;
   result += (0.5*Sqr(g2)) * tmp_2294;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,2) + 1.4142135623730951*Conj(UM(
      gI2,1))*ZN(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(ZN(gI1,2))*UP(gI2,0)) + 0.7071067811865475*g2*Conj(ZN(gI1
      ,4))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2296;
   std::complex<double> tmp_2297;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2297 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2296 += tmp_2297;
   result += (0.7071067811865475*g2) * tmp_2296;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVPHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.3872983346207417*g1*g2*Cos(ThetaW())*(vd*ZP(gI2,0) - vu*ZP(gI2,1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVZHpm(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*g2*(vd*(0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*
      gp*QHd*Sin(ThetaWp()))*ZP(gI2,0) + vu*(-0.7745966692414834*g1*Cos(ThetaWp())
      *Sin(ThetaW()) + 2*gp*QHu*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVZpHpm(unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*g2*(vd*(2*gp*QHd*Cos(ThetaWp()) - 0.7745966692414834*g1*Sin(
      ThetaW())*Sin(ThetaWp()))*ZP(gI2,0) + vu*(2*gp*QHu*Cos(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVWmhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1));

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

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZVZ3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(ThetaWp()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZpVZp1() const
{
   double result = 0.0;

   result = -2*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZpVZp2() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpVWmconjVWmVZpVZp3() const
{
   double result = 0.0;

   result = Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(ThetaWp()));

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

std::complex<double> CLASSNAME::CpUChiconjHpmChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = -(g2*Conj(UM(gI2,0))*KroneckerDelta(3,gO2)*ZP(gI1,0)) + Conj(UM(gI2
      ,1))*(-1.4142135623730951*gp*QHd*KroneckerDelta(0,gO2)*ZP(gI1,0) +
      0.5477225575051661*g1*KroneckerDelta(1,gO2)*ZP(gI1,0) + 0.7071067811865475*
      g2*KroneckerDelta(2,gO2)*ZP(gI1,0) - KroneckerDelta(5,gO2)*Lambdax*ZP(gI1,1)
      );

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHpmChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -(Conj(Lambdax)*KroneckerDelta(5,gO1)*UP(gI2,1)*ZP(gI1,0)) - 0.1*(
      10*g2*KroneckerDelta(4,gO1)*UP(gI2,0) + 1.4142135623730951*(10*gp*QHu*
      KroneckerDelta(0,gO1) + 3.872983346207417*g1*KroneckerDelta(1,gO1) + 5*g2*
      KroneckerDelta(2,gO1))*UP(gI2,1))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSvFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -1.4142135623730951*gp*Ql*KroneckerDelta(0,gO2)*ZV(gI1,gI2);
   }
   if (gI2 < 3) {
      result += 0.5477225575051661*g1*KroneckerDelta(1,gO2)*ZV(gI1,gI2);
   }
   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(2,gO2)*ZV(gI1,gI2);
   }

   return result;
}

double CLASSNAME::CpUChiconjSvFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(3.872983346207417*g1*Conj(ZN(gI2,1))*KroneckerDelta(3,gO2)*ZH(
      gI1,0) - 5*g2*Conj(ZN(gI2,2))*KroneckerDelta(3,gO2)*ZH(gI1,0) +
      7.0710678118654755*Conj(ZN(gI2,5))*KroneckerDelta(4,gO2)*Lambdax*ZH(gI1,0) +
      7.0710678118654755*Conj(ZN(gI2,4))*KroneckerDelta(5,gO2)*Lambdax*ZH(gI1,0)
      - 10*gp*QHu*Conj(ZN(gI2,4))*KroneckerDelta(0,gO2)*ZH(gI1,1) -
      3.872983346207417*g1*Conj(ZN(gI2,4))*KroneckerDelta(1,gO2)*ZH(gI1,1) + 5*g2*
      Conj(ZN(gI2,4))*KroneckerDelta(2,gO2)*ZH(gI1,1) - 3.872983346207417*g1*Conj(
      ZN(gI2,1))*KroneckerDelta(4,gO2)*ZH(gI1,1) + 5*g2*Conj(ZN(gI2,2))*
      KroneckerDelta(4,gO2)*ZH(gI1,1) + 7.0710678118654755*Conj(ZN(gI2,5))*
      KroneckerDelta(3,gO2)*Lambdax*ZH(gI1,1) - 10*gp*Qs*Conj(ZN(gI2,5))*
      KroneckerDelta(0,gO2)*ZH(gI1,2) + 7.0710678118654755*Conj(ZN(gI2,4))*
      KroneckerDelta(3,gO2)*Lambdax*ZH(gI1,2) - 10*gp*Conj(ZN(gI2,0))*(QHd*
      KroneckerDelta(3,gO2)*ZH(gI1,0) + QHu*KroneckerDelta(4,gO2)*ZH(gI1,1) + Qs*
      KroneckerDelta(5,gO2)*ZH(gI1,2)) + Conj(ZN(gI2,3))*(-10*gp*QHd*
      KroneckerDelta(0,gO2)*ZH(gI1,0) + 3.872983346207417*g1*KroneckerDelta(1,gO2)
      *ZH(gI1,0) - 5*g2*KroneckerDelta(2,gO2)*ZH(gI1,0) + 7.0710678118654755*
      KroneckerDelta(5,gO2)*Lambdax*ZH(gI1,1) + 7.0710678118654755*KroneckerDelta(
      4,gO2)*Lambdax*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(-10*gp*Qs*KroneckerDelta(5,gO1)*ZH(gI1,2)*ZN(gI2,0) - 10*gp*
      QHd*KroneckerDelta(0,gO1)*ZH(gI1,0)*ZN(gI2,3) + 3.872983346207417*g1*
      KroneckerDelta(1,gO1)*ZH(gI1,0)*ZN(gI2,3) - 5*g2*KroneckerDelta(2,gO1)*ZH(
      gI1,0)*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*KroneckerDelta(5,gO1)*ZH
      (gI1,1)*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*KroneckerDelta(5,gO1)*
      ZH(gI1,0)*ZN(gI2,4) - 10*gp*QHu*KroneckerDelta(0,gO1)*ZH(gI1,1)*ZN(gI2,4) -
      3.872983346207417*g1*KroneckerDelta(1,gO1)*ZH(gI1,1)*ZN(gI2,4) + 5*g2*
      KroneckerDelta(2,gO1)*ZH(gI1,1)*ZN(gI2,4) - 10*gp*Qs*KroneckerDelta(0,gO1)*
      ZH(gI1,2)*ZN(gI2,5) + KroneckerDelta(4,gO1)*(-(ZH(gI1,1)*(10*gp*QHu*ZN(gI2,0
      ) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(gI2,2))) + 7.0710678118654755*
      Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,3) + ZH(gI1,0)*ZN(gI2,5))) + KroneckerDelta(
      3,gO1)*(ZH(gI1,0)*(-10*gp*QHd*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5
      *g2*ZN(gI2,2)) + 7.0710678118654755*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,4) + ZH(
      gI1,1)*ZN(gI2,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*Conj(ZA(gI2,2))*(2*gp*Qs*Conj(ZN(gI1
      ,5))*KroneckerDelta(0,gO2) + 2*gp*Qs*Conj(ZN(gI1,0))*KroneckerDelta(5,gO2) +
      1.4142135623730951*Conj(ZN(gI1,4))*KroneckerDelta(3,gO2)*Lambdax +
      1.4142135623730951*Conj(ZN(gI1,3))*KroneckerDelta(4,gO2)*Lambdax) + Conj(ZA(
      gI2,1))*(Conj(ZN(gI1,4))*(10*gp*QHu*KroneckerDelta(0,gO2) +
      3.872983346207417*g1*KroneckerDelta(1,gO2) - 5*g2*KroneckerDelta(2,gO2)) +
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

std::complex<double> CLASSNAME::CpUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*Conj(ZA(gI2,2))*(2*gp*Qs*
      KroneckerDelta(5,gO1)*ZN(gI1,0) + 1.4142135623730951*Conj(Lambdax)*(
      KroneckerDelta(4,gO1)*ZN(gI1,3) + KroneckerDelta(3,gO1)*ZN(gI1,4)) + 2*gp*Qs
      *KroneckerDelta(0,gO1)*ZN(gI1,5)) + Conj(ZA(gI2,0))*(KroneckerDelta(3,gO1)*(
      10*gp*QHd*ZN(gI1,0) - 3.872983346207417*g1*ZN(gI1,1) + 5*g2*ZN(gI1,2)) + 10*
      gp*QHd*KroneckerDelta(0,gO1)*ZN(gI1,3) - 3.872983346207417*g1*KroneckerDelta
      (1,gO1)*ZN(gI1,3) + 5*g2*KroneckerDelta(2,gO1)*ZN(gI1,3) +
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(5,gO1)*ZN(gI1,4) +
      7.0710678118654755*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZN(gI1,5)) + Conj(ZA(
      gI2,1))*(KroneckerDelta(4,gO1)*(10*gp*QHu*ZN(gI1,0) + 3.872983346207417*g1*
      ZN(gI1,1) - 5*g2*ZN(gI1,2)) + (10*gp*QHu*KroneckerDelta(0,gO1) +
      3.872983346207417*g1*KroneckerDelta(1,gO1) - 5*g2*KroneckerDelta(2,gO1))*ZN(
      gI1,4) + 7.0710678118654755*Conj(Lambdax)*(KroneckerDelta(5,gO1)*ZN(gI1,3) +
      KroneckerDelta(3,gO1)*ZN(gI1,5))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_2298;
   std::complex<double> tmp_2299;
   std::complex<double> tmp_2300;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2300 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2299 += tmp_2300;
   tmp_2298 += (std::complex<double>(0,-1.4142135623730951)*gp*Qq*
      KroneckerDelta(0,gO2)) * tmp_2299;
   std::complex<double> tmp_2301;
   std::complex<double> tmp_2302;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2302 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2301 += tmp_2302;
   tmp_2298 += (std::complex<double>(0,-0.18257418583505536)*g1*KroneckerDelta(
      1,gO2)) * tmp_2301;
   std::complex<double> tmp_2303;
   std::complex<double> tmp_2304;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2304 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2303 += tmp_2304;
   tmp_2298 += (std::complex<double>(0,0.7071067811865475)*g2*KroneckerDelta(2,
      gO2)) * tmp_2303;
   std::complex<double> tmp_2305;
   std::complex<double> tmp_2306;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2307;
      std::complex<double> tmp_2308;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2308 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2307 += tmp_2308;
      tmp_2306 += (Conj(ZDL(gI2,j2))) * tmp_2307;
   }
   tmp_2305 += tmp_2306;
   tmp_2298 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_2305;
   result += (std::complex<double>(0,-1)) * tmp_2298;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2309;
   std::complex<double> tmp_2310;
   std::complex<double> tmp_2311;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2311 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_2310 += tmp_2311;
   tmp_2309 += (std::complex<double>(0,-1.4142135623730951)*gp*Qd*
      KroneckerDelta(0,gO1)) * tmp_2310;
   std::complex<double> tmp_2312;
   std::complex<double> tmp_2313;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2313 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_2312 += tmp_2313;
   tmp_2309 += (std::complex<double>(0,-0.3651483716701107)*g1*KroneckerDelta(1
      ,gO1)) * tmp_2312;
   std::complex<double> tmp_2314;
   std::complex<double> tmp_2315;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2316;
      std::complex<double> tmp_2317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2317 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2316 += tmp_2317;
      tmp_2315 += (ZD(gI1,j2)) * tmp_2316;
   }
   tmp_2314 += tmp_2315;
   tmp_2309 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_2314;
   result += (std::complex<double>(0,-1)) * tmp_2309;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_2318;
   std::complex<double> tmp_2319;
   std::complex<double> tmp_2320;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2320 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2319 += tmp_2320;
   tmp_2318 += (std::complex<double>(0,-1.4142135623730951)*gp*Ql*
      KroneckerDelta(0,gO2)) * tmp_2319;
   std::complex<double> tmp_2321;
   std::complex<double> tmp_2322;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2322 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2321 += tmp_2322;
   tmp_2318 += (std::complex<double>(0,0.5477225575051661)*g1*KroneckerDelta(1,
      gO2)) * tmp_2321;
   std::complex<double> tmp_2323;
   std::complex<double> tmp_2324;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2324 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2323 += tmp_2324;
   tmp_2318 += (std::complex<double>(0,0.7071067811865475)*g2*KroneckerDelta(2,
      gO2)) * tmp_2323;
   std::complex<double> tmp_2325;
   std::complex<double> tmp_2326;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2327;
      std::complex<double> tmp_2328;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2328 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2327 += tmp_2328;
      tmp_2326 += (Conj(ZEL(gI2,j2))) * tmp_2327;
   }
   tmp_2325 += tmp_2326;
   tmp_2318 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_2325;
   result += (std::complex<double>(0,-1)) * tmp_2318;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2329;
   std::complex<double> tmp_2330;
   std::complex<double> tmp_2331;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2331 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_2330 += tmp_2331;
   tmp_2329 += (std::complex<double>(0,-1.4142135623730951)*gp*Qe*
      KroneckerDelta(0,gO1)) * tmp_2330;
   std::complex<double> tmp_2332;
   std::complex<double> tmp_2333;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2333 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_2332 += tmp_2333;
   tmp_2329 += (std::complex<double>(0,-1.0954451150103321)*g1*KroneckerDelta(1
      ,gO1)) * tmp_2332;
   std::complex<double> tmp_2334;
   std::complex<double> tmp_2335;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2336;
      std::complex<double> tmp_2337;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2337 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2336 += tmp_2337;
      tmp_2335 += (ZE(gI1,j2)) * tmp_2336;
   }
   tmp_2334 += tmp_2335;
   tmp_2329 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_2334;
   result += (std::complex<double>(0,-1)) * tmp_2329;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_2338;
   std::complex<double> tmp_2339;
   std::complex<double> tmp_2340;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2340 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2339 += tmp_2340;
   tmp_2338 += (std::complex<double>(0,-1.4142135623730951)*gp*Qq*
      KroneckerDelta(0,gO2)) * tmp_2339;
   std::complex<double> tmp_2341;
   std::complex<double> tmp_2342;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2342 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2341 += tmp_2342;
   tmp_2338 += (std::complex<double>(0,-0.18257418583505536)*g1*KroneckerDelta(
      1,gO2)) * tmp_2341;
   std::complex<double> tmp_2343;
   std::complex<double> tmp_2344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2344 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2343 += tmp_2344;
   tmp_2338 += (std::complex<double>(0,-0.7071067811865475)*g2*KroneckerDelta(2
      ,gO2)) * tmp_2343;
   std::complex<double> tmp_2345;
   std::complex<double> tmp_2346;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2347;
      std::complex<double> tmp_2348;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2348 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2347 += tmp_2348;
      tmp_2346 += (Conj(ZUL(gI2,j2))) * tmp_2347;
   }
   tmp_2345 += tmp_2346;
   tmp_2338 += (std::complex<double>(0,-1)*KroneckerDelta(4,gO2)) * tmp_2345;
   result += (std::complex<double>(0,-1)) * tmp_2338;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2349;
   std::complex<double> tmp_2350;
   std::complex<double> tmp_2351;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2351 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_2350 += tmp_2351;
   tmp_2349 += (std::complex<double>(0,-1.4142135623730951)*gp*Qu*
      KroneckerDelta(0,gO1)) * tmp_2350;
   std::complex<double> tmp_2352;
   std::complex<double> tmp_2353;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2353 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_2352 += tmp_2353;
   tmp_2349 += (std::complex<double>(0,0.7302967433402214)*g1*KroneckerDelta(1,
      gO1)) * tmp_2352;
   std::complex<double> tmp_2354;
   std::complex<double> tmp_2355;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2356;
      std::complex<double> tmp_2357;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2357 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2356 += tmp_2357;
      tmp_2355 += (ZU(gI1,j2)) * tmp_2356;
   }
   tmp_2354 += tmp_2355;
   tmp_2349 += (std::complex<double>(0,-1)*KroneckerDelta(4,gO1)) * tmp_2354;
   result += (std::complex<double>(0,-1)) * tmp_2349;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(2,gO2)*UP(gI2,0)) + 0.7071067811865475*g2*
      KroneckerDelta(4,gO2)*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(2,gO1) +
      1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(3,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZChiPR(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(KroneckerDelta(3,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())
      )*ZN(gI2,3) - KroneckerDelta(4,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) +
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())
      )*ZN(gI2,4) + 10*gp*Qs*KroneckerDelta(5,gO2)*Sin(ThetaWp())*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZChiPL(unsigned gO1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(-10*gp*Qs*Conj(ZN(gI2,5))*KroneckerDelta(5,gO1)*Sin(ThetaWp())
      - Conj(ZN(gI2,3))*KroneckerDelta(3,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp())
      + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp(
      ))) + Conj(ZN(gI2,4))*KroneckerDelta(4,gO1)*(5*g2*Cos(ThetaW())*Cos(ThetaWp(
      )) + 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHu*Sin(
      ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZpChiPR(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(KroneckerDelta(3,gO2)*(10*gp*QHd*Cos(ThetaWp()) - (5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*ZN(gI2,3) +
      KroneckerDelta(4,gO2)*(10*gp*QHu*Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(
      ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(ThetaWp()))*ZN(gI2,4) +
      10*gp*Qs*Cos(ThetaWp())*KroneckerDelta(5,gO2)*ZN(gI2,5));

   return result;
}

std::complex<double> CLASSNAME::CpUChiVZpChiPL(unsigned gO1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(-10*gp*Qs*Conj(ZN(gI2,5))*Cos(ThetaWp())*KroneckerDelta(5,gO1)
      - Conj(ZN(gI2,4))*KroneckerDelta(4,gO1)*(10*gp*QHu*Cos(ThetaWp()) + 5*g2*
      Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*Sin(ThetaW())*Sin(
      ThetaWp())) + Conj(ZN(gI2,3))*KroneckerDelta(3,gO1)*(-10*gp*QHd*Cos(ThetaWp(
      )) + (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()
      )));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.7071067811865475)*(g2*Conj(UM(gI1,0))*Conj
      (ZA(gI2,1))*KroneckerDelta(1,gO2) + Conj(UM(gI1,1))*(g2*Conj(ZA(gI2,0))*
      KroneckerDelta(0,gO2) - Conj(ZA(gI2,2))*KroneckerDelta(1,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.7071067811865475)*(g2*Conj(ZA(gI2,0))*
      KroneckerDelta(1,gO1)*UP(gI1,0) + (g2*Conj(ZA(gI2,1))*KroneckerDelta(0,gO1)
      - Conj(Lambdax)*Conj(ZA(gI2,2))*KroneckerDelta(1,gO1))*UP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -(Conj(ZN(gI2,5))*KroneckerDelta(1,gO2)*Lambdax*ZP(gI1,0)) - 0.1*(
      10*g2*Conj(ZN(gI2,4))*KroneckerDelta(0,gO2) + 1.4142135623730951*(10*gp*QHu*
      Conj(ZN(gI2,0)) + 3.872983346207417*g1*Conj(ZN(gI2,1)) + 5*g2*Conj(ZN(gI2,2)
      ))*KroneckerDelta(1,gO2))*ZP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO1)*ZN(gI2,3)*ZP(gI1,0)) + KroneckerDelta(1,
      gO1)*(-1.4142135623730951*gp*QHd*ZN(gI2,0)*ZP(gI1,0) + 0.5477225575051661*g1
      *ZN(gI2,1)*ZP(gI1,0) + 0.7071067811865475*g2*ZN(gI2,2)*ZP(gI1,0) - Conj(
      Lambdax)*ZN(gI2,5)*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*KroneckerDelta(1,gO2)*ZH(
      gI1,1) + Conj(UM(gI2,1))*(g2*KroneckerDelta(0,gO2)*ZH(gI1,0) +
      KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*KroneckerDelta(0,gO1)*UP(gI2,1)*ZH(gI1,1) +
      KroneckerDelta(1,gO1)*(g2*UP(gI2,0)*ZH(gI1,0) + Conj(Lambdax)*UP(gI2,1)*ZH(
      gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2358;
   std::complex<double> tmp_2359;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2359 += Conj(ZEL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2358 += tmp_2359;
   result += (-(g2*KroneckerDelta(0,gO2))) * tmp_2358;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2360;
   std::complex<double> tmp_2361;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2362;
      std::complex<double> tmp_2363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2363 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2362 += tmp_2363;
      tmp_2361 += (ZV(gI1,j2)) * tmp_2362;
   }
   tmp_2360 += tmp_2361;
   result += (KroneckerDelta(1,gO1)) * tmp_2360;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2364;
   std::complex<double> tmp_2365;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2366;
      std::complex<double> tmp_2367;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2367 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_2366 += tmp_2367;
      tmp_2365 += (Conj(ZD(gI2,j2))) * tmp_2366;
   }
   tmp_2364 += tmp_2365;
   result += (KroneckerDelta(1,gO2)) * tmp_2364;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2368;
   std::complex<double> tmp_2369;
   std::complex<double> tmp_2370;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2370 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_2369 += tmp_2370;
   tmp_2368 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_2369
      ;
   std::complex<double> tmp_2371;
   std::complex<double> tmp_2372;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2373;
      std::complex<double> tmp_2374;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2374 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2373 += tmp_2374;
      tmp_2372 += (ZUL(gI1,j2)) * tmp_2373;
   }
   tmp_2371 += tmp_2372;
   tmp_2368 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_2371;
   result += (std::complex<double>(0,-1)) * tmp_2368;

   return result;
}

double CLASSNAME::CpbarUChabarFvSePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2375;
   std::complex<double> tmp_2376;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2376 += Conj(Ye(j1,gI1))*Conj(ZE(gI2,3 + j1));
   }
   tmp_2375 += tmp_2376;
   result += (KroneckerDelta(1,gO1)) * tmp_2375;
   if (gI1 < 3) {
      result += -(g2*Conj(ZE(gI2,gI1))*KroneckerDelta(0,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2377;
   std::complex<double> tmp_2378;
   std::complex<double> tmp_2379;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2379 += Conj(ZDL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2378 += tmp_2379;
   tmp_2377 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_2378
      ;
   std::complex<double> tmp_2380;
   std::complex<double> tmp_2381;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2382;
      std::complex<double> tmp_2383;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2383 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2382 += tmp_2383;
      tmp_2381 += (Conj(ZDL(gI2,j2))) * tmp_2382;
   }
   tmp_2380 += tmp_2381;
   tmp_2377 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_2380;
   result += (std::complex<double>(0,-1)) * tmp_2377;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2384;
   std::complex<double> tmp_2385;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2386;
      std::complex<double> tmp_2387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2387 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2386 += tmp_2387;
      tmp_2385 += (ZU(gI1,j2)) * tmp_2386;
   }
   tmp_2384 += tmp_2385;
   result += (KroneckerDelta(1,gO1)) * tmp_2384;

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
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = g2*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(0,gO2)*UP(gI2,0) +
      0.1*KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) -
      3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) + 10*gp*QHu*Sin(ThetaWp())
      )*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPL(unsigned gO1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(0,
      gO1) + 0.1*Conj(UM(gI2,1))*KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW())*Cos(
      ThetaWp()) - 3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW()) - 10*gp*QHd*
      Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZpChaPR(unsigned gO2, unsigned gI2) const
{
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = -(g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*Sin(ThetaWp())*UP(gI2,0))
      + 0.1*KroneckerDelta(1,gO2)*(10*gp*QHu*Cos(ThetaWp()) + (-5*g2*Cos(ThetaW())
      + 3.872983346207417*g1*Sin(ThetaW()))*Sin(ThetaWp()))*UP(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZpChaPL(unsigned gO1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = -(g2*Conj(UM(gI2,0))*Cos(ThetaW())*KroneckerDelta(0,gO1)*Sin(
      ThetaWp())) - 0.1*Conj(UM(gI2,1))*KroneckerDelta(1,gO1)*(10*gp*QHd*Cos(
      ThetaWp()) + (5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*Sin(
      ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,2)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZN(gI2,4);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += Ye(gO2,gI2)*ZP(gI1,0);
   }

   return result;
}

double CLASSNAME::CpbarUFeHpmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2388;
      std::complex<double> tmp_2389;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2389 += Conj(ZV(gI1,j2))*Ye(gO2,j2);
      }
      tmp_2388 += tmp_2389;
      result += (Conj(UM(gI2,1))) * tmp_2388;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZV(gI1,gO1))*UP(gI2,0));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2390;
      std::complex<double> tmp_2391;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2391 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_2390 += tmp_2391;
      result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,0))
         ) * tmp_2390;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2392;
      std::complex<double> tmp_2393;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2393 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_2392 += tmp_2393;
      result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,0)))
         * tmp_2392;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2394;
      std::complex<double> tmp_2395;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2395 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_2394 += tmp_2395;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2394;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2396;
      std::complex<double> tmp_2397;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2397 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_2396 += tmp_2397;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2396;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.4142135623730951*gp*Qe*Conj(ZE(gI1,3 + gO2))*Conj(ZN(gI2,
         0));
   }
   if (gO2 < 3) {
      result += -1.0954451150103321*g1*Conj(ZE(gI1,3 + gO2))*Conj(ZN(gI2,1))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_2398;
      std::complex<double> tmp_2399;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2399 += Conj(ZE(gI1,j2))*Ye(gO2,j2);
      }
      tmp_2398 += tmp_2399;
      result += (-Conj(ZN(gI2,3))) * tmp_2398;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Ql*Conj(ZE(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZE(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZE(gI1,gO1))*ZN(gI2,2);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_2400;
      std::complex<double> tmp_2401;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2401 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_2400 += tmp_2401;
      result += (-ZN(gI2,3)) * tmp_2400;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Cos(ThetaW())*ZER(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVPFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZEL(gI2,gO1))*Sin(ThetaW());
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
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())*ZER(gI2,
         gO2);
   }
   if (gI2 < 3) {
      result += gp*Qe*Sin(ThetaWp())*ZER(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZEL(gI2,gO1))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Cos(ThetaWp())*Sin
         (ThetaW());
   }
   if (gI2 < 3) {
      result += -(gp*Ql*Conj(ZEL(gI2,gO1))*Sin(ThetaWp()));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZpFePR(unsigned gO2, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   if (gI2 < 3) {
      result += gp*Qe*Cos(ThetaWp())*ZER(gI2,gO2);
   }
   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())*ZER(gI2,
         gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZpFePL(unsigned gO1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -(gp*Ql*Conj(ZEL(gI2,gO1))*Cos(ThetaWp()));
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZEL(gI2,gO1))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2402;
      std::complex<double> tmp_2403;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2403 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_2402 += tmp_2403;
      result += (ZP(gI1,0)) * tmp_2402;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2404;
      std::complex<double> tmp_2405;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2405 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_2404 += tmp_2405;
      result += (ZP(gI1,1)) * tmp_2404;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2406;
      std::complex<double> tmp_2407;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2407 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_2406 += tmp_2407;
      result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,0))
         ) * tmp_2406;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2408;
      std::complex<double> tmp_2409;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2409 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_2408 += tmp_2409;
      result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,0)))
         * tmp_2408;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2410;
      std::complex<double> tmp_2411;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2411 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_2410 += tmp_2411;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2410;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2412;
      std::complex<double> tmp_2413;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2413 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_2412 += tmp_2413;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2412;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2414;
      std::complex<double> tmp_2415;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2415 += Conj(ZU(gI1,j2))*Yd(gO2,j2);
      }
      tmp_2414 += tmp_2415;
      result += (Conj(UM(gI2,1))) * tmp_2414;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZU(gI1,gO1))*UP(gI2,0));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_2416;
      std::complex<double> tmp_2417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2417 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2416 += tmp_2417;
      result += (UP(gI2,1)) * tmp_2416;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.4142135623730951*gp*Qd*Conj(ZD(gI1,3 + gO2))*Conj(ZN(gI2,
         0));
   }
   if (gO2 < 3) {
      result += -0.3651483716701107*g1*Conj(ZD(gI1,3 + gO2))*Conj(ZN(gI2,1))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_2418;
      std::complex<double> tmp_2419;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2419 += Conj(ZD(gI1,j2))*Yd(gO2,j2);
      }
      tmp_2418 += tmp_2419;
      result += (-Conj(ZN(gI2,3))) * tmp_2418;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Qq*Conj(ZD(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZD(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZD(gI1,gO1))*ZN(gI2,2);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_2420;
      std::complex<double> tmp_2421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2421 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2420 += tmp_2421;
      result += (-ZN(gI2,3)) * tmp_2420;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 1.4142135623730951*g3*PhaseGlu*Conj(ZD(gI1,3 + gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(PhaseGlu)*Conj(ZD(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*ZDR(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVGFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.2581988897471611*g1*Cos(ThetaW())*ZDR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVPFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZDL(gI2,gO1))*Sin(ThetaW());
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
      result += -0.7071067811865475*g2*Conj(ZUL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.2581988897471611*g1*Cos(ThetaWp())*Sin(ThetaW())*ZDR(gI2,
         gO2);
   }
   if (gI2 < 3) {
      result += gp*Qd*Sin(ThetaWp())*ZDR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZDL(gI2,gO1))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Cos(ThetaWp())*Sin
         (ThetaW());
   }
   if (gI2 < 3) {
      result += -(gp*Qq*Conj(ZDL(gI2,gO1))*Sin(ThetaWp()));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZpFdPR(unsigned gO2, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   if (gI2 < 3) {
      result += gp*Qd*Cos(ThetaWp())*ZDR(gI2,gO2);
   }
   if (gI2 < 3) {
      result += 0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp())*ZDR(gI2,
         gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZpFdPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -(gp*Qq*Conj(ZDL(gI2,gO1))*Cos(ThetaWp()));
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZDL(gI2,gO1))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Sin(ThetaW())*Sin
         (ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2422;
      std::complex<double> tmp_2423;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2423 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_2422 += tmp_2423;
      result += (ZP(gI1,1)) * tmp_2422;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2424;
      std::complex<double> tmp_2425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2425 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_2424 += tmp_2425;
      result += (ZP(gI1,0)) * tmp_2424;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2426;
      std::complex<double> tmp_2427;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2427 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_2426 += tmp_2427;
      result += (Conj(UP(gI1,1))) * tmp_2426;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZD(gI2,gO1))*UM(gI1,0));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_2428;
      std::complex<double> tmp_2429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2429 += Conj(Yd(j1,gO1))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2428 += tmp_2429;
      result += (UM(gI1,1)) * tmp_2428;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2430;
      std::complex<double> tmp_2431;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2431 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_2430 += tmp_2431;
      result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,1))
         ) * tmp_2430;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2432;
      std::complex<double> tmp_2433;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2433 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_2432 += tmp_2433;
      result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,1)))
         * tmp_2432;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_2434;
      std::complex<double> tmp_2435;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2435 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_2434 += tmp_2435;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2434;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_2436;
      std::complex<double> tmp_2437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2437 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_2436 += tmp_2437;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2436;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.4142135623730951*gp*Qu*Conj(ZN(gI2,0))*Conj(ZU(gI1,3 +
         gO2));
   }
   if (gO2 < 3) {
      result += 0.7302967433402214*g1*Conj(ZN(gI2,1))*Conj(ZU(gI1,3 + gO2));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_2438;
      std::complex<double> tmp_2439;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_2439 += Conj(ZU(gI1,j2))*Yu(gO2,j2);
      }
      tmp_2438 += tmp_2439;
      result += (-Conj(ZN(gI2,4))) * tmp_2438;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Qq*Conj(ZU(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZU(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZU(gI1,gO1))*ZN(gI2,2);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_2440;
      std::complex<double> tmp_2441;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2441 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2440 += tmp_2441;
      result += (-ZN(gI2,4)) * tmp_2440;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 1.4142135623730951*g3*PhaseGlu*Conj(ZU(gI1,3 + gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(PhaseGlu)*Conj(ZU(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*ZUR(gI2,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVGFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -(g3*Conj(ZUL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5163977794943222*g1*Cos(ThetaW())*ZUR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVPFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZUL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPR(unsigned gO2, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5163977794943222*g1*Cos(ThetaWp())*Sin(ThetaW())*ZUR(gI2,
         gO2);
   }
   if (gI2 < 3) {
      result += gp*Qu*Sin(ThetaWp())*ZUR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZUL(gI2,gO1))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Cos(ThetaWp())*Sin
         (ThetaW());
   }
   if (gI2 < 3) {
      result += -(gp*Qq*Conj(ZUL(gI2,gO1))*Sin(ThetaWp()));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZpFuPR(unsigned gO2, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   if (gI2 < 3) {
      result += gp*Qu*Cos(ThetaWp())*ZUR(gI2,gO2);
   }
   if (gI2 < 3) {
      result += -0.5163977794943222*g1*Sin(ThetaW())*Sin(ThetaWp())*ZUR(gI2,
         gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZpFuPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -(gp*Qq*Conj(ZUL(gI2,gO1))*Cos(ThetaWp()));
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZUL(gI2,gO1))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gI2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Sin(ThetaW())*Sin
         (ThetaWp());
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
      result += -0.7071067811865475*g2*Conj(ZDL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2442;
   std::complex<double> tmp_2443;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2443 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2442 += tmp_2443;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_2442;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2444;
   std::complex<double> tmp_2445;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2445 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_2444 += tmp_2445;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2444;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2446;
   std::complex<double> tmp_2447;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2447 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2446 += tmp_2447;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_2446;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2448;
   std::complex<double> tmp_2449;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2449 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_2448 += tmp_2449;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2448;

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

double CLASSNAME::CpconjVWmbarVWmVZp() const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpconjVWmbarVZpVWm() const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2450;
   std::complex<double> tmp_2451;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2451 += Conj(ZER(gO2,j1))*Ye(j1,gI2);
   }
   tmp_2450 += tmp_2451;
   result += (ZP(gI1,0)) * tmp_2450;

   return result;
}

double CLASSNAME::CpbarFeHpmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2452;
   std::complex<double> tmp_2453;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2454;
      std::complex<double> tmp_2455;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2455 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2454 += tmp_2455;
      tmp_2453 += (Conj(ZV(gI1,j2))) * tmp_2454;
   }
   tmp_2452 += tmp_2453;
   result += (Conj(UM(gI2,1))) * tmp_2452;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2456;
   std::complex<double> tmp_2457;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2457 += Conj(ZV(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_2456 += tmp_2457;
   result += (-(g2*UP(gI2,0))) * tmp_2456;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2458;
   std::complex<double> tmp_2459;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2460;
      std::complex<double> tmp_2461;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2461 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2460 += tmp_2461;
      tmp_2459 += (Conj(ZEL(gI1,j2))) * tmp_2460;
   }
   tmp_2458 += tmp_2459;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_2458;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2462;
   std::complex<double> tmp_2463;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2464;
      std::complex<double> tmp_2465;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2465 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_2464 += tmp_2465;
      tmp_2463 += (ZEL(gO1,j2)) * tmp_2464;
   }
   tmp_2462 += tmp_2463;
   result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_2462;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2466;
   std::complex<double> tmp_2467;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2468;
      std::complex<double> tmp_2469;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2469 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2468 += tmp_2469;
      tmp_2467 += (Conj(ZEL(gI2,j2))) * tmp_2468;
   }
   tmp_2466 += tmp_2467;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2466;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2470;
   std::complex<double> tmp_2471;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2472;
      std::complex<double> tmp_2473;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2473 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2472 += tmp_2473;
      tmp_2471 += (ZEL(gO1,j2)) * tmp_2472;
   }
   tmp_2470 += tmp_2471;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2470;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2474;
   std::complex<double> tmp_2475;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2475 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_2474 += tmp_2475;
   result += (-1.4142135623730951*gp*Qe*Conj(ZN(gI2,0))) * tmp_2474;
   std::complex<double> tmp_2476;
   std::complex<double> tmp_2477;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2477 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_2476 += tmp_2477;
   result += (-1.0954451150103321*g1*Conj(ZN(gI2,1))) * tmp_2476;
   std::complex<double> tmp_2478;
   std::complex<double> tmp_2479;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2480;
      std::complex<double> tmp_2481;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2481 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2480 += tmp_2481;
      tmp_2479 += (Conj(ZE(gI1,j2))) * tmp_2480;
   }
   tmp_2478 += tmp_2479;
   result += (-Conj(ZN(gI2,3))) * tmp_2478;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_2482;
   std::complex<double> tmp_2483;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2483 += Conj(ZE(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_2482 += tmp_2483;
   result += (-0.7071067811865475*(2*gp*Ql*ZN(gI2,0) - 0.7745966692414834*g1*ZN
      (gI2,1) - g2*ZN(gI2,2))) * tmp_2482;
   std::complex<double> tmp_2484;
   std::complex<double> tmp_2485;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2486;
      std::complex<double> tmp_2487;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2487 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_2486 += tmp_2487;
      tmp_2485 += (ZEL(gO1,j2)) * tmp_2486;
   }
   tmp_2484 += tmp_2485;
   result += (-ZN(gI2,3)) * tmp_2484;

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
      result += -0.7071067811865475*g2*ZEL(gO1,gI2);
   }

   return result;
}

double CLASSNAME::CpbarFeVZFePR(unsigned gO2, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   double result = 0.0;

   result = -(KroneckerDelta(gI2,gO2)*(0.7745966692414834*g1*Cos(ThetaWp())*Sin
      (ThetaW()) - gp*Qe*Sin(ThetaWp())));

   return result;
}

double CLASSNAME::CpbarFeVZFePL(unsigned gO1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = 0.5*KroneckerDelta(gI2,gO1)*(g2*Cos(ThetaW())*Cos(ThetaWp()) -
      0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*Ql*Sin(ThetaWp()))
      ;

   return result;
}

double CLASSNAME::CpbarFeVZpFePR(unsigned gO2, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   double result = 0.0;

   result = KroneckerDelta(gI2,gO2)*(gp*Qe*Cos(ThetaWp()) + 0.7745966692414834*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFeVZpFePL(unsigned gO1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   double result = 0.0;

   result = -0.5*KroneckerDelta(gI2,gO1)*(2*gp*Ql*Cos(ThetaWp()) + g2*Cos(
      ThetaW())*Sin(ThetaWp()) - 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp()
      ));

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2488;
   std::complex<double> tmp_2489;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2490;
      std::complex<double> tmp_2491;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2491 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2490 += tmp_2491;
      tmp_2489 += (Conj(ZUL(gI2,j2))) * tmp_2490;
   }
   tmp_2488 += tmp_2489;
   result += (ZP(gI1,0)) * tmp_2488;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2492;
   std::complex<double> tmp_2493;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2494;
      std::complex<double> tmp_2495;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2495 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2494 += tmp_2495;
      tmp_2493 += (ZDL(gO1,j2)) * tmp_2494;
   }
   tmp_2492 += tmp_2493;
   result += (ZP(gI1,1)) * tmp_2492;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2496;
   std::complex<double> tmp_2497;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2498;
      std::complex<double> tmp_2499;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2499 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2498 += tmp_2499;
      tmp_2497 += (Conj(ZDL(gI1,j2))) * tmp_2498;
   }
   tmp_2496 += tmp_2497;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_2496;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2500;
   std::complex<double> tmp_2501;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2502;
      std::complex<double> tmp_2503;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2503 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_2502 += tmp_2503;
      tmp_2501 += (ZDL(gO1,j2)) * tmp_2502;
   }
   tmp_2500 += tmp_2501;
   result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_2500;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2504;
   std::complex<double> tmp_2505;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2506;
      std::complex<double> tmp_2507;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2507 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2506 += tmp_2507;
      tmp_2505 += (Conj(ZDL(gI2,j2))) * tmp_2506;
   }
   tmp_2504 += tmp_2505;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2504;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2508;
   std::complex<double> tmp_2509;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2510;
      std::complex<double> tmp_2511;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2511 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2510 += tmp_2511;
      tmp_2509 += (ZDL(gO1,j2)) * tmp_2510;
   }
   tmp_2508 += tmp_2509;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2508;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2512;
   std::complex<double> tmp_2513;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2514;
      std::complex<double> tmp_2515;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2515 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2514 += tmp_2515;
      tmp_2513 += (Conj(ZU(gI1,j2))) * tmp_2514;
   }
   tmp_2512 += tmp_2513;
   result += (Conj(UM(gI2,1))) * tmp_2512;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2516;
   std::complex<double> tmp_2517;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2517 += Conj(ZU(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2516 += tmp_2517;
   result += (-(g2*UP(gI2,0))) * tmp_2516;
   std::complex<double> tmp_2518;
   std::complex<double> tmp_2519;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2520;
      std::complex<double> tmp_2521;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2521 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2520 += tmp_2521;
      tmp_2519 += (ZDL(gO1,j2)) * tmp_2520;
   }
   tmp_2518 += tmp_2519;
   result += (UP(gI2,1)) * tmp_2518;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2522;
   std::complex<double> tmp_2523;
   std::complex<double> tmp_2524;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2524 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2523 += tmp_2524;
   tmp_2522 += (-4.242640687119286*gp*Qd*Conj(ZN(gI2,0))) * tmp_2523;
   std::complex<double> tmp_2525;
   std::complex<double> tmp_2526;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2526 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2525 += tmp_2526;
   tmp_2522 += (-1.0954451150103321*g1*Conj(ZN(gI2,1))) * tmp_2525;
   std::complex<double> tmp_2527;
   std::complex<double> tmp_2528;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2529;
      std::complex<double> tmp_2530;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2530 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2529 += tmp_2530;
      tmp_2528 += (Conj(ZD(gI1,j2))) * tmp_2529;
   }
   tmp_2527 += tmp_2528;
   tmp_2522 += (-3*Conj(ZN(gI2,3))) * tmp_2527;
   result += (0.3333333333333333) * tmp_2522;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_2531;
   std::complex<double> tmp_2532;
   std::complex<double> tmp_2533;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2533 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2532 += tmp_2533;
   tmp_2531 += (-1.4142135623730951*(6*gp*Qq*ZN(gI2,0) + 0.7745966692414834*g1*
      ZN(gI2,1) - 3*g2*ZN(gI2,2))) * tmp_2532;
   std::complex<double> tmp_2534;
   std::complex<double> tmp_2535;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2536;
      std::complex<double> tmp_2537;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2537 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2536 += tmp_2537;
      tmp_2535 += (ZDL(gO1,j2)) * tmp_2536;
   }
   tmp_2534 += tmp_2535;
   tmp_2531 += (-6*ZN(gI2,3)) * tmp_2534;
   result += (0.16666666666666666) * tmp_2531;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2538;
   std::complex<double> tmp_2539;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2539 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2538 += tmp_2539;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2538;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2540;
   std::complex<double> tmp_2541;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2541 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2540 += tmp_2541;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2540;

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

   std::complex<double> tmp_2542;
   std::complex<double> tmp_2543;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2543 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2542 += tmp_2543;
   result += (-0.7071067811865475*g2) * tmp_2542;

   return result;
}

double CLASSNAME::CpbarFdVZFdPR(unsigned gO2, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   double result = 0.0;

   result = KroneckerDelta(gI2,gO2)*(-0.2581988897471611*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + gp*Qd*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = 0.16666666666666666*KroneckerDelta(gI2,gO1)*(3*g2*Cos(ThetaW())*Cos
      (ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 6*gp*Qq*
      Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdVZpFdPR(unsigned gO2, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   double result = 0.0;

   result = KroneckerDelta(gI2,gO2)*(gp*Qd*Cos(ThetaWp()) + 0.2581988897471611*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFdVZpFdPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(6*gp*Qq*Cos(ThetaWp()
      ) + 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2544;
   std::complex<double> tmp_2545;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2546;
      std::complex<double> tmp_2547;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2547 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2546 += tmp_2547;
      tmp_2545 += (Conj(ZDL(gI2,j2))) * tmp_2546;
   }
   tmp_2544 += tmp_2545;
   result += (ZP(gI1,1)) * tmp_2544;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2548;
   std::complex<double> tmp_2549;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2550;
      std::complex<double> tmp_2551;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2551 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2550 += tmp_2551;
      tmp_2549 += (ZUL(gO1,j2)) * tmp_2550;
   }
   tmp_2548 += tmp_2549;
   result += (ZP(gI1,0)) * tmp_2548;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2552;
   std::complex<double> tmp_2553;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2554;
      std::complex<double> tmp_2555;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2555 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2554 += tmp_2555;
      tmp_2553 += (Conj(ZD(gI2,j2))) * tmp_2554;
   }
   tmp_2552 += tmp_2553;
   result += (Conj(UP(gI1,1))) * tmp_2552;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2556;
   std::complex<double> tmp_2557;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2557 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2556 += tmp_2557;
   result += (-(g2*UM(gI1,0))) * tmp_2556;
   std::complex<double> tmp_2558;
   std::complex<double> tmp_2559;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2560;
      std::complex<double> tmp_2561;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2561 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2560 += tmp_2561;
      tmp_2559 += (ZUL(gO1,j2)) * tmp_2560;
   }
   tmp_2558 += tmp_2559;
   result += (UM(gI1,1)) * tmp_2558;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2562;
   std::complex<double> tmp_2563;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2564;
      std::complex<double> tmp_2565;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2565 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2564 += tmp_2565;
      tmp_2563 += (Conj(ZUL(gI1,j2))) * tmp_2564;
   }
   tmp_2562 += tmp_2563;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_2562;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2566;
   std::complex<double> tmp_2567;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2568;
      std::complex<double> tmp_2569;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2569 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_2568 += tmp_2569;
      tmp_2567 += (ZUL(gO1,j2)) * tmp_2568;
   }
   tmp_2566 += tmp_2567;
   result += (std::complex<double>(0,0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_2566;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2570;
   std::complex<double> tmp_2571;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2572;
      std::complex<double> tmp_2573;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2573 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2572 += tmp_2573;
      tmp_2571 += (Conj(ZUL(gI2,j2))) * tmp_2572;
   }
   tmp_2570 += tmp_2571;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2570;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2574;
   std::complex<double> tmp_2575;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2576;
      std::complex<double> tmp_2577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2577 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2576 += tmp_2577;
      tmp_2575 += (ZUL(gO1,j2)) * tmp_2576;
   }
   tmp_2574 += tmp_2575;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2574;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2578;
   std::complex<double> tmp_2579;
   std::complex<double> tmp_2580;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2580 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2579 += tmp_2580;
   tmp_2578 += (-4.242640687119286*gp*Qu*Conj(ZN(gI2,0))) * tmp_2579;
   std::complex<double> tmp_2581;
   std::complex<double> tmp_2582;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2582 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2581 += tmp_2582;
   tmp_2578 += (2.1908902300206643*g1*Conj(ZN(gI2,1))) * tmp_2581;
   std::complex<double> tmp_2583;
   std::complex<double> tmp_2584;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2585;
      std::complex<double> tmp_2586;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2586 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2585 += tmp_2586;
      tmp_2584 += (Conj(ZU(gI1,j2))) * tmp_2585;
   }
   tmp_2583 += tmp_2584;
   tmp_2578 += (-3*Conj(ZN(gI2,4))) * tmp_2583;
   result += (0.3333333333333333) * tmp_2578;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_2587;
   std::complex<double> tmp_2588;
   std::complex<double> tmp_2589;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2589 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2588 += tmp_2589;
   tmp_2587 += (-1.4142135623730951*(6*gp*Qq*ZN(gI2,0) + 0.7745966692414834*g1*
      ZN(gI2,1) + 3*g2*ZN(gI2,2))) * tmp_2588;
   std::complex<double> tmp_2590;
   std::complex<double> tmp_2591;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2592;
      std::complex<double> tmp_2593;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2593 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2592 += tmp_2593;
      tmp_2591 += (ZUL(gO1,j2)) * tmp_2592;
   }
   tmp_2590 += tmp_2591;
   tmp_2587 += (-6*ZN(gI2,4)) * tmp_2590;
   result += (0.16666666666666666) * tmp_2587;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2594;
   std::complex<double> tmp_2595;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2595 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2594 += tmp_2595;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2594;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2596;
   std::complex<double> tmp_2597;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2597 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2596 += tmp_2597;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2596;

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
   const auto Qu = LOCALINPUT(Qu);

   double result = 0.0;

   result = KroneckerDelta(gI2,gO2)*(0.5163977794943222*g1*Cos(ThetaWp())*Sin(
      ThetaW()) + gp*Qu*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(3*g2*Cos(ThetaW())*
      Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 6*gp*
      Qq*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuVZpFuPR(unsigned gO2, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   double result = 0.0;

   result = KroneckerDelta(gI2,gO2)*(gp*Qu*Cos(ThetaWp()) - 0.5163977794943222*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

double CLASSNAME::CpbarFuVZpFuPL(unsigned gO1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   double result = 0.0;

   result = -0.16666666666666666*KroneckerDelta(gI2,gO1)*(6*gp*Qq*Cos(ThetaWp()
      ) - 3*g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()));

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

   std::complex<double> tmp_2598;
   std::complex<double> tmp_2599;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2599 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2598 += tmp_2599;
   result += (-0.7071067811865475*g2) * tmp_2598;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSdconjUSdVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_2600;
   std::complex<double> tmp_2601;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2601 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2600 += tmp_2601;
   result += (-1) * tmp_2600;
   std::complex<double> tmp_2602;
   std::complex<double> tmp_2603;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2603 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2602 += tmp_2603;
   result += (-0.5) * tmp_2602;
   std::complex<double> tmp_2604;
   std::complex<double> tmp_2605;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2605 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2604 += tmp_2605;
   result += (-1) * tmp_2604;
   std::complex<double> tmp_2606;
   std::complex<double> tmp_2607;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2607 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2606 += tmp_2607;
   result += (-0.5) * tmp_2606;
   std::complex<double> tmp_2608;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2609;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2609 += (Conj(CpconjUSdFuChaPL(gO2,gI1,gI2))*
            CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPR(gO2,gI1,gI2))*
            CpconjUSdFuChaPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha(gI2));
      }
      tmp_2608 += tmp_2609;
   }
   result += tmp_2608;
   std::complex<double> tmp_2610;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2611;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2611 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_2610 += tmp_2611;
   }
   result += tmp_2610;
   std::complex<double> tmp_2612;
   std::complex<double> tmp_2613;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2614;
      std::complex<double> tmp_2615;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2615 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2614 += tmp_2615;
      tmp_2613 += (MFd(gI1)) * tmp_2614;
   }
   tmp_2612 += tmp_2613;
   result += (-2) * tmp_2612;
   std::complex<double> tmp_2616;
   std::complex<double> tmp_2617;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2618;
      std::complex<double> tmp_2619;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2619 += B0(p,MFu(gI1),MCha(gI2))*(Conj(CpconjUSdFuChaPR(gO2,
            gI1,gI2))*CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPL(gO2,
            gI1,gI2))*CpconjUSdFuChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2618 += tmp_2619;
      tmp_2617 += (MFu(gI1)) * tmp_2618;
   }
   tmp_2616 += tmp_2617;
   result += (-2) * tmp_2616;
   std::complex<double> tmp_2620;
   std::complex<double> tmp_2621;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2621 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2620 += tmp_2621;
   result += (-1) * tmp_2620;
   std::complex<double> tmp_2622;
   std::complex<double> tmp_2623;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2623 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2622 += tmp_2623;
   result += (-1) * tmp_2622;
   std::complex<double> tmp_2624;
   std::complex<double> tmp_2625;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2625 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2624 += tmp_2625;
   result += (-1) * tmp_2624;
   std::complex<double> tmp_2626;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2627;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2627 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_2626 += tmp_2627;
   }
   result += tmp_2626;
   std::complex<double> tmp_2628;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2629;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2629 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_2628 += tmp_2629;
   }
   result += tmp_2628;
   std::complex<double> tmp_2630;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2631;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2631 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_2630 += tmp_2631;
   }
   result += tmp_2630;
   std::complex<double> tmp_2632;
   std::complex<double> tmp_2633;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2633 += (Conj(CpconjUSdGluFdPL(gO2,1,gI2))*CpconjUSdGluFdPL(gO1,1,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,1,gI2))*CpconjUSdGluFdPR(gO1,1,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_2632 += tmp_2633;
   result += (1.3333333333333333) * tmp_2632;
   std::complex<double> tmp_2634;
   std::complex<double> tmp_2635;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2635 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_2634 += tmp_2635;
   result += (1.3333333333333333) * tmp_2634;
   std::complex<double> tmp_2636;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2636 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_2636;
   std::complex<double> tmp_2637;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2637 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_2637;
   std::complex<double> tmp_2638;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2638 += Conj(CpconjUSdVZpSd(gO2,gI2))*CpconjUSdVZpSd(gO1,gI2)*F0(p
         ,MSd(gI2),MVZp);
   }
   result += tmp_2638;
   std::complex<double> tmp_2639;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2639 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_2639;
   std::complex<double> tmp_2640;
   std::complex<double> tmp_2641;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2641 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,1,gI2))*
         CpconjUSdGluFdPL(gO1,1,gI2) + Conj(CpconjUSdGluFdPL(gO2,1,gI2))*
         CpconjUSdGluFdPR(gO1,1,gI2))*MFd(gI2);
   }
   tmp_2640 += tmp_2641;
   result += (-2.6666666666666665*MGlu) * tmp_2640;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSvconjUSvVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_2642;
   std::complex<double> tmp_2643;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2643 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2642 += tmp_2643;
   result += (-1) * tmp_2642;
   std::complex<double> tmp_2644;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2645;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2645 += (Conj(CpconjUSvbarChaFePL(gO2,gI1,gI2))*
            CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(CpconjUSvbarChaFePR(gO2,gI1,
            gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFe(gI2));
      }
      tmp_2644 += tmp_2645;
   }
   result += tmp_2644;
   std::complex<double> tmp_2646;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2647;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2647 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_2646 += tmp_2647;
   }
   result += tmp_2646;
   std::complex<double> tmp_2648;
   std::complex<double> tmp_2649;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2650;
      std::complex<double> tmp_2651;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2651 += B0(p,MCha(gI1),MFe(gI2))*(Conj(CpconjUSvbarChaFePR(
            gO2,gI1,gI2))*CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(
            CpconjUSvbarChaFePL(gO2,gI1,gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2650 += tmp_2651;
      tmp_2649 += (MCha(gI1)) * tmp_2650;
   }
   tmp_2648 += tmp_2649;
   result += (-2) * tmp_2648;
   std::complex<double> tmp_2652;
   std::complex<double> tmp_2653;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2653 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2652 += tmp_2653;
   result += (-0.5) * tmp_2652;
   std::complex<double> tmp_2654;
   std::complex<double> tmp_2655;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2655 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2654 += tmp_2655;
   result += (-1) * tmp_2654;
   std::complex<double> tmp_2656;
   std::complex<double> tmp_2657;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2657 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2656 += tmp_2657;
   result += (-0.5) * tmp_2656;
   std::complex<double> tmp_2658;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2659;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2659 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_2658 += tmp_2659;
   }
   result += tmp_2658;
   std::complex<double> tmp_2660;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2661;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2661 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_2660 += tmp_2661;
   }
   result += tmp_2660;
   std::complex<double> tmp_2662;
   std::complex<double> tmp_2663;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2664;
      std::complex<double> tmp_2665;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2665 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2664 += tmp_2665;
      tmp_2663 += (MFv(gI1)) * tmp_2664;
   }
   tmp_2662 += tmp_2663;
   result += (-2) * tmp_2662;
   std::complex<double> tmp_2666;
   std::complex<double> tmp_2667;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2667 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2666 += tmp_2667;
   result += (-3) * tmp_2666;
   std::complex<double> tmp_2668;
   std::complex<double> tmp_2669;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2669 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2668 += tmp_2669;
   result += (-1) * tmp_2668;
   std::complex<double> tmp_2670;
   std::complex<double> tmp_2671;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2671 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2670 += tmp_2671;
   result += (-3) * tmp_2670;
   std::complex<double> tmp_2672;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2672 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_2672;
   std::complex<double> tmp_2673;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2673 += Conj(CpconjUSvVZpSv(gO2,gI2))*CpconjUSvVZpSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVZp);
   }
   result += tmp_2673;
   std::complex<double> tmp_2674;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2674 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_2674;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSuconjUSuVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_2675;
   std::complex<double> tmp_2676;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2676 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2675 += tmp_2676;
   result += (-1) * tmp_2675;
   std::complex<double> tmp_2677;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2678;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2678 += (Conj(CpconjUSubarChaFdPL(gO2,gI1,gI2))*
            CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(CpconjUSubarChaFdPR(gO2,gI1,
            gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFd(gI2));
      }
      tmp_2677 += tmp_2678;
   }
   result += tmp_2677;
   std::complex<double> tmp_2679;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2680;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2680 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_2679 += tmp_2680;
   }
   result += tmp_2679;
   std::complex<double> tmp_2681;
   std::complex<double> tmp_2682;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2683;
      std::complex<double> tmp_2684;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2684 += B0(p,MCha(gI1),MFd(gI2))*(Conj(CpconjUSubarChaFdPR(
            gO2,gI1,gI2))*CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChaFdPL(gO2,gI1,gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2683 += tmp_2684;
      tmp_2682 += (MCha(gI1)) * tmp_2683;
   }
   tmp_2681 += tmp_2682;
   result += (-2) * tmp_2681;
   std::complex<double> tmp_2685;
   std::complex<double> tmp_2686;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2686 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2685 += tmp_2686;
   result += (-0.5) * tmp_2685;
   std::complex<double> tmp_2687;
   std::complex<double> tmp_2688;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2688 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2687 += tmp_2688;
   result += (-1) * tmp_2687;
   std::complex<double> tmp_2689;
   std::complex<double> tmp_2690;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2690 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2689 += tmp_2690;
   result += (-0.5) * tmp_2689;
   std::complex<double> tmp_2691;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2692;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2692 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_2691 += tmp_2692;
   }
   result += tmp_2691;
   std::complex<double> tmp_2693;
   std::complex<double> tmp_2694;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2695;
      std::complex<double> tmp_2696;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2696 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2695 += tmp_2696;
      tmp_2694 += (MFu(gI1)) * tmp_2695;
   }
   tmp_2693 += tmp_2694;
   result += (-2) * tmp_2693;
   std::complex<double> tmp_2697;
   std::complex<double> tmp_2698;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2698 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2697 += tmp_2698;
   result += (-1) * tmp_2697;
   std::complex<double> tmp_2699;
   std::complex<double> tmp_2700;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2700 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2699 += tmp_2700;
   result += (-1) * tmp_2699;
   std::complex<double> tmp_2701;
   std::complex<double> tmp_2702;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2702 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2701 += tmp_2702;
   result += (-1) * tmp_2701;
   std::complex<double> tmp_2703;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2704;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2704 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_2703 += tmp_2704;
   }
   result += tmp_2703;
   std::complex<double> tmp_2705;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2706;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2706 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_2705 += tmp_2706;
   }
   result += tmp_2705;
   std::complex<double> tmp_2707;
   std::complex<double> tmp_2708;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2708 += (Conj(CpconjUSuGluFuPL(gO2,1,gI2))*CpconjUSuGluFuPL(gO1,1,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,1,gI2))*CpconjUSuGluFuPR(gO1,1,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_2707 += tmp_2708;
   result += (1.3333333333333333) * tmp_2707;
   std::complex<double> tmp_2709;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2709 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_2709;
   std::complex<double> tmp_2710;
   std::complex<double> tmp_2711;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2711 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_2710 += tmp_2711;
   result += (1.3333333333333333) * tmp_2710;
   std::complex<double> tmp_2712;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2712 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_2712;
   std::complex<double> tmp_2713;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2713 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_2713;
   std::complex<double> tmp_2714;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2714 += Conj(CpconjUSuVZpSu(gO2,gI2))*CpconjUSuVZpSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVZp);
   }
   result += tmp_2714;
   std::complex<double> tmp_2715;
   std::complex<double> tmp_2716;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2716 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,1,gI2))*
         CpconjUSuGluFuPL(gO1,1,gI2) + Conj(CpconjUSuGluFuPL(gO2,1,gI2))*
         CpconjUSuGluFuPR(gO1,1,gI2))*MFu(gI2);
   }
   tmp_2715 += tmp_2716;
   result += (-2.6666666666666665*MGlu) * tmp_2715;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSeconjUSeVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_2717;
   std::complex<double> tmp_2718;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2718 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2717 += tmp_2718;
   result += (-1) * tmp_2717;
   std::complex<double> tmp_2719;
   std::complex<double> tmp_2720;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2720 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2719 += tmp_2720;
   result += (-0.5) * tmp_2719;
   std::complex<double> tmp_2721;
   std::complex<double> tmp_2722;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2722 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2721 += tmp_2722;
   result += (-1) * tmp_2721;
   std::complex<double> tmp_2723;
   std::complex<double> tmp_2724;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2724 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2723 += tmp_2724;
   result += (-0.5) * tmp_2723;
   std::complex<double> tmp_2725;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2726;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2726 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeSvHpm(gO2,gI1
            ,gI2))*CpconjUSeSvHpm(gO1,gI1,gI2);
      }
      tmp_2725 += tmp_2726;
   }
   result += tmp_2725;
   std::complex<double> tmp_2727;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2728;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2728 += (Conj(CpconjUSeFvChaPL(gO2,gI1,gI2))*
            CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPR(gO2,gI1,gI2))*
            CpconjUSeFvChaPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha(gI2));
      }
      tmp_2727 += tmp_2728;
   }
   result += tmp_2727;
   std::complex<double> tmp_2729;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2730;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2730 += (Conj(CpconjUSeFeChiPL(gO2,gI1,gI2))*
            CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPR(gO2,gI1,gI2))*
            CpconjUSeFeChiPR(gO1,gI1,gI2))*G0(p,MFe(gI1),MChi(gI2));
      }
      tmp_2729 += tmp_2730;
   }
   result += tmp_2729;
   std::complex<double> tmp_2731;
   std::complex<double> tmp_2732;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2733;
      std::complex<double> tmp_2734;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2734 += B0(p,MFe(gI1),MChi(gI2))*(Conj(CpconjUSeFeChiPR(gO2,
            gI1,gI2))*CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPL(gO2,
            gI1,gI2))*CpconjUSeFeChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2733 += tmp_2734;
      tmp_2732 += (MFe(gI1)) * tmp_2733;
   }
   tmp_2731 += tmp_2732;
   result += (-2) * tmp_2731;
   std::complex<double> tmp_2735;
   std::complex<double> tmp_2736;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2737;
      std::complex<double> tmp_2738;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2738 += B0(p,MFv(gI1),MCha(gI2))*(Conj(CpconjUSeFvChaPR(gO2,
            gI1,gI2))*CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPL(gO2,
            gI1,gI2))*CpconjUSeFvChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2737 += tmp_2738;
      tmp_2736 += (MFv(gI1)) * tmp_2737;
   }
   tmp_2735 += tmp_2736;
   result += (-2) * tmp_2735;
   std::complex<double> tmp_2739;
   std::complex<double> tmp_2740;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2740 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2739 += tmp_2740;
   result += (-3) * tmp_2739;
   std::complex<double> tmp_2741;
   std::complex<double> tmp_2742;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2742 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2741 += tmp_2742;
   result += (-1) * tmp_2741;
   std::complex<double> tmp_2743;
   std::complex<double> tmp_2744;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2744 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2743 += tmp_2744;
   result += (-3) * tmp_2743;
   std::complex<double> tmp_2745;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2746;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2746 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_2745 += tmp_2746;
   }
   result += tmp_2745;
   std::complex<double> tmp_2747;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2748;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2748 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_2747 += tmp_2748;
   }
   result += tmp_2747;
   std::complex<double> tmp_2749;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2749 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_2749;
   std::complex<double> tmp_2750;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2750 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_2750;
   std::complex<double> tmp_2751;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2751 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_2751;
   std::complex<double> tmp_2752;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2752 += Conj(CpconjUSeVZpSe(gO2,gI2))*CpconjUSeVZpSe(gO1,gI2)*F0(p
         ,MSe(gI2),MVZp);
   }
   result += tmp_2752;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmCgWmC(gO1)*CpUhhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmgWm(gO1)*CpUhhbargWmgWm(gO2));
   result += -(B0(p,MVZ,MVZ)*CpUhhbargZgZ(gO1)*CpUhhbargZgZ(gO2));
   result += -2*B0(p,MVZ,MVZp)*CpUhhbargZpgZ(gO1)*CpUhhbargZpgZ(gO2);
   result += -(B0(p,MVZp,MVZp)*CpUhhbargZpgZp(gO1)*CpUhhbargZpgZp(gO2));
   result += 4*B0(p,MVWm,MVWm)*Conj(CpUhhconjVWmVWm(gO2))*CpUhhconjVWmVWm(gO1);
   result += 4*A0(MVWm)*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUhhUhhVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUhhUhhVZVZ(gO1,gO2);
   result += 4*B0(p,MVZ,MVZp)*Conj(CpUhhVZpVZ(gO2))*CpUhhVZpVZ(gO1);
   result += 2*B0(p,MVZp,MVZp)*Conj(CpUhhVZpVZp(gO2))*CpUhhVZpVZp(gO1);
   result += 2*B0(p,MVZ,MVZ)*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
   std::complex<double> tmp_2753;
   std::complex<double> tmp_2754;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2754 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2753 += tmp_2754;
   result += (-1) * tmp_2753;
   std::complex<double> tmp_2755;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2756;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2756 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2755 += tmp_2756;
   }
   result += tmp_2755;
   std::complex<double> tmp_2757;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2758;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2758 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*
            CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*
            CpUhhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2757 += tmp_2758;
   }
   result += tmp_2757;
   std::complex<double> tmp_2759;
   std::complex<double> tmp_2760;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2761;
      std::complex<double> tmp_2762;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2762 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2
            ,gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2761 += tmp_2762;
      tmp_2760 += (MCha(gI1)) * tmp_2761;
   }
   tmp_2759 += tmp_2760;
   result += (-2) * tmp_2759;
   std::complex<double> tmp_2763;
   std::complex<double> tmp_2764;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2764 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2763 += tmp_2764;
   result += (-0.5) * tmp_2763;
   std::complex<double> tmp_2765;
   std::complex<double> tmp_2766;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2766 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2765 += tmp_2766;
   result += (-1) * tmp_2765;
   std::complex<double> tmp_2767;
   std::complex<double> tmp_2768;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2768 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2767 += tmp_2768;
   result += (-0.5) * tmp_2767;
   std::complex<double> tmp_2769;
   std::complex<double> tmp_2770;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2771;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2771 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_2770 += tmp_2771;
   }
   tmp_2769 += tmp_2770;
   result += (0.5) * tmp_2769;
   std::complex<double> tmp_2772;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2773;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2773 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2772 += tmp_2773;
   }
   result += tmp_2772;
   std::complex<double> tmp_2774;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2775;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2775 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))
            *CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_2774 += tmp_2775;
   }
   result += tmp_2774;
   std::complex<double> tmp_2776;
   std::complex<double> tmp_2777;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2778;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2778 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_2777 += tmp_2778;
   }
   tmp_2776 += tmp_2777;
   result += (0.5) * tmp_2776;
   std::complex<double> tmp_2779;
   std::complex<double> tmp_2780;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2781;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2781 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2780 += tmp_2781;
   }
   tmp_2779 += tmp_2780;
   result += (3) * tmp_2779;
   std::complex<double> tmp_2782;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2783;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2783 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2782 += tmp_2783;
   }
   result += tmp_2782;
   std::complex<double> tmp_2784;
   std::complex<double> tmp_2785;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2786;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2786 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2785 += tmp_2786;
   }
   tmp_2784 += tmp_2785;
   result += (3) * tmp_2784;
   std::complex<double> tmp_2787;
   std::complex<double> tmp_2788;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2789;
      std::complex<double> tmp_2790;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2790 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2789 += tmp_2790;
      tmp_2788 += (MFd(gI1)) * tmp_2789;
   }
   tmp_2787 += tmp_2788;
   result += (-6) * tmp_2787;
   std::complex<double> tmp_2791;
   std::complex<double> tmp_2792;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2793;
      std::complex<double> tmp_2794;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2794 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2793 += tmp_2794;
      tmp_2792 += (MFe(gI1)) * tmp_2793;
   }
   tmp_2791 += tmp_2792;
   result += (-2) * tmp_2791;
   std::complex<double> tmp_2795;
   std::complex<double> tmp_2796;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2797;
      std::complex<double> tmp_2798;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2798 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2797 += tmp_2798;
      tmp_2796 += (MFu(gI1)) * tmp_2797;
   }
   tmp_2795 += tmp_2796;
   result += (-6) * tmp_2795;
   std::complex<double> tmp_2799;
   std::complex<double> tmp_2800;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2800 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2799 += tmp_2800;
   result += (-3) * tmp_2799;
   std::complex<double> tmp_2801;
   std::complex<double> tmp_2802;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2802 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2801 += tmp_2802;
   result += (-1) * tmp_2801;
   std::complex<double> tmp_2803;
   std::complex<double> tmp_2804;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2804 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2803 += tmp_2804;
   result += (-3) * tmp_2803;
   std::complex<double> tmp_2805;
   std::complex<double> tmp_2806;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2807;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2807 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2806 += tmp_2807;
   }
   tmp_2805 += tmp_2806;
   result += (3) * tmp_2805;
   std::complex<double> tmp_2808;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2809;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2809 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2808 += tmp_2809;
   }
   result += tmp_2808;
   std::complex<double> tmp_2810;
   std::complex<double> tmp_2811;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2812;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2812 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2811 += tmp_2812;
   }
   tmp_2810 += tmp_2811;
   result += (3) * tmp_2810;
   std::complex<double> tmp_2813;
   std::complex<double> tmp_2814;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2815;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2815 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2814 += tmp_2815;
   }
   tmp_2813 += tmp_2814;
   result += (0.5) * tmp_2813;
   std::complex<double> tmp_2816;
   std::complex<double> tmp_2817;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2818;
      std::complex<double> tmp_2819;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2819 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,
            gI1,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2)
            )*CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2818 += tmp_2819;
      tmp_2817 += (MChi(gI1)) * tmp_2818;
   }
   tmp_2816 += tmp_2817;
   result += (-1) * tmp_2816;
   std::complex<double> tmp_2820;
   std::complex<double> tmp_2821;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2821 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2820 += tmp_2821;
   result += (2) * tmp_2820;
   std::complex<double> tmp_2822;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2822 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_2822;
   std::complex<double> tmp_2823;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2823 += Conj(CpUhhVZpAh(gO2,gI2))*CpUhhVZpAh(gO1,gI2)*F0(p,MAh(gI2
         ),MVZp);
   }
   result += tmp_2823;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   result += 4*A0(MVWm)*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUAhUAhVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUAhUAhVZVZ(gO1,gO2);
   std::complex<double> tmp_2824;
   std::complex<double> tmp_2825;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2825 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2824 += tmp_2825;
   result += (-1) * tmp_2824;
   std::complex<double> tmp_2826;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2827;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2827 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2826 += tmp_2827;
   }
   result += tmp_2826;
   std::complex<double> tmp_2828;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2829;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2829 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*
            CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*
            CpUAhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2828 += tmp_2829;
   }
   result += tmp_2828;
   std::complex<double> tmp_2830;
   std::complex<double> tmp_2831;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2832;
      std::complex<double> tmp_2833;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2833 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2
            ,gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2832 += tmp_2833;
      tmp_2831 += (MCha(gI1)) * tmp_2832;
   }
   tmp_2830 += tmp_2831;
   result += (-2) * tmp_2830;
   std::complex<double> tmp_2834;
   std::complex<double> tmp_2835;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2835 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2834 += tmp_2835;
   result += (-0.5) * tmp_2834;
   std::complex<double> tmp_2836;
   std::complex<double> tmp_2837;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2837 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2836 += tmp_2837;
   result += (-1) * tmp_2836;
   std::complex<double> tmp_2838;
   std::complex<double> tmp_2839;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2839 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2838 += tmp_2839;
   result += (-0.5) * tmp_2838;
   std::complex<double> tmp_2840;
   std::complex<double> tmp_2841;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2842;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2842 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))
            *CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_2841 += tmp_2842;
   }
   tmp_2840 += tmp_2841;
   result += (0.5) * tmp_2840;
   std::complex<double> tmp_2843;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2844;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2844 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_2843 += tmp_2844;
   }
   result += tmp_2843;
   std::complex<double> tmp_2845;
   std::complex<double> tmp_2846;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2847;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2847 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))
            *CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_2846 += tmp_2847;
   }
   tmp_2845 += tmp_2846;
   result += (0.5) * tmp_2845;
   std::complex<double> tmp_2848;
   std::complex<double> tmp_2849;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2850;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2850 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2849 += tmp_2850;
   }
   tmp_2848 += tmp_2849;
   result += (3) * tmp_2848;
   std::complex<double> tmp_2851;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2852;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2852 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2851 += tmp_2852;
   }
   result += tmp_2851;
   std::complex<double> tmp_2853;
   std::complex<double> tmp_2854;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2855;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2855 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2854 += tmp_2855;
   }
   tmp_2853 += tmp_2854;
   result += (3) * tmp_2853;
   std::complex<double> tmp_2856;
   std::complex<double> tmp_2857;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2858;
      std::complex<double> tmp_2859;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2859 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2858 += tmp_2859;
      tmp_2857 += (MFd(gI1)) * tmp_2858;
   }
   tmp_2856 += tmp_2857;
   result += (-6) * tmp_2856;
   std::complex<double> tmp_2860;
   std::complex<double> tmp_2861;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2862;
      std::complex<double> tmp_2863;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2863 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2862 += tmp_2863;
      tmp_2861 += (MFe(gI1)) * tmp_2862;
   }
   tmp_2860 += tmp_2861;
   result += (-2) * tmp_2860;
   std::complex<double> tmp_2864;
   std::complex<double> tmp_2865;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2866;
      std::complex<double> tmp_2867;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2867 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2866 += tmp_2867;
      tmp_2865 += (MFu(gI1)) * tmp_2866;
   }
   tmp_2864 += tmp_2865;
   result += (-6) * tmp_2864;
   std::complex<double> tmp_2868;
   std::complex<double> tmp_2869;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2869 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2868 += tmp_2869;
   result += (-3) * tmp_2868;
   std::complex<double> tmp_2870;
   std::complex<double> tmp_2871;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2871 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2870 += tmp_2871;
   result += (-1) * tmp_2870;
   std::complex<double> tmp_2872;
   std::complex<double> tmp_2873;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2873 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2872 += tmp_2873;
   result += (-3) * tmp_2872;
   std::complex<double> tmp_2874;
   std::complex<double> tmp_2875;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2876;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2876 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2875 += tmp_2876;
   }
   tmp_2874 += tmp_2875;
   result += (3) * tmp_2874;
   std::complex<double> tmp_2877;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2878;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2878 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2877 += tmp_2878;
   }
   result += tmp_2877;
   std::complex<double> tmp_2879;
   std::complex<double> tmp_2880;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2881;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2881 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2880 += tmp_2881;
   }
   tmp_2879 += tmp_2880;
   result += (3) * tmp_2879;
   std::complex<double> tmp_2882;
   std::complex<double> tmp_2883;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2884;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2884 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2883 += tmp_2884;
   }
   tmp_2882 += tmp_2883;
   result += (0.5) * tmp_2882;
   std::complex<double> tmp_2885;
   std::complex<double> tmp_2886;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2887;
      std::complex<double> tmp_2888;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2888 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,
            gI1,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2)
            )*CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2887 += tmp_2888;
      tmp_2886 += (MChi(gI1)) * tmp_2887;
   }
   tmp_2885 += tmp_2886;
   result += (-1) * tmp_2885;
   std::complex<double> tmp_2889;
   std::complex<double> tmp_2890;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2890 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2889 += tmp_2890;
   result += (2) * tmp_2889;
   std::complex<double> tmp_2891;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2891 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_2891;
   std::complex<double> tmp_2892;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2892 += Conj(CpUAhVZphh(gO2,gI2))*CpUAhVZphh(gO1,gI2)*F0(p,Mhh(gI2
         ),MVZp);
   }
   result += tmp_2892;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Hpm(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*B0(p,0,MVWm)*Conj(CpconjUHpmVWmVP(gO2))*CpconjUHpmVWmVP(gO1);
   result += 4*B0(p,MVWm,MVZp)*Conj(CpconjUHpmVZpVWm(gO2))*CpconjUHpmVZpVWm(gO1
      );
   result += 4*B0(p,MVWm,MVZ)*Conj(CpconjUHpmVZVWm(gO2))*CpconjUHpmVZVWm(gO1);
   result += 4*A0(MVWm)*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUHpmconjUHpmVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUHpmconjUHpmVZVZ(gO1,gO2);
   result += -(B0(p,MVZ,MVWm)*CpconjUHpmbargWmCgZ(gO1)*CpUHpmgWmCbargZ(gO2));
   result += -(B0(p,MVZp,MVWm)*CpconjUHpmbargWmCgZp(gO1)*CpUHpmgWmCbargZp(gO2))
      ;
   result += -(B0(p,MVWm,MVZ)*CpconjUHpmbargZgWm(gO1)*CpUHpmgZbargWm(gO2));
   result += -(B0(p,MVWm,MVZp)*CpconjUHpmbargZpgWm(gO1)*CpUHpmgZpbargWm(gO2));
   std::complex<double> tmp_2893;
   std::complex<double> tmp_2894;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2894 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2893 += tmp_2894;
   result += (-1) * tmp_2893;
   std::complex<double> tmp_2895;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2896;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2896 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2895 += tmp_2896;
   }
   result += tmp_2895;
   std::complex<double> tmp_2897;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2898;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2898 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2897 += tmp_2898;
   }
   result += tmp_2897;
   std::complex<double> tmp_2899;
   std::complex<double> tmp_2900;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2900 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2899 += tmp_2900;
   result += (-0.5) * tmp_2899;
   std::complex<double> tmp_2901;
   std::complex<double> tmp_2902;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2902 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2901 += tmp_2902;
   result += (-1) * tmp_2901;
   std::complex<double> tmp_2903;
   std::complex<double> tmp_2904;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2904 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2903 += tmp_2904;
   result += (-0.5) * tmp_2903;
   std::complex<double> tmp_2905;
   std::complex<double> tmp_2906;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2907;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2907 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_2906 += tmp_2907;
   }
   tmp_2905 += tmp_2906;
   result += (3) * tmp_2905;
   std::complex<double> tmp_2908;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2909;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2909 += (Conj(CpconjUHpmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFvFePR(gO2,gI1,
            gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_2908 += tmp_2909;
   }
   result += tmp_2908;
   std::complex<double> tmp_2910;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2911;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2911 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2910 += tmp_2911;
   }
   result += tmp_2910;
   std::complex<double> tmp_2912;
   std::complex<double> tmp_2913;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2914;
      std::complex<double> tmp_2915;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2915 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2914 += tmp_2915;
      tmp_2913 += (MFu(gI1)) * tmp_2914;
   }
   tmp_2912 += tmp_2913;
   result += (-6) * tmp_2912;
   std::complex<double> tmp_2916;
   std::complex<double> tmp_2917;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2918;
      std::complex<double> tmp_2919;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2919 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmbarFvFePR(
            gO2,gI1,gI2))*CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFvFePL(gO2,gI1,gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2918 += tmp_2919;
      tmp_2917 += (MFv(gI1)) * tmp_2918;
   }
   tmp_2916 += tmp_2917;
   result += (-2) * tmp_2916;
   std::complex<double> tmp_2920;
   std::complex<double> tmp_2921;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2921 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2920 += tmp_2921;
   result += (-3) * tmp_2920;
   std::complex<double> tmp_2922;
   std::complex<double> tmp_2923;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2923 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2922 += tmp_2923;
   result += (-1) * tmp_2922;
   std::complex<double> tmp_2924;
   std::complex<double> tmp_2925;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2925 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2924 += tmp_2925;
   result += (-3) * tmp_2924;
   std::complex<double> tmp_2926;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2927;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2927 += (Conj(CpconjUHpmChiChaPL(gO2,gI1,gI2))*
            CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHpmChiChaPR(gO2,gI1,gI2)
            )*CpconjUHpmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_2926 += tmp_2927;
   }
   result += tmp_2926;
   std::complex<double> tmp_2928;
   std::complex<double> tmp_2929;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2930;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2930 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2929 += tmp_2930;
   }
   tmp_2928 += tmp_2929;
   result += (3) * tmp_2928;
   std::complex<double> tmp_2931;
   std::complex<double> tmp_2932;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2933;
      std::complex<double> tmp_2934;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2934 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHpmChiChaPR(
            gO2,gI1,gI2))*CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmChiChaPL(gO2,gI1,gI2))*CpconjUHpmChiChaPR(gO1,gI1,gI2))*MCha
            (gI2);
      }
      tmp_2933 += tmp_2934;
      tmp_2932 += (MChi(gI1)) * tmp_2933;
   }
   tmp_2931 += tmp_2932;
   result += (-2) * tmp_2931;
   std::complex<double> tmp_2935;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2935 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_2935;
   std::complex<double> tmp_2936;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2936 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_2936;
   std::complex<double> tmp_2937;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2937 += Conj(CpconjUHpmVZpHpm(gO2,gI2))*CpconjUHpmVZpHpm(gO1,gI2)*
         F0(p,MHpm(gI2),MVZp);
   }
   result += tmp_2937;
   std::complex<double> tmp_2938;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2938 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_2938;
   std::complex<double> tmp_2939;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2939 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_2939;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_2940;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2940 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_2940;
   std::complex<double> tmp_2941;
   std::complex<double> tmp_2942;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2943;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2943 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2942 += tmp_2943;
   }
   tmp_2941 += tmp_2942;
   result += (-4) * tmp_2941;
   std::complex<double> tmp_2944;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2945;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2945 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2945 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_2944 += tmp_2945;
   }
   result += tmp_2944;
   std::complex<double> tmp_2946;
   std::complex<double> tmp_2947;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2947 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_2946 += tmp_2947;
   result += (0.5) * tmp_2946;
   std::complex<double> tmp_2948;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2948 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2948;
   std::complex<double> tmp_2949;
   std::complex<double> tmp_2950;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2950 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2949 += tmp_2950;
   result += (0.5) * tmp_2949;
   std::complex<double> tmp_2951;
   std::complex<double> tmp_2952;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2953;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2953 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2952 += tmp_2953;
   }
   tmp_2951 += tmp_2952;
   result += (-4) * tmp_2951;
   std::complex<double> tmp_2954;
   std::complex<double> tmp_2955;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2956;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2956 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_2955 += tmp_2956;
   }
   tmp_2954 += tmp_2955;
   result += (-4) * tmp_2954;
   std::complex<double> tmp_2957;
   std::complex<double> tmp_2958;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2959;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2959 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2959 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_2958 += tmp_2959;
   }
   tmp_2957 += tmp_2958;
   result += (3) * tmp_2957;
   std::complex<double> tmp_2960;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2961;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2961 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2961 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_2960 += tmp_2961;
   }
   result += tmp_2960;
   std::complex<double> tmp_2962;
   std::complex<double> tmp_2963;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2964;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2964 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2964 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_2963 += tmp_2964;
   }
   tmp_2962 += tmp_2963;
   result += (3) * tmp_2962;
   std::complex<double> tmp_2965;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2966;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2966 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_2966 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_2965 += tmp_2966;
   }
   result += tmp_2965;
   std::complex<double> tmp_2967;
   std::complex<double> tmp_2968;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2968 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2967 += tmp_2968;
   result += (3) * tmp_2967;
   std::complex<double> tmp_2969;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2969 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2969;
   std::complex<double> tmp_2970;
   std::complex<double> tmp_2971;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2971 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2970 += tmp_2971;
   result += (3) * tmp_2970;
   std::complex<double> tmp_2972;
   std::complex<double> tmp_2973;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2974;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2974 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2973 += tmp_2974;
   }
   tmp_2972 += tmp_2973;
   result += (-12) * tmp_2972;
   std::complex<double> tmp_2975;
   std::complex<double> tmp_2976;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2977;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2977 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2976 += tmp_2977;
   }
   tmp_2975 += tmp_2976;
   result += (-4) * tmp_2975;
   std::complex<double> tmp_2978;
   std::complex<double> tmp_2979;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2980;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2980 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2979 += tmp_2980;
   }
   tmp_2978 += tmp_2979;
   result += (-12) * tmp_2978;
   std::complex<double> tmp_2981;
   std::complex<double> tmp_2982;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2983;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2983 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2983 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_2982 += tmp_2983;
   }
   tmp_2981 += tmp_2982;
   result += (0.5) * tmp_2981;
   std::complex<double> tmp_2984;
   std::complex<double> tmp_2985;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2985 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2984 += tmp_2985;
   result += (2) * tmp_2984;
   std::complex<double> tmp_2986;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2986 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2986;
   std::complex<double> tmp_2987;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2987 += AbsSqr(CpVZVZphh(gI2))*B0(p,MVZp,Mhh(gI2));
   }
   result += tmp_2987;
   result += -(AbsSqr(CpVZconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) + B0(
      p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZp(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZpbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZpbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZpVZpconjVWmVWm1() + CpVZpVZpconjVWmVWm2() +
      CpVZpVZpconjVWmVWm3()));
   std::complex<double> tmp_2988;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2988 += A0(MHpm(gI1))*CpVZpVZpconjHpmHpm(gI1,gI1);
   }
   result += tmp_2988;
   std::complex<double> tmp_2989;
   std::complex<double> tmp_2990;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2991;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2991 += AbsSqr(CpVZpconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),
            MHpm(gI2));
      }
      tmp_2990 += tmp_2991;
   }
   tmp_2989 += tmp_2990;
   result += (-4) * tmp_2989;
   std::complex<double> tmp_2992;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2993;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2993 += (AbsSqr(CpVZpbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZpbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2993 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZpbarChaChaPL(gI1,gI2))*CpVZpbarChaChaPR(gI1,gI2));
      }
      tmp_2992 += tmp_2993;
   }
   result += tmp_2992;
   std::complex<double> tmp_2994;
   std::complex<double> tmp_2995;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2995 += A0(MAh(gI1))*CpVZpVZpAhAh(gI1,gI1);
   }
   tmp_2994 += tmp_2995;
   result += (0.5) * tmp_2994;
   std::complex<double> tmp_2996;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2996 += A0(MSv(gI1))*CpVZpVZpconjSvSv(gI1,gI1);
   }
   result += tmp_2996;
   std::complex<double> tmp_2997;
   std::complex<double> tmp_2998;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2998 += A0(Mhh(gI1))*CpVZpVZphhhh(gI1,gI1);
   }
   tmp_2997 += tmp_2998;
   result += (0.5) * tmp_2997;
   std::complex<double> tmp_2999;
   std::complex<double> tmp_3000;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3001;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3001 += AbsSqr(CpVZpconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(
            gI2));
      }
      tmp_3000 += tmp_3001;
   }
   tmp_2999 += tmp_3000;
   result += (-4) * tmp_2999;
   std::complex<double> tmp_3002;
   std::complex<double> tmp_3003;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3004;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3004 += AbsSqr(CpVZphhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_3003 += tmp_3004;
   }
   tmp_3002 += tmp_3003;
   result += (-4) * tmp_3002;
   std::complex<double> tmp_3005;
   std::complex<double> tmp_3006;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3007;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3007 += (AbsSqr(CpVZpbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZpbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_3007 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZpbarFdFdPL(gI1,gI2))*CpVZpbarFdFdPR(gI1,gI2));
      }
      tmp_3006 += tmp_3007;
   }
   tmp_3005 += tmp_3006;
   result += (3) * tmp_3005;
   std::complex<double> tmp_3008;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3009;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3009 += (AbsSqr(CpVZpbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZpbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_3009 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZpbarFeFePL(gI1,gI2))*CpVZpbarFeFePR(gI1,gI2));
      }
      tmp_3008 += tmp_3009;
   }
   result += tmp_3008;
   std::complex<double> tmp_3010;
   std::complex<double> tmp_3011;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3012;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3012 += (AbsSqr(CpVZpbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZpbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_3012 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZpbarFuFuPL(gI1,gI2))*CpVZpbarFuFuPR(gI1,gI2));
      }
      tmp_3011 += tmp_3012;
   }
   tmp_3010 += tmp_3011;
   result += (3) * tmp_3010;
   std::complex<double> tmp_3013;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3014;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3014 += (AbsSqr(CpVZpbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZpbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_3014 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZpbarFvFvPL(gI1,gI2))*CpVZpbarFvFvPR(gI1,gI2));
      }
      tmp_3013 += tmp_3014;
   }
   result += tmp_3013;
   std::complex<double> tmp_3015;
   std::complex<double> tmp_3016;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3016 += A0(MSd(gI1))*CpVZpVZpconjSdSd(gI1,gI1);
   }
   tmp_3015 += tmp_3016;
   result += (3) * tmp_3015;
   std::complex<double> tmp_3017;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3017 += A0(MSe(gI1))*CpVZpVZpconjSeSe(gI1,gI1);
   }
   result += tmp_3017;
   std::complex<double> tmp_3018;
   std::complex<double> tmp_3019;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3019 += A0(MSu(gI1))*CpVZpVZpconjSuSu(gI1,gI1);
   }
   tmp_3018 += tmp_3019;
   result += (3) * tmp_3018;
   std::complex<double> tmp_3020;
   std::complex<double> tmp_3021;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3022;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3022 += AbsSqr(CpVZpconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(
            gI2));
      }
      tmp_3021 += tmp_3022;
   }
   tmp_3020 += tmp_3021;
   result += (-12) * tmp_3020;
   std::complex<double> tmp_3023;
   std::complex<double> tmp_3024;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3025;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3025 += AbsSqr(CpVZpconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(
            gI2));
      }
      tmp_3024 += tmp_3025;
   }
   tmp_3023 += tmp_3024;
   result += (-4) * tmp_3023;
   std::complex<double> tmp_3026;
   std::complex<double> tmp_3027;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3028;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3028 += AbsSqr(CpVZpconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(
            gI2));
      }
      tmp_3027 += tmp_3028;
   }
   tmp_3026 += tmp_3027;
   result += (-12) * tmp_3026;
   std::complex<double> tmp_3029;
   std::complex<double> tmp_3030;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3031;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3031 += (AbsSqr(CpVZpChiChiPL(gI1,gI2)) + AbsSqr(
            CpVZpChiChiPR(gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_3031 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZpChiChiPL(gI1,gI2))*CpVZpChiChiPR(gI1,gI2));
      }
      tmp_3030 += tmp_3031;
   }
   tmp_3029 += tmp_3030;
   result += (0.5) * tmp_3029;
   std::complex<double> tmp_3032;
   std::complex<double> tmp_3033;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3033 += AbsSqr(CpVZpconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_3032 += tmp_3033;
   result += (2) * tmp_3032;
   std::complex<double> tmp_3034;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3034 += AbsSqr(CpVZpVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_3034;
   std::complex<double> tmp_3035;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3035 += AbsSqr(CpVZpVZphh(gI2))*B0(p,MVZp,Mhh(gI2));
   }
   result += tmp_3035;
   result += -(AbsSqr(CpVZpconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) + B0
      (p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWmbargPgWm())*B00(p,MVWm,MVP);
   result += AbsSqr(CpconjVWmbargWmCgP())*B00(p,MVP,MVWm);
   result += AbsSqr(CpconjVWmbargWmCgZ())*B00(p,MVZ,MVWm);
   result += AbsSqr(CpconjVWmbargWmCgZp())*B00(p,MVZp,MVWm);
   result += AbsSqr(CpconjVWmbargZgWm())*B00(p,MVWm,MVZ);
   result += AbsSqr(CpconjVWmbargZpgWm())*B00(p,MVWm,MVZp);
   result += -(A0(MVWm)*(4*CpVWmconjVWmconjVWmVWm1() + CpVWmconjVWmconjVWmVWm2(
      ) + CpVWmconjVWmconjVWmVWm3()));
   result += 0;
   result += -0.5*A0(MVZp)*(4*CpVWmconjVWmVZpVZp1() + CpVWmconjVWmVZpVZp2() +
      CpVWmconjVWmVZpVZp3());
   result += -0.5*A0(MVZ)*(4*CpVWmconjVWmVZVZ1() + CpVWmconjVWmVZVZ2() +
      CpVWmconjVWmVZVZ3());
   std::complex<double> tmp_3036;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3036 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_3036;
   std::complex<double> tmp_3037;
   std::complex<double> tmp_3038;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3039;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3039 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_3038 += tmp_3039;
   }
   tmp_3037 += tmp_3038;
   result += (-4) * tmp_3037;
   std::complex<double> tmp_3040;
   std::complex<double> tmp_3041;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3042;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3042 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_3041 += tmp_3042;
   }
   tmp_3040 += tmp_3041;
   result += (-4) * tmp_3040;
   std::complex<double> tmp_3043;
   std::complex<double> tmp_3044;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3044 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_3043 += tmp_3044;
   result += (0.5) * tmp_3043;
   std::complex<double> tmp_3045;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3045 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_3045;
   std::complex<double> tmp_3046;
   std::complex<double> tmp_3047;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3047 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_3046 += tmp_3047;
   result += (0.5) * tmp_3046;
   std::complex<double> tmp_3048;
   std::complex<double> tmp_3049;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3050;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3050 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_3050 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_3049 += tmp_3050;
   }
   tmp_3048 += tmp_3049;
   result += (3) * tmp_3048;
   std::complex<double> tmp_3051;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3052;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3052 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_3052 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_3051 += tmp_3052;
   }
   result += tmp_3051;
   std::complex<double> tmp_3053;
   std::complex<double> tmp_3054;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3055;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3055 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_3054 += tmp_3055;
   }
   tmp_3053 += tmp_3054;
   result += (-4) * tmp_3053;
   std::complex<double> tmp_3056;
   std::complex<double> tmp_3057;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3057 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_3056 += tmp_3057;
   result += (3) * tmp_3056;
   std::complex<double> tmp_3058;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3058 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_3058;
   std::complex<double> tmp_3059;
   std::complex<double> tmp_3060;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3060 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_3059 += tmp_3060;
   result += (3) * tmp_3059;
   std::complex<double> tmp_3061;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3062;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3062 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_3062 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_3061 += tmp_3062;
   }
   result += tmp_3061;
   std::complex<double> tmp_3063;
   std::complex<double> tmp_3064;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3065;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3065 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_3064 += tmp_3065;
   }
   tmp_3063 += tmp_3064;
   result += (-12) * tmp_3063;
   std::complex<double> tmp_3066;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3066 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_3066;
   std::complex<double> tmp_3067;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3067 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_3067;
   std::complex<double> tmp_3068;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3068 += AbsSqr(CpconjVWmVZpHpm(gI2))*B0(p,MVZp,MHpm(gI2));
   }
   result += tmp_3068;
   std::complex<double> tmp_3069;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3069 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_3069;
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) + B0(p,
      MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) + B0(p,MVZ,MVWm)*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZpVWm())*(A0(MVWm) + A0(MVZp) + 10*B00(p,MVZp,
      MVWm) + B0(p,MVZp,MVWm)*(Sqr(MVWm) + Sqr(MVZp) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3070;
   std::complex<double> tmp_3071;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3072;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3072 += B0(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3071 += tmp_3072;
   }
   tmp_3070 += tmp_3071;
   result += (2) * tmp_3070;
   std::complex<double> tmp_3073;
   std::complex<double> tmp_3074;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3075;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3075 += B0(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3074 += tmp_3075;
   }
   tmp_3073 += tmp_3074;
   result += (2) * tmp_3073;
   std::complex<double> tmp_3076;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3077;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3077 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3076 += tmp_3077;
   }
   result += tmp_3076;
   std::complex<double> tmp_3078;
   std::complex<double> tmp_3079;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3080;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3080 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3079 += tmp_3080;
   }
   tmp_3078 += tmp_3079;
   result += (6) * tmp_3078;
   std::complex<double> tmp_3081;
   std::complex<double> tmp_3082;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3083;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3083 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3082 += tmp_3083;
   }
   tmp_3081 += tmp_3082;
   result += (2) * tmp_3081;
   std::complex<double> tmp_3084;
   std::complex<double> tmp_3085;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3086;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3086 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3085 += tmp_3086;
   }
   tmp_3084 += tmp_3085;
   result += (6) * tmp_3084;
   std::complex<double> tmp_3087;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3088;
      std::complex<double> tmp_3089;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3089 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_3088 += tmp_3089;
      tmp_3087 += (MChi(gI1)) * tmp_3088;
   }
   result += tmp_3087;
   std::complex<double> tmp_3090;
   std::complex<double> tmp_3091;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3091 += B0(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_3090 += tmp_3091;
   result += (-8) * tmp_3090;
   std::complex<double> tmp_3092;
   std::complex<double> tmp_3093;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3093 += B0(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_3092 += tmp_3093;
   result += (-4) * tmp_3092;
   std::complex<double> tmp_3094;
   std::complex<double> tmp_3095;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3095 += B0(p,MChi(gI2),MVZp)*Conj(CpUChiVZpChiPR(gO2,gI2))*
         CpUChiVZpChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_3094 += tmp_3095;
   result += (-4) * tmp_3094;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3096;
   std::complex<double> tmp_3097;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3098;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3098 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPR(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2);
      }
      tmp_3097 += tmp_3098;
   }
   tmp_3096 += tmp_3097;
   result += (-1) * tmp_3096;
   std::complex<double> tmp_3099;
   std::complex<double> tmp_3100;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3101;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3101 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPR(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2);
      }
      tmp_3100 += tmp_3101;
   }
   tmp_3099 += tmp_3100;
   result += (-1) * tmp_3099;
   std::complex<double> tmp_3102;
   std::complex<double> tmp_3103;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3104;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3104 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_3103 += tmp_3104;
   }
   tmp_3102 += tmp_3103;
   result += (-0.5) * tmp_3102;
   std::complex<double> tmp_3105;
   std::complex<double> tmp_3106;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3107;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3107 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_3106 += tmp_3107;
   }
   tmp_3105 += tmp_3106;
   result += (-0.5) * tmp_3105;
   std::complex<double> tmp_3108;
   std::complex<double> tmp_3109;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3110;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3110 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPR(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_3109 += tmp_3110;
   }
   tmp_3108 += tmp_3109;
   result += (-3) * tmp_3108;
   std::complex<double> tmp_3111;
   std::complex<double> tmp_3112;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3113;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3113 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePR(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_3112 += tmp_3113;
   }
   tmp_3111 += tmp_3112;
   result += (-1) * tmp_3111;
   std::complex<double> tmp_3114;
   std::complex<double> tmp_3115;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3116;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3116 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPR(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_3115 += tmp_3116;
   }
   tmp_3114 += tmp_3115;
   result += (-3) * tmp_3114;
   std::complex<double> tmp_3117;
   std::complex<double> tmp_3118;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3118 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPL(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_3117 += tmp_3118;
   result += (-2) * tmp_3117;
   std::complex<double> tmp_3119;
   std::complex<double> tmp_3120;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3120 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_3119 += tmp_3120;
   result += (-1) * tmp_3119;
   std::complex<double> tmp_3121;
   std::complex<double> tmp_3122;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3122 += B1(p,MChi(gI2),MVZp)*Conj(CpUChiVZpChiPL(gO2,gI2))*
         CpUChiVZpChiPL(gO1,gI2);
   }
   tmp_3121 += tmp_3122;
   result += (-1) * tmp_3121;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3123;
   std::complex<double> tmp_3124;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3125;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3125 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPL(gO1,gI1,gI2);
      }
      tmp_3124 += tmp_3125;
   }
   tmp_3123 += tmp_3124;
   result += (-1) * tmp_3123;
   std::complex<double> tmp_3126;
   std::complex<double> tmp_3127;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3128;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3128 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPL(gO1,gI1,gI2);
      }
      tmp_3127 += tmp_3128;
   }
   tmp_3126 += tmp_3127;
   result += (-1) * tmp_3126;
   std::complex<double> tmp_3129;
   std::complex<double> tmp_3130;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3131;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3131 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_3130 += tmp_3131;
   }
   tmp_3129 += tmp_3130;
   result += (-0.5) * tmp_3129;
   std::complex<double> tmp_3132;
   std::complex<double> tmp_3133;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3134;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3134 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_3133 += tmp_3134;
   }
   tmp_3132 += tmp_3133;
   result += (-0.5) * tmp_3132;
   std::complex<double> tmp_3135;
   std::complex<double> tmp_3136;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3137;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3137 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_3136 += tmp_3137;
   }
   tmp_3135 += tmp_3136;
   result += (-3) * tmp_3135;
   std::complex<double> tmp_3138;
   std::complex<double> tmp_3139;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3140;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3140 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_3139 += tmp_3140;
   }
   tmp_3138 += tmp_3139;
   result += (-1) * tmp_3138;
   std::complex<double> tmp_3141;
   std::complex<double> tmp_3142;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3143;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3143 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_3142 += tmp_3143;
   }
   tmp_3141 += tmp_3142;
   result += (-3) * tmp_3141;
   std::complex<double> tmp_3144;
   std::complex<double> tmp_3145;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3145 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_3144 += tmp_3145;
   result += (-2) * tmp_3144;
   std::complex<double> tmp_3146;
   std::complex<double> tmp_3147;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3147 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_3146 += tmp_3147;
   result += (-1) * tmp_3146;
   std::complex<double> tmp_3148;
   std::complex<double> tmp_3149;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3149 += B1(p,MChi(gI2),MVZp)*Conj(CpUChiVZpChiPR(gO2,gI2))*
         CpUChiVZpChiPR(gO1,gI2);
   }
   tmp_3148 += tmp_3149;
   result += (-1) * tmp_3148;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3150;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3151;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3151 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3150 += tmp_3151;
   }
   result += tmp_3150;
   std::complex<double> tmp_3152;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3153;
      std::complex<double> tmp_3154;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3154 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_3153 += tmp_3154;
      tmp_3152 += (MCha(gI1)) * tmp_3153;
   }
   result += tmp_3152;
   std::complex<double> tmp_3155;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3156;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3156 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3155 += tmp_3156;
   }
   result += tmp_3155;
   std::complex<double> tmp_3157;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3158;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3158 += B0(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3157 += tmp_3158;
   }
   result += tmp_3157;
   std::complex<double> tmp_3159;
   std::complex<double> tmp_3160;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3161;
      std::complex<double> tmp_3162;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3162 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_3161 += tmp_3162;
      tmp_3160 += (MFu(gI1)) * tmp_3161;
   }
   tmp_3159 += tmp_3160;
   result += (3) * tmp_3159;
   std::complex<double> tmp_3163;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3164;
      std::complex<double> tmp_3165;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3165 += B0(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_3164 += tmp_3165;
      tmp_3163 += (MFv(gI1)) * tmp_3164;
   }
   result += tmp_3163;
   std::complex<double> tmp_3166;
   std::complex<double> tmp_3167;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3168;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3168 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3167 += tmp_3168;
   }
   tmp_3166 += tmp_3167;
   result += (3) * tmp_3166;
   std::complex<double> tmp_3169;
   std::complex<double> tmp_3170;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3170 += B0(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_3169 += tmp_3170;
   result += (-4) * tmp_3169;
   std::complex<double> tmp_3171;
   std::complex<double> tmp_3172;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3172 += B0(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_3171 += tmp_3172;
   result += (-4) * tmp_3171;
   std::complex<double> tmp_3173;
   std::complex<double> tmp_3174;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3174 += B0(p,MCha(gI2),MVZp)*Conj(CpbarUChaVZpChaPR(gO2,gI2))*
         CpbarUChaVZpChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_3173 += tmp_3174;
   result += (-4) * tmp_3173;
   std::complex<double> tmp_3175;
   std::complex<double> tmp_3176;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3176 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_3175 += tmp_3176;
   result += (-4) * tmp_3175;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3177;
   std::complex<double> tmp_3178;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3179;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3179 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_3178 += tmp_3179;
   }
   tmp_3177 += tmp_3178;
   result += (-0.5) * tmp_3177;
   std::complex<double> tmp_3180;
   std::complex<double> tmp_3181;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3182;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3182 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPR(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2);
      }
      tmp_3181 += tmp_3182;
   }
   tmp_3180 += tmp_3181;
   result += (-0.5) * tmp_3180;
   std::complex<double> tmp_3183;
   std::complex<double> tmp_3184;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3185;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3185 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_3184 += tmp_3185;
   }
   tmp_3183 += tmp_3184;
   result += (-0.5) * tmp_3183;
   std::complex<double> tmp_3186;
   std::complex<double> tmp_3187;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3188;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3188 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePR(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2);
      }
      tmp_3187 += tmp_3188;
   }
   tmp_3186 += tmp_3187;
   result += (-0.5) * tmp_3186;
   std::complex<double> tmp_3189;
   std::complex<double> tmp_3190;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3191;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3191 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPR(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_3190 += tmp_3191;
   }
   tmp_3189 += tmp_3190;
   result += (-1.5) * tmp_3189;
   std::complex<double> tmp_3192;
   std::complex<double> tmp_3193;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3194;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3194 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePR(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_3193 += tmp_3194;
   }
   tmp_3192 += tmp_3193;
   result += (-0.5) * tmp_3192;
   std::complex<double> tmp_3195;
   std::complex<double> tmp_3196;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3197;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3197 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPR(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2);
      }
      tmp_3196 += tmp_3197;
   }
   tmp_3195 += tmp_3196;
   result += (-1.5) * tmp_3195;
   std::complex<double> tmp_3198;
   std::complex<double> tmp_3199;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3199 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_3198 += tmp_3199;
   result += (-1) * tmp_3198;
   std::complex<double> tmp_3200;
   std::complex<double> tmp_3201;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3201 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPL(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_3200 += tmp_3201;
   result += (-1) * tmp_3200;
   std::complex<double> tmp_3202;
   std::complex<double> tmp_3203;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3203 += B1(p,MCha(gI2),MVZp)*Conj(CpbarUChaVZpChaPL(gO2,gI2))*
         CpbarUChaVZpChaPL(gO1,gI2);
   }
   tmp_3202 += tmp_3203;
   result += (-1) * tmp_3202;
   std::complex<double> tmp_3204;
   std::complex<double> tmp_3205;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3205 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPL(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_3204 += tmp_3205;
   result += (-1) * tmp_3204;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3206;
   std::complex<double> tmp_3207;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3208;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3208 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_3207 += tmp_3208;
   }
   tmp_3206 += tmp_3207;
   result += (-0.5) * tmp_3206;
   std::complex<double> tmp_3209;
   std::complex<double> tmp_3210;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3211;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3211 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPL(gO1,gI1,gI2);
      }
      tmp_3210 += tmp_3211;
   }
   tmp_3209 += tmp_3210;
   result += (-0.5) * tmp_3209;
   std::complex<double> tmp_3212;
   std::complex<double> tmp_3213;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3214;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3214 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_3213 += tmp_3214;
   }
   tmp_3212 += tmp_3213;
   result += (-0.5) * tmp_3212;
   std::complex<double> tmp_3215;
   std::complex<double> tmp_3216;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3217;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3217 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePL(gO1,gI1,gI2);
      }
      tmp_3216 += tmp_3217;
   }
   tmp_3215 += tmp_3216;
   result += (-0.5) * tmp_3215;
   std::complex<double> tmp_3218;
   std::complex<double> tmp_3219;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3220;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3220 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2);
      }
      tmp_3219 += tmp_3220;
   }
   tmp_3218 += tmp_3219;
   result += (-1.5) * tmp_3218;
   std::complex<double> tmp_3221;
   std::complex<double> tmp_3222;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3223;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3223 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2);
      }
      tmp_3222 += tmp_3223;
   }
   tmp_3221 += tmp_3222;
   result += (-0.5) * tmp_3221;
   std::complex<double> tmp_3224;
   std::complex<double> tmp_3225;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3226;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3226 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPL(gO1,gI1,gI2);
      }
      tmp_3225 += tmp_3226;
   }
   tmp_3224 += tmp_3225;
   result += (-1.5) * tmp_3224;
   std::complex<double> tmp_3227;
   std::complex<double> tmp_3228;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3228 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_3227 += tmp_3228;
   result += (-1) * tmp_3227;
   std::complex<double> tmp_3229;
   std::complex<double> tmp_3230;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3230 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_3229 += tmp_3230;
   result += (-1) * tmp_3229;
   std::complex<double> tmp_3231;
   std::complex<double> tmp_3232;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3232 += B1(p,MCha(gI2),MVZp)*Conj(CpbarUChaVZpChaPR(gO2,gI2))*
         CpbarUChaVZpChaPR(gO1,gI2);
   }
   tmp_3231 += tmp_3232;
   result += (-1) * tmp_3231;
   std::complex<double> tmp_3233;
   std::complex<double> tmp_3234;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3234 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_3233 += tmp_3234;
   result += (-1) * tmp_3233;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3235;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3236;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3236 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3235 += tmp_3236;
   }
   result += tmp_3235;
   std::complex<double> tmp_3237;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3238;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3238 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3237 += tmp_3238;
   }
   result += tmp_3237;
   std::complex<double> tmp_3239;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3240;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3240 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3239 += tmp_3240;
   }
   result += tmp_3239;
   std::complex<double> tmp_3241;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3242;
      std::complex<double> tmp_3243;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3243 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3242 += tmp_3243;
      tmp_3241 += (MFe(gI1)) * tmp_3242;
   }
   result += tmp_3241;
   std::complex<double> tmp_3244;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3245;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3245 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3244 += tmp_3245;
   }
   result += tmp_3244;
   std::complex<double> tmp_3246;
   std::complex<double> tmp_3247;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3247 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3246 += tmp_3247;
   result += (-4) * tmp_3246;
   std::complex<double> tmp_3248;
   std::complex<double> tmp_3249;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3249 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3248 += tmp_3249;
   result += (-4) * tmp_3248;
   std::complex<double> tmp_3250;
   std::complex<double> tmp_3251;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3251 += B0(p,MFe(gI2),MVZp)*Conj(CpbarUFeVZpFePR(gO2,gI2))*
         CpbarUFeVZpFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3250 += tmp_3251;
   result += (-4) * tmp_3250;
   std::complex<double> tmp_3252;
   std::complex<double> tmp_3253;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3253 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3252 += tmp_3253;
   result += (-4) * tmp_3252;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3254;
   std::complex<double> tmp_3255;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3256;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3256 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3255 += tmp_3256;
   }
   tmp_3254 += tmp_3255;
   result += (-0.5) * tmp_3254;
   std::complex<double> tmp_3257;
   std::complex<double> tmp_3258;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3259;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3259 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPR(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_3258 += tmp_3259;
   }
   tmp_3257 += tmp_3258;
   result += (-0.5) * tmp_3257;
   std::complex<double> tmp_3260;
   std::complex<double> tmp_3261;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3262;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3262 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3261 += tmp_3262;
   }
   tmp_3260 += tmp_3261;
   result += (-0.5) * tmp_3260;
   std::complex<double> tmp_3263;
   std::complex<double> tmp_3264;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3265;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3265 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_3264 += tmp_3265;
   }
   tmp_3263 += tmp_3264;
   result += (-0.5) * tmp_3263;
   std::complex<double> tmp_3266;
   std::complex<double> tmp_3267;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3268;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3268 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3267 += tmp_3268;
   }
   tmp_3266 += tmp_3267;
   result += (-0.5) * tmp_3266;
   std::complex<double> tmp_3269;
   std::complex<double> tmp_3270;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3270 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_3269 += tmp_3270;
   result += (-1) * tmp_3269;
   std::complex<double> tmp_3271;
   std::complex<double> tmp_3272;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3272 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_3271 += tmp_3272;
   result += (-1) * tmp_3271;
   std::complex<double> tmp_3273;
   std::complex<double> tmp_3274;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3274 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_3273 += tmp_3274;
   result += (-1) * tmp_3273;
   std::complex<double> tmp_3275;
   std::complex<double> tmp_3276;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3276 += B1(p,MFe(gI2),MVZp)*Conj(CpbarUFeVZpFePL(gO2,gI2))*
         CpbarUFeVZpFePL(gO1,gI2);
   }
   tmp_3275 += tmp_3276;
   result += (-1) * tmp_3275;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3277;
   std::complex<double> tmp_3278;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3279;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3279 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3278 += tmp_3279;
   }
   tmp_3277 += tmp_3278;
   result += (-0.5) * tmp_3277;
   std::complex<double> tmp_3280;
   std::complex<double> tmp_3281;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3282;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3282 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_3281 += tmp_3282;
   }
   tmp_3280 += tmp_3281;
   result += (-0.5) * tmp_3280;
   std::complex<double> tmp_3283;
   std::complex<double> tmp_3284;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3285;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3285 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3284 += tmp_3285;
   }
   tmp_3283 += tmp_3284;
   result += (-0.5) * tmp_3283;
   std::complex<double> tmp_3286;
   std::complex<double> tmp_3287;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3288;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3288 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_3287 += tmp_3288;
   }
   tmp_3286 += tmp_3287;
   result += (-0.5) * tmp_3286;
   std::complex<double> tmp_3289;
   std::complex<double> tmp_3290;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3291;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3291 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3290 += tmp_3291;
   }
   tmp_3289 += tmp_3290;
   result += (-0.5) * tmp_3289;
   std::complex<double> tmp_3292;
   std::complex<double> tmp_3293;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3293 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_3292 += tmp_3293;
   result += (-1) * tmp_3292;
   std::complex<double> tmp_3294;
   std::complex<double> tmp_3295;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3295 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_3294 += tmp_3295;
   result += (-1) * tmp_3294;
   std::complex<double> tmp_3296;
   std::complex<double> tmp_3297;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3297 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_3296 += tmp_3297;
   result += (-1) * tmp_3296;
   std::complex<double> tmp_3298;
   std::complex<double> tmp_3299;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3299 += B1(p,MFe(gI2),MVZp)*Conj(CpbarUFeVZpFePR(gO2,gI2))*
         CpbarUFeVZpFePR(gO1,gI2);
   }
   tmp_3298 += tmp_3299;
   result += (-1) * tmp_3298;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3300;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3301;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3301 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3300 += tmp_3301;
   }
   result += tmp_3300;
   std::complex<double> tmp_3302;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3303;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3303 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3302 += tmp_3303;
   }
   result += tmp_3302;
   std::complex<double> tmp_3304;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3305;
      std::complex<double> tmp_3306;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3306 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3305 += tmp_3306;
      tmp_3304 += (MFd(gI1)) * tmp_3305;
   }
   result += tmp_3304;
   std::complex<double> tmp_3307;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3308;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3308 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3307 += tmp_3308;
   }
   result += tmp_3307;
   std::complex<double> tmp_3309;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3310;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3310 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3309 += tmp_3310;
   }
   result += tmp_3309;
   std::complex<double> tmp_3311;
   std::complex<double> tmp_3312;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3312 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3311 += tmp_3312;
   result += (-5.333333333333333) * tmp_3311;
   std::complex<double> tmp_3313;
   std::complex<double> tmp_3314;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3314 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3313 += tmp_3314;
   result += (-4) * tmp_3313;
   std::complex<double> tmp_3315;
   std::complex<double> tmp_3316;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3316 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3315 += tmp_3316;
   result += (-4) * tmp_3315;
   std::complex<double> tmp_3317;
   std::complex<double> tmp_3318;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3318 += B0(p,MFd(gI2),MVZp)*Conj(CpbarUFdVZpFdPR(gO2,gI2))*
         CpbarUFdVZpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3317 += tmp_3318;
   result += (-4) * tmp_3317;
   std::complex<double> tmp_3319;
   std::complex<double> tmp_3320;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3320 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3319 += tmp_3320;
   result += (-4) * tmp_3319;
   std::complex<double> tmp_3321;
   std::complex<double> tmp_3322;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3322 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_3321 += tmp_3322;
   result += (1.3333333333333333*MGlu) * tmp_3321;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3323;
   std::complex<double> tmp_3324;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3325;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3325 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3324 += tmp_3325;
   }
   tmp_3323 += tmp_3324;
   result += (-0.5) * tmp_3323;
   std::complex<double> tmp_3326;
   std::complex<double> tmp_3327;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3328;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3328 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3327 += tmp_3328;
   }
   tmp_3326 += tmp_3327;
   result += (-0.5) * tmp_3326;
   std::complex<double> tmp_3329;
   std::complex<double> tmp_3330;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3331;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3331 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3330 += tmp_3331;
   }
   tmp_3329 += tmp_3330;
   result += (-0.5) * tmp_3329;
   std::complex<double> tmp_3332;
   std::complex<double> tmp_3333;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3333 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_3332 += tmp_3333;
   result += (-0.6666666666666666) * tmp_3332;
   std::complex<double> tmp_3334;
   std::complex<double> tmp_3335;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3336;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3336 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPR(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_3335 += tmp_3336;
   }
   tmp_3334 += tmp_3335;
   result += (-0.5) * tmp_3334;
   std::complex<double> tmp_3337;
   std::complex<double> tmp_3338;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3339;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3339 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3338 += tmp_3339;
   }
   tmp_3337 += tmp_3338;
   result += (-0.5) * tmp_3337;
   std::complex<double> tmp_3340;
   std::complex<double> tmp_3341;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3341 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_3340 += tmp_3341;
   result += (-1.3333333333333333) * tmp_3340;
   std::complex<double> tmp_3342;
   std::complex<double> tmp_3343;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3343 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_3342 += tmp_3343;
   result += (-1) * tmp_3342;
   std::complex<double> tmp_3344;
   std::complex<double> tmp_3345;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3345 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_3344 += tmp_3345;
   result += (-1) * tmp_3344;
   std::complex<double> tmp_3346;
   std::complex<double> tmp_3347;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3347 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_3346 += tmp_3347;
   result += (-1) * tmp_3346;
   std::complex<double> tmp_3348;
   std::complex<double> tmp_3349;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3349 += B1(p,MFd(gI2),MVZp)*Conj(CpbarUFdVZpFdPL(gO2,gI2))*
         CpbarUFdVZpFdPL(gO1,gI2);
   }
   tmp_3348 += tmp_3349;
   result += (-1) * tmp_3348;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3350;
   std::complex<double> tmp_3351;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3352;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3352 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3351 += tmp_3352;
   }
   tmp_3350 += tmp_3351;
   result += (-0.5) * tmp_3350;
   std::complex<double> tmp_3353;
   std::complex<double> tmp_3354;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3355;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3355 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3354 += tmp_3355;
   }
   tmp_3353 += tmp_3354;
   result += (-0.5) * tmp_3353;
   std::complex<double> tmp_3356;
   std::complex<double> tmp_3357;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3358;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3358 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3357 += tmp_3358;
   }
   tmp_3356 += tmp_3357;
   result += (-0.5) * tmp_3356;
   std::complex<double> tmp_3359;
   std::complex<double> tmp_3360;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3360 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPL(gO1,gI1,1);
   }
   tmp_3359 += tmp_3360;
   result += (-0.6666666666666666) * tmp_3359;
   std::complex<double> tmp_3361;
   std::complex<double> tmp_3362;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3363;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3363 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_3362 += tmp_3363;
   }
   tmp_3361 += tmp_3362;
   result += (-0.5) * tmp_3361;
   std::complex<double> tmp_3364;
   std::complex<double> tmp_3365;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3366;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3366 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3365 += tmp_3366;
   }
   tmp_3364 += tmp_3365;
   result += (-0.5) * tmp_3364;
   std::complex<double> tmp_3367;
   std::complex<double> tmp_3368;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3368 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_3367 += tmp_3368;
   result += (-1.3333333333333333) * tmp_3367;
   std::complex<double> tmp_3369;
   std::complex<double> tmp_3370;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3370 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_3369 += tmp_3370;
   result += (-1) * tmp_3369;
   std::complex<double> tmp_3371;
   std::complex<double> tmp_3372;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3372 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_3371 += tmp_3372;
   result += (-1) * tmp_3371;
   std::complex<double> tmp_3373;
   std::complex<double> tmp_3374;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3374 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_3373 += tmp_3374;
   result += (-1) * tmp_3373;
   std::complex<double> tmp_3375;
   std::complex<double> tmp_3376;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3376 += B1(p,MFd(gI2),MVZp)*Conj(CpbarUFdVZpFdPR(gO2,gI2))*
         CpbarUFdVZpFdPR(gO1,gI2);
   }
   tmp_3375 += tmp_3376;
   result += (-1) * tmp_3375;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3377;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3378;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3378 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3377 += tmp_3378;
   }
   result += tmp_3377;
   std::complex<double> tmp_3379;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3380;
      std::complex<double> tmp_3381;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3381 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3380 += tmp_3381;
      tmp_3379 += (MCha(gI1)) * tmp_3380;
   }
   result += tmp_3379;
   std::complex<double> tmp_3382;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3383;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3383 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3382 += tmp_3383;
   }
   result += tmp_3382;
   std::complex<double> tmp_3384;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3385;
      std::complex<double> tmp_3386;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3386 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3385 += tmp_3386;
      tmp_3384 += (MFu(gI1)) * tmp_3385;
   }
   result += tmp_3384;
   std::complex<double> tmp_3387;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3388;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3388 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3387 += tmp_3388;
   }
   result += tmp_3387;
   std::complex<double> tmp_3389;
   std::complex<double> tmp_3390;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3390 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3389 += tmp_3390;
   result += (-4) * tmp_3389;
   std::complex<double> tmp_3391;
   std::complex<double> tmp_3392;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3392 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3391 += tmp_3392;
   result += (-5.333333333333333) * tmp_3391;
   std::complex<double> tmp_3393;
   std::complex<double> tmp_3394;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3394 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3393 += tmp_3394;
   result += (-4) * tmp_3393;
   std::complex<double> tmp_3395;
   std::complex<double> tmp_3396;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3396 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3395 += tmp_3396;
   result += (-4) * tmp_3395;
   std::complex<double> tmp_3397;
   std::complex<double> tmp_3398;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3398 += B0(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3397 += tmp_3398;
   result += (-4) * tmp_3397;
   std::complex<double> tmp_3399;
   std::complex<double> tmp_3400;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3400 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3399 += tmp_3400;
   result += (1.3333333333333333*MGlu) * tmp_3399;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3401;
   std::complex<double> tmp_3402;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3403;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3403 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3402 += tmp_3403;
   }
   tmp_3401 += tmp_3402;
   result += (-0.5) * tmp_3401;
   std::complex<double> tmp_3404;
   std::complex<double> tmp_3405;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3406;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3406 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3405 += tmp_3406;
   }
   tmp_3404 += tmp_3405;
   result += (-0.5) * tmp_3404;
   std::complex<double> tmp_3407;
   std::complex<double> tmp_3408;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3409;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3409 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3408 += tmp_3409;
   }
   tmp_3407 += tmp_3408;
   result += (-0.5) * tmp_3407;
   std::complex<double> tmp_3410;
   std::complex<double> tmp_3411;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3412;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3412 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3411 += tmp_3412;
   }
   tmp_3410 += tmp_3411;
   result += (-0.5) * tmp_3410;
   std::complex<double> tmp_3413;
   std::complex<double> tmp_3414;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3414 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3413 += tmp_3414;
   result += (-0.6666666666666666) * tmp_3413;
   std::complex<double> tmp_3415;
   std::complex<double> tmp_3416;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3417;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3417 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3416 += tmp_3417;
   }
   tmp_3415 += tmp_3416;
   result += (-0.5) * tmp_3415;
   std::complex<double> tmp_3418;
   std::complex<double> tmp_3419;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3419 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3418 += tmp_3419;
   result += (-1) * tmp_3418;
   std::complex<double> tmp_3420;
   std::complex<double> tmp_3421;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3421 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_3420 += tmp_3421;
   result += (-1.3333333333333333) * tmp_3420;
   std::complex<double> tmp_3422;
   std::complex<double> tmp_3423;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3423 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3422 += tmp_3423;
   result += (-1) * tmp_3422;
   std::complex<double> tmp_3424;
   std::complex<double> tmp_3425;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3425 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3424 += tmp_3425;
   result += (-1) * tmp_3424;
   std::complex<double> tmp_3426;
   std::complex<double> tmp_3427;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3427 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPL(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2);
   }
   tmp_3426 += tmp_3427;
   result += (-1) * tmp_3426;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3428;
   std::complex<double> tmp_3429;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3430;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3430 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3429 += tmp_3430;
   }
   tmp_3428 += tmp_3429;
   result += (-0.5) * tmp_3428;
   std::complex<double> tmp_3431;
   std::complex<double> tmp_3432;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3433;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3433 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3432 += tmp_3433;
   }
   tmp_3431 += tmp_3432;
   result += (-0.5) * tmp_3431;
   std::complex<double> tmp_3434;
   std::complex<double> tmp_3435;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3436;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3436 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3435 += tmp_3436;
   }
   tmp_3434 += tmp_3435;
   result += (-0.5) * tmp_3434;
   std::complex<double> tmp_3437;
   std::complex<double> tmp_3438;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3439;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3439 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3438 += tmp_3439;
   }
   tmp_3437 += tmp_3438;
   result += (-0.5) * tmp_3437;
   std::complex<double> tmp_3440;
   std::complex<double> tmp_3441;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3441 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_3440 += tmp_3441;
   result += (-0.6666666666666666) * tmp_3440;
   std::complex<double> tmp_3442;
   std::complex<double> tmp_3443;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3444;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3444 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3443 += tmp_3444;
   }
   tmp_3442 += tmp_3443;
   result += (-0.5) * tmp_3442;
   std::complex<double> tmp_3445;
   std::complex<double> tmp_3446;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3446 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3445 += tmp_3446;
   result += (-1) * tmp_3445;
   std::complex<double> tmp_3447;
   std::complex<double> tmp_3448;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3448 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_3447 += tmp_3448;
   result += (-1.3333333333333333) * tmp_3447;
   std::complex<double> tmp_3449;
   std::complex<double> tmp_3450;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3450 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3449 += tmp_3450;
   result += (-1) * tmp_3449;
   std::complex<double> tmp_3451;
   std::complex<double> tmp_3452;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3452 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3451 += tmp_3452;
   result += (-1) * tmp_3451;
   std::complex<double> tmp_3453;
   std::complex<double> tmp_3454;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3454 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPR(gO1,gI2);
   }
   tmp_3453 += tmp_3454;
   result += (-1) * tmp_3453;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3455;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3456;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3456 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpGluconjSdFdPL(gI1,gI2
            ))*CpGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_3455 += tmp_3456;
   }
   result += tmp_3455;
   std::complex<double> tmp_3457;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3458;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3458 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpGluconjSuFuPL(gI1,gI2
            ))*CpGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_3457 += tmp_3458;
   }
   result += tmp_3457;
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpGluVGGluPR())*CpGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_3459;
   std::complex<double> tmp_3460;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3461;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3461 += AbsSqr(CpGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_3460 += tmp_3461;
   }
   tmp_3459 += tmp_3460;
   result += (-0.5) * tmp_3459;
   std::complex<double> tmp_3462;
   std::complex<double> tmp_3463;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3464;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3464 += AbsSqr(CpGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_3463 += tmp_3464;
   }
   tmp_3462 += tmp_3463;
   result += (-0.5) * tmp_3462;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_3465;
   std::complex<double> tmp_3466;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3467;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3467 += AbsSqr(CpGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_3466 += tmp_3467;
   }
   tmp_3465 += tmp_3466;
   result += (-0.5) * tmp_3465;
   std::complex<double> tmp_3468;
   std::complex<double> tmp_3469;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3470;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3470 += AbsSqr(CpGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_3469 += tmp_3470;
   }
   tmp_3468 += tmp_3469;
   result += (-0.5) * tmp_3468;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3471;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3472;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3472 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_3472 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_3471 += tmp_3472;
   }
   result += tmp_3471;
   std::complex<double> tmp_3473;
   std::complex<double> tmp_3474;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3474 += AbsSqr(CpVZhhAh(gI1,2))*B00(p,MAh(2),Mhh(gI1));
   }
   tmp_3473 += tmp_3474;
   result += (-4) * tmp_3473;
   std::complex<double> tmp_3475;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3475 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_3475;
   std::complex<double> tmp_3476;
   std::complex<double> tmp_3477;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3477 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_3476 += tmp_3477;
   result += (0.5) * tmp_3476;
   std::complex<double> tmp_3478;
   std::complex<double> tmp_3479;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3480;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3480 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_3479 += tmp_3480;
   }
   tmp_3478 += tmp_3479;
   result += (-4) * tmp_3478;
   std::complex<double> tmp_3481;
   std::complex<double> tmp_3482;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3482 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_3481 += tmp_3482;
   result += (3) * tmp_3481;
   std::complex<double> tmp_3483;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3483 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_3483;
   std::complex<double> tmp_3484;
   std::complex<double> tmp_3485;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3485 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_3484 += tmp_3485;
   result += (3) * tmp_3484;
   std::complex<double> tmp_3486;
   std::complex<double> tmp_3487;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3488;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3488 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_3487 += tmp_3488;
   }
   tmp_3486 += tmp_3487;
   result += (-12) * tmp_3486;
   std::complex<double> tmp_3489;
   std::complex<double> tmp_3490;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3491;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3491 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_3490 += tmp_3491;
   }
   tmp_3489 += tmp_3490;
   result += (-4) * tmp_3489;
   std::complex<double> tmp_3492;
   std::complex<double> tmp_3493;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3494;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3494 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_3493 += tmp_3494;
   }
   tmp_3492 += tmp_3493;
   result += (-12) * tmp_3492;
   std::complex<double> tmp_3495;
   std::complex<double> tmp_3496;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3497;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3497 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_3497 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_3496 += tmp_3497;
   }
   tmp_3495 += tmp_3496;
   result += (0.5) * tmp_3495;
   std::complex<double> tmp_3498;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3498 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_3498;
   std::complex<double> tmp_3499;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3499 += AbsSqr(CpVZVZphh(gI2))*B0(p,MVZp,Mhh(gI2));
   }
   result += tmp_3499;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpconjVWmbarVWmVZp())*B00(p,MVZp,MVWm);
   result += AbsSqr(CpconjVWmbarVZpVWm())*B00(p,MVWm,MVZp);
   result += AbsSqr(CpconjVWmVZpHpm(1))*B0(p,MVZp,MHpm(1));
   result += -0.5*A0(MVZp)*(4*CpVWmconjVWmVZpVZp1() + CpVWmconjVWmVZpVZp2() +
      CpVWmconjVWmVZpVZp3());
   std::complex<double> tmp_3500;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3500 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_3500;
   std::complex<double> tmp_3501;
   std::complex<double> tmp_3502;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3502 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_3501 += tmp_3502;
   result += (0.5) * tmp_3501;
   std::complex<double> tmp_3503;
   std::complex<double> tmp_3504;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3505;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3505 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_3504 += tmp_3505;
   }
   tmp_3503 += tmp_3504;
   result += (-4) * tmp_3503;
   std::complex<double> tmp_3506;
   std::complex<double> tmp_3507;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3507 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_3506 += tmp_3507;
   result += (3) * tmp_3506;
   std::complex<double> tmp_3508;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3508 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_3508;
   std::complex<double> tmp_3509;
   std::complex<double> tmp_3510;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3510 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_3509 += tmp_3510;
   result += (3) * tmp_3509;
   std::complex<double> tmp_3511;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3512;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3512 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_3512 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_3511 += tmp_3512;
   }
   result += tmp_3511;
   std::complex<double> tmp_3513;
   std::complex<double> tmp_3514;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3515;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3515 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_3514 += tmp_3515;
   }
   tmp_3513 += tmp_3514;
   result += (-12) * tmp_3513;
   std::complex<double> tmp_3516;
   std::complex<double> tmp_3517;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3517 += AbsSqr(CpconjVWmHpmhh(1,gI2))*B00(p,Mhh(gI2),MHpm(1));
   }
   tmp_3516 += tmp_3517;
   result += (-4) * tmp_3516;
   std::complex<double> tmp_3518;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3518 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_3518;
   result += -(AbsSqr(CpconjVWmVZpVWm())*(A0(MVWm) + A0(MVZp) + 10*B00(p,MVZp,
      MVWm) + B0(p,MVZp,MVWm)*(Sqr(MVWm) + Sqr(MVZp) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3519;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3520;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3520 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3519 += tmp_3520;
   }
   result += tmp_3519;
   std::complex<double> tmp_3521;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3522;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3522 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3521 += tmp_3522;
   }
   result += tmp_3521;
   std::complex<double> tmp_3523;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3524;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3524 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3523 += tmp_3524;
   }
   result += tmp_3523;
   std::complex<double> tmp_3525;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3526;
      std::complex<double> tmp_3527;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3527 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3526 += tmp_3527;
      tmp_3525 += (MFe(gI1)) * tmp_3526;
   }
   result += tmp_3525;
   std::complex<double> tmp_3528;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3529;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3529 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3528 += tmp_3529;
   }
   result += tmp_3528;
   std::complex<double> tmp_3530;
   std::complex<double> tmp_3531;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3531 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3530 += tmp_3531;
   result += (-4) * tmp_3530;
   std::complex<double> tmp_3532;
   std::complex<double> tmp_3533;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3533 += B0(p,MFe(gI2),MVZp)*Conj(CpbarFeVZpFePR(gO2,gI2))*
         CpbarFeVZpFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3532 += tmp_3533;
   result += (-4) * tmp_3532;
   std::complex<double> tmp_3534;
   std::complex<double> tmp_3535;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3535 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3534 += tmp_3535;
   result += (-4) * tmp_3534;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3536;
   std::complex<double> tmp_3537;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3538;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3538 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3537 += tmp_3538;
   }
   tmp_3536 += tmp_3537;
   result += (-0.5) * tmp_3536;
   std::complex<double> tmp_3539;
   std::complex<double> tmp_3540;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3541;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3541 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPR(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_3540 += tmp_3541;
   }
   tmp_3539 += tmp_3540;
   result += (-0.5) * tmp_3539;
   std::complex<double> tmp_3542;
   std::complex<double> tmp_3543;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3544;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3544 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3543 += tmp_3544;
   }
   tmp_3542 += tmp_3543;
   result += (-0.5) * tmp_3542;
   std::complex<double> tmp_3545;
   std::complex<double> tmp_3546;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3547;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3547 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_3546 += tmp_3547;
   }
   tmp_3545 += tmp_3546;
   result += (-0.5) * tmp_3545;
   std::complex<double> tmp_3548;
   std::complex<double> tmp_3549;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3550;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3550 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3549 += tmp_3550;
   }
   tmp_3548 += tmp_3549;
   result += (-0.5) * tmp_3548;
   std::complex<double> tmp_3551;
   std::complex<double> tmp_3552;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3552 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_3551 += tmp_3552;
   result += (-1) * tmp_3551;
   std::complex<double> tmp_3553;
   std::complex<double> tmp_3554;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3554 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_3553 += tmp_3554;
   result += (-1) * tmp_3553;
   std::complex<double> tmp_3555;
   std::complex<double> tmp_3556;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3556 += B1(p,MFe(gI2),MVZp)*Conj(CpbarFeVZpFePL(gO2,gI2))*
         CpbarFeVZpFePL(gO1,gI2);
   }
   tmp_3555 += tmp_3556;
   result += (-1) * tmp_3555;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3557;
   std::complex<double> tmp_3558;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3559;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3559 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3558 += tmp_3559;
   }
   tmp_3557 += tmp_3558;
   result += (-0.5) * tmp_3557;
   std::complex<double> tmp_3560;
   std::complex<double> tmp_3561;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3562;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3562 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_3561 += tmp_3562;
   }
   tmp_3560 += tmp_3561;
   result += (-0.5) * tmp_3560;
   std::complex<double> tmp_3563;
   std::complex<double> tmp_3564;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3565;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3565 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3564 += tmp_3565;
   }
   tmp_3563 += tmp_3564;
   result += (-0.5) * tmp_3563;
   std::complex<double> tmp_3566;
   std::complex<double> tmp_3567;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3568;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3568 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_3567 += tmp_3568;
   }
   tmp_3566 += tmp_3567;
   result += (-0.5) * tmp_3566;
   std::complex<double> tmp_3569;
   std::complex<double> tmp_3570;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3571;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3571 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3570 += tmp_3571;
   }
   tmp_3569 += tmp_3570;
   result += (-0.5) * tmp_3569;
   std::complex<double> tmp_3572;
   std::complex<double> tmp_3573;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3573 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_3572 += tmp_3573;
   result += (-1) * tmp_3572;
   std::complex<double> tmp_3574;
   std::complex<double> tmp_3575;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3575 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_3574 += tmp_3575;
   result += (-1) * tmp_3574;
   std::complex<double> tmp_3576;
   std::complex<double> tmp_3577;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3577 += B1(p,MFe(gI2),MVZp)*Conj(CpbarFeVZpFePR(gO2,gI2))*
         CpbarFeVZpFePR(gO1,gI2);
   }
   tmp_3576 += tmp_3577;
   result += (-1) * tmp_3576;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3578;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3579;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3579 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3578 += tmp_3579;
   }
   result += tmp_3578;
   std::complex<double> tmp_3580;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3581;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3581 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3580 += tmp_3581;
   }
   result += tmp_3580;
   std::complex<double> tmp_3582;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3583;
      std::complex<double> tmp_3584;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3584 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3583 += tmp_3584;
      tmp_3582 += (MFd(gI1)) * tmp_3583;
   }
   result += tmp_3582;
   std::complex<double> tmp_3585;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3586;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3586 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3585 += tmp_3586;
   }
   result += tmp_3585;
   std::complex<double> tmp_3587;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3588;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3588 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3587 += tmp_3588;
   }
   result += tmp_3587;
   std::complex<double> tmp_3589;
   std::complex<double> tmp_3590;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3590 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3589 += tmp_3590;
   result += (-4) * tmp_3589;
   std::complex<double> tmp_3591;
   std::complex<double> tmp_3592;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3592 += B0(p,MFd(gI2),MVZp)*Conj(CpbarFdVZpFdPR(gO2,gI2))*
         CpbarFdVZpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3591 += tmp_3592;
   result += (-4) * tmp_3591;
   std::complex<double> tmp_3593;
   std::complex<double> tmp_3594;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3594 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3593 += tmp_3594;
   result += (-4) * tmp_3593;
   std::complex<double> tmp_3595;
   std::complex<double> tmp_3596;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3596 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3595 += tmp_3596;
   result += (1.3333333333333333*MGlu) * tmp_3595;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3597;
   std::complex<double> tmp_3598;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3599;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3599 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3598 += tmp_3599;
   }
   tmp_3597 += tmp_3598;
   result += (-0.5) * tmp_3597;
   std::complex<double> tmp_3600;
   std::complex<double> tmp_3601;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3602;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3602 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3601 += tmp_3602;
   }
   tmp_3600 += tmp_3601;
   result += (-0.5) * tmp_3600;
   std::complex<double> tmp_3603;
   std::complex<double> tmp_3604;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3605;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3605 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3604 += tmp_3605;
   }
   tmp_3603 += tmp_3604;
   result += (-0.5) * tmp_3603;
   std::complex<double> tmp_3606;
   std::complex<double> tmp_3607;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3607 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3606 += tmp_3607;
   result += (-0.6666666666666666) * tmp_3606;
   std::complex<double> tmp_3608;
   std::complex<double> tmp_3609;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3610;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3610 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPR(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_3609 += tmp_3610;
   }
   tmp_3608 += tmp_3609;
   result += (-0.5) * tmp_3608;
   std::complex<double> tmp_3611;
   std::complex<double> tmp_3612;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3613;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3613 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3612 += tmp_3613;
   }
   tmp_3611 += tmp_3612;
   result += (-0.5) * tmp_3611;
   std::complex<double> tmp_3614;
   std::complex<double> tmp_3615;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3615 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_3614 += tmp_3615;
   result += (-1) * tmp_3614;
   std::complex<double> tmp_3616;
   std::complex<double> tmp_3617;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3617 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_3616 += tmp_3617;
   result += (-1) * tmp_3616;
   std::complex<double> tmp_3618;
   std::complex<double> tmp_3619;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3619 += B1(p,MFd(gI2),MVZp)*Conj(CpbarFdVZpFdPL(gO2,gI2))*
         CpbarFdVZpFdPL(gO1,gI2);
   }
   tmp_3618 += tmp_3619;
   result += (-1) * tmp_3618;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3620;
   std::complex<double> tmp_3621;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3622;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3622 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3621 += tmp_3622;
   }
   tmp_3620 += tmp_3621;
   result += (-0.5) * tmp_3620;
   std::complex<double> tmp_3623;
   std::complex<double> tmp_3624;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3625;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3625 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3624 += tmp_3625;
   }
   tmp_3623 += tmp_3624;
   result += (-0.5) * tmp_3623;
   std::complex<double> tmp_3626;
   std::complex<double> tmp_3627;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3628;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3628 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3627 += tmp_3628;
   }
   tmp_3626 += tmp_3627;
   result += (-0.5) * tmp_3626;
   std::complex<double> tmp_3629;
   std::complex<double> tmp_3630;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3630 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPL(gO1,gI1,1);
   }
   tmp_3629 += tmp_3630;
   result += (-0.6666666666666666) * tmp_3629;
   std::complex<double> tmp_3631;
   std::complex<double> tmp_3632;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3633;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3633 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_3632 += tmp_3633;
   }
   tmp_3631 += tmp_3632;
   result += (-0.5) * tmp_3631;
   std::complex<double> tmp_3634;
   std::complex<double> tmp_3635;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3636;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3636 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3635 += tmp_3636;
   }
   tmp_3634 += tmp_3635;
   result += (-0.5) * tmp_3634;
   std::complex<double> tmp_3637;
   std::complex<double> tmp_3638;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3638 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_3637 += tmp_3638;
   result += (-1) * tmp_3637;
   std::complex<double> tmp_3639;
   std::complex<double> tmp_3640;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3640 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_3639 += tmp_3640;
   result += (-1) * tmp_3639;
   std::complex<double> tmp_3641;
   std::complex<double> tmp_3642;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3642 += B1(p,MFd(gI2),MVZp)*Conj(CpbarFdVZpFdPR(gO2,gI2))*
         CpbarFdVZpFdPR(gO1,gI2);
   }
   tmp_3641 += tmp_3642;
   result += (-1) * tmp_3641;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3643;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3644;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3644 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3643 += tmp_3644;
   }
   result += tmp_3643;
   std::complex<double> tmp_3645;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3646;
      std::complex<double> tmp_3647;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3647 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3646 += tmp_3647;
      tmp_3645 += (MCha(gI1)) * tmp_3646;
   }
   result += tmp_3645;
   std::complex<double> tmp_3648;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3649;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3649 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3648 += tmp_3649;
   }
   result += tmp_3648;
   std::complex<double> tmp_3650;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3651;
      std::complex<double> tmp_3652;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3652 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3651 += tmp_3652;
      tmp_3650 += (MFu(gI1)) * tmp_3651;
   }
   result += tmp_3650;
   std::complex<double> tmp_3653;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3654;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3654 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3653 += tmp_3654;
   }
   result += tmp_3653;
   std::complex<double> tmp_3655;
   std::complex<double> tmp_3656;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3656 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3655 += tmp_3656;
   result += (-4) * tmp_3655;
   std::complex<double> tmp_3657;
   std::complex<double> tmp_3658;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3658 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3657 += tmp_3658;
   result += (-4) * tmp_3657;
   std::complex<double> tmp_3659;
   std::complex<double> tmp_3660;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3660 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3659 += tmp_3660;
   result += (-4) * tmp_3659;
   std::complex<double> tmp_3661;
   std::complex<double> tmp_3662;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3662 += B0(p,MFu(gI2),MVZp)*Conj(CpbarFuVZpFuPR(gO2,gI2))*
         CpbarFuVZpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3661 += tmp_3662;
   result += (-4) * tmp_3661;
   std::complex<double> tmp_3663;
   std::complex<double> tmp_3664;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3664 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3663 += tmp_3664;
   result += (1.3333333333333333*MGlu) * tmp_3663;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3665;
   std::complex<double> tmp_3666;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3667;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3667 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3666 += tmp_3667;
   }
   tmp_3665 += tmp_3666;
   result += (-0.5) * tmp_3665;
   std::complex<double> tmp_3668;
   std::complex<double> tmp_3669;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3670;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3670 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPR(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3669 += tmp_3670;
   }
   tmp_3668 += tmp_3669;
   result += (-0.5) * tmp_3668;
   std::complex<double> tmp_3671;
   std::complex<double> tmp_3672;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3673;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3673 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3672 += tmp_3673;
   }
   tmp_3671 += tmp_3672;
   result += (-0.5) * tmp_3671;
   std::complex<double> tmp_3674;
   std::complex<double> tmp_3675;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3676;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3676 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3675 += tmp_3676;
   }
   tmp_3674 += tmp_3675;
   result += (-0.5) * tmp_3674;
   std::complex<double> tmp_3677;
   std::complex<double> tmp_3678;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3678 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3677 += tmp_3678;
   result += (-0.6666666666666666) * tmp_3677;
   std::complex<double> tmp_3679;
   std::complex<double> tmp_3680;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3681;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3681 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3680 += tmp_3681;
   }
   tmp_3679 += tmp_3680;
   result += (-0.5) * tmp_3679;
   std::complex<double> tmp_3682;
   std::complex<double> tmp_3683;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3683 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3682 += tmp_3683;
   result += (-1) * tmp_3682;
   std::complex<double> tmp_3684;
   std::complex<double> tmp_3685;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3685 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_3684 += tmp_3685;
   result += (-1) * tmp_3684;
   std::complex<double> tmp_3686;
   std::complex<double> tmp_3687;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3687 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_3686 += tmp_3687;
   result += (-1) * tmp_3686;
   std::complex<double> tmp_3688;
   std::complex<double> tmp_3689;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3689 += B1(p,MFu(gI2),MVZp)*Conj(CpbarFuVZpFuPL(gO2,gI2))*
         CpbarFuVZpFuPL(gO1,gI2);
   }
   tmp_3688 += tmp_3689;
   result += (-1) * tmp_3688;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3690;
   std::complex<double> tmp_3691;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3692;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3692 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3691 += tmp_3692;
   }
   tmp_3690 += tmp_3691;
   result += (-0.5) * tmp_3690;
   std::complex<double> tmp_3693;
   std::complex<double> tmp_3694;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3695;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3695 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3694 += tmp_3695;
   }
   tmp_3693 += tmp_3694;
   result += (-0.5) * tmp_3693;
   std::complex<double> tmp_3696;
   std::complex<double> tmp_3697;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3698;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3698 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3697 += tmp_3698;
   }
   tmp_3696 += tmp_3697;
   result += (-0.5) * tmp_3696;
   std::complex<double> tmp_3699;
   std::complex<double> tmp_3700;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3701;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3701 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3700 += tmp_3701;
   }
   tmp_3699 += tmp_3700;
   result += (-0.5) * tmp_3699;
   std::complex<double> tmp_3702;
   std::complex<double> tmp_3703;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3703 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPL(gO1,gI1,1);
   }
   tmp_3702 += tmp_3703;
   result += (-0.6666666666666666) * tmp_3702;
   std::complex<double> tmp_3704;
   std::complex<double> tmp_3705;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3706;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3706 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3705 += tmp_3706;
   }
   tmp_3704 += tmp_3705;
   result += (-0.5) * tmp_3704;
   std::complex<double> tmp_3707;
   std::complex<double> tmp_3708;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3708 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3707 += tmp_3708;
   result += (-1) * tmp_3707;
   std::complex<double> tmp_3709;
   std::complex<double> tmp_3710;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3710 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_3709 += tmp_3710;
   result += (-1) * tmp_3709;
   std::complex<double> tmp_3711;
   std::complex<double> tmp_3712;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3712 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_3711 += tmp_3712;
   result += (-1) * tmp_3711;
   std::complex<double> tmp_3713;
   std::complex<double> tmp_3714;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3714 += B1(p,MFu(gI2),MVZp)*Conj(CpbarFuVZpFuPR(gO2,gI2))*
         CpbarFuVZpFuPR(gO1,gI2);
   }
   tmp_3713 += tmp_3714;
   result += (-1) * tmp_3713;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3715;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3716;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3716 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3715 += tmp_3716;
   }
   result += tmp_3715;
   std::complex<double> tmp_3717;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3718;
      std::complex<double> tmp_3719;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3719 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3718 += tmp_3719;
      tmp_3717 += (MCha(gI1)) * tmp_3718;
   }
   result += tmp_3717;
   std::complex<double> tmp_3720;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3721;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3721 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3720 += tmp_3721;
   }
   result += tmp_3720;
   std::complex<double> tmp_3722;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3723;
      std::complex<double> tmp_3724;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3724 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3723 += tmp_3724;
      tmp_3722 += (MFu(gI1)) * tmp_3723;
   }
   result += tmp_3722;
   std::complex<double> tmp_3725;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3726;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3726 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3725 += tmp_3726;
   }
   result += tmp_3725;
   std::complex<double> tmp_3727;
   std::complex<double> tmp_3728;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3728 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3727 += tmp_3728;
   result += (-4) * tmp_3727;
   std::complex<double> tmp_3729;
   std::complex<double> tmp_3730;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3730 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3729 += tmp_3730;
   result += (-4) * tmp_3729;
   std::complex<double> tmp_3731;
   std::complex<double> tmp_3732;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3732 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3731 += tmp_3732;
   result += (-4) * tmp_3731;
   std::complex<double> tmp_3733;
   std::complex<double> tmp_3734;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3734 += B0(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3733 += tmp_3734;
   result += (-4) * tmp_3733;
   std::complex<double> tmp_3735;
   std::complex<double> tmp_3736;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3736 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3735 += tmp_3736;
   result += (1.3333333333333333*MGlu) * tmp_3735;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3737;
   std::complex<double> tmp_3738;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3739;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3739 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3738 += tmp_3739;
   }
   tmp_3737 += tmp_3738;
   result += (-0.5) * tmp_3737;
   std::complex<double> tmp_3740;
   std::complex<double> tmp_3741;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3742;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3742 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3741 += tmp_3742;
   }
   tmp_3740 += tmp_3741;
   result += (-0.5) * tmp_3740;
   std::complex<double> tmp_3743;
   std::complex<double> tmp_3744;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3745;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3745 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3744 += tmp_3745;
   }
   tmp_3743 += tmp_3744;
   result += (-0.5) * tmp_3743;
   std::complex<double> tmp_3746;
   std::complex<double> tmp_3747;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3748;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3748 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3747 += tmp_3748;
   }
   tmp_3746 += tmp_3747;
   result += (-0.5) * tmp_3746;
   std::complex<double> tmp_3749;
   std::complex<double> tmp_3750;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3750 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3749 += tmp_3750;
   result += (-0.6666666666666666) * tmp_3749;
   std::complex<double> tmp_3751;
   std::complex<double> tmp_3752;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3753;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3753 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3752 += tmp_3753;
   }
   tmp_3751 += tmp_3752;
   result += (-0.5) * tmp_3751;
   std::complex<double> tmp_3754;
   std::complex<double> tmp_3755;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3755 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3754 += tmp_3755;
   result += (-1) * tmp_3754;
   std::complex<double> tmp_3756;
   std::complex<double> tmp_3757;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3757 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3756 += tmp_3757;
   result += (-1) * tmp_3756;
   std::complex<double> tmp_3758;
   std::complex<double> tmp_3759;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3759 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3758 += tmp_3759;
   result += (-1) * tmp_3758;
   std::complex<double> tmp_3760;
   std::complex<double> tmp_3761;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3761 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPL(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2);
   }
   tmp_3760 += tmp_3761;
   result += (-1) * tmp_3760;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3762;
   std::complex<double> tmp_3763;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3764;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3764 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3763 += tmp_3764;
   }
   tmp_3762 += tmp_3763;
   result += (-0.5) * tmp_3762;
   std::complex<double> tmp_3765;
   std::complex<double> tmp_3766;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3767;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3767 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3766 += tmp_3767;
   }
   tmp_3765 += tmp_3766;
   result += (-0.5) * tmp_3765;
   std::complex<double> tmp_3768;
   std::complex<double> tmp_3769;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3770;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3770 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3769 += tmp_3770;
   }
   tmp_3768 += tmp_3769;
   result += (-0.5) * tmp_3768;
   std::complex<double> tmp_3771;
   std::complex<double> tmp_3772;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3773;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3773 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3772 += tmp_3773;
   }
   tmp_3771 += tmp_3772;
   result += (-0.5) * tmp_3771;
   std::complex<double> tmp_3774;
   std::complex<double> tmp_3775;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3775 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_3774 += tmp_3775;
   result += (-0.6666666666666666) * tmp_3774;
   std::complex<double> tmp_3776;
   std::complex<double> tmp_3777;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3778;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3778 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3777 += tmp_3778;
   }
   tmp_3776 += tmp_3777;
   result += (-0.5) * tmp_3776;
   std::complex<double> tmp_3779;
   std::complex<double> tmp_3780;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3780 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3779 += tmp_3780;
   result += (-1) * tmp_3779;
   std::complex<double> tmp_3781;
   std::complex<double> tmp_3782;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3782 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3781 += tmp_3782;
   result += (-1) * tmp_3781;
   std::complex<double> tmp_3783;
   std::complex<double> tmp_3784;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3784 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3783 += tmp_3784;
   result += (-1) * tmp_3783;
   std::complex<double> tmp_3785;
   std::complex<double> tmp_3786;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3786 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPR(gO1,gI2);
   }
   tmp_3785 += tmp_3786;
   result += (-1) * tmp_3785;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh(unsigned gO1) const
{
   std::complex<double> result;

   result += A0(MVWm)*CpUhhbargWmCgWmC(gO1);
   result += A0(MVWm)*CpUhhbargWmgWm(gO1);
   result += A0(MVZ)*CpUhhbargZgZ(gO1);
   result += A0(MVZp)*CpUhhbargZpgZp(gO1);
   result += 4*A0(MVWm)*CpUhhconjVWmVWm(gO1);
   result += 2*A0(MVZp)*CpUhhVZpVZp(gO1);
   result += 2*A0(MVZ)*CpUhhVZVZ(gO1);
   std::complex<double> tmp_3787;
   std::complex<double> tmp_3788;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3788 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_3787 += tmp_3788;
   result += (-1) * tmp_3787;
   std::complex<double> tmp_3789;
   std::complex<double> tmp_3790;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3790 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_3789 += tmp_3790;
   result += (2) * tmp_3789;
   std::complex<double> tmp_3791;
   std::complex<double> tmp_3792;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3792 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_3791 += tmp_3792;
   result += (-0.5) * tmp_3791;
   std::complex<double> tmp_3793;
   std::complex<double> tmp_3794;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3794 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_3793 += tmp_3794;
   result += (-1) * tmp_3793;
   std::complex<double> tmp_3795;
   std::complex<double> tmp_3796;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3796 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_3795 += tmp_3796;
   result += (-0.5) * tmp_3795;
   std::complex<double> tmp_3797;
   std::complex<double> tmp_3798;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3798 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_3797 += tmp_3798;
   result += (6) * tmp_3797;
   std::complex<double> tmp_3799;
   std::complex<double> tmp_3800;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3800 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_3799 += tmp_3800;
   result += (2) * tmp_3799;
   std::complex<double> tmp_3801;
   std::complex<double> tmp_3802;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3802 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_3801 += tmp_3802;
   result += (6) * tmp_3801;
   std::complex<double> tmp_3803;
   std::complex<double> tmp_3804;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3804 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_3803 += tmp_3804;
   result += (-3) * tmp_3803;
   std::complex<double> tmp_3805;
   std::complex<double> tmp_3806;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3806 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_3805 += tmp_3806;
   result += (-1) * tmp_3805;
   std::complex<double> tmp_3807;
   std::complex<double> tmp_3808;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3808 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_3807 += tmp_3808;
   result += (-3) * tmp_3807;
   std::complex<double> tmp_3809;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3809 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_3809;

   return result * oneOver16PiSqr;

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

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
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

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
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

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
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

   Eigen::Array<double,2,1> msf;

   theta = sfermions::diagonalize_sfermions_2x2(sf_data, msf);
   msf1  = msf(0);
   msf2  = msf(1);
}


void CLASSNAME::self_energy_hh_2loop(double result[6]) const
{
   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   double msnusq = Sqr(msnu_2);
   double sxt = Sin(theta_t), cxt = Cos(theta_t);
   double sxb = Sin(theta_b), cxb = Cos(theta_b);
   double sintau = Sin(theta_tau), costau = Cos(theta_tau);

   double gs = g3;
   double as = Sqr(gs) / (4.0 * Pi);
   double rmt = MFu(2);
   double rmtsq = Sqr(rmt);
   double scalesq = Sqr(get_scale());
   double vev2 = Sqr(vd) + Sqr(vu);
   double vev = Sqrt(Sqr(vd) + Sqr(vu));
   double tanb = vu/vd;
   const double tanb2 = Sqr(tanb);
   const double sinb = tanb / Sqrt(1. + tanb2);
   const double cosb = 1. / Sqrt(1. + tanb2);
   double amu = Re(-0.7071067811865475*vS*Lambdax);
   double mg = MGlu;
   double mAsq = (0.7071067811865475*vS*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
   double cotb = 1.0 / tanb;
   double rmb = MFd(2);
   double rmbsq = Sqr(rmb);
   double rmtausq = Sqr(MFe(2));
   double fmasq = Abs(mAsq);
   double lamS = Re(Lambdax);
   static const double root2 = Sqrt(2.0);
   double vevS =  vev / root2;
   double svevS = vS / root2;
   int loop = 2;
   int scheme = 0; // selects DR-bar scheme

   double s11w = 0., s12w = 0., s22w = 0.;
   double s11tau = 0., s12tau = 0., s22tau = 0.;
   double p2w = 0., p2tau = 0.;

   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};
   double DMSB[3][3] = {{ 0. }}, DMPB[3][3] = {{ 0. }};

   LOCK_MUTEX();

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy_higgs_2loop_at_as_nmssm(
         &loop, &rmt, &mg, &mst1sq, &mst2sq, &sxt, &cxt,
         &scalesq, &tanb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy_higgs_2loop_ab_as_nmssm(
         &loop, &rmb, &mg, &msb1sq, &msb2sq, &sxb, &cxb,
         &scalesq, &cotb, &vevS, &lamS, &svevS, &as, &DMSB, &DMPB);
   }

   // Corrections as in MSSM, not corrected for NMSSM,
   // should be OK for MSSM states when S state is close to decoupled

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy_higgs_2loop_at_at_mssm(
         &rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq,
         &msb2sq, &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb,
         &vev2, &s11w, &s12w, &s22w);
      self_energy_pseudoscalar_2loop_at_at_mssm(
         &rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
         &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &p2w);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy_higgs_2loop_atau_atau_mssm(
         &rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
         &costau, &scalesq, &amu, &tanb, &vev2, &scheme, &s11tau,
         &s22tau, &s12tau);
      self_energy_pseudoscalar_2loop_atau_atau_mssm(
         &rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
         &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
   }

   UNLOCK_MUTEX();

   // Make appropriate substitutions for elements following 0907.4682
   // bottom of page 9
   std::swap(DMSB[0][0], DMSB[1][1]);
   std::swap(DMSB[0][2], DMSB[1][2]);

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
         DMS[i][j] += DMSB[i][j];
      }
   }

   const double dMA = p2w + p2tau;

   // subtract two-loop tadpoles
   double tadpole[3];
   tadpole_hh_2loop(tadpole);

   DMS[0][0] += s11w + s11tau + dMA * Sqr(sinb) - tadpole[0] / vd;
   DMS[0][1] += s12w + s12tau - dMA * sinb * cosb;
   DMS[1][1] += s22w + s22tau + dMA * Sqr(cosb) - tadpole[1] / vu;
   DMS[2][2] += - tadpole[2] / vS;

   result[0] = - DMS[0][0]; // 1,1 element
   result[1] = - DMS[0][1]; // 1,2 element
   result[2] = - DMS[0][2]; // 1,3 element
   result[3] = - DMS[1][1]; // 2,2 element
   result[4] = - DMS[1][2]; // 2,3 element
   result[5] = - DMS[2][2]; // 3,3 element

}

void CLASSNAME::self_energy_Ah_2loop(double result[6]) const
{
   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   double msnusq = Sqr(msnu_2);
   double sxt = Sin(theta_t), cxt = Cos(theta_t);
   double sxb = Sin(theta_b), cxb = Cos(theta_b);
   double sintau = Sin(theta_tau), costau = Cos(theta_tau);

   double gs = g3;
   double as = Sqr(gs) / (4.0 * Pi);
   double rmt = MFu(2);
   double rmtsq = Sqr(rmt);
   double scalesq = Sqr(get_scale());
   double vev2 = Sqr(vd) + Sqr(vu);
   double vev = Sqrt(Sqr(vd) + Sqr(vu));
   double tanb = vu/vd;
   const double tanb2 = Sqr(tanb);
   const double sinb = tanb / Sqrt(1. + tanb2);
   const double cosb = 1. / Sqrt(1. + tanb2);
   const double sinb2 = Sqr(sinb);
   const double cosb2 = Sqr(cosb);
   double amu = Re(-0.7071067811865475*vS*Lambdax);
   double mg = MGlu;
   double mAsq = (0.7071067811865475*vS*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
   double cotb = 1.0 / tanb;
   double rmb = MFd(2);
   double rmbsq = Sqr(rmb);
   double rmtausq = Sqr(MFe(2));
   double fmasq = Abs(mAsq);
   double lamS = Re(Lambdax);
   static const double root2 = Sqrt(2.0);
   double vevS =  vev / root2;
   double svevS = vS / root2;
   int loop = 2;

   double p2w = 0., p2tau = 0.;

   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};
   double DMSB[3][3] = {{ 0. }}, DMPB[3][3] = {{ 0. }};

   LOCK_MUTEX();

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy_higgs_2loop_at_as_nmssm(
         &loop, &rmt, &mg, &mst1sq, &mst2sq, &sxt, &cxt,
         &scalesq, &tanb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy_higgs_2loop_ab_as_nmssm(
         &loop, &rmb, &mg, &msb1sq, &msb2sq, &sxb, &cxb,
         &scalesq, &cotb, &vevS, &lamS, &svevS, &as, &DMSB, &DMPB);
   }

   // Corrections as in MSSM, not corrected for NMSSM,
   // should be OK for MSSM states when S state is close to decoupled

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy_pseudoscalar_2loop_at_at_mssm(
         &rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
         &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &p2w);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy_pseudoscalar_2loop_atau_atau_mssm(
         &rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
         &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
   }

   UNLOCK_MUTEX();

   // Make appropriate substitutions for elements following 0907.4682
   // bottom of page 9
   std::swap(DMPB[0][0], DMPB[1][1]);
   std::swap(DMPB[0][2], DMPB[1][2]);

   for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
         DMP[i][j] += DMPB[i][j];
      }
   }

   const double dMA = p2w + p2tau;

   DMP[0][0] += dMA * sinb2;
   DMP[0][1] += dMA * sinb * cosb;
   DMP[1][1] += dMA * cosb2;

   // subtract two-loop tadpoles
   double tadpole[3];
   tadpole_hh_2loop(tadpole);

   DMP[0][0] += - tadpole[0] / vd;
   DMP[1][1] += - tadpole[1] / vu;
   DMP[2][2] += - tadpole[2] / vS;

   result[0] = - DMP[0][0]; // 1,1 element
   result[1] = - DMP[0][1]; // 1,2 element
   result[2] = - DMP[0][2]; // 1,3 element
   result[3] = - DMP[1][1]; // 2,2 element
   result[4] = - DMP[1][2]; // 2,3 element
   result[5] = - DMP[2][2]; // 3,3 element

}



void CLASSNAME::tadpole_hh_2loop(double result[3]) const
{
   // calculate 3rd generation sfermion masses and mixing angles
   double mst_1, mst_2, theta_t;
   double msb_1, msb_2, theta_b;
   double mstau_1, mstau_2, theta_tau;
   double msnu_1, msnu_2, theta_nu;

   calculate_MSu_3rd_generation(mst_1, mst_2, theta_t);
   calculate_MSd_3rd_generation(msb_1, msb_2, theta_b);
   calculate_MSe_3rd_generation(mstau_1, mstau_2, theta_tau);
   calculate_MSv_3rd_generation(msnu_1, msnu_2, theta_nu);

   double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
   double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
   double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
   double msnusq = Sqr(msnu_2);
   double sxt = Sin(theta_t), cxt = Cos(theta_t);
   double sxb = Sin(theta_b), cxb = Cos(theta_b);
   double sintau = Sin(theta_tau), costau = Cos(theta_tau);

   double gs = g3;
   double rmtsq = Sqr(MFu(2));
   double scalesq = Sqr(get_scale());
   double vev2 = Sqr(vd) + Sqr(vu);
   const double vev = Sqrt(vev2);
   double tanb = vu/vd;
   const double tanb2 = Sqr(tanb);
   const double sinb = tanb / Sqrt(1. + tanb2);
   const double cosb = 1. / Sqrt(1. + tanb2);
   double amu = Re(-0.7071067811865475*vS*Lambdax);
   double mg = MGlu;
   double mAsq = (0.7071067811865475*vS*(Sqr(vd) + Sqr(vu))*TLambdax)/(vd*vu);
   double cotbeta = 1.0 / tanb;
   double rmbsq = Sqr(MFd(2));
   double rmtausq = Sqr(MFe(2));

   double s1s = 0., s2s = 0., s1t = 0., s2t = 0.;
   double s1b = 0., s2b = 0., s1tau = 0., s2tau = 0.;

   LOCK_MUTEX();

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      tadpole_higgs_2loop_at_as_mssm(
         &rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq,
         &amu, &tanb, &vev2, &gs, &s1s, &s2s);
   }

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      tadpole_higgs_2loop_at_at_mssm(
         &rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
         &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s1t, &s2t);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      tadpole_higgs_2loop_ab_as_mssm(
         &rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq,
         &amu, &cotbeta, &vev2, &gs, &s2b, &s1b);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      tadpole_higgs_2loop_atau_atau_mssm(
         &rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
         &costau, &scalesq, &amu, &tanb, &vev2, &s1tau, &s2tau);
   }

   UNLOCK_MUTEX();

   // rescale T1 to get TS
   const double sss = s1s * vev * cosb / vS;
   const double ssb = s1b * vev * sinb / vS;

   if (!std::isnan(s1s * s1t * s1b * s1tau * s2s * s2t * s2b * s2tau
                   * sss * ssb)) {
      result[0] = (- s1s - s1t - s1b - s1tau) * vd;
      result[1] = (- s2s - s2t - s2b - s2tau) * vu;
      result[2] = (- sss - ssb) * vS;
   } else {
      result[0] = 0.;
      result[1] = 0.;
      result[2] = 0.;
   }

}


void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MGlu_pole()
{
   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_Glu());
   const double p = MGlu;
   const double self_energy_1  = Re(self_energy_Glu_1(p));
   const double self_energy_PL = Re(self_energy_Glu_PL(p));
   const double self_energy_PR = Re(self_energy_Glu_PR(p));
   const auto M_1loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MGlu) = calculate_singlet_mass(M_1loop);
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
   if (!force_output && problems.is_tachyon(VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_VZ());
   const double p = MVZ;
   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VZ);

   PHYSICAL(MVZ) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MVZp_pole()
{
   if (!force_output && problems.is_tachyon(VZp))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_VZp());
   const double p = MVZp;
   const double self_energy = Re(self_energy_VZp(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VZp);

   PHYSICAL(MVZp) = AbsSqrt(mass_sqr);
}

void CLASSNAME::calculate_MSd_pole()
{
   if (!force_output && problems.is_tachyon(Sd))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,6,6> self_energy;
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sd());

   for (unsigned es = 0; es < 6; ++es) {

      const double p = Abs(MSd(es));
      for (unsigned i1 = 0; i1 < 6; ++i1) {
         for (unsigned i2 = i1; i2 < 6; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Sd(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,6,6> M_1loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZD;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZD,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Sd, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZD);
      #endif

      if (eigen_values(es) < 0.)
         problems.flag_tachyon(Sd);

      PHYSICAL(MSd(es)) = AbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZD) = mix_ZD;
   }
}

void CLASSNAME::calculate_MSv_pole()
{
   if (!force_output && problems.is_tachyon(Sv))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Sv());

   for (unsigned es = 0; es < 3; ++es) {

      const double p = Abs(MSv(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = i1; i2 < 3; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Sv(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree - self_energy);
      Eigen::Array<double,3,1> eigen_values;
      Eigen::Matrix<double,3,3> mix_ZV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZV,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Sv, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZV);
      #endif

      if (eigen_values(es) < 0.)
         problems.flag_tachyon(Sv);

      PHYSICAL(MSv(es)) = AbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZV) = mix_ZV;
   }
}

void CLASSNAME::calculate_MSu_pole()
{
   if (!force_output && problems.is_tachyon(Su))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,6,6> self_energy;
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Su());

   for (unsigned es = 0; es < 6; ++es) {

      const double p = Abs(MSu(es));
      for (unsigned i1 = 0; i1 < 6; ++i1) {
         for (unsigned i2 = i1; i2 < 6; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Su(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,6,6> M_1loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZU;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZU,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Su, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZU);
      #endif

      if (eigen_values(es) < 0.)
         problems.flag_tachyon(Su);

      PHYSICAL(MSu(es)) = AbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZU) = mix_ZU;
   }
}

void CLASSNAME::calculate_MSe_pole()
{
   if (!force_output && problems.is_tachyon(Se))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,6,6> self_energy;
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Se());

   for (unsigned es = 0; es < 6; ++es) {

      const double p = Abs(MSe(es));
      for (unsigned i1 = 0; i1 < 6; ++i1) {
         for (unsigned i2 = i1; i2 < 6; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Se(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,6,6> M_1loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZE;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZE,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Se, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZE);
      #endif

      if (eigen_values(es) < 0.)
         problems.flag_tachyon(Se);

      PHYSICAL(MSe(es)) = AbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZE) = mix_ZE;
   }
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
      Eigen::Matrix<double,3,3> self_energy;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_hh());

      // two-loop Higgs self-energy contributions
      double two_loop[6] = { 0. };
      if (pole_mass_loop_order > 1)
         self_energy_hh_2loop(two_loop);
         for (unsigned i = 0; i < 6; i++) {
            if (!std::isfinite(two_loop[i])) {
               two_loop[i] = 0.;
               problems.flag_bad_mass(UMSSM_info::hh);
            }
         }

      for (unsigned es = 0; es < 3; ++es) {

         const double p = Abs(old_Mhh(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_hh(p,i1,i2
                  ));
            }
         }

         self_energy(0, 0) += two_loop[0];
         self_energy(0, 1) += two_loop[1];
         self_energy(0, 2) += two_loop[2];
         self_energy(1, 1) += two_loop[3];
         self_energy(1, 2) += two_loop[4];
         self_energy(2, 2) += two_loop[5];

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(UMSSM_info::hh,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZH);
         #endif

         if (eigen_values(es) < 0.)
            problems.flag_tachyon(hh);

         PHYSICAL(Mhh(es)) = AbsSqrt(eigen_values(es));
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
   if (!force_output && problems.is_tachyon(Ah))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MAh) old_MAh(MAh), new_MAh(MAh);

   do {
      Eigen::Matrix<double,3,3> self_energy;
      const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Ah());

      // two-loop Higgs self-energy contributions
      double two_loop[6] = { 0. };
      if (pole_mass_loop_order > 1)
         self_energy_Ah_2loop(two_loop);
         for (unsigned i = 0; i < 6; i++) {
            if (!std::isfinite(two_loop[i])) {
               two_loop[i] = 0.;
               problems.flag_bad_mass(UMSSM_info::Ah);
            }
         }

      for (unsigned es = 0; es < 3; ++es) {
         // skip goldstone bosons
         if (is_equal_rel(MAh(es), MVZ, 1e-10)) {
            PHYSICAL(MAh(es)) = MVZ;
            continue;
         }
         if (is_equal_rel(MAh(es), MVZp, 1e-10)) {
            PHYSICAL(MAh(es)) = MVZp;
            continue;
         }

         const double p = Abs(old_MAh(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Ah(p,i1,i2
                  ));
            }
         }

         self_energy(0, 0) += two_loop[0];
         self_energy(0, 1) += two_loop[1];
         self_energy(0, 2) += two_loop[2];
         self_energy(1, 1) += two_loop[3];
         self_energy(1, 2) += two_loop[4];
         self_energy(2, 2) += two_loop[5];

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(UMSSM_info::Ah,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZA);
         #endif

         if (eigen_values(es) < 0.)
            problems.flag_tachyon(Ah);

         PHYSICAL(MAh(es)) = AbsSqrt(eigen_values(es));
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
   if (!force_output && problems.is_tachyon(Hpm))
      return;

   // diagonalization with high precision
   unsigned iteration = 0;
   double diff = 0.0;
   decltype(MHpm) old_MHpm(MHpm), new_MHpm(MHpm);

   do {
      Eigen::Matrix<double,2,2> self_energy;
      const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Hpm());

      for (unsigned es = 0; es < 2; ++es) {
         // skip goldstone bosons
         if (is_equal_rel(MHpm(es), MVWm, 1e-10)) {
            PHYSICAL(MHpm(es)) = MVWm;
            continue;
         }

         const double p = Abs(old_MHpm(es));
         for (unsigned i1 = 0; i1 < 2; ++i1) {
            for (unsigned i2 = i1; i2 < 2; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Hpm(p,i1,
                  i2));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,2,2> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,2,1> eigen_values;
         Eigen::Matrix<double,2,2> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZP, eigenvalue_error);
            problems.flag_bad_mass(UMSSM_info::Hpm,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZP);
         #endif

         if (eigen_values(es) < 0.)
            problems.flag_tachyon(Hpm);

         PHYSICAL(MHpm(es)) = AbsSqrt(eigen_values(es));
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
   Eigen::Matrix<double,6,6> self_energy_1;
   Eigen::Matrix<double,6,6> self_energy_PL;
   Eigen::Matrix<double,6,6> self_energy_PR;
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Chi());
   for (unsigned es = 0; es < 6; ++es) {
      const double p = Abs(MChi(es));
      for (unsigned i1 = 0; i1 < 6; ++i1) {
         for (unsigned i2 = 0; i2 < 6; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Chi_1(p,i1,i2
               ));
            self_energy_PL(i1,i2) = Re(self_energy_Chi_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Chi_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,6,6> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,6,6> M_1loop(M_tree + 0.5 * (delta_M
         + delta_M.transpose()));
      Eigen::Array<double,6,1> eigen_values;
      decltype(ZN) mix_ZN;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_1loop, eigen_values, mix_ZN,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Chi, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_symmetric(M_1loop, eigen_values, mix_ZN);
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
      const Eigen::Matrix<double,2,2> M_1loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM) mix_UM;
      decltype(UP) mix_UP;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_UM, mix_UP, eigenvalue_error);
      problems.flag_bad_mass(UMSSM_info::Cha, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_UM, mix_UP);
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
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZEL) mix_ZEL;
      decltype(ZER) mix_ZER;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZEL, mix_ZER, eigenvalue_error
         );
      problems.flag_bad_mass(UMSSM_info::Fe, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_ZEL, mix_ZER);
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
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZDL) mix_ZDL;
      decltype(ZDR) mix_ZDR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZDL, mix_ZDR, eigenvalue_error
         );
      problems.flag_bad_mass(UMSSM_info::Fd, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_ZDL, mix_ZDR);
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
   const bool add_2loop_corrections = pole_mass_loop_order > 1 &&
      TOP_2LOOP_CORRECTION_QCD;
   const double currentScale = get_scale();

   const double qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(2))
      /Sqr(currentScale)))*Sqr(g3);

   double qcd_2l = 0.;

   if (add_2loop_corrections) {
      qcd_2l = -0.005191204615668296*Power(g3,4) +
         0.0032883224409535764*Power(g3,4)*Log(Sqr(MFu(2))/Sqr(currentScale)) -
         0.0008822328500119351*Power(g3,4)*Sqr(Log(Power(MFu(2),2)/Sqr(
         currentScale)));
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
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZUL) mix_ZUL;
      decltype(ZUR) mix_ZUR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZUL, mix_ZUR, eigenvalue_error
         );
      problems.flag_bad_mass(UMSSM_info::Fu, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_ZUL, mix_ZUR);
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
   if (!force_output && problems.is_tachyon(VWm))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_VWm());
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
   const double mass_sqr = get_mass_matrix_VWm() - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VWm);

   return AbsSqrt(mass_sqr);
}

double CLASSNAME::calculate_MVZ_pole(double p)
{
   if (!force_output && problems.is_tachyon(VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ(p));
   const double mass_sqr = get_mass_matrix_VZ() - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(VZ);

   return AbsSqrt(mass_sqr);
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
   const double qcd_1l = -0.008443431970194815*(5. - 3.*Log(Sqr(MFu(2))
      /Sqr(currentScale)))*Sqr(g3);
   const double qcd_2l = -0.003408916029785599*Power(g3,4) +
      0.0011495761378943394*Power(g3,4)*Log(Sqr(MFu(2))/Sqr(currentScale)) -
      0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(MFu(2),2)/Sqr(
      currentScale)));

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l);

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
   const double m_tree = MFd(2);
   const double drbar_conversion = 1 - 0.00020496318737651018*Power(g3,4)
      + 0.0006860288475783287*Sqr(g1) + 0.0023747152416172916*Sqr(g2) -
      0.008443431970194815*Sqr(g3);
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar / (1.0 - self_energy_1/m_tree -
      self_energy_PL - self_energy_PR);

   return m_susy_drbar;
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
   const double drbar_conversion = 1 - 0.0023747152416172916*(0.6*Sqr(g1)
      - Sqr(g2));
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;

   const double m_susy_drbar = m_sm_drbar + self_energy_1 + m_sm_drbar *
      (self_energy_PL + self_energy_PR);

   return m_susy_drbar;
}

double CLASSNAME::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
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


double CLASSNAME::ThetaW() const
{
   return ArcTan((0.7745966692414834*g1)/g2);
}

double CLASSNAME::v() const
{
   return 2*Sqrt(Sqr(MVWm)/Sqr(g2));
}

double CLASSNAME::ThetaWp() const
{
   return 0;
}


std::ostream& operator<<(std::ostream& ostr, const UMSSM_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
