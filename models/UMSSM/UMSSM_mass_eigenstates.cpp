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

// File generated at Fri 8 Jan 2016 15:18:34

/**
 * @file UMSSM_mass_eigenstates.cpp
 * @brief implementation of the UMSSM model class
 *
 * Contains the definition of the UMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Fri 8 Jan 2016 15:18:34 with FlexibleSUSY
 * 1.3.1 (git commit: v1.3.1) and SARAH 4.6.0 .
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
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS  1

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
   , MVG(0), MGlu(0), MVP(0), MVZ(0), MVZp(0), MSd(Eigen::Array<double,6,1>
      ::Zero()), MSv(Eigen::Array<double,6,1>::Zero()), MSu(Eigen::Array<double,6,
      1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,
      3,1>::Zero()), MAh(Eigen::Array<double,3,1>::Zero()), MHpm(Eigen::Array<
      double,2,1>::Zero()), MChi(Eigen::Array<double,6,1>::Zero()), MFv(
      Eigen::Array<double,3,1>::Zero()), MCha(Eigen::Array<double,2,1>::Zero()),
      MFe(Eigen::Array<double,3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero())
      , MFu(Eigen::Array<double,3,1>::Zero()), MVWm(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,6,6>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,6,6>::Zero()), ZVL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZVR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), UM(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3,3>::Zero())

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

   std::for_each(solvers, solvers + number_of_solvers, Delete_object());

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
   ostr << "MFv = " << MFv.transpose() << '\n';
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
   calculate_MSd();
   calculate_MSv();
   calculate_MSu();
   calculate_MSe();
   calculate_Mhh();
   calculate_MAh();
   calculate_MHpm();
   calculate_MChi();
   calculate_MFv();
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
      std::thread thread_MVP(Thread(this, &CLASSNAME::calculate_MVP_pole));
      std::thread thread_MVZ(Thread(this, &CLASSNAME::calculate_MVZ_pole));
      std::thread thread_MFv(Thread(this, &CLASSNAME::calculate_MFv_pole));
      std::thread thread_MFe(Thread(this, &CLASSNAME::calculate_MFe_pole));
      std::thread thread_MFd(Thread(this, &CLASSNAME::calculate_MFd_pole));
      std::thread thread_MFu(Thread(this, &CLASSNAME::calculate_MFu_pole));
      std::thread thread_MVWm(Thread(this, &CLASSNAME::calculate_MVWm_pole));
      thread_MVG.join();
      thread_MVP.join();
      thread_MVZ.join();
      thread_MFv.join();
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
   if (PHYSICAL(MSd).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Sd);
   if (PHYSICAL(MSv).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Sv);
   if (PHYSICAL(MSu).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Su);
   if (PHYSICAL(MSe).tail<6>().minCoeff() < 0.) problems.flag_tachyon(Se);
   if (PHYSICAL(Mhh).tail<3>().minCoeff() < 0.) problems.flag_tachyon(hh);
   if (PHYSICAL(MAh).tail<1>().minCoeff() < 0.) problems.flag_tachyon(Ah);
   if (PHYSICAL(MHpm).tail<1>().minCoeff() < 0.) problems.flag_tachyon(Hpm);

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
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;
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

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sv() const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qv = LOCALINPUT(Qv);

   Eigen::Matrix<double,6,6> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(0,0)) + AbsSqr(Yv(0,1)) + AbsSqr(Yv(0,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(1,0) +
      Conj(Yv(0,1))*Yv(1,1) + Conj(Yv(0,2))*Yv(1,2));
   mass_matrix_Sv(0,2) = ml2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,0))*Yv(2,0) +
      Conj(Yv(0,1))*Yv(2,1) + Conj(Yv(0,2))*Yv(2,2));
   mass_matrix_Sv(0,3) = 0.7071067811865475*vu*Conj(TYv(0,0)) - 0.5*vd*vS
      *Conj(Yv(0,0))*Lambdax;
   mass_matrix_Sv(0,4) = 0.7071067811865475*vu*Conj(TYv(0,1)) - 0.5*vd*vS
      *Conj(Yv(0,1))*Lambdax;
   mass_matrix_Sv(0,5) = 0.7071067811865475*vu*Conj(TYv(0,2)) - 0.5*vd*vS
      *Conj(Yv(0,2))*Lambdax;
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(1,0)) + AbsSqr(Yv(1,1)) + AbsSqr(Yv(1,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(1,0))*Yv(2,0) +
      Conj(Yv(1,1))*Yv(2,1) + Conj(Yv(1,2))*Yv(2,2));
   mass_matrix_Sv(1,3) = 0.7071067811865475*vu*Conj(TYv(1,0)) - 0.5*vd*vS
      *Conj(Yv(1,0))*Lambdax;
   mass_matrix_Sv(1,4) = 0.7071067811865475*vu*Conj(TYv(1,1)) - 0.5*vd*vS
      *Conj(Yv(1,1))*Lambdax;
   mass_matrix_Sv(1,5) = 0.7071067811865475*vu*Conj(TYv(1,2)) - 0.5*vd*vS
      *Conj(Yv(1,2))*Lambdax;
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*QHd*Ql*Sqr(gp)*Sqr(vd) + 0.5*Ql*Qs*Sqr(gp)*Sqr(vS) + 0.5*(
      AbsSqr(Yv(2,0)) + AbsSqr(Yv(2,1)) + AbsSqr(Yv(2,2)))*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.5*QHu*Ql*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(2,3) = 0.7071067811865475*vu*Conj(TYv(2,0)) - 0.5*vd*vS
      *Conj(Yv(2,0))*Lambdax;
   mass_matrix_Sv(2,4) = 0.7071067811865475*vu*Conj(TYv(2,1)) - 0.5*vd*vS
      *Conj(Yv(2,1))*Lambdax;
   mass_matrix_Sv(2,5) = 0.7071067811865475*vu*Conj(TYv(2,2)) - 0.5*vd*vS
      *Conj(Yv(2,2))*Lambdax;
   mass_matrix_Sv(3,3) = mvR2(0,0) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qv*
      Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,0)) +
      AbsSqr(Yv(1,0)) + AbsSqr(Yv(2,0)))*Sqr(vu) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*
      QHu*Qv*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(3,4) = mvR2(0,1) + 0.5*Sqr(vu)*(Conj(Yv(0,1))*Yv(0,0) +
      Conj(Yv(1,1))*Yv(1,0) + Conj(Yv(2,1))*Yv(2,0));
   mass_matrix_Sv(3,5) = mvR2(0,2) + 0.5*Sqr(vu)*(Conj(Yv(0,2))*Yv(0,0) +
      Conj(Yv(1,2))*Yv(1,0) + Conj(Yv(2,2))*Yv(2,0));
   mass_matrix_Sv(4,4) = mvR2(1,1) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qv*
      Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,1)) +
      AbsSqr(Yv(1,1)) + AbsSqr(Yv(2,1)))*Sqr(vu) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*
      QHu*Qv*Sqr(gp)*Sqr(vu);
   mass_matrix_Sv(4,5) = mvR2(1,2) + 0.5*Sqr(vu)*(Conj(Yv(0,2))*Yv(0,1) +
      Conj(Yv(1,2))*Yv(1,1) + Conj(Yv(2,2))*Yv(2,1));
   mass_matrix_Sv(5,5) = mvR2(2,2) - 0.15*Sqr(g1)*Sqr(vd) + 0.5*QHd*Qv*
      Sqr(gp)*Sqr(vd) + 0.5*Qs*Qv*Sqr(gp)*Sqr(vS) + 0.5*(AbsSqr(Yv(0,2)) +
      AbsSqr(Yv(1,2)) + AbsSqr(Yv(2,2)))*Sqr(vu) + 0.15*Sqr(g1)*Sqr(vu) + 0.5*
      QHu*Qv*Sqr(gp)*Sqr(vu);

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
   problems.flag_bad_mass(UMSSM_info::Fv, eigenvalue_error > precision *
      Abs(MFv(0)));
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

std::complex<double> CLASSNAME::CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_54;
   std::complex<double> tmp_55;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_55 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_54 += tmp_55;
   result += (0.1*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_54;
   std::complex<double> tmp_56;
   std::complex<double> tmp_57;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_57 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_56 += tmp_57;
   result += (-(Qd*QHd*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_56;
   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_59 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_58 += tmp_59;
   result += (-0.1*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_58;
   std::complex<double> tmp_60;
   std::complex<double> tmp_61;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_61 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_60 += tmp_61;
   result += (-(Qd*QHu*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_60;
   std::complex<double> tmp_62;
   std::complex<double> tmp_63;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_63 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_62 += tmp_63;
   result += (-(Qd*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_62;
   std::complex<double> tmp_64;
   std::complex<double> tmp_65;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_66;
      std::complex<double> tmp_67;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_68;
         std::complex<double> tmp_69;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_69 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_68 += tmp_69;
         tmp_67 += (KroneckerDelta(gO2,3 + j2)) * tmp_68;
      }
      tmp_66 += tmp_67;
      tmp_65 += (KroneckerDelta(gO1,3 + j3)) * tmp_66;
   }
   tmp_64 += tmp_65;
   result += (-(Conj(ZH(gI1,0))*Conj(ZH(gI2,0)))) * tmp_64;
   if (gO1 < 3) {
      result += 0.05*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -0.05*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Qq*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_70;
      std::complex<double> tmp_71;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_71 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_70 += tmp_71;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))) * tmp_70
         ;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_72;
      std::complex<double> tmp_73;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_73 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_72 += tmp_73;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,1))*Conj(ZH(gI2,2))) * tmp_72
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_74;
      std::complex<double> tmp_75;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_75 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_74 += tmp_75;
      result += (0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))*Lambdax) * tmp_74;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_76;
      std::complex<double> tmp_77;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_77 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_76 += tmp_77;
      result += (0.5*Conj(ZH(gI1,1))*Conj(ZH(gI2,2))*Lambdax) * tmp_76;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_78;
      std::complex<double> tmp_79;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_79 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_78 += tmp_79;
      result += (-(Conj(ZH(gI1,0))*Conj(ZH(gI2,0)))) * tmp_78;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_80;
      std::complex<double> tmp_81;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_81 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_80 += tmp_81;
      result += (UP(gI2,1)) * tmp_80;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_82;
   std::complex<double> tmp_83;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_84;
      std::complex<double> tmp_85;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_85 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_84 += tmp_85;
      tmp_83 += (Conj(ZUL(gI1,j2))) * tmp_84;
   }
   tmp_82 += tmp_83;
   result += (Conj(UM(gI2,1))) * tmp_82;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(ZUL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_86;
   std::complex<double> tmp_87;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_87 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_86 += tmp_87;
   result += (-1.4142135623730951*gp*Qd*ZN(gI2,0)) * tmp_86;
   std::complex<double> tmp_88;
   std::complex<double> tmp_89;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_89 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_88 += tmp_89;
   result += (-0.3651483716701107*g1*ZN(gI2,1)) * tmp_88;
   if (gO2 < 3) {
      std::complex<double> tmp_90;
      std::complex<double> tmp_91;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_91 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_90 += tmp_91;
      result += (-ZN(gI2,3)) * tmp_90;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_92;
   std::complex<double> tmp_93;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_94;
      std::complex<double> tmp_95;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_95 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_94 += tmp_95;
      tmp_93 += (Conj(ZDL(gI1,j2))) * tmp_94;
   }
   tmp_92 += tmp_93;
   result += (-Conj(ZN(gI2,3))) * tmp_92;
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

   std::complex<double> tmp_96;
   std::complex<double> tmp_98;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_98 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_96 += tmp_98;
   std::complex<double> tmp_97;
   std::complex<double> tmp_99;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_99 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_97 += tmp_99;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_96 * tmp_97;
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
   result += (-0.6666666666666666*Sqr(g3)) * tmp_100 * tmp_101;
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
   result += (-0.5*Sqr(gp)*Sqr(Qd)) * tmp_104 * tmp_105;
   std::complex<double> tmp_108;
   std::complex<double> tmp_110;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_110 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_108 += tmp_110;
   std::complex<double> tmp_109;
   std::complex<double> tmp_111;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_111 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_109 += tmp_111;
   result += (-0.05*Sqr(g1)) * tmp_108 * tmp_109;
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
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_112 * tmp_113;
   std::complex<double> tmp_116;
   std::complex<double> tmp_118;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_118 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_116 += tmp_118;
   std::complex<double> tmp_117;
   std::complex<double> tmp_119;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_119 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_117 += tmp_119;
   result += (-0.1*Sqr(g1)) * tmp_116 * tmp_117;
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
   result += (-1.5*Sqr(gp)*Sqr(Qd)) * tmp_120 * tmp_121;
   std::complex<double> tmp_124;
   std::complex<double> tmp_126;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_126 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_124 += tmp_126;
   std::complex<double> tmp_125;
   std::complex<double> tmp_127;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_127 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_125 += tmp_127;
   result += (-0.05*Sqr(g1)) * tmp_124 * tmp_125;
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
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_128 * tmp_129;
   std::complex<double> tmp_132;
   std::complex<double> tmp_134;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_134 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_132 += tmp_134;
   std::complex<double> tmp_133;
   std::complex<double> tmp_135;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_135 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_133 += tmp_135;
   result += (-0.1*Sqr(g1)) * tmp_132 * tmp_133;
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
   result += (-1.5*Sqr(gp)*Sqr(Qd)) * tmp_136 * tmp_137;
   std::complex<double> tmp_140;
   std::complex<double> tmp_142;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_142 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_140 += tmp_142;
   std::complex<double> tmp_141;
   std::complex<double> tmp_143;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_143 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_141 += tmp_143;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_140 * tmp_141;
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
   result += (-0.6666666666666666*Sqr(g3)) * tmp_144 * tmp_145;
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
   result += (-0.5*Sqr(gp)*Sqr(Qd)) * tmp_148 * tmp_149;
   std::complex<double> tmp_152;
   std::complex<double> tmp_154;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_155;
      std::complex<double> tmp_156;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_156 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_155 += tmp_156;
      tmp_154 += (Conj(ZD(gI2,j2))) * tmp_155;
   }
   tmp_152 += tmp_154;
   std::complex<double> tmp_153;
   std::complex<double> tmp_157;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_158;
      std::complex<double> tmp_159;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_159 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_158 += tmp_159;
      tmp_157 += (ZD(gI1,j4)) * tmp_158;
   }
   tmp_153 += tmp_157;
   result += (-1) * tmp_152 * tmp_153;
   if (gO1 < 3) {
      std::complex<double> tmp_160;
      std::complex<double> tmp_161;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_161 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_160 += tmp_161;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_160;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_162;
      std::complex<double> tmp_163;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_163 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_162 += tmp_163;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_162;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_164;
      std::complex<double> tmp_165;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_165 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_164 += tmp_165;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_164;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_166;
      std::complex<double> tmp_167;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_167 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_166 += tmp_167;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_166;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_168;
      std::complex<double> tmp_169;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_169 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_168 += tmp_169;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_168;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_170;
      std::complex<double> tmp_171;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_171 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_170 += tmp_171;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_170;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_172;
      std::complex<double> tmp_173;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_173 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_172 += tmp_173;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_172;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_174;
      std::complex<double> tmp_175;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_175 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_174 += tmp_175;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_174;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_176;
      std::complex<double> tmp_177;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_177 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_176 += tmp_177;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_176;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_178;
      std::complex<double> tmp_179;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_179 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_178 += tmp_179;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_178;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_182;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_182 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_180 += tmp_182;
      std::complex<double> tmp_181;
      std::complex<double> tmp_183;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_184;
         std::complex<double> tmp_185;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_185 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_184 += tmp_185;
         tmp_183 += (ZD(gI1,j4)) * tmp_184;
      }
      tmp_181 += tmp_183;
      result += (-3) * tmp_180 * tmp_181;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_186;
      std::complex<double> tmp_187;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_187 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_186 += tmp_187;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_186;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_188;
      std::complex<double> tmp_189;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_189 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_188 += tmp_189;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_188;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_190;
      std::complex<double> tmp_191;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_191 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_190 += tmp_191;
      result += (-0.5*Qd*Qq*Conj(ZD(gI2,gO2))*Sqr(gp)) * tmp_190;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_192;
      std::complex<double> tmp_193;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_193 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_192 += tmp_193;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_192;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_194;
      std::complex<double> tmp_195;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_195 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_194 += tmp_195;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_194;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_196;
      std::complex<double> tmp_197;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_197 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_196 += tmp_197;
      result += (-0.5*Qd*Qq*Conj(ZD(gI2,gO2))*Sqr(gp)) * tmp_196;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_198;
      std::complex<double> tmp_200;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_201;
         std::complex<double> tmp_202;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_202 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_201 += tmp_202;
         tmp_200 += (Conj(ZD(gI2,j2))) * tmp_201;
      }
      tmp_198 += tmp_200;
      std::complex<double> tmp_199;
      std::complex<double> tmp_203;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_203 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_199 += tmp_203;
      result += (-3) * tmp_198 * tmp_199;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_204;
      std::complex<double> tmp_206;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_206 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_204 += tmp_206;
      std::complex<double> tmp_205;
      std::complex<double> tmp_207;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_207 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_205 += tmp_207;
      result += (-1) * tmp_204 * tmp_205;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_208;
      std::complex<double> tmp_209;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_209 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_208 += tmp_209;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_208;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_210;
      std::complex<double> tmp_211;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_211 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_210 += tmp_211;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_210;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_212;
      std::complex<double> tmp_213;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_213 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_212 += tmp_213;
      result += (-0.5*Qd*Qq*Sqr(gp)*ZD(gI1,gO1)) * tmp_212;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_214;
      std::complex<double> tmp_215;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_215 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_214 += tmp_215;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_214;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_216;
      std::complex<double> tmp_217;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_217 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_216 += tmp_217;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_216;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_218;
      std::complex<double> tmp_219;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_219 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_218 += tmp_219;
      result += (-0.5*Qd*Qq*Sqr(gp)*ZD(gI1,gO1)) * tmp_218;
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

   std::complex<double> tmp_220;
   std::complex<double> tmp_222;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_222 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_220 += tmp_222;
   std::complex<double> tmp_221;
   std::complex<double> tmp_223;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_223 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_221 += tmp_223;
   result += (0.05*Sqr(g1)) * tmp_220 * tmp_221;
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
   result += (-0.5*Qd*Ql*Sqr(gp)) * tmp_224 * tmp_225;
   std::complex<double> tmp_228;
   std::complex<double> tmp_230;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_230 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_228 += tmp_230;
   std::complex<double> tmp_229;
   std::complex<double> tmp_231;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_231 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_229 += tmp_231;
   result += (-0.1*Sqr(g1)) * tmp_228 * tmp_229;
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
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_232 * tmp_233;
   std::complex<double> tmp_236;
   std::complex<double> tmp_238;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_238 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_236 += tmp_238;
   std::complex<double> tmp_237;
   std::complex<double> tmp_239;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_239 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_237 += tmp_239;
   result += (0.05*Sqr(g1)) * tmp_236 * tmp_237;
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
   result += (-0.5*Qd*Ql*Sqr(gp)) * tmp_240 * tmp_241;
   std::complex<double> tmp_244;
   std::complex<double> tmp_246;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_246 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_244 += tmp_246;
   std::complex<double> tmp_245;
   std::complex<double> tmp_247;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_247 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_245 += tmp_247;
   result += (-0.1*Sqr(g1)) * tmp_244 * tmp_245;
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
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_248 * tmp_249;
   if (gO1 < 3) {
      std::complex<double> tmp_252;
      std::complex<double> tmp_253;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_253 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_252 += tmp_253;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_252;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_254;
      std::complex<double> tmp_255;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_255 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_254 += tmp_255;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_254;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_256;
      std::complex<double> tmp_257;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_257 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_256 += tmp_257;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_256;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_258;
      std::complex<double> tmp_259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_259 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_258 += tmp_259;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_258;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_260;
      std::complex<double> tmp_261;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_261 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_260 += tmp_261;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_260;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_262;
      std::complex<double> tmp_263;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_263 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_262 += tmp_263;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_262;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_264;
      std::complex<double> tmp_265;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_265 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_264 += tmp_265;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_264;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_266;
      std::complex<double> tmp_267;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_267 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_266 += tmp_267;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_266;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_268;
      std::complex<double> tmp_269;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_269 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_268 += tmp_269;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_268;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_270;
      std::complex<double> tmp_271;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_271 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_270 += tmp_271;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_270;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_272;
      std::complex<double> tmp_274;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_274 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_272 += tmp_274;
      std::complex<double> tmp_273;
      std::complex<double> tmp_275;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_276;
         std::complex<double> tmp_277;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_277 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_276 += tmp_277;
         tmp_275 += (ZE(gI1,j4)) * tmp_276;
      }
      tmp_273 += tmp_275;
      result += (-1) * tmp_272 * tmp_273;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_278;
      std::complex<double> tmp_280;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_281;
         std::complex<double> tmp_282;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_282 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_281 += tmp_282;
         tmp_280 += (Conj(ZE(gI2,j2))) * tmp_281;
      }
      tmp_278 += tmp_280;
      std::complex<double> tmp_279;
      std::complex<double> tmp_283;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_283 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_279 += tmp_283;
      result += (-1) * tmp_278 * tmp_279;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_284;
   std::complex<double> tmp_286;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_286 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_284 += tmp_286;
   std::complex<double> tmp_285;
   std::complex<double> tmp_287;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_287 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_285 += tmp_287;
   result += (-0.05*Sqr(g1)) * tmp_284 * tmp_285;
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
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_288 * tmp_289;
   std::complex<double> tmp_292;
   std::complex<double> tmp_294;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_294 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_292 += tmp_294;
   std::complex<double> tmp_293;
   std::complex<double> tmp_295;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_295 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_293 += tmp_295;
   result += (0.2*Sqr(g1)) * tmp_292 * tmp_293;
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
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_296 * tmp_297;
   std::complex<double> tmp_300;
   std::complex<double> tmp_302;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_302 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_300 += tmp_302;
   std::complex<double> tmp_301;
   std::complex<double> tmp_303;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_303 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_301 += tmp_303;
   result += (-0.05*Sqr(g1)) * tmp_300 * tmp_301;
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
   result += (-1.5*Qd*Qq*Sqr(gp)) * tmp_304 * tmp_305;
   std::complex<double> tmp_308;
   std::complex<double> tmp_310;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_310 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_308 += tmp_310;
   std::complex<double> tmp_309;
   std::complex<double> tmp_311;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_311 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_309 += tmp_311;
   result += (0.2*Sqr(g1)) * tmp_308 * tmp_309;
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
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_312 * tmp_313;
   std::complex<double> tmp_316;
   std::complex<double> tmp_318;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_319;
      std::complex<double> tmp_320;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_320 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_319 += tmp_320;
      tmp_318 += (Conj(ZU(gI2,j2))) * tmp_319;
   }
   tmp_316 += tmp_318;
   std::complex<double> tmp_317;
   std::complex<double> tmp_321;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_322;
      std::complex<double> tmp_323;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_323 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_322 += tmp_323;
      tmp_321 += (ZU(gI1,j4)) * tmp_322;
   }
   tmp_317 += tmp_321;
   result += (-1) * tmp_316 * tmp_317;
   if (gO1 < 3) {
      std::complex<double> tmp_324;
      std::complex<double> tmp_325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_325 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_324 += tmp_325;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_324;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_326;
      std::complex<double> tmp_327;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_327 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_326 += tmp_327;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_326;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_328;
      std::complex<double> tmp_329;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_329 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_328 += tmp_329;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_328;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_330;
      std::complex<double> tmp_331;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_331 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_330 += tmp_331;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_330;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_332;
      std::complex<double> tmp_333;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_333 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_332 += tmp_333;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_332;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_334;
      std::complex<double> tmp_335;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_335 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_334 += tmp_335;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_334;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_336;
      std::complex<double> tmp_337;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_337 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_336 += tmp_337;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_336;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_338;
      std::complex<double> tmp_339;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_339 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_338 += tmp_339;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_338;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_340;
      std::complex<double> tmp_341;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_341 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_340 += tmp_341;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_340;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_342;
      std::complex<double> tmp_343;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_343 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_342 += tmp_343;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_342;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_344;
      std::complex<double> tmp_346;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_346 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_344 += tmp_346;
      std::complex<double> tmp_345;
      std::complex<double> tmp_347;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_347 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_345 += tmp_347;
      result += (-1) * tmp_344 * tmp_345;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_348;
   std::complex<double> tmp_350;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_350 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_348 += tmp_350;
   std::complex<double> tmp_349;
   std::complex<double> tmp_351;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_351 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_349 += tmp_351;
   result += (0.05*Sqr(g1)) * tmp_348 * tmp_349;
   std::complex<double> tmp_352;
   std::complex<double> tmp_354;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_354 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_352 += tmp_354;
   std::complex<double> tmp_353;
   std::complex<double> tmp_355;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_355 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_353 += tmp_355;
   result += (-0.5*Qd*Ql*Sqr(gp)) * tmp_352 * tmp_353;
   std::complex<double> tmp_356;
   std::complex<double> tmp_358;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_358 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_356 += tmp_358;
   std::complex<double> tmp_357;
   std::complex<double> tmp_359;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_359 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_357 += tmp_359;
   result += (-0.1*Sqr(g1)) * tmp_356 * tmp_357;
   std::complex<double> tmp_360;
   std::complex<double> tmp_362;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_362 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_360 += tmp_362;
   std::complex<double> tmp_361;
   std::complex<double> tmp_363;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_363 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_361 += tmp_363;
   result += (-0.5*Qd*Qv*Sqr(gp)) * tmp_360 * tmp_361;
   std::complex<double> tmp_364;
   std::complex<double> tmp_366;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_366 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_364 += tmp_366;
   std::complex<double> tmp_365;
   std::complex<double> tmp_367;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_367 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_365 += tmp_367;
   result += (0.05*Sqr(g1)) * tmp_364 * tmp_365;
   std::complex<double> tmp_368;
   std::complex<double> tmp_370;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_370 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_368 += tmp_370;
   std::complex<double> tmp_369;
   std::complex<double> tmp_371;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_371 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_369 += tmp_371;
   result += (-0.5*Qd*Ql*Sqr(gp)) * tmp_368 * tmp_369;
   std::complex<double> tmp_372;
   std::complex<double> tmp_374;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_374 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_372 += tmp_374;
   std::complex<double> tmp_373;
   std::complex<double> tmp_375;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_375 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_373 += tmp_375;
   result += (-0.1*Sqr(g1)) * tmp_372 * tmp_373;
   std::complex<double> tmp_376;
   std::complex<double> tmp_378;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_378 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_376 += tmp_378;
   std::complex<double> tmp_377;
   std::complex<double> tmp_379;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_379 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_377 += tmp_379;
   result += (-0.5*Qd*Qv*Sqr(gp)) * tmp_376 * tmp_377;
   if (gO1 < 3) {
      std::complex<double> tmp_380;
      std::complex<double> tmp_381;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_381 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_380 += tmp_381;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_380;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_382;
      std::complex<double> tmp_383;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_383 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_382 += tmp_383;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_382;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_384;
      std::complex<double> tmp_385;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_385 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_384 += tmp_385;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_384;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_386;
      std::complex<double> tmp_387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_387 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_386 += tmp_387;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_386;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_388;
      std::complex<double> tmp_389;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_389 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_388 += tmp_389;
      result += (-0.5*Qq*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_388;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_390;
      std::complex<double> tmp_391;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_391 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_390 += tmp_391;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_390;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_392;
      std::complex<double> tmp_393;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_393 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_392 += tmp_393;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_392;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_394;
      std::complex<double> tmp_395;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_395 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_394 += tmp_395;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_394;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_396;
      std::complex<double> tmp_397;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_397 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_396 += tmp_397;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_396;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_398;
      std::complex<double> tmp_399;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_399 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_398 += tmp_399;
      result += (-0.5*Qq*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_398;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_400;
   std::complex<double> tmp_401;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_402;
      std::complex<double> tmp_403;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_403 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_402 += tmp_403;
      tmp_401 += (Conj(ZU(gI1,j2))) * tmp_402;
   }
   tmp_400 += tmp_401;
   result += (ZP(gI2,0)) * tmp_400;
   std::complex<double> tmp_404;
   std::complex<double> tmp_405;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_406;
      std::complex<double> tmp_407;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_408;
         std::complex<double> tmp_409;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_409 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_408 += tmp_409;
         tmp_407 += (KroneckerDelta(gO2,3 + j2)) * tmp_408;
      }
      tmp_406 += tmp_407;
      tmp_405 += (Conj(ZU(gI1,3 + j3))) * tmp_406;
   }
   tmp_404 += tmp_405;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_404;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_410;
      std::complex<double> tmp_411;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_411 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_410 += tmp_411;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI2,0)) * tmp_410;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_412;
      std::complex<double> tmp_413;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_414;
         std::complex<double> tmp_415;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_415 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_414 += tmp_415;
         tmp_413 += (Conj(ZU(gI1,j2))) * tmp_414;
      }
      tmp_412 += tmp_413;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_412;
   }
   std::complex<double> tmp_416;
   std::complex<double> tmp_417;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_418;
      std::complex<double> tmp_419;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_419 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_418 += tmp_419;
      tmp_417 += (Conj(ZU(gI1,j2))) * tmp_418;
   }
   tmp_416 += tmp_417;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI2,1)) * tmp_416;
   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_422;
      std::complex<double> tmp_423;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_424;
         std::complex<double> tmp_425;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_425 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_424 += tmp_425;
         tmp_423 += (KroneckerDelta(gO2,3 + j2)) * tmp_424;
      }
      tmp_422 += tmp_423;
      tmp_421 += (Conj(ZU(gI1,3 + j3))) * tmp_422;
   }
   tmp_420 += tmp_421;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_420;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_426;
      std::complex<double> tmp_427;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_427 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_426 += tmp_427;
      result += (ZP(gI2,1)) * tmp_426;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_428;
      std::complex<double> tmp_429;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_430;
         std::complex<double> tmp_431;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_431 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_430 += tmp_431;
         tmp_429 += (Conj(ZU(gI1,j2))) * tmp_430;
      }
      tmp_428 += tmp_429;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_428;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_432;
   std::complex<double> tmp_433;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_434;
      std::complex<double> tmp_435;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_435 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_434 += tmp_435;
      tmp_433 += (Conj(ZD(gI1,j2))) * tmp_434;
   }
   tmp_432 += tmp_433;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,1))) *
      tmp_432;
   std::complex<double> tmp_436;
   std::complex<double> tmp_437;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_438;
      std::complex<double> tmp_439;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_439 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_438 += tmp_439;
      tmp_437 += (Conj(ZD(gI1,j2))) * tmp_438;
   }
   tmp_436 += tmp_437;
   result += (std::complex<double>(0,-0.5)*vu*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_436;
   std::complex<double> tmp_440;
   std::complex<double> tmp_441;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_442;
      std::complex<double> tmp_443;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_443 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_442 += tmp_443;
      tmp_441 += (Conj(ZD(gI1,j2))) * tmp_442;
   }
   tmp_440 += tmp_441;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_440;
   if (gO2 < 3) {
      std::complex<double> tmp_444;
      std::complex<double> tmp_445;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_445 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_444 += tmp_445;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,1))*Lambdax) *
         tmp_444;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_446;
      std::complex<double> tmp_447;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_447 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_446 += tmp_447;
      result += (std::complex<double>(0,0.5)*vu*Conj(ZA(gI2,2))*Lambdax) *
         tmp_446;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_448;
      std::complex<double> tmp_449;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_449 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_448 += tmp_449;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))
         ) * tmp_448;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_450;
   std::complex<double> tmp_451;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_451 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_450 += tmp_451;
   result += (0.1*vd*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_450;
   std::complex<double> tmp_452;
   std::complex<double> tmp_453;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_453 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_452 += tmp_453;
   result += (-(Qd*QHd*vd*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_452;
   std::complex<double> tmp_454;
   std::complex<double> tmp_455;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_455 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_454 += tmp_455;
   result += (-0.1*vu*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_454;
   std::complex<double> tmp_456;
   std::complex<double> tmp_457;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_457 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_456 += tmp_457;
   result += (-(Qd*QHu*vu*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_456;
   std::complex<double> tmp_458;
   std::complex<double> tmp_459;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_459 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_458 += tmp_459;
   result += (-(Qd*Qs*vS*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_458;
   std::complex<double> tmp_460;
   std::complex<double> tmp_461;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_462;
      std::complex<double> tmp_463;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_463 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_462 += tmp_463;
      tmp_461 += (Conj(ZD(gI1,j2))) * tmp_462;
   }
   tmp_460 += tmp_461;
   result += (0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,1))) * tmp_460;
   std::complex<double> tmp_464;
   std::complex<double> tmp_465;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_466;
      std::complex<double> tmp_467;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_467 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_466 += tmp_467;
      tmp_465 += (Conj(ZD(gI1,j2))) * tmp_466;
   }
   tmp_464 += tmp_465;
   result += (0.5*vu*Conj(Lambdax)*Conj(ZH(gI2,2))) * tmp_464;
   std::complex<double> tmp_468;
   std::complex<double> tmp_469;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_470;
      std::complex<double> tmp_471;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_471 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_470 += tmp_471;
      tmp_469 += (Conj(ZD(gI1,j2))) * tmp_470;
   }
   tmp_468 += tmp_469;
   result += (-0.7071067811865475*Conj(ZH(gI2,0))) * tmp_468;
   std::complex<double> tmp_472;
   std::complex<double> tmp_473;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_474;
      std::complex<double> tmp_475;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_476;
         std::complex<double> tmp_477;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_477 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_476 += tmp_477;
         tmp_475 += (KroneckerDelta(gO2,3 + j2)) * tmp_476;
      }
      tmp_474 += tmp_475;
      tmp_473 += (Conj(ZD(gI1,3 + j3))) * tmp_474;
   }
   tmp_472 += tmp_473;
   result += (-(vd*Conj(ZH(gI2,0)))) * tmp_472;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHd*Qq*vd*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHu*Qq*vu*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += -(Qq*Qs*vS*Conj(ZD(gI1,gO2))*Conj(ZH(gI2,2))*Sqr(gp));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_478;
      std::complex<double> tmp_479;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_479 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_478 += tmp_479;
      result += (0.5*vS*Conj(ZH(gI2,1))*Lambdax) * tmp_478;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_480;
      std::complex<double> tmp_481;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_481 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_480 += tmp_481;
      result += (0.5*vu*Conj(ZH(gI2,2))*Lambdax) * tmp_480;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_482;
      std::complex<double> tmp_483;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_483 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_482 += tmp_483;
      result += (-0.7071067811865475*Conj(ZH(gI2,0))) * tmp_482;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_484;
      std::complex<double> tmp_485;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_486;
         std::complex<double> tmp_487;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_487 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_486 += tmp_487;
         tmp_485 += (Conj(ZD(gI1,j2))) * tmp_486;
      }
      tmp_484 += tmp_485;
      result += (-(vd*Conj(ZH(gI2,0)))) * tmp_484;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_488;
   std::complex<double> tmp_489;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_489 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_488 += tmp_489;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_488;

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

   std::complex<double> tmp_490;
   std::complex<double> tmp_491;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_491 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_490 += tmp_491;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_490;
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

   std::complex<double> tmp_492;
   std::complex<double> tmp_493;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_493 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_492 += tmp_493;
   result += (0.2581988897471611*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_492;
   std::complex<double> tmp_494;
   std::complex<double> tmp_495;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_495 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_494 += tmp_495;
   result += (-(gp*Qd*Sin(ThetaWp()))) * tmp_494;
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

   std::complex<double> tmp_496;
   std::complex<double> tmp_497;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_497 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_496 += tmp_497;
   result += (-(gp*Qd*Cos(ThetaWp()))) * tmp_496;
   std::complex<double> tmp_498;
   std::complex<double> tmp_499;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_499 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_498 += tmp_499;
   result += (-0.2581988897471611*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_498;
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
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_500;
   std::complex<double> tmp_501;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_501 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_500 += tmp_501;
   result += (1.2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_500;
   std::complex<double> tmp_502;
   std::complex<double> tmp_503;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_503 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_502 += tmp_503;
   result += (-3.0983866769659336*g1*gp*Qv*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_502;
   std::complex<double> tmp_504;
   std::complex<double> tmp_505;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_505 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_504 += tmp_505;
   result += (2*Sqr(gp)*Sqr(Qv)*Sqr(Sin(ThetaWp()))) * tmp_504;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(
         Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Cos(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(
         Sin(ThetaW()));
   }
   if (gO1 < 3) {
      result += 2*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
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

std::complex<double> CLASSNAME::CpUSvconjUSvVZpVZp(unsigned gO1, unsigned gO2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_506;
   std::complex<double> tmp_507;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_507 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_506 += tmp_507;
   result += (2*Sqr(gp)*Sqr(Qv)*Sqr(Cos(ThetaWp()))) * tmp_506;
   std::complex<double> tmp_508;
   std::complex<double> tmp_509;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_509 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_508 += tmp_509;
   result += (3.0983866769659336*g1*gp*Qv*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_508;
   std::complex<double> tmp_510;
   std::complex<double> tmp_511;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_511 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_510 += tmp_511;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_510;
   if (gO1 < 3) {
      result += 2*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)*Sqr(Cos(ThetaWp())
         );
   }
   if (gO1 < 3) {
      result += -2*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())*KroneckerDelta(gO1,
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
      result += 0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW())*Sqr(Sin(ThetaWp()));
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(
         Sin(ThetaWp()));
   }

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   std::complex<double> tmp_512;
   std::complex<double> tmp_513;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_513 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_512 += tmp_513;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_512;
   std::complex<double> tmp_514;
   std::complex<double> tmp_515;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_515 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_514 += tmp_515;
   result += (-(QHd*Qv*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0))) * tmp_514;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -(QHd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0)
         );
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_516;
      std::complex<double> tmp_517;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_517 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
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
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_518;
   std::complex<double> tmp_520;
   std::complex<double> tmp_521;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_521 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_520 += tmp_521;
   result += (-(QHu*Qv*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1))) * tmp_520;
   std::complex<double> tmp_522;
   std::complex<double> tmp_523;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_524;
      std::complex<double> tmp_525;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_526;
         std::complex<double> tmp_527;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_527 += Conj(Yv(j1,j3))*Yv(j1,j2);
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
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -(QHu*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1)
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_528;
      std::complex<double> tmp_529;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_529 += Conj(Ye(j1,gO2))*ZER(gI2,j1);
      }
      tmp_528 += tmp_529;
      result += (UM(gI1,1)) * tmp_528;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_530;
   std::complex<double> tmp_531;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_532;
      std::complex<double> tmp_533;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_533 += Conj(ZEL(gI2,j1))*Yv(j1,j2);
      }
      tmp_532 += tmp_533;
      tmp_531 += (KroneckerDelta(gO1,3 + j2)) * tmp_532;
   }
   tmp_530 += tmp_531;
   result += (Conj(UP(gI1,1))) * tmp_530;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZEL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_534;
   std::complex<double> tmp_535;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_536;
      std::complex<double> tmp_537;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_537 += Conj(ZE(gI2,j1))*Yv(j1,j2);
      }
      tmp_536 += tmp_537;
      tmp_535 += (KroneckerDelta(gO2,3 + j2)) * tmp_536;
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
            tmp_543 += Conj(Ye(j3,j1))*Yv(j1,j2);
         }
         tmp_542 += tmp_543;
         tmp_541 += (KroneckerDelta(gO2,3 + j2)) * tmp_542;
      }
      tmp_540 += tmp_541;
      tmp_539 += (Conj(ZE(gI2,3 + j3))) * tmp_540;
   }
   tmp_538 += tmp_539;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_538;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_544;
      std::complex<double> tmp_545;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_545 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,gO2));
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
            tmp_549 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_548 += tmp_549;
         tmp_547 += (Conj(ZE(gI2,j2))) * tmp_548;
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
         tmp_553 += Conj(ZE(gI2,j1))*TYv(j1,j2);
      }
      tmp_552 += tmp_553;
      tmp_551 += (KroneckerDelta(gO2,3 + j2)) * tmp_552;
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
            tmp_559 += Conj(Ye(j3,j1))*Yv(j1,j2);
         }
         tmp_558 += tmp_559;
         tmp_557 += (KroneckerDelta(gO2,3 + j2)) * tmp_558;
      }
      tmp_556 += tmp_557;
      tmp_555 += (Conj(ZE(gI2,3 + j3))) * tmp_556;
   }
   tmp_554 += tmp_555;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_554;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_560;
      std::complex<double> tmp_561;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_561 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
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
            tmp_565 += Conj(Yv(gO2,j1))*Yv(j2,j1);
         }
         tmp_564 += tmp_565;
         tmp_563 += (Conj(ZE(gI2,j2))) * tmp_564;
      }
      tmp_562 += tmp_563;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_562;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_566;
   std::complex<double> tmp_567;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_567 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_566 += tmp_567;
   result += (0.3*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)) * tmp_566;
   std::complex<double> tmp_568;
   std::complex<double> tmp_569;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_569 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_568 += tmp_569;
   result += (-(QHd*Qv*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(gp))) * tmp_568;
   std::complex<double> tmp_570;
   std::complex<double> tmp_571;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_571 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_570 += tmp_571;
   result += (-0.3*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)) * tmp_570;
   std::complex<double> tmp_572;
   std::complex<double> tmp_573;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_573 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_572 += tmp_573;
   result += (-(QHu*Qv*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp))) * tmp_572;
   std::complex<double> tmp_574;
   std::complex<double> tmp_575;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_575 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_574 += tmp_575;
   result += (-(Qs*Qv*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp))) * tmp_574;
   std::complex<double> tmp_576;
   std::complex<double> tmp_577;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_578;
      std::complex<double> tmp_579;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_580;
         std::complex<double> tmp_581;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_581 += Conj(Yv(j1,j3))*Yv(j1,j2);
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
      result += -0.15*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
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
      result += 0.25*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
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
      std::complex<double> tmp_582;
      std::complex<double> tmp_583;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_583 += KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2);
      }
      tmp_582 += tmp_583;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))) *
         tmp_582;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_584;
      std::complex<double> tmp_585;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_585 += KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2);
      }
      tmp_584 += tmp_585;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))) *
         tmp_584;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_586;
      std::complex<double> tmp_587;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_587 += Conj(Yv(gO2,j2))*KroneckerDelta(gO1,3 + j2);
      }
      tmp_586 += tmp_587;
      result += (-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*Lambdax) * tmp_586;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_588;
      std::complex<double> tmp_589;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_589 += Conj(Yv(gO2,j2))*KroneckerDelta(gO1,3 + j2);
      }
      tmp_588 += tmp_589;
      result += (-0.5*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))*Lambdax) * tmp_588;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_590;
      std::complex<double> tmp_591;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_591 += Conj(Yv(gO2,j1))*Yv(gO1,j1);
      }
      tmp_590 += tmp_591;
      result += (-(Conj(ZA(gI1,1))*Conj(ZA(gI2,1)))) * tmp_590;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_592;
   std::complex<double> tmp_593;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_593 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_592 += tmp_593;
   result += (0.3*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_592;
   std::complex<double> tmp_594;
   std::complex<double> tmp_595;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_595 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_594 += tmp_595;
   result += (-(QHd*Qv*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_594;
   std::complex<double> tmp_596;
   std::complex<double> tmp_597;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_597 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_596 += tmp_597;
   result += (-0.3*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_596;
   std::complex<double> tmp_598;
   std::complex<double> tmp_599;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_599 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_598 += tmp_599;
   result += (-(QHu*Qv*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_598;
   std::complex<double> tmp_600;
   std::complex<double> tmp_601;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_601 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_600 += tmp_601;
   result += (-(Qs*Qv*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_600;
   std::complex<double> tmp_602;
   std::complex<double> tmp_603;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_604;
      std::complex<double> tmp_605;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_606;
         std::complex<double> tmp_607;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_607 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_606 += tmp_607;
         tmp_605 += (KroneckerDelta(gO2,3 + j2)) * tmp_606;
      }
      tmp_604 += tmp_605;
      tmp_603 += (KroneckerDelta(gO1,3 + j3)) * tmp_604;
   }
   tmp_602 += tmp_603;
   result += (-(Conj(ZH(gI1,1))*Conj(ZH(gI2,1)))) * tmp_602;
   if (gO1 < 3) {
      result += -0.15*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Ql*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += 0.15*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Ql*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Ql*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_608;
      std::complex<double> tmp_609;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_609 += KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2);
      }
      tmp_608 += tmp_609;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))) *
         tmp_608;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_610;
      std::complex<double> tmp_611;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_611 += KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2);
      }
      tmp_610 += tmp_611;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,0))*Conj(ZH(gI2,2))) *
         tmp_610;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_612;
      std::complex<double> tmp_613;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_613 += Conj(Yv(gO2,j2))*KroneckerDelta(gO1,3 + j2);
      }
      tmp_612 += tmp_613;
      result += (0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))*Lambdax) * tmp_612;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_614;
      std::complex<double> tmp_615;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_615 += Conj(Yv(gO2,j2))*KroneckerDelta(gO1,3 + j2);
      }
      tmp_614 += tmp_615;
      result += (0.5*Conj(ZH(gI1,0))*Conj(ZH(gI2,2))*Lambdax) * tmp_614;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_616;
      std::complex<double> tmp_617;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_617 += Conj(Yv(gO2,j1))*Yv(gO1,j1);
      }
      tmp_616 += tmp_617;
      result += (-(Conj(ZH(gI1,1))*Conj(ZH(gI2,1)))) * tmp_616;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_618;
   std::complex<double> tmp_619;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_619 += KroneckerDelta(gO2,3 + j1)*ZVR(gI1,j1);
   }
   tmp_618 += tmp_619;
   result += (-1.4142135623730951*gp*Qv*ZN(gI2,0)) * tmp_618;
   std::complex<double> tmp_620;
   std::complex<double> tmp_621;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_621 += KroneckerDelta(gO2,3 + j1)*ZVR(gI1,j1);
   }
   tmp_620 += tmp_621;
   result += (-1.0954451150103321*g1*ZN(gI2,1)) * tmp_620;
   if (gO2 < 3) {
      std::complex<double> tmp_622;
      std::complex<double> tmp_623;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_623 += Conj(Yv(gO2,j2))*ZVR(gI1,j2);
      }
      tmp_622 += tmp_623;
      result += (-ZN(gI2,4)) * tmp_622;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_624;
   std::complex<double> tmp_625;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_627 += Conj(ZVL(gI1,j1))*Yv(j1,j2);
      }
      tmp_626 += tmp_627;
      tmp_625 += (KroneckerDelta(gO1,3 + j2)) * tmp_626;
   }
   tmp_624 += tmp_625;
   result += (-Conj(ZN(gI2,4))) * tmp_624;
   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Ql*Conj(ZN(gI2,0))*Conj(ZVL(gI1,gO1))
         ;
   }
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZN(gI2,1))*Conj(ZVL(gI1,gO1));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,2))*Conj(ZVL(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qd = LOCALINPUT(Qd);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_628;
   std::complex<double> tmp_630;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_630 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_628 += tmp_630;
   std::complex<double> tmp_629;
   std::complex<double> tmp_631;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_631 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_629 += tmp_631;
   result += (-0.05*Sqr(g1)) * tmp_628 * tmp_629;
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
   result += (-0.5*Qq*Qv*Sqr(gp)) * tmp_632 * tmp_633;
   std::complex<double> tmp_636;
   std::complex<double> tmp_638;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_638 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_636 += tmp_638;
   std::complex<double> tmp_637;
   std::complex<double> tmp_639;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_639 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_637 += tmp_639;
   result += (-0.1*Sqr(g1)) * tmp_636 * tmp_637;
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
   result += (-0.5*Qd*Qv*Sqr(gp)) * tmp_640 * tmp_641;
   std::complex<double> tmp_644;
   std::complex<double> tmp_646;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_646 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_644 += tmp_646;
   std::complex<double> tmp_645;
   std::complex<double> tmp_647;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_647 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_645 += tmp_647;
   result += (-0.05*Sqr(g1)) * tmp_644 * tmp_645;
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
   result += (-0.5*Qq*Qv*Sqr(gp)) * tmp_648 * tmp_649;
   std::complex<double> tmp_652;
   std::complex<double> tmp_654;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_654 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_652 += tmp_654;
   std::complex<double> tmp_653;
   std::complex<double> tmp_655;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_655 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_653 += tmp_655;
   result += (-0.1*Sqr(g1)) * tmp_652 * tmp_653;
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
   result += (-0.5*Qd*Qv*Sqr(gp)) * tmp_656 * tmp_657;
   if (gO1 < 3) {
      std::complex<double> tmp_660;
      std::complex<double> tmp_661;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_661 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_660 += tmp_661;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_660;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_662;
      std::complex<double> tmp_663;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_663 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_662 += tmp_663;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_662;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_664;
      std::complex<double> tmp_665;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_665 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_664 += tmp_665;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_664;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_666;
      std::complex<double> tmp_667;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_667 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_666 += tmp_667;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_666;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_668;
      std::complex<double> tmp_669;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_669 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_668 += tmp_669;
      result += (-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_668;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_670;
      std::complex<double> tmp_671;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_671 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_670 += tmp_671;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_670;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_672;
      std::complex<double> tmp_673;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_673 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_672 += tmp_673;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_672;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_674;
      std::complex<double> tmp_675;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_675 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_674 += tmp_675;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_674;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_676;
      std::complex<double> tmp_677;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_677 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_676 += tmp_677;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_676;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_678;
      std::complex<double> tmp_679;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_679 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_678 += tmp_679;
      result += (-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_678;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_680;
   std::complex<double> tmp_682;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_682 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_680 += tmp_682;
   std::complex<double> tmp_681;
   std::complex<double> tmp_683;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_683 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_681 += tmp_683;
   result += (0.15*Sqr(g1)) * tmp_680 * tmp_681;
   std::complex<double> tmp_684;
   std::complex<double> tmp_686;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_686 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_684 += tmp_686;
   std::complex<double> tmp_685;
   std::complex<double> tmp_687;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_687 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_685 += tmp_687;
   result += (-0.5*Ql*Qv*Sqr(gp)) * tmp_684 * tmp_685;
   std::complex<double> tmp_688;
   std::complex<double> tmp_690;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_690 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_688 += tmp_690;
   std::complex<double> tmp_689;
   std::complex<double> tmp_691;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_691 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_689 += tmp_691;
   result += (-0.3*Sqr(g1)) * tmp_688 * tmp_689;
   std::complex<double> tmp_692;
   std::complex<double> tmp_694;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_694 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_692 += tmp_694;
   std::complex<double> tmp_693;
   std::complex<double> tmp_695;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_695 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_693 += tmp_695;
   result += (-0.5*Qe*Qv*Sqr(gp)) * tmp_692 * tmp_693;
   std::complex<double> tmp_696;
   std::complex<double> tmp_698;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_698 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_696 += tmp_698;
   std::complex<double> tmp_697;
   std::complex<double> tmp_699;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_699 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_697 += tmp_699;
   result += (0.15*Sqr(g1)) * tmp_696 * tmp_697;
   std::complex<double> tmp_700;
   std::complex<double> tmp_702;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_702 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_700 += tmp_702;
   std::complex<double> tmp_701;
   std::complex<double> tmp_703;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_703 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_701 += tmp_703;
   result += (-0.5*Ql*Qv*Sqr(gp)) * tmp_700 * tmp_701;
   std::complex<double> tmp_704;
   std::complex<double> tmp_706;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_706 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_704 += tmp_706;
   std::complex<double> tmp_705;
   std::complex<double> tmp_707;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_707 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_705 += tmp_707;
   result += (-0.3*Sqr(g1)) * tmp_704 * tmp_705;
   std::complex<double> tmp_708;
   std::complex<double> tmp_710;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_710 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_708 += tmp_710;
   std::complex<double> tmp_709;
   std::complex<double> tmp_711;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_711 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_709 += tmp_711;
   result += (-0.5*Qe*Qv*Sqr(gp)) * tmp_708 * tmp_709;
   std::complex<double> tmp_712;
   std::complex<double> tmp_714;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_715;
      std::complex<double> tmp_716;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_716 += Conj(ZE(gI2,j1))*Yv(j1,j2);
      }
      tmp_715 += tmp_716;
      tmp_714 += (KroneckerDelta(gO2,3 + j2)) * tmp_715;
   }
   tmp_712 += tmp_714;
   std::complex<double> tmp_713;
   std::complex<double> tmp_717;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_718;
      std::complex<double> tmp_719;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_719 += Conj(Yv(j3,j4))*ZE(gI1,j3);
      }
      tmp_718 += tmp_719;
      tmp_717 += (KroneckerDelta(gO1,3 + j4)) * tmp_718;
   }
   tmp_713 += tmp_717;
   result += (-1) * tmp_712 * tmp_713;
   if (gO1 < 3) {
      std::complex<double> tmp_720;
      std::complex<double> tmp_721;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_721 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_720 += tmp_721;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_720;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_722;
      std::complex<double> tmp_723;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_723 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_722 += tmp_723;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_722;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_724;
      std::complex<double> tmp_725;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_725 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_724 += tmp_725;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_724;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_726;
      std::complex<double> tmp_727;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_727 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_726 += tmp_727;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_726;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_728;
      std::complex<double> tmp_729;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_729 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_728 += tmp_729;
      result += (-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_728;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_730;
      std::complex<double> tmp_731;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_731 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_730 += tmp_731;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_730;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_732;
      std::complex<double> tmp_733;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_733 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_732 += tmp_733;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_732;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_734;
      std::complex<double> tmp_735;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_735 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_734 += tmp_735;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_734;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_736;
      std::complex<double> tmp_737;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_737 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_736 += tmp_737;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_736;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_738;
      std::complex<double> tmp_739;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_739 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_738 += tmp_739;
      result += (-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_738;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_740;
      std::complex<double> tmp_742;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_742 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_740 += tmp_742;
      std::complex<double> tmp_741;
      std::complex<double> tmp_743;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_743 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_741 += tmp_743;
      result += (-1) * tmp_740 * tmp_741;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qu = LOCALINPUT(Qu);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_744;
   std::complex<double> tmp_746;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_746 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_744 += tmp_746;
   std::complex<double> tmp_745;
   std::complex<double> tmp_747;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_747 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_745 += tmp_747;
   result += (-0.05*Sqr(g1)) * tmp_744 * tmp_745;
   std::complex<double> tmp_748;
   std::complex<double> tmp_750;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_750 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_748 += tmp_750;
   std::complex<double> tmp_749;
   std::complex<double> tmp_751;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_751 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_749 += tmp_751;
   result += (-0.5*Qq*Qv*Sqr(gp)) * tmp_748 * tmp_749;
   std::complex<double> tmp_752;
   std::complex<double> tmp_754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_754 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_752 += tmp_754;
   std::complex<double> tmp_753;
   std::complex<double> tmp_755;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_755 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_753 += tmp_755;
   result += (0.2*Sqr(g1)) * tmp_752 * tmp_753;
   std::complex<double> tmp_756;
   std::complex<double> tmp_758;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_758 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_756 += tmp_758;
   std::complex<double> tmp_757;
   std::complex<double> tmp_759;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_759 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_757 += tmp_759;
   result += (-0.5*Qu*Qv*Sqr(gp)) * tmp_756 * tmp_757;
   std::complex<double> tmp_760;
   std::complex<double> tmp_762;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_762 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_760 += tmp_762;
   std::complex<double> tmp_761;
   std::complex<double> tmp_763;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_763 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_761 += tmp_763;
   result += (-0.05*Sqr(g1)) * tmp_760 * tmp_761;
   std::complex<double> tmp_764;
   std::complex<double> tmp_766;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_766 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_764 += tmp_766;
   std::complex<double> tmp_765;
   std::complex<double> tmp_767;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_767 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_765 += tmp_767;
   result += (-0.5*Qq*Qv*Sqr(gp)) * tmp_764 * tmp_765;
   std::complex<double> tmp_768;
   std::complex<double> tmp_770;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_770 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_768 += tmp_770;
   std::complex<double> tmp_769;
   std::complex<double> tmp_771;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_771 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_769 += tmp_771;
   result += (0.2*Sqr(g1)) * tmp_768 * tmp_769;
   std::complex<double> tmp_772;
   std::complex<double> tmp_774;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_774 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_772 += tmp_774;
   std::complex<double> tmp_773;
   std::complex<double> tmp_775;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_775 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_773 += tmp_775;
   result += (-0.5*Qu*Qv*Sqr(gp)) * tmp_772 * tmp_773;
   if (gO1 < 3) {
      std::complex<double> tmp_776;
      std::complex<double> tmp_777;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_777 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_776 += tmp_777;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_776;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_778;
      std::complex<double> tmp_779;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_779 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_778 += tmp_779;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_778;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_780;
      std::complex<double> tmp_781;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_781 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_780 += tmp_781;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_780;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_782;
      std::complex<double> tmp_783;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_783 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_782 += tmp_783;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_782;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_784;
      std::complex<double> tmp_785;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_785 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_784 += tmp_785;
      result += (-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_784;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_786;
      std::complex<double> tmp_787;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_787 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_786 += tmp_787;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_786;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_788;
      std::complex<double> tmp_789;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_789 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_788 += tmp_789;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_788;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_790;
      std::complex<double> tmp_791;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_791 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_790 += tmp_791;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_790;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_792;
      std::complex<double> tmp_793;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_793 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_792 += tmp_793;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_792;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_794;
      std::complex<double> tmp_795;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_795 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_794 += tmp_795;
      result += (-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_794;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_796;
      std::complex<double> tmp_798;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_798 += KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2);
      }
      tmp_796 += tmp_798;
      std::complex<double> tmp_797;
      std::complex<double> tmp_799;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_800;
         std::complex<double> tmp_801;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_801 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_800 += tmp_801;
         tmp_799 += (ZU(gI1,j4)) * tmp_800;
      }
      tmp_797 += tmp_799;
      result += (-1) * tmp_796 * tmp_797;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_802;
      std::complex<double> tmp_804;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_805;
         std::complex<double> tmp_806;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_806 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_805 += tmp_806;
         tmp_804 += (Conj(ZU(gI2,j2))) * tmp_805;
      }
      tmp_802 += tmp_804;
      std::complex<double> tmp_803;
      std::complex<double> tmp_807;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         tmp_807 += Conj(Yv(gO2,j4))*KroneckerDelta(gO1,3 + j4);
      }
      tmp_803 += tmp_807;
      result += (-1) * tmp_802 * tmp_803;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_808;
   std::complex<double> tmp_810;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_810 += KroneckerDelta(gO1,3 + j1)*ZV(gI1,3 + j1);
   }
   tmp_808 += tmp_810;
   std::complex<double> tmp_809;
   std::complex<double> tmp_811;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_811 += Conj(ZV(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_809 += tmp_811;
   result += (-0.3*Sqr(g1)) * tmp_808 * tmp_809;
   std::complex<double> tmp_812;
   std::complex<double> tmp_814;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_814 += KroneckerDelta(gO1,3 + j1)*ZV(gI1,3 + j1);
   }
   tmp_812 += tmp_814;
   std::complex<double> tmp_813;
   std::complex<double> tmp_815;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_815 += Conj(ZV(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_813 += tmp_815;
   result += (-0.5*Sqr(gp)*Sqr(Qv)) * tmp_812 * tmp_813;
   std::complex<double> tmp_816;
   std::complex<double> tmp_818;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_818 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_816 += tmp_818;
   std::complex<double> tmp_817;
   std::complex<double> tmp_819;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_819 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_817 += tmp_819;
   result += (0.15*Sqr(g1)) * tmp_816 * tmp_817;
   std::complex<double> tmp_820;
   std::complex<double> tmp_822;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_822 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_820 += tmp_822;
   std::complex<double> tmp_821;
   std::complex<double> tmp_823;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_823 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_821 += tmp_823;
   result += (-0.5*Ql*Qv*Sqr(gp)) * tmp_820 * tmp_821;
   std::complex<double> tmp_824;
   std::complex<double> tmp_826;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_826 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_824 += tmp_826;
   std::complex<double> tmp_825;
   std::complex<double> tmp_827;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_827 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_825 += tmp_827;
   result += (-0.3*Sqr(g1)) * tmp_824 * tmp_825;
   std::complex<double> tmp_828;
   std::complex<double> tmp_830;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_830 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_828 += tmp_830;
   std::complex<double> tmp_829;
   std::complex<double> tmp_831;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_831 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_829 += tmp_831;
   result += (-0.5*Sqr(gp)*Sqr(Qv)) * tmp_828 * tmp_829;
   std::complex<double> tmp_832;
   std::complex<double> tmp_834;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_834 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_832 += tmp_834;
   std::complex<double> tmp_833;
   std::complex<double> tmp_835;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_835 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_833 += tmp_835;
   result += (0.15*Sqr(g1)) * tmp_832 * tmp_833;
   std::complex<double> tmp_836;
   std::complex<double> tmp_838;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_838 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_836 += tmp_838;
   std::complex<double> tmp_837;
   std::complex<double> tmp_839;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_839 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_837 += tmp_839;
   result += (-0.5*Ql*Qv*Sqr(gp)) * tmp_836 * tmp_837;
   std::complex<double> tmp_840;
   std::complex<double> tmp_842;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_842 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_840 += tmp_842;
   std::complex<double> tmp_841;
   std::complex<double> tmp_843;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_843 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_841 += tmp_843;
   result += (-0.3*Sqr(g1)) * tmp_840 * tmp_841;
   std::complex<double> tmp_844;
   std::complex<double> tmp_846;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_846 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_844 += tmp_846;
   std::complex<double> tmp_845;
   std::complex<double> tmp_847;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_847 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_845 += tmp_847;
   result += (-0.5*Sqr(gp)*Sqr(Qv)) * tmp_844 * tmp_845;
   std::complex<double> tmp_848;
   std::complex<double> tmp_850;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_850 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_848 += tmp_850;
   std::complex<double> tmp_849;
   std::complex<double> tmp_851;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_851 += KroneckerDelta(gO1,3 + j2)*ZV(gI1,3 + j2);
   }
   tmp_849 += tmp_851;
   result += (-0.3*Sqr(g1)) * tmp_848 * tmp_849;
   std::complex<double> tmp_852;
   std::complex<double> tmp_854;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_854 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_852 += tmp_854;
   std::complex<double> tmp_853;
   std::complex<double> tmp_855;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_855 += KroneckerDelta(gO1,3 + j2)*ZV(gI1,3 + j2);
   }
   tmp_853 += tmp_855;
   result += (-0.5*Sqr(gp)*Sqr(Qv)) * tmp_852 * tmp_853;
   std::complex<double> tmp_856;
   std::complex<double> tmp_858;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_859;
      std::complex<double> tmp_860;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_860 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_859 += tmp_860;
      tmp_858 += (KroneckerDelta(gO2,3 + j2)) * tmp_859;
   }
   tmp_856 += tmp_858;
   std::complex<double> tmp_857;
   std::complex<double> tmp_861;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_862;
      std::complex<double> tmp_863;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_863 += Conj(Yv(j3,j4))*ZV(gI1,j3);
      }
      tmp_862 += tmp_863;
      tmp_861 += (KroneckerDelta(gO1,3 + j4)) * tmp_862;
   }
   tmp_857 += tmp_861;
   result += (-1) * tmp_856 * tmp_857;
   if (gO1 < 3) {
      std::complex<double> tmp_864;
      std::complex<double> tmp_865;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_865 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_864 += tmp_865;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_864;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_866;
      std::complex<double> tmp_867;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_867 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_866 += tmp_867;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_866;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_868;
      std::complex<double> tmp_869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_869 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_868 += tmp_869;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_868;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_870;
      std::complex<double> tmp_871;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_871 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_870 += tmp_871;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_870;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_872;
      std::complex<double> tmp_873;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_873 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_872 += tmp_873;
      result += (-0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_872;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_874;
      std::complex<double> tmp_875;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_875 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_874 += tmp_875;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_874;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_876;
      std::complex<double> tmp_877;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_877 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_876 += tmp_877;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_876;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_878;
      std::complex<double> tmp_879;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_879 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_878 += tmp_879;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_878;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_880;
      std::complex<double> tmp_881;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_881 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_880 += tmp_881;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_880;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_882;
      std::complex<double> tmp_883;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_883 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_882 += tmp_883;
      result += (-0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_882;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_884;
      std::complex<double> tmp_886;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_886 += KroneckerDelta(gO2,3 + j2)*Yv(gO1,j2);
      }
      tmp_884 += tmp_886;
      std::complex<double> tmp_885;
      std::complex<double> tmp_887;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_888;
         std::complex<double> tmp_889;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_889 += Conj(Yv(j3,j4))*ZV(gI1,j3);
         }
         tmp_888 += tmp_889;
         tmp_887 += (Conj(ZV(gI2,3 + j4))) * tmp_888;
      }
      tmp_885 += tmp_887;
      result += (-1) * tmp_884 * tmp_885;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_890;
      std::complex<double> tmp_891;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_891 += KroneckerDelta(gO1,3 + j1)*ZV(gI1,3 + j1);
      }
      tmp_890 += tmp_891;
      result += (0.15*Conj(ZV(gI2,gO2))*Sqr(g1)) * tmp_890;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_892;
      std::complex<double> tmp_893;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_893 += KroneckerDelta(gO1,3 + j1)*ZV(gI1,3 + j1);
      }
      tmp_892 += tmp_893;
      result += (-0.5*Ql*Qv*Conj(ZV(gI2,gO2))*Sqr(gp)) * tmp_892;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_894;
      std::complex<double> tmp_895;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_895 += KroneckerDelta(gO1,3 + j2)*ZV(gI1,3 + j2);
      }
      tmp_894 += tmp_895;
      result += (0.15*Conj(ZV(gI2,gO2))*Sqr(g1)) * tmp_894;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_896;
      std::complex<double> tmp_897;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_897 += KroneckerDelta(gO1,3 + j2)*ZV(gI1,3 + j2);
      }
      tmp_896 += tmp_897;
      result += (-0.5*Ql*Qv*Conj(ZV(gI2,gO2))*Sqr(gp)) * tmp_896;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_898;
      std::complex<double> tmp_900;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_901;
         std::complex<double> tmp_902;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_902 += Conj(ZV(gI2,j1))*Yv(j1,j2);
         }
         tmp_901 += tmp_902;
         tmp_900 += (ZV(gI1,3 + j2)) * tmp_901;
      }
      tmp_898 += tmp_900;
      std::complex<double> tmp_899;
      std::complex<double> tmp_903;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         tmp_903 += Conj(Yv(gO2,j4))*KroneckerDelta(gO1,3 + j4);
      }
      tmp_899 += tmp_903;
      result += (-1) * tmp_898 * tmp_899;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_904;
      std::complex<double> tmp_906;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_906 += Yv(gO1,j2)*ZV(gI1,3 + j2);
      }
      tmp_904 += tmp_906;
      std::complex<double> tmp_905;
      std::complex<double> tmp_907;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         tmp_907 += Conj(Yv(gO2,j4))*Conj(ZV(gI2,3 + j4));
      }
      tmp_905 += tmp_907;
      result += (-1) * tmp_904 * tmp_905;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_908;
      std::complex<double> tmp_909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_909 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_908 += tmp_909;
      result += (0.15*Sqr(g1)*ZV(gI1,gO1)) * tmp_908;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_910;
      std::complex<double> tmp_911;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_911 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_910 += tmp_911;
      result += (-0.5*Ql*Qv*Sqr(gp)*ZV(gI1,gO1)) * tmp_910;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_912;
      std::complex<double> tmp_913;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_913 += Conj(ZV(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_912 += tmp_913;
      result += (0.15*Sqr(g1)*ZV(gI1,gO1)) * tmp_912;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_914;
      std::complex<double> tmp_915;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_915 += Conj(ZV(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_914 += tmp_915;
      result += (-0.5*Ql*Qv*Sqr(gp)*ZV(gI1,gO1)) * tmp_914;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.15*Conj(ZV(gI2,gO2))*Sqr(g1)*ZV(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -(Conj(ZV(gI2,gO2))*Sqr(gp)*Sqr(Ql)*ZV(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_916;
   std::complex<double> tmp_917;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_918;
      std::complex<double> tmp_919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_919 += Conj(ZV(gI1,j1))*Yv(j1,j2);
      }
      tmp_918 += tmp_919;
      tmp_917 += (KroneckerDelta(gO2,3 + j2)) * tmp_918;
   }
   tmp_916 += tmp_917;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,0))) *
      tmp_916;
   std::complex<double> tmp_920;
   std::complex<double> tmp_921;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_922;
      std::complex<double> tmp_923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_923 += Conj(ZV(gI1,j1))*Yv(j1,j2);
      }
      tmp_922 += tmp_923;
      tmp_921 += (KroneckerDelta(gO2,3 + j2)) * tmp_922;
   }
   tmp_920 += tmp_921;
   result += (std::complex<double>(0,-0.5)*vd*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_920;
   std::complex<double> tmp_924;
   std::complex<double> tmp_925;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_926;
      std::complex<double> tmp_927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_927 += Conj(ZV(gI1,j1))*TYv(j1,j2);
      }
      tmp_926 += tmp_927;
      tmp_925 += (KroneckerDelta(gO2,3 + j2)) * tmp_926;
   }
   tmp_924 += tmp_925;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_924;
   if (gO2 < 3) {
      std::complex<double> tmp_928;
      std::complex<double> tmp_929;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_929 += Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_928 += tmp_929;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,0))*Lambdax) *
         tmp_928;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_930;
      std::complex<double> tmp_931;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_931 += Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_930 += tmp_931;
      result += (std::complex<double>(0,0.5)*vd*Conj(ZA(gI2,2))*Lambdax) *
         tmp_930;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_932;
      std::complex<double> tmp_933;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_933 += Conj(ZV(gI1,3 + j2))*Conj(TYv(gO2,j2));
      }
      tmp_932 += tmp_933;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))
         ) * tmp_932;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qv = LOCALINPUT(Qv);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_934;
   std::complex<double> tmp_935;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_935 += Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_934 += tmp_935;
   result += (0.3*vd*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_934;
   std::complex<double> tmp_936;
   std::complex<double> tmp_937;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_937 += Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_936 += tmp_937;
   result += (-(QHd*Qv*vd*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_936;
   std::complex<double> tmp_938;
   std::complex<double> tmp_939;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_939 += Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_938 += tmp_939;
   result += (-0.3*vu*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_938;
   std::complex<double> tmp_940;
   std::complex<double> tmp_941;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_941 += Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_940 += tmp_941;
   result += (-(QHu*Qv*vu*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_940;
   std::complex<double> tmp_942;
   std::complex<double> tmp_943;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_943 += Conj(ZV(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_942 += tmp_943;
   result += (-(Qs*Qv*vS*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_942;
   std::complex<double> tmp_944;
   std::complex<double> tmp_945;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_946;
      std::complex<double> tmp_947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_947 += Conj(ZV(gI1,j1))*Yv(j1,j2);
      }
      tmp_946 += tmp_947;
      tmp_945 += (KroneckerDelta(gO2,3 + j2)) * tmp_946;
   }
   tmp_944 += tmp_945;
   result += (0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,0))) * tmp_944;
   std::complex<double> tmp_948;
   std::complex<double> tmp_949;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_950;
      std::complex<double> tmp_951;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_951 += Conj(ZV(gI1,j1))*Yv(j1,j2);
      }
      tmp_950 += tmp_951;
      tmp_949 += (KroneckerDelta(gO2,3 + j2)) * tmp_950;
   }
   tmp_948 += tmp_949;
   result += (0.5*vd*Conj(Lambdax)*Conj(ZH(gI2,2))) * tmp_948;
   std::complex<double> tmp_952;
   std::complex<double> tmp_953;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_954;
      std::complex<double> tmp_955;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_955 += Conj(ZV(gI1,j1))*TYv(j1,j2);
      }
      tmp_954 += tmp_955;
      tmp_953 += (KroneckerDelta(gO2,3 + j2)) * tmp_954;
   }
   tmp_952 += tmp_953;
   result += (-0.7071067811865475*Conj(ZH(gI2,1))) * tmp_952;
   std::complex<double> tmp_956;
   std::complex<double> tmp_957;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_958;
      std::complex<double> tmp_959;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_960;
         std::complex<double> tmp_961;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_961 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_960 += tmp_961;
         tmp_959 += (KroneckerDelta(gO2,3 + j2)) * tmp_960;
      }
      tmp_958 += tmp_959;
      tmp_957 += (Conj(ZV(gI1,3 + j3))) * tmp_958;
   }
   tmp_956 += tmp_957;
   result += (-(vu*Conj(ZH(gI2,1)))) * tmp_956;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZH(gI2,0))*Conj(ZV(gI1,gO2))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZH(gI2,0))*Conj(ZV(gI1,gO2))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHd*Ql*vd*Conj(ZH(gI2,0))*Conj(ZV(gI1,gO2))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZH(gI2,1))*Conj(ZV(gI1,gO2))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZH(gI2,1))*Conj(ZV(gI1,gO2))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHu*Ql*vu*Conj(ZH(gI2,1))*Conj(ZV(gI1,gO2))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += -(Ql*Qs*vS*Conj(ZH(gI2,2))*Conj(ZV(gI1,gO2))*Sqr(gp));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_962;
      std::complex<double> tmp_963;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_963 += Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_962 += tmp_963;
      result += (0.5*vS*Conj(ZH(gI2,0))*Lambdax) * tmp_962;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_964;
      std::complex<double> tmp_965;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_965 += Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_964 += tmp_965;
      result += (0.5*vd*Conj(ZH(gI2,2))*Lambdax) * tmp_964;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_966;
      std::complex<double> tmp_967;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_967 += Conj(ZV(gI1,3 + j2))*Conj(TYv(gO2,j2));
      }
      tmp_966 += tmp_967;
      result += (-0.7071067811865475*Conj(ZH(gI2,1))) * tmp_966;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_968;
      std::complex<double> tmp_969;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_970;
         std::complex<double> tmp_971;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_971 += Conj(Yv(gO2,j1))*Yv(j2,j1);
         }
         tmp_970 += tmp_971;
         tmp_969 += (Conj(ZV(gI1,j2))) * tmp_970;
      }
      tmp_968 += tmp_969;
      result += (-(vu*Conj(ZH(gI2,1)))) * tmp_968;
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

std::complex<double> CLASSNAME::CpconjUSvVPSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_972;
   std::complex<double> tmp_973;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_973 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_972 += tmp_973;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_972;
   if (gO2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZV(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvVZSv(unsigned gO2, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_974;
   std::complex<double> tmp_975;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_975 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_974 += tmp_975;
   result += (0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_974;
   std::complex<double> tmp_976;
   std::complex<double> tmp_977;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_977 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_976 += tmp_977;
   result += (-(gp*Qv*Sin(ThetaWp()))) * tmp_976;
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Cos(ThetaWp())*Sin(
         ThetaW());
   }
   if (gO2 < 3) {
      result += gp*Ql*Conj(ZV(gI2,gO2))*Sin(ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvVZpSv(unsigned gO2, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_978;
   std::complex<double> tmp_979;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_979 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_978 += tmp_979;
   result += (-(gp*Qv*Cos(ThetaWp()))) * tmp_978;
   std::complex<double> tmp_980;
   std::complex<double> tmp_981;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_981 += Conj(ZV(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_980 += tmp_981;
   result += (-0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_980;
   if (gO2 < 3) {
      result += gp*Ql*Conj(ZV(gI2,gO2))*Cos(ThetaWp());
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gO2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuVZVZ(unsigned gO1, unsigned gO2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_982;
   std::complex<double> tmp_983;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_983 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_982 += tmp_983;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))
      ) * tmp_982;
   std::complex<double> tmp_984;
   std::complex<double> tmp_985;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_985 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_984 += tmp_985;
   result += (2.065591117977289*g1*gp*Qu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_984;
   std::complex<double> tmp_986;
   std::complex<double> tmp_987;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_987 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_986 += tmp_987;
   result += (2*Sqr(gp)*Sqr(Qu)*Sqr(Sin(ThetaWp()))) * tmp_986;
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

   std::complex<double> tmp_988;
   std::complex<double> tmp_989;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_989 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_988 += tmp_989;
   result += (2*Sqr(gp)*Sqr(Qu)*Sqr(Cos(ThetaWp()))) * tmp_988;
   std::complex<double> tmp_990;
   std::complex<double> tmp_991;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_991 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_990 += tmp_991;
   result += (-2.065591117977289*g1*gp*Qu*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_990;
   std::complex<double> tmp_992;
   std::complex<double> tmp_993;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_993 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_992 += tmp_993;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))
      ) * tmp_992;
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

   std::complex<double> tmp_994;
   std::complex<double> tmp_995;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_995 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_994 += tmp_995;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_994;
   std::complex<double> tmp_996;
   std::complex<double> tmp_997;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_997 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_996 += tmp_997;
   result += (-(QHd*Qu*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0))) * tmp_996;
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
      std::complex<double> tmp_998;
      std::complex<double> tmp_999;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_999 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_998 += tmp_999;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_998;
   }
   std::complex<double> tmp_1000;
   std::complex<double> tmp_1001;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1001 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1000 += tmp_1001;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_1000;
   std::complex<double> tmp_1002;
   std::complex<double> tmp_1003;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1003 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1002 += tmp_1003;
   result += (-(QHu*Qu*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1))) * tmp_1002;
   std::complex<double> tmp_1004;
   std::complex<double> tmp_1005;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1006;
      std::complex<double> tmp_1007;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1008;
         std::complex<double> tmp_1009;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1009 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1008 += tmp_1009;
         tmp_1007 += (KroneckerDelta(gO2,3 + j2)) * tmp_1008;
      }
      tmp_1006 += tmp_1007;
      tmp_1005 += (KroneckerDelta(gO1,3 + j3)) * tmp_1006;
   }
   tmp_1004 += tmp_1005;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_1004;
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
      std::complex<double> tmp_1010;
      std::complex<double> tmp_1011;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1011 += Conj(Yd(j1,gO2))*ZDR(gI2,j1);
      }
      tmp_1010 += tmp_1011;
      result += (UM(gI1,1)) * tmp_1010;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1012;
   std::complex<double> tmp_1013;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1014;
      std::complex<double> tmp_1015;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1015 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_1014 += tmp_1015;
      tmp_1013 += (Conj(ZDL(gI2,j2))) * tmp_1014;
   }
   tmp_1012 += tmp_1013;
   result += (Conj(UP(gI1,1))) * tmp_1012;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1016;
   std::complex<double> tmp_1017;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1018;
      std::complex<double> tmp_1019;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1019 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1018 += tmp_1019;
      tmp_1017 += (Conj(ZD(gI2,j2))) * tmp_1018;
   }
   tmp_1016 += tmp_1017;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI1,0)) * tmp_1016;
   std::complex<double> tmp_1020;
   std::complex<double> tmp_1021;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1022;
      std::complex<double> tmp_1023;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1024;
         std::complex<double> tmp_1025;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1025 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1024 += tmp_1025;
         tmp_1023 += (KroneckerDelta(gO2,3 + j2)) * tmp_1024;
      }
      tmp_1022 += tmp_1023;
      tmp_1021 += (Conj(ZD(gI2,3 + j3))) * tmp_1022;
   }
   tmp_1020 += tmp_1021;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_1020;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1026;
      std::complex<double> tmp_1027;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1027 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_1026 += tmp_1027;
      result += (ZP(gI1,0)) * tmp_1026;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1028;
      std::complex<double> tmp_1029;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1030;
         std::complex<double> tmp_1031;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1031 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_1030 += tmp_1031;
         tmp_1029 += (Conj(ZD(gI2,j2))) * tmp_1030;
      }
      tmp_1028 += tmp_1029;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_1028;
   }
   std::complex<double> tmp_1032;
   std::complex<double> tmp_1033;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1034;
      std::complex<double> tmp_1035;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1035 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_1034 += tmp_1035;
      tmp_1033 += (Conj(ZD(gI2,j2))) * tmp_1034;
   }
   tmp_1032 += tmp_1033;
   result += (ZP(gI1,1)) * tmp_1032;
   std::complex<double> tmp_1036;
   std::complex<double> tmp_1037;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1038;
      std::complex<double> tmp_1039;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1040;
         std::complex<double> tmp_1041;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1041 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1040 += tmp_1041;
         tmp_1039 += (KroneckerDelta(gO2,3 + j2)) * tmp_1040;
      }
      tmp_1038 += tmp_1039;
      tmp_1037 += (Conj(ZD(gI2,3 + j3))) * tmp_1038;
   }
   tmp_1036 += tmp_1037;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_1036;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1042;
      std::complex<double> tmp_1043;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1043 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1042 += tmp_1043;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI1,1)) * tmp_1042;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1044;
      std::complex<double> tmp_1045;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1046;
         std::complex<double> tmp_1047;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1047 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_1046 += tmp_1047;
         tmp_1045 += (Conj(ZD(gI2,j2))) * tmp_1046;
      }
      tmp_1044 += tmp_1045;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_1044;
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

   std::complex<double> tmp_1048;
   std::complex<double> tmp_1049;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1049 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1048 += tmp_1049;
   result += (-0.2*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)) * tmp_1048;
   std::complex<double> tmp_1050;
   std::complex<double> tmp_1051;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1051 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1050 += tmp_1051;
   result += (-(QHd*Qu*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(gp))) * tmp_1050;
   std::complex<double> tmp_1052;
   std::complex<double> tmp_1053;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1053 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1052 += tmp_1053;
   result += (0.2*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)) * tmp_1052;
   std::complex<double> tmp_1054;
   std::complex<double> tmp_1055;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1055 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1054 += tmp_1055;
   result += (-(QHu*Qu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp))) * tmp_1054;
   std::complex<double> tmp_1056;
   std::complex<double> tmp_1057;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1057 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1056 += tmp_1057;
   result += (-(Qs*Qu*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp))) * tmp_1056;
   std::complex<double> tmp_1058;
   std::complex<double> tmp_1059;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1060;
      std::complex<double> tmp_1061;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1062;
         std::complex<double> tmp_1063;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1063 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1062 += tmp_1063;
         tmp_1061 += (KroneckerDelta(gO2,3 + j2)) * tmp_1062;
      }
      tmp_1060 += tmp_1061;
      tmp_1059 += (KroneckerDelta(gO1,3 + j3)) * tmp_1060;
   }
   tmp_1058 += tmp_1059;
   result += (-(Conj(ZA(gI1,1))*Conj(ZA(gI2,1)))) * tmp_1058;
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
      std::complex<double> tmp_1064;
      std::complex<double> tmp_1065;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1065 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_1064 += tmp_1065;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))) *
         tmp_1064;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1066;
      std::complex<double> tmp_1067;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1067 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_1066 += tmp_1067;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))) *
         tmp_1066;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1068;
      std::complex<double> tmp_1069;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1069 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1068 += tmp_1069;
      result += (-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,0))*Lambdax) * tmp_1068;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1070;
      std::complex<double> tmp_1071;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1071 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1070 += tmp_1071;
      result += (-0.5*Conj(ZA(gI1,0))*Conj(ZA(gI2,2))*Lambdax) * tmp_1070;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1072;
      std::complex<double> tmp_1073;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1073 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_1072 += tmp_1073;
      result += (-(Conj(ZA(gI1,1))*Conj(ZA(gI2,1)))) * tmp_1072;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_1074;
   std::complex<double> tmp_1075;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1075 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1074 += tmp_1075;
   result += (-0.2*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_1074;
   std::complex<double> tmp_1076;
   std::complex<double> tmp_1077;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1077 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1076 += tmp_1077;
   result += (-(QHd*Qu*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_1076;
   std::complex<double> tmp_1078;
   std::complex<double> tmp_1079;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1079 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1078 += tmp_1079;
   result += (0.2*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_1078;
   std::complex<double> tmp_1080;
   std::complex<double> tmp_1081;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1081 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1080 += tmp_1081;
   result += (-(QHu*Qu*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_1080;
   std::complex<double> tmp_1082;
   std::complex<double> tmp_1083;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1083 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1082 += tmp_1083;
   result += (-(Qs*Qu*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_1082;
   std::complex<double> tmp_1084;
   std::complex<double> tmp_1085;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1086;
      std::complex<double> tmp_1087;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1088;
         std::complex<double> tmp_1089;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1089 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1088 += tmp_1089;
         tmp_1087 += (KroneckerDelta(gO2,3 + j2)) * tmp_1088;
      }
      tmp_1086 += tmp_1087;
      tmp_1085 += (KroneckerDelta(gO1,3 + j3)) * tmp_1086;
   }
   tmp_1084 += tmp_1085;
   result += (-(Conj(ZH(gI1,1))*Conj(ZH(gI2,1)))) * tmp_1084;
   if (gO1 < 3) {
      result += 0.05*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Qq*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -0.05*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Qq*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Qq*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1090;
      std::complex<double> tmp_1091;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1091 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_1090 += tmp_1091;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))) *
         tmp_1090;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1092;
      std::complex<double> tmp_1093;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1093 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_1092 += tmp_1093;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,0))*Conj(ZH(gI2,2))) *
         tmp_1092;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1094;
      std::complex<double> tmp_1095;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1095 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1094 += tmp_1095;
      result += (0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,0))*Lambdax) * tmp_1094;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1096;
      std::complex<double> tmp_1097;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1097 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1096 += tmp_1097;
      result += (0.5*Conj(ZH(gI1,0))*Conj(ZH(gI2,2))*Lambdax) * tmp_1096;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1098;
      std::complex<double> tmp_1099;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1099 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_1098 += tmp_1099;
      result += (-(Conj(ZH(gI1,1))*Conj(ZH(gI2,1)))) * tmp_1098;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_1100;
   std::complex<double> tmp_1101;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1101 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_1100 += tmp_1101;
   result += (-1.4142135623730951*gp*Qu*ZN(gI2,0)) * tmp_1100;
   std::complex<double> tmp_1102;
   std::complex<double> tmp_1103;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1103 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_1102 += tmp_1103;
   result += (0.7302967433402214*g1*ZN(gI2,1)) * tmp_1102;
   if (gO2 < 3) {
      std::complex<double> tmp_1104;
      std::complex<double> tmp_1105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1105 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_1104 += tmp_1105;
      result += (-ZN(gI2,4)) * tmp_1104;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_1106;
   std::complex<double> tmp_1107;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1108;
      std::complex<double> tmp_1109;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1109 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_1108 += tmp_1109;
      tmp_1107 += (Conj(ZUL(gI1,j2))) * tmp_1108;
   }
   tmp_1106 += tmp_1107;
   result += (-Conj(ZN(gI2,4))) * tmp_1106;
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

   std::complex<double> tmp_1110;
   std::complex<double> tmp_1112;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1112 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1110 += tmp_1112;
   std::complex<double> tmp_1111;
   std::complex<double> tmp_1113;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1113 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1111 += tmp_1113;
   result += (0.1*Sqr(g1)) * tmp_1110 * tmp_1111;
   std::complex<double> tmp_1114;
   std::complex<double> tmp_1116;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1116 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1114 += tmp_1116;
   std::complex<double> tmp_1115;
   std::complex<double> tmp_1117;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1117 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1115 += tmp_1117;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_1114 * tmp_1115;
   std::complex<double> tmp_1118;
   std::complex<double> tmp_1120;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1120 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1118 += tmp_1120;
   std::complex<double> tmp_1119;
   std::complex<double> tmp_1121;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1121 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1119 += tmp_1121;
   result += (0.2*Sqr(g1)) * tmp_1118 * tmp_1119;
   std::complex<double> tmp_1122;
   std::complex<double> tmp_1124;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1124 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1122 += tmp_1124;
   std::complex<double> tmp_1123;
   std::complex<double> tmp_1125;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1125 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1123 += tmp_1125;
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_1122 * tmp_1123;
   std::complex<double> tmp_1126;
   std::complex<double> tmp_1128;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1128 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1126 += tmp_1128;
   std::complex<double> tmp_1127;
   std::complex<double> tmp_1129;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1129 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_1127 += tmp_1129;
   result += (0.1*Sqr(g1)) * tmp_1126 * tmp_1127;
   std::complex<double> tmp_1130;
   std::complex<double> tmp_1132;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1132 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1130 += tmp_1132;
   std::complex<double> tmp_1131;
   std::complex<double> tmp_1133;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1133 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_1131 += tmp_1133;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_1130 * tmp_1131;
   std::complex<double> tmp_1134;
   std::complex<double> tmp_1136;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1136 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1134 += tmp_1136;
   std::complex<double> tmp_1135;
   std::complex<double> tmp_1137;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1137 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_1135 += tmp_1137;
   result += (0.2*Sqr(g1)) * tmp_1134 * tmp_1135;
   std::complex<double> tmp_1138;
   std::complex<double> tmp_1140;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1140 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1138 += tmp_1140;
   std::complex<double> tmp_1139;
   std::complex<double> tmp_1141;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1141 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_1139 += tmp_1141;
   result += (-1.5*Qd*Qu*Sqr(gp)) * tmp_1138 * tmp_1139;
   std::complex<double> tmp_1142;
   std::complex<double> tmp_1144;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1145;
      std::complex<double> tmp_1146;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1146 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1145 += tmp_1146;
      tmp_1144 += (Conj(ZD(gI2,j2))) * tmp_1145;
   }
   tmp_1142 += tmp_1144;
   std::complex<double> tmp_1143;
   std::complex<double> tmp_1147;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1148;
      std::complex<double> tmp_1149;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1149 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1148 += tmp_1149;
      tmp_1147 += (ZD(gI1,j4)) * tmp_1148;
   }
   tmp_1143 += tmp_1147;
   result += (-1) * tmp_1142 * tmp_1143;
   if (gO1 < 3) {
      std::complex<double> tmp_1150;
      std::complex<double> tmp_1151;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1151 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1150 += tmp_1151;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1150;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1152;
      std::complex<double> tmp_1153;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1153 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1152 += tmp_1153;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1152;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1154;
      std::complex<double> tmp_1155;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1155 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1154 += tmp_1155;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_1154;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1156;
      std::complex<double> tmp_1157;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1157 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_1156 += tmp_1157;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1156;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1158;
      std::complex<double> tmp_1159;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1159 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_1158 += tmp_1159;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1158;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1160;
      std::complex<double> tmp_1161;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1161 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1160 += tmp_1161;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1160;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1162;
      std::complex<double> tmp_1163;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1163 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1162 += tmp_1163;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1162;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1164;
      std::complex<double> tmp_1165;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1165 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1164 += tmp_1165;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_1164;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1166;
      std::complex<double> tmp_1167;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1167 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_1166 += tmp_1167;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1166;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1168;
      std::complex<double> tmp_1169;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1169 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_1168 += tmp_1169;
      result += (-1.5*Qd*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1168;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1170;
      std::complex<double> tmp_1172;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1172 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_1170 += tmp_1172;
      std::complex<double> tmp_1171;
      std::complex<double> tmp_1173;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1173 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_1171 += tmp_1173;
      result += (-1) * tmp_1170 * tmp_1171;
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

   std::complex<double> tmp_1174;
   std::complex<double> tmp_1176;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1176 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1174 += tmp_1176;
   std::complex<double> tmp_1175;
   std::complex<double> tmp_1177;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1177 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1175 += tmp_1177;
   result += (-0.1*Sqr(g1)) * tmp_1174 * tmp_1175;
   std::complex<double> tmp_1178;
   std::complex<double> tmp_1180;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1180 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1178 += tmp_1180;
   std::complex<double> tmp_1179;
   std::complex<double> tmp_1181;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1181 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1179 += tmp_1181;
   result += (-0.5*Ql*Qu*Sqr(gp)) * tmp_1178 * tmp_1179;
   std::complex<double> tmp_1182;
   std::complex<double> tmp_1184;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1184 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1182 += tmp_1184;
   std::complex<double> tmp_1183;
   std::complex<double> tmp_1185;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1185 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1183 += tmp_1185;
   result += (0.2*Sqr(g1)) * tmp_1182 * tmp_1183;
   std::complex<double> tmp_1186;
   std::complex<double> tmp_1188;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1188 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1186 += tmp_1188;
   std::complex<double> tmp_1187;
   std::complex<double> tmp_1189;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1189 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1187 += tmp_1189;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_1186 * tmp_1187;
   std::complex<double> tmp_1190;
   std::complex<double> tmp_1192;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1192 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1190 += tmp_1192;
   std::complex<double> tmp_1191;
   std::complex<double> tmp_1193;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1193 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_1191 += tmp_1193;
   result += (-0.1*Sqr(g1)) * tmp_1190 * tmp_1191;
   std::complex<double> tmp_1194;
   std::complex<double> tmp_1196;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1196 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1194 += tmp_1196;
   std::complex<double> tmp_1195;
   std::complex<double> tmp_1197;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1197 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_1195 += tmp_1197;
   result += (-0.5*Ql*Qu*Sqr(gp)) * tmp_1194 * tmp_1195;
   std::complex<double> tmp_1198;
   std::complex<double> tmp_1200;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1200 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1198 += tmp_1200;
   std::complex<double> tmp_1199;
   std::complex<double> tmp_1201;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1201 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_1199 += tmp_1201;
   result += (0.2*Sqr(g1)) * tmp_1198 * tmp_1199;
   std::complex<double> tmp_1202;
   std::complex<double> tmp_1204;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1204 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1202 += tmp_1204;
   std::complex<double> tmp_1203;
   std::complex<double> tmp_1205;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1205 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_1203 += tmp_1205;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_1202 * tmp_1203;
   if (gO1 < 3) {
      std::complex<double> tmp_1206;
      std::complex<double> tmp_1207;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1207 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1206 += tmp_1207;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1206;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1208;
      std::complex<double> tmp_1209;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1209 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1208 += tmp_1209;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1208;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1210;
      std::complex<double> tmp_1211;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1211 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1210 += tmp_1211;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1210;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1212;
      std::complex<double> tmp_1213;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1213 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1212 += tmp_1213;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1212;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1214;
      std::complex<double> tmp_1215;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1215 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1214 += tmp_1215;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1214;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1216;
      std::complex<double> tmp_1217;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1217 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1216 += tmp_1217;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1216;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1218;
      std::complex<double> tmp_1219;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1219 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1218 += tmp_1219;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1218;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1220;
      std::complex<double> tmp_1221;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1221 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1220 += tmp_1221;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1220;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1222;
      std::complex<double> tmp_1223;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1223 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1222 += tmp_1223;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1222;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1224;
      std::complex<double> tmp_1225;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1225 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1224 += tmp_1225;
      result += (-0.5*Qe*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1224;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_1226;
   std::complex<double> tmp_1228;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1228 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_1226 += tmp_1228;
   std::complex<double> tmp_1227;
   std::complex<double> tmp_1229;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1229 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1227 += tmp_1229;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_1226 * tmp_1227;
   std::complex<double> tmp_1230;
   std::complex<double> tmp_1232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1232 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_1230 += tmp_1232;
   std::complex<double> tmp_1231;
   std::complex<double> tmp_1233;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1233 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1231 += tmp_1233;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_1230 * tmp_1231;
   std::complex<double> tmp_1234;
   std::complex<double> tmp_1236;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1236 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_1234 += tmp_1236;
   std::complex<double> tmp_1235;
   std::complex<double> tmp_1237;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1237 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1235 += tmp_1237;
   result += (-0.5*Sqr(gp)*Sqr(Qu)) * tmp_1234 * tmp_1235;
   std::complex<double> tmp_1238;
   std::complex<double> tmp_1240;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1240 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1238 += tmp_1240;
   std::complex<double> tmp_1239;
   std::complex<double> tmp_1241;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1241 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1239 += tmp_1241;
   result += (0.1*Sqr(g1)) * tmp_1238 * tmp_1239;
   std::complex<double> tmp_1242;
   std::complex<double> tmp_1244;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1244 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1242 += tmp_1244;
   std::complex<double> tmp_1243;
   std::complex<double> tmp_1245;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1245 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1243 += tmp_1245;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_1242 * tmp_1243;
   std::complex<double> tmp_1246;
   std::complex<double> tmp_1248;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1248 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1246 += tmp_1248;
   std::complex<double> tmp_1247;
   std::complex<double> tmp_1249;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1249 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1247 += tmp_1249;
   result += (-0.4*Sqr(g1)) * tmp_1246 * tmp_1247;
   std::complex<double> tmp_1250;
   std::complex<double> tmp_1252;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1252 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1250 += tmp_1252;
   std::complex<double> tmp_1251;
   std::complex<double> tmp_1253;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1253 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1251 += tmp_1253;
   result += (-1.5*Sqr(gp)*Sqr(Qu)) * tmp_1250 * tmp_1251;
   std::complex<double> tmp_1254;
   std::complex<double> tmp_1256;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1256 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1254 += tmp_1256;
   std::complex<double> tmp_1255;
   std::complex<double> tmp_1257;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1257 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_1255 += tmp_1257;
   result += (0.1*Sqr(g1)) * tmp_1254 * tmp_1255;
   std::complex<double> tmp_1258;
   std::complex<double> tmp_1260;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1260 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1258 += tmp_1260;
   std::complex<double> tmp_1259;
   std::complex<double> tmp_1261;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1261 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_1259 += tmp_1261;
   result += (-1.5*Qq*Qu*Sqr(gp)) * tmp_1258 * tmp_1259;
   std::complex<double> tmp_1262;
   std::complex<double> tmp_1264;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1264 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1262 += tmp_1264;
   std::complex<double> tmp_1263;
   std::complex<double> tmp_1265;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1265 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_1263 += tmp_1265;
   result += (-0.4*Sqr(g1)) * tmp_1262 * tmp_1263;
   std::complex<double> tmp_1266;
   std::complex<double> tmp_1268;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1268 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1266 += tmp_1268;
   std::complex<double> tmp_1267;
   std::complex<double> tmp_1269;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1269 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_1267 += tmp_1269;
   result += (-1.5*Sqr(gp)*Sqr(Qu)) * tmp_1266 * tmp_1267;
   std::complex<double> tmp_1270;
   std::complex<double> tmp_1272;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1272 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1270 += tmp_1272;
   std::complex<double> tmp_1271;
   std::complex<double> tmp_1273;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1273 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_1271 += tmp_1273;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_1270 * tmp_1271;
   std::complex<double> tmp_1274;
   std::complex<double> tmp_1276;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1276 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1274 += tmp_1276;
   std::complex<double> tmp_1275;
   std::complex<double> tmp_1277;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1277 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_1275 += tmp_1277;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_1274 * tmp_1275;
   std::complex<double> tmp_1278;
   std::complex<double> tmp_1280;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1280 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1278 += tmp_1280;
   std::complex<double> tmp_1279;
   std::complex<double> tmp_1281;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1281 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_1279 += tmp_1281;
   result += (-0.5*Sqr(gp)*Sqr(Qu)) * tmp_1278 * tmp_1279;
   std::complex<double> tmp_1282;
   std::complex<double> tmp_1284;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1285;
      std::complex<double> tmp_1286;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1286 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1285 += tmp_1286;
      tmp_1284 += (Conj(ZU(gI2,j2))) * tmp_1285;
   }
   tmp_1282 += tmp_1284;
   std::complex<double> tmp_1283;
   std::complex<double> tmp_1287;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1288;
      std::complex<double> tmp_1289;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1289 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1288 += tmp_1289;
      tmp_1287 += (ZU(gI1,j4)) * tmp_1288;
   }
   tmp_1283 += tmp_1287;
   result += (-1) * tmp_1282 * tmp_1283;
   if (gO1 < 3) {
      std::complex<double> tmp_1290;
      std::complex<double> tmp_1291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1291 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1290 += tmp_1291;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1290;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1292;
      std::complex<double> tmp_1293;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1293 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1292 += tmp_1293;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1292;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1294;
      std::complex<double> tmp_1295;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1295 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1294 += tmp_1295;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_1294;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1296;
      std::complex<double> tmp_1297;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1297 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1296 += tmp_1297;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1296;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1298;
      std::complex<double> tmp_1299;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1299 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1298 += tmp_1299;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1298;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1300;
      std::complex<double> tmp_1301;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1301 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1300 += tmp_1301;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1300;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1302;
      std::complex<double> tmp_1303;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1303 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1302 += tmp_1303;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1302;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1304;
      std::complex<double> tmp_1305;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1305 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1304 += tmp_1305;
      result += (-1.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Qq)) * tmp_1304;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1306;
      std::complex<double> tmp_1307;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1307 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1306 += tmp_1307;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1306;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1308;
      std::complex<double> tmp_1309;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1309 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1308 += tmp_1309;
      result += (-1.5*Qq*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1308;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1310;
      std::complex<double> tmp_1312;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1312 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_1310 += tmp_1312;
      std::complex<double> tmp_1311;
      std::complex<double> tmp_1313;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1314;
         std::complex<double> tmp_1315;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1315 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_1314 += tmp_1315;
         tmp_1313 += (ZU(gI1,j4)) * tmp_1314;
      }
      tmp_1311 += tmp_1313;
      result += (-3) * tmp_1310 * tmp_1311;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1316;
      std::complex<double> tmp_1317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1317 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_1316 += tmp_1317;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_1316;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1318;
      std::complex<double> tmp_1319;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1319 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_1318 += tmp_1319;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_1318;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1320;
      std::complex<double> tmp_1321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1321 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_1320 += tmp_1321;
      result += (-0.5*Qq*Qu*Conj(ZU(gI2,gO2))*Sqr(gp)) * tmp_1320;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1322;
      std::complex<double> tmp_1323;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1323 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_1322 += tmp_1323;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_1322;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1324;
      std::complex<double> tmp_1325;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1325 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_1324 += tmp_1325;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_1324;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1326;
      std::complex<double> tmp_1327;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1327 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_1326 += tmp_1327;
      result += (-0.5*Qq*Qu*Conj(ZU(gI2,gO2))*Sqr(gp)) * tmp_1326;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1328;
      std::complex<double> tmp_1330;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1331;
         std::complex<double> tmp_1332;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1332 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_1331 += tmp_1332;
         tmp_1330 += (Conj(ZU(gI2,j2))) * tmp_1331;
      }
      tmp_1328 += tmp_1330;
      std::complex<double> tmp_1329;
      std::complex<double> tmp_1333;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1333 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1329 += tmp_1333;
      result += (-3) * tmp_1328 * tmp_1329;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1334;
      std::complex<double> tmp_1336;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1336 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_1334 += tmp_1336;
      std::complex<double> tmp_1335;
      std::complex<double> tmp_1337;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1337 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_1335 += tmp_1337;
      result += (-1) * tmp_1334 * tmp_1335;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1338;
      std::complex<double> tmp_1339;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1339 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1338 += tmp_1339;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_1338;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1340;
      std::complex<double> tmp_1341;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1341 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1340 += tmp_1341;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_1340;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1342;
      std::complex<double> tmp_1343;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1343 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1342 += tmp_1343;
      result += (-0.5*Qq*Qu*Sqr(gp)*ZU(gI1,gO1)) * tmp_1342;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1344;
      std::complex<double> tmp_1345;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1345 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1344 += tmp_1345;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_1344;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1346;
      std::complex<double> tmp_1347;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1347 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1346 += tmp_1347;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_1346;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1348;
      std::complex<double> tmp_1349;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1349 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1348 += tmp_1349;
      result += (-0.5*Qq*Qu*Sqr(gp)*ZU(gI1,gO1)) * tmp_1348;
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

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qu = LOCALINPUT(Qu);
   const auto Qv = LOCALINPUT(Qv);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_1350;
   std::complex<double> tmp_1352;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1352 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1350 += tmp_1352;
   std::complex<double> tmp_1351;
   std::complex<double> tmp_1353;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1353 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1351 += tmp_1353;
   result += (-0.1*Sqr(g1)) * tmp_1350 * tmp_1351;
   std::complex<double> tmp_1354;
   std::complex<double> tmp_1356;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1356 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1354 += tmp_1356;
   std::complex<double> tmp_1355;
   std::complex<double> tmp_1357;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1357 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1355 += tmp_1357;
   result += (-0.5*Ql*Qu*Sqr(gp)) * tmp_1354 * tmp_1355;
   std::complex<double> tmp_1358;
   std::complex<double> tmp_1360;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1360 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_1358 += tmp_1360;
   std::complex<double> tmp_1359;
   std::complex<double> tmp_1361;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1361 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1359 += tmp_1361;
   result += (0.2*Sqr(g1)) * tmp_1358 * tmp_1359;
   std::complex<double> tmp_1362;
   std::complex<double> tmp_1364;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1364 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_1362 += tmp_1364;
   std::complex<double> tmp_1363;
   std::complex<double> tmp_1365;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1365 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1363 += tmp_1365;
   result += (-0.5*Qu*Qv*Sqr(gp)) * tmp_1362 * tmp_1363;
   std::complex<double> tmp_1366;
   std::complex<double> tmp_1368;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1368 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1366 += tmp_1368;
   std::complex<double> tmp_1367;
   std::complex<double> tmp_1369;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1369 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_1367 += tmp_1369;
   result += (-0.1*Sqr(g1)) * tmp_1366 * tmp_1367;
   std::complex<double> tmp_1370;
   std::complex<double> tmp_1372;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1372 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1370 += tmp_1372;
   std::complex<double> tmp_1371;
   std::complex<double> tmp_1373;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1373 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_1371 += tmp_1373;
   result += (-0.5*Ql*Qu*Sqr(gp)) * tmp_1370 * tmp_1371;
   std::complex<double> tmp_1374;
   std::complex<double> tmp_1376;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1376 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1374 += tmp_1376;
   std::complex<double> tmp_1375;
   std::complex<double> tmp_1377;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1377 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_1375 += tmp_1377;
   result += (0.2*Sqr(g1)) * tmp_1374 * tmp_1375;
   std::complex<double> tmp_1378;
   std::complex<double> tmp_1380;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1380 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1378 += tmp_1380;
   std::complex<double> tmp_1379;
   std::complex<double> tmp_1381;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1381 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_1379 += tmp_1381;
   result += (-0.5*Qu*Qv*Sqr(gp)) * tmp_1378 * tmp_1379;
   if (gO1 < 3) {
      std::complex<double> tmp_1382;
      std::complex<double> tmp_1383;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1383 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1382 += tmp_1383;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1382;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1384;
      std::complex<double> tmp_1385;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1385 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1384 += tmp_1385;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1384;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1386;
      std::complex<double> tmp_1387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1387 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1386 += tmp_1387;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1386;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1388;
      std::complex<double> tmp_1389;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1389 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_1388 += tmp_1389;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1388;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1390;
      std::complex<double> tmp_1391;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1391 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_1390 += tmp_1391;
      result += (-0.5*Qq*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1390;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1392;
      std::complex<double> tmp_1393;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1393 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1392 += tmp_1393;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1392;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1394;
      std::complex<double> tmp_1395;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1395 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1394 += tmp_1395;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1394;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1396;
      std::complex<double> tmp_1397;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1397 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1396 += tmp_1397;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1396;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1398;
      std::complex<double> tmp_1399;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1399 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_1398 += tmp_1399;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1398;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1400;
      std::complex<double> tmp_1401;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1401 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_1400 += tmp_1401;
      result += (-0.5*Qq*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1400;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1402;
      std::complex<double> tmp_1404;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1404 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_1402 += tmp_1404;
      std::complex<double> tmp_1403;
      std::complex<double> tmp_1405;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1406;
         std::complex<double> tmp_1407;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1407 += Conj(Yv(j3,j4))*ZV(gI1,j3);
         }
         tmp_1406 += tmp_1407;
         tmp_1405 += (Conj(ZV(gI2,3 + j4))) * tmp_1406;
      }
      tmp_1403 += tmp_1405;
      result += (-1) * tmp_1402 * tmp_1403;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1408;
      std::complex<double> tmp_1410;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1411;
         std::complex<double> tmp_1412;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1412 += Conj(ZV(gI2,j1))*Yv(j1,j2);
         }
         tmp_1411 += tmp_1412;
         tmp_1410 += (ZV(gI1,3 + j2)) * tmp_1411;
      }
      tmp_1408 += tmp_1410;
      std::complex<double> tmp_1409;
      std::complex<double> tmp_1413;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1413 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1409 += tmp_1413;
      result += (-1) * tmp_1408 * tmp_1409;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1414;
   std::complex<double> tmp_1415;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1416;
      std::complex<double> tmp_1417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1417 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1416 += tmp_1417;
      tmp_1415 += (Conj(ZU(gI1,j2))) * tmp_1416;
   }
   tmp_1414 += tmp_1415;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,0))) *
      tmp_1414;
   std::complex<double> tmp_1418;
   std::complex<double> tmp_1419;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1420;
      std::complex<double> tmp_1421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1421 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1420 += tmp_1421;
      tmp_1419 += (Conj(ZU(gI1,j2))) * tmp_1420;
   }
   tmp_1418 += tmp_1419;
   result += (std::complex<double>(0,-0.5)*vd*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_1418;
   std::complex<double> tmp_1422;
   std::complex<double> tmp_1423;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1424;
      std::complex<double> tmp_1425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1425 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_1424 += tmp_1425;
      tmp_1423 += (Conj(ZU(gI1,j2))) * tmp_1424;
   }
   tmp_1422 += tmp_1423;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_1422;
   if (gO2 < 3) {
      std::complex<double> tmp_1426;
      std::complex<double> tmp_1427;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1427 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1426 += tmp_1427;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,0))*Lambdax) *
         tmp_1426;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1428;
      std::complex<double> tmp_1429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1429 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1428 += tmp_1429;
      result += (std::complex<double>(0,0.5)*vd*Conj(ZA(gI2,2))*Lambdax) *
         tmp_1428;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1430;
      std::complex<double> tmp_1431;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1431 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_1430 += tmp_1431;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))
         ) * tmp_1430;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qu = LOCALINPUT(Qu);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_1432;
   std::complex<double> tmp_1433;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1433 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1432 += tmp_1433;
   result += (-0.2*vd*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_1432;
   std::complex<double> tmp_1434;
   std::complex<double> tmp_1435;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1435 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1434 += tmp_1435;
   result += (-(QHd*Qu*vd*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_1434;
   std::complex<double> tmp_1436;
   std::complex<double> tmp_1437;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1437 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1436 += tmp_1437;
   result += (0.2*vu*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_1436;
   std::complex<double> tmp_1438;
   std::complex<double> tmp_1439;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1439 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1438 += tmp_1439;
   result += (-(QHu*Qu*vu*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_1438;
   std::complex<double> tmp_1440;
   std::complex<double> tmp_1441;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1441 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1440 += tmp_1441;
   result += (-(Qs*Qu*vS*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_1440;
   std::complex<double> tmp_1442;
   std::complex<double> tmp_1443;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1444;
      std::complex<double> tmp_1445;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1445 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1444 += tmp_1445;
      tmp_1443 += (Conj(ZU(gI1,j2))) * tmp_1444;
   }
   tmp_1442 += tmp_1443;
   result += (0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,0))) * tmp_1442;
   std::complex<double> tmp_1446;
   std::complex<double> tmp_1447;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1448;
      std::complex<double> tmp_1449;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1449 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_1448 += tmp_1449;
      tmp_1447 += (Conj(ZU(gI1,j2))) * tmp_1448;
   }
   tmp_1446 += tmp_1447;
   result += (0.5*vd*Conj(Lambdax)*Conj(ZH(gI2,2))) * tmp_1446;
   std::complex<double> tmp_1450;
   std::complex<double> tmp_1451;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1452;
      std::complex<double> tmp_1453;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1453 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_1452 += tmp_1453;
      tmp_1451 += (Conj(ZU(gI1,j2))) * tmp_1452;
   }
   tmp_1450 += tmp_1451;
   result += (-0.7071067811865475*Conj(ZH(gI2,1))) * tmp_1450;
   std::complex<double> tmp_1454;
   std::complex<double> tmp_1455;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1456;
      std::complex<double> tmp_1457;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1458;
         std::complex<double> tmp_1459;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1459 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1458 += tmp_1459;
         tmp_1457 += (KroneckerDelta(gO2,3 + j2)) * tmp_1458;
      }
      tmp_1456 += tmp_1457;
      tmp_1455 += (Conj(ZU(gI1,3 + j3))) * tmp_1456;
   }
   tmp_1454 += tmp_1455;
   result += (-(vu*Conj(ZH(gI2,1)))) * tmp_1454;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZH(gI2,0))*Conj(ZU(gI1,gO2))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZH(gI2,0))*Conj(ZU(gI1,gO2))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHd*Qq*vd*Conj(ZH(gI2,0))*Conj(ZU(gI1,gO2))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZH(gI2,1))*Conj(ZU(gI1,gO2))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZH(gI2,1))*Conj(ZU(gI1,gO2))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHu*Qq*vu*Conj(ZH(gI2,1))*Conj(ZU(gI1,gO2))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += -(Qq*Qs*vS*Conj(ZH(gI2,2))*Conj(ZU(gI1,gO2))*Sqr(gp));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1460;
      std::complex<double> tmp_1461;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1461 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1460 += tmp_1461;
      result += (0.5*vS*Conj(ZH(gI2,0))*Lambdax) * tmp_1460;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1462;
      std::complex<double> tmp_1463;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1463 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1462 += tmp_1463;
      result += (0.5*vd*Conj(ZH(gI2,2))*Lambdax) * tmp_1462;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1464;
      std::complex<double> tmp_1465;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1465 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_1464 += tmp_1465;
      result += (-0.7071067811865475*Conj(ZH(gI2,1))) * tmp_1464;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1466;
      std::complex<double> tmp_1467;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1468;
         std::complex<double> tmp_1469;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1469 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_1468 += tmp_1469;
         tmp_1467 += (Conj(ZU(gI1,j2))) * tmp_1468;
      }
      tmp_1466 += tmp_1467;
      result += (-(vu*Conj(ZH(gI2,1)))) * tmp_1466;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1470;
   std::complex<double> tmp_1471;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1471 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1470 += tmp_1471;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1470;

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

   std::complex<double> tmp_1472;
   std::complex<double> tmp_1473;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1473 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1472 += tmp_1473;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_1472;
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

   std::complex<double> tmp_1474;
   std::complex<double> tmp_1475;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1475 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1474 += tmp_1475;
   result += (-0.5163977794943222*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_1474;
   std::complex<double> tmp_1476;
   std::complex<double> tmp_1477;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1477 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1476 += tmp_1477;
   result += (-(gp*Qu*Sin(ThetaWp()))) * tmp_1476;
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

   std::complex<double> tmp_1478;
   std::complex<double> tmp_1479;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1479 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1478 += tmp_1479;
   result += (-(gp*Qu*Cos(ThetaWp()))) * tmp_1478;
   std::complex<double> tmp_1480;
   std::complex<double> tmp_1481;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1481 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1480 += tmp_1481;
   result += (0.5163977794943222*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_1480;
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

   std::complex<double> tmp_1482;
   std::complex<double> tmp_1483;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1483 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1482 += tmp_1483;
   result += (1.2*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_1482;
   std::complex<double> tmp_1484;
   std::complex<double> tmp_1485;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1485 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1484 += tmp_1485;
   result += (-3.0983866769659336*g1*gp*Qe*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_1484;
   std::complex<double> tmp_1486;
   std::complex<double> tmp_1487;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1487 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1486 += tmp_1487;
   result += (2*Sqr(gp)*Sqr(Qe)*Sqr(Sin(ThetaWp()))) * tmp_1486;
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

   std::complex<double> tmp_1488;
   std::complex<double> tmp_1489;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1489 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1488 += tmp_1489;
   result += (2*Sqr(gp)*Sqr(Qe)*Sqr(Cos(ThetaWp()))) * tmp_1488;
   std::complex<double> tmp_1490;
   std::complex<double> tmp_1491;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1491 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1490 += tmp_1491;
   result += (3.0983866769659336*g1*gp*Qe*Cos(ThetaWp())*Sin(ThetaW())*Sin(
      ThetaWp())) * tmp_1490;
   std::complex<double> tmp_1492;
   std::complex<double> tmp_1493;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1493 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1492 += tmp_1493;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_1492;
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

   std::complex<double> tmp_1494;
   std::complex<double> tmp_1495;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1495 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1494 += tmp_1495;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_1494;
   std::complex<double> tmp_1496;
   std::complex<double> tmp_1497;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1497 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1496 += tmp_1497;
   result += (-(Qe*QHd*Sqr(gp)*ZP(gI1,0)*ZP(gI2,0))) * tmp_1496;
   std::complex<double> tmp_1498;
   std::complex<double> tmp_1499;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1500;
      std::complex<double> tmp_1501;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1502;
         std::complex<double> tmp_1503;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1503 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1502 += tmp_1503;
         tmp_1501 += (KroneckerDelta(gO2,3 + j2)) * tmp_1502;
      }
      tmp_1500 += tmp_1501;
      tmp_1499 += (KroneckerDelta(gO1,3 + j3)) * tmp_1500;
   }
   tmp_1498 += tmp_1499;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_1498;
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
   std::complex<double> tmp_1504;
   std::complex<double> tmp_1505;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1505 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1504 += tmp_1505;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_1504;
   std::complex<double> tmp_1506;
   std::complex<double> tmp_1507;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1507 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1506 += tmp_1507;
   result += (-(Qe*QHu*Sqr(gp)*ZP(gI1,1)*ZP(gI2,1))) * tmp_1506;
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
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1508;
      std::complex<double> tmp_1509;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1509 += Conj(Yv(gO2,j1))*Yv(gO1,j1);
      }
      tmp_1508 += tmp_1509;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_1508;
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

   std::complex<double> tmp_1510;
   std::complex<double> tmp_1511;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1511 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1510 += tmp_1511;
   result += (0.3*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(g1)) * tmp_1510;
   std::complex<double> tmp_1512;
   std::complex<double> tmp_1513;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1513 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1512 += tmp_1513;
   result += (-(Qe*QHd*Conj(ZA(gI1,0))*Conj(ZA(gI2,0))*Sqr(gp))) * tmp_1512;
   std::complex<double> tmp_1514;
   std::complex<double> tmp_1515;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1515 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1514 += tmp_1515;
   result += (-0.3*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(g1)) * tmp_1514;
   std::complex<double> tmp_1516;
   std::complex<double> tmp_1517;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1517 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1516 += tmp_1517;
   result += (-(Qe*QHu*Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*Sqr(gp))) * tmp_1516;
   std::complex<double> tmp_1518;
   std::complex<double> tmp_1519;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1519 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1518 += tmp_1519;
   result += (-(Qe*Qs*Conj(ZA(gI1,2))*Conj(ZA(gI2,2))*Sqr(gp))) * tmp_1518;
   std::complex<double> tmp_1520;
   std::complex<double> tmp_1521;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1522;
      std::complex<double> tmp_1523;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1524;
         std::complex<double> tmp_1525;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1525 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1524 += tmp_1525;
         tmp_1523 += (KroneckerDelta(gO2,3 + j2)) * tmp_1524;
      }
      tmp_1522 += tmp_1523;
      tmp_1521 += (KroneckerDelta(gO1,3 + j3)) * tmp_1522;
   }
   tmp_1520 += tmp_1521;
   result += (-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)))) * tmp_1520;
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
      std::complex<double> tmp_1526;
      std::complex<double> tmp_1527;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1527 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1526 += tmp_1527;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))) *
         tmp_1526;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1528;
      std::complex<double> tmp_1529;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1529 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1528 += tmp_1529;
      result += (-0.5*Conj(Lambdax)*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))) *
         tmp_1528;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1530;
      std::complex<double> tmp_1531;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1531 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1530 += tmp_1531;
      result += (-0.5*Conj(ZA(gI1,2))*Conj(ZA(gI2,1))*Lambdax) * tmp_1530;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1532;
      std::complex<double> tmp_1533;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1533 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1532 += tmp_1533;
      result += (-0.5*Conj(ZA(gI1,1))*Conj(ZA(gI2,2))*Lambdax) * tmp_1532;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1534;
      std::complex<double> tmp_1535;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1535 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_1534 += tmp_1535;
      result += (-(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)))) * tmp_1534;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1536;
   std::complex<double> tmp_1537;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1537 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1536 += tmp_1537;
   result += (0.3*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_1536;
   std::complex<double> tmp_1538;
   std::complex<double> tmp_1539;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1539 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1538 += tmp_1539;
   result += (-(Qe*QHd*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_1538;
   std::complex<double> tmp_1540;
   std::complex<double> tmp_1541;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1541 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1540 += tmp_1541;
   result += (-0.3*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_1540;
   std::complex<double> tmp_1542;
   std::complex<double> tmp_1543;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1543 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1542 += tmp_1543;
   result += (-(Qe*QHu*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_1542;
   std::complex<double> tmp_1544;
   std::complex<double> tmp_1545;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1545 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1544 += tmp_1545;
   result += (-(Qe*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_1544;
   std::complex<double> tmp_1546;
   std::complex<double> tmp_1547;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1548;
      std::complex<double> tmp_1549;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1550;
         std::complex<double> tmp_1551;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1551 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1550 += tmp_1551;
         tmp_1549 += (KroneckerDelta(gO2,3 + j2)) * tmp_1550;
      }
      tmp_1548 += tmp_1549;
      tmp_1547 += (KroneckerDelta(gO1,3 + j3)) * tmp_1548;
   }
   tmp_1546 += tmp_1547;
   result += (-(Conj(ZH(gI1,0))*Conj(ZH(gI2,0)))) * tmp_1546;
   if (gO1 < 3) {
      result += -0.15*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2
         )*Sqr(g1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,gO2)
         *Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHd*Ql*Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += 0.15*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2)
         *Sqr(g1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,gO2
         )*Sqr(g2);
   }
   if (gO1 < 3) {
      result += -(QHu*Ql*Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      result += -(Ql*Qs*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*KroneckerDelta(gO1,
         gO2)*Sqr(gp));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1552;
      std::complex<double> tmp_1553;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1553 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1552 += tmp_1553;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))) *
         tmp_1552;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1554;
      std::complex<double> tmp_1555;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1555 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1554 += tmp_1555;
      result += (0.5*Conj(Lambdax)*Conj(ZH(gI1,1))*Conj(ZH(gI2,2))) *
         tmp_1554;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1556;
      std::complex<double> tmp_1557;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1557 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1556 += tmp_1557;
      result += (0.5*Conj(ZH(gI1,2))*Conj(ZH(gI2,1))*Lambdax) * tmp_1556;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1558;
      std::complex<double> tmp_1559;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1559 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_1558 += tmp_1559;
      result += (0.5*Conj(ZH(gI1,1))*Conj(ZH(gI2,2))*Lambdax) * tmp_1558;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1560;
      std::complex<double> tmp_1561;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1561 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_1560 += tmp_1561;
      result += (-(Conj(ZH(gI1,0))*Conj(ZH(gI2,0)))) * tmp_1560;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1562;
      std::complex<double> tmp_1563;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1563 += Conj(Yv(gO2,j2))*ZVR(gI1,j2);
      }
      tmp_1562 += tmp_1563;
      result += (UP(gI2,1)) * tmp_1562;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1564;
   std::complex<double> tmp_1565;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1566;
      std::complex<double> tmp_1567;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1567 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_1566 += tmp_1567;
      tmp_1565 += (Conj(ZVL(gI1,j2))) * tmp_1566;
   }
   tmp_1564 += tmp_1565;
   result += (Conj(UM(gI2,1))) * tmp_1564;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(ZVL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_1568;
   std::complex<double> tmp_1569;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1569 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_1568 += tmp_1569;
   result += (-1.4142135623730951*gp*Qe*ZN(gI2,0)) * tmp_1568;
   std::complex<double> tmp_1570;
   std::complex<double> tmp_1571;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1571 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_1570 += tmp_1571;
   result += (-1.0954451150103321*g1*ZN(gI2,1)) * tmp_1570;
   if (gO2 < 3) {
      std::complex<double> tmp_1572;
      std::complex<double> tmp_1573;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1573 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_1572 += tmp_1573;
      result += (-ZN(gI2,3)) * tmp_1572;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1574;
   std::complex<double> tmp_1575;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1576;
      std::complex<double> tmp_1577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1577 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_1576 += tmp_1577;
      tmp_1575 += (Conj(ZEL(gI1,j2))) * tmp_1576;
   }
   tmp_1574 += tmp_1575;
   result += (-Conj(ZN(gI2,3))) * tmp_1574;
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

   std::complex<double> tmp_1578;
   std::complex<double> tmp_1580;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1580 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1578 += tmp_1580;
   std::complex<double> tmp_1579;
   std::complex<double> tmp_1581;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1581 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1579 += tmp_1581;
   result += (-0.05*Sqr(g1)) * tmp_1578 * tmp_1579;
   std::complex<double> tmp_1582;
   std::complex<double> tmp_1584;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1584 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1582 += tmp_1584;
   std::complex<double> tmp_1583;
   std::complex<double> tmp_1585;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1585 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1583 += tmp_1585;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1582 * tmp_1583;
   std::complex<double> tmp_1586;
   std::complex<double> tmp_1588;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1588 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1586 += tmp_1588;
   std::complex<double> tmp_1587;
   std::complex<double> tmp_1589;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1589 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1587 += tmp_1589;
   result += (-0.1*Sqr(g1)) * tmp_1586 * tmp_1587;
   std::complex<double> tmp_1590;
   std::complex<double> tmp_1592;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1592 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1590 += tmp_1592;
   std::complex<double> tmp_1591;
   std::complex<double> tmp_1593;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1593 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1591 += tmp_1593;
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_1590 * tmp_1591;
   std::complex<double> tmp_1594;
   std::complex<double> tmp_1596;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1596 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1594 += tmp_1596;
   std::complex<double> tmp_1595;
   std::complex<double> tmp_1597;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1597 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_1595 += tmp_1597;
   result += (-0.05*Sqr(g1)) * tmp_1594 * tmp_1595;
   std::complex<double> tmp_1598;
   std::complex<double> tmp_1600;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1600 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1598 += tmp_1600;
   std::complex<double> tmp_1599;
   std::complex<double> tmp_1601;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1601 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_1599 += tmp_1601;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1598 * tmp_1599;
   std::complex<double> tmp_1602;
   std::complex<double> tmp_1604;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1604 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1602 += tmp_1604;
   std::complex<double> tmp_1603;
   std::complex<double> tmp_1605;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1605 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_1603 += tmp_1605;
   result += (-0.1*Sqr(g1)) * tmp_1602 * tmp_1603;
   std::complex<double> tmp_1606;
   std::complex<double> tmp_1608;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1608 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1606 += tmp_1608;
   std::complex<double> tmp_1607;
   std::complex<double> tmp_1609;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1609 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_1607 += tmp_1609;
   result += (-0.5*Qd*Qe*Sqr(gp)) * tmp_1606 * tmp_1607;
   if (gO1 < 3) {
      std::complex<double> tmp_1610;
      std::complex<double> tmp_1611;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1611 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1610 += tmp_1611;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1610;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1612;
      std::complex<double> tmp_1613;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1613 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1612 += tmp_1613;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1612;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1614;
      std::complex<double> tmp_1615;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1615 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_1614 += tmp_1615;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1614;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1616;
      std::complex<double> tmp_1617;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1617 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_1616 += tmp_1617;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1616;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1618;
      std::complex<double> tmp_1619;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1619 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_1618 += tmp_1619;
      result += (-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1618;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1620;
      std::complex<double> tmp_1621;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1621 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1620 += tmp_1621;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1620;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1622;
      std::complex<double> tmp_1623;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1623 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1622 += tmp_1623;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1622;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1624;
      std::complex<double> tmp_1625;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1625 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_1624 += tmp_1625;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1624;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1626;
      std::complex<double> tmp_1627;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1627 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_1626 += tmp_1627;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1626;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1628;
      std::complex<double> tmp_1629;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1629 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_1628 += tmp_1629;
      result += (-0.5*Qd*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1628;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1630;
      std::complex<double> tmp_1632;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1632 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1630 += tmp_1632;
      std::complex<double> tmp_1631;
      std::complex<double> tmp_1633;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1634;
         std::complex<double> tmp_1635;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1635 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_1634 += tmp_1635;
         tmp_1633 += (ZD(gI1,j4)) * tmp_1634;
      }
      tmp_1631 += tmp_1633;
      result += (-1) * tmp_1630 * tmp_1631;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1636;
      std::complex<double> tmp_1638;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1639;
         std::complex<double> tmp_1640;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1640 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_1639 += tmp_1640;
         tmp_1638 += (Conj(ZD(gI2,j2))) * tmp_1639;
      }
      tmp_1636 += tmp_1638;
      std::complex<double> tmp_1637;
      std::complex<double> tmp_1641;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1641 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1637 += tmp_1641;
      result += (-1) * tmp_1636 * tmp_1637;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1642;
   std::complex<double> tmp_1644;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1644 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_1642 += tmp_1644;
   std::complex<double> tmp_1643;
   std::complex<double> tmp_1645;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1645 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1643 += tmp_1645;
   result += (-0.3*Sqr(g1)) * tmp_1642 * tmp_1643;
   std::complex<double> tmp_1646;
   std::complex<double> tmp_1648;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1648 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_1646 += tmp_1648;
   std::complex<double> tmp_1647;
   std::complex<double> tmp_1649;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1649 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1647 += tmp_1649;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1646 * tmp_1647;
   std::complex<double> tmp_1650;
   std::complex<double> tmp_1652;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1652 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1650 += tmp_1652;
   std::complex<double> tmp_1651;
   std::complex<double> tmp_1653;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1653 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1651 += tmp_1653;
   result += (0.15*Sqr(g1)) * tmp_1650 * tmp_1651;
   std::complex<double> tmp_1654;
   std::complex<double> tmp_1656;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1656 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1654 += tmp_1656;
   std::complex<double> tmp_1655;
   std::complex<double> tmp_1657;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1657 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1655 += tmp_1657;
   result += (-0.5*Qe*Ql*Sqr(gp)) * tmp_1654 * tmp_1655;
   std::complex<double> tmp_1658;
   std::complex<double> tmp_1660;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1660 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1658 += tmp_1660;
   std::complex<double> tmp_1659;
   std::complex<double> tmp_1661;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1661 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1659 += tmp_1661;
   result += (-0.3*Sqr(g1)) * tmp_1658 * tmp_1659;
   std::complex<double> tmp_1662;
   std::complex<double> tmp_1664;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1664 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1662 += tmp_1664;
   std::complex<double> tmp_1663;
   std::complex<double> tmp_1665;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1665 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1663 += tmp_1665;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1662 * tmp_1663;
   std::complex<double> tmp_1666;
   std::complex<double> tmp_1668;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1668 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1666 += tmp_1668;
   std::complex<double> tmp_1667;
   std::complex<double> tmp_1669;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1669 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_1667 += tmp_1669;
   result += (0.15*Sqr(g1)) * tmp_1666 * tmp_1667;
   std::complex<double> tmp_1670;
   std::complex<double> tmp_1672;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1672 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1670 += tmp_1672;
   std::complex<double> tmp_1671;
   std::complex<double> tmp_1673;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1673 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_1671 += tmp_1673;
   result += (-0.5*Qe*Ql*Sqr(gp)) * tmp_1670 * tmp_1671;
   std::complex<double> tmp_1674;
   std::complex<double> tmp_1676;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1676 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1674 += tmp_1676;
   std::complex<double> tmp_1675;
   std::complex<double> tmp_1677;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1677 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_1675 += tmp_1677;
   result += (-0.3*Sqr(g1)) * tmp_1674 * tmp_1675;
   std::complex<double> tmp_1678;
   std::complex<double> tmp_1680;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1680 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1678 += tmp_1680;
   std::complex<double> tmp_1679;
   std::complex<double> tmp_1681;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1681 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_1679 += tmp_1681;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1678 * tmp_1679;
   std::complex<double> tmp_1682;
   std::complex<double> tmp_1684;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1684 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1682 += tmp_1684;
   std::complex<double> tmp_1683;
   std::complex<double> tmp_1685;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1685 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_1683 += tmp_1685;
   result += (-0.3*Sqr(g1)) * tmp_1682 * tmp_1683;
   std::complex<double> tmp_1686;
   std::complex<double> tmp_1688;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1688 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1686 += tmp_1688;
   std::complex<double> tmp_1687;
   std::complex<double> tmp_1689;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1689 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_1687 += tmp_1689;
   result += (-0.5*Sqr(gp)*Sqr(Qe)) * tmp_1686 * tmp_1687;
   std::complex<double> tmp_1690;
   std::complex<double> tmp_1692;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1693;
      std::complex<double> tmp_1694;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1694 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1693 += tmp_1694;
      tmp_1692 += (Conj(ZE(gI2,j2))) * tmp_1693;
   }
   tmp_1690 += tmp_1692;
   std::complex<double> tmp_1691;
   std::complex<double> tmp_1695;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1696;
      std::complex<double> tmp_1697;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1697 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1696 += tmp_1697;
      tmp_1695 += (ZE(gI1,j4)) * tmp_1696;
   }
   tmp_1691 += tmp_1695;
   result += (-1) * tmp_1690 * tmp_1691;
   if (gO1 < 3) {
      std::complex<double> tmp_1698;
      std::complex<double> tmp_1699;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1699 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1698 += tmp_1699;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1698;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1700;
      std::complex<double> tmp_1701;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1701 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1700 += tmp_1701;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1700;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1702;
      std::complex<double> tmp_1703;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1703 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_1702 += tmp_1703;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_1702;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1704;
      std::complex<double> tmp_1705;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1705 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1704 += tmp_1705;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1704;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1706;
      std::complex<double> tmp_1707;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1707 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_1706 += tmp_1707;
      result += (-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1706;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1708;
      std::complex<double> tmp_1709;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1709 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1708 += tmp_1709;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1708;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1710;
      std::complex<double> tmp_1711;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1711 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1710 += tmp_1711;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1710;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1712;
      std::complex<double> tmp_1713;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1713 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_1712 += tmp_1713;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_1712;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1714;
      std::complex<double> tmp_1715;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1715 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1714 += tmp_1715;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1714;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1716;
      std::complex<double> tmp_1717;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1717 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_1716 += tmp_1717;
      result += (-0.5*Qe*Ql*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1716;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1718;
      std::complex<double> tmp_1720;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1720 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_1718 += tmp_1720;
      std::complex<double> tmp_1719;
      std::complex<double> tmp_1721;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_1722;
         std::complex<double> tmp_1723;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_1723 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_1722 += tmp_1723;
         tmp_1721 += (ZE(gI1,j4)) * tmp_1722;
      }
      tmp_1719 += tmp_1721;
      result += (-1) * tmp_1718 * tmp_1719;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1724;
      std::complex<double> tmp_1725;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1725 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_1724 += tmp_1725;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_1724;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1726;
      std::complex<double> tmp_1727;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1727 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_1726 += tmp_1727;
      result += (-0.5*Qe*Ql*Conj(ZE(gI2,gO2))*Sqr(gp)) * tmp_1726;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1728;
      std::complex<double> tmp_1729;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1729 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_1728 += tmp_1729;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_1728;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1730;
      std::complex<double> tmp_1731;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1731 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_1730 += tmp_1731;
      result += (-0.5*Qe*Ql*Conj(ZE(gI2,gO2))*Sqr(gp)) * tmp_1730;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1732;
      std::complex<double> tmp_1734;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1735;
         std::complex<double> tmp_1736;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1736 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_1735 += tmp_1736;
         tmp_1734 += (Conj(ZE(gI2,j2))) * tmp_1735;
      }
      tmp_1732 += tmp_1734;
      std::complex<double> tmp_1733;
      std::complex<double> tmp_1737;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1737 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1733 += tmp_1737;
      result += (-1) * tmp_1732 * tmp_1733;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1738;
      std::complex<double> tmp_1740;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1740 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_1738 += tmp_1740;
      std::complex<double> tmp_1739;
      std::complex<double> tmp_1741;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1741 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_1739 += tmp_1741;
      result += (-1) * tmp_1738 * tmp_1739;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1742;
      std::complex<double> tmp_1743;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1743 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1742 += tmp_1743;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_1742;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1744;
      std::complex<double> tmp_1745;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1745 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_1744 += tmp_1745;
      result += (-0.5*Qe*Ql*Sqr(gp)*ZE(gI1,gO1)) * tmp_1744;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1746;
      std::complex<double> tmp_1747;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1747 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1746 += tmp_1747;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_1746;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1748;
      std::complex<double> tmp_1749;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1749 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_1748 += tmp_1749;
      result += (-0.5*Qe*Ql*Sqr(gp)*ZE(gI1,gO1)) * tmp_1748;
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

   std::complex<double> tmp_1750;
   std::complex<double> tmp_1752;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1752 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1750 += tmp_1752;
   std::complex<double> tmp_1751;
   std::complex<double> tmp_1753;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1753 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1751 += tmp_1753;
   result += (-0.05*Sqr(g1)) * tmp_1750 * tmp_1751;
   std::complex<double> tmp_1754;
   std::complex<double> tmp_1756;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1756 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1754 += tmp_1756;
   std::complex<double> tmp_1755;
   std::complex<double> tmp_1757;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1757 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1755 += tmp_1757;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1754 * tmp_1755;
   std::complex<double> tmp_1758;
   std::complex<double> tmp_1760;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1760 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1758 += tmp_1760;
   std::complex<double> tmp_1759;
   std::complex<double> tmp_1761;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1761 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1759 += tmp_1761;
   result += (0.2*Sqr(g1)) * tmp_1758 * tmp_1759;
   std::complex<double> tmp_1762;
   std::complex<double> tmp_1764;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1764 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1762 += tmp_1764;
   std::complex<double> tmp_1763;
   std::complex<double> tmp_1765;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1765 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1763 += tmp_1765;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_1762 * tmp_1763;
   std::complex<double> tmp_1766;
   std::complex<double> tmp_1768;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1768 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1766 += tmp_1768;
   std::complex<double> tmp_1767;
   std::complex<double> tmp_1769;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1769 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_1767 += tmp_1769;
   result += (-0.05*Sqr(g1)) * tmp_1766 * tmp_1767;
   std::complex<double> tmp_1770;
   std::complex<double> tmp_1772;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1772 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1770 += tmp_1772;
   std::complex<double> tmp_1771;
   std::complex<double> tmp_1773;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1773 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_1771 += tmp_1773;
   result += (-0.5*Qe*Qq*Sqr(gp)) * tmp_1770 * tmp_1771;
   std::complex<double> tmp_1774;
   std::complex<double> tmp_1776;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1776 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1774 += tmp_1776;
   std::complex<double> tmp_1775;
   std::complex<double> tmp_1777;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1777 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_1775 += tmp_1777;
   result += (0.2*Sqr(g1)) * tmp_1774 * tmp_1775;
   std::complex<double> tmp_1778;
   std::complex<double> tmp_1780;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1780 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1778 += tmp_1780;
   std::complex<double> tmp_1779;
   std::complex<double> tmp_1781;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1781 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_1779 += tmp_1781;
   result += (-0.5*Qe*Qu*Sqr(gp)) * tmp_1778 * tmp_1779;
   if (gO1 < 3) {
      std::complex<double> tmp_1782;
      std::complex<double> tmp_1783;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1783 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1782 += tmp_1783;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1782;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1784;
      std::complex<double> tmp_1785;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1785 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1784 += tmp_1785;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1784;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1786;
      std::complex<double> tmp_1787;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1787 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_1786 += tmp_1787;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1786;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1788;
      std::complex<double> tmp_1789;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1789 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1788 += tmp_1789;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1788;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1790;
      std::complex<double> tmp_1791;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1791 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_1790 += tmp_1791;
      result += (-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1790;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1792;
      std::complex<double> tmp_1793;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1793 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1792 += tmp_1793;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1792;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1794;
      std::complex<double> tmp_1795;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1795 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1794 += tmp_1795;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1794;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1796;
      std::complex<double> tmp_1797;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1797 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_1796 += tmp_1797;
      result += (-0.5*Ql*Qq*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1796;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1798;
      std::complex<double> tmp_1799;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1799 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1798 += tmp_1799;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1798;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1800;
      std::complex<double> tmp_1801;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1801 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_1800 += tmp_1801;
      result += (-0.5*Ql*Qu*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1800;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_1802;
   std::complex<double> tmp_1804;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1804 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1802 += tmp_1804;
   std::complex<double> tmp_1803;
   std::complex<double> tmp_1805;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1805 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1803 += tmp_1805;
   result += (0.15*Sqr(g1)) * tmp_1802 * tmp_1803;
   std::complex<double> tmp_1806;
   std::complex<double> tmp_1808;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1808 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1806 += tmp_1808;
   std::complex<double> tmp_1807;
   std::complex<double> tmp_1809;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1809 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1807 += tmp_1809;
   result += (-0.5*Qe*Ql*Sqr(gp)) * tmp_1806 * tmp_1807;
   std::complex<double> tmp_1810;
   std::complex<double> tmp_1812;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1812 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_1810 += tmp_1812;
   std::complex<double> tmp_1811;
   std::complex<double> tmp_1813;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1813 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1811 += tmp_1813;
   result += (-0.3*Sqr(g1)) * tmp_1810 * tmp_1811;
   std::complex<double> tmp_1814;
   std::complex<double> tmp_1816;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1816 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_1814 += tmp_1816;
   std::complex<double> tmp_1815;
   std::complex<double> tmp_1817;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1817 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_1815 += tmp_1817;
   result += (-0.5*Qe*Qv*Sqr(gp)) * tmp_1814 * tmp_1815;
   std::complex<double> tmp_1818;
   std::complex<double> tmp_1820;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1820 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1818 += tmp_1820;
   std::complex<double> tmp_1819;
   std::complex<double> tmp_1821;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1821 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_1819 += tmp_1821;
   result += (0.15*Sqr(g1)) * tmp_1818 * tmp_1819;
   std::complex<double> tmp_1822;
   std::complex<double> tmp_1824;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1824 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1822 += tmp_1824;
   std::complex<double> tmp_1823;
   std::complex<double> tmp_1825;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1825 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
   }
   tmp_1823 += tmp_1825;
   result += (-0.5*Qe*Ql*Sqr(gp)) * tmp_1822 * tmp_1823;
   std::complex<double> tmp_1826;
   std::complex<double> tmp_1828;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1828 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1826 += tmp_1828;
   std::complex<double> tmp_1827;
   std::complex<double> tmp_1829;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1829 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_1827 += tmp_1829;
   result += (-0.3*Sqr(g1)) * tmp_1826 * tmp_1827;
   std::complex<double> tmp_1830;
   std::complex<double> tmp_1832;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1832 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1830 += tmp_1832;
   std::complex<double> tmp_1831;
   std::complex<double> tmp_1833;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_1833 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
   }
   tmp_1831 += tmp_1833;
   result += (-0.5*Qe*Qv*Sqr(gp)) * tmp_1830 * tmp_1831;
   std::complex<double> tmp_1834;
   std::complex<double> tmp_1836;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1837;
      std::complex<double> tmp_1838;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1838 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1837 += tmp_1838;
      tmp_1836 += (Conj(ZV(gI2,j2))) * tmp_1837;
   }
   tmp_1834 += tmp_1836;
   std::complex<double> tmp_1835;
   std::complex<double> tmp_1839;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_1840;
      std::complex<double> tmp_1841;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_1841 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_1840 += tmp_1841;
      tmp_1839 += (ZV(gI1,j4)) * tmp_1840;
   }
   tmp_1835 += tmp_1839;
   result += (-1) * tmp_1834 * tmp_1835;
   if (gO1 < 3) {
      std::complex<double> tmp_1842;
      std::complex<double> tmp_1843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1843 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1842 += tmp_1843;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1842;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1844;
      std::complex<double> tmp_1845;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1845 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1844 += tmp_1845;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1844;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1846;
      std::complex<double> tmp_1847;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1847 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
      }
      tmp_1846 += tmp_1847;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_1846;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1848;
      std::complex<double> tmp_1849;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1849 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_1848 += tmp_1849;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1848;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1850;
      std::complex<double> tmp_1851;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1851 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
      }
      tmp_1850 += tmp_1851;
      result += (-0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1850;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1852;
      std::complex<double> tmp_1853;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1853 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1852 += tmp_1853;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1852;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1854;
      std::complex<double> tmp_1855;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1855 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1854 += tmp_1855;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_1854;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1856;
      std::complex<double> tmp_1857;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1857 += Conj(ZV(gI2,j2))*ZV(gI1,j2);
      }
      tmp_1856 += tmp_1857;
      result += (-0.5*KroneckerDelta(gO1,gO2)*Sqr(gp)*Sqr(Ql)) * tmp_1856;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1858;
      std::complex<double> tmp_1859;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1859 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_1858 += tmp_1859;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_1858;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1860;
      std::complex<double> tmp_1861;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1861 += Conj(ZV(gI2,3 + j2))*ZV(gI1,3 + j2);
      }
      tmp_1860 += tmp_1861;
      result += (-0.5*Ql*Qv*KroneckerDelta(gO1,gO2)*Sqr(gp)) * tmp_1860;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_1862;
      std::complex<double> tmp_1864;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1864 += Yv(gO1,j2)*ZV(gI1,3 + j2);
      }
      tmp_1862 += tmp_1864;
      std::complex<double> tmp_1863;
      std::complex<double> tmp_1865;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         tmp_1865 += Conj(Yv(gO2,j4))*Conj(ZV(gI2,3 + j4));
      }
      tmp_1863 += tmp_1865;
      result += (-1) * tmp_1862 * tmp_1863;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1866;
   std::complex<double> tmp_1867;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1868;
      std::complex<double> tmp_1869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1869 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_1868 += tmp_1869;
      tmp_1867 += (Conj(ZV(gI1,j2))) * tmp_1868;
   }
   tmp_1866 += tmp_1867;
   result += (ZP(gI2,0)) * tmp_1866;
   std::complex<double> tmp_1870;
   std::complex<double> tmp_1871;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1872;
      std::complex<double> tmp_1873;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1874;
         std::complex<double> tmp_1875;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1875 += Conj(Yv(j1,j3))*Ye(j2,j1);
         }
         tmp_1874 += tmp_1875;
         tmp_1873 += (KroneckerDelta(gO2,3 + j2)) * tmp_1874;
      }
      tmp_1872 += tmp_1873;
      tmp_1871 += (Conj(ZV(gI1,3 + j3))) * tmp_1872;
   }
   tmp_1870 += tmp_1871;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_1870;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1876;
      std::complex<double> tmp_1877;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1877 += Conj(Yv(gO2,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_1876 += tmp_1877;
      result += (0.7071067811865475*vS*Lambdax*ZP(gI2,0)) * tmp_1876;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1878;
      std::complex<double> tmp_1879;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1880;
         std::complex<double> tmp_1881;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1881 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_1880 += tmp_1881;
         tmp_1879 += (Conj(ZV(gI1,j2))) * tmp_1880;
      }
      tmp_1878 += tmp_1879;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_1878;
   }
   std::complex<double> tmp_1882;
   std::complex<double> tmp_1883;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1884;
      std::complex<double> tmp_1885;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1885 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1884 += tmp_1885;
      tmp_1883 += (Conj(ZV(gI1,j2))) * tmp_1884;
   }
   tmp_1882 += tmp_1883;
   result += (0.7071067811865475*vS*Conj(Lambdax)*ZP(gI2,1)) * tmp_1882;
   std::complex<double> tmp_1886;
   std::complex<double> tmp_1887;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1888;
      std::complex<double> tmp_1889;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1890;
         std::complex<double> tmp_1891;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1891 += Conj(Yv(j1,j3))*Ye(j2,j1);
         }
         tmp_1890 += tmp_1891;
         tmp_1889 += (KroneckerDelta(gO2,3 + j2)) * tmp_1890;
      }
      tmp_1888 += tmp_1889;
      tmp_1887 += (Conj(ZV(gI1,3 + j3))) * tmp_1888;
   }
   tmp_1886 += tmp_1887;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_1886;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1892;
      std::complex<double> tmp_1893;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1893 += Conj(ZV(gI1,3 + j2))*Conj(TYv(gO2,j2));
      }
      tmp_1892 += tmp_1893;
      result += (ZP(gI2,1)) * tmp_1892;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1894;
      std::complex<double> tmp_1895;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1896;
         std::complex<double> tmp_1897;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1897 += Conj(Yv(gO2,j1))*Yv(j2,j1);
         }
         tmp_1896 += tmp_1897;
         tmp_1895 += (Conj(ZV(gI1,j2))) * tmp_1896;
      }
      tmp_1894 += tmp_1895;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_1894;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1898;
   std::complex<double> tmp_1899;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1900;
      std::complex<double> tmp_1901;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1901 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1900 += tmp_1901;
      tmp_1899 += (Conj(ZE(gI1,j2))) * tmp_1900;
   }
   tmp_1898 += tmp_1899;
   result += (std::complex<double>(0,-0.5)*vS*Conj(Lambdax)*Conj(ZA(gI2,1))) *
      tmp_1898;
   std::complex<double> tmp_1902;
   std::complex<double> tmp_1903;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1904;
      std::complex<double> tmp_1905;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1905 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1904 += tmp_1905;
      tmp_1903 += (Conj(ZE(gI1,j2))) * tmp_1904;
   }
   tmp_1902 += tmp_1903;
   result += (std::complex<double>(0,-0.5)*vu*Conj(Lambdax)*Conj(ZA(gI2,2))) *
      tmp_1902;
   std::complex<double> tmp_1906;
   std::complex<double> tmp_1907;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1908;
      std::complex<double> tmp_1909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1909 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_1908 += tmp_1909;
      tmp_1907 += (Conj(ZE(gI1,j2))) * tmp_1908;
   }
   tmp_1906 += tmp_1907;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_1906;
   if (gO2 < 3) {
      std::complex<double> tmp_1910;
      std::complex<double> tmp_1911;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1911 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1910 += tmp_1911;
      result += (std::complex<double>(0,0.5)*vS*Conj(ZA(gI2,1))*Lambdax) *
         tmp_1910;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1912;
      std::complex<double> tmp_1913;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1913 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1912 += tmp_1913;
      result += (std::complex<double>(0,0.5)*vu*Conj(ZA(gI2,2))*Lambdax) *
         tmp_1912;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1914;
      std::complex<double> tmp_1915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1915 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_1914 += tmp_1915;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))
         ) * tmp_1914;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_1916;
   std::complex<double> tmp_1917;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1917 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1916 += tmp_1917;
   result += (0.3*vd*Conj(ZH(gI2,0))*Sqr(g1)) * tmp_1916;
   std::complex<double> tmp_1918;
   std::complex<double> tmp_1919;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1919 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1918 += tmp_1919;
   result += (-(Qe*QHd*vd*Conj(ZH(gI2,0))*Sqr(gp))) * tmp_1918;
   std::complex<double> tmp_1920;
   std::complex<double> tmp_1921;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1921 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1920 += tmp_1921;
   result += (-0.3*vu*Conj(ZH(gI2,1))*Sqr(g1)) * tmp_1920;
   std::complex<double> tmp_1922;
   std::complex<double> tmp_1923;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1923 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1922 += tmp_1923;
   result += (-(Qe*QHu*vu*Conj(ZH(gI2,1))*Sqr(gp))) * tmp_1922;
   std::complex<double> tmp_1924;
   std::complex<double> tmp_1925;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1925 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1924 += tmp_1925;
   result += (-(Qe*Qs*vS*Conj(ZH(gI2,2))*Sqr(gp))) * tmp_1924;
   std::complex<double> tmp_1926;
   std::complex<double> tmp_1927;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1928;
      std::complex<double> tmp_1929;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1929 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1928 += tmp_1929;
      tmp_1927 += (Conj(ZE(gI1,j2))) * tmp_1928;
   }
   tmp_1926 += tmp_1927;
   result += (0.5*vS*Conj(Lambdax)*Conj(ZH(gI2,1))) * tmp_1926;
   std::complex<double> tmp_1930;
   std::complex<double> tmp_1931;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1932;
      std::complex<double> tmp_1933;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1933 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1932 += tmp_1933;
      tmp_1931 += (Conj(ZE(gI1,j2))) * tmp_1932;
   }
   tmp_1930 += tmp_1931;
   result += (0.5*vu*Conj(Lambdax)*Conj(ZH(gI2,2))) * tmp_1930;
   std::complex<double> tmp_1934;
   std::complex<double> tmp_1935;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1936;
      std::complex<double> tmp_1937;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1937 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_1936 += tmp_1937;
      tmp_1935 += (Conj(ZE(gI1,j2))) * tmp_1936;
   }
   tmp_1934 += tmp_1935;
   result += (-0.7071067811865475*Conj(ZH(gI2,0))) * tmp_1934;
   std::complex<double> tmp_1938;
   std::complex<double> tmp_1939;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1940;
      std::complex<double> tmp_1941;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1942;
         std::complex<double> tmp_1943;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1943 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1942 += tmp_1943;
         tmp_1941 += (KroneckerDelta(gO2,3 + j2)) * tmp_1942;
      }
      tmp_1940 += tmp_1941;
      tmp_1939 += (Conj(ZE(gI1,3 + j3))) * tmp_1940;
   }
   tmp_1938 += tmp_1939;
   result += (-(vd*Conj(ZH(gI2,0)))) * tmp_1938;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHd*Ql*vd*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,0))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(g1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(g2);
   }
   if (gO2 < 3) {
      result += -(QHu*Ql*vu*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,1))*Sqr(gp));
   }
   if (gO2 < 3) {
      result += -(Ql*Qs*vS*Conj(ZE(gI1,gO2))*Conj(ZH(gI2,2))*Sqr(gp));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1944;
      std::complex<double> tmp_1945;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1945 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1944 += tmp_1945;
      result += (0.5*vS*Conj(ZH(gI2,1))*Lambdax) * tmp_1944;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1946;
      std::complex<double> tmp_1947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1947 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1946 += tmp_1947;
      result += (0.5*vu*Conj(ZH(gI2,2))*Lambdax) * tmp_1946;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1948;
      std::complex<double> tmp_1949;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1949 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_1948 += tmp_1949;
      result += (-0.7071067811865475*Conj(ZH(gI2,0))) * tmp_1948;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1950;
      std::complex<double> tmp_1951;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1952;
         std::complex<double> tmp_1953;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1953 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_1952 += tmp_1953;
         tmp_1951 += (Conj(ZE(gI1,j2))) * tmp_1952;
      }
      tmp_1950 += tmp_1951;
      result += (-(vd*Conj(ZH(gI2,0)))) * tmp_1950;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeVPSe(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1954;
   std::complex<double> tmp_1955;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1955 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1954 += tmp_1955;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_1954;
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

   std::complex<double> tmp_1956;
   std::complex<double> tmp_1957;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1957 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1956 += tmp_1957;
   result += (0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())) * tmp_1956;
   std::complex<double> tmp_1958;
   std::complex<double> tmp_1959;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1959 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1958 += tmp_1959;
   result += (-(gp*Qe*Sin(ThetaWp()))) * tmp_1958;
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

   std::complex<double> tmp_1960;
   std::complex<double> tmp_1961;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1961 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1960 += tmp_1961;
   result += (-(gp*Qe*Cos(ThetaWp()))) * tmp_1960;
   std::complex<double> tmp_1962;
   std::complex<double> tmp_1963;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1963 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1962 += tmp_1963;
   result += (-0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())) * tmp_1962;
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

std::complex<double> CLASSNAME::CpconjUSeVWmSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7071067811865475*g2*Conj(ZV(gI2,gO2));
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

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(Conj(ZH(gI1,0))*(20*Conj(ZH(gI2,2))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2))*(AbsSqr
      (Lambdax) + QHd*Qs*Sqr(gp)) - Conj(ZH(gI2,1))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + Conj(ZH(gI2
      ,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*
      Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 3*KroneckerDelta(
      0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))
      )) + Conj(ZH(gI1,1))*(-20*Conj(ZH(gI2,2))*(KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1) + KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2))*(AbsSqr
      (Lambdax) + QHu*Qs*Sqr(gp)) + Conj(ZH(gI2,0))*(KroneckerDelta(0,gO2)*
      KroneckerDelta(1,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*(Sqr(g2) - 4*QHd*QHu*Sqr(gp))) + Conj(ZH(gI2
      ,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*
      Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - 3*KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))))
      - 20*Conj(ZH(gI1,2))*(Conj(ZH(gI2,0))*(KroneckerDelta(0,gO2)*KroneckerDelta
      (2,gO1) + KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) +
      QHd*Qs*Sqr(gp)) + Conj(ZH(gI2,1))*(KroneckerDelta(1,gO2)*KroneckerDelta(2,
      gO1) + KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2))*(AbsSqr(Lambdax) + QHu*
      Qs*Sqr(gp)) + Conj(ZH(gI2,2))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(
      AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + 3*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs))));

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

std::complex<double> CLASSNAME::CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.35355339059327373)*(Conj(ZA(gI2,2))*(
      Conj(ZH(gI1,1))*KroneckerDelta(0,gO2) + Conj(ZH(gI1,0))*KroneckerDelta(1,gO2
      )) + Conj(ZA(gI2,1))*(Conj(ZH(gI1,2))*KroneckerDelta(0,gO2) + Conj(ZH(gI1,0)
      )*KroneckerDelta(2,gO2)) + Conj(ZA(gI2,0))*(Conj(ZH(gI1,2))*KroneckerDelta(1
      ,gO2) + Conj(ZH(gI1,1))*KroneckerDelta(2,gO2)))*(Conj(TLambdax) - TLambdax);

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-5*Conj(ZH(gI1,2))*(-1.4142135623730951*Conj(TLambdax)*Conj(
      ZH(gI2,1))*KroneckerDelta(0,gO2) + 4*vd*AbsSqr(Lambdax)*Conj(ZH(gI2,2))*
      KroneckerDelta(0,gO2) + 4*vS*AbsSqr(Lambdax)*Conj(ZH(gI2,1))*KroneckerDelta(
      1,gO2) + 4*vu*AbsSqr(Lambdax)*Conj(ZH(gI2,2))*KroneckerDelta(1,gO2) + 4*vu*
      AbsSqr(Lambdax)*Conj(ZH(gI2,1))*KroneckerDelta(2,gO2) + 4*QHd*Qs*vd*Conj(ZH(
      gI2,2))*KroneckerDelta(0,gO2)*Sqr(gp) + 4*QHu*Qs*vS*Conj(ZH(gI2,1))*
      KroneckerDelta(1,gO2)*Sqr(gp) + 4*QHu*Qs*vu*Conj(ZH(gI2,2))*KroneckerDelta(1
      ,gO2)*Sqr(gp) + 4*QHu*Qs*vu*Conj(ZH(gI2,1))*KroneckerDelta(2,gO2)*Sqr(gp) +
      12*vS*Conj(ZH(gI2,2))*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs) -
      1.4142135623730951*Conj(ZH(gI2,1))*KroneckerDelta(0,gO2)*TLambdax + Conj(ZH(
      gI2,0))*(-1.4142135623730951*Conj(TLambdax)*KroneckerDelta(1,gO2) + 4*vd*
      AbsSqr(Lambdax)*KroneckerDelta(2,gO2) + 4*QHd*Qs*vd*KroneckerDelta(2,gO2)*
      Sqr(gp) + 4*vS*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) -
      1.4142135623730951*KroneckerDelta(1,gO2)*TLambdax)) + Conj(ZH(gI1,1))*(Conj(
      ZH(gI2,1))*(vd*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*
      Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*vS*KroneckerDelta(2,gO2)*(AbsSqr(Lambdax)
      + QHu*Qs*Sqr(gp)) - 3*vu*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*
      Sqr(gp)*Sqr(QHu))) + 5*Conj(ZH(gI2,2))*(1.4142135623730951*Conj(TLambdax)*
      KroneckerDelta(0,gO2) - 4*vu*AbsSqr(Lambdax)*KroneckerDelta(2,gO2) - 4*QHu*
      Qs*vu*KroneckerDelta(2,gO2)*Sqr(gp) - 4*vS*KroneckerDelta(1,gO2)*(AbsSqr(
      Lambdax) + QHu*Qs*Sqr(gp)) + 1.4142135623730951*KroneckerDelta(0,gO2)*
      TLambdax) + Conj(ZH(gI2,0))*(vu*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) +
      3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) + vd*KroneckerDelta(1,gO2)*(-20
      *AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) +
      7.0710678118654755*KroneckerDelta(2,gO2)*(Conj(TLambdax) + TLambdax))) -
      Conj(ZH(gI1,0))*(Conj(ZH(gI2,0))*(vu*KroneckerDelta(1,gO2)*(20*AbsSqr(
      Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 20*vS*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + 3*vd*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd)))) + 5*
      Conj(ZH(gI2,2))*(-1.4142135623730951*Conj(TLambdax)*KroneckerDelta(1,gO2) +
      4*vd*AbsSqr(Lambdax)*KroneckerDelta(2,gO2) + 4*QHd*Qs*vd*KroneckerDelta(2,
      gO2)*Sqr(gp) + 4*vS*KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp))
      - 1.4142135623730951*KroneckerDelta(1,gO2)*TLambdax) - Conj(ZH(gI2,1))*(vu*
      KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2) - 20*QHd*
      QHu*Sqr(gp)) + vd*KroneckerDelta(1,gO2)*(-20*AbsSqr(Lambdax) + 3*Sqr(g1) + 5
      *Sqr(g2) - 20*QHd*QHu*Sqr(gp)) + 7.0710678118654755*KroneckerDelta(2,gO2)*(
      Conj(TLambdax) + TLambdax))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1964;
   std::complex<double> tmp_1965;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1966;
      std::complex<double> tmp_1967;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1967 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1966 += tmp_1967;
      tmp_1965 += (ZDL(gI1,j2)) * tmp_1966;
   }
   tmp_1964 += tmp_1965;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1964;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1968;
   std::complex<double> tmp_1969;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1970;
      std::complex<double> tmp_1971;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1971 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1970 += tmp_1971;
      tmp_1969 += (Conj(ZDL(gI2,j2))) * tmp_1970;
   }
   tmp_1968 += tmp_1969;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1968;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1972;
   std::complex<double> tmp_1973;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1974;
      std::complex<double> tmp_1975;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1975 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1974 += tmp_1975;
      tmp_1973 += (ZEL(gI1,j2)) * tmp_1974;
   }
   tmp_1972 += tmp_1973;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1972;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1976;
   std::complex<double> tmp_1977;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1978;
      std::complex<double> tmp_1979;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1979 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1978 += tmp_1979;
      tmp_1977 += (Conj(ZEL(gI2,j2))) * tmp_1978;
   }
   tmp_1976 += tmp_1977;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1976;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1980;
   std::complex<double> tmp_1981;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1982;
      std::complex<double> tmp_1983;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1983 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1982 += tmp_1983;
      tmp_1981 += (ZUL(gI1,j2)) * tmp_1982;
   }
   tmp_1980 += tmp_1981;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1980;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1984;
   std::complex<double> tmp_1985;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1986;
      std::complex<double> tmp_1987;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1987 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1986 += tmp_1987;
      tmp_1985 += (Conj(ZUL(gI2,j2))) * tmp_1986;
   }
   tmp_1984 += tmp_1985;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_1984;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFvFvPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1988;
   std::complex<double> tmp_1989;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1990;
      std::complex<double> tmp_1991;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1991 += Conj(Yv(j1,j2))*ZVL(gI1,j1);
      }
      tmp_1990 += tmp_1991;
      tmp_1989 += (ZVR(gI2,j2)) * tmp_1990;
   }
   tmp_1988 += tmp_1989;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1988;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFvFvPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1992;
   std::complex<double> tmp_1993;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1994;
      std::complex<double> tmp_1995;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1995 += Conj(ZVL(gI2,j1))*Yv(j1,j2);
      }
      tmp_1994 += tmp_1995;
      tmp_1993 += (Conj(ZVR(gI1,j2))) * tmp_1994;
   }
   tmp_1992 += tmp_1993;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_1992;

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

   std::complex<double> tmp_1996;
   std::complex<double> tmp_1997;
   std::complex<double> tmp_1998;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1998 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1997 += tmp_1998;
   tmp_1996 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1997;
   std::complex<double> tmp_1999;
   std::complex<double> tmp_2000;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2000 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1999 += tmp_2000;
   tmp_1996 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1999;
   std::complex<double> tmp_2001;
   std::complex<double> tmp_2002;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2002 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2001 += tmp_2002;
   tmp_1996 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2001;
   std::complex<double> tmp_2003;
   std::complex<double> tmp_2004;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2004 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2003 += tmp_2004;
   tmp_1996 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2003;
   std::complex<double> tmp_2005;
   std::complex<double> tmp_2006;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2006 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2005 += tmp_2006;
   tmp_1996 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2005;
   std::complex<double> tmp_2007;
   std::complex<double> tmp_2008;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2008 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2007 += tmp_2008;
   tmp_1996 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2007;
   std::complex<double> tmp_2009;
   std::complex<double> tmp_2010;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2010 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2009 += tmp_2010;
   tmp_1996 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2009;
   std::complex<double> tmp_2011;
   std::complex<double> tmp_2012;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2012 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2011 += tmp_2012;
   tmp_1996 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2011;
   std::complex<double> tmp_2013;
   std::complex<double> tmp_2014;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2014 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2013 += tmp_2014;
   tmp_1996 += (std::complex<double>(0,-1)*Qd*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2013;
   std::complex<double> tmp_2015;
   std::complex<double> tmp_2016;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2016 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2015 += tmp_2016;
   tmp_1996 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2015;
   std::complex<double> tmp_2017;
   std::complex<double> tmp_2018;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2018 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2017 += tmp_2018;
   tmp_1996 += (std::complex<double>(0,-1)*Qd*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2017;
   std::complex<double> tmp_2019;
   std::complex<double> tmp_2020;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2020 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2019 += tmp_2020;
   tmp_1996 += (std::complex<double>(0,-1)*Qd*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2019;
   std::complex<double> tmp_2021;
   std::complex<double> tmp_2022;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2023;
      std::complex<double> tmp_2024;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2024 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2023 += tmp_2024;
      tmp_2022 += (Conj(ZD(gI2,j2))) * tmp_2023;
   }
   tmp_2021 += tmp_2022;
   tmp_1996 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1)) * tmp_2021;
   std::complex<double> tmp_2025;
   std::complex<double> tmp_2026;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2027;
      std::complex<double> tmp_2028;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2028 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2027 += tmp_2028;
      tmp_2026 += (Conj(ZD(gI2,j2))) * tmp_2027;
   }
   tmp_2025 += tmp_2026;
   tmp_1996 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1)
      *KroneckerDelta(2,gO2)) * tmp_2025;
   std::complex<double> tmp_2029;
   std::complex<double> tmp_2030;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2031;
      std::complex<double> tmp_2032;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2032 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2031 += tmp_2032;
      tmp_2030 += (ZD(gI1,j2)) * tmp_2031;
   }
   tmp_2029 += tmp_2030;
   tmp_1996 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2029;
   std::complex<double> tmp_2033;
   std::complex<double> tmp_2034;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2035;
      std::complex<double> tmp_2036;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2036 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2035 += tmp_2036;
      tmp_2034 += (ZD(gI1,j2)) * tmp_2035;
   }
   tmp_2033 += tmp_2034;
   tmp_1996 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2033;
   std::complex<double> tmp_2037;
   std::complex<double> tmp_2038;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2039;
      std::complex<double> tmp_2040;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2041;
         std::complex<double> tmp_2042;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2042 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_2041 += tmp_2042;
         tmp_2040 += (ZD(gI1,3 + j2)) * tmp_2041;
      }
      tmp_2039 += tmp_2040;
      tmp_2038 += (Conj(ZD(gI2,3 + j3))) * tmp_2039;
   }
   tmp_2037 += tmp_2038;
   tmp_1996 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2037;
   std::complex<double> tmp_2043;
   std::complex<double> tmp_2044;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2045;
      std::complex<double> tmp_2046;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2047;
         std::complex<double> tmp_2048;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2048 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2047 += tmp_2048;
         tmp_2046 += (Conj(ZD(gI2,j2))) * tmp_2047;
      }
      tmp_2045 += tmp_2046;
      tmp_2044 += (ZD(gI1,j3)) * tmp_2045;
   }
   tmp_2043 += tmp_2044;
   tmp_1996 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2043;
   result += (std::complex<double>(0,-1)) * tmp_1996;

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

   std::complex<double> tmp_2049;
   std::complex<double> tmp_2050;
   std::complex<double> tmp_2051;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2051 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2050 += tmp_2051;
   tmp_2049 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2050;
   std::complex<double> tmp_2052;
   std::complex<double> tmp_2053;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2053 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2052 += tmp_2053;
   tmp_2049 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2052;
   std::complex<double> tmp_2054;
   std::complex<double> tmp_2055;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2055 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2054 += tmp_2055;
   tmp_2049 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2054;
   std::complex<double> tmp_2056;
   std::complex<double> tmp_2057;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2057 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2056 += tmp_2057;
   tmp_2049 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2056;
   std::complex<double> tmp_2058;
   std::complex<double> tmp_2059;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2059 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2058 += tmp_2059;
   tmp_2049 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2058;
   std::complex<double> tmp_2060;
   std::complex<double> tmp_2061;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2061 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2060 += tmp_2061;
   tmp_2049 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2060;
   std::complex<double> tmp_2062;
   std::complex<double> tmp_2063;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2063 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2062 += tmp_2063;
   tmp_2049 += (std::complex<double>(0,-1)*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2062;
   std::complex<double> tmp_2064;
   std::complex<double> tmp_2065;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2065 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2064 += tmp_2065;
   tmp_2049 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2064;
   std::complex<double> tmp_2066;
   std::complex<double> tmp_2067;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2067 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2066 += tmp_2067;
   tmp_2049 += (std::complex<double>(0,-1)*Qe*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2066;
   std::complex<double> tmp_2068;
   std::complex<double> tmp_2069;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2069 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2068 += tmp_2069;
   tmp_2049 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2068;
   std::complex<double> tmp_2070;
   std::complex<double> tmp_2071;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2071 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2070 += tmp_2071;
   tmp_2049 += (std::complex<double>(0,-1)*Qe*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2070;
   std::complex<double> tmp_2072;
   std::complex<double> tmp_2073;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2073 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2072 += tmp_2073;
   tmp_2049 += (std::complex<double>(0,-1)*Qe*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2072;
   std::complex<double> tmp_2074;
   std::complex<double> tmp_2075;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2076;
      std::complex<double> tmp_2077;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2077 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2076 += tmp_2077;
      tmp_2075 += (Conj(ZE(gI2,j2))) * tmp_2076;
   }
   tmp_2074 += tmp_2075;
   tmp_2049 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2)
      *KroneckerDelta(2,gO1)) * tmp_2074;
   std::complex<double> tmp_2078;
   std::complex<double> tmp_2079;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2080;
      std::complex<double> tmp_2081;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2081 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2080 += tmp_2081;
      tmp_2079 += (Conj(ZE(gI2,j2))) * tmp_2080;
   }
   tmp_2078 += tmp_2079;
   tmp_2049 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1)
      *KroneckerDelta(2,gO2)) * tmp_2078;
   std::complex<double> tmp_2082;
   std::complex<double> tmp_2083;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2084;
      std::complex<double> tmp_2085;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2085 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2084 += tmp_2085;
      tmp_2083 += (ZE(gI1,j2)) * tmp_2084;
   }
   tmp_2082 += tmp_2083;
   tmp_2049 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2082;
   std::complex<double> tmp_2086;
   std::complex<double> tmp_2087;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2088;
      std::complex<double> tmp_2089;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2089 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2088 += tmp_2089;
      tmp_2087 += (ZE(gI1,j2)) * tmp_2088;
   }
   tmp_2086 += tmp_2087;
   tmp_2049 += (std::complex<double>(0,0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2086;
   std::complex<double> tmp_2090;
   std::complex<double> tmp_2091;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2092;
      std::complex<double> tmp_2093;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2094;
         std::complex<double> tmp_2095;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2095 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_2094 += tmp_2095;
         tmp_2093 += (ZE(gI1,3 + j2)) * tmp_2094;
      }
      tmp_2092 += tmp_2093;
      tmp_2091 += (Conj(ZE(gI2,3 + j3))) * tmp_2092;
   }
   tmp_2090 += tmp_2091;
   tmp_2049 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2090;
   std::complex<double> tmp_2096;
   std::complex<double> tmp_2097;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2098;
      std::complex<double> tmp_2099;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2100;
         std::complex<double> tmp_2101;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2101 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_2100 += tmp_2101;
         tmp_2099 += (Conj(ZE(gI2,j2))) * tmp_2100;
      }
      tmp_2098 += tmp_2099;
      tmp_2097 += (ZE(gI1,j3)) * tmp_2098;
   }
   tmp_2096 += tmp_2097;
   tmp_2049 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2096;
   result += (std::complex<double>(0,-1)) * tmp_2049;

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

   std::complex<double> tmp_2102;
   std::complex<double> tmp_2103;
   std::complex<double> tmp_2104;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2104 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2103 += tmp_2104;
   tmp_2102 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2103;
   std::complex<double> tmp_2105;
   std::complex<double> tmp_2106;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2106 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2105 += tmp_2106;
   tmp_2102 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2105;
   std::complex<double> tmp_2107;
   std::complex<double> tmp_2108;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2108 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2107 += tmp_2108;
   tmp_2102 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2107;
   std::complex<double> tmp_2109;
   std::complex<double> tmp_2110;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2110 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2109 += tmp_2110;
   tmp_2102 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2109;
   std::complex<double> tmp_2111;
   std::complex<double> tmp_2112;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2112 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2111 += tmp_2112;
   tmp_2102 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2111;
   std::complex<double> tmp_2113;
   std::complex<double> tmp_2114;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2114 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2113 += tmp_2114;
   tmp_2102 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2113;
   std::complex<double> tmp_2115;
   std::complex<double> tmp_2116;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2116 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2115 += tmp_2116;
   tmp_2102 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2115;
   std::complex<double> tmp_2117;
   std::complex<double> tmp_2118;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2118 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2117 += tmp_2118;
   tmp_2102 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2117;
   std::complex<double> tmp_2119;
   std::complex<double> tmp_2120;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2120 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2119 += tmp_2120;
   tmp_2102 += (std::complex<double>(0,-1)*QHd*Qu*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2119;
   std::complex<double> tmp_2121;
   std::complex<double> tmp_2122;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2122 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2121 += tmp_2122;
   tmp_2102 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2121;
   std::complex<double> tmp_2123;
   std::complex<double> tmp_2124;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2124 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2123 += tmp_2124;
   tmp_2102 += (std::complex<double>(0,-1)*QHu*Qu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2123;
   std::complex<double> tmp_2125;
   std::complex<double> tmp_2126;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2126 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2125 += tmp_2126;
   tmp_2102 += (std::complex<double>(0,-1)*Qs*Qu*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2125;
   std::complex<double> tmp_2127;
   std::complex<double> tmp_2128;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2129;
      std::complex<double> tmp_2130;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2130 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2129 += tmp_2130;
      tmp_2128 += (Conj(ZU(gI2,j2))) * tmp_2129;
   }
   tmp_2127 += tmp_2128;
   tmp_2102 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2)
      *KroneckerDelta(2,gO1)) * tmp_2127;
   std::complex<double> tmp_2131;
   std::complex<double> tmp_2132;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2133;
      std::complex<double> tmp_2134;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2134 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2133 += tmp_2134;
      tmp_2132 += (Conj(ZU(gI2,j2))) * tmp_2133;
   }
   tmp_2131 += tmp_2132;
   tmp_2102 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1)
      *KroneckerDelta(2,gO2)) * tmp_2131;
   std::complex<double> tmp_2135;
   std::complex<double> tmp_2136;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2137;
      std::complex<double> tmp_2138;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2138 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2137 += tmp_2138;
      tmp_2136 += (ZU(gI1,j2)) * tmp_2137;
   }
   tmp_2135 += tmp_2136;
   tmp_2102 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2135;
   std::complex<double> tmp_2139;
   std::complex<double> tmp_2140;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2141;
      std::complex<double> tmp_2142;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2142 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2141 += tmp_2142;
      tmp_2140 += (ZU(gI1,j2)) * tmp_2141;
   }
   tmp_2139 += tmp_2140;
   tmp_2102 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2139;
   std::complex<double> tmp_2143;
   std::complex<double> tmp_2144;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2145;
      std::complex<double> tmp_2146;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2147;
         std::complex<double> tmp_2148;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2148 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_2147 += tmp_2148;
         tmp_2146 += (ZU(gI1,3 + j2)) * tmp_2147;
      }
      tmp_2145 += tmp_2146;
      tmp_2144 += (Conj(ZU(gI2,3 + j3))) * tmp_2145;
   }
   tmp_2143 += tmp_2144;
   tmp_2102 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2143;
   std::complex<double> tmp_2149;
   std::complex<double> tmp_2150;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2151;
      std::complex<double> tmp_2152;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2153;
         std::complex<double> tmp_2154;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2154 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2153 += tmp_2154;
         tmp_2152 += (Conj(ZU(gI2,j2))) * tmp_2153;
      }
      tmp_2151 += tmp_2152;
      tmp_2150 += (ZU(gI1,j3)) * tmp_2151;
   }
   tmp_2149 += tmp_2150;
   tmp_2102 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2149;
   result += (std::complex<double>(0,-1)) * tmp_2102;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_2155;
   std::complex<double> tmp_2156;
   std::complex<double> tmp_2157;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2157 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2156 += tmp_2157;
   tmp_2155 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2156;
   std::complex<double> tmp_2158;
   std::complex<double> tmp_2159;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2159 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2158 += tmp_2159;
   tmp_2155 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2158;
   std::complex<double> tmp_2160;
   std::complex<double> tmp_2161;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2161 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2160 += tmp_2161;
   tmp_2155 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2160;
   std::complex<double> tmp_2162;
   std::complex<double> tmp_2163;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2163 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2162 += tmp_2163;
   tmp_2155 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2162;
   std::complex<double> tmp_2164;
   std::complex<double> tmp_2165;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2165 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2164 += tmp_2165;
   tmp_2155 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2164;
   std::complex<double> tmp_2166;
   std::complex<double> tmp_2167;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2167 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2166 += tmp_2167;
   tmp_2155 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2166;
   std::complex<double> tmp_2168;
   std::complex<double> tmp_2169;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2169 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2168 += tmp_2169;
   tmp_2155 += (std::complex<double>(0,-1)*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2168;
   std::complex<double> tmp_2170;
   std::complex<double> tmp_2171;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2171 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2170 += tmp_2171;
   tmp_2155 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2170;
   std::complex<double> tmp_2172;
   std::complex<double> tmp_2173;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2173 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2172 += tmp_2173;
   tmp_2155 += (std::complex<double>(0,-1)*QHd*Qv*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2172;
   std::complex<double> tmp_2174;
   std::complex<double> tmp_2175;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2175 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2174 += tmp_2175;
   tmp_2155 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2174;
   std::complex<double> tmp_2176;
   std::complex<double> tmp_2177;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2177 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2176 += tmp_2177;
   tmp_2155 += (std::complex<double>(0,-1)*QHu*Qv*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2176;
   std::complex<double> tmp_2178;
   std::complex<double> tmp_2179;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2179 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2178 += tmp_2179;
   tmp_2155 += (std::complex<double>(0,-1)*Qs*Qv*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2178;
   std::complex<double> tmp_2180;
   std::complex<double> tmp_2181;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2182;
      std::complex<double> tmp_2183;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2183 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2182 += tmp_2183;
      tmp_2181 += (Conj(ZV(gI2,3 + j2))) * tmp_2182;
   }
   tmp_2180 += tmp_2181;
   tmp_2155 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2180;
   std::complex<double> tmp_2184;
   std::complex<double> tmp_2185;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2186;
      std::complex<double> tmp_2187;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2187 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2186 += tmp_2187;
      tmp_2185 += (Conj(ZV(gI2,3 + j2))) * tmp_2186;
   }
   tmp_2184 += tmp_2185;
   tmp_2155 += (std::complex<double>(0,0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2184;
   std::complex<double> tmp_2188;
   std::complex<double> tmp_2189;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2190;
      std::complex<double> tmp_2191;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2191 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2190 += tmp_2191;
      tmp_2189 += (ZV(gI1,3 + j2)) * tmp_2190;
   }
   tmp_2188 += tmp_2189;
   tmp_2155 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2)
      *KroneckerDelta(2,gO1)) * tmp_2188;
   std::complex<double> tmp_2192;
   std::complex<double> tmp_2193;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2194;
      std::complex<double> tmp_2195;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2195 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2194 += tmp_2195;
      tmp_2193 += (ZV(gI1,3 + j2)) * tmp_2194;
   }
   tmp_2192 += tmp_2193;
   tmp_2155 += (std::complex<double>(0,0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1)
      *KroneckerDelta(2,gO2)) * tmp_2192;
   std::complex<double> tmp_2196;
   std::complex<double> tmp_2197;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2198;
      std::complex<double> tmp_2199;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2200;
         std::complex<double> tmp_2201;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2201 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_2200 += tmp_2201;
         tmp_2199 += (ZV(gI1,3 + j2)) * tmp_2200;
      }
      tmp_2198 += tmp_2199;
      tmp_2197 += (Conj(ZV(gI2,3 + j3))) * tmp_2198;
   }
   tmp_2196 += tmp_2197;
   tmp_2155 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2196;
   std::complex<double> tmp_2202;
   std::complex<double> tmp_2203;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2204;
      std::complex<double> tmp_2205;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2206;
         std::complex<double> tmp_2207;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2207 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_2206 += tmp_2207;
         tmp_2205 += (Conj(ZV(gI2,j2))) * tmp_2206;
      }
      tmp_2204 += tmp_2205;
      tmp_2203 += (ZV(gI1,j3)) * tmp_2204;
   }
   tmp_2202 += tmp_2203;
   tmp_2155 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2202;
   result += (std::complex<double>(0,-1)) * tmp_2155;

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

   std::complex<double> tmp_2208;
   std::complex<double> tmp_2209;
   std::complex<double> tmp_2210;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2210 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2209 += tmp_2210;
   tmp_2208 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_2209;
   std::complex<double> tmp_2211;
   std::complex<double> tmp_2212;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2212 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2211 += tmp_2212;
   tmp_2208 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_2211;
   std::complex<double> tmp_2213;
   std::complex<double> tmp_2214;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2214 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2213 += tmp_2214;
   tmp_2208 += (std::complex<double>(0,-1)*QHd*Qq*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2213;
   std::complex<double> tmp_2215;
   std::complex<double> tmp_2216;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2216 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2215 += tmp_2216;
   tmp_2208 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2215;
   std::complex<double> tmp_2217;
   std::complex<double> tmp_2218;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2218 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2217 += tmp_2218;
   tmp_2208 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_2217;
   std::complex<double> tmp_2219;
   std::complex<double> tmp_2220;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2220 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2219 += tmp_2220;
   tmp_2208 += (std::complex<double>(0,-1)*QHu*Qq*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2219;
   std::complex<double> tmp_2221;
   std::complex<double> tmp_2222;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2222 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2221 += tmp_2222;
   tmp_2208 += (std::complex<double>(0,-1)*Qq*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2221;
   std::complex<double> tmp_2223;
   std::complex<double> tmp_2224;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2224 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2223 += tmp_2224;
   tmp_2208 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_2223;
   std::complex<double> tmp_2225;
   std::complex<double> tmp_2226;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2226 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2225 += tmp_2226;
   tmp_2208 += (std::complex<double>(0,-1)*Qd*QHd*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2225;
   std::complex<double> tmp_2227;
   std::complex<double> tmp_2228;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2228 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2227 += tmp_2228;
   tmp_2208 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2227;
   std::complex<double> tmp_2229;
   std::complex<double> tmp_2230;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2230 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2229 += tmp_2230;
   tmp_2208 += (std::complex<double>(0,-1)*Qd*QHu*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2229;
   std::complex<double> tmp_2231;
   std::complex<double> tmp_2232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2232 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2231 += tmp_2232;
   tmp_2208 += (std::complex<double>(0,-1)*Qd*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2231;
   std::complex<double> tmp_2233;
   std::complex<double> tmp_2234;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2235;
      std::complex<double> tmp_2236;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2236 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2235 += tmp_2236;
      tmp_2234 += (Conj(ZD(gI2,j2))) * tmp_2235;
   }
   tmp_2233 += tmp_2234;
   tmp_2208 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(1,
      gO2)) * tmp_2233;
   std::complex<double> tmp_2237;
   std::complex<double> tmp_2238;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2239;
      std::complex<double> tmp_2240;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2240 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2239 += tmp_2240;
      tmp_2238 += (Conj(ZD(gI2,j2))) * tmp_2239;
   }
   tmp_2237 += tmp_2238;
   tmp_2208 += (std::complex<double>(0,0.5)*vu*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_2237;
   std::complex<double> tmp_2241;
   std::complex<double> tmp_2242;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2243;
      std::complex<double> tmp_2244;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2244 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_2243 += tmp_2244;
      tmp_2242 += (Conj(ZD(gI2,j2))) * tmp_2243;
   }
   tmp_2241 += tmp_2242;
   tmp_2208 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_2241;
   std::complex<double> tmp_2245;
   std::complex<double> tmp_2246;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2247;
      std::complex<double> tmp_2248;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2248 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2247 += tmp_2248;
      tmp_2246 += (ZD(gI1,j2)) * tmp_2247;
   }
   tmp_2245 += tmp_2246;
   tmp_2208 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_2245;
   std::complex<double> tmp_2249;
   std::complex<double> tmp_2250;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2251;
      std::complex<double> tmp_2252;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2252 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2251 += tmp_2252;
      tmp_2250 += (ZD(gI1,j2)) * tmp_2251;
   }
   tmp_2249 += tmp_2250;
   tmp_2208 += (std::complex<double>(0,0.5)*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_2249;
   std::complex<double> tmp_2253;
   std::complex<double> tmp_2254;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2255;
      std::complex<double> tmp_2256;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2256 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_2255 += tmp_2256;
      tmp_2254 += (ZD(gI1,j2)) * tmp_2255;
   }
   tmp_2253 += tmp_2254;
   tmp_2208 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_2253;
   std::complex<double> tmp_2257;
   std::complex<double> tmp_2258;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2259;
      std::complex<double> tmp_2260;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2261;
         std::complex<double> tmp_2262;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2262 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_2261 += tmp_2262;
         tmp_2260 += (ZD(gI1,3 + j2)) * tmp_2261;
      }
      tmp_2259 += tmp_2260;
      tmp_2258 += (Conj(ZD(gI2,3 + j3))) * tmp_2259;
   }
   tmp_2257 += tmp_2258;
   tmp_2208 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_2257
      ;
   std::complex<double> tmp_2263;
   std::complex<double> tmp_2264;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2265;
      std::complex<double> tmp_2266;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2267;
         std::complex<double> tmp_2268;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2268 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2267 += tmp_2268;
         tmp_2266 += (Conj(ZD(gI2,j2))) * tmp_2267;
      }
      tmp_2265 += tmp_2266;
      tmp_2264 += (ZD(gI1,j3)) * tmp_2265;
   }
   tmp_2263 += tmp_2264;
   tmp_2208 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_2263
      ;
   result += (std::complex<double>(0,-1)) * tmp_2208;

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

   std::complex<double> tmp_2269;
   std::complex<double> tmp_2270;
   std::complex<double> tmp_2271;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2271 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2270 += tmp_2271;
   tmp_2269 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_2270;
   std::complex<double> tmp_2272;
   std::complex<double> tmp_2273;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2273 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2272 += tmp_2273;
   tmp_2269 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_2272;
   std::complex<double> tmp_2274;
   std::complex<double> tmp_2275;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2275 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2274 += tmp_2275;
   tmp_2269 += (std::complex<double>(0,-1)*QHd*Ql*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2274;
   std::complex<double> tmp_2276;
   std::complex<double> tmp_2277;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2277 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2276 += tmp_2277;
   tmp_2269 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2276;
   std::complex<double> tmp_2278;
   std::complex<double> tmp_2279;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2279 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2278 += tmp_2279;
   tmp_2269 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_2278;
   std::complex<double> tmp_2280;
   std::complex<double> tmp_2281;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2281 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2280 += tmp_2281;
   tmp_2269 += (std::complex<double>(0,-1)*QHu*Ql*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2280;
   std::complex<double> tmp_2282;
   std::complex<double> tmp_2283;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2283 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2282 += tmp_2283;
   tmp_2269 += (std::complex<double>(0,-1)*Ql*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2282;
   std::complex<double> tmp_2284;
   std::complex<double> tmp_2285;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2285 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2284 += tmp_2285;
   tmp_2269 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_2284;
   std::complex<double> tmp_2286;
   std::complex<double> tmp_2287;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2287 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2286 += tmp_2287;
   tmp_2269 += (std::complex<double>(0,-1)*Qe*QHd*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2286;
   std::complex<double> tmp_2288;
   std::complex<double> tmp_2289;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2289 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2288 += tmp_2289;
   tmp_2269 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2288;
   std::complex<double> tmp_2290;
   std::complex<double> tmp_2291;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2291 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2290 += tmp_2291;
   tmp_2269 += (std::complex<double>(0,-1)*Qe*QHu*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2290;
   std::complex<double> tmp_2292;
   std::complex<double> tmp_2293;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2293 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2292 += tmp_2293;
   tmp_2269 += (std::complex<double>(0,-1)*Qe*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2292;
   std::complex<double> tmp_2294;
   std::complex<double> tmp_2295;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2296;
      std::complex<double> tmp_2297;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2297 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2296 += tmp_2297;
      tmp_2295 += (Conj(ZE(gI2,j2))) * tmp_2296;
   }
   tmp_2294 += tmp_2295;
   tmp_2269 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(1,
      gO2)) * tmp_2294;
   std::complex<double> tmp_2298;
   std::complex<double> tmp_2299;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2300;
      std::complex<double> tmp_2301;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2301 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2300 += tmp_2301;
      tmp_2299 += (Conj(ZE(gI2,j2))) * tmp_2300;
   }
   tmp_2298 += tmp_2299;
   tmp_2269 += (std::complex<double>(0,0.5)*vu*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_2298;
   std::complex<double> tmp_2302;
   std::complex<double> tmp_2303;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2304;
      std::complex<double> tmp_2305;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2305 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_2304 += tmp_2305;
      tmp_2303 += (Conj(ZE(gI2,j2))) * tmp_2304;
   }
   tmp_2302 += tmp_2303;
   tmp_2269 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_2302;
   std::complex<double> tmp_2306;
   std::complex<double> tmp_2307;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2308;
      std::complex<double> tmp_2309;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2309 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2308 += tmp_2309;
      tmp_2307 += (ZE(gI1,j2)) * tmp_2308;
   }
   tmp_2306 += tmp_2307;
   tmp_2269 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_2306;
   std::complex<double> tmp_2310;
   std::complex<double> tmp_2311;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2312;
      std::complex<double> tmp_2313;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2313 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2312 += tmp_2313;
      tmp_2311 += (ZE(gI1,j2)) * tmp_2312;
   }
   tmp_2310 += tmp_2311;
   tmp_2269 += (std::complex<double>(0,0.5)*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_2310;
   std::complex<double> tmp_2314;
   std::complex<double> tmp_2315;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2316;
      std::complex<double> tmp_2317;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2317 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_2316 += tmp_2317;
      tmp_2315 += (ZE(gI1,j2)) * tmp_2316;
   }
   tmp_2314 += tmp_2315;
   tmp_2269 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_2314;
   std::complex<double> tmp_2318;
   std::complex<double> tmp_2319;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2320;
      std::complex<double> tmp_2321;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2322;
         std::complex<double> tmp_2323;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2323 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_2322 += tmp_2323;
         tmp_2321 += (ZE(gI1,3 + j2)) * tmp_2322;
      }
      tmp_2320 += tmp_2321;
      tmp_2319 += (Conj(ZE(gI2,3 + j3))) * tmp_2320;
   }
   tmp_2318 += tmp_2319;
   tmp_2269 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_2318
      ;
   std::complex<double> tmp_2324;
   std::complex<double> tmp_2325;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2326;
      std::complex<double> tmp_2327;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2328;
         std::complex<double> tmp_2329;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2329 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_2328 += tmp_2329;
         tmp_2327 += (Conj(ZE(gI2,j2))) * tmp_2328;
      }
      tmp_2326 += tmp_2327;
      tmp_2325 += (ZE(gI1,j3)) * tmp_2326;
   }
   tmp_2324 += tmp_2325;
   tmp_2269 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_2324
      ;
   result += (std::complex<double>(0,-1)) * tmp_2269;

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

   std::complex<double> tmp_2330;
   std::complex<double> tmp_2331;
   std::complex<double> tmp_2332;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2332 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2331 += tmp_2332;
   tmp_2330 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_2331;
   std::complex<double> tmp_2333;
   std::complex<double> tmp_2334;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2334 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2333 += tmp_2334;
   tmp_2330 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_2333;
   std::complex<double> tmp_2335;
   std::complex<double> tmp_2336;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2336 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2335 += tmp_2336;
   tmp_2330 += (std::complex<double>(0,-1)*QHd*Qq*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2335;
   std::complex<double> tmp_2337;
   std::complex<double> tmp_2338;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2338 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2337 += tmp_2338;
   tmp_2330 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2337;
   std::complex<double> tmp_2339;
   std::complex<double> tmp_2340;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2340 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2339 += tmp_2340;
   tmp_2330 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_2339;
   std::complex<double> tmp_2341;
   std::complex<double> tmp_2342;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2342 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2341 += tmp_2342;
   tmp_2330 += (std::complex<double>(0,-1)*QHu*Qq*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2341;
   std::complex<double> tmp_2343;
   std::complex<double> tmp_2344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2344 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2343 += tmp_2344;
   tmp_2330 += (std::complex<double>(0,-1)*Qq*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2343;
   std::complex<double> tmp_2345;
   std::complex<double> tmp_2346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2346 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2345 += tmp_2346;
   tmp_2330 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_2345;
   std::complex<double> tmp_2347;
   std::complex<double> tmp_2348;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2348 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2347 += tmp_2348;
   tmp_2330 += (std::complex<double>(0,-1)*QHd*Qu*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2347;
   std::complex<double> tmp_2349;
   std::complex<double> tmp_2350;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2350 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2349 += tmp_2350;
   tmp_2330 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_2349;
   std::complex<double> tmp_2351;
   std::complex<double> tmp_2352;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2352 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2351 += tmp_2352;
   tmp_2330 += (std::complex<double>(0,-1)*QHu*Qu*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2351;
   std::complex<double> tmp_2353;
   std::complex<double> tmp_2354;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2354 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2353 += tmp_2354;
   tmp_2330 += (std::complex<double>(0,-1)*Qs*Qu*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2353;
   std::complex<double> tmp_2355;
   std::complex<double> tmp_2356;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2357;
      std::complex<double> tmp_2358;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2358 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2357 += tmp_2358;
      tmp_2356 += (Conj(ZU(gI2,j2))) * tmp_2357;
   }
   tmp_2355 += tmp_2356;
   tmp_2330 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(0,
      gO2)) * tmp_2355;
   std::complex<double> tmp_2359;
   std::complex<double> tmp_2360;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2361;
      std::complex<double> tmp_2362;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2362 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2361 += tmp_2362;
      tmp_2360 += (Conj(ZU(gI2,j2))) * tmp_2361;
   }
   tmp_2359 += tmp_2360;
   tmp_2330 += (std::complex<double>(0,0.5)*vd*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_2359;
   std::complex<double> tmp_2363;
   std::complex<double> tmp_2364;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2365;
      std::complex<double> tmp_2366;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2366 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_2365 += tmp_2366;
      tmp_2364 += (Conj(ZU(gI2,j2))) * tmp_2365;
   }
   tmp_2363 += tmp_2364;
   tmp_2330 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_2363;
   std::complex<double> tmp_2367;
   std::complex<double> tmp_2368;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2369;
      std::complex<double> tmp_2370;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2370 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2369 += tmp_2370;
      tmp_2368 += (ZU(gI1,j2)) * tmp_2369;
   }
   tmp_2367 += tmp_2368;
   tmp_2330 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(0,gO2)*Lambdax) *
      tmp_2367;
   std::complex<double> tmp_2371;
   std::complex<double> tmp_2372;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2373;
      std::complex<double> tmp_2374;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2374 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2373 += tmp_2374;
      tmp_2372 += (ZU(gI1,j2)) * tmp_2373;
   }
   tmp_2371 += tmp_2372;
   tmp_2330 += (std::complex<double>(0,0.5)*vd*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_2371;
   std::complex<double> tmp_2375;
   std::complex<double> tmp_2376;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2377;
      std::complex<double> tmp_2378;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2378 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_2377 += tmp_2378;
      tmp_2376 += (ZU(gI1,j2)) * tmp_2377;
   }
   tmp_2375 += tmp_2376;
   tmp_2330 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_2375;
   std::complex<double> tmp_2379;
   std::complex<double> tmp_2380;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2381;
      std::complex<double> tmp_2382;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2383;
         std::complex<double> tmp_2384;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2384 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_2383 += tmp_2384;
         tmp_2382 += (ZU(gI1,3 + j2)) * tmp_2383;
      }
      tmp_2381 += tmp_2382;
      tmp_2380 += (Conj(ZU(gI2,3 + j3))) * tmp_2381;
   }
   tmp_2379 += tmp_2380;
   tmp_2330 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_2379
      ;
   std::complex<double> tmp_2385;
   std::complex<double> tmp_2386;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2387;
      std::complex<double> tmp_2388;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2389;
         std::complex<double> tmp_2390;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2390 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2389 += tmp_2390;
         tmp_2388 += (Conj(ZU(gI2,j2))) * tmp_2389;
      }
      tmp_2387 += tmp_2388;
      tmp_2386 += (ZU(gI1,j3)) * tmp_2387;
   }
   tmp_2385 += tmp_2386;
   tmp_2330 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_2385
      ;
   result += (std::complex<double>(0,-1)) * tmp_2330;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_2391;
   std::complex<double> tmp_2392;
   std::complex<double> tmp_2393;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2393 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2392 += tmp_2393;
   tmp_2391 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_2392;
   std::complex<double> tmp_2394;
   std::complex<double> tmp_2395;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2395 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2394 += tmp_2395;
   tmp_2391 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_2394;
   std::complex<double> tmp_2396;
   std::complex<double> tmp_2397;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2397 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2396 += tmp_2397;
   tmp_2391 += (std::complex<double>(0,-1)*QHd*Ql*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2396;
   std::complex<double> tmp_2398;
   std::complex<double> tmp_2399;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2399 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2398 += tmp_2399;
   tmp_2391 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2398;
   std::complex<double> tmp_2400;
   std::complex<double> tmp_2401;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2401 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2400 += tmp_2401;
   tmp_2391 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_2400;
   std::complex<double> tmp_2402;
   std::complex<double> tmp_2403;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2403 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2402 += tmp_2403;
   tmp_2391 += (std::complex<double>(0,-1)*QHu*Ql*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2402;
   std::complex<double> tmp_2404;
   std::complex<double> tmp_2405;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2405 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2404 += tmp_2405;
   tmp_2391 += (std::complex<double>(0,-1)*Ql*Qs*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2404;
   std::complex<double> tmp_2406;
   std::complex<double> tmp_2407;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2407 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2406 += tmp_2407;
   tmp_2391 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_2406;
   std::complex<double> tmp_2408;
   std::complex<double> tmp_2409;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2409 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2408 += tmp_2409;
   tmp_2391 += (std::complex<double>(0,-1)*QHd*Qv*vd*KroneckerDelta(0,gO2)*Sqr(
      gp)) * tmp_2408;
   std::complex<double> tmp_2410;
   std::complex<double> tmp_2411;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2411 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2410 += tmp_2411;
   tmp_2391 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_2410;
   std::complex<double> tmp_2412;
   std::complex<double> tmp_2413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2413 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2412 += tmp_2413;
   tmp_2391 += (std::complex<double>(0,-1)*QHu*Qv*vu*KroneckerDelta(1,gO2)*Sqr(
      gp)) * tmp_2412;
   std::complex<double> tmp_2414;
   std::complex<double> tmp_2415;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2415 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2414 += tmp_2415;
   tmp_2391 += (std::complex<double>(0,-1)*Qs*Qv*vS*KroneckerDelta(2,gO2)*Sqr(
      gp)) * tmp_2414;
   std::complex<double> tmp_2416;
   std::complex<double> tmp_2417;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2418;
      std::complex<double> tmp_2419;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2419 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2418 += tmp_2419;
      tmp_2417 += (Conj(ZV(gI2,3 + j2))) * tmp_2418;
   }
   tmp_2416 += tmp_2417;
   tmp_2391 += (std::complex<double>(0,0.5)*vS*KroneckerDelta(0,gO2)*Lambdax) *
      tmp_2416;
   std::complex<double> tmp_2420;
   std::complex<double> tmp_2421;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2422;
      std::complex<double> tmp_2423;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2423 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2422 += tmp_2423;
      tmp_2421 += (Conj(ZV(gI2,3 + j2))) * tmp_2422;
   }
   tmp_2420 += tmp_2421;
   tmp_2391 += (std::complex<double>(0,0.5)*vd*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_2420;
   std::complex<double> tmp_2424;
   std::complex<double> tmp_2425;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2426;
      std::complex<double> tmp_2427;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2427 += Conj(TYv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2426 += tmp_2427;
      tmp_2425 += (Conj(ZV(gI2,3 + j2))) * tmp_2426;
   }
   tmp_2424 += tmp_2425;
   tmp_2391 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_2424;
   std::complex<double> tmp_2428;
   std::complex<double> tmp_2429;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2430;
      std::complex<double> tmp_2431;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2431 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2430 += tmp_2431;
      tmp_2429 += (ZV(gI1,3 + j2)) * tmp_2430;
   }
   tmp_2428 += tmp_2429;
   tmp_2391 += (std::complex<double>(0,0.5)*vS*Conj(Lambdax)*KroneckerDelta(0,
      gO2)) * tmp_2428;
   std::complex<double> tmp_2432;
   std::complex<double> tmp_2433;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2434;
      std::complex<double> tmp_2435;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2435 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2434 += tmp_2435;
      tmp_2433 += (ZV(gI1,3 + j2)) * tmp_2434;
   }
   tmp_2432 += tmp_2433;
   tmp_2391 += (std::complex<double>(0,0.5)*vd*Conj(Lambdax)*KroneckerDelta(2,
      gO2)) * tmp_2432;
   std::complex<double> tmp_2436;
   std::complex<double> tmp_2437;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2438;
      std::complex<double> tmp_2439;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2439 += Conj(ZV(gI2,j1))*TYv(j1,j2);
      }
      tmp_2438 += tmp_2439;
      tmp_2437 += (ZV(gI1,3 + j2)) * tmp_2438;
   }
   tmp_2436 += tmp_2437;
   tmp_2391 += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_2436;
   std::complex<double> tmp_2440;
   std::complex<double> tmp_2441;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2442;
      std::complex<double> tmp_2443;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2444;
         std::complex<double> tmp_2445;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2445 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_2444 += tmp_2445;
         tmp_2443 += (ZV(gI1,3 + j2)) * tmp_2444;
      }
      tmp_2442 += tmp_2443;
      tmp_2441 += (Conj(ZV(gI2,3 + j3))) * tmp_2442;
   }
   tmp_2440 += tmp_2441;
   tmp_2391 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_2440
      ;
   std::complex<double> tmp_2446;
   std::complex<double> tmp_2447;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2448;
      std::complex<double> tmp_2449;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2450;
         std::complex<double> tmp_2451;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2451 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_2450 += tmp_2451;
         tmp_2449 += (Conj(ZV(gI2,j2))) * tmp_2450;
      }
      tmp_2448 += tmp_2449;
      tmp_2447 += (ZV(gI1,j3)) * tmp_2448;
   }
   tmp_2446 += tmp_2447;
   tmp_2391 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_2446
      ;
   result += (std::complex<double>(0,-1)) * tmp_2391;

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

   result = std::complex<double>(0.,-0.7071067811865475)*(g2*KroneckerDelta(0,
      gO2)*UM(gI1,1)*UP(gI2,0) + (g2*KroneckerDelta(1,gO2)*UM(gI1,0) - Conj(
      Lambdax)*KroneckerDelta(2,gO2)*UM(gI1,1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gI2,0))*
      Conj(UP(gI1,1))*KroneckerDelta(1,gO1) + Conj(UM(gI2,1))*(g2*Conj(UP(gI1,0))*
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

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-(Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2) + 20*QHd*
      QHu*Sqr(gp)) + 20*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*(AbsSqr(
      Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*
      Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))))) + Conj(ZH(gI1,1))*Conj(ZH(gI2,
      1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-20*AbsSqr(Lambdax) + 3*
      Sqr(g1) + 5*Sqr(g2) - 20*QHd*QHu*Sqr(gp)) - 20*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) - KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu))) -
      20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,
      gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp)) + KroneckerDelta(2,
      gO1)*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.35355339059327373)*(Conj(ZA(gI1,2))*(Conj
      (ZA(gI2,1))*KroneckerDelta(0,gO2) + Conj(ZA(gI2,0))*KroneckerDelta(1,gO2)) +
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

   result = 0.05*(-20*vd*AbsSqr(Lambdax)*Conj(ZA(gI2,1))*Conj(ZH(gI1,0))*
      KroneckerDelta(1,gO2) - 20*vS*AbsSqr(Lambdax)*Conj(ZA(gI2,1))*Conj(ZH(gI1,2)
      )*KroneckerDelta(1,gO2) - 20*vd*AbsSqr(Lambdax)*Conj(ZA(gI2,2))*Conj(ZH(gI1,
      0))*KroneckerDelta(2,gO2) - 20*vu*AbsSqr(Lambdax)*Conj(ZA(gI2,2))*Conj(ZH(
      gI1,1))*KroneckerDelta(2,gO2) - 7.0710678118654755*Conj(TLambdax)*(Conj(ZA(
      gI2,2))*(Conj(ZH(gI1,1))*KroneckerDelta(0,gO2) + Conj(ZH(gI1,0))*
      KroneckerDelta(1,gO2)) + Conj(ZA(gI2,1))*(Conj(ZH(gI1,2))*KroneckerDelta(0,
      gO2) + Conj(ZH(gI1,0))*KroneckerDelta(2,gO2))) + 3*vd*Conj(ZA(gI2,1))*Conj(
      ZH(gI1,0))*KroneckerDelta(1,gO2)*Sqr(g1) - 3*vu*Conj(ZA(gI2,1))*Conj(ZH(gI1,
      1))*KroneckerDelta(1,gO2)*Sqr(g1) + 5*vd*Conj(ZA(gI2,1))*Conj(ZH(gI1,0))*
      KroneckerDelta(1,gO2)*Sqr(g2) - 5*vu*Conj(ZA(gI2,1))*Conj(ZH(gI1,1))*
      KroneckerDelta(1,gO2)*Sqr(g2) - 20*QHd*QHu*vd*Conj(ZA(gI2,1))*Conj(ZH(gI1,0)
      )*KroneckerDelta(1,gO2)*Sqr(gp) - 20*QHu*Qs*vS*Conj(ZA(gI2,1))*Conj(ZH(gI1,2
      ))*KroneckerDelta(1,gO2)*Sqr(gp) - 20*QHd*Qs*vd*Conj(ZA(gI2,2))*Conj(ZH(gI1,
      0))*KroneckerDelta(2,gO2)*Sqr(gp) - 20*QHu*Qs*vu*Conj(ZA(gI2,2))*Conj(ZH(gI1
      ,1))*KroneckerDelta(2,gO2)*Sqr(gp) - 20*vu*Conj(ZA(gI2,1))*Conj(ZH(gI1,1))*
      KroneckerDelta(1,gO2)*Sqr(gp)*Sqr(QHu) - 20*vS*Conj(ZA(gI2,2))*Conj(ZH(gI1,2
      ))*KroneckerDelta(2,gO2)*Sqr(gp)*Sqr(Qs) - 7.0710678118654755*Conj(ZA(gI2,2)
      )*Conj(ZH(gI1,1))*KroneckerDelta(0,gO2)*TLambdax - 7.0710678118654755*Conj(
      ZA(gI2,1))*Conj(ZH(gI1,2))*KroneckerDelta(0,gO2)*TLambdax -
      7.0710678118654755*Conj(ZA(gI2,2))*Conj(ZH(gI1,0))*KroneckerDelta(1,gO2)*
      TLambdax - 7.0710678118654755*Conj(ZA(gI2,1))*Conj(ZH(gI1,0))*KroneckerDelta
      (2,gO2)*TLambdax - Conj(ZA(gI2,0))*(vd*Conj(ZH(gI1,0))*KroneckerDelta(0,gO2)
      *(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))) + 5*Conj(ZH(gI1,2))*(4*vS*
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) +
      1.4142135623730951*KroneckerDelta(1,gO2)*(Conj(TLambdax) + TLambdax)) + Conj
      (ZH(gI1,1))*(vu*KroneckerDelta(0,gO2)*(20*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*
      Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + 7.0710678118654755*KroneckerDelta(2,gO2)*(
      Conj(TLambdax) + TLambdax))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.35355339059327373)*(Conj(ZH(gI1,2))*(
      Conj(ZH(gI2,1))*KroneckerDelta(0,gO2) + Conj(ZH(gI2,0))*KroneckerDelta(1,gO2
      )) + Conj(ZH(gI1,1))*(Conj(ZH(gI2,2))*KroneckerDelta(0,gO2) + Conj(ZH(gI2,0)
      )*KroneckerDelta(2,gO2)) + Conj(ZH(gI1,0))*(Conj(ZH(gI2,2))*KroneckerDelta(1
      ,gO2) + Conj(ZH(gI2,1))*KroneckerDelta(2,gO2)))*(Conj(TLambdax) - TLambdax);

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2452;
   std::complex<double> tmp_2453;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2454;
      std::complex<double> tmp_2455;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2455 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2454 += tmp_2455;
      tmp_2453 += (ZDL(gI1,j2)) * tmp_2454;
   }
   tmp_2452 += tmp_2453;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_2452;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2456;
   std::complex<double> tmp_2457;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2458;
      std::complex<double> tmp_2459;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2459 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_2458 += tmp_2459;
      tmp_2457 += (Conj(ZDL(gI2,j2))) * tmp_2458;
   }
   tmp_2456 += tmp_2457;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_2456;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2460;
   std::complex<double> tmp_2461;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2462;
      std::complex<double> tmp_2463;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2463 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2462 += tmp_2463;
      tmp_2461 += (ZEL(gI1,j2)) * tmp_2462;
   }
   tmp_2460 += tmp_2461;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(0,gO2)
      ) * tmp_2460;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2464;
   std::complex<double> tmp_2465;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2466;
      std::complex<double> tmp_2467;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2467 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_2466 += tmp_2467;
      tmp_2465 += (Conj(ZEL(gI2,j2))) * tmp_2466;
   }
   tmp_2464 += tmp_2465;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(0,gO1
      )) * tmp_2464;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2468;
   std::complex<double> tmp_2469;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2470;
      std::complex<double> tmp_2471;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2471 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2470 += tmp_2471;
      tmp_2469 += (ZUL(gI1,j2)) * tmp_2470;
   }
   tmp_2468 += tmp_2469;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO2)
      ) * tmp_2468;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2472;
   std::complex<double> tmp_2473;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2474;
      std::complex<double> tmp_2475;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2475 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_2474 += tmp_2475;
      tmp_2473 += (Conj(ZUL(gI2,j2))) * tmp_2474;
   }
   tmp_2472 += tmp_2473;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO1
      )) * tmp_2472;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFvFvPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2476;
   std::complex<double> tmp_2477;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2478;
      std::complex<double> tmp_2479;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2479 += Conj(Yv(j1,j2))*ZVL(gI1,j1);
      }
      tmp_2478 += tmp_2479;
      tmp_2477 += (ZVR(gI2,j2)) * tmp_2478;
   }
   tmp_2476 += tmp_2477;
   result += (std::complex<double>(0.,0.7071067811865475)*KroneckerDelta(1,gO2)
      ) * tmp_2476;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFvFvPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2480;
   std::complex<double> tmp_2481;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2482;
      std::complex<double> tmp_2483;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2483 += Conj(ZVL(gI2,j1))*Yv(j1,j2);
      }
      tmp_2482 += tmp_2483;
      tmp_2481 += (Conj(ZVR(gI1,j2))) * tmp_2482;
   }
   tmp_2480 += tmp_2481;
   result += (std::complex<double>(0.,-0.7071067811865475)*KroneckerDelta(1,gO1
      )) * tmp_2480;

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

   std::complex<double> tmp_2484;
   std::complex<double> tmp_2485;
   std::complex<double> tmp_2486;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2486 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2485 += tmp_2486;
   tmp_2484 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2485;
   std::complex<double> tmp_2487;
   std::complex<double> tmp_2488;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2488 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2487 += tmp_2488;
   tmp_2484 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2487;
   std::complex<double> tmp_2489;
   std::complex<double> tmp_2490;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2490 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2489 += tmp_2490;
   tmp_2484 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2489;
   std::complex<double> tmp_2491;
   std::complex<double> tmp_2492;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2492 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2491 += tmp_2492;
   tmp_2484 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2491;
   std::complex<double> tmp_2493;
   std::complex<double> tmp_2494;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2494 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2493 += tmp_2494;
   tmp_2484 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2493;
   std::complex<double> tmp_2495;
   std::complex<double> tmp_2496;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2496 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2495 += tmp_2496;
   tmp_2484 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2495;
   std::complex<double> tmp_2497;
   std::complex<double> tmp_2498;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2498 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2497 += tmp_2498;
   tmp_2484 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2497;
   std::complex<double> tmp_2499;
   std::complex<double> tmp_2500;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2500 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2499 += tmp_2500;
   tmp_2484 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2499;
   std::complex<double> tmp_2501;
   std::complex<double> tmp_2502;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2502 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2501 += tmp_2502;
   tmp_2484 += (std::complex<double>(0,-1)*Qd*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2501;
   std::complex<double> tmp_2503;
   std::complex<double> tmp_2504;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2504 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2503 += tmp_2504;
   tmp_2484 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2503;
   std::complex<double> tmp_2505;
   std::complex<double> tmp_2506;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2506 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2505 += tmp_2506;
   tmp_2484 += (std::complex<double>(0,-1)*Qd*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2505;
   std::complex<double> tmp_2507;
   std::complex<double> tmp_2508;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2508 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2507 += tmp_2508;
   tmp_2484 += (std::complex<double>(0,-1)*Qd*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2507;
   std::complex<double> tmp_2509;
   std::complex<double> tmp_2510;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2511;
      std::complex<double> tmp_2512;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2512 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2511 += tmp_2512;
      tmp_2510 += (Conj(ZD(gI2,j2))) * tmp_2511;
   }
   tmp_2509 += tmp_2510;
   tmp_2484 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)) * tmp_2509;
   std::complex<double> tmp_2513;
   std::complex<double> tmp_2514;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2515;
      std::complex<double> tmp_2516;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2516 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2515 += tmp_2516;
      tmp_2514 += (Conj(ZD(gI2,j2))) * tmp_2515;
   }
   tmp_2513 += tmp_2514;
   tmp_2484 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1
      )*KroneckerDelta(2,gO2)) * tmp_2513;
   std::complex<double> tmp_2517;
   std::complex<double> tmp_2518;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2519;
      std::complex<double> tmp_2520;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2520 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2519 += tmp_2520;
      tmp_2518 += (ZD(gI1,j2)) * tmp_2519;
   }
   tmp_2517 += tmp_2518;
   tmp_2484 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2517;
   std::complex<double> tmp_2521;
   std::complex<double> tmp_2522;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2523;
      std::complex<double> tmp_2524;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2524 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2523 += tmp_2524;
      tmp_2522 += (ZD(gI1,j2)) * tmp_2523;
   }
   tmp_2521 += tmp_2522;
   tmp_2484 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2521;
   std::complex<double> tmp_2525;
   std::complex<double> tmp_2526;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2527;
      std::complex<double> tmp_2528;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2529;
         std::complex<double> tmp_2530;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2530 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_2529 += tmp_2530;
         tmp_2528 += (ZD(gI1,3 + j2)) * tmp_2529;
      }
      tmp_2527 += tmp_2528;
      tmp_2526 += (Conj(ZD(gI2,3 + j3))) * tmp_2527;
   }
   tmp_2525 += tmp_2526;
   tmp_2484 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2525;
   std::complex<double> tmp_2531;
   std::complex<double> tmp_2532;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2533;
      std::complex<double> tmp_2534;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2535;
         std::complex<double> tmp_2536;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2536 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2535 += tmp_2536;
         tmp_2534 += (Conj(ZD(gI2,j2))) * tmp_2535;
      }
      tmp_2533 += tmp_2534;
      tmp_2532 += (ZD(gI1,j3)) * tmp_2533;
   }
   tmp_2531 += tmp_2532;
   tmp_2484 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2531;
   result += (std::complex<double>(0,-1)) * tmp_2484;

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

   std::complex<double> tmp_2537;
   std::complex<double> tmp_2538;
   std::complex<double> tmp_2539;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2539 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2538 += tmp_2539;
   tmp_2537 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2538;
   std::complex<double> tmp_2540;
   std::complex<double> tmp_2541;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2541 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2540 += tmp_2541;
   tmp_2537 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2540;
   std::complex<double> tmp_2542;
   std::complex<double> tmp_2543;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2543 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2542 += tmp_2543;
   tmp_2537 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2542;
   std::complex<double> tmp_2544;
   std::complex<double> tmp_2545;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2545 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2544 += tmp_2545;
   tmp_2537 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2544;
   std::complex<double> tmp_2546;
   std::complex<double> tmp_2547;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2547 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2546 += tmp_2547;
   tmp_2537 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2546;
   std::complex<double> tmp_2548;
   std::complex<double> tmp_2549;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2549 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2548 += tmp_2549;
   tmp_2537 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2548;
   std::complex<double> tmp_2550;
   std::complex<double> tmp_2551;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2551 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2550 += tmp_2551;
   tmp_2537 += (std::complex<double>(0,-1)*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2550;
   std::complex<double> tmp_2552;
   std::complex<double> tmp_2553;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2553 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2552 += tmp_2553;
   tmp_2537 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2552;
   std::complex<double> tmp_2554;
   std::complex<double> tmp_2555;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2555 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2554 += tmp_2555;
   tmp_2537 += (std::complex<double>(0,-1)*Qe*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2554;
   std::complex<double> tmp_2556;
   std::complex<double> tmp_2557;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2557 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2556 += tmp_2557;
   tmp_2537 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2556;
   std::complex<double> tmp_2558;
   std::complex<double> tmp_2559;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2559 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2558 += tmp_2559;
   tmp_2537 += (std::complex<double>(0,-1)*Qe*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2558;
   std::complex<double> tmp_2560;
   std::complex<double> tmp_2561;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2561 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2560 += tmp_2561;
   tmp_2537 += (std::complex<double>(0,-1)*Qe*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2560;
   std::complex<double> tmp_2562;
   std::complex<double> tmp_2563;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2564;
      std::complex<double> tmp_2565;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2565 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2564 += tmp_2565;
      tmp_2563 += (Conj(ZE(gI2,j2))) * tmp_2564;
   }
   tmp_2562 += tmp_2563;
   tmp_2537 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)) * tmp_2562;
   std::complex<double> tmp_2566;
   std::complex<double> tmp_2567;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2568;
      std::complex<double> tmp_2569;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2569 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2568 += tmp_2569;
      tmp_2567 += (Conj(ZE(gI2,j2))) * tmp_2568;
   }
   tmp_2566 += tmp_2567;
   tmp_2537 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(1,gO1
      )*KroneckerDelta(2,gO2)) * tmp_2566;
   std::complex<double> tmp_2570;
   std::complex<double> tmp_2571;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2572;
      std::complex<double> tmp_2573;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2573 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2572 += tmp_2573;
      tmp_2571 += (ZE(gI1,j2)) * tmp_2572;
   }
   tmp_2570 += tmp_2571;
   tmp_2537 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2570;
   std::complex<double> tmp_2574;
   std::complex<double> tmp_2575;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2576;
      std::complex<double> tmp_2577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2577 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2576 += tmp_2577;
      tmp_2575 += (ZE(gI1,j2)) * tmp_2576;
   }
   tmp_2574 += tmp_2575;
   tmp_2537 += (std::complex<double>(0,-0.5)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2574;
   std::complex<double> tmp_2578;
   std::complex<double> tmp_2579;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2580;
      std::complex<double> tmp_2581;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2582;
         std::complex<double> tmp_2583;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2583 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_2582 += tmp_2583;
         tmp_2581 += (ZE(gI1,3 + j2)) * tmp_2582;
      }
      tmp_2580 += tmp_2581;
      tmp_2579 += (Conj(ZE(gI2,3 + j3))) * tmp_2580;
   }
   tmp_2578 += tmp_2579;
   tmp_2537 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2578;
   std::complex<double> tmp_2584;
   std::complex<double> tmp_2585;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2586;
      std::complex<double> tmp_2587;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2588;
         std::complex<double> tmp_2589;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2589 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_2588 += tmp_2589;
         tmp_2587 += (Conj(ZE(gI2,j2))) * tmp_2588;
      }
      tmp_2586 += tmp_2587;
      tmp_2585 += (ZE(gI1,j3)) * tmp_2586;
   }
   tmp_2584 += tmp_2585;
   tmp_2537 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2584;
   result += (std::complex<double>(0,-1)) * tmp_2537;

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

   std::complex<double> tmp_2590;
   std::complex<double> tmp_2591;
   std::complex<double> tmp_2592;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2592 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2591 += tmp_2592;
   tmp_2590 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2591;
   std::complex<double> tmp_2593;
   std::complex<double> tmp_2594;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2594 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2593 += tmp_2594;
   tmp_2590 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2593;
   std::complex<double> tmp_2595;
   std::complex<double> tmp_2596;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2596 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2595 += tmp_2596;
   tmp_2590 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2595;
   std::complex<double> tmp_2597;
   std::complex<double> tmp_2598;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2598 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2597 += tmp_2598;
   tmp_2590 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2597;
   std::complex<double> tmp_2599;
   std::complex<double> tmp_2600;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2600 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2599 += tmp_2600;
   tmp_2590 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2599;
   std::complex<double> tmp_2601;
   std::complex<double> tmp_2602;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2602 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2601 += tmp_2602;
   tmp_2590 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2601;
   std::complex<double> tmp_2603;
   std::complex<double> tmp_2604;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2604 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2603 += tmp_2604;
   tmp_2590 += (std::complex<double>(0,-1)*Qq*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2603;
   std::complex<double> tmp_2605;
   std::complex<double> tmp_2606;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2606 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2605 += tmp_2606;
   tmp_2590 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2605;
   std::complex<double> tmp_2607;
   std::complex<double> tmp_2608;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2608 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2607 += tmp_2608;
   tmp_2590 += (std::complex<double>(0,-1)*QHd*Qu*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2607;
   std::complex<double> tmp_2609;
   std::complex<double> tmp_2610;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2610 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2609 += tmp_2610;
   tmp_2590 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2609;
   std::complex<double> tmp_2611;
   std::complex<double> tmp_2612;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2612 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2611 += tmp_2612;
   tmp_2590 += (std::complex<double>(0,-1)*QHu*Qu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2611;
   std::complex<double> tmp_2613;
   std::complex<double> tmp_2614;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2614 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2613 += tmp_2614;
   tmp_2590 += (std::complex<double>(0,-1)*Qs*Qu*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2613;
   std::complex<double> tmp_2615;
   std::complex<double> tmp_2616;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2617;
      std::complex<double> tmp_2618;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2618 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2617 += tmp_2618;
      tmp_2616 += (Conj(ZU(gI2,j2))) * tmp_2617;
   }
   tmp_2615 += tmp_2616;
   tmp_2590 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2
      )*KroneckerDelta(2,gO1)) * tmp_2615;
   std::complex<double> tmp_2619;
   std::complex<double> tmp_2620;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2621;
      std::complex<double> tmp_2622;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2622 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2621 += tmp_2622;
      tmp_2620 += (Conj(ZU(gI2,j2))) * tmp_2621;
   }
   tmp_2619 += tmp_2620;
   tmp_2590 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1
      )*KroneckerDelta(2,gO2)) * tmp_2619;
   std::complex<double> tmp_2623;
   std::complex<double> tmp_2624;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2625;
      std::complex<double> tmp_2626;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2626 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2625 += tmp_2626;
      tmp_2624 += (ZU(gI1,j2)) * tmp_2625;
   }
   tmp_2623 += tmp_2624;
   tmp_2590 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2623;
   std::complex<double> tmp_2627;
   std::complex<double> tmp_2628;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2629;
      std::complex<double> tmp_2630;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2630 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2629 += tmp_2630;
      tmp_2628 += (ZU(gI1,j2)) * tmp_2629;
   }
   tmp_2627 += tmp_2628;
   tmp_2590 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2627;
   std::complex<double> tmp_2631;
   std::complex<double> tmp_2632;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2633;
      std::complex<double> tmp_2634;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2635;
         std::complex<double> tmp_2636;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2636 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_2635 += tmp_2636;
         tmp_2634 += (ZU(gI1,3 + j2)) * tmp_2635;
      }
      tmp_2633 += tmp_2634;
      tmp_2632 += (Conj(ZU(gI2,3 + j3))) * tmp_2633;
   }
   tmp_2631 += tmp_2632;
   tmp_2590 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2631;
   std::complex<double> tmp_2637;
   std::complex<double> tmp_2638;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2639;
      std::complex<double> tmp_2640;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2641;
         std::complex<double> tmp_2642;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2642 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2641 += tmp_2642;
         tmp_2640 += (Conj(ZU(gI2,j2))) * tmp_2641;
      }
      tmp_2639 += tmp_2640;
      tmp_2638 += (ZU(gI1,j3)) * tmp_2639;
   }
   tmp_2637 += tmp_2638;
   tmp_2590 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2637;
   result += (std::complex<double>(0,-1)) * tmp_2590;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_2643;
   std::complex<double> tmp_2644;
   std::complex<double> tmp_2645;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2645 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2644 += tmp_2645;
   tmp_2643 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2644;
   std::complex<double> tmp_2646;
   std::complex<double> tmp_2647;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2647 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2646 += tmp_2647;
   tmp_2643 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2646;
   std::complex<double> tmp_2648;
   std::complex<double> tmp_2649;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2649 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2648 += tmp_2649;
   tmp_2643 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2648;
   std::complex<double> tmp_2650;
   std::complex<double> tmp_2651;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2651 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2650 += tmp_2651;
   tmp_2643 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2650;
   std::complex<double> tmp_2652;
   std::complex<double> tmp_2653;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2653 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2652 += tmp_2653;
   tmp_2643 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2652;
   std::complex<double> tmp_2654;
   std::complex<double> tmp_2655;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2655 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2654 += tmp_2655;
   tmp_2643 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2654;
   std::complex<double> tmp_2656;
   std::complex<double> tmp_2657;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2657 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2656 += tmp_2657;
   tmp_2643 += (std::complex<double>(0,-1)*Ql*Qs*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2656;
   std::complex<double> tmp_2658;
   std::complex<double> tmp_2659;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2659 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2658 += tmp_2659;
   tmp_2643 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2658;
   std::complex<double> tmp_2660;
   std::complex<double> tmp_2661;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2661 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2660 += tmp_2661;
   tmp_2643 += (std::complex<double>(0,-1)*QHd*Qv*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2660;
   std::complex<double> tmp_2662;
   std::complex<double> tmp_2663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2663 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2662 += tmp_2663;
   tmp_2643 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2662;
   std::complex<double> tmp_2664;
   std::complex<double> tmp_2665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2665 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2664 += tmp_2665;
   tmp_2643 += (std::complex<double>(0,-1)*QHu*Qv*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2664;
   std::complex<double> tmp_2666;
   std::complex<double> tmp_2667;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2667 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2666 += tmp_2667;
   tmp_2643 += (std::complex<double>(0,-1)*Qs*Qv*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(gp)) * tmp_2666;
   std::complex<double> tmp_2668;
   std::complex<double> tmp_2669;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2670;
      std::complex<double> tmp_2671;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2671 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2670 += tmp_2671;
      tmp_2669 += (Conj(ZV(gI2,3 + j2))) * tmp_2670;
   }
   tmp_2668 += tmp_2669;
   tmp_2643 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO2)*
      KroneckerDelta(2,gO1)*Lambdax) * tmp_2668;
   std::complex<double> tmp_2672;
   std::complex<double> tmp_2673;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2674;
      std::complex<double> tmp_2675;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2675 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2674 += tmp_2675;
      tmp_2673 += (Conj(ZV(gI2,3 + j2))) * tmp_2674;
   }
   tmp_2672 += tmp_2673;
   tmp_2643 += (std::complex<double>(0,-0.5)*KroneckerDelta(0,gO1)*
      KroneckerDelta(2,gO2)*Lambdax) * tmp_2672;
   std::complex<double> tmp_2676;
   std::complex<double> tmp_2677;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2678;
      std::complex<double> tmp_2679;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2679 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2678 += tmp_2679;
      tmp_2677 += (ZV(gI1,3 + j2)) * tmp_2678;
   }
   tmp_2676 += tmp_2677;
   tmp_2643 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO2
      )*KroneckerDelta(2,gO1)) * tmp_2676;
   std::complex<double> tmp_2680;
   std::complex<double> tmp_2681;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2682;
      std::complex<double> tmp_2683;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2683 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2682 += tmp_2683;
      tmp_2681 += (ZV(gI1,3 + j2)) * tmp_2682;
   }
   tmp_2680 += tmp_2681;
   tmp_2643 += (std::complex<double>(0,-0.5)*Conj(Lambdax)*KroneckerDelta(0,gO1
      )*KroneckerDelta(2,gO2)) * tmp_2680;
   std::complex<double> tmp_2684;
   std::complex<double> tmp_2685;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2686;
      std::complex<double> tmp_2687;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2688;
         std::complex<double> tmp_2689;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2689 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_2688 += tmp_2689;
         tmp_2687 += (ZV(gI1,3 + j2)) * tmp_2688;
      }
      tmp_2686 += tmp_2687;
      tmp_2685 += (Conj(ZV(gI2,3 + j3))) * tmp_2686;
   }
   tmp_2684 += tmp_2685;
   tmp_2643 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2684;
   std::complex<double> tmp_2690;
   std::complex<double> tmp_2691;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2692;
      std::complex<double> tmp_2693;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2694;
         std::complex<double> tmp_2695;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2695 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_2694 += tmp_2695;
         tmp_2693 += (Conj(ZV(gI2,j2))) * tmp_2694;
      }
      tmp_2692 += tmp_2693;
      tmp_2691 += (ZV(gI1,j3)) * tmp_2692;
   }
   tmp_2690 += tmp_2691;
   tmp_2643 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2690;
   result += (std::complex<double>(0,-1)) * tmp_2643;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2696;
   std::complex<double> tmp_2697;
   std::complex<double> tmp_2698;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2699;
      std::complex<double> tmp_2700;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2700 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2699 += tmp_2700;
      tmp_2698 += (Conj(ZD(gI2,j2))) * tmp_2699;
   }
   tmp_2697 += tmp_2698;
   tmp_2696 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(1,gO2)) * tmp_2697;
   std::complex<double> tmp_2701;
   std::complex<double> tmp_2702;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2703;
      std::complex<double> tmp_2704;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2704 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_2703 += tmp_2704;
      tmp_2702 += (Conj(ZD(gI2,j2))) * tmp_2703;
   }
   tmp_2701 += tmp_2702;
   tmp_2696 += (0.5*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_2701;
   std::complex<double> tmp_2705;
   std::complex<double> tmp_2706;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2707;
      std::complex<double> tmp_2708;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2708 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_2707 += tmp_2708;
      tmp_2706 += (Conj(ZD(gI2,j2))) * tmp_2707;
   }
   tmp_2705 += tmp_2706;
   tmp_2696 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_2705;
   std::complex<double> tmp_2709;
   std::complex<double> tmp_2710;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2711;
      std::complex<double> tmp_2712;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2712 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2711 += tmp_2712;
      tmp_2710 += (ZD(gI1,j2)) * tmp_2711;
   }
   tmp_2709 += tmp_2710;
   tmp_2696 += (-0.5*vS*KroneckerDelta(1,gO2)*Lambdax) * tmp_2709;
   std::complex<double> tmp_2713;
   std::complex<double> tmp_2714;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2715;
      std::complex<double> tmp_2716;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2716 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2715 += tmp_2716;
      tmp_2714 += (ZD(gI1,j2)) * tmp_2715;
   }
   tmp_2713 += tmp_2714;
   tmp_2696 += (-0.5*vu*KroneckerDelta(2,gO2)*Lambdax) * tmp_2713;
   std::complex<double> tmp_2717;
   std::complex<double> tmp_2718;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2719;
      std::complex<double> tmp_2720;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2720 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_2719 += tmp_2720;
      tmp_2718 += (ZD(gI1,j2)) * tmp_2719;
   }
   tmp_2717 += tmp_2718;
   tmp_2696 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_2717;
   result += (std::complex<double>(0,-1)) * tmp_2696;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2721;
   std::complex<double> tmp_2722;
   std::complex<double> tmp_2723;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2724;
      std::complex<double> tmp_2725;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2725 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2724 += tmp_2725;
      tmp_2723 += (Conj(ZE(gI2,j2))) * tmp_2724;
   }
   tmp_2722 += tmp_2723;
   tmp_2721 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(1,gO2)) * tmp_2722;
   std::complex<double> tmp_2726;
   std::complex<double> tmp_2727;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2728;
      std::complex<double> tmp_2729;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2729 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_2728 += tmp_2729;
      tmp_2727 += (Conj(ZE(gI2,j2))) * tmp_2728;
   }
   tmp_2726 += tmp_2727;
   tmp_2721 += (0.5*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_2726;
   std::complex<double> tmp_2730;
   std::complex<double> tmp_2731;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2732;
      std::complex<double> tmp_2733;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2733 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_2732 += tmp_2733;
      tmp_2731 += (Conj(ZE(gI2,j2))) * tmp_2732;
   }
   tmp_2730 += tmp_2731;
   tmp_2721 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_2730;
   std::complex<double> tmp_2734;
   std::complex<double> tmp_2735;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2736;
      std::complex<double> tmp_2737;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2737 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2736 += tmp_2737;
      tmp_2735 += (ZE(gI1,j2)) * tmp_2736;
   }
   tmp_2734 += tmp_2735;
   tmp_2721 += (-0.5*vS*KroneckerDelta(1,gO2)*Lambdax) * tmp_2734;
   std::complex<double> tmp_2738;
   std::complex<double> tmp_2739;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2740;
      std::complex<double> tmp_2741;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2741 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2740 += tmp_2741;
      tmp_2739 += (ZE(gI1,j2)) * tmp_2740;
   }
   tmp_2738 += tmp_2739;
   tmp_2721 += (-0.5*vu*KroneckerDelta(2,gO2)*Lambdax) * tmp_2738;
   std::complex<double> tmp_2742;
   std::complex<double> tmp_2743;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2744;
      std::complex<double> tmp_2745;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2745 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_2744 += tmp_2745;
      tmp_2743 += (ZE(gI1,j2)) * tmp_2744;
   }
   tmp_2742 += tmp_2743;
   tmp_2721 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_2742;
   result += (std::complex<double>(0,-1)) * tmp_2721;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2746;
   std::complex<double> tmp_2747;
   std::complex<double> tmp_2748;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2749;
      std::complex<double> tmp_2750;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2750 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2749 += tmp_2750;
      tmp_2748 += (Conj(ZU(gI2,j2))) * tmp_2749;
   }
   tmp_2747 += tmp_2748;
   tmp_2746 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(0,gO2)) * tmp_2747;
   std::complex<double> tmp_2751;
   std::complex<double> tmp_2752;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2753;
      std::complex<double> tmp_2754;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2754 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2753 += tmp_2754;
      tmp_2752 += (Conj(ZU(gI2,j2))) * tmp_2753;
   }
   tmp_2751 += tmp_2752;
   tmp_2746 += (0.5*vd*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_2751;
   std::complex<double> tmp_2755;
   std::complex<double> tmp_2756;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2757;
      std::complex<double> tmp_2758;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2758 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_2757 += tmp_2758;
      tmp_2756 += (Conj(ZU(gI2,j2))) * tmp_2757;
   }
   tmp_2755 += tmp_2756;
   tmp_2746 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_2755;
   std::complex<double> tmp_2759;
   std::complex<double> tmp_2760;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2761;
      std::complex<double> tmp_2762;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2762 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2761 += tmp_2762;
      tmp_2760 += (ZU(gI1,j2)) * tmp_2761;
   }
   tmp_2759 += tmp_2760;
   tmp_2746 += (-0.5*vS*KroneckerDelta(0,gO2)*Lambdax) * tmp_2759;
   std::complex<double> tmp_2763;
   std::complex<double> tmp_2764;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2765;
      std::complex<double> tmp_2766;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2766 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_2765 += tmp_2766;
      tmp_2764 += (ZU(gI1,j2)) * tmp_2765;
   }
   tmp_2763 += tmp_2764;
   tmp_2746 += (-0.5*vd*KroneckerDelta(2,gO2)*Lambdax) * tmp_2763;
   std::complex<double> tmp_2767;
   std::complex<double> tmp_2768;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2769;
      std::complex<double> tmp_2770;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2770 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_2769 += tmp_2770;
      tmp_2768 += (ZU(gI1,j2)) * tmp_2769;
   }
   tmp_2767 += tmp_2768;
   tmp_2746 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_2767;
   result += (std::complex<double>(0,-1)) * tmp_2746;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2771;
   std::complex<double> tmp_2772;
   std::complex<double> tmp_2773;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2774;
      std::complex<double> tmp_2775;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2775 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2774 += tmp_2775;
      tmp_2773 += (Conj(ZV(gI2,3 + j2))) * tmp_2774;
   }
   tmp_2772 += tmp_2773;
   tmp_2771 += (-0.5*vS*KroneckerDelta(0,gO2)*Lambdax) * tmp_2772;
   std::complex<double> tmp_2776;
   std::complex<double> tmp_2777;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2778;
      std::complex<double> tmp_2779;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2779 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2778 += tmp_2779;
      tmp_2777 += (Conj(ZV(gI2,3 + j2))) * tmp_2778;
   }
   tmp_2776 += tmp_2777;
   tmp_2771 += (-0.5*vd*KroneckerDelta(2,gO2)*Lambdax) * tmp_2776;
   std::complex<double> tmp_2780;
   std::complex<double> tmp_2781;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2782;
      std::complex<double> tmp_2783;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2783 += Conj(TYv(j1,j2))*ZV(gI1,j1);
      }
      tmp_2782 += tmp_2783;
      tmp_2781 += (Conj(ZV(gI2,3 + j2))) * tmp_2782;
   }
   tmp_2780 += tmp_2781;
   tmp_2771 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_2780;
   std::complex<double> tmp_2784;
   std::complex<double> tmp_2785;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2786;
      std::complex<double> tmp_2787;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2787 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2786 += tmp_2787;
      tmp_2785 += (ZV(gI1,3 + j2)) * tmp_2786;
   }
   tmp_2784 += tmp_2785;
   tmp_2771 += (0.5*vS*Conj(Lambdax)*KroneckerDelta(0,gO2)) * tmp_2784;
   std::complex<double> tmp_2788;
   std::complex<double> tmp_2789;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2790;
      std::complex<double> tmp_2791;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2791 += Conj(ZV(gI2,j1))*Yv(j1,j2);
      }
      tmp_2790 += tmp_2791;
      tmp_2789 += (ZV(gI1,3 + j2)) * tmp_2790;
   }
   tmp_2788 += tmp_2789;
   tmp_2771 += (0.5*vd*Conj(Lambdax)*KroneckerDelta(2,gO2)) * tmp_2788;
   std::complex<double> tmp_2792;
   std::complex<double> tmp_2793;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2794;
      std::complex<double> tmp_2795;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2795 += Conj(ZV(gI2,j1))*TYv(j1,j2);
      }
      tmp_2794 += tmp_2795;
      tmp_2793 += (ZV(gI1,3 + j2)) * tmp_2794;
   }
   tmp_2792 += tmp_2793;
   tmp_2771 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_2792;
   result += (std::complex<double>(0,-1)) * tmp_2771;

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

   result = std::complex<double>(0,-0.1)*(10*gp*Qs*Conj(ZH(gI2,2))*
      KroneckerDelta(2,gO2)*Sin(ThetaWp()) + Conj(ZH(gI2,0))*KroneckerDelta(0,gO2)
      *(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp())*
      Sin(ThetaW()) + 10*gp*QHd*Sin(ThetaWp())) - Conj(ZH(gI2,1))*KroneckerDelta(1
      ,gO2)*(5*g2*Cos(ThetaW())*Cos(ThetaWp()) + 3.872983346207417*g1*Cos(ThetaWp(
      ))*Sin(ThetaW()) - 10*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpUAhVZphh(unsigned gO2, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(10*gp*Qs*Conj(ZH(gI2,2))*Cos(ThetaWp(
      ))*KroneckerDelta(2,gO2) + Conj(ZH(gI2,1))*KroneckerDelta(1,gO2)*(10*gp*QHu*
      Cos(ThetaWp()) + 5*g2*Cos(ThetaW())*Sin(ThetaWp()) + 3.872983346207417*g1*
      Sin(ThetaW())*Sin(ThetaWp())) + Conj(ZH(gI2,0))*KroneckerDelta(0,gO2)*(10*gp
      *QHd*Cos(ThetaWp()) - (5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW(
      )))*Sin(ThetaWp())));

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

   result = 0.05*(-(Conj(ZH(gI2,0))*(KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) + 5*
      Sqr(g2) + 20*Sqr(gp)*Sqr(QHd))*ZP(gI1,0) + 5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2
      ))*ZP(gI1,1)) + KroneckerDelta(1,gO2)*(5*vu*(-2*AbsSqr(Lambdax) + Sqr(g2))*
      ZP(gI1,0) + vd*(-3*Sqr(g1) + 5*Sqr(g2) + 20*QHd*QHu*Sqr(gp))*ZP(gI1,1)))) -
      10*Conj(ZH(gI2,2))*(KroneckerDelta(0,gO2)*(2*vS*(AbsSqr(Lambdax) + QHd*Qs*
      Sqr(gp))*ZP(gI1,0) + 1.4142135623730951*Conj(TLambdax)*ZP(gI1,1)) +
      KroneckerDelta(1,gO2)*(1.4142135623730951*TLambdax*ZP(gI1,0) + 2*vS*(AbsSqr(
      Lambdax) + QHu*Qs*Sqr(gp))*ZP(gI1,1))) + Conj(ZH(gI2,1))*(KroneckerDelta(0,
      gO2)*(vu*(3*Sqr(g1) - 5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp)))*ZP(gI1,0) - 5*vd*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)) - KroneckerDelta(1,gO2)*(5*vd*(-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) + vu*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp
      )*Sqr(QHu))*ZP(gI1,1))));

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

std::complex<double> CLASSNAME::CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.05*(-20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*(KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + QHd*Qs*Sqr(gp)) + KroneckerDelta(1,
      gO1)*KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + QHu*Qs*Sqr(gp))) - Conj(ZH(gI1
      ,0))*(5*Conj(ZH(gI2,1))*(KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) +
      KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*AbsSqr(Lambdax) + Sqr(g2))
      + Conj(ZH(gI2,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-3*Sqr(g1) +
      5*Sqr(g2) + 20*QHd*QHu*Sqr(gp)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,
      gO2)*(3*Sqr(g1) + 5*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd)))) + Conj(ZH(gI1,1))*(-5*
      Conj(ZH(gI2,0))*(KroneckerDelta(0,gO2)*KroneckerDelta(1,gO1) +
      KroneckerDelta(0,gO1)*KroneckerDelta(1,gO2))*(-2*AbsSqr(Lambdax) + Sqr(g2))
      + Conj(ZH(gI2,1))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) -
      5*(Sqr(g2) + 4*QHd*QHu*Sqr(gp))) - KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(3*Sqr(g1) + 5*(Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2796;
   std::complex<double> tmp_2797;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2798;
      std::complex<double> tmp_2799;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2799 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2798 += tmp_2799;
      tmp_2797 += (ZUL(gI1,j2)) * tmp_2798;
   }
   tmp_2796 += tmp_2797;
   result += (KroneckerDelta(0,gO2)) * tmp_2796;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2800;
   std::complex<double> tmp_2801;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2802;
      std::complex<double> tmp_2803;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2803 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_2802 += tmp_2803;
      tmp_2801 += (Conj(ZDL(gI2,j2))) * tmp_2802;
   }
   tmp_2800 += tmp_2801;
   result += (KroneckerDelta(1,gO1)) * tmp_2800;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2804;
   std::complex<double> tmp_2805;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2806;
      std::complex<double> tmp_2807;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2807 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2806 += tmp_2807;
      tmp_2805 += (ZVL(gI1,j2)) * tmp_2806;
   }
   tmp_2804 += tmp_2805;
   result += (KroneckerDelta(0,gO2)) * tmp_2804;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2808;
   std::complex<double> tmp_2809;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2810;
      std::complex<double> tmp_2811;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2811 += Conj(ZEL(gI2,j1))*Yv(j1,j2);
      }
      tmp_2810 += tmp_2811;
      tmp_2809 += (Conj(ZVR(gI1,j2))) * tmp_2810;
   }
   tmp_2808 += tmp_2809;
   result += (KroneckerDelta(1,gO1)) * tmp_2808;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_2812;
   std::complex<double> tmp_2813;
   std::complex<double> tmp_2814;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2814 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2813 += tmp_2814;
   tmp_2812 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2813;
   std::complex<double> tmp_2815;
   std::complex<double> tmp_2816;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2816 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2815 += tmp_2816;
   tmp_2812 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2815;
   std::complex<double> tmp_2817;
   std::complex<double> tmp_2818;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2818 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2817 += tmp_2818;
   tmp_2812 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2817;
   std::complex<double> tmp_2819;
   std::complex<double> tmp_2820;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2820 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2819 += tmp_2820;
   tmp_2812 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2819;
   std::complex<double> tmp_2821;
   std::complex<double> tmp_2822;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2822 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2821 += tmp_2822;
   tmp_2812 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2821;
   std::complex<double> tmp_2823;
   std::complex<double> tmp_2824;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2824 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_2823 += tmp_2824;
   tmp_2812 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2823;
   std::complex<double> tmp_2825;
   std::complex<double> tmp_2826;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2826 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2825 += tmp_2826;
   tmp_2812 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2825;
   std::complex<double> tmp_2827;
   std::complex<double> tmp_2828;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2828 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2827 += tmp_2828;
   tmp_2812 += (std::complex<double>(0,-1)*Qd*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2827;
   std::complex<double> tmp_2829;
   std::complex<double> tmp_2830;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2830 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2829 += tmp_2830;
   tmp_2812 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2829;
   std::complex<double> tmp_2831;
   std::complex<double> tmp_2832;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2832 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_2831 += tmp_2832;
   tmp_2812 += (std::complex<double>(0,-1)*Qd*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2831;
   std::complex<double> tmp_2833;
   std::complex<double> tmp_2834;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2835;
      std::complex<double> tmp_2836;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2837;
         std::complex<double> tmp_2838;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2838 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_2837 += tmp_2838;
         tmp_2836 += (ZD(gI1,3 + j2)) * tmp_2837;
      }
      tmp_2835 += tmp_2836;
      tmp_2834 += (Conj(ZD(gI2,3 + j3))) * tmp_2835;
   }
   tmp_2833 += tmp_2834;
   tmp_2812 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2833;
   std::complex<double> tmp_2839;
   std::complex<double> tmp_2840;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2841;
      std::complex<double> tmp_2842;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2843;
         std::complex<double> tmp_2844;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2844 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2843 += tmp_2844;
         tmp_2842 += (Conj(ZD(gI2,j2))) * tmp_2843;
      }
      tmp_2841 += tmp_2842;
      tmp_2840 += (ZD(gI1,j3)) * tmp_2841;
   }
   tmp_2839 += tmp_2840;
   tmp_2812 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2839;
   result += (std::complex<double>(0,-1)) * tmp_2812;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_2845;
   std::complex<double> tmp_2846;
   std::complex<double> tmp_2847;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2847 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2846 += tmp_2847;
   tmp_2845 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2846;
   std::complex<double> tmp_2848;
   std::complex<double> tmp_2849;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2849 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2848 += tmp_2849;
   tmp_2845 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2848;
   std::complex<double> tmp_2850;
   std::complex<double> tmp_2851;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2851 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2850 += tmp_2851;
   tmp_2845 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2850;
   std::complex<double> tmp_2852;
   std::complex<double> tmp_2853;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2853 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2852 += tmp_2853;
   tmp_2845 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2852;
   std::complex<double> tmp_2854;
   std::complex<double> tmp_2855;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2855 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2854 += tmp_2855;
   tmp_2845 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2854;
   std::complex<double> tmp_2856;
   std::complex<double> tmp_2857;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2857 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_2856 += tmp_2857;
   tmp_2845 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2856;
   std::complex<double> tmp_2858;
   std::complex<double> tmp_2859;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2859 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2858 += tmp_2859;
   tmp_2845 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2858;
   std::complex<double> tmp_2860;
   std::complex<double> tmp_2861;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2861 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2860 += tmp_2861;
   tmp_2845 += (std::complex<double>(0,-1)*Qe*QHd*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2860;
   std::complex<double> tmp_2862;
   std::complex<double> tmp_2863;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2863 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2862 += tmp_2863;
   tmp_2845 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2862;
   std::complex<double> tmp_2864;
   std::complex<double> tmp_2865;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2865 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_2864 += tmp_2865;
   tmp_2845 += (std::complex<double>(0,-1)*Qe*QHu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2864;
   std::complex<double> tmp_2866;
   std::complex<double> tmp_2867;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2868;
      std::complex<double> tmp_2869;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2870;
         std::complex<double> tmp_2871;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2871 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_2870 += tmp_2871;
         tmp_2869 += (ZE(gI1,3 + j2)) * tmp_2870;
      }
      tmp_2868 += tmp_2869;
      tmp_2867 += (Conj(ZE(gI2,3 + j3))) * tmp_2868;
   }
   tmp_2866 += tmp_2867;
   tmp_2845 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2866;
   std::complex<double> tmp_2872;
   std::complex<double> tmp_2873;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2874;
      std::complex<double> tmp_2875;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2876;
         std::complex<double> tmp_2877;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2877 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_2876 += tmp_2877;
         tmp_2875 += (Conj(ZE(gI2,j2))) * tmp_2876;
      }
      tmp_2874 += tmp_2875;
      tmp_2873 += (ZE(gI1,j3)) * tmp_2874;
   }
   tmp_2872 += tmp_2873;
   tmp_2845 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2872;
   result += (std::complex<double>(0,-1)) * tmp_2845;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Qq = LOCALINPUT(Qq);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_2878;
   std::complex<double> tmp_2879;
   std::complex<double> tmp_2880;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2880 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2879 += tmp_2880;
   tmp_2878 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2879;
   std::complex<double> tmp_2881;
   std::complex<double> tmp_2882;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2882 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2881 += tmp_2882;
   tmp_2878 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2881;
   std::complex<double> tmp_2883;
   std::complex<double> tmp_2884;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2884 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2883 += tmp_2884;
   tmp_2878 += (std::complex<double>(0,-1)*QHd*Qq*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2883;
   std::complex<double> tmp_2885;
   std::complex<double> tmp_2886;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2886 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2885 += tmp_2886;
   tmp_2878 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2885;
   std::complex<double> tmp_2887;
   std::complex<double> tmp_2888;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2888 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2887 += tmp_2888;
   tmp_2878 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2887;
   std::complex<double> tmp_2889;
   std::complex<double> tmp_2890;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2890 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2889 += tmp_2890;
   tmp_2878 += (std::complex<double>(0,-1)*QHu*Qq*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2889;
   std::complex<double> tmp_2891;
   std::complex<double> tmp_2892;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2892 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2891 += tmp_2892;
   tmp_2878 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2891;
   std::complex<double> tmp_2893;
   std::complex<double> tmp_2894;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2894 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2893 += tmp_2894;
   tmp_2878 += (std::complex<double>(0,-1)*QHd*Qu*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2893;
   std::complex<double> tmp_2895;
   std::complex<double> tmp_2896;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2896 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2895 += tmp_2896;
   tmp_2878 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2895;
   std::complex<double> tmp_2897;
   std::complex<double> tmp_2898;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2898 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_2897 += tmp_2898;
   tmp_2878 += (std::complex<double>(0,-1)*QHu*Qu*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2897;
   std::complex<double> tmp_2899;
   std::complex<double> tmp_2900;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2901;
      std::complex<double> tmp_2902;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2903;
         std::complex<double> tmp_2904;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2904 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_2903 += tmp_2904;
         tmp_2902 += (ZU(gI1,3 + j2)) * tmp_2903;
      }
      tmp_2901 += tmp_2902;
      tmp_2900 += (Conj(ZU(gI2,3 + j3))) * tmp_2901;
   }
   tmp_2899 += tmp_2900;
   tmp_2878 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2899;
   std::complex<double> tmp_2905;
   std::complex<double> tmp_2906;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2907;
      std::complex<double> tmp_2908;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2909;
         std::complex<double> tmp_2910;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2910 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2909 += tmp_2910;
         tmp_2908 += (Conj(ZU(gI2,j2))) * tmp_2909;
      }
      tmp_2907 += tmp_2908;
      tmp_2906 += (ZU(gI1,j3)) * tmp_2907;
   }
   tmp_2905 += tmp_2906;
   tmp_2878 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2905;
   result += (std::complex<double>(0,-1)) * tmp_2878;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto Ql = LOCALINPUT(Ql);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_2911;
   std::complex<double> tmp_2912;
   std::complex<double> tmp_2913;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2913 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2912 += tmp_2913;
   tmp_2911 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2912;
   std::complex<double> tmp_2914;
   std::complex<double> tmp_2915;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2915 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2914 += tmp_2915;
   tmp_2911 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_2914;
   std::complex<double> tmp_2916;
   std::complex<double> tmp_2917;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2917 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2916 += tmp_2917;
   tmp_2911 += (std::complex<double>(0,-1)*QHd*Ql*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2916;
   std::complex<double> tmp_2918;
   std::complex<double> tmp_2919;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2919 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2918 += tmp_2919;
   tmp_2911 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2918;
   std::complex<double> tmp_2920;
   std::complex<double> tmp_2921;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2921 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2920 += tmp_2921;
   tmp_2911 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_2920;
   std::complex<double> tmp_2922;
   std::complex<double> tmp_2923;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2923 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2922 += tmp_2923;
   tmp_2911 += (std::complex<double>(0,-1)*QHu*Ql*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2922;
   std::complex<double> tmp_2924;
   std::complex<double> tmp_2925;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2925 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2924 += tmp_2925;
   tmp_2911 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_2924;
   std::complex<double> tmp_2926;
   std::complex<double> tmp_2927;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2927 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2926 += tmp_2927;
   tmp_2911 += (std::complex<double>(0,-1)*QHd*Qv*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(gp)) * tmp_2926;
   std::complex<double> tmp_2928;
   std::complex<double> tmp_2929;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2929 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2928 += tmp_2929;
   tmp_2911 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_2928;
   std::complex<double> tmp_2930;
   std::complex<double> tmp_2931;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2931 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_2930 += tmp_2931;
   tmp_2911 += (std::complex<double>(0,-1)*QHu*Qv*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(gp)) * tmp_2930;
   std::complex<double> tmp_2932;
   std::complex<double> tmp_2933;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2934;
      std::complex<double> tmp_2935;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2936;
         std::complex<double> tmp_2937;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2937 += Conj(Yv(j1,j3))*Yv(j1,j2);
         }
         tmp_2936 += tmp_2937;
         tmp_2935 += (ZV(gI1,3 + j2)) * tmp_2936;
      }
      tmp_2934 += tmp_2935;
      tmp_2933 += (Conj(ZV(gI2,3 + j3))) * tmp_2934;
   }
   tmp_2932 += tmp_2933;
   tmp_2911 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_2932;
   std::complex<double> tmp_2938;
   std::complex<double> tmp_2939;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2940;
      std::complex<double> tmp_2941;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2942;
         std::complex<double> tmp_2943;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2943 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_2942 += tmp_2943;
         tmp_2941 += (Conj(ZV(gI2,j2))) * tmp_2942;
      }
      tmp_2940 += tmp_2941;
      tmp_2939 += (ZV(gI1,j3)) * tmp_2940;
   }
   tmp_2938 += tmp_2939;
   tmp_2911 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_2938;
   result += (std::complex<double>(0,-1)) * tmp_2911;

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

   std::complex<double> tmp_2944;
   std::complex<double> tmp_2945;
   std::complex<double> tmp_2946;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2946 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2945 += tmp_2946;
   tmp_2944 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_2945;
   std::complex<double> tmp_2947;
   std::complex<double> tmp_2948;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2948 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_2947 += tmp_2948;
   tmp_2944 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_2947;
   std::complex<double> tmp_2949;
   std::complex<double> tmp_2950;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2951;
      std::complex<double> tmp_2952;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2952 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_2951 += tmp_2952;
      tmp_2950 += (Conj(ZD(gI2,j2))) * tmp_2951;
   }
   tmp_2949 += tmp_2950;
   tmp_2944 += (std::complex<double>(0.,0.7071067811865475)*vS*Conj(Lambdax)*
      KroneckerDelta(0,gO2)) * tmp_2949;
   std::complex<double> tmp_2953;
   std::complex<double> tmp_2954;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2955;
      std::complex<double> tmp_2956;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2956 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_2955 += tmp_2956;
      tmp_2954 += (Conj(ZD(gI2,j2))) * tmp_2955;
   }
   tmp_2953 += tmp_2954;
   tmp_2944 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_2953;
   std::complex<double> tmp_2957;
   std::complex<double> tmp_2958;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2959;
      std::complex<double> tmp_2960;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2960 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2959 += tmp_2960;
      tmp_2958 += (ZU(gI1,j2)) * tmp_2959;
   }
   tmp_2957 += tmp_2958;
   tmp_2944 += (std::complex<double>(0.,0.7071067811865475)*vS*KroneckerDelta(1
      ,gO2)*Lambdax) * tmp_2957;
   std::complex<double> tmp_2961;
   std::complex<double> tmp_2962;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2963;
      std::complex<double> tmp_2964;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2964 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_2963 += tmp_2964;
      tmp_2962 += (ZU(gI1,j2)) * tmp_2963;
   }
   tmp_2961 += tmp_2962;
   tmp_2944 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_2961;
   std::complex<double> tmp_2965;
   std::complex<double> tmp_2966;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2967;
      std::complex<double> tmp_2968;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2969;
         std::complex<double> tmp_2970;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2970 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_2969 += tmp_2970;
         tmp_2968 += (ZU(gI1,3 + j2)) * tmp_2969;
      }
      tmp_2967 += tmp_2968;
      tmp_2966 += (Conj(ZD(gI2,3 + j3))) * tmp_2967;
   }
   tmp_2965 += tmp_2966;
   tmp_2944 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(0
      ,gO2)) * tmp_2965;
   std::complex<double> tmp_2971;
   std::complex<double> tmp_2972;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2973;
      std::complex<double> tmp_2974;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2975;
         std::complex<double> tmp_2976;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2976 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_2975 += tmp_2976;
         tmp_2974 += (ZU(gI1,3 + j2)) * tmp_2975;
      }
      tmp_2973 += tmp_2974;
      tmp_2972 += (Conj(ZD(gI2,3 + j3))) * tmp_2973;
   }
   tmp_2971 += tmp_2972;
   tmp_2944 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(1
      ,gO2)) * tmp_2971;
   std::complex<double> tmp_2977;
   std::complex<double> tmp_2978;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2979;
      std::complex<double> tmp_2980;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2981;
         std::complex<double> tmp_2982;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2982 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_2981 += tmp_2982;
         tmp_2980 += (Conj(ZD(gI2,j2))) * tmp_2981;
      }
      tmp_2979 += tmp_2980;
      tmp_2978 += (ZU(gI1,j3)) * tmp_2979;
   }
   tmp_2977 += tmp_2978;
   tmp_2944 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_2977;
   std::complex<double> tmp_2983;
   std::complex<double> tmp_2984;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_2985;
      std::complex<double> tmp_2986;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_2987;
         std::complex<double> tmp_2988;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_2988 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_2987 += tmp_2988;
         tmp_2986 += (Conj(ZD(gI2,j2))) * tmp_2987;
      }
      tmp_2985 += tmp_2986;
      tmp_2984 += (ZU(gI1,j3)) * tmp_2985;
   }
   tmp_2983 += tmp_2984;
   tmp_2944 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(1
      ,gO2)) * tmp_2983;
   result += (std::complex<double>(0,-1)) * tmp_2944;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2989;
   std::complex<double> tmp_2990;
   std::complex<double> tmp_2991;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2991 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2990 += tmp_2991;
   tmp_2989 += (std::complex<double>(0.,-0.35355339059327373)*vd*KroneckerDelta
      (0,gO2)*Sqr(g2)) * tmp_2990;
   std::complex<double> tmp_2992;
   std::complex<double> tmp_2993;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2993 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_2992 += tmp_2993;
   tmp_2989 += (std::complex<double>(0.,-0.35355339059327373)*vu*KroneckerDelta
      (1,gO2)*Sqr(g2)) * tmp_2992;
   std::complex<double> tmp_2994;
   std::complex<double> tmp_2995;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2996;
      std::complex<double> tmp_2997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2997 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_2996 += tmp_2997;
      tmp_2995 += (ZV(gI1,j2)) * tmp_2996;
   }
   tmp_2994 += tmp_2995;
   tmp_2989 += (std::complex<double>(0.,0.7071067811865475)*vS*KroneckerDelta(1
      ,gO2)*Lambdax) * tmp_2994;
   std::complex<double> tmp_2998;
   std::complex<double> tmp_2999;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3000;
      std::complex<double> tmp_3001;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3001 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_3000 += tmp_3001;
      tmp_2999 += (ZV(gI1,j2)) * tmp_3000;
   }
   tmp_2998 += tmp_2999;
   tmp_2989 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_2998;
   std::complex<double> tmp_3002;
   std::complex<double> tmp_3003;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3004;
      std::complex<double> tmp_3005;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3005 += Conj(ZE(gI2,j1))*Yv(j1,j2);
      }
      tmp_3004 += tmp_3005;
      tmp_3003 += (ZV(gI1,3 + j2)) * tmp_3004;
   }
   tmp_3002 += tmp_3003;
   tmp_2989 += (std::complex<double>(0.,0.7071067811865475)*vS*Conj(Lambdax)*
      KroneckerDelta(0,gO2)) * tmp_3002;
   std::complex<double> tmp_3006;
   std::complex<double> tmp_3007;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3008;
      std::complex<double> tmp_3009;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3009 += Conj(ZE(gI2,j1))*TYv(j1,j2);
      }
      tmp_3008 += tmp_3009;
      tmp_3007 += (ZV(gI1,3 + j2)) * tmp_3008;
   }
   tmp_3006 += tmp_3007;
   tmp_2989 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_3006;
   std::complex<double> tmp_3010;
   std::complex<double> tmp_3011;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3012;
      std::complex<double> tmp_3013;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3014;
         std::complex<double> tmp_3015;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3015 += Conj(Ye(j3,j1))*Yv(j1,j2);
         }
         tmp_3014 += tmp_3015;
         tmp_3013 += (ZV(gI1,3 + j2)) * tmp_3014;
      }
      tmp_3012 += tmp_3013;
      tmp_3011 += (Conj(ZE(gI2,3 + j3))) * tmp_3012;
   }
   tmp_3010 += tmp_3011;
   tmp_2989 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(0
      ,gO2)) * tmp_3010;
   std::complex<double> tmp_3016;
   std::complex<double> tmp_3017;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3018;
      std::complex<double> tmp_3019;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3020;
         std::complex<double> tmp_3021;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3021 += Conj(Ye(j3,j1))*Yv(j1,j2);
         }
         tmp_3020 += tmp_3021;
         tmp_3019 += (ZV(gI1,3 + j2)) * tmp_3020;
      }
      tmp_3018 += tmp_3019;
      tmp_3017 += (Conj(ZE(gI2,3 + j3))) * tmp_3018;
   }
   tmp_3016 += tmp_3017;
   tmp_2989 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(1
      ,gO2)) * tmp_3016;
   std::complex<double> tmp_3022;
   std::complex<double> tmp_3023;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3024;
      std::complex<double> tmp_3025;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3026;
         std::complex<double> tmp_3027;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3027 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_3026 += tmp_3027;
         tmp_3025 += (Conj(ZE(gI2,j2))) * tmp_3026;
      }
      tmp_3024 += tmp_3025;
      tmp_3023 += (ZV(gI1,j3)) * tmp_3024;
   }
   tmp_3022 += tmp_3023;
   tmp_2989 += (std::complex<double>(0.,0.7071067811865475)*vd*KroneckerDelta(0
      ,gO2)) * tmp_3022;
   std::complex<double> tmp_3028;
   std::complex<double> tmp_3029;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_3030;
      std::complex<double> tmp_3031;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_3032;
         std::complex<double> tmp_3033;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_3033 += Conj(Yv(j3,j1))*Yv(j2,j1);
         }
         tmp_3032 += tmp_3033;
         tmp_3031 += (Conj(ZE(gI2,j2))) * tmp_3032;
      }
      tmp_3030 += tmp_3031;
      tmp_3029 += (ZV(gI1,j3)) * tmp_3030;
   }
   tmp_3028 += tmp_3029;
   tmp_2989 += (std::complex<double>(0.,0.7071067811865475)*vu*KroneckerDelta(1
      ,gO2)) * tmp_3028;
   result += (std::complex<double>(0,-1)) * tmp_2989;

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

   result = 0.5*g2*(Conj(ZH(gI2,0))*KroneckerDelta(0,gO2) - Conj(ZH(gI2,1))*
      KroneckerDelta(1,gO2));

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

std::complex<double> CLASSNAME::CpVZVZhhhh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(
      ThetaWp())) + Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*(20*g2*gp*QHd*Cos(ThetaW())*
      Cos(ThetaWp())*Sin(ThetaWp()) + 15.491933384829668*g1*gp*QHd*Cos(ThetaWp())*
      Sin(ThetaW())*Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(
      ThetaW()) + 3*g1*Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW()))*Sqr(Cos(ThetaWp())) + 20*Sqr(gp)*Sqr(QHd)*Sqr(Sin(ThetaWp()))) +
      Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*(-20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp())*
      Sin(ThetaWp()) - 15.491933384829668*g1*gp*QHu*Cos(ThetaWp())*Sin(ThetaW())*
      Sin(ThetaWp()) + g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()) + 3*g1
      *Sin(ThetaW()))*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp())) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Sin(ThetaWp()))));

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh(unsigned gI1, unsigned gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(2*gp*Qs*Conj(ZA(gI2,2))*Conj(ZH(gI1,2
      ))*Sin(ThetaWp()) + Conj(ZA(gI2,0))*Conj(ZH(gI1,0))*(g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*
      Sin(ThetaWp())) - Conj(ZA(gI2,1))*Conj(ZH(gI1,1))*(g2*Cos(ThetaW())*Cos(
      ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 2*gp*QHu*
      Sin(ThetaWp())));

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

double CLASSNAME::CpVZbarFvFvPR(unsigned gI1, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   double result = 0.0;

   result = -(KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos(ThetaWp())*Sin
      (ThetaW()) - gp*Qv*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_3034;
   std::complex<double> tmp_3035;
   std::complex<double> tmp_3036;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3036 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3035 += tmp_3036;
   tmp_3034 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_3035;
   std::complex<double> tmp_3037;
   std::complex<double> tmp_3038;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3038 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3037 += tmp_3038;
   tmp_3034 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_3037;
   std::complex<double> tmp_3039;
   std::complex<double> tmp_3040;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3040 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3039 += tmp_3040;
   tmp_3034 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_3039;
   std::complex<double> tmp_3041;
   std::complex<double> tmp_3042;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3042 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3041 += tmp_3042;
   tmp_3034 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_3041;
   std::complex<double> tmp_3043;
   std::complex<double> tmp_3044;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3044 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3043 += tmp_3044;
   tmp_3034 += (std::complex<double>(0,-2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_3043;
   std::complex<double> tmp_3045;
   std::complex<double> tmp_3046;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3046 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3045 += tmp_3046;
   tmp_3034 += (std::complex<double>(0.,-0.5163977794943222)*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3045;
   std::complex<double> tmp_3047;
   std::complex<double> tmp_3048;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3048 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3047 += tmp_3048;
   tmp_3034 += (std::complex<double>(0.,-1.0327955589886444)*g1*gp*Qd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3047;
   std::complex<double> tmp_3049;
   std::complex<double> tmp_3050;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3050 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3049 += tmp_3050;
   tmp_3034 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())))
      * tmp_3049;
   std::complex<double> tmp_3051;
   std::complex<double> tmp_3052;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3052 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3051 += tmp_3052;
   tmp_3034 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qd)*Sqr(Sin(ThetaWp())))
      * tmp_3051;
   result += (std::complex<double>(0,-1)) * tmp_3034;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_3053;
   std::complex<double> tmp_3054;
   std::complex<double> tmp_3055;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3055 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3054 += tmp_3055;
   tmp_3053 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_3054;
   std::complex<double> tmp_3056;
   std::complex<double> tmp_3057;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3057 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3056 += tmp_3057;
   tmp_3053 += (std::complex<double>(0.,-0.7745966692414834)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_3056;
   std::complex<double> tmp_3058;
   std::complex<double> tmp_3059;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3059 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3058 += tmp_3059;
   tmp_3053 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin
      (ThetaW()))) * tmp_3058;
   std::complex<double> tmp_3060;
   std::complex<double> tmp_3061;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3061 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3060 += tmp_3061;
   tmp_3053 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin
      (ThetaW()))) * tmp_3060;
   std::complex<double> tmp_3062;
   std::complex<double> tmp_3063;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3063 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3062 += tmp_3063;
   tmp_3053 += (std::complex<double>(0,-2)*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_3062;
   std::complex<double> tmp_3064;
   std::complex<double> tmp_3065;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3065 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3064 += tmp_3065;
   tmp_3053 += (std::complex<double>(0.,1.5491933384829668)*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3064;
   std::complex<double> tmp_3066;
   std::complex<double> tmp_3067;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3067 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3066 += tmp_3067;
   tmp_3053 += (std::complex<double>(0.,-3.0983866769659336)*g1*gp*Qe*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3066;
   std::complex<double> tmp_3068;
   std::complex<double> tmp_3069;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3069 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3068 += tmp_3069;
   tmp_3053 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())))
      * tmp_3068;
   std::complex<double> tmp_3070;
   std::complex<double> tmp_3071;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3071 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3070 += tmp_3071;
   tmp_3053 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qe)*Sqr(Sin(ThetaWp())))
      * tmp_3070;
   result += (std::complex<double>(0,-1)) * tmp_3053;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_3072;
   std::complex<double> tmp_3073;
   std::complex<double> tmp_3074;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3074 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3073 += tmp_3074;
   tmp_3072 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_3073;
   std::complex<double> tmp_3075;
   std::complex<double> tmp_3076;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3076 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3075 += tmp_3076;
   tmp_3072 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_3075;
   std::complex<double> tmp_3077;
   std::complex<double> tmp_3078;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3078 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3077 += tmp_3078;
   tmp_3072 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_3077;
   std::complex<double> tmp_3079;
   std::complex<double> tmp_3080;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3080 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3079 += tmp_3080;
   tmp_3072 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Cos(
      ThetaWp()))*Sqr(Sin(ThetaW()))) * tmp_3079;
   std::complex<double> tmp_3081;
   std::complex<double> tmp_3082;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3082 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3081 += tmp_3082;
   tmp_3072 += (std::complex<double>(0,2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_3081;
   std::complex<double> tmp_3083;
   std::complex<double> tmp_3084;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3084 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3083 += tmp_3084;
   tmp_3072 += (std::complex<double>(0.,-0.5163977794943222)*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3083;
   std::complex<double> tmp_3085;
   std::complex<double> tmp_3086;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3086 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3085 += tmp_3086;
   tmp_3072 += (std::complex<double>(0.,2.065591117977289)*g1*gp*Qu*Cos(ThetaWp
      ())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3085;
   std::complex<double> tmp_3087;
   std::complex<double> tmp_3088;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3088 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3087 += tmp_3088;
   tmp_3072 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Sin(ThetaWp())))
      * tmp_3087;
   std::complex<double> tmp_3089;
   std::complex<double> tmp_3090;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3090 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3089 += tmp_3090;
   tmp_3072 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qu)*Sqr(Sin(ThetaWp())))
      * tmp_3089;
   result += (std::complex<double>(0,-1)) * tmp_3072;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_3091;
   std::complex<double> tmp_3092;
   std::complex<double> tmp_3093;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3093 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3092 += tmp_3093;
   tmp_3091 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Cos(
      ThetaWp()))) * tmp_3092;
   std::complex<double> tmp_3094;
   std::complex<double> tmp_3095;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3095 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3094 += tmp_3095;
   tmp_3091 += (std::complex<double>(0.,0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Cos(ThetaWp()))) * tmp_3094;
   std::complex<double> tmp_3096;
   std::complex<double> tmp_3097;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3097 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3096 += tmp_3097;
   tmp_3091 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin
      (ThetaW()))) * tmp_3096;
   std::complex<double> tmp_3098;
   std::complex<double> tmp_3099;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3099 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3098 += tmp_3099;
   tmp_3091 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Cos(ThetaWp()))*Sqr(Sin
      (ThetaW()))) * tmp_3098;
   std::complex<double> tmp_3100;
   std::complex<double> tmp_3101;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3101 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3100 += tmp_3101;
   tmp_3091 += (std::complex<double>(0,2)*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_3100;
   std::complex<double> tmp_3102;
   std::complex<double> tmp_3103;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3103 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3102 += tmp_3103;
   tmp_3091 += (std::complex<double>(0.,1.5491933384829668)*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3102;
   std::complex<double> tmp_3104;
   std::complex<double> tmp_3105;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3105 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3104 += tmp_3105;
   tmp_3091 += (std::complex<double>(0.,-3.0983866769659336)*g1*gp*Qv*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3104;
   std::complex<double> tmp_3106;
   std::complex<double> tmp_3107;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3107 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3106 += tmp_3107;
   tmp_3091 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Ql)*Sqr(Sin(ThetaWp())))
      * tmp_3106;
   std::complex<double> tmp_3108;
   std::complex<double> tmp_3109;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3109 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3108 += tmp_3109;
   tmp_3091 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qv)*Sqr(Sin(ThetaWp())))
      * tmp_3108;
   result += (std::complex<double>(0,-1)) * tmp_3091;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_3110;
   std::complex<double> tmp_3111;
   std::complex<double> tmp_3112;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3112 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3111 += tmp_3112;
   tmp_3110 += (2*(0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - 3*gp*Qd
      *Sin(ThetaWp()))) * tmp_3111;
   std::complex<double> tmp_3113;
   std::complex<double> tmp_3114;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3114 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3113 += tmp_3114;
   tmp_3110 += (-3*g2*Cos(ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 6*gp*Qq*Sin(ThetaWp())) * tmp_3113;
   result += (0.16666666666666666) * tmp_3110;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_3115;
   std::complex<double> tmp_3116;
   std::complex<double> tmp_3117;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3117 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3116 += tmp_3117;
   tmp_3115 += (-2*(-0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) + gp*Qe
      *Sin(ThetaWp()))) * tmp_3116;
   std::complex<double> tmp_3118;
   std::complex<double> tmp_3119;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3119 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3118 += tmp_3119;
   tmp_3115 += (-(g2*Cos(ThetaW())*Cos(ThetaWp())) + 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*Sin(ThetaWp())) * tmp_3118;
   result += (0.5) * tmp_3115;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_3120;
   std::complex<double> tmp_3121;
   std::complex<double> tmp_3122;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3122 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3121 += tmp_3122;
   tmp_3120 += (3*g2*Cos(ThetaW())*Cos(ThetaWp()) - 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 6*gp*Qq*Sin(ThetaWp())) * tmp_3121;
   std::complex<double> tmp_3123;
   std::complex<double> tmp_3124;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3124 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3123 += tmp_3124;
   tmp_3120 += (-2*(1.5491933384829668*g1*Cos(ThetaWp())*Sin(ThetaW()) + 3*gp*
      Qu*Sin(ThetaWp()))) * tmp_3123;
   result += (0.16666666666666666) * tmp_3120;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_3125;
   std::complex<double> tmp_3126;
   std::complex<double> tmp_3127;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3127 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3126 += tmp_3127;
   tmp_3125 += (g2*Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(
      ThetaWp())*Sin(ThetaW()) + 2*gp*Ql*Sin(ThetaWp())) * tmp_3126;
   std::complex<double> tmp_3128;
   std::complex<double> tmp_3129;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3129 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3128 += tmp_3129;
   tmp_3125 += (2*(0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW()) - gp*Qv*
      Sin(ThetaWp()))) * tmp_3128;
   result += (0.5) * tmp_3125;

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
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*(4*vS*Conj(ZH(gI2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Sin(ThetaWp())) + vd*
      Conj(ZH(gI2,0))*Sqr(g2*Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*
      Cos(ThetaWp())*Sin(ThetaW()) + 2*gp*QHd*Sin(ThetaWp())) + vu*Conj(ZH(gI2,1))
      *Sqr(g2*Cos(ThetaW())*Cos(ThetaWp()) + 0.7745966692414834*g1*Cos(ThetaWp())*
      Sin(ThetaW()) - 2*gp*QHu*Sin(ThetaWp())));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZphh(unsigned gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(20*vS*Conj(ZH(gI2,2))*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(gp)*
      Sqr(Qs) - vd*Conj(ZH(gI2,0))*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(
      Cos(ThetaW())) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Cos(ThetaWp()
      )) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin
      (ThetaW()))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Sin(ThetaWp()))
      + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin
      (ThetaWp()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(ThetaWp()))))
      - vu*Conj(ZH(gI2,1))*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW
      ())) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(
      ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()
      ))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*
      Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp
      ()) + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(Sin(ThetaWp())))));

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
      + Conj(ZA(gI1,1))*Conj(ZA(gI2,1))*(20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(
      Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*(7.745966692414834*gp*QHu*Sin(ThetaW
      ())*Sin(2*ThetaWp()) + 3.872983346207417*g2*Sin(2*ThetaW())*Sqr(Sin(ThetaWp(
      ))) + 3*g1*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZphhhh(unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(20*Conj(ZH(gI1,2))*Conj(ZH(gI2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(
      ThetaWp())) + Conj(ZH(gI1,0))*Conj(ZH(gI2,0))*(-2*gp*QHd*(5*g2*Cos(ThetaW())
      + 3.872983346207417*g1*Sin(ThetaW()))*Sin(2*ThetaWp()) + 20*Sqr(gp)*Sqr(QHd
      )*Sqr(Cos(ThetaWp())) + (g1*Sin(ThetaW())*(7.745966692414834*g2*Cos(ThetaW()
      ) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())))*Sqr(Sin(ThetaWp())))
      + Conj(ZH(gI1,1))*Conj(ZH(gI2,1))*(20*g2*gp*QHu*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp()) + 20*Sqr(gp)*Sqr(QHu)*Sqr(Cos(ThetaWp())) + 5*Sqr(g2)*Sqr(
      Cos(ThetaW()))*Sqr(Sin(ThetaWp())) + g1*(7.745966692414834*gp*QHu*Sin(ThetaW
      ())*Sin(2*ThetaWp()) + 3.872983346207417*g2*Sin(2*ThetaW())*Sqr(Sin(ThetaWp(
      ))) + 3*g1*Sqr(Sin(ThetaW()))*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpVZphhAh(unsigned gI1, unsigned gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHu = LOCALINPUT(QHu);
   const auto QHd = LOCALINPUT(QHd);

   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(2*gp*Qs*Conj(ZA(gI2,2))*Conj(ZH(gI1,2
      ))*Cos(ThetaWp()) + Conj(ZA(gI2,1))*Conj(ZH(gI1,1))*(2*gp*QHu*Cos(ThetaWp())
      + g2*Cos(ThetaW())*Sin(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*Sin
      (ThetaWp())) + Conj(ZA(gI2,0))*Conj(ZH(gI1,0))*(2*gp*QHd*Cos(ThetaWp()) - (
      g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp())));

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

double CLASSNAME::CpVZpbarFvFvPR(unsigned gI1, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   double result = 0.0;

   result = KroneckerDelta(gI1,gI2)*(gp*Qv*Cos(ThetaWp()) + 0.7745966692414834*
      g1*Sin(ThetaW())*Sin(ThetaWp()));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_3130;
   std::complex<double> tmp_3131;
   std::complex<double> tmp_3132;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3132 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3131 += tmp_3132;
   tmp_3130 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())))
      * tmp_3131;
   std::complex<double> tmp_3133;
   std::complex<double> tmp_3134;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3134 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3133 += tmp_3134;
   tmp_3130 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qd)*Sqr(Cos(ThetaWp())))
      * tmp_3133;
   std::complex<double> tmp_3135;
   std::complex<double> tmp_3136;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3136 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3135 += tmp_3136;
   tmp_3130 += (std::complex<double>(0,2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_3135;
   std::complex<double> tmp_3137;
   std::complex<double> tmp_3138;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3138 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3137 += tmp_3138;
   tmp_3130 += (std::complex<double>(0.,0.5163977794943222)*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3137;
   std::complex<double> tmp_3139;
   std::complex<double> tmp_3140;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3140 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3139 += tmp_3140;
   tmp_3130 += (std::complex<double>(0.,1.0327955589886444)*g1*gp*Qd*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3139;
   std::complex<double> tmp_3141;
   std::complex<double> tmp_3142;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3142 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3141 += tmp_3142;
   tmp_3130 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3141;
   std::complex<double> tmp_3143;
   std::complex<double> tmp_3144;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3144 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3143 += tmp_3144;
   tmp_3130 += (std::complex<double>(0.,0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_3143;
   std::complex<double> tmp_3145;
   std::complex<double> tmp_3146;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3146 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3145 += tmp_3146;
   tmp_3130 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_3145;
   std::complex<double> tmp_3147;
   std::complex<double> tmp_3148;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3148 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3147 += tmp_3148;
   tmp_3130 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_3147;
   result += (std::complex<double>(0,-1)) * tmp_3130;

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_3149;
   std::complex<double> tmp_3150;
   std::complex<double> tmp_3151;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3151 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3150 += tmp_3151;
   tmp_3149 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Ql)*Sqr(Cos(ThetaWp())))
      * tmp_3150;
   std::complex<double> tmp_3152;
   std::complex<double> tmp_3153;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3153 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3152 += tmp_3153;
   tmp_3149 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qe)*Sqr(Cos(ThetaWp())))
      * tmp_3152;
   std::complex<double> tmp_3154;
   std::complex<double> tmp_3155;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3155 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3154 += tmp_3155;
   tmp_3149 += (std::complex<double>(0,2)*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp())
      *Sin(ThetaWp())) * tmp_3154;
   std::complex<double> tmp_3156;
   std::complex<double> tmp_3157;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3157 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3156 += tmp_3157;
   tmp_3149 += (std::complex<double>(0.,-1.5491933384829668)*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3156;
   std::complex<double> tmp_3158;
   std::complex<double> tmp_3159;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3159 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3158 += tmp_3159;
   tmp_3149 += (std::complex<double>(0.,3.0983866769659336)*g1*gp*Qe*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3158;
   std::complex<double> tmp_3160;
   std::complex<double> tmp_3161;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3161 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3160 += tmp_3161;
   tmp_3149 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3160;
   std::complex<double> tmp_3162;
   std::complex<double> tmp_3163;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3163 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3162 += tmp_3163;
   tmp_3149 += (std::complex<double>(0.,-0.7745966692414834)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_3162;
   std::complex<double> tmp_3164;
   std::complex<double> tmp_3165;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3165 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3164 += tmp_3165;
   tmp_3149 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3164;
   std::complex<double> tmp_3166;
   std::complex<double> tmp_3167;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3167 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3166 += tmp_3167;
   tmp_3149 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3166;
   result += (std::complex<double>(0,-1)) * tmp_3149;

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_3168;
   std::complex<double> tmp_3169;
   std::complex<double> tmp_3170;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3170 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3169 += tmp_3170;
   tmp_3168 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qq)*Sqr(Cos(ThetaWp())))
      * tmp_3169;
   std::complex<double> tmp_3171;
   std::complex<double> tmp_3172;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3172 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3171 += tmp_3172;
   tmp_3168 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qu)*Sqr(Cos(ThetaWp())))
      * tmp_3171;
   std::complex<double> tmp_3173;
   std::complex<double> tmp_3174;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3174 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3173 += tmp_3174;
   tmp_3168 += (std::complex<double>(0,-2)*g2*gp*Qq*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_3173;
   std::complex<double> tmp_3175;
   std::complex<double> tmp_3176;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3176 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3175 += tmp_3176;
   tmp_3168 += (std::complex<double>(0.,0.5163977794943222)*g1*gp*Qq*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3175;
   std::complex<double> tmp_3177;
   std::complex<double> tmp_3178;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3178 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3177 += tmp_3178;
   tmp_3168 += (std::complex<double>(0.,-2.065591117977289)*g1*gp*Qu*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3177;
   std::complex<double> tmp_3179;
   std::complex<double> tmp_3180;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3180 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3179 += tmp_3180;
   tmp_3168 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3179;
   std::complex<double> tmp_3181;
   std::complex<double> tmp_3182;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3182 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3181 += tmp_3182;
   tmp_3168 += (std::complex<double>(0.,-0.2581988897471611)*g1*g2*Cos(ThetaW()
      )*Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_3181;
   std::complex<double> tmp_3183;
   std::complex<double> tmp_3184;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3184 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3183 += tmp_3184;
   tmp_3168 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_3183;
   std::complex<double> tmp_3185;
   std::complex<double> tmp_3186;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3186 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3185 += tmp_3186;
   tmp_3168 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))*Sqr(Sin(ThetaWp()))) * tmp_3185;
   result += (std::complex<double>(0,-1)) * tmp_3168;

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZpconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_3187;
   std::complex<double> tmp_3188;
   std::complex<double> tmp_3189;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3189 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3188 += tmp_3189;
   tmp_3187 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Ql)*Sqr(Cos(ThetaWp())))
      * tmp_3188;
   std::complex<double> tmp_3190;
   std::complex<double> tmp_3191;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3191 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3190 += tmp_3191;
   tmp_3187 += (std::complex<double>(0,2)*Sqr(gp)*Sqr(Qv)*Sqr(Cos(ThetaWp())))
      * tmp_3190;
   std::complex<double> tmp_3192;
   std::complex<double> tmp_3193;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3193 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3192 += tmp_3193;
   tmp_3187 += (std::complex<double>(0,-2)*g2*gp*Ql*Cos(ThetaW())*Cos(ThetaWp()
      )*Sin(ThetaWp())) * tmp_3192;
   std::complex<double> tmp_3194;
   std::complex<double> tmp_3195;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3195 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3194 += tmp_3195;
   tmp_3187 += (std::complex<double>(0.,-1.5491933384829668)*g1*gp*Ql*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3194;
   std::complex<double> tmp_3196;
   std::complex<double> tmp_3197;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3197 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3196 += tmp_3197;
   tmp_3187 += (std::complex<double>(0.,3.0983866769659336)*g1*gp*Qv*Cos(
      ThetaWp())*Sin(ThetaW())*Sin(ThetaWp())) * tmp_3196;
   std::complex<double> tmp_3198;
   std::complex<double> tmp_3199;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3199 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3198 += tmp_3199;
   tmp_3187 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3198;
   std::complex<double> tmp_3200;
   std::complex<double> tmp_3201;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3201 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3200 += tmp_3201;
   tmp_3187 += (std::complex<double>(0.,0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())*Sqr(Sin(ThetaWp()))) * tmp_3200;
   std::complex<double> tmp_3202;
   std::complex<double> tmp_3203;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3203 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3202 += tmp_3203;
   tmp_3187 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3202;
   std::complex<double> tmp_3204;
   std::complex<double> tmp_3205;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3205 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3204 += tmp_3205;
   tmp_3187 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))*Sqr(Sin(
      ThetaWp()))) * tmp_3204;
   result += (std::complex<double>(0,-1)) * tmp_3187;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSdSd(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_3206;
   std::complex<double> tmp_3207;
   std::complex<double> tmp_3208;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3208 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_3207 += tmp_3208;
   tmp_3206 += (-2*(3*gp*Qd*Cos(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW()
      )*Sin(ThetaWp()))) * tmp_3207;
   std::complex<double> tmp_3209;
   std::complex<double> tmp_3210;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3210 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3209 += tmp_3210;
   tmp_3206 += (6*gp*Qq*Cos(ThetaWp()) + (3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_3209;
   result += (0.16666666666666666) * tmp_3206;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSeSe(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_3211;
   std::complex<double> tmp_3212;
   std::complex<double> tmp_3213;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3213 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_3212 += tmp_3213;
   tmp_3211 += (-2*(gp*Qe*Cos(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()))) * tmp_3212;
   std::complex<double> tmp_3214;
   std::complex<double> tmp_3215;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3215 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3214 += tmp_3215;
   tmp_3211 += (2*gp*Ql*Cos(ThetaWp()) + (g2*Cos(ThetaW()) - 0.7745966692414834
      *g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_3214;
   result += (0.5) * tmp_3211;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSuSu(unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_3216;
   std::complex<double> tmp_3217;
   std::complex<double> tmp_3218;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3218 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_3217 += tmp_3218;
   tmp_3216 += (2*(-3*gp*Qu*Cos(ThetaWp()) + 1.5491933384829668*g1*Sin(ThetaW()
      )*Sin(ThetaWp()))) * tmp_3217;
   std::complex<double> tmp_3219;
   std::complex<double> tmp_3220;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3220 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3219 += tmp_3220;
   tmp_3216 += (6*gp*Qq*Cos(ThetaWp()) + (-3*g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_3219;
   result += (0.16666666666666666) * tmp_3216;

   return result;
}

std::complex<double> CLASSNAME::CpVZpconjSvSv(unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_3221;
   std::complex<double> tmp_3222;
   std::complex<double> tmp_3223;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3223 += Conj(ZV(gI2,3 + j1))*ZV(gI1,3 + j1);
   }
   tmp_3222 += tmp_3223;
   tmp_3221 += (-2*(gp*Qv*Cos(ThetaWp()) + 0.7745966692414834*g1*Sin(ThetaW())*
      Sin(ThetaWp()))) * tmp_3222;
   std::complex<double> tmp_3224;
   std::complex<double> tmp_3225;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3225 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3224 += tmp_3225;
   tmp_3221 += (2*gp*Ql*Cos(ThetaWp()) - (g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*Sin(ThetaWp())) * tmp_3224;
   result += (0.5) * tmp_3221;

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
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.1*(20*vS*Conj(ZH(gI2,2))*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(gp)*
      Sqr(Qs) - vd*Conj(ZH(gI2,0))*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(
      Cos(ThetaW())) - 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Cos(ThetaWp()
      )) + Cos(ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHd) + 3*Sqr(g1)*Sqr(Sin
      (ThetaW()))) + 7.745966692414834*g1*gp*QHd*Sin(ThetaW())*Sqr(Sin(ThetaWp()))
      + 2*g2*Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin
      (ThetaWp()) - 5*gp*QHd*Sqr(Cos(ThetaWp())) + 5*gp*QHd*Sqr(Sin(ThetaWp()))))
      - vu*Conj(ZH(gI2,1))*(5*Cos(ThetaWp())*Sin(ThetaWp())*Sqr(g2)*Sqr(Cos(ThetaW
      ())) + 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Cos(ThetaWp())) + Cos(
      ThetaWp())*Sin(ThetaWp())*(-20*Sqr(gp)*Sqr(QHu) + 3*Sqr(g1)*Sqr(Sin(ThetaW()
      ))) - 7.745966692414834*g1*gp*QHu*Sin(ThetaW())*Sqr(Sin(ThetaWp())) + 2*g2*
      Cos(ThetaW())*(3.872983346207417*g1*Cos(ThetaWp())*Sin(ThetaW())*Sin(ThetaWp
      ()) + 5*gp*QHu*Sqr(Cos(ThetaWp())) - 5*gp*QHu*Sqr(Sin(ThetaWp())))));

   return result;
}

std::complex<double> CLASSNAME::CpVZpVZphh(unsigned gI2) const
{
   const auto Qs = LOCALINPUT(Qs);
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);

   std::complex<double> result;

   result = 0.5*(4*vS*Conj(ZH(gI2,2))*Sqr(gp)*Sqr(Qs)*Sqr(Cos(ThetaWp())) + vd*
      Conj(ZH(gI2,0))*Sqr(-2*gp*QHd*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp(
      )) + 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())) + vu*Conj(ZH(gI2,1)
      )*Sqr(2*gp*QHu*Cos(ThetaWp()) + g2*Cos(ThetaW())*Sin(ThetaWp()) +
      0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())));

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

   result = 0.5*g2*(Conj(ZH(gI2,0))*ZP(gI1,0) - Conj(ZH(gI2,1))*ZP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZA(gI1,0))*Conj(ZA(gI2,0)) + Conj(ZA(gI1,1))*Conj(ZA(gI2,
      1)))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(ZH(gI1,0))*Conj(ZH(gI2,0)) + Conj(ZH(gI1,1))*Conj(ZH(gI2,
      1)))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3226;
   std::complex<double> tmp_3227;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3227 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_3226 += tmp_3227;
   result += (-0.7071067811865475*g2) * tmp_3226;

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

   std::complex<double> tmp_3228;
   std::complex<double> tmp_3229;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3229 += Conj(ZEL(gI2,j1))*ZVL(gI1,j1);
   }
   tmp_3228 += tmp_3229;
   result += (-0.7071067811865475*g2) * tmp_3228;

   return result;
}

double CLASSNAME::CpconjVWmbarFvFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3230;
   std::complex<double> tmp_3231;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3231 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3230 += tmp_3231;
   result += (0.5*Sqr(g2)) * tmp_3230;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3232;
   std::complex<double> tmp_3233;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3233 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3232 += tmp_3233;
   result += (0.5*Sqr(g2)) * tmp_3232;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3234;
   std::complex<double> tmp_3235;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3235 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3234 += tmp_3235;
   result += (0.5*Sqr(g2)) * tmp_3234;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3236;
   std::complex<double> tmp_3237;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3237 += Conj(ZV(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3236 += tmp_3237;
   result += (0.5*Sqr(g2)) * tmp_3236;

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

   std::complex<double> tmp_3238;
   std::complex<double> tmp_3239;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3239 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3238 += tmp_3239;
   result += (0.7071067811865475*g2) * tmp_3238;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSvSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3240;
   std::complex<double> tmp_3241;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3241 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3240 += tmp_3241;
   result += (0.7071067811865475*g2) * tmp_3240;

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

   result = 0.5*(vd*Conj(ZH(gI2,0)) + vu*Conj(ZH(gI2,1)))*Sqr(g2);

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

std::complex<double> CLASSNAME::CpUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(5*Conj(ZH(gI1,2))*(-2*gp*Qs*Conj(ZN(gI2,5))*KroneckerDelta(0,
      gO2) - 2*gp*Qs*Conj(ZN(gI2,0))*KroneckerDelta(5,gO2) + 1.4142135623730951*
      Conj(ZN(gI2,4))*KroneckerDelta(3,gO2)*Lambdax + 1.4142135623730951*Conj(ZN(
      gI2,3))*KroneckerDelta(4,gO2)*Lambdax) - Conj(ZH(gI1,1))*(Conj(ZN(gI2,4))*(
      10*gp*QHu*KroneckerDelta(0,gO2) + 3.872983346207417*g1*KroneckerDelta(1,gO2)
      - 5*g2*KroneckerDelta(2,gO2)) + 10*gp*QHu*Conj(ZN(gI2,0))*KroneckerDelta(4,
      gO2) + 3.872983346207417*g1*Conj(ZN(gI2,1))*KroneckerDelta(4,gO2) - 5*g2*
      Conj(ZN(gI2,2))*KroneckerDelta(4,gO2) - 7.0710678118654755*Conj(ZN(gI2,5))*
      KroneckerDelta(3,gO2)*Lambdax - 7.0710678118654755*Conj(ZN(gI2,3))*
      KroneckerDelta(5,gO2)*Lambdax) + Conj(ZH(gI1,0))*(Conj(ZN(gI2,3))*(-10*gp*
      QHd*KroneckerDelta(0,gO2) + 3.872983346207417*g1*KroneckerDelta(1,gO2) - 5*
      g2*KroneckerDelta(2,gO2)) - 10*gp*QHd*Conj(ZN(gI2,0))*KroneckerDelta(3,gO2)
      + 3.872983346207417*g1*Conj(ZN(gI2,1))*KroneckerDelta(3,gO2) - 5*g2*Conj(ZN(
      gI2,2))*KroneckerDelta(3,gO2) + 7.0710678118654755*Conj(ZN(gI2,5))*
      KroneckerDelta(4,gO2)*Lambdax + 7.0710678118654755*Conj(ZN(gI2,4))*
      KroneckerDelta(5,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto QHd = LOCALINPUT(QHd);
   const auto QHu = LOCALINPUT(QHu);
   const auto Qs = LOCALINPUT(Qs);

   std::complex<double> result;

   result = 0.1*(5*Conj(ZH(gI1,2))*(-2*gp*Qs*KroneckerDelta(5,gO1)*ZN(gI2,0) +
      1.4142135623730951*Conj(Lambdax)*(KroneckerDelta(4,gO1)*ZN(gI2,3) +
      KroneckerDelta(3,gO1)*ZN(gI2,4)) - 2*gp*Qs*KroneckerDelta(0,gO1)*ZN(gI2,5))
      + Conj(ZH(gI1,0))*(KroneckerDelta(3,gO1)*(-10*gp*QHd*ZN(gI2,0) +
      3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(gI2,2)) - 10*gp*QHd*KroneckerDelta(
      0,gO1)*ZN(gI2,3) + 3.872983346207417*g1*KroneckerDelta(1,gO1)*ZN(gI2,3) - 5*
      g2*KroneckerDelta(2,gO1)*ZN(gI2,3) + 7.0710678118654755*Conj(Lambdax)*
      KroneckerDelta(5,gO1)*ZN(gI2,4) + 7.0710678118654755*Conj(Lambdax)*
      KroneckerDelta(4,gO1)*ZN(gI2,5)) - Conj(ZH(gI1,1))*(KroneckerDelta(4,gO1)*(
      10*gp*QHu*ZN(gI2,0) + 3.872983346207417*g1*ZN(gI2,1) - 5*g2*ZN(gI2,2)) + (10
      *gp*QHu*KroneckerDelta(0,gO1) + 3.872983346207417*g1*KroneckerDelta(1,gO1) -
      5*g2*KroneckerDelta(2,gO1))*ZN(gI2,4) - 7.0710678118654755*Conj(Lambdax)*(
      KroneckerDelta(5,gO1)*ZN(gI2,3) + KroneckerDelta(3,gO1)*ZN(gI2,5))));

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

   std::complex<double> tmp_3242;
   std::complex<double> tmp_3243;
   std::complex<double> tmp_3244;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3244 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3243 += tmp_3244;
   tmp_3242 += (std::complex<double>(0.,-1.4142135623730951)*gp*Qq*
      KroneckerDelta(0,gO2)) * tmp_3243;
   std::complex<double> tmp_3245;
   std::complex<double> tmp_3246;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3246 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3245 += tmp_3246;
   tmp_3242 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (1,gO2)) * tmp_3245;
   std::complex<double> tmp_3247;
   std::complex<double> tmp_3248;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3248 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3247 += tmp_3248;
   tmp_3242 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(2
      ,gO2)) * tmp_3247;
   std::complex<double> tmp_3249;
   std::complex<double> tmp_3250;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3251;
      std::complex<double> tmp_3252;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3252 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_3251 += tmp_3252;
      tmp_3250 += (Conj(ZDL(gI2,j2))) * tmp_3251;
   }
   tmp_3249 += tmp_3250;
   tmp_3242 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_3249;
   result += (std::complex<double>(0,-1)) * tmp_3242;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_3253;
   std::complex<double> tmp_3254;
   std::complex<double> tmp_3255;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3255 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_3254 += tmp_3255;
   tmp_3253 += (std::complex<double>(0.,-1.4142135623730951)*gp*Qd*
      KroneckerDelta(0,gO1)) * tmp_3254;
   std::complex<double> tmp_3256;
   std::complex<double> tmp_3257;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3257 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_3256 += tmp_3257;
   tmp_3253 += (std::complex<double>(0.,-0.3651483716701107)*g1*KroneckerDelta(
      1,gO1)) * tmp_3256;
   std::complex<double> tmp_3258;
   std::complex<double> tmp_3259;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3260;
      std::complex<double> tmp_3261;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3261 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_3260 += tmp_3261;
      tmp_3259 += (ZD(gI1,j2)) * tmp_3260;
   }
   tmp_3258 += tmp_3259;
   tmp_3253 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_3258;
   result += (std::complex<double>(0,-1)) * tmp_3253;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_3262;
   std::complex<double> tmp_3263;
   std::complex<double> tmp_3264;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3264 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3263 += tmp_3264;
   tmp_3262 += (std::complex<double>(0.,-1.4142135623730951)*gp*Ql*
      KroneckerDelta(0,gO2)) * tmp_3263;
   std::complex<double> tmp_3265;
   std::complex<double> tmp_3266;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3266 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3265 += tmp_3266;
   tmp_3262 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(1
      ,gO2)) * tmp_3265;
   std::complex<double> tmp_3267;
   std::complex<double> tmp_3268;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3268 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_3267 += tmp_3268;
   tmp_3262 += (std::complex<double>(0.,0.7071067811865475)*g2*KroneckerDelta(2
      ,gO2)) * tmp_3267;
   std::complex<double> tmp_3269;
   std::complex<double> tmp_3270;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3271;
      std::complex<double> tmp_3272;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3272 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_3271 += tmp_3272;
      tmp_3270 += (Conj(ZEL(gI2,j2))) * tmp_3271;
   }
   tmp_3269 += tmp_3270;
   tmp_3262 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_3269;
   result += (std::complex<double>(0,-1)) * tmp_3262;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_3273;
   std::complex<double> tmp_3274;
   std::complex<double> tmp_3275;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3275 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_3274 += tmp_3275;
   tmp_3273 += (std::complex<double>(0.,-1.4142135623730951)*gp*Qe*
      KroneckerDelta(0,gO1)) * tmp_3274;
   std::complex<double> tmp_3276;
   std::complex<double> tmp_3277;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3277 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_3276 += tmp_3277;
   tmp_3273 += (std::complex<double>(0.,-1.0954451150103321)*g1*KroneckerDelta(
      1,gO1)) * tmp_3276;
   std::complex<double> tmp_3278;
   std::complex<double> tmp_3279;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3280;
      std::complex<double> tmp_3281;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3281 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_3280 += tmp_3281;
      tmp_3279 += (ZE(gI1,j2)) * tmp_3280;
   }
   tmp_3278 += tmp_3279;
   tmp_3273 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_3278;
   result += (std::complex<double>(0,-1)) * tmp_3273;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_3282;
   std::complex<double> tmp_3283;
   std::complex<double> tmp_3284;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3284 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3283 += tmp_3284;
   tmp_3282 += (std::complex<double>(0.,-1.4142135623730951)*gp*Qq*
      KroneckerDelta(0,gO2)) * tmp_3283;
   std::complex<double> tmp_3285;
   std::complex<double> tmp_3286;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3286 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3285 += tmp_3286;
   tmp_3282 += (std::complex<double>(0.,-0.18257418583505536)*g1*KroneckerDelta
      (1,gO2)) * tmp_3285;
   std::complex<double> tmp_3287;
   std::complex<double> tmp_3288;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3288 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3287 += tmp_3288;
   tmp_3282 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      2,gO2)) * tmp_3287;
   std::complex<double> tmp_3289;
   std::complex<double> tmp_3290;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3291;
      std::complex<double> tmp_3292;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3292 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_3291 += tmp_3292;
      tmp_3290 += (Conj(ZUL(gI2,j2))) * tmp_3291;
   }
   tmp_3289 += tmp_3290;
   tmp_3282 += (std::complex<double>(0,-1)*KroneckerDelta(4,gO2)) * tmp_3289;
   result += (std::complex<double>(0,-1)) * tmp_3282;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_3293;
   std::complex<double> tmp_3294;
   std::complex<double> tmp_3295;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3295 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_3294 += tmp_3295;
   tmp_3293 += (std::complex<double>(0.,-1.4142135623730951)*gp*Qu*
      KroneckerDelta(0,gO1)) * tmp_3294;
   std::complex<double> tmp_3296;
   std::complex<double> tmp_3297;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3297 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_3296 += tmp_3297;
   tmp_3293 += (std::complex<double>(0.,0.7302967433402214)*g1*KroneckerDelta(1
      ,gO1)) * tmp_3296;
   std::complex<double> tmp_3298;
   std::complex<double> tmp_3299;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3300;
      std::complex<double> tmp_3301;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3301 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_3300 += tmp_3301;
      tmp_3299 += (ZU(gI1,j2)) * tmp_3300;
   }
   tmp_3298 += tmp_3299;
   tmp_3293 += (std::complex<double>(0,-1)*KroneckerDelta(4,gO1)) * tmp_3298;
   result += (std::complex<double>(0,-1)) * tmp_3293;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSvFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_3302;
   std::complex<double> tmp_3303;
   std::complex<double> tmp_3304;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3304 += Conj(ZVL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3303 += tmp_3304;
   tmp_3302 += (std::complex<double>(0.,-1.4142135623730951)*gp*Ql*
      KroneckerDelta(0,gO2)) * tmp_3303;
   std::complex<double> tmp_3305;
   std::complex<double> tmp_3306;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3306 += Conj(ZVL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3305 += tmp_3306;
   tmp_3302 += (std::complex<double>(0.,0.5477225575051661)*g1*KroneckerDelta(1
      ,gO2)) * tmp_3305;
   std::complex<double> tmp_3307;
   std::complex<double> tmp_3308;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3308 += Conj(ZVL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3307 += tmp_3308;
   tmp_3302 += (std::complex<double>(0.,-0.7071067811865475)*g2*KroneckerDelta(
      2,gO2)) * tmp_3307;
   std::complex<double> tmp_3309;
   std::complex<double> tmp_3310;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3311;
      std::complex<double> tmp_3312;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3312 += Conj(ZVL(gI2,j1))*Yv(j1,j2);
      }
      tmp_3311 += tmp_3312;
      tmp_3310 += (ZV(gI1,3 + j2)) * tmp_3311;
   }
   tmp_3309 += tmp_3310;
   tmp_3302 += (std::complex<double>(0,-1)*KroneckerDelta(4,gO2)) * tmp_3309;
   result += (std::complex<double>(0,-1)) * tmp_3302;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSvFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   std::complex<double> tmp_3313;
   std::complex<double> tmp_3314;
   std::complex<double> tmp_3315;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3315 += ZV(gI1,3 + j1)*ZVR(gI2,j1);
   }
   tmp_3314 += tmp_3315;
   tmp_3313 += (std::complex<double>(0.,-1.4142135623730951)*gp*Qv*
      KroneckerDelta(0,gO1)) * tmp_3314;
   std::complex<double> tmp_3316;
   std::complex<double> tmp_3317;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3317 += ZV(gI1,3 + j1)*ZVR(gI2,j1);
   }
   tmp_3316 += tmp_3317;
   tmp_3313 += (std::complex<double>(0.,-1.0954451150103321)*g1*KroneckerDelta(
      1,gO1)) * tmp_3316;
   std::complex<double> tmp_3318;
   std::complex<double> tmp_3319;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3320;
      std::complex<double> tmp_3321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3321 += Conj(Yv(j1,j2))*ZV(gI1,j1);
      }
      tmp_3320 += tmp_3321;
      tmp_3319 += (ZVR(gI2,j2)) * tmp_3320;
   }
   tmp_3318 += tmp_3319;
   tmp_3313 += (std::complex<double>(0,-1)*KroneckerDelta(4,gO1)) * tmp_3318;
   result += (std::complex<double>(0,-1)) * tmp_3313;

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

std::complex<double> CLASSNAME::CpbarUFvconjHpmFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3322;
      std::complex<double> tmp_3323;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3323 += Conj(ZEL(gI2,j1))*Yv(j1,gO2);
      }
      tmp_3322 += tmp_3323;
      result += (ZP(gI1,1)) * tmp_3322;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvconjHpmFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3324;
      std::complex<double> tmp_3325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3325 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_3324 += tmp_3325;
      result += (ZP(gI1,0)) * tmp_3324;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvbarChaSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3326;
      std::complex<double> tmp_3327;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3327 += Conj(ZE(gI2,j1))*Yv(j1,gO2);
      }
      tmp_3326 += tmp_3327;
      result += (Conj(UP(gI1,1))) * tmp_3326;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvbarChaSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZE(gI2,gO1))*UM(gI1,0));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_3328;
      std::complex<double> tmp_3329;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3329 += Conj(Ye(j1,gO1))*Conj(ZE(gI2,3 + j1));
      }
      tmp_3328 += tmp_3329;
      result += (UM(gI1,1)) * tmp_3328;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3330;
      std::complex<double> tmp_3331;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3331 += Conj(ZVL(gI1,j1))*Yv(j1,gO2);
      }
      tmp_3330 += tmp_3331;
      result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,1)
         )) * tmp_3330;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvFvAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3332;
      std::complex<double> tmp_3333;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3333 += Conj(Yv(gO1,j2))*ZVR(gI1,j2);
      }
      tmp_3332 += tmp_3333;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))
         ) * tmp_3332;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvhhFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3334;
      std::complex<double> tmp_3335;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3335 += Conj(ZVL(gI2,j1))*Yv(j1,gO2);
      }
      tmp_3334 += tmp_3335;
      result += (-0.7071067811865475*Conj(ZH(gI1,1))) * tmp_3334;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvhhFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3336;
      std::complex<double> tmp_3337;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3337 += Conj(Yv(gO1,j2))*ZVR(gI2,j2);
      }
      tmp_3336 += tmp_3337;
      result += (-0.7071067811865475*Conj(ZH(gI1,1))) * tmp_3336;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvSvChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.4142135623730951*gp*Qv*Conj(ZN(gI2,0))*Conj(ZV(gI1,3 +
         gO2));
   }
   if (gO2 < 3) {
      result += -1.0954451150103321*g1*Conj(ZN(gI2,1))*Conj(ZV(gI1,3 + gO2))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_3338;
      std::complex<double> tmp_3339;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3339 += Conj(ZV(gI1,j1))*Yv(j1,gO2);
      }
      tmp_3338 += tmp_3339;
      result += (-Conj(ZN(gI2,4))) * tmp_3338;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvSvChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*gp*Ql*Conj(ZV(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZV(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZV(gI1,gO1))*ZN(gI2,2);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_3340;
      std::complex<double> tmp_3341;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3341 += Conj(Yv(gO1,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_3340 += tmp_3341;
      result += (-ZN(gI2,4)) * tmp_3340;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvVPFvPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Cos(ThetaW())*ZVR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvVPFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZVL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZVL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvVZFvPR(unsigned gO2, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Cos(ThetaWp())*Sin(ThetaW())*ZVR(gI2,
         gO2);
   }
   if (gI2 < 3) {
      result += gp*Qv*Sin(ThetaWp())*ZVR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvVZFvPL(unsigned gO1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZVL(gI2,gO1))*Cos(ThetaW())*Cos(ThetaWp());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZVL(gI2,gO1))*Cos(ThetaWp())*Sin
         (ThetaW());
   }
   if (gI2 < 3) {
      result += -(gp*Ql*Conj(ZVL(gI2,gO1))*Sin(ThetaWp()));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvVZpFvPR(unsigned gO2, unsigned gI2) const
{
   const auto Qv = LOCALINPUT(Qv);

   std::complex<double> result;

   if (gI2 < 3) {
      result += gp*Qv*Cos(ThetaWp())*ZVR(gI2,gO2);
   }
   if (gI2 < 3) {
      result += 0.7745966692414834*g1*Sin(ThetaW())*Sin(ThetaWp())*ZVR(gI2,
         gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvVZpFvPL(unsigned gO1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   if (gI2 < 3) {
      result += -(gp*Ql*Conj(ZVL(gI2,gO1))*Cos(ThetaWp()));
   }
   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZVL(gI2,gO1))*Cos(ThetaW())*Sin(ThetaWp());
   }
   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZVL(gI2,gO1))*Sin(ThetaW())*Sin(
         ThetaWp());
   }

   return result;
}

double CLASSNAME::CpbarUFvconjVWmFePR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFvconjVWmFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZEL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,0.7071067811865475)*(g2*Conj(UM(gI1,0))*
      Conj(ZA(gI2,1))*KroneckerDelta(1,gO2) + Conj(UM(gI1,1))*(g2*Conj(ZA(gI2,0))*
      KroneckerDelta(0,gO2) - Conj(ZA(gI2,2))*KroneckerDelta(1,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0.,-0.7071067811865475)*(g2*Conj(ZA(gI2,0))*
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

   result = -0.7071067811865475*(g2*Conj(UM(gI2,0))*Conj(ZH(gI1,1))*
      KroneckerDelta(1,gO2) + Conj(UM(gI2,1))*(g2*Conj(ZH(gI1,0))*KroneckerDelta(0
      ,gO2) + Conj(ZH(gI1,2))*KroneckerDelta(1,gO2)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*(g2*Conj(ZH(gI1,0))*KroneckerDelta(1,gO1)*UP(
      gI2,0) + (g2*Conj(ZH(gI1,1))*KroneckerDelta(0,gO1) + Conj(Lambdax)*Conj(ZH(
      gI1,2))*KroneckerDelta(1,gO1))*UP(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3342;
   std::complex<double> tmp_3343;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3344;
      std::complex<double> tmp_3345;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3345 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_3344 += tmp_3345;
      tmp_3343 += (Conj(ZD(gI2,j2))) * tmp_3344;
   }
   tmp_3342 += tmp_3343;
   result += (KroneckerDelta(1,gO2)) * tmp_3342;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3346;
   std::complex<double> tmp_3347;
   std::complex<double> tmp_3348;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3348 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_3347 += tmp_3348;
   tmp_3346 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_3347
      ;
   std::complex<double> tmp_3349;
   std::complex<double> tmp_3350;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3351;
      std::complex<double> tmp_3352;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3352 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_3351 += tmp_3352;
      tmp_3350 += (ZUL(gI1,j2)) * tmp_3351;
   }
   tmp_3349 += tmp_3350;
   tmp_3346 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_3349;
   result += (std::complex<double>(0,-1)) * tmp_3346;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3353;
   std::complex<double> tmp_3354;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3355;
      std::complex<double> tmp_3356;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3356 += Conj(ZE(gI2,j1))*Yv(j1,j2);
      }
      tmp_3355 += tmp_3356;
      tmp_3354 += (Conj(ZVR(gI1,j2))) * tmp_3355;
   }
   tmp_3353 += tmp_3354;
   result += (KroneckerDelta(1,gO2)) * tmp_3353;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFvSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3357;
   std::complex<double> tmp_3358;
   std::complex<double> tmp_3359;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3359 += Conj(ZE(gI2,j1))*ZVL(gI1,j1);
   }
   tmp_3358 += tmp_3359;
   tmp_3357 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_3358
      ;
   std::complex<double> tmp_3360;
   std::complex<double> tmp_3361;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3362;
      std::complex<double> tmp_3363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3363 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_3362 += tmp_3363;
      tmp_3361 += (ZVL(gI1,j2)) * tmp_3362;
   }
   tmp_3360 += tmp_3361;
   tmp_3357 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_3360;
   result += (std::complex<double>(0,-1)) * tmp_3357;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3364;
   std::complex<double> tmp_3365;
   std::complex<double> tmp_3366;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3366 += Conj(ZDL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3365 += tmp_3366;
   tmp_3364 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_3365
      ;
   std::complex<double> tmp_3367;
   std::complex<double> tmp_3368;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3369;
      std::complex<double> tmp_3370;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3370 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_3369 += tmp_3370;
      tmp_3368 += (Conj(ZDL(gI2,j2))) * tmp_3369;
   }
   tmp_3367 += tmp_3368;
   tmp_3364 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_3367;
   result += (std::complex<double>(0,-1)) * tmp_3364;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3371;
   std::complex<double> tmp_3372;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3373;
      std::complex<double> tmp_3374;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3374 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_3373 += tmp_3374;
      tmp_3372 += (ZU(gI1,j2)) * tmp_3373;
   }
   tmp_3371 += tmp_3372;
   result += (KroneckerDelta(1,gO1)) * tmp_3371;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3375;
   std::complex<double> tmp_3376;
   std::complex<double> tmp_3377;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3377 += Conj(ZEL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_3376 += tmp_3377;
   tmp_3375 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_3376
      ;
   std::complex<double> tmp_3378;
   std::complex<double> tmp_3379;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3380;
      std::complex<double> tmp_3381;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3381 += Conj(ZEL(gI2,j1))*Yv(j1,j2);
      }
      tmp_3380 += tmp_3381;
      tmp_3379 += (ZV(gI1,3 + j2)) * tmp_3380;
   }
   tmp_3378 += tmp_3379;
   tmp_3375 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_3378;
   result += (std::complex<double>(0,-1)) * tmp_3375;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3382;
   std::complex<double> tmp_3383;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3384;
      std::complex<double> tmp_3385;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3385 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_3384 += tmp_3385;
      tmp_3383 += (ZV(gI1,j2)) * tmp_3384;
   }
   tmp_3382 += tmp_3383;
   result += (KroneckerDelta(1,gO1)) * tmp_3382;

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
      std::complex<double> tmp_3386;
      std::complex<double> tmp_3387;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3387 += Conj(ZVL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_3386 += tmp_3387;
      result += (ZP(gI1,0)) * tmp_3386;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeHpmFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3388;
      std::complex<double> tmp_3389;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3389 += Conj(Yv(gO1,j2))*ZVR(gI2,j2);
      }
      tmp_3388 += tmp_3389;
      result += (ZP(gI1,1)) * tmp_3388;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3390;
      std::complex<double> tmp_3391;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3391 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_3390 += tmp_3391;
      result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,0)
         )) * tmp_3390;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3392;
      std::complex<double> tmp_3393;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3393 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_3392 += tmp_3393;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))
         ) * tmp_3392;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3394;
      std::complex<double> tmp_3395;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3395 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_3394 += tmp_3395;
      result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3394;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3396;
      std::complex<double> tmp_3397;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3397 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_3396 += tmp_3397;
      result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3396;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3398;
      std::complex<double> tmp_3399;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3399 += Conj(ZV(gI1,j2))*Ye(gO2,j2);
      }
      tmp_3398 += tmp_3399;
      result += (Conj(UM(gI2,1))) * tmp_3398;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZV(gI1,gO1))*UP(gI2,0));
   }
   if (gO1 < 3) {
      std::complex<double> tmp_3400;
      std::complex<double> tmp_3401;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3401 += Conj(Yv(gO1,j2))*Conj(ZV(gI1,3 + j2));
      }
      tmp_3400 += tmp_3401;
      result += (UP(gI2,1)) * tmp_3400;
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
      std::complex<double> tmp_3402;
      std::complex<double> tmp_3403;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3403 += Conj(ZE(gI1,j2))*Ye(gO2,j2);
      }
      tmp_3402 += tmp_3403;
      result += (-Conj(ZN(gI2,3))) * tmp_3402;
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
      std::complex<double> tmp_3404;
      std::complex<double> tmp_3405;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3405 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_3404 += tmp_3405;
      result += (-ZN(gI2,3)) * tmp_3404;
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

std::complex<double> CLASSNAME::CpbarUFeVWmFvPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZVL(gI2,gO1));
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
      std::complex<double> tmp_3406;
      std::complex<double> tmp_3407;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3407 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_3406 += tmp_3407;
      result += (ZP(gI1,0)) * tmp_3406;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3408;
      std::complex<double> tmp_3409;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3409 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_3408 += tmp_3409;
      result += (ZP(gI1,1)) * tmp_3408;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3410;
      std::complex<double> tmp_3411;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3411 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_3410 += tmp_3411;
      result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,0)
         )) * tmp_3410;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3412;
      std::complex<double> tmp_3413;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3413 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_3412 += tmp_3413;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))
         ) * tmp_3412;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3414;
      std::complex<double> tmp_3415;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3415 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_3414 += tmp_3415;
      result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3414;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3416;
      std::complex<double> tmp_3417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3417 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_3416 += tmp_3417;
      result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3416;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3418;
      std::complex<double> tmp_3419;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3419 += Conj(ZU(gI1,j2))*Yd(gO2,j2);
      }
      tmp_3418 += tmp_3419;
      result += (Conj(UM(gI2,1))) * tmp_3418;
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
      std::complex<double> tmp_3420;
      std::complex<double> tmp_3421;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3421 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_3420 += tmp_3421;
      result += (UP(gI2,1)) * tmp_3420;
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
      std::complex<double> tmp_3422;
      std::complex<double> tmp_3423;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3423 += Conj(ZD(gI1,j2))*Yd(gO2,j2);
      }
      tmp_3422 += tmp_3423;
      result += (-Conj(ZN(gI2,3))) * tmp_3422;
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
      std::complex<double> tmp_3424;
      std::complex<double> tmp_3425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3425 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_3424 += tmp_3425;
      result += (-ZN(gI2,3)) * tmp_3424;
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
      std::complex<double> tmp_3426;
      std::complex<double> tmp_3427;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3427 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_3426 += tmp_3427;
      result += (ZP(gI1,1)) * tmp_3426;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3428;
      std::complex<double> tmp_3429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3429 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_3428 += tmp_3429;
      result += (ZP(gI1,0)) * tmp_3428;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3430;
      std::complex<double> tmp_3431;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3431 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_3430 += tmp_3431;
      result += (Conj(UP(gI1,1))) * tmp_3430;
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
      std::complex<double> tmp_3432;
      std::complex<double> tmp_3433;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3433 += Conj(Yd(j1,gO1))*Conj(ZD(gI2,3 + j1));
      }
      tmp_3432 += tmp_3433;
      result += (UM(gI1,1)) * tmp_3432;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3434;
      std::complex<double> tmp_3435;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3435 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_3434 += tmp_3435;
      result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,1)
         )) * tmp_3434;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3436;
      std::complex<double> tmp_3437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3437 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_3436 += tmp_3437;
      result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))
         ) * tmp_3436;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_3438;
      std::complex<double> tmp_3439;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3439 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_3438 += tmp_3439;
      result += (-0.7071067811865475*Conj(ZH(gI1,1))) * tmp_3438;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_3440;
      std::complex<double> tmp_3441;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3441 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_3440 += tmp_3441;
      result += (-0.7071067811865475*Conj(ZH(gI1,1))) * tmp_3440;
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
      std::complex<double> tmp_3442;
      std::complex<double> tmp_3443;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_3443 += Conj(ZU(gI1,j2))*Yu(gO2,j2);
      }
      tmp_3442 += tmp_3443;
      result += (-Conj(ZN(gI2,4))) * tmp_3442;
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
      std::complex<double> tmp_3444;
      std::complex<double> tmp_3445;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3445 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_3444 += tmp_3445;
      result += (-ZN(gI2,4)) * tmp_3444;
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

   std::complex<double> tmp_3446;
   std::complex<double> tmp_3447;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3447 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_3446 += tmp_3447;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_3446;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3448;
   std::complex<double> tmp_3449;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3449 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_3448 += tmp_3449;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_3448;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3450;
   std::complex<double> tmp_3451;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3451 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_3450 += tmp_3451;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_3450;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3452;
   std::complex<double> tmp_3453;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3453 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_3452 += tmp_3453;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_3452;

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

   std::complex<double> tmp_3454;
   std::complex<double> tmp_3455;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3456;
      std::complex<double> tmp_3457;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3457 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_3456 += tmp_3457;
      tmp_3455 += (Conj(ZVL(gI2,j2))) * tmp_3456;
   }
   tmp_3454 += tmp_3455;
   result += (ZP(gI1,0)) * tmp_3454;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3458;
   std::complex<double> tmp_3459;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3460;
      std::complex<double> tmp_3461;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3461 += Conj(Yv(j1,j2))*ZEL(gO1,j1);
      }
      tmp_3460 += tmp_3461;
      tmp_3459 += (ZVR(gI2,j2)) * tmp_3460;
   }
   tmp_3458 += tmp_3459;
   result += (ZP(gI1,1)) * tmp_3458;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3462;
   std::complex<double> tmp_3463;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3464;
      std::complex<double> tmp_3465;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3465 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_3464 += tmp_3465;
      tmp_3463 += (Conj(ZEL(gI1,j2))) * tmp_3464;
   }
   tmp_3462 += tmp_3463;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_3462;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3466;
   std::complex<double> tmp_3467;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3468;
      std::complex<double> tmp_3469;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3469 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_3468 += tmp_3469;
      tmp_3467 += (ZEL(gO1,j2)) * tmp_3468;
   }
   tmp_3466 += tmp_3467;
   result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_3466;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3470;
   std::complex<double> tmp_3471;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3472;
      std::complex<double> tmp_3473;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3473 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_3472 += tmp_3473;
      tmp_3471 += (Conj(ZEL(gI2,j2))) * tmp_3472;
   }
   tmp_3470 += tmp_3471;
   result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3470;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3474;
   std::complex<double> tmp_3475;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3476;
      std::complex<double> tmp_3477;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3477 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_3476 += tmp_3477;
      tmp_3475 += (ZEL(gO1,j2)) * tmp_3476;
   }
   tmp_3474 += tmp_3475;
   result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3474;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3478;
   std::complex<double> tmp_3479;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3480;
      std::complex<double> tmp_3481;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3481 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_3480 += tmp_3481;
      tmp_3479 += (Conj(ZV(gI1,j2))) * tmp_3480;
   }
   tmp_3478 += tmp_3479;
   result += (Conj(UM(gI2,1))) * tmp_3478;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3482;
   std::complex<double> tmp_3483;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3483 += Conj(ZV(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_3482 += tmp_3483;
   result += (-(g2*UP(gI2,0))) * tmp_3482;
   std::complex<double> tmp_3484;
   std::complex<double> tmp_3485;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3486;
      std::complex<double> tmp_3487;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3487 += Conj(Yv(j1,j2))*ZEL(gO1,j1);
      }
      tmp_3486 += tmp_3487;
      tmp_3485 += (Conj(ZV(gI1,3 + j2))) * tmp_3486;
   }
   tmp_3484 += tmp_3485;
   result += (UP(gI2,1)) * tmp_3484;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qe = LOCALINPUT(Qe);

   std::complex<double> result;

   std::complex<double> tmp_3488;
   std::complex<double> tmp_3489;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3489 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_3488 += tmp_3489;
   result += (-1.4142135623730951*gp*Qe*Conj(ZN(gI2,0))) * tmp_3488;
   std::complex<double> tmp_3490;
   std::complex<double> tmp_3491;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3491 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_3490 += tmp_3491;
   result += (-1.0954451150103321*g1*Conj(ZN(gI2,1))) * tmp_3490;
   std::complex<double> tmp_3492;
   std::complex<double> tmp_3493;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3494;
      std::complex<double> tmp_3495;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3495 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_3494 += tmp_3495;
      tmp_3493 += (Conj(ZE(gI1,j2))) * tmp_3494;
   }
   tmp_3492 += tmp_3493;
   result += (-Conj(ZN(gI2,3))) * tmp_3492;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Ql = LOCALINPUT(Ql);

   std::complex<double> result;

   std::complex<double> tmp_3496;
   std::complex<double> tmp_3497;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3497 += Conj(ZE(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_3496 += tmp_3497;
   result += (-0.7071067811865475*(2*gp*Ql*ZN(gI2,0) - 0.7745966692414834*g1*ZN
      (gI2,1) - g2*ZN(gI2,2))) * tmp_3496;
   std::complex<double> tmp_3498;
   std::complex<double> tmp_3499;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3500;
      std::complex<double> tmp_3501;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3501 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_3500 += tmp_3501;
      tmp_3499 += (ZEL(gO1,j2)) * tmp_3500;
   }
   tmp_3498 += tmp_3499;
   result += (-ZN(gI2,3)) * tmp_3498;

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

   std::complex<double> tmp_3502;
   std::complex<double> tmp_3503;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3503 += Conj(ZVL(gI2,j1))*ZEL(gO1,j1);
   }
   tmp_3502 += tmp_3503;
   result += (-0.7071067811865475*g2) * tmp_3502;

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

   std::complex<double> tmp_3504;
   std::complex<double> tmp_3505;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3506;
      std::complex<double> tmp_3507;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3507 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_3506 += tmp_3507;
      tmp_3505 += (Conj(ZUL(gI2,j2))) * tmp_3506;
   }
   tmp_3504 += tmp_3505;
   result += (ZP(gI1,0)) * tmp_3504;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3508;
   std::complex<double> tmp_3509;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3510;
      std::complex<double> tmp_3511;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3511 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_3510 += tmp_3511;
      tmp_3509 += (ZDL(gO1,j2)) * tmp_3510;
   }
   tmp_3508 += tmp_3509;
   result += (ZP(gI1,1)) * tmp_3508;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3512;
   std::complex<double> tmp_3513;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3514;
      std::complex<double> tmp_3515;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3515 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_3514 += tmp_3515;
      tmp_3513 += (Conj(ZDL(gI1,j2))) * tmp_3514;
   }
   tmp_3512 += tmp_3513;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_3512;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3516;
   std::complex<double> tmp_3517;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3518;
      std::complex<double> tmp_3519;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3519 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_3518 += tmp_3519;
      tmp_3517 += (ZDL(gO1,j2)) * tmp_3518;
   }
   tmp_3516 += tmp_3517;
   result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,0))) *
      tmp_3516;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3520;
   std::complex<double> tmp_3521;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3522;
      std::complex<double> tmp_3523;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3523 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_3522 += tmp_3523;
      tmp_3521 += (Conj(ZDL(gI2,j2))) * tmp_3522;
   }
   tmp_3520 += tmp_3521;
   result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3520;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3524;
   std::complex<double> tmp_3525;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3526;
      std::complex<double> tmp_3527;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3527 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_3526 += tmp_3527;
      tmp_3525 += (ZDL(gO1,j2)) * tmp_3526;
   }
   tmp_3524 += tmp_3525;
   result += (-0.7071067811865475*Conj(ZH(gI1,0))) * tmp_3524;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3528;
   std::complex<double> tmp_3529;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3530;
      std::complex<double> tmp_3531;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3531 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_3530 += tmp_3531;
      tmp_3529 += (Conj(ZU(gI1,j2))) * tmp_3530;
   }
   tmp_3528 += tmp_3529;
   result += (Conj(UM(gI2,1))) * tmp_3528;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3532;
   std::complex<double> tmp_3533;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3533 += Conj(ZU(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_3532 += tmp_3533;
   result += (-(g2*UP(gI2,0))) * tmp_3532;
   std::complex<double> tmp_3534;
   std::complex<double> tmp_3535;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3536;
      std::complex<double> tmp_3537;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3537 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_3536 += tmp_3537;
      tmp_3535 += (ZDL(gO1,j2)) * tmp_3536;
   }
   tmp_3534 += tmp_3535;
   result += (UP(gI2,1)) * tmp_3534;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qd = LOCALINPUT(Qd);

   std::complex<double> result;

   std::complex<double> tmp_3538;
   std::complex<double> tmp_3539;
   std::complex<double> tmp_3540;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3540 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_3539 += tmp_3540;
   tmp_3538 += (-4.242640687119286*gp*Qd*Conj(ZN(gI2,0))) * tmp_3539;
   std::complex<double> tmp_3541;
   std::complex<double> tmp_3542;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3542 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_3541 += tmp_3542;
   tmp_3538 += (-1.0954451150103321*g1*Conj(ZN(gI2,1))) * tmp_3541;
   std::complex<double> tmp_3543;
   std::complex<double> tmp_3544;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3545;
      std::complex<double> tmp_3546;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3546 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_3545 += tmp_3546;
      tmp_3544 += (Conj(ZD(gI1,j2))) * tmp_3545;
   }
   tmp_3543 += tmp_3544;
   tmp_3538 += (-3*Conj(ZN(gI2,3))) * tmp_3543;
   result += (0.3333333333333333) * tmp_3538;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_3547;
   std::complex<double> tmp_3548;
   std::complex<double> tmp_3549;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3549 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_3548 += tmp_3549;
   tmp_3547 += (-1.4142135623730951*(6*gp*Qq*ZN(gI2,0) + 0.7745966692414834*g1*
      ZN(gI2,1) - 3*g2*ZN(gI2,2))) * tmp_3548;
   std::complex<double> tmp_3550;
   std::complex<double> tmp_3551;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3552;
      std::complex<double> tmp_3553;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3553 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_3552 += tmp_3553;
      tmp_3551 += (ZDL(gO1,j2)) * tmp_3552;
   }
   tmp_3550 += tmp_3551;
   tmp_3547 += (-6*ZN(gI2,3)) * tmp_3550;
   result += (0.16666666666666666) * tmp_3547;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3554;
   std::complex<double> tmp_3555;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3555 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_3554 += tmp_3555;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_3554;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3556;
   std::complex<double> tmp_3557;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3557 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_3556 += tmp_3557;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_3556;

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

   std::complex<double> tmp_3558;
   std::complex<double> tmp_3559;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3559 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_3558 += tmp_3559;
   result += (-0.7071067811865475*g2) * tmp_3558;

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

   std::complex<double> tmp_3560;
   std::complex<double> tmp_3561;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3562;
      std::complex<double> tmp_3563;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3563 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_3562 += tmp_3563;
      tmp_3561 += (Conj(ZDL(gI2,j2))) * tmp_3562;
   }
   tmp_3560 += tmp_3561;
   result += (ZP(gI1,1)) * tmp_3560;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3564;
   std::complex<double> tmp_3565;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3566;
      std::complex<double> tmp_3567;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3567 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_3566 += tmp_3567;
      tmp_3565 += (ZUL(gO1,j2)) * tmp_3566;
   }
   tmp_3564 += tmp_3565;
   result += (ZP(gI1,0)) * tmp_3564;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3568;
   std::complex<double> tmp_3569;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3570;
      std::complex<double> tmp_3571;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3571 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_3570 += tmp_3571;
      tmp_3569 += (Conj(ZD(gI2,j2))) * tmp_3570;
   }
   tmp_3568 += tmp_3569;
   result += (Conj(UP(gI1,1))) * tmp_3568;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3572;
   std::complex<double> tmp_3573;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3573 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_3572 += tmp_3573;
   result += (-(g2*UM(gI1,0))) * tmp_3572;
   std::complex<double> tmp_3574;
   std::complex<double> tmp_3575;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3576;
      std::complex<double> tmp_3577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3577 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_3576 += tmp_3577;
      tmp_3575 += (ZUL(gO1,j2)) * tmp_3576;
   }
   tmp_3574 += tmp_3575;
   result += (UM(gI1,1)) * tmp_3574;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3578;
   std::complex<double> tmp_3579;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3580;
      std::complex<double> tmp_3581;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3581 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_3580 += tmp_3581;
      tmp_3579 += (Conj(ZUL(gI1,j2))) * tmp_3580;
   }
   tmp_3578 += tmp_3579;
   result += (std::complex<double>(0.,-0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_3578;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3582;
   std::complex<double> tmp_3583;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3584;
      std::complex<double> tmp_3585;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3585 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_3584 += tmp_3585;
      tmp_3583 += (ZUL(gO1,j2)) * tmp_3584;
   }
   tmp_3582 += tmp_3583;
   result += (std::complex<double>(0.,0.7071067811865475)*Conj(ZA(gI2,1))) *
      tmp_3582;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3586;
   std::complex<double> tmp_3587;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3588;
      std::complex<double> tmp_3589;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3589 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_3588 += tmp_3589;
      tmp_3587 += (Conj(ZUL(gI2,j2))) * tmp_3588;
   }
   tmp_3586 += tmp_3587;
   result += (-0.7071067811865475*Conj(ZH(gI1,1))) * tmp_3586;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3590;
   std::complex<double> tmp_3591;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3592;
      std::complex<double> tmp_3593;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3593 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_3592 += tmp_3593;
      tmp_3591 += (ZUL(gO1,j2)) * tmp_3592;
   }
   tmp_3590 += tmp_3591;
   result += (-0.7071067811865475*Conj(ZH(gI1,1))) * tmp_3590;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   const auto Qu = LOCALINPUT(Qu);

   std::complex<double> result;

   std::complex<double> tmp_3594;
   std::complex<double> tmp_3595;
   std::complex<double> tmp_3596;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3596 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_3595 += tmp_3596;
   tmp_3594 += (-4.242640687119286*gp*Qu*Conj(ZN(gI2,0))) * tmp_3595;
   std::complex<double> tmp_3597;
   std::complex<double> tmp_3598;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3598 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_3597 += tmp_3598;
   tmp_3594 += (2.1908902300206643*g1*Conj(ZN(gI2,1))) * tmp_3597;
   std::complex<double> tmp_3599;
   std::complex<double> tmp_3600;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3601;
      std::complex<double> tmp_3602;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3602 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_3601 += tmp_3602;
      tmp_3600 += (Conj(ZU(gI1,j2))) * tmp_3601;
   }
   tmp_3599 += tmp_3600;
   tmp_3594 += (-3*Conj(ZN(gI2,4))) * tmp_3599;
   result += (0.3333333333333333) * tmp_3594;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   const auto Qq = LOCALINPUT(Qq);

   std::complex<double> result;

   std::complex<double> tmp_3603;
   std::complex<double> tmp_3604;
   std::complex<double> tmp_3605;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3605 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_3604 += tmp_3605;
   tmp_3603 += (-1.4142135623730951*(6*gp*Qq*ZN(gI2,0) + 0.7745966692414834*g1*
      ZN(gI2,1) + 3*g2*ZN(gI2,2))) * tmp_3604;
   std::complex<double> tmp_3606;
   std::complex<double> tmp_3607;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_3608;
      std::complex<double> tmp_3609;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_3609 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_3608 += tmp_3609;
      tmp_3607 += (ZUL(gO1,j2)) * tmp_3608;
   }
   tmp_3606 += tmp_3607;
   tmp_3603 += (-6*ZN(gI2,4)) * tmp_3606;
   result += (0.16666666666666666) * tmp_3603;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3610;
   std::complex<double> tmp_3611;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3611 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_3610 += tmp_3611;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_3610;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3612;
   std::complex<double> tmp_3613;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3613 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_3612 += tmp_3613;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_3612;

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

   std::complex<double> tmp_3614;
   std::complex<double> tmp_3615;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3615 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_3614 += tmp_3615;
   result += (-0.7071067811865475*g2) * tmp_3614;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSdconjUSdVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_3616;
   std::complex<double> tmp_3617;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3617 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3616 += tmp_3617;
   result += (-1) * tmp_3616;
   std::complex<double> tmp_3618;
   std::complex<double> tmp_3619;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3619 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3618 += tmp_3619;
   result += (-0.5) * tmp_3618;
   std::complex<double> tmp_3620;
   std::complex<double> tmp_3621;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3621 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3620 += tmp_3621;
   result += (-0.5) * tmp_3620;
   std::complex<double> tmp_3622;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3623;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3623 += (Conj(CpconjUSdFuChaPL(gO2,gI1,gI2))*
            CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPR(gO2,gI1,gI2))*
            CpconjUSdFuChaPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha(gI2));
      }
      tmp_3622 += tmp_3623;
   }
   result += tmp_3622;
   std::complex<double> tmp_3624;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3625;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3625 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_3624 += tmp_3625;
   }
   result += tmp_3624;
   std::complex<double> tmp_3626;
   std::complex<double> tmp_3627;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3628;
      std::complex<double> tmp_3629;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3629 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_3628 += tmp_3629;
      tmp_3627 += (MFd(gI1)) * tmp_3628;
   }
   tmp_3626 += tmp_3627;
   result += (-2) * tmp_3626;
   std::complex<double> tmp_3630;
   std::complex<double> tmp_3631;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3632;
      std::complex<double> tmp_3633;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3633 += B0(p,MFu(gI1),MCha(gI2))*(Conj(CpconjUSdFuChaPR(gO2,
            gI1,gI2))*CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPL(gO2,
            gI1,gI2))*CpconjUSdFuChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_3632 += tmp_3633;
      tmp_3631 += (MFu(gI1)) * tmp_3632;
   }
   tmp_3630 += tmp_3631;
   result += (-2) * tmp_3630;
   std::complex<double> tmp_3634;
   std::complex<double> tmp_3635;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3635 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3634 += tmp_3635;
   result += (-1) * tmp_3634;
   std::complex<double> tmp_3636;
   std::complex<double> tmp_3637;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3637 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3636 += tmp_3637;
   result += (-1) * tmp_3636;
   std::complex<double> tmp_3638;
   std::complex<double> tmp_3639;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3639 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3638 += tmp_3639;
   result += (-1) * tmp_3638;
   std::complex<double> tmp_3640;
   std::complex<double> tmp_3641;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3641 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3640 += tmp_3641;
   result += (-1) * tmp_3640;
   std::complex<double> tmp_3642;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3643;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3643 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_3642 += tmp_3643;
   }
   result += tmp_3642;
   std::complex<double> tmp_3644;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3645;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3645 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_3644 += tmp_3645;
   }
   result += tmp_3644;
   std::complex<double> tmp_3646;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3647;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3647 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_3646 += tmp_3647;
   }
   result += tmp_3646;
   std::complex<double> tmp_3648;
   std::complex<double> tmp_3649;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3649 += (Conj(CpconjUSdGluFdPL(gO2,1,gI2))*CpconjUSdGluFdPL(gO1,1,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,1,gI2))*CpconjUSdGluFdPR(gO1,1,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_3648 += tmp_3649;
   result += (1.3333333333333333) * tmp_3648;
   std::complex<double> tmp_3650;
   std::complex<double> tmp_3651;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3651 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_3650 += tmp_3651;
   result += (1.3333333333333333) * tmp_3650;
   std::complex<double> tmp_3652;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3652 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_3652;
   std::complex<double> tmp_3653;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3653 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_3653;
   std::complex<double> tmp_3654;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3654 += Conj(CpconjUSdVZpSd(gO2,gI2))*CpconjUSdVZpSd(gO1,gI2)*F0(p
         ,MSd(gI2),MVZp);
   }
   result += tmp_3654;
   std::complex<double> tmp_3655;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3655 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_3655;
   std::complex<double> tmp_3656;
   std::complex<double> tmp_3657;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3657 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,1,gI2))*
         CpconjUSdGluFdPL(gO1,1,gI2) + Conj(CpconjUSdGluFdPL(gO2,1,gI2))*
         CpconjUSdGluFdPR(gO1,1,gI2))*MFd(gI2);
   }
   tmp_3656 += tmp_3657;
   result += (-2.6666666666666665*MGlu) * tmp_3656;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSvconjUSvVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_3658;
   std::complex<double> tmp_3659;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3659 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3658 += tmp_3659;
   result += (-1) * tmp_3658;
   std::complex<double> tmp_3660;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3661;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3661 += (Conj(CpconjUSvbarChaFePL(gO2,gI1,gI2))*
            CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(CpconjUSvbarChaFePR(gO2,gI1,
            gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFe(gI2));
      }
      tmp_3660 += tmp_3661;
   }
   result += tmp_3660;
   std::complex<double> tmp_3662;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3663;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3663 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_3662 += tmp_3663;
   }
   result += tmp_3662;
   std::complex<double> tmp_3664;
   std::complex<double> tmp_3665;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3666;
      std::complex<double> tmp_3667;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3667 += B0(p,MCha(gI1),MFe(gI2))*(Conj(CpconjUSvbarChaFePR(
            gO2,gI1,gI2))*CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(
            CpconjUSvbarChaFePL(gO2,gI1,gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_3666 += tmp_3667;
      tmp_3665 += (MCha(gI1)) * tmp_3666;
   }
   tmp_3664 += tmp_3665;
   result += (-2) * tmp_3664;
   std::complex<double> tmp_3668;
   std::complex<double> tmp_3669;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3669 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3668 += tmp_3669;
   result += (-0.5) * tmp_3668;
   std::complex<double> tmp_3670;
   std::complex<double> tmp_3671;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3671 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3670 += tmp_3671;
   result += (-0.5) * tmp_3670;
   std::complex<double> tmp_3672;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3673;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3673 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_3672 += tmp_3673;
   }
   result += tmp_3672;
   std::complex<double> tmp_3674;
   std::complex<double> tmp_3675;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3676;
      std::complex<double> tmp_3677;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3677 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_3676 += tmp_3677;
      tmp_3675 += (MFv(gI1)) * tmp_3676;
   }
   tmp_3674 += tmp_3675;
   result += (-2) * tmp_3674;
   std::complex<double> tmp_3678;
   std::complex<double> tmp_3679;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3679 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3678 += tmp_3679;
   result += (-3) * tmp_3678;
   std::complex<double> tmp_3680;
   std::complex<double> tmp_3681;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3681 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3680 += tmp_3681;
   result += (-1) * tmp_3680;
   std::complex<double> tmp_3682;
   std::complex<double> tmp_3683;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3683 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3682 += tmp_3683;
   result += (-3) * tmp_3682;
   std::complex<double> tmp_3684;
   std::complex<double> tmp_3685;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3685 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3684 += tmp_3685;
   result += (-1) * tmp_3684;
   std::complex<double> tmp_3686;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3687;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3687 += B0(p,MSv(gI1),MAh(gI2))*Conj(CpconjUSvSvAh(gO2,gI1,
            gI2))*CpconjUSvSvAh(gO1,gI1,gI2);
      }
      tmp_3686 += tmp_3687;
   }
   result += tmp_3686;
   std::complex<double> tmp_3688;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3689;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3689 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_3688 += tmp_3689;
   }
   result += tmp_3688;
   std::complex<double> tmp_3690;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3690 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_3690;
   std::complex<double> tmp_3691;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3691 += Conj(CpconjUSvVPSv(gO2,gI2))*CpconjUSvVPSv(gO1,gI2)*F0(p,
         MSv(gI2),0);
   }
   result += tmp_3691;
   std::complex<double> tmp_3692;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3692 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_3692;
   std::complex<double> tmp_3693;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3693 += Conj(CpconjUSvVZpSv(gO2,gI2))*CpconjUSvVZpSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVZp);
   }
   result += tmp_3693;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSuconjUSuVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_3694;
   std::complex<double> tmp_3695;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3695 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3694 += tmp_3695;
   result += (-1) * tmp_3694;
   std::complex<double> tmp_3696;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3697;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3697 += (Conj(CpconjUSubarChaFdPL(gO2,gI1,gI2))*
            CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(CpconjUSubarChaFdPR(gO2,gI1,
            gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFd(gI2));
      }
      tmp_3696 += tmp_3697;
   }
   result += tmp_3696;
   std::complex<double> tmp_3698;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3699;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3699 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_3698 += tmp_3699;
   }
   result += tmp_3698;
   std::complex<double> tmp_3700;
   std::complex<double> tmp_3701;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3702;
      std::complex<double> tmp_3703;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3703 += B0(p,MCha(gI1),MFd(gI2))*(Conj(CpconjUSubarChaFdPR(
            gO2,gI1,gI2))*CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChaFdPL(gO2,gI1,gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_3702 += tmp_3703;
      tmp_3701 += (MCha(gI1)) * tmp_3702;
   }
   tmp_3700 += tmp_3701;
   result += (-2) * tmp_3700;
   std::complex<double> tmp_3704;
   std::complex<double> tmp_3705;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3705 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3704 += tmp_3705;
   result += (-0.5) * tmp_3704;
   std::complex<double> tmp_3706;
   std::complex<double> tmp_3707;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3707 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3706 += tmp_3707;
   result += (-0.5) * tmp_3706;
   std::complex<double> tmp_3708;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3709;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3709 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_3708 += tmp_3709;
   }
   result += tmp_3708;
   std::complex<double> tmp_3710;
   std::complex<double> tmp_3711;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3712;
      std::complex<double> tmp_3713;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3713 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_3712 += tmp_3713;
      tmp_3711 += (MFu(gI1)) * tmp_3712;
   }
   tmp_3710 += tmp_3711;
   result += (-2) * tmp_3710;
   std::complex<double> tmp_3714;
   std::complex<double> tmp_3715;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3715 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3714 += tmp_3715;
   result += (-1) * tmp_3714;
   std::complex<double> tmp_3716;
   std::complex<double> tmp_3717;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3717 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3716 += tmp_3717;
   result += (-1) * tmp_3716;
   std::complex<double> tmp_3718;
   std::complex<double> tmp_3719;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3719 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3718 += tmp_3719;
   result += (-1) * tmp_3718;
   std::complex<double> tmp_3720;
   std::complex<double> tmp_3721;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3721 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3720 += tmp_3721;
   result += (-1) * tmp_3720;
   std::complex<double> tmp_3722;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3723;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3723 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_3722 += tmp_3723;
   }
   result += tmp_3722;
   std::complex<double> tmp_3724;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3725;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3725 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_3724 += tmp_3725;
   }
   result += tmp_3724;
   std::complex<double> tmp_3726;
   std::complex<double> tmp_3727;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3727 += (Conj(CpconjUSuGluFuPL(gO2,1,gI2))*CpconjUSuGluFuPL(gO1,1,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,1,gI2))*CpconjUSuGluFuPR(gO1,1,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_3726 += tmp_3727;
   result += (1.3333333333333333) * tmp_3726;
   std::complex<double> tmp_3728;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3728 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_3728;
   std::complex<double> tmp_3729;
   std::complex<double> tmp_3730;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3730 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_3729 += tmp_3730;
   result += (1.3333333333333333) * tmp_3729;
   std::complex<double> tmp_3731;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3731 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_3731;
   std::complex<double> tmp_3732;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3732 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_3732;
   std::complex<double> tmp_3733;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3733 += Conj(CpconjUSuVZpSu(gO2,gI2))*CpconjUSuVZpSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVZp);
   }
   result += tmp_3733;
   std::complex<double> tmp_3734;
   std::complex<double> tmp_3735;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3735 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,1,gI2))*
         CpconjUSuGluFuPL(gO1,1,gI2) + Conj(CpconjUSuGluFuPL(gO2,1,gI2))*
         CpconjUSuGluFuPR(gO1,1,gI2))*MFu(gI2);
   }
   tmp_3734 += tmp_3735;
   result += (-2.6666666666666665*MGlu) * tmp_3734;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZp)*CpUSeconjUSeVZpVZp(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_3736;
   std::complex<double> tmp_3737;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3737 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3736 += tmp_3737;
   result += (-1) * tmp_3736;
   std::complex<double> tmp_3738;
   std::complex<double> tmp_3739;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3739 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3738 += tmp_3739;
   result += (-0.5) * tmp_3738;
   std::complex<double> tmp_3740;
   std::complex<double> tmp_3741;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3741 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3740 += tmp_3741;
   result += (-0.5) * tmp_3740;
   std::complex<double> tmp_3742;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3743;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3743 += (Conj(CpconjUSeFvChaPL(gO2,gI1,gI2))*
            CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPR(gO2,gI1,gI2))*
            CpconjUSeFvChaPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha(gI2));
      }
      tmp_3742 += tmp_3743;
   }
   result += tmp_3742;
   std::complex<double> tmp_3744;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3745;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3745 += (Conj(CpconjUSeFeChiPL(gO2,gI1,gI2))*
            CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPR(gO2,gI1,gI2))*
            CpconjUSeFeChiPR(gO1,gI1,gI2))*G0(p,MFe(gI1),MChi(gI2));
      }
      tmp_3744 += tmp_3745;
   }
   result += tmp_3744;
   std::complex<double> tmp_3746;
   std::complex<double> tmp_3747;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3748;
      std::complex<double> tmp_3749;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3749 += B0(p,MFe(gI1),MChi(gI2))*(Conj(CpconjUSeFeChiPR(gO2,
            gI1,gI2))*CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPL(gO2,
            gI1,gI2))*CpconjUSeFeChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_3748 += tmp_3749;
      tmp_3747 += (MFe(gI1)) * tmp_3748;
   }
   tmp_3746 += tmp_3747;
   result += (-2) * tmp_3746;
   std::complex<double> tmp_3750;
   std::complex<double> tmp_3751;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3752;
      std::complex<double> tmp_3753;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3753 += B0(p,MFv(gI1),MCha(gI2))*(Conj(CpconjUSeFvChaPR(gO2,
            gI1,gI2))*CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPL(gO2,
            gI1,gI2))*CpconjUSeFvChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_3752 += tmp_3753;
      tmp_3751 += (MFv(gI1)) * tmp_3752;
   }
   tmp_3750 += tmp_3751;
   result += (-2) * tmp_3750;
   std::complex<double> tmp_3754;
   std::complex<double> tmp_3755;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3755 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3754 += tmp_3755;
   result += (-3) * tmp_3754;
   std::complex<double> tmp_3756;
   std::complex<double> tmp_3757;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3757 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3756 += tmp_3757;
   result += (-1) * tmp_3756;
   std::complex<double> tmp_3758;
   std::complex<double> tmp_3759;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3759 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3758 += tmp_3759;
   result += (-3) * tmp_3758;
   std::complex<double> tmp_3760;
   std::complex<double> tmp_3761;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3761 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3760 += tmp_3761;
   result += (-1) * tmp_3760;
   std::complex<double> tmp_3762;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3763;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3763 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeSvHpm(gO2,gI1
            ,gI2))*CpconjUSeSvHpm(gO1,gI1,gI2);
      }
      tmp_3762 += tmp_3763;
   }
   result += tmp_3762;
   std::complex<double> tmp_3764;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3765;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3765 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_3764 += tmp_3765;
   }
   result += tmp_3764;
   std::complex<double> tmp_3766;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3767;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3767 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_3766 += tmp_3767;
   }
   result += tmp_3766;
   std::complex<double> tmp_3768;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3768 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_3768;
   std::complex<double> tmp_3769;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3769 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_3769;
   std::complex<double> tmp_3770;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3770 += Conj(CpconjUSeVZpSe(gO2,gI2))*CpconjUSeVZpSe(gO1,gI2)*F0(p
         ,MSe(gI2),MVZp);
   }
   result += tmp_3770;
   std::complex<double> tmp_3771;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3771 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_3771;

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
   std::complex<double> tmp_3772;
   std::complex<double> tmp_3773;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3773 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3772 += tmp_3773;
   result += (-1) * tmp_3772;
   std::complex<double> tmp_3774;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3775;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3775 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_3774 += tmp_3775;
   }
   result += tmp_3774;
   std::complex<double> tmp_3776;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3777;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3777 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*
            CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*
            CpUhhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_3776 += tmp_3777;
   }
   result += tmp_3776;
   std::complex<double> tmp_3778;
   std::complex<double> tmp_3779;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3780;
      std::complex<double> tmp_3781;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3781 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2
            ,gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_3780 += tmp_3781;
      tmp_3779 += (MCha(gI1)) * tmp_3780;
   }
   tmp_3778 += tmp_3779;
   result += (-2) * tmp_3778;
   std::complex<double> tmp_3782;
   std::complex<double> tmp_3783;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3783 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3782 += tmp_3783;
   result += (-0.5) * tmp_3782;
   std::complex<double> tmp_3784;
   std::complex<double> tmp_3785;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3785 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3784 += tmp_3785;
   result += (-0.5) * tmp_3784;
   std::complex<double> tmp_3786;
   std::complex<double> tmp_3787;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3788;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3788 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_3787 += tmp_3788;
   }
   tmp_3786 += tmp_3787;
   result += (0.5) * tmp_3786;
   std::complex<double> tmp_3789;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3790;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3790 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))
            *CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_3789 += tmp_3790;
   }
   result += tmp_3789;
   std::complex<double> tmp_3791;
   std::complex<double> tmp_3792;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3793;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3793 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_3792 += tmp_3793;
   }
   tmp_3791 += tmp_3792;
   result += (0.5) * tmp_3791;
   std::complex<double> tmp_3794;
   std::complex<double> tmp_3795;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3796;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3796 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_3795 += tmp_3796;
   }
   tmp_3794 += tmp_3795;
   result += (3) * tmp_3794;
   std::complex<double> tmp_3797;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3798;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3798 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_3797 += tmp_3798;
   }
   result += tmp_3797;
   std::complex<double> tmp_3799;
   std::complex<double> tmp_3800;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3801;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3801 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_3800 += tmp_3801;
   }
   tmp_3799 += tmp_3800;
   result += (3) * tmp_3799;
   std::complex<double> tmp_3802;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3803;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3803 += (Conj(CpUhhbarFvFvPL(gO2,gI1,gI2))*CpUhhbarFvFvPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFvFvPR(gO2,gI1,gI2))*CpUhhbarFvFvPR(gO1,
            gI1,gI2))*G0(p,MFv(gI1),MFv(gI2));
      }
      tmp_3802 += tmp_3803;
   }
   result += tmp_3802;
   std::complex<double> tmp_3804;
   std::complex<double> tmp_3805;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3806;
      std::complex<double> tmp_3807;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3807 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_3806 += tmp_3807;
      tmp_3805 += (MFd(gI1)) * tmp_3806;
   }
   tmp_3804 += tmp_3805;
   result += (-6) * tmp_3804;
   std::complex<double> tmp_3808;
   std::complex<double> tmp_3809;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3810;
      std::complex<double> tmp_3811;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3811 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_3810 += tmp_3811;
      tmp_3809 += (MFe(gI1)) * tmp_3810;
   }
   tmp_3808 += tmp_3809;
   result += (-2) * tmp_3808;
   std::complex<double> tmp_3812;
   std::complex<double> tmp_3813;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3814;
      std::complex<double> tmp_3815;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3815 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_3814 += tmp_3815;
      tmp_3813 += (MFu(gI1)) * tmp_3814;
   }
   tmp_3812 += tmp_3813;
   result += (-6) * tmp_3812;
   std::complex<double> tmp_3816;
   std::complex<double> tmp_3817;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3818;
      std::complex<double> tmp_3819;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3819 += B0(p,MFv(gI1),MFv(gI2))*(Conj(CpUhhbarFvFvPR(gO2,gI1
            ,gI2))*CpUhhbarFvFvPL(gO1,gI1,gI2) + Conj(CpUhhbarFvFvPL(gO2,gI1,gI2))
            *CpUhhbarFvFvPR(gO1,gI1,gI2))*MFv(gI2);
      }
      tmp_3818 += tmp_3819;
      tmp_3817 += (MFv(gI1)) * tmp_3818;
   }
   tmp_3816 += tmp_3817;
   result += (-2) * tmp_3816;
   std::complex<double> tmp_3820;
   std::complex<double> tmp_3821;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3821 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3820 += tmp_3821;
   result += (-3) * tmp_3820;
   std::complex<double> tmp_3822;
   std::complex<double> tmp_3823;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3823 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3822 += tmp_3823;
   result += (-1) * tmp_3822;
   std::complex<double> tmp_3824;
   std::complex<double> tmp_3825;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3825 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3824 += tmp_3825;
   result += (-3) * tmp_3824;
   std::complex<double> tmp_3826;
   std::complex<double> tmp_3827;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3827 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3826 += tmp_3827;
   result += (-1) * tmp_3826;
   std::complex<double> tmp_3828;
   std::complex<double> tmp_3829;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3830;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3830 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_3829 += tmp_3830;
   }
   tmp_3828 += tmp_3829;
   result += (3) * tmp_3828;
   std::complex<double> tmp_3831;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3832;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3832 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_3831 += tmp_3832;
   }
   result += tmp_3831;
   std::complex<double> tmp_3833;
   std::complex<double> tmp_3834;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3835;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3835 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_3834 += tmp_3835;
   }
   tmp_3833 += tmp_3834;
   result += (3) * tmp_3833;
   std::complex<double> tmp_3836;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3837;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3837 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_3836 += tmp_3837;
   }
   result += tmp_3836;
   std::complex<double> tmp_3838;
   std::complex<double> tmp_3839;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3840;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3840 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_3839 += tmp_3840;
   }
   tmp_3838 += tmp_3839;
   result += (0.5) * tmp_3838;
   std::complex<double> tmp_3841;
   std::complex<double> tmp_3842;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3843;
      std::complex<double> tmp_3844;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3844 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,
            gI1,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2)
            )*CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_3843 += tmp_3844;
      tmp_3842 += (MChi(gI1)) * tmp_3843;
   }
   tmp_3841 += tmp_3842;
   result += (-1) * tmp_3841;
   std::complex<double> tmp_3845;
   std::complex<double> tmp_3846;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3846 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_3845 += tmp_3846;
   result += (2) * tmp_3845;
   std::complex<double> tmp_3847;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3847 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_3847;
   std::complex<double> tmp_3848;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3848 += Conj(CpUhhVZpAh(gO2,gI2))*CpUhhVZpAh(gO1,gI2)*F0(p,MAh(gI2
         ),MVZp);
   }
   result += tmp_3848;

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
   std::complex<double> tmp_3849;
   std::complex<double> tmp_3850;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3850 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3849 += tmp_3850;
   result += (-1) * tmp_3849;
   std::complex<double> tmp_3851;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3852;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3852 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_3851 += tmp_3852;
   }
   result += tmp_3851;
   std::complex<double> tmp_3853;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3854;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3854 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*
            CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*
            CpUAhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_3853 += tmp_3854;
   }
   result += tmp_3853;
   std::complex<double> tmp_3855;
   std::complex<double> tmp_3856;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3857;
      std::complex<double> tmp_3858;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3858 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2
            ,gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_3857 += tmp_3858;
      tmp_3856 += (MCha(gI1)) * tmp_3857;
   }
   tmp_3855 += tmp_3856;
   result += (-2) * tmp_3855;
   std::complex<double> tmp_3859;
   std::complex<double> tmp_3860;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3860 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3859 += tmp_3860;
   result += (-0.5) * tmp_3859;
   std::complex<double> tmp_3861;
   std::complex<double> tmp_3862;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3862 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3861 += tmp_3862;
   result += (-0.5) * tmp_3861;
   std::complex<double> tmp_3863;
   std::complex<double> tmp_3864;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3865;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3865 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))
            *CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_3864 += tmp_3865;
   }
   tmp_3863 += tmp_3864;
   result += (0.5) * tmp_3863;
   std::complex<double> tmp_3866;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3867;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3867 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_3866 += tmp_3867;
   }
   result += tmp_3866;
   std::complex<double> tmp_3868;
   std::complex<double> tmp_3869;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3870;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3870 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))
            *CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_3869 += tmp_3870;
   }
   tmp_3868 += tmp_3869;
   result += (0.5) * tmp_3868;
   std::complex<double> tmp_3871;
   std::complex<double> tmp_3872;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3873;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3873 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_3872 += tmp_3873;
   }
   tmp_3871 += tmp_3872;
   result += (3) * tmp_3871;
   std::complex<double> tmp_3874;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3875;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3875 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_3874 += tmp_3875;
   }
   result += tmp_3874;
   std::complex<double> tmp_3876;
   std::complex<double> tmp_3877;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3878;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3878 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_3877 += tmp_3878;
   }
   tmp_3876 += tmp_3877;
   result += (3) * tmp_3876;
   std::complex<double> tmp_3879;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3880;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3880 += (Conj(CpUAhbarFvFvPL(gO2,gI1,gI2))*CpUAhbarFvFvPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFvFvPR(gO2,gI1,gI2))*CpUAhbarFvFvPR(gO1,
            gI1,gI2))*G0(p,MFv(gI1),MFv(gI2));
      }
      tmp_3879 += tmp_3880;
   }
   result += tmp_3879;
   std::complex<double> tmp_3881;
   std::complex<double> tmp_3882;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3883;
      std::complex<double> tmp_3884;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3884 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_3883 += tmp_3884;
      tmp_3882 += (MFd(gI1)) * tmp_3883;
   }
   tmp_3881 += tmp_3882;
   result += (-6) * tmp_3881;
   std::complex<double> tmp_3885;
   std::complex<double> tmp_3886;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3887;
      std::complex<double> tmp_3888;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3888 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_3887 += tmp_3888;
      tmp_3886 += (MFe(gI1)) * tmp_3887;
   }
   tmp_3885 += tmp_3886;
   result += (-2) * tmp_3885;
   std::complex<double> tmp_3889;
   std::complex<double> tmp_3890;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3891;
      std::complex<double> tmp_3892;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3892 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_3891 += tmp_3892;
      tmp_3890 += (MFu(gI1)) * tmp_3891;
   }
   tmp_3889 += tmp_3890;
   result += (-6) * tmp_3889;
   std::complex<double> tmp_3893;
   std::complex<double> tmp_3894;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3895;
      std::complex<double> tmp_3896;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3896 += B0(p,MFv(gI1),MFv(gI2))*(Conj(CpUAhbarFvFvPR(gO2,gI1
            ,gI2))*CpUAhbarFvFvPL(gO1,gI1,gI2) + Conj(CpUAhbarFvFvPL(gO2,gI1,gI2))
            *CpUAhbarFvFvPR(gO1,gI1,gI2))*MFv(gI2);
      }
      tmp_3895 += tmp_3896;
      tmp_3894 += (MFv(gI1)) * tmp_3895;
   }
   tmp_3893 += tmp_3894;
   result += (-2) * tmp_3893;
   std::complex<double> tmp_3897;
   std::complex<double> tmp_3898;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3898 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3897 += tmp_3898;
   result += (-3) * tmp_3897;
   std::complex<double> tmp_3899;
   std::complex<double> tmp_3900;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3900 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3899 += tmp_3900;
   result += (-1) * tmp_3899;
   std::complex<double> tmp_3901;
   std::complex<double> tmp_3902;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3902 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3901 += tmp_3902;
   result += (-3) * tmp_3901;
   std::complex<double> tmp_3903;
   std::complex<double> tmp_3904;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3904 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3903 += tmp_3904;
   result += (-1) * tmp_3903;
   std::complex<double> tmp_3905;
   std::complex<double> tmp_3906;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3907;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3907 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_3906 += tmp_3907;
   }
   tmp_3905 += tmp_3906;
   result += (3) * tmp_3905;
   std::complex<double> tmp_3908;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3909;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3909 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_3908 += tmp_3909;
   }
   result += tmp_3908;
   std::complex<double> tmp_3910;
   std::complex<double> tmp_3911;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3912;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3912 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_3911 += tmp_3912;
   }
   tmp_3910 += tmp_3911;
   result += (3) * tmp_3910;
   std::complex<double> tmp_3913;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3914;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3914 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUAhconjSvSv(gO2,gI1,
            gI2))*CpUAhconjSvSv(gO1,gI1,gI2);
      }
      tmp_3913 += tmp_3914;
   }
   result += tmp_3913;
   std::complex<double> tmp_3915;
   std::complex<double> tmp_3916;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3917;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3917 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_3916 += tmp_3917;
   }
   tmp_3915 += tmp_3916;
   result += (0.5) * tmp_3915;
   std::complex<double> tmp_3918;
   std::complex<double> tmp_3919;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3920;
      std::complex<double> tmp_3921;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3921 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,
            gI1,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2)
            )*CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_3920 += tmp_3921;
      tmp_3919 += (MChi(gI1)) * tmp_3920;
   }
   tmp_3918 += tmp_3919;
   result += (-1) * tmp_3918;
   std::complex<double> tmp_3922;
   std::complex<double> tmp_3923;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3923 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_3922 += tmp_3923;
   result += (2) * tmp_3922;
   std::complex<double> tmp_3924;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3924 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_3924;
   std::complex<double> tmp_3925;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3925 += Conj(CpUAhVZphh(gO2,gI2))*CpUAhVZphh(gO1,gI2)*F0(p,Mhh(gI2
         ),MVZp);
   }
   result += tmp_3925;

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
   std::complex<double> tmp_3926;
   std::complex<double> tmp_3927;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3927 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_3926 += tmp_3927;
   result += (-1) * tmp_3926;
   std::complex<double> tmp_3928;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3929;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3929 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_3928 += tmp_3929;
   }
   result += tmp_3928;
   std::complex<double> tmp_3930;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3931;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3931 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_3930 += tmp_3931;
   }
   result += tmp_3930;
   std::complex<double> tmp_3932;
   std::complex<double> tmp_3933;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3933 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_3932 += tmp_3933;
   result += (-0.5) * tmp_3932;
   std::complex<double> tmp_3934;
   std::complex<double> tmp_3935;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3935 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_3934 += tmp_3935;
   result += (-0.5) * tmp_3934;
   std::complex<double> tmp_3936;
   std::complex<double> tmp_3937;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3938;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3938 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_3937 += tmp_3938;
   }
   tmp_3936 += tmp_3937;
   result += (3) * tmp_3936;
   std::complex<double> tmp_3939;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3940;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3940 += (Conj(CpconjUHpmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFvFePR(gO2,gI1,
            gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_3939 += tmp_3940;
   }
   result += tmp_3939;
   std::complex<double> tmp_3941;
   std::complex<double> tmp_3942;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3943;
      std::complex<double> tmp_3944;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3944 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_3943 += tmp_3944;
      tmp_3942 += (MFu(gI1)) * tmp_3943;
   }
   tmp_3941 += tmp_3942;
   result += (-6) * tmp_3941;
   std::complex<double> tmp_3945;
   std::complex<double> tmp_3946;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3947;
      std::complex<double> tmp_3948;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3948 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmbarFvFePR(
            gO2,gI1,gI2))*CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFvFePL(gO2,gI1,gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_3947 += tmp_3948;
      tmp_3946 += (MFv(gI1)) * tmp_3947;
   }
   tmp_3945 += tmp_3946;
   result += (-2) * tmp_3945;
   std::complex<double> tmp_3949;
   std::complex<double> tmp_3950;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3950 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_3949 += tmp_3950;
   result += (-3) * tmp_3949;
   std::complex<double> tmp_3951;
   std::complex<double> tmp_3952;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3952 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_3951 += tmp_3952;
   result += (-1) * tmp_3951;
   std::complex<double> tmp_3953;
   std::complex<double> tmp_3954;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3954 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_3953 += tmp_3954;
   result += (-3) * tmp_3953;
   std::complex<double> tmp_3955;
   std::complex<double> tmp_3956;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3956 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_3955 += tmp_3956;
   result += (-1) * tmp_3955;
   std::complex<double> tmp_3957;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3958;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3958 += (Conj(CpconjUHpmChiChaPL(gO2,gI1,gI2))*
            CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHpmChiChaPR(gO2,gI1,gI2)
            )*CpconjUHpmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_3957 += tmp_3958;
   }
   result += tmp_3957;
   std::complex<double> tmp_3959;
   std::complex<double> tmp_3960;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3961;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3961 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_3960 += tmp_3961;
   }
   tmp_3959 += tmp_3960;
   result += (3) * tmp_3959;
   std::complex<double> tmp_3962;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3963;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3963 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_3962 += tmp_3963;
   }
   result += tmp_3962;
   std::complex<double> tmp_3964;
   std::complex<double> tmp_3965;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3966;
      std::complex<double> tmp_3967;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3967 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHpmChiChaPR(
            gO2,gI1,gI2))*CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmChiChaPL(gO2,gI1,gI2))*CpconjUHpmChiChaPR(gO1,gI1,gI2))*MCha
            (gI2);
      }
      tmp_3966 += tmp_3967;
      tmp_3965 += (MChi(gI1)) * tmp_3966;
   }
   tmp_3964 += tmp_3965;
   result += (-2) * tmp_3964;
   std::complex<double> tmp_3968;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3968 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_3968;
   std::complex<double> tmp_3969;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3969 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_3969;
   std::complex<double> tmp_3970;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3970 += Conj(CpconjUHpmVZpHpm(gO2,gI2))*CpconjUHpmVZpHpm(gO1,gI2)*
         F0(p,MHpm(gI2),MVZp);
   }
   result += tmp_3970;
   std::complex<double> tmp_3971;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3971 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_3971;
   std::complex<double> tmp_3972;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3972 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_3972;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_3973;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3973 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_3973;
   std::complex<double> tmp_3974;
   std::complex<double> tmp_3975;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3976;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3976 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_3975 += tmp_3976;
   }
   tmp_3974 += tmp_3975;
   result += (-4) * tmp_3974;
   std::complex<double> tmp_3977;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3978;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3978 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_3978 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_3977 += tmp_3978;
   }
   result += tmp_3977;
   std::complex<double> tmp_3979;
   std::complex<double> tmp_3980;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3980 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_3979 += tmp_3980;
   result += (0.5) * tmp_3979;
   std::complex<double> tmp_3981;
   std::complex<double> tmp_3982;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3982 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_3981 += tmp_3982;
   result += (0.5) * tmp_3981;
   std::complex<double> tmp_3983;
   std::complex<double> tmp_3984;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3985;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3985 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_3984 += tmp_3985;
   }
   tmp_3983 += tmp_3984;
   result += (-4) * tmp_3983;
   std::complex<double> tmp_3986;
   std::complex<double> tmp_3987;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3988;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3988 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_3988 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_3987 += tmp_3988;
   }
   tmp_3986 += tmp_3987;
   result += (3) * tmp_3986;
   std::complex<double> tmp_3989;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3990;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3990 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_3990 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_3989 += tmp_3990;
   }
   result += tmp_3989;
   std::complex<double> tmp_3991;
   std::complex<double> tmp_3992;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3993;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3993 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_3993 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_3992 += tmp_3993;
   }
   tmp_3991 += tmp_3992;
   result += (3) * tmp_3991;
   std::complex<double> tmp_3994;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3995;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3995 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_3995 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_3994 += tmp_3995;
   }
   result += tmp_3994;
   std::complex<double> tmp_3996;
   std::complex<double> tmp_3997;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3997 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_3996 += tmp_3997;
   result += (3) * tmp_3996;
   std::complex<double> tmp_3998;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3998 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_3998;
   std::complex<double> tmp_3999;
   std::complex<double> tmp_4000;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4000 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_3999 += tmp_4000;
   result += (3) * tmp_3999;
   std::complex<double> tmp_4001;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4001 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_4001;
   std::complex<double> tmp_4002;
   std::complex<double> tmp_4003;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4004;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4004 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_4003 += tmp_4004;
   }
   tmp_4002 += tmp_4003;
   result += (-12) * tmp_4002;
   std::complex<double> tmp_4005;
   std::complex<double> tmp_4006;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4007;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4007 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_4006 += tmp_4007;
   }
   tmp_4005 += tmp_4006;
   result += (-4) * tmp_4005;
   std::complex<double> tmp_4008;
   std::complex<double> tmp_4009;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4010;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4010 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_4009 += tmp_4010;
   }
   tmp_4008 += tmp_4009;
   result += (-12) * tmp_4008;
   std::complex<double> tmp_4011;
   std::complex<double> tmp_4012;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4013;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4013 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_4012 += tmp_4013;
   }
   tmp_4011 += tmp_4012;
   result += (-4) * tmp_4011;
   std::complex<double> tmp_4014;
   std::complex<double> tmp_4015;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4016;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4016 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_4016 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_4015 += tmp_4016;
   }
   tmp_4014 += tmp_4015;
   result += (0.5) * tmp_4014;
   std::complex<double> tmp_4017;
   std::complex<double> tmp_4018;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4018 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_4017 += tmp_4018;
   result += (2) * tmp_4017;
   std::complex<double> tmp_4019;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4019 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_4019;
   std::complex<double> tmp_4020;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4020 += AbsSqr(CpVZVZphh(gI2))*B0(p,MVZp,Mhh(gI2));
   }
   result += tmp_4020;
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
   std::complex<double> tmp_4021;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_4021 += A0(MHpm(gI1))*CpVZpVZpconjHpmHpm(gI1,gI1);
   }
   result += tmp_4021;
   std::complex<double> tmp_4022;
   std::complex<double> tmp_4023;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4024;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4024 += AbsSqr(CpVZpconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),
            MHpm(gI2));
      }
      tmp_4023 += tmp_4024;
   }
   tmp_4022 += tmp_4023;
   result += (-4) * tmp_4022;
   std::complex<double> tmp_4025;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4026;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4026 += (AbsSqr(CpVZpbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZpbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_4026 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZpbarChaChaPL(gI1,gI2))*CpVZpbarChaChaPR(gI1,gI2));
      }
      tmp_4025 += tmp_4026;
   }
   result += tmp_4025;
   std::complex<double> tmp_4027;
   std::complex<double> tmp_4028;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4028 += A0(MAh(gI1))*CpVZpVZpAhAh(gI1,gI1);
   }
   tmp_4027 += tmp_4028;
   result += (0.5) * tmp_4027;
   std::complex<double> tmp_4029;
   std::complex<double> tmp_4030;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4030 += A0(Mhh(gI1))*CpVZpVZphhhh(gI1,gI1);
   }
   tmp_4029 += tmp_4030;
   result += (0.5) * tmp_4029;
   std::complex<double> tmp_4031;
   std::complex<double> tmp_4032;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4033;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4033 += AbsSqr(CpVZphhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_4032 += tmp_4033;
   }
   tmp_4031 += tmp_4032;
   result += (-4) * tmp_4031;
   std::complex<double> tmp_4034;
   std::complex<double> tmp_4035;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4036;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4036 += (AbsSqr(CpVZpbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZpbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_4036 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZpbarFdFdPL(gI1,gI2))*CpVZpbarFdFdPR(gI1,gI2));
      }
      tmp_4035 += tmp_4036;
   }
   tmp_4034 += tmp_4035;
   result += (3) * tmp_4034;
   std::complex<double> tmp_4037;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4038;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4038 += (AbsSqr(CpVZpbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZpbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_4038 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZpbarFeFePL(gI1,gI2))*CpVZpbarFeFePR(gI1,gI2));
      }
      tmp_4037 += tmp_4038;
   }
   result += tmp_4037;
   std::complex<double> tmp_4039;
   std::complex<double> tmp_4040;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4041;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4041 += (AbsSqr(CpVZpbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZpbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_4041 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZpbarFuFuPL(gI1,gI2))*CpVZpbarFuFuPR(gI1,gI2));
      }
      tmp_4040 += tmp_4041;
   }
   tmp_4039 += tmp_4040;
   result += (3) * tmp_4039;
   std::complex<double> tmp_4042;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4043;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4043 += (AbsSqr(CpVZpbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZpbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_4043 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZpbarFvFvPL(gI1,gI2))*CpVZpbarFvFvPR(gI1,gI2));
      }
      tmp_4042 += tmp_4043;
   }
   result += tmp_4042;
   std::complex<double> tmp_4044;
   std::complex<double> tmp_4045;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4045 += A0(MSd(gI1))*CpVZpVZpconjSdSd(gI1,gI1);
   }
   tmp_4044 += tmp_4045;
   result += (3) * tmp_4044;
   std::complex<double> tmp_4046;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4046 += A0(MSe(gI1))*CpVZpVZpconjSeSe(gI1,gI1);
   }
   result += tmp_4046;
   std::complex<double> tmp_4047;
   std::complex<double> tmp_4048;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4048 += A0(MSu(gI1))*CpVZpVZpconjSuSu(gI1,gI1);
   }
   tmp_4047 += tmp_4048;
   result += (3) * tmp_4047;
   std::complex<double> tmp_4049;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4049 += A0(MSv(gI1))*CpVZpVZpconjSvSv(gI1,gI1);
   }
   result += tmp_4049;
   std::complex<double> tmp_4050;
   std::complex<double> tmp_4051;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4052;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4052 += AbsSqr(CpVZpconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(
            gI2));
      }
      tmp_4051 += tmp_4052;
   }
   tmp_4050 += tmp_4051;
   result += (-12) * tmp_4050;
   std::complex<double> tmp_4053;
   std::complex<double> tmp_4054;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4055;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4055 += AbsSqr(CpVZpconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(
            gI2));
      }
      tmp_4054 += tmp_4055;
   }
   tmp_4053 += tmp_4054;
   result += (-4) * tmp_4053;
   std::complex<double> tmp_4056;
   std::complex<double> tmp_4057;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4058;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4058 += AbsSqr(CpVZpconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(
            gI2));
      }
      tmp_4057 += tmp_4058;
   }
   tmp_4056 += tmp_4057;
   result += (-12) * tmp_4056;
   std::complex<double> tmp_4059;
   std::complex<double> tmp_4060;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4061;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4061 += AbsSqr(CpVZpconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(
            gI2));
      }
      tmp_4060 += tmp_4061;
   }
   tmp_4059 += tmp_4060;
   result += (-4) * tmp_4059;
   std::complex<double> tmp_4062;
   std::complex<double> tmp_4063;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4064;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4064 += (AbsSqr(CpVZpChiChiPL(gI1,gI2)) + AbsSqr(
            CpVZpChiChiPR(gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_4064 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZpChiChiPL(gI1,gI2))*CpVZpChiChiPR(gI1,gI2));
      }
      tmp_4063 += tmp_4064;
   }
   tmp_4062 += tmp_4063;
   result += (0.5) * tmp_4062;
   std::complex<double> tmp_4065;
   std::complex<double> tmp_4066;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4066 += AbsSqr(CpVZpconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_4065 += tmp_4066;
   result += (2) * tmp_4065;
   std::complex<double> tmp_4067;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4067 += AbsSqr(CpVZpVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_4067;
   std::complex<double> tmp_4068;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4068 += AbsSqr(CpVZpVZphh(gI2))*B0(p,MVZp,Mhh(gI2));
   }
   result += tmp_4068;
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
   std::complex<double> tmp_4069;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_4069 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_4069;
   std::complex<double> tmp_4070;
   std::complex<double> tmp_4071;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4072;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4072 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_4071 += tmp_4072;
   }
   tmp_4070 += tmp_4071;
   result += (-4) * tmp_4070;
   std::complex<double> tmp_4073;
   std::complex<double> tmp_4074;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4075;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4075 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_4074 += tmp_4075;
   }
   tmp_4073 += tmp_4074;
   result += (-4) * tmp_4073;
   std::complex<double> tmp_4076;
   std::complex<double> tmp_4077;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4077 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_4076 += tmp_4077;
   result += (0.5) * tmp_4076;
   std::complex<double> tmp_4078;
   std::complex<double> tmp_4079;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4079 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_4078 += tmp_4079;
   result += (0.5) * tmp_4078;
   std::complex<double> tmp_4080;
   std::complex<double> tmp_4081;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4082;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4082 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_4082 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_4081 += tmp_4082;
   }
   tmp_4080 += tmp_4081;
   result += (3) * tmp_4080;
   std::complex<double> tmp_4083;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4084;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4084 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_4084 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_4083 += tmp_4084;
   }
   result += tmp_4083;
   std::complex<double> tmp_4085;
   std::complex<double> tmp_4086;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4086 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_4085 += tmp_4086;
   result += (3) * tmp_4085;
   std::complex<double> tmp_4087;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4087 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_4087;
   std::complex<double> tmp_4088;
   std::complex<double> tmp_4089;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4089 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_4088 += tmp_4089;
   result += (3) * tmp_4088;
   std::complex<double> tmp_4090;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4090 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_4090;
   std::complex<double> tmp_4091;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4092;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4092 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_4092 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_4091 += tmp_4092;
   }
   result += tmp_4091;
   std::complex<double> tmp_4093;
   std::complex<double> tmp_4094;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4095;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4095 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_4094 += tmp_4095;
   }
   tmp_4093 += tmp_4094;
   result += (-12) * tmp_4093;
   std::complex<double> tmp_4096;
   std::complex<double> tmp_4097;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4098;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4098 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_4097 += tmp_4098;
   }
   tmp_4096 += tmp_4097;
   result += (-4) * tmp_4096;
   std::complex<double> tmp_4099;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4099 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_4099;
   std::complex<double> tmp_4100;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4100 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_4100;
   std::complex<double> tmp_4101;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4101 += AbsSqr(CpconjVWmVZpHpm(gI2))*B0(p,MVZp,MHpm(gI2));
   }
   result += tmp_4101;
   std::complex<double> tmp_4102;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4102 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_4102;
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

   std::complex<double> tmp_4103;
   std::complex<double> tmp_4104;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4105;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4105 += B0(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_4104 += tmp_4105;
   }
   tmp_4103 += tmp_4104;
   result += (2) * tmp_4103;
   std::complex<double> tmp_4106;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4107;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4107 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4106 += tmp_4107;
   }
   result += tmp_4106;
   std::complex<double> tmp_4108;
   std::complex<double> tmp_4109;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4110;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4110 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4109 += tmp_4110;
   }
   tmp_4108 += tmp_4109;
   result += (6) * tmp_4108;
   std::complex<double> tmp_4111;
   std::complex<double> tmp_4112;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4113;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4113 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_4112 += tmp_4113;
   }
   tmp_4111 += tmp_4112;
   result += (2) * tmp_4111;
   std::complex<double> tmp_4114;
   std::complex<double> tmp_4115;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4116;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4116 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_4115 += tmp_4116;
   }
   tmp_4114 += tmp_4115;
   result += (6) * tmp_4114;
   std::complex<double> tmp_4117;
   std::complex<double> tmp_4118;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4119;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4119 += B0(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_4118 += tmp_4119;
   }
   tmp_4117 += tmp_4118;
   result += (2) * tmp_4117;
   std::complex<double> tmp_4120;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4121;
      std::complex<double> tmp_4122;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4122 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_4121 += tmp_4122;
      tmp_4120 += (MChi(gI1)) * tmp_4121;
   }
   result += tmp_4120;
   std::complex<double> tmp_4123;
   std::complex<double> tmp_4124;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4124 += B0(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_4123 += tmp_4124;
   result += (-8) * tmp_4123;
   std::complex<double> tmp_4125;
   std::complex<double> tmp_4126;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4126 += B0(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_4125 += tmp_4126;
   result += (-4) * tmp_4125;
   std::complex<double> tmp_4127;
   std::complex<double> tmp_4128;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4128 += B0(p,MChi(gI2),MVZp)*Conj(CpUChiVZpChiPR(gO2,gI2))*
         CpUChiVZpChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_4127 += tmp_4128;
   result += (-4) * tmp_4127;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4129;
   std::complex<double> tmp_4130;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4131;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4131 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPR(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2);
      }
      tmp_4130 += tmp_4131;
   }
   tmp_4129 += tmp_4130;
   result += (-1) * tmp_4129;
   std::complex<double> tmp_4132;
   std::complex<double> tmp_4133;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4134;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4134 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_4133 += tmp_4134;
   }
   tmp_4132 += tmp_4133;
   result += (-0.5) * tmp_4132;
   std::complex<double> tmp_4135;
   std::complex<double> tmp_4136;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4137;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4137 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_4136 += tmp_4137;
   }
   tmp_4135 += tmp_4136;
   result += (-0.5) * tmp_4135;
   std::complex<double> tmp_4138;
   std::complex<double> tmp_4139;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4140;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4140 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPR(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_4139 += tmp_4140;
   }
   tmp_4138 += tmp_4139;
   result += (-3) * tmp_4138;
   std::complex<double> tmp_4141;
   std::complex<double> tmp_4142;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4143;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4143 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePR(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_4142 += tmp_4143;
   }
   tmp_4141 += tmp_4142;
   result += (-1) * tmp_4141;
   std::complex<double> tmp_4144;
   std::complex<double> tmp_4145;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4146;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4146 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPR(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_4145 += tmp_4146;
   }
   tmp_4144 += tmp_4145;
   result += (-3) * tmp_4144;
   std::complex<double> tmp_4147;
   std::complex<double> tmp_4148;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4149;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4149 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPR(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2);
      }
      tmp_4148 += tmp_4149;
   }
   tmp_4147 += tmp_4148;
   result += (-1) * tmp_4147;
   std::complex<double> tmp_4150;
   std::complex<double> tmp_4151;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4151 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPL(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_4150 += tmp_4151;
   result += (-2) * tmp_4150;
   std::complex<double> tmp_4152;
   std::complex<double> tmp_4153;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4153 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_4152 += tmp_4153;
   result += (-1) * tmp_4152;
   std::complex<double> tmp_4154;
   std::complex<double> tmp_4155;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4155 += B1(p,MChi(gI2),MVZp)*Conj(CpUChiVZpChiPL(gO2,gI2))*
         CpUChiVZpChiPL(gO1,gI2);
   }
   tmp_4154 += tmp_4155;
   result += (-1) * tmp_4154;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4156;
   std::complex<double> tmp_4157;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4158;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4158 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPL(gO1,gI1,gI2);
      }
      tmp_4157 += tmp_4158;
   }
   tmp_4156 += tmp_4157;
   result += (-1) * tmp_4156;
   std::complex<double> tmp_4159;
   std::complex<double> tmp_4160;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4161;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4161 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_4160 += tmp_4161;
   }
   tmp_4159 += tmp_4160;
   result += (-0.5) * tmp_4159;
   std::complex<double> tmp_4162;
   std::complex<double> tmp_4163;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4164;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4164 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_4163 += tmp_4164;
   }
   tmp_4162 += tmp_4163;
   result += (-0.5) * tmp_4162;
   std::complex<double> tmp_4165;
   std::complex<double> tmp_4166;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4167;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4167 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_4166 += tmp_4167;
   }
   tmp_4165 += tmp_4166;
   result += (-3) * tmp_4165;
   std::complex<double> tmp_4168;
   std::complex<double> tmp_4169;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4170;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4170 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_4169 += tmp_4170;
   }
   tmp_4168 += tmp_4169;
   result += (-1) * tmp_4168;
   std::complex<double> tmp_4171;
   std::complex<double> tmp_4172;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4173;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4173 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_4172 += tmp_4173;
   }
   tmp_4171 += tmp_4172;
   result += (-3) * tmp_4171;
   std::complex<double> tmp_4174;
   std::complex<double> tmp_4175;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4176;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4176 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPL(gO1,gI1,gI2);
      }
      tmp_4175 += tmp_4176;
   }
   tmp_4174 += tmp_4175;
   result += (-1) * tmp_4174;
   std::complex<double> tmp_4177;
   std::complex<double> tmp_4178;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4178 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_4177 += tmp_4178;
   result += (-2) * tmp_4177;
   std::complex<double> tmp_4179;
   std::complex<double> tmp_4180;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4180 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_4179 += tmp_4180;
   result += (-1) * tmp_4179;
   std::complex<double> tmp_4181;
   std::complex<double> tmp_4182;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4182 += B1(p,MChi(gI2),MVZp)*Conj(CpUChiVZpChiPR(gO2,gI2))*
         CpUChiVZpChiPR(gO1,gI2);
   }
   tmp_4181 += tmp_4182;
   result += (-1) * tmp_4181;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4183;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4184;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4184 += B0(p,MFe(gI2),MHpm(gI1))*Conj(CpbarUFvconjHpmFePL(
            gO2,gI1,gI2))*CpbarUFvconjHpmFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_4183 += tmp_4184;
   }
   result += tmp_4183;
   std::complex<double> tmp_4185;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4186;
      std::complex<double> tmp_4187;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4187 += B0(p,MCha(gI1),MSe(gI2))*Conj(CpbarUFvbarChaSePL(gO2
            ,gI1,gI2))*CpbarUFvbarChaSePR(gO1,gI1,gI2);
      }
      tmp_4186 += tmp_4187;
      tmp_4185 += (MCha(gI1)) * tmp_4186;
   }
   result += tmp_4185;
   std::complex<double> tmp_4188;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4189;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4189 += B0(p,MFv(gI2),Mhh(gI1))*Conj(CpbarUFvhhFvPL(gO2,gI1,
            gI2))*CpbarUFvhhFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_4188 += tmp_4189;
   }
   result += tmp_4188;
   std::complex<double> tmp_4190;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4191;
      std::complex<double> tmp_4192;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4192 += B0(p,MFv(gI1),MAh(gI2))*Conj(CpbarUFvFvAhPL(gO2,gI1,
            gI2))*CpbarUFvFvAhPR(gO1,gI1,gI2);
      }
      tmp_4191 += tmp_4192;
      tmp_4190 += (MFv(gI1)) * tmp_4191;
   }
   result += tmp_4190;
   std::complex<double> tmp_4193;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4194;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4194 += B0(p,MChi(gI2),MSv(gI1))*Conj(CpbarUFvSvChiPL(gO2,
            gI1,gI2))*CpbarUFvSvChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4193 += tmp_4194;
   }
   result += tmp_4193;
   std::complex<double> tmp_4195;
   std::complex<double> tmp_4196;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4196 += B0(p,MFe(gI2),MVWm)*Conj(CpbarUFvconjVWmFePR(gO2,gI2))*
         CpbarUFvconjVWmFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_4195 += tmp_4196;
   result += (-4) * tmp_4195;
   std::complex<double> tmp_4197;
   std::complex<double> tmp_4198;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4198 += B0(p,MFv(gI2),0)*Conj(CpbarUFvVPFvPR(gO2,gI2))*
         CpbarUFvVPFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_4197 += tmp_4198;
   result += (-4) * tmp_4197;
   std::complex<double> tmp_4199;
   std::complex<double> tmp_4200;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4200 += B0(p,MFv(gI2),MVZ)*Conj(CpbarUFvVZFvPR(gO2,gI2))*
         CpbarUFvVZFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_4199 += tmp_4200;
   result += (-4) * tmp_4199;
   std::complex<double> tmp_4201;
   std::complex<double> tmp_4202;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4202 += B0(p,MFv(gI2),MVZp)*Conj(CpbarUFvVZpFvPR(gO2,gI2))*
         CpbarUFvVZpFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_4201 += tmp_4202;
   result += (-4) * tmp_4201;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4203;
   std::complex<double> tmp_4204;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4205;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4205 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpbarUFvconjHpmFePR(
            gO2,gI1,gI2))*CpbarUFvconjHpmFePR(gO1,gI1,gI2);
      }
      tmp_4204 += tmp_4205;
   }
   tmp_4203 += tmp_4204;
   result += (-0.5) * tmp_4203;
   std::complex<double> tmp_4206;
   std::complex<double> tmp_4207;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4208;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4208 += B1(p,MCha(gI1),MSe(gI2))*Conj(CpbarUFvbarChaSePR(gO2
            ,gI1,gI2))*CpbarUFvbarChaSePR(gO1,gI1,gI2);
      }
      tmp_4207 += tmp_4208;
   }
   tmp_4206 += tmp_4207;
   result += (-0.5) * tmp_4206;
   std::complex<double> tmp_4209;
   std::complex<double> tmp_4210;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4211;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4211 += B1(p,MFv(gI1),MAh(gI2))*Conj(CpbarUFvFvAhPR(gO2,gI1,
            gI2))*CpbarUFvFvAhPR(gO1,gI1,gI2);
      }
      tmp_4210 += tmp_4211;
   }
   tmp_4209 += tmp_4210;
   result += (-0.5) * tmp_4209;
   std::complex<double> tmp_4212;
   std::complex<double> tmp_4213;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4214;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4214 += B1(p,MFv(gI2),Mhh(gI1))*Conj(CpbarUFvhhFvPR(gO2,gI1,
            gI2))*CpbarUFvhhFvPR(gO1,gI1,gI2);
      }
      tmp_4213 += tmp_4214;
   }
   tmp_4212 += tmp_4213;
   result += (-0.5) * tmp_4212;
   std::complex<double> tmp_4215;
   std::complex<double> tmp_4216;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4217;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4217 += B1(p,MChi(gI2),MSv(gI1))*Conj(CpbarUFvSvChiPR(gO2,
            gI1,gI2))*CpbarUFvSvChiPR(gO1,gI1,gI2);
      }
      tmp_4216 += tmp_4217;
   }
   tmp_4215 += tmp_4216;
   result += (-0.5) * tmp_4215;
   std::complex<double> tmp_4218;
   std::complex<double> tmp_4219;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4219 += B1(p,MFe(gI2),MVWm)*Conj(CpbarUFvconjVWmFePL(gO2,gI2))*
         CpbarUFvconjVWmFePL(gO1,gI2);
   }
   tmp_4218 += tmp_4219;
   result += (-1) * tmp_4218;
   std::complex<double> tmp_4220;
   std::complex<double> tmp_4221;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4221 += B1(p,MFv(gI2),0)*Conj(CpbarUFvVPFvPL(gO2,gI2))*
         CpbarUFvVPFvPL(gO1,gI2);
   }
   tmp_4220 += tmp_4221;
   result += (-1) * tmp_4220;
   std::complex<double> tmp_4222;
   std::complex<double> tmp_4223;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4223 += B1(p,MFv(gI2),MVZ)*Conj(CpbarUFvVZFvPL(gO2,gI2))*
         CpbarUFvVZFvPL(gO1,gI2);
   }
   tmp_4222 += tmp_4223;
   result += (-1) * tmp_4222;
   std::complex<double> tmp_4224;
   std::complex<double> tmp_4225;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4225 += B1(p,MFv(gI2),MVZp)*Conj(CpbarUFvVZpFvPL(gO2,gI2))*
         CpbarUFvVZpFvPL(gO1,gI2);
   }
   tmp_4224 += tmp_4225;
   result += (-1) * tmp_4224;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fv_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4226;
   std::complex<double> tmp_4227;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4228;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4228 += B1(p,MFe(gI2),MHpm(gI1))*Conj(CpbarUFvconjHpmFePL(
            gO2,gI1,gI2))*CpbarUFvconjHpmFePL(gO1,gI1,gI2);
      }
      tmp_4227 += tmp_4228;
   }
   tmp_4226 += tmp_4227;
   result += (-0.5) * tmp_4226;
   std::complex<double> tmp_4229;
   std::complex<double> tmp_4230;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4231;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4231 += B1(p,MCha(gI1),MSe(gI2))*Conj(CpbarUFvbarChaSePL(gO2
            ,gI1,gI2))*CpbarUFvbarChaSePL(gO1,gI1,gI2);
      }
      tmp_4230 += tmp_4231;
   }
   tmp_4229 += tmp_4230;
   result += (-0.5) * tmp_4229;
   std::complex<double> tmp_4232;
   std::complex<double> tmp_4233;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4234;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4234 += B1(p,MFv(gI1),MAh(gI2))*Conj(CpbarUFvFvAhPL(gO2,gI1,
            gI2))*CpbarUFvFvAhPL(gO1,gI1,gI2);
      }
      tmp_4233 += tmp_4234;
   }
   tmp_4232 += tmp_4233;
   result += (-0.5) * tmp_4232;
   std::complex<double> tmp_4235;
   std::complex<double> tmp_4236;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4237;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4237 += B1(p,MFv(gI2),Mhh(gI1))*Conj(CpbarUFvhhFvPL(gO2,gI1,
            gI2))*CpbarUFvhhFvPL(gO1,gI1,gI2);
      }
      tmp_4236 += tmp_4237;
   }
   tmp_4235 += tmp_4236;
   result += (-0.5) * tmp_4235;
   std::complex<double> tmp_4238;
   std::complex<double> tmp_4239;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4240;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4240 += B1(p,MChi(gI2),MSv(gI1))*Conj(CpbarUFvSvChiPL(gO2,
            gI1,gI2))*CpbarUFvSvChiPL(gO1,gI1,gI2);
      }
      tmp_4239 += tmp_4240;
   }
   tmp_4238 += tmp_4239;
   result += (-0.5) * tmp_4238;
   std::complex<double> tmp_4241;
   std::complex<double> tmp_4242;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4242 += B1(p,MFe(gI2),MVWm)*Conj(CpbarUFvconjVWmFePR(gO2,gI2))*
         CpbarUFvconjVWmFePR(gO1,gI2);
   }
   tmp_4241 += tmp_4242;
   result += (-1) * tmp_4241;
   std::complex<double> tmp_4243;
   std::complex<double> tmp_4244;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4244 += B1(p,MFv(gI2),0)*Conj(CpbarUFvVPFvPR(gO2,gI2))*
         CpbarUFvVPFvPR(gO1,gI2);
   }
   tmp_4243 += tmp_4244;
   result += (-1) * tmp_4243;
   std::complex<double> tmp_4245;
   std::complex<double> tmp_4246;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4246 += B1(p,MFv(gI2),MVZ)*Conj(CpbarUFvVZFvPR(gO2,gI2))*
         CpbarUFvVZFvPR(gO1,gI2);
   }
   tmp_4245 += tmp_4246;
   result += (-1) * tmp_4245;
   std::complex<double> tmp_4247;
   std::complex<double> tmp_4248;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4248 += B1(p,MFv(gI2),MVZp)*Conj(CpbarUFvVZpFvPR(gO2,gI2))*
         CpbarUFvVZpFvPR(gO1,gI2);
   }
   tmp_4247 += tmp_4248;
   result += (-1) * tmp_4247;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4249;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4250;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4250 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4249 += tmp_4250;
   }
   result += tmp_4249;
   std::complex<double> tmp_4251;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4252;
      std::complex<double> tmp_4253;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4253 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_4252 += tmp_4253;
      tmp_4251 += (MCha(gI1)) * tmp_4252;
   }
   result += tmp_4251;
   std::complex<double> tmp_4254;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4255;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4255 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_4254 += tmp_4255;
   }
   result += tmp_4254;
   std::complex<double> tmp_4256;
   std::complex<double> tmp_4257;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4258;
      std::complex<double> tmp_4259;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4259 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_4258 += tmp_4259;
      tmp_4257 += (MFu(gI1)) * tmp_4258;
   }
   tmp_4256 += tmp_4257;
   result += (3) * tmp_4256;
   std::complex<double> tmp_4260;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4261;
      std::complex<double> tmp_4262;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4262 += B0(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_4261 += tmp_4262;
      tmp_4260 += (MFv(gI1)) * tmp_4261;
   }
   result += tmp_4260;
   std::complex<double> tmp_4263;
   std::complex<double> tmp_4264;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4265;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4265 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4264 += tmp_4265;
   }
   tmp_4263 += tmp_4264;
   result += (3) * tmp_4263;
   std::complex<double> tmp_4266;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4267;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4267 += B0(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_4266 += tmp_4267;
   }
   result += tmp_4266;
   std::complex<double> tmp_4268;
   std::complex<double> tmp_4269;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4269 += B0(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_4268 += tmp_4269;
   result += (-4) * tmp_4268;
   std::complex<double> tmp_4270;
   std::complex<double> tmp_4271;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4271 += B0(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_4270 += tmp_4271;
   result += (-4) * tmp_4270;
   std::complex<double> tmp_4272;
   std::complex<double> tmp_4273;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4273 += B0(p,MCha(gI2),MVZp)*Conj(CpbarUChaVZpChaPR(gO2,gI2))*
         CpbarUChaVZpChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_4272 += tmp_4273;
   result += (-4) * tmp_4272;
   std::complex<double> tmp_4274;
   std::complex<double> tmp_4275;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4275 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_4274 += tmp_4275;
   result += (-4) * tmp_4274;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4276;
   std::complex<double> tmp_4277;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4278;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4278 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_4277 += tmp_4278;
   }
   tmp_4276 += tmp_4277;
   result += (-0.5) * tmp_4276;
   std::complex<double> tmp_4279;
   std::complex<double> tmp_4280;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4281;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4281 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPR(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2);
      }
      tmp_4280 += tmp_4281;
   }
   tmp_4279 += tmp_4280;
   result += (-0.5) * tmp_4279;
   std::complex<double> tmp_4282;
   std::complex<double> tmp_4283;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4284;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4284 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_4283 += tmp_4284;
   }
   tmp_4282 += tmp_4283;
   result += (-0.5) * tmp_4282;
   std::complex<double> tmp_4285;
   std::complex<double> tmp_4286;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4287;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4287 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPR(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_4286 += tmp_4287;
   }
   tmp_4285 += tmp_4286;
   result += (-1.5) * tmp_4285;
   std::complex<double> tmp_4288;
   std::complex<double> tmp_4289;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4290;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4290 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePR(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_4289 += tmp_4290;
   }
   tmp_4288 += tmp_4289;
   result += (-0.5) * tmp_4288;
   std::complex<double> tmp_4291;
   std::complex<double> tmp_4292;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4293;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4293 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPR(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2);
      }
      tmp_4292 += tmp_4293;
   }
   tmp_4291 += tmp_4292;
   result += (-1.5) * tmp_4291;
   std::complex<double> tmp_4294;
   std::complex<double> tmp_4295;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4296;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4296 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePR(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2);
      }
      tmp_4295 += tmp_4296;
   }
   tmp_4294 += tmp_4295;
   result += (-0.5) * tmp_4294;
   std::complex<double> tmp_4297;
   std::complex<double> tmp_4298;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4298 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_4297 += tmp_4298;
   result += (-1) * tmp_4297;
   std::complex<double> tmp_4299;
   std::complex<double> tmp_4300;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4300 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPL(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_4299 += tmp_4300;
   result += (-1) * tmp_4299;
   std::complex<double> tmp_4301;
   std::complex<double> tmp_4302;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4302 += B1(p,MCha(gI2),MVZp)*Conj(CpbarUChaVZpChaPL(gO2,gI2))*
         CpbarUChaVZpChaPL(gO1,gI2);
   }
   tmp_4301 += tmp_4302;
   result += (-1) * tmp_4301;
   std::complex<double> tmp_4303;
   std::complex<double> tmp_4304;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4304 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPL(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_4303 += tmp_4304;
   result += (-1) * tmp_4303;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4305;
   std::complex<double> tmp_4306;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4307;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4307 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_4306 += tmp_4307;
   }
   tmp_4305 += tmp_4306;
   result += (-0.5) * tmp_4305;
   std::complex<double> tmp_4308;
   std::complex<double> tmp_4309;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4310;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4310 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPL(gO1,gI1,gI2);
      }
      tmp_4309 += tmp_4310;
   }
   tmp_4308 += tmp_4309;
   result += (-0.5) * tmp_4308;
   std::complex<double> tmp_4311;
   std::complex<double> tmp_4312;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4313;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4313 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_4312 += tmp_4313;
   }
   tmp_4311 += tmp_4312;
   result += (-0.5) * tmp_4311;
   std::complex<double> tmp_4314;
   std::complex<double> tmp_4315;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4316;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4316 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2);
      }
      tmp_4315 += tmp_4316;
   }
   tmp_4314 += tmp_4315;
   result += (-1.5) * tmp_4314;
   std::complex<double> tmp_4317;
   std::complex<double> tmp_4318;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4319;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4319 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2);
      }
      tmp_4318 += tmp_4319;
   }
   tmp_4317 += tmp_4318;
   result += (-0.5) * tmp_4317;
   std::complex<double> tmp_4320;
   std::complex<double> tmp_4321;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4322;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4322 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPL(gO1,gI1,gI2);
      }
      tmp_4321 += tmp_4322;
   }
   tmp_4320 += tmp_4321;
   result += (-1.5) * tmp_4320;
   std::complex<double> tmp_4323;
   std::complex<double> tmp_4324;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4325;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4325 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePL(gO1,gI1,gI2);
      }
      tmp_4324 += tmp_4325;
   }
   tmp_4323 += tmp_4324;
   result += (-0.5) * tmp_4323;
   std::complex<double> tmp_4326;
   std::complex<double> tmp_4327;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4327 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_4326 += tmp_4327;
   result += (-1) * tmp_4326;
   std::complex<double> tmp_4328;
   std::complex<double> tmp_4329;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4329 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_4328 += tmp_4329;
   result += (-1) * tmp_4328;
   std::complex<double> tmp_4330;
   std::complex<double> tmp_4331;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_4331 += B1(p,MCha(gI2),MVZp)*Conj(CpbarUChaVZpChaPR(gO2,gI2))*
         CpbarUChaVZpChaPR(gO1,gI2);
   }
   tmp_4330 += tmp_4331;
   result += (-1) * tmp_4330;
   std::complex<double> tmp_4332;
   std::complex<double> tmp_4333;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_4333 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_4332 += tmp_4333;
   result += (-1) * tmp_4332;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4334;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4335;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4335 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_4334 += tmp_4335;
   }
   result += tmp_4334;
   std::complex<double> tmp_4336;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4337;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4337 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_4336 += tmp_4337;
   }
   result += tmp_4336;
   std::complex<double> tmp_4338;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4339;
      std::complex<double> tmp_4340;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4340 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_4339 += tmp_4340;
      tmp_4338 += (MFe(gI1)) * tmp_4339;
   }
   result += tmp_4338;
   std::complex<double> tmp_4341;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4342;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4342 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_4341 += tmp_4342;
   }
   result += tmp_4341;
   std::complex<double> tmp_4343;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4344;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4344 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4343 += tmp_4344;
   }
   result += tmp_4343;
   std::complex<double> tmp_4345;
   std::complex<double> tmp_4346;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4346 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_4345 += tmp_4346;
   result += (-4) * tmp_4345;
   std::complex<double> tmp_4347;
   std::complex<double> tmp_4348;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4348 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_4347 += tmp_4348;
   result += (-4) * tmp_4347;
   std::complex<double> tmp_4349;
   std::complex<double> tmp_4350;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4350 += B0(p,MFe(gI2),MVZp)*Conj(CpbarUFeVZpFePR(gO2,gI2))*
         CpbarUFeVZpFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_4349 += tmp_4350;
   result += (-4) * tmp_4349;
   std::complex<double> tmp_4351;
   std::complex<double> tmp_4352;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4352 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_4351 += tmp_4352;
   result += (-4) * tmp_4351;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4353;
   std::complex<double> tmp_4354;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4355;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4355 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_4354 += tmp_4355;
   }
   tmp_4353 += tmp_4354;
   result += (-0.5) * tmp_4353;
   std::complex<double> tmp_4356;
   std::complex<double> tmp_4357;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4358;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4358 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_4357 += tmp_4358;
   }
   tmp_4356 += tmp_4357;
   result += (-0.5) * tmp_4356;
   std::complex<double> tmp_4359;
   std::complex<double> tmp_4360;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4361;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4361 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_4360 += tmp_4361;
   }
   tmp_4359 += tmp_4360;
   result += (-0.5) * tmp_4359;
   std::complex<double> tmp_4362;
   std::complex<double> tmp_4363;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4364;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4364 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPR(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_4363 += tmp_4364;
   }
   tmp_4362 += tmp_4363;
   result += (-0.5) * tmp_4362;
   std::complex<double> tmp_4365;
   std::complex<double> tmp_4366;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4367;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4367 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_4366 += tmp_4367;
   }
   tmp_4365 += tmp_4366;
   result += (-0.5) * tmp_4365;
   std::complex<double> tmp_4368;
   std::complex<double> tmp_4369;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4369 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_4368 += tmp_4369;
   result += (-1) * tmp_4368;
   std::complex<double> tmp_4370;
   std::complex<double> tmp_4371;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4371 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_4370 += tmp_4371;
   result += (-1) * tmp_4370;
   std::complex<double> tmp_4372;
   std::complex<double> tmp_4373;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4373 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_4372 += tmp_4373;
   result += (-1) * tmp_4372;
   std::complex<double> tmp_4374;
   std::complex<double> tmp_4375;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4375 += B1(p,MFe(gI2),MVZp)*Conj(CpbarUFeVZpFePL(gO2,gI2))*
         CpbarUFeVZpFePL(gO1,gI2);
   }
   tmp_4374 += tmp_4375;
   result += (-1) * tmp_4374;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4376;
   std::complex<double> tmp_4377;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4378;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4378 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_4377 += tmp_4378;
   }
   tmp_4376 += tmp_4377;
   result += (-0.5) * tmp_4376;
   std::complex<double> tmp_4379;
   std::complex<double> tmp_4380;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4381;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4381 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_4380 += tmp_4381;
   }
   tmp_4379 += tmp_4380;
   result += (-0.5) * tmp_4379;
   std::complex<double> tmp_4382;
   std::complex<double> tmp_4383;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4384;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4384 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_4383 += tmp_4384;
   }
   tmp_4382 += tmp_4383;
   result += (-0.5) * tmp_4382;
   std::complex<double> tmp_4385;
   std::complex<double> tmp_4386;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4387;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4387 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_4386 += tmp_4387;
   }
   tmp_4385 += tmp_4386;
   result += (-0.5) * tmp_4385;
   std::complex<double> tmp_4388;
   std::complex<double> tmp_4389;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4390;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4390 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_4389 += tmp_4390;
   }
   tmp_4388 += tmp_4389;
   result += (-0.5) * tmp_4388;
   std::complex<double> tmp_4391;
   std::complex<double> tmp_4392;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4392 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_4391 += tmp_4392;
   result += (-1) * tmp_4391;
   std::complex<double> tmp_4393;
   std::complex<double> tmp_4394;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4394 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_4393 += tmp_4394;
   result += (-1) * tmp_4393;
   std::complex<double> tmp_4395;
   std::complex<double> tmp_4396;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4396 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_4395 += tmp_4396;
   result += (-1) * tmp_4395;
   std::complex<double> tmp_4397;
   std::complex<double> tmp_4398;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4398 += B1(p,MFe(gI2),MVZp)*Conj(CpbarUFeVZpFePR(gO2,gI2))*
         CpbarUFeVZpFePR(gO1,gI2);
   }
   tmp_4397 += tmp_4398;
   result += (-1) * tmp_4397;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4399;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4400;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4400 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_4399 += tmp_4400;
   }
   result += tmp_4399;
   std::complex<double> tmp_4401;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4402;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4402 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4401 += tmp_4402;
   }
   result += tmp_4401;
   std::complex<double> tmp_4403;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4404;
      std::complex<double> tmp_4405;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4405 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_4404 += tmp_4405;
      tmp_4403 += (MFd(gI1)) * tmp_4404;
   }
   result += tmp_4403;
   std::complex<double> tmp_4406;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4407;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4407 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_4406 += tmp_4407;
   }
   result += tmp_4406;
   std::complex<double> tmp_4408;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4409;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4409 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4408 += tmp_4409;
   }
   result += tmp_4408;
   std::complex<double> tmp_4410;
   std::complex<double> tmp_4411;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4411 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4410 += tmp_4411;
   result += (-5.333333333333333) * tmp_4410;
   std::complex<double> tmp_4412;
   std::complex<double> tmp_4413;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4413 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4412 += tmp_4413;
   result += (-4) * tmp_4412;
   std::complex<double> tmp_4414;
   std::complex<double> tmp_4415;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4415 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4414 += tmp_4415;
   result += (-4) * tmp_4414;
   std::complex<double> tmp_4416;
   std::complex<double> tmp_4417;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4417 += B0(p,MFd(gI2),MVZp)*Conj(CpbarUFdVZpFdPR(gO2,gI2))*
         CpbarUFdVZpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4416 += tmp_4417;
   result += (-4) * tmp_4416;
   std::complex<double> tmp_4418;
   std::complex<double> tmp_4419;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4419 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4418 += tmp_4419;
   result += (-4) * tmp_4418;
   std::complex<double> tmp_4420;
   std::complex<double> tmp_4421;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4421 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_4420 += tmp_4421;
   result += (1.3333333333333333*MGlu) * tmp_4420;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4422;
   std::complex<double> tmp_4423;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4424;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4424 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_4423 += tmp_4424;
   }
   tmp_4422 += tmp_4423;
   result += (-0.5) * tmp_4422;
   std::complex<double> tmp_4425;
   std::complex<double> tmp_4426;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4427;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4427 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_4426 += tmp_4427;
   }
   tmp_4425 += tmp_4426;
   result += (-0.5) * tmp_4425;
   std::complex<double> tmp_4428;
   std::complex<double> tmp_4429;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4430;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4430 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_4429 += tmp_4430;
   }
   tmp_4428 += tmp_4429;
   result += (-0.5) * tmp_4428;
   std::complex<double> tmp_4431;
   std::complex<double> tmp_4432;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4432 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_4431 += tmp_4432;
   result += (-0.6666666666666666) * tmp_4431;
   std::complex<double> tmp_4433;
   std::complex<double> tmp_4434;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4435;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4435 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPR(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_4434 += tmp_4435;
   }
   tmp_4433 += tmp_4434;
   result += (-0.5) * tmp_4433;
   std::complex<double> tmp_4436;
   std::complex<double> tmp_4437;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4438;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4438 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_4437 += tmp_4438;
   }
   tmp_4436 += tmp_4437;
   result += (-0.5) * tmp_4436;
   std::complex<double> tmp_4439;
   std::complex<double> tmp_4440;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4440 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_4439 += tmp_4440;
   result += (-1.3333333333333333) * tmp_4439;
   std::complex<double> tmp_4441;
   std::complex<double> tmp_4442;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4442 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_4441 += tmp_4442;
   result += (-1) * tmp_4441;
   std::complex<double> tmp_4443;
   std::complex<double> tmp_4444;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4444 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_4443 += tmp_4444;
   result += (-1) * tmp_4443;
   std::complex<double> tmp_4445;
   std::complex<double> tmp_4446;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4446 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_4445 += tmp_4446;
   result += (-1) * tmp_4445;
   std::complex<double> tmp_4447;
   std::complex<double> tmp_4448;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4448 += B1(p,MFd(gI2),MVZp)*Conj(CpbarUFdVZpFdPL(gO2,gI2))*
         CpbarUFdVZpFdPL(gO1,gI2);
   }
   tmp_4447 += tmp_4448;
   result += (-1) * tmp_4447;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4449;
   std::complex<double> tmp_4450;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4451;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4451 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_4450 += tmp_4451;
   }
   tmp_4449 += tmp_4450;
   result += (-0.5) * tmp_4449;
   std::complex<double> tmp_4452;
   std::complex<double> tmp_4453;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4454;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4454 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_4453 += tmp_4454;
   }
   tmp_4452 += tmp_4453;
   result += (-0.5) * tmp_4452;
   std::complex<double> tmp_4455;
   std::complex<double> tmp_4456;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4457;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4457 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_4456 += tmp_4457;
   }
   tmp_4455 += tmp_4456;
   result += (-0.5) * tmp_4455;
   std::complex<double> tmp_4458;
   std::complex<double> tmp_4459;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4459 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPL(gO1,gI1,1);
   }
   tmp_4458 += tmp_4459;
   result += (-0.6666666666666666) * tmp_4458;
   std::complex<double> tmp_4460;
   std::complex<double> tmp_4461;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4462;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4462 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_4461 += tmp_4462;
   }
   tmp_4460 += tmp_4461;
   result += (-0.5) * tmp_4460;
   std::complex<double> tmp_4463;
   std::complex<double> tmp_4464;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4465;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4465 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_4464 += tmp_4465;
   }
   tmp_4463 += tmp_4464;
   result += (-0.5) * tmp_4463;
   std::complex<double> tmp_4466;
   std::complex<double> tmp_4467;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4467 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_4466 += tmp_4467;
   result += (-1.3333333333333333) * tmp_4466;
   std::complex<double> tmp_4468;
   std::complex<double> tmp_4469;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4469 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_4468 += tmp_4469;
   result += (-1) * tmp_4468;
   std::complex<double> tmp_4470;
   std::complex<double> tmp_4471;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4471 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_4470 += tmp_4471;
   result += (-1) * tmp_4470;
   std::complex<double> tmp_4472;
   std::complex<double> tmp_4473;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4473 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_4472 += tmp_4473;
   result += (-1) * tmp_4472;
   std::complex<double> tmp_4474;
   std::complex<double> tmp_4475;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4475 += B1(p,MFd(gI2),MVZp)*Conj(CpbarUFdVZpFdPR(gO2,gI2))*
         CpbarUFdVZpFdPR(gO1,gI2);
   }
   tmp_4474 += tmp_4475;
   result += (-1) * tmp_4474;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4476;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4477;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4477 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4476 += tmp_4477;
   }
   result += tmp_4476;
   std::complex<double> tmp_4478;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4479;
      std::complex<double> tmp_4480;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4480 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_4479 += tmp_4480;
      tmp_4478 += (MCha(gI1)) * tmp_4479;
   }
   result += tmp_4478;
   std::complex<double> tmp_4481;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4482;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4482 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_4481 += tmp_4482;
   }
   result += tmp_4481;
   std::complex<double> tmp_4483;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4484;
      std::complex<double> tmp_4485;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4485 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_4484 += tmp_4485;
      tmp_4483 += (MFu(gI1)) * tmp_4484;
   }
   result += tmp_4483;
   std::complex<double> tmp_4486;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4487;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4487 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4486 += tmp_4487;
   }
   result += tmp_4486;
   std::complex<double> tmp_4488;
   std::complex<double> tmp_4489;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4489 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4488 += tmp_4489;
   result += (-4) * tmp_4488;
   std::complex<double> tmp_4490;
   std::complex<double> tmp_4491;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4491 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4490 += tmp_4491;
   result += (-5.333333333333333) * tmp_4490;
   std::complex<double> tmp_4492;
   std::complex<double> tmp_4493;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4493 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4492 += tmp_4493;
   result += (-4) * tmp_4492;
   std::complex<double> tmp_4494;
   std::complex<double> tmp_4495;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4495 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4494 += tmp_4495;
   result += (-4) * tmp_4494;
   std::complex<double> tmp_4496;
   std::complex<double> tmp_4497;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4497 += B0(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4496 += tmp_4497;
   result += (-4) * tmp_4496;
   std::complex<double> tmp_4498;
   std::complex<double> tmp_4499;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4499 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_4498 += tmp_4499;
   result += (1.3333333333333333*MGlu) * tmp_4498;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4500;
   std::complex<double> tmp_4501;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4502;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4502 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_4501 += tmp_4502;
   }
   tmp_4500 += tmp_4501;
   result += (-0.5) * tmp_4500;
   std::complex<double> tmp_4503;
   std::complex<double> tmp_4504;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4505;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4505 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_4504 += tmp_4505;
   }
   tmp_4503 += tmp_4504;
   result += (-0.5) * tmp_4503;
   std::complex<double> tmp_4506;
   std::complex<double> tmp_4507;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4508;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4508 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_4507 += tmp_4508;
   }
   tmp_4506 += tmp_4507;
   result += (-0.5) * tmp_4506;
   std::complex<double> tmp_4509;
   std::complex<double> tmp_4510;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4511;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4511 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_4510 += tmp_4511;
   }
   tmp_4509 += tmp_4510;
   result += (-0.5) * tmp_4509;
   std::complex<double> tmp_4512;
   std::complex<double> tmp_4513;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4513 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_4512 += tmp_4513;
   result += (-0.6666666666666666) * tmp_4512;
   std::complex<double> tmp_4514;
   std::complex<double> tmp_4515;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4516;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4516 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_4515 += tmp_4516;
   }
   tmp_4514 += tmp_4515;
   result += (-0.5) * tmp_4514;
   std::complex<double> tmp_4517;
   std::complex<double> tmp_4518;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4518 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_4517 += tmp_4518;
   result += (-1) * tmp_4517;
   std::complex<double> tmp_4519;
   std::complex<double> tmp_4520;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4520 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_4519 += tmp_4520;
   result += (-1.3333333333333333) * tmp_4519;
   std::complex<double> tmp_4521;
   std::complex<double> tmp_4522;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4522 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_4521 += tmp_4522;
   result += (-1) * tmp_4521;
   std::complex<double> tmp_4523;
   std::complex<double> tmp_4524;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4524 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_4523 += tmp_4524;
   result += (-1) * tmp_4523;
   std::complex<double> tmp_4525;
   std::complex<double> tmp_4526;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4526 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPL(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2);
   }
   tmp_4525 += tmp_4526;
   result += (-1) * tmp_4525;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4527;
   std::complex<double> tmp_4528;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4529;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4529 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_4528 += tmp_4529;
   }
   tmp_4527 += tmp_4528;
   result += (-0.5) * tmp_4527;
   std::complex<double> tmp_4530;
   std::complex<double> tmp_4531;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4532;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4532 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_4531 += tmp_4532;
   }
   tmp_4530 += tmp_4531;
   result += (-0.5) * tmp_4530;
   std::complex<double> tmp_4533;
   std::complex<double> tmp_4534;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4535;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4535 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_4534 += tmp_4535;
   }
   tmp_4533 += tmp_4534;
   result += (-0.5) * tmp_4533;
   std::complex<double> tmp_4536;
   std::complex<double> tmp_4537;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4538;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4538 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_4537 += tmp_4538;
   }
   tmp_4536 += tmp_4537;
   result += (-0.5) * tmp_4536;
   std::complex<double> tmp_4539;
   std::complex<double> tmp_4540;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4540 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_4539 += tmp_4540;
   result += (-0.6666666666666666) * tmp_4539;
   std::complex<double> tmp_4541;
   std::complex<double> tmp_4542;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4543;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4543 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_4542 += tmp_4543;
   }
   tmp_4541 += tmp_4542;
   result += (-0.5) * tmp_4541;
   std::complex<double> tmp_4544;
   std::complex<double> tmp_4545;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4545 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_4544 += tmp_4545;
   result += (-1) * tmp_4544;
   std::complex<double> tmp_4546;
   std::complex<double> tmp_4547;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4547 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_4546 += tmp_4547;
   result += (-1.3333333333333333) * tmp_4546;
   std::complex<double> tmp_4548;
   std::complex<double> tmp_4549;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4549 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_4548 += tmp_4549;
   result += (-1) * tmp_4548;
   std::complex<double> tmp_4550;
   std::complex<double> tmp_4551;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4551 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_4550 += tmp_4551;
   result += (-1) * tmp_4550;
   std::complex<double> tmp_4552;
   std::complex<double> tmp_4553;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4553 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPR(gO1,gI2);
   }
   tmp_4552 += tmp_4553;
   result += (-1) * tmp_4552;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_4554;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4555;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4555 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpGluconjSdFdPL(gI1,gI2
            ))*CpGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_4554 += tmp_4555;
   }
   result += tmp_4554;
   std::complex<double> tmp_4556;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4557;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4557 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpGluconjSuFuPL(gI1,gI2
            ))*CpGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_4556 += tmp_4557;
   }
   result += tmp_4556;
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpGluVGGluPR())*CpGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_4558;
   std::complex<double> tmp_4559;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4560;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4560 += AbsSqr(CpGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_4559 += tmp_4560;
   }
   tmp_4558 += tmp_4559;
   result += (-0.5) * tmp_4558;
   std::complex<double> tmp_4561;
   std::complex<double> tmp_4562;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4563;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4563 += AbsSqr(CpGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_4562 += tmp_4563;
   }
   tmp_4561 += tmp_4562;
   result += (-0.5) * tmp_4561;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_4564;
   std::complex<double> tmp_4565;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4566;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4566 += AbsSqr(CpGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_4565 += tmp_4566;
   }
   tmp_4564 += tmp_4565;
   result += (-0.5) * tmp_4564;
   std::complex<double> tmp_4567;
   std::complex<double> tmp_4568;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4569;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4569 += AbsSqr(CpGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_4568 += tmp_4569;
   }
   tmp_4567 += tmp_4568;
   result += (-0.5) * tmp_4567;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_4570;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4571;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4571 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_4571 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_4570 += tmp_4571;
   }
   result += tmp_4570;
   std::complex<double> tmp_4572;
   std::complex<double> tmp_4573;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4573 += AbsSqr(CpVZhhAh(gI1,2))*B00(p,MAh(2),Mhh(gI1));
   }
   tmp_4572 += tmp_4573;
   result += (-4) * tmp_4572;
   std::complex<double> tmp_4574;
   std::complex<double> tmp_4575;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4575 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_4574 += tmp_4575;
   result += (0.5) * tmp_4574;
   std::complex<double> tmp_4576;
   std::complex<double> tmp_4577;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4577 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_4576 += tmp_4577;
   result += (3) * tmp_4576;
   std::complex<double> tmp_4578;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4578 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_4578;
   std::complex<double> tmp_4579;
   std::complex<double> tmp_4580;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4580 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_4579 += tmp_4580;
   result += (3) * tmp_4579;
   std::complex<double> tmp_4581;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4581 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_4581;
   std::complex<double> tmp_4582;
   std::complex<double> tmp_4583;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4584;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4584 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_4583 += tmp_4584;
   }
   tmp_4582 += tmp_4583;
   result += (-12) * tmp_4582;
   std::complex<double> tmp_4585;
   std::complex<double> tmp_4586;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4587;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4587 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_4586 += tmp_4587;
   }
   tmp_4585 += tmp_4586;
   result += (-4) * tmp_4585;
   std::complex<double> tmp_4588;
   std::complex<double> tmp_4589;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4590;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4590 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_4589 += tmp_4590;
   }
   tmp_4588 += tmp_4589;
   result += (-12) * tmp_4588;
   std::complex<double> tmp_4591;
   std::complex<double> tmp_4592;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4593;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4593 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_4592 += tmp_4593;
   }
   tmp_4591 += tmp_4592;
   result += (-4) * tmp_4591;
   std::complex<double> tmp_4594;
   std::complex<double> tmp_4595;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4596;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4596 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_4596 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_4595 += tmp_4596;
   }
   tmp_4594 += tmp_4595;
   result += (0.5) * tmp_4594;
   std::complex<double> tmp_4597;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4597 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_4597;
   std::complex<double> tmp_4598;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4598 += AbsSqr(CpVZVZphh(gI2))*B0(p,MVZp,Mhh(gI2));
   }
   result += tmp_4598;

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
   std::complex<double> tmp_4599;
   std::complex<double> tmp_4600;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4600 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_4599 += tmp_4600;
   result += (0.5) * tmp_4599;
   std::complex<double> tmp_4601;
   std::complex<double> tmp_4602;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4602 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_4601 += tmp_4602;
   result += (3) * tmp_4601;
   std::complex<double> tmp_4603;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4603 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_4603;
   std::complex<double> tmp_4604;
   std::complex<double> tmp_4605;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4605 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_4604 += tmp_4605;
   result += (3) * tmp_4604;
   std::complex<double> tmp_4606;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4606 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_4606;
   std::complex<double> tmp_4607;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4608;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4608 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_4608 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_4607 += tmp_4608;
   }
   result += tmp_4607;
   std::complex<double> tmp_4609;
   std::complex<double> tmp_4610;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4611;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4611 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_4610 += tmp_4611;
   }
   tmp_4609 += tmp_4610;
   result += (-12) * tmp_4609;
   std::complex<double> tmp_4612;
   std::complex<double> tmp_4613;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4614;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4614 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_4613 += tmp_4614;
   }
   tmp_4612 += tmp_4613;
   result += (-4) * tmp_4612;
   std::complex<double> tmp_4615;
   std::complex<double> tmp_4616;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4616 += AbsSqr(CpconjVWmHpmhh(1,gI2))*B00(p,Mhh(gI2),MHpm(1));
   }
   tmp_4615 += tmp_4616;
   result += (-4) * tmp_4615;
   std::complex<double> tmp_4617;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4617 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_4617;
   result += -(AbsSqr(CpconjVWmVZpVWm())*(A0(MVWm) + A0(MVZp) + 10*B00(p,MVZp,
      MVWm) + B0(p,MVZp,MVWm)*(Sqr(MVWm) + Sqr(MVZp) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4618;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4619;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4619 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_4618 += tmp_4619;
   }
   result += tmp_4618;
   std::complex<double> tmp_4620;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4621;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4621 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_4620 += tmp_4621;
   }
   result += tmp_4620;
   std::complex<double> tmp_4622;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4623;
      std::complex<double> tmp_4624;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4624 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_4623 += tmp_4624;
      tmp_4622 += (MFe(gI1)) * tmp_4623;
   }
   result += tmp_4622;
   std::complex<double> tmp_4625;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4626;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4626 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_4625 += tmp_4626;
   }
   result += tmp_4625;
   std::complex<double> tmp_4627;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4628;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4628 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4627 += tmp_4628;
   }
   result += tmp_4627;
   std::complex<double> tmp_4629;
   std::complex<double> tmp_4630;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4630 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_4629 += tmp_4630;
   result += (-4) * tmp_4629;
   std::complex<double> tmp_4631;
   std::complex<double> tmp_4632;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4632 += B0(p,MFe(gI2),MVZp)*Conj(CpbarFeVZpFePR(gO2,gI2))*
         CpbarFeVZpFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_4631 += tmp_4632;
   result += (-4) * tmp_4631;
   std::complex<double> tmp_4633;
   std::complex<double> tmp_4634;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4634 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_4633 += tmp_4634;
   result += (-4) * tmp_4633;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4635;
   std::complex<double> tmp_4636;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4637;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4637 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_4636 += tmp_4637;
   }
   tmp_4635 += tmp_4636;
   result += (-0.5) * tmp_4635;
   std::complex<double> tmp_4638;
   std::complex<double> tmp_4639;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4640;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4640 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_4639 += tmp_4640;
   }
   tmp_4638 += tmp_4639;
   result += (-0.5) * tmp_4638;
   std::complex<double> tmp_4641;
   std::complex<double> tmp_4642;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4643;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4643 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_4642 += tmp_4643;
   }
   tmp_4641 += tmp_4642;
   result += (-0.5) * tmp_4641;
   std::complex<double> tmp_4644;
   std::complex<double> tmp_4645;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4646;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4646 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPR(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_4645 += tmp_4646;
   }
   tmp_4644 += tmp_4645;
   result += (-0.5) * tmp_4644;
   std::complex<double> tmp_4647;
   std::complex<double> tmp_4648;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4649;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4649 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_4648 += tmp_4649;
   }
   tmp_4647 += tmp_4648;
   result += (-0.5) * tmp_4647;
   std::complex<double> tmp_4650;
   std::complex<double> tmp_4651;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4651 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_4650 += tmp_4651;
   result += (-1) * tmp_4650;
   std::complex<double> tmp_4652;
   std::complex<double> tmp_4653;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4653 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_4652 += tmp_4653;
   result += (-1) * tmp_4652;
   std::complex<double> tmp_4654;
   std::complex<double> tmp_4655;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4655 += B1(p,MFe(gI2),MVZp)*Conj(CpbarFeVZpFePL(gO2,gI2))*
         CpbarFeVZpFePL(gO1,gI2);
   }
   tmp_4654 += tmp_4655;
   result += (-1) * tmp_4654;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4656;
   std::complex<double> tmp_4657;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4658;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4658 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_4657 += tmp_4658;
   }
   tmp_4656 += tmp_4657;
   result += (-0.5) * tmp_4656;
   std::complex<double> tmp_4659;
   std::complex<double> tmp_4660;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4661;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4661 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_4660 += tmp_4661;
   }
   tmp_4659 += tmp_4660;
   result += (-0.5) * tmp_4659;
   std::complex<double> tmp_4662;
   std::complex<double> tmp_4663;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4664;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4664 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_4663 += tmp_4664;
   }
   tmp_4662 += tmp_4663;
   result += (-0.5) * tmp_4662;
   std::complex<double> tmp_4665;
   std::complex<double> tmp_4666;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4667;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4667 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_4666 += tmp_4667;
   }
   tmp_4665 += tmp_4666;
   result += (-0.5) * tmp_4665;
   std::complex<double> tmp_4668;
   std::complex<double> tmp_4669;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4670;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4670 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_4669 += tmp_4670;
   }
   tmp_4668 += tmp_4669;
   result += (-0.5) * tmp_4668;
   std::complex<double> tmp_4671;
   std::complex<double> tmp_4672;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4672 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_4671 += tmp_4672;
   result += (-1) * tmp_4671;
   std::complex<double> tmp_4673;
   std::complex<double> tmp_4674;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4674 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_4673 += tmp_4674;
   result += (-1) * tmp_4673;
   std::complex<double> tmp_4675;
   std::complex<double> tmp_4676;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4676 += B1(p,MFe(gI2),MVZp)*Conj(CpbarFeVZpFePR(gO2,gI2))*
         CpbarFeVZpFePR(gO1,gI2);
   }
   tmp_4675 += tmp_4676;
   result += (-1) * tmp_4675;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4677;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4678;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4678 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_4677 += tmp_4678;
   }
   result += tmp_4677;
   std::complex<double> tmp_4679;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4680;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4680 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4679 += tmp_4680;
   }
   result += tmp_4679;
   std::complex<double> tmp_4681;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4682;
      std::complex<double> tmp_4683;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4683 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_4682 += tmp_4683;
      tmp_4681 += (MFd(gI1)) * tmp_4682;
   }
   result += tmp_4681;
   std::complex<double> tmp_4684;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4685;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4685 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_4684 += tmp_4685;
   }
   result += tmp_4684;
   std::complex<double> tmp_4686;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4687;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4687 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4686 += tmp_4687;
   }
   result += tmp_4686;
   std::complex<double> tmp_4688;
   std::complex<double> tmp_4689;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4689 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4688 += tmp_4689;
   result += (-4) * tmp_4688;
   std::complex<double> tmp_4690;
   std::complex<double> tmp_4691;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4691 += B0(p,MFd(gI2),MVZp)*Conj(CpbarFdVZpFdPR(gO2,gI2))*
         CpbarFdVZpFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4690 += tmp_4691;
   result += (-4) * tmp_4690;
   std::complex<double> tmp_4692;
   std::complex<double> tmp_4693;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4693 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4692 += tmp_4693;
   result += (-4) * tmp_4692;
   std::complex<double> tmp_4694;
   std::complex<double> tmp_4695;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4695 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_4694 += tmp_4695;
   result += (1.3333333333333333*MGlu) * tmp_4694;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4696;
   std::complex<double> tmp_4697;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4698;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4698 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_4697 += tmp_4698;
   }
   tmp_4696 += tmp_4697;
   result += (-0.5) * tmp_4696;
   std::complex<double> tmp_4699;
   std::complex<double> tmp_4700;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4701;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4701 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_4700 += tmp_4701;
   }
   tmp_4699 += tmp_4700;
   result += (-0.5) * tmp_4699;
   std::complex<double> tmp_4702;
   std::complex<double> tmp_4703;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4704;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4704 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_4703 += tmp_4704;
   }
   tmp_4702 += tmp_4703;
   result += (-0.5) * tmp_4702;
   std::complex<double> tmp_4705;
   std::complex<double> tmp_4706;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4706 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_4705 += tmp_4706;
   result += (-0.6666666666666666) * tmp_4705;
   std::complex<double> tmp_4707;
   std::complex<double> tmp_4708;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4709;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4709 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPR(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_4708 += tmp_4709;
   }
   tmp_4707 += tmp_4708;
   result += (-0.5) * tmp_4707;
   std::complex<double> tmp_4710;
   std::complex<double> tmp_4711;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4712;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4712 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_4711 += tmp_4712;
   }
   tmp_4710 += tmp_4711;
   result += (-0.5) * tmp_4710;
   std::complex<double> tmp_4713;
   std::complex<double> tmp_4714;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4714 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_4713 += tmp_4714;
   result += (-1) * tmp_4713;
   std::complex<double> tmp_4715;
   std::complex<double> tmp_4716;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4716 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_4715 += tmp_4716;
   result += (-1) * tmp_4715;
   std::complex<double> tmp_4717;
   std::complex<double> tmp_4718;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4718 += B1(p,MFd(gI2),MVZp)*Conj(CpbarFdVZpFdPL(gO2,gI2))*
         CpbarFdVZpFdPL(gO1,gI2);
   }
   tmp_4717 += tmp_4718;
   result += (-1) * tmp_4717;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4719;
   std::complex<double> tmp_4720;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4721;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4721 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_4720 += tmp_4721;
   }
   tmp_4719 += tmp_4720;
   result += (-0.5) * tmp_4719;
   std::complex<double> tmp_4722;
   std::complex<double> tmp_4723;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4724;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4724 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_4723 += tmp_4724;
   }
   tmp_4722 += tmp_4723;
   result += (-0.5) * tmp_4722;
   std::complex<double> tmp_4725;
   std::complex<double> tmp_4726;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4727;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4727 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_4726 += tmp_4727;
   }
   tmp_4725 += tmp_4726;
   result += (-0.5) * tmp_4725;
   std::complex<double> tmp_4728;
   std::complex<double> tmp_4729;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4729 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPL(gO1,gI1,1);
   }
   tmp_4728 += tmp_4729;
   result += (-0.6666666666666666) * tmp_4728;
   std::complex<double> tmp_4730;
   std::complex<double> tmp_4731;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4732;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_4732 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_4731 += tmp_4732;
   }
   tmp_4730 += tmp_4731;
   result += (-0.5) * tmp_4730;
   std::complex<double> tmp_4733;
   std::complex<double> tmp_4734;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4735;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4735 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_4734 += tmp_4735;
   }
   tmp_4733 += tmp_4734;
   result += (-0.5) * tmp_4733;
   std::complex<double> tmp_4736;
   std::complex<double> tmp_4737;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4737 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_4736 += tmp_4737;
   result += (-1) * tmp_4736;
   std::complex<double> tmp_4738;
   std::complex<double> tmp_4739;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4739 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_4738 += tmp_4739;
   result += (-1) * tmp_4738;
   std::complex<double> tmp_4740;
   std::complex<double> tmp_4741;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4741 += B1(p,MFd(gI2),MVZp)*Conj(CpbarFdVZpFdPR(gO2,gI2))*
         CpbarFdVZpFdPR(gO1,gI2);
   }
   tmp_4740 += tmp_4741;
   result += (-1) * tmp_4740;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4742;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4743;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4743 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4742 += tmp_4743;
   }
   result += tmp_4742;
   std::complex<double> tmp_4744;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4745;
      std::complex<double> tmp_4746;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4746 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_4745 += tmp_4746;
      tmp_4744 += (MCha(gI1)) * tmp_4745;
   }
   result += tmp_4744;
   std::complex<double> tmp_4747;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4748;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4748 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_4747 += tmp_4748;
   }
   result += tmp_4747;
   std::complex<double> tmp_4749;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4750;
      std::complex<double> tmp_4751;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4751 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_4750 += tmp_4751;
      tmp_4749 += (MFu(gI1)) * tmp_4750;
   }
   result += tmp_4749;
   std::complex<double> tmp_4752;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4753;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4753 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4752 += tmp_4753;
   }
   result += tmp_4752;
   std::complex<double> tmp_4754;
   std::complex<double> tmp_4755;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4755 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4754 += tmp_4755;
   result += (-4) * tmp_4754;
   std::complex<double> tmp_4756;
   std::complex<double> tmp_4757;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4757 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4756 += tmp_4757;
   result += (-4) * tmp_4756;
   std::complex<double> tmp_4758;
   std::complex<double> tmp_4759;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4759 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4758 += tmp_4759;
   result += (-4) * tmp_4758;
   std::complex<double> tmp_4760;
   std::complex<double> tmp_4761;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4761 += B0(p,MFu(gI2),MVZp)*Conj(CpbarFuVZpFuPR(gO2,gI2))*
         CpbarFuVZpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4760 += tmp_4761;
   result += (-4) * tmp_4760;
   std::complex<double> tmp_4762;
   std::complex<double> tmp_4763;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4763 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_4762 += tmp_4763;
   result += (1.3333333333333333*MGlu) * tmp_4762;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4764;
   std::complex<double> tmp_4765;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4766;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4766 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_4765 += tmp_4766;
   }
   tmp_4764 += tmp_4765;
   result += (-0.5) * tmp_4764;
   std::complex<double> tmp_4767;
   std::complex<double> tmp_4768;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4769;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4769 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPR(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_4768 += tmp_4769;
   }
   tmp_4767 += tmp_4768;
   result += (-0.5) * tmp_4767;
   std::complex<double> tmp_4770;
   std::complex<double> tmp_4771;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4772;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4772 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_4771 += tmp_4772;
   }
   tmp_4770 += tmp_4771;
   result += (-0.5) * tmp_4770;
   std::complex<double> tmp_4773;
   std::complex<double> tmp_4774;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4775;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4775 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_4774 += tmp_4775;
   }
   tmp_4773 += tmp_4774;
   result += (-0.5) * tmp_4773;
   std::complex<double> tmp_4776;
   std::complex<double> tmp_4777;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4777 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_4776 += tmp_4777;
   result += (-0.6666666666666666) * tmp_4776;
   std::complex<double> tmp_4778;
   std::complex<double> tmp_4779;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4780;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4780 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_4779 += tmp_4780;
   }
   tmp_4778 += tmp_4779;
   result += (-0.5) * tmp_4778;
   std::complex<double> tmp_4781;
   std::complex<double> tmp_4782;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4782 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_4781 += tmp_4782;
   result += (-1) * tmp_4781;
   std::complex<double> tmp_4783;
   std::complex<double> tmp_4784;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4784 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_4783 += tmp_4784;
   result += (-1) * tmp_4783;
   std::complex<double> tmp_4785;
   std::complex<double> tmp_4786;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4786 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_4785 += tmp_4786;
   result += (-1) * tmp_4785;
   std::complex<double> tmp_4787;
   std::complex<double> tmp_4788;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4788 += B1(p,MFu(gI2),MVZp)*Conj(CpbarFuVZpFuPL(gO2,gI2))*
         CpbarFuVZpFuPL(gO1,gI2);
   }
   tmp_4787 += tmp_4788;
   result += (-1) * tmp_4787;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4789;
   std::complex<double> tmp_4790;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4791;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4791 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_4790 += tmp_4791;
   }
   tmp_4789 += tmp_4790;
   result += (-0.5) * tmp_4789;
   std::complex<double> tmp_4792;
   std::complex<double> tmp_4793;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4794;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4794 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_4793 += tmp_4794;
   }
   tmp_4792 += tmp_4793;
   result += (-0.5) * tmp_4792;
   std::complex<double> tmp_4795;
   std::complex<double> tmp_4796;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4797;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4797 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_4796 += tmp_4797;
   }
   tmp_4795 += tmp_4796;
   result += (-0.5) * tmp_4795;
   std::complex<double> tmp_4798;
   std::complex<double> tmp_4799;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4800;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4800 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_4799 += tmp_4800;
   }
   tmp_4798 += tmp_4799;
   result += (-0.5) * tmp_4798;
   std::complex<double> tmp_4801;
   std::complex<double> tmp_4802;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4802 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPL(gO1,gI1,1);
   }
   tmp_4801 += tmp_4802;
   result += (-0.6666666666666666) * tmp_4801;
   std::complex<double> tmp_4803;
   std::complex<double> tmp_4804;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4805;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4805 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_4804 += tmp_4805;
   }
   tmp_4803 += tmp_4804;
   result += (-0.5) * tmp_4803;
   std::complex<double> tmp_4806;
   std::complex<double> tmp_4807;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4807 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_4806 += tmp_4807;
   result += (-1) * tmp_4806;
   std::complex<double> tmp_4808;
   std::complex<double> tmp_4809;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4809 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_4808 += tmp_4809;
   result += (-1) * tmp_4808;
   std::complex<double> tmp_4810;
   std::complex<double> tmp_4811;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4811 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_4810 += tmp_4811;
   result += (-1) * tmp_4810;
   std::complex<double> tmp_4812;
   std::complex<double> tmp_4813;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4813 += B1(p,MFu(gI2),MVZp)*Conj(CpbarFuVZpFuPR(gO2,gI2))*
         CpbarFuVZpFuPR(gO1,gI2);
   }
   tmp_4812 += tmp_4813;
   result += (-1) * tmp_4812;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4814;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4815;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4815 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_4814 += tmp_4815;
   }
   result += tmp_4814;
   std::complex<double> tmp_4816;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4817;
      std::complex<double> tmp_4818;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4818 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_4817 += tmp_4818;
      tmp_4816 += (MCha(gI1)) * tmp_4817;
   }
   result += tmp_4816;
   std::complex<double> tmp_4819;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4820;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4820 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_4819 += tmp_4820;
   }
   result += tmp_4819;
   std::complex<double> tmp_4821;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4822;
      std::complex<double> tmp_4823;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4823 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_4822 += tmp_4823;
      tmp_4821 += (MFu(gI1)) * tmp_4822;
   }
   result += tmp_4821;
   std::complex<double> tmp_4824;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4825;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4825 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_4824 += tmp_4825;
   }
   result += tmp_4824;
   std::complex<double> tmp_4826;
   std::complex<double> tmp_4827;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4827 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_4826 += tmp_4827;
   result += (-4) * tmp_4826;
   std::complex<double> tmp_4828;
   std::complex<double> tmp_4829;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4829 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4828 += tmp_4829;
   result += (-4) * tmp_4828;
   std::complex<double> tmp_4830;
   std::complex<double> tmp_4831;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4831 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4830 += tmp_4831;
   result += (-4) * tmp_4830;
   std::complex<double> tmp_4832;
   std::complex<double> tmp_4833;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4833 += B0(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_4832 += tmp_4833;
   result += (-4) * tmp_4832;
   std::complex<double> tmp_4834;
   std::complex<double> tmp_4835;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4835 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_4834 += tmp_4835;
   result += (1.3333333333333333*MGlu) * tmp_4834;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4836;
   std::complex<double> tmp_4837;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4838;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4838 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_4837 += tmp_4838;
   }
   tmp_4836 += tmp_4837;
   result += (-0.5) * tmp_4836;
   std::complex<double> tmp_4839;
   std::complex<double> tmp_4840;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4841;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4841 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_4840 += tmp_4841;
   }
   tmp_4839 += tmp_4840;
   result += (-0.5) * tmp_4839;
   std::complex<double> tmp_4842;
   std::complex<double> tmp_4843;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4844;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4844 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_4843 += tmp_4844;
   }
   tmp_4842 += tmp_4843;
   result += (-0.5) * tmp_4842;
   std::complex<double> tmp_4845;
   std::complex<double> tmp_4846;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4847;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4847 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_4846 += tmp_4847;
   }
   tmp_4845 += tmp_4846;
   result += (-0.5) * tmp_4845;
   std::complex<double> tmp_4848;
   std::complex<double> tmp_4849;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4849 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_4848 += tmp_4849;
   result += (-0.6666666666666666) * tmp_4848;
   std::complex<double> tmp_4850;
   std::complex<double> tmp_4851;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4852;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4852 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_4851 += tmp_4852;
   }
   tmp_4850 += tmp_4851;
   result += (-0.5) * tmp_4850;
   std::complex<double> tmp_4853;
   std::complex<double> tmp_4854;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4854 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_4853 += tmp_4854;
   result += (-1) * tmp_4853;
   std::complex<double> tmp_4855;
   std::complex<double> tmp_4856;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4856 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_4855 += tmp_4856;
   result += (-1) * tmp_4855;
   std::complex<double> tmp_4857;
   std::complex<double> tmp_4858;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4858 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_4857 += tmp_4858;
   result += (-1) * tmp_4857;
   std::complex<double> tmp_4859;
   std::complex<double> tmp_4860;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4860 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPL(gO2,gI2))*
         CpbarUFuVZpFuPL(gO1,gI2);
   }
   tmp_4859 += tmp_4860;
   result += (-1) * tmp_4859;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_4861;
   std::complex<double> tmp_4862;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4863;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4863 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_4862 += tmp_4863;
   }
   tmp_4861 += tmp_4862;
   result += (-0.5) * tmp_4861;
   std::complex<double> tmp_4864;
   std::complex<double> tmp_4865;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_4866;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4866 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_4865 += tmp_4866;
   }
   tmp_4864 += tmp_4865;
   result += (-0.5) * tmp_4864;
   std::complex<double> tmp_4867;
   std::complex<double> tmp_4868;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4869;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4869 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_4868 += tmp_4869;
   }
   tmp_4867 += tmp_4868;
   result += (-0.5) * tmp_4867;
   std::complex<double> tmp_4870;
   std::complex<double> tmp_4871;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_4872;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_4872 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_4871 += tmp_4872;
   }
   tmp_4870 += tmp_4871;
   result += (-0.5) * tmp_4870;
   std::complex<double> tmp_4873;
   std::complex<double> tmp_4874;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4874 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_4873 += tmp_4874;
   result += (-0.6666666666666666) * tmp_4873;
   std::complex<double> tmp_4875;
   std::complex<double> tmp_4876;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_4877;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_4877 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_4876 += tmp_4877;
   }
   tmp_4875 += tmp_4876;
   result += (-0.5) * tmp_4875;
   std::complex<double> tmp_4878;
   std::complex<double> tmp_4879;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4879 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_4878 += tmp_4879;
   result += (-1) * tmp_4878;
   std::complex<double> tmp_4880;
   std::complex<double> tmp_4881;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4881 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_4880 += tmp_4881;
   result += (-1) * tmp_4880;
   std::complex<double> tmp_4882;
   std::complex<double> tmp_4883;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4883 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_4882 += tmp_4883;
   result += (-1) * tmp_4882;
   std::complex<double> tmp_4884;
   std::complex<double> tmp_4885;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_4885 += B1(p,MFu(gI2),MVZp)*Conj(CpbarUFuVZpFuPR(gO2,gI2))*
         CpbarUFuVZpFuPR(gO1,gI2);
   }
   tmp_4884 += tmp_4885;
   result += (-1) * tmp_4884;

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
   std::complex<double> tmp_4886;
   std::complex<double> tmp_4887;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_4887 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_4886 += tmp_4887;
   result += (-1) * tmp_4886;
   std::complex<double> tmp_4888;
   std::complex<double> tmp_4889;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_4889 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_4888 += tmp_4889;
   result += (2) * tmp_4888;
   std::complex<double> tmp_4890;
   std::complex<double> tmp_4891;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4891 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_4890 += tmp_4891;
   result += (-0.5) * tmp_4890;
   std::complex<double> tmp_4892;
   std::complex<double> tmp_4893;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4893 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_4892 += tmp_4893;
   result += (-0.5) * tmp_4892;
   std::complex<double> tmp_4894;
   std::complex<double> tmp_4895;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4895 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_4894 += tmp_4895;
   result += (6) * tmp_4894;
   std::complex<double> tmp_4896;
   std::complex<double> tmp_4897;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4897 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_4896 += tmp_4897;
   result += (2) * tmp_4896;
   std::complex<double> tmp_4898;
   std::complex<double> tmp_4899;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4899 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_4898 += tmp_4899;
   result += (6) * tmp_4898;
   std::complex<double> tmp_4900;
   std::complex<double> tmp_4901;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_4901 += A0(MFv(gI1))*(CpUhhbarFvFvPL(gO1,gI1,gI1) + CpUhhbarFvFvPR
         (gO1,gI1,gI1))*MFv(gI1);
   }
   tmp_4900 += tmp_4901;
   result += (2) * tmp_4900;
   std::complex<double> tmp_4902;
   std::complex<double> tmp_4903;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4903 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_4902 += tmp_4903;
   result += (-3) * tmp_4902;
   std::complex<double> tmp_4904;
   std::complex<double> tmp_4905;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4905 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_4904 += tmp_4905;
   result += (-1) * tmp_4904;
   std::complex<double> tmp_4906;
   std::complex<double> tmp_4907;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4907 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_4906 += tmp_4907;
   result += (-3) * tmp_4906;
   std::complex<double> tmp_4908;
   std::complex<double> tmp_4909;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4909 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_4908 += tmp_4909;
   result += (-1) * tmp_4908;
   std::complex<double> tmp_4910;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_4910 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_4910;

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
   const double M_tree(MGlu);
   const double p = MGlu;
   const double self_energy_1  = Re(self_energy_Glu_1(p));
   const double self_energy_PL = Re(self_energy_Glu_PL(p));
   const double self_energy_PR = Re(self_energy_Glu_PR(p));
   const auto M_1loop = M_tree - self_energy_1 - M_tree * (self_energy_PL
      + self_energy_PR);
   PHYSICAL(MGlu) = calculate_singlet_mass(M_1loop);
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

      PHYSICAL(MSd(es)) = SignedAbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZD) = mix_ZD;
   }
}

void CLASSNAME::calculate_MSv_pole()
{
   if (!force_output && problems.is_tachyon(Sv))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,6,6> self_energy;
   const Eigen::Matrix<double,6,6> M_tree(get_mass_matrix_Sv());

   for (unsigned es = 0; es < 6; ++es) {
      const double p = Abs(MSv(es));
      for (unsigned i1 = 0; i1 < 6; ++i1) {
         for (unsigned i2 = i1; i2 < 6; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Sv(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,6,6> M_1loop(M_tree - self_energy);
      Eigen::Array<double,6,1> eigen_values;
      Eigen::Matrix<double,6,6> mix_ZV;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZV,
            eigenvalue_error);
         problems.flag_bad_mass(UMSSM_info::Sv, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZV);
      #endif

      PHYSICAL(MSv(es)) = SignedAbsSqrt(eigen_values(es));
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

      PHYSICAL(MSu(es)) = SignedAbsSqrt(eigen_values(es));
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

      PHYSICAL(MSe(es)) = SignedAbsSqrt(eigen_values(es));
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

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fv());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFv(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fv_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fv_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fv_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZVL) mix_ZVL;
      decltype(ZVR) mix_ZVR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZVL, mix_ZVR, eigenvalue_error
         );
      problems.flag_bad_mass(UMSSM_info::Fv, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_ZVL, mix_ZVR);
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

   const double qcd_1l = 0.025330295910584444*(-1.6666666666666667 + 1.*
      Log(Sqr(MFu(2))/Sqr(currentScale)))*Sqr(g3);

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


double CLASSNAME::calculate_MFv_DRbar(double m_pole, int idx) const
{
   const double p = m_pole;
   const double self_energy_1  = Re(self_energy_Fv_1(p, idx, idx));
   const double self_energy_PL = Re(self_energy_Fv_PL(p, idx, idx));
   const double self_energy_PR = Re(self_energy_Fv_PR(p, idx, idx));

   const double m_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR);

   return m_drbar;
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
   const double qcd_1l = 0.025330295910584444*(-1.6666666666666667 + 1.*
      Log(Sqr(MFu(idx))/Sqr(currentScale)))*Sqr(g3);
   const double qcd_2l = -0.003408916029785599*Power(g3,4) +
      0.0011495761378943394*Power(g3,4)*Log(Sqr(MFu(idx))/Sqr(currentScale)) -
      0.00024060895909416413*Power(g3,4)*Sqr(Log(Power(MFu(idx),2)/Sqr(
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
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1 - 0.00020496318737651018*Power(g3,4)
      + 0.0006860288475783287*Sqr(g1) + 0.0023747152416172916*Sqr(g2) -
      0.008443431970194815*Sqr(g3);
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
