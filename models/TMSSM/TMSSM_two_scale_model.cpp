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

// File generated at Mon 23 Feb 2015 12:32:43

/**
 * @file TMSSM_two_scale_model.cpp
 * @brief implementation of the TMSSM model class
 *
 * Contains the definition of the TMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 23 Feb 2015 12:32:43 with FlexibleSUSY
 * 1.0.4 (git commit: v1.0.4-341-gb865fd3) and SARAH 4.4.6 .
 */

#include "TMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "gsl_utils.hpp"
#include "config.h"
#include "pv.hpp"



#include <cmath>
#include <iostream>

#ifdef ENABLE_THREADS
#include <thread>
#endif

#include <gsl/gsl_multiroots.h>

namespace flexiblesusy {

using namespace TMSSM_info;

#define CLASSNAME TMSSM<Two_scale>

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

#define HIGGS_2LOOP_CORRECTION_AT_AS     higgs_2loop_corrections.at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS     higgs_2loop_corrections.ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT     higgs_2loop_corrections.at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU higgs_2loop_corrections.atau_atau

#ifdef ENABLE_THREADS
   std::mutex CLASSNAME::mtx_fortran;
   #define LOCK_MUTEX() mtx_fortran.lock()
   #define UNLOCK_MUTEX() mtx_fortran.unlock()
#else
   #define LOCK_MUTEX()
   #define UNLOCK_MUTEX()
#endif

CLASSNAME::TMSSM(const TMSSM_input_parameters& input_)
   : Two_scale_model()
   , TMSSM_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(TMSSM_info::particle_names)
   , higgs_2loop_corrections()
#ifdef ENABLE_THREADS
   , thread_exception()
#endif
   , MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MVZ(0), MSd(Eigen::Array<
      double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>::Zero()), MSu(
      Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()),
      Mhh(Eigen::Array<double,3,1>::Zero()), MAh(Eigen::Array<double,3,1>::Zero())
      , MHpm(Eigen::Array<double,4,1>::Zero()), MChi(Eigen::Array<double,5,1>
      ::Zero()), MCha(Eigen::Array<double,3,1>::Zero()), MFe(Eigen::Array<double,3
      ,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double
      ,3,1>::Zero()), MVG(0), MVP(0), MVWm(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,4,4>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,5,5>::Zero()), UM(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), UP(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZEL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero())

   , PhaseGlu(1,0)

{
}

CLASSNAME::~TMSSM()
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

void CLASSNAME::set_higgs_2loop_corrections(const Higgs_2loop_corrections& higgs_2loop_corrections_)
{
   higgs_2loop_corrections = higgs_2loop_corrections_;
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

const Problems<TMSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<TMSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
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
   if (contains_nan(x, number_of_ewsb_equations)) {
      for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
         gsl_vector_set(f, i, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   TMSSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_BMu(gsl_vector_get(x, 0));
   model->set_Mu(gsl_vector_get(x, 1));
   model->set_mT2(gsl_vector_get(x, 2));


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double tadpole[number_of_ewsb_equations] = { 0. };

   model->tadpole_equations(tadpole);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, tadpole[i]);

   bool is_finite = true;

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      is_finite = is_finite && std::isfinite(tadpole[i]);

   return (is_finite ? GSL_SUCCESS : GSL_EDOM);
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

   double x_init[number_of_ewsb_equations];
   ewsb_initial_guess(x_init);

#ifdef ENABLE_VERBOSE
   std::cout << "Solving EWSB equations ...\n"
      "\tInitial guess: x_init =";
   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      std::cout << " " << x_init[i];
   std::cout << '\n';
#endif

   int status;
   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i) {
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

   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i)
      delete solvers[i];

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

   BMu = solver->get_solution(0);
   Mu = solver->get_solution(1);
   mT2 = solver->get_solution(2);


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

   const double old_BMu = BMu;
   const double old_Mu = Mu;
   const double old_mT2 = mT2;

   BMu = (0.05*(-20*mHd2*vd*vu + 20*mHu2*vd*vu + 5*Power(vd,3)*vu*AbsSqr(
      Lambdax) - 5*vd*Power(vu,3)*AbsSqr(Lambdax) - 3*Power(vd,3)*vu*Sqr(g1) + 3*
      vd*Power(vu,3)*Sqr(g1) - 5*Power(vd,3)*vu*Sqr(g2) + 5*vd*Power(vu,3)*Sqr(g2)
      - 10*MT*vT*Conj(Lambdax)*Sqr(vd) - 5*vT*Conj(TLambdax)*Sqr(vd) - 10*vT*Conj
      (MT)*Lambdax*Sqr(vd) + 10*MT*vT*Conj(Lambdax)*Sqr(vu) + 5*vT*Conj(TLambdax)*
      Sqr(vu) + 10*vT*Conj(MT)*Lambdax*Sqr(vu) - 5*vT*Sqr(vd)*TLambdax + 5*vT*Sqr(
      vu)*TLambdax))/(Sqr(vd) - Sqr(vu));
   Mu = (0.0125*(-20*vd*vT*Conj(Lambdax) - 20*vd*vT*Lambdax + LOCALINPUT(SignMu
      )*Sqrt(Sqr(20*vd*vT*Conj(Lambdax) + 20*vd*vT*Lambdax) - 160*vd*(40*mHd2*vd -
      40*vu*BMu - 20*MT*vT*vu*Conj(Lambdax) - 10*vT*vu*Conj(TLambdax) - 20*vT*vu*
      Conj(MT)*Lambdax + 3*Power(vd,3)*Sqr(g1) + 5*Power(vd,3)*Sqr(g2) + 10*vd*
      AbsSqr(Lambdax)*Sqr(vT) + 10*vd*AbsSqr(Lambdax)*Sqr(vu) - 3*vd*Sqr(g1)*Sqr(
      vu) - 5*vd*Sqr(g2)*Sqr(vu) - 10*vT*vu*TLambdax))))/vd;
   mT2 = (0.25*(-16*vT*AbsSqr(MT) - 4*vT*BMT - 4*vT*Conj(BMT) + 2*MT*vd*vu*Conj
      (Lambdax) + vd*vu*Conj(TLambdax) + 2*vd*vu*Conj(MT)*Lambdax - vT*AbsSqr(
      Lambdax)*Sqr(vd) - Conj(Lambdax)*Mu*Sqr(vd) - Lambdax*Mu*Sqr(vd) - vT*AbsSqr
      (Lambdax)*Sqr(vu) - Conj(Lambdax)*Mu*Sqr(vu) - Lambdax*Mu*Sqr(vu) + vd*vu*
      TLambdax))/vT;

   const bool is_finite = std::isfinite(BMu) && std::isfinite(Mu) &&
      std::isfinite(mT2);

   if (!is_finite) {
      BMu = old_BMu;
      Mu = old_Mu;
      mT2 = old_mT2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_via_soft_higgs_masses()
{
   int error = 0;

   const double new_mHd2 = (0.025*(-40*vd*AbsSqr(Mu) + 20*vu*BMu + 20*vu*Conj(
      BMu) + 20*MT*vT*vu*Conj(Lambdax) + 10*vT*vu*Conj(TLambdax) + 20*vT*vu*Conj(
      MT)*Lambdax - 20*vd*vT*Conj(Mu)*Lambdax - 20*vd*vT*Conj(Lambdax)*Mu - 3*
      Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 10*vd*AbsSqr(Lambdax)*Sqr(vT)
      - 10*vd*AbsSqr(Lambdax)*Sqr(vu) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*Sqr(g2)*Sqr(vu
      ) + 10*vT*vu*TLambdax))/vd;
   const double new_mHu2 = (0.025*(-40*vu*AbsSqr(Mu) + 20*vd*BMu + 20*vd*Conj(
      BMu) + 20*MT*vd*vT*Conj(Lambdax) + 10*vd*vT*Conj(TLambdax) + 20*vd*vT*Conj(
      MT)*Lambdax - 20*vT*vu*Conj(Mu)*Lambdax - 20*vT*vu*Conj(Lambdax)*Mu - 3*
      Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) - 10*vu*AbsSqr(Lambdax)*Sqr(vd)
      + 3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 10*vu*AbsSqr(Lambdax)*Sqr(vT
      ) + 10*vd*vT*TLambdax))/vu;
   const double new_mT2 = (0.25*(-16*vT*AbsSqr(MT) - 4*vT*BMT - 4*vT*Conj(BMT)
      + 2*MT*vd*vu*Conj(Lambdax) + vd*vu*Conj(TLambdax) + 2*vd*vu*Conj(MT)*Lambdax
      - vT*AbsSqr(Lambdax)*Sqr(vd) - Conj(Mu)*Lambdax*Sqr(vd) - Conj(Lambdax)*Mu*
      Sqr(vd) - vT*AbsSqr(Lambdax)*Sqr(vu) - Conj(Mu)*Lambdax*Sqr(vu) - Conj(
      Lambdax)*Mu*Sqr(vu) + vd*vu*TLambdax))/vT;

   if (std::isfinite(new_mHd2))
      mHd2 = new_mHd2;
   else
      error = 1;

   if (std::isfinite(new_mHu2))
      mHu2 = new_mHu2;
   else
      error = 1;

   if (std::isfinite(new_mT2))
      mT2 = new_mT2;
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
   x_init[0] = BMu;
   x_init[1] = Mu;
   x_init[2] = mT2;

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

      }
   }

   double BMu;
   double Mu;
   double mT2;

   BMu = (0.05*(-20*mHd2*vd*vu + 20*mHu2*vd*vu + 5*Power(vd,3)*vu*AbsSqr(
      Lambdax) - 5*vd*Power(vu,3)*AbsSqr(Lambdax) + 20*vu*tadpole[0] - 20*vd*
      tadpole[1] - 3*Power(vd,3)*vu*Sqr(g1) + 3*vd*Power(vu,3)*Sqr(g1) - 5*Power(
      vd,3)*vu*Sqr(g2) + 5*vd*Power(vu,3)*Sqr(g2) - 10*MT*vT*Conj(Lambdax)*Sqr(vd)
      - 5*vT*Conj(TLambdax)*Sqr(vd) - 10*vT*Conj(MT)*Lambdax*Sqr(vd) + 10*MT*vT*
      Conj(Lambdax)*Sqr(vu) + 5*vT*Conj(TLambdax)*Sqr(vu) + 10*vT*Conj(MT)*Lambdax
      *Sqr(vu) - 5*vT*Sqr(vd)*TLambdax + 5*vT*Sqr(vu)*TLambdax))/(Sqr(vd) - Sqr(vu
      ));
   Mu = (0.0125*(-20*vd*vT*Conj(Lambdax) - 20*vd*vT*Lambdax + LOCALINPUT(SignMu
      )*Sqrt(Sqr(20*vd*vT*Conj(Lambdax) + 20*vd*vT*Lambdax) - 160*vd*(40*mHd2*vd -
      40*vu*BMu - 20*MT*vT*vu*Conj(Lambdax) - 10*vT*vu*Conj(TLambdax) - 20*vT*vu*
      Conj(MT)*Lambdax - 40*tadpole[0] + 3*Power(vd,3)*Sqr(g1) + 5*Power(vd,3)*Sqr
      (g2) + 10*vd*AbsSqr(Lambdax)*Sqr(vT) + 10*vd*AbsSqr(Lambdax)*Sqr(vu) - 3*vd*
      Sqr(g1)*Sqr(vu) - 5*vd*Sqr(g2)*Sqr(vu) - 10*vT*vu*TLambdax))))/vd;
   mT2 = (0.25*(-16*vT*AbsSqr(MT) - 4*vT*BMT - 4*vT*Conj(BMT) + 2*MT*vd*vu*Conj
      (Lambdax) + vd*vu*Conj(TLambdax) + 2*vd*vu*Conj(MT)*Lambdax + 4*tadpole[2] -
      vT*AbsSqr(Lambdax)*Sqr(vd) - Conj(Lambdax)*Mu*Sqr(vd) - Lambdax*Mu*Sqr(vd)
      - vT*AbsSqr(Lambdax)*Sqr(vu) - Conj(Lambdax)*Mu*Sqr(vu) - Lambdax*Mu*Sqr(vu)
      + vd*vu*TLambdax))/vT;

   const bool is_finite = std::isfinite(BMu) && std::isfinite(Mu) &&
      std::isfinite(mT2);


   if (is_finite) {
      error = GSL_SUCCESS;
      ewsb_parameters[0] = BMu;
      ewsb_parameters[1] = Mu;
      ewsb_parameters[2] = mT2;

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
   if (contains_nan(x, number_of_ewsb_equations)) {
      for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
         gsl_vector_set(f, i, std::numeric_limits<double>::max());
      return GSL_EDOM;
   }

   const CLASSNAME::EWSB_args* ewsb_args
      = static_cast<CLASSNAME::EWSB_args*>(params);
   TMSSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double BMu = gsl_vector_get(x, 0);
   const double Mu = gsl_vector_get(x, 1);
   const double mT2 = gsl_vector_get(x, 2);

   model->set_BMu(BMu);
   model->set_Mu(Mu);
   model->set_mT2(mT2);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { BMu, Mu, mT2 };

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
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
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MVZ = " << MVZ << '\n';
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
   ostr << "MVG = " << MVG << '\n';
   ostr << "MVP = " << MVP << '\n';
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
   const auto old_mT2 = mT2;

   solve_ewsb_tree_level_via_soft_higgs_masses();

   calculate_MVZ();
   calculate_MVG();
   calculate_MVP();
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
   mT2 = old_mT2;

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

   if (calculate_sm_pole_masses) {
      std::thread thread_MFv(Thread(this, &CLASSNAME::calculate_MFv_pole));
      std::thread thread_MVZ(Thread(this, &CLASSNAME::calculate_MVZ_pole));
      std::thread thread_MFe(Thread(this, &CLASSNAME::calculate_MFe_pole));
      std::thread thread_MFd(Thread(this, &CLASSNAME::calculate_MFd_pole));
      std::thread thread_MFu(Thread(this, &CLASSNAME::calculate_MFu_pole));
      std::thread thread_MVG(Thread(this, &CLASSNAME::calculate_MVG_pole));
      std::thread thread_MVP(Thread(this, &CLASSNAME::calculate_MVP_pole));
      std::thread thread_MVWm(Thread(this, &CLASSNAME::calculate_MVWm_pole));
      thread_MFv.join();
      thread_MVZ.join();
      thread_MFe.join();
      thread_MFd.join();
      thread_MFu.join();
      thread_MVG.join();
      thread_MVP.join();
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

   if (calculate_sm_pole_masses) {
      calculate_MFv_pole();
      calculate_MVZ_pole();
      calculate_MFe_pole();
      calculate_MFd_pole();
      calculate_MFu_pole();
      calculate_MVG_pole();
      calculate_MVP_pole();
      calculate_MVWm_pole();
   }

#endif
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(MVZ) = MVZ;
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
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MVP) = MVP;
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
   move_goldstone_to(0, MVWm, MHpm, ZP);

}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associuated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{
   move_goldstone_to(0, MVZ, PHYSICAL(MAh), PHYSICAL(ZA));
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
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MVZ = 0.;
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
   MVG = 0.;
   MVP = 0.;
   MVWm = 0.;

   PhaseGlu = std::complex<double>(1.,0.);

}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   TMSSM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
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

Eigen::Array<double,3,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,3,1> MHpm_ChargedHiggs;
   Eigen::Array<double,1,1> MHpm_goldstone;

   MHpm_goldstone(0) = MVWm;

   remove_if_equal(MHpm, MHpm_goldstone, MHpm_ChargedHiggs);

   return MHpm_ChargedHiggs;
}

Eigen::Array<double,2,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,2,1> MAh_PseudoscalarHiggs;
   Eigen::Array<double,1,1> MAh_goldstone;

   MAh_goldstone(0) = MVZ;

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


double CLASSNAME::get_mass_matrix_Glu() const
{
   const double mass_matrix_Glu = MassG;

   return mass_matrix_Glu;
}

void CLASSNAME::calculate_MGlu()
{
   MGlu = get_mass_matrix_Glu();

   if (MGlu < 0.) {
      MGlu *= -1;
      PhaseGlu = std::complex<double>(0., 1.);
   }
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Fv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,0) = mass_matrix_Fv(0,1);
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,0) = mass_matrix_Fv(0,2);
   mass_matrix_Fv(2,1) = mass_matrix_Fv(1,2);
   mass_matrix_Fv(2,2) = 0;

   return mass_matrix_Fv;
}

void CLASSNAME::calculate_MFv()
{
   MFv.setConstant(0);
}

double CLASSNAME::get_mass_matrix_VZ() const
{
   const double mass_matrix_VZ = 0.25*(Sqr(vd) + Sqr(vu))*Sqr(g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return mass_matrix_VZ;
}

void CLASSNAME::calculate_MVZ()
{
   MVZ = get_mass_matrix_VZ();

   if (MVZ < 0.)
      problems.flag_tachyon(TMSSM_info::VZ);

   MVZ = AbsSqrt(MVZ);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Sd() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(1,0)
      ) + AbsSqr(Yd(2,0)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(0,1) = mq2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,1) +
      Conj(Yd(1,0))*Yd(1,1) + Conj(Yd(2,0))*Yd(2,1));
   mass_matrix_Sd(0,2) = mq2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(0,0))*Yd(0,2) +
      Conj(Yd(1,0))*Yd(1,2) + Conj(Yd(2,0))*Yd(2,2));
   mass_matrix_Sd(0,3) = 0.7071067811865475*vd*Conj(TYd(0,0)) -
      0.35355339059327373*vT*vu*Conj(Yd(0,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(0,0))*Mu;
   mass_matrix_Sd(0,4) = 0.7071067811865475*vd*Conj(TYd(1,0)) -
      0.35355339059327373*vT*vu*Conj(Yd(1,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(1,0))*Mu;
   mass_matrix_Sd(0,5) = 0.7071067811865475*vd*Conj(TYd(2,0)) -
      0.35355339059327373*vT*vu*Conj(Yd(2,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(2,0))*Mu;
   mass_matrix_Sd(1,0) = Conj(mass_matrix_Sd(0,1));
   mass_matrix_Sd(1,1) = mq2(1,1) + 0.5*(AbsSqr(Yd(0,1)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(2,1)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(1,2) = mq2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(0,1))*Yd(0,2) +
      Conj(Yd(1,1))*Yd(1,2) + Conj(Yd(2,1))*Yd(2,2));
   mass_matrix_Sd(1,3) = 0.7071067811865475*vd*Conj(TYd(0,1)) -
      0.35355339059327373*vT*vu*Conj(Yd(0,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(0,1))*Mu;
   mass_matrix_Sd(1,4) = 0.7071067811865475*vd*Conj(TYd(1,1)) -
      0.35355339059327373*vT*vu*Conj(Yd(1,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(1,1))*Mu;
   mass_matrix_Sd(1,5) = 0.7071067811865475*vd*Conj(TYd(2,1)) -
      0.35355339059327373*vT*vu*Conj(Yd(2,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(2,1))*Mu;
   mass_matrix_Sd(2,0) = Conj(mass_matrix_Sd(0,2));
   mass_matrix_Sd(2,1) = Conj(mass_matrix_Sd(1,2));
   mass_matrix_Sd(2,2) = mq2(2,2) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) + 0.025*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = 0.7071067811865475*vd*Conj(TYd(0,2)) -
      0.35355339059327373*vT*vu*Conj(Yd(0,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = 0.7071067811865475*vd*Conj(TYd(1,2)) -
      0.35355339059327373*vT*vu*Conj(Yd(1,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = 0.7071067811865475*vd*Conj(TYd(2,2)) -
      0.35355339059327373*vT*vu*Conj(Yd(2,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,0) = Conj(mass_matrix_Sd(0,3));
   mass_matrix_Sd(3,1) = Conj(mass_matrix_Sd(1,3));
   mass_matrix_Sd(3,2) = Conj(mass_matrix_Sd(2,3));
   mass_matrix_Sd(3,3) = md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) + AbsSqr(Yd(0,1)
      ) + AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) +
      Conj(Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) +
      Conj(Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,0) = Conj(mass_matrix_Sd(0,4));
   mass_matrix_Sd(4,1) = Conj(mass_matrix_Sd(1,4));
   mass_matrix_Sd(4,2) = Conj(mass_matrix_Sd(2,4));
   mass_matrix_Sd(4,3) = Conj(mass_matrix_Sd(3,4));
   mass_matrix_Sd(4,4) = md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) + AbsSqr(Yd(1,1)
      ) + AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) +
      Conj(Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,0) = Conj(mass_matrix_Sd(0,5));
   mass_matrix_Sd(5,1) = Conj(mass_matrix_Sd(1,5));
   mass_matrix_Sd(5,2) = Conj(mass_matrix_Sd(2,5));
   mass_matrix_Sd(5,3) = Conj(mass_matrix_Sd(3,5));
   mass_matrix_Sd(5,4) = Conj(mass_matrix_Sd(4,5));
   mass_matrix_Sd(5,5) = md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) + AbsSqr(Yd(2,1)
      ) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*Sqr(g1)*Sqr(vu
      );

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Sd, eigenvalue_error > precision *
      Abs(MSd(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
#endif

   if (MSd.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::Sd);

   MSd = AbsSqrt(MSd);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Sv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Sv;

   mass_matrix_Sv(0,0) = ml2(0,0) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(0,1) = ml2(0,1);
   mass_matrix_Sv(0,2) = ml2(0,2);
   mass_matrix_Sv(1,0) = Conj(mass_matrix_Sv(0,1));
   mass_matrix_Sv(1,1) = ml2(1,1) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,0) = Conj(mass_matrix_Sv(0,2));
   mass_matrix_Sv(2,1) = Conj(mass_matrix_Sv(1,2));
   mass_matrix_Sv(2,2) = ml2(2,2) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);

   return mass_matrix_Sv;
}

void CLASSNAME::calculate_MSv()
{
   const auto mass_matrix_Sv(get_mass_matrix_Sv());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Sv, eigenvalue_error > precision *
      Abs(MSv(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Sv, MSv, ZV);
#endif

   if (MSv.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::Sv);

   MSv = AbsSqrt(MSv);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Su() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(1,0)) + AbsSqr(Yu(2,0)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(0,1) = mq2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,1) +
      Conj(Yu(1,0))*Yu(1,1) + Conj(Yu(2,0))*Yu(2,1));
   mass_matrix_Su(0,2) = mq2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(0,0))*Yu(0,2) +
      Conj(Yu(1,0))*Yu(1,2) + Conj(Yu(2,0))*Yu(2,2));
   mass_matrix_Su(0,3) = 0.7071067811865475*vu*Conj(TYu(0,0)) -
      0.35355339059327373*vd*vT*Conj(Yu(0,0))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(0,0))*Mu;
   mass_matrix_Su(0,4) = 0.7071067811865475*vu*Conj(TYu(1,0)) -
      0.35355339059327373*vd*vT*Conj(Yu(1,0))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(1,0))*Mu;
   mass_matrix_Su(0,5) = 0.7071067811865475*vu*Conj(TYu(2,0)) -
      0.35355339059327373*vd*vT*Conj(Yu(2,0))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(2,0))*Mu;
   mass_matrix_Su(1,0) = Conj(mass_matrix_Su(0,1));
   mass_matrix_Su(1,1) = mq2(1,1) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,1)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(2,1)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(1,2) = mq2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(0,1))*Yu(0,2) +
      Conj(Yu(1,1))*Yu(1,2) + Conj(Yu(2,1))*Yu(2,2));
   mass_matrix_Su(1,3) = 0.7071067811865475*vu*Conj(TYu(0,1)) -
      0.35355339059327373*vd*vT*Conj(Yu(0,1))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(0,1))*Mu;
   mass_matrix_Su(1,4) = 0.7071067811865475*vu*Conj(TYu(1,1)) -
      0.35355339059327373*vd*vT*Conj(Yu(1,1))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(1,1))*Mu;
   mass_matrix_Su(1,5) = 0.7071067811865475*vu*Conj(TYu(2,1)) -
      0.35355339059327373*vd*vT*Conj(Yu(2,1))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(2,1))*Mu;
   mass_matrix_Su(2,0) = Conj(mass_matrix_Su(0,2));
   mass_matrix_Su(2,1) = Conj(mass_matrix_Su(1,2));
   mass_matrix_Su(2,2) = mq2(2,2) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)
      *Sqr(vd) + 0.5*(AbsSqr(Yu(0,2)) + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(
      vu) + 0.025*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = 0.7071067811865475*vu*Conj(TYu(0,2)) -
      0.35355339059327373*vd*vT*Conj(Yu(0,2))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = 0.7071067811865475*vu*Conj(TYu(1,2)) -
      0.35355339059327373*vd*vT*Conj(Yu(1,2))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = 0.7071067811865475*vu*Conj(TYu(2,2)) -
      0.35355339059327373*vd*vT*Conj(Yu(2,2))*Lambdax - 0.7071067811865475*vd*
      Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,0) = Conj(mass_matrix_Su(0,3));
   mass_matrix_Su(3,1) = Conj(mass_matrix_Su(1,3));
   mass_matrix_Su(3,2) = Conj(mass_matrix_Su(2,3));
   mass_matrix_Su(3,3) = mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      0,0)) + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) +
      Conj(Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) +
      Conj(Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,0) = Conj(mass_matrix_Su(0,4));
   mass_matrix_Su(4,1) = Conj(mass_matrix_Su(1,4));
   mass_matrix_Su(4,2) = Conj(mass_matrix_Su(2,4));
   mass_matrix_Su(4,3) = Conj(mass_matrix_Su(3,4));
   mass_matrix_Su(4,4) = mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      1,0)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) +
      Conj(Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,0) = Conj(mass_matrix_Su(0,5));
   mass_matrix_Su(5,1) = Conj(mass_matrix_Su(1,5));
   mass_matrix_Su(5,2) = Conj(mass_matrix_Su(2,5));
   mass_matrix_Su(5,3) = Conj(mass_matrix_Su(3,5));
   mass_matrix_Su(5,4) = Conj(mass_matrix_Su(4,5));
   mass_matrix_Su(5,5) = mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) + 0.5*(AbsSqr(Yu(
      2,0)) + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*Sqr(g1)*Sqr(vu);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Su, eigenvalue_error > precision *
      Abs(MSu(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
#endif

   if (MSu.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::Su);

   MSu = AbsSqrt(MSu);
}

Eigen::Matrix<double,6,6> CLASSNAME::get_mass_matrix_Se() const
{
   Eigen::Matrix<double,6,6> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(1,0)
      ) + AbsSqr(Ye(2,0)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(0,1) = ml2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,1) +
      Conj(Ye(1,0))*Ye(1,1) + Conj(Ye(2,0))*Ye(2,1));
   mass_matrix_Se(0,2) = ml2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(0,0))*Ye(0,2) +
      Conj(Ye(1,0))*Ye(1,2) + Conj(Ye(2,0))*Ye(2,2));
   mass_matrix_Se(0,3) = 0.7071067811865475*vd*Conj(TYe(0,0)) -
      0.35355339059327373*vT*vu*Conj(Ye(0,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(0,0))*Mu;
   mass_matrix_Se(0,4) = 0.7071067811865475*vd*Conj(TYe(1,0)) -
      0.35355339059327373*vT*vu*Conj(Ye(1,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(1,0))*Mu;
   mass_matrix_Se(0,5) = 0.7071067811865475*vd*Conj(TYe(2,0)) -
      0.35355339059327373*vT*vu*Conj(Ye(2,0))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(2,0))*Mu;
   mass_matrix_Se(1,0) = Conj(mass_matrix_Se(0,1));
   mass_matrix_Se(1,1) = ml2(1,1) + 0.5*(AbsSqr(Ye(0,1)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(2,1)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(1,2) = ml2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(0,1))*Ye(0,2) +
      Conj(Ye(1,1))*Ye(1,2) + Conj(Ye(2,1))*Ye(2,2));
   mass_matrix_Se(1,3) = 0.7071067811865475*vd*Conj(TYe(0,1)) -
      0.35355339059327373*vT*vu*Conj(Ye(0,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(0,1))*Mu;
   mass_matrix_Se(1,4) = 0.7071067811865475*vd*Conj(TYe(1,1)) -
      0.35355339059327373*vT*vu*Conj(Ye(1,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(1,1))*Mu;
   mass_matrix_Se(1,5) = 0.7071067811865475*vd*Conj(TYe(2,1)) -
      0.35355339059327373*vT*vu*Conj(Ye(2,1))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(2,1))*Mu;
   mass_matrix_Se(2,0) = Conj(mass_matrix_Se(0,2));
   mass_matrix_Se(2,1) = Conj(mass_matrix_Se(1,2));
   mass_matrix_Se(2,2) = ml2(2,2) + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(
      vd) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = 0.7071067811865475*vd*Conj(TYe(0,2)) -
      0.35355339059327373*vT*vu*Conj(Ye(0,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = 0.7071067811865475*vd*Conj(TYe(1,2)) -
      0.35355339059327373*vT*vu*Conj(Ye(1,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = 0.7071067811865475*vd*Conj(TYe(2,2)) -
      0.35355339059327373*vT*vu*Conj(Ye(2,2))*Lambdax - 0.7071067811865475*vu*
      Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,0) = Conj(mass_matrix_Se(0,3));
   mass_matrix_Se(3,1) = Conj(mass_matrix_Se(1,3));
   mass_matrix_Se(3,2) = Conj(mass_matrix_Se(2,3));
   mass_matrix_Se(3,3) = me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) + AbsSqr(Ye(0,1)
      ) + AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) +
      Conj(Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) +
      Conj(Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,0) = Conj(mass_matrix_Se(0,4));
   mass_matrix_Se(4,1) = Conj(mass_matrix_Se(1,4));
   mass_matrix_Se(4,2) = Conj(mass_matrix_Se(2,4));
   mass_matrix_Se(4,3) = Conj(mass_matrix_Se(3,4));
   mass_matrix_Se(4,4) = me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) + AbsSqr(Ye(1,1)
      ) + AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) +
      Conj(Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,0) = Conj(mass_matrix_Se(0,5));
   mass_matrix_Se(5,1) = Conj(mass_matrix_Se(1,5));
   mass_matrix_Se(5,2) = Conj(mass_matrix_Se(2,5));
   mass_matrix_Se(5,3) = Conj(mass_matrix_Se(3,5));
   mass_matrix_Se(5,4) = Conj(mass_matrix_Se(4,5));
   mass_matrix_Se(5,5) = me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) + AbsSqr(Ye(2,1)
      ) + AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*Sqr(g1)*Sqr(vu
      );

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Se, eigenvalue_error > precision *
      Abs(MSe(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
#endif

   if (MSe.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::Se);

   MSe = AbsSqrt(MSe);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,3,3> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax +
      0.5*vT*Conj(Lambdax)*Mu + 0.225*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.25*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(0,1) = 0.5*vd*vu*AbsSqr(Lambdax) - 0.5*BMu - 0.5*Conj(
      BMu) - 0.5*MT*vT*Conj(Lambdax) - 0.25*vT*Conj(TLambdax) - 0.5*vT*Conj(MT)
      *Lambdax - 0.15*vd*vu*Sqr(g1) - 0.25*vd*vu*Sqr(g2) - 0.25*vT*TLambdax;
   mass_matrix_hh(0,2) = 0.5*vd*vT*AbsSqr(Lambdax) - 0.5*MT*vu*Conj(
      Lambdax) - 0.25*vu*Conj(TLambdax) - 0.5*vu*Conj(MT)*Lambdax + 0.5*vd*Conj
      (Mu)*Lambdax + 0.5*vd*Conj(Lambdax)*Mu - 0.25*vu*TLambdax;
   mass_matrix_hh(1,0) = mass_matrix_hh(0,1);
   mass_matrix_hh(1,1) = mHu2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax +
      0.5*vT*Conj(Lambdax)*Mu + 0.25*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*
      Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.225*
      Sqr(g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_hh(1,2) = 0.5*vT*vu*AbsSqr(Lambdax) - 0.5*MT*vd*Conj(
      Lambdax) - 0.25*vd*Conj(TLambdax) - 0.5*vd*Conj(MT)*Lambdax + 0.5*vu*Conj
      (Mu)*Lambdax + 0.5*vu*Conj(Lambdax)*Mu - 0.25*vd*TLambdax;
   mass_matrix_hh(2,0) = mass_matrix_hh(0,2);
   mass_matrix_hh(2,1) = mass_matrix_hh(1,2);
   mass_matrix_hh(2,2) = mT2 + 4*AbsSqr(MT) + BMT + Conj(BMT) + 0.25*
      AbsSqr(Lambdax)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vu);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::hh, eigenvalue_error > precision *
      Abs(Mhh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);
#endif

   if (Mhh.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::hh);

   Mhh = AbsSqrt(Mhh);
}

Eigen::Matrix<double,3,3> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax +
      0.5*vT*Conj(Lambdax)*Mu + 0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(
      ThetaW())*Sqr(vd) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.25*
      AbsSqr(Lambdax)*Sqr(vT) + 0.25*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(g1)*
      Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vd)*Sqr(Cos(ThetaW()))
      + 0.15*Sqr(g1)*Sqr(vd)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*Conj(BMu) + 0.5*MT*vT*Conj(Lambdax
      ) + 0.25*vT*Conj(TLambdax) + 0.5*vT*Conj(MT)*Lambdax - 0.3872983346207417
      *g1*g2*vd*vu*Cos(ThetaW())*Sin(ThetaW()) - 0.25*vd*vu*Sqr(g2)*Sqr(Cos(
      ThetaW())) - 0.15*vd*vu*Sqr(g1)*Sqr(Sin(ThetaW())) + 0.25*vT*TLambdax;
   mass_matrix_Ah(0,2) = -0.5*MT*vu*Conj(Lambdax) + 0.25*vu*Conj(TLambdax
      ) - 0.5*vu*Conj(MT)*Lambdax + 0.25*vu*TLambdax;
   mass_matrix_Ah(1,0) = mass_matrix_Ah(0,1);
   mass_matrix_Ah(1,1) = mHu2 + AbsSqr(Mu) + 0.5*vT*Conj(Mu)*Lambdax +
      0.5*vT*Conj(Lambdax)*Mu + 0.25*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*
      Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vT) +
      0.3872983346207417*g1*g2*Cos(ThetaW())*Sin(ThetaW())*Sqr(vu) + 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu) + 0.25*Sqr(g2)*Sqr(vu)*Sqr(Cos(ThetaW
      ())) + 0.15*Sqr(g1)*Sqr(vu)*Sqr(Sin(ThetaW()));
   mass_matrix_Ah(1,2) = -0.5*MT*vd*Conj(Lambdax) + 0.25*vd*Conj(TLambdax
      ) - 0.5*vd*Conj(MT)*Lambdax + 0.25*vd*TLambdax;
   mass_matrix_Ah(2,0) = mass_matrix_Ah(0,2);
   mass_matrix_Ah(2,1) = mass_matrix_Ah(1,2);
   mass_matrix_Ah(2,2) = mT2 + 4*AbsSqr(MT) - BMT - Conj(BMT) + 0.25*
      AbsSqr(Lambdax)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vu);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Ah, eigenvalue_error > precision *
      Abs(MAh(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);
#endif

   if (MAh.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::Ah);

   MAh = AbsSqrt(MAh);
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Hpm() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + AbsSqr(Mu) - 0.5*vT*Conj(Mu)*Lambdax -
      0.5*vT*Conj(Lambdax)*Mu + 0.075*Sqr(g1)*Sqr(vd) + 0.375*Sqr(g2)*Sqr(vd) +
      0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.075*Sqr(
      g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(0,1) = 0.25*vd*vu*AbsSqr(Lambdax) + Conj(BMu) - MT*vT*
      Conj(Lambdax) - 0.5*vT*Conj(TLambdax);
   mass_matrix_Hpm(0,2) = -0.35355339059327373*vd*vT*AbsSqr(Lambdax) -
      1.4142135623730951*MT*vu*Conj(Lambdax) + 0.7071067811865475*vd*Conj(Mu)*
      Lambdax + 0.7071067811865475*vd*vT*Sqr(g2);
   mass_matrix_Hpm(0,3) = 0.35355339059327373*vd*vT*AbsSqr(Lambdax) -
      0.7071067811865475*vu*Conj(TLambdax) + 0.7071067811865475*vd*Conj(Lambdax
      )*Mu;
   mass_matrix_Hpm(1,0) = Conj(mass_matrix_Hpm(0,1));
   mass_matrix_Hpm(1,1) = mHu2 + AbsSqr(Mu) - 0.5*vT*Conj(Mu)*Lambdax -
      0.5*vT*Conj(Lambdax)*Mu + 0.5*AbsSqr(Lambdax)*Sqr(vd) - 0.075*Sqr(g1)*Sqr
      (vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.25*AbsSqr(Lambdax)*Sqr(vT) + 0.075*Sqr(
      g1)*Sqr(vu) + 0.375*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(1,2) = -0.35355339059327373*vT*vu*AbsSqr(Lambdax) -
      0.7071067811865475*vu*Conj(Mu)*Lambdax + 0.7071067811865475*vd*TLambdax;
   mass_matrix_Hpm(1,3) = 0.35355339059327373*vT*vu*AbsSqr(Lambdax) +
      1.4142135623730951*vd*Conj(MT)*Lambdax - 0.7071067811865475*vu*Conj(
      Lambdax)*Mu - 0.7071067811865475*vT*vu*Sqr(g2);
   mass_matrix_Hpm(2,0) = Conj(mass_matrix_Hpm(0,2));
   mass_matrix_Hpm(2,1) = Conj(mass_matrix_Hpm(1,2));
   mass_matrix_Hpm(2,2) = mT2 + 4*AbsSqr(MT) + 0.5*AbsSqr(Lambdax)*Sqr(vd
      ) - 0.25*Sqr(g2)*Sqr(vd) + Sqr(g2)*Sqr(vT) + 0.25*Sqr(g2)*Sqr(vu);
   mass_matrix_Hpm(2,3) = 2*Conj(BMT);
   mass_matrix_Hpm(3,0) = Conj(mass_matrix_Hpm(0,3));
   mass_matrix_Hpm(3,1) = Conj(mass_matrix_Hpm(1,3));
   mass_matrix_Hpm(3,2) = Conj(mass_matrix_Hpm(2,3));
   mass_matrix_Hpm(3,3) = mT2 + 4*AbsSqr(MT) + 0.25*Sqr(g2)*Sqr(vd) + Sqr
      (g2)*Sqr(vT) + 0.5*AbsSqr(Lambdax)*Sqr(vu) - 0.25*Sqr(g2)*Sqr(vu);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Hpm, eigenvalue_error > precision *
      Abs(MHpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);
#endif

   if (MHpm.minCoeff() < 0.)
      problems.flag_tachyon(TMSSM_info::Hpm);

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
   mass_matrix_Chi(1,0) = mass_matrix_Chi(0,1);
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(1,4) = 0;
   mass_matrix_Chi(2,0) = mass_matrix_Chi(0,2);
   mass_matrix_Chi(2,1) = mass_matrix_Chi(1,2);
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -0.5*vT*Lambdax - Mu;
   mass_matrix_Chi(2,4) = -0.5*vu*Lambdax;
   mass_matrix_Chi(3,0) = mass_matrix_Chi(0,3);
   mass_matrix_Chi(3,1) = mass_matrix_Chi(1,3);
   mass_matrix_Chi(3,2) = mass_matrix_Chi(2,3);
   mass_matrix_Chi(3,3) = 0;
   mass_matrix_Chi(3,4) = -0.5*vd*Lambdax;
   mass_matrix_Chi(4,0) = mass_matrix_Chi(0,4);
   mass_matrix_Chi(4,1) = mass_matrix_Chi(1,4);
   mass_matrix_Chi(4,2) = mass_matrix_Chi(2,4);
   mass_matrix_Chi(4,3) = mass_matrix_Chi(3,4);
   mass_matrix_Chi(4,4) = 2*MT;

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN, eigenvalue_error);
   problems.flag_bad_mass(TMSSM_info::Chi, eigenvalue_error > precision *
      Abs(MChi(0)));
#else
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
#endif
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
   problems.flag_bad_mass(TMSSM_info::Cha, eigenvalue_error > precision *
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
   problems.flag_bad_mass(TMSSM_info::Fe, eigenvalue_error > precision *
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
   problems.flag_bad_mass(TMSSM_info::Fd, eigenvalue_error > precision *
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
   problems.flag_bad_mass(TMSSM_info::Fu, eigenvalue_error > precision *
      Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, ZUL, ZUR);
#endif
}

double CLASSNAME::get_mass_matrix_VG() const
{
   const double mass_matrix_VG = 0;

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{
   MVG = get_mass_matrix_VG();
}

double CLASSNAME::get_mass_matrix_VP() const
{
   const double mass_matrix_VP = 0;

   return mass_matrix_VP;
}

void CLASSNAME::calculate_MVP()
{
   MVP = get_mass_matrix_VP();
}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = 0.25*Sqr(g2)*(Sqr(vd) + 4*Sqr(vT) + Sqr
      (vu));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   MVWm = get_mass_matrix_VWm();

   if (MVWm < 0.)
      problems.flag_tachyon(TMSSM_info::VWm);

   MVWm = AbsSqrt(MVWm);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = mHd2*vd + vd*AbsSqr(Mu) - 0.5*vu*BMu - 0.5*vu*Conj(BMu) -
      0.5*MT*vT*vu*Conj(Lambdax) - 0.25*vT*vu*Conj(TLambdax) - 0.5*vT*vu*Conj(MT)*
      Lambdax + 0.5*vd*vT*Conj(Mu)*Lambdax + 0.5*vd*vT*Conj(Lambdax)*Mu + 0.075*
      Power(vd,3)*Sqr(g1) + 0.125*Power(vd,3)*Sqr(g2) + 0.25*vd*AbsSqr(Lambdax)*
      Sqr(vT) + 0.25*vd*AbsSqr(Lambdax)*Sqr(vu) - 0.075*vd*Sqr(g1)*Sqr(vu) - 0.125
      *vd*Sqr(g2)*Sqr(vu) - 0.25*vT*vu*TLambdax;

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = mHu2*vu + vu*AbsSqr(Mu) - 0.5*vd*BMu - 0.5*vd*Conj(BMu) -
      0.5*MT*vd*vT*Conj(Lambdax) - 0.25*vd*vT*Conj(TLambdax) - 0.5*vd*vT*Conj(MT)*
      Lambdax + 0.5*vT*vu*Conj(Mu)*Lambdax + 0.5*vT*vu*Conj(Lambdax)*Mu + 0.075*
      Power(vu,3)*Sqr(g1) + 0.125*Power(vu,3)*Sqr(g2) + 0.25*vu*AbsSqr(Lambdax)*
      Sqr(vd) - 0.075*vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.25*vu*
      AbsSqr(Lambdax)*Sqr(vT) - 0.25*vd*vT*TLambdax;

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   double result = mT2*vT + 4*vT*AbsSqr(MT) + vT*BMT + vT*Conj(BMT) - 0.5*MT*vd
      *vu*Conj(Lambdax) - 0.25*vd*vu*Conj(TLambdax) - 0.5*vd*vu*Conj(MT)*Lambdax +
      0.25*vT*AbsSqr(Lambdax)*Sqr(vd) + 0.25*Conj(Mu)*Lambdax*Sqr(vd) + 0.25*Conj
      (Lambdax)*Mu*Sqr(vd) + 0.25*vT*AbsSqr(Lambdax)*Sqr(vu) + 0.25*Conj(Mu)*
      Lambdax*Sqr(vu) + 0.25*Conj(Lambdax)*Mu*Sqr(vu) - 0.25*vd*vu*TLambdax;

   return result;
}



std::complex<double> CLASSNAME::CpUSdconjUSdVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_0;
   std::complex<double> tmp_1;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_0 += tmp_1;
   result += (0.13333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_0;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += 0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
         ThetaW()));
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

std::complex<double> CLASSNAME::CpUSdconjUSdAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2;
   std::complex<double> tmp_3;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_2 += tmp_3;
   result += (0.1*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_2;
   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_6;
      std::complex<double> tmp_7;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_8;
         std::complex<double> tmp_9;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_9 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_8 += tmp_9;
         tmp_7 += (KroneckerDelta(gO2,3 + j2)) * tmp_8;
      }
      tmp_6 += tmp_7;
      tmp_5 += (KroneckerDelta(gO1,3 + j3)) * tmp_6;
   }
   tmp_4 += tmp_5;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_4;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_10;
      std::complex<double> tmp_11;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_11 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_10 += tmp_11;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_10;
   }
   std::complex<double> tmp_12;
   std::complex<double> tmp_13;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_13 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_12 += tmp_13;
   result += (-0.1*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_12;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_14;
      std::complex<double> tmp_15;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_15 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_14 += tmp_15;
      result += (-0.35355339059327373*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,1)) *
         tmp_14;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_16;
      std::complex<double> tmp_17;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_17 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_16 += tmp_17;
      result += (-0.35355339059327373*Lambdax*ZA(gI1,2)*ZA(gI2,1)) * tmp_16;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_18;
      std::complex<double> tmp_19;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_19 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_18 += tmp_19;
      result += (-0.35355339059327373*Conj(Lambdax)*ZA(gI1,1)*ZA(gI2,2)) *
         tmp_18;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_20;
      std::complex<double> tmp_21;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_21 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_20 += tmp_21;
      result += (-0.35355339059327373*Lambdax*ZA(gI1,1)*ZA(gI2,2)) * tmp_20;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_22;
   std::complex<double> tmp_23;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_23 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_22 += tmp_23;
   result += (0.1*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_22;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1)
         ;
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2)
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_24;
   std::complex<double> tmp_25;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_25 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_24 += tmp_25;
   result += (0.1*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_24;
   std::complex<double> tmp_26;
   std::complex<double> tmp_27;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_28;
      std::complex<double> tmp_29;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_30;
         std::complex<double> tmp_31;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_31 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_30 += tmp_31;
         tmp_29 += (KroneckerDelta(gO2,3 + j2)) * tmp_30;
      }
      tmp_28 += tmp_29;
      tmp_27 += (KroneckerDelta(gO1,3 + j3)) * tmp_28;
   }
   tmp_26 += tmp_27;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_26;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_32;
      std::complex<double> tmp_33;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_33 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_32 += tmp_33;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_32;
   }
   std::complex<double> tmp_34;
   std::complex<double> tmp_35;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_35 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_34 += tmp_35;
   result += (-0.1*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_34;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_36;
      std::complex<double> tmp_37;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_37 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_36 += tmp_37;
      result += (0.35355339059327373*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,1)) *
         tmp_36;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_38;
      std::complex<double> tmp_39;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_39 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_38 += tmp_39;
      result += (0.35355339059327373*Lambdax*ZH(gI1,2)*ZH(gI2,1)) * tmp_38;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_40;
      std::complex<double> tmp_41;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_41 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_40 += tmp_41;
      result += (0.35355339059327373*Conj(Lambdax)*ZH(gI1,1)*ZH(gI2,2)) *
         tmp_40;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_42;
      std::complex<double> tmp_43;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_43 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_42 += tmp_43;
      result += (0.35355339059327373*Lambdax*ZH(gI1,1)*ZH(gI2,2)) * tmp_42;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_44;
      std::complex<double> tmp_45;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_45 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_44 += tmp_45;
      result += (UP(gI2,1)) * tmp_44;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_46;
   std::complex<double> tmp_47;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_48;
      std::complex<double> tmp_49;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_49 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_48 += tmp_49;
      tmp_47 += (Conj(ZUL(gI1,j2))) * tmp_48;
   }
   tmp_46 += tmp_47;
   result += (Conj(UM(gI2,1))) * tmp_46;
   if (gO1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*Conj(ZUL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_50;
   std::complex<double> tmp_51;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_51 += KroneckerDelta(gO2,3 + j1)*ZDR(gI1,j1);
   }
   tmp_50 += tmp_51;
   result += (-0.3651483716701107*g1*ZN(gI2,0)) * tmp_50;
   if (gO2 < 3) {
      std::complex<double> tmp_52;
      std::complex<double> tmp_53;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_53 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_52 += tmp_53;
      result += (-ZN(gI2,2)) * tmp_52;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_54;
   std::complex<double> tmp_55;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_56;
      std::complex<double> tmp_57;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_57 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_56 += tmp_57;
      tmp_55 += (Conj(ZDL(gI1,j2))) * tmp_56;
   }
   tmp_54 += tmp_55;
   result += (-Conj(ZN(gI2,2))) * tmp_54;
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZDL(gI1,gO1))*Conj(ZN(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_59 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_58 += tmp_59;
   result += (0.1*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_58;
   std::complex<double> tmp_60;
   std::complex<double> tmp_61;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_62;
      std::complex<double> tmp_63;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_64;
         std::complex<double> tmp_65;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_65 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_64 += tmp_65;
         tmp_63 += (KroneckerDelta(gO2,3 + j2)) * tmp_64;
      }
      tmp_62 += tmp_63;
      tmp_61 += (KroneckerDelta(gO1,3 + j3)) * tmp_62;
   }
   tmp_60 += tmp_61;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_60;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_66;
   std::complex<double> tmp_67;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_67 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_66 += tmp_67;
   result += (-0.1*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_66;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_68;
      std::complex<double> tmp_69;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_69 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_68 += tmp_69;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_68;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_70;
      std::complex<double> tmp_71;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_71 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_70 += tmp_71;
      result += (-(Conj(Lambdax)*ZP(gI1,2)*ZP(gI2,1))) * tmp_70;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_72;
      std::complex<double> tmp_73;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_73 += Conj(Yd(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_72 += tmp_73;
      result += (-(Lambdax*ZP(gI1,1)*ZP(gI2,2))) * tmp_72;
   }
   if (gO1 < 3) {
      result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_74;
   std::complex<double> tmp_76;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_76 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_74 += tmp_76;
   std::complex<double> tmp_75;
   std::complex<double> tmp_77;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_77 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_75 += tmp_77;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_74 * tmp_75;
   std::complex<double> tmp_78;
   std::complex<double> tmp_80;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_80 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_78 += tmp_80;
   std::complex<double> tmp_79;
   std::complex<double> tmp_81;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_81 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_79 += tmp_81;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_78 * tmp_79;
   std::complex<double> tmp_82;
   std::complex<double> tmp_84;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_84 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_82 += tmp_84;
   std::complex<double> tmp_83;
   std::complex<double> tmp_85;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_85 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_83 += tmp_85;
   result += (-0.05*Sqr(g1)) * tmp_82 * tmp_83;
   std::complex<double> tmp_86;
   std::complex<double> tmp_88;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_88 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_86 += tmp_88;
   std::complex<double> tmp_87;
   std::complex<double> tmp_89;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_89 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_87 += tmp_89;
   result += (-0.1*Sqr(g1)) * tmp_86 * tmp_87;
   std::complex<double> tmp_90;
   std::complex<double> tmp_92;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_92 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_90 += tmp_92;
   std::complex<double> tmp_91;
   std::complex<double> tmp_93;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_93 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_91 += tmp_93;
   result += (-0.05*Sqr(g1)) * tmp_90 * tmp_91;
   std::complex<double> tmp_94;
   std::complex<double> tmp_96;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_96 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_94 += tmp_96;
   std::complex<double> tmp_95;
   std::complex<double> tmp_97;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_97 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_95 += tmp_97;
   result += (-0.1*Sqr(g1)) * tmp_94 * tmp_95;
   std::complex<double> tmp_98;
   std::complex<double> tmp_100;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_100 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_98 += tmp_100;
   std::complex<double> tmp_99;
   std::complex<double> tmp_101;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_101 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_99 += tmp_101;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_98 * tmp_99;
   std::complex<double> tmp_102;
   std::complex<double> tmp_104;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_104 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_102 += tmp_104;
   std::complex<double> tmp_103;
   std::complex<double> tmp_105;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_105 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_103 += tmp_105;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_102 * tmp_103;
   std::complex<double> tmp_106;
   std::complex<double> tmp_108;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_109;
      std::complex<double> tmp_110;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_110 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_109 += tmp_110;
      tmp_108 += (Conj(ZD(gI2,j2))) * tmp_109;
   }
   tmp_106 += tmp_108;
   std::complex<double> tmp_107;
   std::complex<double> tmp_111;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_113;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_113 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_112 += tmp_113;
      tmp_111 += (ZD(gI1,j4)) * tmp_112;
   }
   tmp_107 += tmp_111;
   result += (-1) * tmp_106 * tmp_107;
   if (gO1 < 3) {
      std::complex<double> tmp_114;
      std::complex<double> tmp_115;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_115 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_114 += tmp_115;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_114;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_116;
      std::complex<double> tmp_117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_117 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_116 += tmp_117;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_116;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_118;
      std::complex<double> tmp_119;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_119 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_118 += tmp_119;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_118;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_120;
      std::complex<double> tmp_121;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_121 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_120 += tmp_121;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_120;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_122;
      std::complex<double> tmp_123;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_123 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_122 += tmp_123;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_122;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_124;
      std::complex<double> tmp_125;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_125 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_124 += tmp_125;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_124;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_126;
      std::complex<double> tmp_128;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_128 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_126 += tmp_128;
      std::complex<double> tmp_127;
      std::complex<double> tmp_129;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_130;
         std::complex<double> tmp_131;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_131 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_130 += tmp_131;
         tmp_129 += (ZD(gI1,j4)) * tmp_130;
      }
      tmp_127 += tmp_129;
      result += (-3) * tmp_126 * tmp_127;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_132;
      std::complex<double> tmp_133;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_133 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_132 += tmp_133;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_132;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_135;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_135 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_134 += tmp_135;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_134;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_136;
      std::complex<double> tmp_137;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_137 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_136 += tmp_137;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_136;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_138;
      std::complex<double> tmp_139;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_139 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_138 += tmp_139;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_138;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_140;
      std::complex<double> tmp_142;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_143;
         std::complex<double> tmp_144;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_144 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_143 += tmp_144;
         tmp_142 += (Conj(ZD(gI2,j2))) * tmp_143;
      }
      tmp_140 += tmp_142;
      std::complex<double> tmp_141;
      std::complex<double> tmp_145;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_145 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_141 += tmp_145;
      result += (-3) * tmp_140 * tmp_141;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_146;
      std::complex<double> tmp_148;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_148 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_146 += tmp_148;
      std::complex<double> tmp_147;
      std::complex<double> tmp_149;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_149 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_147 += tmp_149;
      result += (-1) * tmp_146 * tmp_147;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_150;
      std::complex<double> tmp_151;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_151 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_150 += tmp_151;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_150;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_152;
      std::complex<double> tmp_153;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_153 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_152 += tmp_153;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_152;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_154;
      std::complex<double> tmp_155;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_155 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_154 += tmp_155;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_154;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_156;
      std::complex<double> tmp_157;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_157 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_156 += tmp_157;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_156;
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

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_158;
   std::complex<double> tmp_160;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_160 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_158 += tmp_160;
   std::complex<double> tmp_159;
   std::complex<double> tmp_161;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_161 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_159 += tmp_161;
   result += (0.05*Sqr(g1)) * tmp_158 * tmp_159;
   std::complex<double> tmp_162;
   std::complex<double> tmp_164;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_164 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_162 += tmp_164;
   std::complex<double> tmp_163;
   std::complex<double> tmp_165;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_165 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_163 += tmp_165;
   result += (-0.1*Sqr(g1)) * tmp_162 * tmp_163;
   std::complex<double> tmp_166;
   std::complex<double> tmp_168;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_168 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_166 += tmp_168;
   std::complex<double> tmp_167;
   std::complex<double> tmp_169;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_169 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_167 += tmp_169;
   result += (0.05*Sqr(g1)) * tmp_166 * tmp_167;
   std::complex<double> tmp_170;
   std::complex<double> tmp_172;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_172 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_170 += tmp_172;
   std::complex<double> tmp_171;
   std::complex<double> tmp_173;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_173 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_171 += tmp_173;
   result += (-0.1*Sqr(g1)) * tmp_170 * tmp_171;
   if (gO1 < 3) {
      std::complex<double> tmp_174;
      std::complex<double> tmp_175;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_175 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_174 += tmp_175;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_174;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_176;
      std::complex<double> tmp_177;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_177 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_176 += tmp_177;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_176;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_178;
      std::complex<double> tmp_179;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_179 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_178 += tmp_179;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_178;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_181;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_181 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_180 += tmp_181;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_180;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_182;
      std::complex<double> tmp_183;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_183 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_182 += tmp_183;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_182;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_184;
      std::complex<double> tmp_185;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_185 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_184 += tmp_185;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_184;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_186;
      std::complex<double> tmp_188;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_188 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_186 += tmp_188;
      std::complex<double> tmp_187;
      std::complex<double> tmp_189;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_190;
         std::complex<double> tmp_191;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_191 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_190 += tmp_191;
         tmp_189 += (ZE(gI1,j4)) * tmp_190;
      }
      tmp_187 += tmp_189;
      result += (-1) * tmp_186 * tmp_187;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_192;
      std::complex<double> tmp_194;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_195;
         std::complex<double> tmp_196;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_196 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_195 += tmp_196;
         tmp_194 += (Conj(ZE(gI2,j2))) * tmp_195;
      }
      tmp_192 += tmp_194;
      std::complex<double> tmp_193;
      std::complex<double> tmp_197;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_197 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_193 += tmp_197;
      result += (-1) * tmp_192 * tmp_193;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_198;
   std::complex<double> tmp_200;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_200 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_198 += tmp_200;
   std::complex<double> tmp_199;
   std::complex<double> tmp_201;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_201 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_199 += tmp_201;
   result += (-0.05*Sqr(g1)) * tmp_198 * tmp_199;
   std::complex<double> tmp_202;
   std::complex<double> tmp_204;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_204 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_202 += tmp_204;
   std::complex<double> tmp_203;
   std::complex<double> tmp_205;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_205 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_203 += tmp_205;
   result += (0.2*Sqr(g1)) * tmp_202 * tmp_203;
   std::complex<double> tmp_206;
   std::complex<double> tmp_208;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_208 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_206 += tmp_208;
   std::complex<double> tmp_207;
   std::complex<double> tmp_209;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_209 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_207 += tmp_209;
   result += (-0.05*Sqr(g1)) * tmp_206 * tmp_207;
   std::complex<double> tmp_210;
   std::complex<double> tmp_212;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_212 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_210 += tmp_212;
   std::complex<double> tmp_211;
   std::complex<double> tmp_213;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_213 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_211 += tmp_213;
   result += (0.2*Sqr(g1)) * tmp_210 * tmp_211;
   std::complex<double> tmp_214;
   std::complex<double> tmp_216;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_217;
      std::complex<double> tmp_218;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_218 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_217 += tmp_218;
      tmp_216 += (Conj(ZU(gI2,j2))) * tmp_217;
   }
   tmp_214 += tmp_216;
   std::complex<double> tmp_215;
   std::complex<double> tmp_219;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_220;
      std::complex<double> tmp_221;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_221 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_220 += tmp_221;
      tmp_219 += (ZU(gI1,j4)) * tmp_220;
   }
   tmp_215 += tmp_219;
   result += (-1) * tmp_214 * tmp_215;
   if (gO1 < 3) {
      std::complex<double> tmp_222;
      std::complex<double> tmp_223;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_223 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_222 += tmp_223;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_222;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_224;
      std::complex<double> tmp_225;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_225 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_224 += tmp_225;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_224;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_226;
      std::complex<double> tmp_227;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_227 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_226 += tmp_227;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_226;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_228;
      std::complex<double> tmp_229;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_229 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_228 += tmp_229;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_228;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_230;
      std::complex<double> tmp_231;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_231 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_230 += tmp_231;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_230;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_232;
      std::complex<double> tmp_233;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_233 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_232 += tmp_233;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_232;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_234;
      std::complex<double> tmp_236;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_236 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_234 += tmp_236;
      std::complex<double> tmp_235;
      std::complex<double> tmp_237;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_237 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_235 += tmp_237;
      result += (-1) * tmp_234 * tmp_235;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_238;
   std::complex<double> tmp_239;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_240;
      std::complex<double> tmp_241;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_241 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_240 += tmp_241;
      tmp_239 += (Conj(ZD(gI1,j2))) * tmp_240;
   }
   tmp_238 += tmp_239;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) * tmp_238;
   if (gO2 < 3) {
      std::complex<double> tmp_242;
      std::complex<double> tmp_243;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_243 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_242 += tmp_243;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) *
         tmp_242;
   }
   std::complex<double> tmp_244;
   std::complex<double> tmp_245;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_246;
      std::complex<double> tmp_247;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_247 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_246 += tmp_247;
      tmp_245 += (Conj(ZD(gI1,j2))) * tmp_246;
   }
   tmp_244 += tmp_245;
   result += (std::complex<double>(0,-0.35355339059327373)*vT*Conj(Lambdax)*ZA(
      gI2,1)) * tmp_244;
   std::complex<double> tmp_248;
   std::complex<double> tmp_249;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_250;
      std::complex<double> tmp_251;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_251 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_250 += tmp_251;
      tmp_249 += (Conj(ZD(gI1,j2))) * tmp_250;
   }
   tmp_248 += tmp_249;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1)) *
      tmp_248;
   if (gO2 < 3) {
      std::complex<double> tmp_252;
      std::complex<double> tmp_253;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_253 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_252 += tmp_253;
      result += (std::complex<double>(0,0.35355339059327373)*vT*Lambdax*ZA(
         gI2,1)) * tmp_252;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_254;
      std::complex<double> tmp_255;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_255 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_254 += tmp_255;
      result += (std::complex<double>(0,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_254;
   }
   std::complex<double> tmp_256;
   std::complex<double> tmp_257;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_258;
      std::complex<double> tmp_259;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_259 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_258 += tmp_259;
      tmp_257 += (Conj(ZD(gI1,j2))) * tmp_258;
   }
   tmp_256 += tmp_257;
   result += (std::complex<double>(0,-0.35355339059327373)*vu*Conj(Lambdax)*ZA(
      gI2,2)) * tmp_256;
   if (gO2 < 3) {
      std::complex<double> tmp_260;
      std::complex<double> tmp_261;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_261 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_260 += tmp_261;
      result += (std::complex<double>(0,0.35355339059327373)*vu*Lambdax*ZA(
         gI2,2)) * tmp_260;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_262;
   std::complex<double> tmp_263;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_263 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_262 += tmp_263;
   result += (0.1*vd*Sqr(g1)*ZH(gI2,0)) * tmp_262;
   std::complex<double> tmp_264;
   std::complex<double> tmp_265;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_266;
      std::complex<double> tmp_267;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_267 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_266 += tmp_267;
      tmp_265 += (Conj(ZD(gI1,j2))) * tmp_266;
   }
   tmp_264 += tmp_265;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_264;
   std::complex<double> tmp_268;
   std::complex<double> tmp_269;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_270;
      std::complex<double> tmp_271;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_272;
         std::complex<double> tmp_273;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_273 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_272 += tmp_273;
         tmp_271 += (KroneckerDelta(gO2,3 + j2)) * tmp_272;
      }
      tmp_270 += tmp_271;
      tmp_269 += (Conj(ZD(gI1,3 + j3))) * tmp_270;
   }
   tmp_268 += tmp_269;
   result += (-(vd*ZH(gI2,0))) * tmp_268;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_274;
      std::complex<double> tmp_275;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_275 += Conj(ZD(gI1,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_274 += tmp_275;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_274;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_276;
      std::complex<double> tmp_277;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_278;
         std::complex<double> tmp_279;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_279 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_278 += tmp_279;
         tmp_277 += (Conj(ZD(gI1,j2))) * tmp_278;
      }
      tmp_276 += tmp_277;
      result += (-(vd*ZH(gI2,0))) * tmp_276;
   }
   std::complex<double> tmp_280;
   std::complex<double> tmp_281;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_281 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_280 += tmp_281;
   result += (-0.1*vu*Sqr(g1)*ZH(gI2,1)) * tmp_280;
   std::complex<double> tmp_282;
   std::complex<double> tmp_283;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_284;
      std::complex<double> tmp_285;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_285 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_284 += tmp_285;
      tmp_283 += (Conj(ZD(gI1,j2))) * tmp_284;
   }
   tmp_282 += tmp_283;
   result += (0.35355339059327373*vT*Conj(Lambdax)*ZH(gI2,1)) * tmp_282;
   std::complex<double> tmp_286;
   std::complex<double> tmp_287;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_288;
      std::complex<double> tmp_289;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_289 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_288 += tmp_289;
      tmp_287 += (Conj(ZD(gI1,j2))) * tmp_288;
   }
   tmp_286 += tmp_287;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_286;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_290;
      std::complex<double> tmp_291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_291 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_290 += tmp_291;
      result += (0.35355339059327373*vT*Lambdax*ZH(gI2,1)) * tmp_290;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_292;
      std::complex<double> tmp_293;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_293 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_292 += tmp_293;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_292;
   }
   std::complex<double> tmp_294;
   std::complex<double> tmp_295;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_296;
      std::complex<double> tmp_297;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_297 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_296 += tmp_297;
      tmp_295 += (Conj(ZD(gI1,j2))) * tmp_296;
   }
   tmp_294 += tmp_295;
   result += (0.35355339059327373*vu*Conj(Lambdax)*ZH(gI2,2)) * tmp_294;
   if (gO2 < 3) {
      std::complex<double> tmp_298;
      std::complex<double> tmp_299;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_299 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_298 += tmp_299;
      result += (0.35355339059327373*vu*Lambdax*ZH(gI2,2)) * tmp_298;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_300;
   std::complex<double> tmp_301;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_302;
      std::complex<double> tmp_303;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_303 += KroneckerDelta(gO2,3 + j1)*TYd(j1,j2);
      }
      tmp_302 += tmp_303;
      tmp_301 += (Conj(ZU(gI1,j2))) * tmp_302;
   }
   tmp_300 += tmp_301;
   result += (ZP(gI2,0)) * tmp_300;
   std::complex<double> tmp_304;
   std::complex<double> tmp_305;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_306;
      std::complex<double> tmp_307;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_308;
         std::complex<double> tmp_309;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_309 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_308 += tmp_309;
         tmp_307 += (KroneckerDelta(gO2,3 + j2)) * tmp_308;
      }
      tmp_306 += tmp_307;
      tmp_305 += (Conj(ZU(gI1,3 + j3))) * tmp_306;
   }
   tmp_304 += tmp_305;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_304;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_310;
      std::complex<double> tmp_311;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_311 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_310 += tmp_311;
      result += (-0.5*vT*Lambdax*ZP(gI2,0)) * tmp_310;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_312;
      std::complex<double> tmp_313;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_313 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_312 += tmp_313;
      result += (Mu*ZP(gI2,0)) * tmp_312;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_314;
      std::complex<double> tmp_315;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_316;
         std::complex<double> tmp_317;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_317 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_316 += tmp_317;
         tmp_315 += (Conj(ZU(gI1,j2))) * tmp_316;
      }
      tmp_314 += tmp_315;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_314;
   }
   std::complex<double> tmp_318;
   std::complex<double> tmp_319;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_320;
      std::complex<double> tmp_321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_321 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_320 += tmp_321;
      tmp_319 += (Conj(ZU(gI1,j2))) * tmp_320;
   }
   tmp_318 += tmp_319;
   result += (-0.5*vT*Conj(Lambdax)*ZP(gI2,1)) * tmp_318;
   std::complex<double> tmp_322;
   std::complex<double> tmp_323;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_324;
      std::complex<double> tmp_325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_325 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_324 += tmp_325;
      tmp_323 += (Conj(ZU(gI1,j2))) * tmp_324;
   }
   tmp_322 += tmp_323;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_322;
   std::complex<double> tmp_326;
   std::complex<double> tmp_327;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_328;
      std::complex<double> tmp_329;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_330;
         std::complex<double> tmp_331;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_331 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_330 += tmp_331;
         tmp_329 += (KroneckerDelta(gO2,3 + j2)) * tmp_330;
      }
      tmp_328 += tmp_329;
      tmp_327 += (Conj(ZU(gI1,3 + j3))) * tmp_328;
   }
   tmp_326 += tmp_327;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_326;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_332;
      std::complex<double> tmp_333;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_333 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_332 += tmp_333;
      result += (ZP(gI2,1)) * tmp_332;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_334;
      std::complex<double> tmp_335;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_336;
         std::complex<double> tmp_337;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_337 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_336 += tmp_337;
         tmp_335 += (Conj(ZU(gI1,j2))) * tmp_336;
      }
      tmp_334 += tmp_335;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_334;
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,2);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_338;
      std::complex<double> tmp_339;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_339 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_338 += tmp_339;
      result += (0.7071067811865475*vd*Lambdax*ZP(gI2,2)) * tmp_338;
   }
   std::complex<double> tmp_340;
   std::complex<double> tmp_341;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_342;
      std::complex<double> tmp_343;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_343 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_342 += tmp_343;
      tmp_341 += (Conj(ZU(gI1,j2))) * tmp_342;
   }
   tmp_340 += tmp_341;
   result += (-0.7071067811865475*vu*Conj(Lambdax)*ZP(gI2,3)) * tmp_340;
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_344;
   std::complex<double> tmp_345;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_345 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_344 += tmp_345;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_344;

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

   std::complex<double> tmp_346;
   std::complex<double> tmp_347;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_347 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_346 += tmp_347;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_346;
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
   std::complex<double> result;

   std::complex<double> tmp_348;
   std::complex<double> tmp_349;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_349 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_348 += tmp_349;
   result += (0.2581988897471611*g1*Sin(ThetaW())) * tmp_348;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZD(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZD(gI2,gO2))*Sin(ThetaW());
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
   std::complex<double> result;

   result = 0.1*KroneckerDelta(gO1,gO2)*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpUSvconjUSvconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,0)*ZA
      (gI2,0) - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result += -0.15*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1);
   result += -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2);
   if (gI1 < 3 && gI2 < 3) {
      result += -0.15*Conj(ZV(gI2,gO2))*Sqr(g1)*ZV(gI1,gO1);
   }
   if (gI1 < 3 && gI2 < 3) {
      result += -0.25*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,0)*ZH
      (gI2,0) - ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += -0.15*vd*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gI1 < 3) {
      result += -0.25*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gI1 < 3) {
      result += 0.15*vu*Conj(ZV(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gI1 < 3) {
      result += 0.25*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvbarChaFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_350;
      std::complex<double> tmp_351;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_351 += Conj(Ye(j1,gO2))*ZER(gI2,j1);
      }
      tmp_350 += tmp_351;
      result += (UM(gI1,1)) * tmp_350;
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

double CLASSNAME::CpconjUSvFvChiPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFvChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZN(gI2,0))*KroneckerDelta(gI1,gO1
         );
   }
   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,1))*KroneckerDelta(gI1,
         gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_352;
      std::complex<double> tmp_353;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_353 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_352 += tmp_353;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_352;
   }
   result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_354;
      std::complex<double> tmp_355;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_355 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_354 += tmp_355;
      result += (ZP(gI1,0)) * tmp_354;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_356;
      std::complex<double> tmp_357;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_358;
         std::complex<double> tmp_359;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_359 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_358 += tmp_359;
         tmp_357 += (Conj(ZE(gI2,j2))) * tmp_358;
      }
      tmp_356 += tmp_357;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_356;
   }
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_360;
      std::complex<double> tmp_361;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_361 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_360 += tmp_361;
      result += (-0.5*vT*Lambdax*ZP(gI1,1)) * tmp_360;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_362;
      std::complex<double> tmp_363;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_363 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_362 += tmp_363;
      result += (Mu*ZP(gI1,1)) * tmp_362;
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,2);
   }
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,3);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_364;
      std::complex<double> tmp_365;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_365 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_364 += tmp_365;
      result += (-0.7071067811865475*vu*Lambdax*ZP(gI1,3)) * tmp_364;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_366;
   std::complex<double> tmp_367;
   std::complex<double> tmp_368;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_368 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_367 += tmp_368;
   tmp_366 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_367;
   std::complex<double> tmp_369;
   std::complex<double> tmp_370;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_370 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_369 += tmp_370;
   tmp_366 += (std::complex<double>(0,0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_369;
   std::complex<double> tmp_371;
   std::complex<double> tmp_372;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_372 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_371 += tmp_372;
   tmp_366 += (std::complex<double>(0,0.1)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_371;
   result += (std::complex<double>(0,-1)) * tmp_366;

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_373;
   std::complex<double> tmp_374;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_374 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_373 += tmp_374;
   result += (-0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_373;
   std::complex<double> tmp_375;
   std::complex<double> tmp_376;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_376 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_375 += tmp_376;
   result += (0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_375;
   std::complex<double> tmp_377;
   std::complex<double> tmp_378;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_378 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_377 += tmp_378;
   result += (0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_377;
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_379;
      std::complex<double> tmp_381;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_381 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_379 += tmp_381;
      std::complex<double> tmp_380;
      std::complex<double> tmp_382;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_382 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_380 += tmp_382;
      result += (-1) * tmp_379 * tmp_380;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_383;
   std::complex<double> tmp_384;
   std::complex<double> tmp_385;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_385 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_384 += tmp_385;
   tmp_383 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_384;
   std::complex<double> tmp_386;
   std::complex<double> tmp_387;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_387 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_386 += tmp_387;
   tmp_383 += (std::complex<double>(0,-0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_386;
   std::complex<double> tmp_388;
   std::complex<double> tmp_389;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_389 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_388 += tmp_389;
   tmp_383 += (std::complex<double>(0,-0.2)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_388;
   result += (std::complex<double>(0,-1)) * tmp_383;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvVZSv(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZV(gI2,gO2))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZV(gI2,gO2))*Sin(ThetaW());
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
   std::complex<double> result;

   std::complex<double> tmp_390;
   std::complex<double> tmp_391;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_391 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_390 += tmp_391;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_390;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += -0.2581988897471611*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.03333333333333333*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(
         ThetaW()));
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

std::complex<double> CLASSNAME::CpUSuconjUSuAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_392;
   std::complex<double> tmp_393;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_393 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_392 += tmp_393;
   result += (-0.2*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_392;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_394;
      std::complex<double> tmp_395;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_395 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_394 += tmp_395;
      result += (-0.35355339059327373*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,0)) *
         tmp_394;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_396;
      std::complex<double> tmp_397;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_397 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_396 += tmp_397;
      result += (-0.35355339059327373*Lambdax*ZA(gI1,2)*ZA(gI2,0)) * tmp_396
         ;
   }
   std::complex<double> tmp_398;
   std::complex<double> tmp_399;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_399 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_398 += tmp_399;
   result += (0.2*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_398;
   std::complex<double> tmp_400;
   std::complex<double> tmp_401;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_402;
      std::complex<double> tmp_403;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_404;
         std::complex<double> tmp_405;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_405 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_404 += tmp_405;
         tmp_403 += (KroneckerDelta(gO2,3 + j2)) * tmp_404;
      }
      tmp_402 += tmp_403;
      tmp_401 += (KroneckerDelta(gO1,3 + j3)) * tmp_402;
   }
   tmp_400 += tmp_401;
   result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_400;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_406;
      std::complex<double> tmp_407;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_407 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_406 += tmp_407;
      result += (-(ZA(gI1,1)*ZA(gI2,1))) * tmp_406;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_408;
      std::complex<double> tmp_409;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_409 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_408 += tmp_409;
      result += (-0.35355339059327373*Conj(Lambdax)*ZA(gI1,0)*ZA(gI2,2)) *
         tmp_408;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_410;
      std::complex<double> tmp_411;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_411 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_410 += tmp_411;
      result += (-0.35355339059327373*Lambdax*ZA(gI1,0)*ZA(gI2,2)) * tmp_410
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_412;
   std::complex<double> tmp_413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_413 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_412 += tmp_413;
   result += (-0.2*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_412;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1)
         ;
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_415 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_414 += tmp_415;
   result += (-0.2*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_414;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_416;
      std::complex<double> tmp_417;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_417 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_416 += tmp_417;
      result += (0.35355339059327373*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,0)) *
         tmp_416;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_418;
      std::complex<double> tmp_419;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_419 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_418 += tmp_419;
      result += (0.35355339059327373*Lambdax*ZH(gI1,2)*ZH(gI2,0)) * tmp_418;
   }
   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_421 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_420 += tmp_421;
   result += (0.2*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_420;
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_424;
      std::complex<double> tmp_425;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_426;
         std::complex<double> tmp_427;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_427 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_426 += tmp_427;
         tmp_425 += (KroneckerDelta(gO2,3 + j2)) * tmp_426;
      }
      tmp_424 += tmp_425;
      tmp_423 += (KroneckerDelta(gO1,3 + j3)) * tmp_424;
   }
   tmp_422 += tmp_423;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_422;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_428;
      std::complex<double> tmp_429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_429 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_428 += tmp_429;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_428;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_430;
      std::complex<double> tmp_431;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_431 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_430 += tmp_431;
      result += (0.35355339059327373*Conj(Lambdax)*ZH(gI1,0)*ZH(gI2,2)) *
         tmp_430;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_432;
      std::complex<double> tmp_433;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_433 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_432 += tmp_433;
      result += (0.35355339059327373*Lambdax*ZH(gI1,0)*ZH(gI2,2)) * tmp_432;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_434;
      std::complex<double> tmp_435;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_435 += Conj(Yd(j1,gO2))*ZDR(gI2,j1);
      }
      tmp_434 += tmp_435;
      result += (UM(gI1,1)) * tmp_434;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChaFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_436;
   std::complex<double> tmp_437;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_438;
      std::complex<double> tmp_439;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_439 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_438 += tmp_439;
      tmp_437 += (Conj(ZDL(gI2,j2))) * tmp_438;
   }
   tmp_436 += tmp_437;
   result += (Conj(UP(gI1,1))) * tmp_436;
   if (gO1 < 3) {
      result += -(g2*Conj(UP(gI1,0))*Conj(ZDL(gI2,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_440;
   std::complex<double> tmp_441;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_441 += KroneckerDelta(gO2,3 + j1)*ZUR(gI1,j1);
   }
   tmp_440 += tmp_441;
   result += (0.7302967433402214*g1*ZN(gI2,0)) * tmp_440;
   if (gO2 < 3) {
      std::complex<double> tmp_442;
      std::complex<double> tmp_443;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_443 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_442 += tmp_443;
      result += (-ZN(gI2,3)) * tmp_442;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_444;
   std::complex<double> tmp_445;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_446;
      std::complex<double> tmp_447;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_447 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_446 += tmp_447;
      tmp_445 += (Conj(ZUL(gI1,j2))) * tmp_446;
   }
   tmp_444 += tmp_445;
   result += (-Conj(ZN(gI2,3))) * tmp_444;
   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZN(gI2,0))*Conj(ZUL(gI1,gO1));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN(gI2,1))*Conj(ZUL(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_448;
   std::complex<double> tmp_449;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_449 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_448 += tmp_449;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_448;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_450;
      std::complex<double> tmp_451;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_451 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_450 += tmp_451;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_450;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_452;
      std::complex<double> tmp_453;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_453 += Conj(Yu(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_452 += tmp_453;
      result += (Lambdax*ZP(gI1,3)*ZP(gI2,0)) * tmp_452;
   }
   std::complex<double> tmp_454;
   std::complex<double> tmp_455;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_455 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_454 += tmp_455;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_454;
   std::complex<double> tmp_456;
   std::complex<double> tmp_457;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_458;
      std::complex<double> tmp_459;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_460;
         std::complex<double> tmp_461;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_461 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_460 += tmp_461;
         tmp_459 += (KroneckerDelta(gO2,3 + j2)) * tmp_460;
      }
      tmp_458 += tmp_459;
      tmp_457 += (KroneckerDelta(gO1,3 + j3)) * tmp_458;
   }
   tmp_456 += tmp_457;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_456;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_462;
      std::complex<double> tmp_463;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_463 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_462 += tmp_463;
      result += (Conj(Lambdax)*ZP(gI1,0)*ZP(gI2,3)) * tmp_462;
   }
   if (gO1 < 3) {
      result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_464;
   std::complex<double> tmp_465;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_466;
      std::complex<double> tmp_467;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_467 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_466 += tmp_467;
      tmp_465 += (Conj(ZD(gI2,j2))) * tmp_466;
   }
   tmp_464 += tmp_465;
   result += (-0.5*vT*Conj(Lambdax)*ZP(gI1,0)) * tmp_464;
   std::complex<double> tmp_468;
   std::complex<double> tmp_469;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_470;
      std::complex<double> tmp_471;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_471 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_470 += tmp_471;
      tmp_469 += (Conj(ZD(gI2,j2))) * tmp_470;
   }
   tmp_468 += tmp_469;
   result += (Conj(Mu)*ZP(gI1,0)) * tmp_468;
   std::complex<double> tmp_472;
   std::complex<double> tmp_473;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_474;
      std::complex<double> tmp_475;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_476;
         std::complex<double> tmp_477;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_477 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_476 += tmp_477;
         tmp_475 += (KroneckerDelta(gO2,3 + j2)) * tmp_476;
      }
      tmp_474 += tmp_475;
      tmp_473 += (Conj(ZD(gI2,3 + j3))) * tmp_474;
   }
   tmp_472 += tmp_473;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_472;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_478;
      std::complex<double> tmp_479;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_479 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,gO2));
      }
      tmp_478 += tmp_479;
      result += (ZP(gI1,0)) * tmp_478;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_480;
      std::complex<double> tmp_481;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_482;
         std::complex<double> tmp_483;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_483 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_482 += tmp_483;
         tmp_481 += (Conj(ZD(gI2,j2))) * tmp_482;
      }
      tmp_480 += tmp_481;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_480;
   }
   std::complex<double> tmp_484;
   std::complex<double> tmp_485;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_486;
      std::complex<double> tmp_487;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_487 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_486 += tmp_487;
      tmp_485 += (Conj(ZD(gI2,j2))) * tmp_486;
   }
   tmp_484 += tmp_485;
   result += (ZP(gI1,1)) * tmp_484;
   std::complex<double> tmp_488;
   std::complex<double> tmp_489;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_490;
      std::complex<double> tmp_491;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_492;
         std::complex<double> tmp_493;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_493 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_492 += tmp_493;
         tmp_491 += (KroneckerDelta(gO2,3 + j2)) * tmp_492;
      }
      tmp_490 += tmp_491;
      tmp_489 += (Conj(ZD(gI2,3 + j3))) * tmp_490;
   }
   tmp_488 += tmp_489;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_488;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_494;
      std::complex<double> tmp_495;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_495 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_494 += tmp_495;
      result += (-0.5*vT*Lambdax*ZP(gI1,1)) * tmp_494;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_496;
      std::complex<double> tmp_497;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_497 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_496 += tmp_497;
      result += (Mu*ZP(gI1,1)) * tmp_496;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_498;
      std::complex<double> tmp_499;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_500;
         std::complex<double> tmp_501;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_501 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_500 += tmp_501;
         tmp_499 += (Conj(ZD(gI2,j2))) * tmp_500;
      }
      tmp_498 += tmp_499;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_498;
   }
   std::complex<double> tmp_502;
   std::complex<double> tmp_503;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_504;
      std::complex<double> tmp_505;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_505 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_504 += tmp_505;
      tmp_503 += (Conj(ZD(gI2,j2))) * tmp_504;
   }
   tmp_502 += tmp_503;
   result += (0.7071067811865475*vd*Conj(Lambdax)*ZP(gI1,2)) * tmp_502;
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,2);
   }
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,3);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_506;
      std::complex<double> tmp_507;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_507 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_506 += tmp_507;
      result += (-0.7071067811865475*vu*Lambdax*ZP(gI1,3)) * tmp_506;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_508;
   std::complex<double> tmp_510;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_510 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_508 += tmp_510;
   std::complex<double> tmp_509;
   std::complex<double> tmp_511;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_511 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_509 += tmp_511;
   result += (0.1*Sqr(g1)) * tmp_508 * tmp_509;
   std::complex<double> tmp_512;
   std::complex<double> tmp_514;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_514 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_512 += tmp_514;
   std::complex<double> tmp_513;
   std::complex<double> tmp_515;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_515 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_513 += tmp_515;
   result += (0.2*Sqr(g1)) * tmp_512 * tmp_513;
   std::complex<double> tmp_516;
   std::complex<double> tmp_518;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_518 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_516 += tmp_518;
   std::complex<double> tmp_517;
   std::complex<double> tmp_519;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_519 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_517 += tmp_519;
   result += (0.1*Sqr(g1)) * tmp_516 * tmp_517;
   std::complex<double> tmp_520;
   std::complex<double> tmp_522;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_522 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_520 += tmp_522;
   std::complex<double> tmp_521;
   std::complex<double> tmp_523;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_523 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_521 += tmp_523;
   result += (0.2*Sqr(g1)) * tmp_520 * tmp_521;
   std::complex<double> tmp_524;
   std::complex<double> tmp_526;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_527;
      std::complex<double> tmp_528;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_528 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_527 += tmp_528;
      tmp_526 += (Conj(ZD(gI2,j2))) * tmp_527;
   }
   tmp_524 += tmp_526;
   std::complex<double> tmp_525;
   std::complex<double> tmp_529;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_530;
      std::complex<double> tmp_531;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_531 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_530 += tmp_531;
      tmp_529 += (ZD(gI1,j4)) * tmp_530;
   }
   tmp_525 += tmp_529;
   result += (-1) * tmp_524 * tmp_525;
   if (gO1 < 3) {
      std::complex<double> tmp_532;
      std::complex<double> tmp_533;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_533 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_532 += tmp_533;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_532;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_534;
      std::complex<double> tmp_535;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_535 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_534 += tmp_535;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_534;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_536;
      std::complex<double> tmp_537;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_537 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_536 += tmp_537;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_536;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_538;
      std::complex<double> tmp_539;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_539 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_538 += tmp_539;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_538;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_540;
      std::complex<double> tmp_541;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_541 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_540 += tmp_541;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_540;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_542;
      std::complex<double> tmp_543;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_543 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_542 += tmp_543;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_542;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_544;
      std::complex<double> tmp_546;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_546 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_544 += tmp_546;
      std::complex<double> tmp_545;
      std::complex<double> tmp_547;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_547 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_545 += tmp_547;
      result += (-1) * tmp_544 * tmp_545;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_548;
   std::complex<double> tmp_550;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_550 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_548 += tmp_550;
   std::complex<double> tmp_549;
   std::complex<double> tmp_551;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_551 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_549 += tmp_551;
   result += (-0.1*Sqr(g1)) * tmp_548 * tmp_549;
   std::complex<double> tmp_552;
   std::complex<double> tmp_554;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_554 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_552 += tmp_554;
   std::complex<double> tmp_553;
   std::complex<double> tmp_555;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_555 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_553 += tmp_555;
   result += (0.2*Sqr(g1)) * tmp_552 * tmp_553;
   std::complex<double> tmp_556;
   std::complex<double> tmp_558;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_558 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_556 += tmp_558;
   std::complex<double> tmp_557;
   std::complex<double> tmp_559;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_559 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_557 += tmp_559;
   result += (-0.1*Sqr(g1)) * tmp_556 * tmp_557;
   std::complex<double> tmp_560;
   std::complex<double> tmp_562;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_562 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_560 += tmp_562;
   std::complex<double> tmp_561;
   std::complex<double> tmp_563;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_563 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_561 += tmp_563;
   result += (0.2*Sqr(g1)) * tmp_560 * tmp_561;
   if (gO1 < 3) {
      std::complex<double> tmp_564;
      std::complex<double> tmp_565;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_565 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_564 += tmp_565;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_564;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_566;
      std::complex<double> tmp_567;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_567 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_566 += tmp_567;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_566;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_568;
      std::complex<double> tmp_569;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_569 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_568 += tmp_569;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_568;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_570;
      std::complex<double> tmp_571;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_571 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_570 += tmp_571;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_570;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_572;
      std::complex<double> tmp_573;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_573 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_572 += tmp_573;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_572;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_574;
      std::complex<double> tmp_575;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_575 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_574 += tmp_575;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_574;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_576;
   std::complex<double> tmp_578;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_578 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_576 += tmp_578;
   std::complex<double> tmp_577;
   std::complex<double> tmp_579;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_579 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_577 += tmp_579;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_576 * tmp_577;
   std::complex<double> tmp_580;
   std::complex<double> tmp_582;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_582 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_580 += tmp_582;
   std::complex<double> tmp_581;
   std::complex<double> tmp_583;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_583 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_581 += tmp_583;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_580 * tmp_581;
   std::complex<double> tmp_584;
   std::complex<double> tmp_586;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_586 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_584 += tmp_586;
   std::complex<double> tmp_585;
   std::complex<double> tmp_587;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_587 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_585 += tmp_587;
   result += (0.1*Sqr(g1)) * tmp_584 * tmp_585;
   std::complex<double> tmp_588;
   std::complex<double> tmp_590;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_590 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_588 += tmp_590;
   std::complex<double> tmp_589;
   std::complex<double> tmp_591;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_591 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_589 += tmp_591;
   result += (-0.4*Sqr(g1)) * tmp_588 * tmp_589;
   std::complex<double> tmp_592;
   std::complex<double> tmp_594;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_594 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_592 += tmp_594;
   std::complex<double> tmp_593;
   std::complex<double> tmp_595;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_595 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_593 += tmp_595;
   result += (0.1*Sqr(g1)) * tmp_592 * tmp_593;
   std::complex<double> tmp_596;
   std::complex<double> tmp_598;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_598 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_596 += tmp_598;
   std::complex<double> tmp_597;
   std::complex<double> tmp_599;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_599 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_597 += tmp_599;
   result += (-0.4*Sqr(g1)) * tmp_596 * tmp_597;
   std::complex<double> tmp_600;
   std::complex<double> tmp_602;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_602 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_600 += tmp_602;
   std::complex<double> tmp_601;
   std::complex<double> tmp_603;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_603 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_601 += tmp_603;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_600 * tmp_601;
   std::complex<double> tmp_604;
   std::complex<double> tmp_606;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_606 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_604 += tmp_606;
   std::complex<double> tmp_605;
   std::complex<double> tmp_607;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_607 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_605 += tmp_607;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_604 * tmp_605;
   std::complex<double> tmp_608;
   std::complex<double> tmp_610;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_611;
      std::complex<double> tmp_612;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_612 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_611 += tmp_612;
      tmp_610 += (Conj(ZU(gI2,j2))) * tmp_611;
   }
   tmp_608 += tmp_610;
   std::complex<double> tmp_609;
   std::complex<double> tmp_613;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_614;
      std::complex<double> tmp_615;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_615 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_614 += tmp_615;
      tmp_613 += (ZU(gI1,j4)) * tmp_614;
   }
   tmp_609 += tmp_613;
   result += (-1) * tmp_608 * tmp_609;
   if (gO1 < 3) {
      std::complex<double> tmp_616;
      std::complex<double> tmp_617;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_617 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_616 += tmp_617;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_616;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_618;
      std::complex<double> tmp_619;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_619 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_618 += tmp_619;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_618;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_620;
      std::complex<double> tmp_621;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_621 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_620 += tmp_621;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_620;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_622;
      std::complex<double> tmp_623;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_623 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_622 += tmp_623;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_622;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_624;
      std::complex<double> tmp_625;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_625 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_624 += tmp_625;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_624;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_627 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_626 += tmp_627;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_626;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_628;
      std::complex<double> tmp_630;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_630 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_628 += tmp_630;
      std::complex<double> tmp_629;
      std::complex<double> tmp_631;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_632;
         std::complex<double> tmp_633;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_633 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_632 += tmp_633;
         tmp_631 += (ZU(gI1,j4)) * tmp_632;
      }
      tmp_629 += tmp_631;
      result += (-3) * tmp_628 * tmp_629;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_634;
      std::complex<double> tmp_635;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_635 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_634 += tmp_635;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_634;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_636;
      std::complex<double> tmp_637;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_637 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_636 += tmp_637;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_636;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_638;
      std::complex<double> tmp_639;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_639 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_638 += tmp_639;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_638;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_640;
      std::complex<double> tmp_641;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_641 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_640 += tmp_641;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_640;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_642;
      std::complex<double> tmp_644;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_645;
         std::complex<double> tmp_646;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_646 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_645 += tmp_646;
         tmp_644 += (Conj(ZU(gI2,j2))) * tmp_645;
      }
      tmp_642 += tmp_644;
      std::complex<double> tmp_643;
      std::complex<double> tmp_647;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_647 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_643 += tmp_647;
      result += (-3) * tmp_642 * tmp_643;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_648;
      std::complex<double> tmp_650;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_650 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_648 += tmp_650;
      std::complex<double> tmp_649;
      std::complex<double> tmp_651;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_651 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_649 += tmp_651;
      result += (-1) * tmp_648 * tmp_649;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_652;
      std::complex<double> tmp_653;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_653 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_652 += tmp_653;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_652;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_654;
      std::complex<double> tmp_655;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_655 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_654 += tmp_655;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_654;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_656;
      std::complex<double> tmp_657;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_657 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_656 += tmp_657;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_656;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_658;
      std::complex<double> tmp_659;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_659 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_658 += tmp_659;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_658;
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

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_660;
   std::complex<double> tmp_661;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_662;
      std::complex<double> tmp_663;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_663 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_662 += tmp_663;
      tmp_661 += (Conj(ZU(gI1,j2))) * tmp_662;
   }
   tmp_660 += tmp_661;
   result += (std::complex<double>(0,-0.35355339059327373)*vT*Conj(Lambdax)*ZA(
      gI2,0)) * tmp_660;
   std::complex<double> tmp_664;
   std::complex<double> tmp_665;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_666;
      std::complex<double> tmp_667;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_667 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_666 += tmp_667;
      tmp_665 += (Conj(ZU(gI1,j2))) * tmp_666;
   }
   tmp_664 += tmp_665;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(Mu)*ZA(gI2,0)) *
      tmp_664;
   if (gO2 < 3) {
      std::complex<double> tmp_668;
      std::complex<double> tmp_669;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_669 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_668 += tmp_669;
      result += (std::complex<double>(0,0.35355339059327373)*vT*Lambdax*ZA(
         gI2,0)) * tmp_668;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_670;
      std::complex<double> tmp_671;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_671 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_670 += tmp_671;
      result += (std::complex<double>(0,0.7071067811865475)*Mu*ZA(gI2,0)) *
         tmp_670;
   }
   std::complex<double> tmp_672;
   std::complex<double> tmp_673;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_674;
      std::complex<double> tmp_675;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_675 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_674 += tmp_675;
      tmp_673 += (Conj(ZU(gI1,j2))) * tmp_674;
   }
   tmp_672 += tmp_673;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,1)) * tmp_672;
   if (gO2 < 3) {
      std::complex<double> tmp_676;
      std::complex<double> tmp_677;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_677 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_676 += tmp_677;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,1)) *
         tmp_676;
   }
   std::complex<double> tmp_678;
   std::complex<double> tmp_679;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_680;
      std::complex<double> tmp_681;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_681 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_680 += tmp_681;
      tmp_679 += (Conj(ZU(gI1,j2))) * tmp_680;
   }
   tmp_678 += tmp_679;
   result += (std::complex<double>(0,-0.35355339059327373)*vd*Conj(Lambdax)*ZA(
      gI2,2)) * tmp_678;
   if (gO2 < 3) {
      std::complex<double> tmp_682;
      std::complex<double> tmp_683;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_683 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_682 += tmp_683;
      result += (std::complex<double>(0,0.35355339059327373)*vd*Lambdax*ZA(
         gI2,2)) * tmp_682;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_684;
   std::complex<double> tmp_685;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_685 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_684 += tmp_685;
   result += (-0.2*vd*Sqr(g1)*ZH(gI2,0)) * tmp_684;
   std::complex<double> tmp_686;
   std::complex<double> tmp_687;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_688;
      std::complex<double> tmp_689;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_689 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_688 += tmp_689;
      tmp_687 += (Conj(ZU(gI1,j2))) * tmp_688;
   }
   tmp_686 += tmp_687;
   result += (0.35355339059327373*vT*Conj(Lambdax)*ZH(gI2,0)) * tmp_686;
   std::complex<double> tmp_690;
   std::complex<double> tmp_691;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_692;
      std::complex<double> tmp_693;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_693 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_692 += tmp_693;
      tmp_691 += (Conj(ZU(gI1,j2))) * tmp_692;
   }
   tmp_690 += tmp_691;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,0)) * tmp_690;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_694;
      std::complex<double> tmp_695;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_695 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_694 += tmp_695;
      result += (0.35355339059327373*vT*Lambdax*ZH(gI2,0)) * tmp_694;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_696;
      std::complex<double> tmp_697;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_697 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_696 += tmp_697;
      result += (0.7071067811865475*Mu*ZH(gI2,0)) * tmp_696;
   }
   std::complex<double> tmp_698;
   std::complex<double> tmp_699;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_699 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_698 += tmp_699;
   result += (0.2*vu*Sqr(g1)*ZH(gI2,1)) * tmp_698;
   std::complex<double> tmp_700;
   std::complex<double> tmp_701;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_702;
      std::complex<double> tmp_703;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_703 += KroneckerDelta(gO2,3 + j1)*TYu(j1,j2);
      }
      tmp_702 += tmp_703;
      tmp_701 += (Conj(ZU(gI1,j2))) * tmp_702;
   }
   tmp_700 += tmp_701;
   result += (-0.7071067811865475*ZH(gI2,1)) * tmp_700;
   std::complex<double> tmp_704;
   std::complex<double> tmp_705;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_706;
      std::complex<double> tmp_707;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_708;
         std::complex<double> tmp_709;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_709 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_708 += tmp_709;
         tmp_707 += (KroneckerDelta(gO2,3 + j2)) * tmp_708;
      }
      tmp_706 += tmp_707;
      tmp_705 += (Conj(ZU(gI1,3 + j3))) * tmp_706;
   }
   tmp_704 += tmp_705;
   result += (-(vu*ZH(gI2,1))) * tmp_704;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_710;
      std::complex<double> tmp_711;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_711 += Conj(ZU(gI1,3 + j1))*Conj(TYu(j1,gO2));
      }
      tmp_710 += tmp_711;
      result += (-0.7071067811865475*ZH(gI2,1)) * tmp_710;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_712;
      std::complex<double> tmp_713;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_714;
         std::complex<double> tmp_715;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_715 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_714 += tmp_715;
         tmp_713 += (Conj(ZU(gI1,j2))) * tmp_714;
      }
      tmp_712 += tmp_713;
      result += (-(vu*ZH(gI2,1))) * tmp_712;
   }
   std::complex<double> tmp_716;
   std::complex<double> tmp_717;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_718;
      std::complex<double> tmp_719;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_719 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_718 += tmp_719;
      tmp_717 += (Conj(ZU(gI1,j2))) * tmp_718;
   }
   tmp_716 += tmp_717;
   result += (0.35355339059327373*vd*Conj(Lambdax)*ZH(gI2,2)) * tmp_716;
   if (gO2 < 3) {
      std::complex<double> tmp_720;
      std::complex<double> tmp_721;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_721 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_720 += tmp_721;
      result += (0.35355339059327373*vd*Lambdax*ZH(gI2,2)) * tmp_720;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_722;
   std::complex<double> tmp_723;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_723 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_722 += tmp_723;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_722;

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

   std::complex<double> tmp_724;
   std::complex<double> tmp_725;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_725 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_724 += tmp_725;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_724;
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
   std::complex<double> result;

   std::complex<double> tmp_726;
   std::complex<double> tmp_727;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_727 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_726 += tmp_727;
   result += (-0.5163977794943222*g1*Sin(ThetaW())) * tmp_726;
   if (gO2 < 3) {
      result += 0.5*g2*Conj(ZU(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(ZU(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_728;
   std::complex<double> tmp_729;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_729 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_728 += tmp_729;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_728;
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*Sqr(Cos(ThetaW()));
   }
   if (gO1 < 3) {
      result += -0.7745966692414834*g1*g2*Cos(ThetaW())*KroneckerDelta(gO1,
         gO2)*Sin(ThetaW());
   }
   if (gO1 < 3) {
      result += 0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)*Sqr(Sin(ThetaW()));
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

std::complex<double> CLASSNAME::CpUSeconjUSeAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_730;
   std::complex<double> tmp_731;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_731 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_730 += tmp_731;
   result += (0.3*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_730;
   std::complex<double> tmp_732;
   std::complex<double> tmp_733;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_734;
      std::complex<double> tmp_735;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_736;
         std::complex<double> tmp_737;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_737 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_736 += tmp_737;
         tmp_735 += (KroneckerDelta(gO2,3 + j2)) * tmp_736;
      }
      tmp_734 += tmp_735;
      tmp_733 += (KroneckerDelta(gO1,3 + j3)) * tmp_734;
   }
   tmp_732 += tmp_733;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_732;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_738;
      std::complex<double> tmp_739;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_739 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_738 += tmp_739;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_738;
   }
   std::complex<double> tmp_740;
   std::complex<double> tmp_741;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_741 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_740 += tmp_741;
   result += (-0.3*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_740;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_742;
      std::complex<double> tmp_743;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_743 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_742 += tmp_743;
      result += (-0.35355339059327373*Conj(Lambdax)*ZA(gI1,2)*ZA(gI2,1)) *
         tmp_742;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_744;
      std::complex<double> tmp_745;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_745 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_744 += tmp_745;
      result += (-0.35355339059327373*Lambdax*ZA(gI1,2)*ZA(gI2,1)) * tmp_744
         ;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_746;
      std::complex<double> tmp_747;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_747 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_746 += tmp_747;
      result += (-0.35355339059327373*Conj(Lambdax)*ZA(gI1,1)*ZA(gI2,2)) *
         tmp_746;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_748;
      std::complex<double> tmp_749;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_749 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_748 += tmp_749;
      result += (-0.35355339059327373*Lambdax*ZA(gI1,1)*ZA(gI2,2)) * tmp_748
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_750;
   std::complex<double> tmp_751;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_751 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_750 += tmp_751;
   result += (0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_750;
   std::complex<double> tmp_752;
   std::complex<double> tmp_754;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_755;
      std::complex<double> tmp_756;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_756 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_755 += tmp_756;
      tmp_754 += (Conj(ZV(gI2,j2))) * tmp_755;
   }
   tmp_752 += tmp_754;
   std::complex<double> tmp_753;
   std::complex<double> tmp_757;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_758;
      std::complex<double> tmp_759;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_759 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_758 += tmp_759;
      tmp_757 += (ZV(gI1,j4)) * tmp_758;
   }
   tmp_753 += tmp_757;
   result += (-1) * tmp_752 * tmp_753;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g1
         );
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gI1,gI2)*KroneckerDelta(gO1,gO2)*Sqr(g2)
         ;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZV(gI2,gO2))*Sqr(g2)*ZV(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_760;
   std::complex<double> tmp_761;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_761 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_760 += tmp_761;
   result += (0.3*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_760;
   std::complex<double> tmp_762;
   std::complex<double> tmp_763;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_764;
      std::complex<double> tmp_765;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_766;
         std::complex<double> tmp_767;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_767 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_766 += tmp_767;
         tmp_765 += (KroneckerDelta(gO2,3 + j2)) * tmp_766;
      }
      tmp_764 += tmp_765;
      tmp_763 += (KroneckerDelta(gO1,3 + j3)) * tmp_764;
   }
   tmp_762 += tmp_763;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_762;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_768;
      std::complex<double> tmp_769;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_769 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_768 += tmp_769;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_768;
   }
   std::complex<double> tmp_770;
   std::complex<double> tmp_771;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_771 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_770 += tmp_771;
   result += (-0.3*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_770;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_772;
      std::complex<double> tmp_773;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_773 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_772 += tmp_773;
      result += (0.35355339059327373*Conj(Lambdax)*ZH(gI1,2)*ZH(gI2,1)) *
         tmp_772;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_774;
      std::complex<double> tmp_775;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_775 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_774 += tmp_775;
      result += (0.35355339059327373*Lambdax*ZH(gI1,2)*ZH(gI2,1)) * tmp_774;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_776;
      std::complex<double> tmp_777;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_777 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_776 += tmp_777;
      result += (0.35355339059327373*Conj(Lambdax)*ZH(gI1,1)*ZH(gI2,2)) *
         tmp_776;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_778;
      std::complex<double> tmp_779;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_779 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_778 += tmp_779;
      result += (0.35355339059327373*Lambdax*ZH(gI1,1)*ZH(gI2,2)) * tmp_778;
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

   std::complex<double> tmp_780;
   std::complex<double> tmp_781;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_781 += KroneckerDelta(gO1,3 + j1)*Ye(j1,gI1);
   }
   tmp_780 += tmp_781;
   result += (Conj(UM(gI2,1))) * tmp_780;
   if (gI1 < 3) {
      result += -(g2*Conj(UM(gI2,0))*KroneckerDelta(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_782;
   std::complex<double> tmp_783;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_784;
      std::complex<double> tmp_785;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_785 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_784 += tmp_785;
      tmp_783 += (Conj(ZV(gI1,j2))) * tmp_784;
   }
   tmp_782 += tmp_783;
   result += (ZP(gI2,0)) * tmp_782;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_786;
      std::complex<double> tmp_787;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_788;
         std::complex<double> tmp_789;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_789 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_788 += tmp_789;
         tmp_787 += (Conj(ZV(gI1,j2))) * tmp_788;
      }
      tmp_786 += tmp_787;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_786;
   }
   std::complex<double> tmp_790;
   std::complex<double> tmp_791;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_792;
      std::complex<double> tmp_793;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_793 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_792 += tmp_793;
      tmp_791 += (Conj(ZV(gI1,j2))) * tmp_792;
   }
   tmp_790 += tmp_791;
   result += (-0.5*vT*Conj(Lambdax)*ZP(gI2,1)) * tmp_790;
   std::complex<double> tmp_794;
   std::complex<double> tmp_795;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_796;
      std::complex<double> tmp_797;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_797 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_796 += tmp_797;
      tmp_795 += (Conj(ZV(gI1,j2))) * tmp_796;
   }
   tmp_794 += tmp_795;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_794;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,2);
   }
   std::complex<double> tmp_798;
   std::complex<double> tmp_799;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_800;
      std::complex<double> tmp_801;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_801 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_800 += tmp_801;
      tmp_799 += (Conj(ZV(gI1,j2))) * tmp_800;
   }
   tmp_798 += tmp_799;
   result += (-0.7071067811865475*vu*Conj(Lambdax)*ZP(gI2,3)) * tmp_798;
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_802;
   std::complex<double> tmp_803;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_803 += KroneckerDelta(gO2,3 + j1)*ZER(gI1,j1);
   }
   tmp_802 += tmp_803;
   result += (-1.0954451150103321*g1*ZN(gI2,0)) * tmp_802;
   if (gO2 < 3) {
      std::complex<double> tmp_804;
      std::complex<double> tmp_805;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_805 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_804 += tmp_805;
      result += (-ZN(gI2,2)) * tmp_804;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_806;
   std::complex<double> tmp_807;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_808;
      std::complex<double> tmp_809;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_809 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_808 += tmp_809;
      tmp_807 += (Conj(ZEL(gI1,j2))) * tmp_808;
   }
   tmp_806 += tmp_807;
   result += (-Conj(ZN(gI2,2))) * tmp_806;
   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZEL(gI1,gO1))*Conj(ZN(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_810;
   std::complex<double> tmp_811;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_811 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_810 += tmp_811;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_810;
   std::complex<double> tmp_812;
   std::complex<double> tmp_813;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_814;
      std::complex<double> tmp_815;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_816;
         std::complex<double> tmp_817;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_817 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_816 += tmp_817;
         tmp_815 += (KroneckerDelta(gO2,3 + j2)) * tmp_816;
      }
      tmp_814 += tmp_815;
      tmp_813 += (KroneckerDelta(gO1,3 + j3)) * tmp_814;
   }
   tmp_812 += tmp_813;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_812;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_818;
   std::complex<double> tmp_819;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_819 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_818 += tmp_819;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_818;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_820;
      std::complex<double> tmp_821;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_821 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_820 += tmp_821;
      result += (-(Conj(Lambdax)*ZP(gI1,2)*ZP(gI2,1))) * tmp_820;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_822;
      std::complex<double> tmp_823;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_823 += Conj(Ye(j1,gO2))*KroneckerDelta(gO1,3 + j1);
      }
      tmp_822 += tmp_823;
      result += (-(Lambdax*ZP(gI1,1)*ZP(gI2,2))) * tmp_822;
   }
   if (gO1 < 3) {
      result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_824;
   std::complex<double> tmp_826;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_826 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_824 += tmp_826;
   std::complex<double> tmp_825;
   std::complex<double> tmp_827;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_827 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_825 += tmp_827;
   result += (-0.05*Sqr(g1)) * tmp_824 * tmp_825;
   std::complex<double> tmp_828;
   std::complex<double> tmp_830;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_830 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_828 += tmp_830;
   std::complex<double> tmp_829;
   std::complex<double> tmp_831;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_831 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_829 += tmp_831;
   result += (-0.1*Sqr(g1)) * tmp_828 * tmp_829;
   std::complex<double> tmp_832;
   std::complex<double> tmp_834;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_834 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_832 += tmp_834;
   std::complex<double> tmp_833;
   std::complex<double> tmp_835;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_835 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_833 += tmp_835;
   result += (-0.05*Sqr(g1)) * tmp_832 * tmp_833;
   std::complex<double> tmp_836;
   std::complex<double> tmp_838;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_838 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_836 += tmp_838;
   std::complex<double> tmp_837;
   std::complex<double> tmp_839;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_839 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_837 += tmp_839;
   result += (-0.1*Sqr(g1)) * tmp_836 * tmp_837;
   if (gO1 < 3) {
      std::complex<double> tmp_840;
      std::complex<double> tmp_841;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_841 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_840 += tmp_841;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_840;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_842;
      std::complex<double> tmp_843;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_843 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_842 += tmp_843;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_842;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_844;
      std::complex<double> tmp_845;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_845 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_844 += tmp_845;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_844;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_846;
      std::complex<double> tmp_847;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_847 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_846 += tmp_847;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_846;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_848;
      std::complex<double> tmp_849;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_849 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_848 += tmp_849;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_848;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_850;
      std::complex<double> tmp_851;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_851 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_850 += tmp_851;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_850;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_852;
      std::complex<double> tmp_854;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_854 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_852 += tmp_854;
      std::complex<double> tmp_853;
      std::complex<double> tmp_855;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_856;
         std::complex<double> tmp_857;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_857 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_856 += tmp_857;
         tmp_855 += (ZD(gI1,j4)) * tmp_856;
      }
      tmp_853 += tmp_855;
      result += (-1) * tmp_852 * tmp_853;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_858;
      std::complex<double> tmp_860;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_861;
         std::complex<double> tmp_862;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_862 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_861 += tmp_862;
         tmp_860 += (Conj(ZD(gI2,j2))) * tmp_861;
      }
      tmp_858 += tmp_860;
      std::complex<double> tmp_859;
      std::complex<double> tmp_863;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_863 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_859 += tmp_863;
      result += (-1) * tmp_858 * tmp_859;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_864;
   std::complex<double> tmp_866;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_866 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_864 += tmp_866;
   std::complex<double> tmp_865;
   std::complex<double> tmp_867;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_867 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_865 += tmp_867;
   result += (-0.3*Sqr(g1)) * tmp_864 * tmp_865;
   std::complex<double> tmp_868;
   std::complex<double> tmp_870;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_870 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_868 += tmp_870;
   std::complex<double> tmp_869;
   std::complex<double> tmp_871;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_871 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_869 += tmp_871;
   result += (0.15*Sqr(g1)) * tmp_868 * tmp_869;
   std::complex<double> tmp_872;
   std::complex<double> tmp_874;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_874 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_872 += tmp_874;
   std::complex<double> tmp_873;
   std::complex<double> tmp_875;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_875 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_873 += tmp_875;
   result += (-0.3*Sqr(g1)) * tmp_872 * tmp_873;
   std::complex<double> tmp_876;
   std::complex<double> tmp_878;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_878 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_876 += tmp_878;
   std::complex<double> tmp_877;
   std::complex<double> tmp_879;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_879 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_877 += tmp_879;
   result += (0.15*Sqr(g1)) * tmp_876 * tmp_877;
   std::complex<double> tmp_880;
   std::complex<double> tmp_882;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_882 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_880 += tmp_882;
   std::complex<double> tmp_881;
   std::complex<double> tmp_883;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_883 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_881 += tmp_883;
   result += (-0.3*Sqr(g1)) * tmp_880 * tmp_881;
   std::complex<double> tmp_884;
   std::complex<double> tmp_886;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_886 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_884 += tmp_886;
   std::complex<double> tmp_885;
   std::complex<double> tmp_887;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_887 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_885 += tmp_887;
   result += (-0.3*Sqr(g1)) * tmp_884 * tmp_885;
   std::complex<double> tmp_888;
   std::complex<double> tmp_890;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_891;
      std::complex<double> tmp_892;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_892 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_891 += tmp_892;
      tmp_890 += (Conj(ZE(gI2,j2))) * tmp_891;
   }
   tmp_888 += tmp_890;
   std::complex<double> tmp_889;
   std::complex<double> tmp_893;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_894;
      std::complex<double> tmp_895;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_895 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_894 += tmp_895;
      tmp_893 += (ZE(gI1,j4)) * tmp_894;
   }
   tmp_889 += tmp_893;
   result += (-1) * tmp_888 * tmp_889;
   if (gO1 < 3) {
      std::complex<double> tmp_896;
      std::complex<double> tmp_897;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_897 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_896 += tmp_897;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_896;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_898;
      std::complex<double> tmp_899;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_899 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_898 += tmp_899;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_898;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_900;
      std::complex<double> tmp_901;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_901 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_900 += tmp_901;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_900;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_902;
      std::complex<double> tmp_903;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_903 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_902 += tmp_903;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_902;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_904;
      std::complex<double> tmp_905;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_905 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_904 += tmp_905;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_904;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_906;
      std::complex<double> tmp_907;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_907 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_906 += tmp_907;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_906;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_908;
      std::complex<double> tmp_910;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_910 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_908 += tmp_910;
      std::complex<double> tmp_909;
      std::complex<double> tmp_911;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_912;
         std::complex<double> tmp_913;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_913 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_912 += tmp_913;
         tmp_911 += (ZE(gI1,j4)) * tmp_912;
      }
      tmp_909 += tmp_911;
      result += (-1) * tmp_908 * tmp_909;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_914;
      std::complex<double> tmp_915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_915 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_914 += tmp_915;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_914;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_916;
      std::complex<double> tmp_917;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_917 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_916 += tmp_917;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_916;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_918;
      std::complex<double> tmp_920;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_921;
         std::complex<double> tmp_922;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_922 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_921 += tmp_922;
         tmp_920 += (Conj(ZE(gI2,j2))) * tmp_921;
      }
      tmp_918 += tmp_920;
      std::complex<double> tmp_919;
      std::complex<double> tmp_923;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_923 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_919 += tmp_923;
      result += (-1) * tmp_918 * tmp_919;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_924;
      std::complex<double> tmp_926;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_926 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_924 += tmp_926;
      std::complex<double> tmp_925;
      std::complex<double> tmp_927;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_927 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_925 += tmp_927;
      result += (-1) * tmp_924 * tmp_925;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_928;
      std::complex<double> tmp_929;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_929 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_928 += tmp_929;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_928;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_930;
      std::complex<double> tmp_931;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_931 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_930 += tmp_931;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_930;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.15*Conj(ZE(gI2,gO2))*Sqr(g1)*ZE(gI1,gO1);
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.25*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_932;
   std::complex<double> tmp_934;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_934 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_932 += tmp_934;
   std::complex<double> tmp_933;
   std::complex<double> tmp_935;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_935 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_933 += tmp_935;
   result += (-0.05*Sqr(g1)) * tmp_932 * tmp_933;
   std::complex<double> tmp_936;
   std::complex<double> tmp_938;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_938 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_936 += tmp_938;
   std::complex<double> tmp_937;
   std::complex<double> tmp_939;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_939 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_937 += tmp_939;
   result += (0.2*Sqr(g1)) * tmp_936 * tmp_937;
   std::complex<double> tmp_940;
   std::complex<double> tmp_942;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_942 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_940 += tmp_942;
   std::complex<double> tmp_941;
   std::complex<double> tmp_943;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_943 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_941 += tmp_943;
   result += (-0.05*Sqr(g1)) * tmp_940 * tmp_941;
   std::complex<double> tmp_944;
   std::complex<double> tmp_946;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_946 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_944 += tmp_946;
   std::complex<double> tmp_945;
   std::complex<double> tmp_947;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_947 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_945 += tmp_947;
   result += (0.2*Sqr(g1)) * tmp_944 * tmp_945;
   if (gO1 < 3) {
      std::complex<double> tmp_948;
      std::complex<double> tmp_949;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_949 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_948 += tmp_949;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_948;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_950;
      std::complex<double> tmp_951;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_951 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_950 += tmp_951;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_950;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_952;
      std::complex<double> tmp_953;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_953 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_952 += tmp_953;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_952;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_954;
      std::complex<double> tmp_955;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_955 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_954 += tmp_955;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_954;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_956;
      std::complex<double> tmp_957;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_957 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_956 += tmp_957;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_956;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_958;
      std::complex<double> tmp_959;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_959 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_958 += tmp_959;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_958;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_960;
   std::complex<double> tmp_961;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_962;
      std::complex<double> tmp_963;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_963 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_962 += tmp_963;
      tmp_961 += (Conj(ZE(gI1,j2))) * tmp_962;
   }
   tmp_960 += tmp_961;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) * tmp_960;
   if (gO2 < 3) {
      std::complex<double> tmp_964;
      std::complex<double> tmp_965;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_965 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_964 += tmp_965;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) *
         tmp_964;
   }
   std::complex<double> tmp_966;
   std::complex<double> tmp_967;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_968;
      std::complex<double> tmp_969;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_969 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_968 += tmp_969;
      tmp_967 += (Conj(ZE(gI1,j2))) * tmp_968;
   }
   tmp_966 += tmp_967;
   result += (std::complex<double>(0,-0.35355339059327373)*vT*Conj(Lambdax)*ZA(
      gI2,1)) * tmp_966;
   std::complex<double> tmp_970;
   std::complex<double> tmp_971;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_972;
      std::complex<double> tmp_973;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_973 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_972 += tmp_973;
      tmp_971 += (Conj(ZE(gI1,j2))) * tmp_972;
   }
   tmp_970 += tmp_971;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1)) *
      tmp_970;
   if (gO2 < 3) {
      std::complex<double> tmp_974;
      std::complex<double> tmp_975;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_975 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_974 += tmp_975;
      result += (std::complex<double>(0,0.35355339059327373)*vT*Lambdax*ZA(
         gI2,1)) * tmp_974;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_976;
      std::complex<double> tmp_977;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_977 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_976 += tmp_977;
      result += (std::complex<double>(0,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_976;
   }
   std::complex<double> tmp_978;
   std::complex<double> tmp_979;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_980;
      std::complex<double> tmp_981;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_981 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_980 += tmp_981;
      tmp_979 += (Conj(ZE(gI1,j2))) * tmp_980;
   }
   tmp_978 += tmp_979;
   result += (std::complex<double>(0,-0.35355339059327373)*vu*Conj(Lambdax)*ZA(
      gI2,2)) * tmp_978;
   if (gO2 < 3) {
      std::complex<double> tmp_982;
      std::complex<double> tmp_983;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_983 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_982 += tmp_983;
      result += (std::complex<double>(0,0.35355339059327373)*vu*Lambdax*ZA(
         gI2,2)) * tmp_982;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_984;
   std::complex<double> tmp_985;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_985 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_984 += tmp_985;
   result += (0.3*vd*Sqr(g1)*ZH(gI2,0)) * tmp_984;
   std::complex<double> tmp_986;
   std::complex<double> tmp_987;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_988;
      std::complex<double> tmp_989;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_989 += KroneckerDelta(gO2,3 + j1)*TYe(j1,j2);
      }
      tmp_988 += tmp_989;
      tmp_987 += (Conj(ZE(gI1,j2))) * tmp_988;
   }
   tmp_986 += tmp_987;
   result += (-0.7071067811865475*ZH(gI2,0)) * tmp_986;
   std::complex<double> tmp_990;
   std::complex<double> tmp_991;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_992;
      std::complex<double> tmp_993;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_994;
         std::complex<double> tmp_995;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_995 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_994 += tmp_995;
         tmp_993 += (KroneckerDelta(gO2,3 + j2)) * tmp_994;
      }
      tmp_992 += tmp_993;
      tmp_991 += (Conj(ZE(gI1,3 + j3))) * tmp_992;
   }
   tmp_990 += tmp_991;
   result += (-(vd*ZH(gI2,0))) * tmp_990;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_996;
      std::complex<double> tmp_997;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_997 += Conj(ZE(gI1,3 + j1))*Conj(TYe(j1,gO2));
      }
      tmp_996 += tmp_997;
      result += (-0.7071067811865475*ZH(gI2,0)) * tmp_996;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_998;
      std::complex<double> tmp_999;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1000;
         std::complex<double> tmp_1001;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1001 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_1000 += tmp_1001;
         tmp_999 += (Conj(ZE(gI1,j2))) * tmp_1000;
      }
      tmp_998 += tmp_999;
      result += (-(vd*ZH(gI2,0))) * tmp_998;
   }
   std::complex<double> tmp_1002;
   std::complex<double> tmp_1003;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1003 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1002 += tmp_1003;
   result += (-0.3*vu*Sqr(g1)*ZH(gI2,1)) * tmp_1002;
   std::complex<double> tmp_1004;
   std::complex<double> tmp_1005;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1006;
      std::complex<double> tmp_1007;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1007 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1006 += tmp_1007;
      tmp_1005 += (Conj(ZE(gI1,j2))) * tmp_1006;
   }
   tmp_1004 += tmp_1005;
   result += (0.35355339059327373*vT*Conj(Lambdax)*ZH(gI2,1)) * tmp_1004;
   std::complex<double> tmp_1008;
   std::complex<double> tmp_1009;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1010;
      std::complex<double> tmp_1011;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1011 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1010 += tmp_1011;
      tmp_1009 += (Conj(ZE(gI1,j2))) * tmp_1010;
   }
   tmp_1008 += tmp_1009;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_1008;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1012;
      std::complex<double> tmp_1013;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1013 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1012 += tmp_1013;
      result += (0.35355339059327373*vT*Lambdax*ZH(gI2,1)) * tmp_1012;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1014;
      std::complex<double> tmp_1015;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1015 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1014 += tmp_1015;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_1014;
   }
   std::complex<double> tmp_1016;
   std::complex<double> tmp_1017;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1018;
      std::complex<double> tmp_1019;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1019 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_1018 += tmp_1019;
      tmp_1017 += (Conj(ZE(gI1,j2))) * tmp_1018;
   }
   tmp_1016 += tmp_1017;
   result += (0.35355339059327373*vu*Conj(Lambdax)*ZH(gI2,2)) * tmp_1016;
   if (gO2 < 3) {
      std::complex<double> tmp_1020;
      std::complex<double> tmp_1021;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1021 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1020 += tmp_1021;
      result += (0.35355339059327373*vu*Lambdax*ZH(gI2,2)) * tmp_1020;
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

   std::complex<double> tmp_1022;
   std::complex<double> tmp_1023;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1023 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1022 += tmp_1023;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_1022;
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
   std::complex<double> result;

   std::complex<double> tmp_1024;
   std::complex<double> tmp_1025;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1025 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_1024 += tmp_1025;
   result += (0.7745966692414834*g1*Sin(ThetaW())) * tmp_1024;
   if (gO2 < 3) {
      result += -0.5*g2*Conj(ZE(gI2,gO2))*Cos(ThetaW());
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(ZE(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZVZ(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(vd*KroneckerDelta(0,gO2) + vu*KroneckerDelta(1,gO2))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(vd*KroneckerDelta(0,gO2) + vu*KroneckerDelta(1,gO2) + 4*vT*
      KroneckerDelta(2,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1) + 4*vT*
      KroneckerDelta(2,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1) + 4*vT*
      KroneckerDelta(2,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargZgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.025*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1))*(
      7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(g2) + Cos(2*
      ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

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

std::complex<double> CLASSNAME::CpUhhUhhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2) + 4*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2))
      *Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,gO1)*KroneckerDelta(2,
      gO2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,
      0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) - 10*AbsSqr(
      Lambdax)*ZA(gI1,2)*ZA(gI2,2)) - KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*
      ((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (10*AbsSqr(Lambdax) - 3*Sqr(
      g1) - 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2
      )));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,gO1)*(KroneckerDelta(2,
      gO2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) + KroneckerDelta(0,gO2)*(ZH
      (gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*(ZH(gI1,2)*
      ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2))) + KroneckerDelta(1,gO1)*(KroneckerDelta(0,
      gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH
      (gI1,0)*ZH(gI2,1)) - 10*AbsSqr(Lambdax)*KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH(
      gI2,1) + ZH(gI1,1)*ZH(gI2,2)) + KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax)
      + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZH(
      gI1,1)*ZH(gI2,1) - 10*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2))) - KroneckerDelta
      (0,gO1)*(-(KroneckerDelta(1,gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2
      ))*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1))) + 10*AbsSqr(Lambdax)*
      KroneckerDelta(2,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) +
      KroneckerDelta(0,gO2)*(3*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (10*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 10*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*KroneckerDelta(2,gO2)*((Conj(TLambdax) + 2*Conj(MT)*
      Lambdax + TLambdax)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + 2*Conj(Mu)
      *Lambdax*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1)) + 2*Conj(Lambdax)*(ZA(
      gI1,0)*((vT*Lambdax + Mu)*ZA(gI2,0) + MT*ZA(gI2,1)) + ZA(gI1,1)*(MT*ZA(gI2,0
      ) + (vT*Lambdax + Mu)*ZA(gI2,1)))) + KroneckerDelta(1,gO2)*(10*MT*Conj(
      Lambdax)*ZA(gI1,2)*ZA(gI2,0) - 5*Conj(TLambdax)*ZA(gI1,2)*ZA(gI2,0) + 10*
      Conj(MT)*Lambdax*ZA(gI1,2)*ZA(gI2,0) - 5*TLambdax*ZA(gI1,2)*ZA(gI2,0) - 3*vu
      *Sqr(g1)*ZA(gI1,1)*ZA(gI2,1) - 5*vu*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1) - 10*vu*
      AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2) + ZA(gI1,0)*(vu*(-10*AbsSqr(Lambdax) + 3
      *Sqr(g1) + 5*Sqr(g2))*ZA(gI2,0) + 5*(2*MT*Conj(Lambdax) - Conj(TLambdax) + 2
      *Conj(MT)*Lambdax - TLambdax)*ZA(gI2,2))) - KroneckerDelta(0,gO2)*(vd*(3*Sqr
      (g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*(vd*(10*AbsSqr(Lambdax) -
      3*Sqr(g1) - 5*Sqr(g2))*ZA(gI2,1) + 5*(-2*MT*Conj(Lambdax) + Conj(TLambdax) -
      2*Conj(MT)*Lambdax + TLambdax)*ZA(gI2,2)) + 5*ZA(gI1,2)*(-2*Conj(MT)*
      Lambdax*ZA(gI2,1) + (Conj(TLambdax) + TLambdax)*ZA(gI2,1) + Conj(Lambdax)*(
      -2*MT*ZA(gI2,1) + 2*vd*Lambdax*ZA(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*
      KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(-(Conj(TLambdax)*KroneckerDelta(2,gO2
      )*ZA(gI2,1)*ZH(gI1,0)) + KroneckerDelta(2,gO2)*TLambdax*ZA(gI2,1)*ZH(gI1,0)
      - Conj(TLambdax)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,0) - 2*Conj(Mu)*
      KroneckerDelta(0,gO2)*Lambdax*ZA(gI2,2)*ZH(gI1,0) + KroneckerDelta(1,gO2)*
      TLambdax*ZA(gI2,2)*ZH(gI1,0) - Conj(TLambdax)*KroneckerDelta(2,gO2)*ZA(gI2,0
      )*ZH(gI1,1) + KroneckerDelta(2,gO2)*TLambdax*ZA(gI2,0)*ZH(gI1,1) - Conj(
      TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZH(gI1,1) - 2*Conj(Mu)*
      KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,2)*ZH(gI1,1) + KroneckerDelta(0,gO2)*
      TLambdax*ZA(gI2,2)*ZH(gI1,1) - Conj(TLambdax)*KroneckerDelta(1,gO2)*ZA(gI2,0
      )*ZH(gI1,2) + KroneckerDelta(1,gO2)*TLambdax*ZA(gI2,0)*ZH(gI1,2) - Conj(
      TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,1)*ZH(gI1,2) + KroneckerDelta(0,gO2)*
      TLambdax*ZA(gI2,1)*ZH(gI1,2) + 2*Conj(MT)*Lambdax*(KroneckerDelta(2,gO2)*(ZA
      (gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) + KroneckerDelta(1,gO2)*(-(ZA(gI2,2
      )*ZH(gI1,0)) + ZA(gI2,0)*ZH(gI1,2)) + KroneckerDelta(0,gO2)*(-(ZA(gI2,2)*ZH(
      gI1,1)) + ZA(gI2,1)*ZH(gI1,2))) - 2*Conj(Lambdax)*(MT*KroneckerDelta(2,gO2)*
      (ZA(gI2,1)*ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,1)) - KroneckerDelta(1,gO2)*(ZA(gI2,
      2)*(MT*ZH(gI1,0) + Mu*ZH(gI1,1)) - MT*ZA(gI2,0)*ZH(gI1,2)) - KroneckerDelta(
      0,gO2)*(ZA(gI2,2)*(Mu*ZH(gI1,0) + MT*ZH(gI1,1)) - MT*ZA(gI2,1)*ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO2)*(-(ZH(gI1,0)*(3*vd*(3*Sqr(g1) + 5*Sqr(
      g2))*ZH(gI2,0) + vu*(10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI2,1) +
      10*(Conj(Mu)*Lambdax + Conj(Lambdax)*(vT*Lambdax + Mu))*ZH(gI2,2))) + ZH(
      gI1,1)*(vu*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,0) + vd*(-10
      *AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) + 5*(2*MT*Conj(Lambdax)
      + Conj(TLambdax) + 2*Conj(MT)*Lambdax + TLambdax)*ZH(gI2,2)) - 5*ZH(gI1,2)*(
      2*Conj(Mu)*Lambdax*ZH(gI2,0) - (Conj(TLambdax) + 2*Conj(MT)*Lambdax +
      TLambdax)*ZH(gI2,1) + 2*Conj(Lambdax)*((vT*Lambdax + Mu)*ZH(gI2,0) - MT*ZH(
      gI2,1) + vd*Lambdax*ZH(gI2,2)))) + KroneckerDelta(1,gO2)*(ZH(gI1,1)*(vd*(-10
      *AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,0) - 3*vu*(3*Sqr(g1) + 5*
      Sqr(g2))*ZH(gI2,1) - 10*(Conj(Mu)*Lambdax + Conj(Lambdax)*(vT*Lambdax + Mu))
      *ZH(gI2,2)) + ZH(gI1,0)*(vu*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH
      (gI2,0) + vd*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) + 5*(2*
      MT*Conj(Lambdax) + Conj(TLambdax) + 2*Conj(MT)*Lambdax + TLambdax)*ZH(gI2,2)
      ) + 5*ZH(gI1,2)*(Conj(TLambdax)*ZH(gI2,0) + 2*Conj(MT)*Lambdax*ZH(gI2,0) +
      TLambdax*ZH(gI2,0) - 2*Conj(Mu)*Lambdax*ZH(gI2,1) + 2*Conj(Lambdax)*(MT*ZH(
      gI2,0) - (vT*Lambdax + Mu)*ZH(gI2,1) - vu*Lambdax*ZH(gI2,2)))) - 5*
      KroneckerDelta(2,gO2)*(-((Conj(TLambdax) + 2*Conj(MT)*Lambdax + TLambdax)*(
      ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1))) + 2*Conj(Mu)*Lambdax*(ZH(gI1,0)*
      ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) + 2*Conj(Lambdax)*(Lambdax*ZH(gI1,2)*(vd*ZH
      (gI2,0) + vu*ZH(gI2,1)) + ZH(gI1,0)*((vT*Lambdax + Mu)*ZH(gI2,0) - MT*ZH(gI2
      ,1) + vd*Lambdax*ZH(gI2,2)) + ZH(gI1,1)*(-(MT*ZH(gI2,0)) + (vT*Lambdax + Mu)
      *ZH(gI2,1) + vu*Lambdax*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-1.4142135623730951*KroneckerDelta(0,gO2)*(g2*UM(gI1,1)*UP(gI2
      ,0) + Conj(Lambdax)*UM(gI1,2)*UP(gI2,1)) + KroneckerDelta(2,gO2)*(-2*g2*UM(
      gI1,2)*UP(gI2,0) + Conj(Lambdax)*UM(gI1,1)*UP(gI2,1) + 2*g2*UM(gI1,0)*UP(gI2
      ,2)) + 1.4142135623730951*KroneckerDelta(1,gO2)*(-(g2*UM(gI1,0)*UP(gI2,1)) +
      Conj(Lambdax)*UM(gI1,1)*UP(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Conj(UM(gI2,0))*(-1.4142135623730951*Conj(UP(gI1,1))*
      KroneckerDelta(1,gO1) + 2*Conj(UP(gI1,2))*KroneckerDelta(2,gO1)) - Conj(UM(
      gI2,2))*(2*g2*Conj(UP(gI1,0))*KroneckerDelta(2,gO1) + 1.4142135623730951*
      Conj(UP(gI1,1))*KroneckerDelta(0,gO1)*Lambdax) + Conj(UM(gI2,1))*(
      -1.4142135623730951*g2*Conj(UP(gI1,0))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(UP(gI1,2))*KroneckerDelta(1,gO1)*Lambdax + Conj(UP(
      gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1026;
   std::complex<double> tmp_1027;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1028;
      std::complex<double> tmp_1029;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1029 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1028 += tmp_1029;
      tmp_1027 += (ZDL(gI1,j2)) * tmp_1028;
   }
   tmp_1026 += tmp_1027;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1026;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1030;
   std::complex<double> tmp_1031;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1032;
      std::complex<double> tmp_1033;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1033 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1032 += tmp_1033;
      tmp_1031 += (Conj(ZDL(gI2,j2))) * tmp_1032;
   }
   tmp_1030 += tmp_1031;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1030;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1034;
   std::complex<double> tmp_1035;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1036;
      std::complex<double> tmp_1037;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1037 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1036 += tmp_1037;
      tmp_1035 += (ZEL(gI1,j2)) * tmp_1036;
   }
   tmp_1034 += tmp_1035;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1034;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1038;
   std::complex<double> tmp_1039;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1040;
      std::complex<double> tmp_1041;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1041 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1040 += tmp_1041;
      tmp_1039 += (Conj(ZEL(gI2,j2))) * tmp_1040;
   }
   tmp_1038 += tmp_1039;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_1038;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1042;
   std::complex<double> tmp_1043;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1044;
      std::complex<double> tmp_1045;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1045 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1044 += tmp_1045;
      tmp_1043 += (ZUL(gI1,j2)) * tmp_1044;
   }
   tmp_1042 += tmp_1043;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1042;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1046;
   std::complex<double> tmp_1047;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1048;
      std::complex<double> tmp_1049;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1049 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1048 += tmp_1049;
      tmp_1047 += (Conj(ZUL(gI2,j2))) * tmp_1048;
   }
   tmp_1046 += tmp_1047;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_1046;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(5*(KroneckerDelta(1,gO2)*(AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      1.4142135623730951*KroneckerDelta(2,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZP(
      gI1,2)*ZP(gI2,0) - ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) - ZP(gI2,3))))
      + KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (20*
      AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*
      AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3
      ))))) + KroneckerDelta(1,gO1)*(-5*(KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) +
      Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + 1.4142135623730951*
      KroneckerDelta(2,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,2)*ZP(gI2,1) - ZP
      (gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) - ZP(gI2,3)))) + KroneckerDelta(1,
      gO2)*((-20*AbsSqr(Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - (3
      *Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2)
      - (-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3)))) - 5*KroneckerDelta(2
      ,gO1)*(Sqr(g2)*(1.4142135623730951*KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0
      ) - ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) - ZP(gI2,3))) +
      1.4142135623730951*KroneckerDelta(1,gO2)*(ZP(gI1,2)*ZP(gI2,1) - ZP(gI1,3)*ZP
      (gI2,1) + ZP(gI1,1)*(ZP(gI2,2) - ZP(gI2,3))) + 4*KroneckerDelta(2,gO2)*(ZP(
      gI1,2) - ZP(gI1,3))*(ZP(gI2,2) - ZP(gI2,3))) + AbsSqr(Lambdax)*(2*
      KroneckerDelta(2,gO2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1)) +
      1.4142135623730951*(KroneckerDelta(0,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,3
      )*ZP(gI2,0) + ZP(gI1,0)*(-ZP(gI2,2) + ZP(gI2,3))) + KroneckerDelta(1,gO2)*(-
      (ZP(gI1,2)*ZP(gI2,1)) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(-ZP(gI2,2) + ZP(gI2
      ,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*KroneckerDelta(2,gO2)*(-2*Conj(Mu)*Lambdax*ZP(gI1,0)*ZP(
      gI2,0) - 4*Conj(MT)*Lambdax*ZP(gI1,1)*ZP(gI2,0) - 2*TLambdax*ZP(gI1,1)*ZP(
      gI2,0) + 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) -
      1.4142135623730951*vd*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) - 2*Conj(TLambdax)*ZP(gI1,
      0)*ZP(gI2,1) - 2*Conj(Mu)*Lambdax*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*
      vu*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,3)*ZP(
      gI2,1) + 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) + 4*vT*Sqr(g2)*ZP(gI1,2)*
      ZP(gI2,2) - 4*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) - 1.4142135623730951*vd*Sqr(g2)
      *ZP(gI1,0)*ZP(gI2,3) - 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) - 4
      *vT*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) + 4*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3) + Conj(
      Lambdax)*(1.4142135623730951*vd*Lambdax*ZP(gI1,3)*ZP(gI2,0) + 2*vT*Lambdax*
      ZP(gI1,1)*ZP(gI2,1) - 2*Mu*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*vu*
      Lambdax*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*Lambdax*ZP(gI1,2)*(vd*ZP(
      gI2,0) + vu*ZP(gI2,1)) - 1.4142135623730951*vu*Lambdax*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*vu*Lambdax*ZP(gI1,1)*ZP(gI2,3) + ZP(gI1,0)*(2*(vT*
      Lambdax - Mu)*ZP(gI2,0) - 4*MT*ZP(gI2,1) + 1.4142135623730951*vd*Lambdax*(
      -ZP(gI2,2) + ZP(gI2,3))))) + KroneckerDelta(1,gO2)*(5*(5.656854249492381*
      Conj(MT)*Lambdax*ZP(gI1,2)*ZP(gI2,0) + 2.8284271247461903*TLambdax*ZP(gI1,3)
      *ZP(gI2,0) + 1.4142135623730951*vT*AbsSqr(Lambdax)*ZP(gI1,2)*ZP(gI2,1) +
      2.8284271247461903*Conj(Lambdax)*Mu*ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951
      *vT*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951*vT*AbsSqr(Lambdax)*ZP(
      gI1,3)*ZP(gI2,1) + 2.8284271247461903*Conj(Mu)*Lambdax*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) - 2*vu*Sqr(g2)*ZP(gI1,2)*
      ZP(gI2,2) - 4*vu*AbsSqr(Lambdax)*ZP(gI1,3)*ZP(gI2,3) + 2*vu*Sqr(g2)*ZP(gI1,3
      )*ZP(gI2,3)) - ZP(gI1,1)*(5*vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,0) + vu*(3
      *Sqr(g1) + 5*Sqr(g2))*ZP(gI2,1) + 7.0710678118654755*((-(vT*AbsSqr(Lambdax))
      - 2*Conj(Mu)*Lambdax + vT*Sqr(g2))*ZP(gI2,2) + (vT*AbsSqr(Lambdax) - 2*Conj
      (Lambdax)*Mu - vT*Sqr(g2))*ZP(gI2,3))) + ZP(gI1,0)*(vu*(-20*AbsSqr(Lambdax)
      + 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI2,0) - 5*(vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(
      gI2,1) - 2.8284271247461903*(2*MT*Conj(Lambdax)*ZP(gI2,2) + Conj(TLambdax)*
      ZP(gI2,3))))) - KroneckerDelta(0,gO2)*(ZP(gI1,1)*(5*vu*(AbsSqr(Lambdax) +
      Sqr(g2))*ZP(gI2,0) + vd*(20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,
      1) + 14.142135623730951*(TLambdax*ZP(gI2,2) + 2*Conj(MT)*Lambdax*ZP(gI2,3)))
      + 5*(ZP(gI1,2)*(1.4142135623730951*(Conj(Lambdax)*(-(vT*Lambdax) + 2*Mu) +
      vT*Sqr(g2))*ZP(gI2,0) + 2.8284271247461903*Conj(TLambdax)*ZP(gI2,1) - 2*vd*(
      -2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,2)) + ZP(gI1,3)*(1.4142135623730951*(vT
      *AbsSqr(Lambdax) + 2*Conj(Mu)*Lambdax - vT*Sqr(g2))*ZP(gI2,0) +
      5.656854249492381*MT*Conj(Lambdax)*ZP(gI2,1) + 2*vd*Sqr(g2)*ZP(gI2,3))) + ZP
      (gI1,0)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,0) + 5*(vu*(AbsSqr(Lambdax) + Sqr
      (g2))*ZP(gI2,1) + 1.4142135623730951*((-(vT*AbsSqr(Lambdax)) + 2*Conj(Mu)*
      Lambdax + vT*Sqr(g2))*ZP(gI2,2) + (vT*AbsSqr(Lambdax) + 2*Conj(Lambdax)*Mu -
      vT*Sqr(g2))*ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(5*Conj(Lambdax)*KroneckerDelta(2,gO2)*(ZN(gI1,3)*ZN(gI2,2) +
      ZN(gI1,2)*ZN(gI2,3)) + KroneckerDelta(0,gO2)*(ZN(gI1,2)*(3.872983346207417*
      g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + 3.872983346207417*g1*ZN(gI1,0)*ZN(gI2,2) -
      5*g2*ZN(gI1,1)*ZN(gI2,2) + 5*Conj(Lambdax)*ZN(gI1,4)*ZN(gI2,3) + 5*Conj(
      Lambdax)*ZN(gI1,3)*ZN(gI2,4)) + KroneckerDelta(1,gO2)*(ZN(gI1,3)*(
      -3.872983346207417*g1*ZN(gI2,0) + 5*g2*ZN(gI2,1)) + (-3.872983346207417*g1*
      ZN(gI1,0) + 5*g2*ZN(gI1,1))*ZN(gI2,3) + 5*Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2,2)
      + ZN(gI1,2)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) -
      3.872983346207417*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))*KroneckerDelta(1,gO1) +
      5*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta(1,gO1) + 5*g2*Conj(ZN(
      gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN
      (gI1,0))*(Conj(ZN(gI2,2))*KroneckerDelta(0,gO1) - Conj(ZN(gI2,3))*
      KroneckerDelta(1,gO1)) + 5*Conj(ZN(gI1,4))*Conj(ZN(gI2,3))*KroneckerDelta(0,
      gO1)*Lambdax + 5*Conj(ZN(gI1,3))*Conj(ZN(gI2,4))*KroneckerDelta(0,gO1)*
      Lambdax + 5*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1)*Lambdax +
      5*Conj(ZN(gI1,3))*Conj(ZN(gI2,2))*KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(
      gI1,2))*(3.872983346207417*g1*Conj(ZN(gI2,0))*KroneckerDelta(0,gO1) - 5*g2*
      Conj(ZN(gI2,1))*KroneckerDelta(0,gO1) + 5*Conj(ZN(gI2,4))*KroneckerDelta(1,
      gO1)*Lambdax + 5*Conj(ZN(gI2,3))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1050;
   std::complex<double> tmp_1051;
   std::complex<double> tmp_1052;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1052 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1051 += tmp_1052;
   tmp_1050 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1051;
   std::complex<double> tmp_1053;
   std::complex<double> tmp_1054;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1054 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1053 += tmp_1054;
   tmp_1050 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1053;
   std::complex<double> tmp_1055;
   std::complex<double> tmp_1056;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1056 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1055 += tmp_1056;
   tmp_1050 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1055;
   std::complex<double> tmp_1057;
   std::complex<double> tmp_1058;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1058 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1057 += tmp_1058;
   tmp_1050 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1057;
   std::complex<double> tmp_1059;
   std::complex<double> tmp_1060;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1060 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1059 += tmp_1060;
   tmp_1050 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1059;
   std::complex<double> tmp_1061;
   std::complex<double> tmp_1062;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1062 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1061 += tmp_1062;
   tmp_1050 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1061;
   std::complex<double> tmp_1063;
   std::complex<double> tmp_1064;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1065;
      std::complex<double> tmp_1066;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1066 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1065 += tmp_1066;
      tmp_1064 += (Conj(ZD(gI2,j2))) * tmp_1065;
   }
   tmp_1063 += tmp_1064;
   tmp_1050 += (std::complex<double>(0,0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)) * tmp_1063;
   std::complex<double> tmp_1067;
   std::complex<double> tmp_1068;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1069;
      std::complex<double> tmp_1070;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1070 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1069 += tmp_1070;
      tmp_1068 += (Conj(ZD(gI2,j2))) * tmp_1069;
   }
   tmp_1067 += tmp_1068;
   tmp_1050 += (std::complex<double>(0,0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2)) * tmp_1067;
   std::complex<double> tmp_1071;
   std::complex<double> tmp_1072;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1073;
      std::complex<double> tmp_1074;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1074 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1073 += tmp_1074;
      tmp_1072 += (ZD(gI1,j2)) * tmp_1073;
   }
   tmp_1071 += tmp_1072;
   tmp_1050 += (std::complex<double>(0,0.35355339059327373)*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*Lambdax) * tmp_1071;
   std::complex<double> tmp_1075;
   std::complex<double> tmp_1076;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1077;
      std::complex<double> tmp_1078;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1078 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1077 += tmp_1078;
      tmp_1076 += (ZD(gI1,j2)) * tmp_1077;
   }
   tmp_1075 += tmp_1076;
   tmp_1050 += (std::complex<double>(0,0.35355339059327373)*KroneckerDelta(1,
      gO1)*KroneckerDelta(2,gO2)*Lambdax) * tmp_1075;
   std::complex<double> tmp_1079;
   std::complex<double> tmp_1080;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1081;
      std::complex<double> tmp_1082;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1083;
         std::complex<double> tmp_1084;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1084 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1083 += tmp_1084;
         tmp_1082 += (ZD(gI1,3 + j2)) * tmp_1083;
      }
      tmp_1081 += tmp_1082;
      tmp_1080 += (Conj(ZD(gI2,3 + j3))) * tmp_1081;
   }
   tmp_1079 += tmp_1080;
   tmp_1050 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1079;
   std::complex<double> tmp_1085;
   std::complex<double> tmp_1086;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1087;
      std::complex<double> tmp_1088;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1089;
         std::complex<double> tmp_1090;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1090 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1089 += tmp_1090;
         tmp_1088 += (Conj(ZD(gI2,j2))) * tmp_1089;
      }
      tmp_1087 += tmp_1088;
      tmp_1086 += (ZD(gI1,j3)) * tmp_1087;
   }
   tmp_1085 += tmp_1086;
   tmp_1050 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1085;
   result += (std::complex<double>(0,-1)) * tmp_1050;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1091;
   std::complex<double> tmp_1092;
   std::complex<double> tmp_1093;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1093 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1092 += tmp_1093;
   tmp_1091 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1092;
   std::complex<double> tmp_1094;
   std::complex<double> tmp_1095;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1095 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1094 += tmp_1095;
   tmp_1091 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1094;
   std::complex<double> tmp_1096;
   std::complex<double> tmp_1097;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1097 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1096 += tmp_1097;
   tmp_1091 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1096;
   std::complex<double> tmp_1098;
   std::complex<double> tmp_1099;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1099 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1098 += tmp_1099;
   tmp_1091 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1098;
   std::complex<double> tmp_1100;
   std::complex<double> tmp_1101;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1101 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1100 += tmp_1101;
   tmp_1091 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1100;
   std::complex<double> tmp_1102;
   std::complex<double> tmp_1103;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1103 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1102 += tmp_1103;
   tmp_1091 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1102;
   std::complex<double> tmp_1104;
   std::complex<double> tmp_1105;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1106;
      std::complex<double> tmp_1107;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1107 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1106 += tmp_1107;
      tmp_1105 += (Conj(ZE(gI2,j2))) * tmp_1106;
   }
   tmp_1104 += tmp_1105;
   tmp_1091 += (std::complex<double>(0,0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)) * tmp_1104;
   std::complex<double> tmp_1108;
   std::complex<double> tmp_1109;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1110;
      std::complex<double> tmp_1111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1111 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1110 += tmp_1111;
      tmp_1109 += (Conj(ZE(gI2,j2))) * tmp_1110;
   }
   tmp_1108 += tmp_1109;
   tmp_1091 += (std::complex<double>(0,0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2)) * tmp_1108;
   std::complex<double> tmp_1112;
   std::complex<double> tmp_1113;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1114;
      std::complex<double> tmp_1115;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1115 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1114 += tmp_1115;
      tmp_1113 += (ZE(gI1,j2)) * tmp_1114;
   }
   tmp_1112 += tmp_1113;
   tmp_1091 += (std::complex<double>(0,0.35355339059327373)*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*Lambdax) * tmp_1112;
   std::complex<double> tmp_1116;
   std::complex<double> tmp_1117;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1118;
      std::complex<double> tmp_1119;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1119 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1118 += tmp_1119;
      tmp_1117 += (ZE(gI1,j2)) * tmp_1118;
   }
   tmp_1116 += tmp_1117;
   tmp_1091 += (std::complex<double>(0,0.35355339059327373)*KroneckerDelta(1,
      gO1)*KroneckerDelta(2,gO2)*Lambdax) * tmp_1116;
   std::complex<double> tmp_1120;
   std::complex<double> tmp_1121;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1122;
      std::complex<double> tmp_1123;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1124;
         std::complex<double> tmp_1125;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1125 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1124 += tmp_1125;
         tmp_1123 += (ZE(gI1,3 + j2)) * tmp_1124;
      }
      tmp_1122 += tmp_1123;
      tmp_1121 += (Conj(ZE(gI2,3 + j3))) * tmp_1122;
   }
   tmp_1120 += tmp_1121;
   tmp_1091 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1120;
   std::complex<double> tmp_1126;
   std::complex<double> tmp_1127;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1128;
      std::complex<double> tmp_1129;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1130;
         std::complex<double> tmp_1131;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1131 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1130 += tmp_1131;
         tmp_1129 += (Conj(ZE(gI2,j2))) * tmp_1130;
      }
      tmp_1128 += tmp_1129;
      tmp_1127 += (ZE(gI1,j3)) * tmp_1128;
   }
   tmp_1126 += tmp_1127;
   tmp_1091 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1126;
   result += (std::complex<double>(0,-1)) * tmp_1091;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1132;
   std::complex<double> tmp_1133;
   std::complex<double> tmp_1134;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1134 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1133 += tmp_1134;
   tmp_1132 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1133;
   std::complex<double> tmp_1135;
   std::complex<double> tmp_1136;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1136 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1135 += tmp_1136;
   tmp_1132 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1135;
   std::complex<double> tmp_1137;
   std::complex<double> tmp_1138;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1138 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1137 += tmp_1138;
   tmp_1132 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1137;
   std::complex<double> tmp_1139;
   std::complex<double> tmp_1140;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1140 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1139 += tmp_1140;
   tmp_1132 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1139;
   std::complex<double> tmp_1141;
   std::complex<double> tmp_1142;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1142 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1141 += tmp_1142;
   tmp_1132 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1141;
   std::complex<double> tmp_1143;
   std::complex<double> tmp_1144;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1144 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1143 += tmp_1144;
   tmp_1132 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1143;
   std::complex<double> tmp_1145;
   std::complex<double> tmp_1146;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1147;
      std::complex<double> tmp_1148;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1148 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1147 += tmp_1148;
      tmp_1146 += (Conj(ZU(gI2,j2))) * tmp_1147;
   }
   tmp_1145 += tmp_1146;
   tmp_1132 += (std::complex<double>(0,0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(0,gO2)*KroneckerDelta(2,gO1)) * tmp_1145;
   std::complex<double> tmp_1149;
   std::complex<double> tmp_1150;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1151;
      std::complex<double> tmp_1152;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1152 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1151 += tmp_1152;
      tmp_1150 += (Conj(ZU(gI2,j2))) * tmp_1151;
   }
   tmp_1149 += tmp_1150;
   tmp_1132 += (std::complex<double>(0,0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2)) * tmp_1149;
   std::complex<double> tmp_1153;
   std::complex<double> tmp_1154;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1155;
      std::complex<double> tmp_1156;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1156 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1155 += tmp_1156;
      tmp_1154 += (ZU(gI1,j2)) * tmp_1155;
   }
   tmp_1153 += tmp_1154;
   tmp_1132 += (std::complex<double>(0,0.35355339059327373)*KroneckerDelta(0,
      gO2)*KroneckerDelta(2,gO1)*Lambdax) * tmp_1153;
   std::complex<double> tmp_1157;
   std::complex<double> tmp_1158;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1159;
      std::complex<double> tmp_1160;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1160 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1159 += tmp_1160;
      tmp_1158 += (ZU(gI1,j2)) * tmp_1159;
   }
   tmp_1157 += tmp_1158;
   tmp_1132 += (std::complex<double>(0,0.35355339059327373)*KroneckerDelta(0,
      gO1)*KroneckerDelta(2,gO2)*Lambdax) * tmp_1157;
   std::complex<double> tmp_1161;
   std::complex<double> tmp_1162;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1163;
      std::complex<double> tmp_1164;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1165;
         std::complex<double> tmp_1166;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1166 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1165 += tmp_1166;
         tmp_1164 += (ZU(gI1,3 + j2)) * tmp_1165;
      }
      tmp_1163 += tmp_1164;
      tmp_1162 += (Conj(ZU(gI2,3 + j3))) * tmp_1163;
   }
   tmp_1161 += tmp_1162;
   tmp_1132 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1161;
   std::complex<double> tmp_1167;
   std::complex<double> tmp_1168;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1169;
      std::complex<double> tmp_1170;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1171;
         std::complex<double> tmp_1172;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1172 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1171 += tmp_1172;
         tmp_1170 += (Conj(ZU(gI2,j2))) * tmp_1171;
      }
      tmp_1169 += tmp_1170;
      tmp_1168 += (ZU(gI1,j3)) * tmp_1169;
   }
   tmp_1167 += tmp_1168;
   tmp_1132 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1167;
   result += (std::complex<double>(0,-1)) * tmp_1132;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1173;
   std::complex<double> tmp_1174;
   std::complex<double> tmp_1175;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1175 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1174 += tmp_1175;
   tmp_1173 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1174;
   std::complex<double> tmp_1176;
   std::complex<double> tmp_1177;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1177 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1176 += tmp_1177;
   tmp_1173 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1176;
   std::complex<double> tmp_1178;
   std::complex<double> tmp_1179;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1179 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1178 += tmp_1179;
   tmp_1173 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1178;
   std::complex<double> tmp_1180;
   std::complex<double> tmp_1181;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1181 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1180 += tmp_1181;
   tmp_1173 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1180;
   std::complex<double> tmp_1182;
   std::complex<double> tmp_1183;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1183 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1182 += tmp_1183;
   tmp_1173 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1182;
   std::complex<double> tmp_1184;
   std::complex<double> tmp_1185;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1185 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1184 += tmp_1185;
   tmp_1173 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1184;
   std::complex<double> tmp_1186;
   std::complex<double> tmp_1187;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1188;
      std::complex<double> tmp_1189;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1189 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1188 += tmp_1189;
      tmp_1187 += (Conj(ZD(gI2,j2))) * tmp_1188;
   }
   tmp_1186 += tmp_1187;
   tmp_1173 += (std::complex<double>(0,0.35355339059327373)*vT*Conj(Lambdax)*
      KroneckerDelta(1,gO2)) * tmp_1186;
   std::complex<double> tmp_1190;
   std::complex<double> tmp_1191;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1192;
      std::complex<double> tmp_1193;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1193 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1192 += tmp_1193;
      tmp_1191 += (Conj(ZD(gI2,j2))) * tmp_1192;
   }
   tmp_1190 += tmp_1191;
   tmp_1173 += (std::complex<double>(0,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1190;
   std::complex<double> tmp_1194;
   std::complex<double> tmp_1195;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1196;
      std::complex<double> tmp_1197;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1197 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1196 += tmp_1197;
      tmp_1195 += (Conj(ZD(gI2,j2))) * tmp_1196;
   }
   tmp_1194 += tmp_1195;
   tmp_1173 += (std::complex<double>(0,0.35355339059327373)*vu*Conj(Lambdax)*
      KroneckerDelta(2,gO2)) * tmp_1194;
   std::complex<double> tmp_1198;
   std::complex<double> tmp_1199;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1200;
      std::complex<double> tmp_1201;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1201 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1200 += tmp_1201;
      tmp_1199 += (Conj(ZD(gI2,j2))) * tmp_1200;
   }
   tmp_1198 += tmp_1199;
   tmp_1173 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1198;
   std::complex<double> tmp_1202;
   std::complex<double> tmp_1203;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1204;
      std::complex<double> tmp_1205;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1205 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1204 += tmp_1205;
      tmp_1203 += (ZD(gI1,j2)) * tmp_1204;
   }
   tmp_1202 += tmp_1203;
   tmp_1173 += (std::complex<double>(0,0.35355339059327373)*vT*KroneckerDelta(1
      ,gO2)*Lambdax) * tmp_1202;
   std::complex<double> tmp_1206;
   std::complex<double> tmp_1207;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1208;
      std::complex<double> tmp_1209;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1209 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1208 += tmp_1209;
      tmp_1207 += (ZD(gI1,j2)) * tmp_1208;
   }
   tmp_1206 += tmp_1207;
   tmp_1173 += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2
      )*Mu) * tmp_1206;
   std::complex<double> tmp_1210;
   std::complex<double> tmp_1211;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1212;
      std::complex<double> tmp_1213;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1213 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1212 += tmp_1213;
      tmp_1211 += (ZD(gI1,j2)) * tmp_1212;
   }
   tmp_1210 += tmp_1211;
   tmp_1173 += (std::complex<double>(0,0.35355339059327373)*vu*KroneckerDelta(2
      ,gO2)*Lambdax) * tmp_1210;
   std::complex<double> tmp_1214;
   std::complex<double> tmp_1215;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1216;
      std::complex<double> tmp_1217;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1217 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1216 += tmp_1217;
      tmp_1215 += (ZD(gI1,j2)) * tmp_1216;
   }
   tmp_1214 += tmp_1215;
   tmp_1173 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1214;
   std::complex<double> tmp_1218;
   std::complex<double> tmp_1219;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1220;
      std::complex<double> tmp_1221;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1222;
         std::complex<double> tmp_1223;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1223 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1222 += tmp_1223;
         tmp_1221 += (ZD(gI1,3 + j2)) * tmp_1222;
      }
      tmp_1220 += tmp_1221;
      tmp_1219 += (Conj(ZD(gI2,3 + j3))) * tmp_1220;
   }
   tmp_1218 += tmp_1219;
   tmp_1173 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1218
      ;
   std::complex<double> tmp_1224;
   std::complex<double> tmp_1225;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1226;
      std::complex<double> tmp_1227;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1228;
         std::complex<double> tmp_1229;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1229 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1228 += tmp_1229;
         tmp_1227 += (Conj(ZD(gI2,j2))) * tmp_1228;
      }
      tmp_1226 += tmp_1227;
      tmp_1225 += (ZD(gI1,j3)) * tmp_1226;
   }
   tmp_1224 += tmp_1225;
   tmp_1173 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1224
      ;
   result += (std::complex<double>(0,-1)) * tmp_1173;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1230;
   std::complex<double> tmp_1231;
   std::complex<double> tmp_1232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1232 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1231 += tmp_1232;
   tmp_1230 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1231;
   std::complex<double> tmp_1233;
   std::complex<double> tmp_1234;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1234 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1233 += tmp_1234;
   tmp_1230 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1233;
   std::complex<double> tmp_1235;
   std::complex<double> tmp_1236;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1236 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1235 += tmp_1236;
   tmp_1230 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1235;
   std::complex<double> tmp_1237;
   std::complex<double> tmp_1238;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1238 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1237 += tmp_1238;
   tmp_1230 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1237;
   std::complex<double> tmp_1239;
   std::complex<double> tmp_1240;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1240 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1239 += tmp_1240;
   tmp_1230 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1239;
   std::complex<double> tmp_1241;
   std::complex<double> tmp_1242;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1242 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1241 += tmp_1242;
   tmp_1230 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1241;
   std::complex<double> tmp_1243;
   std::complex<double> tmp_1244;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1245;
      std::complex<double> tmp_1246;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1246 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1245 += tmp_1246;
      tmp_1244 += (Conj(ZE(gI2,j2))) * tmp_1245;
   }
   tmp_1243 += tmp_1244;
   tmp_1230 += (std::complex<double>(0,0.35355339059327373)*vT*Conj(Lambdax)*
      KroneckerDelta(1,gO2)) * tmp_1243;
   std::complex<double> tmp_1247;
   std::complex<double> tmp_1248;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1249;
      std::complex<double> tmp_1250;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1250 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1249 += tmp_1250;
      tmp_1248 += (Conj(ZE(gI2,j2))) * tmp_1249;
   }
   tmp_1247 += tmp_1248;
   tmp_1230 += (std::complex<double>(0,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1247;
   std::complex<double> tmp_1251;
   std::complex<double> tmp_1252;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1253;
      std::complex<double> tmp_1254;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1254 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1253 += tmp_1254;
      tmp_1252 += (Conj(ZE(gI2,j2))) * tmp_1253;
   }
   tmp_1251 += tmp_1252;
   tmp_1230 += (std::complex<double>(0,0.35355339059327373)*vu*Conj(Lambdax)*
      KroneckerDelta(2,gO2)) * tmp_1251;
   std::complex<double> tmp_1255;
   std::complex<double> tmp_1256;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1257;
      std::complex<double> tmp_1258;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1258 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1257 += tmp_1258;
      tmp_1256 += (Conj(ZE(gI2,j2))) * tmp_1257;
   }
   tmp_1255 += tmp_1256;
   tmp_1230 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1255;
   std::complex<double> tmp_1259;
   std::complex<double> tmp_1260;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1261;
      std::complex<double> tmp_1262;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1262 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1261 += tmp_1262;
      tmp_1260 += (ZE(gI1,j2)) * tmp_1261;
   }
   tmp_1259 += tmp_1260;
   tmp_1230 += (std::complex<double>(0,0.35355339059327373)*vT*KroneckerDelta(1
      ,gO2)*Lambdax) * tmp_1259;
   std::complex<double> tmp_1263;
   std::complex<double> tmp_1264;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1265;
      std::complex<double> tmp_1266;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1266 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1265 += tmp_1266;
      tmp_1264 += (ZE(gI1,j2)) * tmp_1265;
   }
   tmp_1263 += tmp_1264;
   tmp_1230 += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2
      )*Mu) * tmp_1263;
   std::complex<double> tmp_1267;
   std::complex<double> tmp_1268;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1269;
      std::complex<double> tmp_1270;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1270 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1269 += tmp_1270;
      tmp_1268 += (ZE(gI1,j2)) * tmp_1269;
   }
   tmp_1267 += tmp_1268;
   tmp_1230 += (std::complex<double>(0,0.35355339059327373)*vu*KroneckerDelta(2
      ,gO2)*Lambdax) * tmp_1267;
   std::complex<double> tmp_1271;
   std::complex<double> tmp_1272;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1273;
      std::complex<double> tmp_1274;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1274 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1273 += tmp_1274;
      tmp_1272 += (ZE(gI1,j2)) * tmp_1273;
   }
   tmp_1271 += tmp_1272;
   tmp_1230 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,
      gO2)) * tmp_1271;
   std::complex<double> tmp_1275;
   std::complex<double> tmp_1276;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1277;
      std::complex<double> tmp_1278;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1279;
         std::complex<double> tmp_1280;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1280 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1279 += tmp_1280;
         tmp_1278 += (ZE(gI1,3 + j2)) * tmp_1279;
      }
      tmp_1277 += tmp_1278;
      tmp_1276 += (Conj(ZE(gI2,3 + j3))) * tmp_1277;
   }
   tmp_1275 += tmp_1276;
   tmp_1230 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1275
      ;
   std::complex<double> tmp_1281;
   std::complex<double> tmp_1282;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1283;
      std::complex<double> tmp_1284;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1285;
         std::complex<double> tmp_1286;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1286 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1285 += tmp_1286;
         tmp_1284 += (Conj(ZE(gI2,j2))) * tmp_1285;
      }
      tmp_1283 += tmp_1284;
      tmp_1282 += (ZE(gI1,j3)) * tmp_1283;
   }
   tmp_1281 += tmp_1282;
   tmp_1230 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1281
      ;
   result += (std::complex<double>(0,-1)) * tmp_1230;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1287;
   std::complex<double> tmp_1288;
   std::complex<double> tmp_1289;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1289 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1288 += tmp_1289;
   tmp_1287 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1288;
   std::complex<double> tmp_1290;
   std::complex<double> tmp_1291;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1291 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1290 += tmp_1291;
   tmp_1287 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1290;
   std::complex<double> tmp_1292;
   std::complex<double> tmp_1293;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1293 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1292 += tmp_1293;
   tmp_1287 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1292;
   std::complex<double> tmp_1294;
   std::complex<double> tmp_1295;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1295 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1294 += tmp_1295;
   tmp_1287 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1294;
   std::complex<double> tmp_1296;
   std::complex<double> tmp_1297;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1297 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1296 += tmp_1297;
   tmp_1287 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1296;
   std::complex<double> tmp_1298;
   std::complex<double> tmp_1299;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1299 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1298 += tmp_1299;
   tmp_1287 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_1298;
   std::complex<double> tmp_1300;
   std::complex<double> tmp_1301;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1302;
      std::complex<double> tmp_1303;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1303 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1302 += tmp_1303;
      tmp_1301 += (Conj(ZU(gI2,j2))) * tmp_1302;
   }
   tmp_1300 += tmp_1301;
   tmp_1287 += (std::complex<double>(0,0.35355339059327373)*vT*Conj(Lambdax)*
      KroneckerDelta(0,gO2)) * tmp_1300;
   std::complex<double> tmp_1304;
   std::complex<double> tmp_1305;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1306;
      std::complex<double> tmp_1307;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1307 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1306 += tmp_1307;
      tmp_1305 += (Conj(ZU(gI2,j2))) * tmp_1306;
   }
   tmp_1304 += tmp_1305;
   tmp_1287 += (std::complex<double>(0,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(0,gO2)) * tmp_1304;
   std::complex<double> tmp_1308;
   std::complex<double> tmp_1309;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1310;
      std::complex<double> tmp_1311;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1311 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1310 += tmp_1311;
      tmp_1309 += (Conj(ZU(gI2,j2))) * tmp_1310;
   }
   tmp_1308 += tmp_1309;
   tmp_1287 += (std::complex<double>(0,0.35355339059327373)*vd*Conj(Lambdax)*
      KroneckerDelta(2,gO2)) * tmp_1308;
   std::complex<double> tmp_1312;
   std::complex<double> tmp_1313;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1314;
      std::complex<double> tmp_1315;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1315 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1314 += tmp_1315;
      tmp_1313 += (Conj(ZU(gI2,j2))) * tmp_1314;
   }
   tmp_1312 += tmp_1313;
   tmp_1287 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1312;
   std::complex<double> tmp_1316;
   std::complex<double> tmp_1317;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1318;
      std::complex<double> tmp_1319;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1319 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1318 += tmp_1319;
      tmp_1317 += (ZU(gI1,j2)) * tmp_1318;
   }
   tmp_1316 += tmp_1317;
   tmp_1287 += (std::complex<double>(0,0.35355339059327373)*vT*KroneckerDelta(0
      ,gO2)*Lambdax) * tmp_1316;
   std::complex<double> tmp_1320;
   std::complex<double> tmp_1321;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1322;
      std::complex<double> tmp_1323;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1323 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1322 += tmp_1323;
      tmp_1321 += (ZU(gI1,j2)) * tmp_1322;
   }
   tmp_1320 += tmp_1321;
   tmp_1287 += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2
      )*Mu) * tmp_1320;
   std::complex<double> tmp_1324;
   std::complex<double> tmp_1325;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1326;
      std::complex<double> tmp_1327;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1327 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1326 += tmp_1327;
      tmp_1325 += (ZU(gI1,j2)) * tmp_1326;
   }
   tmp_1324 += tmp_1325;
   tmp_1287 += (std::complex<double>(0,0.35355339059327373)*vd*KroneckerDelta(2
      ,gO2)*Lambdax) * tmp_1324;
   std::complex<double> tmp_1328;
   std::complex<double> tmp_1329;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1330;
      std::complex<double> tmp_1331;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1331 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1330 += tmp_1331;
      tmp_1329 += (ZU(gI1,j2)) * tmp_1330;
   }
   tmp_1328 += tmp_1329;
   tmp_1287 += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,
      gO2)) * tmp_1328;
   std::complex<double> tmp_1332;
   std::complex<double> tmp_1333;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1334;
      std::complex<double> tmp_1335;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1336;
         std::complex<double> tmp_1337;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1337 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1336 += tmp_1337;
         tmp_1335 += (ZU(gI1,3 + j2)) * tmp_1336;
      }
      tmp_1334 += tmp_1335;
      tmp_1333 += (Conj(ZU(gI2,3 + j3))) * tmp_1334;
   }
   tmp_1332 += tmp_1333;
   tmp_1287 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1332
      ;
   std::complex<double> tmp_1338;
   std::complex<double> tmp_1339;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1340;
      std::complex<double> tmp_1341;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1342;
         std::complex<double> tmp_1343;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1343 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1342 += tmp_1343;
         tmp_1341 += (Conj(ZU(gI2,j2))) * tmp_1342;
      }
      tmp_1340 += tmp_1341;
      tmp_1339 += (ZU(gI1,j3)) * tmp_1340;
   }
   tmp_1338 += tmp_1339;
   tmp_1287 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1338
      ;
   result += (std::complex<double>(0,-1)) * tmp_1287;

   return result;
}

std::complex<double> CLASSNAME::CpUhhVZAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZA(gI2,0) -
      KroneckerDelta(1,gO2)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjVWmHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) - KroneckerDelta(1,gO2)*ZP(
      gI2,1) + 1.4142135623730951*KroneckerDelta(2,gO2)*(ZP(gI2,2) + ZP(gI2,3)));

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
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2) + 4*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2))
      *Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,gO1)*(KroneckerDelta(2,
      gO2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1)) + KroneckerDelta(0,gO2)*(ZA
      (gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*
      ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2))) + KroneckerDelta(1,gO1)*(KroneckerDelta(0,
      gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,1)) - 10*AbsSqr(Lambdax)*KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA(
      gI2,1) + ZA(gI1,1)*ZA(gI2,2)) + KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax)
      + 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - 3*(3*Sqr(g1) + 5*Sqr(g2))*ZA(
      gI1,1)*ZA(gI2,1) - 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2))) - KroneckerDelta
      (0,gO1)*(-(KroneckerDelta(1,gO2)*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2
      ))*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1))) + 10*AbsSqr(Lambdax)*
      KroneckerDelta(2,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) +
      KroneckerDelta(0,gO2)*(3*(3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (10*
      AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 10*AbsSqr(
      Lambdax)*ZA(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2
      ));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-10*AbsSqr(Lambdax)*KroneckerDelta(2,gO1)*KroneckerDelta(2,
      gO2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,
      0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) - 10*AbsSqr(
      Lambdax)*ZH(gI1,2)*ZH(gI2,2)) - KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*
      ((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (10*AbsSqr(Lambdax) - 3*Sqr(
      g1) - 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 10*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2
      )));

   return result;
}

std::complex<double> CLASSNAME::CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*((Conj(TLambdax) + 2*Conj(MT)*Lambdax
      - TLambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1
      )) + KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) +
      KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2))) - 2*Conj(
      Mu)*Lambdax*(KroneckerDelta(2,gO2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1
      )) + KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,2)) +
      KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2))) + 2*Conj(
      Lambdax)*(KroneckerDelta(2,gO2)*(ZA(gI1,0)*(Mu*ZA(gI2,0) - MT*ZA(gI2,1)) +
      ZA(gI1,1)*(-(MT*ZA(gI2,0)) + Mu*ZA(gI2,1))) + KroneckerDelta(0,gO2)*(ZA(gI1,
      2)*(Mu*ZA(gI2,0) - MT*ZA(gI2,1)) + (Mu*ZA(gI1,0) - MT*ZA(gI1,1))*ZA(gI2,2))
      + KroneckerDelta(1,gO2)*(ZA(gI1,2)*(-(MT*ZA(gI2,0)) + Mu*ZA(gI2,1)) + (-(MT*
      ZA(gI1,0)) + Mu*ZA(gI1,1))*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(5*KroneckerDelta(2,gO2)*(2*Conj(MT)*Lambdax*(ZA(gI2,1)*ZH(gI1
      ,0) + ZA(gI2,0)*ZH(gI1,1)) - (Conj(TLambdax) + TLambdax)*(ZA(gI2,1)*ZH(gI1,0
      ) + ZA(gI2,0)*ZH(gI1,1)) + 2*Conj(Lambdax)*(MT*ZA(gI2,1)*ZH(gI1,0) + MT*ZA(
      gI2,0)*ZH(gI1,1) - Lambdax*ZA(gI2,2)*(vd*ZH(gI1,0) + vu*ZH(gI1,1)))) +
      KroneckerDelta(1,gO2)*(ZA(gI2,1)*(vd*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*
      Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1) - 10*(Conj(Mu)*
      Lambdax + Conj(Lambdax)*(vT*Lambdax + Mu))*ZH(gI1,2)) + 5*(2*MT*Conj(Lambdax
      )*(ZA(gI2,2)*ZH(gI1,0) - ZA(gI2,0)*ZH(gI1,2)) + 2*Conj(MT)*Lambdax*(ZA(gI2,2
      )*ZH(gI1,0) - ZA(gI2,0)*ZH(gI1,2)) - (Conj(TLambdax) + TLambdax)*(ZA(gI2,2)*
      ZH(gI1,0) + ZA(gI2,0)*ZH(gI1,2)))) - KroneckerDelta(0,gO2)*(ZA(gI2,0)*(vd*(3
      *Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) + vu*(10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr
      (g2))*ZH(gI1,1) + 10*(Conj(Mu)*Lambdax + Conj(Lambdax)*(vT*Lambdax + Mu))*ZH
      (gI1,2)) + 5*((Conj(TLambdax) + TLambdax)*(ZA(gI2,2)*ZH(gI1,1) + ZA(gI2,1)*
      ZH(gI1,2)) + Conj(Lambdax)*(-2*MT*ZA(gI2,2)*ZH(gI1,1) + 2*MT*ZA(gI2,1)*ZH(
      gI1,2)) + Conj(MT)*(-2*Lambdax*ZA(gI2,2)*ZH(gI1,1) + 2*Lambdax*ZA(gI2,1)*ZH(
      gI1,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(-(Conj(TLambdax)*KroneckerDelta(2,gO2
      )*ZH(gI1,1)*ZH(gI2,0)) - 2*Conj(MT)*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1)*
      ZH(gI2,0) + KroneckerDelta(2,gO2)*TLambdax*ZH(gI1,1)*ZH(gI2,0) - Conj(
      TLambdax)*KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,0) + 2*Conj(MT)*
      KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,2)*ZH(gI2,0) + KroneckerDelta(1,gO2)*
      TLambdax*ZH(gI1,2)*ZH(gI2,0) - Conj(TLambdax)*KroneckerDelta(2,gO2)*ZH(gI1,0
      )*ZH(gI2,1) - 2*Conj(MT)*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,0)*ZH(gI2,1) +
      KroneckerDelta(2,gO2)*TLambdax*ZH(gI1,0)*ZH(gI2,1) - Conj(TLambdax)*
      KroneckerDelta(0,gO2)*ZH(gI1,2)*ZH(gI2,1) + 2*Conj(MT)*KroneckerDelta(0,gO2)
      *Lambdax*ZH(gI1,2)*ZH(gI2,1) + KroneckerDelta(0,gO2)*TLambdax*ZH(gI1,2)*ZH(
      gI2,1) - 2*Conj(Mu)*KroneckerDelta(2,gO2)*Lambdax*(ZH(gI1,0)*ZH(gI2,0) + ZH(
      gI1,1)*ZH(gI2,1)) - Conj(TLambdax)*KroneckerDelta(1,gO2)*ZH(gI1,0)*ZH(gI2,2)
      + 2*Conj(MT)*KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,0)*ZH(gI2,2) +
      KroneckerDelta(1,gO2)*TLambdax*ZH(gI1,0)*ZH(gI2,2) - Conj(TLambdax)*
      KroneckerDelta(0,gO2)*ZH(gI1,1)*ZH(gI2,2) + 2*Conj(MT)*KroneckerDelta(0,gO2)
      *Lambdax*ZH(gI1,1)*ZH(gI2,2) + KroneckerDelta(0,gO2)*TLambdax*ZH(gI1,1)*ZH(
      gI2,2) + 2*Conj(Lambdax)*(KroneckerDelta(2,gO2)*(ZH(gI1,0)*(Mu*ZH(gI2,0) +
      MT*ZH(gI2,1)) + ZH(gI1,1)*(MT*ZH(gI2,0) + Mu*ZH(gI2,1))) - MT*(
      KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) +
      KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(1.4142135623730951*KroneckerDelta(0,
      gO2)*(g2*UM(gI1,1)*UP(gI2,0) - Conj(Lambdax)*UM(gI1,2)*UP(gI2,1)) +
      KroneckerDelta(2,gO2)*(2*g2*UM(gI1,2)*UP(gI2,0) + Conj(Lambdax)*UM(gI1,1)*UP
      (gI2,1) - 2*g2*UM(gI1,0)*UP(gI2,2)) + 1.4142135623730951*KroneckerDelta(1,
      gO2)*(g2*UM(gI1,0)*UP(gI2,1) + Conj(Lambdax)*UM(gI1,1)*UP(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChaChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(g2*Conj(UM(gI2,0))*(1.4142135623730951
      *Conj(UP(gI1,1))*KroneckerDelta(1,gO1) - 2*Conj(UP(gI1,2))*KroneckerDelta(2,
      gO1)) + Conj(UM(gI2,2))*(2*g2*Conj(UP(gI1,0))*KroneckerDelta(2,gO1) -
      1.4142135623730951*Conj(UP(gI1,1))*KroneckerDelta(0,gO1)*Lambdax) + Conj(UM(
      gI2,1))*(1.4142135623730951*g2*Conj(UP(gI1,0))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(UP(gI1,2))*KroneckerDelta(1,gO1)*Lambdax + Conj(UP(
      gI1,1))*KroneckerDelta(2,gO1)*Lambdax));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
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
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2))
      * tmp_1344;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
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
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_1348;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
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
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2))
      * tmp_1352;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
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
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_1356;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
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
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2))
      * tmp_1360;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
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
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,gO1)
      ) * tmp_1364;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*((3*Sqr(g1) +
      5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2)
      )*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2)*ZP(
      gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)))) + 5*(KroneckerDelta(1,gO2)*(AbsSqr(
      Lambdax) + Sqr(g2))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) -
      1.4142135623730951*KroneckerDelta(2,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZP(
      gI1,2)*ZP(gI2,0) + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2,3))))
      ) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*AbsSqr(Lambdax) + 3*
      Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)
      *ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - (-2*AbsSqr(Lambdax) + Sqr(g2)
      )*ZP(gI1,3)*ZP(gI2,3))) + 5*(KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + Sqr(g2
      ))*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) + 1.4142135623730951*
      KroneckerDelta(2,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZP(gI1,2)*ZP(gI2,1) + ZP
      (gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3))))) - 5*KroneckerDelta(
      2,gO1)*(Sqr(g2)*(4*KroneckerDelta(2,gO2)*(ZP(gI1,2) + ZP(gI1,3))*(ZP(gI2,2)
      + ZP(gI2,3)) + 1.4142135623730951*KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0)
      + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2,3))) -
      1.4142135623730951*KroneckerDelta(1,gO2)*(ZP(gI1,2)*ZP(gI2,1) + ZP(gI1,3)*ZP
      (gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3)))) + AbsSqr(Lambdax)*(2*
      KroneckerDelta(2,gO2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1)) +
      1.4142135623730951*(-(KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,3)
      *ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2,3)))) + KroneckerDelta(1,gO2)*(ZP
      (gI1,2)*ZP(gI2,1) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3)))
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*(4*Conj(MT)*KroneckerDelta(2,gO2)*
      Lambdax*ZP(gI1,1)*ZP(gI2,0) - vu*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(
      gI2,0) - vd*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,0) - 2*
      KroneckerDelta(2,gO2)*TLambdax*ZP(gI1,1)*ZP(gI2,0) - 5.656854249492381*Conj(
      MT)*KroneckerDelta(1,gO2)*Lambdax*ZP(gI1,2)*ZP(gI2,0) - 1.4142135623730951*
      vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 1.4142135623730951*vd
      *KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 1.4142135623730951*vT*
      KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) + 1.4142135623730951*vd*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) - 2.8284271247461903*
      KroneckerDelta(1,gO2)*TLambdax*ZP(gI1,3)*ZP(gI2,0) + 2*Conj(TLambdax)*
      KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,1) + vu*KroneckerDelta(0,gO2)*Sqr(g2)
      *ZP(gI1,0)*ZP(gI2,1) + vd*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,1)
      - 2.8284271247461903*Conj(TLambdax)*KroneckerDelta(0,gO2)*ZP(gI1,2)*ZP(gI2,1
      ) + 1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1)
      + 1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) -
      1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*vd*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) -
      1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) +
      2.8284271247461903*KroneckerDelta(0,gO2)*TLambdax*ZP(gI1,1)*ZP(gI2,2) + 4*vT
      *KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) - 2*Conj(Mu)*Lambdax*(
      KroneckerDelta(2,gO2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1)) +
      1.4142135623730951*(KroneckerDelta(0,gO2)*(ZP(gI1,3)*ZP(gI2,0) - ZP(gI1,0)*
      ZP(gI2,2)) + KroneckerDelta(1,gO2)*(ZP(gI1,3)*ZP(gI2,1) - ZP(gI1,1)*ZP(gI2,2
      )))) + 2.8284271247461903*Conj(TLambdax)*KroneckerDelta(1,gO2)*ZP(gI1,0)*ZP(
      gI2,3) - 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(
      gI2,3) - 1.4142135623730951*vd*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(
      gI2,3) + 5.656854249492381*Conj(MT)*KroneckerDelta(0,gO2)*Lambdax*ZP(gI1,1)*
      ZP(gI2,3) + 1.4142135623730951*vT*KroneckerDelta(1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP
      (gI2,3) - 1.4142135623730951*vu*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(
      gI2,3) - 4*vT*KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,3) + Conj(
      Lambdax)*(KroneckerDelta(0,gO2)*(-(vu*Lambdax*ZP(gI1,1)*ZP(gI2,0)) +
      1.4142135623730951*(vT*Lambdax - 2*Mu)*ZP(gI1,2)*ZP(gI2,0) -
      1.4142135623730951*vT*Lambdax*ZP(gI1,3)*ZP(gI2,0) + vu*Lambdax*ZP(gI1,0)*ZP(
      gI2,1) - 5.656854249492381*MT*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*vT*
      Lambdax*ZP(gI1,0)*ZP(gI2,2) + 1.4142135623730951*vT*Lambdax*ZP(gI1,0)*ZP(gI2
      ,3) + 2.8284271247461903*Mu*ZP(gI1,0)*ZP(gI2,3)) + KroneckerDelta(2,gO2)*(
      -1.4142135623730951*vd*Lambdax*ZP(gI1,3)*ZP(gI2,0) + 2*Mu*ZP(gI1,1)*ZP(gI2,1
      ) - 1.4142135623730951*vu*Lambdax*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*
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

std::complex<double> CLASSNAME::CpUAhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(-5*Conj(Lambdax)*KroneckerDelta(2,gO2)
      *(ZN(gI1,3)*ZN(gI2,2) + ZN(gI1,2)*ZN(gI2,3)) - KroneckerDelta(1,gO2)*(ZN(gI1
      ,3)*(3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + (3.872983346207417*
      g1*ZN(gI1,0) - 5*g2*ZN(gI1,1))*ZN(gI2,3) + 5*Conj(Lambdax)*(ZN(gI1,4)*ZN(gI2
      ,2) + ZN(gI1,2)*ZN(gI2,4))) + KroneckerDelta(0,gO2)*(ZN(gI1,2)*(
      3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(gI2,1)) + 3.872983346207417*g1*ZN(
      gI1,0)*ZN(gI2,2) - 5*(g2*ZN(gI1,1)*ZN(gI2,2) + Conj(Lambdax)*ZN(gI1,4)*ZN(
      gI2,3) + Conj(Lambdax)*ZN(gI1,3)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,2))*
      KroneckerDelta(0,gO1) + 3.872983346207417*g1*Conj(ZN(gI1,3))*Conj(ZN(gI2,0))
      *KroneckerDelta(1,gO1) - 5*g2*Conj(ZN(gI1,3))*Conj(ZN(gI2,1))*KroneckerDelta
      (1,gO1) - 5*g2*Conj(ZN(gI1,1))*Conj(ZN(gI2,3))*KroneckerDelta(1,gO1) +
      3.872983346207417*g1*Conj(ZN(gI1,0))*(-(Conj(ZN(gI2,2))*KroneckerDelta(0,gO1
      )) + Conj(ZN(gI2,3))*KroneckerDelta(1,gO1)) + 5*Conj(ZN(gI1,4))*Conj(ZN(gI2,
      3))*KroneckerDelta(0,gO1)*Lambdax + 5*Conj(ZN(gI1,3))*Conj(ZN(gI2,4))*
      KroneckerDelta(0,gO1)*Lambdax + 5*Conj(ZN(gI1,4))*Conj(ZN(gI2,2))*
      KroneckerDelta(1,gO1)*Lambdax + 5*Conj(ZN(gI1,3))*Conj(ZN(gI2,2))*
      KroneckerDelta(2,gO1)*Lambdax + Conj(ZN(gI1,2))*(-3.872983346207417*g1*Conj(
      ZN(gI2,0))*KroneckerDelta(0,gO1) + 5*(g2*Conj(ZN(gI2,1))*KroneckerDelta(0,
      gO1) + Conj(ZN(gI2,4))*KroneckerDelta(1,gO1)*Lambdax + Conj(ZN(gI2,3))*
      KroneckerDelta(2,gO1)*Lambdax)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
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
   tmp_1368 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1373;
   std::complex<double> tmp_1375;
   std::complex<double> tmp_1376;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1376 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1375 += tmp_1376;
   tmp_1368 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1375;
   std::complex<double> tmp_1377;
   std::complex<double> tmp_1378;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1378 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1377 += tmp_1378;
   tmp_1368 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1377;
   std::complex<double> tmp_1379;
   std::complex<double> tmp_1380;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1380 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1379 += tmp_1380;
   tmp_1368 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1379;
   std::complex<double> tmp_1381;
   std::complex<double> tmp_1382;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1383;
      std::complex<double> tmp_1384;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1384 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1383 += tmp_1384;
      tmp_1382 += (Conj(ZD(gI2,j2))) * tmp_1383;
   }
   tmp_1381 += tmp_1382;
   tmp_1368 += (std::complex<double>(0,-0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)) * tmp_1381;
   std::complex<double> tmp_1385;
   std::complex<double> tmp_1386;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1387;
      std::complex<double> tmp_1388;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1388 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1387 += tmp_1388;
      tmp_1386 += (Conj(ZD(gI2,j2))) * tmp_1387;
   }
   tmp_1385 += tmp_1386;
   tmp_1368 += (std::complex<double>(0,-0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2)) * tmp_1385;
   std::complex<double> tmp_1389;
   std::complex<double> tmp_1390;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1391;
      std::complex<double> tmp_1392;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1392 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1391 += tmp_1392;
      tmp_1390 += (ZD(gI1,j2)) * tmp_1391;
   }
   tmp_1389 += tmp_1390;
   tmp_1368 += (std::complex<double>(0,-0.35355339059327373)*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*Lambdax) * tmp_1389;
   std::complex<double> tmp_1393;
   std::complex<double> tmp_1394;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1395;
      std::complex<double> tmp_1396;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1396 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1395 += tmp_1396;
      tmp_1394 += (ZD(gI1,j2)) * tmp_1395;
   }
   tmp_1393 += tmp_1394;
   tmp_1368 += (std::complex<double>(0,-0.35355339059327373)*KroneckerDelta(1,
      gO1)*KroneckerDelta(2,gO2)*Lambdax) * tmp_1393;
   std::complex<double> tmp_1397;
   std::complex<double> tmp_1398;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1399;
      std::complex<double> tmp_1400;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1401;
         std::complex<double> tmp_1402;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1402 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1401 += tmp_1402;
         tmp_1400 += (ZD(gI1,3 + j2)) * tmp_1401;
      }
      tmp_1399 += tmp_1400;
      tmp_1398 += (Conj(ZD(gI2,3 + j3))) * tmp_1399;
   }
   tmp_1397 += tmp_1398;
   tmp_1368 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1397;
   std::complex<double> tmp_1403;
   std::complex<double> tmp_1404;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1405;
      std::complex<double> tmp_1406;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1407;
         std::complex<double> tmp_1408;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1408 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1407 += tmp_1408;
         tmp_1406 += (Conj(ZD(gI2,j2))) * tmp_1407;
      }
      tmp_1405 += tmp_1406;
      tmp_1404 += (ZD(gI1,j3)) * tmp_1405;
   }
   tmp_1403 += tmp_1404;
   tmp_1368 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1403;
   result += (std::complex<double>(0,-1)) * tmp_1368;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1409;
   std::complex<double> tmp_1410;
   std::complex<double> tmp_1411;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1411 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1410 += tmp_1411;
   tmp_1409 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1410;
   std::complex<double> tmp_1412;
   std::complex<double> tmp_1413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1413 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1412 += tmp_1413;
   tmp_1409 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1412;
   std::complex<double> tmp_1414;
   std::complex<double> tmp_1415;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1415 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1414 += tmp_1415;
   tmp_1409 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1414;
   std::complex<double> tmp_1416;
   std::complex<double> tmp_1417;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1417 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1416 += tmp_1417;
   tmp_1409 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1416;
   std::complex<double> tmp_1418;
   std::complex<double> tmp_1419;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1419 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1418 += tmp_1419;
   tmp_1409 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1418;
   std::complex<double> tmp_1420;
   std::complex<double> tmp_1421;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1421 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1420 += tmp_1421;
   tmp_1409 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1420;
   std::complex<double> tmp_1422;
   std::complex<double> tmp_1423;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1424;
      std::complex<double> tmp_1425;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1425 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1424 += tmp_1425;
      tmp_1423 += (Conj(ZE(gI2,j2))) * tmp_1424;
   }
   tmp_1422 += tmp_1423;
   tmp_1409 += (std::complex<double>(0,-0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)) * tmp_1422;
   std::complex<double> tmp_1426;
   std::complex<double> tmp_1427;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1428;
      std::complex<double> tmp_1429;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1429 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1428 += tmp_1429;
      tmp_1427 += (Conj(ZE(gI2,j2))) * tmp_1428;
   }
   tmp_1426 += tmp_1427;
   tmp_1409 += (std::complex<double>(0,-0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(1,gO1)*KroneckerDelta(2,gO2)) * tmp_1426;
   std::complex<double> tmp_1430;
   std::complex<double> tmp_1431;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1432;
      std::complex<double> tmp_1433;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1433 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1432 += tmp_1433;
      tmp_1431 += (ZE(gI1,j2)) * tmp_1432;
   }
   tmp_1430 += tmp_1431;
   tmp_1409 += (std::complex<double>(0,-0.35355339059327373)*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*Lambdax) * tmp_1430;
   std::complex<double> tmp_1434;
   std::complex<double> tmp_1435;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1436;
      std::complex<double> tmp_1437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1437 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1436 += tmp_1437;
      tmp_1435 += (ZE(gI1,j2)) * tmp_1436;
   }
   tmp_1434 += tmp_1435;
   tmp_1409 += (std::complex<double>(0,-0.35355339059327373)*KroneckerDelta(1,
      gO1)*KroneckerDelta(2,gO2)*Lambdax) * tmp_1434;
   std::complex<double> tmp_1438;
   std::complex<double> tmp_1439;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1440;
      std::complex<double> tmp_1441;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1442;
         std::complex<double> tmp_1443;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1443 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1442 += tmp_1443;
         tmp_1441 += (ZE(gI1,3 + j2)) * tmp_1442;
      }
      tmp_1440 += tmp_1441;
      tmp_1439 += (Conj(ZE(gI2,3 + j3))) * tmp_1440;
   }
   tmp_1438 += tmp_1439;
   tmp_1409 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1438;
   std::complex<double> tmp_1444;
   std::complex<double> tmp_1445;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1446;
      std::complex<double> tmp_1447;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1448;
         std::complex<double> tmp_1449;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1449 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1448 += tmp_1449;
         tmp_1447 += (Conj(ZE(gI2,j2))) * tmp_1448;
      }
      tmp_1446 += tmp_1447;
      tmp_1445 += (ZE(gI1,j3)) * tmp_1446;
   }
   tmp_1444 += tmp_1445;
   tmp_1409 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1444;
   result += (std::complex<double>(0,-1)) * tmp_1409;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1450;
   std::complex<double> tmp_1451;
   std::complex<double> tmp_1452;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1452 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1451 += tmp_1452;
   tmp_1450 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1451;
   std::complex<double> tmp_1453;
   std::complex<double> tmp_1454;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1454 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1453 += tmp_1454;
   tmp_1450 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1453;
   std::complex<double> tmp_1455;
   std::complex<double> tmp_1456;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1456 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1455 += tmp_1456;
   tmp_1450 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1455;
   std::complex<double> tmp_1457;
   std::complex<double> tmp_1458;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1458 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1457 += tmp_1458;
   tmp_1450 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1457;
   std::complex<double> tmp_1459;
   std::complex<double> tmp_1460;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1460 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1459 += tmp_1460;
   tmp_1450 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1459;
   std::complex<double> tmp_1461;
   std::complex<double> tmp_1462;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1462 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1461 += tmp_1462;
   tmp_1450 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1461;
   std::complex<double> tmp_1463;
   std::complex<double> tmp_1464;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1465;
      std::complex<double> tmp_1466;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1466 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1465 += tmp_1466;
      tmp_1464 += (Conj(ZU(gI2,j2))) * tmp_1465;
   }
   tmp_1463 += tmp_1464;
   tmp_1450 += (std::complex<double>(0,-0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(0,gO2)*KroneckerDelta(2,gO1)) * tmp_1463;
   std::complex<double> tmp_1467;
   std::complex<double> tmp_1468;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1469;
      std::complex<double> tmp_1470;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1470 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1469 += tmp_1470;
      tmp_1468 += (Conj(ZU(gI2,j2))) * tmp_1469;
   }
   tmp_1467 += tmp_1468;
   tmp_1450 += (std::complex<double>(0,-0.35355339059327373)*Conj(Lambdax)*
      KroneckerDelta(0,gO1)*KroneckerDelta(2,gO2)) * tmp_1467;
   std::complex<double> tmp_1471;
   std::complex<double> tmp_1472;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1473;
      std::complex<double> tmp_1474;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1474 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1473 += tmp_1474;
      tmp_1472 += (ZU(gI1,j2)) * tmp_1473;
   }
   tmp_1471 += tmp_1472;
   tmp_1450 += (std::complex<double>(0,-0.35355339059327373)*KroneckerDelta(0,
      gO2)*KroneckerDelta(2,gO1)*Lambdax) * tmp_1471;
   std::complex<double> tmp_1475;
   std::complex<double> tmp_1476;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1477;
      std::complex<double> tmp_1478;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1478 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1477 += tmp_1478;
      tmp_1476 += (ZU(gI1,j2)) * tmp_1477;
   }
   tmp_1475 += tmp_1476;
   tmp_1450 += (std::complex<double>(0,-0.35355339059327373)*KroneckerDelta(0,
      gO1)*KroneckerDelta(2,gO2)*Lambdax) * tmp_1475;
   std::complex<double> tmp_1479;
   std::complex<double> tmp_1480;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1481;
      std::complex<double> tmp_1482;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1483;
         std::complex<double> tmp_1484;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1484 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1483 += tmp_1484;
         tmp_1482 += (ZU(gI1,3 + j2)) * tmp_1483;
      }
      tmp_1481 += tmp_1482;
      tmp_1480 += (Conj(ZU(gI2,3 + j3))) * tmp_1481;
   }
   tmp_1479 += tmp_1480;
   tmp_1450 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1479;
   std::complex<double> tmp_1485;
   std::complex<double> tmp_1486;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1487;
      std::complex<double> tmp_1488;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1489;
         std::complex<double> tmp_1490;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1490 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1489 += tmp_1490;
         tmp_1488 += (Conj(ZU(gI2,j2))) * tmp_1489;
      }
      tmp_1487 += tmp_1488;
      tmp_1486 += (ZU(gI1,j3)) * tmp_1487;
   }
   tmp_1485 += tmp_1486;
   tmp_1450 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1485;
   result += (std::complex<double>(0,-1)) * tmp_1450;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1491;
   std::complex<double> tmp_1492;
   std::complex<double> tmp_1493;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1494;
      std::complex<double> tmp_1495;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1495 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1494 += tmp_1495;
      tmp_1493 += (Conj(ZD(gI2,j2))) * tmp_1494;
   }
   tmp_1492 += tmp_1493;
   tmp_1491 += (0.35355339059327373*vT*Conj(Lambdax)*KroneckerDelta(1,gO2)) *
      tmp_1492;
   std::complex<double> tmp_1496;
   std::complex<double> tmp_1497;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1498;
      std::complex<double> tmp_1499;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1499 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1498 += tmp_1499;
      tmp_1497 += (Conj(ZD(gI2,j2))) * tmp_1498;
   }
   tmp_1496 += tmp_1497;
   tmp_1491 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1496;
   std::complex<double> tmp_1500;
   std::complex<double> tmp_1501;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1502;
      std::complex<double> tmp_1503;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1503 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1502 += tmp_1503;
      tmp_1501 += (Conj(ZD(gI2,j2))) * tmp_1502;
   }
   tmp_1500 += tmp_1501;
   tmp_1491 += (0.35355339059327373*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) *
      tmp_1500;
   std::complex<double> tmp_1504;
   std::complex<double> tmp_1505;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1506;
      std::complex<double> tmp_1507;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1507 += ZD(gI1,3 + j1)*TYd(j1,j2);
      }
      tmp_1506 += tmp_1507;
      tmp_1505 += (Conj(ZD(gI2,j2))) * tmp_1506;
   }
   tmp_1504 += tmp_1505;
   tmp_1491 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1504;
   std::complex<double> tmp_1508;
   std::complex<double> tmp_1509;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1510;
      std::complex<double> tmp_1511;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1511 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1510 += tmp_1511;
      tmp_1509 += (ZD(gI1,j2)) * tmp_1510;
   }
   tmp_1508 += tmp_1509;
   tmp_1491 += (-0.35355339059327373*vT*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_1508;
   std::complex<double> tmp_1512;
   std::complex<double> tmp_1513;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1514;
      std::complex<double> tmp_1515;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1515 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1514 += tmp_1515;
      tmp_1513 += (ZD(gI1,j2)) * tmp_1514;
   }
   tmp_1512 += tmp_1513;
   tmp_1491 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1512;
   std::complex<double> tmp_1516;
   std::complex<double> tmp_1517;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1518;
      std::complex<double> tmp_1519;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1519 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1518 += tmp_1519;
      tmp_1517 += (ZD(gI1,j2)) * tmp_1518;
   }
   tmp_1516 += tmp_1517;
   tmp_1491 += (-0.35355339059327373*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1516;
   std::complex<double> tmp_1520;
   std::complex<double> tmp_1521;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1522;
      std::complex<double> tmp_1523;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1523 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1522 += tmp_1523;
      tmp_1521 += (ZD(gI1,j2)) * tmp_1522;
   }
   tmp_1520 += tmp_1521;
   tmp_1491 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1520;
   result += (std::complex<double>(0,-1)) * tmp_1491;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1524;
   std::complex<double> tmp_1525;
   std::complex<double> tmp_1526;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1527;
      std::complex<double> tmp_1528;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1528 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1527 += tmp_1528;
      tmp_1526 += (Conj(ZE(gI2,j2))) * tmp_1527;
   }
   tmp_1525 += tmp_1526;
   tmp_1524 += (0.35355339059327373*vT*Conj(Lambdax)*KroneckerDelta(1,gO2)) *
      tmp_1525;
   std::complex<double> tmp_1529;
   std::complex<double> tmp_1530;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1531;
      std::complex<double> tmp_1532;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1532 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1531 += tmp_1532;
      tmp_1530 += (Conj(ZE(gI2,j2))) * tmp_1531;
   }
   tmp_1529 += tmp_1530;
   tmp_1524 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1529;
   std::complex<double> tmp_1533;
   std::complex<double> tmp_1534;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1535;
      std::complex<double> tmp_1536;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1536 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1535 += tmp_1536;
      tmp_1534 += (Conj(ZE(gI2,j2))) * tmp_1535;
   }
   tmp_1533 += tmp_1534;
   tmp_1524 += (0.35355339059327373*vu*Conj(Lambdax)*KroneckerDelta(2,gO2)) *
      tmp_1533;
   std::complex<double> tmp_1537;
   std::complex<double> tmp_1538;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1539;
      std::complex<double> tmp_1540;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1540 += ZE(gI1,3 + j1)*TYe(j1,j2);
      }
      tmp_1539 += tmp_1540;
      tmp_1538 += (Conj(ZE(gI2,j2))) * tmp_1539;
   }
   tmp_1537 += tmp_1538;
   tmp_1524 += (0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1537;
   std::complex<double> tmp_1541;
   std::complex<double> tmp_1542;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1543;
      std::complex<double> tmp_1544;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1544 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1543 += tmp_1544;
      tmp_1542 += (ZE(gI1,j2)) * tmp_1543;
   }
   tmp_1541 += tmp_1542;
   tmp_1524 += (-0.35355339059327373*vT*KroneckerDelta(1,gO2)*Lambdax) *
      tmp_1541;
   std::complex<double> tmp_1545;
   std::complex<double> tmp_1546;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1547;
      std::complex<double> tmp_1548;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1548 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1547 += tmp_1548;
      tmp_1546 += (ZE(gI1,j2)) * tmp_1547;
   }
   tmp_1545 += tmp_1546;
   tmp_1524 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1545;
   std::complex<double> tmp_1549;
   std::complex<double> tmp_1550;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1551;
      std::complex<double> tmp_1552;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1552 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1551 += tmp_1552;
      tmp_1550 += (ZE(gI1,j2)) * tmp_1551;
   }
   tmp_1549 += tmp_1550;
   tmp_1524 += (-0.35355339059327373*vu*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1549;
   std::complex<double> tmp_1553;
   std::complex<double> tmp_1554;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1555;
      std::complex<double> tmp_1556;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1556 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1555 += tmp_1556;
      tmp_1554 += (ZE(gI1,j2)) * tmp_1555;
   }
   tmp_1553 += tmp_1554;
   tmp_1524 += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_1553;
   result += (std::complex<double>(0,-1)) * tmp_1524;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1557;
   std::complex<double> tmp_1558;
   std::complex<double> tmp_1559;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1560;
      std::complex<double> tmp_1561;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1561 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1560 += tmp_1561;
      tmp_1559 += (Conj(ZU(gI2,j2))) * tmp_1560;
   }
   tmp_1558 += tmp_1559;
   tmp_1557 += (0.35355339059327373*vT*Conj(Lambdax)*KroneckerDelta(0,gO2)) *
      tmp_1558;
   std::complex<double> tmp_1562;
   std::complex<double> tmp_1563;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1564;
      std::complex<double> tmp_1565;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1565 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1564 += tmp_1565;
      tmp_1563 += (Conj(ZU(gI2,j2))) * tmp_1564;
   }
   tmp_1562 += tmp_1563;
   tmp_1557 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(0,gO2)) * tmp_1562;
   std::complex<double> tmp_1566;
   std::complex<double> tmp_1567;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1568;
      std::complex<double> tmp_1569;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1569 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1568 += tmp_1569;
      tmp_1567 += (Conj(ZU(gI2,j2))) * tmp_1568;
   }
   tmp_1566 += tmp_1567;
   tmp_1557 += (0.35355339059327373*vd*Conj(Lambdax)*KroneckerDelta(2,gO2)) *
      tmp_1566;
   std::complex<double> tmp_1570;
   std::complex<double> tmp_1571;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1572;
      std::complex<double> tmp_1573;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1573 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1572 += tmp_1573;
      tmp_1571 += (Conj(ZU(gI2,j2))) * tmp_1572;
   }
   tmp_1570 += tmp_1571;
   tmp_1557 += (0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1570;
   std::complex<double> tmp_1574;
   std::complex<double> tmp_1575;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1576;
      std::complex<double> tmp_1577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1577 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1576 += tmp_1577;
      tmp_1575 += (ZU(gI1,j2)) * tmp_1576;
   }
   tmp_1574 += tmp_1575;
   tmp_1557 += (-0.35355339059327373*vT*KroneckerDelta(0,gO2)*Lambdax) *
      tmp_1574;
   std::complex<double> tmp_1578;
   std::complex<double> tmp_1579;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1580;
      std::complex<double> tmp_1581;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1581 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1580 += tmp_1581;
      tmp_1579 += (ZU(gI1,j2)) * tmp_1580;
   }
   tmp_1578 += tmp_1579;
   tmp_1557 += (-0.7071067811865475*KroneckerDelta(0,gO2)*Mu) * tmp_1578;
   std::complex<double> tmp_1582;
   std::complex<double> tmp_1583;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1584;
      std::complex<double> tmp_1585;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1585 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1584 += tmp_1585;
      tmp_1583 += (ZU(gI1,j2)) * tmp_1584;
   }
   tmp_1582 += tmp_1583;
   tmp_1557 += (-0.35355339059327373*vd*KroneckerDelta(2,gO2)*Lambdax) *
      tmp_1582;
   std::complex<double> tmp_1586;
   std::complex<double> tmp_1587;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1588;
      std::complex<double> tmp_1589;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1589 += Conj(ZU(gI2,3 + j1))*Conj(TYu(j1,j2));
      }
      tmp_1588 += tmp_1589;
      tmp_1587 += (ZU(gI1,j2)) * tmp_1588;
   }
   tmp_1586 += tmp_1587;
   tmp_1557 += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_1586;
   result += (std::complex<double>(0,-1)) * tmp_1557;

   return result;
}

std::complex<double> CLASSNAME::CpUAhVZhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()))*(KroneckerDelta(0,gO2)*ZH(gI2,0) -
      KroneckerDelta(1,gO2)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjVWmHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(2,gO2)*(
      ZP(gI2,2) - ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmVP(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.1*g2*(3.872983346207417*g1*vd*Cos(ThetaW())*KroneckerDelta(0,gO2
      ) - 3.872983346207417*g1*vu*Cos(ThetaW())*KroneckerDelta(1,gO2) +
      7.0710678118654755*g2*vT*(KroneckerDelta(2,gO2) + KroneckerDelta(3,gO2))*Sin
      (ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZVWm(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.1*g2*(7.0710678118654755*g2*vT*Cos(ThetaW())*KroneckerDelta(2,
      gO2) + 7.0710678118654755*g2*vT*Cos(ThetaW())*KroneckerDelta(3,gO2) +
      3.872983346207417*g1*(-(vd*KroneckerDelta(0,gO2)) + vu*KroneckerDelta(1,gO2)
      )*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargWmCgZ(unsigned gO1) const
{
   std::complex<double> result;

   result = 0.05*g2*(14.142135623730951*g2*vT*Cos(ThetaW())*(KroneckerDelta(2,
      gO1) + KroneckerDelta(3,gO1)) + vd*KroneckerDelta(0,gO1)*(5*g2*Cos(ThetaW())
      - 3.872983346207417*g1*Sin(ThetaW())) + KroneckerDelta(1,gO1)*(-5*g2*vu*Cos
      (ThetaW()) + 3.872983346207417*g1*vu*Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgWmCbargZ(unsigned gO2) const
{
   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO2) - vu*KroneckerDelta(1,gO2))*(5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbargZgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.05*g2*(vd*KroneckerDelta(0,gO1) - vu*KroneckerDelta(1,gO1))*(5*
      g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmgZbargWm(unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*g2*(14.142135623730951*g2*vT*Cos(ThetaW())*(KroneckerDelta(2,
      gO2) + KroneckerDelta(3,gO2)) + vd*KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW())
      - 3.872983346207417*g1*Sin(ThetaW())) + KroneckerDelta(1,gO2)*(-5*g2*vu*Cos
      (ThetaW()) + 3.872983346207417*g1*vu*Sin(ThetaW())));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.1*(20*(KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2) +
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2)*Sqr(Cos(ThetaW())) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-7.745966692414834*g1*g2*Cos(
      ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(
      ThetaW()))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(
      -7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW()))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2) + 2*KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)
      + 2*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(-(KroneckerDelta(0,gO2)*((3*Sqr(g1) +
      5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2)
      )*ZA(gI1,1)*ZA(gI2,1) + 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2))) + 5*(
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) + ZA(
      gI1,0)*ZA(gI2,1)) - 1.4142135623730951*(KroneckerDelta(2,gO2) +
      KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,2)))) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20*
      AbsSqr(Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) +
      5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) - 10*AbsSqr(Lambdax)*ZA(gI1,2)*ZA(gI2,2)) + 5
      *(KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,1)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,1)) + 1.4142135623730951*(KroneckerDelta(2,gO2) +
      KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,1) +
      ZA(gI1,1)*ZA(gI2,2)))) + 5*(-(KroneckerDelta(3,gO1)*(1.4142135623730951*
      AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,1) +
      1.4142135623730951*AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2)
      - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,2) + 4*
      KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,2) + 1.4142135623730951*
      KroneckerDelta(0,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,2)) + 2*KroneckerDelta(3,gO2)*(Sqr(g2)*ZA(gI1,0)*ZA(gI2,0) -
      (-2*AbsSqr(Lambdax) + Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 2*Sqr(g2)*ZA(gI1,2)*ZA(
      gI2,2)))) + KroneckerDelta(2,gO1)*(-1.4142135623730951*AbsSqr(Lambdax)*
      KroneckerDelta(1,gO2)*ZA(gI1,2)*ZA(gI2,1) + 1.4142135623730951*
      KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,1) - 1.4142135623730951*
      AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2) +
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,2) - 4*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,2)*ZA(gI2,2) - 1.4142135623730951*
      KroneckerDelta(0,gO2)*(-AbsSqr(Lambdax) + Sqr(g2))*(ZA(gI1,2)*ZA(gI2,0) + ZA
      (gI1,0)*ZA(gI2,2)) + 2*KroneckerDelta(2,gO2)*((-2*AbsSqr(Lambdax) + Sqr(g2))
      *ZA(gI1,0)*ZA(gI2,0) - Sqr(g2)*(ZA(gI1,1)*ZA(gI2,1) + 2*ZA(gI1,2)*ZA(gI2,2))
      ))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1590;
   tmp_1590 += std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1590 += std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1590 += std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1590 += std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1590 += std::complex<double>(0,0.5)*KroneckerDelta(2,gO1)*KroneckerDelta
      (2,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1590 += std::complex<double>(0,-0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   std::complex<double> tmp_1591;
   std::complex<double> tmp_1592;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1593;
      std::complex<double> tmp_1594;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1595;
         std::complex<double> tmp_1596;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1596 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1595 += tmp_1596;
         tmp_1594 += (Conj(ZV(gI2,j2))) * tmp_1595;
      }
      tmp_1593 += tmp_1594;
      tmp_1592 += (ZV(gI1,j3)) * tmp_1593;
   }
   tmp_1591 += tmp_1592;
   tmp_1590 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1591;
   result += (std::complex<double>(0,-1)) * tmp_1590;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO1)*(KroneckerDelta(0,gO2)*((3*Sqr(g1) +
      5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2)
      )*ZH(gI1,1)*ZH(gI2,1) + 10*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2)) + 5*(
      KroneckerDelta(1,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) + ZH(
      gI1,0)*ZH(gI2,1)) + 1.4142135623730951*(KroneckerDelta(2,gO2) -
      KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,2)*ZH(gI2,0) +
      ZH(gI1,0)*ZH(gI2,2))))) + KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((-20
      *AbsSqr(Lambdax) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) +
      5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) - 10*AbsSqr(Lambdax)*ZH(gI1,2)*ZH(gI2,2)) -
      5*(KroneckerDelta(0,gO2)*(AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,1)*ZH(gI2,0) +
      ZH(gI1,0)*ZH(gI2,1)) + 1.4142135623730951*(KroneckerDelta(2,gO2) -
      KroneckerDelta(3,gO2))*(-AbsSqr(Lambdax) + Sqr(g2))*(ZH(gI1,2)*ZH(gI2,1) +
      ZH(gI1,1)*ZH(gI2,2)))) + 5*(KroneckerDelta(3,gO1)*(-1.4142135623730951*
      AbsSqr(Lambdax)*KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,1) +
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,2)*ZH(gI2,1) -
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

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1597;
   std::complex<double> tmp_1598;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1599;
      std::complex<double> tmp_1600;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1600 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1599 += tmp_1600;
      tmp_1598 += (ZUL(gI1,j2)) * tmp_1599;
   }
   tmp_1597 += tmp_1598;
   result += (KroneckerDelta(0,gO2)) * tmp_1597;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1601;
   std::complex<double> tmp_1602;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1603;
      std::complex<double> tmp_1604;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1604 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1603 += tmp_1604;
      tmp_1602 += (Conj(ZDL(gI2,j2))) * tmp_1603;
   }
   tmp_1601 += tmp_1602;
   result += (KroneckerDelta(1,gO1)) * tmp_1601;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1605;
   std::complex<double> tmp_1606;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1606 += Conj(Ye(j1,gI1))*ZER(gI2,j1);
   }
   tmp_1605 += tmp_1606;
   result += (KroneckerDelta(0,gO2)) * tmp_1605;

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

   std::complex<double> tmp_1607;
   std::complex<double> tmp_1608;
   std::complex<double> tmp_1609;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1609 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1608 += tmp_1609;
   tmp_1607 += (std::complex<double>(0,-0.35355339059327373)*vd*KroneckerDelta(
      0,gO2)*Sqr(g2)) * tmp_1608;
   std::complex<double> tmp_1610;
   std::complex<double> tmp_1611;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1611 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1610 += tmp_1611;
   tmp_1607 += (std::complex<double>(0,-0.35355339059327373)*vu*KroneckerDelta(
      1,gO2)*Sqr(g2)) * tmp_1610;
   std::complex<double> tmp_1612;
   std::complex<double> tmp_1613;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1613 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1612 += tmp_1613;
   tmp_1607 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(2,gO2)*Sqr(g2))
      * tmp_1612;
   std::complex<double> tmp_1614;
   std::complex<double> tmp_1615;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1615 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1614 += tmp_1615;
   tmp_1607 += (std::complex<double>(0,0.5)*vT*KroneckerDelta(3,gO2)*Sqr(g2)) *
      tmp_1614;
   std::complex<double> tmp_1616;
   std::complex<double> tmp_1617;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1618;
      std::complex<double> tmp_1619;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1619 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1618 += tmp_1619;
      tmp_1617 += (ZV(gI1,j2)) * tmp_1618;
   }
   tmp_1616 += tmp_1617;
   tmp_1607 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(1,gO2)*Lambdax)
      * tmp_1616;
   std::complex<double> tmp_1620;
   std::complex<double> tmp_1621;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1622;
      std::complex<double> tmp_1623;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1623 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1622 += tmp_1623;
      tmp_1621 += (ZV(gI1,j2)) * tmp_1622;
   }
   tmp_1620 += tmp_1621;
   tmp_1607 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1620;
   std::complex<double> tmp_1624;
   std::complex<double> tmp_1625;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1626;
      std::complex<double> tmp_1627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1627 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1626 += tmp_1627;
      tmp_1625 += (ZV(gI1,j2)) * tmp_1626;
   }
   tmp_1624 += tmp_1625;
   tmp_1607 += (std::complex<double>(0,-0.7071067811865475)*vu*KroneckerDelta(3
      ,gO2)*Lambdax) * tmp_1624;
   std::complex<double> tmp_1628;
   std::complex<double> tmp_1629;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1630;
      std::complex<double> tmp_1631;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1631 += Conj(ZE(gI2,3 + j1))*Conj(TYe(j1,j2));
      }
      tmp_1630 += tmp_1631;
      tmp_1629 += (ZV(gI1,j2)) * tmp_1630;
   }
   tmp_1628 += tmp_1629;
   tmp_1607 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1628;
   std::complex<double> tmp_1632;
   std::complex<double> tmp_1633;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1634;
      std::complex<double> tmp_1635;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1636;
         std::complex<double> tmp_1637;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1637 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1636 += tmp_1637;
         tmp_1635 += (Conj(ZE(gI2,j2))) * tmp_1636;
      }
      tmp_1634 += tmp_1635;
      tmp_1633 += (ZV(gI1,j3)) * tmp_1634;
   }
   tmp_1632 += tmp_1633;
   tmp_1607 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(0,
      gO2)) * tmp_1632;
   result += (std::complex<double>(0,-1)) * tmp_1607;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(ZP(gI1,1)*(KroneckerDelta(0,gO2)*(-10*
      AbsSqr(Lambdax) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,0) + 10*(KroneckerDelta(2,
      gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI2,2) - KroneckerDelta(3,gO2)*Sqr(g2
      )*ZP(gI2,3))) + KroneckerDelta(1,gO2)*((-10*AbsSqr(Lambdax) + 3*Sqr(g1) + 5*
      Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - 2*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1
      ) - 5*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2) + 5*Sqr(g2)*ZP(gI1,
      3)*ZP(gI2,3)))) - KroneckerDelta(0,gO1)*(ZP(gI1,0)*(KroneckerDelta(1,gO2)*(
      10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI2,1) + 10*(KroneckerDelta(2
      ,gO2)*Sqr(g2)*ZP(gI2,2) - KroneckerDelta(3,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2
      ))*ZP(gI2,3))) + KroneckerDelta(0,gO2)*(2*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*
      ZP(gI2,0) + (10*AbsSqr(Lambdax) - 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1)
      + 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - (-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,3
      )*ZP(gI2,3)))) - 10*(KroneckerDelta(2,gO1)*(ZP(gI1,2)*(KroneckerDelta(0,gO2)
      *Sqr(g2)*ZP(gI2,0) - KroneckerDelta(1,gO2)*(-2*AbsSqr(Lambdax) + Sqr(g2))*ZP
      (gI2,1) - 2*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI2,3)) + KroneckerDelta(2,gO2)
      *(Sqr(g2)*ZP(gI1,0)*ZP(gI2,0) - (-2*AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1)*ZP(
      gI2,1) + 2*Sqr(g2)*(2*ZP(gI1,2)*ZP(gI2,2) - ZP(gI1,3)*ZP(gI2,3)))) +
      KroneckerDelta(3,gO1)*(ZP(gI1,3)*(-(KroneckerDelta(0,gO2)*(-2*AbsSqr(Lambdax
      ) + Sqr(g2))*ZP(gI2,0)) + Sqr(g2)*(KroneckerDelta(1,gO2)*ZP(gI2,1) - 2*
      KroneckerDelta(2,gO2)*ZP(gI2,2))) + KroneckerDelta(3,gO2)*(-((-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gI1,0)*ZP(gI2,0)) + Sqr(g2)*(ZP(gI1,1)*ZP(gI2,1) - 2*
      ZP(gI1,2)*ZP(gI2,2) + 4*ZP(gI1,3)*ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*(1.4142135623730951*vT*AbsSqr(Lambdax)
      *KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) + 2.8284271247461903*Conj(Mu)*
      KroneckerDelta(3,gO2)*Lambdax*ZA(gI2,0)*ZP(gI1,0) - 1.4142135623730951*vT*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,0) + 2.8284271247461903*
      KroneckerDelta(3,gO2)*TLambdax*ZA(gI2,1)*ZP(gI1,0) + 1.4142135623730951*vd*
      AbsSqr(Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZP(gI1,0) + 2*Conj(Mu)*
      KroneckerDelta(0,gO2)*Lambdax*ZA(gI2,2)*ZP(gI1,0) - 2*Conj(Lambdax)*
      KroneckerDelta(0,gO2)*Mu*ZA(gI2,2)*ZP(gI1,0) - 1.4142135623730951*vd*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,0) - vu*AbsSqr(Lambdax)*
      KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(gI1,1) + 5.656854249492381*MT*Conj(
      Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,1) - vu*KroneckerDelta(0,gO2
      )*Sqr(g2)*ZA(gI2,0)*ZP(gI1,1) - vd*AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA(
      gI2,1)*ZP(gI1,1) - 1.4142135623730951*vT*AbsSqr(Lambdax)*KroneckerDelta(3,
      gO2)*ZA(gI2,1)*ZP(gI1,1) + 2.8284271247461903*Conj(Mu)*KroneckerDelta(3,gO2)
      *Lambdax*ZA(gI2,1)*ZP(gI1,1) - vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,1)*ZP
      (gI1,1) + 1.4142135623730951*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,1)*ZP(
      gI1,1) + 4*MT*Conj(Lambdax)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,1) - 2*
      Conj(TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,1) +
      1.4142135623730951*vu*AbsSqr(Lambdax)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZP(gI1
      ,1) - 1.4142135623730951*vu*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,1
      ) + 1.4142135623730951*vT*AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP
      (gI1,2) - 2.8284271247461903*Conj(Mu)*KroneckerDelta(0,gO2)*Lambdax*ZA(gI2,0
      )*ZP(gI1,2) - 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*
      ZP(gI1,2) - 5.656854249492381*MT*Conj(Lambdax)*KroneckerDelta(0,gO2)*ZA(gI2,
      1)*ZP(gI1,2) - 1.4142135623730951*vd*AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*
      ZA(gI2,2)*ZP(gI1,2) + 1.4142135623730951*vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA
      (gI2,2)*ZP(gI1,2) - 4*vT*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,2)*ZP(gI1,2) -
      1.4142135623730951*vT*AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZP(
      gI1,3) - 2.8284271247461903*Conj(Lambdax)*KroneckerDelta(0,gO2)*Mu*ZA(gI2,0)
      *ZP(gI1,3) + 1.4142135623730951*vT*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI2,0)*
      ZP(gI1,3) - 2.8284271247461903*Conj(TLambdax)*KroneckerDelta(0,gO2)*ZA(gI2,1
      )*ZP(gI1,3) - 1.4142135623730951*vd*AbsSqr(Lambdax)*KroneckerDelta(0,gO2)*ZA
      (gI2,2)*ZP(gI1,3) + 1.4142135623730951*vd*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(
      gI2,2)*ZP(gI1,3) + KroneckerDelta(2,gO2)*(5.656854249492381*Conj(MT)*Lambdax
      *ZA(gI2,1)*ZP(gI1,0) + 1.4142135623730951*vd*AbsSqr(Lambdax)*ZA(gI2,2)*ZP(
      gI1,0) - 1.4142135623730951*vd*Sqr(g2)*ZA(gI2,2)*ZP(gI1,0) +
      1.4142135623730951*vT*AbsSqr(Lambdax)*ZA(gI2,1)*ZP(gI1,1) +
      2.8284271247461903*Conj(Lambdax)*Mu*ZA(gI2,1)*ZP(gI1,1) - 1.4142135623730951
      *vT*Sqr(g2)*ZA(gI2,1)*ZP(gI1,1) + 1.4142135623730951*vu*AbsSqr(Lambdax)*ZA(
      gI2,2)*ZP(gI1,1) - 1.4142135623730951*vu*Sqr(g2)*ZA(gI2,2)*ZP(gI1,1) +
      1.4142135623730951*ZA(gI2,0)*((-(vT*AbsSqr(Lambdax)) + 2*Conj(Lambdax)*Mu +
      vT*Sqr(g2))*ZP(gI1,0) + 2*Conj(TLambdax)*ZP(gI1,1)) + 4*vT*Sqr(g2)*ZA(gI2,2)
      *ZP(gI1,3)) + KroneckerDelta(1,gO2)*(ZA(gI2,2)*(-4*Conj(MT)*Lambdax*ZP(gI1,0
      ) + 2*TLambdax*ZP(gI1,0) + 2*Conj(Mu)*Lambdax*ZP(gI1,1) - 2*Conj(Lambdax)*Mu
      *ZP(gI1,1) - 1.4142135623730951*vu*AbsSqr(Lambdax)*ZP(gI1,2) +
      1.4142135623730951*vu*Sqr(g2)*ZP(gI1,2) - 1.4142135623730951*vu*AbsSqr(
      Lambdax)*ZP(gI1,3) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI1,3)) + ZA(gI2,0)*(
      vu*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,0) - 2.8284271247461903*(TLambdax*ZP(
      gI1,2) + 2*Conj(MT)*Lambdax*ZP(gI1,3))) + ZA(gI2,1)*(vd*(AbsSqr(Lambdax) +
      Sqr(g2))*ZP(gI1,0) + 1.4142135623730951*((-(vT*AbsSqr(Lambdax)) - 2*Conj(Mu)
      *Lambdax + vT*Sqr(g2))*ZP(gI1,2) + (vT*AbsSqr(Lambdax) - 2*Conj(Lambdax)*Mu
      - vT*Sqr(g2))*ZP(gI1,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(KroneckerDelta(2,gO2)*(-5.656854249492381*Conj(MT)*
      Lambdax*ZH(gI2,1)*ZP(gI1,0) - 1.4142135623730951*vd*AbsSqr(Lambdax)*ZH(gI2,2
      )*ZP(gI1,0) + 1.4142135623730951*vd*Sqr(g2)*ZH(gI2,2)*ZP(gI1,0) -
      1.4142135623730951*vT*AbsSqr(Lambdax)*ZH(gI2,1)*ZP(gI1,1) -
      2.8284271247461903*Conj(Lambdax)*Mu*ZH(gI2,1)*ZP(gI1,1) + 1.4142135623730951
      *vT*Sqr(g2)*ZH(gI2,1)*ZP(gI1,1) - 1.4142135623730951*vu*AbsSqr(Lambdax)*ZH(
      gI2,2)*ZP(gI1,1) + 1.4142135623730951*vu*Sqr(g2)*ZH(gI2,2)*ZP(gI1,1) + 2*vu*
      Sqr(g2)*ZH(gI2,1)*ZP(gI1,2) + 4*vT*Sqr(g2)*ZH(gI2,2)*ZP(gI1,2) + ZH(gI2,0)*(
      1.4142135623730951*(Conj(Lambdax)*(-(vT*Lambdax) + 2*Mu) + vT*Sqr(g2))*ZP(
      gI1,0) + 2.8284271247461903*Conj(TLambdax)*ZP(gI1,1) - 2*vd*(-2*AbsSqr(
      Lambdax) + Sqr(g2))*ZP(gI1,2)) - 4*vT*Sqr(g2)*ZH(gI2,2)*ZP(gI1,3)) +
      KroneckerDelta(3,gO2)*(-2.8284271247461903*TLambdax*ZH(gI2,1)*ZP(gI1,0) +
      1.4142135623730951*vd*AbsSqr(Lambdax)*ZH(gI2,2)*ZP(gI1,0) -
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
      gI1,3)))) - KroneckerDelta(0,gO2)*(5*ZH(gI2,2)*(-2*Conj(Mu)*Lambdax*ZP(gI1,0
      ) - 2*Conj(TLambdax)*ZP(gI1,1) + 1.4142135623730951*vd*Sqr(g2)*ZP(gI1,2) -
      1.4142135623730951*vd*Sqr(g2)*ZP(gI1,3) + Conj(Lambdax)*(2*(vT*Lambdax - Mu)
      *ZP(gI1,0) - 4*MT*ZP(gI1,1) + 1.4142135623730951*vd*Lambdax*(-ZP(gI1,2) + ZP
      (gI1,3)))) + ZH(gI2,1)*(vu*(20*AbsSqr(Lambdax) - 3*Sqr(g1) + 5*Sqr(g2))*ZP(
      gI1,0) + 5*(vd*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1) - 2.8284271247461903*(2
      *MT*Conj(Lambdax)*ZP(gI1,2) + Conj(TLambdax)*ZP(gI1,3)))) + ZH(gI2,0)*(vd*(3
      *Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0) + 5*(vu*(AbsSqr(Lambdax) + Sqr(g2))*ZP(gI1,1
      ) + 1.4142135623730951*((-(vT*AbsSqr(Lambdax)) + 2*Conj(Mu)*Lambdax + vT*Sqr
      (g2))*ZP(gI1,2) + (vT*AbsSqr(Lambdax) + 2*Conj(Lambdax)*Mu - vT*Sqr(g2))*ZP(
      gI1,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(Lambdax)*KroneckerDelta(2,gO2)*UP(gI2,1)*ZN(gI1,2)) + Conj(
      Lambdax)*KroneckerDelta(0,gO2)*UP(gI2,2)*ZN(gI1,3) - 0.1*KroneckerDelta(1,
      gO2)*(UP(gI2,1)*(5.477225575051661*g1*ZN(gI1,0) + 7.0710678118654755*g2*ZN(
      gI1,1)) + 10*g2*UP(gI2,0)*ZN(gI1,3)) + 0.7071067811865475*Conj(Lambdax)*
      KroneckerDelta(0,gO2)*UP(gI2,1)*ZN(gI1,4) + 1.4142135623730951*g2*
      KroneckerDelta(3,gO2)*(-(UP(gI2,2)*ZN(gI1,1)) + UP(gI2,0)*ZN(gI1,4));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmChiChaPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UM(gI2,0))*(Conj(ZN(gI1,2))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN(gI1,4))*KroneckerDelta(2,gO1))) + Conj(UM(gI2,2))
      *(1.4142135623730951*g2*Conj(ZN(gI1,1))*KroneckerDelta(2,gO1) - Conj(ZN(gI1,
      2))*KroneckerDelta(1,gO1)*Lambdax) + Conj(UM(gI2,1))*(0.5477225575051661*g1*
      Conj(ZN(gI1,0))*KroneckerDelta(0,gO1) + 0.7071067811865475*g2*Conj(ZN(gI1,1)
      )*KroneckerDelta(0,gO1) + 0.7071067811865475*Conj(ZN(gI1,4))*KroneckerDelta(
      1,gO1)*Lambdax + Conj(ZN(gI1,3))*KroneckerDelta(3,gO1)*Lambdax);

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1638;
   std::complex<double> tmp_1639;
   std::complex<double> tmp_1640;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1640 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1639 += tmp_1640;
   tmp_1638 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1639;
   std::complex<double> tmp_1641;
   std::complex<double> tmp_1642;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1642 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1641 += tmp_1642;
   tmp_1638 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1641;
   std::complex<double> tmp_1643;
   std::complex<double> tmp_1644;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1644 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1643 += tmp_1644;
   tmp_1638 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1643;
   std::complex<double> tmp_1645;
   std::complex<double> tmp_1646;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1646 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1645 += tmp_1646;
   tmp_1638 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1645;
   std::complex<double> tmp_1647;
   std::complex<double> tmp_1648;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1648 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1647 += tmp_1648;
   tmp_1638 += (std::complex<double>(0,-0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1647;
   std::complex<double> tmp_1649;
   std::complex<double> tmp_1650;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1650 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1649 += tmp_1650;
   tmp_1638 += (std::complex<double>(0,0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1649;
   std::complex<double> tmp_1651;
   std::complex<double> tmp_1652;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1652 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1651 += tmp_1652;
   tmp_1638 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1651;
   std::complex<double> tmp_1653;
   std::complex<double> tmp_1654;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1654 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1653 += tmp_1654;
   tmp_1638 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1653;
   std::complex<double> tmp_1655;
   std::complex<double> tmp_1656;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1657;
      std::complex<double> tmp_1658;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1658 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1657 += tmp_1658;
      tmp_1656 += (Conj(ZD(gI2,j2))) * tmp_1657;
   }
   tmp_1655 += tmp_1656;
   tmp_1638 += (std::complex<double>(0,-1)*Conj(Lambdax)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)) * tmp_1655;
   std::complex<double> tmp_1659;
   std::complex<double> tmp_1660;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1661;
      std::complex<double> tmp_1662;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1662 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1661 += tmp_1662;
      tmp_1660 += (ZD(gI1,j2)) * tmp_1661;
   }
   tmp_1659 += tmp_1660;
   tmp_1638 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO2)*KroneckerDelta
      (2,gO1)*Lambdax) * tmp_1659;
   std::complex<double> tmp_1663;
   std::complex<double> tmp_1664;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1665;
      std::complex<double> tmp_1666;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1667;
         std::complex<double> tmp_1668;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1668 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1667 += tmp_1668;
         tmp_1666 += (ZD(gI1,3 + j2)) * tmp_1667;
      }
      tmp_1665 += tmp_1666;
      tmp_1664 += (Conj(ZD(gI2,3 + j3))) * tmp_1665;
   }
   tmp_1663 += tmp_1664;
   tmp_1638 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1663;
   std::complex<double> tmp_1669;
   std::complex<double> tmp_1670;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1671;
      std::complex<double> tmp_1672;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1673;
         std::complex<double> tmp_1674;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1674 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1673 += tmp_1674;
         tmp_1672 += (Conj(ZD(gI2,j2))) * tmp_1673;
      }
      tmp_1671 += tmp_1672;
      tmp_1670 += (ZD(gI1,j3)) * tmp_1671;
   }
   tmp_1669 += tmp_1670;
   tmp_1638 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1669;
   result += (std::complex<double>(0,-1)) * tmp_1638;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1675;
   std::complex<double> tmp_1676;
   std::complex<double> tmp_1677;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1677 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1676 += tmp_1677;
   tmp_1675 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1676;
   std::complex<double> tmp_1678;
   std::complex<double> tmp_1679;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1679 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1678 += tmp_1679;
   tmp_1675 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1678;
   std::complex<double> tmp_1680;
   std::complex<double> tmp_1681;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1681 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1680 += tmp_1681;
   tmp_1675 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1680;
   std::complex<double> tmp_1682;
   std::complex<double> tmp_1683;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1683 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1682 += tmp_1683;
   tmp_1675 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1682;
   std::complex<double> tmp_1684;
   std::complex<double> tmp_1685;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1685 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1684 += tmp_1685;
   tmp_1675 += (std::complex<double>(0,-0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1684;
   std::complex<double> tmp_1686;
   std::complex<double> tmp_1687;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1687 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1686 += tmp_1687;
   tmp_1675 += (std::complex<double>(0,0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1686;
   std::complex<double> tmp_1688;
   std::complex<double> tmp_1689;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1689 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1688 += tmp_1689;
   tmp_1675 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1688;
   std::complex<double> tmp_1690;
   std::complex<double> tmp_1691;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1691 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1690 += tmp_1691;
   tmp_1675 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1690;
   std::complex<double> tmp_1692;
   std::complex<double> tmp_1693;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1694;
      std::complex<double> tmp_1695;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1695 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1694 += tmp_1695;
      tmp_1693 += (Conj(ZE(gI2,j2))) * tmp_1694;
   }
   tmp_1692 += tmp_1693;
   tmp_1675 += (std::complex<double>(0,-1)*Conj(Lambdax)*KroneckerDelta(1,gO1)*
      KroneckerDelta(2,gO2)) * tmp_1692;
   std::complex<double> tmp_1696;
   std::complex<double> tmp_1697;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1698;
      std::complex<double> tmp_1699;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1699 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1698 += tmp_1699;
      tmp_1697 += (ZE(gI1,j2)) * tmp_1698;
   }
   tmp_1696 += tmp_1697;
   tmp_1675 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO2)*KroneckerDelta
      (2,gO1)*Lambdax) * tmp_1696;
   std::complex<double> tmp_1700;
   std::complex<double> tmp_1701;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1702;
      std::complex<double> tmp_1703;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1704;
         std::complex<double> tmp_1705;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1705 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1704 += tmp_1705;
         tmp_1703 += (ZE(gI1,3 + j2)) * tmp_1704;
      }
      tmp_1702 += tmp_1703;
      tmp_1701 += (Conj(ZE(gI2,3 + j3))) * tmp_1702;
   }
   tmp_1700 += tmp_1701;
   tmp_1675 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1700;
   result += (std::complex<double>(0,-1)) * tmp_1675;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1706;
   std::complex<double> tmp_1707;
   std::complex<double> tmp_1708;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1708 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1707 += tmp_1708;
   tmp_1706 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1707;
   std::complex<double> tmp_1709;
   std::complex<double> tmp_1710;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1710 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1709 += tmp_1710;
   tmp_1706 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1709;
   std::complex<double> tmp_1711;
   std::complex<double> tmp_1712;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1712 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1711 += tmp_1712;
   tmp_1706 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1711;
   std::complex<double> tmp_1713;
   std::complex<double> tmp_1714;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1714 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1713 += tmp_1714;
   tmp_1706 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1713;
   std::complex<double> tmp_1715;
   std::complex<double> tmp_1716;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1716 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1715 += tmp_1716;
   tmp_1706 += (std::complex<double>(0,0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1715;
   std::complex<double> tmp_1717;
   std::complex<double> tmp_1718;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1718 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1717 += tmp_1718;
   tmp_1706 += (std::complex<double>(0,-0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1717;
   std::complex<double> tmp_1719;
   std::complex<double> tmp_1720;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1720 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1719 += tmp_1720;
   tmp_1706 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1719;
   std::complex<double> tmp_1721;
   std::complex<double> tmp_1722;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1722 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1721 += tmp_1722;
   tmp_1706 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1721;
   std::complex<double> tmp_1723;
   std::complex<double> tmp_1724;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1725;
      std::complex<double> tmp_1726;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1726 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1725 += tmp_1726;
      tmp_1724 += (Conj(ZU(gI2,j2))) * tmp_1725;
   }
   tmp_1723 += tmp_1724;
   tmp_1706 += (std::complex<double>(0,1)*Conj(Lambdax)*KroneckerDelta(0,gO2)*
      KroneckerDelta(3,gO1)) * tmp_1723;
   std::complex<double> tmp_1727;
   std::complex<double> tmp_1728;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1729;
      std::complex<double> tmp_1730;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1730 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1729 += tmp_1730;
      tmp_1728 += (ZU(gI1,j2)) * tmp_1729;
   }
   tmp_1727 += tmp_1728;
   tmp_1706 += (std::complex<double>(0,1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      3,gO2)*Lambdax) * tmp_1727;
   std::complex<double> tmp_1731;
   std::complex<double> tmp_1732;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1733;
      std::complex<double> tmp_1734;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1735;
         std::complex<double> tmp_1736;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1736 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1735 += tmp_1736;
         tmp_1734 += (ZU(gI1,3 + j2)) * tmp_1735;
      }
      tmp_1733 += tmp_1734;
      tmp_1732 += (Conj(ZU(gI2,3 + j3))) * tmp_1733;
   }
   tmp_1731 += tmp_1732;
   tmp_1706 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1731;
   std::complex<double> tmp_1737;
   std::complex<double> tmp_1738;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1739;
      std::complex<double> tmp_1740;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1741;
         std::complex<double> tmp_1742;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1742 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1741 += tmp_1742;
         tmp_1740 += (Conj(ZU(gI2,j2))) * tmp_1741;
      }
      tmp_1739 += tmp_1740;
      tmp_1738 += (ZU(gI1,j3)) * tmp_1739;
   }
   tmp_1737 += tmp_1738;
   tmp_1706 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1737;
   result += (std::complex<double>(0,-1)) * tmp_1706;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1743;
   std::complex<double> tmp_1744;
   std::complex<double> tmp_1745;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1745 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1744 += tmp_1745;
   tmp_1743 += (std::complex<double>(0,-0.35355339059327373)*vd*KroneckerDelta(
      0,gO2)*Sqr(g2)) * tmp_1744;
   std::complex<double> tmp_1746;
   std::complex<double> tmp_1747;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1747 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1746 += tmp_1747;
   tmp_1743 += (std::complex<double>(0,-0.35355339059327373)*vu*KroneckerDelta(
      1,gO2)*Sqr(g2)) * tmp_1746;
   std::complex<double> tmp_1748;
   std::complex<double> tmp_1749;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1749 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1748 += tmp_1749;
   tmp_1743 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(2,gO2)*Sqr(g2))
      * tmp_1748;
   std::complex<double> tmp_1750;
   std::complex<double> tmp_1751;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1751 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1750 += tmp_1751;
   tmp_1743 += (std::complex<double>(0,0.5)*vT*KroneckerDelta(3,gO2)*Sqr(g2)) *
      tmp_1750;
   std::complex<double> tmp_1752;
   std::complex<double> tmp_1753;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1754;
      std::complex<double> tmp_1755;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1755 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1754 += tmp_1755;
      tmp_1753 += (Conj(ZD(gI2,j2))) * tmp_1754;
   }
   tmp_1752 += tmp_1753;
   tmp_1743 += (std::complex<double>(0,-0.5)*vT*Conj(Lambdax)*KroneckerDelta(0,
      gO2)) * tmp_1752;
   std::complex<double> tmp_1756;
   std::complex<double> tmp_1757;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1758;
      std::complex<double> tmp_1759;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1759 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1758 += tmp_1759;
      tmp_1757 += (Conj(ZD(gI2,j2))) * tmp_1758;
   }
   tmp_1756 += tmp_1757;
   tmp_1743 += (std::complex<double>(0,1)*Conj(Mu)*KroneckerDelta(0,gO2)) *
      tmp_1756;
   std::complex<double> tmp_1760;
   std::complex<double> tmp_1761;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1762;
      std::complex<double> tmp_1763;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1763 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1762 += tmp_1763;
      tmp_1761 += (Conj(ZD(gI2,j2))) * tmp_1762;
   }
   tmp_1760 += tmp_1761;
   tmp_1743 += (std::complex<double>(0,0.7071067811865475)*vd*Conj(Lambdax)*
      KroneckerDelta(2,gO2)) * tmp_1760;
   std::complex<double> tmp_1764;
   std::complex<double> tmp_1765;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1766;
      std::complex<double> tmp_1767;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1767 += ZU(gI1,3 + j1)*TYu(j1,j2);
      }
      tmp_1766 += tmp_1767;
      tmp_1765 += (Conj(ZD(gI2,j2))) * tmp_1766;
   }
   tmp_1764 += tmp_1765;
   tmp_1743 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_1764;
   std::complex<double> tmp_1768;
   std::complex<double> tmp_1769;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1770;
      std::complex<double> tmp_1771;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1771 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1770 += tmp_1771;
      tmp_1769 += (ZU(gI1,j2)) * tmp_1770;
   }
   tmp_1768 += tmp_1769;
   tmp_1743 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(1,gO2)*Lambdax)
      * tmp_1768;
   std::complex<double> tmp_1772;
   std::complex<double> tmp_1773;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1774;
      std::complex<double> tmp_1775;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1775 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1774 += tmp_1775;
      tmp_1773 += (ZU(gI1,j2)) * tmp_1774;
   }
   tmp_1772 += tmp_1773;
   tmp_1743 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1772;
   std::complex<double> tmp_1776;
   std::complex<double> tmp_1777;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1778;
      std::complex<double> tmp_1779;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1779 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1778 += tmp_1779;
      tmp_1777 += (ZU(gI1,j2)) * tmp_1778;
   }
   tmp_1776 += tmp_1777;
   tmp_1743 += (std::complex<double>(0,-0.7071067811865475)*vu*KroneckerDelta(3
      ,gO2)*Lambdax) * tmp_1776;
   std::complex<double> tmp_1780;
   std::complex<double> tmp_1781;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1782;
      std::complex<double> tmp_1783;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1783 += Conj(ZD(gI2,3 + j1))*Conj(TYd(j1,j2));
      }
      tmp_1782 += tmp_1783;
      tmp_1781 += (ZU(gI1,j2)) * tmp_1782;
   }
   tmp_1780 += tmp_1781;
   tmp_1743 += (std::complex<double>(0,1)*KroneckerDelta(0,gO2)) * tmp_1780;
   std::complex<double> tmp_1784;
   std::complex<double> tmp_1785;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1786;
      std::complex<double> tmp_1787;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1788;
         std::complex<double> tmp_1789;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1789 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1788 += tmp_1789;
         tmp_1787 += (ZU(gI1,3 + j2)) * tmp_1788;
      }
      tmp_1786 += tmp_1787;
      tmp_1785 += (Conj(ZD(gI2,3 + j3))) * tmp_1786;
   }
   tmp_1784 += tmp_1785;
   tmp_1743 += (std::complex<double>(0,0.7071067811865475)*vu*KroneckerDelta(0,
      gO2)) * tmp_1784;
   std::complex<double> tmp_1790;
   std::complex<double> tmp_1791;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1792;
      std::complex<double> tmp_1793;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1794;
         std::complex<double> tmp_1795;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1795 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1794 += tmp_1795;
         tmp_1793 += (ZU(gI1,3 + j2)) * tmp_1794;
      }
      tmp_1792 += tmp_1793;
      tmp_1791 += (Conj(ZD(gI2,3 + j3))) * tmp_1792;
   }
   tmp_1790 += tmp_1791;
   tmp_1743 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(1,
      gO2)) * tmp_1790;
   std::complex<double> tmp_1796;
   std::complex<double> tmp_1797;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1798;
      std::complex<double> tmp_1799;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1800;
         std::complex<double> tmp_1801;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1801 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1800 += tmp_1801;
         tmp_1799 += (Conj(ZD(gI2,j2))) * tmp_1800;
      }
      tmp_1798 += tmp_1799;
      tmp_1797 += (ZU(gI1,j3)) * tmp_1798;
   }
   tmp_1796 += tmp_1797;
   tmp_1743 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(0,
      gO2)) * tmp_1796;
   std::complex<double> tmp_1802;
   std::complex<double> tmp_1803;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1804;
      std::complex<double> tmp_1805;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1806;
         std::complex<double> tmp_1807;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1807 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1806 += tmp_1807;
         tmp_1805 += (Conj(ZD(gI2,j2))) * tmp_1806;
      }
      tmp_1804 += tmp_1805;
      tmp_1803 += (ZU(gI1,j3)) * tmp_1804;
   }
   tmp_1802 += tmp_1803;
   tmp_1743 += (std::complex<double>(0,0.7071067811865475)*vu*KroneckerDelta(1,
      gO2)) * tmp_1802;
   result += (std::complex<double>(0,-1)) * tmp_1743;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZA(gI2,0) +
      KroneckerDelta(1,gO2)*ZA(gI2,1) + 1.4142135623730951*(KroneckerDelta(2,gO2)
      - KroneckerDelta(3,gO2))*ZA(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) - KroneckerDelta(1,gO2)*ZH(
      gI2,1) + 1.4142135623730951*(KroneckerDelta(2,gO2) + KroneckerDelta(3,gO2))*
      ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVPHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-(KroneckerDelta(0,gO2)*(3.872983346207417*g1*Cos(ThetaW()) +
      5*g2*Sin(ThetaW()))*ZP(gI2,0)) - KroneckerDelta(1,gO2)*(3.872983346207417*g1
      *Cos(ThetaW()) + 5*g2*Sin(ThetaW()))*ZP(gI2,1) - 10*g2*Sin(ThetaW())*(
      KroneckerDelta(2,gO2)*ZP(gI2,2) + KroneckerDelta(3,gO2)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVZHpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-(KroneckerDelta(0,gO2)*(5*g2*Cos(ThetaW()) -
      3.872983346207417*g1*Sin(ThetaW()))*ZP(gI2,0)) - KroneckerDelta(1,gO2)*(5*g2
      *Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW()))*ZP(gI2,1) - 10*g2*Cos(
      ThetaW())*(KroneckerDelta(2,gO2)*ZP(gI2,2) + KroneckerDelta(3,gO2)*ZP(gI2,3)
      ));

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

std::complex<double> CLASSNAME::CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*KroneckerDelta(gI1,gI2)*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

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

double CLASSNAME::CpVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834*
      g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZhhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) - ZA(gI2,1)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM(gI2,0))*Cos(ThetaW())*UM(gI1,0) + Conj(UM(gI2,1))
      *(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UM(gI1,1) + 2*g2*
      Conj(UM(gI2,2))*Cos(ThetaW())*UM(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChaChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP(gI1,0))*Cos(ThetaW())*UP(gI2,0) + Conj(UP(gI1,1))
      *(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UP(gI2,1) + 2*g2*
      Conj(UP(gI1,2))*Cos(ThetaW())*UP(gI2,2));

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

std::complex<double> CLASSNAME::CpVZVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*((-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(
      g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*ZP(gI1,0)*ZP(gI2,0) +
      (-7.745966692414834*g1*g2*Cos(ThetaW())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(
      ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW())))*ZP(gI1,1)*ZP(gI2,1) + 20*Sqr(g2)*
      Sqr(Cos(ThetaW()))*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpVZconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*((-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*ZP(
      gI1,0)*ZP(gI2,0) + (-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW())
      )*ZP(gI1,1)*ZP(gI2,1) - 10*g2*Cos(ThetaW())*(ZP(gI1,2)*ZP(gI2,2) + ZP(gI1,3)
      *ZP(gI2,3)));

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

std::complex<double> CLASSNAME::CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1808;
   std::complex<double> tmp_1809;
   std::complex<double> tmp_1810;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1810 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1809 += tmp_1810;
   tmp_1808 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1809;
   std::complex<double> tmp_1811;
   std::complex<double> tmp_1812;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1812 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1811 += tmp_1812;
   tmp_1808 += (std::complex<double>(0,0.2581988897471611)*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())) * tmp_1811;
   std::complex<double> tmp_1813;
   std::complex<double> tmp_1814;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1814 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1813 += tmp_1814;
   tmp_1808 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1813;
   std::complex<double> tmp_1815;
   std::complex<double> tmp_1816;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1816 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1815 += tmp_1816;
   tmp_1808 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1815;
   result += (std::complex<double>(0,-1)) * tmp_1808;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1817;
   std::complex<double> tmp_1818;
   std::complex<double> tmp_1819;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1819 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1818 += tmp_1819;
   tmp_1817 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1818;
   std::complex<double> tmp_1820;
   std::complex<double> tmp_1821;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1821 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1820 += tmp_1821;
   tmp_1817 += (std::complex<double>(0,-0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1820;
   std::complex<double> tmp_1822;
   std::complex<double> tmp_1823;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1823 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1822 += tmp_1823;
   tmp_1817 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1822;
   std::complex<double> tmp_1824;
   std::complex<double> tmp_1825;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1825 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1824 += tmp_1825;
   tmp_1817 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1824;
   result += (std::complex<double>(0,-1)) * tmp_1817;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1826;
   std::complex<double> tmp_1827;
   std::complex<double> tmp_1828;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1828 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1827 += tmp_1828;
   tmp_1826 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1827;
   std::complex<double> tmp_1829;
   std::complex<double> tmp_1830;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1830 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1829 += tmp_1830;
   tmp_1826 += (std::complex<double>(0,-0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1829;
   std::complex<double> tmp_1831;
   std::complex<double> tmp_1832;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1832 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1831 += tmp_1832;
   tmp_1826 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1831;
   std::complex<double> tmp_1833;
   std::complex<double> tmp_1834;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1834 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1833 += tmp_1834;
   tmp_1826 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1833;
   result += (std::complex<double>(0,-1)) * tmp_1826;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1835;
   std::complex<double> tmp_1836;
   std::complex<double> tmp_1837;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1837 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1836 += tmp_1837;
   tmp_1835 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1836;
   std::complex<double> tmp_1838;
   std::complex<double> tmp_1839;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1839 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1838 += tmp_1839;
   tmp_1835 += (-3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1838;
   result += (0.16666666666666666) * tmp_1835;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1840;
   std::complex<double> tmp_1841;
   std::complex<double> tmp_1842;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1842 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1841 += tmp_1842;
   tmp_1840 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1841;
   std::complex<double> tmp_1843;
   std::complex<double> tmp_1844;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1844 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1843 += tmp_1844;
   tmp_1840 += (-(g2*Cos(ThetaW())) + 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1843;
   result += (0.5) * tmp_1840;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1845;
   std::complex<double> tmp_1846;
   std::complex<double> tmp_1847;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1847 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1846 += tmp_1847;
   tmp_1845 += (-3.0983866769659336*g1*Sin(ThetaW())) * tmp_1846;
   std::complex<double> tmp_1848;
   std::complex<double> tmp_1849;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1849 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1848 += tmp_1849;
   tmp_1845 += (3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1848;
   result += (0.16666666666666666) * tmp_1845;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(vd
      *ZH(gI2,0) + vu*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZconjVWmHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(-0.7745966692414834*g1*vd*Sin(ThetaW())*ZP(gI2,0) +
      0.7745966692414834*g1*vu*Sin(ThetaW())*ZP(gI2,1) + 1.4142135623730951*g2*vT*
      Cos(ThetaW())*(ZP(gI2,2) + ZP(gI2,3)));

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

   result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1) + 4*ZA(gI1,2
      )*ZA(gI2,2));

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

   result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1) + 4*ZH(gI1,2
      )*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1850;
   std::complex<double> tmp_1851;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1851 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1850 += tmp_1851;
   result += (-0.7071067811865475*g2) * tmp_1850;

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

   std::complex<double> tmp_1852;
   std::complex<double> tmp_1853;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1853 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1852 += tmp_1853;
   result += (0.7071067811865475*g2) * tmp_1852;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1) + 2*ZP(gI1,2
      )*ZP(gI2,2) + 2*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(ZA(gI2,0)*ZP(gI1,0) + ZA(gI2,1)*ZP
      (gI1,1) + 1.4142135623730951*ZA(gI2,2)*(ZP(gI1,2) - ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP(gI1,1) +
      1.4142135623730951*ZH(gI2,2)*(ZP(gI1,2) + ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*ZN(gI1,1) + 1.4142135623730951*Conj(UM(
      gI2,1))*ZN(gI1,2) + 2*Conj(UM(gI2,2))*ZN(gI1,4));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmChiChaPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(ZN(gI1,1))*UP(gI2,0)) + 0.7071067811865475*g2*Conj(ZN(gI1
      ,3))*UP(gI2,1) - g2*Conj(ZN(gI1,4))*UP(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1854;
   std::complex<double> tmp_1855;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1855 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1854 += tmp_1855;
   result += (0.5*Sqr(g2)) * tmp_1854;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1856;
   std::complex<double> tmp_1857;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1857 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1856 += tmp_1857;
   result += (0.5*Sqr(g2)) * tmp_1856;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1858;
   std::complex<double> tmp_1859;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1859 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1858 += tmp_1859;
   result += (0.5*Sqr(g2)) * tmp_1858;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1860;
   std::complex<double> tmp_1861;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1861 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1860 += tmp_1861;
   result += (0.7071067811865475*g2) * tmp_1860;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVWmhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1) + 4*vT*ZH(gI2,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVPHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(0.7745966692414834*g1*vd*Cos(ThetaW())*ZP(gI2,0) -
      0.7745966692414834*g1*vu*Cos(ThetaW())*ZP(gI2,1) + 1.4142135623730951*g2*vT*
      Sin(ThetaW())*(ZP(gI2,2) + ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmVZHpm(unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(-0.7745966692414834*g1*vd*Sin(ThetaW())*ZP(gI2,0) +
      0.7745966692414834*g1*vu*Sin(ThetaW())*ZP(gI2,1) + 1.4142135623730951*g2*vT*
      Cos(ThetaW())*(ZP(gI2,2) + ZP(gI2,3)));

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

std::complex<double> CLASSNAME::CpUChiconjSvFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5477225575051661*g1*KroneckerDelta(0,gO2)*ZV(gI1,gI2);
   }
   if (gI2 < 3) {
      result += -0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZV(gI1,gI2);
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
   std::complex<double> result;

   result = 0.1*(-5*g2*Conj(ZN(gI2,1))*KroneckerDelta(2,gO2)*ZH(gI1,0) + 5*Conj
      (ZN(gI2,4))*KroneckerDelta(3,gO2)*Lambdax*ZH(gI1,0) + 5*Conj(ZN(gI2,3))*
      KroneckerDelta(4,gO2)*Lambdax*ZH(gI1,0) - 3.872983346207417*g1*Conj(ZN(gI2,3
      ))*KroneckerDelta(0,gO2)*ZH(gI1,1) + 5*g2*Conj(ZN(gI2,3))*KroneckerDelta(1,
      gO2)*ZH(gI1,1) + 5*g2*Conj(ZN(gI2,1))*KroneckerDelta(3,gO2)*ZH(gI1,1) + 5*
      Conj(ZN(gI2,4))*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1) + 3.872983346207417*
      g1*Conj(ZN(gI2,0))*(KroneckerDelta(2,gO2)*ZH(gI1,0) - KroneckerDelta(3,gO2)*
      ZH(gI1,1)) + 5*Conj(ZN(gI2,3))*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,2) +
      Conj(ZN(gI2,2))*(3.872983346207417*g1*KroneckerDelta(0,gO2)*ZH(gI1,0) - 5*g2
      *KroneckerDelta(1,gO2)*ZH(gI1,0) + 5*KroneckerDelta(4,gO2)*Lambdax*ZH(gI1,1)
      + 5*KroneckerDelta(3,gO2)*Lambdax*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(3.872983346207417*g1*KroneckerDelta(0,gO1)*ZH(gI1,0)*ZN(gI2,2)
      - 5*g2*KroneckerDelta(1,gO1)*ZH(gI1,0)*ZN(gI2,2) + 5*Conj(Lambdax)*
      KroneckerDelta(4,gO1)*ZH(gI1,1)*ZN(gI2,2) + 5*Conj(Lambdax)*KroneckerDelta(4
      ,gO1)*ZH(gI1,0)*ZN(gI2,3) - 3.872983346207417*g1*KroneckerDelta(0,gO1)*ZH(
      gI1,1)*ZN(gI2,3) + 5*g2*KroneckerDelta(1,gO1)*ZH(gI1,1)*ZN(gI2,3) +
      KroneckerDelta(3,gO1)*(ZH(gI1,1)*(-3.872983346207417*g1*ZN(gI2,0) + 5*g2*ZN(
      gI2,1)) + 5*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,2) + ZH(gI1,0)*ZN(gI2,4))) +
      KroneckerDelta(2,gO1)*(ZH(gI1,0)*(3.872983346207417*g1*ZN(gI2,0) - 5*g2*ZN(
      gI2,1)) + 5*Conj(Lambdax)*(ZH(gI1,2)*ZN(gI2,3) + ZH(gI1,1)*ZN(gI2,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHpmChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(UM(gI2,2))*(-(KroneckerDelta(2,gO2)*Lambdax*ZP(gI1,1)) +
      1.4142135623730951*g2*KroneckerDelta(1,gO2)*ZP(gI1,2)) - g2*Conj(UM(gI2,0))*
      (KroneckerDelta(2,gO2)*ZP(gI1,0) + 1.4142135623730951*KroneckerDelta(4,gO2)*
      ZP(gI1,2)) + Conj(UM(gI2,1))*(0.5477225575051661*g1*KroneckerDelta(0,gO2)*ZP
      (gI1,0) + 0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZP(gI1,0) +
      0.7071067811865475*KroneckerDelta(4,gO2)*Lambdax*ZP(gI1,1) + KroneckerDelta(
      3,gO2)*Lambdax*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjHpmChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(3,gO1)*UP(gI2,0)*ZP(gI1,1)) + Conj(Lambdax)*(
      0.7071067811865475*KroneckerDelta(4,gO1)*UP(gI2,1)*ZP(gI1,0) +
      KroneckerDelta(3,gO1)*UP(gI2,2)*ZP(gI1,0) - KroneckerDelta(2,gO1)*UP(gI2,1)*
      ZP(gI1,2)) + 0.1414213562373095*(-3.872983346207417*g1*KroneckerDelta(0,gO1)
      *UP(gI2,1)*ZP(gI1,1) - 5*g2*(-2*KroneckerDelta(4,gO1)*UP(gI2,0)*ZP(gI1,3) +
      KroneckerDelta(1,gO1)*(UP(gI2,1)*ZP(gI1,1) + 2*UP(gI2,2)*ZP(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(5*g2*Conj(ZN(gI1,1))*KroneckerDelta(2,
      gO2)*ZA(gI2,0) + 5*Conj(ZN(gI1,4))*KroneckerDelta(3,gO2)*Lambdax*ZA(gI2,0) +
      5*Conj(ZN(gI1,3))*KroneckerDelta(4,gO2)*Lambdax*ZA(gI2,0) +
      3.872983346207417*g1*Conj(ZN(gI1,3))*KroneckerDelta(0,gO2)*ZA(gI2,1) - 5*g2*
      Conj(ZN(gI1,3))*KroneckerDelta(1,gO2)*ZA(gI2,1) - 5*g2*Conj(ZN(gI1,1))*
      KroneckerDelta(3,gO2)*ZA(gI2,1) + 5*Conj(ZN(gI1,4))*KroneckerDelta(2,gO2)*
      Lambdax*ZA(gI2,1) + 3.872983346207417*g1*Conj(ZN(gI1,0))*(-(KroneckerDelta(2
      ,gO2)*ZA(gI2,0)) + KroneckerDelta(3,gO2)*ZA(gI2,1)) + 5*Conj(ZN(gI1,3))*
      KroneckerDelta(2,gO2)*Lambdax*ZA(gI2,2) + Conj(ZN(gI1,2))*(
      -3.872983346207417*g1*KroneckerDelta(0,gO2)*ZA(gI2,0) + 5*(g2*KroneckerDelta
      (1,gO2)*ZA(gI2,0) + KroneckerDelta(4,gO2)*Lambdax*ZA(gI2,1) + KroneckerDelta
      (3,gO2)*Lambdax*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(3.872983346207417*g1*KroneckerDelta(0,
      gO1)*ZA(gI2,0)*ZN(gI1,2) - 5*g2*KroneckerDelta(1,gO1)*ZA(gI2,0)*ZN(gI1,2) -
      5*Conj(Lambdax)*KroneckerDelta(4,gO1)*ZA(gI2,1)*ZN(gI1,2) - 5*Conj(Lambdax)*
      KroneckerDelta(4,gO1)*ZA(gI2,0)*ZN(gI1,3) - 3.872983346207417*g1*
      KroneckerDelta(0,gO1)*ZA(gI2,1)*ZN(gI1,3) + 5*g2*KroneckerDelta(1,gO1)*ZA(
      gI2,1)*ZN(gI1,3) - KroneckerDelta(3,gO1)*(ZA(gI2,1)*(3.872983346207417*g1*ZN
      (gI1,0) - 5*g2*ZN(gI1,1)) + 5*Conj(Lambdax)*(ZA(gI2,2)*ZN(gI1,2) + ZA(gI2,0)
      *ZN(gI1,4))) + KroneckerDelta(2,gO1)*(ZA(gI2,0)*(3.872983346207417*g1*ZN(gI1
      ,0) - 5*g2*ZN(gI1,1)) - 5*Conj(Lambdax)*(ZA(gI2,2)*ZN(gI1,3) + ZA(gI2,1)*ZN(
      gI1,4))));

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1862;
   std::complex<double> tmp_1863;
   std::complex<double> tmp_1864;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1864 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1863 += tmp_1864;
   tmp_1862 += (std::complex<double>(0,-0.18257418583505536)*g1*KroneckerDelta(
      0,gO2)) * tmp_1863;
   std::complex<double> tmp_1865;
   std::complex<double> tmp_1866;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1866 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1865 += tmp_1866;
   tmp_1862 += (std::complex<double>(0,0.7071067811865475)*g2*KroneckerDelta(1,
      gO2)) * tmp_1865;
   std::complex<double> tmp_1867;
   std::complex<double> tmp_1868;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1869;
      std::complex<double> tmp_1870;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1870 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1869 += tmp_1870;
      tmp_1868 += (Conj(ZDL(gI2,j2))) * tmp_1869;
   }
   tmp_1867 += tmp_1868;
   tmp_1862 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1867;
   result += (std::complex<double>(0,-1)) * tmp_1862;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1871;
   std::complex<double> tmp_1872;
   std::complex<double> tmp_1873;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1873 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1872 += tmp_1873;
   tmp_1871 += (std::complex<double>(0,-0.3651483716701107)*g1*KroneckerDelta(0
      ,gO1)) * tmp_1872;
   std::complex<double> tmp_1874;
   std::complex<double> tmp_1875;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1876;
      std::complex<double> tmp_1877;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1877 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1876 += tmp_1877;
      tmp_1875 += (ZD(gI1,j2)) * tmp_1876;
   }
   tmp_1874 += tmp_1875;
   tmp_1871 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1874;
   result += (std::complex<double>(0,-1)) * tmp_1871;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1878;
   std::complex<double> tmp_1879;
   std::complex<double> tmp_1880;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1880 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1879 += tmp_1880;
   tmp_1878 += (std::complex<double>(0,0.5477225575051661)*g1*KroneckerDelta(0,
      gO2)) * tmp_1879;
   std::complex<double> tmp_1881;
   std::complex<double> tmp_1882;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1882 += Conj(ZEL(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1881 += tmp_1882;
   tmp_1878 += (std::complex<double>(0,0.7071067811865475)*g2*KroneckerDelta(1,
      gO2)) * tmp_1881;
   std::complex<double> tmp_1883;
   std::complex<double> tmp_1884;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1885;
      std::complex<double> tmp_1886;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1886 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1885 += tmp_1886;
      tmp_1884 += (Conj(ZEL(gI2,j2))) * tmp_1885;
   }
   tmp_1883 += tmp_1884;
   tmp_1878 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO2)) * tmp_1883;
   result += (std::complex<double>(0,-1)) * tmp_1878;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1887;
   std::complex<double> tmp_1888;
   std::complex<double> tmp_1889;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1889 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_1888 += tmp_1889;
   tmp_1887 += (std::complex<double>(0,-1.0954451150103321)*g1*KroneckerDelta(0
      ,gO1)) * tmp_1888;
   std::complex<double> tmp_1890;
   std::complex<double> tmp_1891;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1892;
      std::complex<double> tmp_1893;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1893 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1892 += tmp_1893;
      tmp_1891 += (ZE(gI1,j2)) * tmp_1892;
   }
   tmp_1890 += tmp_1891;
   tmp_1887 += (std::complex<double>(0,-1)*KroneckerDelta(2,gO1)) * tmp_1890;
   result += (std::complex<double>(0,-1)) * tmp_1887;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1894;
   std::complex<double> tmp_1895;
   std::complex<double> tmp_1896;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1896 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1895 += tmp_1896;
   tmp_1894 += (std::complex<double>(0,-0.18257418583505536)*g1*KroneckerDelta(
      0,gO2)) * tmp_1895;
   std::complex<double> tmp_1897;
   std::complex<double> tmp_1898;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1898 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1897 += tmp_1898;
   tmp_1894 += (std::complex<double>(0,-0.7071067811865475)*g2*KroneckerDelta(1
      ,gO2)) * tmp_1897;
   std::complex<double> tmp_1899;
   std::complex<double> tmp_1900;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1901;
      std::complex<double> tmp_1902;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1902 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1901 += tmp_1902;
      tmp_1900 += (Conj(ZUL(gI2,j2))) * tmp_1901;
   }
   tmp_1899 += tmp_1900;
   tmp_1894 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO2)) * tmp_1899;
   result += (std::complex<double>(0,-1)) * tmp_1894;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1903;
   std::complex<double> tmp_1904;
   std::complex<double> tmp_1905;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1905 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1904 += tmp_1905;
   tmp_1903 += (std::complex<double>(0,0.7302967433402214)*g1*KroneckerDelta(0,
      gO1)) * tmp_1904;
   std::complex<double> tmp_1906;
   std::complex<double> tmp_1907;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1908;
      std::complex<double> tmp_1909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1909 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1908 += tmp_1909;
      tmp_1907 += (ZU(gI1,j2)) * tmp_1908;
   }
   tmp_1906 += tmp_1907;
   tmp_1903 += (std::complex<double>(0,-1)*KroneckerDelta(3,gO1)) * tmp_1906;
   result += (std::complex<double>(0,-1)) * tmp_1903;

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(1,gO2)*UP(gI2,0)) + 0.7071067811865475*g2*
      KroneckerDelta(3,gO2)*UP(gI2,1) - g2*KroneckerDelta(4,gO2)*UP(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpUChiconjVWmChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM(gI2,0))*KroneckerDelta(1,gO1) +
      1.4142135623730951*Conj(UM(gI2,1))*KroneckerDelta(2,gO1) + 2*Conj(UM(gI2,2))
      *KroneckerDelta(4,gO1));

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

   result = std::complex<double>(0,0.5)*(Conj(UM(gI1,2))*(-1.4142135623730951*
      KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,0) + 2*g2*KroneckerDelta(0,gO2)*ZA(gI2,
      2)) + g2*Conj(UM(gI1,0))*(1.4142135623730951*KroneckerDelta(1,gO2)*ZA(gI2,1)
      - 2*KroneckerDelta(2,gO2)*ZA(gI2,2)) + Conj(UM(gI1,1))*(1.4142135623730951*
      g2*KroneckerDelta(0,gO2)*ZA(gI2,0) + 1.4142135623730951*KroneckerDelta(2,gO2
      )*Lambdax*ZA(gI2,1) + KroneckerDelta(1,gO2)*Lambdax*ZA(gI2,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaChaAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(-1.4142135623730951*Conj(Lambdax)*
      KroneckerDelta(2,gO1)*UP(gI1,1)*ZA(gI2,0) + KroneckerDelta(1,gO1)*(
      1.4142135623730951*g2*UP(gI1,0)*ZA(gI2,0) + Conj(Lambdax)*(
      1.4142135623730951*UP(gI1,2)*ZA(gI2,1) + UP(gI1,1)*ZA(gI2,2))) + g2*(2*
      KroneckerDelta(2,gO1)*UP(gI1,0)*ZA(gI2,2) + KroneckerDelta(0,gO1)*(
      1.4142135623730951*UP(gI1,1)*ZA(gI2,1) - 2*UP(gI1,2)*ZA(gI2,2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1910;
   std::complex<double> tmp_1911;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1911 += Conj(ZEL(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1910 += tmp_1911;
   result += (-(g2*KroneckerDelta(0,gO2))) * tmp_1910;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSvFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1912;
   std::complex<double> tmp_1913;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1914;
      std::complex<double> tmp_1915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1915 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1914 += tmp_1915;
      tmp_1913 += (ZV(gI1,j2)) * tmp_1914;
   }
   tmp_1912 += tmp_1913;
   result += (KroneckerDelta(1,gO1)) * tmp_1912;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(Conj(UM(gI2,2))*(1.4142135623730951*KroneckerDelta(1,gO2)*
      Lambdax*ZH(gI1,0) + 2*g2*KroneckerDelta(0,gO2)*ZH(gI1,2))) + g2*Conj(UM(gI2,
      0))*(-1.4142135623730951*KroneckerDelta(1,gO2)*ZH(gI1,1) + 2*KroneckerDelta(
      2,gO2)*ZH(gI1,2)) + Conj(UM(gI2,1))*(-1.4142135623730951*g2*KroneckerDelta(0
      ,gO2)*ZH(gI1,0) + 1.4142135623730951*KroneckerDelta(2,gO2)*Lambdax*ZH(gI1,1)
      + KroneckerDelta(1,gO2)*Lambdax*ZH(gI1,2)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChahhChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-1.4142135623730951*Conj(Lambdax)*KroneckerDelta(2,gO1)*UP(gI2
      ,1)*ZH(gI1,0) + KroneckerDelta(1,gO1)*(-1.4142135623730951*g2*UP(gI2,0)*ZH(
      gI1,0) + Conj(Lambdax)*(1.4142135623730951*UP(gI2,2)*ZH(gI1,1) + UP(gI2,1)*
      ZH(gI1,2))) - g2*(2*KroneckerDelta(2,gO1)*UP(gI2,0)*ZH(gI1,2) +
      KroneckerDelta(0,gO1)*(1.4142135623730951*UP(gI2,1)*ZH(gI1,1) - 2*UP(gI2,2)*
      ZH(gI1,2))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1916;
   std::complex<double> tmp_1917;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1918;
      std::complex<double> tmp_1919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1919 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1918 += tmp_1919;
      tmp_1917 += (Conj(ZD(gI2,j2))) * tmp_1918;
   }
   tmp_1916 += tmp_1917;
   result += (KroneckerDelta(1,gO2)) * tmp_1916;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChabarFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1920;
   std::complex<double> tmp_1921;
   std::complex<double> tmp_1922;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1922 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1921 += tmp_1922;
   tmp_1920 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO1)) * tmp_1921
      ;
   std::complex<double> tmp_1923;
   std::complex<double> tmp_1924;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1925;
      std::complex<double> tmp_1926;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1926 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1925 += tmp_1926;
      tmp_1924 += (ZUL(gI1,j2)) * tmp_1925;
   }
   tmp_1923 += tmp_1924;
   tmp_1920 += (std::complex<double>(0,1)*KroneckerDelta(1,gO1)) * tmp_1923;
   result += (std::complex<double>(0,-1)) * tmp_1920;

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

   std::complex<double> tmp_1927;
   std::complex<double> tmp_1928;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1928 += Conj(Ye(j1,gI1))*Conj(ZE(gI2,3 + j1));
   }
   tmp_1927 += tmp_1928;
   result += (KroneckerDelta(1,gO1)) * tmp_1927;
   if (gI1 < 3) {
      result += -(g2*Conj(ZE(gI2,gI1))*KroneckerDelta(0,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5477225575051661*g1*Conj(ZN(gI2,0))*KroneckerDelta(1,gO2)*ZP(gI1
      ,1) - 0.7071067811865475*g2*Conj(ZN(gI2,1))*KroneckerDelta(1,gO2)*ZP(gI1,1)
      + Conj(ZN(gI2,3))*(KroneckerDelta(2,gO2)*Lambdax*ZP(gI1,0) - g2*
      KroneckerDelta(0,gO2)*ZP(gI1,1)) - Conj(ZN(gI2,2))*KroneckerDelta(1,gO2)*
      Lambdax*ZP(gI1,2) - 1.4142135623730951*g2*Conj(ZN(gI2,1))*KroneckerDelta(2,
      gO2)*ZP(gI1,3) + 0.7071067811865475*Conj(ZN(gI2,4))*(KroneckerDelta(1,gO2)*
      Lambdax*ZP(gI1,0) + 2*g2*KroneckerDelta(0,gO2)*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = KroneckerDelta(2,gO1)*(-(Conj(Lambdax)*ZN(gI2,2)*ZP(gI1,1)) +
      1.4142135623730951*g2*ZN(gI2,1)*ZP(gI1,2)) - g2*KroneckerDelta(0,gO1)*(ZN(
      gI2,2)*ZP(gI1,0) + 1.4142135623730951*ZN(gI2,4)*ZP(gI1,2)) + KroneckerDelta(
      1,gO1)*(0.5477225575051661*g1*ZN(gI2,0)*ZP(gI1,0) + 0.7071067811865475*g2*ZN
      (gI2,1)*ZP(gI1,0) + 0.7071067811865475*Conj(Lambdax)*ZN(gI2,4)*ZP(gI1,1) +
      Conj(Lambdax)*ZN(gI2,3)*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1929;
   std::complex<double> tmp_1930;
   std::complex<double> tmp_1931;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1931 += Conj(ZDL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1930 += tmp_1931;
   tmp_1929 += (std::complex<double>(0,-1)*g2*KroneckerDelta(0,gO2)) * tmp_1930
      ;
   std::complex<double> tmp_1932;
   std::complex<double> tmp_1933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1934;
      std::complex<double> tmp_1935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1935 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1934 += tmp_1935;
      tmp_1933 += (Conj(ZDL(gI2,j2))) * tmp_1934;
   }
   tmp_1932 += tmp_1933;
   tmp_1929 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)) * tmp_1932;
   result += (std::complex<double>(0,-1)) * tmp_1929;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaconjSuFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1936;
   std::complex<double> tmp_1937;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1938;
      std::complex<double> tmp_1939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1939 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1938 += tmp_1939;
      tmp_1937 += (ZU(gI1,j2)) * tmp_1938;
   }
   tmp_1936 += tmp_1937;
   result += (KroneckerDelta(1,gO1)) * tmp_1936;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVPChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*UP(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW(
      )))*UP(gI2,1) + g2*KroneckerDelta(2,gO2)*Sin(ThetaW())*UP(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVPChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*(Conj(UM(gI2,0))*KroneckerDelta(0,gO1) + Conj(UM(gI2,2))*
      KroneckerDelta(2,gO1))*Sin(ThetaW()) + 0.1*Conj(UM(gI2,1))*KroneckerDelta(1,
      gO1)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*UP(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      )))*UP(gI2,1) + g2*Cos(ThetaW())*KroneckerDelta(2,gO2)*UP(gI2,2);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVZChaPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM(gI2,0))*Cos(ThetaW())*KroneckerDelta(0,gO1) + g2*Conj(UM
      (gI2,2))*Cos(ThetaW())*KroneckerDelta(2,gO1) + 0.1*Conj(UM(gI2,1))*
      KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      )));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO2)*ZN(gI2,1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZN(gI2,3) - g2*KroneckerDelta(2,gO2)*ZN(gI2,4);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChaVWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(ZN(gI2,1))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN(gI2,2))*KroneckerDelta(1,gO1) + 2*Conj(ZN(gI2,4))
      *KroneckerDelta(2,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1940;
      std::complex<double> tmp_1941;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1941 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1940 += tmp_1941;
      result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1940;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1942;
      std::complex<double> tmp_1943;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1943 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_1942 += tmp_1943;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1942;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1944;
      std::complex<double> tmp_1945;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1945 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1944 += tmp_1945;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1944;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1946;
      std::complex<double> tmp_1947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1947 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_1946 += tmp_1947;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1946;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1948;
      std::complex<double> tmp_1949;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1949 += Conj(ZV(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1948 += tmp_1949;
      result += (Conj(UM(gI2,1))) * tmp_1948;
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

std::complex<double> CLASSNAME::CpbarUFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.0954451150103321*g1*Conj(ZE(gI1,3 + gO2))*Conj(ZN(gI2,0))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1950;
      std::complex<double> tmp_1951;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1951 += Conj(ZE(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1950 += tmp_1951;
      result += (-Conj(ZN(gI2,2))) * tmp_1950;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZE(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZE(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1952;
      std::complex<double> tmp_1953;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1953 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1952 += tmp_1953;
      result += (-ZN(gI2,2)) * tmp_1952;
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
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.7745966692414834*g1*Sin(ThetaW())*ZER(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeVZFePL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZEL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += -0.3872983346207417*g1*Conj(ZEL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1954;
      std::complex<double> tmp_1955;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1955 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1954 += tmp_1955;
      result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1954;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1956;
      std::complex<double> tmp_1957;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1957 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_1956 += tmp_1957;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1956;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1958;
      std::complex<double> tmp_1959;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1959 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1958 += tmp_1959;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1958;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1960;
      std::complex<double> tmp_1961;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1961 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1960 += tmp_1961;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1960;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1962;
      std::complex<double> tmp_1963;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1963 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1962 += tmp_1963;
      result += (ZP(gI1,0)) * tmp_1962;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1964;
      std::complex<double> tmp_1965;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1965 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1964 += tmp_1965;
      result += (ZP(gI1,1)) * tmp_1964;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1966;
      std::complex<double> tmp_1967;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1967 += Conj(ZU(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1966 += tmp_1967;
      result += (Conj(UM(gI2,1))) * tmp_1966;
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
      std::complex<double> tmp_1968;
      std::complex<double> tmp_1969;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1969 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1968 += tmp_1969;
      result += (UP(gI2,1)) * tmp_1968;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.3651483716701107*g1*Conj(ZD(gI1,3 + gO2))*Conj(ZN(gI2,0))
         ;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1970;
      std::complex<double> tmp_1971;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1971 += Conj(ZD(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1970 += tmp_1971;
      result += (-Conj(ZN(gI2,2))) * tmp_1970;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZD(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZD(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1972;
      std::complex<double> tmp_1973;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1973 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1972 += tmp_1973;
      result += (-ZN(gI2,2)) * tmp_1972;
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
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.2581988897471611*g1*Sin(ThetaW())*ZDR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdVZFdPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5*g2*Conj(ZDL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZDL(gI2,gO1))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1974;
      std::complex<double> tmp_1975;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1975 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_1974 += tmp_1975;
      result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_1974;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1976;
      std::complex<double> tmp_1977;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1977 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_1976 += tmp_1977;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,1)) *
         tmp_1976;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1978;
      std::complex<double> tmp_1979;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1979 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1978 += tmp_1979;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1978;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1980;
      std::complex<double> tmp_1981;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1981 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1980 += tmp_1981;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1980;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1982;
      std::complex<double> tmp_1983;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1983 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1982 += tmp_1983;
      result += (Conj(UP(gI1,1))) * tmp_1982;
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
      std::complex<double> tmp_1984;
      std::complex<double> tmp_1985;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1985 += Conj(Yd(j1,gO1))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1984 += tmp_1985;
      result += (UM(gI1,1)) * tmp_1984;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1986;
      std::complex<double> tmp_1987;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1987 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1986 += tmp_1987;
      result += (ZP(gI1,1)) * tmp_1986;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1988;
      std::complex<double> tmp_1989;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1989 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1988 += tmp_1989;
      result += (ZP(gI1,0)) * tmp_1988;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7302967433402214*g1*Conj(ZN(gI2,0))*Conj(ZU(gI1,3 + gO2));
   }
   if (gO2 < 3) {
      std::complex<double> tmp_1990;
      std::complex<double> tmp_1991;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1991 += Conj(ZU(gI1,j2))*Yu(gO2,j2);
      }
      tmp_1990 += tmp_1991;
      result += (-Conj(ZN(gI2,3))) * tmp_1990;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZU(gI1,gO1))*ZN(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZU(gI1,gO1))*ZN(gI2,1);
   }
   if (gO1 < 3) {
      std::complex<double> tmp_1992;
      std::complex<double> tmp_1993;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1993 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1992 += tmp_1993;
      result += (-ZN(gI2,3)) * tmp_1992;
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
   std::complex<double> result;

   if (gI2 < 3) {
      result += 0.5163977794943222*g1*Sin(ThetaW())*ZUR(gI2,gO2);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuVZFuPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 3) {
      result += -0.5*g2*Conj(ZUL(gI2,gO1))*Cos(ThetaW());
   }
   if (gI2 < 3) {
      result += 0.12909944487358055*g1*Conj(ZUL(gI2,gO1))*Sin(ThetaW());
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

   std::complex<double> tmp_1994;
   std::complex<double> tmp_1995;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1995 += Conj(ZDL(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1994 += tmp_1995;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_1994;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1996;
   std::complex<double> tmp_1997;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1997 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1996 += tmp_1997;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_1996;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1998;
   std::complex<double> tmp_1999;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1999 += Conj(ZUL(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1998 += tmp_1999;
   result += (-1.4142135623730951*g3*PhaseGlu) * tmp_1998;

   return result;
}

std::complex<double> CLASSNAME::CpGluconjSuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2000;
   std::complex<double> tmp_2001;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2001 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_2000 += tmp_2001;
   result += (1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2000;

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

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2002;
   std::complex<double> tmp_2003;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2004;
      std::complex<double> tmp_2005;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2005 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2004 += tmp_2005;
      tmp_2003 += (Conj(ZEL(gI1,j2))) * tmp_2004;
   }
   tmp_2002 += tmp_2003;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) * tmp_2002
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2006;
   std::complex<double> tmp_2007;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2008;
      std::complex<double> tmp_2009;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2009 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_2008 += tmp_2009;
      tmp_2007 += (ZEL(gO1,j2)) * tmp_2008;
   }
   tmp_2006 += tmp_2007;
   result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) * tmp_2006;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2010;
   std::complex<double> tmp_2011;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2012;
      std::complex<double> tmp_2013;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2013 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2012 += tmp_2013;
      tmp_2011 += (Conj(ZEL(gI2,j2))) * tmp_2012;
   }
   tmp_2010 += tmp_2011;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2010;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2014;
   std::complex<double> tmp_2015;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2016;
      std::complex<double> tmp_2017;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2017 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_2016 += tmp_2017;
      tmp_2015 += (ZEL(gO1,j2)) * tmp_2016;
   }
   tmp_2014 += tmp_2015;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2014;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2018;
   std::complex<double> tmp_2019;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2020;
      std::complex<double> tmp_2021;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2021 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2020 += tmp_2021;
      tmp_2019 += (Conj(ZV(gI1,j2))) * tmp_2020;
   }
   tmp_2018 += tmp_2019;
   result += (Conj(UM(gI2,1))) * tmp_2018;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSvChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2022;
   std::complex<double> tmp_2023;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2023 += Conj(ZV(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_2022 += tmp_2023;
   result += (-(g2*UP(gI2,0))) * tmp_2022;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2024;
   std::complex<double> tmp_2025;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2025 += Conj(ZER(gO2,j1))*Ye(j1,gI2);
   }
   tmp_2024 += tmp_2025;
   result += (ZP(gI1,0)) * tmp_2024;

   return result;
}

double CLASSNAME::CpbarFeHpmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2026;
   std::complex<double> tmp_2027;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2027 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_2026 += tmp_2027;
   result += (-1.0954451150103321*g1*Conj(ZN(gI2,0))) * tmp_2026;
   std::complex<double> tmp_2028;
   std::complex<double> tmp_2029;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2030;
      std::complex<double> tmp_2031;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2031 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_2030 += tmp_2031;
      tmp_2029 += (Conj(ZE(gI1,j2))) * tmp_2030;
   }
   tmp_2028 += tmp_2029;
   result += (-Conj(ZN(gI2,2))) * tmp_2028;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2032;
   std::complex<double> tmp_2033;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2033 += Conj(ZE(gI1,j1))*ZEL(gO1,j1);
   }
   tmp_2032 += tmp_2033;
   result += (0.7071067811865475*(0.7745966692414834*g1*ZN(gI2,0) + g2*ZN(gI2,1
      ))) * tmp_2032;
   std::complex<double> tmp_2034;
   std::complex<double> tmp_2035;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2036;
      std::complex<double> tmp_2037;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2037 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_2036 += tmp_2037;
      tmp_2035 += (ZEL(gO1,j2)) * tmp_2036;
   }
   tmp_2034 += tmp_2035;
   result += (-ZN(gI2,2)) * tmp_2034;

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

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2038;
   std::complex<double> tmp_2039;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2040;
      std::complex<double> tmp_2041;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2041 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2040 += tmp_2041;
      tmp_2039 += (Conj(ZDL(gI1,j2))) * tmp_2040;
   }
   tmp_2038 += tmp_2039;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) * tmp_2038
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2042;
   std::complex<double> tmp_2043;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2044;
      std::complex<double> tmp_2045;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2045 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_2044 += tmp_2045;
      tmp_2043 += (ZDL(gO1,j2)) * tmp_2044;
   }
   tmp_2042 += tmp_2043;
   result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) * tmp_2042;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2046;
   std::complex<double> tmp_2047;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2048;
      std::complex<double> tmp_2049;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2049 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2048 += tmp_2049;
      tmp_2047 += (Conj(ZDL(gI2,j2))) * tmp_2048;
   }
   tmp_2046 += tmp_2047;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2046;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2050;
   std::complex<double> tmp_2051;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2052;
      std::complex<double> tmp_2053;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2053 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2052 += tmp_2053;
      tmp_2051 += (ZDL(gO1,j2)) * tmp_2052;
   }
   tmp_2050 += tmp_2051;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_2050;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2054;
   std::complex<double> tmp_2055;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2056;
      std::complex<double> tmp_2057;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2057 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2056 += tmp_2057;
      tmp_2055 += (Conj(ZUL(gI2,j2))) * tmp_2056;
   }
   tmp_2054 += tmp_2055;
   result += (ZP(gI1,0)) * tmp_2054;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2058;
   std::complex<double> tmp_2059;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2060;
      std::complex<double> tmp_2061;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2061 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2060 += tmp_2061;
      tmp_2059 += (ZDL(gO1,j2)) * tmp_2060;
   }
   tmp_2058 += tmp_2059;
   result += (ZP(gI1,1)) * tmp_2058;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2062;
   std::complex<double> tmp_2063;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2064;
      std::complex<double> tmp_2065;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2065 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2064 += tmp_2065;
      tmp_2063 += (Conj(ZU(gI1,j2))) * tmp_2064;
   }
   tmp_2062 += tmp_2063;
   result += (Conj(UM(gI2,1))) * tmp_2062;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuChaPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2066;
   std::complex<double> tmp_2067;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2067 += Conj(ZU(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2066 += tmp_2067;
   result += (-(g2*UP(gI2,0))) * tmp_2066;
   std::complex<double> tmp_2068;
   std::complex<double> tmp_2069;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2070;
      std::complex<double> tmp_2071;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2071 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2070 += tmp_2071;
      tmp_2069 += (ZDL(gO1,j2)) * tmp_2070;
   }
   tmp_2068 += tmp_2069;
   result += (UP(gI2,1)) * tmp_2068;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2072;
   std::complex<double> tmp_2073;
   std::complex<double> tmp_2074;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2074 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2073 += tmp_2074;
   tmp_2072 += (-1.0954451150103321*g1*Conj(ZN(gI2,0))) * tmp_2073;
   std::complex<double> tmp_2075;
   std::complex<double> tmp_2076;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2077;
      std::complex<double> tmp_2078;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2078 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2077 += tmp_2078;
      tmp_2076 += (Conj(ZD(gI1,j2))) * tmp_2077;
   }
   tmp_2075 += tmp_2076;
   tmp_2072 += (-3*Conj(ZN(gI2,2))) * tmp_2075;
   result += (0.3333333333333333) * tmp_2072;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2079;
   std::complex<double> tmp_2080;
   std::complex<double> tmp_2081;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2081 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2080 += tmp_2081;
   tmp_2079 += (-1.4142135623730951*(0.7745966692414834*g1*ZN(gI2,0) - 3*g2*ZN(
      gI2,1))) * tmp_2080;
   std::complex<double> tmp_2082;
   std::complex<double> tmp_2083;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2084;
      std::complex<double> tmp_2085;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2085 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2084 += tmp_2085;
      tmp_2083 += (ZDL(gO1,j2)) * tmp_2084;
   }
   tmp_2082 += tmp_2083;
   tmp_2079 += (-6*ZN(gI2,2)) * tmp_2082;
   result += (0.16666666666666666) * tmp_2079;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2086;
   std::complex<double> tmp_2087;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2087 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2086 += tmp_2087;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2086;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2088;
   std::complex<double> tmp_2089;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2089 += Conj(ZD(gI1,j1))*ZDL(gO1,j1);
   }
   tmp_2088 += tmp_2089;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2088;

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

   std::complex<double> tmp_2090;
   std::complex<double> tmp_2091;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2091 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2090 += tmp_2091;
   result += (-0.7071067811865475*g2) * tmp_2090;

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

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2092;
   std::complex<double> tmp_2093;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2094;
      std::complex<double> tmp_2095;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2095 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2094 += tmp_2095;
      tmp_2093 += (Conj(ZUL(gI1,j2))) * tmp_2094;
   }
   tmp_2092 += tmp_2093;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,1)) * tmp_2092
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2096;
   std::complex<double> tmp_2097;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2098;
      std::complex<double> tmp_2099;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2099 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_2098 += tmp_2099;
      tmp_2097 += (ZUL(gO1,j2)) * tmp_2098;
   }
   tmp_2096 += tmp_2097;
   result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,1)) * tmp_2096;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2100;
   std::complex<double> tmp_2101;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2102;
      std::complex<double> tmp_2103;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2103 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2102 += tmp_2103;
      tmp_2101 += (Conj(ZUL(gI2,j2))) * tmp_2102;
   }
   tmp_2100 += tmp_2101;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2100;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2104;
   std::complex<double> tmp_2105;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2106;
      std::complex<double> tmp_2107;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2107 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2106 += tmp_2107;
      tmp_2105 += (ZUL(gO1,j2)) * tmp_2106;
   }
   tmp_2104 += tmp_2105;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2104;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2108;
   std::complex<double> tmp_2109;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2110;
      std::complex<double> tmp_2111;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2111 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2110 += tmp_2111;
      tmp_2109 += (Conj(ZD(gI2,j2))) * tmp_2110;
   }
   tmp_2108 += tmp_2109;
   result += (Conj(UP(gI1,1))) * tmp_2108;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChaSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2112;
   std::complex<double> tmp_2113;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2113 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2112 += tmp_2113;
   result += (-(g2*UM(gI1,0))) * tmp_2112;
   std::complex<double> tmp_2114;
   std::complex<double> tmp_2115;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2116;
      std::complex<double> tmp_2117;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2117 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_2116 += tmp_2117;
      tmp_2115 += (ZUL(gO1,j2)) * tmp_2116;
   }
   tmp_2114 += tmp_2115;
   result += (UM(gI1,1)) * tmp_2114;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2118;
   std::complex<double> tmp_2119;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2120;
      std::complex<double> tmp_2121;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2121 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2120 += tmp_2121;
      tmp_2119 += (Conj(ZDL(gI2,j2))) * tmp_2120;
   }
   tmp_2118 += tmp_2119;
   result += (ZP(gI1,1)) * tmp_2118;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2122;
   std::complex<double> tmp_2123;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2124;
      std::complex<double> tmp_2125;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2125 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2124 += tmp_2125;
      tmp_2123 += (ZUL(gO1,j2)) * tmp_2124;
   }
   tmp_2122 += tmp_2123;
   result += (ZP(gI1,0)) * tmp_2122;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2126;
   std::complex<double> tmp_2127;
   std::complex<double> tmp_2128;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2128 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2127 += tmp_2128;
   tmp_2126 += (2.1908902300206643*g1*Conj(ZN(gI2,0))) * tmp_2127;
   std::complex<double> tmp_2129;
   std::complex<double> tmp_2130;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2131;
      std::complex<double> tmp_2132;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2132 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2131 += tmp_2132;
      tmp_2130 += (Conj(ZU(gI1,j2))) * tmp_2131;
   }
   tmp_2129 += tmp_2130;
   tmp_2126 += (-3*Conj(ZN(gI2,3))) * tmp_2129;
   result += (0.3333333333333333) * tmp_2126;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2133;
   std::complex<double> tmp_2134;
   std::complex<double> tmp_2135;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2135 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2134 += tmp_2135;
   tmp_2133 += (-1.4142135623730951*(0.7745966692414834*g1*ZN(gI2,0) + 3*g2*ZN(
      gI2,1))) * tmp_2134;
   std::complex<double> tmp_2136;
   std::complex<double> tmp_2137;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2138;
      std::complex<double> tmp_2139;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2139 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2138 += tmp_2139;
      tmp_2137 += (ZUL(gO1,j2)) * tmp_2138;
   }
   tmp_2136 += tmp_2137;
   tmp_2133 += (-6*ZN(gI2,3)) * tmp_2136;
   result += (0.16666666666666666) * tmp_2133;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2140;
   std::complex<double> tmp_2141;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2141 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2140 += tmp_2141;
   result += (1.4142135623730951*g3*PhaseGlu) * tmp_2140;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPR(unsigned gO1, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2142;
   std::complex<double> tmp_2143;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2143 += Conj(ZU(gI1,j1))*ZUL(gO1,j1);
   }
   tmp_2142 += tmp_2143;
   result += (-1.4142135623730951*g3*Conj(PhaseGlu)) * tmp_2142;

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

   std::complex<double> tmp_2144;
   std::complex<double> tmp_2145;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2145 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2144 += tmp_2145;
   result += (-0.7071067811865475*g2) * tmp_2144;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_2146;
   std::complex<double> tmp_2147;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2147 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2146 += tmp_2147;
   result += (-0.5) * tmp_2146;
   std::complex<double> tmp_2148;
   std::complex<double> tmp_2149;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2149 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2148 += tmp_2149;
   result += (-1) * tmp_2148;
   std::complex<double> tmp_2150;
   std::complex<double> tmp_2151;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2151 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2150 += tmp_2151;
   result += (-0.5) * tmp_2150;
   std::complex<double> tmp_2152;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2153;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2153 += (Conj(CpconjUSdFuChaPL(gO2,gI1,gI2))*
            CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPR(gO2,gI1,gI2))*
            CpconjUSdFuChaPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha(gI2));
      }
      tmp_2152 += tmp_2153;
   }
   result += tmp_2152;
   std::complex<double> tmp_2154;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2155;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2155 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_2154 += tmp_2155;
   }
   result += tmp_2154;
   std::complex<double> tmp_2156;
   std::complex<double> tmp_2157;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2158;
      std::complex<double> tmp_2159;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2159 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2158 += tmp_2159;
      tmp_2157 += (MFd(gI1)) * tmp_2158;
   }
   tmp_2156 += tmp_2157;
   result += (-2) * tmp_2156;
   std::complex<double> tmp_2160;
   std::complex<double> tmp_2161;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2162;
      std::complex<double> tmp_2163;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2163 += B0(p,MFu(gI1),MCha(gI2))*(Conj(CpconjUSdFuChaPR(gO2,
            gI1,gI2))*CpconjUSdFuChaPL(gO1,gI1,gI2) + Conj(CpconjUSdFuChaPL(gO2,
            gI1,gI2))*CpconjUSdFuChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2162 += tmp_2163;
      tmp_2161 += (MFu(gI1)) * tmp_2162;
   }
   tmp_2160 += tmp_2161;
   result += (-2) * tmp_2160;
   std::complex<double> tmp_2164;
   std::complex<double> tmp_2165;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2165 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2164 += tmp_2165;
   result += (-1) * tmp_2164;
   std::complex<double> tmp_2166;
   std::complex<double> tmp_2167;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2167 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2166 += tmp_2167;
   result += (-1) * tmp_2166;
   std::complex<double> tmp_2168;
   std::complex<double> tmp_2169;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2169 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2168 += tmp_2169;
   result += (-1) * tmp_2168;
   std::complex<double> tmp_2170;
   std::complex<double> tmp_2171;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2171 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2170 += tmp_2171;
   result += (-1) * tmp_2170;
   std::complex<double> tmp_2172;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2173;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2173 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_2172 += tmp_2173;
   }
   result += tmp_2172;
   std::complex<double> tmp_2174;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2175;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2175 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_2174 += tmp_2175;
   }
   result += tmp_2174;
   std::complex<double> tmp_2176;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2177;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2177 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_2176 += tmp_2177;
   }
   result += tmp_2176;
   std::complex<double> tmp_2178;
   std::complex<double> tmp_2179;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2179 += (Conj(CpconjUSdGluFdPL(gO2,1,gI2))*CpconjUSdGluFdPL(gO1,1,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,1,gI2))*CpconjUSdGluFdPR(gO1,1,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_2178 += tmp_2179;
   result += (1.3333333333333333) * tmp_2178;
   std::complex<double> tmp_2180;
   std::complex<double> tmp_2181;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2181 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_2180 += tmp_2181;
   result += (1.3333333333333333) * tmp_2180;
   std::complex<double> tmp_2182;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2182 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_2182;
   std::complex<double> tmp_2183;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2183 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_2183;
   std::complex<double> tmp_2184;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2184 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_2184;
   std::complex<double> tmp_2185;
   std::complex<double> tmp_2186;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2186 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,1,gI2))*
         CpconjUSdGluFdPL(gO1,1,gI2) + Conj(CpconjUSdGluFdPL(gO2,1,gI2))*
         CpconjUSdGluFdPR(gO1,1,gI2))*MFd(gI2);
   }
   tmp_2185 += tmp_2186;
   result += (-2.6666666666666665*MGlu) * tmp_2185;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_2187;
   std::complex<double> tmp_2188;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2188 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2187 += tmp_2188;
   result += (-0.5) * tmp_2187;
   std::complex<double> tmp_2189;
   std::complex<double> tmp_2190;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2190 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2189 += tmp_2190;
   result += (-1) * tmp_2189;
   std::complex<double> tmp_2191;
   std::complex<double> tmp_2192;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2192 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2191 += tmp_2192;
   result += (-0.5) * tmp_2191;
   std::complex<double> tmp_2193;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2194;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2194 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_2193 += tmp_2194;
   }
   result += tmp_2193;
   std::complex<double> tmp_2195;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2196;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2196 += (Conj(CpconjUSvbarChaFePL(gO2,gI1,gI2))*
            CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(CpconjUSvbarChaFePR(gO2,gI1,
            gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFe(gI2));
      }
      tmp_2195 += tmp_2196;
   }
   result += tmp_2195;
   std::complex<double> tmp_2197;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2198;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2198 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_2197 += tmp_2198;
   }
   result += tmp_2197;
   std::complex<double> tmp_2199;
   std::complex<double> tmp_2200;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2201;
      std::complex<double> tmp_2202;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2202 += B0(p,MCha(gI1),MFe(gI2))*(Conj(CpconjUSvbarChaFePR(
            gO2,gI1,gI2))*CpconjUSvbarChaFePL(gO1,gI1,gI2) + Conj(
            CpconjUSvbarChaFePL(gO2,gI1,gI2))*CpconjUSvbarChaFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2201 += tmp_2202;
      tmp_2200 += (MCha(gI1)) * tmp_2201;
   }
   tmp_2199 += tmp_2200;
   result += (-2) * tmp_2199;
   std::complex<double> tmp_2203;
   std::complex<double> tmp_2204;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2205;
      std::complex<double> tmp_2206;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2206 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2205 += tmp_2206;
      tmp_2204 += (MFv(gI1)) * tmp_2205;
   }
   tmp_2203 += tmp_2204;
   result += (-2) * tmp_2203;
   std::complex<double> tmp_2207;
   std::complex<double> tmp_2208;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2208 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2207 += tmp_2208;
   result += (-1) * tmp_2207;
   std::complex<double> tmp_2209;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2210;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2210 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_2209 += tmp_2210;
   }
   result += tmp_2209;
   std::complex<double> tmp_2211;
   std::complex<double> tmp_2212;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2212 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2211 += tmp_2212;
   result += (-3) * tmp_2211;
   std::complex<double> tmp_2213;
   std::complex<double> tmp_2214;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2214 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2213 += tmp_2214;
   result += (-1) * tmp_2213;
   std::complex<double> tmp_2215;
   std::complex<double> tmp_2216;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2216 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2215 += tmp_2216;
   result += (-3) * tmp_2215;
   std::complex<double> tmp_2217;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2217 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_2217;
   std::complex<double> tmp_2218;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2218 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_2218;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_2219;
   std::complex<double> tmp_2220;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2220 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2219 += tmp_2220;
   result += (-0.5) * tmp_2219;
   std::complex<double> tmp_2221;
   std::complex<double> tmp_2222;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2222 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2221 += tmp_2222;
   result += (-1) * tmp_2221;
   std::complex<double> tmp_2223;
   std::complex<double> tmp_2224;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2224 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2223 += tmp_2224;
   result += (-0.5) * tmp_2223;
   std::complex<double> tmp_2225;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2226;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2226 += (Conj(CpconjUSubarChaFdPL(gO2,gI1,gI2))*
            CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(CpconjUSubarChaFdPR(gO2,gI1,
            gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MFd(gI2));
      }
      tmp_2225 += tmp_2226;
   }
   result += tmp_2225;
   std::complex<double> tmp_2227;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2228;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2228 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_2227 += tmp_2228;
   }
   result += tmp_2227;
   std::complex<double> tmp_2229;
   std::complex<double> tmp_2230;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2231;
      std::complex<double> tmp_2232;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2232 += B0(p,MCha(gI1),MFd(gI2))*(Conj(CpconjUSubarChaFdPR(
            gO2,gI1,gI2))*CpconjUSubarChaFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChaFdPL(gO2,gI1,gI2))*CpconjUSubarChaFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2231 += tmp_2232;
      tmp_2230 += (MCha(gI1)) * tmp_2231;
   }
   tmp_2229 += tmp_2230;
   result += (-2) * tmp_2229;
   std::complex<double> tmp_2233;
   std::complex<double> tmp_2234;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2235;
      std::complex<double> tmp_2236;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2236 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2235 += tmp_2236;
      tmp_2234 += (MFu(gI1)) * tmp_2235;
   }
   tmp_2233 += tmp_2234;
   result += (-2) * tmp_2233;
   std::complex<double> tmp_2237;
   std::complex<double> tmp_2238;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2238 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2237 += tmp_2238;
   result += (-1) * tmp_2237;
   std::complex<double> tmp_2239;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2240;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2240 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_2239 += tmp_2240;
   }
   result += tmp_2239;
   std::complex<double> tmp_2241;
   std::complex<double> tmp_2242;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2242 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2241 += tmp_2242;
   result += (-1) * tmp_2241;
   std::complex<double> tmp_2243;
   std::complex<double> tmp_2244;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2244 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2243 += tmp_2244;
   result += (-1) * tmp_2243;
   std::complex<double> tmp_2245;
   std::complex<double> tmp_2246;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2246 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2245 += tmp_2246;
   result += (-1) * tmp_2245;
   std::complex<double> tmp_2247;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2248;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2248 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_2247 += tmp_2248;
   }
   result += tmp_2247;
   std::complex<double> tmp_2249;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2250;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2250 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_2249 += tmp_2250;
   }
   result += tmp_2249;
   std::complex<double> tmp_2251;
   std::complex<double> tmp_2252;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2252 += (Conj(CpconjUSuGluFuPL(gO2,1,gI2))*CpconjUSuGluFuPL(gO1,1,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,1,gI2))*CpconjUSuGluFuPR(gO1,1,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_2251 += tmp_2252;
   result += (1.3333333333333333) * tmp_2251;
   std::complex<double> tmp_2253;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2253 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_2253;
   std::complex<double> tmp_2254;
   std::complex<double> tmp_2255;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2255 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_2254 += tmp_2255;
   result += (1.3333333333333333) * tmp_2254;
   std::complex<double> tmp_2256;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2256 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_2256;
   std::complex<double> tmp_2257;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2257 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_2257;
   std::complex<double> tmp_2258;
   std::complex<double> tmp_2259;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2259 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,1,gI2))*
         CpconjUSuGluFuPL(gO1,1,gI2) + Conj(CpconjUSuGluFuPL(gO2,1,gI2))*
         CpconjUSuGluFuPR(gO1,1,gI2))*MFu(gI2);
   }
   tmp_2258 += tmp_2259;
   result += (-2.6666666666666665*MGlu) * tmp_2258;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_2260;
   std::complex<double> tmp_2261;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2261 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2260 += tmp_2261;
   result += (-0.5) * tmp_2260;
   std::complex<double> tmp_2262;
   std::complex<double> tmp_2263;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2263 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2262 += tmp_2263;
   result += (-1) * tmp_2262;
   std::complex<double> tmp_2264;
   std::complex<double> tmp_2265;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2265 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2264 += tmp_2265;
   result += (-0.5) * tmp_2264;
   std::complex<double> tmp_2266;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2267;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2267 += (Conj(CpconjUSeFvChaPL(gO2,gI1,gI2))*
            CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPR(gO2,gI1,gI2))*
            CpconjUSeFvChaPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha(gI2));
      }
      tmp_2266 += tmp_2267;
   }
   result += tmp_2266;
   std::complex<double> tmp_2268;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2269;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2269 += B0(p,MSv(gI1),MHpm(gI2))*Conj(CpconjUSeSvHpm(gO2,gI1
            ,gI2))*CpconjUSeSvHpm(gO1,gI1,gI2);
      }
      tmp_2268 += tmp_2269;
   }
   result += tmp_2268;
   std::complex<double> tmp_2270;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2271;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2271 += (Conj(CpconjUSeFeChiPL(gO2,gI1,gI2))*
            CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPR(gO2,gI1,gI2))*
            CpconjUSeFeChiPR(gO1,gI1,gI2))*G0(p,MFe(gI1),MChi(gI2));
      }
      tmp_2270 += tmp_2271;
   }
   result += tmp_2270;
   std::complex<double> tmp_2272;
   std::complex<double> tmp_2273;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2274;
      std::complex<double> tmp_2275;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2275 += B0(p,MFe(gI1),MChi(gI2))*(Conj(CpconjUSeFeChiPR(gO2,
            gI1,gI2))*CpconjUSeFeChiPL(gO1,gI1,gI2) + Conj(CpconjUSeFeChiPL(gO2,
            gI1,gI2))*CpconjUSeFeChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2274 += tmp_2275;
      tmp_2273 += (MFe(gI1)) * tmp_2274;
   }
   tmp_2272 += tmp_2273;
   result += (-2) * tmp_2272;
   std::complex<double> tmp_2276;
   std::complex<double> tmp_2277;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2278;
      std::complex<double> tmp_2279;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2279 += B0(p,MFv(gI1),MCha(gI2))*(Conj(CpconjUSeFvChaPR(gO2,
            gI1,gI2))*CpconjUSeFvChaPL(gO1,gI1,gI2) + Conj(CpconjUSeFvChaPL(gO2,
            gI1,gI2))*CpconjUSeFvChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2278 += tmp_2279;
      tmp_2277 += (MFv(gI1)) * tmp_2278;
   }
   tmp_2276 += tmp_2277;
   result += (-2) * tmp_2276;
   std::complex<double> tmp_2280;
   std::complex<double> tmp_2281;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2281 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2280 += tmp_2281;
   result += (-1) * tmp_2280;
   std::complex<double> tmp_2282;
   std::complex<double> tmp_2283;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2283 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2282 += tmp_2283;
   result += (-3) * tmp_2282;
   std::complex<double> tmp_2284;
   std::complex<double> tmp_2285;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2285 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2284 += tmp_2285;
   result += (-1) * tmp_2284;
   std::complex<double> tmp_2286;
   std::complex<double> tmp_2287;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2287 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2286 += tmp_2287;
   result += (-3) * tmp_2286;
   std::complex<double> tmp_2288;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2289;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2289 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_2288 += tmp_2289;
   }
   result += tmp_2288;
   std::complex<double> tmp_2290;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2291;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2291 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_2290 += tmp_2291;
   }
   result += tmp_2290;
   std::complex<double> tmp_2292;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2292 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_2292;
   std::complex<double> tmp_2293;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2293 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_2293;
   std::complex<double> tmp_2294;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2294 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_2294;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_hh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmCgWmC(gO1)*CpUhhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUhhbargWmgWm(gO1)*CpUhhbargWmgWm(gO2));
   result += -(B0(p,MVZ,MVZ)*CpUhhbargZgZ(gO1)*CpUhhbargZgZ(gO2));
   result += 4*B0(p,MVWm,MVWm)*Conj(CpUhhconjVWmVWm(gO2))*CpUhhconjVWmVWm(gO1);
   result += 4*A0(MVWm)*CpUhhUhhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUhhUhhVZVZ(gO1,gO2);
   result += 2*B0(p,MVZ,MVZ)*Conj(CpUhhVZVZ(gO2))*CpUhhVZVZ(gO1);
   std::complex<double> tmp_2295;
   std::complex<double> tmp_2296;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2296 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2295 += tmp_2296;
   result += (-0.5) * tmp_2295;
   std::complex<double> tmp_2297;
   std::complex<double> tmp_2298;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2298 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2297 += tmp_2298;
   result += (-1) * tmp_2297;
   std::complex<double> tmp_2299;
   std::complex<double> tmp_2300;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2300 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2299 += tmp_2300;
   result += (-0.5) * tmp_2299;
   std::complex<double> tmp_2301;
   std::complex<double> tmp_2302;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2303;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2303 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_2302 += tmp_2303;
   }
   tmp_2301 += tmp_2302;
   result += (0.5) * tmp_2301;
   std::complex<double> tmp_2304;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2305;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2305 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2304 += tmp_2305;
   }
   result += tmp_2304;
   std::complex<double> tmp_2306;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2307;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2307 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))
            *CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_2306 += tmp_2307;
   }
   result += tmp_2306;
   std::complex<double> tmp_2308;
   std::complex<double> tmp_2309;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2310;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2310 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_2309 += tmp_2310;
   }
   tmp_2308 += tmp_2309;
   result += (0.5) * tmp_2308;
   std::complex<double> tmp_2311;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2312;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2312 += (Conj(CpUhhbarChaChaPL(gO2,gI1,gI2))*
            CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPR(gO2,gI1,gI2))*
            CpUhhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2311 += tmp_2312;
   }
   result += tmp_2311;
   std::complex<double> tmp_2313;
   std::complex<double> tmp_2314;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2315;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2315 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2314 += tmp_2315;
   }
   tmp_2313 += tmp_2314;
   result += (3) * tmp_2313;
   std::complex<double> tmp_2316;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2317;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2317 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2316 += tmp_2317;
   }
   result += tmp_2316;
   std::complex<double> tmp_2318;
   std::complex<double> tmp_2319;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2320;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2320 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2319 += tmp_2320;
   }
   tmp_2318 += tmp_2319;
   result += (3) * tmp_2318;
   std::complex<double> tmp_2321;
   std::complex<double> tmp_2322;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2323;
      std::complex<double> tmp_2324;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2324 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUhhbarChaChaPR(gO2
            ,gI1,gI2))*CpUhhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUhhbarChaChaPL(gO2,
            gI1,gI2))*CpUhhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2323 += tmp_2324;
      tmp_2322 += (MCha(gI1)) * tmp_2323;
   }
   tmp_2321 += tmp_2322;
   result += (-2) * tmp_2321;
   std::complex<double> tmp_2325;
   std::complex<double> tmp_2326;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2327;
      std::complex<double> tmp_2328;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2328 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2327 += tmp_2328;
      tmp_2326 += (MFd(gI1)) * tmp_2327;
   }
   tmp_2325 += tmp_2326;
   result += (-6) * tmp_2325;
   std::complex<double> tmp_2329;
   std::complex<double> tmp_2330;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2331;
      std::complex<double> tmp_2332;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2332 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2331 += tmp_2332;
      tmp_2330 += (MFe(gI1)) * tmp_2331;
   }
   tmp_2329 += tmp_2330;
   result += (-2) * tmp_2329;
   std::complex<double> tmp_2333;
   std::complex<double> tmp_2334;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2335;
      std::complex<double> tmp_2336;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2336 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2335 += tmp_2336;
      tmp_2334 += (MFu(gI1)) * tmp_2335;
   }
   tmp_2333 += tmp_2334;
   result += (-6) * tmp_2333;
   std::complex<double> tmp_2337;
   std::complex<double> tmp_2338;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2338 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2337 += tmp_2338;
   result += (-1) * tmp_2337;
   std::complex<double> tmp_2339;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2340;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2340 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2339 += tmp_2340;
   }
   result += tmp_2339;
   std::complex<double> tmp_2341;
   std::complex<double> tmp_2342;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2343;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2343 += (Conj(CpUhhChiChiPL(gO2,gI1,gI2))*CpUhhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUhhChiChiPR(gO2,gI1,gI2))*CpUhhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2342 += tmp_2343;
   }
   tmp_2341 += tmp_2342;
   result += (0.5) * tmp_2341;
   std::complex<double> tmp_2344;
   std::complex<double> tmp_2345;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2346;
      std::complex<double> tmp_2347;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2347 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhChiChiPR(gO2,
            gI1,gI2))*CpUhhChiChiPL(gO1,gI1,gI2) + Conj(CpUhhChiChiPL(gO2,gI1,gI2)
            )*CpUhhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2346 += tmp_2347;
      tmp_2345 += (MChi(gI1)) * tmp_2346;
   }
   tmp_2344 += tmp_2345;
   result += (-1) * tmp_2344;
   std::complex<double> tmp_2348;
   std::complex<double> tmp_2349;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2349 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2348 += tmp_2349;
   result += (-3) * tmp_2348;
   std::complex<double> tmp_2350;
   std::complex<double> tmp_2351;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2351 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2350 += tmp_2351;
   result += (-1) * tmp_2350;
   std::complex<double> tmp_2352;
   std::complex<double> tmp_2353;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2353 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2352 += tmp_2353;
   result += (-3) * tmp_2352;
   std::complex<double> tmp_2354;
   std::complex<double> tmp_2355;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2356;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2356 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2355 += tmp_2356;
   }
   tmp_2354 += tmp_2355;
   result += (3) * tmp_2354;
   std::complex<double> tmp_2357;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2358;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2358 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2357 += tmp_2358;
   }
   result += tmp_2357;
   std::complex<double> tmp_2359;
   std::complex<double> tmp_2360;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2361;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2361 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2360 += tmp_2361;
   }
   tmp_2359 += tmp_2360;
   result += (3) * tmp_2359;
   std::complex<double> tmp_2362;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2362 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_2362;
   std::complex<double> tmp_2363;
   std::complex<double> tmp_2364;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2364 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2363 += tmp_2364;
   result += (2) * tmp_2363;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   result += 4*A0(MVWm)*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUAhUAhVZVZ(gO1,gO2);
   std::complex<double> tmp_2365;
   std::complex<double> tmp_2366;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2366 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2365 += tmp_2366;
   result += (-0.5) * tmp_2365;
   std::complex<double> tmp_2367;
   std::complex<double> tmp_2368;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2368 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2367 += tmp_2368;
   result += (-1) * tmp_2367;
   std::complex<double> tmp_2369;
   std::complex<double> tmp_2370;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2370 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2369 += tmp_2370;
   result += (-0.5) * tmp_2369;
   std::complex<double> tmp_2371;
   std::complex<double> tmp_2372;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2373;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2373 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))
            *CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_2372 += tmp_2373;
   }
   tmp_2371 += tmp_2372;
   result += (0.5) * tmp_2371;
   std::complex<double> tmp_2374;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2375;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2375 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_2374 += tmp_2375;
   }
   result += tmp_2374;
   std::complex<double> tmp_2376;
   std::complex<double> tmp_2377;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2378;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2378 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))
            *CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_2377 += tmp_2378;
   }
   tmp_2376 += tmp_2377;
   result += (0.5) * tmp_2376;
   std::complex<double> tmp_2379;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2380;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2380 += (Conj(CpUAhbarChaChaPL(gO2,gI1,gI2))*
            CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPR(gO2,gI1,gI2))*
            CpUAhbarChaChaPR(gO1,gI1,gI2))*G0(p,MCha(gI1),MCha(gI2));
      }
      tmp_2379 += tmp_2380;
   }
   result += tmp_2379;
   std::complex<double> tmp_2381;
   std::complex<double> tmp_2382;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2383;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2383 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2382 += tmp_2383;
   }
   tmp_2381 += tmp_2382;
   result += (3) * tmp_2381;
   std::complex<double> tmp_2384;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2385;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2385 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2384 += tmp_2385;
   }
   result += tmp_2384;
   std::complex<double> tmp_2386;
   std::complex<double> tmp_2387;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2388;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2388 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2387 += tmp_2388;
   }
   tmp_2386 += tmp_2387;
   result += (3) * tmp_2386;
   std::complex<double> tmp_2389;
   std::complex<double> tmp_2390;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2391;
      std::complex<double> tmp_2392;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2392 += B0(p,MCha(gI1),MCha(gI2))*(Conj(CpUAhbarChaChaPR(gO2
            ,gI1,gI2))*CpUAhbarChaChaPL(gO1,gI1,gI2) + Conj(CpUAhbarChaChaPL(gO2,
            gI1,gI2))*CpUAhbarChaChaPR(gO1,gI1,gI2))*MCha(gI2);
      }
      tmp_2391 += tmp_2392;
      tmp_2390 += (MCha(gI1)) * tmp_2391;
   }
   tmp_2389 += tmp_2390;
   result += (-2) * tmp_2389;
   std::complex<double> tmp_2393;
   std::complex<double> tmp_2394;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2395;
      std::complex<double> tmp_2396;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2396 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2395 += tmp_2396;
      tmp_2394 += (MFd(gI1)) * tmp_2395;
   }
   tmp_2393 += tmp_2394;
   result += (-6) * tmp_2393;
   std::complex<double> tmp_2397;
   std::complex<double> tmp_2398;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2399;
      std::complex<double> tmp_2400;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2400 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2399 += tmp_2400;
      tmp_2398 += (MFe(gI1)) * tmp_2399;
   }
   tmp_2397 += tmp_2398;
   result += (-2) * tmp_2397;
   std::complex<double> tmp_2401;
   std::complex<double> tmp_2402;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2403;
      std::complex<double> tmp_2404;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2404 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2403 += tmp_2404;
      tmp_2402 += (MFu(gI1)) * tmp_2403;
   }
   tmp_2401 += tmp_2402;
   result += (-6) * tmp_2401;
   std::complex<double> tmp_2405;
   std::complex<double> tmp_2406;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2406 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2405 += tmp_2406;
   result += (-1) * tmp_2405;
   std::complex<double> tmp_2407;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2408;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2408 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2407 += tmp_2408;
   }
   result += tmp_2407;
   std::complex<double> tmp_2409;
   std::complex<double> tmp_2410;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2411;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2411 += (Conj(CpUAhChiChiPL(gO2,gI1,gI2))*CpUAhChiChiPL(gO1,
            gI1,gI2) + Conj(CpUAhChiChiPR(gO2,gI1,gI2))*CpUAhChiChiPR(gO1,gI1,gI2)
            )*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2410 += tmp_2411;
   }
   tmp_2409 += tmp_2410;
   result += (0.5) * tmp_2409;
   std::complex<double> tmp_2412;
   std::complex<double> tmp_2413;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2414;
      std::complex<double> tmp_2415;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2415 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhChiChiPR(gO2,
            gI1,gI2))*CpUAhChiChiPL(gO1,gI1,gI2) + Conj(CpUAhChiChiPL(gO2,gI1,gI2)
            )*CpUAhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2414 += tmp_2415;
      tmp_2413 += (MChi(gI1)) * tmp_2414;
   }
   tmp_2412 += tmp_2413;
   result += (-1) * tmp_2412;
   std::complex<double> tmp_2416;
   std::complex<double> tmp_2417;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2417 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2416 += tmp_2417;
   result += (-3) * tmp_2416;
   std::complex<double> tmp_2418;
   std::complex<double> tmp_2419;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2419 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2418 += tmp_2419;
   result += (-1) * tmp_2418;
   std::complex<double> tmp_2420;
   std::complex<double> tmp_2421;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2421 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2420 += tmp_2421;
   result += (-3) * tmp_2420;
   std::complex<double> tmp_2422;
   std::complex<double> tmp_2423;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2424;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2424 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2423 += tmp_2424;
   }
   tmp_2422 += tmp_2423;
   result += (3) * tmp_2422;
   std::complex<double> tmp_2425;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2426;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2426 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2425 += tmp_2426;
   }
   result += tmp_2425;
   std::complex<double> tmp_2427;
   std::complex<double> tmp_2428;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2429;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2429 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2428 += tmp_2429;
   }
   tmp_2427 += tmp_2428;
   result += (3) * tmp_2427;
   std::complex<double> tmp_2430;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2430 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_2430;
   std::complex<double> tmp_2431;
   std::complex<double> tmp_2432;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2432 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2431 += tmp_2432;
   result += (2) * tmp_2431;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Hpm(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*B0(p,0,MVWm)*Conj(CpconjUHpmVWmVP(gO2))*CpconjUHpmVWmVP(gO1);
   result += 4*B0(p,MVWm,MVZ)*Conj(CpconjUHpmVZVWm(gO2))*CpconjUHpmVZVWm(gO1);
   result += 4*A0(MVWm)*CpUHpmconjUHpmconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUHpmconjUHpmVZVZ(gO1,gO2);
   result += -(B0(p,MVZ,MVWm)*CpconjUHpmbargWmCgZ(gO1)*CpUHpmgWmCbargZ(gO2));
   result += -(B0(p,MVWm,MVZ)*CpconjUHpmbargZgWm(gO1)*CpUHpmgZbargWm(gO2));
   std::complex<double> tmp_2433;
   std::complex<double> tmp_2434;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2434 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2433 += tmp_2434;
   result += (-0.5) * tmp_2433;
   std::complex<double> tmp_2435;
   std::complex<double> tmp_2436;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2436 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2435 += tmp_2436;
   result += (-1) * tmp_2435;
   std::complex<double> tmp_2437;
   std::complex<double> tmp_2438;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2438 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2437 += tmp_2438;
   result += (-0.5) * tmp_2437;
   std::complex<double> tmp_2439;
   std::complex<double> tmp_2440;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2441;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2441 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_2440 += tmp_2441;
   }
   tmp_2439 += tmp_2440;
   result += (3) * tmp_2439;
   std::complex<double> tmp_2442;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2443;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2443 += (Conj(CpconjUHpmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFvFePR(gO2,gI1,
            gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_2442 += tmp_2443;
   }
   result += tmp_2442;
   std::complex<double> tmp_2444;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2445;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2445 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2444 += tmp_2445;
   }
   result += tmp_2444;
   std::complex<double> tmp_2446;
   std::complex<double> tmp_2447;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2448;
      std::complex<double> tmp_2449;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2449 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2448 += tmp_2449;
      tmp_2447 += (MFu(gI1)) * tmp_2448;
   }
   tmp_2446 += tmp_2447;
   result += (-6) * tmp_2446;
   std::complex<double> tmp_2450;
   std::complex<double> tmp_2451;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2452;
      std::complex<double> tmp_2453;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2453 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmbarFvFePR(
            gO2,gI1,gI2))*CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFvFePL(gO2,gI1,gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2452 += tmp_2453;
      tmp_2451 += (MFv(gI1)) * tmp_2452;
   }
   tmp_2450 += tmp_2451;
   result += (-2) * tmp_2450;
   std::complex<double> tmp_2454;
   std::complex<double> tmp_2455;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2455 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2454 += tmp_2455;
   result += (-1) * tmp_2454;
   std::complex<double> tmp_2456;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2457;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2457 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2456 += tmp_2457;
   }
   result += tmp_2456;
   std::complex<double> tmp_2458;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2459;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2459 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2458 += tmp_2459;
   }
   result += tmp_2458;
   std::complex<double> tmp_2460;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2461;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2461 += (Conj(CpconjUHpmChiChaPL(gO2,gI1,gI2))*
            CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(CpconjUHpmChiChaPR(gO2,gI1,gI2)
            )*CpconjUHpmChiChaPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha(gI2));
      }
      tmp_2460 += tmp_2461;
   }
   result += tmp_2460;
   std::complex<double> tmp_2462;
   std::complex<double> tmp_2463;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2464;
      std::complex<double> tmp_2465;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2465 += B0(p,MChi(gI1),MCha(gI2))*(Conj(CpconjUHpmChiChaPR(
            gO2,gI1,gI2))*CpconjUHpmChiChaPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmChiChaPL(gO2,gI1,gI2))*CpconjUHpmChiChaPR(gO1,gI1,gI2))*MCha
            (gI2);
      }
      tmp_2464 += tmp_2465;
      tmp_2463 += (MChi(gI1)) * tmp_2464;
   }
   tmp_2462 += tmp_2463;
   result += (-2) * tmp_2462;
   std::complex<double> tmp_2466;
   std::complex<double> tmp_2467;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2467 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2466 += tmp_2467;
   result += (-3) * tmp_2466;
   std::complex<double> tmp_2468;
   std::complex<double> tmp_2469;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2469 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2468 += tmp_2469;
   result += (-1) * tmp_2468;
   std::complex<double> tmp_2470;
   std::complex<double> tmp_2471;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2471 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2470 += tmp_2471;
   result += (-3) * tmp_2470;
   std::complex<double> tmp_2472;
   std::complex<double> tmp_2473;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2474;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2474 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2473 += tmp_2474;
   }
   tmp_2472 += tmp_2473;
   result += (3) * tmp_2472;
   std::complex<double> tmp_2475;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2475 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_2475;
   std::complex<double> tmp_2476;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2476 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_2476;
   std::complex<double> tmp_2477;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2477 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_2477;
   std::complex<double> tmp_2478;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2478 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_2478;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_2479;
   std::complex<double> tmp_2480;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2480 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_2479 += tmp_2480;
   result += (0.5) * tmp_2479;
   std::complex<double> tmp_2481;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2481 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2481;
   std::complex<double> tmp_2482;
   std::complex<double> tmp_2483;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2483 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2482 += tmp_2483;
   result += (0.5) * tmp_2482;
   std::complex<double> tmp_2484;
   std::complex<double> tmp_2485;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2486;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2486 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2485 += tmp_2486;
   }
   tmp_2484 += tmp_2485;
   result += (-4) * tmp_2484;
   std::complex<double> tmp_2487;
   std::complex<double> tmp_2488;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2489;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2489 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_2488 += tmp_2489;
   }
   tmp_2487 += tmp_2488;
   result += (-4) * tmp_2487;
   std::complex<double> tmp_2490;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2491;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2491 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2491 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_2490 += tmp_2491;
   }
   result += tmp_2490;
   std::complex<double> tmp_2492;
   std::complex<double> tmp_2493;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2494;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2494 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2494 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_2493 += tmp_2494;
   }
   tmp_2492 += tmp_2493;
   result += (3) * tmp_2492;
   std::complex<double> tmp_2495;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2496;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2496 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2496 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_2495 += tmp_2496;
   }
   result += tmp_2495;
   std::complex<double> tmp_2497;
   std::complex<double> tmp_2498;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2499;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2499 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2499 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_2498 += tmp_2499;
   }
   tmp_2497 += tmp_2498;
   result += (3) * tmp_2497;
   std::complex<double> tmp_2500;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2501;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2501 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_2501 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_2500 += tmp_2501;
   }
   result += tmp_2500;
   std::complex<double> tmp_2502;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2502 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_2502;
   std::complex<double> tmp_2503;
   std::complex<double> tmp_2504;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2505;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2505 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2504 += tmp_2505;
   }
   tmp_2503 += tmp_2504;
   result += (-4) * tmp_2503;
   std::complex<double> tmp_2506;
   std::complex<double> tmp_2507;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2508;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2508 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2508 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_2507 += tmp_2508;
   }
   tmp_2506 += tmp_2507;
   result += (0.5) * tmp_2506;
   std::complex<double> tmp_2509;
   std::complex<double> tmp_2510;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2510 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2509 += tmp_2510;
   result += (3) * tmp_2509;
   std::complex<double> tmp_2511;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2511 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2511;
   std::complex<double> tmp_2512;
   std::complex<double> tmp_2513;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2513 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2512 += tmp_2513;
   result += (3) * tmp_2512;
   std::complex<double> tmp_2514;
   std::complex<double> tmp_2515;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2516;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2516 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2515 += tmp_2516;
   }
   tmp_2514 += tmp_2515;
   result += (-12) * tmp_2514;
   std::complex<double> tmp_2517;
   std::complex<double> tmp_2518;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2519;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2519 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2518 += tmp_2519;
   }
   tmp_2517 += tmp_2518;
   result += (-4) * tmp_2517;
   std::complex<double> tmp_2520;
   std::complex<double> tmp_2521;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2522;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2522 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2521 += tmp_2522;
   }
   tmp_2520 += tmp_2521;
   result += (-12) * tmp_2520;
   std::complex<double> tmp_2523;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2523 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2523;
   std::complex<double> tmp_2524;
   std::complex<double> tmp_2525;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2525 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2524 += tmp_2525;
   result += (2) * tmp_2524;
   result += -(AbsSqr(CpVZconjVWmVWm())*(2*A0(MVWm) + 10*B00(p,MVWm,MVWm) + B0(
      p,MVWm,MVWm)*(2*Sqr(MVWm) + 4*Sqr(p))));

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
   result += -0.5*A0(MVZ)*(4*CpVWmconjVWmVZVZ1() + CpVWmconjVWmVZVZ2() +
      CpVWmconjVWmVZVZ3());
   std::complex<double> tmp_2526;
   std::complex<double> tmp_2527;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2527 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_2526 += tmp_2527;
   result += (0.5) * tmp_2526;
   std::complex<double> tmp_2528;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2528 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_2528;
   std::complex<double> tmp_2529;
   std::complex<double> tmp_2530;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2530 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_2529 += tmp_2530;
   result += (0.5) * tmp_2529;
   std::complex<double> tmp_2531;
   std::complex<double> tmp_2532;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2533;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2533 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_2533 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_2532 += tmp_2533;
   }
   tmp_2531 += tmp_2532;
   result += (3) * tmp_2531;
   std::complex<double> tmp_2534;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2535;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2535 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_2535 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_2534 += tmp_2535;
   }
   result += tmp_2534;
   std::complex<double> tmp_2536;
   std::complex<double> tmp_2537;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2538;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2538 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_2537 += tmp_2538;
   }
   tmp_2536 += tmp_2537;
   result += (-4) * tmp_2536;
   std::complex<double> tmp_2539;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2539 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_2539;
   std::complex<double> tmp_2540;
   std::complex<double> tmp_2541;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2542;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2542 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_2541 += tmp_2542;
   }
   tmp_2540 += tmp_2541;
   result += (-4) * tmp_2540;
   std::complex<double> tmp_2543;
   std::complex<double> tmp_2544;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2545;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2545 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_2544 += tmp_2545;
   }
   tmp_2543 += tmp_2544;
   result += (-4) * tmp_2543;
   std::complex<double> tmp_2546;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2547;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2547 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_2547 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_2546 += tmp_2547;
   }
   result += tmp_2546;
   std::complex<double> tmp_2548;
   std::complex<double> tmp_2549;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2549 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_2548 += tmp_2549;
   result += (3) * tmp_2548;
   std::complex<double> tmp_2550;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2550 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_2550;
   std::complex<double> tmp_2551;
   std::complex<double> tmp_2552;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2552 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_2551 += tmp_2552;
   result += (3) * tmp_2551;
   std::complex<double> tmp_2553;
   std::complex<double> tmp_2554;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2555;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2555 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_2554 += tmp_2555;
   }
   tmp_2553 += tmp_2554;
   result += (-12) * tmp_2553;
   std::complex<double> tmp_2556;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2556 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_2556;
   std::complex<double> tmp_2557;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2557 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_2557;
   std::complex<double> tmp_2558;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2558 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_2558;
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) + B0(p,
      MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) + B0(p,MVZ,MVWm)*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2559;
   std::complex<double> tmp_2560;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2561;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2561 += B0(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2560 += tmp_2561;
   }
   tmp_2559 += tmp_2560;
   result += (2) * tmp_2559;
   std::complex<double> tmp_2562;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2563;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2563 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2562 += tmp_2563;
   }
   result += tmp_2562;
   std::complex<double> tmp_2564;
   std::complex<double> tmp_2565;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2566;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2566 += B0(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2565 += tmp_2566;
   }
   tmp_2564 += tmp_2565;
   result += (2) * tmp_2564;
   std::complex<double> tmp_2567;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2568;
      std::complex<double> tmp_2569;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2569 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2568 += tmp_2569;
      tmp_2567 += (MChi(gI1)) * tmp_2568;
   }
   result += tmp_2567;
   std::complex<double> tmp_2570;
   std::complex<double> tmp_2571;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2572;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2572 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2571 += tmp_2572;
   }
   tmp_2570 += tmp_2571;
   result += (6) * tmp_2570;
   std::complex<double> tmp_2573;
   std::complex<double> tmp_2574;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2575;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2575 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2574 += tmp_2575;
   }
   tmp_2573 += tmp_2574;
   result += (2) * tmp_2573;
   std::complex<double> tmp_2576;
   std::complex<double> tmp_2577;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2578;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2578 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2577 += tmp_2578;
   }
   tmp_2576 += tmp_2577;
   result += (6) * tmp_2576;
   std::complex<double> tmp_2579;
   std::complex<double> tmp_2580;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2580 += B0(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2579 += tmp_2580;
   result += (-8) * tmp_2579;
   std::complex<double> tmp_2581;
   std::complex<double> tmp_2582;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2582 += B0(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2581 += tmp_2582;
   result += (-4) * tmp_2581;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2583;
   std::complex<double> tmp_2584;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2585;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2585 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPR(gO2,
            gI1,gI2))*CpUChiconjSvFvPR(gO1,gI1,gI2);
      }
      tmp_2584 += tmp_2585;
   }
   tmp_2583 += tmp_2584;
   result += (-1) * tmp_2583;
   std::complex<double> tmp_2586;
   std::complex<double> tmp_2587;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2588;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2588 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPR(gO2,gI1,
            gI2))*CpUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_2587 += tmp_2588;
   }
   tmp_2586 += tmp_2587;
   result += (-0.5) * tmp_2586;
   std::complex<double> tmp_2589;
   std::complex<double> tmp_2590;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2591;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2591 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPR(
            gO2,gI1,gI2))*CpUChiconjHpmChaPR(gO1,gI1,gI2);
      }
      tmp_2590 += tmp_2591;
   }
   tmp_2589 += tmp_2590;
   result += (-1) * tmp_2589;
   std::complex<double> tmp_2592;
   std::complex<double> tmp_2593;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2594;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2594 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPR(gO2,gI1,
            gI2))*CpUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2593 += tmp_2594;
   }
   tmp_2592 += tmp_2593;
   result += (-0.5) * tmp_2592;
   std::complex<double> tmp_2595;
   std::complex<double> tmp_2596;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2597;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2597 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPR(gO2,
            gI1,gI2))*CpUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_2596 += tmp_2597;
   }
   tmp_2595 += tmp_2596;
   result += (-3) * tmp_2595;
   std::complex<double> tmp_2598;
   std::complex<double> tmp_2599;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2600;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2600 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePR(gO2,
            gI1,gI2))*CpUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_2599 += tmp_2600;
   }
   tmp_2598 += tmp_2599;
   result += (-1) * tmp_2598;
   std::complex<double> tmp_2601;
   std::complex<double> tmp_2602;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2603;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2603 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPR(gO2,
            gI1,gI2))*CpUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_2602 += tmp_2603;
   }
   tmp_2601 += tmp_2602;
   result += (-3) * tmp_2601;
   std::complex<double> tmp_2604;
   std::complex<double> tmp_2605;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2605 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPL(gO2,gI2))*
         CpUChiconjVWmChaPL(gO1,gI2);
   }
   tmp_2604 += tmp_2605;
   result += (-2) * tmp_2604;
   std::complex<double> tmp_2606;
   std::complex<double> tmp_2607;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2607 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPL(gO2,gI2))*
         CpUChiVZChiPL(gO1,gI2);
   }
   tmp_2606 += tmp_2607;
   result += (-1) * tmp_2606;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2608;
   std::complex<double> tmp_2609;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2610;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2610 += B1(p,MFv(gI2),MSv(gI1))*Conj(CpUChiconjSvFvPL(gO2,
            gI1,gI2))*CpUChiconjSvFvPL(gO1,gI1,gI2);
      }
      tmp_2609 += tmp_2610;
   }
   tmp_2608 += tmp_2609;
   result += (-1) * tmp_2608;
   std::complex<double> tmp_2611;
   std::complex<double> tmp_2612;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2613;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2613 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpUChihhChiPL(gO2,gI1,
            gI2))*CpUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_2612 += tmp_2613;
   }
   tmp_2611 += tmp_2612;
   result += (-0.5) * tmp_2611;
   std::complex<double> tmp_2614;
   std::complex<double> tmp_2615;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2616;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2616 += B1(p,MCha(gI2),MHpm(gI1))*Conj(CpUChiconjHpmChaPL(
            gO2,gI1,gI2))*CpUChiconjHpmChaPL(gO1,gI1,gI2);
      }
      tmp_2615 += tmp_2616;
   }
   tmp_2614 += tmp_2615;
   result += (-1) * tmp_2614;
   std::complex<double> tmp_2617;
   std::complex<double> tmp_2618;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2619;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2619 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpUChiChiAhPL(gO2,gI1,
            gI2))*CpUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_2618 += tmp_2619;
   }
   tmp_2617 += tmp_2618;
   result += (-0.5) * tmp_2617;
   std::complex<double> tmp_2620;
   std::complex<double> tmp_2621;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2622;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2622 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpUChiconjSdFdPL(gO2,
            gI1,gI2))*CpUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_2621 += tmp_2622;
   }
   tmp_2620 += tmp_2621;
   result += (-3) * tmp_2620;
   std::complex<double> tmp_2623;
   std::complex<double> tmp_2624;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2625;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2625 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpUChiconjSeFePL(gO2,
            gI1,gI2))*CpUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_2624 += tmp_2625;
   }
   tmp_2623 += tmp_2624;
   result += (-1) * tmp_2623;
   std::complex<double> tmp_2626;
   std::complex<double> tmp_2627;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2628;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2628 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpUChiconjSuFuPL(gO2,
            gI1,gI2))*CpUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_2627 += tmp_2628;
   }
   tmp_2626 += tmp_2627;
   result += (-3) * tmp_2626;
   std::complex<double> tmp_2629;
   std::complex<double> tmp_2630;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2630 += B1(p,MCha(gI2),MVWm)*Conj(CpUChiconjVWmChaPR(gO2,gI2))*
         CpUChiconjVWmChaPR(gO1,gI2);
   }
   tmp_2629 += tmp_2630;
   result += (-2) * tmp_2629;
   std::complex<double> tmp_2631;
   std::complex<double> tmp_2632;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2632 += B1(p,MChi(gI2),MVZ)*Conj(CpUChiVZChiPR(gO2,gI2))*
         CpUChiVZChiPR(gO1,gI2);
   }
   tmp_2631 += tmp_2632;
   result += (-1) * tmp_2631;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2633;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2634;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2634 += B0(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2633 += tmp_2634;
   }
   result += tmp_2633;
   std::complex<double> tmp_2635;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2636;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2636 += B0(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2635 += tmp_2636;
   }
   result += tmp_2635;
   std::complex<double> tmp_2637;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2638;
      std::complex<double> tmp_2639;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2639 += B0(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_2638 += tmp_2639;
      tmp_2637 += (MCha(gI1)) * tmp_2638;
   }
   result += tmp_2637;
   std::complex<double> tmp_2640;
   std::complex<double> tmp_2641;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2642;
      std::complex<double> tmp_2643;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2643 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_2642 += tmp_2643;
      tmp_2641 += (MFu(gI1)) * tmp_2642;
   }
   tmp_2640 += tmp_2641;
   result += (3) * tmp_2640;
   std::complex<double> tmp_2644;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2645;
      std::complex<double> tmp_2646;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2646 += B0(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_2645 += tmp_2646;
      tmp_2644 += (MFv(gI1)) * tmp_2645;
   }
   result += tmp_2644;
   std::complex<double> tmp_2647;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2648;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2648 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2647 += tmp_2648;
   }
   result += tmp_2647;
   std::complex<double> tmp_2649;
   std::complex<double> tmp_2650;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2651;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2651 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2650 += tmp_2651;
   }
   tmp_2649 += tmp_2650;
   result += (3) * tmp_2649;
   std::complex<double> tmp_2652;
   std::complex<double> tmp_2653;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2653 += B0(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2652 += tmp_2653;
   result += (-4) * tmp_2652;
   std::complex<double> tmp_2654;
   std::complex<double> tmp_2655;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2655 += B0(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2)*MCha(gI2);
   }
   tmp_2654 += tmp_2655;
   result += (-4) * tmp_2654;
   std::complex<double> tmp_2656;
   std::complex<double> tmp_2657;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2657 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2656 += tmp_2657;
   result += (-4) * tmp_2656;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2658;
   std::complex<double> tmp_2659;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2660;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2660 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPR(gO2,
            gI1,gI2))*CpbarUChaChaAhPR(gO1,gI1,gI2);
      }
      tmp_2659 += tmp_2660;
   }
   tmp_2658 += tmp_2659;
   result += (-0.5) * tmp_2658;
   std::complex<double> tmp_2661;
   std::complex<double> tmp_2662;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2663;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2663 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePR(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePR(gO1,gI1,gI2);
      }
      tmp_2662 += tmp_2663;
   }
   tmp_2661 += tmp_2662;
   result += (-0.5) * tmp_2661;
   std::complex<double> tmp_2664;
   std::complex<double> tmp_2665;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2666;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2666 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPR(gO2,
            gI1,gI2))*CpbarUChahhChaPR(gO1,gI1,gI2);
      }
      tmp_2665 += tmp_2666;
   }
   tmp_2664 += tmp_2665;
   result += (-0.5) * tmp_2664;
   std::complex<double> tmp_2667;
   std::complex<double> tmp_2668;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2669;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2669 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPR(gO2,
            gI1,gI2))*CpbarUChabarFuSdPR(gO1,gI1,gI2);
      }
      tmp_2668 += tmp_2669;
   }
   tmp_2667 += tmp_2668;
   result += (-1.5) * tmp_2667;
   std::complex<double> tmp_2670;
   std::complex<double> tmp_2671;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2672;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2672 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePR(gO2,
            gI1,gI2))*CpbarUChabarFvSePR(gO1,gI1,gI2);
      }
      tmp_2671 += tmp_2672;
   }
   tmp_2670 += tmp_2671;
   result += (-0.5) * tmp_2670;
   std::complex<double> tmp_2673;
   std::complex<double> tmp_2674;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2675;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2675 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPR(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPR(gO1,gI1,gI2);
      }
      tmp_2674 += tmp_2675;
   }
   tmp_2673 += tmp_2674;
   result += (-0.5) * tmp_2673;
   std::complex<double> tmp_2676;
   std::complex<double> tmp_2677;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2678;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2678 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPR(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPR(gO1,gI1,gI2);
      }
      tmp_2677 += tmp_2678;
   }
   tmp_2676 += tmp_2677;
   result += (-1.5) * tmp_2676;
   std::complex<double> tmp_2679;
   std::complex<double> tmp_2680;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2680 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPL(gO2,gI2))*
         CpbarUChaVPChaPL(gO1,gI2);
   }
   tmp_2679 += tmp_2680;
   result += (-1) * tmp_2679;
   std::complex<double> tmp_2681;
   std::complex<double> tmp_2682;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2682 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPL(gO2,gI2))*
         CpbarUChaVZChaPL(gO1,gI2);
   }
   tmp_2681 += tmp_2682;
   result += (-1) * tmp_2681;
   std::complex<double> tmp_2683;
   std::complex<double> tmp_2684;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2684 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPL(gO2,gI2))*
         CpbarUChaVWmChiPL(gO1,gI2);
   }
   tmp_2683 += tmp_2684;
   result += (-1) * tmp_2683;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2685;
   std::complex<double> tmp_2686;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2687;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2687 += B1(p,MCha(gI1),MAh(gI2))*Conj(CpbarUChaChaAhPL(gO2,
            gI1,gI2))*CpbarUChaChaAhPL(gO1,gI1,gI2);
      }
      tmp_2686 += tmp_2687;
   }
   tmp_2685 += tmp_2686;
   result += (-0.5) * tmp_2685;
   std::complex<double> tmp_2688;
   std::complex<double> tmp_2689;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2690;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2690 += B1(p,MFe(gI2),MSv(gI1))*Conj(CpbarUChaconjSvFePL(gO2
            ,gI1,gI2))*CpbarUChaconjSvFePL(gO1,gI1,gI2);
      }
      tmp_2689 += tmp_2690;
   }
   tmp_2688 += tmp_2689;
   result += (-0.5) * tmp_2688;
   std::complex<double> tmp_2691;
   std::complex<double> tmp_2692;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2693;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2693 += B1(p,MCha(gI2),Mhh(gI1))*Conj(CpbarUChahhChaPL(gO2,
            gI1,gI2))*CpbarUChahhChaPL(gO1,gI1,gI2);
      }
      tmp_2692 += tmp_2693;
   }
   tmp_2691 += tmp_2692;
   result += (-0.5) * tmp_2691;
   std::complex<double> tmp_2694;
   std::complex<double> tmp_2695;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2696;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2696 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUChabarFuSdPL(gO2,
            gI1,gI2))*CpbarUChabarFuSdPL(gO1,gI1,gI2);
      }
      tmp_2695 += tmp_2696;
   }
   tmp_2694 += tmp_2695;
   result += (-1.5) * tmp_2694;
   std::complex<double> tmp_2697;
   std::complex<double> tmp_2698;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2699;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2699 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUChabarFvSePL(gO2,
            gI1,gI2))*CpbarUChabarFvSePL(gO1,gI1,gI2);
      }
      tmp_2698 += tmp_2699;
   }
   tmp_2697 += tmp_2698;
   result += (-0.5) * tmp_2697;
   std::complex<double> tmp_2700;
   std::complex<double> tmp_2701;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2702;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2702 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUChaHpmChiPL(gO2
            ,gI1,gI2))*CpbarUChaHpmChiPL(gO1,gI1,gI2);
      }
      tmp_2701 += tmp_2702;
   }
   tmp_2700 += tmp_2701;
   result += (-0.5) * tmp_2700;
   std::complex<double> tmp_2703;
   std::complex<double> tmp_2704;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2705;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2705 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUChaconjSuFdPL(gO2
            ,gI1,gI2))*CpbarUChaconjSuFdPL(gO1,gI1,gI2);
      }
      tmp_2704 += tmp_2705;
   }
   tmp_2703 += tmp_2704;
   result += (-1.5) * tmp_2703;
   std::complex<double> tmp_2706;
   std::complex<double> tmp_2707;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2707 += B1(p,MCha(gI2),0)*Conj(CpbarUChaVPChaPR(gO2,gI2))*
         CpbarUChaVPChaPR(gO1,gI2);
   }
   tmp_2706 += tmp_2707;
   result += (-1) * tmp_2706;
   std::complex<double> tmp_2708;
   std::complex<double> tmp_2709;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2709 += B1(p,MCha(gI2),MVZ)*Conj(CpbarUChaVZChaPR(gO2,gI2))*
         CpbarUChaVZChaPR(gO1,gI2);
   }
   tmp_2708 += tmp_2709;
   result += (-1) * tmp_2708;
   std::complex<double> tmp_2710;
   std::complex<double> tmp_2711;
   for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
      tmp_2711 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUChaVWmChiPR(gO2,gI2))*
         CpbarUChaVWmChiPR(gO1,gI2);
   }
   tmp_2710 += tmp_2711;
   result += (-1) * tmp_2710;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2712;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2713;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2713 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2712 += tmp_2713;
   }
   result += tmp_2712;
   std::complex<double> tmp_2714;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2715;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2715 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2714 += tmp_2715;
   }
   result += tmp_2714;
   std::complex<double> tmp_2716;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2717;
      std::complex<double> tmp_2718;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2718 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_2717 += tmp_2718;
      tmp_2716 += (MFe(gI1)) * tmp_2717;
   }
   result += tmp_2716;
   std::complex<double> tmp_2719;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2720;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2720 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2719 += tmp_2720;
   }
   result += tmp_2719;
   std::complex<double> tmp_2721;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2722;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2722 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2721 += tmp_2722;
   }
   result += tmp_2721;
   std::complex<double> tmp_2723;
   std::complex<double> tmp_2724;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2724 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2723 += tmp_2724;
   result += (-4) * tmp_2723;
   std::complex<double> tmp_2725;
   std::complex<double> tmp_2726;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2726 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2725 += tmp_2726;
   result += (-4) * tmp_2725;
   std::complex<double> tmp_2727;
   std::complex<double> tmp_2728;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2728 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_2727 += tmp_2728;
   result += (-4) * tmp_2727;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2729;
   std::complex<double> tmp_2730;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2731;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2731 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_2730 += tmp_2731;
   }
   tmp_2729 += tmp_2730;
   result += (-0.5) * tmp_2729;
   std::complex<double> tmp_2732;
   std::complex<double> tmp_2733;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2734;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2734 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_2733 += tmp_2734;
   }
   tmp_2732 += tmp_2733;
   result += (-0.5) * tmp_2732;
   std::complex<double> tmp_2735;
   std::complex<double> tmp_2736;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2737;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2737 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPR(gO2,
            gI1,gI2))*CpbarUFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_2736 += tmp_2737;
   }
   tmp_2735 += tmp_2736;
   result += (-0.5) * tmp_2735;
   std::complex<double> tmp_2738;
   std::complex<double> tmp_2739;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2740;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2740 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_2739 += tmp_2740;
   }
   tmp_2738 += tmp_2739;
   result += (-0.5) * tmp_2738;
   std::complex<double> tmp_2741;
   std::complex<double> tmp_2742;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2743;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2743 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_2742 += tmp_2743;
   }
   tmp_2741 += tmp_2742;
   result += (-0.5) * tmp_2741;
   std::complex<double> tmp_2744;
   std::complex<double> tmp_2745;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2745 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_2744 += tmp_2745;
   result += (-1) * tmp_2744;
   std::complex<double> tmp_2746;
   std::complex<double> tmp_2747;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2747 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_2746 += tmp_2747;
   result += (-1) * tmp_2746;
   std::complex<double> tmp_2748;
   std::complex<double> tmp_2749;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2749 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_2748 += tmp_2749;
   result += (-1) * tmp_2748;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2750;
   std::complex<double> tmp_2751;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2752;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2752 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_2751 += tmp_2752;
   }
   tmp_2750 += tmp_2751;
   result += (-0.5) * tmp_2750;
   std::complex<double> tmp_2753;
   std::complex<double> tmp_2754;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2755;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2755 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_2754 += tmp_2755;
   }
   tmp_2753 += tmp_2754;
   result += (-0.5) * tmp_2753;
   std::complex<double> tmp_2756;
   std::complex<double> tmp_2757;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2758;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2758 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarUFeSvChaPL(gO2,
            gI1,gI2))*CpbarUFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_2757 += tmp_2758;
   }
   tmp_2756 += tmp_2757;
   result += (-0.5) * tmp_2756;
   std::complex<double> tmp_2759;
   std::complex<double> tmp_2760;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2761;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2761 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_2760 += tmp_2761;
   }
   tmp_2759 += tmp_2760;
   result += (-0.5) * tmp_2759;
   std::complex<double> tmp_2762;
   std::complex<double> tmp_2763;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2764;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2764 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_2763 += tmp_2764;
   }
   tmp_2762 += tmp_2763;
   result += (-0.5) * tmp_2762;
   std::complex<double> tmp_2765;
   std::complex<double> tmp_2766;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2766 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_2765 += tmp_2766;
   result += (-1) * tmp_2765;
   std::complex<double> tmp_2767;
   std::complex<double> tmp_2768;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2768 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_2767 += tmp_2768;
   result += (-1) * tmp_2767;
   std::complex<double> tmp_2769;
   std::complex<double> tmp_2770;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2770 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_2769 += tmp_2770;
   result += (-1) * tmp_2769;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2771;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2772;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2772 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2771 += tmp_2772;
   }
   result += tmp_2771;
   std::complex<double> tmp_2773;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2774;
      std::complex<double> tmp_2775;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2775 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_2774 += tmp_2775;
      tmp_2773 += (MFd(gI1)) * tmp_2774;
   }
   result += tmp_2773;
   std::complex<double> tmp_2776;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2777;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2777 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2776 += tmp_2777;
   }
   result += tmp_2776;
   std::complex<double> tmp_2778;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2779;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2779 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2778 += tmp_2779;
   }
   result += tmp_2778;
   std::complex<double> tmp_2780;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2781;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2781 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2780 += tmp_2781;
   }
   result += tmp_2780;
   std::complex<double> tmp_2782;
   std::complex<double> tmp_2783;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2783 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2782 += tmp_2783;
   result += (-5.333333333333333) * tmp_2782;
   std::complex<double> tmp_2784;
   std::complex<double> tmp_2785;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2785 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2784 += tmp_2785;
   result += (-4) * tmp_2784;
   std::complex<double> tmp_2786;
   std::complex<double> tmp_2787;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2787 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2786 += tmp_2787;
   result += (-4) * tmp_2786;
   std::complex<double> tmp_2788;
   std::complex<double> tmp_2789;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2789 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2788 += tmp_2789;
   result += (-4) * tmp_2788;
   std::complex<double> tmp_2790;
   std::complex<double> tmp_2791;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2791 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_2790 += tmp_2791;
   result += (1.3333333333333333*MGlu) * tmp_2790;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2792;
   std::complex<double> tmp_2793;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2794;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2794 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_2793 += tmp_2794;
   }
   tmp_2792 += tmp_2793;
   result += (-0.5) * tmp_2792;
   std::complex<double> tmp_2795;
   std::complex<double> tmp_2796;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2797;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2797 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_2796 += tmp_2797;
   }
   tmp_2795 += tmp_2796;
   result += (-0.5) * tmp_2795;
   std::complex<double> tmp_2798;
   std::complex<double> tmp_2799;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2800;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2800 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_2799 += tmp_2800;
   }
   tmp_2798 += tmp_2799;
   result += (-0.5) * tmp_2798;
   std::complex<double> tmp_2801;
   std::complex<double> tmp_2802;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2802 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_2801 += tmp_2802;
   result += (-0.6666666666666666) * tmp_2801;
   std::complex<double> tmp_2803;
   std::complex<double> tmp_2804;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2805;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2805 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPR(gO2,
            gI1,gI2))*CpbarUFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_2804 += tmp_2805;
   }
   tmp_2803 += tmp_2804;
   result += (-0.5) * tmp_2803;
   std::complex<double> tmp_2806;
   std::complex<double> tmp_2807;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2808;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2808 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_2807 += tmp_2808;
   }
   tmp_2806 += tmp_2807;
   result += (-0.5) * tmp_2806;
   std::complex<double> tmp_2809;
   std::complex<double> tmp_2810;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2810 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_2809 += tmp_2810;
   result += (-1.3333333333333333) * tmp_2809;
   std::complex<double> tmp_2811;
   std::complex<double> tmp_2812;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2812 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_2811 += tmp_2812;
   result += (-1) * tmp_2811;
   std::complex<double> tmp_2813;
   std::complex<double> tmp_2814;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2814 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_2813 += tmp_2814;
   result += (-1) * tmp_2813;
   std::complex<double> tmp_2815;
   std::complex<double> tmp_2816;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2816 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_2815 += tmp_2816;
   result += (-1) * tmp_2815;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2817;
   std::complex<double> tmp_2818;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2819;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2819 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_2818 += tmp_2819;
   }
   tmp_2817 += tmp_2818;
   result += (-0.5) * tmp_2817;
   std::complex<double> tmp_2820;
   std::complex<double> tmp_2821;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2822;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2822 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_2821 += tmp_2822;
   }
   tmp_2820 += tmp_2821;
   result += (-0.5) * tmp_2820;
   std::complex<double> tmp_2823;
   std::complex<double> tmp_2824;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2825;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2825 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_2824 += tmp_2825;
   }
   tmp_2823 += tmp_2824;
   result += (-0.5) * tmp_2823;
   std::complex<double> tmp_2826;
   std::complex<double> tmp_2827;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2827 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPL(gO1,gI1,1);
   }
   tmp_2826 += tmp_2827;
   result += (-0.6666666666666666) * tmp_2826;
   std::complex<double> tmp_2828;
   std::complex<double> tmp_2829;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2830;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2830 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarUFdSuChaPL(gO2,
            gI1,gI2))*CpbarUFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_2829 += tmp_2830;
   }
   tmp_2828 += tmp_2829;
   result += (-0.5) * tmp_2828;
   std::complex<double> tmp_2831;
   std::complex<double> tmp_2832;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2833;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2833 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_2832 += tmp_2833;
   }
   tmp_2831 += tmp_2832;
   result += (-0.5) * tmp_2831;
   std::complex<double> tmp_2834;
   std::complex<double> tmp_2835;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2835 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_2834 += tmp_2835;
   result += (-1.3333333333333333) * tmp_2834;
   std::complex<double> tmp_2836;
   std::complex<double> tmp_2837;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2837 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_2836 += tmp_2837;
   result += (-1) * tmp_2836;
   std::complex<double> tmp_2838;
   std::complex<double> tmp_2839;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2839 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_2838 += tmp_2839;
   result += (-1) * tmp_2838;
   std::complex<double> tmp_2840;
   std::complex<double> tmp_2841;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2841 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_2840 += tmp_2841;
   result += (-1) * tmp_2840;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2842;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2843;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2843 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2842 += tmp_2843;
   }
   result += tmp_2842;
   std::complex<double> tmp_2844;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2845;
      std::complex<double> tmp_2846;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2846 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_2845 += tmp_2846;
      tmp_2844 += (MCha(gI1)) * tmp_2845;
   }
   result += tmp_2844;
   std::complex<double> tmp_2847;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2848;
      std::complex<double> tmp_2849;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2849 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_2848 += tmp_2849;
      tmp_2847 += (MFu(gI1)) * tmp_2848;
   }
   result += tmp_2847;
   std::complex<double> tmp_2850;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2851;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2851 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2850 += tmp_2851;
   }
   result += tmp_2850;
   std::complex<double> tmp_2852;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2853;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2853 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2852 += tmp_2853;
   }
   result += tmp_2852;
   std::complex<double> tmp_2854;
   std::complex<double> tmp_2855;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2855 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_2854 += tmp_2855;
   result += (-4) * tmp_2854;
   std::complex<double> tmp_2856;
   std::complex<double> tmp_2857;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2857 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2856 += tmp_2857;
   result += (-5.333333333333333) * tmp_2856;
   std::complex<double> tmp_2858;
   std::complex<double> tmp_2859;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2859 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2858 += tmp_2859;
   result += (-4) * tmp_2858;
   std::complex<double> tmp_2860;
   std::complex<double> tmp_2861;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2861 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_2860 += tmp_2861;
   result += (-4) * tmp_2860;
   std::complex<double> tmp_2862;
   std::complex<double> tmp_2863;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2863 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_2862 += tmp_2863;
   result += (1.3333333333333333*MGlu) * tmp_2862;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2864;
   std::complex<double> tmp_2865;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2866;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2866 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_2865 += tmp_2866;
   }
   tmp_2864 += tmp_2865;
   result += (-0.5) * tmp_2864;
   std::complex<double> tmp_2867;
   std::complex<double> tmp_2868;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2869;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2869 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_2868 += tmp_2869;
   }
   tmp_2867 += tmp_2868;
   result += (-0.5) * tmp_2867;
   std::complex<double> tmp_2870;
   std::complex<double> tmp_2871;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2872;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2872 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPR(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_2871 += tmp_2872;
   }
   tmp_2870 += tmp_2871;
   result += (-0.5) * tmp_2870;
   std::complex<double> tmp_2873;
   std::complex<double> tmp_2874;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2875;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2875 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_2874 += tmp_2875;
   }
   tmp_2873 += tmp_2874;
   result += (-0.5) * tmp_2873;
   std::complex<double> tmp_2876;
   std::complex<double> tmp_2877;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2877 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_2876 += tmp_2877;
   result += (-0.6666666666666666) * tmp_2876;
   std::complex<double> tmp_2878;
   std::complex<double> tmp_2879;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2880;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2880 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_2879 += tmp_2880;
   }
   tmp_2878 += tmp_2879;
   result += (-0.5) * tmp_2878;
   std::complex<double> tmp_2881;
   std::complex<double> tmp_2882;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2882 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_2881 += tmp_2882;
   result += (-1) * tmp_2881;
   std::complex<double> tmp_2883;
   std::complex<double> tmp_2884;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2884 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_2883 += tmp_2884;
   result += (-1.3333333333333333) * tmp_2883;
   std::complex<double> tmp_2885;
   std::complex<double> tmp_2886;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2886 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_2885 += tmp_2886;
   result += (-1) * tmp_2885;
   std::complex<double> tmp_2887;
   std::complex<double> tmp_2888;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2888 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_2887 += tmp_2888;
   result += (-1) * tmp_2887;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2889;
   std::complex<double> tmp_2890;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2891;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2891 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_2890 += tmp_2891;
   }
   tmp_2889 += tmp_2890;
   result += (-0.5) * tmp_2889;
   std::complex<double> tmp_2892;
   std::complex<double> tmp_2893;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2894;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2894 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_2893 += tmp_2894;
   }
   tmp_2892 += tmp_2893;
   result += (-0.5) * tmp_2892;
   std::complex<double> tmp_2895;
   std::complex<double> tmp_2896;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2897;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2897 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarUFubarChaSdPL(gO2
            ,gI1,gI2))*CpbarUFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_2896 += tmp_2897;
   }
   tmp_2895 += tmp_2896;
   result += (-0.5) * tmp_2895;
   std::complex<double> tmp_2898;
   std::complex<double> tmp_2899;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2900;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2900 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_2899 += tmp_2900;
   }
   tmp_2898 += tmp_2899;
   result += (-0.5) * tmp_2898;
   std::complex<double> tmp_2901;
   std::complex<double> tmp_2902;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2902 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_2901 += tmp_2902;
   result += (-0.6666666666666666) * tmp_2901;
   std::complex<double> tmp_2903;
   std::complex<double> tmp_2904;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2905;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2905 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_2904 += tmp_2905;
   }
   tmp_2903 += tmp_2904;
   result += (-0.5) * tmp_2903;
   std::complex<double> tmp_2906;
   std::complex<double> tmp_2907;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2907 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_2906 += tmp_2907;
   result += (-1) * tmp_2906;
   std::complex<double> tmp_2908;
   std::complex<double> tmp_2909;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2909 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_2908 += tmp_2909;
   result += (-1.3333333333333333) * tmp_2908;
   std::complex<double> tmp_2910;
   std::complex<double> tmp_2911;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2911 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_2910 += tmp_2911;
   result += (-1) * tmp_2910;
   std::complex<double> tmp_2912;
   std::complex<double> tmp_2913;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2913 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_2912 += tmp_2913;
   result += (-1) * tmp_2912;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2914;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2915;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2915 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpGluconjSdFdPL(gI1,gI2
            ))*CpGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_2914 += tmp_2915;
   }
   result += tmp_2914;
   std::complex<double> tmp_2916;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2917;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2917 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpGluconjSuFuPL(gI1,gI2
            ))*CpGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_2916 += tmp_2917;
   }
   result += tmp_2916;
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpGluVGGluPR())*CpGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_2918;
   std::complex<double> tmp_2919;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2920;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2920 += AbsSqr(CpGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_2919 += tmp_2920;
   }
   tmp_2918 += tmp_2919;
   result += (-0.5) * tmp_2918;
   std::complex<double> tmp_2921;
   std::complex<double> tmp_2922;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2923;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2923 += AbsSqr(CpGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_2922 += tmp_2923;
   }
   tmp_2921 += tmp_2922;
   result += (-0.5) * tmp_2921;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += -3*AbsSqr(CpGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_2924;
   std::complex<double> tmp_2925;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2926;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2926 += AbsSqr(CpGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),MSd(
            gI1));
      }
      tmp_2925 += tmp_2926;
   }
   tmp_2924 += tmp_2925;
   result += (-0.5) * tmp_2924;
   std::complex<double> tmp_2927;
   std::complex<double> tmp_2928;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2929;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2929 += AbsSqr(CpGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),MSu(
            gI1));
      }
      tmp_2928 += tmp_2929;
   }
   tmp_2927 += tmp_2928;
   result += (-0.5) * tmp_2927;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2930;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2930 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2930;
   std::complex<double> tmp_2931;
   std::complex<double> tmp_2932;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2932 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2931 += tmp_2932;
   result += (0.5) * tmp_2931;
   std::complex<double> tmp_2933;
   std::complex<double> tmp_2934;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2935;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2935 += AbsSqr(CpVZhhAh(gI1,1 + gI2))*B00(p,MAh(1 + gI2),Mhh
            (gI1));
      }
      tmp_2934 += tmp_2935;
   }
   tmp_2933 += tmp_2934;
   result += (-4) * tmp_2933;
   std::complex<double> tmp_2936;
   std::complex<double> tmp_2937;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2938;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2938 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2937 += tmp_2938;
   }
   tmp_2936 += tmp_2937;
   result += (-4) * tmp_2936;
   std::complex<double> tmp_2939;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2940;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2940 += (AbsSqr(CpVZbarChaChaPL(gI1,gI2)) + AbsSqr(
            CpVZbarChaChaPR(gI1,gI2)))*H0(p,MCha(gI1),MCha(gI2));
         tmp_2940 += 4*B0(p,MCha(gI1),MCha(gI2))*MCha(gI1)*MCha(gI2)*Re(
            Conj(CpVZbarChaChaPL(gI1,gI2))*CpVZbarChaChaPR(gI1,gI2));
      }
      tmp_2939 += tmp_2940;
   }
   result += tmp_2939;
   std::complex<double> tmp_2941;
   std::complex<double> tmp_2942;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2943;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2943 += (AbsSqr(CpVZChiChiPL(gI1,gI2)) + AbsSqr(CpVZChiChiPR
            (gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2943 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZChiChiPL(gI1,gI2))*CpVZChiChiPR(gI1,gI2));
      }
      tmp_2942 += tmp_2943;
   }
   tmp_2941 += tmp_2942;
   result += (0.5) * tmp_2941;
   std::complex<double> tmp_2944;
   std::complex<double> tmp_2945;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2945 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2944 += tmp_2945;
   result += (3) * tmp_2944;
   std::complex<double> tmp_2946;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2946 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2946;
   std::complex<double> tmp_2947;
   std::complex<double> tmp_2948;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2948 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2947 += tmp_2948;
   result += (3) * tmp_2947;
   std::complex<double> tmp_2949;
   std::complex<double> tmp_2950;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2951;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2951 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2950 += tmp_2951;
   }
   tmp_2949 += tmp_2950;
   result += (-12) * tmp_2949;
   std::complex<double> tmp_2952;
   std::complex<double> tmp_2953;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2954;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2954 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2953 += tmp_2954;
   }
   tmp_2952 += tmp_2953;
   result += (-4) * tmp_2952;
   std::complex<double> tmp_2955;
   std::complex<double> tmp_2956;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2957;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2957 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2956 += tmp_2957;
   }
   tmp_2955 += tmp_2956;
   result += (-12) * tmp_2955;
   std::complex<double> tmp_2958;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2958 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2958;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2959;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2959 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_2959;
   std::complex<double> tmp_2960;
   std::complex<double> tmp_2961;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2961 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_2960 += tmp_2961;
   result += (0.5) * tmp_2960;
   std::complex<double> tmp_2962;
   std::complex<double> tmp_2963;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2964;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2964 += AbsSqr(CpconjVWmHpmhh(1 + gI1,gI2))*B00(p,Mhh(gI2),
            MHpm(1 + gI1));
      }
      tmp_2963 += tmp_2964;
   }
   tmp_2962 += tmp_2963;
   result += (-4) * tmp_2962;
   std::complex<double> tmp_2965;
   std::complex<double> tmp_2966;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2967;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2967 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_2966 += tmp_2967;
   }
   tmp_2965 += tmp_2966;
   result += (-4) * tmp_2965;
   std::complex<double> tmp_2968;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      std::complex<double> tmp_2969;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2969 += (AbsSqr(CpconjVWmChiChaPL(gI1,gI2)) + AbsSqr(
            CpconjVWmChiChaPR(gI1,gI2)))*H0(p,MChi(gI1),MCha(gI2));
         tmp_2969 += 4*B0(p,MChi(gI1),MCha(gI2))*MCha(gI2)*MChi(gI1)*Re(
            Conj(CpconjVWmChiChaPL(gI1,gI2))*CpconjVWmChiChaPR(gI1,gI2));
      }
      tmp_2968 += tmp_2969;
   }
   result += tmp_2968;
   std::complex<double> tmp_2970;
   std::complex<double> tmp_2971;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2971 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_2970 += tmp_2971;
   result += (3) * tmp_2970;
   std::complex<double> tmp_2972;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2972 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_2972;
   std::complex<double> tmp_2973;
   std::complex<double> tmp_2974;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2974 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_2973 += tmp_2974;
   result += (3) * tmp_2973;
   std::complex<double> tmp_2975;
   std::complex<double> tmp_2976;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2977;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2977 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_2976 += tmp_2977;
   }
   tmp_2975 += tmp_2976;
   result += (-12) * tmp_2975;
   std::complex<double> tmp_2978;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2978 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_2978;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2979;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2980;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2980 += B0(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_2979 += tmp_2980;
   }
   result += tmp_2979;
   std::complex<double> tmp_2981;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2982;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2982 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2981 += tmp_2982;
   }
   result += tmp_2981;
   std::complex<double> tmp_2983;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2984;
      std::complex<double> tmp_2985;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2985 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_2984 += tmp_2985;
      tmp_2983 += (MFe(gI1)) * tmp_2984;
   }
   result += tmp_2983;
   std::complex<double> tmp_2986;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2987;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2987 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2986 += tmp_2987;
   }
   result += tmp_2986;
   std::complex<double> tmp_2988;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2989;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_2989 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2988 += tmp_2989;
   }
   result += tmp_2988;
   std::complex<double> tmp_2990;
   std::complex<double> tmp_2991;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2991 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_2990 += tmp_2991;
   result += (-4) * tmp_2990;
   std::complex<double> tmp_2992;
   std::complex<double> tmp_2993;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2993 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_2992 += tmp_2993;
   result += (-4) * tmp_2992;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2994;
   std::complex<double> tmp_2995;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2996;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2996 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_2995 += tmp_2996;
   }
   tmp_2994 += tmp_2995;
   result += (-0.5) * tmp_2994;
   std::complex<double> tmp_2997;
   std::complex<double> tmp_2998;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2999;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2999 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_2998 += tmp_2999;
   }
   tmp_2997 += tmp_2998;
   result += (-0.5) * tmp_2997;
   std::complex<double> tmp_3000;
   std::complex<double> tmp_3001;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3002;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3002 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPR(gO2,gI1
            ,gI2))*CpbarFeSvChaPR(gO1,gI1,gI2);
      }
      tmp_3001 += tmp_3002;
   }
   tmp_3000 += tmp_3001;
   result += (-0.5) * tmp_3000;
   std::complex<double> tmp_3003;
   std::complex<double> tmp_3004;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3005;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3005 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3004 += tmp_3005;
   }
   tmp_3003 += tmp_3004;
   result += (-0.5) * tmp_3003;
   std::complex<double> tmp_3006;
   std::complex<double> tmp_3007;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3008;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3008 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3007 += tmp_3008;
   }
   tmp_3006 += tmp_3007;
   result += (-0.5) * tmp_3006;
   std::complex<double> tmp_3009;
   std::complex<double> tmp_3010;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3010 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_3009 += tmp_3010;
   result += (-1) * tmp_3009;
   std::complex<double> tmp_3011;
   std::complex<double> tmp_3012;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3012 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_3011 += tmp_3012;
   result += (-1) * tmp_3011;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3013;
   std::complex<double> tmp_3014;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3015;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3015 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3014 += tmp_3015;
   }
   tmp_3013 += tmp_3014;
   result += (-0.5) * tmp_3013;
   std::complex<double> tmp_3016;
   std::complex<double> tmp_3017;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3018;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3018 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_3017 += tmp_3018;
   }
   tmp_3016 += tmp_3017;
   result += (-0.5) * tmp_3016;
   std::complex<double> tmp_3019;
   std::complex<double> tmp_3020;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3021;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3021 += B1(p,MCha(gI2),MSv(gI1))*Conj(CpbarFeSvChaPL(gO2,gI1
            ,gI2))*CpbarFeSvChaPL(gO1,gI1,gI2);
      }
      tmp_3020 += tmp_3021;
   }
   tmp_3019 += tmp_3020;
   result += (-0.5) * tmp_3019;
   std::complex<double> tmp_3022;
   std::complex<double> tmp_3023;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3024;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3024 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3023 += tmp_3024;
   }
   tmp_3022 += tmp_3023;
   result += (-0.5) * tmp_3022;
   std::complex<double> tmp_3025;
   std::complex<double> tmp_3026;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3027;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3027 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3026 += tmp_3027;
   }
   tmp_3025 += tmp_3026;
   result += (-0.5) * tmp_3025;
   std::complex<double> tmp_3028;
   std::complex<double> tmp_3029;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3029 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_3028 += tmp_3029;
   result += (-1) * tmp_3028;
   std::complex<double> tmp_3030;
   std::complex<double> tmp_3031;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3031 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_3030 += tmp_3031;
   result += (-1) * tmp_3030;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3032;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3033;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3033 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3032 += tmp_3033;
   }
   result += tmp_3032;
   std::complex<double> tmp_3034;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3035;
      std::complex<double> tmp_3036;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3036 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3035 += tmp_3036;
      tmp_3034 += (MFd(gI1)) * tmp_3035;
   }
   result += tmp_3034;
   std::complex<double> tmp_3037;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3038;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3038 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3037 += tmp_3038;
   }
   result += tmp_3037;
   std::complex<double> tmp_3039;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3040;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3040 += B0(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2)*MCha(gI2);
      }
      tmp_3039 += tmp_3040;
   }
   result += tmp_3039;
   std::complex<double> tmp_3041;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3042;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3042 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3041 += tmp_3042;
   }
   result += tmp_3041;
   std::complex<double> tmp_3043;
   std::complex<double> tmp_3044;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3044 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3043 += tmp_3044;
   result += (-4) * tmp_3043;
   std::complex<double> tmp_3045;
   std::complex<double> tmp_3046;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3046 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3045 += tmp_3046;
   result += (-4) * tmp_3045;
   std::complex<double> tmp_3047;
   std::complex<double> tmp_3048;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3048 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3047 += tmp_3048;
   result += (1.3333333333333333*MGlu) * tmp_3047;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3049;
   std::complex<double> tmp_3050;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3051;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3051 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3050 += tmp_3051;
   }
   tmp_3049 += tmp_3050;
   result += (-0.5) * tmp_3049;
   std::complex<double> tmp_3052;
   std::complex<double> tmp_3053;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3054;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3054 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3053 += tmp_3054;
   }
   tmp_3052 += tmp_3053;
   result += (-0.5) * tmp_3052;
   std::complex<double> tmp_3055;
   std::complex<double> tmp_3056;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3057;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3057 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3056 += tmp_3057;
   }
   tmp_3055 += tmp_3056;
   result += (-0.5) * tmp_3055;
   std::complex<double> tmp_3058;
   std::complex<double> tmp_3059;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3059 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3058 += tmp_3059;
   result += (-0.6666666666666666) * tmp_3058;
   std::complex<double> tmp_3060;
   std::complex<double> tmp_3061;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3062;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3062 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPR(gO2,gI1
            ,gI2))*CpbarFdSuChaPR(gO1,gI1,gI2);
      }
      tmp_3061 += tmp_3062;
   }
   tmp_3060 += tmp_3061;
   result += (-0.5) * tmp_3060;
   std::complex<double> tmp_3063;
   std::complex<double> tmp_3064;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3065;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3065 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3064 += tmp_3065;
   }
   tmp_3063 += tmp_3064;
   result += (-0.5) * tmp_3063;
   std::complex<double> tmp_3066;
   std::complex<double> tmp_3067;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3067 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_3066 += tmp_3067;
   result += (-1) * tmp_3066;
   std::complex<double> tmp_3068;
   std::complex<double> tmp_3069;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3069 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_3068 += tmp_3069;
   result += (-1) * tmp_3068;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3070;
   std::complex<double> tmp_3071;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3072;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3072 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3071 += tmp_3072;
   }
   tmp_3070 += tmp_3071;
   result += (-0.5) * tmp_3070;
   std::complex<double> tmp_3073;
   std::complex<double> tmp_3074;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3075;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3075 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3074 += tmp_3075;
   }
   tmp_3073 += tmp_3074;
   result += (-0.5) * tmp_3073;
   std::complex<double> tmp_3076;
   std::complex<double> tmp_3077;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3078;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3078 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3077 += tmp_3078;
   }
   tmp_3076 += tmp_3077;
   result += (-0.5) * tmp_3076;
   std::complex<double> tmp_3079;
   std::complex<double> tmp_3080;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3080 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPL(gO1,gI1,1);
   }
   tmp_3079 += tmp_3080;
   result += (-0.6666666666666666) * tmp_3079;
   std::complex<double> tmp_3081;
   std::complex<double> tmp_3082;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3083;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3083 += B1(p,MCha(gI2),MSu(gI1))*Conj(CpbarFdSuChaPL(gO2,gI1
            ,gI2))*CpbarFdSuChaPL(gO1,gI1,gI2);
      }
      tmp_3082 += tmp_3083;
   }
   tmp_3081 += tmp_3082;
   result += (-0.5) * tmp_3081;
   std::complex<double> tmp_3084;
   std::complex<double> tmp_3085;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3086;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3086 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3085 += tmp_3086;
   }
   tmp_3084 += tmp_3085;
   result += (-0.5) * tmp_3084;
   std::complex<double> tmp_3087;
   std::complex<double> tmp_3088;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3088 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_3087 += tmp_3088;
   result += (-1) * tmp_3087;
   std::complex<double> tmp_3089;
   std::complex<double> tmp_3090;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3090 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_3089 += tmp_3090;
   result += (-1) * tmp_3089;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3091;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3092;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3092 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3091 += tmp_3092;
   }
   result += tmp_3091;
   std::complex<double> tmp_3093;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3094;
      std::complex<double> tmp_3095;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3095 += B0(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3094 += tmp_3095;
      tmp_3093 += (MCha(gI1)) * tmp_3094;
   }
   result += tmp_3093;
   std::complex<double> tmp_3096;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3097;
      std::complex<double> tmp_3098;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3098 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3097 += tmp_3098;
      tmp_3096 += (MFu(gI1)) * tmp_3097;
   }
   result += tmp_3096;
   std::complex<double> tmp_3099;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3100;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3100 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3099 += tmp_3100;
   }
   result += tmp_3099;
   std::complex<double> tmp_3101;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3102;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3102 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3101 += tmp_3102;
   }
   result += tmp_3101;
   std::complex<double> tmp_3103;
   std::complex<double> tmp_3104;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3104 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3103 += tmp_3104;
   result += (-4) * tmp_3103;
   std::complex<double> tmp_3105;
   std::complex<double> tmp_3106;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3106 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3105 += tmp_3106;
   result += (-4) * tmp_3105;
   std::complex<double> tmp_3107;
   std::complex<double> tmp_3108;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3108 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3107 += tmp_3108;
   result += (-4) * tmp_3107;
   std::complex<double> tmp_3109;
   std::complex<double> tmp_3110;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3110 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3109 += tmp_3110;
   result += (1.3333333333333333*MGlu) * tmp_3109;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3111;
   std::complex<double> tmp_3112;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3113;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3113 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3112 += tmp_3113;
   }
   tmp_3111 += tmp_3112;
   result += (-0.5) * tmp_3111;
   std::complex<double> tmp_3114;
   std::complex<double> tmp_3115;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3116;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3116 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3115 += tmp_3116;
   }
   tmp_3114 += tmp_3115;
   result += (-0.5) * tmp_3114;
   std::complex<double> tmp_3117;
   std::complex<double> tmp_3118;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3119;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3119 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPR(gO2,
            gI1,gI2))*CpbarFubarChaSdPR(gO1,gI1,gI2);
      }
      tmp_3118 += tmp_3119;
   }
   tmp_3117 += tmp_3118;
   result += (-0.5) * tmp_3117;
   std::complex<double> tmp_3120;
   std::complex<double> tmp_3121;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3122;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3122 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3121 += tmp_3122;
   }
   tmp_3120 += tmp_3121;
   result += (-0.5) * tmp_3120;
   std::complex<double> tmp_3123;
   std::complex<double> tmp_3124;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3124 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3123 += tmp_3124;
   result += (-0.6666666666666666) * tmp_3123;
   std::complex<double> tmp_3125;
   std::complex<double> tmp_3126;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3127;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3127 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3126 += tmp_3127;
   }
   tmp_3125 += tmp_3126;
   result += (-0.5) * tmp_3125;
   std::complex<double> tmp_3128;
   std::complex<double> tmp_3129;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3129 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3128 += tmp_3129;
   result += (-1) * tmp_3128;
   std::complex<double> tmp_3130;
   std::complex<double> tmp_3131;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3131 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_3130 += tmp_3131;
   result += (-1) * tmp_3130;
   std::complex<double> tmp_3132;
   std::complex<double> tmp_3133;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3133 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_3132 += tmp_3133;
   result += (-1) * tmp_3132;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3134;
   std::complex<double> tmp_3135;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3136;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3136 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3135 += tmp_3136;
   }
   tmp_3134 += tmp_3135;
   result += (-0.5) * tmp_3134;
   std::complex<double> tmp_3137;
   std::complex<double> tmp_3138;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3139;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3139 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3138 += tmp_3139;
   }
   tmp_3137 += tmp_3138;
   result += (-0.5) * tmp_3137;
   std::complex<double> tmp_3140;
   std::complex<double> tmp_3141;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3142;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3142 += B1(p,MCha(gI1),MSd(gI2))*Conj(CpbarFubarChaSdPL(gO2,
            gI1,gI2))*CpbarFubarChaSdPL(gO1,gI1,gI2);
      }
      tmp_3141 += tmp_3142;
   }
   tmp_3140 += tmp_3141;
   result += (-0.5) * tmp_3140;
   std::complex<double> tmp_3143;
   std::complex<double> tmp_3144;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3145;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3145 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3144 += tmp_3145;
   }
   tmp_3143 += tmp_3144;
   result += (-0.5) * tmp_3143;
   std::complex<double> tmp_3146;
   std::complex<double> tmp_3147;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3147 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPL(gO1,gI1,1);
   }
   tmp_3146 += tmp_3147;
   result += (-0.6666666666666666) * tmp_3146;
   std::complex<double> tmp_3148;
   std::complex<double> tmp_3149;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3150;
      for (unsigned gI2 = 0; gI2 < 5; ++gI2) {
         tmp_3150 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3149 += tmp_3150;
   }
   tmp_3148 += tmp_3149;
   result += (-0.5) * tmp_3148;
   std::complex<double> tmp_3151;
   std::complex<double> tmp_3152;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3152 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3151 += tmp_3152;
   result += (-1) * tmp_3151;
   std::complex<double> tmp_3153;
   std::complex<double> tmp_3154;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3154 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_3153 += tmp_3154;
   result += (-1) * tmp_3153;
   std::complex<double> tmp_3155;
   std::complex<double> tmp_3156;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3156 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_3155 += tmp_3156;
   result += (-1) * tmp_3155;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::tadpole_hh(unsigned gO1) const
{
   std::complex<double> result;

   result += A0(MVWm)*CpUhhbargWmCgWmC(gO1);
   result += A0(MVWm)*CpUhhbargWmgWm(gO1);
   result += A0(MVZ)*CpUhhbargZgZ(gO1);
   result += 4*A0(MVWm)*CpUhhconjVWmVWm(gO1);
   result += 2*A0(MVZ)*CpUhhVZVZ(gO1);
   std::complex<double> tmp_3157;
   std::complex<double> tmp_3158;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3158 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_3157 += tmp_3158;
   result += (-0.5) * tmp_3157;
   std::complex<double> tmp_3159;
   std::complex<double> tmp_3160;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3160 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_3159 += tmp_3160;
   result += (-1) * tmp_3159;
   std::complex<double> tmp_3161;
   std::complex<double> tmp_3162;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3162 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_3161 += tmp_3162;
   result += (-0.5) * tmp_3161;
   std::complex<double> tmp_3163;
   std::complex<double> tmp_3164;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3164 += A0(MCha(gI1))*(CpUhhbarChaChaPL(gO1,gI1,gI1) +
         CpUhhbarChaChaPR(gO1,gI1,gI1))*MCha(gI1);
   }
   tmp_3163 += tmp_3164;
   result += (2) * tmp_3163;
   std::complex<double> tmp_3165;
   std::complex<double> tmp_3166;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3166 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_3165 += tmp_3166;
   result += (6) * tmp_3165;
   std::complex<double> tmp_3167;
   std::complex<double> tmp_3168;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3168 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_3167 += tmp_3168;
   result += (2) * tmp_3167;
   std::complex<double> tmp_3169;
   std::complex<double> tmp_3170;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3170 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_3169 += tmp_3170;
   result += (6) * tmp_3169;
   std::complex<double> tmp_3171;
   std::complex<double> tmp_3172;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3172 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_3171 += tmp_3172;
   result += (-1) * tmp_3171;
   std::complex<double> tmp_3173;
   for (unsigned gI1 = 0; gI1 < 5; ++gI1) {
      tmp_3173 += A0(MChi(gI1))*(CpUhhChiChiPL(gO1,gI1,gI1) + CpUhhChiChiPR(
         gO1,gI1,gI1))*MChi(gI1);
   }
   result += tmp_3173;
   std::complex<double> tmp_3174;
   std::complex<double> tmp_3175;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3175 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_3174 += tmp_3175;
   result += (-3) * tmp_3174;
   std::complex<double> tmp_3176;
   std::complex<double> tmp_3177;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3177 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_3176 += tmp_3177;
   result += (-1) * tmp_3176;
   std::complex<double> tmp_3178;
   std::complex<double> tmp_3179;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3179 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_3178 += tmp_3179;
   result += (-3) * tmp_3178;

   return result * oneOver16PiSqr;

}








void CLASSNAME::calculate_MGlu_pole()
{
   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_Glu());
   const double p = MGlu;
   const double self_energy_1  = Re(self_energy_Glu_1(p));
   const double self_energy_PL = Re(self_energy_Glu_PL(p));
   const double self_energy_PR = Re(self_energy_Glu_PR(p));
   PHYSICAL(MGlu) = M_tree - self_energy_1 - M_tree * (self_energy_PL +
      self_energy_PR);
}

void CLASSNAME::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
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
         problems.flag_bad_mass(TMSSM_info::Sd, eigenvalue_error >
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
         problems.flag_bad_mass(TMSSM_info::Sv, eigenvalue_error >
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
         problems.flag_bad_mass(TMSSM_info::Su, eigenvalue_error >
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
         problems.flag_bad_mass(TMSSM_info::Se, eigenvalue_error >
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

      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_Mhh(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_hh(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(TMSSM_info::hh,
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

      for (unsigned es = 0; es < 3; ++es) {
         const double p = Abs(old_MAh(es));
         for (unsigned i1 = 0; i1 < 3; ++i1) {
            for (unsigned i2 = i1; i2 < 3; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Ah(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,3,3> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,3,1> eigen_values;
         Eigen::Matrix<double,3,3> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(TMSSM_info::Ah,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZA);
         #endif

         if (eigen_values(es) < 0.)
            problems.flag_tachyon(Ah);

         PHYSICAL(MAh(es)) = AbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZA) = mix_ZA;
      }

      new_MAh = PHYSICAL(MAh);
      diff = MaxRelDiff(new_MAh, old_MAh);
      old_MAh = new_MAh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);
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
      Eigen::Matrix<double,4,4> self_energy;
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Hpm());

      for (unsigned es = 0; es < 4; ++es) {
         const double p = Abs(old_MHpm(es));
         for (unsigned i1 = 0; i1 < 4; ++i1) {
            for (unsigned i2 = i1; i2 < 4; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Hpm(p,i1,
                  i2));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,4,4> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZP;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZP, eigenvalue_error);
            problems.flag_bad_mass(TMSSM_info::Hpm,
               eigenvalue_error > precision * Abs(eigen_values(0)));
         #else
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZP);
         #endif

         if (eigen_values(es) < 0.)
            problems.flag_tachyon(Hpm);

         PHYSICAL(MHpm(es)) = AbsSqrt(eigen_values(es));
         if (es == 0)
            PHYSICAL(ZP) = mix_ZP;
      }

      new_MHpm = PHYSICAL(MHpm);
      diff = MaxRelDiff(new_MHpm, old_MHpm);
      old_MHpm = new_MHpm;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);
}

void CLASSNAME::calculate_MChi_pole()
{
   // diagonalization with medium precision
   Eigen::Matrix<double,5,5> self_energy_1;
   Eigen::Matrix<double,5,5> self_energy_PL;
   Eigen::Matrix<double,5,5> self_energy_PR;
   const Eigen::Matrix<double,5,5> M_tree(get_mass_matrix_Chi());
   for (unsigned es = 0; es < 5; ++es) {
      const double p = Abs(MChi(es));
      for (unsigned i1 = 0; i1 < 5; ++i1) {
         for (unsigned i2 = 0; i2 < 5; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Chi_1(p,i1,i2
               ));
            self_energy_PL(i1,i2) = Re(self_energy_Chi_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Chi_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,5,5> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,5,5> M_1loop(M_tree + 0.5 * (delta_M
         + delta_M.transpose()));
      Eigen::Array<double,5,1> eigen_values;
      decltype(ZN) mix_ZN;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_symmetric(M_1loop, eigen_values, mix_ZN,
            eigenvalue_error);
         problems.flag_bad_mass(TMSSM_info::Chi, eigenvalue_error >
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
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Cha());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MCha(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Cha_1(p,i1,i2
               ));
            self_energy_PL(i1,i2) = Re(self_energy_Cha_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Cha_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(UM) mix_UM;
      decltype(UP) mix_UP;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_UM, mix_UP, eigenvalue_error);
      problems.flag_bad_mass(TMSSM_info::Cha, eigenvalue_error >
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
      problems.flag_bad_mass(TMSSM_info::Fe, eigenvalue_error >
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
      problems.flag_bad_mass(TMSSM_info::Fd, eigenvalue_error >
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
   Eigen::Matrix<double,3,3> self_energy_1;
   Eigen::Matrix<double,3,3> self_energy_PL;
   Eigen::Matrix<double,3,3> self_energy_PR;
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (unsigned es = 0; es < 3; ++es) {
      const double p = Abs(MFu(es));
      for (unsigned i1 = 0; i1 < 3; ++i1) {
         for (unsigned i2 = 0; i2 < 3; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Fu_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fu_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fu_PR(p,i1,i2
               ));
         }
      }
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZUL) mix_ZUL;
      decltype(ZUR) mix_ZUR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZUL, mix_ZUR, eigenvalue_error
         );
      problems.flag_bad_mass(TMSSM_info::Fu, eigenvalue_error >
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

void CLASSNAME::calculate_MVG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVG) = 0.;
}

void CLASSNAME::calculate_MVP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MVP) = 0.;
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
   const double m_sm_drbar = m_sm_msbar * (1 - 0.00020496318737651018*
      Power(g3,4) + 0.0006860288475783287*Sqr(g1) + 0.0023747152416172916*Sqr(
      g2) - 0.008443431970194815*Sqr(g3));

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
   const double m_sm_drbar = m_sm_msbar * (1 - 0.0023747152416172916*(0.6
      *Sqr(g1) - Sqr(g2)));

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


std::ostream& operator<<(std::ostream& ostr, const TMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
