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

// File generated at Tue 24 Feb 2015 17:34:04

/**
 * @file MRSSM_two_scale_model.cpp
 * @brief implementation of the MRSSM model class
 *
 * Contains the definition of the MRSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 24 Feb 2015 17:34:04 with FlexibleSUSY
 * 1.0.4 (git commit: v1.0.4-355-g539c238) and SARAH 4.4.6 .
 */

#include "MRSSM_two_scale_model.hpp"
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

using namespace MRSSM_info;

#define CLASSNAME MRSSM<Two_scale>

#define PHYSICAL(parameter) physical.parameter
#define INPUT(parameter) model->get_input().parameter
#define LOCALINPUT(parameter) input.parameter

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

CLASSNAME::MRSSM(const MRSSM_input_parameters& input_)
   : Two_scale_model()
   , MRSSM_soft_parameters(input_)
   , number_of_ewsb_iterations(100)
   , number_of_mass_iterations(20)
   , ewsb_loop_order(2)
   , pole_mass_loop_order(2)
   , calculate_sm_pole_masses(false)
   , force_output(false)
   , precision(1.0e-3)
   , ewsb_iteration_precision(1.0e-5)
   , physical()
   , problems(MRSSM_info::particle_names)
   , two_loop_corrections()
#ifdef ENABLE_THREADS
   , thread_exception()
#endif
   , MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MSOc(0), MVP(0),
      MVZ(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>
      ::Zero()), MSu(Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,
      1>::Zero()), Mhh(Eigen::Array<double,4,1>::Zero()), MAh(Eigen::Array<double,
      4,1>::Zero()), MRh(Eigen::Array<double,2,1>::Zero()), MHpm(Eigen::Array<
      double,4,1>::Zero()), MRpm(Eigen::Array<double,2,1>::Zero()), MChi(
      Eigen::Array<double,4,1>::Zero()), MCha1(Eigen::Array<double,2,1>::Zero()),
      MCha2(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero(
      )), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double,3,1>
      ::Zero()), MVWm(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,4,4>::Zero()), ZA(Eigen::Matrix<double,4,
      4>::Zero()), ZHR(Eigen::Matrix<double,2,2>::Zero()), ZP(Eigen::Matrix<double
      ,4,4>::Zero()), ZRP(Eigen::Matrix<double,2,2>::Zero()), ZN1(Eigen::Matrix<
      std::complex<double>,4,4>::Zero()), ZN2(Eigen::Matrix<std::complex<double>,4
      ,4>::Zero()), UM1(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP1(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), UM2(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), UP2(Eigen::Matrix<std::complex<double>,2
      ,2>::Zero()), ZEL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZUL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(
      Eigen::Matrix<std::complex<double>,3,3>::Zero())


{
}

CLASSNAME::~MRSSM()
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

const Problems<MRSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems() const
{
   return problems;
}

Problems<MRSSM_info::NUMBER_OF_PARTICLES>& CLASSNAME::get_problems()
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
   tadpole[3] = get_ewsb_eq_hh_4();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh(0));
      tadpole[1] -= Re(tadpole_hh(1));
      tadpole[2] -= Re(tadpole_hh(3));
      tadpole[3] -= Re(tadpole_hh(2));

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
   MRSSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   model->set_mHd2(gsl_vector_get(x, 0));
   model->set_mHu2(gsl_vector_get(x, 1));
   model->set_mS2(gsl_vector_get(x, 2));
   model->set_mT2(gsl_vector_get(x, 3));


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

   mHd2 = solver->get_solution(0);
   mHu2 = solver->get_solution(1);
   mS2 = solver->get_solution(2);
   mT2 = solver->get_solution(3);


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

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;
   const double old_mT2 = mT2;
   const double old_mS2 = mS2;

   mHd2 = (0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*vd*
      AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS*
      Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*AbsSqr(
      LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) + 5*vd*
      Sqr(g2)*Sqr(vu)))/vd;
   mHu2 = (0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40*vu
      *AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*vu*
      Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*vu*
      Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*Sqr(g1)*Sqr
      (vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu*AbsSqr(
      LamTU)*Sqr(vT)))/vu;
   mT2 = (0.125*(-32*vT*Sqr(MDWBT) - 2*g2*MDWBT*Sqr(vd) - 2*vT*AbsSqr(LamTD)*
      Sqr(vd) - 1.4142135623730951*LamTD*vS*Conj(LamSD)*Sqr(vd) - 2*MuD*Conj(LamTD
      )*Sqr(vd) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*Sqr(vd) - 2*g2*Conj(
      MDWBT)*Sqr(vd) - 2*LamTD*Conj(MuD)*Sqr(vd) + 2*g2*MDWBT*Sqr(vu) - 2*vT*
      AbsSqr(LamTU)*Sqr(vu) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*Sqr(vu) + 2*
      MuU*Conj(LamTU)*Sqr(vu) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*Sqr(vu) +
      2*g2*Conj(MDWBT)*Sqr(vu) + 2*LamTU*Conj(MuU)*Sqr(vu)))/vT;
   mS2 = (0.025*(-160*vS*Sqr(MDBS) + 7.745966692414834*g1*MDBS*Sqr(vd) - 20*vS*
      AbsSqr(LamSD)*Sqr(vd) - 14.142135623730951*MuD*Conj(LamSD)*Sqr(vd) -
      7.0710678118654755*LamTD*vT*Conj(LamSD)*Sqr(vd) - 7.0710678118654755*LamSD*
      vT*Conj(LamTD)*Sqr(vd) + 7.745966692414834*g1*Conj(MDBS)*Sqr(vd) -
      14.142135623730951*LamSD*Conj(MuD)*Sqr(vd) - 7.745966692414834*g1*MDBS*Sqr(
      vu) - 20*vS*AbsSqr(LamSU)*Sqr(vu) - 14.142135623730951*MuU*Conj(LamSU)*Sqr(
      vu) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*Sqr(vu) + 7.0710678118654755*
      LamSU*vT*Conj(LamTU)*Sqr(vu) - 7.745966692414834*g1*Conj(MDBS)*Sqr(vu) -
      14.142135623730951*LamSU*Conj(MuU)*Sqr(vu)))/vS;

   const bool is_finite = std::isfinite(mHd2) && std::isfinite(mHu2) &&
      std::isfinite(mT2) && std::isfinite(mS2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      mT2 = old_mT2;
      mS2 = old_mS2;
      error = 1;
   }


   return error;
}

int CLASSNAME::solve_ewsb_tree_level_via_soft_higgs_masses()
{
   int error = 0;

   const double new_mHd2 = (0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*
      MDWBT*vd*vT - 40*vd*AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu -
      28.284271247461902*MuD*vd*vS*Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT
      *Conj(LamSD) - 20*MuD*vd*vT*Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*
      Conj(LamTD) + 15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(
      MDWBT) - 28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD)
      + 20*vu*Conj(BMu) - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) - 20*vd*
      AbsSqr(LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr(vu) +
      5*vd*Sqr(g2)*Sqr(vu)))/vd;
   const double new_mHu2 = (0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*
      MDWBT*vT*vu - 40*vu*AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu -
      28.284271247461902*MuU*vS*vu*Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu
      *Conj(LamSU) + 20*MuU*vT*vu*Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*
      Conj(LamTU) - 15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(
      MDWBT) - 28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU)
      + 20*vd*Conj(BMu) - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) + 3*vu*
      Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) - 10*vu
      *AbsSqr(LamTU)*Sqr(vT)))/vu;
   const double new_mT2 = (0.125*(-32*vT*Sqr(MDWBT) - 2*g2*MDWBT*Sqr(vd) - 2*vT
      *AbsSqr(LamTD)*Sqr(vd) - 1.4142135623730951*LamTD*vS*Conj(LamSD)*Sqr(vd) - 2
      *MuD*Conj(LamTD)*Sqr(vd) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*Sqr(vd) -
      2*g2*Conj(MDWBT)*Sqr(vd) - 2*LamTD*Conj(MuD)*Sqr(vd) + 2*g2*MDWBT*Sqr(vu) -
      2*vT*AbsSqr(LamTU)*Sqr(vu) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*Sqr(vu
      ) + 2*MuU*Conj(LamTU)*Sqr(vu) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*Sqr(
      vu) + 2*g2*Conj(MDWBT)*Sqr(vu) + 2*LamTU*Conj(MuU)*Sqr(vu)))/vT;
   const double new_mS2 = (0.025*(-160*vS*Sqr(MDBS) + 7.745966692414834*g1*MDBS
      *Sqr(vd) - 20*vS*AbsSqr(LamSD)*Sqr(vd) - 14.142135623730951*MuD*Conj(LamSD)*
      Sqr(vd) - 7.0710678118654755*LamTD*vT*Conj(LamSD)*Sqr(vd) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*Sqr(vd) + 7.745966692414834*g1*Conj(
      MDBS)*Sqr(vd) - 14.142135623730951*LamSD*Conj(MuD)*Sqr(vd) -
      7.745966692414834*g1*MDBS*Sqr(vu) - 20*vS*AbsSqr(LamSU)*Sqr(vu) -
      14.142135623730951*MuU*Conj(LamSU)*Sqr(vu) + 7.0710678118654755*LamTU*vT*
      Conj(LamSU)*Sqr(vu) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*Sqr(vu) -
      7.745966692414834*g1*Conj(MDBS)*Sqr(vu) - 14.142135623730951*LamSU*Conj(MuU)
      *Sqr(vu)))/vS;

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

   if (std::isfinite(new_mS2))
      mS2 = new_mS2;
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
   x_init[2] = mS2;
   x_init[3] = mT2;

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
      tadpole[2] += Re(tadpole_hh(3));
      tadpole[3] += Re(tadpole_hh(2));

      if (ewsb_loop_order > 1) {

      }
   }

   double mHd2;
   double mHu2;
   double mT2;
   double mS2;

   mHd2 = (0.025*(15.491933384829668*g1*MDBS*vd*vS - 20*g2*MDWBT*vd*vT - 40*vd*
      AbsSqr(MuD) - 40*vd*AbsSqr(Mu) + 20*vu*BMu - 28.284271247461902*MuD*vd*vS*
      Conj(LamSD) - 14.142135623730951*LamTD*vd*vS*vT*Conj(LamSD) - 20*MuD*vd*vT*
      Conj(LamTD) - 14.142135623730951*LamSD*vd*vS*vT*Conj(LamTD) +
      15.491933384829668*g1*vd*vS*Conj(MDBS) - 20*g2*vd*vT*Conj(MDWBT) -
      28.284271247461902*LamSD*vd*vS*Conj(MuD) - 20*LamTD*vd*vT*Conj(MuD) + 20*vu*
      Conj(BMu) + 40*tadpole[0] - 3*Power(vd,3)*Sqr(g1) - 5*Power(vd,3)*Sqr(g2) -
      20*vd*AbsSqr(LamSD)*Sqr(vS) - 10*vd*AbsSqr(LamTD)*Sqr(vT) + 3*vd*Sqr(g1)*Sqr
      (vu) + 5*vd*Sqr(g2)*Sqr(vu)))/vd;
   mHu2 = (0.025*(-15.491933384829668*g1*MDBS*vS*vu + 20*g2*MDWBT*vT*vu - 40*vu
      *AbsSqr(MuU) - 40*vu*AbsSqr(Mu) + 20*vd*BMu - 28.284271247461902*MuU*vS*vu*
      Conj(LamSU) + 14.142135623730951*LamTU*vS*vT*vu*Conj(LamSU) + 20*MuU*vT*vu*
      Conj(LamTU) + 14.142135623730951*LamSU*vS*vT*vu*Conj(LamTU) -
      15.491933384829668*g1*vS*vu*Conj(MDBS) + 20*g2*vT*vu*Conj(MDWBT) -
      28.284271247461902*LamSU*vS*vu*Conj(MuU) + 20*LamTU*vT*vu*Conj(MuU) + 20*vd*
      Conj(BMu) + 40*tadpole[1] - 3*Power(vu,3)*Sqr(g1) - 5*Power(vu,3)*Sqr(g2) +
      3*vu*Sqr(g1)*Sqr(vd) + 5*vu*Sqr(g2)*Sqr(vd) - 20*vu*AbsSqr(LamSU)*Sqr(vS) -
      10*vu*AbsSqr(LamTU)*Sqr(vT)))/vu;
   mT2 = (0.125*(8*tadpole[2] - 32*vT*Sqr(MDWBT) - 2*g2*MDWBT*Sqr(vd) - 2*vT*
      AbsSqr(LamTD)*Sqr(vd) - 1.4142135623730951*LamTD*vS*Conj(LamSD)*Sqr(vd) - 2*
      MuD*Conj(LamTD)*Sqr(vd) - 1.4142135623730951*LamSD*vS*Conj(LamTD)*Sqr(vd) -
      2*g2*Conj(MDWBT)*Sqr(vd) - 2*LamTD*Conj(MuD)*Sqr(vd) + 2*g2*MDWBT*Sqr(vu) -
      2*vT*AbsSqr(LamTU)*Sqr(vu) + 1.4142135623730951*LamTU*vS*Conj(LamSU)*Sqr(vu)
      + 2*MuU*Conj(LamTU)*Sqr(vu) + 1.4142135623730951*LamSU*vS*Conj(LamTU)*Sqr(
      vu) + 2*g2*Conj(MDWBT)*Sqr(vu) + 2*LamTU*Conj(MuU)*Sqr(vu)))/vT;
   mS2 = (0.025*(40*tadpole[3] - 160*vS*Sqr(MDBS) + 7.745966692414834*g1*MDBS*
      Sqr(vd) - 20*vS*AbsSqr(LamSD)*Sqr(vd) - 14.142135623730951*MuD*Conj(LamSD)*
      Sqr(vd) - 7.0710678118654755*LamTD*vT*Conj(LamSD)*Sqr(vd) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*Sqr(vd) + 7.745966692414834*g1*Conj(
      MDBS)*Sqr(vd) - 14.142135623730951*LamSD*Conj(MuD)*Sqr(vd) -
      7.745966692414834*g1*MDBS*Sqr(vu) - 20*vS*AbsSqr(LamSU)*Sqr(vu) -
      14.142135623730951*MuU*Conj(LamSU)*Sqr(vu) + 7.0710678118654755*LamTU*vT*
      Conj(LamSU)*Sqr(vu) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*Sqr(vu) -
      7.745966692414834*g1*Conj(MDBS)*Sqr(vu) - 14.142135623730951*LamSU*Conj(MuU)
      *Sqr(vu)))/vS;

   const bool is_finite = std::isfinite(mHd2) && std::isfinite(mHu2) &&
      std::isfinite(mT2) && std::isfinite(mS2);


   if (is_finite) {
      error = GSL_SUCCESS;
      ewsb_parameters[0] = mHd2;
      ewsb_parameters[1] = mHu2;
      ewsb_parameters[2] = mS2;
      ewsb_parameters[3] = mT2;

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
   MRSSM* model = ewsb_args->model;
   const unsigned ewsb_loop_order = ewsb_args->ewsb_loop_order;

   const double mHd2 = gsl_vector_get(x, 0);
   const double mHu2 = gsl_vector_get(x, 1);
   const double mS2 = gsl_vector_get(x, 2);
   const double mT2 = gsl_vector_get(x, 3);

   model->set_mHd2(mHd2);
   model->set_mHu2(mHu2);
   model->set_mS2(mS2);
   model->set_mT2(mT2);


   if (ewsb_loop_order > 0)
      model->calculate_DRbar_masses();

   double ewsb_parameters[number_of_ewsb_equations] =
      { mHd2, mHu2, mS2, mT2 };

   const int status = model->ewsb_step(ewsb_parameters);

   for (std::size_t i = 0; i < number_of_ewsb_equations; ++i)
      gsl_vector_set(f, i, ewsb_parameters[i]);

   return status;
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
   ostr << "MSOc = " << MSOc << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MRh = " << MRh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MRpm = " << MRpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha1 = " << MCha1.transpose() << '\n';
   ostr << "MCha2 = " << MCha2.transpose() << '\n';
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
   ostr << "ZHR = " << ZHR << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZRP = " << ZRP << '\n';
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
   const auto old_mS2 = mS2;
   const auto old_mT2 = mT2;

   solve_ewsb_tree_level_via_soft_higgs_masses();

   calculate_MVG();
   calculate_MVP();
   calculate_MVZ();
   calculate_MVWm();
   calculate_MGlu();
   calculate_MFv();
   calculate_MSOc();
   calculate_MSd();
   calculate_MSv();
   calculate_MSu();
   calculate_MSe();
   calculate_Mhh();
   calculate_MAh();
   calculate_MRh();
   calculate_MHpm();
   calculate_MRpm();
   calculate_MChi();
   calculate_MCha1();
   calculate_MCha2();
   calculate_MFe();
   calculate_MFd();
   calculate_MFu();

   mHd2 = old_mHd2;
   mHu2 = old_mHu2;
   mS2 = old_mS2;
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
   std::thread thread_MCha1(Thread(this, &CLASSNAME::calculate_MCha1_pole));
   std::thread thread_MCha2(Thread(this, &CLASSNAME::calculate_MCha2_pole));
   std::thread thread_MChi(Thread(this, &CLASSNAME::calculate_MChi_pole));
   std::thread thread_MGlu(Thread(this, &CLASSNAME::calculate_MGlu_pole));
   std::thread thread_Mhh(Thread(this, &CLASSNAME::calculate_Mhh_pole));
   std::thread thread_MHpm(Thread(this, &CLASSNAME::calculate_MHpm_pole));
   std::thread thread_MRh(Thread(this, &CLASSNAME::calculate_MRh_pole));
   std::thread thread_MRpm(Thread(this, &CLASSNAME::calculate_MRpm_pole));
   std::thread thread_MSd(Thread(this, &CLASSNAME::calculate_MSd_pole));
   std::thread thread_MSe(Thread(this, &CLASSNAME::calculate_MSe_pole));
   std::thread thread_MSOc(Thread(this, &CLASSNAME::calculate_MSOc_pole));
   std::thread thread_MSu(Thread(this, &CLASSNAME::calculate_MSu_pole));
   std::thread thread_MSv(Thread(this, &CLASSNAME::calculate_MSv_pole));

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
   thread_MCha1.join();
   thread_MCha2.join();
   thread_MChi.join();
   thread_MGlu.join();
   thread_Mhh.join();
   thread_MHpm.join();
   thread_MRh.join();
   thread_MRpm.join();
   thread_MSd.join();
   thread_MSe.join();
   thread_MSOc.join();
   thread_MSu.join();
   thread_MSv.join();


   if (thread_exception != 0)
      std::rethrow_exception(thread_exception);
#else
   calculate_MAh_pole();
   calculate_MCha1_pole();
   calculate_MCha2_pole();
   calculate_MChi_pole();
   calculate_MGlu_pole();
   calculate_Mhh_pole();
   calculate_MHpm_pole();
   calculate_MRh_pole();
   calculate_MRpm_pole();
   calculate_MSd_pole();
   calculate_MSe_pole();
   calculate_MSOc_pole();
   calculate_MSu_pole();
   calculate_MSv_pole();

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
   PHYSICAL(MSOc) = MSOc;
   PHYSICAL(MVP) = MVP;
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
   PHYSICAL(MRh) = MRh;
   PHYSICAL(ZHR) = ZHR;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MRpm) = MRpm;
   PHYSICAL(ZRP) = ZRP;
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
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MSOc = 0.;
   MVP = 0.;
   MVZ = 0.;
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
   MRpm = Eigen::Matrix<double,2,1>::Zero();
   ZRP = Eigen::Matrix<double,2,2>::Zero();
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


}

void CLASSNAME::clear_problems()
{
   problems.unflag_all_tachyons();
}

void CLASSNAME::clear()
{
   MRSSM_soft_parameters::clear();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
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
   Eigen::Array<double,3,1> MHpm_ChargedHiggs;
   Eigen::Array<double,1,1> MHpm_goldstone;

   MHpm_goldstone(0) = MVWm;

   remove_if_equal(MHpm, MHpm_goldstone, MHpm_ChargedHiggs);

   return MHpm_ChargedHiggs;
}

Eigen::Array<double,3,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,3,1> MAh_PseudoscalarHiggs;
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
   const double mass_matrix_VG = 0;

   return mass_matrix_VG;
}

void CLASSNAME::calculate_MVG()
{
   MVG = get_mass_matrix_VG();
}

double CLASSNAME::get_mass_matrix_Glu() const
{
   const double mass_matrix_Glu = MDGoc;

   return mass_matrix_Glu;
}

void CLASSNAME::calculate_MGlu()
{
   MGlu = get_mass_matrix_Glu();
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

double CLASSNAME::get_mass_matrix_SOc() const
{
   const double mass_matrix_SOc = moc2 + 2*Sqr(MDGoc);

   return mass_matrix_SOc;
}

void CLASSNAME::calculate_MSOc()
{
   MSOc = get_mass_matrix_SOc();

   if (MSOc < 0.)
      problems.flag_tachyon(MRSSM_info::SOc);

   MSOc = AbsSqrt(MSOc);
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
      problems.flag_tachyon(MRSSM_info::VZ);

   MVZ = AbsSqrt(MVZ);
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
   mass_matrix_Sd(1,0) = Conj(mass_matrix_Sd(0,1));
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
   mass_matrix_Sd(2,0) = Conj(mass_matrix_Sd(0,2));
   mass_matrix_Sd(2,1) = Conj(mass_matrix_Sd(1,2));
   mass_matrix_Sd(2,2) = 0.12909944487358055*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + mq2(2,2
      ) + 0.5*(AbsSqr(Yd(0,2)) + AbsSqr(Yd(1,2)) + AbsSqr(Yd(2,2)))*Sqr(vd) -
      0.025*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) + 0.025*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sd(2,3) = -0.7071067811865475*vu*Conj(Yd(0,2))*Mu;
   mass_matrix_Sd(2,4) = -0.7071067811865475*vu*Conj(Yd(1,2))*Mu;
   mass_matrix_Sd(2,5) = -0.7071067811865475*vu*Conj(Yd(2,2))*Mu;
   mass_matrix_Sd(3,0) = Conj(mass_matrix_Sd(0,3));
   mass_matrix_Sd(3,1) = Conj(mass_matrix_Sd(1,3));
   mass_matrix_Sd(3,2) = Conj(mass_matrix_Sd(2,3));
   mass_matrix_Sd(3,3) = 0.2581988897471611*g1*MDBS*vS +
      0.2581988897471611*g1*vS*Conj(MDBS) + md2(0,0) + 0.5*(AbsSqr(Yd(0,0)) +
      AbsSqr(Yd(0,1)) + AbsSqr(Yd(0,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(3,4) = md2(0,1) + 0.5*Sqr(vd)*(Conj(Yd(1,0))*Yd(0,0) +
      Conj(Yd(1,1))*Yd(0,1) + Conj(Yd(1,2))*Yd(0,2));
   mass_matrix_Sd(3,5) = md2(0,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(0,0) +
      Conj(Yd(2,1))*Yd(0,1) + Conj(Yd(2,2))*Yd(0,2));
   mass_matrix_Sd(4,0) = Conj(mass_matrix_Sd(0,4));
   mass_matrix_Sd(4,1) = Conj(mass_matrix_Sd(1,4));
   mass_matrix_Sd(4,2) = Conj(mass_matrix_Sd(2,4));
   mass_matrix_Sd(4,3) = Conj(mass_matrix_Sd(3,4));
   mass_matrix_Sd(4,4) = 0.2581988897471611*g1*MDBS*vS +
      0.2581988897471611*g1*vS*Conj(MDBS) + md2(1,1) + 0.5*(AbsSqr(Yd(1,0)) +
      AbsSqr(Yd(1,1)) + AbsSqr(Yd(1,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Sd(4,5) = md2(1,2) + 0.5*Sqr(vd)*(Conj(Yd(2,0))*Yd(1,0) +
      Conj(Yd(2,1))*Yd(1,1) + Conj(Yd(2,2))*Yd(1,2));
   mass_matrix_Sd(5,0) = Conj(mass_matrix_Sd(0,5));
   mass_matrix_Sd(5,1) = Conj(mass_matrix_Sd(1,5));
   mass_matrix_Sd(5,2) = Conj(mass_matrix_Sd(2,5));
   mass_matrix_Sd(5,3) = Conj(mass_matrix_Sd(3,5));
   mass_matrix_Sd(5,4) = Conj(mass_matrix_Sd(4,5));
   mass_matrix_Sd(5,5) = 0.2581988897471611*g1*MDBS*vS +
      0.2581988897471611*g1*vS*Conj(MDBS) + md2(2,2) + 0.5*(AbsSqr(Yd(2,0)) +
      AbsSqr(Yd(2,1)) + AbsSqr(Yd(2,2)))*Sqr(vd) - 0.05*Sqr(g1)*Sqr(vd) + 0.05*
      Sqr(g1)*Sqr(vu);

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

   if (MSd.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Sd);

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
   mass_matrix_Sv(1,0) = Conj(mass_matrix_Sv(0,1));
   mass_matrix_Sv(1,1) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(1,1)
      + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu)
      - 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Sv(1,2) = ml2(1,2);
   mass_matrix_Sv(2,0) = Conj(mass_matrix_Sv(0,2));
   mass_matrix_Sv(2,1) = Conj(mass_matrix_Sv(1,2));
   mass_matrix_Sv(2,2) = -0.3872983346207417*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + ml2(2,2)
      + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu)
      - 0.125*Sqr(g2)*Sqr(vu);

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

   if (MSv.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Sv);

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
   mass_matrix_Su(1,0) = Conj(mass_matrix_Su(0,1));
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
   mass_matrix_Su(2,0) = Conj(mass_matrix_Su(0,2));
   mass_matrix_Su(2,1) = Conj(mass_matrix_Su(1,2));
   mass_matrix_Su(2,2) = 0.12909944487358055*g1*MDBS*vS + 0.5*g2*MDWBT*vT
      + 0.12909944487358055*g1*vS*Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + mq2(2,2
      ) - 0.025*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*(AbsSqr(Yu(0,2))
      + AbsSqr(Yu(1,2)) + AbsSqr(Yu(2,2)))*Sqr(vu) + 0.025*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Su(2,3) = -0.7071067811865475*vd*Conj(Yu(0,2))*Mu;
   mass_matrix_Su(2,4) = -0.7071067811865475*vd*Conj(Yu(1,2))*Mu;
   mass_matrix_Su(2,5) = -0.7071067811865475*vd*Conj(Yu(2,2))*Mu;
   mass_matrix_Su(3,0) = Conj(mass_matrix_Su(0,3));
   mass_matrix_Su(3,1) = Conj(mass_matrix_Su(1,3));
   mass_matrix_Su(3,2) = Conj(mass_matrix_Su(2,3));
   mass_matrix_Su(3,3) = -0.5163977794943222*g1*MDBS*vS -
      0.5163977794943222*g1*vS*Conj(MDBS) + mu2(0,0) + 0.1*Sqr(g1)*Sqr(vd) +
      0.5*(AbsSqr(Yu(0,0)) + AbsSqr(Yu(0,1)) + AbsSqr(Yu(0,2)))*Sqr(vu) - 0.1*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Su(3,4) = mu2(0,1) + 0.5*Sqr(vu)*(Conj(Yu(1,0))*Yu(0,0) +
      Conj(Yu(1,1))*Yu(0,1) + Conj(Yu(1,2))*Yu(0,2));
   mass_matrix_Su(3,5) = mu2(0,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(0,0) +
      Conj(Yu(2,1))*Yu(0,1) + Conj(Yu(2,2))*Yu(0,2));
   mass_matrix_Su(4,0) = Conj(mass_matrix_Su(0,4));
   mass_matrix_Su(4,1) = Conj(mass_matrix_Su(1,4));
   mass_matrix_Su(4,2) = Conj(mass_matrix_Su(2,4));
   mass_matrix_Su(4,3) = Conj(mass_matrix_Su(3,4));
   mass_matrix_Su(4,4) = -0.5163977794943222*g1*MDBS*vS -
      0.5163977794943222*g1*vS*Conj(MDBS) + mu2(1,1) + 0.1*Sqr(g1)*Sqr(vd) +
      0.5*(AbsSqr(Yu(1,0)) + AbsSqr(Yu(1,1)) + AbsSqr(Yu(1,2)))*Sqr(vu) - 0.1*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Su(4,5) = mu2(1,2) + 0.5*Sqr(vu)*(Conj(Yu(2,0))*Yu(1,0) +
      Conj(Yu(2,1))*Yu(1,1) + Conj(Yu(2,2))*Yu(1,2));
   mass_matrix_Su(5,0) = Conj(mass_matrix_Su(0,5));
   mass_matrix_Su(5,1) = Conj(mass_matrix_Su(1,5));
   mass_matrix_Su(5,2) = Conj(mass_matrix_Su(2,5));
   mass_matrix_Su(5,3) = Conj(mass_matrix_Su(3,5));
   mass_matrix_Su(5,4) = Conj(mass_matrix_Su(4,5));
   mass_matrix_Su(5,5) = -0.5163977794943222*g1*MDBS*vS -
      0.5163977794943222*g1*vS*Conj(MDBS) + mu2(2,2) + 0.1*Sqr(g1)*Sqr(vd) +
      0.5*(AbsSqr(Yu(2,0)) + AbsSqr(Yu(2,1)) + AbsSqr(Yu(2,2)))*Sqr(vu) - 0.1*
      Sqr(g1)*Sqr(vu);

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

   if (MSu.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Su);

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
   mass_matrix_Se(1,0) = Conj(mass_matrix_Se(0,1));
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
   mass_matrix_Se(2,0) = Conj(mass_matrix_Se(0,2));
   mass_matrix_Se(2,1) = Conj(mass_matrix_Se(1,2));
   mass_matrix_Se(2,2) = -0.3872983346207417*g1*MDBS*vS - 0.5*g2*MDWBT*vT
      - 0.3872983346207417*g1*vS*Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + ml2(2,2)
      + 0.5*(AbsSqr(Ye(0,2)) + AbsSqr(Ye(1,2)) + AbsSqr(Ye(2,2)))*Sqr(vd) +
      0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd) - 0.075*Sqr(g1)*Sqr(vu) +
      0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Se(2,3) = -0.7071067811865475*vu*Conj(Ye(0,2))*Mu;
   mass_matrix_Se(2,4) = -0.7071067811865475*vu*Conj(Ye(1,2))*Mu;
   mass_matrix_Se(2,5) = -0.7071067811865475*vu*Conj(Ye(2,2))*Mu;
   mass_matrix_Se(3,0) = Conj(mass_matrix_Se(0,3));
   mass_matrix_Se(3,1) = Conj(mass_matrix_Se(1,3));
   mass_matrix_Se(3,2) = Conj(mass_matrix_Se(2,3));
   mass_matrix_Se(3,3) = 0.7745966692414834*g1*MDBS*vS +
      0.7745966692414834*g1*vS*Conj(MDBS) + me2(0,0) + 0.5*(AbsSqr(Ye(0,0)) +
      AbsSqr(Ye(0,1)) + AbsSqr(Ye(0,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Se(3,4) = me2(0,1) + 0.5*Sqr(vd)*(Conj(Ye(1,0))*Ye(0,0) +
      Conj(Ye(1,1))*Ye(0,1) + Conj(Ye(1,2))*Ye(0,2));
   mass_matrix_Se(3,5) = me2(0,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(0,0) +
      Conj(Ye(2,1))*Ye(0,1) + Conj(Ye(2,2))*Ye(0,2));
   mass_matrix_Se(4,0) = Conj(mass_matrix_Se(0,4));
   mass_matrix_Se(4,1) = Conj(mass_matrix_Se(1,4));
   mass_matrix_Se(4,2) = Conj(mass_matrix_Se(2,4));
   mass_matrix_Se(4,3) = Conj(mass_matrix_Se(3,4));
   mass_matrix_Se(4,4) = 0.7745966692414834*g1*MDBS*vS +
      0.7745966692414834*g1*vS*Conj(MDBS) + me2(1,1) + 0.5*(AbsSqr(Ye(1,0)) +
      AbsSqr(Ye(1,1)) + AbsSqr(Ye(1,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*
      Sqr(g1)*Sqr(vu);
   mass_matrix_Se(4,5) = me2(1,2) + 0.5*Sqr(vd)*(Conj(Ye(2,0))*Ye(1,0) +
      Conj(Ye(2,1))*Ye(1,1) + Conj(Ye(2,2))*Ye(1,2));
   mass_matrix_Se(5,0) = Conj(mass_matrix_Se(0,5));
   mass_matrix_Se(5,1) = Conj(mass_matrix_Se(1,5));
   mass_matrix_Se(5,2) = Conj(mass_matrix_Se(2,5));
   mass_matrix_Se(5,3) = Conj(mass_matrix_Se(3,5));
   mass_matrix_Se(5,4) = Conj(mass_matrix_Se(4,5));
   mass_matrix_Se(5,5) = 0.7745966692414834*g1*MDBS*vS +
      0.7745966692414834*g1*vS*Conj(MDBS) + me2(2,2) + 0.5*(AbsSqr(Ye(2,0)) +
      AbsSqr(Ye(2,1)) + AbsSqr(Ye(2,2)))*Sqr(vd) - 0.15*Sqr(g1)*Sqr(vd) + 0.15*
      Sqr(g1)*Sqr(vu);

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

   if (MSe.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Se);

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
   mass_matrix_hh(1,0) = mass_matrix_hh(0,1);
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
   mass_matrix_hh(2,0) = mass_matrix_hh(0,2);
   mass_matrix_hh(2,1) = mass_matrix_hh(1,2);
   mass_matrix_hh(2,2) = mS2 + 4*Sqr(MDBS) + 0.5*AbsSqr(LamSD)*Sqr(vd) +
      0.5*AbsSqr(LamSU)*Sqr(vu);
   mass_matrix_hh(2,3) = 0.17677669529663687*LamTD*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*Conj(LamTD)*Sqr(vd) - 0.17677669529663687*LamTU
      *Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*Conj(LamTU)*Sqr(vu);
   mass_matrix_hh(3,0) = mass_matrix_hh(0,3);
   mass_matrix_hh(3,1) = mass_matrix_hh(1,3);
   mass_matrix_hh(3,2) = mass_matrix_hh(2,3);
   mass_matrix_hh(3,3) = mT2 + 4*Sqr(MDWBT) + 0.25*AbsSqr(LamTD)*Sqr(vd)
      + 0.25*AbsSqr(LamTU)*Sqr(vu);

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

   if (Mhh.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::hh);

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
   mass_matrix_Ah(1,0) = mass_matrix_Ah(0,1);
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
   mass_matrix_Ah(2,0) = mass_matrix_Ah(0,2);
   mass_matrix_Ah(2,1) = mass_matrix_Ah(1,2);
   mass_matrix_Ah(2,2) = mS2 + 0.5*AbsSqr(LamSD)*Sqr(vd) + 0.5*AbsSqr(
      LamSU)*Sqr(vu);
   mass_matrix_Ah(2,3) = 0.17677669529663687*LamTD*Conj(LamSD)*Sqr(vd) +
      0.17677669529663687*LamSD*Conj(LamTD)*Sqr(vd) - 0.17677669529663687*LamTU
      *Conj(LamSU)*Sqr(vu) - 0.17677669529663687*LamSU*Conj(LamTU)*Sqr(vu);
   mass_matrix_Ah(3,0) = mass_matrix_Ah(0,3);
   mass_matrix_Ah(3,1) = mass_matrix_Ah(1,3);
   mass_matrix_Ah(3,2) = mass_matrix_Ah(2,3);
   mass_matrix_Ah(3,3) = mT2 + 0.25*AbsSqr(LamTD)*Sqr(vd) + 0.25*AbsSqr(
      LamTU)*Sqr(vu);

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

   if (MAh.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Ah);

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
   mass_matrix_Rh(1,0) = Conj(mass_matrix_Rh(0,1));
   mass_matrix_Rh(1,1) = mRu2 - 0.3872983346207417*g1*MDBS*vS + 0.5*g2*
      MDWBT*vT + AbsSqr(MuU) + 0.7071067811865475*MuU*vS*Conj(LamSU) -
      0.35355339059327373*LamTU*vS*vT*Conj(LamSU) - 0.5*MuU*vT*Conj(LamTU) -
      0.35355339059327373*LamSU*vS*vT*Conj(LamTU) - 0.3872983346207417*g1*vS*
      Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*Conj(MuU
      ) - 0.5*LamTU*vT*Conj(MuU) + 0.075*Sqr(g1)*Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) + 0.5*AbsSqr(
      LamSU)*Sqr(vu) + 0.25*AbsSqr(LamTU)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) -
      0.125*Sqr(g2)*Sqr(vu);

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

   if (MRh.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Rh);

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

   if (MHpm.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Hpm);

   MHpm = AbsSqrt(MHpm);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Rpm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Rpm;

   mass_matrix_Rpm(0,0) = mRu2 - 0.3872983346207417*g1*MDBS*vS - 0.5*g2*
      MDWBT*vT + AbsSqr(MuU) + 0.7071067811865475*MuU*vS*Conj(LamSU) +
      0.35355339059327373*LamTU*vS*vT*Conj(LamSU) + 0.5*MuU*vT*Conj(LamTU) +
      0.35355339059327373*LamSU*vS*vT*Conj(LamTU) - 0.3872983346207417*g1*vS*
      Conj(MDBS) - 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSU*vS*Conj(MuU
      ) + 0.5*LamTU*vT*Conj(MuU) + 0.075*Sqr(g1)*Sqr(vd) - 0.125*Sqr(g2)*Sqr(vd
      ) + 0.5*AbsSqr(LamSU)*Sqr(vS) + 0.25*AbsSqr(LamTU)*Sqr(vT) + 0.5*AbsSqr(
      LamTU)*Sqr(vu) - 0.075*Sqr(g1)*Sqr(vu) + 0.125*Sqr(g2)*Sqr(vu);
   mass_matrix_Rpm(0,1) = 0;
   mass_matrix_Rpm(1,0) = mass_matrix_Rpm(0,1);
   mass_matrix_Rpm(1,1) = mRd2 + 0.3872983346207417*g1*MDBS*vS + 0.5*g2*
      MDWBT*vT + AbsSqr(MuD) + 0.7071067811865475*MuD*vS*Conj(LamSD) -
      0.35355339059327373*LamTD*vS*vT*Conj(LamSD) - 0.5*MuD*vT*Conj(LamTD) -
      0.35355339059327373*LamSD*vS*vT*Conj(LamTD) + 0.3872983346207417*g1*vS*
      Conj(MDBS) + 0.5*g2*vT*Conj(MDWBT) + 0.7071067811865475*LamSD*vS*Conj(MuD
      ) - 0.5*LamTD*vT*Conj(MuD) + 0.5*AbsSqr(LamTD)*Sqr(vd) - 0.075*Sqr(g1)*
      Sqr(vd) + 0.125*Sqr(g2)*Sqr(vd) + 0.5*AbsSqr(LamSD)*Sqr(vS) + 0.25*AbsSqr
      (LamTD)*Sqr(vT) + 0.075*Sqr(g1)*Sqr(vu) - 0.125*Sqr(g2)*Sqr(vu);

   return mass_matrix_Rpm;
}

void CLASSNAME::calculate_MRpm()
{
   const auto mass_matrix_Rpm(get_mass_matrix_Rpm());

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_Rpm, MRpm, ZRP, eigenvalue_error)
      ;
   problems.flag_bad_mass(MRSSM_info::Rpm, eigenvalue_error > precision *
      Abs(MRpm(0)));
#else
   fs_diagonalize_hermitian(mass_matrix_Rpm, MRpm, ZRP);
#endif

   if (MRpm.minCoeff() < 0.)
      problems.flag_tachyon(MRSSM_info::Rpm);

   MRpm = AbsSqrt(MRpm);
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
   const double mass_matrix_VWm = 0.25*Sqr(g2)*(Sqr(vd) + 4*Sqr(vT) + Sqr
      (vu));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   MVWm = get_mass_matrix_VWm();

   if (MVWm < 0.)
      problems.flag_tachyon(MRSSM_info::VWm);

   MVWm = AbsSqrt(MVWm);
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = mHd2*vd - 0.3872983346207417*g1*MDBS*vd*vS + 0.5*g2*MDWBT*vd
      *vT + vd*AbsSqr(MuD) + vd*AbsSqr(Mu) - 0.5*vu*BMu + 0.7071067811865475*MuD*
      vd*vS*Conj(LamSD) + 0.35355339059327373*LamTD*vd*vS*vT*Conj(LamSD) + 0.5*MuD
      *vd*vT*Conj(LamTD) + 0.35355339059327373*LamSD*vd*vS*vT*Conj(LamTD) -
      0.3872983346207417*g1*vd*vS*Conj(MDBS) + 0.5*g2*vd*vT*Conj(MDWBT) +
      0.7071067811865475*LamSD*vd*vS*Conj(MuD) + 0.5*LamTD*vd*vT*Conj(MuD) - 0.5*
      vu*Conj(BMu) + 0.075*Power(vd,3)*Sqr(g1) + 0.125*Power(vd,3)*Sqr(g2) + 0.5*
      vd*AbsSqr(LamSD)*Sqr(vS) + 0.25*vd*AbsSqr(LamTD)*Sqr(vT) - 0.075*vd*Sqr(g1)*
      Sqr(vu) - 0.125*vd*Sqr(g2)*Sqr(vu);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = mHu2*vu + 0.3872983346207417*g1*MDBS*vS*vu - 0.5*g2*MDWBT*vT
      *vu + vu*AbsSqr(MuU) + vu*AbsSqr(Mu) - 0.5*vd*BMu + 0.7071067811865475*MuU*
      vS*vu*Conj(LamSU) - 0.35355339059327373*LamTU*vS*vT*vu*Conj(LamSU) - 0.5*MuU
      *vT*vu*Conj(LamTU) - 0.35355339059327373*LamSU*vS*vT*vu*Conj(LamTU) +
      0.3872983346207417*g1*vS*vu*Conj(MDBS) - 0.5*g2*vT*vu*Conj(MDWBT) +
      0.7071067811865475*LamSU*vS*vu*Conj(MuU) - 0.5*LamTU*vT*vu*Conj(MuU) - 0.5*
      vd*Conj(BMu) + 0.075*Power(vu,3)*Sqr(g1) + 0.125*Power(vu,3)*Sqr(g2) - 0.075
      *vu*Sqr(g1)*Sqr(vd) - 0.125*vu*Sqr(g2)*Sqr(vd) + 0.5*vu*AbsSqr(LamSU)*Sqr(vS
      ) + 0.25*vu*AbsSqr(LamTU)*Sqr(vT);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_3() const
{
   double result = mT2*vT + 4*vT*Sqr(MDWBT) + 0.25*g2*MDWBT*Sqr(vd) + 0.25*vT*
      AbsSqr(LamTD)*Sqr(vd) + 0.17677669529663687*LamTD*vS*Conj(LamSD)*Sqr(vd) +
      0.25*MuD*Conj(LamTD)*Sqr(vd) + 0.17677669529663687*LamSD*vS*Conj(LamTD)*Sqr(
      vd) + 0.25*g2*Conj(MDWBT)*Sqr(vd) + 0.25*LamTD*Conj(MuD)*Sqr(vd) - 0.25*g2*
      MDWBT*Sqr(vu) + 0.25*vT*AbsSqr(LamTU)*Sqr(vu) - 0.17677669529663687*LamTU*vS
      *Conj(LamSU)*Sqr(vu) - 0.25*MuU*Conj(LamTU)*Sqr(vu) - 0.17677669529663687*
      LamSU*vS*Conj(LamTU)*Sqr(vu) - 0.25*g2*Conj(MDWBT)*Sqr(vu) - 0.25*LamTU*Conj
      (MuU)*Sqr(vu);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_4() const
{
   double result = mS2*vS + 4*vS*Sqr(MDBS) - 0.19364916731037085*g1*MDBS*Sqr(vd
      ) + 0.5*vS*AbsSqr(LamSD)*Sqr(vd) + 0.35355339059327373*MuD*Conj(LamSD)*Sqr(
      vd) + 0.17677669529663687*LamTD*vT*Conj(LamSD)*Sqr(vd) + 0.17677669529663687
      *LamSD*vT*Conj(LamTD)*Sqr(vd) - 0.19364916731037085*g1*Conj(MDBS)*Sqr(vd) +
      0.35355339059327373*LamSD*Conj(MuD)*Sqr(vd) + 0.19364916731037085*g1*MDBS*
      Sqr(vu) + 0.5*vS*AbsSqr(LamSU)*Sqr(vu) + 0.35355339059327373*MuU*Conj(LamSU)
      *Sqr(vu) - 0.17677669529663687*LamTU*vT*Conj(LamSU)*Sqr(vu) -
      0.17677669529663687*LamSU*vT*Conj(LamTU)*Sqr(vu) + 0.19364916731037085*g1*
      Conj(MDBS)*Sqr(vu) + 0.35355339059327373*LamSU*Conj(MuU)*Sqr(vu);

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

double CLASSNAME::CpUSdconjUSdconjSOcSOc(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

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

std::complex<double> CLASSNAME::CpUSdconjUSdconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2;
   std::complex<double> tmp_3;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_3 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_2 += tmp_3;
   result += (-0.1*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0)) * tmp_2;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   std::complex<double> tmp_4;
   std::complex<double> tmp_5;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_5 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_4 += tmp_5;
   result += (0.1*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1)) * tmp_4;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_6;
   std::complex<double> tmp_7;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_7 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_6 += tmp_7;
   result += (0.1*Conj(ZRP(gI2,0))*Sqr(g1)*ZRP(gI1,0)) * tmp_6;
   if (gO1 < 3) {
      result += 0.05*Conj(ZRP(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZRP(
         gI1,0);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZRP(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZRP(
         gI1,0);
   }
   std::complex<double> tmp_8;
   std::complex<double> tmp_9;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_9 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_8 += tmp_9;
   result += (-0.1*Conj(ZRP(gI2,1))*Sqr(g1)*ZRP(gI1,1)) * tmp_8;
   if (gO1 < 3) {
      result += -0.05*Conj(ZRP(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZRP(
         gI1,1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZRP(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZRP(
         gI1,1);
   }

   return result;
}

double CLASSNAME::CpconjUSdbarCha1FuPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdbarCha1FuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_10;
   std::complex<double> tmp_11;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_12;
      std::complex<double> tmp_13;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_13 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_12 += tmp_13;
      tmp_11 += (Conj(ZUL(gI2,j2))) * tmp_12;
   }
   tmp_10 += tmp_11;
   result += (Conj(UM1(gI1,1))) * tmp_10;

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_14;
   std::complex<double> tmp_15;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_15 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_14 += tmp_15;
   result += (0.1*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_14;
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

std::complex<double> CLASSNAME::CpconjUSdFuCha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_16;
      std::complex<double> tmp_17;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_17 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_16 += tmp_17;
      result += (UP2(gI2,1)) * tmp_16;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFuCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(UM2(gI2,0))*Conj(ZUL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_18;
      std::complex<double> tmp_19;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_19 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_18 += tmp_19;
      result += (-ZN2(gI2,2)) * tmp_18;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdFdChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZDL(gI1,gO1))*Conj(ZN1(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZDL(gI1,gO1))*Conj(ZN1(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_20;
   std::complex<double> tmp_21;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_21 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_20 += tmp_21;
   result += (0.1*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_20;
   std::complex<double> tmp_22;
   std::complex<double> tmp_23;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_24;
      std::complex<double> tmp_25;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_26;
         std::complex<double> tmp_27;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_27 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_26 += tmp_27;
         tmp_25 += (KroneckerDelta(gO2,3 + j2)) * tmp_26;
      }
      tmp_24 += tmp_25;
      tmp_23 += (KroneckerDelta(gO1,3 + j3)) * tmp_24;
   }
   tmp_22 += tmp_23;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_22;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_28;
      std::complex<double> tmp_29;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_29 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_28 += tmp_29;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_28;
   }
   std::complex<double> tmp_30;
   std::complex<double> tmp_31;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_31 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_30 += tmp_31;
   result += (-0.1*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_30;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_32;
   std::complex<double> tmp_33;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_33 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_32 += tmp_33;
   result += (0.1*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_32;
   std::complex<double> tmp_34;
   std::complex<double> tmp_35;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_36;
      std::complex<double> tmp_37;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_38;
         std::complex<double> tmp_39;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_39 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_38 += tmp_39;
         tmp_37 += (KroneckerDelta(gO2,3 + j2)) * tmp_38;
      }
      tmp_36 += tmp_37;
      tmp_35 += (KroneckerDelta(gO1,3 + j3)) * tmp_36;
   }
   tmp_34 += tmp_35;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_34;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_40;
   std::complex<double> tmp_41;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_41 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_40 += tmp_41;
   result += (-0.1*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_40;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_42;
      std::complex<double> tmp_43;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_43 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_42 += tmp_43;
      result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_42;
   }
   if (gO1 < 3) {
      result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_44;
   std::complex<double> tmp_45;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_45 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_44 += tmp_45;
   result += (0.1*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_44;
   std::complex<double> tmp_46;
   std::complex<double> tmp_47;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_48;
      std::complex<double> tmp_49;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_50;
         std::complex<double> tmp_51;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_51 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_50 += tmp_51;
         tmp_49 += (KroneckerDelta(gO2,3 + j2)) * tmp_50;
      }
      tmp_48 += tmp_49;
      tmp_47 += (KroneckerDelta(gO1,3 + j3)) * tmp_48;
   }
   tmp_46 += tmp_47;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_46;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_52;
      std::complex<double> tmp_53;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_53 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_52 += tmp_53;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_52;
   }
   std::complex<double> tmp_54;
   std::complex<double> tmp_55;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_55 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_54 += tmp_55;
   result += (-0.1*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_54;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdbarChiFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_56;
   std::complex<double> tmp_57;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_57 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_56 += tmp_57;
   result += (-0.3651483716701107*g1*ZN1(gI1,0)) * tmp_56;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdbarChiFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_58;
   std::complex<double> tmp_59;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_60;
      std::complex<double> tmp_61;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_61 += KroneckerDelta(gO1,3 + j1)*Yd(j1,j2);
      }
      tmp_60 += tmp_61;
      tmp_59 += (Conj(ZDL(gI2,j2))) * tmp_60;
   }
   tmp_58 += tmp_59;
   result += (-Conj(ZN2(gI1,2))) * tmp_58;

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_62;
   std::complex<double> tmp_64;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_64 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_62 += tmp_64;
   std::complex<double> tmp_63;
   std::complex<double> tmp_65;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_65 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_63 += tmp_65;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_62 * tmp_63;
   std::complex<double> tmp_66;
   std::complex<double> tmp_68;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_68 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
   }
   tmp_66 += tmp_68;
   std::complex<double> tmp_67;
   std::complex<double> tmp_69;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_69 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_67 += tmp_69;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_66 * tmp_67;
   std::complex<double> tmp_70;
   std::complex<double> tmp_72;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_72 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_70 += tmp_72;
   std::complex<double> tmp_71;
   std::complex<double> tmp_73;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_73 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_71 += tmp_73;
   result += (-0.05*Sqr(g1)) * tmp_70 * tmp_71;
   std::complex<double> tmp_74;
   std::complex<double> tmp_76;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_76 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_74 += tmp_76;
   std::complex<double> tmp_75;
   std::complex<double> tmp_77;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_77 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_75 += tmp_77;
   result += (-0.1*Sqr(g1)) * tmp_74 * tmp_75;
   std::complex<double> tmp_78;
   std::complex<double> tmp_80;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_80 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_78 += tmp_80;
   std::complex<double> tmp_79;
   std::complex<double> tmp_81;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_81 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_79 += tmp_81;
   result += (-0.05*Sqr(g1)) * tmp_78 * tmp_79;
   std::complex<double> tmp_82;
   std::complex<double> tmp_84;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_84 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_82 += tmp_84;
   std::complex<double> tmp_83;
   std::complex<double> tmp_85;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_85 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_83 += tmp_85;
   result += (-0.1*Sqr(g1)) * tmp_82 * tmp_83;
   std::complex<double> tmp_86;
   std::complex<double> tmp_88;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_88 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_86 += tmp_88;
   std::complex<double> tmp_87;
   std::complex<double> tmp_89;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_89 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_87 += tmp_89;
   result += (-0.03333333333333333*Sqr(g1)) * tmp_86 * tmp_87;
   std::complex<double> tmp_90;
   std::complex<double> tmp_92;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_92 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_90 += tmp_92;
   std::complex<double> tmp_91;
   std::complex<double> tmp_93;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_93 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
   }
   tmp_91 += tmp_93;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_90 * tmp_91;
   std::complex<double> tmp_94;
   std::complex<double> tmp_96;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_97;
      std::complex<double> tmp_98;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_98 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_97 += tmp_98;
      tmp_96 += (Conj(ZD(gI2,j2))) * tmp_97;
   }
   tmp_94 += tmp_96;
   std::complex<double> tmp_95;
   std::complex<double> tmp_99;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_100;
      std::complex<double> tmp_101;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_101 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_100 += tmp_101;
      tmp_99 += (ZD(gI1,j4)) * tmp_100;
   }
   tmp_95 += tmp_99;
   result += (-1) * tmp_94 * tmp_95;
   if (gO1 < 3) {
      std::complex<double> tmp_102;
      std::complex<double> tmp_103;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_103 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_102 += tmp_103;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_102;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_104;
      std::complex<double> tmp_105;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_105 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_104 += tmp_105;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_104;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_106;
      std::complex<double> tmp_107;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_107 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_106 += tmp_107;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_106;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_108;
      std::complex<double> tmp_109;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_109 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_108 += tmp_109;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_108;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_110;
      std::complex<double> tmp_111;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_111 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_110 += tmp_111;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_110;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_112;
      std::complex<double> tmp_113;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_113 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_112 += tmp_113;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_112;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_114;
      std::complex<double> tmp_116;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_116 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_114 += tmp_116;
      std::complex<double> tmp_115;
      std::complex<double> tmp_117;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_118;
         std::complex<double> tmp_119;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_119 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_118 += tmp_119;
         tmp_117 += (ZD(gI1,j4)) * tmp_118;
      }
      tmp_115 += tmp_117;
      result += (-3) * tmp_114 * tmp_115;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_120;
      std::complex<double> tmp_121;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_121 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_120 += tmp_121;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_120;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_122;
      std::complex<double> tmp_123;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_123 += KroneckerDelta(gO1,3 + j1)*ZD(gI1,3 + j1);
      }
      tmp_122 += tmp_123;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_122;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_124;
      std::complex<double> tmp_125;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_125 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_124 += tmp_125;
      result += (-0.016666666666666666*Conj(ZD(gI2,gO2))*Sqr(g1)) * tmp_124;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_126;
      std::complex<double> tmp_127;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_127 += KroneckerDelta(gO1,3 + j2)*ZD(gI1,3 + j2);
      }
      tmp_126 += tmp_127;
      result += (0.6666666666666666*Conj(ZD(gI2,gO2))*Sqr(g3)) * tmp_126;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_128;
      std::complex<double> tmp_130;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_131;
         std::complex<double> tmp_132;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_132 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_131 += tmp_132;
         tmp_130 += (Conj(ZD(gI2,j2))) * tmp_131;
      }
      tmp_128 += tmp_130;
      std::complex<double> tmp_129;
      std::complex<double> tmp_133;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_133 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_129 += tmp_133;
      result += (-3) * tmp_128 * tmp_129;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_134;
      std::complex<double> tmp_136;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_136 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_134 += tmp_136;
      std::complex<double> tmp_135;
      std::complex<double> tmp_137;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_137 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_135 += tmp_137;
      result += (-1) * tmp_134 * tmp_135;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_138;
      std::complex<double> tmp_139;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_139 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_138 += tmp_139;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_138;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_140;
      std::complex<double> tmp_141;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_141 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_140 += tmp_141;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_140;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_142;
      std::complex<double> tmp_143;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_143 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_142 += tmp_143;
      result += (-0.016666666666666666*Sqr(g1)*ZD(gI1,gO1)) * tmp_142;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_144;
      std::complex<double> tmp_145;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_145 += Conj(ZD(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_144 += tmp_145;
      result += (0.6666666666666666*Sqr(g3)*ZD(gI1,gO1)) * tmp_144;
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

   std::complex<double> tmp_146;
   std::complex<double> tmp_148;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_148 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_146 += tmp_148;
   std::complex<double> tmp_147;
   std::complex<double> tmp_149;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_149 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_147 += tmp_149;
   result += (0.05*Sqr(g1)) * tmp_146 * tmp_147;
   std::complex<double> tmp_150;
   std::complex<double> tmp_152;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_152 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_150 += tmp_152;
   std::complex<double> tmp_151;
   std::complex<double> tmp_153;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_153 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_151 += tmp_153;
   result += (-0.1*Sqr(g1)) * tmp_150 * tmp_151;
   std::complex<double> tmp_154;
   std::complex<double> tmp_156;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_156 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_154 += tmp_156;
   std::complex<double> tmp_155;
   std::complex<double> tmp_157;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_157 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_155 += tmp_157;
   result += (0.05*Sqr(g1)) * tmp_154 * tmp_155;
   std::complex<double> tmp_158;
   std::complex<double> tmp_160;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_160 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_158 += tmp_160;
   std::complex<double> tmp_159;
   std::complex<double> tmp_161;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_161 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_159 += tmp_161;
   result += (-0.1*Sqr(g1)) * tmp_158 * tmp_159;
   if (gO1 < 3) {
      std::complex<double> tmp_162;
      std::complex<double> tmp_163;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_163 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_162 += tmp_163;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_162;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_164;
      std::complex<double> tmp_165;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_165 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_164 += tmp_165;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_164;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_166;
      std::complex<double> tmp_167;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_167 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_166 += tmp_167;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_166;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_168;
      std::complex<double> tmp_169;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_169 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_168 += tmp_169;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_168;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_170;
      std::complex<double> tmp_171;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_171 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_170 += tmp_171;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_170;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_172;
      std::complex<double> tmp_173;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_173 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_172 += tmp_173;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_172;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_174;
      std::complex<double> tmp_176;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_176 += KroneckerDelta(gO2,3 + j1)*Yd(j1,gO1);
      }
      tmp_174 += tmp_176;
      std::complex<double> tmp_175;
      std::complex<double> tmp_177;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_178;
         std::complex<double> tmp_179;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_179 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_178 += tmp_179;
         tmp_177 += (ZE(gI1,j4)) * tmp_178;
      }
      tmp_175 += tmp_177;
      result += (-1) * tmp_174 * tmp_175;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_180;
      std::complex<double> tmp_182;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_183;
         std::complex<double> tmp_184;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_184 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_183 += tmp_184;
         tmp_182 += (Conj(ZE(gI2,j2))) * tmp_183;
      }
      tmp_180 += tmp_182;
      std::complex<double> tmp_181;
      std::complex<double> tmp_185;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_185 += Conj(Yd(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_181 += tmp_185;
      result += (-1) * tmp_180 * tmp_181;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSdconjUSdconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_186;
   std::complex<double> tmp_188;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_188 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_186 += tmp_188;
   std::complex<double> tmp_187;
   std::complex<double> tmp_189;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_189 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_187 += tmp_189;
   result += (-0.05*Sqr(g1)) * tmp_186 * tmp_187;
   std::complex<double> tmp_190;
   std::complex<double> tmp_192;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_192 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_190 += tmp_192;
   std::complex<double> tmp_191;
   std::complex<double> tmp_193;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_193 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_191 += tmp_193;
   result += (0.2*Sqr(g1)) * tmp_190 * tmp_191;
   std::complex<double> tmp_194;
   std::complex<double> tmp_196;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_196 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_194 += tmp_196;
   std::complex<double> tmp_195;
   std::complex<double> tmp_197;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_197 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_195 += tmp_197;
   result += (-0.05*Sqr(g1)) * tmp_194 * tmp_195;
   std::complex<double> tmp_198;
   std::complex<double> tmp_200;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_200 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_198 += tmp_200;
   std::complex<double> tmp_199;
   std::complex<double> tmp_201;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_201 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_199 += tmp_201;
   result += (0.2*Sqr(g1)) * tmp_198 * tmp_199;
   std::complex<double> tmp_202;
   std::complex<double> tmp_204;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_205;
      std::complex<double> tmp_206;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_206 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_205 += tmp_206;
      tmp_204 += (Conj(ZU(gI2,j2))) * tmp_205;
   }
   tmp_202 += tmp_204;
   std::complex<double> tmp_203;
   std::complex<double> tmp_207;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_208;
      std::complex<double> tmp_209;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_209 += Conj(Yd(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_208 += tmp_209;
      tmp_207 += (ZU(gI1,j4)) * tmp_208;
   }
   tmp_203 += tmp_207;
   result += (-1) * tmp_202 * tmp_203;
   if (gO1 < 3) {
      std::complex<double> tmp_210;
      std::complex<double> tmp_211;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_211 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_210 += tmp_211;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_210;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_212;
      std::complex<double> tmp_213;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_213 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_212 += tmp_213;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_212;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_214;
      std::complex<double> tmp_215;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_215 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_214 += tmp_215;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_214;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_216;
      std::complex<double> tmp_217;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_217 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_216 += tmp_217;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_216;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_218;
      std::complex<double> tmp_219;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_219 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_218 += tmp_219;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_218;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_220;
      std::complex<double> tmp_221;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_221 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_220 += tmp_221;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_220;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_222;
      std::complex<double> tmp_224;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_224 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_222 += tmp_224;
      std::complex<double> tmp_223;
      std::complex<double> tmp_225;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_225 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_223 += tmp_225;
      result += (-1) * tmp_222 * tmp_223;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZU(gI2,gO2))*Sqr(g2)*ZU(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_226;
      std::complex<double> tmp_227;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_227 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_226 += tmp_227;
      result += (MuD*ZHR(gI2,0)) * tmp_226;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_228;
      std::complex<double> tmp_229;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_229 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_228 += tmp_229;
      result += (0.7071067811865475*LamSD*vS*ZHR(gI2,0)) * tmp_228;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_230;
      std::complex<double> tmp_231;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_231 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_230 += tmp_231;
      result += (0.5*LamTD*vT*ZHR(gI2,0)) * tmp_230;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuRpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_232;
   std::complex<double> tmp_233;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_234;
      std::complex<double> tmp_235;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_235 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_234 += tmp_235;
      tmp_233 += (Conj(ZU(gI1,j2))) * tmp_234;
   }
   tmp_232 += tmp_233;
   result += (0.7071067811865475*vS*Conj(LamSD)*Conj(ZRP(gI2,1))) * tmp_232;
   std::complex<double> tmp_236;
   std::complex<double> tmp_237;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_238;
      std::complex<double> tmp_239;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_239 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_238 += tmp_239;
      tmp_237 += (Conj(ZU(gI1,j2))) * tmp_238;
   }
   tmp_236 += tmp_237;
   result += (-0.5*vT*Conj(LamTD)*Conj(ZRP(gI2,1))) * tmp_236;
   std::complex<double> tmp_240;
   std::complex<double> tmp_241;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_242;
      std::complex<double> tmp_243;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_243 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_242 += tmp_243;
      tmp_241 += (Conj(ZU(gI1,j2))) * tmp_242;
   }
   tmp_240 += tmp_241;
   result += (Conj(MuD)*Conj(ZRP(gI2,1))) * tmp_240;
   if (gO2 < 3) {
      std::complex<double> tmp_244;
      std::complex<double> tmp_245;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_245 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_244 += tmp_245;
      result += (-(MuU*Conj(ZRP(gI2,0)))) * tmp_244;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_246;
      std::complex<double> tmp_247;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_247 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_246 += tmp_247;
      result += (-0.7071067811865475*LamSU*vS*Conj(ZRP(gI2,0))) * tmp_246;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_248;
      std::complex<double> tmp_249;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_249 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_248 += tmp_249;
      result += (-0.5*LamTU*vT*Conj(ZRP(gI2,0))) * tmp_248;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_250;
   std::complex<double> tmp_251;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_252;
      std::complex<double> tmp_253;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_253 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_252 += tmp_253;
      tmp_251 += (Conj(ZD(gI1,j2))) * tmp_252;
   }
   tmp_250 += tmp_251;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1)) *
      tmp_250;
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
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_257 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_256 += tmp_257;
   result += (std::complex<double>(0,-0.2581988897471611)*g1*MDBS*ZA(gI2,2)) *
      tmp_256;
   std::complex<double> tmp_258;
   std::complex<double> tmp_259;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_259 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_258 += tmp_259;
   result += (std::complex<double>(0,0.2581988897471611)*g1*Conj(MDBS)*ZA(gI2,2
      )) * tmp_258;
   if (gO2 < 3) {
      result += std::complex<double>(0,-0.12909944487358055)*g1*MDBS*Conj(ZD
         (gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,0.12909944487358055)*g1*Conj(MDBS)*
         Conj(ZD(gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,0.5)*g2*MDWBT*Conj(ZD(gI1,gO2))*ZA(
         gI2,3);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,-0.5)*g2*Conj(MDWBT)*Conj(ZD(gI1,gO2)
         )*ZA(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSdhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_260;
   std::complex<double> tmp_261;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_261 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_260 += tmp_261;
   result += (0.1*vd*Sqr(g1)*ZH(gI2,0)) * tmp_260;
   std::complex<double> tmp_262;
   std::complex<double> tmp_263;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_264;
      std::complex<double> tmp_265;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_266;
         std::complex<double> tmp_267;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_267 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_266 += tmp_267;
         tmp_265 += (KroneckerDelta(gO2,3 + j2)) * tmp_266;
      }
      tmp_264 += tmp_265;
      tmp_263 += (Conj(ZD(gI1,3 + j3))) * tmp_264;
   }
   tmp_262 += tmp_263;
   result += (-(vd*ZH(gI2,0))) * tmp_262;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_268;
      std::complex<double> tmp_269;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_270;
         std::complex<double> tmp_271;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_271 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_270 += tmp_271;
         tmp_269 += (Conj(ZD(gI1,j2))) * tmp_270;
      }
      tmp_268 += tmp_269;
      result += (-(vd*ZH(gI2,0))) * tmp_268;
   }
   std::complex<double> tmp_272;
   std::complex<double> tmp_273;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_273 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_272 += tmp_273;
   result += (-0.1*vu*Sqr(g1)*ZH(gI2,1)) * tmp_272;
   std::complex<double> tmp_274;
   std::complex<double> tmp_275;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_276;
      std::complex<double> tmp_277;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_277 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_276 += tmp_277;
      tmp_275 += (Conj(ZD(gI1,j2))) * tmp_276;
   }
   tmp_274 += tmp_275;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_274;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZD(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZD(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_278;
      std::complex<double> tmp_279;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_279 += Conj(Yd(j1,gO2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_278 += tmp_279;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_278;
   }
   std::complex<double> tmp_280;
   std::complex<double> tmp_281;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_281 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_280 += tmp_281;
   result += (-0.2581988897471611*g1*MDBS*ZH(gI2,2)) * tmp_280;
   std::complex<double> tmp_282;
   std::complex<double> tmp_283;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_283 += Conj(ZD(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_282 += tmp_283;
   result += (-0.2581988897471611*g1*Conj(MDBS)*ZH(gI2,2)) * tmp_282;
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*MDBS*Conj(ZD(gI1,gO2))*ZH(gI2,2);
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(MDBS)*Conj(ZD(gI1,gO2))*ZH(gI2,
         2);
   }
   if (gO2 < 3) {
      result += 0.5*g2*MDWBT*Conj(ZD(gI1,gO2))*ZH(gI2,3);
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(MDWBT)*Conj(ZD(gI1,gO2))*ZH(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdSuHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_284;
   std::complex<double> tmp_285;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_286;
      std::complex<double> tmp_287;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_288;
         std::complex<double> tmp_289;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_289 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_288 += tmp_289;
         tmp_287 += (KroneckerDelta(gO2,3 + j2)) * tmp_288;
      }
      tmp_286 += tmp_287;
      tmp_285 += (Conj(ZU(gI1,3 + j3))) * tmp_286;
   }
   tmp_284 += tmp_285;
   result += (0.7071067811865475*vu*ZP(gI2,0)) * tmp_284;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_290;
      std::complex<double> tmp_291;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_291 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_290 += tmp_291;
      result += (Mu*ZP(gI2,0)) * tmp_290;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_292;
      std::complex<double> tmp_293;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_294;
         std::complex<double> tmp_295;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_295 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_294 += tmp_295;
         tmp_293 += (Conj(ZU(gI1,j2))) * tmp_294;
      }
      tmp_292 += tmp_293;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_292;
   }
   std::complex<double> tmp_296;
   std::complex<double> tmp_297;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_298;
      std::complex<double> tmp_299;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_299 += KroneckerDelta(gO2,3 + j1)*Yd(j1,j2);
      }
      tmp_298 += tmp_299;
      tmp_297 += (Conj(ZU(gI1,j2))) * tmp_298;
   }
   tmp_296 += tmp_297;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_296;
   std::complex<double> tmp_300;
   std::complex<double> tmp_301;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_302;
      std::complex<double> tmp_303;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_304;
         std::complex<double> tmp_305;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_305 += Conj(Yu(j3,j1))*Yd(j2,j1);
         }
         tmp_304 += tmp_305;
         tmp_303 += (KroneckerDelta(gO2,3 + j2)) * tmp_304;
      }
      tmp_302 += tmp_303;
      tmp_301 += (Conj(ZU(gI1,3 + j3))) * tmp_302;
   }
   tmp_300 += tmp_301;
   result += (0.7071067811865475*vd*ZP(gI2,1)) * tmp_300;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_306;
      std::complex<double> tmp_307;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_308;
         std::complex<double> tmp_309;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_309 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_308 += tmp_309;
         tmp_307 += (Conj(ZU(gI1,j2))) * tmp_308;
      }
      tmp_306 += tmp_307;
      result += (0.7071067811865475*vu*ZP(gI2,1)) * tmp_306;
   }
   if (gO2 < 3) {
      result += -(g2*MDWBT*Conj(ZU(gI1,gO2))*ZP(gI2,2));
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,2);
   }
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZU(gI1,gO2))*Sqr(g2)*ZP(gI2,3);
   }
   if (gO2 < 3) {
      result += -(g2*Conj(MDWBT)*Conj(ZU(gI1,gO2))*ZP(gI2,3));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdbarGluFdPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_310;
   std::complex<double> tmp_311;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_311 += KroneckerDelta(gO2,3 + j1)*ZDR(gI2,j1);
   }
   tmp_310 += tmp_311;
   result += (1.4142135623730951*g3) * tmp_310;

   return result;
}

double CLASSNAME::CpconjUSdbarGluFdPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpconjUSdGluFdPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdGluFdPL(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(ZDL(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSdconjSOcSd(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_312;
   std::complex<double> tmp_313;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_313 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_312 += tmp_313;
   result += (1.4142135623730951*g3*Conj(MDGoc)) * tmp_312;
   if (gO2 < 3) {
      result += -1.4142135623730951*g3*Conj(MDGoc)*Conj(ZD(gI2,gO2));
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

   std::complex<double> tmp_314;
   std::complex<double> tmp_315;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_315 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_314 += tmp_315;
   result += (-0.2581988897471611*g1*Cos(ThetaW())) * tmp_314;
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

   std::complex<double> tmp_316;
   std::complex<double> tmp_317;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_317 += Conj(ZD(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_316 += tmp_317;
   result += (0.2581988897471611*g1*Sin(ThetaW())) * tmp_316;
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

std::complex<double> CLASSNAME::CpUSvconjUSvconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZHR(gI1,0)*
      ZHR(gI2,0) - ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) - 5*Sqr(g2))*(Conj(ZRP(gI2
      ,0))*ZRP(gI1,0) - Conj(ZRP(gI2,1))*ZRP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjRpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_318;
      std::complex<double> tmp_319;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_319 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_318 += tmp_319;
      result += (MuD*ZRP(gI1,1)) * tmp_318;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_320;
      std::complex<double> tmp_321;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_321 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_320 += tmp_321;
      result += (0.7071067811865475*LamSD*vS*ZRP(gI1,1)) * tmp_320;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_322;
      std::complex<double> tmp_323;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_323 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_322 += tmp_323;
      result += (-0.5*LamTD*vT*ZRP(gI1,1)) * tmp_322;
   }

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

std::complex<double> CLASSNAME::CpconjUSvFeCha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_324;
      std::complex<double> tmp_325;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_325 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_324 += tmp_325;
      result += (UM1(gI2,1)) * tmp_324;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvFeCha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(UP1(gI2,0))*Conj(ZEL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvSvAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += std::complex<double>(0,0.3872983346207417)*g1*MDBS*Conj(ZV(
         gI1,gO2))*ZA(gI2,2);
   }
   if (gI1 < 3) {
      result += std::complex<double>(0,-0.3872983346207417)*g1*Conj(MDBS)*
         Conj(ZV(gI1,gO2))*ZA(gI2,2);
   }
   if (gI1 < 3) {
      result += std::complex<double>(0,-0.5)*g2*MDWBT*Conj(ZV(gI1,gO2))*ZA(
         gI2,3);
   }
   if (gI1 < 3) {
      result += std::complex<double>(0,0.5)*g2*Conj(MDWBT)*Conj(ZV(gI1,gO2))
         *ZA(gI2,3);
   }

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
   if (gI1 < 3) {
      result += 0.3872983346207417*g1*MDBS*Conj(ZV(gI1,gO2))*ZH(gI2,2);
   }
   if (gI1 < 3) {
      result += 0.3872983346207417*g1*Conj(MDBS)*Conj(ZV(gI1,gO2))*ZH(gI2,2)
         ;
   }
   if (gI1 < 3) {
      result += -0.5*g2*MDWBT*Conj(ZV(gI1,gO2))*ZH(gI2,3);
   }
   if (gI1 < 3) {
      result += -0.5*g2*Conj(MDWBT)*Conj(ZV(gI1,gO2))*ZH(gI2,3);
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
      result += 0.5477225575051661*g1*Conj(ZN1(gI2,0))*KroneckerDelta(gI1,
         gO1);
   }
   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN1(gI2,1))*KroneckerDelta(gI1,
         gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZA(gI1,0)*ZA
      (gI2,0) - ZA(gI1,1)*ZA(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_326;
      std::complex<double> tmp_327;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_327 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_326 += tmp_327;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_326;
   }
   result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*KroneckerDelta(gO1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*(ZH(gI1,0)*ZH
      (gI2,0) - ZH(gI1,1)*ZH(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUSvconjHpmSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_328;
      std::complex<double> tmp_329;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_330;
         std::complex<double> tmp_331;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_331 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_330 += tmp_331;
         tmp_329 += (Conj(ZE(gI2,j2))) * tmp_330;
      }
      tmp_328 += tmp_329;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_328;
   }
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_332;
      std::complex<double> tmp_333;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_333 += Conj(Ye(j1,gO2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_332 += tmp_333;
      result += (Mu*ZP(gI1,1)) * tmp_332;
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,2);
   }
   if (gO2 < 3) {
      result += -(g2*Conj(MDWBT)*Conj(ZE(gI2,gO2))*ZP(gI1,2));
   }
   if (gO2 < 3) {
      result += -(g2*MDWBT*Conj(ZE(gI2,gO2))*ZP(gI1,3));
   }
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZE(gI2,gO2))*Sqr(g2)*ZP(gI1,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_334;
   std::complex<double> tmp_335;
   std::complex<double> tmp_336;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_336 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_335 += tmp_336;
   tmp_334 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_335;
   std::complex<double> tmp_337;
   std::complex<double> tmp_338;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_338 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_337 += tmp_338;
   tmp_334 += (std::complex<double>(0,0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_337;
   std::complex<double> tmp_339;
   std::complex<double> tmp_340;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_340 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_339 += tmp_340;
   tmp_334 += (std::complex<double>(0,0.1)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_339;
   result += (std::complex<double>(0,-1)) * tmp_334;

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_341;
   std::complex<double> tmp_342;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_342 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_341 += tmp_342;
   result += (-0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_341;
   std::complex<double> tmp_343;
   std::complex<double> tmp_344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_344 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_343 += tmp_344;
   result += (0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_343;
   std::complex<double> tmp_345;
   std::complex<double> tmp_346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_346 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_345 += tmp_346;
   result += (0.3*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_345;
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_347;
      std::complex<double> tmp_349;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_349 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_347 += tmp_349;
      std::complex<double> tmp_348;
      std::complex<double> tmp_350;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_350 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_348 += tmp_350;
      result += (-1) * tmp_347 * tmp_348;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZE(gI2,gO2))*Sqr(g2)*ZE(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSvconjUSvconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_351;
   std::complex<double> tmp_352;
   std::complex<double> tmp_353;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_353 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_352 += tmp_353;
   tmp_351 += (std::complex<double>(0,0.05)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_352;
   std::complex<double> tmp_354;
   std::complex<double> tmp_355;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_355 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_354 += tmp_355;
   tmp_351 += (std::complex<double>(0,-0.25)*KroneckerDelta(gO1,gO2)*Sqr(g2)) *
      tmp_354;
   std::complex<double> tmp_356;
   std::complex<double> tmp_357;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_357 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_356 += tmp_357;
   tmp_351 += (std::complex<double>(0,-0.2)*KroneckerDelta(gO1,gO2)*Sqr(g1)) *
      tmp_356;
   result += (std::complex<double>(0,-1)) * tmp_351;

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

   std::complex<double> tmp_358;
   std::complex<double> tmp_359;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_359 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_358 += tmp_359;
   result += (0.5333333333333333*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_358;
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

double CLASSNAME::CpUSuconjUSuconjSOcSOc(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

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

std::complex<double> CLASSNAME::CpUSuconjUSuconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_360;
   std::complex<double> tmp_361;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_361 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_360 += tmp_361;
   result += (0.2*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0)) * tmp_360;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   std::complex<double> tmp_362;
   std::complex<double> tmp_363;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_363 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_362 += tmp_363;
   result += (-0.2*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1)) * tmp_362;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_364;
   std::complex<double> tmp_365;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_365 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_364 += tmp_365;
   result += (-0.2*Conj(ZRP(gI2,0))*Sqr(g1)*ZRP(gI1,0)) * tmp_364;
   if (gO1 < 3) {
      result += 0.05*Conj(ZRP(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZRP(
         gI1,0);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZRP(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZRP(
         gI1,0);
   }
   std::complex<double> tmp_366;
   std::complex<double> tmp_367;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_367 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_366 += tmp_367;
   result += (0.2*Conj(ZRP(gI2,1))*Sqr(g1)*ZRP(gI1,1)) * tmp_366;
   if (gO1 < 3) {
      result += -0.05*Conj(ZRP(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZRP(
         gI1,1);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZRP(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZRP(
         gI1,1);
   }

   return result;
}

double CLASSNAME::CpconjUSubarCha2FdPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarCha2FdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_368;
   std::complex<double> tmp_369;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_370;
      std::complex<double> tmp_371;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_371 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_370 += tmp_371;
      tmp_369 += (Conj(ZDL(gI2,j2))) * tmp_370;
   }
   tmp_368 += tmp_369;
   result += (Conj(UP2(gI1,1))) * tmp_368;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjRpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_372;
   std::complex<double> tmp_373;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_374;
      std::complex<double> tmp_375;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_375 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_374 += tmp_375;
      tmp_373 += (Conj(ZD(gI2,j2))) * tmp_374;
   }
   tmp_372 += tmp_373;
   result += (-0.7071067811865475*vS*Conj(LamSU)*ZRP(gI1,0)) * tmp_372;
   std::complex<double> tmp_376;
   std::complex<double> tmp_377;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_378;
      std::complex<double> tmp_379;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_379 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_378 += tmp_379;
      tmp_377 += (Conj(ZD(gI2,j2))) * tmp_378;
   }
   tmp_376 += tmp_377;
   result += (-0.5*vT*Conj(LamTU)*ZRP(gI1,0)) * tmp_376;
   std::complex<double> tmp_380;
   std::complex<double> tmp_381;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_382;
      std::complex<double> tmp_383;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_383 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_382 += tmp_383;
      tmp_381 += (Conj(ZD(gI2,j2))) * tmp_382;
   }
   tmp_380 += tmp_381;
   result += (-(Conj(MuU)*ZRP(gI1,0))) * tmp_380;
   if (gO2 < 3) {
      std::complex<double> tmp_384;
      std::complex<double> tmp_385;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_385 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_384 += tmp_385;
      result += (MuD*ZRP(gI1,1)) * tmp_384;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_386;
      std::complex<double> tmp_387;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_387 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_386 += tmp_387;
      result += (0.7071067811865475*LamSD*vS*ZRP(gI1,1)) * tmp_386;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_388;
      std::complex<double> tmp_389;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_389 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_388 += tmp_389;
      result += (-0.5*LamTD*vT*ZRP(gI1,1)) * tmp_388;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_390;
   std::complex<double> tmp_391;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_391 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_390 += tmp_391;
   result += (-0.2*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_390;
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

std::complex<double> CLASSNAME::CpconjUSuFdCha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_392;
      std::complex<double> tmp_393;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_393 += Conj(Yd(j1,gO2))*ZDR(gI1,j1);
      }
      tmp_392 += tmp_393;
      result += (UM1(gI2,1)) * tmp_392;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFdCha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(UP1(gI2,0))*Conj(ZDL(gI1,gO1)));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_394;
      std::complex<double> tmp_395;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_395 += Conj(Yu(j1,gO2))*ZUR(gI1,j1);
      }
      tmp_394 += tmp_395;
      result += (-ZN2(gI2,3)) * tmp_394;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuFuChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZN1(gI2,0))*Conj(ZUL(gI1,gO1));
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZN1(gI2,1))*Conj(ZUL(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_396;
   std::complex<double> tmp_397;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_397 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_396 += tmp_397;
   result += (-0.2*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_396;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
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

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_408;
   std::complex<double> tmp_409;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_409 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_408 += tmp_409;
   result += (-0.2*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_408;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_410;
      std::complex<double> tmp_411;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_411 += Conj(Yd(j1,gO2))*Yd(j1,gO1);
      }
      tmp_410 += tmp_411;
      result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_410;
   }
   std::complex<double> tmp_412;
   std::complex<double> tmp_413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_413 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_412 += tmp_413;
   result += (0.2*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_412;
   std::complex<double> tmp_414;
   std::complex<double> tmp_415;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_416;
      std::complex<double> tmp_417;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_418;
         std::complex<double> tmp_419;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_419 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_418 += tmp_419;
         tmp_417 += (KroneckerDelta(gO2,3 + j2)) * tmp_418;
      }
      tmp_416 += tmp_417;
      tmp_415 += (KroneckerDelta(gO1,3 + j3)) * tmp_416;
   }
   tmp_414 += tmp_415;
   result += (-(ZP(gI1,1)*ZP(gI2,1))) * tmp_414;
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
      result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_420;
   std::complex<double> tmp_421;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_421 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_420 += tmp_421;
   result += (-0.2*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_420;
   if (gO1 < 3) {
      result += 0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   std::complex<double> tmp_422;
   std::complex<double> tmp_423;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_423 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_422 += tmp_423;
   result += (0.2*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_422;
   std::complex<double> tmp_424;
   std::complex<double> tmp_425;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_426;
      std::complex<double> tmp_427;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_428;
         std::complex<double> tmp_429;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_429 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_428 += tmp_429;
         tmp_427 += (KroneckerDelta(gO2,3 + j2)) * tmp_428;
      }
      tmp_426 += tmp_427;
      tmp_425 += (KroneckerDelta(gO1,3 + j3)) * tmp_426;
   }
   tmp_424 += tmp_425;
   result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_424;
   if (gO1 < 3) {
      result += -0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_430;
      std::complex<double> tmp_431;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_431 += Conj(Yu(j1,gO2))*Yu(j1,gO1);
      }
      tmp_430 += tmp_431;
      result += (-(ZH(gI1,1)*ZH(gI2,1))) * tmp_430;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChiFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_432;
   std::complex<double> tmp_433;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_433 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_432 += tmp_433;
   result += (0.7302967433402214*g1*ZN1(gI1,0)) * tmp_432;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarChiFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_434;
   std::complex<double> tmp_435;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_436;
      std::complex<double> tmp_437;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_437 += KroneckerDelta(gO1,3 + j1)*Yu(j1,j2);
      }
      tmp_436 += tmp_437;
      tmp_435 += (Conj(ZUL(gI2,j2))) * tmp_436;
   }
   tmp_434 += tmp_435;
   result += (-Conj(ZN2(gI1,3))) * tmp_434;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuconjHpmSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_438;
   std::complex<double> tmp_439;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_440;
      std::complex<double> tmp_441;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_441 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_440 += tmp_441;
      tmp_439 += (Conj(ZD(gI2,j2))) * tmp_440;
   }
   tmp_438 += tmp_439;
   result += (Conj(Mu)*ZP(gI1,0)) * tmp_438;
   std::complex<double> tmp_442;
   std::complex<double> tmp_443;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_444;
      std::complex<double> tmp_445;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_446;
         std::complex<double> tmp_447;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_447 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_446 += tmp_447;
         tmp_445 += (KroneckerDelta(gO2,3 + j2)) * tmp_446;
      }
      tmp_444 += tmp_445;
      tmp_443 += (Conj(ZD(gI2,3 + j3))) * tmp_444;
   }
   tmp_442 += tmp_443;
   result += (0.7071067811865475*vu*ZP(gI1,0)) * tmp_442;
   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_448;
      std::complex<double> tmp_449;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_450;
         std::complex<double> tmp_451;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_451 += Conj(Yd(j1,gO2))*Yd(j1,j2);
         }
         tmp_450 += tmp_451;
         tmp_449 += (Conj(ZD(gI2,j2))) * tmp_450;
      }
      tmp_448 += tmp_449;
      result += (0.7071067811865475*vd*ZP(gI1,0)) * tmp_448;
   }
   std::complex<double> tmp_452;
   std::complex<double> tmp_453;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_454;
      std::complex<double> tmp_455;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_456;
         std::complex<double> tmp_457;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_457 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_456 += tmp_457;
         tmp_455 += (KroneckerDelta(gO2,3 + j2)) * tmp_456;
      }
      tmp_454 += tmp_455;
      tmp_453 += (Conj(ZD(gI2,3 + j3))) * tmp_454;
   }
   tmp_452 += tmp_453;
   result += (0.7071067811865475*vd*ZP(gI1,1)) * tmp_452;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_458;
      std::complex<double> tmp_459;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_459 += Conj(Yd(j1,gO2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_458 += tmp_459;
      result += (Mu*ZP(gI1,1)) * tmp_458;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_460;
      std::complex<double> tmp_461;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_462;
         std::complex<double> tmp_463;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_463 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_462 += tmp_463;
         tmp_461 += (Conj(ZD(gI2,j2))) * tmp_462;
      }
      tmp_460 += tmp_461;
      result += (0.7071067811865475*vu*ZP(gI1,1)) * tmp_460;
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,2);
   }
   if (gO2 < 3) {
      result += -(g2*Conj(MDWBT)*Conj(ZD(gI2,gO2))*ZP(gI1,2));
   }
   if (gO2 < 3) {
      result += -(g2*MDWBT*Conj(ZD(gI2,gO2))*ZP(gI1,3));
   }
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZD(gI2,gO2))*Sqr(g2)*ZP(gI1,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuSOc(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   std::complex<double> tmp_464;
   std::complex<double> tmp_465;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_465 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_464 += tmp_465;
   result += (1.4142135623730951*g3*MDGoc) * tmp_464;
   if (gO2 < 3) {
      result += -1.4142135623730951*g3*MDGoc*Conj(ZU(gI1,gO2));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_466;
   std::complex<double> tmp_468;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_468 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_466 += tmp_468;
   std::complex<double> tmp_467;
   std::complex<double> tmp_469;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_469 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_467 += tmp_469;
   result += (0.1*Sqr(g1)) * tmp_466 * tmp_467;
   std::complex<double> tmp_470;
   std::complex<double> tmp_472;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_472 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_470 += tmp_472;
   std::complex<double> tmp_471;
   std::complex<double> tmp_473;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_473 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_471 += tmp_473;
   result += (0.2*Sqr(g1)) * tmp_470 * tmp_471;
   std::complex<double> tmp_474;
   std::complex<double> tmp_476;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_476 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_474 += tmp_476;
   std::complex<double> tmp_475;
   std::complex<double> tmp_477;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_477 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_475 += tmp_477;
   result += (0.1*Sqr(g1)) * tmp_474 * tmp_475;
   std::complex<double> tmp_478;
   std::complex<double> tmp_480;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_480 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_478 += tmp_480;
   std::complex<double> tmp_479;
   std::complex<double> tmp_481;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_481 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_479 += tmp_481;
   result += (0.2*Sqr(g1)) * tmp_478 * tmp_479;
   std::complex<double> tmp_482;
   std::complex<double> tmp_484;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_485;
      std::complex<double> tmp_486;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_486 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_485 += tmp_486;
      tmp_484 += (Conj(ZD(gI2,j2))) * tmp_485;
   }
   tmp_482 += tmp_484;
   std::complex<double> tmp_483;
   std::complex<double> tmp_487;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_488;
      std::complex<double> tmp_489;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_489 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_488 += tmp_489;
      tmp_487 += (ZD(gI1,j4)) * tmp_488;
   }
   tmp_483 += tmp_487;
   result += (-1) * tmp_482 * tmp_483;
   if (gO1 < 3) {
      std::complex<double> tmp_490;
      std::complex<double> tmp_491;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_491 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_490 += tmp_491;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_490;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_492;
      std::complex<double> tmp_493;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_493 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_492 += tmp_493;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_492;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_494;
      std::complex<double> tmp_495;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_495 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_494 += tmp_495;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_494;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_496;
      std::complex<double> tmp_497;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_497 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_496 += tmp_497;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_496;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_498;
      std::complex<double> tmp_499;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_499 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_498 += tmp_499;
      result += (0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_498;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_500;
      std::complex<double> tmp_501;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_501 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_500 += tmp_501;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_500;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_502;
      std::complex<double> tmp_504;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_504 += Yd(j1,gO1)*ZD(gI1,3 + j1);
      }
      tmp_502 += tmp_504;
      std::complex<double> tmp_503;
      std::complex<double> tmp_505;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_505 += Conj(Yd(j3,gO2))*Conj(ZD(gI2,3 + j3));
      }
      tmp_503 += tmp_505;
      result += (-1) * tmp_502 * tmp_503;
   }
   if (gO1 < 3 && gO2 < 3) {
      result += -0.5*Conj(ZD(gI2,gO2))*Sqr(g2)*ZD(gI1,gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_506;
   std::complex<double> tmp_508;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_508 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_506 += tmp_508;
   std::complex<double> tmp_507;
   std::complex<double> tmp_509;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_509 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_507 += tmp_509;
   result += (-0.1*Sqr(g1)) * tmp_506 * tmp_507;
   std::complex<double> tmp_510;
   std::complex<double> tmp_512;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_512 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_510 += tmp_512;
   std::complex<double> tmp_511;
   std::complex<double> tmp_513;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_513 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_511 += tmp_513;
   result += (0.2*Sqr(g1)) * tmp_510 * tmp_511;
   std::complex<double> tmp_514;
   std::complex<double> tmp_516;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_516 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_514 += tmp_516;
   std::complex<double> tmp_515;
   std::complex<double> tmp_517;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_517 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_515 += tmp_517;
   result += (-0.1*Sqr(g1)) * tmp_514 * tmp_515;
   std::complex<double> tmp_518;
   std::complex<double> tmp_520;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_520 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_518 += tmp_520;
   std::complex<double> tmp_519;
   std::complex<double> tmp_521;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_521 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_519 += tmp_521;
   result += (0.2*Sqr(g1)) * tmp_518 * tmp_519;
   if (gO1 < 3) {
      std::complex<double> tmp_522;
      std::complex<double> tmp_523;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_523 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_522 += tmp_523;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_522;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_524;
      std::complex<double> tmp_525;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_525 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_524 += tmp_525;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_524;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_526;
      std::complex<double> tmp_527;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_527 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_526 += tmp_527;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_526;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_528;
      std::complex<double> tmp_529;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_529 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_528 += tmp_529;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_528;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_530;
      std::complex<double> tmp_531;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_531 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_530 += tmp_531;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_530;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_532;
      std::complex<double> tmp_533;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_533 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_532 += tmp_533;
      result += (-0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_532;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSuconjUSuconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_534;
   std::complex<double> tmp_536;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_536 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_534 += tmp_536;
   std::complex<double> tmp_535;
   std::complex<double> tmp_537;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_537 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_535 += tmp_537;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_534 * tmp_535;
   std::complex<double> tmp_538;
   std::complex<double> tmp_540;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_540 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
   }
   tmp_538 += tmp_540;
   std::complex<double> tmp_539;
   std::complex<double> tmp_541;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_541 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_539 += tmp_541;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_538 * tmp_539;
   std::complex<double> tmp_542;
   std::complex<double> tmp_544;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_544 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_542 += tmp_544;
   std::complex<double> tmp_543;
   std::complex<double> tmp_545;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_545 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_543 += tmp_545;
   result += (0.1*Sqr(g1)) * tmp_542 * tmp_543;
   std::complex<double> tmp_546;
   std::complex<double> tmp_548;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_548 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_546 += tmp_548;
   std::complex<double> tmp_547;
   std::complex<double> tmp_549;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_549 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_547 += tmp_549;
   result += (-0.4*Sqr(g1)) * tmp_546 * tmp_547;
   std::complex<double> tmp_550;
   std::complex<double> tmp_552;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_552 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_550 += tmp_552;
   std::complex<double> tmp_551;
   std::complex<double> tmp_553;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_553 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_551 += tmp_553;
   result += (0.1*Sqr(g1)) * tmp_550 * tmp_551;
   std::complex<double> tmp_554;
   std::complex<double> tmp_556;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_556 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_554 += tmp_556;
   std::complex<double> tmp_555;
   std::complex<double> tmp_557;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_557 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_555 += tmp_557;
   result += (-0.4*Sqr(g1)) * tmp_554 * tmp_555;
   std::complex<double> tmp_558;
   std::complex<double> tmp_560;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_560 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_558 += tmp_560;
   std::complex<double> tmp_559;
   std::complex<double> tmp_561;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_561 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_559 += tmp_561;
   result += (-0.13333333333333333*Sqr(g1)) * tmp_558 * tmp_559;
   std::complex<double> tmp_562;
   std::complex<double> tmp_564;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_564 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_562 += tmp_564;
   std::complex<double> tmp_563;
   std::complex<double> tmp_565;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_565 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
   }
   tmp_563 += tmp_565;
   result += (-0.6666666666666666*Sqr(g3)) * tmp_562 * tmp_563;
   std::complex<double> tmp_566;
   std::complex<double> tmp_568;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_569;
      std::complex<double> tmp_570;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_570 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_569 += tmp_570;
      tmp_568 += (Conj(ZU(gI2,j2))) * tmp_569;
   }
   tmp_566 += tmp_568;
   std::complex<double> tmp_567;
   std::complex<double> tmp_571;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_572;
      std::complex<double> tmp_573;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_573 += Conj(Yu(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_572 += tmp_573;
      tmp_571 += (ZU(gI1,j4)) * tmp_572;
   }
   tmp_567 += tmp_571;
   result += (-1) * tmp_566 * tmp_567;
   if (gO1 < 3) {
      std::complex<double> tmp_574;
      std::complex<double> tmp_575;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_575 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_574 += tmp_575;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_574;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_576;
      std::complex<double> tmp_577;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_577 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_576 += tmp_577;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_576;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_578;
      std::complex<double> tmp_579;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_579 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_578 += tmp_579;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_578;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_580;
      std::complex<double> tmp_581;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_581 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_580 += tmp_581;
      result += (-0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_580;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_582;
      std::complex<double> tmp_583;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_583 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_582 += tmp_583;
      result += (-0.375*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_582;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_584;
      std::complex<double> tmp_585;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_585 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_584 += tmp_585;
      result += (0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_584;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_586;
      std::complex<double> tmp_588;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_588 += KroneckerDelta(gO2,3 + j1)*Yu(j1,gO1);
      }
      tmp_586 += tmp_588;
      std::complex<double> tmp_587;
      std::complex<double> tmp_589;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_590;
         std::complex<double> tmp_591;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_591 += Conj(Yu(j3,j4))*Conj(ZU(gI2,3 + j3));
         }
         tmp_590 += tmp_591;
         tmp_589 += (ZU(gI1,j4)) * tmp_590;
      }
      tmp_587 += tmp_589;
      result += (-3) * tmp_586 * tmp_587;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_592;
      std::complex<double> tmp_593;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_593 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_592 += tmp_593;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_592;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_594;
      std::complex<double> tmp_595;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_595 += KroneckerDelta(gO1,3 + j1)*ZU(gI1,3 + j1);
      }
      tmp_594 += tmp_595;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_594;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_596;
      std::complex<double> tmp_597;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_597 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_596 += tmp_597;
      result += (0.03333333333333333*Conj(ZU(gI2,gO2))*Sqr(g1)) * tmp_596;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_598;
      std::complex<double> tmp_599;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_599 += KroneckerDelta(gO1,3 + j2)*ZU(gI1,3 + j2);
      }
      tmp_598 += tmp_599;
      result += (0.6666666666666666*Conj(ZU(gI2,gO2))*Sqr(g3)) * tmp_598;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_600;
      std::complex<double> tmp_602;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_603;
         std::complex<double> tmp_604;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_604 += Yu(j1,j2)*ZU(gI1,3 + j1);
         }
         tmp_603 += tmp_604;
         tmp_602 += (Conj(ZU(gI2,j2))) * tmp_603;
      }
      tmp_600 += tmp_602;
      std::complex<double> tmp_601;
      std::complex<double> tmp_605;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_605 += Conj(Yu(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_601 += tmp_605;
      result += (-3) * tmp_600 * tmp_601;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_606;
      std::complex<double> tmp_608;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_608 += Yu(j1,gO1)*ZU(gI1,3 + j1);
      }
      tmp_606 += tmp_608;
      std::complex<double> tmp_607;
      std::complex<double> tmp_609;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_609 += Conj(Yu(j3,gO2))*Conj(ZU(gI2,3 + j3));
      }
      tmp_607 += tmp_609;
      result += (-1) * tmp_606 * tmp_607;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_610;
      std::complex<double> tmp_611;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_611 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_610 += tmp_611;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_610;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_612;
      std::complex<double> tmp_613;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_613 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_612 += tmp_613;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_612;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_614;
      std::complex<double> tmp_615;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_615 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_614 += tmp_615;
      result += (0.03333333333333333*Sqr(g1)*ZU(gI1,gO1)) * tmp_614;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_616;
      std::complex<double> tmp_617;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_617 += Conj(ZU(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_616 += tmp_617;
      result += (0.6666666666666666*Sqr(g3)*ZU(gI1,gO1)) * tmp_616;
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

std::complex<double> CLASSNAME::CpconjUSuSuRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_618;
      std::complex<double> tmp_619;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_619 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_618 += tmp_619;
      result += (-(MuU*ZHR(gI2,1))) * tmp_618;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_620;
      std::complex<double> tmp_621;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_621 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_620 += tmp_621;
      result += (-0.7071067811865475*LamSU*vS*ZHR(gI2,1)) * tmp_620;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_622;
      std::complex<double> tmp_623;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_623 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_622 += tmp_623;
      result += (0.5*LamTU*vT*ZHR(gI2,1)) * tmp_622;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_624;
   std::complex<double> tmp_625;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_626;
      std::complex<double> tmp_627;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_627 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_626 += tmp_627;
      tmp_625 += (Conj(ZU(gI1,j2))) * tmp_626;
   }
   tmp_624 += tmp_625;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(Mu)*ZA(gI2,0)) *
      tmp_624;
   if (gO2 < 3) {
      std::complex<double> tmp_628;
      std::complex<double> tmp_629;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_629 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_628 += tmp_629;
      result += (std::complex<double>(0,0.7071067811865475)*Mu*ZA(gI2,0)) *
         tmp_628;
   }
   std::complex<double> tmp_630;
   std::complex<double> tmp_631;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_631 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_630 += tmp_631;
   result += (std::complex<double>(0,0.5163977794943222)*g1*MDBS*ZA(gI2,2)) *
      tmp_630;
   std::complex<double> tmp_632;
   std::complex<double> tmp_633;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_633 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_632 += tmp_633;
   result += (std::complex<double>(0,-0.5163977794943222)*g1*Conj(MDBS)*ZA(gI2,
      2)) * tmp_632;
   if (gO2 < 3) {
      result += std::complex<double>(0,-0.12909944487358055)*g1*MDBS*Conj(ZU
         (gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,0.12909944487358055)*g1*Conj(MDBS)*
         Conj(ZU(gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,-0.5)*g2*MDWBT*Conj(ZU(gI1,gO2))*ZA(
         gI2,3);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,0.5)*g2*Conj(MDWBT)*Conj(ZU(gI1,gO2))
         *ZA(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuSuhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_634;
   std::complex<double> tmp_635;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_635 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_634 += tmp_635;
   result += (-0.2*vd*Sqr(g1)*ZH(gI2,0)) * tmp_634;
   std::complex<double> tmp_636;
   std::complex<double> tmp_637;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_638;
      std::complex<double> tmp_639;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_639 += KroneckerDelta(gO2,3 + j1)*Yu(j1,j2);
      }
      tmp_638 += tmp_639;
      tmp_637 += (Conj(ZU(gI1,j2))) * tmp_638;
   }
   tmp_636 += tmp_637;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,0)) * tmp_636;
   if (gO2 < 3) {
      result += 0.05*vd*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += -0.25*vd*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_640;
      std::complex<double> tmp_641;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_641 += Conj(Yu(j1,gO2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_640 += tmp_641;
      result += (0.7071067811865475*Mu*ZH(gI2,0)) * tmp_640;
   }
   std::complex<double> tmp_642;
   std::complex<double> tmp_643;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_643 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_642 += tmp_643;
   result += (0.2*vu*Sqr(g1)*ZH(gI2,1)) * tmp_642;
   std::complex<double> tmp_644;
   std::complex<double> tmp_645;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_646;
      std::complex<double> tmp_647;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_648;
         std::complex<double> tmp_649;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_649 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_648 += tmp_649;
         tmp_647 += (KroneckerDelta(gO2,3 + j2)) * tmp_648;
      }
      tmp_646 += tmp_647;
      tmp_645 += (Conj(ZU(gI1,3 + j3))) * tmp_646;
   }
   tmp_644 += tmp_645;
   result += (-(vu*ZH(gI2,1))) * tmp_644;
   if (gO2 < 3) {
      result += -0.05*vu*Conj(ZU(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += 0.25*vu*Conj(ZU(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_650;
      std::complex<double> tmp_651;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_652;
         std::complex<double> tmp_653;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_653 += Conj(Yu(j1,gO2))*Yu(j1,j2);
         }
         tmp_652 += tmp_653;
         tmp_651 += (Conj(ZU(gI1,j2))) * tmp_652;
      }
      tmp_650 += tmp_651;
      result += (-(vu*ZH(gI2,1))) * tmp_650;
   }
   std::complex<double> tmp_654;
   std::complex<double> tmp_655;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_655 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_654 += tmp_655;
   result += (0.5163977794943222*g1*MDBS*ZH(gI2,2)) * tmp_654;
   std::complex<double> tmp_656;
   std::complex<double> tmp_657;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_657 += Conj(ZU(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_656 += tmp_657;
   result += (0.5163977794943222*g1*Conj(MDBS)*ZH(gI2,2)) * tmp_656;
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*MDBS*Conj(ZU(gI1,gO2))*ZH(gI2,2);
   }
   if (gO2 < 3) {
      result += -0.12909944487358055*g1*Conj(MDBS)*Conj(ZU(gI1,gO2))*ZH(gI2,
         2);
   }
   if (gO2 < 3) {
      result += -0.5*g2*MDWBT*Conj(ZU(gI1,gO2))*ZH(gI2,3);
   }
   if (gO2 < 3) {
      result += -0.5*g2*Conj(MDWBT)*Conj(ZU(gI1,gO2))*ZH(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSubarGluFuPR(unsigned gO2, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_658;
   std::complex<double> tmp_659;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_659 += KroneckerDelta(gO2,3 + j1)*ZUR(gI2,j1);
   }
   tmp_658 += tmp_659;
   result += (1.4142135623730951*g3) * tmp_658;

   return result;
}

double CLASSNAME::CpconjUSubarGluFuPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpconjUSuGluFuPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSuGluFuPL(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(ZUL(gI2,gO1));
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

   std::complex<double> tmp_660;
   std::complex<double> tmp_661;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_661 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_660 += tmp_661;
   result += (0.5163977794943222*g1*Cos(ThetaW())) * tmp_660;
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

   std::complex<double> tmp_662;
   std::complex<double> tmp_663;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_663 += Conj(ZU(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_662 += tmp_663;
   result += (-0.5163977794943222*g1*Sin(ThetaW())) * tmp_662;
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

   std::complex<double> tmp_664;
   std::complex<double> tmp_665;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_665 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_664 += tmp_665;
   result += (1.2*Sqr(g1)*Sqr(Sin(ThetaW()))) * tmp_664;
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

std::complex<double> CLASSNAME::CpUSeconjUSeconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_666;
   std::complex<double> tmp_667;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_667 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_666 += tmp_667;
   result += (-0.3*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0)) * tmp_666;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0);
   }
   std::complex<double> tmp_668;
   std::complex<double> tmp_669;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_669 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_668 += tmp_669;
   result += (0.3*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1)) * tmp_668;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZHR(gI1,1)*ZHR(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_670;
   std::complex<double> tmp_671;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_671 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_670 += tmp_671;
   result += (0.3*Conj(ZRP(gI2,0))*Sqr(g1)*ZRP(gI1,0)) * tmp_670;
   if (gO1 < 3) {
      result += -0.15*Conj(ZRP(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZRP(
         gI1,0);
   }
   if (gO1 < 3) {
      result += -0.25*Conj(ZRP(gI2,0))*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZRP(
         gI1,0);
   }
   std::complex<double> tmp_672;
   std::complex<double> tmp_673;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_673 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_672 += tmp_673;
   result += (-0.3*Conj(ZRP(gI2,1))*Sqr(g1)*ZRP(gI1,1)) * tmp_672;
   if (gO1 < 3) {
      result += 0.15*Conj(ZRP(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZRP(
         gI1,1);
   }
   if (gO1 < 3) {
      result += 0.25*Conj(ZRP(gI2,1))*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZRP(
         gI1,1);
   }

   return result;
}

double CLASSNAME::CpconjUSebarCha1FvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSebarCha1FvPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_674;
   std::complex<double> tmp_675;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_675 += KroneckerDelta(gO1,3 + j1)*Ye(j1,gI2);
   }
   tmp_674 += tmp_675;
   result += (Conj(UM1(gI1,1))) * tmp_674;

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_676;
   std::complex<double> tmp_677;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_677 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_676 += tmp_677;
   result += (0.3*KroneckerDelta(gI1,gI2)*Sqr(g1)) * tmp_676;
   std::complex<double> tmp_678;
   std::complex<double> tmp_680;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_681;
      std::complex<double> tmp_682;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_682 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_681 += tmp_682;
      tmp_680 += (Conj(ZV(gI2,j2))) * tmp_681;
   }
   tmp_678 += tmp_680;
   std::complex<double> tmp_679;
   std::complex<double> tmp_683;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_684;
      std::complex<double> tmp_685;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_685 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_684 += tmp_685;
      tmp_683 += (ZV(gI1,j4)) * tmp_684;
   }
   tmp_679 += tmp_683;
   result += (-1) * tmp_678 * tmp_679;
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

std::complex<double> CLASSNAME::CpconjUSeSvRpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_686;
   std::complex<double> tmp_687;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_688;
      std::complex<double> tmp_689;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_689 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_688 += tmp_689;
      tmp_687 += (Conj(ZV(gI1,j2))) * tmp_688;
   }
   tmp_686 += tmp_687;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) - vT*Conj(LamTD) + 2*Conj(
      MuD))*Conj(ZRP(gI2,1))) * tmp_686;

   return result;
}

double CLASSNAME::CpconjUSeFvCha2PR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFvCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += -(g2*Conj(UM2(gI2,0))*KroneckerDelta(gI1,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSvHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.35355339059327373*vd*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_690;
      std::complex<double> tmp_691;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_692;
         std::complex<double> tmp_693;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_693 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_692 += tmp_693;
         tmp_691 += (Conj(ZV(gI1,j2))) * tmp_692;
      }
      tmp_690 += tmp_691;
      result += (0.7071067811865475*vd*ZP(gI2,0)) * tmp_690;
   }
   std::complex<double> tmp_694;
   std::complex<double> tmp_695;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_696;
      std::complex<double> tmp_697;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_697 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_696 += tmp_697;
      tmp_695 += (Conj(ZV(gI1,j2))) * tmp_696;
   }
   tmp_694 += tmp_695;
   result += (Conj(Mu)*ZP(gI2,1)) * tmp_694;
   if (gO2 < 3) {
      result += -0.35355339059327373*vu*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,1);
   }
   if (gO2 < 3) {
      result += -(g2*MDWBT*Conj(ZV(gI1,gO2))*ZP(gI2,2));
   }
   if (gO2 < 3) {
      result += -0.5*vT*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,2);
   }
   if (gO2 < 3) {
      result += 0.5*vT*Conj(ZV(gI1,gO2))*Sqr(g2)*ZP(gI2,3);
   }
   if (gO2 < 3) {
      result += -(g2*Conj(MDWBT)*Conj(ZV(gI1,gO2))*ZP(gI2,3));
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_698;
      std::complex<double> tmp_699;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_699 += Conj(Ye(j1,gO2))*ZER(gI1,j1);
      }
      tmp_698 += tmp_699;
      result += (-ZN2(gI2,2)) * tmp_698;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeFeChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZEL(gI1,gO1))*Conj(ZN1(gI2,0));
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZEL(gI1,gO1))*Conj(ZN1(gI2,1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_700;
   std::complex<double> tmp_701;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_701 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_700 += tmp_701;
   result += (0.3*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0)) * tmp_700;
   std::complex<double> tmp_702;
   std::complex<double> tmp_703;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_704;
      std::complex<double> tmp_705;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_706;
         std::complex<double> tmp_707;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_707 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_706 += tmp_707;
         tmp_705 += (KroneckerDelta(gO2,3 + j2)) * tmp_706;
      }
      tmp_704 += tmp_705;
      tmp_703 += (KroneckerDelta(gO1,3 + j3)) * tmp_704;
   }
   tmp_702 += tmp_703;
   result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_702;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_708;
      std::complex<double> tmp_709;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_709 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_708 += tmp_709;
      result += (-(ZA(gI1,0)*ZA(gI2,0))) * tmp_708;
   }
   std::complex<double> tmp_710;
   std::complex<double> tmp_711;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_711 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_710 += tmp_711;
   result += (-0.3*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1)) * tmp_710;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZA(gI1,1)*ZA(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_712;
   std::complex<double> tmp_713;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_713 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_712 += tmp_713;
   result += (0.3*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0)) * tmp_712;
   std::complex<double> tmp_714;
   std::complex<double> tmp_715;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_716;
      std::complex<double> tmp_717;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_718;
         std::complex<double> tmp_719;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_719 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_718 += tmp_719;
         tmp_717 += (KroneckerDelta(gO2,3 + j2)) * tmp_718;
      }
      tmp_716 += tmp_717;
      tmp_715 += (KroneckerDelta(gO1,3 + j3)) * tmp_716;
   }
   tmp_714 += tmp_715;
   result += (-(ZP(gI1,0)*ZP(gI2,0))) * tmp_714;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,0)*ZP(gI2,0);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,0);
   }
   std::complex<double> tmp_720;
   std::complex<double> tmp_721;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_721 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_720 += tmp_721;
   result += (-0.3*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1)) * tmp_720;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2);
   }
   if (gO1 < 3) {
      result += 0.5*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSehhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_722;
   std::complex<double> tmp_723;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_723 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_722 += tmp_723;
   result += (0.3*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0)) * tmp_722;
   std::complex<double> tmp_724;
   std::complex<double> tmp_725;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_726;
      std::complex<double> tmp_727;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_728;
         std::complex<double> tmp_729;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_729 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_728 += tmp_729;
         tmp_727 += (KroneckerDelta(gO2,3 + j2)) * tmp_728;
      }
      tmp_726 += tmp_727;
      tmp_725 += (KroneckerDelta(gO1,3 + j3)) * tmp_726;
   }
   tmp_724 += tmp_725;
   result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_724;
   if (gO1 < 3) {
      result += -0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3) {
      result += 0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,0);
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_730;
      std::complex<double> tmp_731;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_731 += Conj(Ye(j1,gO2))*Ye(j1,gO1);
      }
      tmp_730 += tmp_731;
      result += (-(ZH(gI1,0)*ZH(gI2,0))) * tmp_730;
   }
   std::complex<double> tmp_732;
   std::complex<double> tmp_733;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_733 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_732 += tmp_733;
   result += (-0.3*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1)) * tmp_732;
   if (gO1 < 3) {
      result += 0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)*ZH(gI1,1)*ZH(gI2,1);
   }
   if (gO1 < 3) {
      result += -0.25*KroneckerDelta(gO1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSebarChiFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_734;
   std::complex<double> tmp_735;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_735 += KroneckerDelta(gO2,3 + j1)*ZER(gI2,j1);
   }
   tmp_734 += tmp_735;
   result += (-1.0954451150103321*g1*ZN1(gI1,0)) * tmp_734;

   return result;
}

std::complex<double> CLASSNAME::CpconjUSebarChiFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_736;
   std::complex<double> tmp_737;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_738;
      std::complex<double> tmp_739;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_739 += KroneckerDelta(gO1,3 + j1)*Ye(j1,j2);
      }
      tmp_738 += tmp_739;
      tmp_737 += (Conj(ZEL(gI2,j2))) * tmp_738;
   }
   tmp_736 += tmp_737;
   result += (-Conj(ZN2(gI1,2))) * tmp_736;

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_740;
   std::complex<double> tmp_742;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_742 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_740 += tmp_742;
   std::complex<double> tmp_741;
   std::complex<double> tmp_743;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_743 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_741 += tmp_743;
   result += (-0.05*Sqr(g1)) * tmp_740 * tmp_741;
   std::complex<double> tmp_744;
   std::complex<double> tmp_746;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_746 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_744 += tmp_746;
   std::complex<double> tmp_745;
   std::complex<double> tmp_747;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_747 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_745 += tmp_747;
   result += (-0.1*Sqr(g1)) * tmp_744 * tmp_745;
   std::complex<double> tmp_748;
   std::complex<double> tmp_750;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_750 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_748 += tmp_750;
   std::complex<double> tmp_749;
   std::complex<double> tmp_751;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_751 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
   }
   tmp_749 += tmp_751;
   result += (-0.05*Sqr(g1)) * tmp_748 * tmp_749;
   std::complex<double> tmp_752;
   std::complex<double> tmp_754;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_754 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_752 += tmp_754;
   std::complex<double> tmp_753;
   std::complex<double> tmp_755;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_755 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
   }
   tmp_753 += tmp_755;
   result += (-0.1*Sqr(g1)) * tmp_752 * tmp_753;
   if (gO1 < 3) {
      std::complex<double> tmp_756;
      std::complex<double> tmp_757;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_757 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_756 += tmp_757;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_756;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_758;
      std::complex<double> tmp_759;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_759 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
      }
      tmp_758 += tmp_759;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_758;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_760;
      std::complex<double> tmp_761;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_761 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
      }
      tmp_760 += tmp_761;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_760;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_762;
      std::complex<double> tmp_763;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_763 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_762 += tmp_763;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_762;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_764;
      std::complex<double> tmp_765;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_765 += Conj(ZD(gI2,j2))*ZD(gI1,j2);
      }
      tmp_764 += tmp_765;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_764;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_766;
      std::complex<double> tmp_767;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_767 += Conj(ZD(gI2,3 + j2))*ZD(gI1,3 + j2);
      }
      tmp_766 += tmp_767;
      result += (0.05*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_766;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_768;
      std::complex<double> tmp_770;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_770 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_768 += tmp_770;
      std::complex<double> tmp_769;
      std::complex<double> tmp_771;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_772;
         std::complex<double> tmp_773;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_773 += Conj(Yd(j3,j4))*Conj(ZD(gI2,3 + j3));
         }
         tmp_772 += tmp_773;
         tmp_771 += (ZD(gI1,j4)) * tmp_772;
      }
      tmp_769 += tmp_771;
      result += (-1) * tmp_768 * tmp_769;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_774;
      std::complex<double> tmp_776;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_777;
         std::complex<double> tmp_778;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_778 += Yd(j1,j2)*ZD(gI1,3 + j1);
         }
         tmp_777 += tmp_778;
         tmp_776 += (Conj(ZD(gI2,j2))) * tmp_777;
      }
      tmp_774 += tmp_776;
      std::complex<double> tmp_775;
      std::complex<double> tmp_779;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_779 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_775 += tmp_779;
      result += (-1) * tmp_774 * tmp_775;
   }

   return result;
}

std::complex<double> CLASSNAME::CpUSeconjUSeconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_780;
   std::complex<double> tmp_782;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_782 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
   }
   tmp_780 += tmp_782;
   std::complex<double> tmp_781;
   std::complex<double> tmp_783;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_783 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
   }
   tmp_781 += tmp_783;
   result += (-0.3*Sqr(g1)) * tmp_780 * tmp_781;
   std::complex<double> tmp_784;
   std::complex<double> tmp_786;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_786 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_784 += tmp_786;
   std::complex<double> tmp_785;
   std::complex<double> tmp_787;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_787 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_785 += tmp_787;
   result += (0.15*Sqr(g1)) * tmp_784 * tmp_785;
   std::complex<double> tmp_788;
   std::complex<double> tmp_790;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_790 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_788 += tmp_790;
   std::complex<double> tmp_789;
   std::complex<double> tmp_791;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_791 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_789 += tmp_791;
   result += (-0.3*Sqr(g1)) * tmp_788 * tmp_789;
   std::complex<double> tmp_792;
   std::complex<double> tmp_794;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_794 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_792 += tmp_794;
   std::complex<double> tmp_793;
   std::complex<double> tmp_795;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_795 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
   }
   tmp_793 += tmp_795;
   result += (0.15*Sqr(g1)) * tmp_792 * tmp_793;
   std::complex<double> tmp_796;
   std::complex<double> tmp_798;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_798 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_796 += tmp_798;
   std::complex<double> tmp_797;
   std::complex<double> tmp_799;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_799 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
   }
   tmp_797 += tmp_799;
   result += (-0.3*Sqr(g1)) * tmp_796 * tmp_797;
   std::complex<double> tmp_800;
   std::complex<double> tmp_802;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_802 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_800 += tmp_802;
   std::complex<double> tmp_801;
   std::complex<double> tmp_803;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_803 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
   }
   tmp_801 += tmp_803;
   result += (-0.3*Sqr(g1)) * tmp_800 * tmp_801;
   std::complex<double> tmp_804;
   std::complex<double> tmp_806;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_807;
      std::complex<double> tmp_808;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_808 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_807 += tmp_808;
      tmp_806 += (Conj(ZE(gI2,j2))) * tmp_807;
   }
   tmp_804 += tmp_806;
   std::complex<double> tmp_805;
   std::complex<double> tmp_809;
   for (unsigned j4 = 0; j4 < 3; ++j4) {
      std::complex<double> tmp_810;
      std::complex<double> tmp_811;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_811 += Conj(Ye(j3,j4))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_810 += tmp_811;
      tmp_809 += (ZE(gI1,j4)) * tmp_810;
   }
   tmp_805 += tmp_809;
   result += (-1) * tmp_804 * tmp_805;
   if (gO1 < 3) {
      std::complex<double> tmp_812;
      std::complex<double> tmp_813;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_813 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_812 += tmp_813;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_812;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_814;
      std::complex<double> tmp_815;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_815 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
      }
      tmp_814 += tmp_815;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_814;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_816;
      std::complex<double> tmp_817;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_817 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
      }
      tmp_816 += tmp_817;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_816;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_818;
      std::complex<double> tmp_819;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_819 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_818 += tmp_819;
      result += (-0.075*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_818;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_820;
      std::complex<double> tmp_821;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_821 += Conj(ZE(gI2,j2))*ZE(gI1,j2);
      }
      tmp_820 += tmp_821;
      result += (-0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_820;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_822;
      std::complex<double> tmp_823;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_823 += Conj(ZE(gI2,3 + j2))*ZE(gI1,3 + j2);
      }
      tmp_822 += tmp_823;
      result += (0.15*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_822;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_824;
      std::complex<double> tmp_826;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_826 += KroneckerDelta(gO2,3 + j1)*Ye(j1,gO1);
      }
      tmp_824 += tmp_826;
      std::complex<double> tmp_825;
      std::complex<double> tmp_827;
      for (unsigned j4 = 0; j4 < 3; ++j4) {
         std::complex<double> tmp_828;
         std::complex<double> tmp_829;
         for (unsigned j3 = 0; j3 < 3; ++j3) {
            tmp_829 += Conj(Ye(j3,j4))*Conj(ZE(gI2,3 + j3));
         }
         tmp_828 += tmp_829;
         tmp_827 += (ZE(gI1,j4)) * tmp_828;
      }
      tmp_825 += tmp_827;
      result += (-1) * tmp_824 * tmp_825;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_830;
      std::complex<double> tmp_831;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_831 += KroneckerDelta(gO1,3 + j1)*ZE(gI1,3 + j1);
      }
      tmp_830 += tmp_831;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_830;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_832;
      std::complex<double> tmp_833;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_833 += KroneckerDelta(gO1,3 + j2)*ZE(gI1,3 + j2);
      }
      tmp_832 += tmp_833;
      result += (0.15*Conj(ZE(gI2,gO2))*Sqr(g1)) * tmp_832;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_834;
      std::complex<double> tmp_836;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_837;
         std::complex<double> tmp_838;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_838 += Ye(j1,j2)*ZE(gI1,3 + j1);
         }
         tmp_837 += tmp_838;
         tmp_836 += (Conj(ZE(gI2,j2))) * tmp_837;
      }
      tmp_834 += tmp_836;
      std::complex<double> tmp_835;
      std::complex<double> tmp_839;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_839 += Conj(Ye(j3,gO2))*KroneckerDelta(gO1,3 + j3);
      }
      tmp_835 += tmp_839;
      result += (-1) * tmp_834 * tmp_835;
   }
   if (gO1 < 3 && gO2 < 3) {
      std::complex<double> tmp_840;
      std::complex<double> tmp_842;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_842 += Ye(j1,gO1)*ZE(gI1,3 + j1);
      }
      tmp_840 += tmp_842;
      std::complex<double> tmp_841;
      std::complex<double> tmp_843;
      for (unsigned j3 = 0; j3 < 3; ++j3) {
         tmp_843 += Conj(Ye(j3,gO2))*Conj(ZE(gI2,3 + j3));
      }
      tmp_841 += tmp_843;
      result += (-1) * tmp_840 * tmp_841;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_844;
      std::complex<double> tmp_845;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_845 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
      }
      tmp_844 += tmp_845;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_844;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_846;
      std::complex<double> tmp_847;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_847 += Conj(ZE(gI2,3 + j2))*KroneckerDelta(gO2,3 + j2);
      }
      tmp_846 += tmp_847;
      result += (0.15*Sqr(g1)*ZE(gI1,gO1)) * tmp_846;
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

   std::complex<double> tmp_848;
   std::complex<double> tmp_850;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_850 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_848 += tmp_850;
   std::complex<double> tmp_849;
   std::complex<double> tmp_851;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_851 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_849 += tmp_851;
   result += (-0.05*Sqr(g1)) * tmp_848 * tmp_849;
   std::complex<double> tmp_852;
   std::complex<double> tmp_854;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_854 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_852 += tmp_854;
   std::complex<double> tmp_853;
   std::complex<double> tmp_855;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_855 += KroneckerDelta(gO1,3 + j2)*KroneckerDelta(gO2,3 + j2);
   }
   tmp_853 += tmp_855;
   result += (0.2*Sqr(g1)) * tmp_852 * tmp_853;
   std::complex<double> tmp_856;
   std::complex<double> tmp_858;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_858 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_856 += tmp_858;
   std::complex<double> tmp_857;
   std::complex<double> tmp_859;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_859 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
   }
   tmp_857 += tmp_859;
   result += (-0.05*Sqr(g1)) * tmp_856 * tmp_857;
   std::complex<double> tmp_860;
   std::complex<double> tmp_862;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_862 += KroneckerDelta(gO1,3 + j1)*KroneckerDelta(gO2,3 + j1);
   }
   tmp_860 += tmp_862;
   std::complex<double> tmp_861;
   std::complex<double> tmp_863;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      tmp_863 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
   }
   tmp_861 += tmp_863;
   result += (0.2*Sqr(g1)) * tmp_860 * tmp_861;
   if (gO1 < 3) {
      std::complex<double> tmp_864;
      std::complex<double> tmp_865;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_865 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_864 += tmp_865;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_864;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_866;
      std::complex<double> tmp_867;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_867 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
      }
      tmp_866 += tmp_867;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_866;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_868;
      std::complex<double> tmp_869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_869 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
      }
      tmp_868 += tmp_869;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_868;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_870;
      std::complex<double> tmp_871;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_871 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_870 += tmp_871;
      result += (0.025*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_870;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_872;
      std::complex<double> tmp_873;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_873 += Conj(ZU(gI2,j2))*ZU(gI1,j2);
      }
      tmp_872 += tmp_873;
      result += (0.125*KroneckerDelta(gO1,gO2)*Sqr(g2)) * tmp_872;
   }
   if (gO1 < 3) {
      std::complex<double> tmp_874;
      std::complex<double> tmp_875;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_875 += Conj(ZU(gI2,3 + j2))*ZU(gI1,3 + j2);
      }
      tmp_874 += tmp_875;
      result += (-0.1*KroneckerDelta(gO1,gO2)*Sqr(g1)) * tmp_874;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_876;
      std::complex<double> tmp_877;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_877 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_876 += tmp_877;
      result += (MuD*ZHR(gI2,0)) * tmp_876;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_878;
      std::complex<double> tmp_879;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_879 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_878 += tmp_879;
      result += (0.7071067811865475*LamSD*vS*ZHR(gI2,0)) * tmp_878;
   }
   if (gO2 < 3) {
      std::complex<double> tmp_880;
      std::complex<double> tmp_881;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_881 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_880 += tmp_881;
      result += (0.5*LamTD*vT*ZHR(gI2,0)) * tmp_880;
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSeAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_882;
   std::complex<double> tmp_883;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_884;
      std::complex<double> tmp_885;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_885 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_884 += tmp_885;
      tmp_883 += (Conj(ZE(gI1,j2))) * tmp_884;
   }
   tmp_882 += tmp_883;
   result += (std::complex<double>(0,-0.7071067811865475)*Conj(Mu)*ZA(gI2,1)) *
      tmp_882;
   if (gO2 < 3) {
      std::complex<double> tmp_886;
      std::complex<double> tmp_887;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_887 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_886 += tmp_887;
      result += (std::complex<double>(0,0.7071067811865475)*Mu*ZA(gI2,1)) *
         tmp_886;
   }
   std::complex<double> tmp_888;
   std::complex<double> tmp_889;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_889 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_888 += tmp_889;
   result += (std::complex<double>(0,-0.7745966692414834)*g1*MDBS*ZA(gI2,2)) *
      tmp_888;
   std::complex<double> tmp_890;
   std::complex<double> tmp_891;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_891 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_890 += tmp_891;
   result += (std::complex<double>(0,0.7745966692414834)*g1*Conj(MDBS)*ZA(gI2,2
      )) * tmp_890;
   if (gO2 < 3) {
      result += std::complex<double>(0,0.3872983346207417)*g1*MDBS*Conj(ZE(
         gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,-0.3872983346207417)*g1*Conj(MDBS)*
         Conj(ZE(gI1,gO2))*ZA(gI2,2);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,0.5)*g2*MDWBT*Conj(ZE(gI1,gO2))*ZA(
         gI2,3);
   }
   if (gO2 < 3) {
      result += std::complex<double>(0,-0.5)*g2*Conj(MDWBT)*Conj(ZE(gI1,gO2)
         )*ZA(gI2,3);
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjUSeSehh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_892;
   std::complex<double> tmp_893;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_893 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_892 += tmp_893;
   result += (0.3*vd*Sqr(g1)*ZH(gI2,0)) * tmp_892;
   std::complex<double> tmp_894;
   std::complex<double> tmp_895;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_896;
      std::complex<double> tmp_897;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_898;
         std::complex<double> tmp_899;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_899 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_898 += tmp_899;
         tmp_897 += (KroneckerDelta(gO2,3 + j2)) * tmp_898;
      }
      tmp_896 += tmp_897;
      tmp_895 += (Conj(ZE(gI1,3 + j3))) * tmp_896;
   }
   tmp_894 += tmp_895;
   result += (-(vd*ZH(gI2,0))) * tmp_894;
   if (gO2 < 3) {
      result += -0.15*vd*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      result += 0.25*vd*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,0);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_900;
      std::complex<double> tmp_901;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_902;
         std::complex<double> tmp_903;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_903 += Conj(Ye(j1,gO2))*Ye(j1,j2);
         }
         tmp_902 += tmp_903;
         tmp_901 += (Conj(ZE(gI1,j2))) * tmp_902;
      }
      tmp_900 += tmp_901;
      result += (-(vd*ZH(gI2,0))) * tmp_900;
   }
   std::complex<double> tmp_904;
   std::complex<double> tmp_905;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_905 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_904 += tmp_905;
   result += (-0.3*vu*Sqr(g1)*ZH(gI2,1)) * tmp_904;
   std::complex<double> tmp_906;
   std::complex<double> tmp_907;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_908;
      std::complex<double> tmp_909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_909 += KroneckerDelta(gO2,3 + j1)*Ye(j1,j2);
      }
      tmp_908 += tmp_909;
      tmp_907 += (Conj(ZE(gI1,j2))) * tmp_908;
   }
   tmp_906 += tmp_907;
   result += (0.7071067811865475*Conj(Mu)*ZH(gI2,1)) * tmp_906;
   if (gO2 < 3) {
      result += 0.15*vu*Conj(ZE(gI1,gO2))*Sqr(g1)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      result += -0.25*vu*Conj(ZE(gI1,gO2))*Sqr(g2)*ZH(gI2,1);
   }
   if (gO2 < 3) {
      std::complex<double> tmp_910;
      std::complex<double> tmp_911;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_911 += Conj(Ye(j1,gO2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_910 += tmp_911;
      result += (0.7071067811865475*Mu*ZH(gI2,1)) * tmp_910;
   }
   std::complex<double> tmp_912;
   std::complex<double> tmp_913;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_913 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_912 += tmp_913;
   result += (-0.7745966692414834*g1*MDBS*ZH(gI2,2)) * tmp_912;
   std::complex<double> tmp_914;
   std::complex<double> tmp_915;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_915 += Conj(ZE(gI1,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_914 += tmp_915;
   result += (-0.7745966692414834*g1*Conj(MDBS)*ZH(gI2,2)) * tmp_914;
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*MDBS*Conj(ZE(gI1,gO2))*ZH(gI2,2);
   }
   if (gO2 < 3) {
      result += 0.3872983346207417*g1*Conj(MDBS)*Conj(ZE(gI1,gO2))*ZH(gI2,2)
         ;
   }
   if (gO2 < 3) {
      result += 0.5*g2*MDWBT*Conj(ZE(gI1,gO2))*ZH(gI2,3);
   }
   if (gO2 < 3) {
      result += 0.5*g2*Conj(MDWBT)*Conj(ZE(gI1,gO2))*ZH(gI2,3);
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

   std::complex<double> tmp_916;
   std::complex<double> tmp_917;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_917 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_916 += tmp_917;
   result += (-0.7745966692414834*g1*Cos(ThetaW())) * tmp_916;
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

   std::complex<double> tmp_918;
   std::complex<double> tmp_919;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_919 += Conj(ZE(gI2,3 + j1))*KroneckerDelta(gO2,3 + j1);
   }
   tmp_918 += tmp_919;
   result += (0.7745966692414834*g1*Sin(ThetaW())) * tmp_918;
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
      KroneckerDelta(3,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmgWm(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1) + 4*vT*
      KroneckerDelta(3,gO1))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhbargWmCgWmC(unsigned gO1) const
{
   std::complex<double> result;

   result = -0.25*(vd*KroneckerDelta(0,gO1) + vu*KroneckerDelta(1,gO1) + 4*vT*
      KroneckerDelta(3,gO1))*Sqr(g2);

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
      ,gO1)*KroneckerDelta(1,gO2) + 4*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))
      *Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*KroneckerDelta(2,gO1)
      + 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))*ZHR(gI1,0)*ZHR(gI2,
      0) + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)
      + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))*ZHR(gI1,0)*ZHR(gI2,0) +
      (-(Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + (1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*
      LamTU*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(
      1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(-4*LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTU*KroneckerDelta(3,gO2))))*ZHR(gI1,1)*ZHR(gI2,1)) - KroneckerDelta(1,gO1
      )*(5*KroneckerDelta(0,gO2)*(-2*LamSD*Conj(LamSU)*ZHR(gI1,1)*ZHR(gI2,0) +
      LamTD*Conj(LamTU)*ZHR(gI1,1)*ZHR(gI2,0) + (-2*LamSU*Conj(LamSD) + LamTU*Conj
      (LamTD))*ZHR(gI1,0)*ZHR(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(
      g2))*ZHR(gI1,0)*ZHR(gI2,0) + (20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1
      ) - 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1))) + KroneckerDelta(0,gO1)*(5*
      KroneckerDelta(1,gO2)*(2*LamSD*Conj(LamSU)*ZHR(gI1,1)*ZHR(gI2,0) - LamTD*
      Conj(LamTU)*ZHR(gI1,1)*ZHR(gI2,0) + (2*LamSU*Conj(LamSD) - LamTU*Conj(LamTD)
      )*ZHR(gI1,0)*ZHR(gI2,1)) + KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) - (3*Sqr(g1) +
      5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(Conj(ZRP(gI2,0))*(5*(Conj(LamTU)*(1.4142135623730951*LamSU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSU*
      KroneckerDelta(2,gO1) + 2*LamTU*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)
      ) + Conj(LamSU)*(1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSU*KroneckerDelta(2,gO2)
      + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1)
      *KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2)))*ZRP(gI1,0)
      ) + Conj(ZRP(gI2,1))*(5*(Conj(LamTD)*(1.4142135623730951*LamSD*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*
      KroneckerDelta(2,gO1) - 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)
      ) + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*LamSD*KroneckerDelta(2,gO2
      ) + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1
      )*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)))*ZRP(
      gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjRhRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZHR(gI1,0)*
      ZHR(gI2,0) - 20*vS*AbsSqr(LamSD)*KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(gI2,0)
      - 14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(
      gI2,0) - 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZHR(
      gI1,0)*ZHR(gI2,0) - 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2
      ,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta
      (2,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*
      KroneckerDelta(2,gO2)*ZHR(gI1,0)*ZHR(gI2,0) + 10*g2*MDWBT*KroneckerDelta(3,
      gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZHR(
      gI1,0)*ZHR(gI2,0) - 7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(3
      ,gO2)*ZHR(gI1,0)*ZHR(gI2,0) - 10*MuD*Conj(LamTD)*KroneckerDelta(3,gO2)*ZHR(
      gI1,0)*ZHR(gI2,0) - 7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(3
      ,gO2)*ZHR(gI1,0)*ZHR(gI2,0) + 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZHR(
      gI1,0)*ZHR(gI2,0) - 10*LamTD*Conj(MuD)*KroneckerDelta(3,gO2)*ZHR(gI1,0)*ZHR(
      gI2,0) + 7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZHR(gI1,1)*ZHR(gI2,
      1) - 20*vS*AbsSqr(LamSU)*KroneckerDelta(2,gO2)*ZHR(gI1,1)*ZHR(gI2,1) -
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

std::complex<double> CLASSNAME::CpUhhconjRpmRpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(Conj(ZRP(gI2,0))*(-7.745966692414834*g1*MDBS*KroneckerDelta
      (2,gO2) + 20*vS*AbsSqr(LamSU)*KroneckerDelta(2,gO2) + 14.142135623730951*MuU
      *Conj(LamSU)*KroneckerDelta(2,gO2) + 7.0710678118654755*LamTU*vT*Conj(LamSU)
      *KroneckerDelta(2,gO2) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*
      KroneckerDelta(2,gO2) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2
      ) + 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(2,gO2) - 10*g2*MDWBT*
      KroneckerDelta(3,gO2) + 10*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2) +
      7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2) + 10*MuU*Conj(
      LamTU)*KroneckerDelta(3,gO2) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*
      KroneckerDelta(3,gO2) - 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2) + 10*LamTU*
      Conj(MuU)*KroneckerDelta(3,gO2) + vd*KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*
      Sqr(g2)) + vu*KroneckerDelta(1,gO2)*(20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2
      )))*ZRP(gI1,0)) + Conj(ZRP(gI2,1))*(-7.745966692414834*g1*MDBS*
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
      KroneckerDelta(1,gO2)*(-3*vu*Sqr(g1) + 5*vu*Sqr(g2)))*ZRP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarCha1Cha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(3,gO2)*(-2*g2*UM1(gI2,0)*UP1(gI1,0) + Conj(
      LamTD)*UM1(gI2,1)*UP1(gI1,1)) - 1.4142135623730951*(Conj(LamSD)*
      KroneckerDelta(2,gO2)*UM1(gI2,1)*UP1(gI1,1) + KroneckerDelta(0,gO2)*(g2*UM1(
      gI2,1)*UP1(gI1,0) + Conj(LamTD)*UM1(gI2,0)*UP1(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarCha1Cha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(-(Conj(UM1(gI1,0))*(1.4142135623730951*LamTD*Conj(UP1(gI2,1))*
      KroneckerDelta(0,gO1) + 2*g2*Conj(UP1(gI2,0))*KroneckerDelta(3,gO1))) - Conj
      (UM1(gI1,1))*(1.4142135623730951*g2*Conj(UP1(gI2,0))*KroneckerDelta(0,gO1) +
      Conj(UP1(gI2,1))*(1.4142135623730951*LamSD*KroneckerDelta(2,gO1) - LamTD*
      KroneckerDelta(3,gO1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarCha2Cha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(3,gO2)*(2*g2*UM2(gI1,0)*UP2(gI2,0) + Conj(LamTU
      )*UM2(gI1,1)*UP2(gI2,1)) + 1.4142135623730951*(Conj(LamTU)*KroneckerDelta(1,
      gO2)*UM2(gI1,1)*UP2(gI2,0) + (-(g2*KroneckerDelta(1,gO2)*UM2(gI1,0)) + Conj(
      LamSU)*KroneckerDelta(2,gO2)*UM2(gI1,1))*UP2(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarCha2Cha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Conj(UM2(gI2,0))*(-1.4142135623730951*Conj(UP2(gI1,1))*
      KroneckerDelta(1,gO1) + 2*Conj(UP2(gI1,0))*KroneckerDelta(3,gO1)) + Conj(UM2
      (gI2,1))*(1.4142135623730951*LamTU*Conj(UP2(gI1,0))*KroneckerDelta(1,gO1) +
      Conj(UP2(gI1,1))*(1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + LamTU*
      KroneckerDelta(3,gO1))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjRhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*Mu*(2*Conj(LamSD)*(KroneckerDelta(2,
      gO2)*ZA(gI2,1) - KroneckerDelta(1,gO2)*ZA(gI2,2))*ZHR(gI1,0) +
      1.4142135623730951*Conj(LamTD)*(KroneckerDelta(3,gO2)*ZA(gI2,1) -
      KroneckerDelta(1,gO2)*ZA(gI2,3))*ZHR(gI1,0) + (Conj(LamSU)*(-2*
      KroneckerDelta(2,gO2)*ZA(gI2,0) + 2*KroneckerDelta(0,gO2)*ZA(gI2,2)) +
      1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZA(gI2,0) -
      KroneckerDelta(0,gO2)*ZA(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjRhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*Mu*(2*Conj(LamSD)*(KroneckerDelta(2,gO2)*ZH(gI2,1) +
      KroneckerDelta(1,gO2)*ZH(gI2,2))*ZHR(gI1,0) + 1.4142135623730951*Conj(LamTD)
      *(KroneckerDelta(3,gO2)*ZH(gI2,1) + KroneckerDelta(1,gO2)*ZH(gI2,3))*ZHR(gI1
      ,0) + (-2*Conj(LamSU)*(KroneckerDelta(2,gO2)*ZH(gI2,0) + KroneckerDelta(0,
      gO2)*ZH(gI2,2)) + 1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZH(
      gI2,0) + KroneckerDelta(0,gO2)*ZH(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjRpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(1.4142135623730951*Conj(LamSU)*KroneckerDelta(2,gO2)*Mu*ZP(gI2
      ,0)*ZRP(gI1,0) + Conj(LamTU)*Mu*(KroneckerDelta(3,gO2)*ZP(gI2,0) -
      1.4142135623730951*KroneckerDelta(0,gO2)*ZP(gI2,3))*ZRP(gI1,0) + Conj(Mu)*(
      -1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*ZP(gI2,1) + LamTD*
      KroneckerDelta(3,gO2)*ZP(gI2,1) + 1.4142135623730951*LamTD*KroneckerDelta(1,
      gO2)*ZP(gI2,2))*ZRP(gI1,1));

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

std::complex<double> CLASSNAME::CpUhhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*KroneckerDelta(gI1,gI2)*(2*(3.872983346207417*g1*(MDBS + Conj(
      MDBS))*KroneckerDelta(2,gO2) - 5*g2*(MDWBT + Conj(MDWBT))*KroneckerDelta(3,
      gO2)) - vd*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2)) + vu*KroneckerDelta
      (1,gO2)*(3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_920;
   std::complex<double> tmp_921;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_922;
      std::complex<double> tmp_923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_923 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_922 += tmp_923;
      tmp_921 += (ZDL(gI1,j2)) * tmp_922;
   }
   tmp_920 += tmp_921;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_920;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_924;
   std::complex<double> tmp_925;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_926;
      std::complex<double> tmp_927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_927 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_926 += tmp_927;
      tmp_925 += (Conj(ZDL(gI2,j2))) * tmp_926;
   }
   tmp_924 += tmp_925;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_924;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_928;
   std::complex<double> tmp_929;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_930;
      std::complex<double> tmp_931;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_931 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_930 += tmp_931;
      tmp_929 += (ZEL(gI1,j2)) * tmp_930;
   }
   tmp_928 += tmp_929;
   result += (-0.7071067811865475*KroneckerDelta(0,gO2)) * tmp_928;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_932;
   std::complex<double> tmp_933;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_934;
      std::complex<double> tmp_935;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_935 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_934 += tmp_935;
      tmp_933 += (Conj(ZEL(gI2,j2))) * tmp_934;
   }
   tmp_932 += tmp_933;
   result += (-0.7071067811865475*KroneckerDelta(0,gO1)) * tmp_932;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_936;
   std::complex<double> tmp_937;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_938;
      std::complex<double> tmp_939;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_939 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_938 += tmp_939;
      tmp_937 += (ZUL(gI1,j2)) * tmp_938;
   }
   tmp_936 += tmp_937;
   result += (-0.7071067811865475*KroneckerDelta(1,gO2)) * tmp_936;

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_940;
   std::complex<double> tmp_941;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_942;
      std::complex<double> tmp_943;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_943 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_942 += tmp_943;
      tmp_941 += (Conj(ZUL(gI2,j2))) * tmp_942;
   }
   tmp_940 += tmp_941;
   result += (-0.7071067811865475*KroneckerDelta(1,gO1)) * tmp_940;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*KroneckerDelta(2,gO1)
      + 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))*ZA(gI1,0)*ZA(gI2,0)
      + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)
      + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))*ZA(gI1,0)*ZA(gI2,0) + (-
      (Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(
      3,gO1) + (1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*LamTU*
      KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(
      1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(-4*LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTU*KroneckerDelta(3,gO2))))*ZA(gI1,1)*ZA(gI2,1)) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(
      g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*
      LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(4*LamSD*ZA(gI2,2) +
      1.4142135623730951*LamTD*ZA(gI2,3))) + Conj(LamTD)*(1.4142135623730951*LamSD
      *ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSD*ZA(gI2,2) + 2*
      LamTD*ZA(gI2,3))))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1
      ) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(
      gI2,1) + 5*(Conj(LamTU)*(1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(
      gI1,3)*(1.4142135623730951*LamSU*ZA(gI2,2) - 2*LamTU*ZA(gI2,3))) + Conj(
      LamSU)*(1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSU*
      ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3
      ,gO1)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) - 1.4142135623730951*KroneckerDelta(0,gO2)
      *KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) + 4*AbsSqr(LamSU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,
      gO1)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamSU*Conj(LamTU)*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamTU*Conj(LamSU)*KroneckerDelta(2,gO1)*KroneckerDelta(3,
      gO2)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamSU*Conj(LamTU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) + 2*AbsSqr(
      LamTU)*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) + 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,2)*ZP(
      gI2,1) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZP(gI1,2)*ZP(gI2,1) + 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) + 2*
      LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,3)*ZP(
      gI2,1) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,2) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      2,gO1)*ZP(gI1,1)*ZP(gI2,2) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta
      (1,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,2) + 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) +
      Conj(LamSD)*(-1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(
      3,gO1)*ZP(gI1,0)*ZP(gI2,0) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2
      ,gO2)*ZP(gI1,0)*ZP(gI2,0) + LamTD*(-1.4142135623730951*KroneckerDelta(3,gO2)
      *ZP(gI1,0)*ZP(gI2,0) + 2*KroneckerDelta(0,gO2)*(ZP(gI1,3)*ZP(gI2,0) + ZP(gI1
      ,0)*ZP(gI2,2))))) - 1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(
      3,gO1)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) + 2*LamSU*Conj(LamTU)*KroneckerDelta(1,
      gO2)*KroneckerDelta(2,gO1)*ZP(gI1,1)*ZP(gI2,3) + 1.4142135623730951*AbsSqr(
      LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,3) -
      1.4142135623730951*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,1)*ZP(gI2,3) - 4*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(
      gI1,2)*ZP(gI2,3) + 4*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,3) + Conj(LamTD)*(-1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1)*ZP(gI1,0)*ZP(gI2,0) + LamSD*KroneckerDelta(2,gO1)
      *(-1.4142135623730951*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 2*
      KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,3))) + LamTD*
      KroneckerDelta(3,gO1)*(2*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) +
      1.4142135623730951*KroneckerDelta(0,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,3)
      *ZP(gI2,0) + ZP(gI1,0)*(-ZP(gI2,2) + ZP(gI2,3)))))) - KroneckerDelta(0,gO1)*
      (KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*AbsSqr(LamTD) + Sqr(g2)
      )*ZP(gI1,2)*ZP(gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3))) + 5*(
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 2*
      LamTD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,0) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      2*LamTD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) + Conj(
      LamTD)*(2*LamSD*KroneckerDelta(2,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,0)*ZP(
      gI2,3)) + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*(-(ZP(gI1,2)*ZP(gI2
      ,0)) + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(-ZP(gI2,2) + ZP(gI2,3)))))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1
      ,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP
      (gI1,2)*ZP(gI2,2) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3))) - 5*(
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) + 2*
      LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,1) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) +
      KroneckerDelta(0,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      2*LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) + Conj(
      LamTU)*(2*LamSU*KroneckerDelta(2,gO2)*(ZP(gI1,2)*ZP(gI2,1) + ZP(gI1,1)*ZP(
      gI2,3)) + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*(-(ZP(gI1,2)*ZP(gI2
      ,1)) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(-ZP(gI2,2) + ZP(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(4*AbsSqr(LamSU)*KroneckerDelta(2,gO1)*KroneckerDelta(2,
      gO2)*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZH(gI1,1)*ZH(gI2,1) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,
      gO1)*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
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

std::complex<double> CLASSNAME::CpUhhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(
      gI2,0) - 20*vS*AbsSqr(LamSD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0) -
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*ZA(gI2,0)
      - 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZA(gI1,0)*
      ZA(gI2,0) - 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZA
      (gI1,0)*ZA(gI2,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZA
      (gI1,0)*ZA(gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)
      *ZA(gI1,0)*ZA(gI2,0) - 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0)
      - 10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) -
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(
      gI2,0) - 10*MuD*Conj(LamTD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) -
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(
      gI2,0) - 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) - 10*
      LamTD*Conj(MuD)*KroneckerDelta(3,gO2)*ZA(gI1,0)*ZA(gI2,0) -
      7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) - 20*vS*
      AbsSqr(LamSU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) - 14.142135623730951
      *MuU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(gI2,1) +
      7.0710678118654755*LamTU*vT*Conj(LamSU)*KroneckerDelta(2,gO2)*ZA(gI1,1)*ZA(
      gI2,1) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(2,gO2)*ZA(
      gI1,1)*ZA(gI2,1) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZA(
      gI1,1)*ZA(gI2,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(2,gO2)*
      ZA(gI1,1)*ZA(gI2,1) + 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1)
      - 10*vT*AbsSqr(LamTU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) +
      7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(
      gI2,1) + 10*MuU*Conj(LamTU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) +
      7.0710678118654755*LamSU*vS*Conj(LamTU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(
      gI2,1) + 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) + 10*
      LamTU*Conj(MuU)*KroneckerDelta(3,gO2)*ZA(gI1,1)*ZA(gI2,1) - vd*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(
      g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*
      LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(4*LamSD*ZA(gI2,2) +
      1.4142135623730951*LamTD*ZA(gI2,3))) + Conj(LamTD)*(1.4142135623730951*LamSD
      *ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSD*ZA(gI2,2) + 2*
      LamTD*ZA(gI2,3))))) + vu*KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZA(
      gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(
      LamTU)*(1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(
      1.4142135623730951*LamSU*ZA(gI2,2) - 2*LamTU*ZA(gI2,3))) + Conj(LamSU)*(
      1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSU*ZA(gI2,2)
      + 1.4142135623730951*LamTU*ZA(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(
      gI2,0) - 20*vS*AbsSqr(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0) -
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,0)
      + 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*
      ZP(gI2,0) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZP
      (gI1,0)*ZP(gI2,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZP
      (gI1,0)*ZP(gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)
      *ZP(gI1,0)*ZP(gI2,0) + 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0)
      - 10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) +
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(
      gI2,0) + 10*MuD*Conj(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(
      gI2,0) + 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 10*
      LamTD*Conj(MuD)*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) - 10*LamSD*vd*Conj
      (LamTD)*KroneckerDelta(2,gO2)*ZP(gI1,2)*ZP(gI2,0) + 7.0710678118654755*vd*
      AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,2)*ZP(gI2,0) - 7.0710678118654755
      *vd*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) - 10*LamTD*vd*Conj(
      LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,0) - 7.0710678118654755*vd*
      AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZP(gI1,3)*ZP(gI2,0) + 7.0710678118654755
      *vd*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) - 7.745966692414834*g1
      *MDBS*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) - 20*vS*AbsSqr(LamSU)*
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

std::complex<double> CLASSNAME::CpUhhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(-7.0710678118654755*LamSD*vd*Conj(
      LamTD)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZH(gI1,0) + 7.0710678118654755*LamSD*
      vd*Conj(LamTD)*KroneckerDelta(2,gO2)*ZA(gI2,3)*ZH(gI1,0) +
      7.0710678118654755*LamTD*vd*Conj(LamSD)*(KroneckerDelta(3,gO2)*ZA(gI2,2) -
      KroneckerDelta(2,gO2)*ZA(gI2,3))*ZH(gI1,0) - 7.745966692414834*g1*MDBS*
      KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,1) + 14.142135623730951*MuU*Conj(
      LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,1) - 7.0710678118654755*LamTU*
      vT*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,1) +
      7.0710678118654755*LamSU*vT*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(
      gI1,1) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(
      gI1,1) - 14.142135623730951*LamSU*Conj(MuU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*
      ZH(gI1,1) - 7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(3,gO2)*ZA
      (gI2,2)*ZH(gI1,1) + 7.0710678118654755*LamSU*vu*Conj(LamTU)*KroneckerDelta(3
      ,gO2)*ZA(gI2,2)*ZH(gI1,1) + 10*g2*MDWBT*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(
      gI1,1) + 7.0710678118654755*LamTU*vS*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(
      gI2,3)*ZH(gI1,1) - 10*MuU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1
      ,1) - 7.0710678118654755*LamSU*vS*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI2,3
      )*ZH(gI1,1) - 10*g2*Conj(MDWBT)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) +
      10*LamTU*Conj(MuU)*KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,1) +
      7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(2,gO2)*ZA(gI2,3)*ZH(
      gI1,1) - 7.0710678118654755*LamSU*vu*Conj(LamTU)*KroneckerDelta(2,gO2)*ZA(
      gI2,3)*ZH(gI1,1) + 7.0710678118654755*LamTU*vu*Conj(LamSU)*KroneckerDelta(1,
      gO2)*ZA(gI2,3)*ZH(gI1,2) - 7.0710678118654755*LamSU*vu*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZA(gI2,3)*ZH(gI1,2) - 7.0710678118654755*LamTU*vu*Conj
      (LamSU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,3) + 7.0710678118654755*LamSU
      *vu*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI2,2)*ZH(gI1,3) + KroneckerDelta(0
      ,gO2)*(-5*ZA(gI2,3)*((2*g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) -
      2*MuD*Conj(LamTD) - 1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(
      MDWBT) + 2*LamTD*Conj(MuD))*ZH(gI1,0) + 1.4142135623730951*vd*(LamTD*Conj(
      LamSD) - LamSD*Conj(LamTD))*ZH(gI1,2)) + ZA(gI2,2)*((7.745966692414834*g1*
      MDBS + 7.0710678118654755*(2*MuD + LamTD*vT)*Conj(LamSD) -
      7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS) -
      14.142135623730951*LamSD*Conj(MuD))*ZH(gI1,0) + 7.0710678118654755*vd*(LamTD
      *Conj(LamSD) - LamSD*Conj(LamTD))*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(
      gI2,0) - 20*vS*AbsSqr(LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0) -
      14.142135623730951*MuD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*ZH(gI2,0)
      - 7.0710678118654755*LamTD*vT*Conj(LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,0)*
      ZH(gI2,0) - 7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(2,gO2)*ZH
      (gI1,0)*ZH(gI2,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)*ZH
      (gI1,0)*ZH(gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(2,gO2)
      *ZH(gI1,0)*ZH(gI2,0) - 10*g2*MDWBT*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0)
      - 10*vT*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) -
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(
      gI2,0) - 10*MuD*Conj(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) -
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(
      gI2,0) - 10*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 10*
      LamTD*Conj(MuD)*KroneckerDelta(3,gO2)*ZH(gI1,0)*ZH(gI2,0) - 20*vd*AbsSqr(
      LamSD)*KroneckerDelta(2,gO2)*ZH(gI1,2)*ZH(gI2,0) - 7.0710678118654755*LamTD*
      vd*Conj(LamSD)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(gI2,0) -
      7.0710678118654755*LamSD*vd*Conj(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,2)*ZH(
      gI2,0) - 7.0710678118654755*LamTD*vd*Conj(LamSD)*KroneckerDelta(2,gO2)*ZH(
      gI1,3)*ZH(gI2,0) - 7.0710678118654755*LamSD*vd*Conj(LamTD)*KroneckerDelta(2,
      gO2)*ZH(gI1,3)*ZH(gI2,0) - 10*vd*AbsSqr(LamTD)*KroneckerDelta(3,gO2)*ZH(gI1,
      3)*ZH(gI2,0) - 7.745966692414834*g1*MDBS*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(
      gI2,1) - 20*vS*AbsSqr(LamSU)*KroneckerDelta(2,gO2)*ZH(gI1,1)*ZH(gI2,1) -
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

std::complex<double> CLASSNAME::CpUhhbarChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(5*Conj(LamTD)*KroneckerDelta(0,gO2)*ZN1(gI1,2)*ZN2(gI2,1) + 5*
      Conj(LamTU)*KroneckerDelta(1,gO2)*ZN1(gI1,3)*ZN2(gI2,1) + 3.872983346207417*
      g1*KroneckerDelta(0,gO2)*ZN1(gI1,0)*ZN2(gI2,2) - 5*g2*KroneckerDelta(0,gO2)*
      ZN1(gI1,1)*ZN2(gI2,2) + 5*Conj(LamTD)*KroneckerDelta(3,gO2)*ZN1(gI1,2)*ZN2(
      gI2,2) + 7.0710678118654755*Conj(LamSD)*ZN1(gI1,2)*(KroneckerDelta(0,gO2)*
      ZN2(gI2,0) + KroneckerDelta(2,gO2)*ZN2(gI2,2)) - 3.872983346207417*g1*
      KroneckerDelta(1,gO2)*ZN1(gI1,0)*ZN2(gI2,3) + 5*g2*KroneckerDelta(1,gO2)*ZN1
      (gI1,1)*ZN2(gI2,3) + 5*Conj(LamTU)*KroneckerDelta(3,gO2)*ZN1(gI1,3)*ZN2(gI2,
      3) - 7.0710678118654755*Conj(LamSU)*ZN1(gI1,3)*(KroneckerDelta(1,gO2)*ZN2(
      gI2,0) + KroneckerDelta(2,gO2)*ZN2(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUhhbarChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(3.872983346207417*g1*Conj(ZN1(gI2,0))*(Conj(ZN2(gI1,2))*
      KroneckerDelta(0,gO1) - Conj(ZN2(gI1,3))*KroneckerDelta(1,gO1)) + 5*Conj(ZN1
      (gI2,2))*(1.4142135623730951*LamSD*Conj(ZN2(gI1,0))*KroneckerDelta(0,gO1) +
      LamTD*Conj(ZN2(gI1,1))*KroneckerDelta(0,gO1) + Conj(ZN2(gI1,2))*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + LamTD*KroneckerDelta(3,gO1)
      )) + 5*(Conj(ZN1(gI2,1))*(-(g2*Conj(ZN2(gI1,2))*KroneckerDelta(0,gO1)) + g2*
      Conj(ZN2(gI1,3))*KroneckerDelta(1,gO1)) + Conj(ZN1(gI2,3))*(
      -1.4142135623730951*LamSU*Conj(ZN2(gI1,0))*KroneckerDelta(1,gO1) + LamTU*
      Conj(ZN2(gI1,1))*KroneckerDelta(1,gO1) + Conj(ZN2(gI1,3))*(
      -1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + LamTU*KroneckerDelta(3,gO1
      )))));

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_944;
   std::complex<double> tmp_945;
   std::complex<double> tmp_946;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_946 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_945 += tmp_946;
   tmp_944 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_945;
   std::complex<double> tmp_947;
   std::complex<double> tmp_948;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_948 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_947 += tmp_948;
   tmp_944 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_947;
   std::complex<double> tmp_949;
   std::complex<double> tmp_950;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_950 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_949 += tmp_950;
   tmp_944 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_949;
   std::complex<double> tmp_951;
   std::complex<double> tmp_952;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_952 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_951 += tmp_952;
   tmp_944 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_951;
   std::complex<double> tmp_953;
   std::complex<double> tmp_954;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_954 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_953 += tmp_954;
   tmp_944 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*Sqr(g1)) * tmp_953;
   std::complex<double> tmp_955;
   std::complex<double> tmp_956;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_956 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_955 += tmp_956;
   tmp_944 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_955;
   std::complex<double> tmp_957;
   std::complex<double> tmp_958;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_959;
      std::complex<double> tmp_960;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_961;
         std::complex<double> tmp_962;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_962 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_961 += tmp_962;
         tmp_960 += (ZD(gI1,3 + j2)) * tmp_961;
      }
      tmp_959 += tmp_960;
      tmp_958 += (Conj(ZD(gI2,3 + j3))) * tmp_959;
   }
   tmp_957 += tmp_958;
   tmp_944 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_957;
   std::complex<double> tmp_963;
   std::complex<double> tmp_964;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_965;
      std::complex<double> tmp_966;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_967;
         std::complex<double> tmp_968;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_968 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_967 += tmp_968;
         tmp_966 += (Conj(ZD(gI2,j2))) * tmp_967;
      }
      tmp_965 += tmp_966;
      tmp_964 += (ZD(gI1,j3)) * tmp_965;
   }
   tmp_963 += tmp_964;
   tmp_944 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_963;
   result += (std::complex<double>(0,-1)) * tmp_944;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_969;
   std::complex<double> tmp_970;
   std::complex<double> tmp_971;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_971 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_970 += tmp_971;
   tmp_969 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_970;
   std::complex<double> tmp_972;
   std::complex<double> tmp_973;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_973 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_972 += tmp_973;
   tmp_969 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_972;
   std::complex<double> tmp_974;
   std::complex<double> tmp_975;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_975 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_974 += tmp_975;
   tmp_969 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_974;
   std::complex<double> tmp_976;
   std::complex<double> tmp_977;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_977 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_976 += tmp_977;
   tmp_969 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_976;
   std::complex<double> tmp_978;
   std::complex<double> tmp_979;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_979 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_978 += tmp_979;
   tmp_969 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)*Sqr(g1)) * tmp_978;
   std::complex<double> tmp_980;
   std::complex<double> tmp_981;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_981 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_980 += tmp_981;
   tmp_969 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_980;
   std::complex<double> tmp_982;
   std::complex<double> tmp_983;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_984;
      std::complex<double> tmp_985;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_986;
         std::complex<double> tmp_987;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_987 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_986 += tmp_987;
         tmp_985 += (ZE(gI1,3 + j2)) * tmp_986;
      }
      tmp_984 += tmp_985;
      tmp_983 += (Conj(ZE(gI2,3 + j3))) * tmp_984;
   }
   tmp_982 += tmp_983;
   tmp_969 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_982;
   std::complex<double> tmp_988;
   std::complex<double> tmp_989;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_990;
      std::complex<double> tmp_991;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_992;
         std::complex<double> tmp_993;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_993 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_992 += tmp_993;
         tmp_991 += (Conj(ZE(gI2,j2))) * tmp_992;
      }
      tmp_990 += tmp_991;
      tmp_989 += (ZE(gI1,j3)) * tmp_990;
   }
   tmp_988 += tmp_989;
   tmp_969 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta(
      0,gO2)) * tmp_988;
   result += (std::complex<double>(0,-1)) * tmp_969;

   return result;
}

std::complex<double> CLASSNAME::CpUhhUhhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_994;
   std::complex<double> tmp_995;
   std::complex<double> tmp_996;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_996 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_995 += tmp_996;
   tmp_994 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_995;
   std::complex<double> tmp_997;
   std::complex<double> tmp_998;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_998 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_997 += tmp_998;
   tmp_994 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_997;
   std::complex<double> tmp_999;
   std::complex<double> tmp_1000;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1000 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_999 += tmp_1000;
   tmp_994 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_999;
   std::complex<double> tmp_1001;
   std::complex<double> tmp_1002;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1002 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1001 += tmp_1002;
   tmp_994 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1001;
   std::complex<double> tmp_1003;
   std::complex<double> tmp_1004;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1004 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1003 += tmp_1004;
   tmp_994 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1003;
   std::complex<double> tmp_1005;
   std::complex<double> tmp_1006;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1006 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1005 += tmp_1006;
   tmp_994 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)*Sqr(g1)) * tmp_1005;
   std::complex<double> tmp_1007;
   std::complex<double> tmp_1008;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1009;
      std::complex<double> tmp_1010;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1011;
         std::complex<double> tmp_1012;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1012 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1011 += tmp_1012;
         tmp_1010 += (ZU(gI1,3 + j2)) * tmp_1011;
      }
      tmp_1009 += tmp_1010;
      tmp_1008 += (Conj(ZU(gI2,3 + j3))) * tmp_1009;
   }
   tmp_1007 += tmp_1008;
   tmp_994 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta(
      1,gO2)) * tmp_1007;
   std::complex<double> tmp_1013;
   std::complex<double> tmp_1014;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1015;
      std::complex<double> tmp_1016;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1017;
         std::complex<double> tmp_1018;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1018 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1017 += tmp_1018;
         tmp_1016 += (Conj(ZU(gI2,j2))) * tmp_1017;
      }
      tmp_1015 += tmp_1016;
      tmp_1014 += (ZU(gI1,j3)) * tmp_1015;
   }
   tmp_1013 += tmp_1014;
   tmp_994 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta(
      1,gO2)) * tmp_1013;
   result += (std::complex<double>(0,-1)) * tmp_994;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1019;
   std::complex<double> tmp_1020;
   std::complex<double> tmp_1021;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1021 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1020 += tmp_1021;
   tmp_1019 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1020;
   std::complex<double> tmp_1022;
   std::complex<double> tmp_1023;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1023 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1022 += tmp_1023;
   tmp_1019 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1022;
   std::complex<double> tmp_1024;
   std::complex<double> tmp_1025;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1025 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1024 += tmp_1025;
   tmp_1019 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1024;
   std::complex<double> tmp_1026;
   std::complex<double> tmp_1027;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1027 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1026 += tmp_1027;
   tmp_1019 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1026;
   std::complex<double> tmp_1028;
   std::complex<double> tmp_1029;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1029 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1028 += tmp_1029;
   tmp_1019 += (std::complex<double>(0,-0.12909944487358055)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1028;
   std::complex<double> tmp_1030;
   std::complex<double> tmp_1031;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1031 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1030 += tmp_1031;
   tmp_1019 += (std::complex<double>(0,-0.12909944487358055)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1030;
   std::complex<double> tmp_1032;
   std::complex<double> tmp_1033;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1033 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1032 += tmp_1033;
   tmp_1019 += (std::complex<double>(0,0.5)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1032;
   std::complex<double> tmp_1034;
   std::complex<double> tmp_1035;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1035 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1034 += tmp_1035;
   tmp_1019 += (std::complex<double>(0,0.5)*g2*Conj(MDWBT)*KroneckerDelta(3,gO2
      )) * tmp_1034;
   std::complex<double> tmp_1036;
   std::complex<double> tmp_1037;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1037 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1036 += tmp_1037;
   tmp_1019 += (std::complex<double>(0,0.1)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1036;
   std::complex<double> tmp_1038;
   std::complex<double> tmp_1039;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1039 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1038 += tmp_1039;
   tmp_1019 += (std::complex<double>(0,-0.1)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1038;
   std::complex<double> tmp_1040;
   std::complex<double> tmp_1041;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1041 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1040 += tmp_1041;
   tmp_1019 += (std::complex<double>(0,-0.2581988897471611)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1040;
   std::complex<double> tmp_1042;
   std::complex<double> tmp_1043;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1043 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1042 += tmp_1043;
   tmp_1019 += (std::complex<double>(0,-0.2581988897471611)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1042;
   std::complex<double> tmp_1044;
   std::complex<double> tmp_1045;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1046;
      std::complex<double> tmp_1047;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1047 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1046 += tmp_1047;
      tmp_1045 += (Conj(ZD(gI2,j2))) * tmp_1046;
   }
   tmp_1044 += tmp_1045;
   tmp_1019 += (std::complex<double>(0,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1044;
   std::complex<double> tmp_1048;
   std::complex<double> tmp_1049;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1050;
      std::complex<double> tmp_1051;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1051 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1050 += tmp_1051;
      tmp_1049 += (ZD(gI1,j2)) * tmp_1050;
   }
   tmp_1048 += tmp_1049;
   tmp_1019 += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2
      )*Mu) * tmp_1048;
   std::complex<double> tmp_1052;
   std::complex<double> tmp_1053;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1054;
      std::complex<double> tmp_1055;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1056;
         std::complex<double> tmp_1057;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1057 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1056 += tmp_1057;
         tmp_1055 += (ZD(gI1,3 + j2)) * tmp_1056;
      }
      tmp_1054 += tmp_1055;
      tmp_1053 += (Conj(ZD(gI2,3 + j3))) * tmp_1054;
   }
   tmp_1052 += tmp_1053;
   tmp_1019 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1052
      ;
   std::complex<double> tmp_1058;
   std::complex<double> tmp_1059;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1060;
      std::complex<double> tmp_1061;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1062;
         std::complex<double> tmp_1063;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1063 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1062 += tmp_1063;
         tmp_1061 += (Conj(ZD(gI2,j2))) * tmp_1062;
      }
      tmp_1060 += tmp_1061;
      tmp_1059 += (ZD(gI1,j3)) * tmp_1060;
   }
   tmp_1058 += tmp_1059;
   tmp_1019 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1058
      ;
   result += (std::complex<double>(0,-1)) * tmp_1019;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1064;
   std::complex<double> tmp_1065;
   std::complex<double> tmp_1066;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1066 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1065 += tmp_1066;
   tmp_1064 += (std::complex<double>(0,-0.15)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1065;
   std::complex<double> tmp_1067;
   std::complex<double> tmp_1068;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1068 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1067 += tmp_1068;
   tmp_1064 += (std::complex<double>(0,0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1067;
   std::complex<double> tmp_1069;
   std::complex<double> tmp_1070;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1070 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1069 += tmp_1070;
   tmp_1064 += (std::complex<double>(0,0.15)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1069;
   std::complex<double> tmp_1071;
   std::complex<double> tmp_1072;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1072 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1071 += tmp_1072;
   tmp_1064 += (std::complex<double>(0,-0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1071;
   std::complex<double> tmp_1073;
   std::complex<double> tmp_1074;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1074 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1073 += tmp_1074;
   tmp_1064 += (std::complex<double>(0,0.3872983346207417)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1073;
   std::complex<double> tmp_1075;
   std::complex<double> tmp_1076;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1076 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1075 += tmp_1076;
   tmp_1064 += (std::complex<double>(0,0.3872983346207417)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1075;
   std::complex<double> tmp_1077;
   std::complex<double> tmp_1078;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1078 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1077 += tmp_1078;
   tmp_1064 += (std::complex<double>(0,0.5)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1077;
   std::complex<double> tmp_1079;
   std::complex<double> tmp_1080;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1080 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1079 += tmp_1080;
   tmp_1064 += (std::complex<double>(0,0.5)*g2*Conj(MDWBT)*KroneckerDelta(3,gO2
      )) * tmp_1079;
   std::complex<double> tmp_1081;
   std::complex<double> tmp_1082;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1082 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1081 += tmp_1082;
   tmp_1064 += (std::complex<double>(0,0.3)*vd*KroneckerDelta(0,gO2)*Sqr(g1)) *
      tmp_1081;
   std::complex<double> tmp_1083;
   std::complex<double> tmp_1084;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1084 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1083 += tmp_1084;
   tmp_1064 += (std::complex<double>(0,-0.3)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1083;
   std::complex<double> tmp_1085;
   std::complex<double> tmp_1086;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1086 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1085 += tmp_1086;
   tmp_1064 += (std::complex<double>(0,-0.7745966692414834)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1085;
   std::complex<double> tmp_1087;
   std::complex<double> tmp_1088;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1088 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1087 += tmp_1088;
   tmp_1064 += (std::complex<double>(0,-0.7745966692414834)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1087;
   std::complex<double> tmp_1089;
   std::complex<double> tmp_1090;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1091;
      std::complex<double> tmp_1092;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1092 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1091 += tmp_1092;
      tmp_1090 += (Conj(ZE(gI2,j2))) * tmp_1091;
   }
   tmp_1089 += tmp_1090;
   tmp_1064 += (std::complex<double>(0,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(1,gO2)) * tmp_1089;
   std::complex<double> tmp_1093;
   std::complex<double> tmp_1094;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1095;
      std::complex<double> tmp_1096;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1096 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1095 += tmp_1096;
      tmp_1094 += (ZE(gI1,j2)) * tmp_1095;
   }
   tmp_1093 += tmp_1094;
   tmp_1064 += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2
      )*Mu) * tmp_1093;
   std::complex<double> tmp_1097;
   std::complex<double> tmp_1098;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1099;
      std::complex<double> tmp_1100;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1101;
         std::complex<double> tmp_1102;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1102 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1101 += tmp_1102;
         tmp_1100 += (ZE(gI1,3 + j2)) * tmp_1101;
      }
      tmp_1099 += tmp_1100;
      tmp_1098 += (Conj(ZE(gI2,3 + j3))) * tmp_1099;
   }
   tmp_1097 += tmp_1098;
   tmp_1064 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1097
      ;
   std::complex<double> tmp_1103;
   std::complex<double> tmp_1104;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1105;
      std::complex<double> tmp_1106;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1107;
         std::complex<double> tmp_1108;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1108 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1107 += tmp_1108;
         tmp_1106 += (Conj(ZE(gI2,j2))) * tmp_1107;
      }
      tmp_1105 += tmp_1106;
      tmp_1104 += (ZE(gI1,j3)) * tmp_1105;
   }
   tmp_1103 += tmp_1104;
   tmp_1064 += (std::complex<double>(0,-1)*vd*KroneckerDelta(0,gO2)) * tmp_1103
      ;
   result += (std::complex<double>(0,-1)) * tmp_1064;

   return result;
}

std::complex<double> CLASSNAME::CpUhhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1109;
   std::complex<double> tmp_1110;
   std::complex<double> tmp_1111;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1111 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1110 += tmp_1111;
   tmp_1109 += (std::complex<double>(0,0.05)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1110;
   std::complex<double> tmp_1112;
   std::complex<double> tmp_1113;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1113 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1112 += tmp_1113;
   tmp_1109 += (std::complex<double>(0,-0.25)*vd*KroneckerDelta(0,gO2)*Sqr(g2))
      * tmp_1112;
   std::complex<double> tmp_1114;
   std::complex<double> tmp_1115;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1115 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1114 += tmp_1115;
   tmp_1109 += (std::complex<double>(0,-0.05)*vu*KroneckerDelta(1,gO2)*Sqr(g1))
      * tmp_1114;
   std::complex<double> tmp_1116;
   std::complex<double> tmp_1117;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1117 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1116 += tmp_1117;
   tmp_1109 += (std::complex<double>(0,0.25)*vu*KroneckerDelta(1,gO2)*Sqr(g2))
      * tmp_1116;
   std::complex<double> tmp_1118;
   std::complex<double> tmp_1119;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1119 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1118 += tmp_1119;
   tmp_1109 += (std::complex<double>(0,-0.12909944487358055)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1118;
   std::complex<double> tmp_1120;
   std::complex<double> tmp_1121;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1121 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1120 += tmp_1121;
   tmp_1109 += (std::complex<double>(0,-0.12909944487358055)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1120;
   std::complex<double> tmp_1122;
   std::complex<double> tmp_1123;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1123 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1122 += tmp_1123;
   tmp_1109 += (std::complex<double>(0,-0.5)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1122;
   std::complex<double> tmp_1124;
   std::complex<double> tmp_1125;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1125 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1124 += tmp_1125;
   tmp_1109 += (std::complex<double>(0,-0.5)*g2*Conj(MDWBT)*KroneckerDelta(3,
      gO2)) * tmp_1124;
   std::complex<double> tmp_1126;
   std::complex<double> tmp_1127;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1127 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1126 += tmp_1127;
   tmp_1109 += (std::complex<double>(0,-0.2)*vd*KroneckerDelta(0,gO2)*Sqr(g1))
      * tmp_1126;
   std::complex<double> tmp_1128;
   std::complex<double> tmp_1129;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1129 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1128 += tmp_1129;
   tmp_1109 += (std::complex<double>(0,0.2)*vu*KroneckerDelta(1,gO2)*Sqr(g1)) *
      tmp_1128;
   std::complex<double> tmp_1130;
   std::complex<double> tmp_1131;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1131 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1130 += tmp_1131;
   tmp_1109 += (std::complex<double>(0,0.5163977794943222)*g1*MDBS*
      KroneckerDelta(2,gO2)) * tmp_1130;
   std::complex<double> tmp_1132;
   std::complex<double> tmp_1133;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1133 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1132 += tmp_1133;
   tmp_1109 += (std::complex<double>(0,0.5163977794943222)*g1*Conj(MDBS)*
      KroneckerDelta(2,gO2)) * tmp_1132;
   std::complex<double> tmp_1134;
   std::complex<double> tmp_1135;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1136;
      std::complex<double> tmp_1137;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1137 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1136 += tmp_1137;
      tmp_1135 += (Conj(ZU(gI2,j2))) * tmp_1136;
   }
   tmp_1134 += tmp_1135;
   tmp_1109 += (std::complex<double>(0,0.7071067811865475)*Conj(Mu)*
      KroneckerDelta(0,gO2)) * tmp_1134;
   std::complex<double> tmp_1138;
   std::complex<double> tmp_1139;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1140;
      std::complex<double> tmp_1141;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1141 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1140 += tmp_1141;
      tmp_1139 += (ZU(gI1,j2)) * tmp_1140;
   }
   tmp_1138 += tmp_1139;
   tmp_1109 += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2
      )*Mu) * tmp_1138;
   std::complex<double> tmp_1142;
   std::complex<double> tmp_1143;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1144;
      std::complex<double> tmp_1145;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1146;
         std::complex<double> tmp_1147;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1147 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1146 += tmp_1147;
         tmp_1145 += (ZU(gI1,3 + j2)) * tmp_1146;
      }
      tmp_1144 += tmp_1145;
      tmp_1143 += (Conj(ZU(gI2,3 + j3))) * tmp_1144;
   }
   tmp_1142 += tmp_1143;
   tmp_1109 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1142
      ;
   std::complex<double> tmp_1148;
   std::complex<double> tmp_1149;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1150;
      std::complex<double> tmp_1151;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1152;
         std::complex<double> tmp_1153;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1153 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1152 += tmp_1153;
         tmp_1151 += (Conj(ZU(gI2,j2))) * tmp_1152;
      }
      tmp_1150 += tmp_1151;
      tmp_1149 += (ZU(gI1,j3)) * tmp_1150;
   }
   tmp_1148 += tmp_1149;
   tmp_1109 += (std::complex<double>(0,-1)*vu*KroneckerDelta(1,gO2)) * tmp_1148
      ;
   result += (std::complex<double>(0,-1)) * tmp_1109;

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
      gI2,1) + 1.4142135623730951*KroneckerDelta(3,gO2)*(ZP(gI2,2) + ZP(gI2,3)));

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
      ,gO1)*KroneckerDelta(1,gO2) + 4*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2))
      *Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*KroneckerDelta(2,gO1)
      + 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))*ZHR(gI1,0)*ZHR(gI2,
      0) + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)
      + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))*ZHR(gI1,0)*ZHR(gI2,0) +
      (-(Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + (1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*
      LamTU*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(
      1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(-4*LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTU*KroneckerDelta(3,gO2))))*ZHR(gI1,1)*ZHR(gI2,1)) - KroneckerDelta(1,gO1
      )*(5*KroneckerDelta(0,gO2)*(-2*LamSD*Conj(LamSU)*ZHR(gI1,1)*ZHR(gI2,0) +
      LamTD*Conj(LamTU)*ZHR(gI1,1)*ZHR(gI2,0) + (-2*LamSU*Conj(LamSD) + LamTU*Conj
      (LamTD))*ZHR(gI1,0)*ZHR(gI2,1)) + KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(
      g2))*ZHR(gI1,0)*ZHR(gI2,0) + (20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1
      ) - 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1))) + KroneckerDelta(0,gO1)*(5*
      KroneckerDelta(1,gO2)*(2*LamSD*Conj(LamSU)*ZHR(gI1,1)*ZHR(gI2,0) - LamTD*
      Conj(LamTU)*ZHR(gI1,1)*ZHR(gI2,0) + (2*LamSU*Conj(LamSD) - LamTU*Conj(LamTD)
      )*ZHR(gI1,0)*ZHR(gI2,1)) + KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) - (3*Sqr(g1) +
      5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(Conj(ZRP(gI2,0))*(5*(Conj(LamTU)*(1.4142135623730951*LamSU*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSU*
      KroneckerDelta(2,gO1) + 2*LamTU*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)
      ) + Conj(LamSU)*(1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSU*KroneckerDelta(2,gO2)
      + 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1)
      *KroneckerDelta(0,gO2)*(3*Sqr(g1) - 5*Sqr(g2)) + KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2)))*ZRP(gI1,0)
      ) + Conj(ZRP(gI2,1))*(5*(Conj(LamTD)*(1.4142135623730951*LamSD*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*
      KroneckerDelta(2,gO1) - 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2)
      ) + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(-4*LamSD*KroneckerDelta(2,gO2
      ) + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))) + KroneckerDelta(0,gO1
      )*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)))*ZRP(
      gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjRhRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.05)*(KroneckerDelta(2,gO2)*((
      7.745966692414834*g1*MDBS - 7.0710678118654755*(2*MuD + LamTD*vT)*Conj(LamSD
      ) + 7.0710678118654755*LamSD*vT*Conj(LamTD) - 7.745966692414834*g1*Conj(MDBS
      ) + 14.142135623730951*LamSD*Conj(MuD))*ZHR(gI1,0)*ZHR(gI2,0) + (
      -7.745966692414834*g1*MDBS + 7.0710678118654755*(-2*MuU + LamTU*vT)*Conj(
      LamSU) - 7.0710678118654755*LamSU*vT*Conj(LamTU) + 7.745966692414834*g1*Conj
      (MDBS) + 14.142135623730951*LamSU*Conj(MuU))*ZHR(gI1,1)*ZHR(gI2,1)) - 5*((vu
      *KroneckerDelta(0,gO2) - vd*KroneckerDelta(1,gO2))*(2*LamSD*Conj(LamSU)*ZHR(
      gI1,1)*ZHR(gI2,0) - LamTD*Conj(LamTU)*ZHR(gI1,1)*ZHR(gI2,0) + (-2*LamSU*Conj
      (LamSD) + LamTU*Conj(LamTD))*ZHR(gI1,0)*ZHR(gI2,1)) + KroneckerDelta(3,gO2)*
      ((2*g2*MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD +
      1.4142135623730951*LamSD*vS)*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(
      MuD))*ZHR(gI1,0)*ZHR(gI2,0) + (-2*g2*MDWBT + 1.4142135623730951*LamTU*vS*
      Conj(LamSU) - (2*MuU + 1.4142135623730951*LamSU*vS)*Conj(LamTU) + 2*g2*Conj(
      MDWBT) + 2*LamTU*Conj(MuU))*ZHR(gI1,1)*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjRpmRpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(Conj(ZRP(gI2,0))*((7.745966692414834*
      g1*MDBS + 7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU) -
      7.0710678118654755*LamSU*vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS) -
      14.142135623730951*LamSU*Conj(MuU))*KroneckerDelta(2,gO2) + 5*(2*g2*MDWBT -
      1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) +
      1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(
      MuU))*KroneckerDelta(3,gO2))*ZRP(gI1,0) + Conj(ZRP(gI2,1))*((
      -7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT)*Conj(
      LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) + 7.745966692414834*g1*Conj
      (MDBS) - 14.142135623730951*LamSD*Conj(MuD))*KroneckerDelta(2,gO2) - 5*(2*g2
      *MDWBT - 1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) +
      1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(
      MuD))*KroneckerDelta(3,gO2))*ZRP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarCha1Cha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(KroneckerDelta(3,gO2)*(2*g2*UM1(gI2,0
      )*UP1(gI1,0) + Conj(LamTD)*UM1(gI2,1)*UP1(gI1,1)) + 1.4142135623730951*(-(
      Conj(LamSD)*KroneckerDelta(2,gO2)*UM1(gI2,1)*UP1(gI1,1)) + KroneckerDelta(0,
      gO2)*(g2*UM1(gI2,1)*UP1(gI1,0) - Conj(LamTD)*UM1(gI2,0)*UP1(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarCha1Cha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(UM1(gI1,0))*(-1.4142135623730951*
      LamTD*Conj(UP1(gI2,1))*KroneckerDelta(0,gO1) + 2*g2*Conj(UP1(gI2,0))*
      KroneckerDelta(3,gO1)) + Conj(UM1(gI1,1))*(1.4142135623730951*g2*Conj(UP1(
      gI2,0))*KroneckerDelta(0,gO1) + Conj(UP1(gI2,1))*(-1.4142135623730951*LamSD*
      KroneckerDelta(2,gO1) + LamTD*KroneckerDelta(3,gO1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarCha2Cha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(KroneckerDelta(3,gO2)*(2*g2*UM2(gI1,0)
      *UP2(gI2,0) - Conj(LamTU)*UM2(gI1,1)*UP2(gI2,1)) - 1.4142135623730951*(Conj(
      LamTU)*KroneckerDelta(1,gO2)*UM2(gI1,1)*UP2(gI2,0) + (g2*KroneckerDelta(1,
      gO2)*UM2(gI1,0) + Conj(LamSU)*KroneckerDelta(2,gO2)*UM2(gI1,1))*UP2(gI2,1)))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarCha2Cha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(g2*Conj(UM2(gI2,0))*(
      1.4142135623730951*Conj(UP2(gI1,1))*KroneckerDelta(1,gO1) - 2*Conj(UP2(gI1,0
      ))*KroneckerDelta(3,gO1)) + Conj(UM2(gI2,1))*(1.4142135623730951*LamTU*Conj(
      UP2(gI1,0))*KroneckerDelta(1,gO1) + Conj(UP2(gI1,1))*(1.4142135623730951*
      LamSU*KroneckerDelta(2,gO1) + LamTU*KroneckerDelta(3,gO1))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjRhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*Mu*(2*Conj(LamSD)*(KroneckerDelta(2,gO2)*ZA(gI2,1) +
      KroneckerDelta(1,gO2)*ZA(gI2,2))*ZHR(gI1,0) + 1.4142135623730951*Conj(LamTD)
      *(KroneckerDelta(3,gO2)*ZA(gI2,1) + KroneckerDelta(1,gO2)*ZA(gI2,3))*ZHR(gI1
      ,0) + (-2*Conj(LamSU)*(KroneckerDelta(2,gO2)*ZA(gI2,0) + KroneckerDelta(0,
      gO2)*ZA(gI2,2)) + 1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZA(
      gI2,0) + KroneckerDelta(0,gO2)*ZA(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjRhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.25)*Mu*(2*Conj(LamSD)*(KroneckerDelta(2,
      gO2)*ZH(gI2,1) - KroneckerDelta(1,gO2)*ZH(gI2,2))*ZHR(gI1,0) +
      1.4142135623730951*Conj(LamTD)*(KroneckerDelta(3,gO2)*ZH(gI2,1) -
      KroneckerDelta(1,gO2)*ZH(gI2,3))*ZHR(gI1,0) + (Conj(LamSU)*(-2*
      KroneckerDelta(2,gO2)*ZH(gI2,0) + 2*KroneckerDelta(0,gO2)*ZH(gI2,2)) +
      1.4142135623730951*Conj(LamTU)*(KroneckerDelta(3,gO2)*ZH(gI2,0) -
      KroneckerDelta(0,gO2)*ZH(gI2,3)))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjRpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(1.4142135623730951*Conj(LamSU)*
      KroneckerDelta(2,gO2)*Mu*ZP(gI2,0)*ZRP(gI1,0) + Conj(LamTU)*Mu*(
      KroneckerDelta(3,gO2)*ZP(gI2,0) + 1.4142135623730951*KroneckerDelta(0,gO2)*
      ZP(gI2,3))*ZRP(gI1,0) + Conj(Mu)*(1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*ZP(gI2,1) - LamTD*KroneckerDelta(3,gO2)*ZP(gI2,1) + 1.4142135623730951*
      LamTD*KroneckerDelta(1,gO2)*ZP(gI2,2))*ZRP(gI1,1));

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

std::complex<double> CLASSNAME::CpUAhconjSvSv(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(3.872983346207417*g1*(MDBS - Conj(MDBS
      ))*KroneckerDelta(2,gO2) + 5*g2*(-MDWBT + Conj(MDWBT))*KroneckerDelta(3,gO2)
      )*KroneckerDelta(gI1,gI2);

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1154;
   std::complex<double> tmp_1155;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1156;
      std::complex<double> tmp_1157;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1157 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1156 += tmp_1157;
      tmp_1155 += (ZDL(gI1,j2)) * tmp_1156;
   }
   tmp_1154 += tmp_1155;
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2))
      * tmp_1154;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFdFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1158;
   std::complex<double> tmp_1159;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1160;
      std::complex<double> tmp_1161;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1161 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1160 += tmp_1161;
      tmp_1159 += (Conj(ZDL(gI2,j2))) * tmp_1160;
   }
   tmp_1158 += tmp_1159;
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_1158;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1162;
   std::complex<double> tmp_1163;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1164;
      std::complex<double> tmp_1165;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1165 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1164 += tmp_1165;
      tmp_1163 += (ZEL(gI1,j2)) * tmp_1164;
   }
   tmp_1162 += tmp_1163;
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(0,gO2))
      * tmp_1162;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFeFePL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1166;
   std::complex<double> tmp_1167;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1168;
      std::complex<double> tmp_1169;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1169 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1168 += tmp_1169;
      tmp_1167 += (Conj(ZEL(gI2,j2))) * tmp_1168;
   }
   tmp_1166 += tmp_1167;
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(0,gO1)
      ) * tmp_1166;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1170;
   std::complex<double> tmp_1171;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1172;
      std::complex<double> tmp_1173;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1173 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_1172 += tmp_1173;
      tmp_1171 += (ZUL(gI1,j2)) * tmp_1172;
   }
   tmp_1170 += tmp_1171;
   result += (std::complex<double>(0,0.7071067811865475)*KroneckerDelta(1,gO2))
      * tmp_1170;

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarFuFuPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1174;
   std::complex<double> tmp_1175;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1176;
      std::complex<double> tmp_1177;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1177 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1176 += tmp_1177;
      tmp_1175 += (Conj(ZUL(gI2,j2))) * tmp_1176;
   }
   tmp_1174 += tmp_1175;
   result += (std::complex<double>(0,-0.7071067811865475)*KroneckerDelta(1,gO1)
      ) * tmp_1174;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(4*AbsSqr(LamSU)*KroneckerDelta(2,gO1)*KroneckerDelta(2,
      gO2)*ZA(gI1,1)*ZA(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZA(gI1,1)*ZA(gI2,1) -
      1.4142135623730951*LamSU*Conj(LamTU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,
      gO1)*ZA(gI1,1)*ZA(gI2,1) - 1.4142135623730951*LamTU*Conj(LamSU)*
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

std::complex<double> CLASSNAME::CpUAhUAhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3
      ,gO1)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 1.4142135623730951*KroneckerDelta(0,gO2)
      *KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) + 4*AbsSqr(LamSU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*KroneckerDelta(3,
      gO1)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamSU*Conj(LamTU)*
      KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,1) +
      1.4142135623730951*LamTU*Conj(LamSU)*KroneckerDelta(2,gO1)*KroneckerDelta(3,
      gO2)*ZP(gI1,1)*ZP(gI2,1) + 1.4142135623730951*LamSU*Conj(LamTU)*
      KroneckerDelta(2,gO1)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) + 2*AbsSqr(
      LamTU)*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*ZP(gI1,1)*ZP(gI2,1) - 2*
      LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,2)*ZP(
      gI2,1) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZP(gI1,2)*ZP(gI2,1) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) + 2*
      LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(2,gO1)*ZP(gI1,3)*ZP(
      gI2,1) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*
      KroneckerDelta(3,gO1)*ZP(gI1,3)*ZP(gI2,1) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,0)*ZP(gI2,2) - 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*KroneckerDelta(
      2,gO1)*ZP(gI1,1)*ZP(gI2,2) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta
      (1,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,2) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) + 4*
      KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,2) +
      Conj(LamSD)*(-1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*KroneckerDelta(
      3,gO1)*ZP(gI1,0)*ZP(gI2,0) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2
      ,gO2)*ZP(gI1,0)*ZP(gI2,0) - LamTD*(1.4142135623730951*KroneckerDelta(3,gO2)*
      ZP(gI1,0)*ZP(gI2,0) + 2*KroneckerDelta(0,gO2)*(ZP(gI1,3)*ZP(gI2,0) - ZP(gI1,
      0)*ZP(gI2,2))))) + 1.4142135623730951*KroneckerDelta(0,gO2)*KroneckerDelta(3
      ,gO1)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) + 2*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2
      )*KroneckerDelta(2,gO1)*ZP(gI1,1)*ZP(gI2,3) + 1.4142135623730951*AbsSqr(
      LamTU)*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*ZP(gI1,1)*ZP(gI2,3) -
      1.4142135623730951*KroneckerDelta(1,gO2)*KroneckerDelta(3,gO1)*Sqr(g2)*ZP(
      gI1,1)*ZP(gI2,3) + 4*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(
      gI1,2)*ZP(gI2,3) + 4*KroneckerDelta(3,gO1)*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(
      gI1,3)*ZP(gI2,3) - Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1)*ZP(gI1,0)*ZP(gI2,0) + LamSD*KroneckerDelta(2,gO1)
      *(1.4142135623730951*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) + 2*
      KroneckerDelta(0,gO2)*(-(ZP(gI1,2)*ZP(gI2,0)) + ZP(gI1,0)*ZP(gI2,3))) +
      LamTD*KroneckerDelta(3,gO1)*(-2*KroneckerDelta(3,gO2)*ZP(gI1,0)*ZP(gI2,0) +
      1.4142135623730951*KroneckerDelta(0,gO2)*(ZP(gI1,2)*ZP(gI2,0) + ZP(gI1,3)*ZP
      (gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2,3)))))) + KroneckerDelta(0,gO1)*(-(
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*Sqr
      (g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*AbsSqr(LamTD) + Sqr(g2))*
      ZP(gI1,2)*ZP(gI2,2)) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)))) + 5*(
      -1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,0) + 2*
      LamTD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,0) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) +
      KroneckerDelta(1,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) -
      2*LamTD*Conj(LamSD)*KroneckerDelta(2,gO2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) -
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,0)*ZP(gI2,3) + Conj(
      LamTD)*(KroneckerDelta(2,gO2)*(-2*LamSD*ZP(gI1,2)*ZP(gI2,0) + 2*LamSD*ZP(gI1
      ,0)*ZP(gI2,3)) + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)*(ZP(gI1,2)*
      ZP(gI2,0) + ZP(gI1,3)*ZP(gI2,0) + ZP(gI1,0)*(ZP(gI2,2) + ZP(gI2,3)))))) +
      KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1
      ,0)*ZP(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP
      (gI1,2)*ZP(gI2,2) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3))) + 5*(
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,2)*ZP(gI2,1) - 2*
      LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI1,3)*ZP(gI2,1) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,3)*ZP(gI2,1) +
      KroneckerDelta(0,gO2)*Sqr(g2)*(ZP(gI1,1)*ZP(gI2,0) + ZP(gI1,0)*ZP(gI2,1)) +
      2*LamTU*Conj(LamSU)*KroneckerDelta(2,gO2)*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,2) +
      1.4142135623730951*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI1,1)*ZP(gI2,3) + Conj(
      LamTU)*(2*LamSU*KroneckerDelta(2,gO2)*(ZP(gI1,2)*ZP(gI2,1) - ZP(gI1,1)*ZP(
      gI2,3)) - 1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*(ZP(gI1,2)*ZP(gI2,1
      ) + ZP(gI1,3)*ZP(gI2,1) + ZP(gI1,1)*(ZP(gI2,2) + ZP(gI2,3)))))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(Conj(LamTD)*(1.4142135623730951*LamSD*KroneckerDelta(2,
      gO2)*KroneckerDelta(3,gO1) + (1.4142135623730951*LamSD*KroneckerDelta(2,gO1)
      + 2*LamTD*KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))*ZH(gI1,0)*ZH(gI2,0)
      + Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(2,gO2)*
      KroneckerDelta(3,gO1) + KroneckerDelta(2,gO1)*(4*LamSD*KroneckerDelta(2,gO2)
      + 1.4142135623730951*LamTD*KroneckerDelta(3,gO2)))*ZH(gI1,0)*ZH(gI2,0) + (-
      (Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*KroneckerDelta(
      3,gO1) + (1.4142135623730951*LamSU*KroneckerDelta(2,gO1) - 2*LamTU*
      KroneckerDelta(3,gO1))*KroneckerDelta(3,gO2))) - Conj(LamSU)*(
      1.4142135623730951*LamTU*KroneckerDelta(2,gO2)*KroneckerDelta(3,gO1) +
      KroneckerDelta(2,gO1)*(-4*LamSU*KroneckerDelta(2,gO2) + 1.4142135623730951*
      LamTU*KroneckerDelta(3,gO2))))*ZH(gI1,1)*ZH(gI2,1)) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(
      g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 5*(Conj(LamSD)*(1.4142135623730951*
      LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(4*LamSD*ZH(gI2,2) +
      1.4142135623730951*LamTD*ZH(gI2,3))) + Conj(LamTD)*(1.4142135623730951*LamSD
      *ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951*LamSD*ZH(gI2,2) + 2*
      LamTD*ZH(gI2,3))))) + KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1
      ) + 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(
      gI2,1) + 5*(Conj(LamTU)*(1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(
      gI1,3)*(1.4142135623730951*LamSU*ZH(gI2,2) - 2*LamTU*ZH(gI2,3))) + Conj(
      LamSU)*(1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(-4*LamSU*
      ZH(gI2,2) + 1.4142135623730951*LamTU*ZH(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(7.745966692414834*g1*MDBS*
      KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(gI2,0) + 14.142135623730951*MuD*Conj(
      LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(gI2,0) + 7.0710678118654755*LamTD*
      vT*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(gI2,0) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(
      gI2,0) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(0,gO2)*ZA(gI1,2)*ZA(
      gI2,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI1,2)*
      ZA(gI2,0) - 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(gI2,0) -
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(
      gI2,0) + 10*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,3)*ZA(gI2,0) +
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

std::complex<double> CLASSNAME::CpUAhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(KroneckerDelta(2,gO2)*(-10*LamSD*vd*
      Conj(LamTD)*ZP(gI1,2)*ZP(gI2,0) + 10*LamTD*vd*Conj(LamSD)*ZP(gI1,3)*ZP(gI2,0
      ) - 7.745966692414834*g1*MDBS*ZP(gI1,1)*ZP(gI2,1) + 14.142135623730951*MuU*
      Conj(LamSU)*ZP(gI1,1)*ZP(gI2,1) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZP
      (gI1,1)*ZP(gI2,1) - 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZP(gI1,1)*ZP(gI2
      ,1) + 7.745966692414834*g1*Conj(MDBS)*ZP(gI1,1)*ZP(gI2,1) -
      14.142135623730951*LamSU*Conj(MuU)*ZP(gI1,1)*ZP(gI2,1) - 10*LamSU*vu*Conj(
      LamTU)*ZP(gI1,2)*ZP(gI2,1) + 10*LamTU*vu*Conj(LamSU)*ZP(gI1,3)*ZP(gI2,1) +
      10*LamTU*vu*Conj(LamSU)*ZP(gI1,1)*ZP(gI2,2) - 10*LamSU*vu*Conj(LamTU)*ZP(gI1
      ,1)*ZP(gI2,3) + ZP(gI1,0)*((7.745966692414834*g1*MDBS + 7.0710678118654755*(
      2*MuD - LamTD*vT)*Conj(LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) -
      7.745966692414834*g1*Conj(MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZP(gI2
      ,0) + 10*vd*(LamTD*Conj(LamSD)*ZP(gI2,2) - LamSD*Conj(LamTD)*ZP(gI2,3)))) +
      5*(KroneckerDelta(0,gO2)*(vu*Sqr(g2)*ZP(gI1,1)*ZP(gI2,0) + ((
      2.8284271247461903*MuD + 2*LamSD*vS - 1.4142135623730951*LamTD*vT)*Conj(
      LamTD) + 1.4142135623730951*g2*(g2*vT + 2*Conj(MDWBT)))*ZP(gI1,2)*ZP(gI2,0)
      + 2.8284271247461903*g2*MDWBT*ZP(gI1,3)*ZP(gI2,0) + 1.4142135623730951*vT*
      AbsSqr(LamTD)*ZP(gI1,3)*ZP(gI2,0) + 2*LamTD*vS*Conj(LamSD)*ZP(gI1,3)*ZP(gI2,
      0) + 2.8284271247461903*LamTD*Conj(MuD)*ZP(gI1,3)*ZP(gI2,0) -
      1.4142135623730951*vT*Sqr(g2)*ZP(gI1,3)*ZP(gI2,0) - vu*Sqr(g2)*ZP(gI1,0)*ZP(
      gI2,1) - 2.8284271247461903*g2*MDWBT*ZP(gI1,0)*ZP(gI2,2) +
      1.4142135623730951*vT*AbsSqr(LamTD)*ZP(gI1,0)*ZP(gI2,2) - 2*LamTD*vS*Conj(
      LamSD)*ZP(gI1,0)*ZP(gI2,2) - 2.8284271247461903*LamTD*Conj(MuD)*ZP(gI1,0)*ZP
      (gI2,2) - 1.4142135623730951*vT*Sqr(g2)*ZP(gI1,0)*ZP(gI2,2) -
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

std::complex<double> CLASSNAME::CpUAhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-5*(vd*Conj(LamSD)*(1.4142135623730951*LamTD*KroneckerDelta(3
      ,gO2)*ZA(gI2,2) + KroneckerDelta(2,gO2)*(4*LamSD*ZA(gI2,2) +
      1.4142135623730951*LamTD*ZA(gI2,3)))*ZH(gI1,0) + vd*Conj(LamTD)*(
      1.4142135623730951*LamSD*KroneckerDelta(2,gO2)*ZA(gI2,3) + KroneckerDelta(3,
      gO2)*(1.4142135623730951*LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3)))*ZH(gI1,0) +
      vu*(-(Conj(LamTU)*(1.4142135623730951*LamSU*KroneckerDelta(2,gO2)*ZA(gI2,3)
      + KroneckerDelta(3,gO2)*(1.4142135623730951*LamSU*ZA(gI2,2) - 2*LamTU*ZA(gI2
      ,3)))) - Conj(LamSU)*(1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*ZA(gI2,
      2) + KroneckerDelta(2,gO2)*(-4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA
      (gI2,3))))*ZH(gI1,1)) - KroneckerDelta(0,gO2)*ZA(gI2,0)*(vd*(3*Sqr(g1) + 5*
      Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1) -
      7.745966692414834*g1*MDBS*ZH(gI1,2) + 20*vS*AbsSqr(LamSD)*ZH(gI1,2) +
      14.142135623730951*MuD*Conj(LamSD)*ZH(gI1,2) + 7.0710678118654755*LamTD*vT*
      Conj(LamSD)*ZH(gI1,2) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gI1,2) -
      7.745966692414834*g1*Conj(MDBS)*ZH(gI1,2) + 14.142135623730951*LamSD*Conj(
      MuD)*ZH(gI1,2) + 10*g2*MDWBT*ZH(gI1,3) + 10*vT*AbsSqr(LamTD)*ZH(gI1,3) +
      7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gI1,3) + 10*MuD*Conj(LamTD)*ZH(
      gI1,3) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gI1,3) + 10*g2*Conj(
      MDWBT)*ZH(gI1,3) + 10*LamTD*Conj(MuD)*ZH(gI1,3)) + KroneckerDelta(1,gO2)*ZA(
      gI2,1)*(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,0) - vu*(3*Sqr(g1) + 5*Sqr(g2))*ZH
      (gI1,1) - 7.745966692414834*g1*MDBS*ZH(gI1,2) - 20*vS*AbsSqr(LamSU)*ZH(gI1,2
      ) - 14.142135623730951*MuU*Conj(LamSU)*ZH(gI1,2) + 7.0710678118654755*LamTU*
      vT*Conj(LamSU)*ZH(gI1,2) + 7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gI1,2)
      - 7.745966692414834*g1*Conj(MDBS)*ZH(gI1,2) - 14.142135623730951*LamSU*Conj
      (MuU)*ZH(gI1,2) + 10*g2*MDWBT*ZH(gI1,3) - 10*vT*AbsSqr(LamTU)*ZH(gI1,3) +
      7.0710678118654755*LamTU*vS*Conj(LamSU)*ZH(gI1,3) + 10*MuU*Conj(LamTU)*ZH(
      gI1,3) + 7.0710678118654755*LamSU*vS*Conj(LamTU)*ZH(gI1,3) + 10*g2*Conj(
      MDWBT)*ZH(gI1,3) + 10*LamTU*Conj(MuU)*ZH(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(-5*KroneckerDelta(3,gO2)*(
      1.4142135623730951*LamTD*vd*Conj(LamSD)*ZH(gI1,2)*ZH(gI2,0) -
      1.4142135623730951*LamSD*vd*Conj(LamTD)*ZH(gI1,2)*ZH(gI2,0) - 2*g2*MDWBT*ZH(
      gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*vS*Conj(LamSU)*ZH(gI1,1)*ZH(gI2,
      1) + 2*MuU*Conj(LamTU)*ZH(gI1,1)*ZH(gI2,1) + 1.4142135623730951*LamSU*vS*
      Conj(LamTU)*ZH(gI1,1)*ZH(gI2,1) + 2*g2*Conj(MDWBT)*ZH(gI1,1)*ZH(gI2,1) - 2*
      LamTU*Conj(MuU)*ZH(gI1,1)*ZH(gI2,1) - 1.4142135623730951*LamTU*vu*Conj(LamSU
      )*ZH(gI1,2)*ZH(gI2,1) + 1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gI1,2)*ZH
      (gI2,1) - 1.4142135623730951*LamTU*vu*Conj(LamSU)*ZH(gI1,1)*ZH(gI2,2) +
      1.4142135623730951*LamSU*vu*Conj(LamTU)*ZH(gI1,1)*ZH(gI2,2) + ZH(gI1,0)*((2*
      g2*MDWBT + 1.4142135623730951*LamTD*vS*Conj(LamSD) - 2*MuD*Conj(LamTD) -
      1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) + 2*LamTD*Conj(
      MuD))*ZH(gI2,0) + 1.4142135623730951*vd*(LamTD*Conj(LamSD) - LamSD*Conj(
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

std::complex<double> CLASSNAME::CpUAhbarChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*Conj(LamTD)*KroneckerDelta(0,gO2)*
      ZN1(gI1,2)*ZN2(gI2,1) + 5*Conj(LamTU)*KroneckerDelta(1,gO2)*ZN1(gI1,3)*ZN2(
      gI2,1) - 3.872983346207417*g1*KroneckerDelta(0,gO2)*ZN1(gI1,0)*ZN2(gI2,2) +
      5*g2*KroneckerDelta(0,gO2)*ZN1(gI1,1)*ZN2(gI2,2) + 5*Conj(LamTD)*
      KroneckerDelta(3,gO2)*ZN1(gI1,2)*ZN2(gI2,2) + 7.0710678118654755*Conj(LamSD)
      *ZN1(gI1,2)*(KroneckerDelta(0,gO2)*ZN2(gI2,0) + KroneckerDelta(2,gO2)*ZN2(
      gI2,2)) + 3.872983346207417*g1*KroneckerDelta(1,gO2)*ZN1(gI1,0)*ZN2(gI2,3) -
      5*g2*KroneckerDelta(1,gO2)*ZN1(gI1,1)*ZN2(gI2,3) + 5*Conj(LamTU)*
      KroneckerDelta(3,gO2)*ZN1(gI1,3)*ZN2(gI2,3) - 7.0710678118654755*Conj(LamSU)
      *ZN1(gI1,3)*(KroneckerDelta(1,gO2)*ZN2(gI2,0) + KroneckerDelta(2,gO2)*ZN2(
      gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpUAhbarChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(3.872983346207417*g1*Conj(ZN1(gI2,0))*
      (-(Conj(ZN2(gI1,2))*KroneckerDelta(0,gO1)) + Conj(ZN2(gI1,3))*KroneckerDelta
      (1,gO1)) + 5*Conj(ZN1(gI2,2))*(1.4142135623730951*LamSD*Conj(ZN2(gI1,0))*
      KroneckerDelta(0,gO1) + LamTD*Conj(ZN2(gI1,1))*KroneckerDelta(0,gO1) + Conj(
      ZN2(gI1,2))*(1.4142135623730951*LamSD*KroneckerDelta(2,gO1) + LamTD*
      KroneckerDelta(3,gO1))) + 5*(g2*Conj(ZN1(gI2,1))*(Conj(ZN2(gI1,2))*
      KroneckerDelta(0,gO1) - Conj(ZN2(gI1,3))*KroneckerDelta(1,gO1)) + Conj(ZN1(
      gI2,3))*(-1.4142135623730951*LamSU*Conj(ZN2(gI1,0))*KroneckerDelta(1,gO1) +
      LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(1,gO1) + Conj(ZN2(gI1,3))*(
      -1.4142135623730951*LamSU*KroneckerDelta(2,gO1) + LamTU*KroneckerDelta(3,gO1
      )))));

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1178;
   std::complex<double> tmp_1179;
   std::complex<double> tmp_1180;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1180 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1179 += tmp_1180;
   tmp_1178 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1179;
   std::complex<double> tmp_1181;
   std::complex<double> tmp_1182;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1182 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1181 += tmp_1182;
   tmp_1178 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1181;
   std::complex<double> tmp_1183;
   std::complex<double> tmp_1184;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1184 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1183 += tmp_1184;
   tmp_1178 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1183;
   std::complex<double> tmp_1185;
   std::complex<double> tmp_1186;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1186 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1185 += tmp_1186;
   tmp_1178 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1185;
   std::complex<double> tmp_1187;
   std::complex<double> tmp_1188;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1188 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1187 += tmp_1188;
   tmp_1178 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1187;
   std::complex<double> tmp_1189;
   std::complex<double> tmp_1190;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1190 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1189 += tmp_1190;
   tmp_1178 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1189;
   std::complex<double> tmp_1191;
   std::complex<double> tmp_1192;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1193;
      std::complex<double> tmp_1194;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1195;
         std::complex<double> tmp_1196;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1196 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1195 += tmp_1196;
         tmp_1194 += (ZD(gI1,3 + j2)) * tmp_1195;
      }
      tmp_1193 += tmp_1194;
      tmp_1192 += (Conj(ZD(gI2,3 + j3))) * tmp_1193;
   }
   tmp_1191 += tmp_1192;
   tmp_1178 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1191;
   std::complex<double> tmp_1197;
   std::complex<double> tmp_1198;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1199;
      std::complex<double> tmp_1200;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1201;
         std::complex<double> tmp_1202;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1202 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1201 += tmp_1202;
         tmp_1200 += (Conj(ZD(gI2,j2))) * tmp_1201;
      }
      tmp_1199 += tmp_1200;
      tmp_1198 += (ZD(gI1,j3)) * tmp_1199;
   }
   tmp_1197 += tmp_1198;
   tmp_1178 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1197;
   result += (std::complex<double>(0,-1)) * tmp_1178;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1203;
   std::complex<double> tmp_1204;
   std::complex<double> tmp_1205;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1205 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1204 += tmp_1205;
   tmp_1203 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1204;
   std::complex<double> tmp_1206;
   std::complex<double> tmp_1207;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1207 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1206 += tmp_1207;
   tmp_1203 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1206;
   std::complex<double> tmp_1208;
   std::complex<double> tmp_1209;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1209 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1208 += tmp_1209;
   tmp_1203 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1208;
   std::complex<double> tmp_1210;
   std::complex<double> tmp_1211;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1211 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1210 += tmp_1211;
   tmp_1203 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1210;
   std::complex<double> tmp_1212;
   std::complex<double> tmp_1213;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1213 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1212 += tmp_1213;
   tmp_1203 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1212;
   std::complex<double> tmp_1214;
   std::complex<double> tmp_1215;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1215 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1214 += tmp_1215;
   tmp_1203 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1214;
   std::complex<double> tmp_1216;
   std::complex<double> tmp_1217;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1218;
      std::complex<double> tmp_1219;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1220;
         std::complex<double> tmp_1221;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1221 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1220 += tmp_1221;
         tmp_1219 += (ZE(gI1,3 + j2)) * tmp_1220;
      }
      tmp_1218 += tmp_1219;
      tmp_1217 += (Conj(ZE(gI2,3 + j3))) * tmp_1218;
   }
   tmp_1216 += tmp_1217;
   tmp_1203 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1216;
   std::complex<double> tmp_1222;
   std::complex<double> tmp_1223;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1224;
      std::complex<double> tmp_1225;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1226;
         std::complex<double> tmp_1227;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1227 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1226 += tmp_1227;
         tmp_1225 += (Conj(ZE(gI2,j2))) * tmp_1226;
      }
      tmp_1224 += tmp_1225;
      tmp_1223 += (ZE(gI1,j3)) * tmp_1224;
   }
   tmp_1222 += tmp_1223;
   tmp_1203 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1222;
   result += (std::complex<double>(0,-1)) * tmp_1203;

   return result;
}

std::complex<double> CLASSNAME::CpUAhUAhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1228;
   std::complex<double> tmp_1229;
   std::complex<double> tmp_1230;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1230 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1229 += tmp_1230;
   tmp_1228 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1229;
   std::complex<double> tmp_1231;
   std::complex<double> tmp_1232;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1232 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1231 += tmp_1232;
   tmp_1228 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1231;
   std::complex<double> tmp_1233;
   std::complex<double> tmp_1234;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1234 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1233 += tmp_1234;
   tmp_1228 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1233;
   std::complex<double> tmp_1235;
   std::complex<double> tmp_1236;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1236 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1235 += tmp_1236;
   tmp_1228 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1235;
   std::complex<double> tmp_1237;
   std::complex<double> tmp_1238;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1238 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1237 += tmp_1238;
   tmp_1228 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1237;
   std::complex<double> tmp_1239;
   std::complex<double> tmp_1240;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1240 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1239 += tmp_1240;
   tmp_1228 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1239;
   std::complex<double> tmp_1241;
   std::complex<double> tmp_1242;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1243;
      std::complex<double> tmp_1244;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1245;
         std::complex<double> tmp_1246;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1246 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1245 += tmp_1246;
         tmp_1244 += (ZU(gI1,3 + j2)) * tmp_1245;
      }
      tmp_1243 += tmp_1244;
      tmp_1242 += (Conj(ZU(gI2,3 + j3))) * tmp_1243;
   }
   tmp_1241 += tmp_1242;
   tmp_1228 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1241;
   std::complex<double> tmp_1247;
   std::complex<double> tmp_1248;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1249;
      std::complex<double> tmp_1250;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1251;
         std::complex<double> tmp_1252;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1252 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1251 += tmp_1252;
         tmp_1250 += (Conj(ZU(gI2,j2))) * tmp_1251;
      }
      tmp_1249 += tmp_1250;
      tmp_1248 += (ZU(gI1,j3)) * tmp_1249;
   }
   tmp_1247 += tmp_1248;
   tmp_1228 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1247;
   result += (std::complex<double>(0,-1)) * tmp_1228;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1253;
   std::complex<double> tmp_1254;
   std::complex<double> tmp_1255;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1255 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1254 += tmp_1255;
   tmp_1253 += (0.12909944487358055*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1254;
   std::complex<double> tmp_1256;
   std::complex<double> tmp_1257;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1257 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1256 += tmp_1257;
   tmp_1253 += (-0.12909944487358055*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1256;
   std::complex<double> tmp_1258;
   std::complex<double> tmp_1259;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1259 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1258 += tmp_1259;
   tmp_1253 += (-0.5*g2*MDWBT*KroneckerDelta(3,gO2)) * tmp_1258;
   std::complex<double> tmp_1260;
   std::complex<double> tmp_1261;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1261 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1260 += tmp_1261;
   tmp_1253 += (0.5*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)) * tmp_1260;
   std::complex<double> tmp_1262;
   std::complex<double> tmp_1263;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1263 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1262 += tmp_1263;
   tmp_1253 += (0.2581988897471611*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1262;
   std::complex<double> tmp_1264;
   std::complex<double> tmp_1265;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1265 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1264 += tmp_1265;
   tmp_1253 += (-0.2581988897471611*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1264;
   std::complex<double> tmp_1266;
   std::complex<double> tmp_1267;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1268;
      std::complex<double> tmp_1269;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1269 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1268 += tmp_1269;
      tmp_1267 += (Conj(ZD(gI2,j2))) * tmp_1268;
   }
   tmp_1266 += tmp_1267;
   tmp_1253 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1266;
   std::complex<double> tmp_1270;
   std::complex<double> tmp_1271;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1272;
      std::complex<double> tmp_1273;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1273 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1272 += tmp_1273;
      tmp_1271 += (ZD(gI1,j2)) * tmp_1272;
   }
   tmp_1270 += tmp_1271;
   tmp_1253 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1270;
   result += (std::complex<double>(0,-1)) * tmp_1253;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1274;
   std::complex<double> tmp_1275;
   std::complex<double> tmp_1276;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1276 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1275 += tmp_1276;
   tmp_1274 += (-0.3872983346207417*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1275;
   std::complex<double> tmp_1277;
   std::complex<double> tmp_1278;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1278 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1277 += tmp_1278;
   tmp_1274 += (0.3872983346207417*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1277;
   std::complex<double> tmp_1279;
   std::complex<double> tmp_1280;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1280 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1279 += tmp_1280;
   tmp_1274 += (-0.5*g2*MDWBT*KroneckerDelta(3,gO2)) * tmp_1279;
   std::complex<double> tmp_1281;
   std::complex<double> tmp_1282;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1282 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1281 += tmp_1282;
   tmp_1274 += (0.5*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)) * tmp_1281;
   std::complex<double> tmp_1283;
   std::complex<double> tmp_1284;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1284 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1283 += tmp_1284;
   tmp_1274 += (0.7745966692414834*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1283;
   std::complex<double> tmp_1285;
   std::complex<double> tmp_1286;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1286 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1285 += tmp_1286;
   tmp_1274 += (-0.7745966692414834*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1285;
   std::complex<double> tmp_1287;
   std::complex<double> tmp_1288;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1289;
      std::complex<double> tmp_1290;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1290 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1289 += tmp_1290;
      tmp_1288 += (Conj(ZE(gI2,j2))) * tmp_1289;
   }
   tmp_1287 += tmp_1288;
   tmp_1274 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(1,gO2)) * tmp_1287;
   std::complex<double> tmp_1291;
   std::complex<double> tmp_1292;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1293;
      std::complex<double> tmp_1294;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1294 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1293 += tmp_1294;
      tmp_1292 += (ZE(gI1,j2)) * tmp_1293;
   }
   tmp_1291 += tmp_1292;
   tmp_1274 += (-0.7071067811865475*KroneckerDelta(1,gO2)*Mu) * tmp_1291;
   result += (std::complex<double>(0,-1)) * tmp_1274;

   return result;
}

std::complex<double> CLASSNAME::CpUAhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1295;
   std::complex<double> tmp_1296;
   std::complex<double> tmp_1297;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1297 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1296 += tmp_1297;
   tmp_1295 += (0.12909944487358055*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1296;
   std::complex<double> tmp_1298;
   std::complex<double> tmp_1299;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1299 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1298 += tmp_1299;
   tmp_1295 += (-0.12909944487358055*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1298;
   std::complex<double> tmp_1300;
   std::complex<double> tmp_1301;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1301 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1300 += tmp_1301;
   tmp_1295 += (0.5*g2*MDWBT*KroneckerDelta(3,gO2)) * tmp_1300;
   std::complex<double> tmp_1302;
   std::complex<double> tmp_1303;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1303 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1302 += tmp_1303;
   tmp_1295 += (-0.5*g2*Conj(MDWBT)*KroneckerDelta(3,gO2)) * tmp_1302;
   std::complex<double> tmp_1304;
   std::complex<double> tmp_1305;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1305 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1304 += tmp_1305;
   tmp_1295 += (-0.5163977794943222*g1*MDBS*KroneckerDelta(2,gO2)) * tmp_1304;
   std::complex<double> tmp_1306;
   std::complex<double> tmp_1307;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1307 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1306 += tmp_1307;
   tmp_1295 += (0.5163977794943222*g1*Conj(MDBS)*KroneckerDelta(2,gO2)) *
      tmp_1306;
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
   tmp_1295 += (0.7071067811865475*Conj(Mu)*KroneckerDelta(0,gO2)) * tmp_1308;
   std::complex<double> tmp_1312;
   std::complex<double> tmp_1313;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1314;
      std::complex<double> tmp_1315;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1315 += Conj(Yu(j1,j2))*Conj(ZU(gI2,3 + j1));
      }
      tmp_1314 += tmp_1315;
      tmp_1313 += (ZU(gI1,j2)) * tmp_1314;
   }
   tmp_1312 += tmp_1313;
   tmp_1295 += (-0.7071067811865475*KroneckerDelta(0,gO2)*Mu) * tmp_1312;
   result += (std::complex<double>(0,-1)) * tmp_1295;

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
      KroneckerDelta(1,gO2)*ZP(gI2,1) + 1.4142135623730951*KroneckerDelta(3,gO2)*(
      ZP(gI2,2) - ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(3*Sqr(g1) + 5*Sqr(g2))*(KroneckerDelta(1,gO1)*(KroneckerDelta
      (0,gO2)*ZHR(gI1,1)*ZHR(gI2,0) + KroneckerDelta(1,gO2)*(ZHR(gI1,0)*ZHR(gI2,0)
      - 2*ZHR(gI1,1)*ZHR(gI2,1))) + KroneckerDelta(0,gO1)*(KroneckerDelta(1,gO2)*
      ZHR(gI1,0)*ZHR(gI2,1) + KroneckerDelta(0,gO2)*(-2*ZHR(gI1,0)*ZHR(gI2,0) +
      ZHR(gI1,1)*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(Conj(ZRP(gI2,0))*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)
      *(3*Sqr(g1) - 5*Sqr(g2)) - KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*
      Sqr(g1) + 5*Sqr(g2)))*ZRP(gI1,0) - Conj(ZRP(gI2,1))*(KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*(-3*Sqr(g1) + 5*Sqr(g2)) + KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2)))*ZRP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjRpmVWm(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*KroneckerDelta(0,gO2)*ZRP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjURhCha2Cha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(LamTD)*KroneckerDelta(0,gO2)*UM1(gI2,1)*UP2(gI1,0) - Conj(
      LamTU)*KroneckerDelta(1,gO2)*UM1(gI2,0)*UP2(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjURhCha2Cha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*(Conj(UM2(gI1,0))*Conj(UP1(gI2,1))*KroneckerDelta(0,gO1) +
      Conj(UM2(gI1,1))*Conj(UP1(gI2,0))*KroneckerDelta(1,gO1)));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjRpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*(1.4142135623730951*KroneckerDelta(1,gO2)*(2*LamSD*vu*Conj(
      LamSU)*ZP(gI2,0) + LamTD*Conj(LamTU)*(vu*ZP(gI2,0) + 2*vd*ZP(gI2,1))) +
      KroneckerDelta(0,gO2)*(1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(
      LamTD) + Sqr(g2))*ZP(gI2,0) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI2,1) + 4*g2
      *MDWBT*ZP(gI2,2) - 2*vT*AbsSqr(LamTD)*ZP(gI2,2) - 2.8284271247461903*LamTD*
      vS*Conj(LamSD)*ZP(gI2,2) - 4*LamTD*Conj(MuD)*ZP(gI2,2) + 2*vT*Sqr(g2)*ZP(gI2
      ,2) + 2*vT*AbsSqr(LamTD)*ZP(gI2,3) - 4*MuD*Conj(LamTD)*ZP(gI2,3) -
      2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gI2,3) + 4*g2*Conj(MDWBT)*ZP(gI2,
      3) - 2*vT*Sqr(g2)*ZP(gI2,3)))*ZRP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjURhRhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(-7.745966692414834*g1*MDBS*
      KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR(gI1,0) + 14.142135623730951*MuD*Conj(
      LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR(gI1,0) + 7.0710678118654755*LamTD
      *vT*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR(gI1,0) -
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR(
      gI1,0) + 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZHR
      (gI1,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI2,2)
      *ZHR(gI1,0) + 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(gI1,0) -
      7.0710678118654755*LamTD*vS*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(
      gI1,0) + 10*MuD*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(gI1,0) +
      7.0710678118654755*LamSD*vS*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(
      gI1,0) - 10*g2*Conj(MDWBT)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(gI1,0) - 10*
      LamTD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZHR(gI1,0) - 10*LamSU*vu*
      Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZHR(gI1,1) + 5*LamTU*vu*Conj(
      LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,0)*ZHR(gI1,1) + 10*LamSU*vd*Conj(LamSD)*
      KroneckerDelta(0,gO2)*ZA(gI2,1)*ZHR(gI1,1) - 5*LamTU*vd*Conj(LamTD)*
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

std::complex<double> CLASSNAME::CpconjURhRhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO2)*(-(((5*(2.8284271247461903*MuD + 4*
      LamSD*vS + 1.4142135623730951*LamTD*vT)*Conj(LamSD) + 7.0710678118654755*
      LamSD*vT*Conj(LamTD) + 2*(3.872983346207417*g1*MDBS + 3.872983346207417*g1*
      Conj(MDBS) + 7.0710678118654755*LamSD*Conj(MuD)))*ZH(gI2,2) + 5*(
      1.4142135623730951*LamTD*vS*Conj(LamSD) + (2*MuD + 1.4142135623730951*LamSD*
      vS + 2*LamTD*vT)*Conj(LamTD) - 2*(g2*MDWBT + g2*Conj(MDWBT) - LamTD*Conj(MuD
      )))*ZH(gI2,3))*ZHR(gI1,0)) + ZH(gI2,0)*(vd*(-20*AbsSqr(LamSD) - 10*AbsSqr(
      LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,0) + 5*vu*(2*LamSU*Conj(LamSD) -
      LamTU*Conj(LamTD))*ZHR(gI1,1)) - ZH(gI2,1)*(vu*(3*Sqr(g1) + 5*Sqr(g2))*ZHR(
      gI1,0) + 5*vd*(-2*LamSU*Conj(LamSD) + LamTU*Conj(LamTD))*ZHR(gI1,1))) +
      KroneckerDelta(1,gO2)*((-(vd*(3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,0)) + vu*(3*Sqr(
      g1) + 5*Sqr(g2))*ZH(gI2,1) + 2*((3.872983346207417*g1*MDBS +
      3.872983346207417*g1*Conj(MDBS) - 7.0710678118654755*LamSU*Conj(MuU))*ZH(gI2
      ,2) - 5*(g2*MDWBT + g2*Conj(MDWBT) - LamTU*Conj(MuU))*ZH(gI2,3)))*ZHR(gI1,1)
      + 5*Conj(LamSU)*(2*LamSD*vu*ZH(gI2,0)*ZHR(gI1,0) + ((-2.8284271247461903*
      MuU - 4*LamSU*vS + 1.4142135623730951*LamTU*vT)*ZH(gI2,2) +
      1.4142135623730951*LamTU*vS*ZH(gI2,3))*ZHR(gI1,1) + 2*ZH(gI2,1)*(LamSD*vd*
      ZHR(gI1,0) - 2*LamSU*vu*ZHR(gI1,1))) - 5*Conj(LamTU)*(LamTD*vu*ZH(gI2,0)*ZHR
      (gI1,0) - (1.4142135623730951*LamSU*vT*ZH(gI2,2) + (2*MuU +
      1.4142135623730951*LamSU*vS - 2*LamTU*vT)*ZH(gI2,3))*ZHR(gI1,1) + ZH(gI2,1)*
      (LamTD*vd*ZHR(gI1,0) + 2*LamTU*vu*ZHR(gI1,1)))));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) + 5*Sqr(g2)
      );

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(5*(2*LamSD*Conj(LamSU) - LamTD*Conj(
      LamTU))*KroneckerDelta(1,gO2)*(ZA(gI1,1)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,1)) +
      KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5
      *Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1)
      - 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(
      4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))) + Conj(LamTD)*(
      1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSD*ZA(gI2,2) + 2*LamTD*ZA(gI2,3)))))) + KroneckerDelta(1,gO1)*(5*(2*
      LamSU*Conj(LamSD) - LamTU*Conj(LamTD))*KroneckerDelta(0,gO2)*(ZA(gI1,1)*ZA(
      gI2,0) + ZA(gI1,0)*ZA(gI2,1)) + KroneckerDelta(1,gO2)*(-((3*Sqr(g1) + 5*Sqr(
      g2))*ZA(gI1,0)*ZA(gI2,0)) + (-20*AbsSqr(LamSU) - 10*AbsSqr(LamTU) + 3*Sqr(g1
      ) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTU)*(1.4142135623730951*
      LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSU*ZA(gI2,2) -
      2*LamTU*ZA(gI2,3))) + Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2
      ,2) + ZA(gI1,2)*(-4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3)))))
      ));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1) - 5
      *Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2))*
      ZP(gI1,1)*ZP(gI2,1) + 10*(-((-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,2)*ZP(gI2,2)
      ) + Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,
      gO2)*((-20*AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (-3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) -
      (-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,3)*ZP(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*(5*(2*LamSD*Conj(LamSU) - LamTD*Conj(
      LamTU))*KroneckerDelta(1,gO2)*(ZH(gI1,1)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,1)) +
      KroneckerDelta(0,gO2)*((-20*AbsSqr(LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5
      *Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1)
      - 5*(Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*(
      4*LamSD*ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3))) + Conj(LamTD)*(
      1.4142135623730951*LamSD*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951
      *LamSD*ZH(gI2,2) + 2*LamTD*ZH(gI2,3)))))) + KroneckerDelta(1,gO1)*(5*(2*
      LamSU*Conj(LamSD) - LamTU*Conj(LamTD))*KroneckerDelta(0,gO2)*(ZH(gI1,1)*ZH(
      gI2,0) + ZH(gI1,0)*ZH(gI2,1)) + KroneckerDelta(1,gO2)*(-((3*Sqr(g1) + 5*Sqr(
      g2))*ZH(gI1,0)*ZH(gI2,0)) + (-20*AbsSqr(LamSU) - 10*AbsSqr(LamTU) + 3*Sqr(g1
      ) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 5*(Conj(LamTU)*(1.4142135623730951*
      LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951*LamSU*ZH(gI2,2) -
      2*LamTU*ZH(gI2,3))) + Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2
      ,2) + ZH(gI1,2)*(-4*LamSU*ZH(gI2,2) + 1.4142135623730951*LamTU*ZH(gI2,3)))))
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjHpmRpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*Conj(ZRP(gI2,0))*(2.8284271247461903*LamSU*vd*Conj(LamSD)*
      KroneckerDelta(0,gO2)*ZP(gI1,1) + 1.4142135623730951*LamTU*Conj(LamTD)*
      KroneckerDelta(0,gO2)*(2*vu*ZP(gI1,0) + vd*ZP(gI1,1)) + KroneckerDelta(1,gO2
      )*(1.4142135623730951*vd*Sqr(g2)*ZP(gI1,0) + 1.4142135623730951*vu*(-2*
      AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,1) + 2*((-((2*MuU +
      1.4142135623730951*LamSU*vS + LamTU*vT)*Conj(LamTU)) + g2*(g2*vT + 2*Conj(
      MDWBT)))*ZP(gI1,2) - (-2*g2*MDWBT - vT*AbsSqr(LamTU) + 1.4142135623730951*
      LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(MuU) + vT*Sqr(g2))*ZP(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhAhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*Mu*(2*Conj(LamSU)*KroneckerDelta(1,gO2)*(ZA(gI1,2)*ZA(gI2,0) +
      ZA(gI1,0)*ZA(gI2,2)) - 1.4142135623730951*Conj(LamTU)*KroneckerDelta(1,gO2)
      *(ZA(gI1,3)*ZA(gI2,0) + ZA(gI1,0)*ZA(gI2,3)) - KroneckerDelta(0,gO2)*(2*Conj
      (LamSD)*(ZA(gI1,2)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,2)) + 1.4142135623730951*
      Conj(LamTD)*(ZA(gI1,3)*ZA(gI2,1) + ZA(gI1,1)*ZA(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjHpmHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(LamTU)*KroneckerDelta(1,gO2)*Mu*ZP(gI1,2)*ZP(gI2,0)) + Conj(
      LamTD)*KroneckerDelta(0,gO2)*Mu*ZP(gI1,1)*ZP(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpconjURhhhAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.25)*Mu*(2*Conj(LamSU)*KroneckerDelta(1,
      gO2)*(ZA(gI2,2)*ZH(gI1,0) - ZA(gI2,0)*ZH(gI1,2)) + 1.4142135623730951*Conj(
      LamTU)*KroneckerDelta(1,gO2)*(-(ZA(gI2,3)*ZH(gI1,0)) + ZA(gI2,0)*ZH(gI1,3))
      + KroneckerDelta(0,gO2)*(Conj(LamSD)*(-2*ZA(gI2,2)*ZH(gI1,1) + 2*ZA(gI2,1)*
      ZH(gI1,2)) + 1.4142135623730951*Conj(LamTD)*(-(ZA(gI2,3)*ZH(gI1,1)) + ZA(gI2
      ,1)*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhhhhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.25*Mu*(2*Conj(LamSU)*KroneckerDelta(1,gO2)*(ZH(gI1,2)*ZH(gI2,0) +
      ZH(gI1,0)*ZH(gI2,2)) - 1.4142135623730951*Conj(LamTU)*KroneckerDelta(1,gO2)
      *(ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,3)) - KroneckerDelta(0,gO2)*(2*Conj
      (LamSD)*(ZH(gI1,2)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,2)) + 1.4142135623730951*
      Conj(LamTD)*(ZH(gI1,3)*ZH(gI2,1) + ZH(gI1,1)*ZH(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhChiChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(LamSD)*KroneckerDelta(0,gO2)*(ZN2(gI1,2)*ZN2(gI2,0) + ZN2(gI1,
      0)*ZN2(gI2,2)) + 0.5*(-2*Conj(LamSU)*KroneckerDelta(1,gO2)*(ZN2(gI1,3)*ZN2(
      gI2,0) + ZN2(gI1,0)*ZN2(gI2,3)) + 1.4142135623730951*(Conj(LamTD)*
      KroneckerDelta(0,gO2)*(ZN2(gI1,2)*ZN2(gI2,1) + ZN2(gI1,1)*ZN2(gI2,2)) + Conj
      (LamTU)*KroneckerDelta(1,gO2)*(ZN2(gI1,3)*ZN2(gI2,1) + ZN2(gI1,1)*ZN2(gI2,3)
      )));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhChiChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1414213562373095*(Conj(ZN1(gI1,2))*(-3.872983346207417*g1*Conj(
      ZN1(gI2,0)) + 5*g2*Conj(ZN1(gI2,1)))*KroneckerDelta(0,gO1) + 5*g2*Conj(ZN1(
      gI1,1))*Conj(ZN1(gI2,2))*KroneckerDelta(0,gO1) + 3.872983346207417*g1*Conj(
      ZN1(gI1,3))*Conj(ZN1(gI2,0))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN1(gI1,3))*
      Conj(ZN1(gI2,1))*KroneckerDelta(1,gO1) - 5*g2*Conj(ZN1(gI1,1))*Conj(ZN1(gI2,
      3))*KroneckerDelta(1,gO1) + 3.872983346207417*g1*Conj(ZN1(gI1,0))*(-(Conj(
      ZN1(gI2,2))*KroneckerDelta(0,gO1)) + Conj(ZN1(gI2,3))*KroneckerDelta(1,gO1))
      );

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1316;
   std::complex<double> tmp_1317;
   std::complex<double> tmp_1318;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1318 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1317 += tmp_1318;
   tmp_1316 += (std::complex<double>(0,-0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1317;
   std::complex<double> tmp_1319;
   std::complex<double> tmp_1320;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1320 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1319 += tmp_1320;
   tmp_1316 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1319;
   std::complex<double> tmp_1321;
   std::complex<double> tmp_1322;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1322 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1321 += tmp_1322;
   tmp_1316 += (std::complex<double>(0,0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1321;
   std::complex<double> tmp_1323;
   std::complex<double> tmp_1324;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1324 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1323 += tmp_1324;
   tmp_1316 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1323;
   std::complex<double> tmp_1325;
   std::complex<double> tmp_1326;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1326 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1325 += tmp_1326;
   tmp_1316 += (std::complex<double>(0,-0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1325;
   std::complex<double> tmp_1327;
   std::complex<double> tmp_1328;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1328 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1327 += tmp_1328;
   tmp_1316 += (std::complex<double>(0,0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1327;
   result += (std::complex<double>(0,-1)) * tmp_1316;

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1329;
   std::complex<double> tmp_1330;
   std::complex<double> tmp_1331;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1331 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1330 += tmp_1331;
   tmp_1329 += (std::complex<double>(0,0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1330;
   std::complex<double> tmp_1332;
   std::complex<double> tmp_1333;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1333 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1332 += tmp_1333;
   tmp_1329 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1332;
   std::complex<double> tmp_1334;
   std::complex<double> tmp_1335;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1335 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1334 += tmp_1335;
   tmp_1329 += (std::complex<double>(0,-0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1334;
   std::complex<double> tmp_1336;
   std::complex<double> tmp_1337;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1337 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1336 += tmp_1337;
   tmp_1329 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1336;
   std::complex<double> tmp_1338;
   std::complex<double> tmp_1339;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1339 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1338 += tmp_1339;
   tmp_1329 += (std::complex<double>(0,-0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1338;
   std::complex<double> tmp_1340;
   std::complex<double> tmp_1341;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1341 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1340 += tmp_1341;
   tmp_1329 += (std::complex<double>(0,0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1340;
   result += (std::complex<double>(0,-1)) * tmp_1329;

   return result;
}

std::complex<double> CLASSNAME::CpURhconjURhconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1342;
   std::complex<double> tmp_1343;
   std::complex<double> tmp_1344;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1344 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1343 += tmp_1344;
   tmp_1342 += (std::complex<double>(0,-0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1343;
   std::complex<double> tmp_1345;
   std::complex<double> tmp_1346;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1346 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1345 += tmp_1346;
   tmp_1342 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1345;
   std::complex<double> tmp_1347;
   std::complex<double> tmp_1348;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1348 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1347 += tmp_1348;
   tmp_1342 += (std::complex<double>(0,0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1347;
   std::complex<double> tmp_1349;
   std::complex<double> tmp_1350;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1350 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1349 += tmp_1350;
   tmp_1342 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1349;
   std::complex<double> tmp_1351;
   std::complex<double> tmp_1352;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1352 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1351 += tmp_1352;
   tmp_1342 += (std::complex<double>(0,0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1351;
   std::complex<double> tmp_1353;
   std::complex<double> tmp_1354;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1354 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1353 += tmp_1354;
   tmp_1342 += (std::complex<double>(0,-0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1353;
   result += (std::complex<double>(0,-1)) * tmp_1342;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjSdSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1355;
   std::complex<double> tmp_1356;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1357;
      std::complex<double> tmp_1358;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1358 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1357 += tmp_1358;
      tmp_1356 += (Conj(ZD(gI2,j2))) * tmp_1357;
   }
   tmp_1355 += tmp_1356;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) + vT*Conj(LamTD) + 2*Conj(
      MuD))*KroneckerDelta(0,gO2)) * tmp_1355;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjSeSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1359;
   std::complex<double> tmp_1360;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1361;
      std::complex<double> tmp_1362;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1362 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1361 += tmp_1362;
      tmp_1360 += (Conj(ZE(gI2,j2))) * tmp_1361;
   }
   tmp_1359 += tmp_1360;
   result += (0.5*(1.4142135623730951*vS*Conj(LamSD) + vT*Conj(LamTD) + 2*Conj(
      MuD))*KroneckerDelta(0,gO2)) * tmp_1359;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjSuSu(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1363;
   std::complex<double> tmp_1364;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1365;
      std::complex<double> tmp_1366;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1366 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1365 += tmp_1366;
      tmp_1364 += (Conj(ZU(gI2,j2))) * tmp_1365;
   }
   tmp_1363 += tmp_1364;
   result += (-0.5*(1.4142135623730951*vS*Conj(LamSU) - vT*Conj(LamTU) + 2*Conj
      (MuU))*KroneckerDelta(1,gO2)) * tmp_1363;

   return result;
}

std::complex<double> CLASSNAME::CpconjURhVZRh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(
      KroneckerDelta(0,gO2)*ZHR(gI2,0) - KroneckerDelta(1,gO2)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjURhconjVWmRpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*Conj(ZRP(gI2,0))*KroneckerDelta(1,gO2);

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

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-3*Sqr(g1) + 5*
      Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) + (-20*AbsSqr(LamTU) + 3*Sqr(g1) - 5*Sqr(g2))
      *ZHR(gI1,1)*ZHR(gI2,1)) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((-20*
      AbsSqr(LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) + (-3*Sqr(g1) +
      5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1)) - 10*(KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*(-((-2*AbsSqr(LamTD) + Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0))
      + Sqr(g2)*ZHR(gI1,1)*ZHR(gI2,1)) + KroneckerDelta(2,gO1)*KroneckerDelta(2,
      gO2)*(Sqr(g2)*ZHR(gI1,0)*ZHR(gI2,0) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZHR(gI1,1
      )*ZHR(gI2,1))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(Conj(ZRP(gI2,0))*(KroneckerDelta(1,gO1)*KroneckerDelta(1,
      gO2)*(20*AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1) - 5*Sqr(g2)) +
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(3*Sqr(g1) + 5*Sqr(g2)) + 10*(
      KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)*Sqr(g2) - KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*(-2*AbsSqr(LamTU) + Sqr(g2))))*ZRP(gI1,0)) + Conj(ZRP(
      gI2,1))*(-(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2
      ))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*(-20*AbsSqr(LamSD) - 10*
      AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2)) + 10*(-(KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) + KroneckerDelta(2,gO1)*KroneckerDelta(2,gO2)
      *(-2*AbsSqr(LamTD) + Sqr(g2))))*ZRP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjRhRpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*Conj(ZRP(gI2,0))*(2.8284271247461903*LamSU*vd*Conj(LamSD)*
      KroneckerDelta(1,gO2)*ZHR(gI1,0) + 1.4142135623730951*LamTU*Conj(LamTD)*(2*
      vu*KroneckerDelta(0,gO2) + vd*KroneckerDelta(1,gO2))*ZHR(gI1,0) + (
      1.4142135623730951*vd*KroneckerDelta(0,gO2)*Sqr(g2) + 1.4142135623730951*vu*
      KroneckerDelta(1,gO2)*(-2*AbsSqr(LamSU) + AbsSqr(LamTU) + Sqr(g2)) + 2*((-((
      2*MuU + 1.4142135623730951*LamSU*vS + LamTU*vT)*Conj(LamTU)) + g2*(g2*vT + 2
      *Conj(MDWBT)))*KroneckerDelta(2,gO2) - KroneckerDelta(3,gO2)*(-2*g2*MDWBT -
      vT*AbsSqr(LamTU) + 1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*LamTU*Conj(
      MuU) + vT*Sqr(g2))))*ZHR(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmRpmRh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*Conj(ZRP(gI1,1))*(2*((-((2*MuD + 1.4142135623730951*LamSD*vS
      + LamTD*vT)*Conj(LamTD)) + g2*(g2*vT + 2*Conj(MDWBT)))*KroneckerDelta(2,gO2)
      - KroneckerDelta(3,gO2)*(-2*g2*MDWBT - vT*AbsSqr(LamTD) +
      1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*LamTD*Conj(MuD) + vT*Sqr(g2)))*
      ZHR(gI2,0) + 1.4142135623730951*KroneckerDelta(1,gO2)*(vu*Sqr(g2)*ZHR(gI2,0)
      + 2*LamTU*vd*Conj(LamTD)*ZHR(gI2,1)) + 1.4142135623730951*KroneckerDelta(0,
      gO2)*(vd*(-2*AbsSqr(LamSD) + AbsSqr(LamTD) + Sqr(g2))*ZHR(gI2,0) + vu*(2*
      LamSU*Conj(LamSD) + LamTU*Conj(LamTD))*ZHR(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjRhHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(LamTU)*KroneckerDelta(2,gO2)*Mu*ZHR(gI1,1)*ZP(gI2,0)) + Conj
      (LamTD)*KroneckerDelta(1,gO2)*Mu*ZHR(gI1,0)*ZP(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmRpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(Conj(Mu)*Conj(ZRP(gI1,0))*(
      1.4142135623730951*LamTU*KroneckerDelta(3,gO2)*ZA(gI2,0) + KroneckerDelta(0,
      gO2)*(1.4142135623730951*LamSU*ZA(gI2,2) + LamTU*ZA(gI2,3))) + Conj(ZRP(gI1,
      1))*Mu*(1.4142135623730951*Conj(LamSD)*KroneckerDelta(1,gO2)*ZA(gI2,2) +
      Conj(LamTD)*(1.4142135623730951*KroneckerDelta(2,gO2)*ZA(gI2,1) -
      KroneckerDelta(1,gO2)*ZA(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmRpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(Mu)*Conj(ZRP(gI1,0))*(-1.4142135623730951*LamTU*
      KroneckerDelta(3,gO2)*ZH(gI2,0) + KroneckerDelta(0,gO2)*(1.4142135623730951*
      LamSU*ZH(gI2,2) + LamTU*ZH(gI2,3))) + Conj(ZRP(gI1,1))*Mu*(
      -1.4142135623730951*Conj(LamSD)*KroneckerDelta(1,gO2)*ZH(gI2,2) + Conj(LamTD
      )*(1.4142135623730951*KroneckerDelta(2,gO2)*ZH(gI2,1) + KroneckerDelta(1,gO2
      )*ZH(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarCha1ChiPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(Conj(LamSD)*KroneckerDelta(0,gO2)*UP1(gI1,1)*ZN2(gI2,0)) +
      1.4142135623730951*g2*KroneckerDelta(3,gO2)*UP1(gI1,0)*ZN2(gI2,1) +
      0.7071067811865475*Conj(LamTD)*KroneckerDelta(0,gO2)*UP1(gI1,1)*ZN2(gI2,1) -
      Conj(LamTD)*KroneckerDelta(2,gO2)*UP1(gI1,1)*ZN2(gI2,2) - g2*KroneckerDelta
      (1,gO2)*UP1(gI1,0)*ZN2(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarCha1ChiPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(UM1(gI1,0))*(-(LamTU*Conj(ZN1(gI2,3))*KroneckerDelta(1,gO1)) +
      1.4142135623730951*g2*Conj(ZN1(gI2,1))*KroneckerDelta(2,gO1)) + Conj(UM1(
      gI1,1))*(0.5477225575051661*g1*Conj(ZN1(gI2,0))*KroneckerDelta(0,gO1) +
      0.7071067811865475*g2*Conj(ZN1(gI2,1))*KroneckerDelta(0,gO1) + LamTD*Conj(
      ZN1(gI2,2))*KroneckerDelta(3,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1367;
   tmp_1367 += std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1367 += std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1367 += std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g1);
   tmp_1367 += std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1367 += std::complex<double>(0,0.5)*KroneckerDelta(2,gO1)*KroneckerDelta
      (2,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   tmp_1367 += std::complex<double>(0,-0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*KroneckerDelta(gI1,gI2)*Sqr(g2);
   std::complex<double> tmp_1368;
   std::complex<double> tmp_1369;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1370;
      std::complex<double> tmp_1371;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1372;
         std::complex<double> tmp_1373;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1373 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1372 += tmp_1373;
         tmp_1371 += (Conj(ZV(gI2,j2))) * tmp_1372;
      }
      tmp_1370 += tmp_1371;
      tmp_1369 += (ZV(gI1,j3)) * tmp_1370;
   }
   tmp_1368 += tmp_1369;
   tmp_1367 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1368;
   result += (std::complex<double>(0,-1)) * tmp_1367;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1374;
   std::complex<double> tmp_1375;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1376;
      std::complex<double> tmp_1377;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1377 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1376 += tmp_1377;
      tmp_1375 += (ZUL(gI1,j2)) * tmp_1376;
   }
   tmp_1374 += tmp_1375;
   result += (KroneckerDelta(0,gO2)) * tmp_1374;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFuFdPL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1378;
   std::complex<double> tmp_1379;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1380;
      std::complex<double> tmp_1381;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1381 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1380 += tmp_1381;
      tmp_1379 += (Conj(ZDL(gI2,j2))) * tmp_1380;
   }
   tmp_1378 += tmp_1379;
   result += (KroneckerDelta(1,gO1)) * tmp_1378;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarFvFePR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1382;
   std::complex<double> tmp_1383;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1383 += Conj(Ye(j1,gI1))*ZER(gI2,j1);
   }
   tmp_1382 += tmp_1383;
   result += (KroneckerDelta(0,gO2)) * tmp_1382;

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

   std::complex<double> tmp_1384;
   std::complex<double> tmp_1385;
   std::complex<double> tmp_1386;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1386 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1385 += tmp_1386;
   tmp_1384 += (std::complex<double>(0,-0.35355339059327373)*vd*KroneckerDelta(
      0,gO2)*Sqr(g2)) * tmp_1385;
   std::complex<double> tmp_1387;
   std::complex<double> tmp_1388;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1388 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1387 += tmp_1388;
   tmp_1384 += (std::complex<double>(0,-0.35355339059327373)*vu*KroneckerDelta(
      1,gO2)*Sqr(g2)) * tmp_1387;
   std::complex<double> tmp_1389;
   std::complex<double> tmp_1390;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1390 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1389 += tmp_1390;
   tmp_1384 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(2,gO2)*Sqr(g2))
      * tmp_1389;
   std::complex<double> tmp_1391;
   std::complex<double> tmp_1392;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1392 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1391 += tmp_1392;
   tmp_1384 += (std::complex<double>(0,-1)*g2*Conj(MDWBT)*KroneckerDelta(2,gO2)
      ) * tmp_1391;
   std::complex<double> tmp_1393;
   std::complex<double> tmp_1394;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1394 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1393 += tmp_1394;
   tmp_1384 += (std::complex<double>(0,-1)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1393;
   std::complex<double> tmp_1395;
   std::complex<double> tmp_1396;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1396 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1395 += tmp_1396;
   tmp_1384 += (std::complex<double>(0,0.5)*vT*KroneckerDelta(3,gO2)*Sqr(g2)) *
      tmp_1395;
   std::complex<double> tmp_1397;
   std::complex<double> tmp_1398;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1399;
      std::complex<double> tmp_1400;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1400 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1399 += tmp_1400;
      tmp_1398 += (ZV(gI1,j2)) * tmp_1399;
   }
   tmp_1397 += tmp_1398;
   tmp_1384 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1397;
   std::complex<double> tmp_1401;
   std::complex<double> tmp_1402;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1403;
      std::complex<double> tmp_1404;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1405;
         std::complex<double> tmp_1406;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1406 += Conj(Ye(j1,j3))*Ye(j1,j2);
         }
         tmp_1405 += tmp_1406;
         tmp_1404 += (Conj(ZE(gI2,j2))) * tmp_1405;
      }
      tmp_1403 += tmp_1404;
      tmp_1402 += (ZV(gI1,j3)) * tmp_1403;
   }
   tmp_1401 += tmp_1402;
   tmp_1384 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(0,
      gO2)) * tmp_1401;
   result += (std::complex<double>(0,-1)) * tmp_1384;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(5*(-(KroneckerDelta(3,gO1)*(1.4142135623730951*KroneckerDelta
      (0,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,0) + 2*LamSU*Conj(LamTU)*KroneckerDelta(1,
      gO2)*ZA(gI1,2)*ZA(gI2,1) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1
      ,gO2)*ZA(gI1,3)*ZA(gI2,1) - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)
      *ZA(gI1,3)*ZA(gI2,1) + 2*LamSU*Conj(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*
      ZA(gI2,2) + 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(
      gI2,3) + 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA
      (gI2,3) - 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,
      3) + 4*KroneckerDelta(2,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,3) + 2*KroneckerDelta(
      3,gO2)*(Sqr(g2)*ZA(gI1,0)*ZA(gI2,0) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZA(gI1,1)
      *ZA(gI2,1) + 2*Sqr(g2)*ZA(gI1,3)*ZA(gI2,3)) - Conj(LamTD)*KroneckerDelta(0,
      gO2)*(2*LamSD*ZA(gI1,2)*ZA(gI2,0) + 1.4142135623730951*LamTD*ZA(gI1,3)*ZA(
      gI2,0) + ZA(gI1,0)*(2*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3)))
      )) + KroneckerDelta(2,gO1)*(1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(
      0,gO2)*ZA(gI1,3)*ZA(gI2,0) - 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2
      )*ZA(gI1,3)*ZA(gI2,0) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,2)*
      ZA(gI2,1) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,3)
      *ZA(gI2,1) + 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(
      gI2,1) + 2*LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA(gI2,2) - 2*
      LamTD*Conj(LamSD)*KroneckerDelta(0,gO2)*(ZA(gI1,2)*ZA(gI2,0) + ZA(gI1,0)*ZA(
      gI2,2)) + 1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZA(gI1,0)*
      ZA(gI2,3) - 1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZA(gI1,0)*ZA(
      gI2,3) - 1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZA(gI1,1)*ZA
      (gI2,3) + 1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZA(gI1,1)*ZA(gI2,
      3) - 4*KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI1,3)*ZA(gI2,3) + 2*KroneckerDelta(
      2,gO2)*((-2*AbsSqr(LamTD) + Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) - Sqr(g2)*(ZA(gI1,1
      )*ZA(gI2,1) + 2*ZA(gI1,3)*ZA(gI2,3))))) + KroneckerDelta(0,gO1)*(
      KroneckerDelta(0,gO2)*(-((3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0)) + (3*
      Sqr(g1) - 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTD)*(
      1.4142135623730951*LamSD*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSD*ZA(gI2,2) - 2*LamTD*ZA(gI2,3))) + Conj(LamSD)*(1.4142135623730951*
      LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*(-4*LamSD*ZA(gI2,2) +
      1.4142135623730951*LamTD*ZA(gI2,3))))) + 5*(2*LamTD*Conj(LamSD)*
      KroneckerDelta(3,gO2)*ZA(gI1,2)*ZA(gI2,0) - 1.4142135623730951*
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

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*(ZP(gI1,1)*(KroneckerDelta(0,gO2)*(3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI2,0) + 10*(KroneckerDelta(2,gO2)*(-2*AbsSqr(LamTU)
      + Sqr(g2))*ZP(gI2,2) - KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI2,3))) +
      KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - 2*((3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 5*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP
      (gI1,2)*ZP(gI2,2) + 5*Sqr(g2)*ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(0,gO1)
      *(ZP(gI1,0)*(KroneckerDelta(1,gO2)*(3*Sqr(g1) + 5*Sqr(g2))*ZP(gI2,1) - 10*(
      KroneckerDelta(2,gO2)*Sqr(g2)*ZP(gI2,2) - KroneckerDelta(3,gO2)*(-2*AbsSqr(
      LamTD) + Sqr(g2))*ZP(gI2,3))) + KroneckerDelta(0,gO2)*(-2*(3*Sqr(g1) + 5*Sqr
      (g2))*ZP(gI1,0)*ZP(gI2,0) + (3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) - 10
      *(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - (-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,3)*ZP(
      gI2,3)))) - 10*(KroneckerDelta(2,gO1)*(ZP(gI1,2)*(KroneckerDelta(0,gO2)*Sqr(
      g2)*ZP(gI2,0) - KroneckerDelta(1,gO2)*(-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI2,1)
      - 2*KroneckerDelta(3,gO2)*Sqr(g2)*ZP(gI2,3)) + KroneckerDelta(2,gO2)*(Sqr(
      g2)*ZP(gI1,0)*ZP(gI2,0) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) +
      2*Sqr(g2)*(2*ZP(gI1,2)*ZP(gI2,2) - ZP(gI1,3)*ZP(gI2,3)))) + KroneckerDelta(
      3,gO1)*(ZP(gI1,3)*(-(KroneckerDelta(0,gO2)*(-2*AbsSqr(LamTD) + Sqr(g2))*ZP(
      gI2,0)) + Sqr(g2)*(KroneckerDelta(1,gO2)*ZP(gI2,1) - 2*KroneckerDelta(2,gO2)
      *ZP(gI2,2))) + KroneckerDelta(3,gO2)*(-((-2*AbsSqr(LamTD) + Sqr(g2))*ZP(gI1,
      0)*ZP(gI2,0)) + Sqr(g2)*(ZP(gI1,1)*ZP(gI2,1) - 2*ZP(gI1,2)*ZP(gI2,2) + 4*ZP(
      gI1,3)*ZP(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(5*(-(KroneckerDelta(3,gO1)*(-1.4142135623730951*
      KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,0) + 2*LamSU*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,1) + 1.4142135623730951*AbsSqr(LamTU)
      *KroneckerDelta(1,gO2)*ZH(gI1,3)*ZH(gI2,1) - 1.4142135623730951*
      KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,1) + 2*LamSU*Conj(LamTU)*
      KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,2) - 1.4142135623730951*
      KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,3) + 1.4142135623730951*
      AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,3) - 1.4142135623730951
      *KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,3) - 4*KroneckerDelta(2,gO2)
      *Sqr(g2)*ZH(gI1,3)*ZH(gI2,3) + 2*KroneckerDelta(3,gO2)*(Sqr(g2)*ZH(gI1,0)*ZH
      (gI2,0) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZH(gI1,1)*ZH(gI2,1) + 2*Sqr(g2)*ZH(
      gI1,3)*ZH(gI2,3)) + Conj(LamTD)*KroneckerDelta(0,gO2)*(2*LamSD*ZH(gI1,2)*ZH(
      gI2,0) + 1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,0) + ZH(gI1,0)*(2*LamSD*
      ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3))))) + KroneckerDelta(2,gO1)*(
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZH(gI1,3)*ZH(gI2,0) -
      1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,0) - 2*
      LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZH(gI1,2)*ZH(gI2,1) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,3)*ZH(gI2,1) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,3)*ZH(gI2,1) - 2*
      LamTU*Conj(LamSU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,2) - 2*LamTD*Conj(
      LamSD)*KroneckerDelta(0,gO2)*(ZH(gI1,2)*ZH(gI2,0) + ZH(gI1,0)*ZH(gI2,2)) +
      1.4142135623730951*AbsSqr(LamTD)*KroneckerDelta(0,gO2)*ZH(gI1,0)*ZH(gI2,3) -
      1.4142135623730951*KroneckerDelta(0,gO2)*Sqr(g2)*ZH(gI1,0)*ZH(gI2,3) +
      1.4142135623730951*AbsSqr(LamTU)*KroneckerDelta(1,gO2)*ZH(gI1,1)*ZH(gI2,3) -
      1.4142135623730951*KroneckerDelta(1,gO2)*Sqr(g2)*ZH(gI1,1)*ZH(gI2,3) + 4*
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

std::complex<double> CLASSNAME::CpconjUHpmbarChiCha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -1.4142135623730951*g2*KroneckerDelta(3,gO2)*UP2(gI2,0)*ZN1(gI1,1)
      - 0.1414213562373095*KroneckerDelta(1,gO2)*UP2(gI2,1)*(3.872983346207417*g1*
      ZN1(gI1,0) + 5*g2*ZN1(gI1,1)) + Conj(LamTD)*KroneckerDelta(0,gO2)*UP2(gI2,0)
      *ZN1(gI1,2) - Conj(LamTU)*KroneckerDelta(2,gO2)*UP2(gI2,1)*ZN1(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmbarChiCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UM2(gI2,0))*(Conj(ZN2(gI1,2))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN2(gI1,1))*KroneckerDelta(2,gO1))) + 0.5*Conj(UM2(
      gI2,1))*(2*LamSU*Conj(ZN2(gI1,0))*KroneckerDelta(1,gO1) + 1.4142135623730951
      *LamTU*Conj(ZN2(gI1,1))*KroneckerDelta(1,gO1) + 2*LamTU*Conj(ZN2(gI1,3))*
      KroneckerDelta(3,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(14.142135623730951*g2*MDWBT*
      KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) + 7.0710678118654755*vT*AbsSqr(
      LamTD)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) + 10*LamTD*vS*Conj(LamSD)*
      KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) + 14.142135623730951*LamTD*Conj(
      MuD)*KroneckerDelta(3,gO2)*ZA(gI2,0)*ZP(gI1,0) - 7.0710678118654755*vT*
      KroneckerDelta(3,gO2)*Sqr(g2)*ZA(gI2,0)*ZP(gI1,0) + 7.745966692414834*g1*
      MDBS*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,0) + 14.142135623730951*MuD*Conj
      (LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,0) - 7.0710678118654755*LamTD
      *vT*Conj(LamSD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(gI1,0) +
      7.0710678118654755*LamSD*vT*Conj(LamTD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(
      gI1,0) - 7.745966692414834*g1*Conj(MDBS)*KroneckerDelta(0,gO2)*ZA(gI2,2)*ZP(
      gI1,0) - 14.142135623730951*LamSD*Conj(MuD)*KroneckerDelta(0,gO2)*ZA(gI2,2)*
      ZP(gI1,0) + 10*LamTD*vd*Conj(LamSD)*KroneckerDelta(3,gO2)*ZA(gI2,2)*ZP(gI1,0
      ) + 10*g2*MDWBT*KroneckerDelta(0,gO2)*ZA(gI2,3)*ZP(gI1,0) +
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

std::complex<double> CLASSNAME::CpconjUHpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(0,gO2)*(-3*vu*Sqr(g1)*ZH(gI2,1)*ZP(gI1,0) +
      5*vu*Sqr(g2)*ZH(gI2,1)*ZP(gI1,0) - 7.745966692414834*g1*MDBS*ZH(gI2,2)*ZP(
      gI1,0) + 20*vS*AbsSqr(LamSD)*ZH(gI2,2)*ZP(gI1,0) + 14.142135623730951*MuD*
      Conj(LamSD)*ZH(gI2,2)*ZP(gI1,0) - 7.0710678118654755*LamTD*vT*Conj(LamSD)*ZH
      (gI2,2)*ZP(gI1,0) - 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gI2,2)*ZP(gI1
      ,0) - 7.745966692414834*g1*Conj(MDBS)*ZH(gI2,2)*ZP(gI1,0) +
      14.142135623730951*LamSD*Conj(MuD)*ZH(gI2,2)*ZP(gI1,0) - 10*g2*MDWBT*ZH(gI2,
      3)*ZP(gI1,0) + 10*vT*AbsSqr(LamTD)*ZH(gI2,3)*ZP(gI1,0) - 7.0710678118654755*
      LamTD*vS*Conj(LamSD)*ZH(gI2,3)*ZP(gI1,0) - 10*MuD*Conj(LamTD)*ZH(gI2,3)*ZP(
      gI1,0) - 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gI2,3)*ZP(gI1,0) - 10*g2
      *Conj(MDWBT)*ZH(gI2,3)*ZP(gI1,0) - 10*LamTD*Conj(MuD)*ZH(gI2,3)*ZP(gI1,0) +
      5*vd*Sqr(g2)*ZH(gI2,1)*ZP(gI1,1) + 10*LamTD*vd*Conj(LamSD)*ZH(gI2,2)*ZP(gI1,
      2) - 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gI2,3)*ZP(gI1,2) +
      7.0710678118654755*vd*Sqr(g2)*ZH(gI2,3)*ZP(gI1,2) + 10*LamSD*vd*Conj(LamTD)*
      ZH(gI2,2)*ZP(gI1,3) + 7.0710678118654755*vd*AbsSqr(LamTD)*ZH(gI2,3)*ZP(gI1,3
      ) - 7.0710678118654755*vd*Sqr(g2)*ZH(gI2,3)*ZP(gI1,3) + ZH(gI2,0)*(vd*(3*Sqr
      (g1) + 5*Sqr(g2))*ZP(gI1,0) + 5*(vu*Sqr(g2)*ZP(gI1,1) + (2*LamTD*vS*Conj(
      LamSD) + 1.4142135623730951*(2*g2*MDWBT - vT*AbsSqr(LamTD) + 2*LamTD*Conj(
      MuD) + vT*Sqr(g2)))*ZP(gI1,2) + ((2.8284271247461903*MuD + 2*LamSD*vS +
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

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1407;
   std::complex<double> tmp_1408;
   std::complex<double> tmp_1409;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1409 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1408 += tmp_1409;
   tmp_1407 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1408;
   std::complex<double> tmp_1410;
   std::complex<double> tmp_1411;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1411 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1410 += tmp_1411;
   tmp_1407 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1410;
   std::complex<double> tmp_1412;
   std::complex<double> tmp_1413;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1413 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1412 += tmp_1413;
   tmp_1407 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1412;
   std::complex<double> tmp_1414;
   std::complex<double> tmp_1415;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1415 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1414 += tmp_1415;
   tmp_1407 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1414;
   std::complex<double> tmp_1416;
   std::complex<double> tmp_1417;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1417 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1416 += tmp_1417;
   tmp_1407 += (std::complex<double>(0,-0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1416;
   std::complex<double> tmp_1418;
   std::complex<double> tmp_1419;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1419 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1418 += tmp_1419;
   tmp_1407 += (std::complex<double>(0,0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1418;
   std::complex<double> tmp_1420;
   std::complex<double> tmp_1421;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1421 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1420 += tmp_1421;
   tmp_1407 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1420;
   std::complex<double> tmp_1422;
   std::complex<double> tmp_1423;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1423 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1422 += tmp_1423;
   tmp_1407 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1422;
   std::complex<double> tmp_1424;
   std::complex<double> tmp_1425;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1426;
      std::complex<double> tmp_1427;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1428;
         std::complex<double> tmp_1429;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1429 += Conj(Yd(j3,j1))*Yd(j2,j1);
         }
         tmp_1428 += tmp_1429;
         tmp_1427 += (ZD(gI1,3 + j2)) * tmp_1428;
      }
      tmp_1426 += tmp_1427;
      tmp_1425 += (Conj(ZD(gI2,3 + j3))) * tmp_1426;
   }
   tmp_1424 += tmp_1425;
   tmp_1407 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1424;
   std::complex<double> tmp_1430;
   std::complex<double> tmp_1431;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1432;
      std::complex<double> tmp_1433;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1434;
         std::complex<double> tmp_1435;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1435 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1434 += tmp_1435;
         tmp_1433 += (Conj(ZD(gI2,j2))) * tmp_1434;
      }
      tmp_1432 += tmp_1433;
      tmp_1431 += (ZD(gI1,j3)) * tmp_1432;
   }
   tmp_1430 += tmp_1431;
   tmp_1407 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1430;
   result += (std::complex<double>(0,-1)) * tmp_1407;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1436;
   std::complex<double> tmp_1437;
   std::complex<double> tmp_1438;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1438 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1437 += tmp_1438;
   tmp_1436 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1437;
   std::complex<double> tmp_1439;
   std::complex<double> tmp_1440;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1440 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1439 += tmp_1440;
   tmp_1436 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1439;
   std::complex<double> tmp_1441;
   std::complex<double> tmp_1442;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1442 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1441 += tmp_1442;
   tmp_1436 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1441;
   std::complex<double> tmp_1443;
   std::complex<double> tmp_1444;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1444 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1443 += tmp_1444;
   tmp_1436 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1443;
   std::complex<double> tmp_1445;
   std::complex<double> tmp_1446;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1446 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1445 += tmp_1446;
   tmp_1436 += (std::complex<double>(0,-0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1445;
   std::complex<double> tmp_1447;
   std::complex<double> tmp_1448;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1448 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1447 += tmp_1448;
   tmp_1436 += (std::complex<double>(0,0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1447;
   std::complex<double> tmp_1449;
   std::complex<double> tmp_1450;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1450 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1449 += tmp_1450;
   tmp_1436 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1449;
   std::complex<double> tmp_1451;
   std::complex<double> tmp_1452;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1452 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1451 += tmp_1452;
   tmp_1436 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1451;
   std::complex<double> tmp_1453;
   std::complex<double> tmp_1454;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1455;
      std::complex<double> tmp_1456;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1457;
         std::complex<double> tmp_1458;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1458 += Conj(Ye(j3,j1))*Ye(j2,j1);
         }
         tmp_1457 += tmp_1458;
         tmp_1456 += (ZE(gI1,3 + j2)) * tmp_1457;
      }
      tmp_1455 += tmp_1456;
      tmp_1454 += (Conj(ZE(gI2,3 + j3))) * tmp_1455;
   }
   tmp_1453 += tmp_1454;
   tmp_1436 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1453;
   result += (std::complex<double>(0,-1)) * tmp_1436;

   return result;
}

std::complex<double> CLASSNAME::CpUHpmconjUHpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1459;
   std::complex<double> tmp_1460;
   std::complex<double> tmp_1461;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1461 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1460 += tmp_1461;
   tmp_1459 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1460;
   std::complex<double> tmp_1462;
   std::complex<double> tmp_1463;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1463 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1462 += tmp_1463;
   tmp_1459 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1462;
   std::complex<double> tmp_1464;
   std::complex<double> tmp_1465;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1465 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1464 += tmp_1465;
   tmp_1459 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1464;
   std::complex<double> tmp_1466;
   std::complex<double> tmp_1467;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1467 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1466 += tmp_1467;
   tmp_1459 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1466;
   std::complex<double> tmp_1468;
   std::complex<double> tmp_1469;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1469 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1468 += tmp_1469;
   tmp_1459 += (std::complex<double>(0,0.5)*KroneckerDelta(2,gO1)*
      KroneckerDelta(2,gO2)*Sqr(g2)) * tmp_1468;
   std::complex<double> tmp_1470;
   std::complex<double> tmp_1471;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1471 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1470 += tmp_1471;
   tmp_1459 += (std::complex<double>(0,-0.5)*KroneckerDelta(3,gO1)*
      KroneckerDelta(3,gO2)*Sqr(g2)) * tmp_1470;
   std::complex<double> tmp_1472;
   std::complex<double> tmp_1473;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1473 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1472 += tmp_1473;
   tmp_1459 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1472;
   std::complex<double> tmp_1474;
   std::complex<double> tmp_1475;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1475 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1474 += tmp_1475;
   tmp_1459 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1474;
   std::complex<double> tmp_1476;
   std::complex<double> tmp_1477;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1478;
      std::complex<double> tmp_1479;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1480;
         std::complex<double> tmp_1481;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1481 += Conj(Yu(j3,j1))*Yu(j2,j1);
         }
         tmp_1480 += tmp_1481;
         tmp_1479 += (ZU(gI1,3 + j2)) * tmp_1480;
      }
      tmp_1478 += tmp_1479;
      tmp_1477 += (Conj(ZU(gI2,3 + j3))) * tmp_1478;
   }
   tmp_1476 += tmp_1477;
   tmp_1459 += (std::complex<double>(0,-1)*KroneckerDelta(1,gO1)*KroneckerDelta
      (1,gO2)) * tmp_1476;
   std::complex<double> tmp_1482;
   std::complex<double> tmp_1483;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1484;
      std::complex<double> tmp_1485;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1486;
         std::complex<double> tmp_1487;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1487 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1486 += tmp_1487;
         tmp_1485 += (Conj(ZU(gI2,j2))) * tmp_1486;
      }
      tmp_1484 += tmp_1485;
      tmp_1483 += (ZU(gI1,j3)) * tmp_1484;
   }
   tmp_1482 += tmp_1483;
   tmp_1459 += (std::complex<double>(0,-1)*KroneckerDelta(0,gO1)*KroneckerDelta
      (0,gO2)) * tmp_1482;
   result += (std::complex<double>(0,-1)) * tmp_1459;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1488;
   std::complex<double> tmp_1489;
   std::complex<double> tmp_1490;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1490 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1489 += tmp_1490;
   tmp_1488 += (std::complex<double>(0,-0.35355339059327373)*vd*KroneckerDelta(
      0,gO2)*Sqr(g2)) * tmp_1489;
   std::complex<double> tmp_1491;
   std::complex<double> tmp_1492;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1492 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1491 += tmp_1492;
   tmp_1488 += (std::complex<double>(0,-0.35355339059327373)*vu*KroneckerDelta(
      1,gO2)*Sqr(g2)) * tmp_1491;
   std::complex<double> tmp_1493;
   std::complex<double> tmp_1494;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1494 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1493 += tmp_1494;
   tmp_1488 += (std::complex<double>(0,-0.5)*vT*KroneckerDelta(2,gO2)*Sqr(g2))
      * tmp_1493;
   std::complex<double> tmp_1495;
   std::complex<double> tmp_1496;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1496 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1495 += tmp_1496;
   tmp_1488 += (std::complex<double>(0,-1)*g2*Conj(MDWBT)*KroneckerDelta(2,gO2)
      ) * tmp_1495;
   std::complex<double> tmp_1497;
   std::complex<double> tmp_1498;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1498 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1497 += tmp_1498;
   tmp_1488 += (std::complex<double>(0,-1)*g2*MDWBT*KroneckerDelta(3,gO2)) *
      tmp_1497;
   std::complex<double> tmp_1499;
   std::complex<double> tmp_1500;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1500 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1499 += tmp_1500;
   tmp_1488 += (std::complex<double>(0,0.5)*vT*KroneckerDelta(3,gO2)*Sqr(g2)) *
      tmp_1499;
   std::complex<double> tmp_1501;
   std::complex<double> tmp_1502;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1503;
      std::complex<double> tmp_1504;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1504 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1503 += tmp_1504;
      tmp_1502 += (Conj(ZD(gI2,j2))) * tmp_1503;
   }
   tmp_1501 += tmp_1502;
   tmp_1488 += (std::complex<double>(0,1)*Conj(Mu)*KroneckerDelta(0,gO2)) *
      tmp_1501;
   std::complex<double> tmp_1505;
   std::complex<double> tmp_1506;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1507;
      std::complex<double> tmp_1508;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1508 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1507 += tmp_1508;
      tmp_1506 += (ZU(gI1,j2)) * tmp_1507;
   }
   tmp_1505 += tmp_1506;
   tmp_1488 += (std::complex<double>(0,1)*KroneckerDelta(1,gO2)*Mu) * tmp_1505;
   std::complex<double> tmp_1509;
   std::complex<double> tmp_1510;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1511;
      std::complex<double> tmp_1512;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1513;
         std::complex<double> tmp_1514;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1514 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1513 += tmp_1514;
         tmp_1512 += (ZU(gI1,3 + j2)) * tmp_1513;
      }
      tmp_1511 += tmp_1512;
      tmp_1510 += (Conj(ZD(gI2,3 + j3))) * tmp_1511;
   }
   tmp_1509 += tmp_1510;
   tmp_1488 += (std::complex<double>(0,0.7071067811865475)*vu*KroneckerDelta(0,
      gO2)) * tmp_1509;
   std::complex<double> tmp_1515;
   std::complex<double> tmp_1516;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1517;
      std::complex<double> tmp_1518;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1519;
         std::complex<double> tmp_1520;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1520 += Conj(Yd(j3,j1))*Yu(j2,j1);
         }
         tmp_1519 += tmp_1520;
         tmp_1518 += (ZU(gI1,3 + j2)) * tmp_1519;
      }
      tmp_1517 += tmp_1518;
      tmp_1516 += (Conj(ZD(gI2,3 + j3))) * tmp_1517;
   }
   tmp_1515 += tmp_1516;
   tmp_1488 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(1,
      gO2)) * tmp_1515;
   std::complex<double> tmp_1521;
   std::complex<double> tmp_1522;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1523;
      std::complex<double> tmp_1524;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1525;
         std::complex<double> tmp_1526;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1526 += Conj(Yd(j1,j3))*Yd(j1,j2);
         }
         tmp_1525 += tmp_1526;
         tmp_1524 += (Conj(ZD(gI2,j2))) * tmp_1525;
      }
      tmp_1523 += tmp_1524;
      tmp_1522 += (ZU(gI1,j3)) * tmp_1523;
   }
   tmp_1521 += tmp_1522;
   tmp_1488 += (std::complex<double>(0,0.7071067811865475)*vd*KroneckerDelta(0,
      gO2)) * tmp_1521;
   std::complex<double> tmp_1527;
   std::complex<double> tmp_1528;
   for (unsigned j3 = 0; j3 < 3; ++j3) {
      std::complex<double> tmp_1529;
      std::complex<double> tmp_1530;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         std::complex<double> tmp_1531;
         std::complex<double> tmp_1532;
         for (unsigned j1 = 0; j1 < 3; ++j1) {
            tmp_1532 += Conj(Yu(j1,j3))*Yu(j1,j2);
         }
         tmp_1531 += tmp_1532;
         tmp_1530 += (Conj(ZD(gI2,j2))) * tmp_1531;
      }
      tmp_1529 += tmp_1530;
      tmp_1528 += (ZU(gI1,j3)) * tmp_1529;
   }
   tmp_1527 += tmp_1528;
   tmp_1488 += (std::complex<double>(0,0.7071067811865475)*vu*KroneckerDelta(1,
      gO2)) * tmp_1527;
   result += (std::complex<double>(0,-1)) * tmp_1488;

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmAh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(KroneckerDelta(0,gO2)*ZA(gI2,0) +
      KroneckerDelta(1,gO2)*ZA(gI2,1) + 1.4142135623730951*(KroneckerDelta(2,gO2)
      - KroneckerDelta(3,gO2))*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjUHpmVWmhh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(KroneckerDelta(0,gO2)*ZH(gI2,0) - KroneckerDelta(1,gO2)*ZH(
      gI2,1) + 1.4142135623730951*(KroneckerDelta(2,gO2) + KroneckerDelta(3,gO2))*
      ZH(gI2,3));

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

std::complex<double> CLASSNAME::CpURpmconjURpmVZVZ(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(
      1,gO1)*KroneckerDelta(1,gO2))*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*
      Sqr(g1) + 5*Sqr(g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)));

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjVWmVWm(unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result = 0.5*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) + KroneckerDelta(1
      ,gO1)*KroneckerDelta(1,gO2))*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjRhRh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((3*Sqr(g1) + 5
      *Sqr(g2))*ZHR(gI1,0)*ZHR(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,1)*ZHR(
      gI2,1))) + KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*Sqr(g1) - 5*Sqr(
      g2))*ZHR(gI1,0)*ZHR(gI2,0) - (3*Sqr(g1) + 5*Sqr(g2))*ZHR(gI1,1)*ZHR(gI2,1)))
      ;

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjRpmRpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(3*Sqr(g1) + 5*Sqr(g2))*(-(Conj(ZRP(gI2,1))*(-2*
      KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*ZRP(gI1,1) + KroneckerDelta(0,
      gO1)*(KroneckerDelta(1,gO2)*ZRP(gI1,0) + KroneckerDelta(0,gO2)*ZRP(gI1,1))))
      + Conj(ZRP(gI2,0))*(2*KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*ZRP(gI1,0
      ) - KroneckerDelta(1,gO1)*(KroneckerDelta(1,gO2)*ZRP(gI1,0) + KroneckerDelta
      (0,gO2)*ZRP(gI1,1))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmconjRhVWm(unsigned gO2, unsigned gI1) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*KroneckerDelta(1,gO2)*ZHR(gI1,0);

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmconjRhHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*KroneckerDelta(1,gO2)*(1.4142135623730951*ZHR(gI1,1)*(2*LamSD
      *vu*Conj(LamSU)*ZP(gI2,0) + LamTD*Conj(LamTU)*(vu*ZP(gI2,0) + 2*vd*ZP(gI2,1)
      )) + ZHR(gI1,0)*(1.4142135623730951*vd*(-2*AbsSqr(LamSD) + AbsSqr(LamTD) +
      Sqr(g2))*ZP(gI2,0) + 1.4142135623730951*vu*Sqr(g2)*ZP(gI2,1) + 4*g2*MDWBT*ZP
      (gI2,2) - 2*vT*AbsSqr(LamTD)*ZP(gI2,2) - 2.8284271247461903*LamTD*vS*Conj(
      LamSD)*ZP(gI2,2) - 4*LamTD*Conj(MuD)*ZP(gI2,2) + 2*vT*Sqr(g2)*ZP(gI2,2) + 2*
      vT*AbsSqr(LamTD)*ZP(gI2,3) - 4*MuD*Conj(LamTD)*ZP(gI2,3) -
      2.8284271247461903*LamSD*vS*Conj(LamTD)*ZP(gI2,3) + 4*g2*Conj(MDWBT)*ZP(gI2,
      3) - 2*vT*Sqr(g2)*ZP(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmRhHpm(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.25*KroneckerDelta(0,gO2)*(2.8284271247461903*LamSD*vd*Conj(LamSU
      )*ZHR(gI1,0)*ZP(gI2,1) + ZHR(gI1,1)*(1.4142135623730951*vd*Sqr(g2)*ZP(gI2,0)
      + 1.4142135623730951*vu*(-2*AbsSqr(LamSU) + Sqr(g2))*ZP(gI2,1) + 4*g2*MDWBT
      *ZP(gI2,2) - 2.8284271247461903*LamTU*vS*Conj(LamSU)*ZP(gI2,2) - 4*LamTU*
      Conj(MuU)*ZP(gI2,2) + 2*vT*Sqr(g2)*ZP(gI2,2) + 4*g2*Conj(MDWBT)*ZP(gI2,3) -
      2*vT*Sqr(g2)*ZP(gI2,3)) + Conj(LamTU)*(1.4142135623730951*LamTD*ZHR(gI1,0)*(
      2*vu*ZP(gI2,0) + vd*ZP(gI2,1)) + ZHR(gI1,1)*(1.4142135623730951*LamTU*vu*ZP(
      gI2,1) - 2*(LamTU*vT*ZP(gI2,2) + (2*MuU + 1.4142135623730951*LamSU*vS -
      LamTU*vT)*ZP(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmRpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.05)*(Conj(ZRP(gI1,1))*KroneckerDelta(1,gO2
      )*((-7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuD - LamTD*vT)*Conj(
      LamSD) + 7.0710678118654755*LamSD*vT*Conj(LamTD) + 7.745966692414834*g1*Conj
      (MDBS) - 14.142135623730951*LamSD*Conj(MuD))*ZA(gI2,2) - 5*(2*g2*MDWBT -
      1.4142135623730951*LamTD*vS*Conj(LamSD) + 2*MuD*Conj(LamTD) +
      1.4142135623730951*LamSD*vS*Conj(LamTD) - 2*g2*Conj(MDWBT) - 2*LamTD*Conj(
      MuD))*ZA(gI2,3)) + Conj(ZRP(gI1,0))*KroneckerDelta(0,gO2)*((
      7.745966692414834*g1*MDBS + 7.0710678118654755*(2*MuU + LamTU*vT)*Conj(LamSU
      ) - 7.0710678118654755*LamSU*vT*Conj(LamTU) - 7.745966692414834*g1*Conj(MDBS
      ) - 14.142135623730951*LamSU*Conj(MuU))*ZA(gI2,2) + 5*(2*g2*MDWBT -
      1.4142135623730951*LamTU*vS*Conj(LamSU) + 2*MuU*Conj(LamTU) +
      1.4142135623730951*LamSU*vS*Conj(LamTU) - 2*g2*Conj(MDWBT) - 2*LamTU*Conj(
      MuU))*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmRpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(Conj(ZRP(gI1,1))*KroneckerDelta(1,gO2)*(vd*(-20*AbsSqr(LamTD)
      + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI2,0) + (-3*vu*Sqr(g1) + 5*vu*Sqr(g2))*ZH(gI2,
      1) - 7.745966692414834*g1*MDBS*ZH(gI2,2) - 20*vS*AbsSqr(LamSD)*ZH(gI2,2) -
      14.142135623730951*MuD*Conj(LamSD)*ZH(gI2,2) + 7.0710678118654755*LamTD*vT*
      Conj(LamSD)*ZH(gI2,2) + 7.0710678118654755*LamSD*vT*Conj(LamTD)*ZH(gI2,2) -
      7.745966692414834*g1*Conj(MDBS)*ZH(gI2,2) - 14.142135623730951*LamSD*Conj(
      MuD)*ZH(gI2,2) - 10*g2*MDWBT*ZH(gI2,3) - 10*vT*AbsSqr(LamTD)*ZH(gI2,3) +
      7.0710678118654755*LamTD*vS*Conj(LamSD)*ZH(gI2,3) + 10*MuD*Conj(LamTD)*ZH(
      gI2,3) + 7.0710678118654755*LamSD*vS*Conj(LamTD)*ZH(gI2,3) - 10*g2*Conj(
      MDWBT)*ZH(gI2,3) + 10*LamTD*Conj(MuD)*ZH(gI2,3)) - Conj(ZRP(gI1,0))*
      KroneckerDelta(0,gO2)*(vd*(3*Sqr(g1) - 5*Sqr(g2))*ZH(gI2,0) + vu*(20*AbsSqr(
      LamTU) - 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI2,1) - 7.745966692414834*g1*MDBS*ZH(gI2
      ,2) + 20*vS*AbsSqr(LamSU)*ZH(gI2,2) + 14.142135623730951*MuU*Conj(LamSU)*ZH(
      gI2,2) + 7.0710678118654755*LamTU*vT*Conj(LamSU)*ZH(gI2,2) +
      7.0710678118654755*LamSU*vT*Conj(LamTU)*ZH(gI2,2) - 7.745966692414834*g1*
      Conj(MDBS)*ZH(gI2,2) + 14.142135623730951*LamSU*Conj(MuU)*ZH(gI2,2) - 10*g2*
      MDWBT*ZH(gI2,3) + 10*vT*AbsSqr(LamTU)*ZH(gI2,3) + 7.0710678118654755*LamTU*
      vS*Conj(LamSU)*ZH(gI2,3) + 10*MuU*Conj(LamTU)*ZH(gI2,3) + 7.0710678118654755
      *LamSU*vS*Conj(LamTU)*ZH(gI2,3) - 10*g2*Conj(MDWBT)*ZH(gI2,3) + 10*LamTU*
      Conj(MuU)*ZH(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjSvSv(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.05*(KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2) - KroneckerDelta
      (1,gO1)*KroneckerDelta(1,gO2))*KroneckerDelta(gI1,gI2)*(3*Sqr(g1) - 5*Sqr(g2
      ));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmconjSvSe(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1533;
   std::complex<double> tmp_1534;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1535;
      std::complex<double> tmp_1536;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1536 += Conj(Ye(j1,j2))*Conj(ZE(gI2,3 + j1));
      }
      tmp_1535 += tmp_1536;
      tmp_1534 += (ZV(gI1,j2)) * tmp_1535;
   }
   tmp_1533 += tmp_1534;
   result += (0.5*(2*MuD + 1.4142135623730951*LamSD*vS - LamTD*vT)*
      KroneckerDelta(1,gO2)) * tmp_1533;

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmAhAh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-20*AbsSqr(
      LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,0)*ZA(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2
      ))*ZA(gI1,1)*ZA(gI2,1) + 5*(Conj(LamTD)*(1.4142135623730951*LamSD*ZA(gI1,2)*
      ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951*LamSD*ZA(gI2,2) - 2*LamTD*ZA(gI2,3
      ))) + Conj(LamSD)*(1.4142135623730951*LamTD*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)*
      (-4*LamSD*ZA(gI2,2) + 1.4142135623730951*LamTD*ZA(gI2,3))))) -
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZA(gI1,
      0)*ZA(gI2,0) + (20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2))*ZA(gI1,1)*ZA(gI2,1
      ) + 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZA(gI1,3)*ZA(gI2,2) + ZA(gI1,2)
      *(4*LamSU*ZA(gI2,2) + 1.4142135623730951*LamTU*ZA(gI2,3))) + Conj(LamTU)*(
      1.4142135623730951*LamSU*ZA(gI1,2)*ZA(gI2,3) + ZA(gI1,3)*(1.4142135623730951
      *LamSU*ZA(gI2,2) + 2*LamTU*ZA(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjHpmHpm(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-20*AbsSqr(
      LamSD) - 10*AbsSqr(LamTD) + 3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) - (3*
      Sqr(g1) + 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1) + 10*((-2*AbsSqr(LamTD) + Sqr(g2))*
      ZP(gI1,2)*ZP(gI2,2) - Sqr(g2)*ZP(gI1,3)*ZP(gI2,3))) - KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*((3*Sqr(g1) + 5*Sqr(g2))*ZP(gI1,0)*ZP(gI2,0) + (20*
      AbsSqr(LamSU) + 10*AbsSqr(LamTU) - 3*Sqr(g1) - 5*Sqr(g2))*ZP(gI1,1)*ZP(gI2,1
      ) + 10*(Sqr(g2)*ZP(gI1,2)*ZP(gI2,2) - (-2*AbsSqr(LamTU) + Sqr(g2))*ZP(gI1,3)
      *ZP(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmhhhh(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(KroneckerDelta(1,gO1)*KroneckerDelta(1,gO2)*((-20*AbsSqr(
      LamTD) + 3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,0)*ZH(gI2,0) + (-3*Sqr(g1) + 5*Sqr(g2
      ))*ZH(gI1,1)*ZH(gI2,1) + 5*(Conj(LamTD)*(1.4142135623730951*LamSD*ZH(gI1,2)*
      ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951*LamSD*ZH(gI2,2) - 2*LamTD*ZH(gI2,3
      ))) + Conj(LamSD)*(1.4142135623730951*LamTD*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)*
      (-4*LamSD*ZH(gI2,2) + 1.4142135623730951*LamTD*ZH(gI2,3))))) -
      KroneckerDelta(0,gO1)*KroneckerDelta(0,gO2)*((3*Sqr(g1) - 5*Sqr(g2))*ZH(gI1,
      0)*ZH(gI2,0) + (20*AbsSqr(LamTU) - 3*Sqr(g1) + 5*Sqr(g2))*ZH(gI1,1)*ZH(gI2,1
      ) + 5*(Conj(LamSU)*(1.4142135623730951*LamTU*ZH(gI1,3)*ZH(gI2,2) + ZH(gI1,2)
      *(4*LamSU*ZH(gI2,2) + 1.4142135623730951*LamTU*ZH(gI2,3))) + Conj(LamTU)*(
      1.4142135623730951*LamSU*ZH(gI1,2)*ZH(gI2,3) + ZH(gI1,3)*(1.4142135623730951
      *LamSU*ZH(gI2,2) + 2*LamTU*ZH(gI2,3))))));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmbarChibarCha1PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*KroneckerDelta(1,gO2)*(UP1(gI2,1)*(1.0954451150103321*g1*ZN1(
      gI1,0) + 1.4142135623730951*g2*ZN1(gI1,1)) + 2*g2*UP1(gI2,0)*ZN1(gI1,2));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmbarChibarCha1PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(Conj(UM1(gI2,1))*(2*LamSD*Conj(ZN2(gI1,0)) -
      1.4142135623730951*LamTD*Conj(ZN2(gI1,1))) + 2*LamTD*Conj(UM1(gI2,0))*Conj(
      ZN2(gI1,2)))*KroneckerDelta(1,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmChiCha2PR(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*KroneckerDelta(0,gO2)*(2*Conj(LamSU)*UP2(gI2,1)*ZN2(gI1,0) +
      Conj(LamTU)*(1.4142135623730951*UP2(gI2,1)*ZN2(gI1,1) + 2*UP2(gI2,0)*ZN2(gI1
      ,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmChiCha2PL(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(Conj(UM2(gI2,1))*(1.0954451150103321*g1*Conj(ZN1(gI1,0)) +
      1.4142135623730951*g2*Conj(ZN1(gI1,1))) - 2*g2*Conj(UM2(gI2,0))*Conj(ZN1(gI1
      ,3)))*KroneckerDelta(0,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmHpmAh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(1.4142135623730951*Conj(LamSU)*
      KroneckerDelta(0,gO2)*Mu*ZA(gI2,2)*ZP(gI1,0) + Conj(Mu)*KroneckerDelta(1,gO2
      )*(1.4142135623730951*LamSD*ZA(gI2,2)*ZP(gI1,1) - LamTD*ZA(gI2,3)*ZP(gI1,1)
      + 1.4142135623730951*LamTD*ZA(gI2,1)*ZP(gI1,2)) + Conj(LamTU)*KroneckerDelta
      (0,gO2)*Mu*(ZA(gI2,3)*ZP(gI1,0) + 1.4142135623730951*ZA(gI2,0)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmHpmhh(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(1.4142135623730951*Conj(LamSU)*KroneckerDelta(0,gO2)*Mu*ZH(gI2
      ,2)*ZP(gI1,0) + Conj(Mu)*KroneckerDelta(1,gO2)*(-1.4142135623730951*LamSD*ZH
      (gI2,2)*ZP(gI1,1) + LamTD*ZH(gI2,3)*ZP(gI1,1) + 1.4142135623730951*LamTD*ZH(
      gI2,1)*ZP(gI1,2)) + Conj(LamTU)*KroneckerDelta(0,gO2)*Mu*(ZH(gI2,3)*ZP(gI1,0
      ) - 1.4142135623730951*ZH(gI2,0)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjSdSd(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1537;
   std::complex<double> tmp_1538;
   std::complex<double> tmp_1539;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1539 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1538 += tmp_1539;
   tmp_1537 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1538;
   std::complex<double> tmp_1540;
   std::complex<double> tmp_1541;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1541 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1540 += tmp_1541;
   tmp_1537 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1540;
   std::complex<double> tmp_1542;
   std::complex<double> tmp_1543;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1543 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1542 += tmp_1543;
   tmp_1537 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1542;
   std::complex<double> tmp_1544;
   std::complex<double> tmp_1545;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1545 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1544 += tmp_1545;
   tmp_1537 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1544;
   std::complex<double> tmp_1546;
   std::complex<double> tmp_1547;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1547 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1546 += tmp_1547;
   tmp_1537 += (std::complex<double>(0,0.1)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1546;
   std::complex<double> tmp_1548;
   std::complex<double> tmp_1549;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1549 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1548 += tmp_1549;
   tmp_1537 += (std::complex<double>(0,-0.1)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1548;
   result += (std::complex<double>(0,-1)) * tmp_1537;

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjSeSe(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1550;
   std::complex<double> tmp_1551;
   std::complex<double> tmp_1552;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1552 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1551 += tmp_1552;
   tmp_1550 += (std::complex<double>(0,-0.15)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1551;
   std::complex<double> tmp_1553;
   std::complex<double> tmp_1554;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1554 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1553 += tmp_1554;
   tmp_1550 += (std::complex<double>(0,-0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1553;
   std::complex<double> tmp_1555;
   std::complex<double> tmp_1556;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1556 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1555 += tmp_1556;
   tmp_1550 += (std::complex<double>(0,0.15)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1555;
   std::complex<double> tmp_1557;
   std::complex<double> tmp_1558;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1558 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1557 += tmp_1558;
   tmp_1550 += (std::complex<double>(0,0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1557;
   std::complex<double> tmp_1559;
   std::complex<double> tmp_1560;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1560 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1559 += tmp_1560;
   tmp_1550 += (std::complex<double>(0,0.3)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1559;
   std::complex<double> tmp_1561;
   std::complex<double> tmp_1562;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1562 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1561 += tmp_1562;
   tmp_1550 += (std::complex<double>(0,-0.3)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1561;
   result += (std::complex<double>(0,-1)) * tmp_1550;

   return result;
}

std::complex<double> CLASSNAME::CpURpmconjURpmconjSuSu(unsigned gO1, unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1563;
   std::complex<double> tmp_1564;
   std::complex<double> tmp_1565;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1565 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1564 += tmp_1565;
   tmp_1563 += (std::complex<double>(0,0.05)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1564;
   std::complex<double> tmp_1566;
   std::complex<double> tmp_1567;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1567 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1566 += tmp_1567;
   tmp_1563 += (std::complex<double>(0,0.25)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g2)) * tmp_1566;
   std::complex<double> tmp_1568;
   std::complex<double> tmp_1569;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1569 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1568 += tmp_1569;
   tmp_1563 += (std::complex<double>(0,-0.05)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1568;
   std::complex<double> tmp_1570;
   std::complex<double> tmp_1571;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1571 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1570 += tmp_1571;
   tmp_1563 += (std::complex<double>(0,-0.25)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g2)) * tmp_1570;
   std::complex<double> tmp_1572;
   std::complex<double> tmp_1573;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1573 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1572 += tmp_1573;
   tmp_1563 += (std::complex<double>(0,-0.2)*KroneckerDelta(0,gO1)*
      KroneckerDelta(0,gO2)*Sqr(g1)) * tmp_1572;
   std::complex<double> tmp_1574;
   std::complex<double> tmp_1575;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1575 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1574 += tmp_1575;
   tmp_1563 += (std::complex<double>(0,0.2)*KroneckerDelta(1,gO1)*
      KroneckerDelta(1,gO2)*Sqr(g1)) * tmp_1574;
   result += (std::complex<double>(0,-1)) * tmp_1563;

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmconjSuSd(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1576;
   std::complex<double> tmp_1577;
   std::complex<double> tmp_1578;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1579;
      std::complex<double> tmp_1580;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1580 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1579 += tmp_1580;
      tmp_1578 += (Conj(ZD(gI2,j2))) * tmp_1579;
   }
   tmp_1577 += tmp_1578;
   tmp_1576 += (std::complex<double>(0,-0.7071067811865475)*vS*Conj(LamSU)*
      KroneckerDelta(0,gO2)) * tmp_1577;
   std::complex<double> tmp_1581;
   std::complex<double> tmp_1582;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1583;
      std::complex<double> tmp_1584;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1584 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1583 += tmp_1584;
      tmp_1582 += (Conj(ZD(gI2,j2))) * tmp_1583;
   }
   tmp_1581 += tmp_1582;
   tmp_1576 += (std::complex<double>(0,-0.5)*vT*Conj(LamTU)*KroneckerDelta(0,
      gO2)) * tmp_1581;
   std::complex<double> tmp_1585;
   std::complex<double> tmp_1586;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1587;
      std::complex<double> tmp_1588;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1588 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1587 += tmp_1588;
      tmp_1586 += (Conj(ZD(gI2,j2))) * tmp_1587;
   }
   tmp_1585 += tmp_1586;
   tmp_1576 += (std::complex<double>(0,-1)*Conj(MuU)*KroneckerDelta(0,gO2)) *
      tmp_1585;
   std::complex<double> tmp_1589;
   std::complex<double> tmp_1590;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1591;
      std::complex<double> tmp_1592;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1592 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1591 += tmp_1592;
      tmp_1590 += (ZU(gI1,j2)) * tmp_1591;
   }
   tmp_1589 += tmp_1590;
   tmp_1576 += (std::complex<double>(0,1)*MuD*KroneckerDelta(1,gO2)) * tmp_1589
      ;
   std::complex<double> tmp_1593;
   std::complex<double> tmp_1594;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1595;
      std::complex<double> tmp_1596;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1596 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1595 += tmp_1596;
      tmp_1594 += (ZU(gI1,j2)) * tmp_1595;
   }
   tmp_1593 += tmp_1594;
   tmp_1576 += (std::complex<double>(0,0.7071067811865475)*LamSD*vS*
      KroneckerDelta(1,gO2)) * tmp_1593;
   std::complex<double> tmp_1597;
   std::complex<double> tmp_1598;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1599;
      std::complex<double> tmp_1600;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1600 += Conj(Yd(j1,j2))*Conj(ZD(gI2,3 + j1));
      }
      tmp_1599 += tmp_1600;
      tmp_1598 += (ZU(gI1,j2)) * tmp_1599;
   }
   tmp_1597 += tmp_1598;
   tmp_1576 += (std::complex<double>(0,-0.5)*LamTD*vT*KroneckerDelta(1,gO2)) *
      tmp_1597;
   result += (std::complex<double>(0,-1)) * tmp_1576;

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmVWmRh(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*KroneckerDelta(0,gO2)*ZHR(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmVPRpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 2) {
      result += -0.3872983346207417*g1*Conj(ZRP(gI2,gO2))*Cos(ThetaW());
   }
   if (gI2 < 2) {
      result += -0.5*g2*Conj(ZRP(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

std::complex<double> CLASSNAME::CpconjURpmVZRpm(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   if (gI2 < 2) {
      result += -0.5*g2*Conj(ZRP(gI2,gO2))*Cos(ThetaW());
   }
   if (gI2 < 2) {
      result += 0.3872983346207417*g1*Conj(ZRP(gI2,gO2))*Sin(ThetaW());
   }

   return result;
}

double CLASSNAME::CpSOcconjSOcconjSOcSOc() const
{
   double result = 0.0;

   result = -3*Sqr(g3);

   return result;
}

std::complex<double> CLASSNAME::CpconjSOcVGSOc() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

double CLASSNAME::CpSOcconjSOcconjSdSd(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

double CLASSNAME::CpSOcconjSOcconjSuSu(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjSOcconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1601;
   std::complex<double> tmp_1602;
   std::complex<double> tmp_1603;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1603 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1602 += tmp_1603;
   tmp_1601 += (-1) * tmp_1602;
   std::complex<double> tmp_1604;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1604 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1601 += tmp_1604;
   result += (1.4142135623730951*g3*Conj(MDGoc)) * tmp_1601;

   return result;
}

std::complex<double> CLASSNAME::CpconjSOcconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1605;
   std::complex<double> tmp_1606;
   std::complex<double> tmp_1607;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1607 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1606 += tmp_1607;
   tmp_1605 += (-1) * tmp_1606;
   std::complex<double> tmp_1608;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1608 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1605 += tmp_1608;
   result += (1.4142135623730951*g3*Conj(MDGoc)) * tmp_1605;

   return result;
}

double CLASSNAME::CpconjSOcbarGluGluPR(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpconjSOcbarGluGluPL(unsigned , unsigned ) const
{
   std::complex<double> result;

   result = std::complex<double>(0,1.4142135623730951)*g3;

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

std::complex<double> CLASSNAME::CpVZVZconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZHR(gI1,0)*ZHR(gI2,0) + ZHR
      (gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjRpmRpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(-7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(Conj(ZRP(gI2,0))*ZRP(gI1,0)
      + Conj(ZRP(gI2,1))*ZRP(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(ZHR(
      gI1,0)*ZHR(gI2,0) - ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

double CLASSNAME::CpVZconjRpmRpm(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.1*KroneckerDelta(gI1,gI2)*(-5*g2*Cos(ThetaW()) +
      3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarCha1Cha1PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-10*g2*Conj(UP1(gI2,0))*Cos(ThetaW())*UP1(gI1,0) + Conj(UP1(
      gI2,1))*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*UP1(gI1,1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarCha1Cha1PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(-10*g2*Conj(UM1(gI1,0))*Cos(ThetaW())*UM1(gI2,0) + Conj(UM1(
      gI1,1))*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*UM1(gI2,1
      ));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarCha2Cha2PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UM2(gI2,0))*Cos(ThetaW())*UM2(gI1,0) + Conj(UM2(gI2,
      1))*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UM2(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarCha2Cha2PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(2*g2*Conj(UP2(gI1,0))*Cos(ThetaW())*UP2(gI2,0) + Conj(UP2(gI1,
      1))*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW()))*UP2(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*KroneckerDelta(gI1,gI2)*(g1*Sin(ThetaW())*(7.745966692414834*g2
      *Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

double CLASSNAME::CpVZconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) + 0.7745966692414834*
      g1*Sin(ThetaW()));

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

std::complex<double> CLASSNAME::CpVZVZAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZA(gI1,0)*ZA(gI2,0) + ZA(
      gI1,1)*ZA(gI2,1));

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

std::complex<double> CLASSNAME::CpVZVZhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.05*(7.745966692414834*g1*g2*Sin(2*ThetaW()) + 3*Sqr(g1) + 5*Sqr(
      g2) + Cos(2*ThetaW())*(-3*Sqr(g1) + 5*Sqr(g2)))*(ZH(gI1,0)*ZH(gI2,0) + ZH(
      gI1,1)*ZH(gI2,1));

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

std::complex<double> CLASSNAME::CpVZhhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(g2*Cos(ThetaW()) + 0.7745966692414834
      *g1*Sin(ThetaW()))*(ZA(gI2,0)*ZH(gI1,0) - ZA(gI2,1)*ZH(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChiChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(Conj(
      ZN1(gI2,2))*ZN1(gI1,2) - Conj(ZN1(gI2,3))*ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpVZbarChiChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))*(Conj(
      ZN2(gI1,2))*ZN2(gI2,2) - Conj(ZN2(gI1,3))*ZN2(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1609;
   std::complex<double> tmp_1610;
   std::complex<double> tmp_1611;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1611 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1610 += tmp_1611;
   tmp_1609 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1610;
   std::complex<double> tmp_1612;
   std::complex<double> tmp_1613;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1613 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1612 += tmp_1613;
   tmp_1609 += (std::complex<double>(0,0.2581988897471611)*g1*g2*Cos(ThetaW())*
      Sin(ThetaW())) * tmp_1612;
   std::complex<double> tmp_1614;
   std::complex<double> tmp_1615;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1615 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1614 += tmp_1615;
   tmp_1609 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1614;
   std::complex<double> tmp_1616;
   std::complex<double> tmp_1617;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1617 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1616 += tmp_1617;
   tmp_1609 += (std::complex<double>(0,0.13333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1616;
   result += (std::complex<double>(0,-1)) * tmp_1609;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1618;
   std::complex<double> tmp_1619;
   std::complex<double> tmp_1620;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1620 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1619 += tmp_1620;
   tmp_1618 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1619;
   std::complex<double> tmp_1621;
   std::complex<double> tmp_1622;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1622 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1621 += tmp_1622;
   tmp_1618 += (std::complex<double>(0,-0.7745966692414834)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1621;
   std::complex<double> tmp_1623;
   std::complex<double> tmp_1624;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1624 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1623 += tmp_1624;
   tmp_1618 += (std::complex<double>(0,0.3)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1623;
   std::complex<double> tmp_1625;
   std::complex<double> tmp_1626;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1626 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1625 += tmp_1626;
   tmp_1618 += (std::complex<double>(0,1.2)*Sqr(g1)*Sqr(Sin(ThetaW()))) *
      tmp_1625;
   result += (std::complex<double>(0,-1)) * tmp_1618;

   return result;
}

std::complex<double> CLASSNAME::CpVZVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1627;
   std::complex<double> tmp_1628;
   std::complex<double> tmp_1629;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1629 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1628 += tmp_1629;
   tmp_1627 += (std::complex<double>(0,0.5)*Sqr(g2)*Sqr(Cos(ThetaW()))) *
      tmp_1628;
   std::complex<double> tmp_1630;
   std::complex<double> tmp_1631;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1631 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1630 += tmp_1631;
   tmp_1627 += (std::complex<double>(0,-0.2581988897471611)*g1*g2*Cos(ThetaW())
      *Sin(ThetaW())) * tmp_1630;
   std::complex<double> tmp_1632;
   std::complex<double> tmp_1633;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1633 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1632 += tmp_1633;
   tmp_1627 += (std::complex<double>(0,0.03333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1632;
   std::complex<double> tmp_1634;
   std::complex<double> tmp_1635;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1635 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1634 += tmp_1635;
   tmp_1627 += (std::complex<double>(0,0.5333333333333333)*Sqr(g1)*Sqr(Sin(
      ThetaW()))) * tmp_1634;
   result += (std::complex<double>(0,-1)) * tmp_1627;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1636;
   std::complex<double> tmp_1637;
   std::complex<double> tmp_1638;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1638 += Conj(ZD(gI2,3 + j1))*ZD(gI1,3 + j1);
   }
   tmp_1637 += tmp_1638;
   tmp_1636 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1637;
   std::complex<double> tmp_1639;
   std::complex<double> tmp_1640;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1640 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1639 += tmp_1640;
   tmp_1636 += (-3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1639;
   result += (0.16666666666666666) * tmp_1636;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1641;
   std::complex<double> tmp_1642;
   std::complex<double> tmp_1643;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1643 += Conj(ZE(gI2,3 + j1))*ZE(gI1,3 + j1);
   }
   tmp_1642 += tmp_1643;
   tmp_1641 += (1.5491933384829668*g1*Sin(ThetaW())) * tmp_1642;
   std::complex<double> tmp_1644;
   std::complex<double> tmp_1645;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1645 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1644 += tmp_1645;
   tmp_1641 += (-(g2*Cos(ThetaW())) + 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1644;
   result += (0.5) * tmp_1641;

   return result;
}

std::complex<double> CLASSNAME::CpVZconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1646;
   std::complex<double> tmp_1647;
   std::complex<double> tmp_1648;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1648 += Conj(ZU(gI2,3 + j1))*ZU(gI1,3 + j1);
   }
   tmp_1647 += tmp_1648;
   tmp_1646 += (-3.0983866769659336*g1*Sin(ThetaW())) * tmp_1647;
   std::complex<double> tmp_1649;
   std::complex<double> tmp_1650;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1650 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1649 += tmp_1650;
   tmp_1646 += (3*g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(ThetaW())) *
      tmp_1649;
   result += (0.16666666666666666) * tmp_1646;

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

std::complex<double> CLASSNAME::CpVWmconjVWmconjRhRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZHR(gI1,0)*ZHR(gI2,0) + ZHR(gI1,1)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjRpmRpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(Conj(ZRP(gI2,0))*ZRP(gI1,0) + Conj(ZRP(gI2,1))*ZRP(gI1
      ,1));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjRhRpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*g2*Conj(ZRP(gI2,0))*ZHR(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmRpmRh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.7071067811865475*g2*Conj(ZRP(gI1,1))*ZHR(gI2,0);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarCha1ChiPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(ZN1(gI2,1))*UP1(gI1,0) - 0.7071067811865475*g2*Conj(ZN1(gI2
      ,2))*UP1(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarCha1ChiPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM1(gI1,0))*ZN2(gI2,1) + 0.7071067811865475*g2*Conj(UM1(gI1
      ,1))*ZN2(gI2,2);

   return result;
}

double CLASSNAME::CpVWmconjVWmconjSvSv(unsigned gI1, unsigned gI2) const
{
   double result = 0.0;

   result = 0.5*KroneckerDelta(gI1,gI2)*Sqr(g2);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarFuFdPL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1651;
   std::complex<double> tmp_1652;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1652 += Conj(ZDL(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1651 += tmp_1652;
   result += (-0.7071067811865475*g2) * tmp_1651;

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

   std::complex<double> tmp_1653;
   std::complex<double> tmp_1654;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1654 += Conj(ZE(gI2,j1))*ZV(gI1,j1);
   }
   tmp_1653 += tmp_1654;
   result += (0.7071067811865475*g2) * tmp_1653;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmAhAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZA(gI1,0)*ZA(gI2,0) + ZA(gI1,1)*ZA(gI2,1) + 4*ZA(gI1,3
      )*ZA(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjHpmHpm(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZP(gI1,0)*ZP(gI2,0) + ZP(gI1,1)*ZP(gI2,1) + 2*ZP(gI1,2
      )*ZP(gI2,2) + 2*ZP(gI1,3)*ZP(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmhhhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(ZH(gI1,0)*ZH(gI2,0) + ZH(gI1,1)*ZH(gI2,1) + 4*ZH(gI1,3
      )*ZH(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarChiCha2PL(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM2(gI2,0))*ZN1(gI1,1) + 1.4142135623730951*Conj(
      UM2(gI2,1))*ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmbarChiCha2PR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(ZN2(gI1,1))*UP2(gI2,0)) + 0.7071067811865475*g2*Conj(ZN2(
      gI1,3))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmAh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*g2*(ZA(gI2,0)*ZP(gI1,0) + ZA(gI2,1)*ZP
      (gI1,1) + 1.4142135623730951*ZA(gI2,3)*(ZP(gI1,2) - ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmHpmhh(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*g2*(ZH(gI2,0)*ZP(gI1,0) - ZH(gI2,1)*ZP(gI1,1) +
      1.4142135623730951*ZH(gI2,3)*(ZP(gI1,2) + ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSdSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1655;
   std::complex<double> tmp_1656;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1656 += Conj(ZD(gI2,j1))*ZD(gI1,j1);
   }
   tmp_1655 += tmp_1656;
   result += (0.5*Sqr(g2)) * tmp_1655;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSeSe(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1657;
   std::complex<double> tmp_1658;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1658 += Conj(ZE(gI2,j1))*ZE(gI1,j1);
   }
   tmp_1657 += tmp_1658;
   result += (0.5*Sqr(g2)) * tmp_1657;

   return result;
}

std::complex<double> CLASSNAME::CpVWmconjVWmconjSuSu(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1659;
   std::complex<double> tmp_1660;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1660 += Conj(ZU(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1659 += tmp_1660;
   result += (0.5*Sqr(g2)) * tmp_1659;

   return result;
}

std::complex<double> CLASSNAME::CpconjVWmconjSuSd(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1661;
   std::complex<double> tmp_1662;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1662 += Conj(ZD(gI2,j1))*ZU(gI1,j1);
   }
   tmp_1661 += tmp_1662;
   result += (0.7071067811865475*g2) * tmp_1661;

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

std::complex<double> CLASSNAME::CpconjVWmVWmhh(unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Sqr(g2)*(vd*ZH(gI2,0) + vu*ZH(gI2,1) + 4*vT*ZH(gI2,3));

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

std::complex<double> CLASSNAME::CpbarUChibarCha2RpmPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Conj(ZRP(gI2,0))*(Conj(UP2(gI1,1))*(2*LamSU*KroneckerDelta(0,
      gO2) + 1.4142135623730951*LamTU*KroneckerDelta(1,gO2)) + 2*LamTU*Conj(UP2(
      gI1,0))*KroneckerDelta(3,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarCha2RpmPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*Conj(ZRP(gI2,0))*(-10*g2*KroneckerDelta(3,gO1)*UM2(gI1,0) +
      1.4142135623730951*(3.872983346207417*g1*KroneckerDelta(0,gO1) + 5*g2*
      KroneckerDelta(1,gO1))*UM2(gI1,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjRpmbarCha1PL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*(Conj(UM1(gI2,1))*(2*LamSD*KroneckerDelta(0,gO2) -
      1.4142135623730951*LamTD*KroneckerDelta(1,gO2)) + 2*LamTD*Conj(UM1(gI2,0))*
      KroneckerDelta(2,gO2))*ZRP(gI1,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjRpmbarCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.1*(10*g2*KroneckerDelta(2,gO1)*UP1(gI2,0) + 1.4142135623730951*(
      3.872983346207417*g1*KroneckerDelta(0,gO1) + 5*g2*KroneckerDelta(1,gO1))*UP1
      (gI2,1))*ZRP(gI1,1);

   return result;
}

double CLASSNAME::CpbarUChibarFvSvPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFvSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZV(gI2,gI1))*KroneckerDelta(0,gO1
         );
   }
   if (gI1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZV(gI2,gI1))*KroneckerDelta(1,
         gO1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFdSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1663;
   std::complex<double> tmp_1664;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1665;
      std::complex<double> tmp_1666;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1666 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1665 += tmp_1666;
      tmp_1664 += (Conj(ZD(gI2,j2))) * tmp_1665;
   }
   tmp_1663 += tmp_1664;
   result += (-KroneckerDelta(2,gO2)) * tmp_1663;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFdSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1667;
   std::complex<double> tmp_1668;
   std::complex<double> tmp_1669;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1669 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1668 += tmp_1669;
   tmp_1667 += (std::complex<double>(0,-0.18257418583505536)*g1*KroneckerDelta(
      0,gO1)) * tmp_1668;
   std::complex<double> tmp_1670;
   std::complex<double> tmp_1671;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1671 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1670 += tmp_1671;
   tmp_1667 += (std::complex<double>(0,0.7071067811865475)*g2*KroneckerDelta(1,
      gO1)) * tmp_1670;
   result += (std::complex<double>(0,-1)) * tmp_1667;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFeSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1672;
   std::complex<double> tmp_1673;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1674;
      std::complex<double> tmp_1675;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1675 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1674 += tmp_1675;
      tmp_1673 += (Conj(ZE(gI2,j2))) * tmp_1674;
   }
   tmp_1672 += tmp_1673;
   result += (-KroneckerDelta(2,gO2)) * tmp_1672;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFeSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1676;
   std::complex<double> tmp_1677;
   std::complex<double> tmp_1678;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1678 += Conj(ZE(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1677 += tmp_1678;
   tmp_1676 += (std::complex<double>(0,0.5477225575051661)*g1*KroneckerDelta(0,
      gO1)) * tmp_1677;
   std::complex<double> tmp_1679;
   std::complex<double> tmp_1680;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1680 += Conj(ZE(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1679 += tmp_1680;
   tmp_1676 += (std::complex<double>(0,0.7071067811865475)*g2*KroneckerDelta(1,
      gO1)) * tmp_1679;
   result += (std::complex<double>(0,-1)) * tmp_1676;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFuSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1681;
   std::complex<double> tmp_1682;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1683;
      std::complex<double> tmp_1684;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1684 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1683 += tmp_1684;
      tmp_1682 += (Conj(ZU(gI2,j2))) * tmp_1683;
   }
   tmp_1681 += tmp_1682;
   result += (-KroneckerDelta(3,gO2)) * tmp_1681;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarFuSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1685;
   std::complex<double> tmp_1686;
   std::complex<double> tmp_1687;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1687 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1686 += tmp_1687;
   tmp_1685 += (std::complex<double>(0,-0.18257418583505536)*g1*KroneckerDelta(
      0,gO1)) * tmp_1686;
   std::complex<double> tmp_1688;
   std::complex<double> tmp_1689;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1689 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1688 += tmp_1689;
   tmp_1685 += (std::complex<double>(0,-0.7071067811865475)*g2*KroneckerDelta(1
      ,gO1)) * tmp_1688;
   result += (std::complex<double>(0,-1)) * tmp_1685;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarChiRhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.7071067811865475*LamTD*Conj(ZN2(gI1,1))*KroneckerDelta(2,gO2)*ZHR
      (gI2,0) + Conj(ZN2(gI1,2))*(LamSD*KroneckerDelta(0,gO2)*ZHR(gI2,0) +
      0.7071067811865475*LamTD*KroneckerDelta(1,gO2)*ZHR(gI2,0)) - LamSU*Conj(ZN2(
      gI1,3))*KroneckerDelta(0,gO2)*ZHR(gI2,1) + 0.7071067811865475*LamTU*Conj(ZN2
      (gI1,3))*KroneckerDelta(1,gO2)*ZHR(gI2,1) + 0.7071067811865475*LamTU*Conj(
      ZN2(gI1,1))*KroneckerDelta(3,gO2)*ZHR(gI2,1) + Conj(ZN2(gI1,0))*(LamSD*
      KroneckerDelta(2,gO2)*ZHR(gI2,0) - LamSU*KroneckerDelta(3,gO2)*ZHR(gI2,1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChibarChiRhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1414213562373095*(KroneckerDelta(3,gO1)*ZHR(gI2,1)*(
      3.872983346207417*g1*ZN1(gI1,0) - 5*g2*ZN1(gI1,1)) + KroneckerDelta(2,gO1)*
      ZHR(gI2,0)*(-3.872983346207417*g1*ZN1(gI1,0) + 5*g2*ZN1(gI1,1)) - (
      3.872983346207417*g1*KroneckerDelta(0,gO1) - 5*g2*KroneckerDelta(1,gO1))*(
      ZHR(gI2,0)*ZN1(gI1,2) - ZHR(gI2,1)*ZN1(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjHpmCha2PL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*Conj(UM2(gI2,0))*(KroneckerDelta(2,gO2)*ZP(gI1,0) +
      1.4142135623730951*KroneckerDelta(1,gO2)*ZP(gI1,2))) + 0.5*Conj(UM2(gI2,1))*
      (2*LamSU*KroneckerDelta(0,gO2)*ZP(gI1,1) + 1.4142135623730951*LamTU*
      KroneckerDelta(1,gO2)*ZP(gI1,1) + 2*LamTU*KroneckerDelta(3,gO2)*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjHpmCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(LamTD)*KroneckerDelta(2,gO1)*UP2(gI2,0)*ZP(gI1,0) -
      0.5477225575051661*g1*KroneckerDelta(0,gO1)*UP2(gI2,1)*ZP(gI1,1) -
      0.7071067811865475*g2*KroneckerDelta(1,gO1)*UP2(gI2,1)*ZP(gI1,1) - Conj(
      LamTU)*KroneckerDelta(3,gO1)*UP2(gI2,1)*ZP(gI1,2) - 1.4142135623730951*g2*
      KroneckerDelta(1,gO1)*UP2(gI2,0)*ZP(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiHpmCha1PL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = Conj(UP1(gI2,1))*(-(LamSD*KroneckerDelta(0,gO2)*ZP(gI1,0)) +
      0.7071067811865475*LamTD*KroneckerDelta(1,gO2)*ZP(gI1,0) - LamTD*
      KroneckerDelta(2,gO2)*ZP(gI1,2)) + g2*Conj(UP1(gI2,0))*(-(KroneckerDelta(3,
      gO2)*ZP(gI1,1)) + 1.4142135623730951*KroneckerDelta(1,gO2)*ZP(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiHpmCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5477225575051661*g1*KroneckerDelta(0,gO1)*UM1(gI2,1)*ZP(gI1,0) -
      Conj(LamTU)*KroneckerDelta(3,gO1)*UM1(gI2,0)*ZP(gI1,1) + 0.7071067811865475*
      g2*KroneckerDelta(1,gO1)*(UM1(gI2,1)*ZP(gI1,0) + 2*UM1(gI2,0)*ZP(gI1,2)) +
      Conj(LamTD)*KroneckerDelta(2,gO1)*UM1(gI2,1)*ZP(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChiAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.1)*(3.872983346207417*g1*Conj(ZN1(gI1,0))*
      (-(KroneckerDelta(2,gO2)*ZA(gI2,0)) + KroneckerDelta(3,gO2)*ZA(gI2,1)) + 5*
      Conj(ZN1(gI1,2))*(1.4142135623730951*LamSD*KroneckerDelta(0,gO2)*ZA(gI2,0) +
      LamTD*KroneckerDelta(1,gO2)*ZA(gI2,0) + KroneckerDelta(2,gO2)*(
      1.4142135623730951*LamSD*ZA(gI2,2) + LamTD*ZA(gI2,3))) + 5*(g2*Conj(ZN1(gI1,
      1))*(KroneckerDelta(2,gO2)*ZA(gI2,0) - KroneckerDelta(3,gO2)*ZA(gI2,1)) +
      Conj(ZN1(gI1,3))*(-1.4142135623730951*LamSU*KroneckerDelta(0,gO2)*ZA(gI2,1)
      + LamTU*KroneckerDelta(1,gO2)*ZA(gI2,1) + KroneckerDelta(3,gO2)*(
      -1.4142135623730951*LamSU*ZA(gI2,2) + LamTU*ZA(gI2,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiChiAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.1)*(5*Conj(LamTD)*KroneckerDelta(2,gO1)*
      ZA(gI2,0)*ZN2(gI1,1) + 5*Conj(LamTU)*KroneckerDelta(3,gO1)*ZA(gI2,1)*ZN2(gI1
      ,1) - 3.872983346207417*g1*KroneckerDelta(0,gO1)*ZA(gI2,0)*ZN2(gI1,2) + 5*g2
      *KroneckerDelta(1,gO1)*ZA(gI2,0)*ZN2(gI1,2) + 5*Conj(LamTD)*KroneckerDelta(2
      ,gO1)*ZA(gI2,3)*ZN2(gI1,2) + 7.0710678118654755*Conj(LamSD)*KroneckerDelta(2
      ,gO1)*(ZA(gI2,0)*ZN2(gI1,0) + ZA(gI2,2)*ZN2(gI1,2)) + 3.872983346207417*g1*
      KroneckerDelta(0,gO1)*ZA(gI2,1)*ZN2(gI1,3) - 5*g2*KroneckerDelta(1,gO1)*ZA(
      gI2,1)*ZN2(gI1,3) + 5*Conj(LamTU)*KroneckerDelta(3,gO1)*ZA(gI2,3)*ZN2(gI1,3)
      - 7.0710678118654755*Conj(LamSU)*KroneckerDelta(3,gO1)*(ZA(gI2,1)*ZN2(gI1,0
      ) + ZA(gI2,2)*ZN2(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChihhChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(3.872983346207417*g1*Conj(ZN1(gI2,0))*(KroneckerDelta(2,gO2)*
      ZH(gI1,0) - KroneckerDelta(3,gO2)*ZH(gI1,1)) + 5*Conj(ZN1(gI2,2))*(
      1.4142135623730951*LamSD*KroneckerDelta(0,gO2)*ZH(gI1,0) + LamTD*
      KroneckerDelta(1,gO2)*ZH(gI1,0) + KroneckerDelta(2,gO2)*(1.4142135623730951*
      LamSD*ZH(gI1,2) + LamTD*ZH(gI1,3))) + 5*(Conj(ZN1(gI2,1))*(-(g2*
      KroneckerDelta(2,gO2)*ZH(gI1,0)) + g2*KroneckerDelta(3,gO2)*ZH(gI1,1)) +
      Conj(ZN1(gI2,3))*(-1.4142135623730951*LamSU*KroneckerDelta(0,gO2)*ZH(gI1,1)
      + LamTU*KroneckerDelta(1,gO2)*ZH(gI1,1) + KroneckerDelta(3,gO2)*(
      -1.4142135623730951*LamSU*ZH(gI1,2) + LamTU*ZH(gI1,3)))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChihhChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(5*Conj(LamTD)*KroneckerDelta(2,gO1)*ZH(gI1,0)*ZN2(gI2,1) + 5*
      Conj(LamTU)*KroneckerDelta(3,gO1)*ZH(gI1,1)*ZN2(gI2,1) + 3.872983346207417*
      g1*KroneckerDelta(0,gO1)*ZH(gI1,0)*ZN2(gI2,2) - 5*g2*KroneckerDelta(1,gO1)*
      ZH(gI1,0)*ZN2(gI2,2) + 5*Conj(LamTD)*KroneckerDelta(2,gO1)*ZH(gI1,3)*ZN2(gI2
      ,2) + 7.0710678118654755*Conj(LamSD)*KroneckerDelta(2,gO1)*(ZH(gI1,0)*ZN2(
      gI2,0) + ZH(gI1,2)*ZN2(gI2,2)) - 3.872983346207417*g1*KroneckerDelta(0,gO1)*
      ZH(gI1,1)*ZN2(gI2,3) + 5*g2*KroneckerDelta(1,gO1)*ZH(gI1,1)*ZN2(gI2,3) + 5*
      Conj(LamTU)*KroneckerDelta(3,gO1)*ZH(gI1,3)*ZN2(gI2,3) - 7.0710678118654755*
      Conj(LamSU)*KroneckerDelta(3,gO1)*(ZH(gI1,1)*ZN2(gI2,0) + ZH(gI1,2)*ZN2(gI2,
      3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSdFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1690;
   std::complex<double> tmp_1691;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1692;
      std::complex<double> tmp_1693;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1693 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1692 += tmp_1693;
      tmp_1691 += (Conj(ZDL(gI2,j2))) * tmp_1692;
   }
   tmp_1690 += tmp_1691;
   result += (-KroneckerDelta(2,gO2)) * tmp_1690;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSdFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1694;
   std::complex<double> tmp_1695;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1695 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1694 += tmp_1695;
   result += (-0.3651483716701107*g1*KroneckerDelta(0,gO1)) * tmp_1694;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSeFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1696;
   std::complex<double> tmp_1697;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1698;
      std::complex<double> tmp_1699;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1699 += Ye(j1,j2)*ZE(gI1,3 + j1);
      }
      tmp_1698 += tmp_1699;
      tmp_1697 += (Conj(ZEL(gI2,j2))) * tmp_1698;
   }
   tmp_1696 += tmp_1697;
   result += (-KroneckerDelta(2,gO2)) * tmp_1696;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSeFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1700;
   std::complex<double> tmp_1701;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1701 += ZE(gI1,3 + j1)*ZER(gI2,j1);
   }
   tmp_1700 += tmp_1701;
   result += (-1.0954451150103321*g1*KroneckerDelta(0,gO1)) * tmp_1700;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSuFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1702;
   std::complex<double> tmp_1703;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1704;
      std::complex<double> tmp_1705;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1705 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1704 += tmp_1705;
      tmp_1703 += (Conj(ZUL(gI2,j2))) * tmp_1704;
   }
   tmp_1702 += tmp_1703;
   result += (-KroneckerDelta(3,gO2)) * tmp_1702;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjSuFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1706;
   std::complex<double> tmp_1707;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1707 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1706 += tmp_1707;
   result += (0.7302967433402214*g1*KroneckerDelta(0,gO1)) * tmp_1706;

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiVWmCha1PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(1,gO2)*UM1(gI2,0) + 0.7071067811865475*g2*
      KroneckerDelta(2,gO2)*UM1(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiVWmCha1PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UP1(gI2,0))*KroneckerDelta(1,gO1) - 0.7071067811865475*g2*
      Conj(UP1(gI2,1))*KroneckerDelta(2,gO1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjVWmCha2PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(1,gO2)*UP2(gI2,0)) + 0.7071067811865475*g2*
      KroneckerDelta(3,gO2)*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiconjVWmCha2PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(UM2(gI2,0))*KroneckerDelta(1,gO1) +
      1.4142135623730951*Conj(UM2(gI2,1))*KroneckerDelta(3,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiVZChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()))*(
      KroneckerDelta(2,gO2)*ZN2(gI2,2) - KroneckerDelta(3,gO2)*ZN2(gI2,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUChiVZChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*(Conj(ZN1(gI2,2))*KroneckerDelta(2,gO1) - Conj(ZN1(gI2,3))*
      KroneckerDelta(3,gO1))*(5*g2*Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW
      ()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barCha2RhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1708;
   std::complex<double> tmp_1709;
   std::complex<double> tmp_1710;
   for (unsigned gl429 = 0; gl429 < 2; ++gl429) {
      tmp_1710 += Conj(UM1(gl429,1))*UP1(gl429,gO2);
   }
   tmp_1709 += tmp_1710;
   tmp_1708 += (std::complex<double>(0,1)*LamTD*Conj(UP2(gI1,0))*ZHR(gI2,0)) *
      tmp_1709;
   std::complex<double> tmp_1711;
   std::complex<double> tmp_1712;
   for (unsigned gl429 = 0; gl429 < 2; ++gl429) {
      tmp_1712 += Conj(UM1(gl429,0))*UP1(gl429,gO2);
   }
   tmp_1711 += tmp_1712;
   tmp_1708 += (std::complex<double>(0,-1)*LamTU*Conj(UP2(gI1,1))*ZHR(gI2,1)) *
      tmp_1711;
   result += (std::complex<double>(0,-1)) * tmp_1708;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barCha2RhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1713;
   std::complex<double> tmp_1714;
   std::complex<double> tmp_1715;
   for (unsigned gl430 = 0; gl430 < 2; ++gl430) {
      tmp_1715 += Conj(UM1(gl430,gO1))*UP1(gl430,1);
   }
   tmp_1714 += tmp_1715;
   tmp_1713 += (std::complex<double>(0,-1)*g2*UM2(gI1,0)*ZHR(gI2,0)) * tmp_1714
      ;
   std::complex<double> tmp_1716;
   std::complex<double> tmp_1717;
   for (unsigned gl430 = 0; gl430 < 2; ++gl430) {
      tmp_1717 += Conj(UM1(gl430,gO1))*UP1(gl430,0);
   }
   tmp_1716 += tmp_1717;
   tmp_1713 += (std::complex<double>(0,-1)*g2*UM2(gI1,1)*ZHR(gI2,1)) * tmp_1716
      ;
   result += (std::complex<double>(0,-1)) * tmp_1713;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1AhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1718;
   std::complex<double> tmp_1719;
   std::complex<double> tmp_1720;
   for (unsigned gl431 = 0; gl431 < 2; ++gl431) {
      tmp_1720 += Conj(UM1(gl431,0))*UP1(gl431,gO2);
   }
   tmp_1719 += tmp_1720;
   tmp_1718 += (0.7071067811865475*LamTD*Conj(UP1(gI1,1))*ZA(gI2,0)) * tmp_1719
      ;
   std::complex<double> tmp_1721;
   std::complex<double> tmp_1722;
   for (unsigned gl431 = 0; gl431 < 2; ++gl431) {
      tmp_1722 += Conj(UM1(gl431,1))*UP1(gl431,gO2);
   }
   tmp_1721 += tmp_1722;
   tmp_1718 += (-0.7071067811865475*g2*Conj(UP1(gI1,0))*ZA(gI2,0)) * tmp_1721;
   std::complex<double> tmp_1723;
   std::complex<double> tmp_1724;
   for (unsigned gl431 = 0; gl431 < 2; ++gl431) {
      tmp_1724 += Conj(UM1(gl431,1))*UP1(gl431,gO2);
   }
   tmp_1723 += tmp_1724;
   tmp_1718 += (0.7071067811865475*LamSD*Conj(UP1(gI1,1))*ZA(gI2,2)) * tmp_1723
      ;
   std::complex<double> tmp_1725;
   std::complex<double> tmp_1726;
   for (unsigned gl431 = 0; gl431 < 2; ++gl431) {
      tmp_1726 += Conj(UM1(gl431,0))*UP1(gl431,gO2);
   }
   tmp_1725 += tmp_1726;
   tmp_1718 += (-(g2*Conj(UP1(gI1,0))*ZA(gI2,3))) * tmp_1725;
   std::complex<double> tmp_1727;
   std::complex<double> tmp_1728;
   for (unsigned gl431 = 0; gl431 < 2; ++gl431) {
      tmp_1728 += Conj(UM1(gl431,1))*UP1(gl431,gO2);
   }
   tmp_1727 += tmp_1728;
   tmp_1718 += (-0.5*LamTD*Conj(UP1(gI1,1))*ZA(gI2,3)) * tmp_1727;
   result += (std::complex<double>(0,-1)) * tmp_1718;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1Cha1AhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1729;
   std::complex<double> tmp_1730;
   std::complex<double> tmp_1731;
   for (unsigned gl432 = 0; gl432 < 2; ++gl432) {
      tmp_1731 += Conj(UM1(gl432,gO1))*UP1(gl432,1);
   }
   tmp_1730 += tmp_1731;
   tmp_1729 += (-0.7071067811865475*Conj(LamTD)*UM1(gI1,0)*ZA(gI2,0)) *
      tmp_1730;
   std::complex<double> tmp_1732;
   std::complex<double> tmp_1733;
   for (unsigned gl432 = 0; gl432 < 2; ++gl432) {
      tmp_1733 += Conj(UM1(gl432,gO1))*UP1(gl432,0);
   }
   tmp_1732 += tmp_1733;
   tmp_1729 += (0.7071067811865475*g2*UM1(gI1,1)*ZA(gI2,0)) * tmp_1732;
   std::complex<double> tmp_1734;
   std::complex<double> tmp_1735;
   for (unsigned gl432 = 0; gl432 < 2; ++gl432) {
      tmp_1735 += Conj(UM1(gl432,gO1))*UP1(gl432,1);
   }
   tmp_1734 += tmp_1735;
   tmp_1729 += (-0.7071067811865475*Conj(LamSD)*UM1(gI1,1)*ZA(gI2,2)) *
      tmp_1734;
   std::complex<double> tmp_1736;
   std::complex<double> tmp_1737;
   for (unsigned gl432 = 0; gl432 < 2; ++gl432) {
      tmp_1737 += Conj(UM1(gl432,gO1))*UP1(gl432,0);
   }
   tmp_1736 += tmp_1737;
   tmp_1729 += (g2*UM1(gI1,0)*ZA(gI2,3)) * tmp_1736;
   std::complex<double> tmp_1738;
   std::complex<double> tmp_1739;
   for (unsigned gl432 = 0; gl432 < 2; ++gl432) {
      tmp_1739 += Conj(UM1(gl432,gO1))*UP1(gl432,1);
   }
   tmp_1738 += tmp_1739;
   tmp_1729 += (0.5*Conj(LamTD)*UM1(gI1,1)*ZA(gI2,3)) * tmp_1738;
   result += (std::complex<double>(0,-1)) * tmp_1729;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjRpmbarChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1740;
   std::complex<double> tmp_1741;
   std::complex<double> tmp_1742;
   for (unsigned gl433 = 0; gl433 < 2; ++gl433) {
      tmp_1742 += Conj(UM1(gl433,0))*UP1(gl433,gO2);
   }
   tmp_1741 += tmp_1742;
   tmp_1740 += (std::complex<double>(0,-1)*LamTD*Conj(ZN2(gI2,2))*ZRP(gI1,1)) *
      tmp_1741;
   std::complex<double> tmp_1743;
   std::complex<double> tmp_1744;
   for (unsigned gl433 = 0; gl433 < 2; ++gl433) {
      tmp_1744 += Conj(UM1(gl433,1))*UP1(gl433,gO2);
   }
   tmp_1743 += tmp_1744;
   tmp_1740 += (std::complex<double>(0,-1)*LamSD*Conj(ZN2(gI2,0))*ZRP(gI1,1)) *
      tmp_1743;
   std::complex<double> tmp_1745;
   std::complex<double> tmp_1746;
   for (unsigned gl433 = 0; gl433 < 2; ++gl433) {
      tmp_1746 += Conj(UM1(gl433,1))*UP1(gl433,gO2);
   }
   tmp_1745 += tmp_1746;
   tmp_1740 += (std::complex<double>(0,0.7071067811865475)*LamTD*Conj(ZN2(gI2,1
      ))*ZRP(gI1,1)) * tmp_1745;
   result += (std::complex<double>(0,-1)) * tmp_1740;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjRpmbarChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1747;
   std::complex<double> tmp_1748;
   std::complex<double> tmp_1749;
   for (unsigned gl434 = 0; gl434 < 2; ++gl434) {
      tmp_1749 += Conj(UM1(gl434,gO1))*UP1(gl434,1);
   }
   tmp_1748 += tmp_1749;
   tmp_1747 += (std::complex<double>(0,-0.5477225575051661)*g1*ZN1(gI2,0)*ZRP(
      gI1,1)) * tmp_1748;
   std::complex<double> tmp_1750;
   std::complex<double> tmp_1751;
   for (unsigned gl434 = 0; gl434 < 2; ++gl434) {
      tmp_1751 += Conj(UM1(gl434,gO1))*UP1(gl434,1);
   }
   tmp_1750 += tmp_1751;
   tmp_1747 += (std::complex<double>(0,-0.7071067811865475)*g2*ZN1(gI2,1)*ZRP(
      gI1,1)) * tmp_1750;
   std::complex<double> tmp_1752;
   std::complex<double> tmp_1753;
   for (unsigned gl434 = 0; gl434 < 2; ++gl434) {
      tmp_1753 += Conj(UM1(gl434,gO1))*UP1(gl434,0);
   }
   tmp_1752 += tmp_1753;
   tmp_1747 += (std::complex<double>(0,-1)*g2*ZN1(gI2,2)*ZRP(gI1,1)) * tmp_1752
      ;
   result += (std::complex<double>(0,-1)) * tmp_1747;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFeSvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1754;
   std::complex<double> tmp_1756;
   for (unsigned gl435 = 0; gl435 < 2; ++gl435) {
      tmp_1756 += Conj(UM1(gl435,1))*UP1(gl435,gO2);
   }
   tmp_1754 += tmp_1756;
   std::complex<double> tmp_1755;
   std::complex<double> tmp_1757;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1758;
      std::complex<double> tmp_1759;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1759 += Conj(ZER(gI1,j1))*Ye(j1,j2);
      }
      tmp_1758 += tmp_1759;
      tmp_1757 += (Conj(ZV(gI2,j2))) * tmp_1758;
   }
   tmp_1755 += tmp_1757;
   result += (1) * tmp_1754 * tmp_1755;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFeSvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1760;
   std::complex<double> tmp_1762;
   for (unsigned gl436 = 0; gl436 < 2; ++gl436) {
      tmp_1762 += Conj(UM1(gl436,gO1))*UP1(gl436,0);
   }
   tmp_1760 += tmp_1762;
   std::complex<double> tmp_1761;
   std::complex<double> tmp_1763;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1763 += Conj(ZV(gI2,j1))*ZEL(gI1,j1);
   }
   tmp_1761 += tmp_1763;
   result += (-g2) * tmp_1760 * tmp_1761;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFdSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1764;
   std::complex<double> tmp_1766;
   for (unsigned gl437 = 0; gl437 < 2; ++gl437) {
      tmp_1766 += Conj(UM1(gl437,1))*UP1(gl437,gO2);
   }
   tmp_1764 += tmp_1766;
   std::complex<double> tmp_1765;
   std::complex<double> tmp_1767;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1768;
      std::complex<double> tmp_1769;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1769 += Conj(ZDR(gI1,j1))*Yd(j1,j2);
      }
      tmp_1768 += tmp_1769;
      tmp_1767 += (Conj(ZU(gI2,j2))) * tmp_1768;
   }
   tmp_1765 += tmp_1767;
   result += (1) * tmp_1764 * tmp_1765;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1barFdSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1770;
   std::complex<double> tmp_1772;
   for (unsigned gl438 = 0; gl438 < 2; ++gl438) {
      tmp_1772 += Conj(UM1(gl438,gO1))*UP1(gl438,0);
   }
   tmp_1770 += tmp_1772;
   std::complex<double> tmp_1771;
   std::complex<double> tmp_1773;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1773 += Conj(ZU(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1771 += tmp_1773;
   result += (-g2) * tmp_1770 * tmp_1771;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1hhCha1PL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1774;
   std::complex<double> tmp_1775;
   std::complex<double> tmp_1776;
   for (unsigned gl439 = 0; gl439 < 2; ++gl439) {
      tmp_1776 += Conj(UM1(gl439,0))*UP1(gl439,gO2);
   }
   tmp_1775 += tmp_1776;
   tmp_1774 += (std::complex<double>(0,-0.7071067811865475)*LamTD*Conj(UP1(gI2,
      1))*ZH(gI1,0)) * tmp_1775;
   std::complex<double> tmp_1777;
   std::complex<double> tmp_1778;
   for (unsigned gl439 = 0; gl439 < 2; ++gl439) {
      tmp_1778 += Conj(UM1(gl439,1))*UP1(gl439,gO2);
   }
   tmp_1777 += tmp_1778;
   tmp_1774 += (std::complex<double>(0,-0.7071067811865475)*g2*Conj(UP1(gI2,0))
      *ZH(gI1,0)) * tmp_1777;
   std::complex<double> tmp_1779;
   std::complex<double> tmp_1780;
   for (unsigned gl439 = 0; gl439 < 2; ++gl439) {
      tmp_1780 += Conj(UM1(gl439,1))*UP1(gl439,gO2);
   }
   tmp_1779 += tmp_1780;
   tmp_1774 += (std::complex<double>(0,-0.7071067811865475)*LamSD*Conj(UP1(gI2,
      1))*ZH(gI1,2)) * tmp_1779;
   std::complex<double> tmp_1781;
   std::complex<double> tmp_1782;
   for (unsigned gl439 = 0; gl439 < 2; ++gl439) {
      tmp_1782 += Conj(UM1(gl439,0))*UP1(gl439,gO2);
   }
   tmp_1781 += tmp_1782;
   tmp_1774 += (std::complex<double>(0,-1)*g2*Conj(UP1(gI2,0))*ZH(gI1,3)) *
      tmp_1781;
   std::complex<double> tmp_1783;
   std::complex<double> tmp_1784;
   for (unsigned gl439 = 0; gl439 < 2; ++gl439) {
      tmp_1784 += Conj(UM1(gl439,1))*UP1(gl439,gO2);
   }
   tmp_1783 += tmp_1784;
   tmp_1774 += (std::complex<double>(0,0.5)*LamTD*Conj(UP1(gI2,1))*ZH(gI1,3)) *
      tmp_1783;
   result += (std::complex<double>(0,-1)) * tmp_1774;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1hhCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1785;
   std::complex<double> tmp_1786;
   std::complex<double> tmp_1787;
   for (unsigned gl440 = 0; gl440 < 2; ++gl440) {
      tmp_1787 += Conj(UM1(gl440,gO1))*UP1(gl440,1);
   }
   tmp_1786 += tmp_1787;
   tmp_1785 += (std::complex<double>(0,-0.7071067811865475)*Conj(LamTD)*UM1(gI2
      ,0)*ZH(gI1,0)) * tmp_1786;
   std::complex<double> tmp_1788;
   std::complex<double> tmp_1789;
   for (unsigned gl440 = 0; gl440 < 2; ++gl440) {
      tmp_1789 += Conj(UM1(gl440,gO1))*UP1(gl440,0);
   }
   tmp_1788 += tmp_1789;
   tmp_1785 += (std::complex<double>(0,-0.7071067811865475)*g2*UM1(gI2,1)*ZH(
      gI1,0)) * tmp_1788;
   std::complex<double> tmp_1790;
   std::complex<double> tmp_1791;
   for (unsigned gl440 = 0; gl440 < 2; ++gl440) {
      tmp_1791 += Conj(UM1(gl440,gO1))*UP1(gl440,1);
   }
   tmp_1790 += tmp_1791;
   tmp_1785 += (std::complex<double>(0,-0.7071067811865475)*Conj(LamSD)*UM1(gI2
      ,1)*ZH(gI1,2)) * tmp_1790;
   std::complex<double> tmp_1792;
   std::complex<double> tmp_1793;
   for (unsigned gl440 = 0; gl440 < 2; ++gl440) {
      tmp_1793 += Conj(UM1(gl440,gO1))*UP1(gl440,0);
   }
   tmp_1792 += tmp_1793;
   tmp_1785 += (std::complex<double>(0,-1)*g2*UM1(gI2,0)*ZH(gI1,3)) * tmp_1792;
   std::complex<double> tmp_1794;
   std::complex<double> tmp_1795;
   for (unsigned gl440 = 0; gl440 < 2; ++gl440) {
      tmp_1795 += Conj(UM1(gl440,gO1))*UP1(gl440,1);
   }
   tmp_1794 += tmp_1795;
   tmp_1785 += (std::complex<double>(0,0.5)*Conj(LamTD)*UM1(gI2,1)*ZH(gI1,3)) *
      tmp_1794;
   result += (std::complex<double>(0,-1)) * tmp_1785;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjHpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1796;
   std::complex<double> tmp_1797;
   std::complex<double> tmp_1798;
   for (unsigned gl441 = 0; gl441 < 2; ++gl441) {
      tmp_1798 += Conj(UM1(gl441,1))*UP1(gl441,gO2);
   }
   tmp_1797 += tmp_1798;
   tmp_1796 += (std::complex<double>(0,0.5477225575051661)*g1*Conj(ZN1(gI2,0))*
      ZP(gI1,0)) * tmp_1797;
   std::complex<double> tmp_1799;
   std::complex<double> tmp_1800;
   for (unsigned gl441 = 0; gl441 < 2; ++gl441) {
      tmp_1800 += Conj(UM1(gl441,1))*UP1(gl441,gO2);
   }
   tmp_1799 += tmp_1800;
   tmp_1796 += (std::complex<double>(0,0.7071067811865475)*g2*Conj(ZN1(gI2,1))*
      ZP(gI1,0)) * tmp_1799;
   std::complex<double> tmp_1801;
   std::complex<double> tmp_1802;
   for (unsigned gl441 = 0; gl441 < 2; ++gl441) {
      tmp_1802 += Conj(UM1(gl441,0))*UP1(gl441,gO2);
   }
   tmp_1801 += tmp_1802;
   tmp_1796 += (std::complex<double>(0,-1)*LamTU*Conj(ZN1(gI2,3))*ZP(gI1,1)) *
      tmp_1801;
   std::complex<double> tmp_1803;
   std::complex<double> tmp_1804;
   for (unsigned gl441 = 0; gl441 < 2; ++gl441) {
      tmp_1804 += Conj(UM1(gl441,0))*UP1(gl441,gO2);
   }
   tmp_1803 += tmp_1804;
   tmp_1796 += (std::complex<double>(0,1.4142135623730951)*g2*Conj(ZN1(gI2,1))*
      ZP(gI1,2)) * tmp_1803;
   std::complex<double> tmp_1805;
   std::complex<double> tmp_1806;
   for (unsigned gl441 = 0; gl441 < 2; ++gl441) {
      tmp_1806 += Conj(UM1(gl441,1))*UP1(gl441,gO2);
   }
   tmp_1805 += tmp_1806;
   tmp_1796 += (std::complex<double>(0,1)*LamTD*Conj(ZN1(gI2,2))*ZP(gI1,3)) *
      tmp_1805;
   result += (std::complex<double>(0,-1)) * tmp_1796;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjHpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1807;
   std::complex<double> tmp_1808;
   std::complex<double> tmp_1809;
   for (unsigned gl442 = 0; gl442 < 2; ++gl442) {
      tmp_1809 += Conj(UM1(gl442,gO1))*UP1(gl442,1);
   }
   tmp_1808 += tmp_1809;
   tmp_1807 += (std::complex<double>(0,-1)*Conj(LamSD)*ZN2(gI2,0)*ZP(gI1,0)) *
      tmp_1808;
   std::complex<double> tmp_1810;
   std::complex<double> tmp_1811;
   for (unsigned gl442 = 0; gl442 < 2; ++gl442) {
      tmp_1811 += Conj(UM1(gl442,gO1))*UP1(gl442,1);
   }
   tmp_1810 += tmp_1811;
   tmp_1807 += (std::complex<double>(0,0.7071067811865475)*Conj(LamTD)*ZN2(gI2,
      1)*ZP(gI1,0)) * tmp_1810;
   std::complex<double> tmp_1812;
   std::complex<double> tmp_1813;
   for (unsigned gl442 = 0; gl442 < 2; ++gl442) {
      tmp_1813 += Conj(UM1(gl442,gO1))*UP1(gl442,0);
   }
   tmp_1812 += tmp_1813;
   tmp_1807 += (std::complex<double>(0,-1)*g2*ZN2(gI2,3)*ZP(gI1,1)) * tmp_1812;
   std::complex<double> tmp_1814;
   std::complex<double> tmp_1815;
   for (unsigned gl442 = 0; gl442 < 2; ++gl442) {
      tmp_1815 += Conj(UM1(gl442,gO1))*UP1(gl442,1);
   }
   tmp_1814 += tmp_1815;
   tmp_1807 += (std::complex<double>(0,-1)*Conj(LamTD)*ZN2(gI2,2)*ZP(gI1,2)) *
      tmp_1814;
   std::complex<double> tmp_1816;
   std::complex<double> tmp_1817;
   for (unsigned gl442 = 0; gl442 < 2; ++gl442) {
      tmp_1817 += Conj(UM1(gl442,gO1))*UP1(gl442,0);
   }
   tmp_1816 += tmp_1817;
   tmp_1807 += (std::complex<double>(0,1.4142135623730951)*g2*ZN2(gI2,1)*ZP(gI1
      ,3)) * tmp_1816;
   result += (std::complex<double>(0,-1)) * tmp_1807;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjSdFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1818;
   std::complex<double> tmp_1820;
   for (unsigned gl443 = 0; gl443 < 2; ++gl443) {
      tmp_1820 += Conj(UM1(gl443,1))*UP1(gl443,gO2);
   }
   tmp_1818 += tmp_1820;
   std::complex<double> tmp_1819;
   std::complex<double> tmp_1821;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1822;
      std::complex<double> tmp_1823;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1823 += Yd(j1,j2)*ZD(gI1,3 + j1);
      }
      tmp_1822 += tmp_1823;
      tmp_1821 += (Conj(ZUL(gI2,j2))) * tmp_1822;
   }
   tmp_1819 += tmp_1821;
   result += (1) * tmp_1818 * tmp_1819;

   return result;
}

double CLASSNAME::CpbarUCha1conjSdFuPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjSeFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1824;
   std::complex<double> tmp_1826;
   for (unsigned gl445 = 0; gl445 < 2; ++gl445) {
      tmp_1826 += Conj(UM1(gl445,1))*UP1(gl445,gO2);
   }
   tmp_1824 += tmp_1826;
   std::complex<double> tmp_1825;
   std::complex<double> tmp_1827;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1827 += Ye(j1,gI2)*ZE(gI1,3 + j1);
   }
   tmp_1825 += tmp_1827;
   result += (1) * tmp_1824 * tmp_1825;

   return result;
}

double CLASSNAME::CpbarUCha1conjSeFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VPCha1PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1828;
   std::complex<double> tmp_1829;
   std::complex<double> tmp_1830;
   for (unsigned gl447 = 0; gl447 < 2; ++gl447) {
      tmp_1830 += Conj(UM1(gl447,0))*UP1(gl447,gO2);
   }
   tmp_1829 += tmp_1830;
   tmp_1828 += (std::complex<double>(0,-1)*g2*Sin(ThetaW())*UM1(gI2,0)) *
      tmp_1829;
   std::complex<double> tmp_1831;
   std::complex<double> tmp_1832;
   for (unsigned gl447 = 0; gl447 < 2; ++gl447) {
      tmp_1832 += Conj(UM1(gl447,1))*UP1(gl447,gO2);
   }
   tmp_1831 += tmp_1832;
   tmp_1828 += (std::complex<double>(0,-0.3872983346207417)*g1*Cos(ThetaW())*
      UM1(gI2,1)) * tmp_1831;
   std::complex<double> tmp_1833;
   std::complex<double> tmp_1834;
   for (unsigned gl447 = 0; gl447 < 2; ++gl447) {
      tmp_1834 += Conj(UM1(gl447,1))*UP1(gl447,gO2);
   }
   tmp_1833 += tmp_1834;
   tmp_1828 += (std::complex<double>(0,-0.5)*g2*Sin(ThetaW())*UM1(gI2,1)) *
      tmp_1833;
   result += (std::complex<double>(0,-1)) * tmp_1828;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VPCha1PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1835;
   std::complex<double> tmp_1836;
   std::complex<double> tmp_1837;
   for (unsigned gl448 = 0; gl448 < 2; ++gl448) {
      tmp_1837 += Conj(UM1(gl448,gO1))*UP1(gl448,1);
   }
   tmp_1836 += tmp_1837;
   tmp_1835 += (std::complex<double>(0,-0.3872983346207417)*g1*Conj(UP1(gI2,1))
      *Cos(ThetaW())) * tmp_1836;
   std::complex<double> tmp_1838;
   std::complex<double> tmp_1839;
   for (unsigned gl448 = 0; gl448 < 2; ++gl448) {
      tmp_1839 += Conj(UM1(gl448,gO1))*UP1(gl448,0);
   }
   tmp_1838 += tmp_1839;
   tmp_1835 += (std::complex<double>(0,-1)*g2*Conj(UP1(gI2,0))*Sin(ThetaW())) *
      tmp_1838;
   std::complex<double> tmp_1840;
   std::complex<double> tmp_1841;
   for (unsigned gl448 = 0; gl448 < 2; ++gl448) {
      tmp_1841 += Conj(UM1(gl448,gO1))*UP1(gl448,1);
   }
   tmp_1840 += tmp_1841;
   tmp_1835 += (std::complex<double>(0,-0.5)*g2*Conj(UP1(gI2,1))*Sin(ThetaW()))
      * tmp_1840;
   result += (std::complex<double>(0,-1)) * tmp_1835;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VZCha1PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1842;
   std::complex<double> tmp_1843;
   std::complex<double> tmp_1844;
   for (unsigned gl449 = 0; gl449 < 2; ++gl449) {
      tmp_1844 += Conj(UM1(gl449,0))*UP1(gl449,gO2);
   }
   tmp_1843 += tmp_1844;
   tmp_1842 += (std::complex<double>(0,-1)*g2*Cos(ThetaW())*UM1(gI2,0)) *
      tmp_1843;
   std::complex<double> tmp_1845;
   std::complex<double> tmp_1846;
   for (unsigned gl449 = 0; gl449 < 2; ++gl449) {
      tmp_1846 += Conj(UM1(gl449,1))*UP1(gl449,gO2);
   }
   tmp_1845 += tmp_1846;
   tmp_1842 += (std::complex<double>(0,-0.5)*g2*Cos(ThetaW())*UM1(gI2,1)) *
      tmp_1845;
   std::complex<double> tmp_1847;
   std::complex<double> tmp_1848;
   for (unsigned gl449 = 0; gl449 < 2; ++gl449) {
      tmp_1848 += Conj(UM1(gl449,1))*UP1(gl449,gO2);
   }
   tmp_1847 += tmp_1848;
   tmp_1842 += (std::complex<double>(0,0.3872983346207417)*g1*Sin(ThetaW())*UM1
      (gI2,1)) * tmp_1847;
   result += (std::complex<double>(0,-1)) * tmp_1842;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1VZCha1PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1849;
   std::complex<double> tmp_1850;
   std::complex<double> tmp_1851;
   for (unsigned gl450 = 0; gl450 < 2; ++gl450) {
      tmp_1851 += Conj(UM1(gl450,gO1))*UP1(gl450,0);
   }
   tmp_1850 += tmp_1851;
   tmp_1849 += (std::complex<double>(0,-1)*g2*Conj(UP1(gI2,0))*Cos(ThetaW())) *
      tmp_1850;
   std::complex<double> tmp_1852;
   std::complex<double> tmp_1853;
   for (unsigned gl450 = 0; gl450 < 2; ++gl450) {
      tmp_1853 += Conj(UM1(gl450,gO1))*UP1(gl450,1);
   }
   tmp_1852 += tmp_1853;
   tmp_1849 += (std::complex<double>(0,-0.5)*g2*Conj(UP1(gI2,1))*Cos(ThetaW()))
      * tmp_1852;
   std::complex<double> tmp_1854;
   std::complex<double> tmp_1855;
   for (unsigned gl450 = 0; gl450 < 2; ++gl450) {
      tmp_1855 += Conj(UM1(gl450,gO1))*UP1(gl450,1);
   }
   tmp_1854 += tmp_1855;
   tmp_1849 += (std::complex<double>(0,0.3872983346207417)*g1*Conj(UP1(gI2,1))*
      Sin(ThetaW())) * tmp_1854;
   result += (std::complex<double>(0,-1)) * tmp_1849;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjVWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1856;
   std::complex<double> tmp_1857;
   std::complex<double> tmp_1858;
   for (unsigned gl451 = 0; gl451 < 2; ++gl451) {
      tmp_1858 += Conj(UM1(gl451,0))*UP1(gl451,gO2);
   }
   tmp_1857 += tmp_1858;
   tmp_1856 += (std::complex<double>(0,1)*g2*ZN2(gI2,1)) * tmp_1857;
   std::complex<double> tmp_1859;
   std::complex<double> tmp_1860;
   for (unsigned gl451 = 0; gl451 < 2; ++gl451) {
      tmp_1860 += Conj(UM1(gl451,1))*UP1(gl451,gO2);
   }
   tmp_1859 += tmp_1860;
   tmp_1856 += (std::complex<double>(0,0.7071067811865475)*g2*ZN2(gI2,2)) *
      tmp_1859;
   result += (std::complex<double>(0,-1)) * tmp_1856;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha1conjVWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1861;
   std::complex<double> tmp_1862;
   std::complex<double> tmp_1863;
   for (unsigned gl452 = 0; gl452 < 2; ++gl452) {
      tmp_1863 += Conj(UM1(gl452,gO1))*UP1(gl452,0);
   }
   tmp_1862 += tmp_1863;
   tmp_1861 += (std::complex<double>(0,1)*g2*Conj(ZN1(gI2,1))) * tmp_1862;
   std::complex<double> tmp_1864;
   std::complex<double> tmp_1865;
   for (unsigned gl452 = 0; gl452 < 2; ++gl452) {
      tmp_1865 += Conj(UM1(gl452,gO1))*UP1(gl452,1);
   }
   tmp_1864 += tmp_1865;
   tmp_1861 += (std::complex<double>(0,-0.7071067811865475)*g2*Conj(ZN1(gI2,2))
      ) * tmp_1864;
   result += (std::complex<double>(0,-1)) * tmp_1861;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barCha1RhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = LamTD*Conj(UM1(gI1,1))*KroneckerDelta(0,gO2)*ZHR(gI2,0) - LamTU*
      Conj(UM1(gI1,0))*KroneckerDelta(1,gO2)*ZHR(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barCha1RhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*(KroneckerDelta(0,gO1)*UP1(gI1,1)*ZHR(gI2,0) + KroneckerDelta(
      1,gO1)*UP1(gI1,0)*ZHR(gI2,1)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2AhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,0.5)*(g2*Conj(UM2(gI1,0))*(
      1.4142135623730951*KroneckerDelta(1,gO2)*ZA(gI2,1) - 2*KroneckerDelta(0,gO2)
      *ZA(gI2,3)) + Conj(UM2(gI1,1))*(1.4142135623730951*LamTU*KroneckerDelta(0,
      gO2)*ZA(gI2,1) + KroneckerDelta(1,gO2)*(1.4142135623730951*LamSU*ZA(gI2,2) +
      LamTU*ZA(gI2,3))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2Cha2AhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = std::complex<double>(0,-0.5)*(1.4142135623730951*Conj(LamSU)*
      KroneckerDelta(1,gO1)*UP2(gI1,1)*ZA(gI2,2) + g2*KroneckerDelta(0,gO1)*(
      1.4142135623730951*UP2(gI1,1)*ZA(gI2,1) - 2*UP2(gI1,0)*ZA(gI2,3)) + Conj(
      LamTU)*KroneckerDelta(1,gO1)*(1.4142135623730951*UP2(gI1,0)*ZA(gI2,1) + UP2(
      gI1,1)*ZA(gI2,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFuSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1866;
   std::complex<double> tmp_1867;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1868;
      std::complex<double> tmp_1869;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1869 += Conj(ZUR(gI1,j1))*Yu(j1,j2);
      }
      tmp_1868 += tmp_1869;
      tmp_1867 += (Conj(ZD(gI2,j2))) * tmp_1868;
   }
   tmp_1866 += tmp_1867;
   result += (KroneckerDelta(1,gO2)) * tmp_1866;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFuSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1870;
   std::complex<double> tmp_1871;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1871 += Conj(ZD(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1870 += tmp_1871;
   result += (-(g2*KroneckerDelta(0,gO1))) * tmp_1870;

   return result;
}

double CLASSNAME::CpbarUCha2barFvSePL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barFvSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gI1 < 3) {
      result += -(g2*Conj(ZE(gI2,gI1))*KroneckerDelta(0,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barChiRpmPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*Conj(ZRP(gI2,0))*(2*LamTU*Conj(ZN2(gI1,3))*KroneckerDelta(0,gO2
      ) + (2*LamSU*Conj(ZN2(gI1,0)) + 1.4142135623730951*LamTU*Conj(ZN2(gI1,1)))*
      KroneckerDelta(1,gO2));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2barChiRpmPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.1*Conj(ZRP(gI2,0))*(KroneckerDelta(1,gO1)*(5.477225575051661*g1*
      ZN1(gI1,0) + 7.0710678118654755*g2*ZN1(gI1,1)) - 10*g2*KroneckerDelta(0,gO1)
      *ZN1(gI1,3));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2hhCha2PL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(g2*Conj(UM2(gI2,0))*(-1.4142135623730951*KroneckerDelta(1,gO2)
      *ZH(gI1,1) + 2*KroneckerDelta(0,gO2)*ZH(gI1,3)) + Conj(UM2(gI2,1))*(
      1.4142135623730951*LamTU*KroneckerDelta(0,gO2)*ZH(gI1,1) + KroneckerDelta(1,
      gO2)*(1.4142135623730951*LamSU*ZH(gI1,2) + LamTU*ZH(gI1,3))));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2hhCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = 0.5*(1.4142135623730951*Conj(LamSU)*KroneckerDelta(1,gO1)*UP2(gI2,1
      )*ZH(gI1,2) + KroneckerDelta(0,gO1)*(-1.4142135623730951*g2*UP2(gI2,1)*ZH(
      gI1,1) + 2*g2*UP2(gI2,0)*ZH(gI1,3)) + Conj(LamTU)*KroneckerDelta(1,gO1)*(
      1.4142135623730951*UP2(gI2,0)*ZH(gI1,1) + UP2(gI2,1)*ZH(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2HpmChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = LamTD*Conj(ZN1(gI2,2))*KroneckerDelta(0,gO2)*ZP(gI1,0) -
      0.5477225575051661*g1*Conj(ZN1(gI2,0))*KroneckerDelta(1,gO2)*ZP(gI1,1) -
      0.7071067811865475*g2*Conj(ZN1(gI2,1))*KroneckerDelta(1,gO2)*ZP(gI1,1) -
      LamTU*Conj(ZN1(gI2,3))*KroneckerDelta(1,gO2)*ZP(gI1,2) - 1.4142135623730951*
      g2*Conj(ZN1(gI2,1))*KroneckerDelta(0,gO2)*ZP(gI1,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2HpmChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO1)*(ZN2(gI2,2)*ZP(gI1,0) +
      1.4142135623730951*ZN2(gI2,1)*ZP(gI1,2))) + 0.5*KroneckerDelta(1,gO1)*(2*
      Conj(LamSU)*ZN2(gI2,0)*ZP(gI1,1) + Conj(LamTU)*(1.4142135623730951*ZN2(gI2,1
      )*ZP(gI1,1) + 2*ZN2(gI2,3)*ZP(gI1,3)));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2conjSuFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1872;
   std::complex<double> tmp_1873;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1874;
      std::complex<double> tmp_1875;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1875 += Yu(j1,j2)*ZU(gI1,3 + j1);
      }
      tmp_1874 += tmp_1875;
      tmp_1873 += (Conj(ZDL(gI2,j2))) * tmp_1874;
   }
   tmp_1872 += tmp_1873;
   result += (KroneckerDelta(1,gO2)) * tmp_1872;

   return result;
}

double CLASSNAME::CpbarUCha2conjSuFdPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2VPCha2PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*KroneckerDelta(0,gO2)*Sin(ThetaW())*UP2(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW(
      )))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2VPCha2PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM2(gI2,0))*KroneckerDelta(0,gO1)*Sin(ThetaW()) + 0.1*Conj(
      UM2(gI2,1))*KroneckerDelta(1,gO1)*(3.872983346207417*g1*Cos(ThetaW()) + 5*g2
      *Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2VZCha2PR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Cos(ThetaW())*KroneckerDelta(0,gO2)*UP2(gI2,0) + 0.1*
      KroneckerDelta(1,gO2)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1*Sin(ThetaW(
      )))*UP2(gI2,1);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2VZCha2PL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = g2*Conj(UM2(gI2,0))*Cos(ThetaW())*KroneckerDelta(0,gO1) + 0.1*Conj(
      UM2(gI2,1))*KroneckerDelta(1,gO1)*(5*g2*Cos(ThetaW()) - 3.872983346207417*g1
      *Sin(ThetaW()));

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2VWmChiPR(unsigned gO2, unsigned gI2) const
{
   std::complex<double> result;

   result = -(g2*KroneckerDelta(0,gO2)*ZN2(gI2,1)) + 0.7071067811865475*g2*
      KroneckerDelta(1,gO2)*ZN2(gI2,3);

   return result;
}

std::complex<double> CLASSNAME::CpbarUCha2VWmChiPL(unsigned gO1, unsigned gI2) const
{
   std::complex<double> result;

   result = -0.5*g2*(2*Conj(ZN1(gI2,1))*KroneckerDelta(0,gO1) +
      1.4142135623730951*Conj(ZN1(gI2,3))*KroneckerDelta(1,gO1));

   return result;
}

std::complex<double> CLASSNAME::CpbarUFebarCha1SvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1876;
      std::complex<double> tmp_1877;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1877 += Conj(ZV(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1876 += tmp_1877;
      result += (Conj(UM1(gI1,1))) * tmp_1876;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFebarCha1SvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZV(gI2,gO1))*UP1(gI1,0));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1878;
      std::complex<double> tmp_1879;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1879 += Conj(ZEL(gI1,j2))*Ye(gO2,j2);
      }
      tmp_1878 += tmp_1879;
      result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1878;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1880;
      std::complex<double> tmp_1881;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1881 += Conj(Ye(j1,gO1))*ZER(gI1,j1);
      }
      tmp_1880 += tmp_1881;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1880;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1882;
      std::complex<double> tmp_1883;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1883 += Conj(ZEL(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1882 += tmp_1883;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1882;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1884;
      std::complex<double> tmp_1885;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1885 += Conj(Ye(j1,gO1))*ZER(gI2,j1);
      }
      tmp_1884 += tmp_1885;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1884;
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

std::complex<double> CLASSNAME::CpbarUFebarChiSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1886;
      std::complex<double> tmp_1887;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1887 += Conj(ZE(gI2,j2))*Ye(gO2,j2);
      }
      tmp_1886 += tmp_1887;
      result += (-Conj(ZN2(gI1,2))) * tmp_1886;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFebarChiSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += 0.5477225575051661*g1*Conj(ZE(gI2,gO1))*ZN1(gI1,0);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZE(gI2,gO1))*ZN1(gI1,1);
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -1.0954451150103321*g1*Conj(ZE(gI1,3 + gO2))*Conj(ZN1(gI2,0)
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1888;
      std::complex<double> tmp_1889;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1889 += Conj(Ye(j1,gO1))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1888 += tmp_1889;
      result += (-ZN2(gI2,2)) * tmp_1888;
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

std::complex<double> CLASSNAME::CpbarUFdbarCha1SuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1890;
      std::complex<double> tmp_1891;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1891 += Conj(ZU(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1890 += tmp_1891;
      result += (Conj(UM1(gI1,1))) * tmp_1890;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdbarCha1SuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZU(gI2,gO1))*UP1(gI1,0));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1892;
      std::complex<double> tmp_1893;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1893 += Conj(ZDL(gI1,j2))*Yd(gO2,j2);
      }
      tmp_1892 += tmp_1893;
      result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) *
         tmp_1892;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1894;
      std::complex<double> tmp_1895;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1895 += Conj(Yd(j1,gO1))*ZDR(gI1,j1);
      }
      tmp_1894 += tmp_1895;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) *
         tmp_1894;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1896;
      std::complex<double> tmp_1897;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1897 += Conj(ZDL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1896 += tmp_1897;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1896;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1898;
      std::complex<double> tmp_1899;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1899 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1898 += tmp_1899;
      result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1898;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1900;
      std::complex<double> tmp_1901;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1901 += Conj(ZUL(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1900 += tmp_1901;
      result += (ZP(gI1,0)) * tmp_1900;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1902;
      std::complex<double> tmp_1903;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1903 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1902 += tmp_1903;
      result += (ZP(gI1,1)) * tmp_1902;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdbarChiSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1904;
      std::complex<double> tmp_1905;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1905 += Conj(ZD(gI2,j2))*Yd(gO2,j2);
      }
      tmp_1904 += tmp_1905;
      result += (-Conj(ZN2(gI1,2))) * tmp_1904;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdbarChiSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZD(gI2,gO1))*ZN1(gI1,0);
   }
   if (gO1 < 3) {
      result += 0.7071067811865475*g2*Conj(ZD(gI2,gO1))*ZN1(gI1,1);
   }

   return result;
}

double CLASSNAME::CpbarUFdSuCha2PL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSuCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1906;
      std::complex<double> tmp_1907;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1907 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1906 += tmp_1907;
      result += (UP2(gI2,1)) * tmp_1906;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += -0.3651483716701107*g1*Conj(ZD(gI1,3 + gO2))*Conj(ZN1(gI2,0)
         );
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1908;
      std::complex<double> tmp_1909;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1909 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1908 += tmp_1909;
      result += (-ZN2(gI2,2)) * tmp_1908;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 1.4142135623730951*g3*Conj(ZD(gI1,3 + gO2));
   }

   return result;
}

double CLASSNAME::CpbarUFdSdGluPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

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

double CLASSNAME::CpbarUFdbarGluSdPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFdbarGluSdPR(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(ZD(gI2,gO1));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarCha2SdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1910;
      std::complex<double> tmp_1911;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1911 += Conj(ZD(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1910 += tmp_1911;
      result += (Conj(UP2(gI1,1))) * tmp_1910;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarCha2SdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -(g2*Conj(ZD(gI2,gO1))*UM2(gI1,0));
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1912;
      std::complex<double> tmp_1913;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1913 += Conj(ZUL(gI1,j2))*Yu(gO2,j2);
      }
      tmp_1912 += tmp_1913;
      result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,1)) *
         tmp_1912;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1914;
      std::complex<double> tmp_1915;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1915 += Conj(Yu(j1,gO1))*ZUR(gI1,j1);
      }
      tmp_1914 += tmp_1915;
      result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,1)) *
         tmp_1914;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1916;
      std::complex<double> tmp_1917;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1917 += Conj(ZDL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1916 += tmp_1917;
      result += (ZP(gI1,1)) * tmp_1916;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1918;
      std::complex<double> tmp_1919;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1919 += Conj(Yd(j1,gO1))*ZDR(gI2,j1);
      }
      tmp_1918 += tmp_1919;
      result += (ZP(gI1,0)) * tmp_1918;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1920;
      std::complex<double> tmp_1921;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1921 += Conj(ZUL(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1920 += tmp_1921;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1920;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1922;
      std::complex<double> tmp_1923;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1923 += Conj(Yu(j1,gO1))*ZUR(gI2,j1);
      }
      tmp_1922 += tmp_1923;
      result += (-0.7071067811865475*ZH(gI1,1)) * tmp_1922;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChiSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      std::complex<double> tmp_1924;
      std::complex<double> tmp_1925;
      for (unsigned j2 = 0; j2 < 3; ++j2) {
         tmp_1925 += Conj(ZU(gI2,j2))*Yu(gO2,j2);
      }
      tmp_1924 += tmp_1925;
      result += (-Conj(ZN2(gI1,3))) * tmp_1924;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarChiSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -0.18257418583505536*g1*Conj(ZU(gI2,gO1))*ZN1(gI1,0);
   }
   if (gO1 < 3) {
      result += -0.7071067811865475*g2*Conj(ZU(gI2,gO1))*ZN1(gI1,1);
   }

   return result;
}

double CLASSNAME::CpbarUFuSdCha1PL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSdCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1926;
      std::complex<double> tmp_1927;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1927 += Conj(Yd(j1,gO1))*Conj(ZD(gI1,3 + j1));
      }
      tmp_1926 += tmp_1927;
      result += (UM1(gI2,1)) * tmp_1926;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 0.7302967433402214*g1*Conj(ZN1(gI2,0))*Conj(ZU(gI1,3 + gO2))
         ;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      std::complex<double> tmp_1928;
      std::complex<double> tmp_1929;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1929 += Conj(Yu(j1,gO1))*Conj(ZU(gI1,3 + j1));
      }
      tmp_1928 += tmp_1929;
      result += (-ZN2(gI2,3)) * tmp_1928;
   }

   return result;
}

std::complex<double> CLASSNAME::CpbarUFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   if (gO2 < 3) {
      result += 1.4142135623730951*g3*Conj(ZU(gI1,3 + gO2));
   }

   return result;
}

double CLASSNAME::CpbarUFuSuGluPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

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

double CLASSNAME::CpbarUFubarGluSuPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarUFubarGluSuPR(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   if (gO1 < 3) {
      result += -1.4142135623730951*g3*Conj(ZU(gI2,gO1));
   }

   return result;
}

double CLASSNAME::CpbarGlubarFdSdPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGlubarFdSdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1930;
   std::complex<double> tmp_1931;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1931 += Conj(ZD(gI2,j1))*ZDL(gI1,j1);
   }
   tmp_1930 += tmp_1931;
   result += (-1.4142135623730951*g3) * tmp_1930;

   return result;
}

double CLASSNAME::CpbarGlubarFuSuPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGlubarFuSuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1932;
   std::complex<double> tmp_1933;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1933 += Conj(ZU(gI2,j1))*ZUL(gI1,j1);
   }
   tmp_1932 += tmp_1933;
   result += (-1.4142135623730951*g3) * tmp_1932;

   return result;
}

double CLASSNAME::CpbarGluconjSdFdPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluconjSdFdPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1934;
   std::complex<double> tmp_1935;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1935 += ZD(gI1,3 + j1)*ZDR(gI2,j1);
   }
   tmp_1934 += tmp_1935;
   result += (1.4142135623730951*g3) * tmp_1934;

   return result;
}

double CLASSNAME::CpbarGluconjSuFuPL(unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluconjSuFuPR(unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1936;
   std::complex<double> tmp_1937;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1937 += ZU(gI1,3 + j1)*ZUR(gI2,j1);
   }
   tmp_1936 += tmp_1937;
   result += (1.4142135623730951*g3) * tmp_1936;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluVGGluPR() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluVGGluPL() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> CLASSNAME::CpbarGluconjSOcGluPL() const
{
   std::complex<double> result;

   result = std::complex<double>(0,-1.4142135623730951)*g3;

   return result;
}

double CLASSNAME::CpbarGluconjSOcGluPR() const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarCha1SvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1938;
   std::complex<double> tmp_1939;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1940;
      std::complex<double> tmp_1941;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1941 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1940 += tmp_1941;
      tmp_1939 += (Conj(ZV(gI2,j2))) * tmp_1940;
   }
   tmp_1938 += tmp_1939;
   result += (Conj(UM1(gI1,1))) * tmp_1938;

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarCha1SvPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1942;
   std::complex<double> tmp_1943;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1943 += Conj(ZV(gI2,j1))*ZEL(gO1,j1);
   }
   tmp_1942 += tmp_1943;
   result += (-(g2*UP1(gI1,0))) * tmp_1942;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1944;
   std::complex<double> tmp_1945;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1946;
      std::complex<double> tmp_1947;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1947 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1946 += tmp_1947;
      tmp_1945 += (Conj(ZEL(gI1,j2))) * tmp_1946;
   }
   tmp_1944 += tmp_1945;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) * tmp_1944
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeFeAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1948;
   std::complex<double> tmp_1949;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1950;
      std::complex<double> tmp_1951;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1951 += Conj(Ye(j1,j2))*ZER(gI1,j1);
      }
      tmp_1950 += tmp_1951;
      tmp_1949 += (ZEL(gO1,j2)) * tmp_1950;
   }
   tmp_1948 += tmp_1949;
   result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) * tmp_1948;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1952;
   std::complex<double> tmp_1953;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1954;
      std::complex<double> tmp_1955;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1955 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1954 += tmp_1955;
      tmp_1953 += (Conj(ZEL(gI2,j2))) * tmp_1954;
   }
   tmp_1952 += tmp_1953;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1952;

   return result;
}

std::complex<double> CLASSNAME::CpbarFehhFePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1956;
   std::complex<double> tmp_1957;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1958;
      std::complex<double> tmp_1959;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1959 += Conj(Ye(j1,j2))*ZER(gI2,j1);
      }
      tmp_1958 += tmp_1959;
      tmp_1957 += (ZEL(gO1,j2)) * tmp_1958;
   }
   tmp_1956 += tmp_1957;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1956;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeHpmFvPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1960;
   std::complex<double> tmp_1961;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1961 += Conj(ZER(gO2,j1))*Ye(j1,gI2);
   }
   tmp_1960 += tmp_1961;
   result += (ZP(gI1,0)) * tmp_1960;

   return result;
}

double CLASSNAME::CpbarFeHpmFvPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarChiSePL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1962;
   std::complex<double> tmp_1963;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1964;
      std::complex<double> tmp_1965;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1965 += Conj(ZER(gO2,j1))*Ye(j1,j2);
      }
      tmp_1964 += tmp_1965;
      tmp_1963 += (Conj(ZE(gI2,j2))) * tmp_1964;
   }
   tmp_1962 += tmp_1963;
   result += (-Conj(ZN2(gI1,2))) * tmp_1962;

   return result;
}

std::complex<double> CLASSNAME::CpbarFebarChiSePR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1966;
   std::complex<double> tmp_1967;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1967 += Conj(ZE(gI2,j1))*ZEL(gO1,j1);
   }
   tmp_1966 += tmp_1967;
   result += (0.7071067811865475*(0.7745966692414834*g1*ZN1(gI1,0) + g2*ZN1(gI1
      ,1))) * tmp_1966;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1968;
   std::complex<double> tmp_1969;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1969 += Conj(ZE(gI1,3 + j1))*Conj(ZER(gO2,j1));
   }
   tmp_1968 += tmp_1969;
   result += (-1.0954451150103321*g1*Conj(ZN1(gI2,0))) * tmp_1968;

   return result;
}

std::complex<double> CLASSNAME::CpbarFeSeChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1970;
   std::complex<double> tmp_1971;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1972;
      std::complex<double> tmp_1973;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1973 += Conj(Ye(j1,j2))*Conj(ZE(gI1,3 + j1));
      }
      tmp_1972 += tmp_1973;
      tmp_1971 += (ZEL(gO1,j2)) * tmp_1972;
   }
   tmp_1970 += tmp_1971;
   result += (-ZN2(gI2,2)) * tmp_1970;

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

std::complex<double> CLASSNAME::CpbarFdbarCha1SuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1974;
   std::complex<double> tmp_1975;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1976;
      std::complex<double> tmp_1977;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1977 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1976 += tmp_1977;
      tmp_1975 += (Conj(ZU(gI2,j2))) * tmp_1976;
   }
   tmp_1974 += tmp_1975;
   result += (Conj(UM1(gI1,1))) * tmp_1974;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarCha1SuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1978;
   std::complex<double> tmp_1979;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_1979 += Conj(ZU(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_1978 += tmp_1979;
   result += (-(g2*UP1(gI1,0))) * tmp_1978;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1980;
   std::complex<double> tmp_1981;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1982;
      std::complex<double> tmp_1983;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1983 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1982 += tmp_1983;
      tmp_1981 += (Conj(ZDL(gI1,j2))) * tmp_1982;
   }
   tmp_1980 += tmp_1981;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,0)) * tmp_1980
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdFdAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1984;
   std::complex<double> tmp_1985;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1986;
      std::complex<double> tmp_1987;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1987 += Conj(Yd(j1,j2))*ZDR(gI1,j1);
      }
      tmp_1986 += tmp_1987;
      tmp_1985 += (ZDL(gO1,j2)) * tmp_1986;
   }
   tmp_1984 += tmp_1985;
   result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,0)) * tmp_1984;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1988;
   std::complex<double> tmp_1989;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1990;
      std::complex<double> tmp_1991;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1991 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1990 += tmp_1991;
      tmp_1989 += (Conj(ZDL(gI2,j2))) * tmp_1990;
   }
   tmp_1988 += tmp_1989;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1988;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdhhFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1992;
   std::complex<double> tmp_1993;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1994;
      std::complex<double> tmp_1995;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1995 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_1994 += tmp_1995;
      tmp_1993 += (ZDL(gO1,j2)) * tmp_1994;
   }
   tmp_1992 += tmp_1993;
   result += (-0.7071067811865475*ZH(gI1,0)) * tmp_1992;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_1996;
   std::complex<double> tmp_1997;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_1998;
      std::complex<double> tmp_1999;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_1999 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_1998 += tmp_1999;
      tmp_1997 += (Conj(ZUL(gI2,j2))) * tmp_1998;
   }
   tmp_1996 += tmp_1997;
   result += (ZP(gI1,0)) * tmp_1996;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdHpmFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2000;
   std::complex<double> tmp_2001;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2002;
      std::complex<double> tmp_2003;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2003 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2002 += tmp_2003;
      tmp_2001 += (ZDL(gO1,j2)) * tmp_2002;
   }
   tmp_2000 += tmp_2001;
   result += (ZP(gI1,1)) * tmp_2000;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarChiSdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2004;
   std::complex<double> tmp_2005;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2006;
      std::complex<double> tmp_2007;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2007 += Conj(ZDR(gO2,j1))*Yd(j1,j2);
      }
      tmp_2006 += tmp_2007;
      tmp_2005 += (Conj(ZD(gI2,j2))) * tmp_2006;
   }
   tmp_2004 += tmp_2005;
   result += (-Conj(ZN2(gI1,2))) * tmp_2004;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarChiSdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2008;
   std::complex<double> tmp_2009;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2009 += Conj(ZD(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2008 += tmp_2009;
   result += (-0.2357022603955158*(0.7745966692414834*g1*ZN1(gI1,0) - 3*g2*ZN1(
      gI1,1))) * tmp_2008;

   return result;
}

double CLASSNAME::CpbarFdSuCha2PL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSuCha2PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2010;
   std::complex<double> tmp_2011;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2012;
      std::complex<double> tmp_2013;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2013 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2012 += tmp_2013;
      tmp_2011 += (ZDL(gO1,j2)) * tmp_2012;
   }
   tmp_2010 += tmp_2011;
   result += (UP2(gI2,1)) * tmp_2010;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2014;
   std::complex<double> tmp_2015;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2015 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2014 += tmp_2015;
   result += (-0.3651483716701107*g1*Conj(ZN1(gI2,0))) * tmp_2014;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2016;
   std::complex<double> tmp_2017;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2018;
      std::complex<double> tmp_2019;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2019 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2018 += tmp_2019;
      tmp_2017 += (ZDL(gO1,j2)) * tmp_2018;
   }
   tmp_2016 += tmp_2017;
   result += (-ZN2(gI2,2)) * tmp_2016;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdSdGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2020;
   std::complex<double> tmp_2021;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2021 += Conj(ZD(gI1,3 + j1))*Conj(ZDR(gO2,j1));
   }
   tmp_2020 += tmp_2021;
   result += (1.4142135623730951*g3) * tmp_2020;

   return result;
}

double CLASSNAME::CpbarFdSdGluPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

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

   std::complex<double> tmp_2022;
   std::complex<double> tmp_2023;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2023 += Conj(ZUL(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2022 += tmp_2023;
   result += (-0.7071067811865475*g2) * tmp_2022;

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

double CLASSNAME::CpbarFdbarGluSdPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFdbarGluSdPR(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2024;
   std::complex<double> tmp_2025;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2025 += Conj(ZD(gI2,j1))*ZDL(gO1,j1);
   }
   tmp_2024 += tmp_2025;
   result += (-1.4142135623730951*g3) * tmp_2024;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarCha2SdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2026;
   std::complex<double> tmp_2027;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2028;
      std::complex<double> tmp_2029;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2029 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2028 += tmp_2029;
      tmp_2027 += (Conj(ZD(gI2,j2))) * tmp_2028;
   }
   tmp_2026 += tmp_2027;
   result += (Conj(UP2(gI1,1))) * tmp_2026;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarCha2SdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2030;
   std::complex<double> tmp_2031;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2031 += Conj(ZD(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2030 += tmp_2031;
   result += (-(g2*UM2(gI1,0))) * tmp_2030;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2032;
   std::complex<double> tmp_2033;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2034;
      std::complex<double> tmp_2035;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2035 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2034 += tmp_2035;
      tmp_2033 += (Conj(ZUL(gI1,j2))) * tmp_2034;
   }
   tmp_2032 += tmp_2033;
   result += (std::complex<double>(0,-0.7071067811865475)*ZA(gI2,1)) * tmp_2032
      ;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuFuAhPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2036;
   std::complex<double> tmp_2037;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2038;
      std::complex<double> tmp_2039;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2039 += Conj(Yu(j1,j2))*ZUR(gI1,j1);
      }
      tmp_2038 += tmp_2039;
      tmp_2037 += (ZUL(gO1,j2)) * tmp_2038;
   }
   tmp_2036 += tmp_2037;
   result += (std::complex<double>(0,0.7071067811865475)*ZA(gI2,1)) * tmp_2036;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2040;
   std::complex<double> tmp_2041;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2042;
      std::complex<double> tmp_2043;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2043 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2042 += tmp_2043;
      tmp_2041 += (Conj(ZDL(gI2,j2))) * tmp_2042;
   }
   tmp_2040 += tmp_2041;
   result += (ZP(gI1,1)) * tmp_2040;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuconjHpmFdPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2044;
   std::complex<double> tmp_2045;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2046;
      std::complex<double> tmp_2047;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2047 += Conj(Yd(j1,j2))*ZDR(gI2,j1);
      }
      tmp_2046 += tmp_2047;
      tmp_2045 += (ZUL(gO1,j2)) * tmp_2046;
   }
   tmp_2044 += tmp_2045;
   result += (ZP(gI1,0)) * tmp_2044;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2048;
   std::complex<double> tmp_2049;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2050;
      std::complex<double> tmp_2051;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2051 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2050 += tmp_2051;
      tmp_2049 += (Conj(ZUL(gI2,j2))) * tmp_2050;
   }
   tmp_2048 += tmp_2049;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2048;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuhhFuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2052;
   std::complex<double> tmp_2053;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2054;
      std::complex<double> tmp_2055;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2055 += Conj(Yu(j1,j2))*ZUR(gI2,j1);
      }
      tmp_2054 += tmp_2055;
      tmp_2053 += (ZUL(gO1,j2)) * tmp_2054;
   }
   tmp_2052 += tmp_2053;
   result += (-0.7071067811865475*ZH(gI1,1)) * tmp_2052;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChiSuPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2056;
   std::complex<double> tmp_2057;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2058;
      std::complex<double> tmp_2059;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2059 += Conj(ZUR(gO2,j1))*Yu(j1,j2);
      }
      tmp_2058 += tmp_2059;
      tmp_2057 += (Conj(ZU(gI2,j2))) * tmp_2058;
   }
   tmp_2056 += tmp_2057;
   result += (-Conj(ZN2(gI1,3))) * tmp_2056;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarChiSuPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2060;
   std::complex<double> tmp_2061;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2061 += Conj(ZU(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2060 += tmp_2061;
   result += (-0.2357022603955158*(0.7745966692414834*g1*ZN1(gI1,0) + 3*g2*ZN1(
      gI1,1))) * tmp_2060;

   return result;
}

double CLASSNAME::CpbarFuSdCha1PL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSdCha1PR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2062;
   std::complex<double> tmp_2063;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2064;
      std::complex<double> tmp_2065;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2065 += Conj(Yd(j1,j2))*Conj(ZD(gI1,3 + j1));
      }
      tmp_2064 += tmp_2065;
      tmp_2063 += (ZUL(gO1,j2)) * tmp_2064;
   }
   tmp_2062 += tmp_2063;
   result += (UM1(gI2,1)) * tmp_2062;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPL(unsigned gO2, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2066;
   std::complex<double> tmp_2067;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2067 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2066 += tmp_2067;
   result += (0.7302967433402214*g1*Conj(ZN1(gI2,0))) * tmp_2066;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuChiPR(unsigned gO1, unsigned gI1, unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2068;
   std::complex<double> tmp_2069;
   for (unsigned j2 = 0; j2 < 3; ++j2) {
      std::complex<double> tmp_2070;
      std::complex<double> tmp_2071;
      for (unsigned j1 = 0; j1 < 3; ++j1) {
         tmp_2071 += Conj(Yu(j1,j2))*Conj(ZU(gI1,3 + j1));
      }
      tmp_2070 += tmp_2071;
      tmp_2069 += (ZUL(gO1,j2)) * tmp_2070;
   }
   tmp_2068 += tmp_2069;
   result += (-ZN2(gI2,3)) * tmp_2068;

   return result;
}

std::complex<double> CLASSNAME::CpbarFuSuGluPL(unsigned gO2, unsigned gI1, unsigned ) const
{
   std::complex<double> result;

   std::complex<double> tmp_2072;
   std::complex<double> tmp_2073;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2073 += Conj(ZU(gI1,3 + j1))*Conj(ZUR(gO2,j1));
   }
   tmp_2072 += tmp_2073;
   result += (1.4142135623730951*g3) * tmp_2072;

   return result;
}

double CLASSNAME::CpbarFuSuGluPR(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

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

   std::complex<double> tmp_2074;
   std::complex<double> tmp_2075;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2075 += Conj(ZDL(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2074 += tmp_2075;
   result += (-0.7071067811865475*g2) * tmp_2074;

   return result;
}

double CLASSNAME::CpbarFubarGluSuPL(unsigned , unsigned , unsigned ) const
{
   double result = 0.0;

   result = 0;

   return result;
}

std::complex<double> CLASSNAME::CpbarFubarGluSuPR(unsigned gO1, unsigned , unsigned gI2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2076;
   std::complex<double> tmp_2077;
   for (unsigned j1 = 0; j1 < 3; ++j1) {
      tmp_2077 += Conj(ZU(gI2,j1))*ZUL(gO1,j1);
   }
   tmp_2076 += tmp_2077;
   result += (-1.4142135623730951*g3) * tmp_2076;

   return result;
}


std::complex<double> CLASSNAME::self_energy_Sd(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSOc)*CpUSdconjUSdconjSOcSOc(gO1,gO2));
   result += 4*A0(MVWm)*CpUSdconjUSdconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSdconjUSdVZVZ(gO1,gO2);
   std::complex<double> tmp_2078;
   std::complex<double> tmp_2079;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2079 += A0(MRh(gI1))*CpUSdconjUSdconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2078 += tmp_2079;
   result += (-1) * tmp_2078;
   std::complex<double> tmp_2080;
   std::complex<double> tmp_2081;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2081 += A0(MRpm(gI1))*CpUSdconjUSdconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2080 += tmp_2081;
   result += (-1) * tmp_2080;
   std::complex<double> tmp_2082;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2083;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2083 += (Conj(CpconjUSdbarCha1FuPL(gO2,gI1,gI2))*
            CpconjUSdbarCha1FuPL(gO1,gI1,gI2) + Conj(CpconjUSdbarCha1FuPR(gO2,gI1,
            gI2))*CpconjUSdbarCha1FuPR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MFu(gI2));
      }
      tmp_2082 += tmp_2083;
   }
   result += tmp_2082;
   std::complex<double> tmp_2084;
   std::complex<double> tmp_2085;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2086;
      std::complex<double> tmp_2087;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2087 += B0(p,MCha1(gI1),MFu(gI2))*(Conj(CpconjUSdbarCha1FuPR
            (gO2,gI1,gI2))*CpconjUSdbarCha1FuPL(gO1,gI1,gI2) + Conj(
            CpconjUSdbarCha1FuPL(gO2,gI1,gI2))*CpconjUSdbarCha1FuPR(gO1,gI1,gI2))*
            MFu(gI2);
      }
      tmp_2086 += tmp_2087;
      tmp_2085 += (MCha1(gI1)) * tmp_2086;
   }
   tmp_2084 += tmp_2085;
   result += (-2) * tmp_2084;
   std::complex<double> tmp_2088;
   std::complex<double> tmp_2089;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2089 += A0(MSv(gI1))*CpUSdconjUSdconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2088 += tmp_2089;
   result += (-1) * tmp_2088;
   std::complex<double> tmp_2090;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2091;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2091 += (Conj(CpconjUSdFuCha2PL(gO2,gI1,gI2))*
            CpconjUSdFuCha2PL(gO1,gI1,gI2) + Conj(CpconjUSdFuCha2PR(gO2,gI1,gI2))*
            CpconjUSdFuCha2PR(gO1,gI1,gI2))*G0(p,MFu(gI1),MCha2(gI2));
      }
      tmp_2090 += tmp_2091;
   }
   result += tmp_2090;
   std::complex<double> tmp_2092;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2093;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2093 += (Conj(CpconjUSdFdChiPL(gO2,gI1,gI2))*
            CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPR(gO2,gI1,gI2))*
            CpconjUSdFdChiPR(gO1,gI1,gI2))*G0(p,MFd(gI1),MChi(gI2));
      }
      tmp_2092 += tmp_2093;
   }
   result += tmp_2092;
   std::complex<double> tmp_2094;
   std::complex<double> tmp_2095;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2096;
      std::complex<double> tmp_2097;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2097 += B0(p,MFd(gI1),MChi(gI2))*(Conj(CpconjUSdFdChiPR(gO2,
            gI1,gI2))*CpconjUSdFdChiPL(gO1,gI1,gI2) + Conj(CpconjUSdFdChiPL(gO2,
            gI1,gI2))*CpconjUSdFdChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2096 += tmp_2097;
      tmp_2095 += (MFd(gI1)) * tmp_2096;
   }
   tmp_2094 += tmp_2095;
   result += (-2) * tmp_2094;
   std::complex<double> tmp_2098;
   std::complex<double> tmp_2099;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2100;
      std::complex<double> tmp_2101;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2101 += B0(p,MFu(gI1),MCha2(gI2))*(Conj(CpconjUSdFuCha2PR(
            gO2,gI1,gI2))*CpconjUSdFuCha2PL(gO1,gI1,gI2) + Conj(CpconjUSdFuCha2PL(
            gO2,gI1,gI2))*CpconjUSdFuCha2PR(gO1,gI1,gI2))*MCha2(gI2);
      }
      tmp_2100 += tmp_2101;
      tmp_2099 += (MFu(gI1)) * tmp_2100;
   }
   tmp_2098 += tmp_2099;
   result += (-2) * tmp_2098;
   std::complex<double> tmp_2102;
   std::complex<double> tmp_2103;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2103 += A0(MAh(gI1))*CpUSdconjUSdAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2102 += tmp_2103;
   result += (-0.5) * tmp_2102;
   std::complex<double> tmp_2104;
   std::complex<double> tmp_2105;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2105 += A0(MHpm(gI1))*CpUSdconjUSdconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2104 += tmp_2105;
   result += (-1) * tmp_2104;
   std::complex<double> tmp_2106;
   std::complex<double> tmp_2107;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2107 += A0(Mhh(gI1))*CpUSdconjUSdhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2106 += tmp_2107;
   result += (-0.5) * tmp_2106;
   std::complex<double> tmp_2108;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2109;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2109 += (Conj(CpconjUSdbarChiFdPL(gO2,gI1,gI2))*
            CpconjUSdbarChiFdPL(gO1,gI1,gI2) + Conj(CpconjUSdbarChiFdPR(gO2,gI1,
            gI2))*CpconjUSdbarChiFdPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MFd(gI2));
      }
      tmp_2108 += tmp_2109;
   }
   result += tmp_2108;
   std::complex<double> tmp_2110;
   std::complex<double> tmp_2111;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2112;
      std::complex<double> tmp_2113;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2113 += B0(p,MChi(gI1),MFd(gI2))*(Conj(CpconjUSdbarChiFdPR(
            gO2,gI1,gI2))*CpconjUSdbarChiFdPL(gO1,gI1,gI2) + Conj(
            CpconjUSdbarChiFdPL(gO2,gI1,gI2))*CpconjUSdbarChiFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2112 += tmp_2113;
      tmp_2111 += (MChi(gI1)) * tmp_2112;
   }
   tmp_2110 += tmp_2111;
   result += (-2) * tmp_2110;
   std::complex<double> tmp_2114;
   std::complex<double> tmp_2115;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2115 += A0(MSd(gI1))*CpUSdconjUSdconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2114 += tmp_2115;
   result += (-1) * tmp_2114;
   std::complex<double> tmp_2116;
   std::complex<double> tmp_2117;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2117 += A0(MSe(gI1))*CpUSdconjUSdconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2116 += tmp_2117;
   result += (-1) * tmp_2116;
   std::complex<double> tmp_2118;
   std::complex<double> tmp_2119;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2119 += A0(MSu(gI1))*CpUSdconjUSdconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2118 += tmp_2119;
   result += (-1) * tmp_2118;
   std::complex<double> tmp_2120;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2121;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2121 += B0(p,MSd(gI1),MRh(gI2))*Conj(CpconjUSdSdRh(gO2,gI1,
            gI2))*CpconjUSdSdRh(gO1,gI1,gI2);
      }
      tmp_2120 += tmp_2121;
   }
   result += tmp_2120;
   std::complex<double> tmp_2122;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2123;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2123 += B0(p,MSu(gI1),MRpm(gI2))*Conj(CpconjUSdSuRpm(gO2,gI1
            ,gI2))*CpconjUSdSuRpm(gO1,gI1,gI2);
      }
      tmp_2122 += tmp_2123;
   }
   result += tmp_2122;
   std::complex<double> tmp_2124;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2125;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2125 += B0(p,MSd(gI1),MAh(gI2))*Conj(CpconjUSdSdAh(gO2,gI1,
            gI2))*CpconjUSdSdAh(gO1,gI1,gI2);
      }
      tmp_2124 += tmp_2125;
   }
   result += tmp_2124;
   std::complex<double> tmp_2126;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2127;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2127 += B0(p,MSd(gI1),Mhh(gI2))*Conj(CpconjUSdSdhh(gO2,gI1,
            gI2))*CpconjUSdSdhh(gO1,gI1,gI2);
      }
      tmp_2126 += tmp_2127;
   }
   result += tmp_2126;
   std::complex<double> tmp_2128;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2129;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2129 += B0(p,MSu(gI1),MHpm(gI2))*Conj(CpconjUSdSuHpm(gO2,gI1
            ,gI2))*CpconjUSdSuHpm(gO1,gI1,gI2);
      }
      tmp_2128 += tmp_2129;
   }
   result += tmp_2128;
   std::complex<double> tmp_2130;
   std::complex<double> tmp_2131;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2131 += (Conj(CpconjUSdbarGluFdPL(gO2,1,gI2))*CpconjUSdbarGluFdPL(
         gO1,1,gI2) + Conj(CpconjUSdbarGluFdPR(gO2,1,gI2))*CpconjUSdbarGluFdPR(gO1
         ,1,gI2))*G0(p,MGlu,MFd(gI2));
   }
   tmp_2130 += tmp_2131;
   result += (1.3333333333333333) * tmp_2130;
   std::complex<double> tmp_2132;
   std::complex<double> tmp_2133;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2133 += (Conj(CpconjUSdGluFdPL(gO2,1,gI2))*CpconjUSdGluFdPL(gO1,1,
         gI2) + Conj(CpconjUSdGluFdPR(gO2,1,gI2))*CpconjUSdGluFdPR(gO1,1,gI2))*G0(
         p,MGlu,MFd(gI2));
   }
   tmp_2132 += tmp_2133;
   result += (1.3333333333333333) * tmp_2132;
   std::complex<double> tmp_2134;
   std::complex<double> tmp_2135;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2135 += B0(p,MSOc,MSd(gI2))*Conj(CpconjUSdconjSOcSd(gO2,gI2))*
         CpconjUSdconjSOcSd(gO1,gI2);
   }
   tmp_2134 += tmp_2135;
   result += (1.3333333333333333) * tmp_2134;
   std::complex<double> tmp_2136;
   std::complex<double> tmp_2137;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2137 += Conj(CpconjUSdVGSd(gO2,gI2))*CpconjUSdVGSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   tmp_2136 += tmp_2137;
   result += (1.3333333333333333) * tmp_2136;
   std::complex<double> tmp_2138;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2138 += Conj(CpconjUSdVPSd(gO2,gI2))*CpconjUSdVPSd(gO1,gI2)*F0(p,
         MSd(gI2),0);
   }
   result += tmp_2138;
   std::complex<double> tmp_2139;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2139 += Conj(CpconjUSdVZSd(gO2,gI2))*CpconjUSdVZSd(gO1,gI2)*F0(p,
         MSd(gI2),MVZ);
   }
   result += tmp_2139;
   std::complex<double> tmp_2140;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2140 += Conj(CpconjUSdVWmSu(gO2,gI2))*CpconjUSdVWmSu(gO1,gI2)*F0(p
         ,MSu(gI2),MVWm);
   }
   result += tmp_2140;
   std::complex<double> tmp_2141;
   std::complex<double> tmp_2142;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2142 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdbarGluFdPR(gO2,1,gI2))*
         CpconjUSdbarGluFdPL(gO1,1,gI2) + Conj(CpconjUSdbarGluFdPL(gO2,1,gI2))*
         CpconjUSdbarGluFdPR(gO1,1,gI2))*MFd(gI2);
   }
   tmp_2141 += tmp_2142;
   result += (-2.6666666666666665*MGlu) * tmp_2141;
   std::complex<double> tmp_2143;
   std::complex<double> tmp_2144;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2144 += B0(p,MGlu,MFd(gI2))*(Conj(CpconjUSdGluFdPR(gO2,1,gI2))*
         CpconjUSdGluFdPL(gO1,1,gI2) + Conj(CpconjUSdGluFdPL(gO2,1,gI2))*
         CpconjUSdGluFdPR(gO1,1,gI2))*MFd(gI2);
   }
   tmp_2143 += tmp_2144;
   result += (-2.6666666666666665*MGlu) * tmp_2143;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Sv(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSvconjUSvconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSvconjUSvVZVZ(gO1,gO2);
   std::complex<double> tmp_2145;
   std::complex<double> tmp_2146;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2146 += A0(MRh(gI1))*CpUSvconjUSvconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2145 += tmp_2146;
   result += (-1) * tmp_2145;
   std::complex<double> tmp_2147;
   std::complex<double> tmp_2148;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2148 += A0(MRpm(gI1))*CpUSvconjUSvconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2147 += tmp_2148;
   result += (-1) * tmp_2147;
   std::complex<double> tmp_2149;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2150;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2150 += B0(p,MRpm(gI1),MSe(gI2))*Conj(CpconjUSvconjRpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjRpmSe(gO1,gI1,gI2);
      }
      tmp_2149 += tmp_2150;
   }
   result += tmp_2149;
   std::complex<double> tmp_2151;
   std::complex<double> tmp_2152;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2152 += A0(MSv(gI1))*CpUSvconjUSvconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2151 += tmp_2152;
   result += (-1) * tmp_2151;
   std::complex<double> tmp_2153;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2154;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2154 += (Conj(CpconjUSvFeCha1PL(gO2,gI1,gI2))*
            CpconjUSvFeCha1PL(gO1,gI1,gI2) + Conj(CpconjUSvFeCha1PR(gO2,gI1,gI2))*
            CpconjUSvFeCha1PR(gO1,gI1,gI2))*G0(p,MFe(gI1),MCha1(gI2));
      }
      tmp_2153 += tmp_2154;
   }
   result += tmp_2153;
   std::complex<double> tmp_2155;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2156;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2156 += B0(p,MSv(gI1),MAh(gI2))*Conj(CpconjUSvSvAh(gO2,gI1,
            gI2))*CpconjUSvSvAh(gO1,gI1,gI2);
      }
      tmp_2155 += tmp_2156;
   }
   result += tmp_2155;
   std::complex<double> tmp_2157;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2158;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2158 += B0(p,MSv(gI1),Mhh(gI2))*Conj(CpconjUSvSvhh(gO2,gI1,
            gI2))*CpconjUSvSvhh(gO1,gI1,gI2);
      }
      tmp_2157 += tmp_2158;
   }
   result += tmp_2157;
   std::complex<double> tmp_2159;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2160;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2160 += (Conj(CpconjUSvFvChiPL(gO2,gI1,gI2))*
            CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPR(gO2,gI1,gI2))*
            CpconjUSvFvChiPR(gO1,gI1,gI2))*G0(p,MFv(gI1),MChi(gI2));
      }
      tmp_2159 += tmp_2160;
   }
   result += tmp_2159;
   std::complex<double> tmp_2161;
   std::complex<double> tmp_2162;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2163;
      std::complex<double> tmp_2164;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2164 += B0(p,MFe(gI1),MCha1(gI2))*(Conj(CpconjUSvFeCha1PR(
            gO2,gI1,gI2))*CpconjUSvFeCha1PL(gO1,gI1,gI2) + Conj(CpconjUSvFeCha1PL(
            gO2,gI1,gI2))*CpconjUSvFeCha1PR(gO1,gI1,gI2))*MCha1(gI2);
      }
      tmp_2163 += tmp_2164;
      tmp_2162 += (MFe(gI1)) * tmp_2163;
   }
   tmp_2161 += tmp_2162;
   result += (-2) * tmp_2161;
   std::complex<double> tmp_2165;
   std::complex<double> tmp_2166;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2167;
      std::complex<double> tmp_2168;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2168 += B0(p,MFv(gI1),MChi(gI2))*(Conj(CpconjUSvFvChiPR(gO2,
            gI1,gI2))*CpconjUSvFvChiPL(gO1,gI1,gI2) + Conj(CpconjUSvFvChiPL(gO2,
            gI1,gI2))*CpconjUSvFvChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2167 += tmp_2168;
      tmp_2166 += (MFv(gI1)) * tmp_2167;
   }
   tmp_2165 += tmp_2166;
   result += (-2) * tmp_2165;
   std::complex<double> tmp_2169;
   std::complex<double> tmp_2170;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2170 += A0(MAh(gI1))*CpUSvconjUSvAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2169 += tmp_2170;
   result += (-0.5) * tmp_2169;
   std::complex<double> tmp_2171;
   std::complex<double> tmp_2172;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2172 += A0(MHpm(gI1))*CpUSvconjUSvconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2171 += tmp_2172;
   result += (-1) * tmp_2171;
   std::complex<double> tmp_2173;
   std::complex<double> tmp_2174;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2174 += A0(Mhh(gI1))*CpUSvconjUSvhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2173 += tmp_2174;
   result += (-0.5) * tmp_2173;
   std::complex<double> tmp_2175;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2176;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2176 += B0(p,MHpm(gI1),MSe(gI2))*Conj(CpconjUSvconjHpmSe(gO2
            ,gI1,gI2))*CpconjUSvconjHpmSe(gO1,gI1,gI2);
      }
      tmp_2175 += tmp_2176;
   }
   result += tmp_2175;
   std::complex<double> tmp_2177;
   std::complex<double> tmp_2178;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2178 += A0(MSd(gI1))*CpUSvconjUSvconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2177 += tmp_2178;
   result += (-3) * tmp_2177;
   std::complex<double> tmp_2179;
   std::complex<double> tmp_2180;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2180 += A0(MSe(gI1))*CpUSvconjUSvconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2179 += tmp_2180;
   result += (-1) * tmp_2179;
   std::complex<double> tmp_2181;
   std::complex<double> tmp_2182;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2182 += A0(MSu(gI1))*CpUSvconjUSvconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2181 += tmp_2182;
   result += (-3) * tmp_2181;
   std::complex<double> tmp_2183;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2183 += Conj(CpconjUSvVZSv(gO2,gI2))*CpconjUSvVZSv(gO1,gI2)*F0(p,
         MSv(gI2),MVZ);
   }
   result += tmp_2183;
   std::complex<double> tmp_2184;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2184 += Conj(CpconjUSvconjVWmSe(gO2,gI2))*CpconjUSvconjVWmSe(gO1,
         gI2)*F0(p,MSe(gI2),MVWm);
   }
   result += tmp_2184;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Su(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(A0(MSOc)*CpUSuconjUSuconjSOcSOc(gO1,gO2));
   result += 4*A0(MVWm)*CpUSuconjUSuconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSuconjUSuVZVZ(gO1,gO2);
   std::complex<double> tmp_2185;
   std::complex<double> tmp_2186;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2186 += A0(MRh(gI1))*CpUSuconjUSuconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2185 += tmp_2186;
   result += (-1) * tmp_2185;
   std::complex<double> tmp_2187;
   std::complex<double> tmp_2188;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2188 += A0(MRpm(gI1))*CpUSuconjUSuconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2187 += tmp_2188;
   result += (-1) * tmp_2187;
   std::complex<double> tmp_2189;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2190;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2190 += (Conj(CpconjUSubarCha2FdPL(gO2,gI1,gI2))*
            CpconjUSubarCha2FdPL(gO1,gI1,gI2) + Conj(CpconjUSubarCha2FdPR(gO2,gI1,
            gI2))*CpconjUSubarCha2FdPR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MFd(gI2));
      }
      tmp_2189 += tmp_2190;
   }
   result += tmp_2189;
   std::complex<double> tmp_2191;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2192;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2192 += B0(p,MRpm(gI1),MSd(gI2))*Conj(CpconjUSuconjRpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjRpmSd(gO1,gI1,gI2);
      }
      tmp_2191 += tmp_2192;
   }
   result += tmp_2191;
   std::complex<double> tmp_2193;
   std::complex<double> tmp_2194;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2195;
      std::complex<double> tmp_2196;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2196 += B0(p,MCha2(gI1),MFd(gI2))*(Conj(CpconjUSubarCha2FdPR
            (gO2,gI1,gI2))*CpconjUSubarCha2FdPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarCha2FdPL(gO2,gI1,gI2))*CpconjUSubarCha2FdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2195 += tmp_2196;
      tmp_2194 += (MCha2(gI1)) * tmp_2195;
   }
   tmp_2193 += tmp_2194;
   result += (-2) * tmp_2193;
   std::complex<double> tmp_2197;
   std::complex<double> tmp_2198;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2198 += A0(MSv(gI1))*CpUSuconjUSuconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2197 += tmp_2198;
   result += (-1) * tmp_2197;
   std::complex<double> tmp_2199;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2200;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2200 += (Conj(CpconjUSuFdCha1PL(gO2,gI1,gI2))*
            CpconjUSuFdCha1PL(gO1,gI1,gI2) + Conj(CpconjUSuFdCha1PR(gO2,gI1,gI2))*
            CpconjUSuFdCha1PR(gO1,gI1,gI2))*G0(p,MFd(gI1),MCha1(gI2));
      }
      tmp_2199 += tmp_2200;
   }
   result += tmp_2199;
   std::complex<double> tmp_2201;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2202;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2202 += (Conj(CpconjUSuFuChiPL(gO2,gI1,gI2))*
            CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPR(gO2,gI1,gI2))*
            CpconjUSuFuChiPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MChi(gI2));
      }
      tmp_2201 += tmp_2202;
   }
   result += tmp_2201;
   std::complex<double> tmp_2203;
   std::complex<double> tmp_2204;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2205;
      std::complex<double> tmp_2206;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2206 += B0(p,MFd(gI1),MCha1(gI2))*(Conj(CpconjUSuFdCha1PR(
            gO2,gI1,gI2))*CpconjUSuFdCha1PL(gO1,gI1,gI2) + Conj(CpconjUSuFdCha1PL(
            gO2,gI1,gI2))*CpconjUSuFdCha1PR(gO1,gI1,gI2))*MCha1(gI2);
      }
      tmp_2205 += tmp_2206;
      tmp_2204 += (MFd(gI1)) * tmp_2205;
   }
   tmp_2203 += tmp_2204;
   result += (-2) * tmp_2203;
   std::complex<double> tmp_2207;
   std::complex<double> tmp_2208;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2209;
      std::complex<double> tmp_2210;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2210 += B0(p,MFu(gI1),MChi(gI2))*(Conj(CpconjUSuFuChiPR(gO2,
            gI1,gI2))*CpconjUSuFuChiPL(gO1,gI1,gI2) + Conj(CpconjUSuFuChiPL(gO2,
            gI1,gI2))*CpconjUSuFuChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2209 += tmp_2210;
      tmp_2208 += (MFu(gI1)) * tmp_2209;
   }
   tmp_2207 += tmp_2208;
   result += (-2) * tmp_2207;
   std::complex<double> tmp_2211;
   std::complex<double> tmp_2212;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2212 += A0(MAh(gI1))*CpUSuconjUSuAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2211 += tmp_2212;
   result += (-0.5) * tmp_2211;
   std::complex<double> tmp_2213;
   std::complex<double> tmp_2214;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2214 += A0(MHpm(gI1))*CpUSuconjUSuconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2213 += tmp_2214;
   result += (-1) * tmp_2213;
   std::complex<double> tmp_2215;
   std::complex<double> tmp_2216;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2216 += A0(Mhh(gI1))*CpUSuconjUSuhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2215 += tmp_2216;
   result += (-0.5) * tmp_2215;
   std::complex<double> tmp_2217;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2218;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2218 += (Conj(CpconjUSubarChiFuPL(gO2,gI1,gI2))*
            CpconjUSubarChiFuPL(gO1,gI1,gI2) + Conj(CpconjUSubarChiFuPR(gO2,gI1,
            gI2))*CpconjUSubarChiFuPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MFu(gI2));
      }
      tmp_2217 += tmp_2218;
   }
   result += tmp_2217;
   std::complex<double> tmp_2219;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2220;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2220 += B0(p,MHpm(gI1),MSd(gI2))*Conj(CpconjUSuconjHpmSd(gO2
            ,gI1,gI2))*CpconjUSuconjHpmSd(gO1,gI1,gI2);
      }
      tmp_2219 += tmp_2220;
   }
   result += tmp_2219;
   std::complex<double> tmp_2221;
   std::complex<double> tmp_2222;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2223;
      std::complex<double> tmp_2224;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2224 += B0(p,MChi(gI1),MFu(gI2))*(Conj(CpconjUSubarChiFuPR(
            gO2,gI1,gI2))*CpconjUSubarChiFuPL(gO1,gI1,gI2) + Conj(
            CpconjUSubarChiFuPL(gO2,gI1,gI2))*CpconjUSubarChiFuPR(gO1,gI1,gI2))*
            MFu(gI2);
      }
      tmp_2223 += tmp_2224;
      tmp_2222 += (MChi(gI1)) * tmp_2223;
   }
   tmp_2221 += tmp_2222;
   result += (-2) * tmp_2221;
   std::complex<double> tmp_2225;
   std::complex<double> tmp_2226;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2226 += B0(p,MSu(gI1),MSOc)*Conj(CpconjUSuSuSOc(gO2,gI1))*
         CpconjUSuSuSOc(gO1,gI1);
   }
   tmp_2225 += tmp_2226;
   result += (1.3333333333333333) * tmp_2225;
   std::complex<double> tmp_2227;
   std::complex<double> tmp_2228;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2228 += A0(MSd(gI1))*CpUSuconjUSuconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2227 += tmp_2228;
   result += (-1) * tmp_2227;
   std::complex<double> tmp_2229;
   std::complex<double> tmp_2230;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2230 += A0(MSe(gI1))*CpUSuconjUSuconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2229 += tmp_2230;
   result += (-1) * tmp_2229;
   std::complex<double> tmp_2231;
   std::complex<double> tmp_2232;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2232 += A0(MSu(gI1))*CpUSuconjUSuconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2231 += tmp_2232;
   result += (-1) * tmp_2231;
   std::complex<double> tmp_2233;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2234;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2234 += B0(p,MSu(gI1),MRh(gI2))*Conj(CpconjUSuSuRh(gO2,gI1,
            gI2))*CpconjUSuSuRh(gO1,gI1,gI2);
      }
      tmp_2233 += tmp_2234;
   }
   result += tmp_2233;
   std::complex<double> tmp_2235;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2236;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2236 += B0(p,MSu(gI1),MAh(gI2))*Conj(CpconjUSuSuAh(gO2,gI1,
            gI2))*CpconjUSuSuAh(gO1,gI1,gI2);
      }
      tmp_2235 += tmp_2236;
   }
   result += tmp_2235;
   std::complex<double> tmp_2237;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2238;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2238 += B0(p,MSu(gI1),Mhh(gI2))*Conj(CpconjUSuSuhh(gO2,gI1,
            gI2))*CpconjUSuSuhh(gO1,gI1,gI2);
      }
      tmp_2237 += tmp_2238;
   }
   result += tmp_2237;
   std::complex<double> tmp_2239;
   std::complex<double> tmp_2240;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2240 += (Conj(CpconjUSubarGluFuPL(gO2,1,gI2))*CpconjUSubarGluFuPL(
         gO1,1,gI2) + Conj(CpconjUSubarGluFuPR(gO2,1,gI2))*CpconjUSubarGluFuPR(gO1
         ,1,gI2))*G0(p,MGlu,MFu(gI2));
   }
   tmp_2239 += tmp_2240;
   result += (1.3333333333333333) * tmp_2239;
   std::complex<double> tmp_2241;
   std::complex<double> tmp_2242;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2242 += (Conj(CpconjUSuGluFuPL(gO2,1,gI2))*CpconjUSuGluFuPL(gO1,1,
         gI2) + Conj(CpconjUSuGluFuPR(gO2,1,gI2))*CpconjUSuGluFuPR(gO1,1,gI2))*G0(
         p,MGlu,MFu(gI2));
   }
   tmp_2241 += tmp_2242;
   result += (1.3333333333333333) * tmp_2241;
   std::complex<double> tmp_2243;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2243 += Conj(CpconjUSuconjVWmSd(gO2,gI2))*CpconjUSuconjVWmSd(gO1,
         gI2)*F0(p,MSd(gI2),MVWm);
   }
   result += tmp_2243;
   std::complex<double> tmp_2244;
   std::complex<double> tmp_2245;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2245 += Conj(CpconjUSuVGSu(gO2,gI2))*CpconjUSuVGSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   tmp_2244 += tmp_2245;
   result += (1.3333333333333333) * tmp_2244;
   std::complex<double> tmp_2246;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2246 += Conj(CpconjUSuVPSu(gO2,gI2))*CpconjUSuVPSu(gO1,gI2)*F0(p,
         MSu(gI2),0);
   }
   result += tmp_2246;
   std::complex<double> tmp_2247;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2247 += Conj(CpconjUSuVZSu(gO2,gI2))*CpconjUSuVZSu(gO1,gI2)*F0(p,
         MSu(gI2),MVZ);
   }
   result += tmp_2247;
   std::complex<double> tmp_2248;
   std::complex<double> tmp_2249;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2249 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSubarGluFuPR(gO2,1,gI2))*
         CpconjUSubarGluFuPL(gO1,1,gI2) + Conj(CpconjUSubarGluFuPL(gO2,1,gI2))*
         CpconjUSubarGluFuPR(gO1,1,gI2))*MFu(gI2);
   }
   tmp_2248 += tmp_2249;
   result += (-2.6666666666666665*MGlu) * tmp_2248;
   std::complex<double> tmp_2250;
   std::complex<double> tmp_2251;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2251 += B0(p,MGlu,MFu(gI2))*(Conj(CpconjUSuGluFuPR(gO2,1,gI2))*
         CpconjUSuGluFuPL(gO1,1,gI2) + Conj(CpconjUSuGluFuPL(gO2,1,gI2))*
         CpconjUSuGluFuPR(gO1,1,gI2))*MFu(gI2);
   }
   tmp_2250 += tmp_2251;
   result += (-2.6666666666666665*MGlu) * tmp_2250;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Se(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpUSeconjUSeconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUSeconjUSeVZVZ(gO1,gO2);
   std::complex<double> tmp_2252;
   std::complex<double> tmp_2253;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2253 += A0(MRh(gI1))*CpUSeconjUSeconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2252 += tmp_2253;
   result += (-1) * tmp_2252;
   std::complex<double> tmp_2254;
   std::complex<double> tmp_2255;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2255 += A0(MRpm(gI1))*CpUSeconjUSeconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2254 += tmp_2255;
   result += (-1) * tmp_2254;
   std::complex<double> tmp_2256;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2257;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2257 += (Conj(CpconjUSebarCha1FvPL(gO2,gI1,gI2))*
            CpconjUSebarCha1FvPL(gO1,gI1,gI2) + Conj(CpconjUSebarCha1FvPR(gO2,gI1,
            gI2))*CpconjUSebarCha1FvPR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MFv(gI2));
      }
      tmp_2256 += tmp_2257;
   }
   result += tmp_2256;
   std::complex<double> tmp_2258;
   std::complex<double> tmp_2259;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2260;
      std::complex<double> tmp_2261;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2261 += B0(p,MCha1(gI1),MFv(gI2))*(Conj(CpconjUSebarCha1FvPR
            (gO2,gI1,gI2))*CpconjUSebarCha1FvPL(gO1,gI1,gI2) + Conj(
            CpconjUSebarCha1FvPL(gO2,gI1,gI2))*CpconjUSebarCha1FvPR(gO1,gI1,gI2))*
            MFv(gI2);
      }
      tmp_2260 += tmp_2261;
      tmp_2259 += (MCha1(gI1)) * tmp_2260;
   }
   tmp_2258 += tmp_2259;
   result += (-2) * tmp_2258;
   std::complex<double> tmp_2262;
   std::complex<double> tmp_2263;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2263 += A0(MSv(gI1))*CpUSeconjUSeconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2262 += tmp_2263;
   result += (-1) * tmp_2262;
   std::complex<double> tmp_2264;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2265;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2265 += B0(p,MSv(gI1),MRpm(gI2))*Conj(CpconjUSeSvRpm(gO2,gI1
            ,gI2))*CpconjUSeSvRpm(gO1,gI1,gI2);
      }
      tmp_2264 += tmp_2265;
   }
   result += tmp_2264;
   std::complex<double> tmp_2266;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2267;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2267 += (Conj(CpconjUSeFvCha2PL(gO2,gI1,gI2))*
            CpconjUSeFvCha2PL(gO1,gI1,gI2) + Conj(CpconjUSeFvCha2PR(gO2,gI1,gI2))*
            CpconjUSeFvCha2PR(gO1,gI1,gI2))*G0(p,MFv(gI1),MCha2(gI2));
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
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
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
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
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
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2279 += B0(p,MFv(gI1),MCha2(gI2))*(Conj(CpconjUSeFvCha2PR(
            gO2,gI1,gI2))*CpconjUSeFvCha2PL(gO1,gI1,gI2) + Conj(CpconjUSeFvCha2PL(
            gO2,gI1,gI2))*CpconjUSeFvCha2PR(gO1,gI1,gI2))*MCha2(gI2);
      }
      tmp_2278 += tmp_2279;
      tmp_2277 += (MFv(gI1)) * tmp_2278;
   }
   tmp_2276 += tmp_2277;
   result += (-2) * tmp_2276;
   std::complex<double> tmp_2280;
   std::complex<double> tmp_2281;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2281 += A0(MAh(gI1))*CpUSeconjUSeAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2280 += tmp_2281;
   result += (-0.5) * tmp_2280;
   std::complex<double> tmp_2282;
   std::complex<double> tmp_2283;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2283 += A0(MHpm(gI1))*CpUSeconjUSeconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2282 += tmp_2283;
   result += (-1) * tmp_2282;
   std::complex<double> tmp_2284;
   std::complex<double> tmp_2285;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2285 += A0(Mhh(gI1))*CpUSeconjUSehhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2284 += tmp_2285;
   result += (-0.5) * tmp_2284;
   std::complex<double> tmp_2286;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2287;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2287 += (Conj(CpconjUSebarChiFePL(gO2,gI1,gI2))*
            CpconjUSebarChiFePL(gO1,gI1,gI2) + Conj(CpconjUSebarChiFePR(gO2,gI1,
            gI2))*CpconjUSebarChiFePR(gO1,gI1,gI2))*G0(p,MChi(gI1),MFe(gI2));
      }
      tmp_2286 += tmp_2287;
   }
   result += tmp_2286;
   std::complex<double> tmp_2288;
   std::complex<double> tmp_2289;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2290;
      std::complex<double> tmp_2291;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2291 += B0(p,MChi(gI1),MFe(gI2))*(Conj(CpconjUSebarChiFePR(
            gO2,gI1,gI2))*CpconjUSebarChiFePL(gO1,gI1,gI2) + Conj(
            CpconjUSebarChiFePL(gO2,gI1,gI2))*CpconjUSebarChiFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2290 += tmp_2291;
      tmp_2289 += (MChi(gI1)) * tmp_2290;
   }
   tmp_2288 += tmp_2289;
   result += (-2) * tmp_2288;
   std::complex<double> tmp_2292;
   std::complex<double> tmp_2293;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2293 += A0(MSd(gI1))*CpUSeconjUSeconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2292 += tmp_2293;
   result += (-3) * tmp_2292;
   std::complex<double> tmp_2294;
   std::complex<double> tmp_2295;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2295 += A0(MSe(gI1))*CpUSeconjUSeconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2294 += tmp_2295;
   result += (-1) * tmp_2294;
   std::complex<double> tmp_2296;
   std::complex<double> tmp_2297;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2297 += A0(MSu(gI1))*CpUSeconjUSeconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2296 += tmp_2297;
   result += (-3) * tmp_2296;
   std::complex<double> tmp_2298;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2299;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2299 += B0(p,MSe(gI1),MRh(gI2))*Conj(CpconjUSeSeRh(gO2,gI1,
            gI2))*CpconjUSeSeRh(gO1,gI1,gI2);
      }
      tmp_2298 += tmp_2299;
   }
   result += tmp_2298;
   std::complex<double> tmp_2300;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2301;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2301 += B0(p,MSe(gI1),MAh(gI2))*Conj(CpconjUSeSeAh(gO2,gI1,
            gI2))*CpconjUSeSeAh(gO1,gI1,gI2);
      }
      tmp_2300 += tmp_2301;
   }
   result += tmp_2300;
   std::complex<double> tmp_2302;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2303;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2303 += B0(p,MSe(gI1),Mhh(gI2))*Conj(CpconjUSeSehh(gO2,gI1,
            gI2))*CpconjUSeSehh(gO1,gI1,gI2);
      }
      tmp_2302 += tmp_2303;
   }
   result += tmp_2302;
   std::complex<double> tmp_2304;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_2304 += Conj(CpconjUSeVWmSv(gO2,gI2))*CpconjUSeVWmSv(gO1,gI2)*F0(p
         ,MSv(gI2),MVWm);
   }
   result += tmp_2304;
   std::complex<double> tmp_2305;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2305 += Conj(CpconjUSeVPSe(gO2,gI2))*CpconjUSeVPSe(gO1,gI2)*F0(p,
         MSe(gI2),0);
   }
   result += tmp_2305;
   std::complex<double> tmp_2306;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_2306 += Conj(CpconjUSeVZSe(gO2,gI2))*CpconjUSeVZSe(gO1,gI2)*F0(p,
         MSe(gI2),MVZ);
   }
   result += tmp_2306;

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
   std::complex<double> tmp_2307;
   std::complex<double> tmp_2308;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2308 += A0(MRh(gI1))*CpUhhUhhconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2307 += tmp_2308;
   result += (-1) * tmp_2307;
   std::complex<double> tmp_2309;
   std::complex<double> tmp_2310;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2310 += A0(MRpm(gI1))*CpUhhUhhconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2309 += tmp_2310;
   result += (-1) * tmp_2309;
   std::complex<double> tmp_2311;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2312;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2312 += B0(p,MRh(gI1),MRh(gI2))*Conj(CpUhhconjRhRh(gO2,gI1,
            gI2))*CpUhhconjRhRh(gO1,gI1,gI2);
      }
      tmp_2311 += tmp_2312;
   }
   result += tmp_2311;
   std::complex<double> tmp_2313;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2314;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2314 += B0(p,MRpm(gI1),MRpm(gI2))*Conj(CpUhhconjRpmRpm(gO2,
            gI1,gI2))*CpUhhconjRpmRpm(gO1,gI1,gI2);
      }
      tmp_2313 += tmp_2314;
   }
   result += tmp_2313;
   std::complex<double> tmp_2315;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2316;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2316 += (Conj(CpUhhbarCha1Cha1PL(gO2,gI1,gI2))*
            CpUhhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(CpUhhbarCha1Cha1PR(gO2,gI1,gI2)
            )*CpUhhbarCha1Cha1PR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MCha1(gI2));
      }
      tmp_2315 += tmp_2316;
   }
   result += tmp_2315;
   std::complex<double> tmp_2317;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2318;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2318 += (Conj(CpUhhbarCha2Cha2PL(gO2,gI1,gI2))*
            CpUhhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(CpUhhbarCha2Cha2PR(gO2,gI1,gI2)
            )*CpUhhbarCha2Cha2PR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MCha2(gI2));
      }
      tmp_2317 += tmp_2318;
   }
   result += tmp_2317;
   std::complex<double> tmp_2319;
   std::complex<double> tmp_2320;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2321;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2321 += B0(p,MRh(gI1),MAh(gI2))*Conj(CpUhhconjRhAh(gO2,gI1,
            gI2))*CpUhhconjRhAh(gO1,gI1,gI2);
      }
      tmp_2320 += tmp_2321;
   }
   tmp_2319 += tmp_2320;
   result += (2) * tmp_2319;
   std::complex<double> tmp_2322;
   std::complex<double> tmp_2323;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2324;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2324 += B0(p,MRh(gI1),Mhh(gI2))*Conj(CpUhhconjRhhh(gO2,gI1,
            gI2))*CpUhhconjRhhh(gO1,gI1,gI2);
      }
      tmp_2323 += tmp_2324;
   }
   tmp_2322 += tmp_2323;
   result += (2) * tmp_2322;
   std::complex<double> tmp_2325;
   std::complex<double> tmp_2326;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2327;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2327 += B0(p,MRpm(gI1),MHpm(gI2))*Conj(CpUhhconjRpmHpm(gO2,
            gI1,gI2))*CpUhhconjRpmHpm(gO1,gI1,gI2);
      }
      tmp_2326 += tmp_2327;
   }
   tmp_2325 += tmp_2326;
   result += (2) * tmp_2325;
   std::complex<double> tmp_2328;
   std::complex<double> tmp_2329;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2330;
      std::complex<double> tmp_2331;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2331 += B0(p,MCha1(gI1),MCha1(gI2))*(Conj(CpUhhbarCha1Cha1PR
            (gO2,gI1,gI2))*CpUhhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(
            CpUhhbarCha1Cha1PL(gO2,gI1,gI2))*CpUhhbarCha1Cha1PR(gO1,gI1,gI2))*
            MCha1(gI2);
      }
      tmp_2330 += tmp_2331;
      tmp_2329 += (MCha1(gI1)) * tmp_2330;
   }
   tmp_2328 += tmp_2329;
   result += (-2) * tmp_2328;
   std::complex<double> tmp_2332;
   std::complex<double> tmp_2333;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2334;
      std::complex<double> tmp_2335;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2335 += B0(p,MCha2(gI1),MCha2(gI2))*(Conj(CpUhhbarCha2Cha2PR
            (gO2,gI1,gI2))*CpUhhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(
            CpUhhbarCha2Cha2PL(gO2,gI1,gI2))*CpUhhbarCha2Cha2PR(gO1,gI1,gI2))*
            MCha2(gI2);
      }
      tmp_2334 += tmp_2335;
      tmp_2333 += (MCha2(gI1)) * tmp_2334;
   }
   tmp_2332 += tmp_2333;
   result += (-2) * tmp_2332;
   std::complex<double> tmp_2336;
   std::complex<double> tmp_2337;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2337 += A0(MSv(gI1))*CpUhhUhhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2336 += tmp_2337;
   result += (-1) * tmp_2336;
   std::complex<double> tmp_2338;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2339;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2339 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUhhconjSvSv(gO2,gI1,
            gI2))*CpUhhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2338 += tmp_2339;
   }
   result += tmp_2338;
   std::complex<double> tmp_2340;
   std::complex<double> tmp_2341;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2342;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2342 += (Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))*CpUhhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFdFdPR(gO2,gI1,gI2))*CpUhhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2341 += tmp_2342;
   }
   tmp_2340 += tmp_2341;
   result += (3) * tmp_2340;
   std::complex<double> tmp_2343;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2344;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2344 += (Conj(CpUhhbarFeFePL(gO2,gI1,gI2))*CpUhhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUhhbarFeFePR(gO2,gI1,gI2))*CpUhhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2343 += tmp_2344;
   }
   result += tmp_2343;
   std::complex<double> tmp_2345;
   std::complex<double> tmp_2346;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2347;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2347 += (Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))*CpUhhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUhhbarFuFuPR(gO2,gI1,gI2))*CpUhhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2346 += tmp_2347;
   }
   tmp_2345 += tmp_2346;
   result += (3) * tmp_2345;
   std::complex<double> tmp_2348;
   std::complex<double> tmp_2349;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2350;
      std::complex<double> tmp_2351;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2351 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUhhbarFdFdPR(gO2,gI1
            ,gI2))*CpUhhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUhhbarFdFdPL(gO2,gI1,gI2))
            *CpUhhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2350 += tmp_2351;
      tmp_2349 += (MFd(gI1)) * tmp_2350;
   }
   tmp_2348 += tmp_2349;
   result += (-6) * tmp_2348;
   std::complex<double> tmp_2352;
   std::complex<double> tmp_2353;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2354;
      std::complex<double> tmp_2355;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2355 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUhhbarFeFePR(gO2,gI1
            ,gI2))*CpUhhbarFeFePL(gO1,gI1,gI2) + Conj(CpUhhbarFeFePL(gO2,gI1,gI2))
            *CpUhhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2354 += tmp_2355;
      tmp_2353 += (MFe(gI1)) * tmp_2354;
   }
   tmp_2352 += tmp_2353;
   result += (-2) * tmp_2352;
   std::complex<double> tmp_2356;
   std::complex<double> tmp_2357;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2358;
      std::complex<double> tmp_2359;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2359 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUhhbarFuFuPR(gO2,gI1
            ,gI2))*CpUhhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUhhbarFuFuPL(gO2,gI1,gI2))
            *CpUhhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2358 += tmp_2359;
      tmp_2357 += (MFu(gI1)) * tmp_2358;
   }
   tmp_2356 += tmp_2357;
   result += (-6) * tmp_2356;
   std::complex<double> tmp_2360;
   std::complex<double> tmp_2361;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2361 += A0(MAh(gI1))*CpUhhUhhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2360 += tmp_2361;
   result += (-0.5) * tmp_2360;
   std::complex<double> tmp_2362;
   std::complex<double> tmp_2363;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2363 += A0(MHpm(gI1))*CpUhhUhhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2362 += tmp_2363;
   result += (-1) * tmp_2362;
   std::complex<double> tmp_2364;
   std::complex<double> tmp_2365;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2365 += A0(Mhh(gI1))*CpUhhUhhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2364 += tmp_2365;
   result += (-0.5) * tmp_2364;
   std::complex<double> tmp_2366;
   std::complex<double> tmp_2367;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2368;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2368 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUhhAhAh(gO2,gI1,gI2))
            *CpUhhAhAh(gO1,gI1,gI2);
      }
      tmp_2367 += tmp_2368;
   }
   tmp_2366 += tmp_2367;
   result += (0.5) * tmp_2366;
   std::complex<double> tmp_2369;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2370;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2370 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUhhconjHpmHpm(gO2,
            gI1,gI2))*CpUhhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2369 += tmp_2370;
   }
   result += tmp_2369;
   std::complex<double> tmp_2371;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2372;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2372 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUhhhhAh(gO2,gI1,gI2))
            *CpUhhhhAh(gO1,gI1,gI2);
      }
      tmp_2371 += tmp_2372;
   }
   result += tmp_2371;
   std::complex<double> tmp_2373;
   std::complex<double> tmp_2374;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2375;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2375 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUhhhhhh(gO2,gI1,gI2))
            *CpUhhhhhh(gO1,gI1,gI2);
      }
      tmp_2374 += tmp_2375;
   }
   tmp_2373 += tmp_2374;
   result += (0.5) * tmp_2373;
   std::complex<double> tmp_2376;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2377;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2377 += (Conj(CpUhhbarChiChiPL(gO2,gI1,gI2))*
            CpUhhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUhhbarChiChiPR(gO2,gI1,gI2))*
            CpUhhbarChiChiPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2376 += tmp_2377;
   }
   result += tmp_2376;
   std::complex<double> tmp_2378;
   std::complex<double> tmp_2379;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2380;
      std::complex<double> tmp_2381;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2381 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUhhbarChiChiPR(gO2
            ,gI1,gI2))*CpUhhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUhhbarChiChiPL(gO2,
            gI1,gI2))*CpUhhbarChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2380 += tmp_2381;
      tmp_2379 += (MChi(gI1)) * tmp_2380;
   }
   tmp_2378 += tmp_2379;
   result += (-2) * tmp_2378;
   std::complex<double> tmp_2382;
   std::complex<double> tmp_2383;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2383 += A0(MSd(gI1))*CpUhhUhhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2382 += tmp_2383;
   result += (-3) * tmp_2382;
   std::complex<double> tmp_2384;
   std::complex<double> tmp_2385;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2385 += A0(MSe(gI1))*CpUhhUhhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2384 += tmp_2385;
   result += (-1) * tmp_2384;
   std::complex<double> tmp_2386;
   std::complex<double> tmp_2387;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2387 += A0(MSu(gI1))*CpUhhUhhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2386 += tmp_2387;
   result += (-3) * tmp_2386;
   std::complex<double> tmp_2388;
   std::complex<double> tmp_2389;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2390;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2390 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUhhconjSdSd(gO2,gI1,
            gI2))*CpUhhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2389 += tmp_2390;
   }
   tmp_2388 += tmp_2389;
   result += (3) * tmp_2388;
   std::complex<double> tmp_2391;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2392;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2392 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUhhconjSeSe(gO2,gI1,
            gI2))*CpUhhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2391 += tmp_2392;
   }
   result += tmp_2391;
   std::complex<double> tmp_2393;
   std::complex<double> tmp_2394;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2395;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2395 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUhhconjSuSu(gO2,gI1,
            gI2))*CpUhhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2394 += tmp_2395;
   }
   tmp_2393 += tmp_2394;
   result += (3) * tmp_2393;
   std::complex<double> tmp_2396;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2396 += Conj(CpUhhVZAh(gO2,gI2))*CpUhhVZAh(gO1,gI2)*F0(p,MAh(gI2),
         MVZ);
   }
   result += tmp_2396;
   std::complex<double> tmp_2397;
   std::complex<double> tmp_2398;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2398 += Conj(CpUhhconjVWmHpm(gO2,gI2))*CpUhhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2397 += tmp_2398;
   result += (2) * tmp_2397;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Ah(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmCgWmC(gO1)*CpUAhbargWmCgWmC(gO2));
   result += -(B0(p,MVWm,MVWm)*CpUAhbargWmgWm(gO1)*CpUAhbargWmgWm(gO2));
   result += 4*A0(MVWm)*CpUAhUAhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpUAhUAhVZVZ(gO1,gO2);
   std::complex<double> tmp_2399;
   std::complex<double> tmp_2400;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2400 += A0(MRh(gI1))*CpUAhUAhconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2399 += tmp_2400;
   result += (-1) * tmp_2399;
   std::complex<double> tmp_2401;
   std::complex<double> tmp_2402;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2402 += A0(MRpm(gI1))*CpUAhUAhconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2401 += tmp_2402;
   result += (-1) * tmp_2401;
   std::complex<double> tmp_2403;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2404;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2404 += B0(p,MRh(gI1),MRh(gI2))*Conj(CpUAhconjRhRh(gO2,gI1,
            gI2))*CpUAhconjRhRh(gO1,gI1,gI2);
      }
      tmp_2403 += tmp_2404;
   }
   result += tmp_2403;
   std::complex<double> tmp_2405;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2406;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2406 += B0(p,MRpm(gI1),MRpm(gI2))*Conj(CpUAhconjRpmRpm(gO2,
            gI1,gI2))*CpUAhconjRpmRpm(gO1,gI1,gI2);
      }
      tmp_2405 += tmp_2406;
   }
   result += tmp_2405;
   std::complex<double> tmp_2407;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2408;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2408 += (Conj(CpUAhbarCha1Cha1PL(gO2,gI1,gI2))*
            CpUAhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(CpUAhbarCha1Cha1PR(gO2,gI1,gI2)
            )*CpUAhbarCha1Cha1PR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MCha1(gI2));
      }
      tmp_2407 += tmp_2408;
   }
   result += tmp_2407;
   std::complex<double> tmp_2409;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2410;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2410 += (Conj(CpUAhbarCha2Cha2PL(gO2,gI1,gI2))*
            CpUAhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(CpUAhbarCha2Cha2PR(gO2,gI1,gI2)
            )*CpUAhbarCha2Cha2PR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MCha2(gI2));
      }
      tmp_2409 += tmp_2410;
   }
   result += tmp_2409;
   std::complex<double> tmp_2411;
   std::complex<double> tmp_2412;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2413;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2413 += B0(p,MRh(gI1),MAh(gI2))*Conj(CpUAhconjRhAh(gO2,gI1,
            gI2))*CpUAhconjRhAh(gO1,gI1,gI2);
      }
      tmp_2412 += tmp_2413;
   }
   tmp_2411 += tmp_2412;
   result += (2) * tmp_2411;
   std::complex<double> tmp_2414;
   std::complex<double> tmp_2415;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2416;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2416 += B0(p,MRh(gI1),Mhh(gI2))*Conj(CpUAhconjRhhh(gO2,gI1,
            gI2))*CpUAhconjRhhh(gO1,gI1,gI2);
      }
      tmp_2415 += tmp_2416;
   }
   tmp_2414 += tmp_2415;
   result += (2) * tmp_2414;
   std::complex<double> tmp_2417;
   std::complex<double> tmp_2418;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2419;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2419 += B0(p,MRpm(gI1),MHpm(gI2))*Conj(CpUAhconjRpmHpm(gO2,
            gI1,gI2))*CpUAhconjRpmHpm(gO1,gI1,gI2);
      }
      tmp_2418 += tmp_2419;
   }
   tmp_2417 += tmp_2418;
   result += (2) * tmp_2417;
   std::complex<double> tmp_2420;
   std::complex<double> tmp_2421;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2422;
      std::complex<double> tmp_2423;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2423 += B0(p,MCha1(gI1),MCha1(gI2))*(Conj(CpUAhbarCha1Cha1PR
            (gO2,gI1,gI2))*CpUAhbarCha1Cha1PL(gO1,gI1,gI2) + Conj(
            CpUAhbarCha1Cha1PL(gO2,gI1,gI2))*CpUAhbarCha1Cha1PR(gO1,gI1,gI2))*
            MCha1(gI2);
      }
      tmp_2422 += tmp_2423;
      tmp_2421 += (MCha1(gI1)) * tmp_2422;
   }
   tmp_2420 += tmp_2421;
   result += (-2) * tmp_2420;
   std::complex<double> tmp_2424;
   std::complex<double> tmp_2425;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2426;
      std::complex<double> tmp_2427;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2427 += B0(p,MCha2(gI1),MCha2(gI2))*(Conj(CpUAhbarCha2Cha2PR
            (gO2,gI1,gI2))*CpUAhbarCha2Cha2PL(gO1,gI1,gI2) + Conj(
            CpUAhbarCha2Cha2PL(gO2,gI1,gI2))*CpUAhbarCha2Cha2PR(gO1,gI1,gI2))*
            MCha2(gI2);
      }
      tmp_2426 += tmp_2427;
      tmp_2425 += (MCha2(gI1)) * tmp_2426;
   }
   tmp_2424 += tmp_2425;
   result += (-2) * tmp_2424;
   std::complex<double> tmp_2428;
   std::complex<double> tmp_2429;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2429 += A0(MSv(gI1))*CpUAhUAhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2428 += tmp_2429;
   result += (-1) * tmp_2428;
   std::complex<double> tmp_2430;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2431;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2431 += B0(p,MSv(gI1),MSv(gI2))*Conj(CpUAhconjSvSv(gO2,gI1,
            gI2))*CpUAhconjSvSv(gO1,gI1,gI2);
      }
      tmp_2430 += tmp_2431;
   }
   result += tmp_2430;
   std::complex<double> tmp_2432;
   std::complex<double> tmp_2433;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2434;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2434 += (Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))*CpUAhbarFdFdPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFdFdPR(gO2,gI1,gI2))*CpUAhbarFdFdPR(gO1,
            gI1,gI2))*G0(p,MFd(gI1),MFd(gI2));
      }
      tmp_2433 += tmp_2434;
   }
   tmp_2432 += tmp_2433;
   result += (3) * tmp_2432;
   std::complex<double> tmp_2435;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2436;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2436 += (Conj(CpUAhbarFeFePL(gO2,gI1,gI2))*CpUAhbarFeFePL(
            gO1,gI1,gI2) + Conj(CpUAhbarFeFePR(gO2,gI1,gI2))*CpUAhbarFeFePR(gO1,
            gI1,gI2))*G0(p,MFe(gI1),MFe(gI2));
      }
      tmp_2435 += tmp_2436;
   }
   result += tmp_2435;
   std::complex<double> tmp_2437;
   std::complex<double> tmp_2438;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2439;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2439 += (Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))*CpUAhbarFuFuPL(
            gO1,gI1,gI2) + Conj(CpUAhbarFuFuPR(gO2,gI1,gI2))*CpUAhbarFuFuPR(gO1,
            gI1,gI2))*G0(p,MFu(gI1),MFu(gI2));
      }
      tmp_2438 += tmp_2439;
   }
   tmp_2437 += tmp_2438;
   result += (3) * tmp_2437;
   std::complex<double> tmp_2440;
   std::complex<double> tmp_2441;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2442;
      std::complex<double> tmp_2443;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2443 += B0(p,MFd(gI1),MFd(gI2))*(Conj(CpUAhbarFdFdPR(gO2,gI1
            ,gI2))*CpUAhbarFdFdPL(gO1,gI1,gI2) + Conj(CpUAhbarFdFdPL(gO2,gI1,gI2))
            *CpUAhbarFdFdPR(gO1,gI1,gI2))*MFd(gI2);
      }
      tmp_2442 += tmp_2443;
      tmp_2441 += (MFd(gI1)) * tmp_2442;
   }
   tmp_2440 += tmp_2441;
   result += (-6) * tmp_2440;
   std::complex<double> tmp_2444;
   std::complex<double> tmp_2445;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2446;
      std::complex<double> tmp_2447;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2447 += B0(p,MFe(gI1),MFe(gI2))*(Conj(CpUAhbarFeFePR(gO2,gI1
            ,gI2))*CpUAhbarFeFePL(gO1,gI1,gI2) + Conj(CpUAhbarFeFePL(gO2,gI1,gI2))
            *CpUAhbarFeFePR(gO1,gI1,gI2))*MFe(gI2);
      }
      tmp_2446 += tmp_2447;
      tmp_2445 += (MFe(gI1)) * tmp_2446;
   }
   tmp_2444 += tmp_2445;
   result += (-2) * tmp_2444;
   std::complex<double> tmp_2448;
   std::complex<double> tmp_2449;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2450;
      std::complex<double> tmp_2451;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2451 += B0(p,MFu(gI1),MFu(gI2))*(Conj(CpUAhbarFuFuPR(gO2,gI1
            ,gI2))*CpUAhbarFuFuPL(gO1,gI1,gI2) + Conj(CpUAhbarFuFuPL(gO2,gI1,gI2))
            *CpUAhbarFuFuPR(gO1,gI1,gI2))*MFu(gI2);
      }
      tmp_2450 += tmp_2451;
      tmp_2449 += (MFu(gI1)) * tmp_2450;
   }
   tmp_2448 += tmp_2449;
   result += (-6) * tmp_2448;
   std::complex<double> tmp_2452;
   std::complex<double> tmp_2453;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2453 += A0(MAh(gI1))*CpUAhUAhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2452 += tmp_2453;
   result += (-0.5) * tmp_2452;
   std::complex<double> tmp_2454;
   std::complex<double> tmp_2455;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2455 += A0(MHpm(gI1))*CpUAhUAhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2454 += tmp_2455;
   result += (-1) * tmp_2454;
   std::complex<double> tmp_2456;
   std::complex<double> tmp_2457;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2457 += A0(Mhh(gI1))*CpUAhUAhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2456 += tmp_2457;
   result += (-0.5) * tmp_2456;
   std::complex<double> tmp_2458;
   std::complex<double> tmp_2459;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2460;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2460 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpUAhAhAh(gO2,gI1,gI2))
            *CpUAhAhAh(gO1,gI1,gI2);
      }
      tmp_2459 += tmp_2460;
   }
   tmp_2458 += tmp_2459;
   result += (0.5) * tmp_2458;
   std::complex<double> tmp_2461;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2462;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2462 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpUAhconjHpmHpm(gO2,
            gI1,gI2))*CpUAhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2461 += tmp_2462;
   }
   result += tmp_2461;
   std::complex<double> tmp_2463;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2464;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2464 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpUAhhhAh(gO2,gI1,gI2))
            *CpUAhhhAh(gO1,gI1,gI2);
      }
      tmp_2463 += tmp_2464;
   }
   result += tmp_2463;
   std::complex<double> tmp_2465;
   std::complex<double> tmp_2466;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2467;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2467 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpUAhhhhh(gO2,gI1,gI2))
            *CpUAhhhhh(gO1,gI1,gI2);
      }
      tmp_2466 += tmp_2467;
   }
   tmp_2465 += tmp_2466;
   result += (0.5) * tmp_2465;
   std::complex<double> tmp_2468;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2469;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2469 += (Conj(CpUAhbarChiChiPL(gO2,gI1,gI2))*
            CpUAhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUAhbarChiChiPR(gO2,gI1,gI2))*
            CpUAhbarChiChiPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2468 += tmp_2469;
   }
   result += tmp_2468;
   std::complex<double> tmp_2470;
   std::complex<double> tmp_2471;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2472;
      std::complex<double> tmp_2473;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2473 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpUAhbarChiChiPR(gO2
            ,gI1,gI2))*CpUAhbarChiChiPL(gO1,gI1,gI2) + Conj(CpUAhbarChiChiPL(gO2,
            gI1,gI2))*CpUAhbarChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2472 += tmp_2473;
      tmp_2471 += (MChi(gI1)) * tmp_2472;
   }
   tmp_2470 += tmp_2471;
   result += (-2) * tmp_2470;
   std::complex<double> tmp_2474;
   std::complex<double> tmp_2475;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2475 += A0(MSd(gI1))*CpUAhUAhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2474 += tmp_2475;
   result += (-3) * tmp_2474;
   std::complex<double> tmp_2476;
   std::complex<double> tmp_2477;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2477 += A0(MSe(gI1))*CpUAhUAhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2476 += tmp_2477;
   result += (-1) * tmp_2476;
   std::complex<double> tmp_2478;
   std::complex<double> tmp_2479;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2479 += A0(MSu(gI1))*CpUAhUAhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2478 += tmp_2479;
   result += (-3) * tmp_2478;
   std::complex<double> tmp_2480;
   std::complex<double> tmp_2481;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2482;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2482 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpUAhconjSdSd(gO2,gI1,
            gI2))*CpUAhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2481 += tmp_2482;
   }
   tmp_2480 += tmp_2481;
   result += (3) * tmp_2480;
   std::complex<double> tmp_2483;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2484;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2484 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpUAhconjSeSe(gO2,gI1,
            gI2))*CpUAhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2483 += tmp_2484;
   }
   result += tmp_2483;
   std::complex<double> tmp_2485;
   std::complex<double> tmp_2486;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2487;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2487 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpUAhconjSuSu(gO2,gI1,
            gI2))*CpUAhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2486 += tmp_2487;
   }
   tmp_2485 += tmp_2486;
   result += (3) * tmp_2485;
   std::complex<double> tmp_2488;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2488 += Conj(CpUAhVZhh(gO2,gI2))*CpUAhVZhh(gO1,gI2)*F0(p,Mhh(gI2),
         MVZ);
   }
   result += tmp_2488;
   std::complex<double> tmp_2489;
   std::complex<double> tmp_2490;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2490 += Conj(CpUAhconjVWmHpm(gO2,gI2))*CpUAhconjVWmHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVWm);
   }
   tmp_2489 += tmp_2490;
   result += (2) * tmp_2489;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Rh(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpURhconjURhconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpURhconjURhVZVZ(gO1,gO2);
   std::complex<double> tmp_2491;
   std::complex<double> tmp_2492;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2492 += A0(MRh(gI1))*CpURhconjURhconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2491 += tmp_2492;
   result += (-1) * tmp_2491;
   std::complex<double> tmp_2493;
   std::complex<double> tmp_2494;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2494 += A0(MRpm(gI1))*CpURhconjURhconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2493 += tmp_2494;
   result += (-1) * tmp_2493;
   std::complex<double> tmp_2495;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2495 += Conj(CpconjURhconjRpmVWm(gO2,gI1))*CpconjURhconjRpmVWm(gO1
         ,gI1)*F0(p,MRpm(gI1),MVWm);
   }
   result += tmp_2495;
   std::complex<double> tmp_2496;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2497;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2497 += (Conj(CpconjURhCha2Cha1PL(gO2,gI1,gI2))*
            CpconjURhCha2Cha1PL(gO1,gI1,gI2) + Conj(CpconjURhCha2Cha1PR(gO2,gI1,
            gI2))*CpconjURhCha2Cha1PR(gO1,gI1,gI2))*G0(p,MCha2(gI1),MCha1(gI2));
      }
      tmp_2496 += tmp_2497;
   }
   result += tmp_2496;
   std::complex<double> tmp_2498;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2499;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2499 += B0(p,MRpm(gI1),MHpm(gI2))*Conj(CpconjURhconjRpmHpm(
            gO2,gI1,gI2))*CpconjURhconjRpmHpm(gO1,gI1,gI2);
      }
      tmp_2498 += tmp_2499;
   }
   result += tmp_2498;
   std::complex<double> tmp_2500;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2501;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2501 += B0(p,MRh(gI1),MAh(gI2))*Conj(CpconjURhRhAh(gO2,gI1,
            gI2))*CpconjURhRhAh(gO1,gI1,gI2);
      }
      tmp_2500 += tmp_2501;
   }
   result += tmp_2500;
   std::complex<double> tmp_2502;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2503;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2503 += B0(p,MRh(gI1),Mhh(gI2))*Conj(CpconjURhRhhh(gO2,gI1,
            gI2))*CpconjURhRhhh(gO1,gI1,gI2);
      }
      tmp_2502 += tmp_2503;
   }
   result += tmp_2502;
   std::complex<double> tmp_2504;
   std::complex<double> tmp_2505;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2506;
      std::complex<double> tmp_2507;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2507 += B0(p,MCha2(gI1),MCha1(gI2))*(Conj(
            CpconjURhCha2Cha1PR(gO2,gI1,gI2))*CpconjURhCha2Cha1PL(gO1,gI1,gI2) +
            Conj(CpconjURhCha2Cha1PL(gO2,gI1,gI2))*CpconjURhCha2Cha1PR(gO1,gI1,gI2
            ))*MCha1(gI2);
      }
      tmp_2506 += tmp_2507;
      tmp_2505 += (MCha2(gI1)) * tmp_2506;
   }
   tmp_2504 += tmp_2505;
   result += (-2) * tmp_2504;
   std::complex<double> tmp_2508;
   std::complex<double> tmp_2509;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2509 += A0(MSv(gI1))*CpURhconjURhconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2508 += tmp_2509;
   result += (-1) * tmp_2508;
   std::complex<double> tmp_2510;
   std::complex<double> tmp_2511;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2511 += A0(MAh(gI1))*CpURhconjURhAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2510 += tmp_2511;
   result += (-0.5) * tmp_2510;
   std::complex<double> tmp_2512;
   std::complex<double> tmp_2513;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2513 += A0(MHpm(gI1))*CpURhconjURhconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2512 += tmp_2513;
   result += (-1) * tmp_2512;
   std::complex<double> tmp_2514;
   std::complex<double> tmp_2515;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2515 += A0(Mhh(gI1))*CpURhconjURhhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2514 += tmp_2515;
   result += (-0.5) * tmp_2514;
   std::complex<double> tmp_2516;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2517;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2517 += B0(p,MHpm(gI1),MRpm(gI2))*Conj(CpconjURhconjHpmRpm(
            gO2,gI1,gI2))*CpconjURhconjHpmRpm(gO1,gI1,gI2);
      }
      tmp_2516 += tmp_2517;
   }
   result += tmp_2516;
   std::complex<double> tmp_2518;
   std::complex<double> tmp_2519;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2520;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2520 += B0(p,MAh(gI1),MAh(gI2))*Conj(CpconjURhAhAh(gO2,gI1,
            gI2))*CpconjURhAhAh(gO1,gI1,gI2);
      }
      tmp_2519 += tmp_2520;
   }
   tmp_2518 += tmp_2519;
   result += (0.25) * tmp_2518;
   std::complex<double> tmp_2521;
   std::complex<double> tmp_2522;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2523;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2523 += B0(p,MHpm(gI1),MHpm(gI2))*Conj(CpconjURhconjHpmHpm(
            gO2,gI1,gI2))*CpconjURhconjHpmHpm(gO1,gI1,gI2);
      }
      tmp_2522 += tmp_2523;
   }
   tmp_2521 += tmp_2522;
   result += (0.5) * tmp_2521;
   std::complex<double> tmp_2524;
   std::complex<double> tmp_2525;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2526;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2526 += B0(p,Mhh(gI1),MAh(gI2))*Conj(CpconjURhhhAh(gO2,gI1,
            gI2))*CpconjURhhhAh(gO1,gI1,gI2);
      }
      tmp_2525 += tmp_2526;
   }
   tmp_2524 += tmp_2525;
   result += (0.5) * tmp_2524;
   std::complex<double> tmp_2527;
   std::complex<double> tmp_2528;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2529;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2529 += B0(p,Mhh(gI1),Mhh(gI2))*Conj(CpconjURhhhhh(gO2,gI1,
            gI2))*CpconjURhhhhh(gO1,gI1,gI2);
      }
      tmp_2528 += tmp_2529;
   }
   tmp_2527 += tmp_2528;
   result += (0.25) * tmp_2527;
   std::complex<double> tmp_2530;
   std::complex<double> tmp_2531;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2532;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2532 += (Conj(CpconjURhChiChiPL(gO2,gI1,gI2))*
            CpconjURhChiChiPL(gO1,gI1,gI2) + Conj(CpconjURhChiChiPR(gO2,gI1,gI2))*
            CpconjURhChiChiPR(gO1,gI1,gI2))*G0(p,MChi(gI1),MChi(gI2));
      }
      tmp_2531 += tmp_2532;
   }
   tmp_2530 += tmp_2531;
   result += (0.5) * tmp_2530;
   std::complex<double> tmp_2533;
   std::complex<double> tmp_2534;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2535;
      std::complex<double> tmp_2536;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2536 += B0(p,MChi(gI1),MChi(gI2))*(Conj(CpconjURhChiChiPR(
            gO2,gI1,gI2))*CpconjURhChiChiPL(gO1,gI1,gI2) + Conj(CpconjURhChiChiPL(
            gO2,gI1,gI2))*CpconjURhChiChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2535 += tmp_2536;
      tmp_2534 += (MChi(gI1)) * tmp_2535;
   }
   tmp_2533 += tmp_2534;
   result += (-1) * tmp_2533;
   std::complex<double> tmp_2537;
   std::complex<double> tmp_2538;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2538 += A0(MSd(gI1))*CpURhconjURhconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2537 += tmp_2538;
   result += (-3) * tmp_2537;
   std::complex<double> tmp_2539;
   std::complex<double> tmp_2540;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2540 += A0(MSe(gI1))*CpURhconjURhconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2539 += tmp_2540;
   result += (-1) * tmp_2539;
   std::complex<double> tmp_2541;
   std::complex<double> tmp_2542;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2542 += A0(MSu(gI1))*CpURhconjURhconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2541 += tmp_2542;
   result += (-3) * tmp_2541;
   std::complex<double> tmp_2543;
   std::complex<double> tmp_2544;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2545;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2545 += B0(p,MSd(gI1),MSd(gI2))*Conj(CpconjURhconjSdSd(gO2,
            gI1,gI2))*CpconjURhconjSdSd(gO1,gI1,gI2);
      }
      tmp_2544 += tmp_2545;
   }
   tmp_2543 += tmp_2544;
   result += (1.5) * tmp_2543;
   std::complex<double> tmp_2546;
   std::complex<double> tmp_2547;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2548;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2548 += B0(p,MSe(gI1),MSe(gI2))*Conj(CpconjURhconjSeSe(gO2,
            gI1,gI2))*CpconjURhconjSeSe(gO1,gI1,gI2);
      }
      tmp_2547 += tmp_2548;
   }
   tmp_2546 += tmp_2547;
   result += (0.5) * tmp_2546;
   std::complex<double> tmp_2549;
   std::complex<double> tmp_2550;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2551;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2551 += B0(p,MSu(gI1),MSu(gI2))*Conj(CpconjURhconjSuSu(gO2,
            gI1,gI2))*CpconjURhconjSuSu(gO1,gI1,gI2);
      }
      tmp_2550 += tmp_2551;
   }
   tmp_2549 += tmp_2550;
   result += (1.5) * tmp_2549;
   std::complex<double> tmp_2552;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2552 += Conj(CpconjURhVZRh(gO2,gI2))*CpconjURhVZRh(gO1,gI2)*F0(p,
         MRh(gI2),MVZ);
   }
   result += tmp_2552;
   std::complex<double> tmp_2553;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2553 += Conj(CpconjURhconjVWmRpm(gO2,gI2))*CpconjURhconjVWmRpm(gO1
         ,gI2)*F0(p,MRpm(gI2),MVWm);
   }
   result += tmp_2553;

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
   std::complex<double> tmp_2554;
   std::complex<double> tmp_2555;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2555 += A0(MRh(gI1))*CpUHpmconjUHpmconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2554 += tmp_2555;
   result += (-1) * tmp_2554;
   std::complex<double> tmp_2556;
   std::complex<double> tmp_2557;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2557 += A0(MRpm(gI1))*CpUHpmconjUHpmconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2556 += tmp_2557;
   result += (-1) * tmp_2556;
   std::complex<double> tmp_2558;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2559;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2559 += B0(p,MRh(gI1),MRpm(gI2))*Conj(CpconjUHpmconjRhRpm(
            gO2,gI1,gI2))*CpconjUHpmconjRhRpm(gO1,gI1,gI2);
      }
      tmp_2558 += tmp_2559;
   }
   result += tmp_2558;
   std::complex<double> tmp_2560;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2561;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2561 += B0(p,MRpm(gI1),MRh(gI2))*Conj(CpconjUHpmRpmRh(gO2,
            gI1,gI2))*CpconjUHpmRpmRh(gO1,gI1,gI2);
      }
      tmp_2560 += tmp_2561;
   }
   result += tmp_2560;
   std::complex<double> tmp_2562;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2563;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2563 += B0(p,MRh(gI1),MHpm(gI2))*Conj(CpconjUHpmconjRhHpm(
            gO2,gI1,gI2))*CpconjUHpmconjRhHpm(gO1,gI1,gI2);
      }
      tmp_2562 += tmp_2563;
   }
   result += tmp_2562;
   std::complex<double> tmp_2564;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2565;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2565 += B0(p,MRpm(gI1),MAh(gI2))*Conj(CpconjUHpmRpmAh(gO2,
            gI1,gI2))*CpconjUHpmRpmAh(gO1,gI1,gI2);
      }
      tmp_2564 += tmp_2565;
   }
   result += tmp_2564;
   std::complex<double> tmp_2566;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2567;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2567 += B0(p,MRpm(gI1),Mhh(gI2))*Conj(CpconjUHpmRpmhh(gO2,
            gI1,gI2))*CpconjUHpmRpmhh(gO1,gI1,gI2);
      }
      tmp_2566 += tmp_2567;
   }
   result += tmp_2566;
   std::complex<double> tmp_2568;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2569;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2569 += (Conj(CpconjUHpmbarCha1ChiPL(gO2,gI1,gI2))*
            CpconjUHpmbarCha1ChiPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarCha1ChiPR(gO2,
            gI1,gI2))*CpconjUHpmbarCha1ChiPR(gO1,gI1,gI2))*G0(p,MCha1(gI1),MChi(
            gI2));
      }
      tmp_2568 += tmp_2569;
   }
   result += tmp_2568;
   std::complex<double> tmp_2570;
   std::complex<double> tmp_2571;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2572;
      std::complex<double> tmp_2573;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2573 += B0(p,MCha1(gI1),MChi(gI2))*(Conj(
            CpconjUHpmbarCha1ChiPR(gO2,gI1,gI2))*CpconjUHpmbarCha1ChiPL(gO1,gI1,
            gI2) + Conj(CpconjUHpmbarCha1ChiPL(gO2,gI1,gI2))*
            CpconjUHpmbarCha1ChiPR(gO1,gI1,gI2))*MChi(gI2);
      }
      tmp_2572 += tmp_2573;
      tmp_2571 += (MCha1(gI1)) * tmp_2572;
   }
   tmp_2570 += tmp_2571;
   result += (-2) * tmp_2570;
   std::complex<double> tmp_2574;
   std::complex<double> tmp_2575;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2575 += A0(MSv(gI1))*CpUHpmconjUHpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2574 += tmp_2575;
   result += (-1) * tmp_2574;
   std::complex<double> tmp_2576;
   std::complex<double> tmp_2577;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2578;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2578 += (Conj(CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*
            CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFuFdPR(gO2,gI1,
            gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*G0(p,MFu(gI1),MFd(gI2));
      }
      tmp_2577 += tmp_2578;
   }
   tmp_2576 += tmp_2577;
   result += (3) * tmp_2576;
   std::complex<double> tmp_2579;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2580;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2580 += (Conj(CpconjUHpmbarFvFePL(gO2,gI1,gI2))*
            CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(CpconjUHpmbarFvFePR(gO2,gI1,
            gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*G0(p,MFv(gI1),MFe(gI2));
      }
      tmp_2579 += tmp_2580;
   }
   result += tmp_2579;
   std::complex<double> tmp_2581;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2582;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2582 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjUHpmconjSvSe(gO2,
            gI1,gI2))*CpconjUHpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2581 += tmp_2582;
   }
   result += tmp_2581;
   std::complex<double> tmp_2583;
   std::complex<double> tmp_2584;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2585;
      std::complex<double> tmp_2586;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2586 += B0(p,MFu(gI1),MFd(gI2))*(Conj(CpconjUHpmbarFuFdPR(
            gO2,gI1,gI2))*CpconjUHpmbarFuFdPL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFuFdPL(gO2,gI1,gI2))*CpconjUHpmbarFuFdPR(gO1,gI1,gI2))*
            MFd(gI2);
      }
      tmp_2585 += tmp_2586;
      tmp_2584 += (MFu(gI1)) * tmp_2585;
   }
   tmp_2583 += tmp_2584;
   result += (-6) * tmp_2583;
   std::complex<double> tmp_2587;
   std::complex<double> tmp_2588;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2589;
      std::complex<double> tmp_2590;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2590 += B0(p,MFv(gI1),MFe(gI2))*(Conj(CpconjUHpmbarFvFePR(
            gO2,gI1,gI2))*CpconjUHpmbarFvFePL(gO1,gI1,gI2) + Conj(
            CpconjUHpmbarFvFePL(gO2,gI1,gI2))*CpconjUHpmbarFvFePR(gO1,gI1,gI2))*
            MFe(gI2);
      }
      tmp_2589 += tmp_2590;
      tmp_2588 += (MFv(gI1)) * tmp_2589;
   }
   tmp_2587 += tmp_2588;
   result += (-2) * tmp_2587;
   std::complex<double> tmp_2591;
   std::complex<double> tmp_2592;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2592 += A0(MAh(gI1))*CpUHpmconjUHpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2591 += tmp_2592;
   result += (-0.5) * tmp_2591;
   std::complex<double> tmp_2593;
   std::complex<double> tmp_2594;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2594 += A0(MHpm(gI1))*CpUHpmconjUHpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2593 += tmp_2594;
   result += (-1) * tmp_2593;
   std::complex<double> tmp_2595;
   std::complex<double> tmp_2596;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2596 += A0(Mhh(gI1))*CpUHpmconjUHpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2595 += tmp_2596;
   result += (-0.5) * tmp_2595;
   std::complex<double> tmp_2597;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2598;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2598 += (Conj(CpconjUHpmbarChiCha2PL(gO2,gI1,gI2))*
            CpconjUHpmbarChiCha2PL(gO1,gI1,gI2) + Conj(CpconjUHpmbarChiCha2PR(gO2,
            gI1,gI2))*CpconjUHpmbarChiCha2PR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha2(
            gI2));
      }
      tmp_2597 += tmp_2598;
   }
   result += tmp_2597;
   std::complex<double> tmp_2599;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2600;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2600 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjUHpmHpmAh(gO2,
            gI1,gI2))*CpconjUHpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2599 += tmp_2600;
   }
   result += tmp_2599;
   std::complex<double> tmp_2601;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2602;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2602 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjUHpmHpmhh(gO2,
            gI1,gI2))*CpconjUHpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2601 += tmp_2602;
   }
   result += tmp_2601;
   std::complex<double> tmp_2603;
   std::complex<double> tmp_2604;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2605;
      std::complex<double> tmp_2606;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2606 += B0(p,MChi(gI1),MCha2(gI2))*(Conj(
            CpconjUHpmbarChiCha2PR(gO2,gI1,gI2))*CpconjUHpmbarChiCha2PL(gO1,gI1,
            gI2) + Conj(CpconjUHpmbarChiCha2PL(gO2,gI1,gI2))*
            CpconjUHpmbarChiCha2PR(gO1,gI1,gI2))*MCha2(gI2);
      }
      tmp_2605 += tmp_2606;
      tmp_2604 += (MChi(gI1)) * tmp_2605;
   }
   tmp_2603 += tmp_2604;
   result += (-2) * tmp_2603;
   std::complex<double> tmp_2607;
   std::complex<double> tmp_2608;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2608 += A0(MSd(gI1))*CpUHpmconjUHpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2607 += tmp_2608;
   result += (-3) * tmp_2607;
   std::complex<double> tmp_2609;
   std::complex<double> tmp_2610;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2610 += A0(MSe(gI1))*CpUHpmconjUHpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2609 += tmp_2610;
   result += (-1) * tmp_2609;
   std::complex<double> tmp_2611;
   std::complex<double> tmp_2612;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2612 += A0(MSu(gI1))*CpUHpmconjUHpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2611 += tmp_2612;
   result += (-3) * tmp_2611;
   std::complex<double> tmp_2613;
   std::complex<double> tmp_2614;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2615;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2615 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjUHpmconjSuSd(gO2,
            gI1,gI2))*CpconjUHpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2614 += tmp_2615;
   }
   tmp_2613 += tmp_2614;
   result += (3) * tmp_2613;
   std::complex<double> tmp_2616;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2616 += Conj(CpconjUHpmVWmAh(gO2,gI2))*CpconjUHpmVWmAh(gO1,gI2)*F0
         (p,MAh(gI2),MVWm);
   }
   result += tmp_2616;
   std::complex<double> tmp_2617;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2617 += Conj(CpconjUHpmVWmhh(gO2,gI2))*CpconjUHpmVWmhh(gO1,gI2)*F0
         (p,Mhh(gI2),MVWm);
   }
   result += tmp_2617;
   std::complex<double> tmp_2618;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2618 += Conj(CpconjUHpmVPHpm(gO2,gI2))*CpconjUHpmVPHpm(gO1,gI2)*F0
         (p,MHpm(gI2),0);
   }
   result += tmp_2618;
   std::complex<double> tmp_2619;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2619 += Conj(CpconjUHpmVZHpm(gO2,gI2))*CpconjUHpmVZHpm(gO1,gI2)*F0
         (p,MHpm(gI2),MVZ);
   }
   result += tmp_2619;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Rpm(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   result += 4*A0(MVWm)*CpURpmconjURpmconjVWmVWm(gO1,gO2);
   result += 2*A0(MVZ)*CpURpmconjURpmVZVZ(gO1,gO2);
   std::complex<double> tmp_2620;
   std::complex<double> tmp_2621;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2621 += A0(MRh(gI1))*CpURpmconjURpmconjRhRh(gO1,gO2,gI1,gI1);
   }
   tmp_2620 += tmp_2621;
   result += (-1) * tmp_2620;
   std::complex<double> tmp_2622;
   std::complex<double> tmp_2623;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2623 += A0(MRpm(gI1))*CpURpmconjURpmconjRpmRpm(gO1,gO2,gI1,gI1);
   }
   tmp_2622 += tmp_2623;
   result += (-1) * tmp_2622;
   std::complex<double> tmp_2624;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2624 += Conj(CpconjURpmconjRhVWm(gO2,gI1))*CpconjURpmconjRhVWm(gO1
         ,gI1)*F0(p,MRh(gI1),MVWm);
   }
   result += tmp_2624;
   std::complex<double> tmp_2625;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2626;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2626 += B0(p,MRh(gI1),MHpm(gI2))*Conj(CpconjURpmconjRhHpm(
            gO2,gI1,gI2))*CpconjURpmconjRhHpm(gO1,gI1,gI2);
      }
      tmp_2625 += tmp_2626;
   }
   result += tmp_2625;
   std::complex<double> tmp_2627;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2628;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2628 += B0(p,MRh(gI1),MHpm(gI2))*Conj(CpconjURpmRhHpm(gO2,
            gI1,gI2))*CpconjURpmRhHpm(gO1,gI1,gI2);
      }
      tmp_2627 += tmp_2628;
   }
   result += tmp_2627;
   std::complex<double> tmp_2629;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2630;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2630 += B0(p,MRpm(gI1),MAh(gI2))*Conj(CpconjURpmRpmAh(gO2,
            gI1,gI2))*CpconjURpmRpmAh(gO1,gI1,gI2);
      }
      tmp_2629 += tmp_2630;
   }
   result += tmp_2629;
   std::complex<double> tmp_2631;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2632;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2632 += B0(p,MRpm(gI1),Mhh(gI2))*Conj(CpconjURpmRpmhh(gO2,
            gI1,gI2))*CpconjURpmRpmhh(gO1,gI1,gI2);
      }
      tmp_2631 += tmp_2632;
   }
   result += tmp_2631;
   std::complex<double> tmp_2633;
   std::complex<double> tmp_2634;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2634 += A0(MSv(gI1))*CpURpmconjURpmconjSvSv(gO1,gO2,gI1,gI1);
   }
   tmp_2633 += tmp_2634;
   result += (-1) * tmp_2633;
   std::complex<double> tmp_2635;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2636;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2636 += B0(p,MSv(gI1),MSe(gI2))*Conj(CpconjURpmconjSvSe(gO2,
            gI1,gI2))*CpconjURpmconjSvSe(gO1,gI1,gI2);
      }
      tmp_2635 += tmp_2636;
   }
   result += tmp_2635;
   std::complex<double> tmp_2637;
   std::complex<double> tmp_2638;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2638 += A0(MAh(gI1))*CpURpmconjURpmAhAh(gO1,gO2,gI1,gI1);
   }
   tmp_2637 += tmp_2638;
   result += (-0.5) * tmp_2637;
   std::complex<double> tmp_2639;
   std::complex<double> tmp_2640;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2640 += A0(MHpm(gI1))*CpURpmconjURpmconjHpmHpm(gO1,gO2,gI1,gI1);
   }
   tmp_2639 += tmp_2640;
   result += (-1) * tmp_2639;
   std::complex<double> tmp_2641;
   std::complex<double> tmp_2642;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2642 += A0(Mhh(gI1))*CpURpmconjURpmhhhh(gO1,gO2,gI1,gI1);
   }
   tmp_2641 += tmp_2642;
   result += (-0.5) * tmp_2641;
   std::complex<double> tmp_2643;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2644;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2644 += (Conj(CpconjURpmbarChibarCha1PL(gO2,gI1,gI2))*
            CpconjURpmbarChibarCha1PL(gO1,gI1,gI2) + Conj(
            CpconjURpmbarChibarCha1PR(gO2,gI1,gI2))*CpconjURpmbarChibarCha1PR(gO1,
            gI1,gI2))*G0(p,MChi(gI1),MCha1(gI2));
      }
      tmp_2643 += tmp_2644;
   }
   result += tmp_2643;
   std::complex<double> tmp_2645;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2646;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2646 += (Conj(CpconjURpmChiCha2PL(gO2,gI1,gI2))*
            CpconjURpmChiCha2PL(gO1,gI1,gI2) + Conj(CpconjURpmChiCha2PR(gO2,gI1,
            gI2))*CpconjURpmChiCha2PR(gO1,gI1,gI2))*G0(p,MChi(gI1),MCha2(gI2));
      }
      tmp_2645 += tmp_2646;
   }
   result += tmp_2645;
   std::complex<double> tmp_2647;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2648;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2648 += B0(p,MHpm(gI1),MAh(gI2))*Conj(CpconjURpmHpmAh(gO2,
            gI1,gI2))*CpconjURpmHpmAh(gO1,gI1,gI2);
      }
      tmp_2647 += tmp_2648;
   }
   result += tmp_2647;
   std::complex<double> tmp_2649;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2650;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2650 += B0(p,MHpm(gI1),Mhh(gI2))*Conj(CpconjURpmHpmhh(gO2,
            gI1,gI2))*CpconjURpmHpmhh(gO1,gI1,gI2);
      }
      tmp_2649 += tmp_2650;
   }
   result += tmp_2649;
   std::complex<double> tmp_2651;
   std::complex<double> tmp_2652;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2653;
      std::complex<double> tmp_2654;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2654 += B0(p,MChi(gI1),MCha1(gI2))*(Conj(
            CpconjURpmbarChibarCha1PR(gO2,gI1,gI2))*CpconjURpmbarChibarCha1PL(gO1,
            gI1,gI2) + Conj(CpconjURpmbarChibarCha1PL(gO2,gI1,gI2))*
            CpconjURpmbarChibarCha1PR(gO1,gI1,gI2))*MCha1(gI2);
      }
      tmp_2653 += tmp_2654;
      tmp_2652 += (MChi(gI1)) * tmp_2653;
   }
   tmp_2651 += tmp_2652;
   result += (-2) * tmp_2651;
   std::complex<double> tmp_2655;
   std::complex<double> tmp_2656;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2657;
      std::complex<double> tmp_2658;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2658 += B0(p,MChi(gI1),MCha2(gI2))*(Conj(CpconjURpmChiCha2PR
            (gO2,gI1,gI2))*CpconjURpmChiCha2PL(gO1,gI1,gI2) + Conj(
            CpconjURpmChiCha2PL(gO2,gI1,gI2))*CpconjURpmChiCha2PR(gO1,gI1,gI2))*
            MCha2(gI2);
      }
      tmp_2657 += tmp_2658;
      tmp_2656 += (MChi(gI1)) * tmp_2657;
   }
   tmp_2655 += tmp_2656;
   result += (-2) * tmp_2655;
   std::complex<double> tmp_2659;
   std::complex<double> tmp_2660;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2660 += A0(MSd(gI1))*CpURpmconjURpmconjSdSd(gO1,gO2,gI1,gI1);
   }
   tmp_2659 += tmp_2660;
   result += (-3) * tmp_2659;
   std::complex<double> tmp_2661;
   std::complex<double> tmp_2662;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2662 += A0(MSe(gI1))*CpURpmconjURpmconjSeSe(gO1,gO2,gI1,gI1);
   }
   tmp_2661 += tmp_2662;
   result += (-1) * tmp_2661;
   std::complex<double> tmp_2663;
   std::complex<double> tmp_2664;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2664 += A0(MSu(gI1))*CpURpmconjURpmconjSuSu(gO1,gO2,gI1,gI1);
   }
   tmp_2663 += tmp_2664;
   result += (-3) * tmp_2663;
   std::complex<double> tmp_2665;
   std::complex<double> tmp_2666;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2667;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2667 += B0(p,MSu(gI1),MSd(gI2))*Conj(CpconjURpmconjSuSd(gO2,
            gI1,gI2))*CpconjURpmconjSuSd(gO1,gI1,gI2);
      }
      tmp_2666 += tmp_2667;
   }
   tmp_2665 += tmp_2666;
   result += (3) * tmp_2665;
   std::complex<double> tmp_2668;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2668 += Conj(CpconjURpmVWmRh(gO2,gI2))*CpconjURpmVWmRh(gO1,gI2)*F0
         (p,MRh(gI2),MVWm);
   }
   result += tmp_2668;
   std::complex<double> tmp_2669;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2669 += Conj(CpconjURpmVPRpm(gO2,gI2))*CpconjURpmVPRpm(gO1,gI2)*F0
         (p,MRpm(gI2),0);
   }
   result += tmp_2669;
   std::complex<double> tmp_2670;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2670 += Conj(CpconjURpmVZRpm(gO2,gI2))*CpconjURpmVZRpm(gO1,gI2)*F0
         (p,MRpm(gI2),MVZ);
   }
   result += tmp_2670;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_SOc(double p ) const
{
   std::complex<double> result;

   result += -(A0(MSOc)*CpSOcconjSOcconjSOcSOc());
   result += 3*AbsSqr(CpconjSOcVGSOc())*F0(p,MSOc,0);
   result += -1.5*(AbsSqr(CpconjSOcbarGluGluPL(1,1)) + AbsSqr(
      CpconjSOcbarGluGluPR(1,1)))*G0(p,MGlu,MGlu);
   std::complex<double> tmp_2671;
   std::complex<double> tmp_2672;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2672 += A0(MSd(gI1))*CpSOcconjSOcconjSdSd(gI1,gI1);
   }
   tmp_2671 += tmp_2672;
   result += (-1) * tmp_2671;
   std::complex<double> tmp_2673;
   std::complex<double> tmp_2674;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2674 += A0(MSu(gI1))*CpSOcconjSOcconjSuSu(gI1,gI1);
   }
   tmp_2673 += tmp_2674;
   result += (-1) * tmp_2673;
   std::complex<double> tmp_2675;
   std::complex<double> tmp_2676;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2677;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2677 += AbsSqr(CpconjSOcconjSdSd(gI1,gI2))*B0(p,MSd(gI1),MSd
            (gI2));
      }
      tmp_2676 += tmp_2677;
   }
   tmp_2675 += tmp_2676;
   result += (0.25) * tmp_2675;
   std::complex<double> tmp_2678;
   std::complex<double> tmp_2679;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2680;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2680 += AbsSqr(CpconjSOcconjSuSu(gI1,gI2))*B0(p,MSu(gI1),MSu
            (gI2));
      }
      tmp_2679 += tmp_2680;
   }
   tmp_2678 += tmp_2679;
   result += (0.25) * tmp_2678;
   result += 3*B0(p,MGlu,MGlu)*(Conj(CpconjSOcbarGluGluPR(1,1))*
      CpconjSOcbarGluGluPL(1,1) + Conj(CpconjSOcbarGluGluPL(1,1))*
      CpconjSOcbarGluGluPR(1,1))*Sqr(MGlu);

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ(double p ) const
{
   std::complex<double> result;

   result += AbsSqr(CpVZbargWmCgWmC())*B00(p,MVWm,MVWm);
   result += AbsSqr(CpVZbargWmgWm())*B00(p,MVWm,MVWm);
   result += -(A0(MVWm)*(4*CpVZVZconjVWmVWm1() + CpVZVZconjVWmVWm2() +
      CpVZVZconjVWmVWm3()));
   std::complex<double> tmp_2681;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2681 += A0(MRh(gI1))*CpVZVZconjRhRh(gI1,gI1);
   }
   result += tmp_2681;
   std::complex<double> tmp_2682;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2682 += A0(MRpm(gI1))*CpVZVZconjRpmRpm(gI1,gI1);
   }
   result += tmp_2682;
   std::complex<double> tmp_2683;
   std::complex<double> tmp_2684;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2685;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2685 += AbsSqr(CpVZconjRhRh(gI1,gI2))*B00(p,MRh(gI1),MRh(gI2
            ));
      }
      tmp_2684 += tmp_2685;
   }
   tmp_2683 += tmp_2684;
   result += (-4) * tmp_2683;
   std::complex<double> tmp_2686;
   std::complex<double> tmp_2687;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2688;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2688 += AbsSqr(CpVZconjRpmRpm(gI1,gI2))*B00(p,MRpm(gI1),MRpm
            (gI2));
      }
      tmp_2687 += tmp_2688;
   }
   tmp_2686 += tmp_2687;
   result += (-4) * tmp_2686;
   std::complex<double> tmp_2689;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2690;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2690 += (AbsSqr(CpVZbarCha1Cha1PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha1Cha1PR(gI1,gI2)))*H0(p,MCha1(gI1),MCha1(gI2));
         tmp_2690 += 4*B0(p,MCha1(gI1),MCha1(gI2))*MCha1(gI1)*MCha1(gI2)*
            Re(Conj(CpVZbarCha1Cha1PL(gI1,gI2))*CpVZbarCha1Cha1PR(gI1,gI2));
      }
      tmp_2689 += tmp_2690;
   }
   result += tmp_2689;
   std::complex<double> tmp_2691;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2692;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2692 += (AbsSqr(CpVZbarCha2Cha2PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha2Cha2PR(gI1,gI2)))*H0(p,MCha2(gI1),MCha2(gI2));
         tmp_2692 += 4*B0(p,MCha2(gI1),MCha2(gI2))*MCha2(gI1)*MCha2(gI2)*
            Re(Conj(CpVZbarCha2Cha2PL(gI1,gI2))*CpVZbarCha2Cha2PR(gI1,gI2));
      }
      tmp_2691 += tmp_2692;
   }
   result += tmp_2691;
   std::complex<double> tmp_2693;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2693 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_2693;
   std::complex<double> tmp_2694;
   std::complex<double> tmp_2695;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2696;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2696 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_2695 += tmp_2696;
   }
   tmp_2694 += tmp_2695;
   result += (-4) * tmp_2694;
   std::complex<double> tmp_2697;
   std::complex<double> tmp_2698;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2699;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2699 += (AbsSqr(CpVZbarFdFdPL(gI1,gI2)) + AbsSqr(
            CpVZbarFdFdPR(gI1,gI2)))*H0(p,MFd(gI1),MFd(gI2));
         tmp_2699 += 4*B0(p,MFd(gI1),MFd(gI2))*MFd(gI1)*MFd(gI2)*Re(Conj(
            CpVZbarFdFdPL(gI1,gI2))*CpVZbarFdFdPR(gI1,gI2));
      }
      tmp_2698 += tmp_2699;
   }
   tmp_2697 += tmp_2698;
   result += (3) * tmp_2697;
   std::complex<double> tmp_2700;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2701;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2701 += (AbsSqr(CpVZbarFeFePL(gI1,gI2)) + AbsSqr(
            CpVZbarFeFePR(gI1,gI2)))*H0(p,MFe(gI1),MFe(gI2));
         tmp_2701 += 4*B0(p,MFe(gI1),MFe(gI2))*MFe(gI1)*MFe(gI2)*Re(Conj(
            CpVZbarFeFePL(gI1,gI2))*CpVZbarFeFePR(gI1,gI2));
      }
      tmp_2700 += tmp_2701;
   }
   result += tmp_2700;
   std::complex<double> tmp_2702;
   std::complex<double> tmp_2703;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2704;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2704 += (AbsSqr(CpVZbarFuFuPL(gI1,gI2)) + AbsSqr(
            CpVZbarFuFuPR(gI1,gI2)))*H0(p,MFu(gI1),MFu(gI2));
         tmp_2704 += 4*B0(p,MFu(gI1),MFu(gI2))*MFu(gI1)*MFu(gI2)*Re(Conj(
            CpVZbarFuFuPL(gI1,gI2))*CpVZbarFuFuPR(gI1,gI2));
      }
      tmp_2703 += tmp_2704;
   }
   tmp_2702 += tmp_2703;
   result += (3) * tmp_2702;
   std::complex<double> tmp_2705;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2706;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2706 += (AbsSqr(CpVZbarFvFvPL(gI1,gI2)) + AbsSqr(
            CpVZbarFvFvPR(gI1,gI2)))*H0(p,MFv(gI1),MFv(gI2));
         tmp_2706 += 4*B0(p,MFv(gI1),MFv(gI2))*MFv(gI1)*MFv(gI2)*Re(Conj(
            CpVZbarFvFvPL(gI1,gI2))*CpVZbarFvFvPR(gI1,gI2));
      }
      tmp_2705 += tmp_2706;
   }
   result += tmp_2705;
   std::complex<double> tmp_2707;
   std::complex<double> tmp_2708;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2708 += A0(MAh(gI1))*CpVZVZAhAh(gI1,gI1);
   }
   tmp_2707 += tmp_2708;
   result += (0.5) * tmp_2707;
   std::complex<double> tmp_2709;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2709 += A0(MHpm(gI1))*CpVZVZconjHpmHpm(gI1,gI1);
   }
   result += tmp_2709;
   std::complex<double> tmp_2710;
   std::complex<double> tmp_2711;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2711 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_2710 += tmp_2711;
   result += (0.5) * tmp_2710;
   std::complex<double> tmp_2712;
   std::complex<double> tmp_2713;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2714;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2714 += AbsSqr(CpVZconjHpmHpm(gI1,gI2))*B00(p,MHpm(gI1),MHpm
            (gI2));
      }
      tmp_2713 += tmp_2714;
   }
   tmp_2712 += tmp_2713;
   result += (-4) * tmp_2712;
   std::complex<double> tmp_2715;
   std::complex<double> tmp_2716;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2717;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2717 += AbsSqr(CpVZhhAh(gI1,gI2))*B00(p,MAh(gI2),Mhh(gI1));
      }
      tmp_2716 += tmp_2717;
   }
   tmp_2715 += tmp_2716;
   result += (-4) * tmp_2715;
   std::complex<double> tmp_2718;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2719;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2719 += (AbsSqr(CpVZbarChiChiPL(gI1,gI2)) + AbsSqr(
            CpVZbarChiChiPR(gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_2719 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZbarChiChiPL(gI1,gI2))*CpVZbarChiChiPR(gI1,gI2));
      }
      tmp_2718 += tmp_2719;
   }
   result += tmp_2718;
   std::complex<double> tmp_2720;
   std::complex<double> tmp_2721;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2721 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_2720 += tmp_2721;
   result += (3) * tmp_2720;
   std::complex<double> tmp_2722;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2722 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_2722;
   std::complex<double> tmp_2723;
   std::complex<double> tmp_2724;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2724 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_2723 += tmp_2724;
   result += (3) * tmp_2723;
   std::complex<double> tmp_2725;
   std::complex<double> tmp_2726;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2727;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2727 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_2726 += tmp_2727;
   }
   tmp_2725 += tmp_2726;
   result += (-12) * tmp_2725;
   std::complex<double> tmp_2728;
   std::complex<double> tmp_2729;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2730;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2730 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_2729 += tmp_2730;
   }
   tmp_2728 += tmp_2729;
   result += (-4) * tmp_2728;
   std::complex<double> tmp_2731;
   std::complex<double> tmp_2732;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2733;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2733 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_2732 += tmp_2733;
   }
   tmp_2731 += tmp_2732;
   result += (-12) * tmp_2731;
   std::complex<double> tmp_2734;
   std::complex<double> tmp_2735;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2735 += AbsSqr(CpVZconjVWmHpm(gI2))*B0(p,MVWm,MHpm(gI2));
   }
   tmp_2734 += tmp_2735;
   result += (2) * tmp_2734;
   std::complex<double> tmp_2736;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2736 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_2736;
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
   std::complex<double> tmp_2737;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2737 += A0(MRh(gI1))*CpVWmconjVWmconjRhRh(gI1,gI1);
   }
   result += tmp_2737;
   std::complex<double> tmp_2738;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_2738 += A0(MRpm(gI1))*CpVWmconjVWmconjRpmRpm(gI1,gI1);
   }
   result += tmp_2738;
   std::complex<double> tmp_2739;
   std::complex<double> tmp_2740;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2741;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2741 += AbsSqr(CpconjVWmconjRhRpm(gI1,gI2))*B00(p,MRpm(gI2),
            MRh(gI1));
      }
      tmp_2740 += tmp_2741;
   }
   tmp_2739 += tmp_2740;
   result += (-4) * tmp_2739;
   std::complex<double> tmp_2742;
   std::complex<double> tmp_2743;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2744;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2744 += AbsSqr(CpconjVWmRpmRh(gI1,gI2))*B00(p,MRh(gI2),MRpm(
            gI1));
      }
      tmp_2743 += tmp_2744;
   }
   tmp_2742 += tmp_2743;
   result += (-4) * tmp_2742;
   std::complex<double> tmp_2745;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2746;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2746 += (AbsSqr(CpconjVWmbarCha1ChiPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarCha1ChiPR(gI1,gI2)))*H0(p,MCha1(gI1),MChi(gI2));
         tmp_2746 += 4*B0(p,MCha1(gI1),MChi(gI2))*MCha1(gI1)*MChi(gI2)*Re
            (Conj(CpconjVWmbarCha1ChiPL(gI1,gI2))*CpconjVWmbarCha1ChiPR(gI1,gI2));
      }
      tmp_2745 += tmp_2746;
   }
   result += tmp_2745;
   std::complex<double> tmp_2747;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_2747 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_2747;
   std::complex<double> tmp_2748;
   std::complex<double> tmp_2749;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2750;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2750 += (AbsSqr(CpconjVWmbarFuFdPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFuFdPR(gI1,gI2)))*H0(p,MFu(gI1),MFd(gI2));
         tmp_2750 += 4*B0(p,MFu(gI1),MFd(gI2))*MFd(gI2)*MFu(gI1)*Re(Conj(
            CpconjVWmbarFuFdPL(gI1,gI2))*CpconjVWmbarFuFdPR(gI1,gI2));
      }
      tmp_2749 += tmp_2750;
   }
   tmp_2748 += tmp_2749;
   result += (3) * tmp_2748;
   std::complex<double> tmp_2751;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2752;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2752 += (AbsSqr(CpconjVWmbarFvFePL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarFvFePR(gI1,gI2)))*H0(p,MFv(gI1),MFe(gI2));
         tmp_2752 += 4*B0(p,MFv(gI1),MFe(gI2))*MFe(gI2)*MFv(gI1)*Re(Conj(
            CpconjVWmbarFvFePL(gI1,gI2))*CpconjVWmbarFvFePR(gI1,gI2));
      }
      tmp_2751 += tmp_2752;
   }
   result += tmp_2751;
   std::complex<double> tmp_2753;
   std::complex<double> tmp_2754;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2755;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2755 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_2754 += tmp_2755;
   }
   tmp_2753 += tmp_2754;
   result += (-4) * tmp_2753;
   std::complex<double> tmp_2756;
   std::complex<double> tmp_2757;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2757 += A0(MAh(gI1))*CpVWmconjVWmAhAh(gI1,gI1);
   }
   tmp_2756 += tmp_2757;
   result += (0.5) * tmp_2756;
   std::complex<double> tmp_2758;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2758 += A0(MHpm(gI1))*CpVWmconjVWmconjHpmHpm(gI1,gI1);
   }
   result += tmp_2758;
   std::complex<double> tmp_2759;
   std::complex<double> tmp_2760;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_2760 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_2759 += tmp_2760;
   result += (0.5) * tmp_2759;
   std::complex<double> tmp_2761;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2762;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2762 += (AbsSqr(CpconjVWmbarChiCha2PL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarChiCha2PR(gI1,gI2)))*H0(p,MChi(gI1),MCha2(gI2));
         tmp_2762 += 4*B0(p,MChi(gI1),MCha2(gI2))*MCha2(gI2)*MChi(gI1)*Re
            (Conj(CpconjVWmbarChiCha2PL(gI1,gI2))*CpconjVWmbarChiCha2PR(gI1,gI2));
      }
      tmp_2761 += tmp_2762;
   }
   result += tmp_2761;
   std::complex<double> tmp_2763;
   std::complex<double> tmp_2764;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2765;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2765 += AbsSqr(CpconjVWmHpmAh(gI1,gI2))*B00(p,MAh(gI2),MHpm(
            gI1));
      }
      tmp_2764 += tmp_2765;
   }
   tmp_2763 += tmp_2764;
   result += (-4) * tmp_2763;
   std::complex<double> tmp_2766;
   std::complex<double> tmp_2767;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2768;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2768 += AbsSqr(CpconjVWmHpmhh(gI1,gI2))*B00(p,Mhh(gI2),MHpm(
            gI1));
      }
      tmp_2767 += tmp_2768;
   }
   tmp_2766 += tmp_2767;
   result += (-4) * tmp_2766;
   std::complex<double> tmp_2769;
   std::complex<double> tmp_2770;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2770 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_2769 += tmp_2770;
   result += (3) * tmp_2769;
   std::complex<double> tmp_2771;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2771 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_2771;
   std::complex<double> tmp_2772;
   std::complex<double> tmp_2773;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_2773 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_2772 += tmp_2773;
   result += (3) * tmp_2772;
   std::complex<double> tmp_2774;
   std::complex<double> tmp_2775;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2776;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2776 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_2775 += tmp_2776;
   }
   tmp_2774 += tmp_2775;
   result += (-12) * tmp_2774;
   std::complex<double> tmp_2777;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2777 += AbsSqr(CpconjVWmVPHpm(gI2))*B0(p,0,MHpm(gI2));
   }
   result += tmp_2777;
   std::complex<double> tmp_2778;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2778 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_2778;
   std::complex<double> tmp_2779;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2779 += AbsSqr(CpconjVWmVZHpm(gI2))*B0(p,MVZ,MHpm(gI2));
   }
   result += tmp_2779;
   result += -(AbsSqr(CpconjVWmVWmVP())*(A0(MVWm) + 10*B00(p,MVWm,0) + B0(p,
      MVWm,0)*(Sqr(MVWm) + 4*Sqr(p))));
   result += -(AbsSqr(CpconjVWmVZVWm())*(A0(MVWm) + A0(MVZ) + 10*B00(p,MVZ,MVWm
      ) + B0(p,MVZ,MVWm)*(Sqr(MVWm) + Sqr(MVZ) + 4*Sqr(p))));

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2780;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2781;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2781 += B0(p,MCha1(gI2),MRpm(gI1))*Conj(
            CpbarUChiconjRpmbarCha1PL(gO2,gI1,gI2))*CpbarUChiconjRpmbarCha1PR(gO1,
            gI1,gI2)*MCha1(gI2);
      }
      tmp_2780 += tmp_2781;
   }
   result += tmp_2780;
   std::complex<double> tmp_2782;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2783;
      std::complex<double> tmp_2784;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2784 += B0(p,MCha2(gI1),MRpm(gI2))*Conj(
            CpbarUChibarCha2RpmPL(gO2,gI1,gI2))*CpbarUChibarCha2RpmPR(gO1,gI1,gI2)
            ;
      }
      tmp_2783 += tmp_2784;
      tmp_2782 += (MCha2(gI1)) * tmp_2783;
   }
   result += tmp_2782;
   std::complex<double> tmp_2785;
   std::complex<double> tmp_2786;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2787;
      std::complex<double> tmp_2788;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2788 += B0(p,MFd(gI1),MSd(gI2))*Conj(CpbarUChibarFdSdPL(gO2,
            gI1,gI2))*CpbarUChibarFdSdPR(gO1,gI1,gI2);
      }
      tmp_2787 += tmp_2788;
      tmp_2786 += (MFd(gI1)) * tmp_2787;
   }
   tmp_2785 += tmp_2786;
   result += (3) * tmp_2785;
   std::complex<double> tmp_2789;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2790;
      std::complex<double> tmp_2791;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2791 += B0(p,MFe(gI1),MSe(gI2))*Conj(CpbarUChibarFeSePL(gO2,
            gI1,gI2))*CpbarUChibarFeSePR(gO1,gI1,gI2);
      }
      tmp_2790 += tmp_2791;
      tmp_2789 += (MFe(gI1)) * tmp_2790;
   }
   result += tmp_2789;
   std::complex<double> tmp_2792;
   std::complex<double> tmp_2793;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2794;
      std::complex<double> tmp_2795;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2795 += B0(p,MFu(gI1),MSu(gI2))*Conj(CpbarUChibarFuSuPL(gO2,
            gI1,gI2))*CpbarUChibarFuSuPR(gO1,gI1,gI2);
      }
      tmp_2794 += tmp_2795;
      tmp_2793 += (MFu(gI1)) * tmp_2794;
   }
   tmp_2792 += tmp_2793;
   result += (3) * tmp_2792;
   std::complex<double> tmp_2796;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2797;
      std::complex<double> tmp_2798;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2798 += B0(p,MFv(gI1),MSv(gI2))*Conj(CpbarUChibarFvSvPL(gO2,
            gI1,gI2))*CpbarUChibarFvSvPR(gO1,gI1,gI2);
      }
      tmp_2797 += tmp_2798;
      tmp_2796 += (MFv(gI1)) * tmp_2797;
   }
   result += tmp_2796;
   std::complex<double> tmp_2799;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2800;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2800 += B0(p,MCha1(gI2),MHpm(gI1))*Conj(CpbarUChiHpmCha1PL(
            gO2,gI1,gI2))*CpbarUChiHpmCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_2799 += tmp_2800;
   }
   result += tmp_2799;
   std::complex<double> tmp_2801;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2802;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2802 += B0(p,MCha2(gI2),MHpm(gI1))*Conj(
            CpbarUChiconjHpmCha2PL(gO2,gI1,gI2))*CpbarUChiconjHpmCha2PR(gO1,gI1,
            gI2)*MCha2(gI2);
      }
      tmp_2801 += tmp_2802;
   }
   result += tmp_2801;
   std::complex<double> tmp_2803;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2804;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2804 += B0(p,MChi(gI2),Mhh(gI1))*Conj(CpbarUChihhChiPL(gO2,
            gI1,gI2))*CpbarUChihhChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_2803 += tmp_2804;
   }
   result += tmp_2803;
   std::complex<double> tmp_2805;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2806;
      std::complex<double> tmp_2807;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2807 += B0(p,MChi(gI1),MRh(gI2))*Conj(CpbarUChibarChiRhPL(
            gO2,gI1,gI2))*CpbarUChibarChiRhPR(gO1,gI1,gI2);
      }
      tmp_2806 += tmp_2807;
      tmp_2805 += (MChi(gI1)) * tmp_2806;
   }
   result += tmp_2805;
   std::complex<double> tmp_2808;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2809;
      std::complex<double> tmp_2810;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2810 += B0(p,MChi(gI1),MAh(gI2))*Conj(CpbarUChiChiAhPL(gO2,
            gI1,gI2))*CpbarUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2809 += tmp_2810;
      tmp_2808 += (MChi(gI1)) * tmp_2809;
   }
   result += tmp_2808;
   std::complex<double> tmp_2811;
   std::complex<double> tmp_2812;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2813;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2813 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpbarUChiconjSdFdPL(gO2
            ,gI1,gI2))*CpbarUChiconjSdFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_2812 += tmp_2813;
   }
   tmp_2811 += tmp_2812;
   result += (3) * tmp_2811;
   std::complex<double> tmp_2814;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2815;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2815 += B0(p,MFe(gI2),MSe(gI1))*Conj(CpbarUChiconjSeFePL(gO2
            ,gI1,gI2))*CpbarUChiconjSeFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_2814 += tmp_2815;
   }
   result += tmp_2814;
   std::complex<double> tmp_2816;
   std::complex<double> tmp_2817;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2818;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2818 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpbarUChiconjSuFuPL(gO2
            ,gI1,gI2))*CpbarUChiconjSuFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2817 += tmp_2818;
   }
   tmp_2816 += tmp_2817;
   result += (3) * tmp_2816;
   std::complex<double> tmp_2819;
   std::complex<double> tmp_2820;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2820 += B0(p,MCha1(gI2),MVWm)*Conj(CpbarUChiVWmCha1PR(gO2,gI2))*
         CpbarUChiVWmCha1PL(gO1,gI2)*MCha1(gI2);
   }
   tmp_2819 += tmp_2820;
   result += (-4) * tmp_2819;
   std::complex<double> tmp_2821;
   std::complex<double> tmp_2822;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2822 += B0(p,MCha2(gI2),MVWm)*Conj(CpbarUChiconjVWmCha2PR(gO2,gI2)
         )*CpbarUChiconjVWmCha2PL(gO1,gI2)*MCha2(gI2);
   }
   tmp_2821 += tmp_2822;
   result += (-4) * tmp_2821;
   std::complex<double> tmp_2823;
   std::complex<double> tmp_2824;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2824 += B0(p,MChi(gI2),MVZ)*Conj(CpbarUChiVZChiPR(gO2,gI2))*
         CpbarUChiVZChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2823 += tmp_2824;
   result += (-4) * tmp_2823;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2825;
   std::complex<double> tmp_2826;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2827;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2827 += B1(p,MCha2(gI1),MRpm(gI2))*Conj(
            CpbarUChibarCha2RpmPR(gO2,gI1,gI2))*CpbarUChibarCha2RpmPR(gO1,gI1,gI2)
            ;
      }
      tmp_2826 += tmp_2827;
   }
   tmp_2825 += tmp_2826;
   result += (-0.5) * tmp_2825;
   std::complex<double> tmp_2828;
   std::complex<double> tmp_2829;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2830;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2830 += B1(p,MCha1(gI2),MRpm(gI1))*Conj(
            CpbarUChiconjRpmbarCha1PR(gO2,gI1,gI2))*CpbarUChiconjRpmbarCha1PR(gO1,
            gI1,gI2);
      }
      tmp_2829 += tmp_2830;
   }
   tmp_2828 += tmp_2829;
   result += (-0.5) * tmp_2828;
   std::complex<double> tmp_2831;
   std::complex<double> tmp_2832;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2833;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2833 += B1(p,MFv(gI1),MSv(gI2))*Conj(CpbarUChibarFvSvPR(gO2,
            gI1,gI2))*CpbarUChibarFvSvPR(gO1,gI1,gI2);
      }
      tmp_2832 += tmp_2833;
   }
   tmp_2831 += tmp_2832;
   result += (-0.5) * tmp_2831;
   std::complex<double> tmp_2834;
   std::complex<double> tmp_2835;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2836;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2836 += B1(p,MFd(gI1),MSd(gI2))*Conj(CpbarUChibarFdSdPR(gO2,
            gI1,gI2))*CpbarUChibarFdSdPR(gO1,gI1,gI2);
      }
      tmp_2835 += tmp_2836;
   }
   tmp_2834 += tmp_2835;
   result += (-1.5) * tmp_2834;
   std::complex<double> tmp_2837;
   std::complex<double> tmp_2838;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2839;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2839 += B1(p,MFe(gI1),MSe(gI2))*Conj(CpbarUChibarFeSePR(gO2,
            gI1,gI2))*CpbarUChibarFeSePR(gO1,gI1,gI2);
      }
      tmp_2838 += tmp_2839;
   }
   tmp_2837 += tmp_2838;
   result += (-0.5) * tmp_2837;
   std::complex<double> tmp_2840;
   std::complex<double> tmp_2841;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2842;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2842 += B1(p,MFu(gI1),MSu(gI2))*Conj(CpbarUChibarFuSuPR(gO2,
            gI1,gI2))*CpbarUChibarFuSuPR(gO1,gI1,gI2);
      }
      tmp_2841 += tmp_2842;
   }
   tmp_2840 += tmp_2841;
   result += (-1.5) * tmp_2840;
   std::complex<double> tmp_2843;
   std::complex<double> tmp_2844;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2845;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2845 += B1(p,MChi(gI1),MRh(gI2))*Conj(CpbarUChibarChiRhPR(
            gO2,gI1,gI2))*CpbarUChibarChiRhPR(gO1,gI1,gI2);
      }
      tmp_2844 += tmp_2845;
   }
   tmp_2843 += tmp_2844;
   result += (-0.5) * tmp_2843;
   std::complex<double> tmp_2846;
   std::complex<double> tmp_2847;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2848;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2848 += B1(p,MCha2(gI2),MHpm(gI1))*Conj(
            CpbarUChiconjHpmCha2PR(gO2,gI1,gI2))*CpbarUChiconjHpmCha2PR(gO1,gI1,
            gI2);
      }
      tmp_2847 += tmp_2848;
   }
   tmp_2846 += tmp_2847;
   result += (-0.5) * tmp_2846;
   std::complex<double> tmp_2849;
   std::complex<double> tmp_2850;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2851;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2851 += B1(p,MCha1(gI2),MHpm(gI1))*Conj(CpbarUChiHpmCha1PR(
            gO2,gI1,gI2))*CpbarUChiHpmCha1PR(gO1,gI1,gI2);
      }
      tmp_2850 += tmp_2851;
   }
   tmp_2849 += tmp_2850;
   result += (-0.5) * tmp_2849;
   std::complex<double> tmp_2852;
   std::complex<double> tmp_2853;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2854;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2854 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpbarUChiChiAhPR(gO2,
            gI1,gI2))*CpbarUChiChiAhPR(gO1,gI1,gI2);
      }
      tmp_2853 += tmp_2854;
   }
   tmp_2852 += tmp_2853;
   result += (-0.5) * tmp_2852;
   std::complex<double> tmp_2855;
   std::complex<double> tmp_2856;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2857;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2857 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpbarUChihhChiPR(gO2,
            gI1,gI2))*CpbarUChihhChiPR(gO1,gI1,gI2);
      }
      tmp_2856 += tmp_2857;
   }
   tmp_2855 += tmp_2856;
   result += (-0.5) * tmp_2855;
   std::complex<double> tmp_2858;
   std::complex<double> tmp_2859;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2860;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2860 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpbarUChiconjSdFdPR(gO2
            ,gI1,gI2))*CpbarUChiconjSdFdPR(gO1,gI1,gI2);
      }
      tmp_2859 += tmp_2860;
   }
   tmp_2858 += tmp_2859;
   result += (-1.5) * tmp_2858;
   std::complex<double> tmp_2861;
   std::complex<double> tmp_2862;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2863;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2863 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpbarUChiconjSeFePR(gO2
            ,gI1,gI2))*CpbarUChiconjSeFePR(gO1,gI1,gI2);
      }
      tmp_2862 += tmp_2863;
   }
   tmp_2861 += tmp_2862;
   result += (-0.5) * tmp_2861;
   std::complex<double> tmp_2864;
   std::complex<double> tmp_2865;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2866;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2866 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpbarUChiconjSuFuPR(gO2
            ,gI1,gI2))*CpbarUChiconjSuFuPR(gO1,gI1,gI2);
      }
      tmp_2865 += tmp_2866;
   }
   tmp_2864 += tmp_2865;
   result += (-1.5) * tmp_2864;
   std::complex<double> tmp_2867;
   std::complex<double> tmp_2868;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2868 += B1(p,MCha2(gI2),MVWm)*Conj(CpbarUChiconjVWmCha2PL(gO2,gI2)
         )*CpbarUChiconjVWmCha2PL(gO1,gI2);
   }
   tmp_2867 += tmp_2868;
   result += (-1) * tmp_2867;
   std::complex<double> tmp_2869;
   std::complex<double> tmp_2870;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2870 += B1(p,MCha1(gI2),MVWm)*Conj(CpbarUChiVWmCha1PL(gO2,gI2))*
         CpbarUChiVWmCha1PL(gO1,gI2);
   }
   tmp_2869 += tmp_2870;
   result += (-1) * tmp_2869;
   std::complex<double> tmp_2871;
   std::complex<double> tmp_2872;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2872 += B1(p,MChi(gI2),MVZ)*Conj(CpbarUChiVZChiPL(gO2,gI2))*
         CpbarUChiVZChiPL(gO1,gI2);
   }
   tmp_2871 += tmp_2872;
   result += (-1) * tmp_2871;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Chi_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2873;
   std::complex<double> tmp_2874;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2875;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2875 += B1(p,MCha2(gI1),MRpm(gI2))*Conj(
            CpbarUChibarCha2RpmPL(gO2,gI1,gI2))*CpbarUChibarCha2RpmPL(gO1,gI1,gI2)
            ;
      }
      tmp_2874 += tmp_2875;
   }
   tmp_2873 += tmp_2874;
   result += (-0.5) * tmp_2873;
   std::complex<double> tmp_2876;
   std::complex<double> tmp_2877;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2878;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2878 += B1(p,MCha1(gI2),MRpm(gI1))*Conj(
            CpbarUChiconjRpmbarCha1PL(gO2,gI1,gI2))*CpbarUChiconjRpmbarCha1PL(gO1,
            gI1,gI2);
      }
      tmp_2877 += tmp_2878;
   }
   tmp_2876 += tmp_2877;
   result += (-0.5) * tmp_2876;
   std::complex<double> tmp_2879;
   std::complex<double> tmp_2880;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2881;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2881 += B1(p,MFv(gI1),MSv(gI2))*Conj(CpbarUChibarFvSvPL(gO2,
            gI1,gI2))*CpbarUChibarFvSvPL(gO1,gI1,gI2);
      }
      tmp_2880 += tmp_2881;
   }
   tmp_2879 += tmp_2880;
   result += (-0.5) * tmp_2879;
   std::complex<double> tmp_2882;
   std::complex<double> tmp_2883;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2884;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2884 += B1(p,MFd(gI1),MSd(gI2))*Conj(CpbarUChibarFdSdPL(gO2,
            gI1,gI2))*CpbarUChibarFdSdPL(gO1,gI1,gI2);
      }
      tmp_2883 += tmp_2884;
   }
   tmp_2882 += tmp_2883;
   result += (-1.5) * tmp_2882;
   std::complex<double> tmp_2885;
   std::complex<double> tmp_2886;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2887;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2887 += B1(p,MFe(gI1),MSe(gI2))*Conj(CpbarUChibarFeSePL(gO2,
            gI1,gI2))*CpbarUChibarFeSePL(gO1,gI1,gI2);
      }
      tmp_2886 += tmp_2887;
   }
   tmp_2885 += tmp_2886;
   result += (-0.5) * tmp_2885;
   std::complex<double> tmp_2888;
   std::complex<double> tmp_2889;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2890;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2890 += B1(p,MFu(gI1),MSu(gI2))*Conj(CpbarUChibarFuSuPL(gO2,
            gI1,gI2))*CpbarUChibarFuSuPL(gO1,gI1,gI2);
      }
      tmp_2889 += tmp_2890;
   }
   tmp_2888 += tmp_2889;
   result += (-1.5) * tmp_2888;
   std::complex<double> tmp_2891;
   std::complex<double> tmp_2892;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2893;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2893 += B1(p,MChi(gI1),MRh(gI2))*Conj(CpbarUChibarChiRhPL(
            gO2,gI1,gI2))*CpbarUChibarChiRhPL(gO1,gI1,gI2);
      }
      tmp_2892 += tmp_2893;
   }
   tmp_2891 += tmp_2892;
   result += (-0.5) * tmp_2891;
   std::complex<double> tmp_2894;
   std::complex<double> tmp_2895;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2896;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2896 += B1(p,MCha2(gI2),MHpm(gI1))*Conj(
            CpbarUChiconjHpmCha2PL(gO2,gI1,gI2))*CpbarUChiconjHpmCha2PL(gO1,gI1,
            gI2);
      }
      tmp_2895 += tmp_2896;
   }
   tmp_2894 += tmp_2895;
   result += (-0.5) * tmp_2894;
   std::complex<double> tmp_2897;
   std::complex<double> tmp_2898;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2899;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2899 += B1(p,MCha1(gI2),MHpm(gI1))*Conj(CpbarUChiHpmCha1PL(
            gO2,gI1,gI2))*CpbarUChiHpmCha1PL(gO1,gI1,gI2);
      }
      tmp_2898 += tmp_2899;
   }
   tmp_2897 += tmp_2898;
   result += (-0.5) * tmp_2897;
   std::complex<double> tmp_2900;
   std::complex<double> tmp_2901;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2902;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2902 += B1(p,MChi(gI1),MAh(gI2))*Conj(CpbarUChiChiAhPL(gO2,
            gI1,gI2))*CpbarUChiChiAhPL(gO1,gI1,gI2);
      }
      tmp_2901 += tmp_2902;
   }
   tmp_2900 += tmp_2901;
   result += (-0.5) * tmp_2900;
   std::complex<double> tmp_2903;
   std::complex<double> tmp_2904;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2905;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2905 += B1(p,MChi(gI2),Mhh(gI1))*Conj(CpbarUChihhChiPL(gO2,
            gI1,gI2))*CpbarUChihhChiPL(gO1,gI1,gI2);
      }
      tmp_2904 += tmp_2905;
   }
   tmp_2903 += tmp_2904;
   result += (-0.5) * tmp_2903;
   std::complex<double> tmp_2906;
   std::complex<double> tmp_2907;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2908;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2908 += B1(p,MFd(gI2),MSd(gI1))*Conj(CpbarUChiconjSdFdPL(gO2
            ,gI1,gI2))*CpbarUChiconjSdFdPL(gO1,gI1,gI2);
      }
      tmp_2907 += tmp_2908;
   }
   tmp_2906 += tmp_2907;
   result += (-1.5) * tmp_2906;
   std::complex<double> tmp_2909;
   std::complex<double> tmp_2910;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2911;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2911 += B1(p,MFe(gI2),MSe(gI1))*Conj(CpbarUChiconjSeFePL(gO2
            ,gI1,gI2))*CpbarUChiconjSeFePL(gO1,gI1,gI2);
      }
      tmp_2910 += tmp_2911;
   }
   tmp_2909 += tmp_2910;
   result += (-0.5) * tmp_2909;
   std::complex<double> tmp_2912;
   std::complex<double> tmp_2913;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2914;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2914 += B1(p,MFu(gI2),MSu(gI1))*Conj(CpbarUChiconjSuFuPL(gO2
            ,gI1,gI2))*CpbarUChiconjSuFuPL(gO1,gI1,gI2);
      }
      tmp_2913 += tmp_2914;
   }
   tmp_2912 += tmp_2913;
   result += (-1.5) * tmp_2912;
   std::complex<double> tmp_2915;
   std::complex<double> tmp_2916;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2916 += B1(p,MCha2(gI2),MVWm)*Conj(CpbarUChiconjVWmCha2PR(gO2,gI2)
         )*CpbarUChiconjVWmCha2PR(gO1,gI2);
   }
   tmp_2915 += tmp_2916;
   result += (-1) * tmp_2915;
   std::complex<double> tmp_2917;
   std::complex<double> tmp_2918;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2918 += B1(p,MCha1(gI2),MVWm)*Conj(CpbarUChiVWmCha1PR(gO2,gI2))*
         CpbarUChiVWmCha1PR(gO1,gI2);
   }
   tmp_2917 += tmp_2918;
   result += (-1) * tmp_2917;
   std::complex<double> tmp_2919;
   std::complex<double> tmp_2920;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2920 += B1(p,MChi(gI2),MVZ)*Conj(CpbarUChiVZChiPR(gO2,gI2))*
         CpbarUChiVZChiPR(gO1,gI2);
   }
   tmp_2919 += tmp_2920;
   result += (-1) * tmp_2919;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha1_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2921;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2922;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2922 += B0(p,MChi(gI2),MRpm(gI1))*Conj(
            CpbarUCha1conjRpmbarChiPL(gO2,gI1,gI2))*CpbarUCha1conjRpmbarChiPR(gO1,
            gI1,gI2)*MChi(gI2);
      }
      tmp_2921 += tmp_2922;
   }
   result += tmp_2921;
   std::complex<double> tmp_2923;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2924;
      std::complex<double> tmp_2925;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2925 += B0(p,MCha1(gI1),MAh(gI2))*Conj(CpbarUCha1Cha1AhPL(
            gO2,gI1,gI2))*CpbarUCha1Cha1AhPR(gO1,gI1,gI2);
      }
      tmp_2924 += tmp_2925;
      tmp_2923 += (MCha1(gI1)) * tmp_2924;
   }
   result += tmp_2923;
   std::complex<double> tmp_2926;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2927;
      std::complex<double> tmp_2928;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2928 += B0(p,MCha2(gI1),MRh(gI2))*Conj(CpbarUCha1barCha2RhPL
            (gO2,gI1,gI2))*CpbarUCha1barCha2RhPR(gO1,gI1,gI2);
      }
      tmp_2927 += tmp_2928;
      tmp_2926 += (MCha2(gI1)) * tmp_2927;
   }
   result += tmp_2926;
   std::complex<double> tmp_2929;
   std::complex<double> tmp_2930;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2931;
      std::complex<double> tmp_2932;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2932 += B0(p,MFd(gI1),MSu(gI2))*Conj(CpbarUCha1barFdSuPL(gO2
            ,gI1,gI2))*CpbarUCha1barFdSuPR(gO1,gI1,gI2);
      }
      tmp_2931 += tmp_2932;
      tmp_2930 += (MFd(gI1)) * tmp_2931;
   }
   tmp_2929 += tmp_2930;
   result += (3) * tmp_2929;
   std::complex<double> tmp_2933;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2934;
      std::complex<double> tmp_2935;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2935 += B0(p,MFe(gI1),MSv(gI2))*Conj(CpbarUCha1barFeSvPL(gO2
            ,gI1,gI2))*CpbarUCha1barFeSvPR(gO1,gI1,gI2);
      }
      tmp_2934 += tmp_2935;
      tmp_2933 += (MFe(gI1)) * tmp_2934;
   }
   result += tmp_2933;
   std::complex<double> tmp_2936;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2937;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2937 += B0(p,MCha1(gI2),Mhh(gI1))*Conj(CpbarUCha1hhCha1PL(
            gO2,gI1,gI2))*CpbarUCha1hhCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_2936 += tmp_2937;
   }
   result += tmp_2936;
   std::complex<double> tmp_2938;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2939;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2939 += B0(p,MChi(gI2),MHpm(gI1))*Conj(
            CpbarUCha1conjHpmChiPL(gO2,gI1,gI2))*CpbarUCha1conjHpmChiPR(gO1,gI1,
            gI2)*MChi(gI2);
      }
      tmp_2938 += tmp_2939;
   }
   result += tmp_2938;
   std::complex<double> tmp_2940;
   std::complex<double> tmp_2941;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2942;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2942 += B0(p,MFu(gI2),MSd(gI1))*Conj(CpbarUCha1conjSdFuPL(
            gO2,gI1,gI2))*CpbarUCha1conjSdFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_2941 += tmp_2942;
   }
   tmp_2940 += tmp_2941;
   result += (3) * tmp_2940;
   std::complex<double> tmp_2943;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2944;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2944 += B0(p,MFv(gI2),MSe(gI1))*Conj(CpbarUCha1conjSeFvPL(
            gO2,gI1,gI2))*CpbarUCha1conjSeFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_2943 += tmp_2944;
   }
   result += tmp_2943;
   std::complex<double> tmp_2945;
   std::complex<double> tmp_2946;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2946 += B0(p,MCha1(gI2),0)*Conj(CpbarUCha1VPCha1PR(gO2,gI2))*
         CpbarUCha1VPCha1PL(gO1,gI2)*MCha1(gI2);
   }
   tmp_2945 += tmp_2946;
   result += (-4) * tmp_2945;
   std::complex<double> tmp_2947;
   std::complex<double> tmp_2948;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2948 += B0(p,MCha1(gI2),MVZ)*Conj(CpbarUCha1VZCha1PR(gO2,gI2))*
         CpbarUCha1VZCha1PL(gO1,gI2)*MCha1(gI2);
   }
   tmp_2947 += tmp_2948;
   result += (-4) * tmp_2947;
   std::complex<double> tmp_2949;
   std::complex<double> tmp_2950;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2950 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUCha1conjVWmChiPR(gO2,gI2))
         *CpbarUCha1conjVWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_2949 += tmp_2950;
   result += (-4) * tmp_2949;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha1_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2951;
   std::complex<double> tmp_2952;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2953;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2953 += B1(p,MCha2(gI1),MRh(gI2))*Conj(CpbarUCha1barCha2RhPR
            (gO2,gI1,gI2))*CpbarUCha1barCha2RhPR(gO1,gI1,gI2);
      }
      tmp_2952 += tmp_2953;
   }
   tmp_2951 += tmp_2952;
   result += (-0.5) * tmp_2951;
   std::complex<double> tmp_2954;
   std::complex<double> tmp_2955;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2956;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2956 += B1(p,MCha1(gI1),MAh(gI2))*Conj(CpbarUCha1Cha1AhPR(
            gO2,gI1,gI2))*CpbarUCha1Cha1AhPR(gO1,gI1,gI2);
      }
      tmp_2955 += tmp_2956;
   }
   tmp_2954 += tmp_2955;
   result += (-0.5) * tmp_2954;
   std::complex<double> tmp_2957;
   std::complex<double> tmp_2958;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2959;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2959 += B1(p,MChi(gI2),MRpm(gI1))*Conj(
            CpbarUCha1conjRpmbarChiPR(gO2,gI1,gI2))*CpbarUCha1conjRpmbarChiPR(gO1,
            gI1,gI2);
      }
      tmp_2958 += tmp_2959;
   }
   tmp_2957 += tmp_2958;
   result += (-0.5) * tmp_2957;
   std::complex<double> tmp_2960;
   std::complex<double> tmp_2961;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2962;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2962 += B1(p,MFe(gI1),MSv(gI2))*Conj(CpbarUCha1barFeSvPR(gO2
            ,gI1,gI2))*CpbarUCha1barFeSvPR(gO1,gI1,gI2);
      }
      tmp_2961 += tmp_2962;
   }
   tmp_2960 += tmp_2961;
   result += (-0.5) * tmp_2960;
   std::complex<double> tmp_2963;
   std::complex<double> tmp_2964;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2965;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2965 += B1(p,MFd(gI1),MSu(gI2))*Conj(CpbarUCha1barFdSuPR(gO2
            ,gI1,gI2))*CpbarUCha1barFdSuPR(gO1,gI1,gI2);
      }
      tmp_2964 += tmp_2965;
   }
   tmp_2963 += tmp_2964;
   result += (-1.5) * tmp_2963;
   std::complex<double> tmp_2966;
   std::complex<double> tmp_2967;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2968;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2968 += B1(p,MCha1(gI2),Mhh(gI1))*Conj(CpbarUCha1hhCha1PR(
            gO2,gI1,gI2))*CpbarUCha1hhCha1PR(gO1,gI1,gI2);
      }
      tmp_2967 += tmp_2968;
   }
   tmp_2966 += tmp_2967;
   result += (-0.5) * tmp_2966;
   std::complex<double> tmp_2969;
   std::complex<double> tmp_2970;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_2971;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2971 += B1(p,MChi(gI2),MHpm(gI1))*Conj(
            CpbarUCha1conjHpmChiPR(gO2,gI1,gI2))*CpbarUCha1conjHpmChiPR(gO1,gI1,
            gI2);
      }
      tmp_2970 += tmp_2971;
   }
   tmp_2969 += tmp_2970;
   result += (-0.5) * tmp_2969;
   std::complex<double> tmp_2972;
   std::complex<double> tmp_2973;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2974;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2974 += B1(p,MFu(gI2),MSd(gI1))*Conj(CpbarUCha1conjSdFuPR(
            gO2,gI1,gI2))*CpbarUCha1conjSdFuPR(gO1,gI1,gI2);
      }
      tmp_2973 += tmp_2974;
   }
   tmp_2972 += tmp_2973;
   result += (-1.5) * tmp_2972;
   std::complex<double> tmp_2975;
   std::complex<double> tmp_2976;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_2977;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2977 += B1(p,MFv(gI2),MSe(gI1))*Conj(CpbarUCha1conjSeFvPR(
            gO2,gI1,gI2))*CpbarUCha1conjSeFvPR(gO1,gI1,gI2);
      }
      tmp_2976 += tmp_2977;
   }
   tmp_2975 += tmp_2976;
   result += (-0.5) * tmp_2975;
   std::complex<double> tmp_2978;
   std::complex<double> tmp_2979;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2979 += B1(p,MCha1(gI2),0)*Conj(CpbarUCha1VPCha1PL(gO2,gI2))*
         CpbarUCha1VPCha1PL(gO1,gI2);
   }
   tmp_2978 += tmp_2979;
   result += (-1) * tmp_2978;
   std::complex<double> tmp_2980;
   std::complex<double> tmp_2981;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_2981 += B1(p,MCha1(gI2),MVZ)*Conj(CpbarUCha1VZCha1PL(gO2,gI2))*
         CpbarUCha1VZCha1PL(gO1,gI2);
   }
   tmp_2980 += tmp_2981;
   result += (-1) * tmp_2980;
   std::complex<double> tmp_2982;
   std::complex<double> tmp_2983;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_2983 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha1conjVWmChiPL(gO2,gI2))
         *CpbarUCha1conjVWmChiPL(gO1,gI2);
   }
   tmp_2982 += tmp_2983;
   result += (-1) * tmp_2982;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha1_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_2984;
   std::complex<double> tmp_2985;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2986;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_2986 += B1(p,MCha2(gI1),MRh(gI2))*Conj(CpbarUCha1barCha2RhPL
            (gO2,gI1,gI2))*CpbarUCha1barCha2RhPL(gO1,gI1,gI2);
      }
      tmp_2985 += tmp_2986;
   }
   tmp_2984 += tmp_2985;
   result += (-0.5) * tmp_2984;
   std::complex<double> tmp_2987;
   std::complex<double> tmp_2988;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2989;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2989 += B1(p,MCha1(gI1),MAh(gI2))*Conj(CpbarUCha1Cha1AhPL(
            gO2,gI1,gI2))*CpbarUCha1Cha1AhPL(gO1,gI1,gI2);
      }
      tmp_2988 += tmp_2989;
   }
   tmp_2987 += tmp_2988;
   result += (-0.5) * tmp_2987;
   std::complex<double> tmp_2990;
   std::complex<double> tmp_2991;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_2992;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_2992 += B1(p,MChi(gI2),MRpm(gI1))*Conj(
            CpbarUCha1conjRpmbarChiPL(gO2,gI1,gI2))*CpbarUCha1conjRpmbarChiPL(gO1,
            gI1,gI2);
      }
      tmp_2991 += tmp_2992;
   }
   tmp_2990 += tmp_2991;
   result += (-0.5) * tmp_2990;
   std::complex<double> tmp_2993;
   std::complex<double> tmp_2994;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2995;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_2995 += B1(p,MFe(gI1),MSv(gI2))*Conj(CpbarUCha1barFeSvPL(gO2
            ,gI1,gI2))*CpbarUCha1barFeSvPL(gO1,gI1,gI2);
      }
      tmp_2994 += tmp_2995;
   }
   tmp_2993 += tmp_2994;
   result += (-0.5) * tmp_2993;
   std::complex<double> tmp_2996;
   std::complex<double> tmp_2997;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_2998;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_2998 += B1(p,MFd(gI1),MSu(gI2))*Conj(CpbarUCha1barFdSuPL(gO2
            ,gI1,gI2))*CpbarUCha1barFdSuPL(gO1,gI1,gI2);
      }
      tmp_2997 += tmp_2998;
   }
   tmp_2996 += tmp_2997;
   result += (-1.5) * tmp_2996;
   std::complex<double> tmp_2999;
   std::complex<double> tmp_3000;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3001;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3001 += B1(p,MCha1(gI2),Mhh(gI1))*Conj(CpbarUCha1hhCha1PL(
            gO2,gI1,gI2))*CpbarUCha1hhCha1PL(gO1,gI1,gI2);
      }
      tmp_3000 += tmp_3001;
   }
   tmp_2999 += tmp_3000;
   result += (-0.5) * tmp_2999;
   std::complex<double> tmp_3002;
   std::complex<double> tmp_3003;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3004;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3004 += B1(p,MChi(gI2),MHpm(gI1))*Conj(
            CpbarUCha1conjHpmChiPL(gO2,gI1,gI2))*CpbarUCha1conjHpmChiPL(gO1,gI1,
            gI2);
      }
      tmp_3003 += tmp_3004;
   }
   tmp_3002 += tmp_3003;
   result += (-0.5) * tmp_3002;
   std::complex<double> tmp_3005;
   std::complex<double> tmp_3006;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3007;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3007 += B1(p,MFu(gI2),MSd(gI1))*Conj(CpbarUCha1conjSdFuPL(
            gO2,gI1,gI2))*CpbarUCha1conjSdFuPL(gO1,gI1,gI2);
      }
      tmp_3006 += tmp_3007;
   }
   tmp_3005 += tmp_3006;
   result += (-1.5) * tmp_3005;
   std::complex<double> tmp_3008;
   std::complex<double> tmp_3009;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3010;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3010 += B1(p,MFv(gI2),MSe(gI1))*Conj(CpbarUCha1conjSeFvPL(
            gO2,gI1,gI2))*CpbarUCha1conjSeFvPL(gO1,gI1,gI2);
      }
      tmp_3009 += tmp_3010;
   }
   tmp_3008 += tmp_3009;
   result += (-0.5) * tmp_3008;
   std::complex<double> tmp_3011;
   std::complex<double> tmp_3012;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3012 += B1(p,MCha1(gI2),0)*Conj(CpbarUCha1VPCha1PR(gO2,gI2))*
         CpbarUCha1VPCha1PR(gO1,gI2);
   }
   tmp_3011 += tmp_3012;
   result += (-1) * tmp_3011;
   std::complex<double> tmp_3013;
   std::complex<double> tmp_3014;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3014 += B1(p,MCha1(gI2),MVZ)*Conj(CpbarUCha1VZCha1PR(gO2,gI2))*
         CpbarUCha1VZCha1PR(gO1,gI2);
   }
   tmp_3013 += tmp_3014;
   result += (-1) * tmp_3013;
   std::complex<double> tmp_3015;
   std::complex<double> tmp_3016;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3016 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha1conjVWmChiPR(gO2,gI2))
         *CpbarUCha1conjVWmChiPR(gO1,gI2);
   }
   tmp_3015 += tmp_3016;
   result += (-1) * tmp_3015;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha2_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3017;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3018;
      std::complex<double> tmp_3019;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3019 += B0(p,MCha1(gI1),MRh(gI2))*Conj(CpbarUCha2barCha1RhPL
            (gO2,gI1,gI2))*CpbarUCha2barCha1RhPR(gO1,gI1,gI2);
      }
      tmp_3018 += tmp_3019;
      tmp_3017 += (MCha1(gI1)) * tmp_3018;
   }
   result += tmp_3017;
   std::complex<double> tmp_3020;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3021;
      std::complex<double> tmp_3022;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3022 += B0(p,MCha2(gI1),MAh(gI2))*Conj(CpbarUCha2Cha2AhPL(
            gO2,gI1,gI2))*CpbarUCha2Cha2AhPR(gO1,gI1,gI2);
      }
      tmp_3021 += tmp_3022;
      tmp_3020 += (MCha2(gI1)) * tmp_3021;
   }
   result += tmp_3020;
   std::complex<double> tmp_3023;
   std::complex<double> tmp_3024;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3025;
      std::complex<double> tmp_3026;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3026 += B0(p,MFu(gI1),MSd(gI2))*Conj(CpbarUCha2barFuSdPL(gO2
            ,gI1,gI2))*CpbarUCha2barFuSdPR(gO1,gI1,gI2);
      }
      tmp_3025 += tmp_3026;
      tmp_3024 += (MFu(gI1)) * tmp_3025;
   }
   tmp_3023 += tmp_3024;
   result += (3) * tmp_3023;
   std::complex<double> tmp_3027;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3028;
      std::complex<double> tmp_3029;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3029 += B0(p,MFv(gI1),MSe(gI2))*Conj(CpbarUCha2barFvSePL(gO2
            ,gI1,gI2))*CpbarUCha2barFvSePR(gO1,gI1,gI2);
      }
      tmp_3028 += tmp_3029;
      tmp_3027 += (MFv(gI1)) * tmp_3028;
   }
   result += tmp_3027;
   std::complex<double> tmp_3030;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3031;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3031 += B0(p,MCha2(gI2),Mhh(gI1))*Conj(CpbarUCha2hhCha2PL(
            gO2,gI1,gI2))*CpbarUCha2hhCha2PR(gO1,gI1,gI2)*MCha2(gI2);
      }
      tmp_3030 += tmp_3031;
   }
   result += tmp_3030;
   std::complex<double> tmp_3032;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3033;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3033 += B0(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUCha2HpmChiPL(
            gO2,gI1,gI2))*CpbarUCha2HpmChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3032 += tmp_3033;
   }
   result += tmp_3032;
   std::complex<double> tmp_3034;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3035;
      std::complex<double> tmp_3036;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3036 += B0(p,MChi(gI1),MRpm(gI2))*Conj(CpbarUCha2barChiRpmPL
            (gO2,gI1,gI2))*CpbarUCha2barChiRpmPR(gO1,gI1,gI2);
      }
      tmp_3035 += tmp_3036;
      tmp_3034 += (MChi(gI1)) * tmp_3035;
   }
   result += tmp_3034;
   std::complex<double> tmp_3037;
   std::complex<double> tmp_3038;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3039;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3039 += B0(p,MFd(gI2),MSu(gI1))*Conj(CpbarUCha2conjSuFdPL(
            gO2,gI1,gI2))*CpbarUCha2conjSuFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3038 += tmp_3039;
   }
   tmp_3037 += tmp_3038;
   result += (3) * tmp_3037;
   std::complex<double> tmp_3040;
   std::complex<double> tmp_3041;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3041 += B0(p,MCha2(gI2),0)*Conj(CpbarUCha2VPCha2PR(gO2,gI2))*
         CpbarUCha2VPCha2PL(gO1,gI2)*MCha2(gI2);
   }
   tmp_3040 += tmp_3041;
   result += (-4) * tmp_3040;
   std::complex<double> tmp_3042;
   std::complex<double> tmp_3043;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3043 += B0(p,MCha2(gI2),MVZ)*Conj(CpbarUCha2VZCha2PR(gO2,gI2))*
         CpbarUCha2VZCha2PL(gO1,gI2)*MCha2(gI2);
   }
   tmp_3042 += tmp_3043;
   result += (-4) * tmp_3042;
   std::complex<double> tmp_3044;
   std::complex<double> tmp_3045;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3045 += B0(p,MChi(gI2),MVWm)*Conj(CpbarUCha2VWmChiPR(gO2,gI2))*
         CpbarUCha2VWmChiPL(gO1,gI2)*MChi(gI2);
   }
   tmp_3044 += tmp_3045;
   result += (-4) * tmp_3044;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha2_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3046;
   std::complex<double> tmp_3047;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3048;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3048 += B1(p,MCha1(gI1),MRh(gI2))*Conj(CpbarUCha2barCha1RhPR
            (gO2,gI1,gI2))*CpbarUCha2barCha1RhPR(gO1,gI1,gI2);
      }
      tmp_3047 += tmp_3048;
   }
   tmp_3046 += tmp_3047;
   result += (-0.5) * tmp_3046;
   std::complex<double> tmp_3049;
   std::complex<double> tmp_3050;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3051;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3051 += B1(p,MCha2(gI1),MAh(gI2))*Conj(CpbarUCha2Cha2AhPR(
            gO2,gI1,gI2))*CpbarUCha2Cha2AhPR(gO1,gI1,gI2);
      }
      tmp_3050 += tmp_3051;
   }
   tmp_3049 += tmp_3050;
   result += (-0.5) * tmp_3049;
   std::complex<double> tmp_3052;
   std::complex<double> tmp_3053;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3054;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3054 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUCha2barFuSdPR(gO2
            ,gI1,gI2))*CpbarUCha2barFuSdPR(gO1,gI1,gI2);
      }
      tmp_3053 += tmp_3054;
   }
   tmp_3052 += tmp_3053;
   result += (-1.5) * tmp_3052;
   std::complex<double> tmp_3055;
   std::complex<double> tmp_3056;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3057;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3057 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUCha2barFvSePR(gO2
            ,gI1,gI2))*CpbarUCha2barFvSePR(gO1,gI1,gI2);
      }
      tmp_3056 += tmp_3057;
   }
   tmp_3055 += tmp_3056;
   result += (-0.5) * tmp_3055;
   std::complex<double> tmp_3058;
   std::complex<double> tmp_3059;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3060;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3060 += B1(p,MChi(gI1),MRpm(gI2))*Conj(CpbarUCha2barChiRpmPR
            (gO2,gI1,gI2))*CpbarUCha2barChiRpmPR(gO1,gI1,gI2);
      }
      tmp_3059 += tmp_3060;
   }
   tmp_3058 += tmp_3059;
   result += (-0.5) * tmp_3058;
   std::complex<double> tmp_3061;
   std::complex<double> tmp_3062;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3063;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3063 += B1(p,MCha2(gI2),Mhh(gI1))*Conj(CpbarUCha2hhCha2PR(
            gO2,gI1,gI2))*CpbarUCha2hhCha2PR(gO1,gI1,gI2);
      }
      tmp_3062 += tmp_3063;
   }
   tmp_3061 += tmp_3062;
   result += (-0.5) * tmp_3061;
   std::complex<double> tmp_3064;
   std::complex<double> tmp_3065;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3066;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3066 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUCha2HpmChiPR(
            gO2,gI1,gI2))*CpbarUCha2HpmChiPR(gO1,gI1,gI2);
      }
      tmp_3065 += tmp_3066;
   }
   tmp_3064 += tmp_3065;
   result += (-0.5) * tmp_3064;
   std::complex<double> tmp_3067;
   std::complex<double> tmp_3068;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3069;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3069 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUCha2conjSuFdPR(
            gO2,gI1,gI2))*CpbarUCha2conjSuFdPR(gO1,gI1,gI2);
      }
      tmp_3068 += tmp_3069;
   }
   tmp_3067 += tmp_3068;
   result += (-1.5) * tmp_3067;
   std::complex<double> tmp_3070;
   std::complex<double> tmp_3071;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3071 += B1(p,MCha2(gI2),0)*Conj(CpbarUCha2VPCha2PL(gO2,gI2))*
         CpbarUCha2VPCha2PL(gO1,gI2);
   }
   tmp_3070 += tmp_3071;
   result += (-1) * tmp_3070;
   std::complex<double> tmp_3072;
   std::complex<double> tmp_3073;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3073 += B1(p,MCha2(gI2),MVZ)*Conj(CpbarUCha2VZCha2PL(gO2,gI2))*
         CpbarUCha2VZCha2PL(gO1,gI2);
   }
   tmp_3072 += tmp_3073;
   result += (-1) * tmp_3072;
   std::complex<double> tmp_3074;
   std::complex<double> tmp_3075;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3075 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha2VWmChiPL(gO2,gI2))*
         CpbarUCha2VWmChiPL(gO1,gI2);
   }
   tmp_3074 += tmp_3075;
   result += (-1) * tmp_3074;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Cha2_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3076;
   std::complex<double> tmp_3077;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3078;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3078 += B1(p,MCha1(gI1),MRh(gI2))*Conj(CpbarUCha2barCha1RhPL
            (gO2,gI1,gI2))*CpbarUCha2barCha1RhPL(gO1,gI1,gI2);
      }
      tmp_3077 += tmp_3078;
   }
   tmp_3076 += tmp_3077;
   result += (-0.5) * tmp_3076;
   std::complex<double> tmp_3079;
   std::complex<double> tmp_3080;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3081;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3081 += B1(p,MCha2(gI1),MAh(gI2))*Conj(CpbarUCha2Cha2AhPL(
            gO2,gI1,gI2))*CpbarUCha2Cha2AhPL(gO1,gI1,gI2);
      }
      tmp_3080 += tmp_3081;
   }
   tmp_3079 += tmp_3080;
   result += (-0.5) * tmp_3079;
   std::complex<double> tmp_3082;
   std::complex<double> tmp_3083;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3084;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3084 += B1(p,MFu(gI1),MSd(gI2))*Conj(CpbarUCha2barFuSdPL(gO2
            ,gI1,gI2))*CpbarUCha2barFuSdPL(gO1,gI1,gI2);
      }
      tmp_3083 += tmp_3084;
   }
   tmp_3082 += tmp_3083;
   result += (-1.5) * tmp_3082;
   std::complex<double> tmp_3085;
   std::complex<double> tmp_3086;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3087;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3087 += B1(p,MFv(gI1),MSe(gI2))*Conj(CpbarUCha2barFvSePL(gO2
            ,gI1,gI2))*CpbarUCha2barFvSePL(gO1,gI1,gI2);
      }
      tmp_3086 += tmp_3087;
   }
   tmp_3085 += tmp_3086;
   result += (-0.5) * tmp_3085;
   std::complex<double> tmp_3088;
   std::complex<double> tmp_3089;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3090;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3090 += B1(p,MChi(gI1),MRpm(gI2))*Conj(CpbarUCha2barChiRpmPL
            (gO2,gI1,gI2))*CpbarUCha2barChiRpmPL(gO1,gI1,gI2);
      }
      tmp_3089 += tmp_3090;
   }
   tmp_3088 += tmp_3089;
   result += (-0.5) * tmp_3088;
   std::complex<double> tmp_3091;
   std::complex<double> tmp_3092;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3093;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3093 += B1(p,MCha2(gI2),Mhh(gI1))*Conj(CpbarUCha2hhCha2PL(
            gO2,gI1,gI2))*CpbarUCha2hhCha2PL(gO1,gI1,gI2);
      }
      tmp_3092 += tmp_3093;
   }
   tmp_3091 += tmp_3092;
   result += (-0.5) * tmp_3091;
   std::complex<double> tmp_3094;
   std::complex<double> tmp_3095;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3096;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3096 += B1(p,MChi(gI2),MHpm(gI1))*Conj(CpbarUCha2HpmChiPL(
            gO2,gI1,gI2))*CpbarUCha2HpmChiPL(gO1,gI1,gI2);
      }
      tmp_3095 += tmp_3096;
   }
   tmp_3094 += tmp_3095;
   result += (-0.5) * tmp_3094;
   std::complex<double> tmp_3097;
   std::complex<double> tmp_3098;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3099;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3099 += B1(p,MFd(gI2),MSu(gI1))*Conj(CpbarUCha2conjSuFdPL(
            gO2,gI1,gI2))*CpbarUCha2conjSuFdPL(gO1,gI1,gI2);
      }
      tmp_3098 += tmp_3099;
   }
   tmp_3097 += tmp_3098;
   result += (-1.5) * tmp_3097;
   std::complex<double> tmp_3100;
   std::complex<double> tmp_3101;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3101 += B1(p,MCha2(gI2),0)*Conj(CpbarUCha2VPCha2PR(gO2,gI2))*
         CpbarUCha2VPCha2PR(gO1,gI2);
   }
   tmp_3100 += tmp_3101;
   result += (-1) * tmp_3100;
   std::complex<double> tmp_3102;
   std::complex<double> tmp_3103;
   for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
      tmp_3103 += B1(p,MCha2(gI2),MVZ)*Conj(CpbarUCha2VZCha2PR(gO2,gI2))*
         CpbarUCha2VZCha2PR(gO1,gI2);
   }
   tmp_3102 += tmp_3103;
   result += (-1) * tmp_3102;
   std::complex<double> tmp_3104;
   std::complex<double> tmp_3105;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3105 += B1(p,MChi(gI2),MVWm)*Conj(CpbarUCha2VWmChiPR(gO2,gI2))*
         CpbarUCha2VWmChiPR(gO1,gI2);
   }
   tmp_3104 += tmp_3105;
   result += (-1) * tmp_3104;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3106;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3107;
      std::complex<double> tmp_3108;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3108 += B0(p,MCha1(gI1),MSv(gI2))*Conj(CpbarUFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarUFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3107 += tmp_3108;
      tmp_3106 += (MCha1(gI1)) * tmp_3107;
   }
   result += tmp_3106;
   std::complex<double> tmp_3109;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3110;
      std::complex<double> tmp_3111;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3111 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3110 += tmp_3111;
      tmp_3109 += (MFe(gI1)) * tmp_3110;
   }
   result += tmp_3109;
   std::complex<double> tmp_3112;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3113;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3113 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3112 += tmp_3113;
   }
   result += tmp_3112;
   std::complex<double> tmp_3114;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3115;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3115 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3114 += tmp_3115;
   }
   result += tmp_3114;
   std::complex<double> tmp_3116;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3117;
      std::complex<double> tmp_3118;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3118 += B0(p,MChi(gI1),MSe(gI2))*Conj(CpbarUFebarChiSePL(gO2
            ,gI1,gI2))*CpbarUFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3117 += tmp_3118;
      tmp_3116 += (MChi(gI1)) * tmp_3117;
   }
   result += tmp_3116;
   std::complex<double> tmp_3119;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3120;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3120 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3119 += tmp_3120;
   }
   result += tmp_3119;
   std::complex<double> tmp_3121;
   std::complex<double> tmp_3122;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3122 += B0(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3121 += tmp_3122;
   result += (-4) * tmp_3121;
   std::complex<double> tmp_3123;
   std::complex<double> tmp_3124;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3124 += B0(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3123 += tmp_3124;
   result += (-4) * tmp_3123;
   std::complex<double> tmp_3125;
   std::complex<double> tmp_3126;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3126 += B0(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3125 += tmp_3126;
   result += (-4) * tmp_3125;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3127;
   std::complex<double> tmp_3128;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3129;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3129 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarUFebarCha1SvPR(
            gO2,gI1,gI2))*CpbarUFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3128 += tmp_3129;
   }
   tmp_3127 += tmp_3128;
   result += (-0.5) * tmp_3127;
   std::complex<double> tmp_3130;
   std::complex<double> tmp_3131;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3132;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3132 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPR(gO2,gI1,
            gI2))*CpbarUFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3131 += tmp_3132;
   }
   tmp_3130 += tmp_3131;
   result += (-0.5) * tmp_3130;
   std::complex<double> tmp_3133;
   std::complex<double> tmp_3134;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3135;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3135 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePR(gO2,gI1,
            gI2))*CpbarUFehhFePR(gO1,gI1,gI2);
      }
      tmp_3134 += tmp_3135;
   }
   tmp_3133 += tmp_3134;
   result += (-0.5) * tmp_3133;
   std::complex<double> tmp_3136;
   std::complex<double> tmp_3137;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3138;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3138 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPR(gO2,
            gI1,gI2))*CpbarUFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3137 += tmp_3138;
   }
   tmp_3136 += tmp_3137;
   result += (-0.5) * tmp_3136;
   std::complex<double> tmp_3139;
   std::complex<double> tmp_3140;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3141;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3141 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarUFebarChiSePR(gO2
            ,gI1,gI2))*CpbarUFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3140 += tmp_3141;
   }
   tmp_3139 += tmp_3140;
   result += (-0.5) * tmp_3139;
   std::complex<double> tmp_3142;
   std::complex<double> tmp_3143;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3144;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3144 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPR(gO2,
            gI1,gI2))*CpbarUFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3143 += tmp_3144;
   }
   tmp_3142 += tmp_3143;
   result += (-0.5) * tmp_3142;
   std::complex<double> tmp_3145;
   std::complex<double> tmp_3146;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3146 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePL(gO2,gI2))*
         CpbarUFeVPFePL(gO1,gI2);
   }
   tmp_3145 += tmp_3146;
   result += (-1) * tmp_3145;
   std::complex<double> tmp_3147;
   std::complex<double> tmp_3148;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3148 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPL(gO2,gI2))*
         CpbarUFeVWmFvPL(gO1,gI2);
   }
   tmp_3147 += tmp_3148;
   result += (-1) * tmp_3147;
   std::complex<double> tmp_3149;
   std::complex<double> tmp_3150;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3150 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePL(gO2,gI2))*
         CpbarUFeVZFePL(gO1,gI2);
   }
   tmp_3149 += tmp_3150;
   result += (-1) * tmp_3149;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3151;
   std::complex<double> tmp_3152;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3153;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3153 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarUFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarUFebarCha1SvPL(gO1,gI1,gI2);
      }
      tmp_3152 += tmp_3153;
   }
   tmp_3151 += tmp_3152;
   result += (-0.5) * tmp_3151;
   std::complex<double> tmp_3154;
   std::complex<double> tmp_3155;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3156;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3156 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarUFeFeAhPL(gO2,gI1,
            gI2))*CpbarUFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3155 += tmp_3156;
   }
   tmp_3154 += tmp_3155;
   result += (-0.5) * tmp_3154;
   std::complex<double> tmp_3157;
   std::complex<double> tmp_3158;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3159;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3159 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarUFehhFePL(gO2,gI1,
            gI2))*CpbarUFehhFePL(gO1,gI1,gI2);
      }
      tmp_3158 += tmp_3159;
   }
   tmp_3157 += tmp_3158;
   result += (-0.5) * tmp_3157;
   std::complex<double> tmp_3160;
   std::complex<double> tmp_3161;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3162;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3162 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarUFeHpmFvPL(gO2,
            gI1,gI2))*CpbarUFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3161 += tmp_3162;
   }
   tmp_3160 += tmp_3161;
   result += (-0.5) * tmp_3160;
   std::complex<double> tmp_3163;
   std::complex<double> tmp_3164;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3165;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3165 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarUFebarChiSePL(gO2
            ,gI1,gI2))*CpbarUFebarChiSePL(gO1,gI1,gI2);
      }
      tmp_3164 += tmp_3165;
   }
   tmp_3163 += tmp_3164;
   result += (-0.5) * tmp_3163;
   std::complex<double> tmp_3166;
   std::complex<double> tmp_3167;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3168;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3168 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarUFeSeChiPL(gO2,
            gI1,gI2))*CpbarUFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3167 += tmp_3168;
   }
   tmp_3166 += tmp_3167;
   result += (-0.5) * tmp_3166;
   std::complex<double> tmp_3169;
   std::complex<double> tmp_3170;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3170 += B1(p,MFe(gI2),0)*Conj(CpbarUFeVPFePR(gO2,gI2))*
         CpbarUFeVPFePR(gO1,gI2);
   }
   tmp_3169 += tmp_3170;
   result += (-1) * tmp_3169;
   std::complex<double> tmp_3171;
   std::complex<double> tmp_3172;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3172 += B1(p,MFv(gI2),MVWm)*Conj(CpbarUFeVWmFvPR(gO2,gI2))*
         CpbarUFeVWmFvPR(gO1,gI2);
   }
   tmp_3171 += tmp_3172;
   result += (-1) * tmp_3171;
   std::complex<double> tmp_3173;
   std::complex<double> tmp_3174;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3174 += B1(p,MFe(gI2),MVZ)*Conj(CpbarUFeVZFePR(gO2,gI2))*
         CpbarUFeVZFePR(gO1,gI2);
   }
   tmp_3173 += tmp_3174;
   result += (-1) * tmp_3173;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3175;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3176;
      std::complex<double> tmp_3177;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3177 += B0(p,MCha1(gI1),MSu(gI2))*Conj(CpbarUFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarUFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3176 += tmp_3177;
      tmp_3175 += (MCha1(gI1)) * tmp_3176;
   }
   result += tmp_3175;
   std::complex<double> tmp_3178;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3179;
      std::complex<double> tmp_3180;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3180 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3179 += tmp_3180;
      tmp_3178 += (MFd(gI1)) * tmp_3179;
   }
   result += tmp_3178;
   std::complex<double> tmp_3181;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3182;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3182 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3181 += tmp_3182;
   }
   result += tmp_3181;
   std::complex<double> tmp_3183;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3184;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3184 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3183 += tmp_3184;
   }
   result += tmp_3183;
   std::complex<double> tmp_3185;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3186;
      std::complex<double> tmp_3187;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3187 += B0(p,MChi(gI1),MSd(gI2))*Conj(CpbarUFdbarChiSdPL(gO2
            ,gI1,gI2))*CpbarUFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3186 += tmp_3187;
      tmp_3185 += (MChi(gI1)) * tmp_3186;
   }
   result += tmp_3185;
   std::complex<double> tmp_3188;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3189;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3189 += B0(p,MCha2(gI2),MSu(gI1))*Conj(CpbarUFdSuCha2PL(gO2,
            gI1,gI2))*CpbarUFdSuCha2PR(gO1,gI1,gI2)*MCha2(gI2);
      }
      tmp_3188 += tmp_3189;
   }
   result += tmp_3188;
   std::complex<double> tmp_3190;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3191;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3191 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3190 += tmp_3191;
   }
   result += tmp_3190;
   std::complex<double> tmp_3192;
   std::complex<double> tmp_3193;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3193 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3192 += tmp_3193;
   result += (-5.333333333333333) * tmp_3192;
   std::complex<double> tmp_3194;
   std::complex<double> tmp_3195;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3195 += B0(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3194 += tmp_3195;
   result += (-4) * tmp_3194;
   std::complex<double> tmp_3196;
   std::complex<double> tmp_3197;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3197 += B0(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3196 += tmp_3197;
   result += (-4) * tmp_3196;
   std::complex<double> tmp_3198;
   std::complex<double> tmp_3199;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3199 += B0(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3198 += tmp_3199;
   result += (-4) * tmp_3198;
   std::complex<double> tmp_3200;
   std::complex<double> tmp_3201;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3201 += B0(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_3200 += tmp_3201;
   result += (1.3333333333333333*MGlu) * tmp_3200;
   std::complex<double> tmp_3202;
   std::complex<double> tmp_3203;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3203 += B0(p,MGlu,MSd(gI2))*Conj(CpbarUFdbarGluSdPL(gO2,1,gI2))*
         CpbarUFdbarGluSdPR(gO1,1,gI2);
   }
   tmp_3202 += tmp_3203;
   result += (1.3333333333333333*MGlu) * tmp_3202;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3204;
   std::complex<double> tmp_3205;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3206;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3206 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarUFdbarCha1SuPR(
            gO2,gI1,gI2))*CpbarUFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3205 += tmp_3206;
   }
   tmp_3204 += tmp_3205;
   result += (-0.5) * tmp_3204;
   std::complex<double> tmp_3207;
   std::complex<double> tmp_3208;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3209;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3209 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPR(gO2,gI1,
            gI2))*CpbarUFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3208 += tmp_3209;
   }
   tmp_3207 += tmp_3208;
   result += (-0.5) * tmp_3207;
   std::complex<double> tmp_3210;
   std::complex<double> tmp_3211;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3212;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3212 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPR(gO2,gI1,
            gI2))*CpbarUFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3211 += tmp_3212;
   }
   tmp_3210 += tmp_3211;
   result += (-0.5) * tmp_3210;
   std::complex<double> tmp_3213;
   std::complex<double> tmp_3214;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3215;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3215 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPR(gO2,
            gI1,gI2))*CpbarUFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3214 += tmp_3215;
   }
   tmp_3213 += tmp_3214;
   result += (-0.5) * tmp_3213;
   std::complex<double> tmp_3216;
   std::complex<double> tmp_3217;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3218;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3218 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarUFdbarChiSdPR(gO2
            ,gI1,gI2))*CpbarUFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3217 += tmp_3218;
   }
   tmp_3216 += tmp_3217;
   result += (-0.5) * tmp_3216;
   std::complex<double> tmp_3219;
   std::complex<double> tmp_3220;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3220 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPR(gO2,gI1,1))*
         CpbarUFdSdGluPR(gO1,gI1,1);
   }
   tmp_3219 += tmp_3220;
   result += (-0.6666666666666666) * tmp_3219;
   std::complex<double> tmp_3221;
   std::complex<double> tmp_3222;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3223;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3223 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarUFdSuCha2PR(gO2,
            gI1,gI2))*CpbarUFdSuCha2PR(gO1,gI1,gI2);
      }
      tmp_3222 += tmp_3223;
   }
   tmp_3221 += tmp_3222;
   result += (-0.5) * tmp_3221;
   std::complex<double> tmp_3224;
   std::complex<double> tmp_3225;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3226;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3226 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPR(gO2,
            gI1,gI2))*CpbarUFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3225 += tmp_3226;
   }
   tmp_3224 += tmp_3225;
   result += (-0.5) * tmp_3224;
   std::complex<double> tmp_3227;
   std::complex<double> tmp_3228;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3228 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPL(gO2,gI2))*
         CpbarUFdVGFdPL(gO1,gI2);
   }
   tmp_3227 += tmp_3228;
   result += (-1.3333333333333333) * tmp_3227;
   std::complex<double> tmp_3229;
   std::complex<double> tmp_3230;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3230 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPL(gO2,gI2))*
         CpbarUFdVPFdPL(gO1,gI2);
   }
   tmp_3229 += tmp_3230;
   result += (-1) * tmp_3229;
   std::complex<double> tmp_3231;
   std::complex<double> tmp_3232;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3232 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPL(gO2,gI2))*
         CpbarUFdVWmFuPL(gO1,gI2);
   }
   tmp_3231 += tmp_3232;
   result += (-1) * tmp_3231;
   std::complex<double> tmp_3233;
   std::complex<double> tmp_3234;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3234 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPL(gO2,gI2))*
         CpbarUFdVZFdPL(gO1,gI2);
   }
   tmp_3233 += tmp_3234;
   result += (-1) * tmp_3233;
   std::complex<double> tmp_3235;
   std::complex<double> tmp_3236;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3236 += B1(p,MGlu,MSd(gI2))*Conj(CpbarUFdbarGluSdPR(gO2,1,gI2))*
         CpbarUFdbarGluSdPR(gO1,1,gI2);
   }
   tmp_3235 += tmp_3236;
   result += (-0.6666666666666666) * tmp_3235;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3237;
   std::complex<double> tmp_3238;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3239;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3239 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarUFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarUFdbarCha1SuPL(gO1,gI1,gI2);
      }
      tmp_3238 += tmp_3239;
   }
   tmp_3237 += tmp_3238;
   result += (-0.5) * tmp_3237;
   std::complex<double> tmp_3240;
   std::complex<double> tmp_3241;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3242;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3242 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarUFdFdAhPL(gO2,gI1,
            gI2))*CpbarUFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3241 += tmp_3242;
   }
   tmp_3240 += tmp_3241;
   result += (-0.5) * tmp_3240;
   std::complex<double> tmp_3243;
   std::complex<double> tmp_3244;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3245;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3245 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarUFdhhFdPL(gO2,gI1,
            gI2))*CpbarUFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3244 += tmp_3245;
   }
   tmp_3243 += tmp_3244;
   result += (-0.5) * tmp_3243;
   std::complex<double> tmp_3246;
   std::complex<double> tmp_3247;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3248;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3248 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarUFdHpmFuPL(gO2,
            gI1,gI2))*CpbarUFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3247 += tmp_3248;
   }
   tmp_3246 += tmp_3247;
   result += (-0.5) * tmp_3246;
   std::complex<double> tmp_3249;
   std::complex<double> tmp_3250;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3251;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3251 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarUFdbarChiSdPL(gO2
            ,gI1,gI2))*CpbarUFdbarChiSdPL(gO1,gI1,gI2);
      }
      tmp_3250 += tmp_3251;
   }
   tmp_3249 += tmp_3250;
   result += (-0.5) * tmp_3249;
   std::complex<double> tmp_3252;
   std::complex<double> tmp_3253;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3253 += B1(p,MGlu,MSd(gI1))*Conj(CpbarUFdSdGluPL(gO2,gI1,1))*
         CpbarUFdSdGluPL(gO1,gI1,1);
   }
   tmp_3252 += tmp_3253;
   result += (-0.6666666666666666) * tmp_3252;
   std::complex<double> tmp_3254;
   std::complex<double> tmp_3255;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3256;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3256 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarUFdSuCha2PL(gO2,
            gI1,gI2))*CpbarUFdSuCha2PL(gO1,gI1,gI2);
      }
      tmp_3255 += tmp_3256;
   }
   tmp_3254 += tmp_3255;
   result += (-0.5) * tmp_3254;
   std::complex<double> tmp_3257;
   std::complex<double> tmp_3258;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3259;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3259 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarUFdSdChiPL(gO2,
            gI1,gI2))*CpbarUFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3258 += tmp_3259;
   }
   tmp_3257 += tmp_3258;
   result += (-0.5) * tmp_3257;
   std::complex<double> tmp_3260;
   std::complex<double> tmp_3261;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3261 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVGFdPR(gO2,gI2))*
         CpbarUFdVGFdPR(gO1,gI2);
   }
   tmp_3260 += tmp_3261;
   result += (-1.3333333333333333) * tmp_3260;
   std::complex<double> tmp_3262;
   std::complex<double> tmp_3263;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3263 += B1(p,MFd(gI2),0)*Conj(CpbarUFdVPFdPR(gO2,gI2))*
         CpbarUFdVPFdPR(gO1,gI2);
   }
   tmp_3262 += tmp_3263;
   result += (-1) * tmp_3262;
   std::complex<double> tmp_3264;
   std::complex<double> tmp_3265;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3265 += B1(p,MFu(gI2),MVWm)*Conj(CpbarUFdVWmFuPR(gO2,gI2))*
         CpbarUFdVWmFuPR(gO1,gI2);
   }
   tmp_3264 += tmp_3265;
   result += (-1) * tmp_3264;
   std::complex<double> tmp_3266;
   std::complex<double> tmp_3267;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3267 += B1(p,MFd(gI2),MVZ)*Conj(CpbarUFdVZFdPR(gO2,gI2))*
         CpbarUFdVZFdPR(gO1,gI2);
   }
   tmp_3266 += tmp_3267;
   result += (-1) * tmp_3266;
   std::complex<double> tmp_3268;
   std::complex<double> tmp_3269;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3269 += B1(p,MGlu,MSd(gI2))*Conj(CpbarUFdbarGluSdPL(gO2,1,gI2))*
         CpbarUFdbarGluSdPL(gO1,1,gI2);
   }
   tmp_3268 += tmp_3269;
   result += (-0.6666666666666666) * tmp_3268;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3270;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3271;
      std::complex<double> tmp_3272;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3272 += B0(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3271 += tmp_3272;
      tmp_3270 += (MCha2(gI1)) * tmp_3271;
   }
   result += tmp_3270;
   std::complex<double> tmp_3273;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3274;
      std::complex<double> tmp_3275;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3275 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3274 += tmp_3275;
      tmp_3273 += (MFu(gI1)) * tmp_3274;
   }
   result += tmp_3273;
   std::complex<double> tmp_3276;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3277;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3277 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3276 += tmp_3277;
   }
   result += tmp_3276;
   std::complex<double> tmp_3278;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3279;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3279 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3278 += tmp_3279;
   }
   result += tmp_3278;
   std::complex<double> tmp_3280;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3281;
      std::complex<double> tmp_3282;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3282 += B0(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3281 += tmp_3282;
      tmp_3280 += (MChi(gI1)) * tmp_3281;
   }
   result += tmp_3280;
   std::complex<double> tmp_3283;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3284;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3284 += B0(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3283 += tmp_3284;
   }
   result += tmp_3283;
   std::complex<double> tmp_3285;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3286;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3286 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3285 += tmp_3286;
   }
   result += tmp_3285;
   std::complex<double> tmp_3287;
   std::complex<double> tmp_3288;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3288 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3287 += tmp_3288;
   result += (-4) * tmp_3287;
   std::complex<double> tmp_3289;
   std::complex<double> tmp_3290;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3290 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3289 += tmp_3290;
   result += (-5.333333333333333) * tmp_3289;
   std::complex<double> tmp_3291;
   std::complex<double> tmp_3292;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3292 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3291 += tmp_3292;
   result += (-4) * tmp_3291;
   std::complex<double> tmp_3293;
   std::complex<double> tmp_3294;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3294 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3293 += tmp_3294;
   result += (-4) * tmp_3293;
   std::complex<double> tmp_3295;
   std::complex<double> tmp_3296;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3296 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3295 += tmp_3296;
   result += (1.3333333333333333*MGlu) * tmp_3295;
   std::complex<double> tmp_3297;
   std::complex<double> tmp_3298;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3298 += B0(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,1,gI2))*
         CpbarUFubarGluSuPR(gO1,1,gI2);
   }
   tmp_3297 += tmp_3298;
   result += (1.3333333333333333*MGlu) * tmp_3297;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3299;
   std::complex<double> tmp_3300;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3301;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3301 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPR(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3300 += tmp_3301;
   }
   tmp_3299 += tmp_3300;
   result += (-0.5) * tmp_3299;
   std::complex<double> tmp_3302;
   std::complex<double> tmp_3303;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3304;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3304 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3303 += tmp_3304;
   }
   tmp_3302 += tmp_3303;
   result += (-0.5) * tmp_3302;
   std::complex<double> tmp_3305;
   std::complex<double> tmp_3306;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3307;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3307 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3306 += tmp_3307;
   }
   tmp_3305 += tmp_3306;
   result += (-0.5) * tmp_3305;
   std::complex<double> tmp_3308;
   std::complex<double> tmp_3309;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3310;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3310 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3309 += tmp_3310;
   }
   tmp_3308 += tmp_3309;
   result += (-0.5) * tmp_3308;
   std::complex<double> tmp_3311;
   std::complex<double> tmp_3312;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3313;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3313 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPR(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3312 += tmp_3313;
   }
   tmp_3311 += tmp_3312;
   result += (-0.5) * tmp_3311;
   std::complex<double> tmp_3314;
   std::complex<double> tmp_3315;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3315 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3314 += tmp_3315;
   result += (-0.6666666666666666) * tmp_3314;
   std::complex<double> tmp_3316;
   std::complex<double> tmp_3317;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3318;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3318 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PR(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2);
      }
      tmp_3317 += tmp_3318;
   }
   tmp_3316 += tmp_3317;
   result += (-0.5) * tmp_3316;
   std::complex<double> tmp_3319;
   std::complex<double> tmp_3320;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3321;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3321 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3320 += tmp_3321;
   }
   tmp_3319 += tmp_3320;
   result += (-0.5) * tmp_3319;
   std::complex<double> tmp_3322;
   std::complex<double> tmp_3323;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3323 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3322 += tmp_3323;
   result += (-1) * tmp_3322;
   std::complex<double> tmp_3324;
   std::complex<double> tmp_3325;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3325 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPL(gO2,gI2))*
         CpbarUFuVGFuPL(gO1,gI2);
   }
   tmp_3324 += tmp_3325;
   result += (-1.3333333333333333) * tmp_3324;
   std::complex<double> tmp_3326;
   std::complex<double> tmp_3327;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3327 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3326 += tmp_3327;
   result += (-1) * tmp_3326;
   std::complex<double> tmp_3328;
   std::complex<double> tmp_3329;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3329 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3328 += tmp_3329;
   result += (-1) * tmp_3328;
   std::complex<double> tmp_3330;
   std::complex<double> tmp_3331;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3331 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPR(gO2,1,gI2))*
         CpbarUFubarGluSuPR(gO1,1,gI2);
   }
   tmp_3330 += tmp_3331;
   result += (-0.6666666666666666) * tmp_3330;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3332;
   std::complex<double> tmp_3333;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3334;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3334 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPL(gO1,gI1,gI2);
      }
      tmp_3333 += tmp_3334;
   }
   tmp_3332 += tmp_3333;
   result += (-0.5) * tmp_3332;
   std::complex<double> tmp_3335;
   std::complex<double> tmp_3336;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3337;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3337 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3336 += tmp_3337;
   }
   tmp_3335 += tmp_3336;
   result += (-0.5) * tmp_3335;
   std::complex<double> tmp_3338;
   std::complex<double> tmp_3339;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3340;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3340 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3339 += tmp_3340;
   }
   tmp_3338 += tmp_3339;
   result += (-0.5) * tmp_3338;
   std::complex<double> tmp_3341;
   std::complex<double> tmp_3342;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3343;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3343 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3342 += tmp_3343;
   }
   tmp_3341 += tmp_3342;
   result += (-0.5) * tmp_3341;
   std::complex<double> tmp_3344;
   std::complex<double> tmp_3345;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3346;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3346 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPL(gO1,gI1,gI2);
      }
      tmp_3345 += tmp_3346;
   }
   tmp_3344 += tmp_3345;
   result += (-0.5) * tmp_3344;
   std::complex<double> tmp_3347;
   std::complex<double> tmp_3348;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3348 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_3347 += tmp_3348;
   result += (-0.6666666666666666) * tmp_3347;
   std::complex<double> tmp_3349;
   std::complex<double> tmp_3350;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3351;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3351 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PL(gO1,gI1,gI2);
      }
      tmp_3350 += tmp_3351;
   }
   tmp_3349 += tmp_3350;
   result += (-0.5) * tmp_3349;
   std::complex<double> tmp_3352;
   std::complex<double> tmp_3353;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3354;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3354 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3353 += tmp_3354;
   }
   tmp_3352 += tmp_3353;
   result += (-0.5) * tmp_3352;
   std::complex<double> tmp_3355;
   std::complex<double> tmp_3356;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3356 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3355 += tmp_3356;
   result += (-1) * tmp_3355;
   std::complex<double> tmp_3357;
   std::complex<double> tmp_3358;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3358 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVGFuPR(gO2,gI2))*
         CpbarUFuVGFuPR(gO1,gI2);
   }
   tmp_3357 += tmp_3358;
   result += (-1.3333333333333333) * tmp_3357;
   std::complex<double> tmp_3359;
   std::complex<double> tmp_3360;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3360 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3359 += tmp_3360;
   result += (-1) * tmp_3359;
   std::complex<double> tmp_3361;
   std::complex<double> tmp_3362;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3362 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3361 += tmp_3362;
   result += (-1) * tmp_3361;
   std::complex<double> tmp_3363;
   std::complex<double> tmp_3364;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3364 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,1,gI2))*
         CpbarUFubarGluSuPL(gO1,1,gI2);
   }
   tmp_3363 += tmp_3364;
   result += (-0.6666666666666666) * tmp_3363;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_1(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3365;
   std::complex<double> tmp_3366;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3367;
      std::complex<double> tmp_3368;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3368 += B0(p,MFd(gI1),MSd(gI2))*Conj(CpbarGlubarFdSdPL(gI1,
            gI2))*CpbarGlubarFdSdPR(gI1,gI2);
      }
      tmp_3367 += tmp_3368;
      tmp_3366 += (MFd(gI1)) * tmp_3367;
   }
   tmp_3365 += tmp_3366;
   result += (0.5) * tmp_3365;
   std::complex<double> tmp_3369;
   std::complex<double> tmp_3370;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3371;
      std::complex<double> tmp_3372;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3372 += B0(p,MFu(gI1),MSu(gI2))*Conj(CpbarGlubarFuSuPL(gI1,
            gI2))*CpbarGlubarFuSuPR(gI1,gI2);
      }
      tmp_3371 += tmp_3372;
      tmp_3370 += (MFu(gI1)) * tmp_3371;
   }
   tmp_3369 += tmp_3370;
   result += (0.5) * tmp_3369;
   std::complex<double> tmp_3373;
   std::complex<double> tmp_3374;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3375;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3375 += B0(p,MFd(gI2),MSd(gI1))*Conj(CpbarGluconjSdFdPL(gI1,
            gI2))*CpbarGluconjSdFdPR(gI1,gI2)*MFd(gI2);
      }
      tmp_3374 += tmp_3375;
   }
   tmp_3373 += tmp_3374;
   result += (0.5) * tmp_3373;
   std::complex<double> tmp_3376;
   std::complex<double> tmp_3377;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3378;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3378 += B0(p,MFu(gI2),MSu(gI1))*Conj(CpbarGluconjSuFuPL(gI1,
            gI2))*CpbarGluconjSuFuPR(gI1,gI2)*MFu(gI2);
      }
      tmp_3377 += tmp_3378;
   }
   tmp_3376 += tmp_3377;
   result += (0.5) * tmp_3376;
   result += -3*MGlu*B0(p,MGlu,MSOc)*Conj(CpbarGluconjSOcGluPL())*
      CpbarGluconjSOcGluPR();
   result += -12*MGlu*B0(p,MGlu,0)*Conj(CpbarGluVGGluPR())*CpbarGluVGGluPL();

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PR(double p ) const
{
   std::complex<double> result;

   result += 1.5*AbsSqr(CpbarGluconjSOcGluPR())*B1(p,MGlu,MSOc);
   result += -3*AbsSqr(CpbarGluVGGluPL())*B1(p,MGlu,0);
   std::complex<double> tmp_3379;
   std::complex<double> tmp_3380;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3381;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3381 += AbsSqr(CpbarGlubarFdSdPR(gI1,gI2))*B1(p,MFd(gI1),MSd
            (gI2));
      }
      tmp_3380 += tmp_3381;
   }
   tmp_3379 += tmp_3380;
   result += (-0.25) * tmp_3379;
   std::complex<double> tmp_3382;
   std::complex<double> tmp_3383;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3384;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3384 += AbsSqr(CpbarGlubarFuSuPR(gI1,gI2))*B1(p,MFu(gI1),MSu
            (gI2));
      }
      tmp_3383 += tmp_3384;
   }
   tmp_3382 += tmp_3383;
   result += (-0.25) * tmp_3382;
   std::complex<double> tmp_3385;
   std::complex<double> tmp_3386;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3387;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3387 += AbsSqr(CpbarGluconjSdFdPR(gI1,gI2))*B1(p,MFd(gI2),
            MSd(gI1));
      }
      tmp_3386 += tmp_3387;
   }
   tmp_3385 += tmp_3386;
   result += (-0.25) * tmp_3385;
   std::complex<double> tmp_3388;
   std::complex<double> tmp_3389;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3390;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3390 += AbsSqr(CpbarGluconjSuFuPR(gI1,gI2))*B1(p,MFu(gI2),
            MSu(gI1));
      }
      tmp_3389 += tmp_3390;
   }
   tmp_3388 += tmp_3389;
   result += (-0.25) * tmp_3388;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Glu_PL(double p ) const
{
   std::complex<double> result;

   result += 1.5*AbsSqr(CpbarGluconjSOcGluPL())*B1(p,MGlu,MSOc);
   result += -3*AbsSqr(CpbarGluVGGluPR())*B1(p,MGlu,0);
   std::complex<double> tmp_3391;
   std::complex<double> tmp_3392;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3393;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3393 += AbsSqr(CpbarGlubarFdSdPL(gI1,gI2))*B1(p,MFd(gI1),MSd
            (gI2));
      }
      tmp_3392 += tmp_3393;
   }
   tmp_3391 += tmp_3392;
   result += (-0.25) * tmp_3391;
   std::complex<double> tmp_3394;
   std::complex<double> tmp_3395;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3396;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3396 += AbsSqr(CpbarGlubarFuSuPL(gI1,gI2))*B1(p,MFu(gI1),MSu
            (gI2));
      }
      tmp_3395 += tmp_3396;
   }
   tmp_3394 += tmp_3395;
   result += (-0.25) * tmp_3394;
   std::complex<double> tmp_3397;
   std::complex<double> tmp_3398;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3399;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3399 += AbsSqr(CpbarGluconjSdFdPL(gI1,gI2))*B1(p,MFd(gI2),
            MSd(gI1));
      }
      tmp_3398 += tmp_3399;
   }
   tmp_3397 += tmp_3398;
   result += (-0.25) * tmp_3397;
   std::complex<double> tmp_3400;
   std::complex<double> tmp_3401;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3402;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3402 += AbsSqr(CpbarGluconjSuFuPL(gI1,gI2))*B1(p,MFu(gI2),
            MSu(gI1));
      }
      tmp_3401 += tmp_3402;
   }
   tmp_3400 += tmp_3401;
   result += (-0.25) * tmp_3400;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VZ_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3403;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3403 += A0(MRh(gI1))*CpVZVZconjRhRh(gI1,gI1);
   }
   result += tmp_3403;
   std::complex<double> tmp_3404;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3404 += A0(MRpm(gI1))*CpVZVZconjRpmRpm(gI1,gI1);
   }
   result += tmp_3404;
   std::complex<double> tmp_3405;
   std::complex<double> tmp_3406;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3407;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3407 += AbsSqr(CpVZconjRhRh(gI1,gI2))*B00(p,MRh(gI1),MRh(gI2
            ));
      }
      tmp_3406 += tmp_3407;
   }
   tmp_3405 += tmp_3406;
   result += (-4) * tmp_3405;
   std::complex<double> tmp_3408;
   std::complex<double> tmp_3409;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3410;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3410 += AbsSqr(CpVZconjRpmRpm(gI1,gI2))*B00(p,MRpm(gI1),MRpm
            (gI2));
      }
      tmp_3409 += tmp_3410;
   }
   tmp_3408 += tmp_3409;
   result += (-4) * tmp_3408;
   std::complex<double> tmp_3411;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3412;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3412 += (AbsSqr(CpVZbarCha1Cha1PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha1Cha1PR(gI1,gI2)))*H0(p,MCha1(gI1),MCha1(gI2));
         tmp_3412 += 4*B0(p,MCha1(gI1),MCha1(gI2))*MCha1(gI1)*MCha1(gI2)*
            Re(Conj(CpVZbarCha1Cha1PL(gI1,gI2))*CpVZbarCha1Cha1PR(gI1,gI2));
      }
      tmp_3411 += tmp_3412;
   }
   result += tmp_3411;
   std::complex<double> tmp_3413;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3414;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3414 += (AbsSqr(CpVZbarCha2Cha2PL(gI1,gI2)) + AbsSqr(
            CpVZbarCha2Cha2PR(gI1,gI2)))*H0(p,MCha2(gI1),MCha2(gI2));
         tmp_3414 += 4*B0(p,MCha2(gI1),MCha2(gI2))*MCha2(gI1)*MCha2(gI2)*
            Re(Conj(CpVZbarCha2Cha2PL(gI1,gI2))*CpVZbarCha2Cha2PR(gI1,gI2));
      }
      tmp_3413 += tmp_3414;
   }
   result += tmp_3413;
   std::complex<double> tmp_3415;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3415 += A0(MSv(gI1))*CpVZVZconjSvSv(gI1,gI1);
   }
   result += tmp_3415;
   std::complex<double> tmp_3416;
   std::complex<double> tmp_3417;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3418;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3418 += AbsSqr(CpVZconjSvSv(gI1,gI2))*B00(p,MSv(gI1),MSv(gI2
            ));
      }
      tmp_3417 += tmp_3418;
   }
   tmp_3416 += tmp_3417;
   result += (-4) * tmp_3416;
   std::complex<double> tmp_3419;
   std::complex<double> tmp_3420;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3420 += A0(Mhh(gI1))*CpVZVZhhhh(gI1,gI1);
   }
   tmp_3419 += tmp_3420;
   result += (0.5) * tmp_3419;
   std::complex<double> tmp_3421;
   std::complex<double> tmp_3422;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3423;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3423 += AbsSqr(CpVZhhAh(gI1,1 + gI2))*B00(p,MAh(1 + gI2),Mhh
            (gI1));
      }
      tmp_3422 += tmp_3423;
   }
   tmp_3421 += tmp_3422;
   result += (-4) * tmp_3421;
   std::complex<double> tmp_3424;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3425;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3425 += (AbsSqr(CpVZbarChiChiPL(gI1,gI2)) + AbsSqr(
            CpVZbarChiChiPR(gI1,gI2)))*H0(p,MChi(gI1),MChi(gI2));
         tmp_3425 += 4*B0(p,MChi(gI1),MChi(gI2))*MChi(gI1)*MChi(gI2)*Re(
            Conj(CpVZbarChiChiPL(gI1,gI2))*CpVZbarChiChiPR(gI1,gI2));
      }
      tmp_3424 += tmp_3425;
   }
   result += tmp_3424;
   std::complex<double> tmp_3426;
   std::complex<double> tmp_3427;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3427 += A0(MSd(gI1))*CpVZVZconjSdSd(gI1,gI1);
   }
   tmp_3426 += tmp_3427;
   result += (3) * tmp_3426;
   std::complex<double> tmp_3428;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3428 += A0(MSe(gI1))*CpVZVZconjSeSe(gI1,gI1);
   }
   result += tmp_3428;
   std::complex<double> tmp_3429;
   std::complex<double> tmp_3430;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3430 += A0(MSu(gI1))*CpVZVZconjSuSu(gI1,gI1);
   }
   tmp_3429 += tmp_3430;
   result += (3) * tmp_3429;
   std::complex<double> tmp_3431;
   std::complex<double> tmp_3432;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3433;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3433 += AbsSqr(CpVZconjSdSd(gI1,gI2))*B00(p,MSd(gI1),MSd(gI2
            ));
      }
      tmp_3432 += tmp_3433;
   }
   tmp_3431 += tmp_3432;
   result += (-12) * tmp_3431;
   std::complex<double> tmp_3434;
   std::complex<double> tmp_3435;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3436;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3436 += AbsSqr(CpVZconjSeSe(gI1,gI2))*B00(p,MSe(gI1),MSe(gI2
            ));
      }
      tmp_3435 += tmp_3436;
   }
   tmp_3434 += tmp_3435;
   result += (-4) * tmp_3434;
   std::complex<double> tmp_3437;
   std::complex<double> tmp_3438;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3439;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3439 += AbsSqr(CpVZconjSuSu(gI1,gI2))*B00(p,MSu(gI1),MSu(gI2
            ));
      }
      tmp_3438 += tmp_3439;
   }
   tmp_3437 += tmp_3438;
   result += (-12) * tmp_3437;
   std::complex<double> tmp_3440;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3440 += AbsSqr(CpVZVZhh(gI2))*B0(p,MVZ,Mhh(gI2));
   }
   result += tmp_3440;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_VWm_heavy(double p ) const
{
   std::complex<double> result;

   std::complex<double> tmp_3441;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3441 += A0(MRh(gI1))*CpVWmconjVWmconjRhRh(gI1,gI1);
   }
   result += tmp_3441;
   std::complex<double> tmp_3442;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3442 += A0(MRpm(gI1))*CpVWmconjVWmconjRpmRpm(gI1,gI1);
   }
   result += tmp_3442;
   std::complex<double> tmp_3443;
   std::complex<double> tmp_3444;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3445;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3445 += AbsSqr(CpconjVWmconjRhRpm(gI1,gI2))*B00(p,MRpm(gI2),
            MRh(gI1));
      }
      tmp_3444 += tmp_3445;
   }
   tmp_3443 += tmp_3444;
   result += (-4) * tmp_3443;
   std::complex<double> tmp_3446;
   std::complex<double> tmp_3447;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3448;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3448 += AbsSqr(CpconjVWmRpmRh(gI1,gI2))*B00(p,MRh(gI2),MRpm(
            gI1));
      }
      tmp_3447 += tmp_3448;
   }
   tmp_3446 += tmp_3447;
   result += (-4) * tmp_3446;
   std::complex<double> tmp_3449;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3450;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3450 += (AbsSqr(CpconjVWmbarCha1ChiPL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarCha1ChiPR(gI1,gI2)))*H0(p,MCha1(gI1),MChi(gI2));
         tmp_3450 += 4*B0(p,MCha1(gI1),MChi(gI2))*MCha1(gI1)*MChi(gI2)*Re
            (Conj(CpconjVWmbarCha1ChiPL(gI1,gI2))*CpconjVWmbarCha1ChiPR(gI1,gI2));
      }
      tmp_3449 += tmp_3450;
   }
   result += tmp_3449;
   std::complex<double> tmp_3451;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3451 += A0(MSv(gI1))*CpVWmconjVWmconjSvSv(gI1,gI1);
   }
   result += tmp_3451;
   std::complex<double> tmp_3452;
   std::complex<double> tmp_3453;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3454;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3454 += AbsSqr(CpconjVWmHpmhh(1 + gI1,gI2))*B00(p,Mhh(gI2),
            MHpm(1 + gI1));
      }
      tmp_3453 += tmp_3454;
   }
   tmp_3452 += tmp_3453;
   result += (-4) * tmp_3452;
   std::complex<double> tmp_3455;
   std::complex<double> tmp_3456;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3457;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3457 += AbsSqr(CpconjVWmconjSvSe(gI1,gI2))*B00(p,MSe(gI2),
            MSv(gI1));
      }
      tmp_3456 += tmp_3457;
   }
   tmp_3455 += tmp_3456;
   result += (-4) * tmp_3455;
   std::complex<double> tmp_3458;
   std::complex<double> tmp_3459;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3459 += A0(Mhh(gI1))*CpVWmconjVWmhhhh(gI1,gI1);
   }
   tmp_3458 += tmp_3459;
   result += (0.5) * tmp_3458;
   std::complex<double> tmp_3460;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3461;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3461 += (AbsSqr(CpconjVWmbarChiCha2PL(gI1,gI2)) + AbsSqr(
            CpconjVWmbarChiCha2PR(gI1,gI2)))*H0(p,MChi(gI1),MCha2(gI2));
         tmp_3461 += 4*B0(p,MChi(gI1),MCha2(gI2))*MCha2(gI2)*MChi(gI1)*Re
            (Conj(CpconjVWmbarChiCha2PL(gI1,gI2))*CpconjVWmbarChiCha2PR(gI1,gI2));
      }
      tmp_3460 += tmp_3461;
   }
   result += tmp_3460;
   std::complex<double> tmp_3462;
   std::complex<double> tmp_3463;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3463 += A0(MSd(gI1))*CpVWmconjVWmconjSdSd(gI1,gI1);
   }
   tmp_3462 += tmp_3463;
   result += (3) * tmp_3462;
   std::complex<double> tmp_3464;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3464 += A0(MSe(gI1))*CpVWmconjVWmconjSeSe(gI1,gI1);
   }
   result += tmp_3464;
   std::complex<double> tmp_3465;
   std::complex<double> tmp_3466;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3466 += A0(MSu(gI1))*CpVWmconjVWmconjSuSu(gI1,gI1);
   }
   tmp_3465 += tmp_3466;
   result += (3) * tmp_3465;
   std::complex<double> tmp_3467;
   std::complex<double> tmp_3468;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3469;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3469 += AbsSqr(CpconjVWmconjSuSd(gI1,gI2))*B00(p,MSd(gI2),
            MSu(gI1));
      }
      tmp_3468 += tmp_3469;
   }
   tmp_3467 += tmp_3468;
   result += (-12) * tmp_3467;
   std::complex<double> tmp_3470;
   for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
      tmp_3470 += AbsSqr(CpconjVWmVWmhh(gI2))*B0(p,MVWm,Mhh(gI2));
   }
   result += tmp_3470;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3471;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3472;
      std::complex<double> tmp_3473;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3473 += B0(p,MCha1(gI1),MSv(gI2))*Conj(CpbarFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3472 += tmp_3473;
      tmp_3471 += (MCha1(gI1)) * tmp_3472;
   }
   result += tmp_3471;
   std::complex<double> tmp_3474;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3475;
      std::complex<double> tmp_3476;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3476 += B0(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3475 += tmp_3476;
      tmp_3474 += (MFe(gI1)) * tmp_3475;
   }
   result += tmp_3474;
   std::complex<double> tmp_3477;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3478;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3478 += B0(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2)*MFe(gI2);
      }
      tmp_3477 += tmp_3478;
   }
   result += tmp_3477;
   std::complex<double> tmp_3479;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3480;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3480 += B0(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2)*MFv(gI2);
      }
      tmp_3479 += tmp_3480;
   }
   result += tmp_3479;
   std::complex<double> tmp_3481;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3482;
      std::complex<double> tmp_3483;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3483 += B0(p,MChi(gI1),MSe(gI2))*Conj(CpbarFebarChiSePL(gO2,
            gI1,gI2))*CpbarFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3482 += tmp_3483;
      tmp_3481 += (MChi(gI1)) * tmp_3482;
   }
   result += tmp_3481;
   std::complex<double> tmp_3484;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3485;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3485 += B0(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3484 += tmp_3485;
   }
   result += tmp_3484;
   std::complex<double> tmp_3486;
   std::complex<double> tmp_3487;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3487 += B0(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2)*MFe(gI2);
   }
   tmp_3486 += tmp_3487;
   result += (-4) * tmp_3486;
   std::complex<double> tmp_3488;
   std::complex<double> tmp_3489;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3489 += B0(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2)*MFv(gI2);
   }
   tmp_3488 += tmp_3489;
   result += (-4) * tmp_3488;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3490;
   std::complex<double> tmp_3491;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3492;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3492 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarFebarCha1SvPR(
            gO2,gI1,gI2))*CpbarFebarCha1SvPR(gO1,gI1,gI2);
      }
      tmp_3491 += tmp_3492;
   }
   tmp_3490 += tmp_3491;
   result += (-0.5) * tmp_3490;
   std::complex<double> tmp_3493;
   std::complex<double> tmp_3494;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3495;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3495 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPR(gO2,gI1,
            gI2))*CpbarFeFeAhPR(gO1,gI1,gI2);
      }
      tmp_3494 += tmp_3495;
   }
   tmp_3493 += tmp_3494;
   result += (-0.5) * tmp_3493;
   std::complex<double> tmp_3496;
   std::complex<double> tmp_3497;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3498;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3498 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePR(gO2,gI1,
            gI2))*CpbarFehhFePR(gO1,gI1,gI2);
      }
      tmp_3497 += tmp_3498;
   }
   tmp_3496 += tmp_3497;
   result += (-0.5) * tmp_3496;
   std::complex<double> tmp_3499;
   std::complex<double> tmp_3500;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3501;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3501 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPR(gO2,gI1
            ,gI2))*CpbarFeHpmFvPR(gO1,gI1,gI2);
      }
      tmp_3500 += tmp_3501;
   }
   tmp_3499 += tmp_3500;
   result += (-0.5) * tmp_3499;
   std::complex<double> tmp_3502;
   std::complex<double> tmp_3503;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3504;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3504 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarFebarChiSePR(gO2,
            gI1,gI2))*CpbarFebarChiSePR(gO1,gI1,gI2);
      }
      tmp_3503 += tmp_3504;
   }
   tmp_3502 += tmp_3503;
   result += (-0.5) * tmp_3502;
   std::complex<double> tmp_3505;
   std::complex<double> tmp_3506;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3507;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3507 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPR(gO2,gI1
            ,gI2))*CpbarFeSeChiPR(gO1,gI1,gI2);
      }
      tmp_3506 += tmp_3507;
   }
   tmp_3505 += tmp_3506;
   result += (-0.5) * tmp_3505;
   std::complex<double> tmp_3508;
   std::complex<double> tmp_3509;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3509 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPL(gO2,gI2))*
         CpbarFeVWmFvPL(gO1,gI2);
   }
   tmp_3508 += tmp_3509;
   result += (-1) * tmp_3508;
   std::complex<double> tmp_3510;
   std::complex<double> tmp_3511;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3511 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePL(gO2,gI2))*
         CpbarFeVZFePL(gO1,gI2);
   }
   tmp_3510 += tmp_3511;
   result += (-1) * tmp_3510;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fe_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3512;
   std::complex<double> tmp_3513;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3514;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3514 += B1(p,MCha1(gI1),MSv(gI2))*Conj(CpbarFebarCha1SvPL(
            gO2,gI1,gI2))*CpbarFebarCha1SvPL(gO1,gI1,gI2);
      }
      tmp_3513 += tmp_3514;
   }
   tmp_3512 += tmp_3513;
   result += (-0.5) * tmp_3512;
   std::complex<double> tmp_3515;
   std::complex<double> tmp_3516;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3517;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3517 += B1(p,MFe(gI1),MAh(gI2))*Conj(CpbarFeFeAhPL(gO2,gI1,
            gI2))*CpbarFeFeAhPL(gO1,gI1,gI2);
      }
      tmp_3516 += tmp_3517;
   }
   tmp_3515 += tmp_3516;
   result += (-0.5) * tmp_3515;
   std::complex<double> tmp_3518;
   std::complex<double> tmp_3519;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3520;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3520 += B1(p,MFe(gI2),Mhh(gI1))*Conj(CpbarFehhFePL(gO2,gI1,
            gI2))*CpbarFehhFePL(gO1,gI1,gI2);
      }
      tmp_3519 += tmp_3520;
   }
   tmp_3518 += tmp_3519;
   result += (-0.5) * tmp_3518;
   std::complex<double> tmp_3521;
   std::complex<double> tmp_3522;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3523;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3523 += B1(p,MFv(gI2),MHpm(gI1))*Conj(CpbarFeHpmFvPL(gO2,gI1
            ,gI2))*CpbarFeHpmFvPL(gO1,gI1,gI2);
      }
      tmp_3522 += tmp_3523;
   }
   tmp_3521 += tmp_3522;
   result += (-0.5) * tmp_3521;
   std::complex<double> tmp_3524;
   std::complex<double> tmp_3525;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3526;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3526 += B1(p,MChi(gI1),MSe(gI2))*Conj(CpbarFebarChiSePL(gO2,
            gI1,gI2))*CpbarFebarChiSePL(gO1,gI1,gI2);
      }
      tmp_3525 += tmp_3526;
   }
   tmp_3524 += tmp_3525;
   result += (-0.5) * tmp_3524;
   std::complex<double> tmp_3527;
   std::complex<double> tmp_3528;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3529;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3529 += B1(p,MChi(gI2),MSe(gI1))*Conj(CpbarFeSeChiPL(gO2,gI1
            ,gI2))*CpbarFeSeChiPL(gO1,gI1,gI2);
      }
      tmp_3528 += tmp_3529;
   }
   tmp_3527 += tmp_3528;
   result += (-0.5) * tmp_3527;
   std::complex<double> tmp_3530;
   std::complex<double> tmp_3531;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3531 += B1(p,MFv(gI2),MVWm)*Conj(CpbarFeVWmFvPR(gO2,gI2))*
         CpbarFeVWmFvPR(gO1,gI2);
   }
   tmp_3530 += tmp_3531;
   result += (-1) * tmp_3530;
   std::complex<double> tmp_3532;
   std::complex<double> tmp_3533;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3533 += B1(p,MFe(gI2),MVZ)*Conj(CpbarFeVZFePR(gO2,gI2))*
         CpbarFeVZFePR(gO1,gI2);
   }
   tmp_3532 += tmp_3533;
   result += (-1) * tmp_3532;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3534;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3535;
      std::complex<double> tmp_3536;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3536 += B0(p,MCha1(gI1),MSu(gI2))*Conj(CpbarFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3535 += tmp_3536;
      tmp_3534 += (MCha1(gI1)) * tmp_3535;
   }
   result += tmp_3534;
   std::complex<double> tmp_3537;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3538;
      std::complex<double> tmp_3539;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3539 += B0(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3538 += tmp_3539;
      tmp_3537 += (MFd(gI1)) * tmp_3538;
   }
   result += tmp_3537;
   std::complex<double> tmp_3540;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3541;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3541 += B0(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3540 += tmp_3541;
   }
   result += tmp_3540;
   std::complex<double> tmp_3542;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3543;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3543 += B0(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3542 += tmp_3543;
   }
   result += tmp_3542;
   std::complex<double> tmp_3544;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3545;
      std::complex<double> tmp_3546;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3546 += B0(p,MChi(gI1),MSd(gI2))*Conj(CpbarFdbarChiSdPL(gO2,
            gI1,gI2))*CpbarFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3545 += tmp_3546;
      tmp_3544 += (MChi(gI1)) * tmp_3545;
   }
   result += tmp_3544;
   std::complex<double> tmp_3547;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3548;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3548 += B0(p,MCha2(gI2),MSu(gI1))*Conj(CpbarFdSuCha2PL(gO2,
            gI1,gI2))*CpbarFdSuCha2PR(gO1,gI1,gI2)*MCha2(gI2);
      }
      tmp_3547 += tmp_3548;
   }
   result += tmp_3547;
   std::complex<double> tmp_3549;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3550;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3550 += B0(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3549 += tmp_3550;
   }
   result += tmp_3549;
   std::complex<double> tmp_3551;
   std::complex<double> tmp_3552;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3552 += B0(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3551 += tmp_3552;
   result += (-4) * tmp_3551;
   std::complex<double> tmp_3553;
   std::complex<double> tmp_3554;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3554 += B0(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3553 += tmp_3554;
   result += (-4) * tmp_3553;
   std::complex<double> tmp_3555;
   std::complex<double> tmp_3556;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3556 += B0(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3555 += tmp_3556;
   result += (1.3333333333333333*MGlu) * tmp_3555;
   std::complex<double> tmp_3557;
   std::complex<double> tmp_3558;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3558 += B0(p,MGlu,MSd(gI2))*Conj(CpbarFdbarGluSdPL(gO2,1,gI2))*
         CpbarFdbarGluSdPR(gO1,1,gI2);
   }
   tmp_3557 += tmp_3558;
   result += (1.3333333333333333*MGlu) * tmp_3557;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3559;
   std::complex<double> tmp_3560;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3561;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3561 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarFdbarCha1SuPR(
            gO2,gI1,gI2))*CpbarFdbarCha1SuPR(gO1,gI1,gI2);
      }
      tmp_3560 += tmp_3561;
   }
   tmp_3559 += tmp_3560;
   result += (-0.5) * tmp_3559;
   std::complex<double> tmp_3562;
   std::complex<double> tmp_3563;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3564;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3564 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPR(gO2,gI1,
            gI2))*CpbarFdFdAhPR(gO1,gI1,gI2);
      }
      tmp_3563 += tmp_3564;
   }
   tmp_3562 += tmp_3563;
   result += (-0.5) * tmp_3562;
   std::complex<double> tmp_3565;
   std::complex<double> tmp_3566;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3567;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3567 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPR(gO2,gI1,
            gI2))*CpbarFdhhFdPR(gO1,gI1,gI2);
      }
      tmp_3566 += tmp_3567;
   }
   tmp_3565 += tmp_3566;
   result += (-0.5) * tmp_3565;
   std::complex<double> tmp_3568;
   std::complex<double> tmp_3569;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3570;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3570 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPR(gO2,gI1
            ,gI2))*CpbarFdHpmFuPR(gO1,gI1,gI2);
      }
      tmp_3569 += tmp_3570;
   }
   tmp_3568 += tmp_3569;
   result += (-0.5) * tmp_3568;
   std::complex<double> tmp_3571;
   std::complex<double> tmp_3572;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3573;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3573 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarFdbarChiSdPR(gO2,
            gI1,gI2))*CpbarFdbarChiSdPR(gO1,gI1,gI2);
      }
      tmp_3572 += tmp_3573;
   }
   tmp_3571 += tmp_3572;
   result += (-0.5) * tmp_3571;
   std::complex<double> tmp_3574;
   std::complex<double> tmp_3575;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3575 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPR(gO2,gI1,1))*
         CpbarFdSdGluPR(gO1,gI1,1);
   }
   tmp_3574 += tmp_3575;
   result += (-0.6666666666666666) * tmp_3574;
   std::complex<double> tmp_3576;
   std::complex<double> tmp_3577;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3578;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3578 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarFdSuCha2PR(gO2,
            gI1,gI2))*CpbarFdSuCha2PR(gO1,gI1,gI2);
      }
      tmp_3577 += tmp_3578;
   }
   tmp_3576 += tmp_3577;
   result += (-0.5) * tmp_3576;
   std::complex<double> tmp_3579;
   std::complex<double> tmp_3580;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3581;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3581 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPR(gO2,gI1
            ,gI2))*CpbarFdSdChiPR(gO1,gI1,gI2);
      }
      tmp_3580 += tmp_3581;
   }
   tmp_3579 += tmp_3580;
   result += (-0.5) * tmp_3579;
   std::complex<double> tmp_3582;
   std::complex<double> tmp_3583;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3583 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPL(gO2,gI2))*
         CpbarFdVWmFuPL(gO1,gI2);
   }
   tmp_3582 += tmp_3583;
   result += (-1) * tmp_3582;
   std::complex<double> tmp_3584;
   std::complex<double> tmp_3585;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3585 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPL(gO2,gI2))*
         CpbarFdVZFdPL(gO1,gI2);
   }
   tmp_3584 += tmp_3585;
   result += (-1) * tmp_3584;
   std::complex<double> tmp_3586;
   std::complex<double> tmp_3587;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3587 += B1(p,MGlu,MSd(gI2))*Conj(CpbarFdbarGluSdPR(gO2,1,gI2))*
         CpbarFdbarGluSdPR(gO1,1,gI2);
   }
   tmp_3586 += tmp_3587;
   result += (-0.6666666666666666) * tmp_3586;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fd_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3588;
   std::complex<double> tmp_3589;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3590;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3590 += B1(p,MCha1(gI1),MSu(gI2))*Conj(CpbarFdbarCha1SuPL(
            gO2,gI1,gI2))*CpbarFdbarCha1SuPL(gO1,gI1,gI2);
      }
      tmp_3589 += tmp_3590;
   }
   tmp_3588 += tmp_3589;
   result += (-0.5) * tmp_3588;
   std::complex<double> tmp_3591;
   std::complex<double> tmp_3592;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3593;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3593 += B1(p,MFd(gI1),MAh(gI2))*Conj(CpbarFdFdAhPL(gO2,gI1,
            gI2))*CpbarFdFdAhPL(gO1,gI1,gI2);
      }
      tmp_3592 += tmp_3593;
   }
   tmp_3591 += tmp_3592;
   result += (-0.5) * tmp_3591;
   std::complex<double> tmp_3594;
   std::complex<double> tmp_3595;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3596;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3596 += B1(p,MFd(gI2),Mhh(gI1))*Conj(CpbarFdhhFdPL(gO2,gI1,
            gI2))*CpbarFdhhFdPL(gO1,gI1,gI2);
      }
      tmp_3595 += tmp_3596;
   }
   tmp_3594 += tmp_3595;
   result += (-0.5) * tmp_3594;
   std::complex<double> tmp_3597;
   std::complex<double> tmp_3598;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3599;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3599 += B1(p,MFu(gI2),MHpm(gI1))*Conj(CpbarFdHpmFuPL(gO2,gI1
            ,gI2))*CpbarFdHpmFuPL(gO1,gI1,gI2);
      }
      tmp_3598 += tmp_3599;
   }
   tmp_3597 += tmp_3598;
   result += (-0.5) * tmp_3597;
   std::complex<double> tmp_3600;
   std::complex<double> tmp_3601;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3602;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3602 += B1(p,MChi(gI1),MSd(gI2))*Conj(CpbarFdbarChiSdPL(gO2,
            gI1,gI2))*CpbarFdbarChiSdPL(gO1,gI1,gI2);
      }
      tmp_3601 += tmp_3602;
   }
   tmp_3600 += tmp_3601;
   result += (-0.5) * tmp_3600;
   std::complex<double> tmp_3603;
   std::complex<double> tmp_3604;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3604 += B1(p,MGlu,MSd(gI1))*Conj(CpbarFdSdGluPL(gO2,gI1,1))*
         CpbarFdSdGluPL(gO1,gI1,1);
   }
   tmp_3603 += tmp_3604;
   result += (-0.6666666666666666) * tmp_3603;
   std::complex<double> tmp_3605;
   std::complex<double> tmp_3606;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3607;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3607 += B1(p,MCha2(gI2),MSu(gI1))*Conj(CpbarFdSuCha2PL(gO2,
            gI1,gI2))*CpbarFdSuCha2PL(gO1,gI1,gI2);
      }
      tmp_3606 += tmp_3607;
   }
   tmp_3605 += tmp_3606;
   result += (-0.5) * tmp_3605;
   std::complex<double> tmp_3608;
   std::complex<double> tmp_3609;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3610;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3610 += B1(p,MChi(gI2),MSd(gI1))*Conj(CpbarFdSdChiPL(gO2,gI1
            ,gI2))*CpbarFdSdChiPL(gO1,gI1,gI2);
      }
      tmp_3609 += tmp_3610;
   }
   tmp_3608 += tmp_3609;
   result += (-0.5) * tmp_3608;
   std::complex<double> tmp_3611;
   std::complex<double> tmp_3612;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3612 += B1(p,MFu(gI2),MVWm)*Conj(CpbarFdVWmFuPR(gO2,gI2))*
         CpbarFdVWmFuPR(gO1,gI2);
   }
   tmp_3611 += tmp_3612;
   result += (-1) * tmp_3611;
   std::complex<double> tmp_3613;
   std::complex<double> tmp_3614;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3614 += B1(p,MFd(gI2),MVZ)*Conj(CpbarFdVZFdPR(gO2,gI2))*
         CpbarFdVZFdPR(gO1,gI2);
   }
   tmp_3613 += tmp_3614;
   result += (-1) * tmp_3613;
   std::complex<double> tmp_3615;
   std::complex<double> tmp_3616;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3616 += B1(p,MGlu,MSd(gI2))*Conj(CpbarFdbarGluSdPL(gO2,1,gI2))*
         CpbarFdbarGluSdPL(gO1,1,gI2);
   }
   tmp_3615 += tmp_3616;
   result += (-0.6666666666666666) * tmp_3615;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3617;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3618;
      std::complex<double> tmp_3619;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3619 += B0(p,MCha2(gI1),MSd(gI2))*Conj(CpbarFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3618 += tmp_3619;
      tmp_3617 += (MCha2(gI1)) * tmp_3618;
   }
   result += tmp_3617;
   std::complex<double> tmp_3620;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3621;
      std::complex<double> tmp_3622;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3622 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3621 += tmp_3622;
      tmp_3620 += (MFu(gI1)) * tmp_3621;
   }
   result += tmp_3620;
   std::complex<double> tmp_3623;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3624;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3624 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3623 += tmp_3624;
   }
   result += tmp_3623;
   std::complex<double> tmp_3625;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3626;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3626 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3625 += tmp_3626;
   }
   result += tmp_3625;
   std::complex<double> tmp_3627;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3628;
      std::complex<double> tmp_3629;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3629 += B0(p,MChi(gI1),MSu(gI2))*Conj(CpbarFubarChiSuPL(gO2,
            gI1,gI2))*CpbarFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3628 += tmp_3629;
      tmp_3627 += (MChi(gI1)) * tmp_3628;
   }
   result += tmp_3627;
   std::complex<double> tmp_3630;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3631;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3631 += B0(p,MCha1(gI2),MSd(gI1))*Conj(CpbarFuSdCha1PL(gO2,
            gI1,gI2))*CpbarFuSdCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3630 += tmp_3631;
   }
   result += tmp_3630;
   std::complex<double> tmp_3632;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3633;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3633 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3632 += tmp_3633;
   }
   result += tmp_3632;
   std::complex<double> tmp_3634;
   std::complex<double> tmp_3635;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3635 += B0(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3634 += tmp_3635;
   result += (-4) * tmp_3634;
   std::complex<double> tmp_3636;
   std::complex<double> tmp_3637;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3637 += B0(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3636 += tmp_3637;
   result += (-4) * tmp_3636;
   std::complex<double> tmp_3638;
   std::complex<double> tmp_3639;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3639 += B0(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3638 += tmp_3639;
   result += (-4) * tmp_3638;
   std::complex<double> tmp_3640;
   std::complex<double> tmp_3641;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3641 += B0(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3640 += tmp_3641;
   result += (1.3333333333333333*MGlu) * tmp_3640;
   std::complex<double> tmp_3642;
   std::complex<double> tmp_3643;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3643 += B0(p,MGlu,MSu(gI2))*Conj(CpbarFubarGluSuPL(gO2,1,gI2))*
         CpbarFubarGluSuPR(gO1,1,gI2);
   }
   tmp_3642 += tmp_3643;
   result += (1.3333333333333333*MGlu) * tmp_3642;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3644;
   std::complex<double> tmp_3645;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3646;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3646 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarFubarCha2SdPR(
            gO2,gI1,gI2))*CpbarFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3645 += tmp_3646;
   }
   tmp_3644 += tmp_3645;
   result += (-0.5) * tmp_3644;
   std::complex<double> tmp_3647;
   std::complex<double> tmp_3648;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3649;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3649 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPR(gO2,gI1,
            gI2))*CpbarFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3648 += tmp_3649;
   }
   tmp_3647 += tmp_3648;
   result += (-0.5) * tmp_3647;
   std::complex<double> tmp_3650;
   std::complex<double> tmp_3651;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3652;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3652 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPR(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3651 += tmp_3652;
   }
   tmp_3650 += tmp_3651;
   result += (-0.5) * tmp_3650;
   std::complex<double> tmp_3653;
   std::complex<double> tmp_3654;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3655;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3655 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPR(gO2,gI1,
            gI2))*CpbarFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3654 += tmp_3655;
   }
   tmp_3653 += tmp_3654;
   result += (-0.5) * tmp_3653;
   std::complex<double> tmp_3656;
   std::complex<double> tmp_3657;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3658;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3658 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarFubarChiSuPR(gO2,
            gI1,gI2))*CpbarFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3657 += tmp_3658;
   }
   tmp_3656 += tmp_3657;
   result += (-0.5) * tmp_3656;
   std::complex<double> tmp_3659;
   std::complex<double> tmp_3660;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3660 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPR(gO2,gI1,1))*
         CpbarFuSuGluPR(gO1,gI1,1);
   }
   tmp_3659 += tmp_3660;
   result += (-0.6666666666666666) * tmp_3659;
   std::complex<double> tmp_3661;
   std::complex<double> tmp_3662;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3663;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3663 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarFuSdCha1PR(gO2,
            gI1,gI2))*CpbarFuSdCha1PR(gO1,gI1,gI2);
      }
      tmp_3662 += tmp_3663;
   }
   tmp_3661 += tmp_3662;
   result += (-0.5) * tmp_3661;
   std::complex<double> tmp_3664;
   std::complex<double> tmp_3665;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3666;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3666 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPR(gO2,gI1
            ,gI2))*CpbarFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3665 += tmp_3666;
   }
   tmp_3664 += tmp_3665;
   result += (-0.5) * tmp_3664;
   std::complex<double> tmp_3667;
   std::complex<double> tmp_3668;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3668 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPL(gO2,gI2))*
         CpbarFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3667 += tmp_3668;
   result += (-1) * tmp_3667;
   std::complex<double> tmp_3669;
   std::complex<double> tmp_3670;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3670 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPL(gO2,gI2))*
         CpbarFuVPFuPL(gO1,gI2);
   }
   tmp_3669 += tmp_3670;
   result += (-1) * tmp_3669;
   std::complex<double> tmp_3671;
   std::complex<double> tmp_3672;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3672 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPL(gO2,gI2))*
         CpbarFuVZFuPL(gO1,gI2);
   }
   tmp_3671 += tmp_3672;
   result += (-1) * tmp_3671;
   std::complex<double> tmp_3673;
   std::complex<double> tmp_3674;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3674 += B1(p,MGlu,MSu(gI2))*Conj(CpbarFubarGluSuPR(gO2,1,gI2))*
         CpbarFubarGluSuPR(gO1,1,gI2);
   }
   tmp_3673 += tmp_3674;
   result += (-0.6666666666666666) * tmp_3673;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy_rotated(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3675;
   std::complex<double> tmp_3676;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3677;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3677 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarFubarCha2SdPL(gO1,gI1,gI2);
      }
      tmp_3676 += tmp_3677;
   }
   tmp_3675 += tmp_3676;
   result += (-0.5) * tmp_3675;
   std::complex<double> tmp_3678;
   std::complex<double> tmp_3679;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3680;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3680 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarFuFuAhPL(gO2,gI1,
            gI2))*CpbarFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3679 += tmp_3680;
   }
   tmp_3678 += tmp_3679;
   result += (-0.5) * tmp_3678;
   std::complex<double> tmp_3681;
   std::complex<double> tmp_3682;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3683;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3683 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarFuconjHpmFdPL(gO2
            ,gI1,gI2))*CpbarFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3682 += tmp_3683;
   }
   tmp_3681 += tmp_3682;
   result += (-0.5) * tmp_3681;
   std::complex<double> tmp_3684;
   std::complex<double> tmp_3685;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3686;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3686 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarFuhhFuPL(gO2,gI1,
            gI2))*CpbarFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3685 += tmp_3686;
   }
   tmp_3684 += tmp_3685;
   result += (-0.5) * tmp_3684;
   std::complex<double> tmp_3687;
   std::complex<double> tmp_3688;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3689;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3689 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarFubarChiSuPL(gO2,
            gI1,gI2))*CpbarFubarChiSuPL(gO1,gI1,gI2);
      }
      tmp_3688 += tmp_3689;
   }
   tmp_3687 += tmp_3688;
   result += (-0.5) * tmp_3687;
   std::complex<double> tmp_3690;
   std::complex<double> tmp_3691;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3691 += B1(p,MGlu,MSu(gI1))*Conj(CpbarFuSuGluPL(gO2,gI1,1))*
         CpbarFuSuGluPL(gO1,gI1,1);
   }
   tmp_3690 += tmp_3691;
   result += (-0.6666666666666666) * tmp_3690;
   std::complex<double> tmp_3692;
   std::complex<double> tmp_3693;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3694;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3694 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarFuSdCha1PL(gO2,
            gI1,gI2))*CpbarFuSdCha1PL(gO1,gI1,gI2);
      }
      tmp_3693 += tmp_3694;
   }
   tmp_3692 += tmp_3693;
   result += (-0.5) * tmp_3692;
   std::complex<double> tmp_3695;
   std::complex<double> tmp_3696;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3697;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3697 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarFuSuChiPL(gO2,gI1
            ,gI2))*CpbarFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3696 += tmp_3697;
   }
   tmp_3695 += tmp_3696;
   result += (-0.5) * tmp_3695;
   std::complex<double> tmp_3698;
   std::complex<double> tmp_3699;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3699 += B1(p,MFd(gI2),MVWm)*Conj(CpbarFuconjVWmFdPR(gO2,gI2))*
         CpbarFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3698 += tmp_3699;
   result += (-1) * tmp_3698;
   std::complex<double> tmp_3700;
   std::complex<double> tmp_3701;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3701 += B1(p,MFu(gI2),0)*Conj(CpbarFuVPFuPR(gO2,gI2))*
         CpbarFuVPFuPR(gO1,gI2);
   }
   tmp_3700 += tmp_3701;
   result += (-1) * tmp_3700;
   std::complex<double> tmp_3702;
   std::complex<double> tmp_3703;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3703 += B1(p,MFu(gI2),MVZ)*Conj(CpbarFuVZFuPR(gO2,gI2))*
         CpbarFuVZFuPR(gO1,gI2);
   }
   tmp_3702 += tmp_3703;
   result += (-1) * tmp_3702;
   std::complex<double> tmp_3704;
   std::complex<double> tmp_3705;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3705 += B1(p,MGlu,MSu(gI2))*Conj(CpbarFubarGluSuPL(gO2,1,gI2))*
         CpbarFubarGluSuPL(gO1,1,gI2);
   }
   tmp_3704 += tmp_3705;
   result += (-0.6666666666666666) * tmp_3704;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_1_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3706;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3707;
      std::complex<double> tmp_3708;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3708 += B0(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3707 += tmp_3708;
      tmp_3706 += (MCha2(gI1)) * tmp_3707;
   }
   result += tmp_3706;
   std::complex<double> tmp_3709;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3710;
      std::complex<double> tmp_3711;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3711 += B0(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3710 += tmp_3711;
      tmp_3709 += (MFu(gI1)) * tmp_3710;
   }
   result += tmp_3709;
   std::complex<double> tmp_3712;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3713;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3713 += B0(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2)*MFd(gI2);
      }
      tmp_3712 += tmp_3713;
   }
   result += tmp_3712;
   std::complex<double> tmp_3714;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3715;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3715 += B0(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2)*MFu(gI2);
      }
      tmp_3714 += tmp_3715;
   }
   result += tmp_3714;
   std::complex<double> tmp_3716;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3717;
      std::complex<double> tmp_3718;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3718 += B0(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3717 += tmp_3718;
      tmp_3716 += (MChi(gI1)) * tmp_3717;
   }
   result += tmp_3716;
   std::complex<double> tmp_3719;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3720;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3720 += B0(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2)*MCha1(gI2);
      }
      tmp_3719 += tmp_3720;
   }
   result += tmp_3719;
   std::complex<double> tmp_3721;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3722;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3722 += B0(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2)*MChi(gI2);
      }
      tmp_3721 += tmp_3722;
   }
   result += tmp_3721;
   std::complex<double> tmp_3723;
   std::complex<double> tmp_3724;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3724 += B0(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2)*MFd(gI2);
   }
   tmp_3723 += tmp_3724;
   result += (-4) * tmp_3723;
   std::complex<double> tmp_3725;
   std::complex<double> tmp_3726;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3726 += B0(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3725 += tmp_3726;
   result += (-4) * tmp_3725;
   std::complex<double> tmp_3727;
   std::complex<double> tmp_3728;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3728 += B0(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2)*MFu(gI2);
   }
   tmp_3727 += tmp_3728;
   result += (-4) * tmp_3727;
   std::complex<double> tmp_3729;
   std::complex<double> tmp_3730;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3730 += B0(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3729 += tmp_3730;
   result += (1.3333333333333333*MGlu) * tmp_3729;
   std::complex<double> tmp_3731;
   std::complex<double> tmp_3732;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3732 += B0(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,1,gI2))*
         CpbarUFubarGluSuPR(gO1,1,gI2);
   }
   tmp_3731 += tmp_3732;
   result += (1.3333333333333333*MGlu) * tmp_3731;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PR_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3733;
   std::complex<double> tmp_3734;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3735;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3735 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPR(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPR(gO1,gI1,gI2);
      }
      tmp_3734 += tmp_3735;
   }
   tmp_3733 += tmp_3734;
   result += (-0.5) * tmp_3733;
   std::complex<double> tmp_3736;
   std::complex<double> tmp_3737;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3738;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3738 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPR(gO2,gI1,
            gI2))*CpbarUFuFuAhPR(gO1,gI1,gI2);
      }
      tmp_3737 += tmp_3738;
   }
   tmp_3736 += tmp_3737;
   result += (-0.5) * tmp_3736;
   std::complex<double> tmp_3739;
   std::complex<double> tmp_3740;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3741;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3741 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPR(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPR(gO1,gI1,gI2);
      }
      tmp_3740 += tmp_3741;
   }
   tmp_3739 += tmp_3740;
   result += (-0.5) * tmp_3739;
   std::complex<double> tmp_3742;
   std::complex<double> tmp_3743;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3744;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3744 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPR(gO2,gI1,
            gI2))*CpbarUFuhhFuPR(gO1,gI1,gI2);
      }
      tmp_3743 += tmp_3744;
   }
   tmp_3742 += tmp_3743;
   result += (-0.5) * tmp_3742;
   std::complex<double> tmp_3745;
   std::complex<double> tmp_3746;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3747;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3747 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPR(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPR(gO1,gI1,gI2);
      }
      tmp_3746 += tmp_3747;
   }
   tmp_3745 += tmp_3746;
   result += (-0.5) * tmp_3745;
   std::complex<double> tmp_3748;
   std::complex<double> tmp_3749;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3749 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPR(gO2,gI1,1))*
         CpbarUFuSuGluPR(gO1,gI1,1);
   }
   tmp_3748 += tmp_3749;
   result += (-0.6666666666666666) * tmp_3748;
   std::complex<double> tmp_3750;
   std::complex<double> tmp_3751;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3752;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3752 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PR(gO2,
            gI1,gI2))*CpbarUFuSdCha1PR(gO1,gI1,gI2);
      }
      tmp_3751 += tmp_3752;
   }
   tmp_3750 += tmp_3751;
   result += (-0.5) * tmp_3750;
   std::complex<double> tmp_3753;
   std::complex<double> tmp_3754;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3755;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3755 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPR(gO2,
            gI1,gI2))*CpbarUFuSuChiPR(gO1,gI1,gI2);
      }
      tmp_3754 += tmp_3755;
   }
   tmp_3753 += tmp_3754;
   result += (-0.5) * tmp_3753;
   std::complex<double> tmp_3756;
   std::complex<double> tmp_3757;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3757 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPL(gO2,gI2))*
         CpbarUFuconjVWmFdPL(gO1,gI2);
   }
   tmp_3756 += tmp_3757;
   result += (-1) * tmp_3756;
   std::complex<double> tmp_3758;
   std::complex<double> tmp_3759;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3759 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPL(gO2,gI2))*
         CpbarUFuVPFuPL(gO1,gI2);
   }
   tmp_3758 += tmp_3759;
   result += (-1) * tmp_3758;
   std::complex<double> tmp_3760;
   std::complex<double> tmp_3761;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3761 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPL(gO2,gI2))*
         CpbarUFuVZFuPL(gO1,gI2);
   }
   tmp_3760 += tmp_3761;
   result += (-1) * tmp_3760;
   std::complex<double> tmp_3762;
   std::complex<double> tmp_3763;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3763 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPR(gO2,1,gI2))*
         CpbarUFubarGluSuPR(gO1,1,gI2);
   }
   tmp_3762 += tmp_3763;
   result += (-0.6666666666666666) * tmp_3762;

   return result * oneOver16PiSqr;

}

std::complex<double> CLASSNAME::self_energy_Fu_PL_heavy(double p , unsigned gO1, unsigned gO2) const
{
   std::complex<double> result;

   std::complex<double> tmp_3764;
   std::complex<double> tmp_3765;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      std::complex<double> tmp_3766;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3766 += B1(p,MCha2(gI1),MSd(gI2))*Conj(CpbarUFubarCha2SdPL(
            gO2,gI1,gI2))*CpbarUFubarCha2SdPL(gO1,gI1,gI2);
      }
      tmp_3765 += tmp_3766;
   }
   tmp_3764 += tmp_3765;
   result += (-0.5) * tmp_3764;
   std::complex<double> tmp_3767;
   std::complex<double> tmp_3768;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      std::complex<double> tmp_3769;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3769 += B1(p,MFu(gI1),MAh(gI2))*Conj(CpbarUFuFuAhPL(gO2,gI1,
            gI2))*CpbarUFuFuAhPL(gO1,gI1,gI2);
      }
      tmp_3768 += tmp_3769;
   }
   tmp_3767 += tmp_3768;
   result += (-0.5) * tmp_3767;
   std::complex<double> tmp_3770;
   std::complex<double> tmp_3771;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3772;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3772 += B1(p,MFd(gI2),MHpm(gI1))*Conj(CpbarUFuconjHpmFdPL(
            gO2,gI1,gI2))*CpbarUFuconjHpmFdPL(gO1,gI1,gI2);
      }
      tmp_3771 += tmp_3772;
   }
   tmp_3770 += tmp_3771;
   result += (-0.5) * tmp_3770;
   std::complex<double> tmp_3773;
   std::complex<double> tmp_3774;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3775;
      for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
         tmp_3775 += B1(p,MFu(gI2),Mhh(gI1))*Conj(CpbarUFuhhFuPL(gO2,gI1,
            gI2))*CpbarUFuhhFuPL(gO1,gI1,gI2);
      }
      tmp_3774 += tmp_3775;
   }
   tmp_3773 += tmp_3774;
   result += (-0.5) * tmp_3773;
   std::complex<double> tmp_3776;
   std::complex<double> tmp_3777;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      std::complex<double> tmp_3778;
      for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
         tmp_3778 += B1(p,MChi(gI1),MSu(gI2))*Conj(CpbarUFubarChiSuPL(gO2
            ,gI1,gI2))*CpbarUFubarChiSuPL(gO1,gI1,gI2);
      }
      tmp_3777 += tmp_3778;
   }
   tmp_3776 += tmp_3777;
   result += (-0.5) * tmp_3776;
   std::complex<double> tmp_3779;
   std::complex<double> tmp_3780;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3780 += B1(p,MGlu,MSu(gI1))*Conj(CpbarUFuSuGluPL(gO2,gI1,1))*
         CpbarUFuSuGluPL(gO1,gI1,1);
   }
   tmp_3779 += tmp_3780;
   result += (-0.6666666666666666) * tmp_3779;
   std::complex<double> tmp_3781;
   std::complex<double> tmp_3782;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3783;
      for (unsigned gI2 = 0; gI2 < 2; ++gI2) {
         tmp_3783 += B1(p,MCha1(gI2),MSd(gI1))*Conj(CpbarUFuSdCha1PL(gO2,
            gI1,gI2))*CpbarUFuSdCha1PL(gO1,gI1,gI2);
      }
      tmp_3782 += tmp_3783;
   }
   tmp_3781 += tmp_3782;
   result += (-0.5) * tmp_3781;
   std::complex<double> tmp_3784;
   std::complex<double> tmp_3785;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      std::complex<double> tmp_3786;
      for (unsigned gI2 = 0; gI2 < 4; ++gI2) {
         tmp_3786 += B1(p,MChi(gI2),MSu(gI1))*Conj(CpbarUFuSuChiPL(gO2,
            gI1,gI2))*CpbarUFuSuChiPL(gO1,gI1,gI2);
      }
      tmp_3785 += tmp_3786;
   }
   tmp_3784 += tmp_3785;
   result += (-0.5) * tmp_3784;
   std::complex<double> tmp_3787;
   std::complex<double> tmp_3788;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3788 += B1(p,MFd(gI2),MVWm)*Conj(CpbarUFuconjVWmFdPR(gO2,gI2))*
         CpbarUFuconjVWmFdPR(gO1,gI2);
   }
   tmp_3787 += tmp_3788;
   result += (-1) * tmp_3787;
   std::complex<double> tmp_3789;
   std::complex<double> tmp_3790;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3790 += B1(p,MFu(gI2),0)*Conj(CpbarUFuVPFuPR(gO2,gI2))*
         CpbarUFuVPFuPR(gO1,gI2);
   }
   tmp_3789 += tmp_3790;
   result += (-1) * tmp_3789;
   std::complex<double> tmp_3791;
   std::complex<double> tmp_3792;
   for (unsigned gI2 = 0; gI2 < 3; ++gI2) {
      tmp_3792 += B1(p,MFu(gI2),MVZ)*Conj(CpbarUFuVZFuPR(gO2,gI2))*
         CpbarUFuVZFuPR(gO1,gI2);
   }
   tmp_3791 += tmp_3792;
   result += (-1) * tmp_3791;
   std::complex<double> tmp_3793;
   std::complex<double> tmp_3794;
   for (unsigned gI2 = 0; gI2 < 6; ++gI2) {
      tmp_3794 += B1(p,MGlu,MSu(gI2))*Conj(CpbarUFubarGluSuPL(gO2,1,gI2))*
         CpbarUFubarGluSuPL(gO1,1,gI2);
   }
   tmp_3793 += tmp_3794;
   result += (-0.6666666666666666) * tmp_3793;

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
   std::complex<double> tmp_3795;
   std::complex<double> tmp_3796;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3796 += A0(MRh(gI1))*CpUhhconjRhRh(gO1,gI1,gI1);
   }
   tmp_3795 += tmp_3796;
   result += (-1) * tmp_3795;
   std::complex<double> tmp_3797;
   std::complex<double> tmp_3798;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3798 += A0(MRpm(gI1))*CpUhhconjRpmRpm(gO1,gI1,gI1);
   }
   tmp_3797 += tmp_3798;
   result += (-1) * tmp_3797;
   std::complex<double> tmp_3799;
   std::complex<double> tmp_3800;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3800 += A0(MCha1(gI1))*(CpUhhbarCha1Cha1PL(gO1,gI1,gI1) +
         CpUhhbarCha1Cha1PR(gO1,gI1,gI1))*MCha1(gI1);
   }
   tmp_3799 += tmp_3800;
   result += (2) * tmp_3799;
   std::complex<double> tmp_3801;
   std::complex<double> tmp_3802;
   for (unsigned gI1 = 0; gI1 < 2; ++gI1) {
      tmp_3802 += A0(MCha2(gI1))*(CpUhhbarCha2Cha2PL(gO1,gI1,gI1) +
         CpUhhbarCha2Cha2PR(gO1,gI1,gI1))*MCha2(gI1);
   }
   tmp_3801 += tmp_3802;
   result += (2) * tmp_3801;
   std::complex<double> tmp_3803;
   std::complex<double> tmp_3804;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3804 += A0(MSv(gI1))*CpUhhconjSvSv(gO1,gI1,gI1);
   }
   tmp_3803 += tmp_3804;
   result += (-1) * tmp_3803;
   std::complex<double> tmp_3805;
   std::complex<double> tmp_3806;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3806 += A0(MFd(gI1))*(CpUhhbarFdFdPL(gO1,gI1,gI1) + CpUhhbarFdFdPR
         (gO1,gI1,gI1))*MFd(gI1);
   }
   tmp_3805 += tmp_3806;
   result += (6) * tmp_3805;
   std::complex<double> tmp_3807;
   std::complex<double> tmp_3808;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3808 += A0(MFe(gI1))*(CpUhhbarFeFePL(gO1,gI1,gI1) + CpUhhbarFeFePR
         (gO1,gI1,gI1))*MFe(gI1);
   }
   tmp_3807 += tmp_3808;
   result += (2) * tmp_3807;
   std::complex<double> tmp_3809;
   std::complex<double> tmp_3810;
   for (unsigned gI1 = 0; gI1 < 3; ++gI1) {
      tmp_3810 += A0(MFu(gI1))*(CpUhhbarFuFuPL(gO1,gI1,gI1) + CpUhhbarFuFuPR
         (gO1,gI1,gI1))*MFu(gI1);
   }
   tmp_3809 += tmp_3810;
   result += (6) * tmp_3809;
   std::complex<double> tmp_3811;
   std::complex<double> tmp_3812;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3812 += A0(MAh(gI1))*CpUhhAhAh(gO1,gI1,gI1);
   }
   tmp_3811 += tmp_3812;
   result += (-0.5) * tmp_3811;
   std::complex<double> tmp_3813;
   std::complex<double> tmp_3814;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3814 += A0(MHpm(gI1))*CpUhhconjHpmHpm(gO1,gI1,gI1);
   }
   tmp_3813 += tmp_3814;
   result += (-1) * tmp_3813;
   std::complex<double> tmp_3815;
   std::complex<double> tmp_3816;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3816 += A0(Mhh(gI1))*CpUhhhhhh(gO1,gI1,gI1);
   }
   tmp_3815 += tmp_3816;
   result += (-0.5) * tmp_3815;
   std::complex<double> tmp_3817;
   std::complex<double> tmp_3818;
   for (unsigned gI1 = 0; gI1 < 4; ++gI1) {
      tmp_3818 += A0(MChi(gI1))*(CpUhhbarChiChiPL(gO1,gI1,gI1) +
         CpUhhbarChiChiPR(gO1,gI1,gI1))*MChi(gI1);
   }
   tmp_3817 += tmp_3818;
   result += (2) * tmp_3817;
   std::complex<double> tmp_3819;
   std::complex<double> tmp_3820;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3820 += A0(MSd(gI1))*CpUhhconjSdSd(gO1,gI1,gI1);
   }
   tmp_3819 += tmp_3820;
   result += (-3) * tmp_3819;
   std::complex<double> tmp_3821;
   std::complex<double> tmp_3822;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3822 += A0(MSe(gI1))*CpUhhconjSeSe(gO1,gI1,gI1);
   }
   tmp_3821 += tmp_3822;
   result += (-1) * tmp_3821;
   std::complex<double> tmp_3823;
   std::complex<double> tmp_3824;
   for (unsigned gI1 = 0; gI1 < 6; ++gI1) {
      tmp_3824 += A0(MSu(gI1))*CpUhhconjSuSu(gO1,gI1,gI1);
   }
   tmp_3823 += tmp_3824;
   result += (-3) * tmp_3823;

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

void CLASSNAME::calculate_MSOc_pole()
{
   if (!force_output && problems.is_tachyon(SOc))
      return;

   // diagonalization with medium precision
   const double M_tree(get_mass_matrix_SOc());
   const double p = MSOc;
   const double self_energy = Re(self_energy_SOc(p));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_tachyon(SOc);

   PHYSICAL(MSOc) = AbsSqrt(mass_sqr);
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
         problems.flag_bad_mass(MRSSM_info::Sd, eigenvalue_error >
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
         problems.flag_bad_mass(MRSSM_info::Sv, eigenvalue_error >
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
         problems.flag_bad_mass(MRSSM_info::Su, eigenvalue_error >
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
         problems.flag_bad_mass(MRSSM_info::Se, eigenvalue_error >
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
      Eigen::Matrix<double,4,4> self_energy;
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_hh());

      for (unsigned es = 0; es < 4; ++es) {
         const double p = Abs(old_Mhh(es));
         for (unsigned i1 = 0; i1 < 4; ++i1) {
            for (unsigned i2 = i1; i2 < 4; ++i2) {
               self_energy(i1,i2) = Re(self_energy_hh(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,4,4> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZH;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZH, eigenvalue_error);
            problems.flag_bad_mass(MRSSM_info::hh,
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
      Eigen::Matrix<double,4,4> self_energy;
      const Eigen::Matrix<double,4,4> M_tree(get_mass_matrix_Ah());

      for (unsigned es = 0; es < 4; ++es) {
         const double p = Abs(old_MAh(es));
         for (unsigned i1 = 0; i1 < 4; ++i1) {
            for (unsigned i2 = i1; i2 < 4; ++i2) {
               self_energy(i1,i2) = Re(self_energy_Ah(p,i1,i2
                  ));
            }
         }

         Symmetrize(self_energy);
         const Eigen::Matrix<double,4,4> M_1loop(M_tree -
            self_energy);
         Eigen::Array<double,4,1> eigen_values;
         Eigen::Matrix<double,4,4> mix_ZA;
         #ifdef CHECK_EIGENVALUE_ERROR
            double eigenvalue_error;
            fs_diagonalize_hermitian(M_1loop, eigen_values,
               mix_ZA, eigenvalue_error);
            problems.flag_bad_mass(MRSSM_info::Ah,
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

void CLASSNAME::calculate_MRh_pole()
{
   if (!force_output && problems.is_tachyon(Rh))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,2,2> self_energy;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Rh());

   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MRh(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = i1; i2 < 2; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Rh(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,2,2> M_1loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZHR;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZHR,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Rh, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZHR);
      #endif

      if (eigen_values(es) < 0.)
         problems.flag_tachyon(Rh);

      PHYSICAL(MRh(es)) = AbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZHR) = mix_ZHR;
   }
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
            problems.flag_bad_mass(MRSSM_info::Hpm,
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

void CLASSNAME::calculate_MRpm_pole()
{
   if (!force_output && problems.is_tachyon(Rpm))
      return;

   // diagonalization with medium precision
   Eigen::Matrix<double,2,2> self_energy;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Rpm());

   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MRpm(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = i1; i2 < 2; ++i2) {
            self_energy(i1,i2) = Re(self_energy_Rpm(p,i1,i2));
         }
      }

      Symmetrize(self_energy);
      const Eigen::Matrix<double,2,2> M_1loop(M_tree - self_energy);
      Eigen::Array<double,2,1> eigen_values;
      Eigen::Matrix<double,2,2> mix_ZRP;
      #ifdef CHECK_EIGENVALUE_ERROR
         double eigenvalue_error;
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZRP,
            eigenvalue_error);
         problems.flag_bad_mass(MRSSM_info::Rpm, eigenvalue_error >
            precision * Abs(eigen_values(0)));
      #else
         fs_diagonalize_hermitian(M_1loop, eigen_values, mix_ZRP);
      #endif

      if (eigen_values(es) < 0.)
         problems.flag_tachyon(Rpm);

      PHYSICAL(MRpm(es)) = AbsSqrt(eigen_values(es));
      if (es == 0)
         PHYSICAL(ZRP) = mix_ZRP;
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
      const Eigen::Matrix<double,4,4> M_1loop(M_tree + delta_M);
      Eigen::Array<double,4,1> eigen_values;
      decltype(ZN1) mix_ZN1;
      decltype(ZN2) mix_ZN2;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZN1, mix_ZN2, eigenvalue_error
         );
      problems.flag_bad_mass(MRSSM_info::Chi, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_ZN1, mix_ZN2);
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
   Eigen::Matrix<double,2,2> self_energy_1;
   Eigen::Matrix<double,2,2> self_energy_PL;
   Eigen::Matrix<double,2,2> self_energy_PR;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha1());
   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MCha1(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = 0; i2 < 2; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Cha1_1(p,i1,
               i2));
            self_energy_PL(i1,i2) = Re(self_energy_Cha1_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Cha1_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_1loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM1) mix_UM1;
      decltype(UP1) mix_UP1;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_UM1, mix_UP1, eigenvalue_error
         );
      problems.flag_bad_mass(MRSSM_info::Cha1, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_UM1, mix_UP1);
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
   Eigen::Matrix<double,2,2> self_energy_1;
   Eigen::Matrix<double,2,2> self_energy_PL;
   Eigen::Matrix<double,2,2> self_energy_PR;
   const Eigen::Matrix<double,2,2> M_tree(get_mass_matrix_Cha2());
   for (unsigned es = 0; es < 2; ++es) {
      const double p = Abs(MCha2(es));
      for (unsigned i1 = 0; i1 < 2; ++i1) {
         for (unsigned i2 = 0; i2 < 2; ++i2) {
            self_energy_1(i1,i2)  = Re(self_energy_Cha2_1(p,i1,
               i2));
            self_energy_PL(i1,i2) = Re(self_energy_Cha2_PL(p,i1,
               i2));
            self_energy_PR(i1,i2) = Re(self_energy_Cha2_PR(p,i1,
               i2));
         }
      }
      const Eigen::Matrix<double,2,2> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,2,2> M_1loop(M_tree + delta_M);
      Eigen::Array<double,2,1> eigen_values;
      decltype(UM2) mix_UM2;
      decltype(UP2) mix_UP2;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_UM2, mix_UP2, eigenvalue_error
         );
      problems.flag_bad_mass(MRSSM_info::Cha2, eigenvalue_error >
         precision * Abs(eigen_values(0)));
   #else
      fs_svd(M_1loop, eigen_values, mix_UM2, mix_UP2);
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
      problems.flag_bad_mass(MRSSM_info::Fe, eigenvalue_error >
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
      problems.flag_bad_mass(MRSSM_info::Fd, eigenvalue_error >
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

   double qcd_2l = 0.;

   if (add_2loop_corrections) {
      const double currentScale = get_scale();
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
            self_energy_1(i1,i2)  = Re(self_energy_Fu_1(p,i1,i2)
               );
            self_energy_PL(i1,i2) = Re(self_energy_Fu_PL(p,i1,i2
               ));
            self_energy_PR(i1,i2) = Re(self_energy_Fu_PR(p,i1,i2
               ));
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * qcd_2l;
      const Eigen::Matrix<double,3,3> M_1loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(ZUL) mix_ZUL;
      decltype(ZUR) mix_ZUR;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_1loop, eigen_values, mix_ZUL, mix_ZUR, eigenvalue_error
         );
      problems.flag_bad_mass(MRSSM_info::Fu, eigenvalue_error >
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


std::ostream& operator<<(std::ostream& ostr, const MRSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
